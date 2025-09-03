# ========================================================================================================================
# ПРАКТИЧЕСКОЕ ЗАНЯТИЕ: МОДЕЛИ ИСТОЩЕНИЯ ЛЕСЛИ И ДЕЛУРИ
# Курс: "Оценка водных биоресурсов в среде R (для начинающих)"
# Автор: Баканев С. В. Дата: 03.09.2025
# ========================================================================================================================

# 1. ПОДГОТОВКА СРЕДЫ ====================================================================================================

# Установка и подключение необходимых пакетов
# FSA - Fisheries Stock Analysis для моделей истощения
# ggplot2 - для продвинутой визуализации
# dplyr - для манипуляций с данными
# tidyverse - современный подход к обработке данных
if (!require("FSA")) install.packages("FSA")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("dplyr")) install.packages("dplyr")
if (!require("broom")) install.packages("broom")
if (!require("tidyverse")) install.packages("tidyverse")

library(FSA)
library(ggplot2)
library(dplyr)
library(broom)
library(tidyverse)

# Установка рабочей директории - укажите путь к папке с данными
setwd("C:/DEPLETION/")

# 2. ЗАГРУЗКА И ПЕРВИЧНЫЙ АНАЛИЗ ДАННЫХ ==================================================================================

# Чтение данных из CSV-файла
# header = TRUE - первая строка содержит названия колонок
# sep = ";" - разделитель точка с запятой (common для European CSV)
LESLIDATA <- read.csv("DATAdep.csv", header = TRUE, sep = ";")

# Проверка структуры данных
# Функция str() показывает:
# - тип объекта (data.frame)
# - количество наблюдений и переменных
# - тип каждой переменной
str(LESLIDATA)

# 3. БАЗОВЫЙ АНАЛИЗ МОДЕЛИ ЛЕСЛИ ДЛЯ ОДНОГО ПЕРИОДА ======================================================================

# Фильтрация данных для 2008 года (2007 < YEAR < 2009)
DATA <- LESLIDATA[LESLIDATA$YEAR > 2007 & LESLIDATA$YEAR < 2009, ]

# Построение модели истощения по методу Лесли
# Модель Лесли: CPUE ~ кумулятивное усилие
# Использует линейную регрессию для оценки начальной численности
lesli <- depletion(DATA$CATCH, DATA$EFFORT, method = "Leslie")

# Визуализация модели
# График показывает зависимость CPUE от кумулятивного усилия
plot(lesli)

# Доверительные интервалы для параметров модели
confint(lesli)

# Сводная информация по модели
# Включает оценки параметров и статистику качества拟合
summary(lesli)

# 4. РАСШИРЕННЫЙ АНАЛИЗ ПО ГОДАМ (2007-2018) =============================================================================

# Создаем функцию для расчета кумулятивных показателей
# Кумулятивные показатели необходимы для построения моделей истощения
calculate_cumulative <- function(data) {
  data %>%
    arrange(WEEK) %>%  # Сортировка по неделям
    mutate(
      cumulative_effort = cumsum(EFFORT),  # Накопленное усилие
      cumulative_catch = cumsum(CATCH)     # Накопленный улов
    )
}

# Инициализация списка для хранения результатов
leslie_models_list <- list()

# Анализ для каждого года в диапазоне 2007-2018
for (year in 2007:2018) {
  # Фильтрация данных по году
  year_data <- LESLIDATA %>% 
    filter(YEAR == year) %>%
    na.omit()  # Удаление пропущенных значений
  
  # Проверка достаточности данных (минимум 3 наблюдения)
  if (nrow(year_data) < 3) {
    message(paste("Недостаточно данных для анализа в", year))
    next  # Переход к следующему году
  }
  
  # Расчет кумулятивных показателей
  year_data <- calculate_cumulative(year_data)
  
  # Добавление CPUE (улов на единицу усилия)
  year_data$CPUE <- year_data$CATCH / year_data$EFFORT
  
  # Построение модели Лесли через линейную регрессию
  leslie_model <- try(lm(CPUE ~ cumulative_effort, data = year_data), silent = TRUE)
  
  if (inherits(leslie_model, "try-error")) {
    # Обработка ошибок моделирования
    year_data$leslie_predicted <- NA
    year_data$leslie_lwr <- NA
    year_data$leslie_upr <- NA
    message(paste("Ошибка в модели Лесли для", year))
  } else {
    # Получение предсказаний с доверительными интервалами
    predictions <- predict(leslie_model, interval = "confidence", level = 0.95)
    year_data$leslie_predicted <- predictions[, "fit"]
    year_data$leslie_lwr <- predictions[, "lwr"]
    year_data$leslie_upr <- predictions[, "upr"]
    
    # Расчет начальной биомассы (No)
    # No = -a/b, где a - интерсепт, b - коэффициент кумулятивного усилия
    a <- coef(leslie_model)[1]
    b <- coef(leslie_model)[2]
    No <- -a / b
    year_data$No <- No
  }
  
  # Сохранение результатов для года
  leslie_models_list[[as.character(year)]] <- year_data
}

# Объединение данных всех лет
all_years_leslie <- bind_rows(leslie_models_list, .id = "Year")

# Преобразование Year в фактор с сохранением порядка
all_years_leslie$Year <- factor(all_years_leslie$Year, levels = as.character(2007:2018))

# 5. ВИЗУАЛИЗАЦИЯ РЕЗУЛЬТАТОВ МОДЕЛИ ЛЕСЛИ ==============================================================================

# Построение фасетного графика для всех лет
leslie_facet_plot <- ggplot(all_years_leslie, aes(x = cumulative_effort)) +
  geom_point(aes(y = CPUE), size = 2, color = "darkblue", alpha = 0.7) +
  geom_ribbon(aes(ymin = leslie_lwr, ymax = leslie_upr), 
              fill = "red", alpha = 0.2) +  # Доверительный интервал
  geom_line(aes(y = leslie_predicted), color = "red", linewidth = 1) +
  facet_wrap(~ Year, scales = "free", ncol = 3) +  # Свободные масштабы для каждого года
  labs(
    title = "Модель Лесли: зависимость CPUE от кумулятивного усилия",
    subtitle = "С доверительными интервалами (95%)",
    x = "Кумулятивное усилие лова",
    y = "CPUE (улов на единицу усилия)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    strip.text = element_text(face = "bold", size = 10)
  )

# Вывод графика
print(leslie_facet_plot)

# Сохранение графика (раскомментируйте для использования)
# ggsave("Leslie_Model_Facets_2007_2018_with_CI.png", plot = leslie_facet_plot, width = 14, height = 10, dpi = 300)

# 6. СРАВНИТЕЛЬНЫЙ АНАЛИЗ МОДЕЛЕЙ ЛЕСЛИ И ДЕЛУРИ ========================================================================

# Инициализация списка для хранения результатов
results_list <- list()

# Анализ для каждого года
for (year in 2007:2018) {
  year_data <- LESLIDATA %>% 
    filter(YEAR == year) %>%
    na.omit()
  
  if (nrow(year_data) < 3) next
  
  # Модель Лесли через FSA
  leslie_model <- try(depletion(year_data$CATCH, year_data$EFFORT, method = "Leslie"), silent = TRUE)
  if (inherits(leslie_model, "try-error")) {
    leslie_no <- leslie_lci <- leslie_uci <- NA
  } else {
    leslie_ci <- confint(leslie_model)
    leslie_no <- coef(leslie_model)["No"]
    leslie_lci <- leslie_ci["No", "95% LCI"]
    leslie_uci <- leslie_ci["No", "95% UCI"]
  }
  
  # Модель Делури через FSA
  delury_model <- try(depletion(year_data$CATCH, year_data$EFFORT, method = "Delury"), silent = TRUE)
  if (inherits(delury_model, "try-error")) {
    delury_no <- delury_lci <- delury_uci <- NA
  } else {
    delury_ci <- confint(delury_model)
    delury_no <- coef(delury_model)["No"]
    delury_lci <- delury_ci["No", "95% LCI"]
    delury_uci <- delury_ci["No", "95% UCI"]
  }
  
  # Сохранение результатов
  results_list[[as.character(year)]] <- data.frame(
    Year = year,
    Model = c("Лесли", "Делури"),
    Initial_Biomass = c(leslie_no, delury_no),
    LCI = c(leslie_lci, delury_lci),
    UCI = c(leslie_uci, delury_uci)
  )
}

# Преобразование списка в dataframe
results_df <- bind_rows(results_list)

# 7. ВИЗУАЛИЗАЦИЯ СРАВНЕНИЯ МОДЕЛЕЙ =====================================================================================

# График динамики начальной биомассы
biomass_plot <- ggplot(results_df, aes(x = Year, y = Initial_Biomass/1000, color = Model)) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = LCI/1000, ymax = UCI/1000), width = 0.2) +
  geom_line(aes(group = Model), linetype = "dashed") +
  labs(
    title = "Динамика начальной биомассы (2007-2018)",
    subtitle = "Модели Лесли и Делури с 95% доверительными интервалами",
    x = "Год",
    y = "Начальная биомасса, тыс. т",
    color = "Модель"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "top"
  ) +
  scale_color_manual(values = c("Лесли" = "blue", "Делури" = "red"))

print(biomass_plot)

# 8. ДЕТАЛИЗИРОВАННЫЙ АНАЛИЗ ПАРАМЕТРОВ =================================================================================

# Создание таблиц с параметрами моделей
leslie_all_years <- data.frame()
delury_all_years <- data.frame()

for (year in 2007:2018) {
  year_data <- LESLIDATA %>% filter(YEAR == year) %>% na.omit()
  if (nrow(year_data) < 3) next
  
  # Анализ для модели Лесли
  tryCatch({
    leslie_model <- depletion(year_data$CATCH, year_data$EFFORT, method = "Leslie")
    leslie_ci <- confint(leslie_model)
    leslie_all_years <- rbind(leslie_all_years, data.frame(
      Модель = "Лесли",
      Год = year,
      B0 = round(leslie_model$est["No", "Estimate"], 2),
      B0_LCI = round(leslie_ci["No", "95% LCI"], 2),
      B0_UCI = round(leslie_ci["No", "95% UCI"], 2),
      q = round(leslie_model$est["q", "Estimate"], 6),
      q_LCI = round(leslie_ci["q", "95% LCI"], 6),
      q_UCI = round(leslie_ci["q", "95% UCI"], 6),
      R2 = round(summary(leslie_model$lm)$r.squared, 4)
    ))
  }, error = function(e) {
    message(paste("Ошибка в модели Лесли для", year, ":", e$message))
  })
  
  # Анализ для модели Делури
  tryCatch({
    delury_model <- depletion(year_data$CATCH, year_data$EFFORT, method = "DeLury")
    delury_ci <- confint(delury_model)
    delury_all_years <- rbind(delury_all_years, data.frame(
      Модель = "Делури",
      Год = year,
      B0 = round(delury_model$est["No", "Estimate"], 2),
      B0_LCI = round(delury_ci["No", "95% LCI"], 2),
      B0_UCI = round(delury_ci["No", "95% UCI"], 2),
      q = round(delury_model$est["q", "Estimate"], 6),
      q_LCI = round(delury_ci["q", "95% LCI"], 6),
      q_UCI = round(delury_ci["q", "95% UCI"], 6),
      R2 = round(summary(delury_model$lm)$r.squared, 4)
    ))
  }, error = function(e) {
    message(paste("Ошибка в модели Делури для", year, ":", e$message))
  })
}

# Вывод результатов
print("Таблица параметров модели Лесли по годам:")
print(leslie_all_years)

print("Таблица параметров модели Делури по годам:")
print(delury_all_years)

# Сохранение результатов
write.csv(leslie_all_years, "Leslie_parameters_all_years.csv", row.names = FALSE)
write.csv(delury_all_years, "Delury_parameters_all_years.csv", row.names = FALSE)

# ========================================================================================================================
# ИНТЕРПРЕТАЦИЯ РЕЗУЛЬТАТОВ:
# 1. Модель Лесли: CPUE = a + b * (кумулятивное усилие)
# 2. Модель Делури: ln(CPUE) = a + b * (кумулятивное усилие)
# Параметры:
# - B0 (No): начальная биомасса/численность
# - q: коэффициент уловистости
# - R2: показатель качества модели (0-1)
# Доверительные интервалы показывают точность оценок
# ========================================================================================================================