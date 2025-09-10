# =============================================================================
# НАСТРОЙКА СРЕДЫ И ЗАГРУЗКА ПАКЕТОВ
# =============================================================================

# Установка рабочей директории - указываем путь к папке с данными
setwd("C:/LBM/")

# Тихая загрузка пакетов (без вывода сообщений)
suppressPackageStartupMessages({
  library(tidyverse)    # Набор пакетов для обработки и визуализации данных
  library(lubridate)    # Работа с датами
  library(ggridges)     # Риджлайн-графики (горные хребты)
  library(TropFishR)    # Пакет для анализа данных рыболовства
})

# =============================================================================
# 1. ЗАГРУЗКА И ПРЕДВАРИТЕЛЬНАЯ ОБРАБОТКА ДАННЫХ
# =============================================================================

# Загрузка данных из CSV-файлов с разделителем ";"
SURVEY  <- read.csv("SURVEYDATA.csv",  sep = ";", stringsAsFactors = FALSE)  # Данные съемки
FISHERY <- read.csv("FISHERYDATA.csv", sep = ";", stringsAsFactors = FALSE)  # Данные промысла

# ОБРАБОТКА ДАННЫХ СЪЕМКИ:
# 1. Добавляем столбец TYPE для идентификации типа данных
# 2. Стандартизируем обозначения пола (M/F)
# 3. Создаем дату (середина месяца)
# 4. Преобразуем размер карапакса в числовой формат
# 5. Фильтруем некорректные значения (отсутствующие, отрицательные, слишком большие)
survey <- SURVEY %>%
  mutate(
    TYPE     = "Съёмка",
    SEX      = if_else(SEX %in% c("M","male","m"), "M",
                       if_else(SEX %in% c("F","female","f"), "F", NA_character_)),
    date     = ymd(paste(YEAR, MONTH, 15, sep = "-")),
    CARAPACE = as.numeric(CARAPACE)
  ) %>%
  filter(!is.na(CARAPACE), CARAPACE > 0, CARAPACE < 220)

# Аналогичная обработка данных промысла
fishery <- FISHERY %>%
  mutate(
    TYPE     = "Промысел",
    SEX      = if_else(SEX %in% c("M","male","m"), "M",
                       if_else(SEX %in% c("F","female","f"), "F", NA_character_)),
    date     = ymd(paste(YEAR, MONTH, 15, sep = "-")),
    CARAPACE = as.numeric(CARAPACE)
  ) %>%
  filter(!is.na(CARAPACE), CARAPACE > 0, CARAPACE < 240)

# =============================================================================
# 2. ВИЗУАЛИЗАЦИЯ РАЗМЕРНЫХ РАСПРЕДЕЛЕНИЙ (РИДЖЛАЙН-ГРАФИКИ)
# =============================================================================

# Объединяем данные съемки и промысла и фильтруем только самцов
df <- bind_rows(survey, fishery) %>%
  filter(!is.na(YEAR), !is.na(TYPE), SEX == "M") %>%
  mutate(
    YEAR   = as.integer(YEAR),
    YEAR_F = factor(YEAR, levels = sort(unique(YEAR)))  # Преобразуем год в фактор
  )

# Определяем границы для оси X (округленные до десятков)
xmin <- floor(min(df$CARAPACE, na.rm = TRUE) / 10) * 10
xmax <- ceiling(max(df$CARAPACE, na.rm = TRUE) / 10) * 10

# Создаем риджлайн-график (горный хребет)
# Риджлайны показывают распределение размеров по годам
p <- ggplot(df, aes(x = CARAPACE, y = YEAR_F, fill = after_stat(x))) +
  stat_density_ridges(
    geom = "density_ridges_gradient",
    calc_ecdf = FALSE,
    rel_min_height = 0.001,
    bandwidth = 3  # Параметр сглаживания
  ) +
  scale_fill_viridis_c(name = "Размер, мм", option = "C") +  # Цветовая схема
  facet_wrap(~ TYPE, ncol = 1, scales = "free_y") +  # Разделяем по типам данных
  coord_cartesian(xlim = c(xmin, xmax)) +  # Устанавливаем границы оси X
  labs(x = "Ширина карапакса, мм", y = "Год",
       title = "Размерные распределения по годам: съёмка и промысел (самцы)") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "right", panel.grid.minor = element_blank())

p  # Отображаем график

# Сохранение графика (раскомментируйте при необходимости)
# ggsave("size_ridges_by_year_type_M.png", p, width = 10, height = 12, dpi = 200)

# =============================================================================
# 3. ОЦЕНКА ПАРАМЕТРОВ РОСТА МЕТОДОМ ELEFAN_SA
# =============================================================================

# Фильтруем данные только по самцам
survey_m <- survey %>% filter(SEX == "M")
fish_m   <- fishery %>% filter(SEX == "M")

# Не ограничиваем анализ определенными годами для устойчивости результатов
survey_m <- survey_m %>% filter(!is.na(YEAR))  # Все доступные данные

 #years_use <- 2010:2019  # 10 лет вместо неограниченного периода
 #survey_m <- survey_m %>% filter(YEAR %in% years_use)


# Создаем объект LFQ (Length-Frequency Data) для анализа
binw <- 2  # Ширина класса размеров (2 мм)
lfq <- TropFishR::lfqCreate(
  data     = survey_m %>% mutate(date = as.Date(date)),
  Lname    = "CARAPACE",  # Столбец с размерами
  Dname    = "date",      # Столбец с датами
  bin_size = binw         # Ширина класса
)

# Инициализация параметра L∞ (максимальный теоретический размер)
# Основана на 99-м процентиле размеров промысла с поправочным коэффициентом
init_Linf2 <- as.numeric(quantile(fish_m$CARAPACE, 0.99, na.rm = TRUE)) / 0.95
init_Linf2 <- min(max(init_Linf2, 155), 165)  # Ограничиваем разумными пределами

# Оценка параметров роста методом имитации отжига (ELEFAN_SA)
set.seed(1)  # Для воспроизводимости результатов
res_sa2 <- TropFishR::ELEFAN_SA(
  lfq,
  seasonalised = TRUE,  # Используем сезонную модель
  # Начальные значения параметров:
  init_par = list(Linf = init_Linf2, K = 0.2, t_anchor = 0.1, C = 0.2, ts = 0.5),
  # Нижние границы параметров:
  low_par  = list(Linf = 155,        K = 0.20, t_anchor = 0.0, C = 0.0, ts = 0.0),
  # Верхние границы параметров:
  up_par   = list(Linf = 165,        K = 0.35, t_anchor = 1.0, C = 0.4, ts = 1.0),
  SA_time  = 60  # Время работы алгоритма (секунды)
)


# Извлекаем оцененные параметры
pars2 <- res_sa2$par

# Рассчитываем дополнительные показатели:
# φ' (фи-prime) - индекс производительности роста
phi_prime <- log10(pars2$K) + 2*log10(pars2$Linf)
# Время достижения 95% от L∞
t95 <- function(K) -log(1 - 0.95)/K

# Вывод результатов
cat(sprintf("ELEFAN_SA (ограничения): L∞=%.1f мм, K=%.3f, C=%.2f, ts=%.2f, phi' = %.3f, t95 = %.1f лет\n",
            pars2$Linf, pars2$K, pars2$C, pars2$ts, phi_prime, t95(pars2$K)))

# =============================================================================
# 4. СРАВНЕНИЕ РАСПРЕДЕЛЕНИЙ РАЗМЕРОВ СЪЕМКИ И ПРОМЫСЛА
# =============================================================================

# График плотности распределений размеров
p_dens <- ggplot() +
  # Плотность распределения для данных съемки
  geom_density(data = survey_m, aes(CARAPACE, colour = "Съёмка"), bw = 3) +
  # Плотность распределения для данных промысла
  geom_density(data = fish_m,   aes(CARAPACE, colour = "Промысел"), bw = 3) +
  # Вертикальная линия, показывающая оценку L∞
  geom_vline(xintercept = pars2$Linf, linetype = 2) +
  # Цветовая схема
  scale_colour_manual(values = c("Съёмка" = "steelblue", "Промысел" = "tomato")) +
  labs(x = "Ширина карапакса, мм", y = "Плотность", colour = NULL,
       title = "Распределения по размерам; штрих — оценка L∞ (ELEFAN)") +
  theme_minimal(base_size = 12)

p_dens

# Сохранение графика
ggsave("densities_survey_fishery_Linf.png", p_dens, width = 8, height = 5, dpi = 200)

# =============================================================================
# 4.1 Визуализация "реструктурированных" распределений по длине (restructured counts)
# с наложенной оптимальной кривой роста, подобранной методом ELEFAN_SA.
# =============================================================================

plot(res_sa2, Fname = "rcounts", image.col = colorRampPalette(c("red","grey100","green"))(21))




# =============================================================================
# 5. ПОСТРОЕНИЕ КРИВОЙ РОСТА С ДОВЕРИТЕЛЬНЫМИ ИНТЕРВАЛАЛАМИ
# =============================================================================

# Функция сезонного уравнения роста Берталанфи
# t - возраст, Linf - асимптотическая длина, K - коэффициент роста,
# t0 - теоретический возраст при нулевой длине, C - амплитуда сезонных колебаний,
# ts - параметр сдвига сезонных колебаний
seasonal_vbgf <- function(t, Linf, K, t0, C, ts) {
  term <- -K*(t - t0) + (C*K/(2*pi)) * sin(2*pi*(t - ts)) - (C*K/(2*pi)) * sin(2*pi*(t0 - ts))
  Lt <- Linf * (1 - exp(term))
  return(Lt)
}

# Функция для предсказания роста с учетом неопределенности параметров
predict_growth <- function(ages, params, n_sim = 1000, alpha = 0.05) {
  # Извлекаем точечные оценки параметров
  Linf_est <- params$Linf
  K_est <- params$K
  t0_est <- params$t_anchor
  C_est <- params$C
  ts_est <- params$ts
  
  # Матрица для хранения симуляций
  simulations <- matrix(NA, nrow = n_sim, ncol = length(ages))
  
  # Симуляция возможных значений параметров на основе их неопределенности
  set.seed(123)
  Linf_sim <- rnorm(n_sim, mean = Linf_est, sd = Linf_est * 0.05)  # CV 5%
  K_sim <- rnorm(n_sim, mean = K_est, sd = K_est * 0.1)           # CV 10%
  t0_sim <- rnorm(n_sim, mean = t0_est, sd = 0.1)                 # SD 0.1
  C_sim <- rnorm(n_sim, mean = C_est, sd = 0.05)                  # SD 0.05
  ts_sim <- rnorm(n_sim, mean = ts_est, sd = 0.05)                # SD 0.05
  
  # Ограничение параметров физиологически реалистичными значениями
  Linf_sim <- pmax(pmin(Linf_sim, 180), 140)
  K_sim <- pmax(pmin(K_sim, 0.5), 0.05)
  C_sim <- pmax(pmin(C_sim, 1), 0)
  ts_sim <- pmax(pmin(ts_sim, 1), 0)
  
  # Расчет кривых роста для каждой симуляции
  for (i in 1:n_sim) {
    simulations[i, ] <- seasonal_vbgf(ages, Linf_sim[i], K_sim[i], 
                                      t0_sim[i], C_sim[i], ts_sim[i])
  }
  
  # Расчет доверительных интервалов
  mean_pred <- apply(simulations, 2, mean)
  lower_ci <- apply(simulations, 2, quantile, probs = alpha/2)
  upper_ci <- apply(simulations, 2, quantile, probs = 1 - alpha/2)
  
  return(data.frame(
    Age = ages,
    Length = mean_pred,
    Lower = lower_ci,
    Upper = upper_ci
  ))
}

# Создаем данные для графика с доверительными интервалами
ages <- seq(0, 20, by = 0.1)  # Возраст от 0 до 20 лет с шагом 0.1
growth_df <- predict_growth(ages, pars2, n_sim = 1000, alpha = 0.05)  # 95% ДИ

# Построение графика кривой роста
p_growth <- ggplot(growth_df, aes(x = Age, y = Length)) +
  # Доверительные интервалы (заполненная область)
  geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.2, fill = "blue") +
  # Средняя кривая роста
  geom_line(color = "blue", linewidth = 1) +
  # Горизонтальная линия, показывающая L∞
  geom_hline(yintercept = pars2$Linf, linetype = "dashed", color = "red") +
  labs(x = "Возраст (лет)", y = "Ширина карапакса (мм)",
       title = "Кривая роста (сезонный ВБГФ) с доверительными интервалами",
       subtitle = sprintf("L∞=%.1f мм, K=%.3f, C=%.2f, ts=%.2f, t95≈%.1f лет",
                          pars2$Linf, pars2$K, pars2$C, pars2$ts, t95(pars2$K))) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5))

print(p_growth)
# Сохранение графика
ggsave("seasonal_VBGF_with_CI.png", p_growth, width = 8, height = 5, dpi = 200)

# =============================================================================
# 6. СОХРАНЕНИЕ РЕЗУЛЬТАТОВ
# =============================================================================

# Сохранение оцененных параметров роста в CSV-файл
write.csv(tibble(
  Linf = pars2$Linf, 
  K = pars2$K, 
  C = pars2$C, 
  ts = pars2$ts,
  t_anchor = pars2$t_anchor, 
  phi_prime = phi_prime, 
  t95_years = t95(pars2$K)
), "ELEFAN_params_constrained.csv", row.names = FALSE)

# Сообщение о завершении анализа
cat("Готово. Файлы: size_ridges_by_year_type_M.png, densities_survey_fishery_Linf.png, seasonal_VBGF_to15y.png, ELEFAN_params_constrained.csv\n")



