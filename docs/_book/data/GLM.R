# ========================================================================================================================
# ПРАКТИЧЕСКОЕ ЗАНЯТИЕ: СТАНДАРТИЗАЦИЯ CPUE С ИСПОЛЬЗОВАНИЕМ GLM, GAM И GAMM
# Курс: "Оценка водных биоресурсов в среде R (для начинающих)"
# Автор: Баканев С.В. 
# Дата: 20.08.2025
# 
# Структура:
# 1) Загрузка пакетов и настройка среды
# 2) Загрузка и предварительная обработка данных
# 3) Вспомогательные функции для расчета индексов
# 4) Моделирование GLM (Gamma с лог-ссылкой)
# 5) Моделирование GAM (обобщенная аддитивная модель)
# 6) Моделирование GAMM (смешанная модель со случайными эффектами)
# 7) Сравнение моделей и финальная визуализация результатов
# ========================================================================================================================


# ==============================================================================
# БЛОК 1: ЗАГРУЗКА ПАКЕТОВ И НАСТРОЙКА СРЕДЫ
# ==============================================================================

# Отключаем вспомогательные сообщения при загрузке пакетов
suppressPackageStartupMessages({
  library(tidyverse)   # Основные пакеты для обработки данных и визуализации
  library(readxl)      # Чтение данных из Excel-файлов
  library(mgcv)        # Обобщенные аддитивные модели (GAM)
  library(gamm4)       # GAM со смешанными эффектами
  library(emmeans)     # Расчет маргинальных средних и контрастов
  library(broom)       # Преобразование результатов моделей в таблицы
  library(broom.mixed) # Поддержка смешанных моделей для broom
  library(DHARMa)      # Диагностика остатков обобщенных моделей
  library(knitr)       # Форматирование таблиц для отчетов
})

# Установка рабочей директории
setwd("C:/GLM/")

# Фиксируем случайное зерно для воспроизводимости результатов
set.seed(42)

# ==============================================================================
# БЛОК 2: ЗАГРУЗКА И ПРЕДОБРАБОТКА ДАННЫХ
# ==============================================================================

# Определяем путь к файлу с данными
DATA_PATH <- "C:/GLM/data/KARTOGRAPHIC.xlsx"

# Чтение данных из листа "FISHERY" и фильтрация осенних месяцев
DATA <- read_excel(DATA_PATH, sheet = "FISHERY") %>%
  as_tibble() %>%  # Преобразуем в современный формат таблицы
  filter(MONTH > 8 & MONTH < 12)  # Сентябрь-ноябрь (осенний сезон)

# Преобразование типов переменных и обработка пропусков
DATA <- DATA %>%
  mutate(
    YEAR = as.factor(YEAR),           # Год как категориальная переменная
    MONTH = as.factor(MONTH),         # Месяц как фактор
    CALL = as.factor(CALL),           # Идентификатор судна
    REGION = as.factor(REGION),       # Рыбохозяйственный район
    VESSELNUMBER = as.factor(VESSELNUMBER),  # Номер судна
    CPUE = as.numeric(CPUE)           # Целевой показатель - улов на усилие
  ) %>%
  filter(!is.na(CPUE))  # Удаление строк с пропусками в CPUE

# Обработка нулевых значений CPUE для Gamma-моделей
if (any(DATA$CPUE <= 0, na.rm = TRUE)) {
  min_pos <- min(DATA$CPUE[DATA$CPUE > 0], na.rm = TRUE)  # Минимальный положительный улов
  offset <- min_pos / 2  # Величина поправки
  DATA <- DATA %>% 
    mutate(CPUE_POS = if_else(CPUE <= 0, CPUE + offset, CPUE))  # Добавляем поправку
} else {
  DATA <- DATA %>% 
    mutate(CPUE_POS = CPUE)  # Исходные данные если нулей нет
}

# Рассчитываем медианные значения CPUE по годам из исходных данных
actual_medians <- DATA %>%
  group_by(YEAR) %>%
  summarise(median_cpue = median(CPUE, na.rm = TRUE))
# Рассчитываем медианные значения CPUE по годам из исходных данных для последующих графиков
actual_medians

# Экспресс-визуализация распределения CPUE по годам
DATA %>%
  ggplot(aes(x = YEAR, y = CPUE)) +
  geom_boxplot(outlier.alpha = 0.2) +
  labs(title = "Распределение CPUE по годам", 
       x = "Год", 
       y = "CPUE (улов на усилие)")

# ==============================================================================
# БЛОК 3: ВСПОМОГАТЕЛЬНЫЕ ФУНКЦИИ ДЛЯ РАСЧЕТА ИНДЕКСОВ
# ==============================================================================

# Функция нормировки индексов
scale_to_index <- function(values, method = c("mean", "first")) {
  method <- match.arg(method)
  if (method == "mean") {
    # Нормировка на среднее значение
    return(as.numeric(values) / mean(as.numeric(values), na.rm = TRUE))
  }
  if (method == "first") {
    # Нормировка на значение первого года
    return(as.numeric(values) / as.numeric(values[1]))
  }
}

# Функция расчета индексов для GLM/GAM через маргинальные средние
emmeans_standardized_index <- function(model, variable = "YEAR") {
  out <- suppressWarnings(
    emmeans(model, 
            specs = as.formula(paste0("~ ", variable)), 
            type = "response")
  )
  df <- as_tibble(out) %>% 
    select(!!sym(variable), response = response, lower.CL, upper.CL)
  colnames(df) <- c("YEAR", "value", "lcl", "ucl")
  df
}

# Функция расчета индексов для GAMM через бутстреп
compute_standardized_index <- function(model, base_data, year_levels, predict_fun,
                                      response_transform = identity, 
                                      n_boot = 200L, 
                                      seed = 7L) {
  set.seed(seed)
  acc <- vector("list", length(year_levels))
  for (i in seq_along(year_levels)) {
    newdata <- base_data
    newdata$YEAR <- factor(year_levels[i], levels = levels(base_data$YEAR))
    preds <- suppressWarnings(predict_fun(model, newdata))
    mu <- mean(response_transform(preds), na.rm = TRUE)
    # Бутстреп для оценки неопределенности
    boot_vals <- replicate(n_boot, {
      idx <- sample.int(nrow(base_data), nrow(base_data), replace = TRUE)
      bd <- newdata[idx, , drop = FALSE]
      p <- suppressWarnings(predict_fun(model, bd))
      mean(response_transform(p), na.rm = TRUE)
    })
    ci <- quantile(boot_vals, c(0.025, 0.975), na.rm = TRUE)
    acc[[i]] <- tibble(YEAR = year_levels[i], value = mu, lcl = ci[[1]], ucl = ci[[2]])
  }
  bind_rows(acc)
}

# ==============================================================================
# БЛОК 4: МОДЕЛИРОВАНИЕ GLM (GAMMA С ЛОГ-ССЫЛКОЙ)
# ==============================================================================

# Подбор модели с фиксированными эффектами
glm_gamma_fit <- glm(
  CPUE_POS ~ YEAR + MONTH + CALL + REGION,  # Формула с факторными предикторами
  family = Gamma(link = "log"),            # Гамма-распределение с логарифмической связью
  data = DATA
)

# Диагностика модели
summary(glm_gamma_fit)  # Стандартная сводка модели

# Таблица коэффициентов в форматированном виде
broom::tidy(glm_gamma_fit) %>%
  mutate(across(estimate:statistic, ~round(.x, 4))) %>%
  kable(caption = "Коэффициенты GLM модели", align = "lrrrr")

# Графики диагностики остатков
par(mfrow = c(2, 2))
plot(glm_gamma_fit)  # Стандартные диагностические графики GLM
par(mfrow = c(1, 1))

# Диагностика остатков GLM с использованием DHARMa
sim_glm <- simulateResiduals(glm_gamma_fit, n = 1000, refit = FALSE)
plot(sim_glm, main = "GLM")

# Расчет и визуализация индексов
idx_glm <- emmeans_standardized_index(glm_gamma_fit) %>%
  mutate(model = "GLM_Gamma",
         index_mean = scale_to_index(value, "mean"),
         index_first = scale_to_index(value, "first"))

# Добавление доверительных интервалов
idx_glm <- idx_glm %>%
  mutate(
    lcl_index_mean = scale_to_index(lcl, "mean"),
    ucl_index_mean = scale_to_index(ucl, "mean")
  )

# Визуализация индексов GLM

idx_glm %>%
  ggplot(aes(x = YEAR, y = value, group = 1)) +
  geom_line() +
  geom_point() +
  geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.3) +
  geom_point(data = actual_medians, 
           aes(x = YEAR, y = median_cpue), 
           shape = 4,  # 4 соответствует крестику (x)
           size = 3, 
           color = "black", 
           inherit.aes = FALSE)+
labs(title = "Индексы CPUE по GLM модели (крестики - факт)", 
       x = "Год", 
       y = "Стандартизированный индекс")

# ==============================================================================
# БЛОК 5: МОДЕЛИРОВАНИЕ GAM
# ==============================================================================

# Подбор обобщенной аддитивной модели
gam_fit <- gam(
  CPUE_POS ~ YEAR + MONTH + CALL + REGION,  # Линейная формула без сглаживания
  family = Gamma(link = "log"),            # Аналогичное GLM распределение
  method = "REML",                         # Метод оптимизации гиперпараметров
  data = DATA
)

summary(gam_fit)  # Сводка модели

# Проверка адекватности модели (графики остатков)
mgcv::gam.check(gam_fit)

# Диагностика остатков GAM с использованием DHARMa
sim_gam <- simulateResiduals(gam_fit, n = 1000, refit = FALSE)
plot(sim_gam, main = "GAM")

# Расчет индексов
idx_gam <- emmeans_standardized_index(gam_fit) %>%
  mutate(model = "GAM",
         index_mean = scale_to_index(value, "mean"),
         index_first = scale_to_index(value, "first"))

# Доверительные интервалы
idx_gam <- idx_gam %>%
  mutate(
    lcl_index_mean = scale_to_index(lcl, "mean"),
    ucl_index_mean = scale_to_index(ucl, "mean")
  )


# Визуализация
idx_gam %>%
  ggplot(aes(x = YEAR, y = value, group = 1)) +
  geom_line() +
  geom_point() +
  geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.3) +
  geom_point(data = actual_medians, 
           aes(x = YEAR, y = median_cpue), 
           shape = 4,  # 4 соответствует крестику (x)
           size = 3, 
           color = "black", 
           inherit.aes = FALSE)+
  labs(title = "Индексы CPUE по GAM модели (крестики - факт", 
       x = "Год", 
       y = "Стандартизированный индекс")

# ==============================================================================
# БЛОК 6: МОДЕЛИРОВАНИЕ GAMM (СМЕШАННАЯ МОДЕЛЬ)
# ==============================================================================

# Подбор модели со смешанными эффектами
gamm_fit <- gamm4::gamm4(
  formula = CPUE_POS ~ YEAR + MONTH + REGION,  # Фиксированные эффекты
  random = ~ (1 | CALL),                       # Случайный эффект для судна
  family = Gamma(link = "log"),               # Распределение
  data = DATA
)

# 1. График остатков от предсказанных значений
plot(fitted(gamm_fit$gam), residuals(gamm_fit$gam, type = "deviance"),
     xlab = "Предсказанные значения", ylab = "Девиансные остатки",
     main = "Остатки GAMM vs. Предсказания")
abline(h = 0, col = "red", lty = 2)

# 2. QQ-plot для остатков
qqnorm(residuals(gamm_fit$gam, type = "deviance"),
       main = "QQ-plot для остатков GAMM")
qqline(residuals(gamm_fit$gam, type = "deviance"), col = "red")

# 3. Диагностика случайных эффектов
cat("\nСлучайные эффекты (CALL):\n")
print(summary(ranef(gamm_fit$mer)$CALL))

# График случайных эффектов
random_effects <- ranef(gamm_fit$mer)$CALL
plot(density(random_effects[,1]), main = "Распределение случайных эффектов",
     xlab = "Случайный эффект", ylab = "Плотность")

# 5. Проверка гетероскедастичности
library(lmtest)
bptest(gamm_fit$gam$y ~ fitted(gamm_fit$gam)) %>% 
  print()

# 6. Сводка по модели
cat("\nСводка GAMM модели:\n")
print(summary(gamm_fit$gam))
print(summary(gamm_fit$mer))

# Создание сетки для предсказания
newdata_grid <- expand.grid(
  YEAR = levels(DATA$YEAR),
  MONTH = levels(DATA$MONTH),
  REGION = levels(DATA$REGION),
  CALL = levels(DATA$CALL)[1]  # Фиксированное значение для случайного эффекта
)

# Предсказание на сетке
newdata_grid$pred <- predict(gamm_fit$gam, 
                            newdata = newdata_grid, 
                            type = "response")

# Усреднение предсказаний по годам
idx_gamm <- newdata_grid %>%
  group_by(YEAR) %>%
  summarise(value = mean(pred, na.rm = TRUE)) %>%
  mutate(
    model = "GAMM (mgcv)",
    index_mean = scale_to_index(value, "mean"),
    index_first = scale_to_index(value, "first")
  )

# Функция расчета доверительных интервалов через бутстреп
compute_gamm_ci <- function(model, newdata, n_boot = 100) {
  boot_means <- replicate(n_boot, {
    boot_data <- newdata[sample(nrow(newdata), replace = TRUE), ]
    preds <- predict(model, newdata = boot_data, type = "response")
    boot_data %>%
      mutate(pred = preds) %>%
      group_by(YEAR) %>%
      summarise(mean_pred = mean(pred, na.rm = TRUE)) %>%
      pull(mean_pred)
  })
  
  ci <- apply(boot_means, 1, function(x) quantile(x, c(0.025, 0.975), na.rm = TRUE))
  return(list(mean = rowMeans(boot_means), lcl = ci[1, ], ucl = ci[2, ]))
}

# Расчет интервалов
gamm_ci <- compute_gamm_ci(gamm_fit$gam, newdata_grid)

# Добавление интервалов к индексам
idx_gamm <- idx_gamm %>%
  mutate(
    lcl = gamm_ci$lcl,
    ucl = gamm_ci$ucl,
    lcl_index_mean = scale_to_index(lcl, "mean"),
    ucl_index_mean = scale_to_index(ucl, "mean")
  )

# Визуализация
idx_gamm %>%
  ggplot(aes(x = YEAR, y = value, group = 1)) +
  geom_line() +
  geom_point() +
  geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.3) +
  geom_point(data = actual_medians, 
           aes(x = YEAR, y = median_cpue), 
           shape = 4,  # 4 соответствует крестику (x)
           size = 3, 
           color = "black", 
           inherit.aes = FALSE)+
labs(title = "Индексы CPUE по GAMM модели (крестики - факт)", 
       x = "Год", 
       y = "Стандартизированный индекс")

# ==============================================================================
# БЛОК 7: СРАВНЕНИЕ МОДЕЛЕЙ И ФИНАЛЬНАЯ ВИЗУАЛИЗАЦИЯ
# ==============================================================================

# Объединение результатов всех моделей
indices_all <- bind_rows(
  idx_glm %>% select(YEAR, value, lcl, ucl, model, index_mean, lcl_index_mean, ucl_index_mean),
  idx_gam %>% select(YEAR, value, lcl, ucl, model, index_mean, lcl_index_mean, ucl_index_mean),
  idx_gamm %>% select(YEAR, value, lcl, ucl, model, index_mean, lcl_index_mean, ucl_index_mean)
) %>% mutate(YEAR = factor(YEAR, levels = levels(DATA$YEAR)))

# Сводная таблица результатов
indices_all %>% 
  kable(caption = "Сравнение индексов CPUE по разным моделям")


indices_all %>%
  ggplot(aes(x = YEAR, y = value, color = model, group = model, fill = model)) +
  geom_line() +
  geom_point() +
  geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.1, linetype = "dashed") +
geom_point(data = actual_medians, 
           aes(x = YEAR, y = median_cpue), 
           shape = 4,  # 4 соответствует крестику (x)
           size = 3, 
           color = "black", 
           inherit.aes = FALSE)+
  labs(title = "Сравнение стандартизированных индексов CPUE (крестики - факт)", 
       x = "Год", 
       y = "Индекс CPUE (кг/ловушку)", 
       color = "Модель", 
       fill = "Модель") +
  theme(legend.position = "bottom")


# ==============================================================================
# БЛОК 9: СРАВНЕНИЕ МОДЕЛЕЙ ПО ИНФОРМАЦИОННЫМ КРИТЕРИЯМ
# ==============================================================================

cat("\n=== СРАВНИТЕЛЬНЫЙ АНАЛИЗ МОДЕЛЕЙ ПО ИНФОРМАЦИОННЫМ КРИТЕРИЯМ ===\n")

# Упрощенная функция для извлечения ключевых критериев из моделей
extract_model_metrics <- function(model, model_name, model_type = "glm") {
  if (model_type == "glm") {
    aic_val <- AIC(model)
    bic_val <- BIC(model)
    loglik_val <- as.numeric(logLik(model))
    df_val <- model$rank
    null_dev <- model$null.deviance
    dev <- model$deviance
  } else if (model_type == "gam") {
    aic_val <- AIC(model)
    bic_val <- BIC(model)
    loglik_val <- as.numeric(logLik(model))
    df_val <- sum(model$edf)
    null_dev <- model$null.deviance
    dev <- model$deviance
  } else if (model_type == "gamm") {
    aic_val <- AIC(model$mer)
    bic_val <- BIC(model$mer)
    loglik_val <- as.numeric(logLik(model$mer))
    df_val <- length(fixef(model$mer)) + 1  # +1 для случайного эффекта
    null_dev <- model$gam$null.deviance
    dev <- model$gam$deviance
  }
  
  # Вычисляем долю объясненной девиации
  deviance_explained <- ifelse(!is.null(null_dev) && !is.null(dev) && null_dev > 0,
                              (null_dev - dev) / null_dev, NA)
  
  data.frame(
    Model = model_name,
    AIC = round(aic_val, 2),
    BIC = round(bic_val, 2),
    LogLik = round(loglik_val, 2),
    DF = round(df_val, 2),
    Deviance_Explained = round(deviance_explained, 4)
  )
}

# Извлекаем метрики для всех моделей
model_metrics <- bind_rows(
  extract_model_metrics(glm_gamma_fit, "GLM (Gamma)", "glm"),
  extract_model_metrics(gam_fit, "GAM", "gam"),
  extract_model_metrics(gamm_fit, "GAMM", "gamm")
)

# Добавляем разницу в AIC относительно наилучшей модели
min_aic <- min(model_metrics$AIC)
model_metrics <- model_metrics %>%
  mutate(Delta_AIC = AIC - min_aic,
         AIC_Weight = exp(-0.5 * Delta_AIC) / sum(exp(-0.5 * Delta_AIC)))

# Форматируем таблицу для вывода
comparison_table <- model_metrics %>%
  mutate(across(where(is.numeric), ~round(., 3))) %>%
  arrange(AIC)  # Сортируем по AIC (лучшая модель первая)

# Выводим таблицу сравнения
cat("\nТАБЛИЦА СРАВНЕНИЯ МОДЕЛЕЙ:\n")
print(comparison_table)


# Выводим итоговые рекомендации
cat("\n=== ИТОГОВЫЕ РЕКОМЕНДАЦИИ ПО ВЫБОРУ МОДЕЛИ ===\n")
best_model <- comparison_table$Model[1]
cat("Наилучшая модель по критерию AIC:", best_model, "\n")
cat("Вес AIC для наилучшей модели:", round(comparison_table$AIC_Weight[1], 3), "\n")

if (nrow(comparison_table) > 1 && comparison_table$Delta_AIC[2] > 2) {
  cat("Наилучшая модель существенно лучше остальных (?AIC > 2).\n")
} else if (nrow(comparison_table) > 1) {
  cat("Несколько моделей имеют сходное качество (?AIC < 2).\n")
}

cat("Доля объясненной девиации наилучшей модели:", 
    round(comparison_table$Deviance_Explained[1], 3), "\n")






