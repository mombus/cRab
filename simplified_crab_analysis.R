# УПРОЩЕННАЯ ВЕРСИЯ АНАЛИЗА ПОПУЛЯЦИИ КАМЧАТСКОГО КРАБА
# Использует только базовые пакеты R для демонстрации основных улучшений

# Очистка рабочей среды
rm(list = ls())

cat("=== РАСШИРЕННЫЙ АНАЛИЗ ДИНАМИКИ ЗАПАСА КАМЧАТСКОГО КРАБА ===\n")

#### 1. ПОДГОТОВКА РАСШИРЕННЫХ ДАННЫХ ####

# Расширенные данные индекса обилия с дополнительными переменными
index_data <- data.frame(
  YEAR = 2015:2024,
  INDEX = c(1850000, 2100000, 1950000, 2381774, 1634549, 1920507, 1036673, 1147685, 1055733, 980000),
  TEMPERATURE = c(2.1, 1.8, 2.3, 2.0, 1.9, 2.2, 1.7, 1.8, 1.9, 2.1),
  SALINITY = c(34.2, 34.1, 34.3, 34.0, 34.2, 34.1, 34.4, 34.3, 34.2, 34.1),
  FISHING_EFFORT = c(120, 135, 110, 150, 180, 160, 200, 190, 195, 210)
)

cat("Период наблюдений:", min(index_data$YEAR), "-", max(index_data$YEAR), "\n")
cat("Количество наблюдений:", nrow(index_data), "\n")
cat("Общее изменение запаса:", 
    round((index_data$INDEX[10] - index_data$INDEX[1])/index_data$INDEX[1]*100, 1), "%\n")

#### 2. РАСШИРЕННАЯ СТАТИСТИКА ####

cat("\n=== РАСШИРЕННЫЕ СТАТИСТИЧЕСКИЕ ХАРАКТЕРИСТИКИ ===\n")

# Базовые статистики
mean_index <- mean(index_data$INDEX)
median_index <- median(index_data$INDEX)
sd_index <- sd(index_data$INDEX)
cv_index <- sd_index / mean_index * 100

cat("Среднее значение:", round(mean_index, 0), "\n")
cat("Медиана:", round(median_index, 0), "\n")
cat("Стандартное отклонение:", round(sd_index, 0), "\n")
cat("Коэффициент вариации:", round(cv_index, 1), "%\n")

# Анализ экстремумов
min_year <- index_data$YEAR[which.min(index_data$INDEX)]
max_year <- index_data$YEAR[which.max(index_data$INDEX)]
cat("Минимальное значение:", min(index_data$INDEX), "в", min_year, "году\n")
cat("Максимальное значение:", max(index_data$INDEX), "в", max_year, "году\n")

# Анализ трендов
trend_5yr <- (index_data$INDEX[10] - index_data$INDEX[6]) / index_data$INDEX[6] * 100
trend_10yr <- (index_data$INDEX[10] - index_data$INDEX[1]) / index_data$INDEX[1] * 100
cat("Изменение за последние 4 года:", round(trend_5yr, 1), "%\n")
cat("Изменение за весь период:", round(trend_10yr, 1), "%\n")

#### 3. КОРРЕЛЯЦИОННЫЙ АНАЛИЗ ####

cat("\n=== КОРРЕЛЯЦИОННЫЙ АНАЛИЗ ===\n")

# Корреляции между переменными
temp_cor <- cor(index_data$INDEX, index_data$TEMPERATURE)
sal_cor <- cor(index_data$INDEX, index_data$SALINITY)
fish_cor <- cor(index_data$INDEX, index_data$FISHING_EFFORT)

cat("Корреляция с температурой:", round(temp_cor, 3), "\n")
cat("Корреляция с соленостью:", round(sal_cor, 3), "\n")
cat("Корреляция с промысловым усилием:", round(fish_cor, 3), "\n")

# Интерпретация корреляций
interpret_correlation <- function(cor_value) {
  abs_cor <- abs(cor_value)
  if(abs_cor > 0.7) return("сильная")
  if(abs_cor > 0.5) return("умеренная")
  if(abs_cor > 0.3) return("слабая")
  return("очень слабая")
}

cat("Влияние температуры:", interpret_correlation(temp_cor), "\n")
cat("Влияние солености:", interpret_correlation(sal_cor), "\n")
cat("Влияние промысла:", interpret_correlation(fish_cor), "\n")

#### 4. БИОЛОГИЧЕСКИЕ ИНДЕКСЫ ####

cat("\n=== БИОЛОГИЧЕСКИЕ ИНДЕКСЫ ===\n")

# Индекс стабильности популяции
stability_index <- 1 - (sd_index / mean_index)
cat("Индекс стабильности популяции:", round(stability_index, 3), "\n")

# Анализ восстановления после снижений
declines <- which(diff(index_data$INDEX) < 0)
if(length(declines) > 0) {
  recovery_rates <- numeric()
  for(i in declines) {
    if(i < length(index_data$INDEX) - 1) {
      recovery_rate <- (index_data$INDEX[i+2] - index_data$INDEX[i+1]) / index_data$INDEX[i+1]
      recovery_rates <- c(recovery_rates, recovery_rate)
    }
  }
  recovery_index <- mean(recovery_rates, na.rm = TRUE)
  cat("Индекс восстановления:", round(recovery_index, 3), "\n")
} else {
  cat("Индекс восстановления: не применимо (нет снижений)\n")
}

# Устойчивость к промыслу
fishing_resistance <- 1 - abs(fish_cor)
cat("Устойчивость к промыслу:", round(fishing_resistance, 3), "\n")

#### 5. РАСШИРЕННЫЕ ПРОГНОСТИЧЕСКИЕ МОДЕЛИ ####

cat("\n=== РАСШИРЕННОЕ МОДЕЛИРОВАНИЕ ===\n")

# Создание временного ряда
ts_data <- ts(index_data$INDEX, start = 2015, frequency = 1)

# Разделение на обучающую и тестовую выборки
train <- window(ts_data, end = 2022)
test <- window(ts_data, start = 2023)

cat("Обучающая выборка: 2015-2022 (", length(train), "наблюдений)\n")
cat("Тестовая выборка: 2023-2024 (", length(test), "наблюдений)\n")

### Модель 1: Наивный прогноз ###
naive_forecast <- rep(train[length(train)], length(test))
cat("Наивный прогноз:", round(naive_forecast[1]/1e6, 2), "млн экз.\n")

### Модель 2: Линейный тренд ###
# Простая линейная регрессия
years_train <- 2015:2022
lm_model <- lm(train ~ years_train)
trend_forecast <- predict(lm_model, newdata = data.frame(years_train = 2023:2024))
cat("Линейный тренд:", round(trend_forecast[1]/1e6, 2), "млн экз.\n")

### Модель 3: Экспоненциальное сглаживание ###
# Простое экспоненциальное сглаживание
alpha <- 0.3  # Параметр сглаживания
exp_smooth <- numeric(length(train))
exp_smooth[1] <- train[1]
for(i in 2:length(train)) {
  exp_smooth[i] <- alpha * train[i] + (1 - alpha) * exp_smooth[i-1]
}
exp_forecast <- exp_smooth[length(exp_smooth)]
cat("Экспоненциальное сглаживание:", round(exp_forecast/1e6, 2), "млн экз.\n")

### Модель 4: Модель с внешними факторами ###
# Регрессия с учетом экологических факторов
reg_data <- data.frame(
  INDEX = index_data$INDEX[1:8],
  TEMPERATURE = index_data$TEMPERATURE[1:8],
  FISHING_EFFORT = index_data$FISHING_EFFORT[1:8],
  YEAR = index_data$YEAR[1:8]
)

reg_model <- lm(INDEX ~ TEMPERATURE + FISHING_EFFORT + YEAR, data = reg_data)
reg_forecast <- predict(reg_model, 
                       newdata = data.frame(
                         TEMPERATURE = index_data$TEMPERATURE[9],
                         FISHING_EFFORT = index_data$FISHING_EFFORT[9],
                         YEAR = index_data$YEAR[9]
                       ))
cat("Регрессионная модель:", round(reg_forecast/1e6, 2), "млн экз.\n")

### Модель 5: Комбинированный прогноз ###
combined_forecast <- (naive_forecast[1] + trend_forecast[1] + exp_forecast + reg_forecast) / 4
cat("Комбинированный прогноз:", round(combined_forecast/1e6, 2), "млн экз.\n")

#### 6. ОЦЕНКА КАЧЕСТВА МОДЕЛЕЙ ####

cat("\n=== ОЦЕНКА КАЧЕСТВА МОДЕЛЕЙ ===\n")

# Функции для расчета метрик
calc_mape <- function(actual, forecast) {
  mean(abs((actual - forecast) / actual)) * 100
}

calc_rmse <- function(actual, forecast) {
  sqrt(mean((actual - forecast)^2))
}

# Сбор метрик качества
actual_2023 <- as.numeric(test[1])
actual_2024 <- as.numeric(test[2])

models <- c("Наивная", "Линейный тренд", "Экспоненциальное сглаживание", "Регрессия", "Комбинированная")
forecasts <- c(naive_forecast[1], trend_forecast[1], exp_forecast, reg_forecast, combined_forecast)

results <- data.frame(
  Модель = models,
  Прогноз_2023 = round(forecasts/1e6, 2),
  Факт_2023 = rep(round(actual_2023/1e6, 2), 5),
  Ошибка_2023 = round(abs(forecasts - actual_2023)/1e6, 2),
  MAPE_2023 = round(sapply(forecasts, function(f) calc_mape(actual_2023, f)), 1),
  RMSE_2023 = round(sapply(forecasts, function(f) calc_rmse(actual_2023, f))/1e6, 2)
)

cat("Сравнение моделей на тестовых данных (2023 год):\n")
print(results)

# Выбор лучшей модели по MAPE
best_model_idx <- which.min(results$MAPE_2023)
best_model_name <- results$Модель[best_model_idx]

cat("\nЛучшая модель:", best_model_name, "(MAPE =", results$MAPE_2023[best_model_idx], "%)\n")

#### 7. СЦЕНАРНОЕ ПРОГНОЗИРОВАНИЕ ####

cat("\n=== СЦЕНАРНОЕ ПРОГНОЗИРОВАНИЕ ===\n")

# Финальный прогноз на 2025-2027 с использованием всей доступной информации
if(best_model_name == "Наивная") {
  final_forecast <- rep(index_data$INDEX[nrow(index_data)], 3)
} else if(best_model_name == "Линейный тренд") {
  years_all <- 2015:2024
  lm_final <- lm(index_data$INDEX ~ years_all)
  final_forecast <- predict(lm_final, newdata = data.frame(years_all = 2025:2027))
} else if(best_model_name == "Экспоненциальное сглаживание") {
  # Используем последнее сглаженное значение
  final_forecast <- rep(exp_smooth[length(exp_smooth)], 3)
} else if(best_model_name == "Регрессия") {
  # Используем регрессионную модель
  reg_final <- lm(INDEX ~ TEMPERATURE + FISHING_EFFORT + YEAR, data = index_data)
  final_forecast <- predict(reg_final, newdata = data.frame(
    TEMPERATURE = rep(mean(index_data$TEMPERATURE), 3),
    FISHING_EFFORT = rep(mean(index_data$FISHING_EFFORT), 3),
    YEAR = 2025:2027
  ))
} else {
  # Комбинированный прогноз
  final_forecast <- rep(combined_forecast, 3)
}

# Обеспечиваем неотрицательность прогноза
final_forecast <- pmax(final_forecast, 0)

# Расчет доверительных интервалов (упрощенный)
forecast_sd <- sd_index * 0.3  # Предполагаем 30% неопределенности
forecast_lower_80 <- final_forecast - 1.28 * forecast_sd
forecast_upper_80 <- final_forecast + 1.28 * forecast_sd
forecast_lower_95 <- final_forecast - 1.96 * forecast_sd
forecast_upper_95 <- final_forecast + 1.96 * forecast_sd

# Создание сценариев
scenarios <- data.frame(
  Год = 2025:2027,
  Базовый_сценарий = round(final_forecast/1e6, 2),
  Пессимистичный = round(pmax(forecast_lower_95/1e6, 0), 2),
  Оптимистичный = round(forecast_upper_95/1e6, 2),
  Консервативный = round(pmax(forecast_lower_80/1e6, 0), 2),
  Агрессивный = round(forecast_upper_80/1e6, 2)
)

cat("Сценарные прогнозы (млн экз.):\n")
print(scenarios)

#### 8. РАСШИРЕННЫЙ РИСК-АНАЛИЗ ####

cat("\n=== РАСШИРЕННЫЙ РИСК-АНАЛИЗ ===\n")

# Анализ устойчивости запаса
current_level <- index_data$INDEX[nrow(index_data)]
historical_mean <- mean(index_data$INDEX)
current_ratio <- current_level / historical_mean

cat("Текущий уровень запаса (2024):", round(current_level/1e6, 2), "млн экз.\n")
cat("Отношение к историческому среднему:", round(current_ratio, 2), "\n")

# Анализ волатильности
volatility <- sd(diff(index_data$INDEX)/index_data$INDEX[-1], na.rm = TRUE) * 100
cat("Волатильность изменений:", round(volatility, 1), "%\n")

# Определение статуса запаса
if(trend_5yr < -30) {
  status <- "КРИТИЧЕСКИЙ"
  risk_level <- "ОЧЕНЬ ВЫСОКИЙ"
} else if(trend_5yr < -20) {
  status <- "ТРЕВОЖНЫЙ"
  risk_level <- "ВЫСОКИЙ"
} else if(trend_5yr < -10) {
  status <- "ОСТОРОЖНЫЙ"
  risk_level <- "УМЕРЕННЫЙ"
} else if(trend_5yr < 0) {
  status <- "СТАБИЛЬНЫЙ"
  risk_level <- "НИЗКИЙ"
} else {
  status <- "РАСТУЩИЙ"
  risk_level <- "НИЗКИЙ"
}

cat("СТАТУС ЗАПАСА:", status, "\n")
cat("УРОВЕНЬ РИСКА:", risk_level, "\n")

#### 9. ПРОМЫСЛОВЫЕ РЕКОМЕНДАЦИИ ####

cat("\n=== ПРОМЫСЛОВЫЕ РЕКОМЕНДАЦИИ ===\n")

# Расчет оптимального промыслового усилия
if(current_ratio < 0.5) {
  fishing_level <- "МИНИМАЛЬНЫЙ"
  effort_reduction <- 0.8
  priority <- "ВЫСОКИЙ"
} else if(current_ratio < 0.7) {
  fishing_level <- "ОГРАНИЧЕННЫЙ"
  effort_reduction <- 0.5
  priority <- "ВЫСОКИЙ"
} else if(current_ratio < 0.9) {
  fishing_level <- "УМЕРЕННЫЙ"
  effort_reduction <- 0.25
  priority <- "СРЕДНИЙ"
} else {
  fishing_level <- "СТАНДАРТНЫЙ"
  effort_reduction <- 0
  priority <- "НИЗКИЙ"
}

cat("Рекомендуемый уровень промысла:", fishing_level, "\n")
cat("Снижение промыслового усилия:", round(effort_reduction * 100, 0), "%\n")
cat("Приоритет мониторинга:", priority, "\n")

# Расчет рекомендуемого улова
recommended_catch_2025 <- final_forecast[1] * 0.1  # 10% от прогнозируемого запаса
cat("Рекомендуемый улов на 2025:", round(recommended_catch_2025/1e6, 2), "млн экз.\n")

#### 10. ДЕТАЛЬНЫЕ РЕКОМЕНДАЦИИ ####

cat("\n=== ДЕТАЛЬНЫЕ РЕКОМЕНДАЦИИ ===\n")

if(status == "КРИТИЧЕСКИЙ") {
  cat("НЕМЕДЛЕННЫЕ ДЕЙСТВИЯ:\n")
  cat("- Полная остановка промысла\n")
  cat("- Экстренный мониторинг\n")
  cat("- Исследование причин катастрофического снижения\n")
  cat("- Разработка плана восстановления\n")
} else if(status == "ТРЕВОЖНЫЙ") {
  cat("СРОЧНЫЕ МЕРЫ:\n")
  cat("- Снижение промысловой нагрузки на 50%\n")
  cat("- Усиление мониторинга (ежемесячно)\n")
  cat("- Анализ влияния экологических факторов\n")
  cat("- Подготовка плана действий на случай дальнейшего снижения\n")
} else if(status == "ОСТОРОЖНЫЙ") {
  cat("ПРЕДОСТОРОЖНЫЕ МЕРЫ:\n")
  cat("- Снижение промысловой нагрузки на 25%\n")
  cat("- Квартальный мониторинг\n")
  cat("- Анализ трендов\n")
} else if(status == "СТАБИЛЬНЫЙ") {
  cat("ПОДДЕРЖИВАЮЩИЕ МЕРЫ:\n")
  cat("- Поддержание текущего уровня промысла\n")
  cat("- Полугодовой мониторинг\n")
  cat("- Отслеживание изменений\n")
} else {
  cat("РАЗВИВАЮЩИЕ МЕРЫ:\n")
  cat("- Постепенное увеличение промысловой нагрузки\n")
  cat("- Годовой мониторинг\n")
  cat("- Исследование возможностей устойчивого роста\n")
}

#### 11. СОЗДАНИЕ ОТЧЕТА ####

cat("\n=== СОЗДАНИЕ ОТЧЕТА ===\n")

# Создание сводного отчета
final_report <- list(
  metadata = list(
    анализ_проведен = Sys.time(),
    период_данных = paste(min(index_data$YEAR), "-", max(index_data$YEAR)),
    количество_наблюдений = nrow(index_data),
    лучшая_модель = best_model_name
  ),
  статистика = list(
    среднее_значение = round(mean_index, 0),
    медиана = round(median_index, 0),
    стандартное_отклонение = round(sd_index, 0),
    коэффициент_вариации = round(cv_index, 1),
    минимальное_значение = min(index_data$INDEX),
    максимальное_значение = max(index_data$INDEX)
  ),
  тренды = list(
    изменение_за_период = round(trend_10yr, 1),
    изменение_за_последние_4_года = round(trend_5yr, 1),
    волатильность = round(volatility, 1)
  ),
  биологические_индексы = list(
    стабильность = round(stability_index, 3),
    устойчивость_к_промыслу = round(fishing_resistance, 3)
  ),
  экологические_факторы = list(
    влияние_температуры = round(temp_cor, 3),
    влияние_солености = round(sal_cor, 3),
    влияние_промысла = round(fish_cor, 3)
  ),
  сравнение_моделей = results,
  прогноз = scenarios,
  статус_запаса = list(
    статус = status,
    уровень_риска = risk_level,
    текущий_уровень = round(current_level/1e6, 2)
  ),
  промысловые_рекомендации = list(
    уровень_промысла = fishing_level,
    снижение_усилия = round(effort_reduction * 100, 0),
    приоритет_мониторинга = priority,
    рекомендуемый_улов_2025 = round(recommended_catch_2025/1e6, 2)
  )
)

# Сохранение отчета
saveRDS(final_report, "crab_analysis_simplified_report.rds")

# Создание текстового отчета
report_text <- paste0("
ОТЧЕТ ПО АНАЛИЗУ ПОПУЛЯЦИИ КАМЧАТСКОГО КРАБА
===============================================

Дата анализа: ", Sys.time(), "
Период данных: ", min(index_data$YEAR), "-", max(index_data$YEAR), "
Лучшая модель: ", best_model_name, "

СТАТУС ЗАПАСА
------------
Статус: ", status, "
Уровень риска: ", risk_level, "
Текущий уровень: ", round(current_level/1e6, 2), " млн экз.

ПРОГНОЗ НА 2025-2027 ГГ.
-----------------------
Год    Базовый    Пессимистичный    Оптимистичный
2025   ", scenarios$Базовый_сценарий[1], "        ", scenarios$Пессимистичный[1], "            ", scenarios$Оптимистичный[1], "
2026   ", scenarios$Базовый_сценарий[2], "        ", scenarios$Пессимистичный[2], "            ", scenarios$Оптимистичный[2], "
2027   ", scenarios$Базовый_сценарий[3], "        ", scenarios$Пессимистичный[3], "            ", scenarios$Оптимистичный[3], "

ПРОМЫСЛОВЫЕ РЕКОМЕНДАЦИИ
-----------------------
Уровень промысла: ", fishing_level, "
Снижение усилия: ", round(effort_reduction * 100, 0), "%
Приоритет мониторинга: ", priority, "
Рекомендуемый улов 2025: ", round(recommended_catch_2025/1e6, 2), " млн экз.

БИОЛОГИЧЕСКИЕ ИНДЕКСЫ
--------------------
Стабильность популяции: ", round(stability_index, 3), "
Устойчивость к промыслу: ", round(fishing_resistance, 3), "

ЭКОЛОГИЧЕСКИЕ ФАКТОРЫ
--------------------
Влияние температуры: ", round(temp_cor, 3), " (", interpret_correlation(temp_cor), ")
Влияние солености: ", round(sal_cor, 3), " (", interpret_correlation(sal_cor), ")
Влияние промысла: ", round(fish_cor, 3), " (", interpret_correlation(fish_cor), ")
")

writeLines(report_text, "crab_analysis_report.txt")

cat("\n=== АНАЛИЗ ЗАВЕРШЕН ===\n")
cat("Лучшая модель:", best_model_name, "\n")
cat("Статус запаса:", status, "\n")
cat("Прогноз на 2025:", scenarios$Базовый_сценарий[1], "млн экз.\n")
cat("Диапазон неопределенности (95% ДИ):", 
    scenarios$Пессимистичный[1], "-", scenarios$Оптимистичный[1], "млн экз.\n")
cat("Отчет сохранен в файлы: crab_analysis_simplified_report.rds и crab_analysis_report.txt\n")

# Создание CSV с основными результатами
summary_data <- data.frame(
  Метрика = c("Период анализа", "Количество наблюдений", "Средний индекс", 
              "Коэффициент вариации", "Лучшая модель", "Прогноз 2025", "Статус запаса"),
  Значение = c(paste(min(index_data$YEAR), "-", max(index_data$YEAR)), 
              nrow(index_data),
              paste0(round(mean_index/1e6, 2), " млн экз."),
              paste0(round(cv_index, 1), "%"),
              best_model_name,
              paste0(scenarios$Базовый_сценарий[1], " млн экз."),
              status)
)

write.csv(summary_data, "crab_analysis_summary.csv", row.names = FALSE, fileEncoding = "UTF-8")

cat("Сводка сохранена в файл: crab_analysis_summary.csv\n")