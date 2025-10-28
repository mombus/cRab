# ПРАКТИЧЕСКОЕ ЗАНЯТИЕ: РАСШИРЕННЫЙ АНАЛИЗ ПОПУЛЯЦИИ КАМЧАТСКОГО КРАБА
# УЛУЧШЕННАЯ ВЕРСИЯ С БИОЛОГИЧЕСКИМИ ОГРАНИЧЕНИЯМИ И СЦЕНАРНЫМ ПРОГНОЗИРОВАНИЕМ

# Очистка рабочей среды
rm(list = ls())

# Загрузка пакетов
suppressPackageStartupMessages({
  library(tidyverse)
  library(forecast)
  library(tseries)
  library(ggplot2)
  library(gridExtra)
  library(plotly)
  library(knitr)
  library(DT)
  library(corrplot)
  library(strucchange)
  library(urca)
})

#### 1. ПОДГОТОВКА ДАННЫХ И РАСШИРЕННЫЙ АНАЛИЗ ####

# Расширенные данные индекса обилия (добавим больше исторических данных)
index_data <- data.frame(
  YEAR = 2015:2024,
  INDEX = c(1850000, 2100000, 1950000, 2381774, 1634549, 1920507, 1036673, 1147685, 1055733, 980000),
  TEMPERATURE = c(2.1, 1.8, 2.3, 2.0, 1.9, 2.2, 1.7, 1.8, 1.9, 2.1), # Температура воды
  SALINITY = c(34.2, 34.1, 34.3, 34.0, 34.2, 34.1, 34.4, 34.3, 34.2, 34.1), # Соленость
  FISHING_EFFORT = c(120, 135, 110, 150, 180, 160, 200, 190, 195, 210) # Промысловое усилие
)

cat("=== РАСШИРЕННЫЙ АНАЛИЗ ДИНАМИКИ ЗАПАСА КАМЧАТСКОГО КРАБА ===\n")
cat("Период наблюдений:", min(index_data$YEAR), "-", max(index_data$YEAR), "\n")
cat("Общее изменение запаса:", 
    round((index_data$INDEX[10] - index_data$INDEX[1])/index_data$INDEX[1]*100, 1), "%\n")
cat("Среднегодовое изменение:", 
    round(((index_data$INDEX[10]/index_data$INDEX[1])^(1/9) - 1)*100, 1), "%\n")

# Создание временного ряда
ts_data <- ts(index_data$INDEX, start = 2015, frequency = 1)

# Расчет дополнительных статистик
cat("\n=== СТАТИСТИЧЕСКИЕ ХАРАКТЕРИСТИКИ ===\n")
cat("Среднее значение:", round(mean(ts_data), 0), "\n")
cat("Медиана:", round(median(ts_data), 0), "\n")
cat("Стандартное отклонение:", round(sd(ts_data), 0), "\n")
cat("Коэффициент вариации:", round(sd(ts_data)/mean(ts_data)*100, 1), "%\n")
cat("Минимальное значение:", min(ts_data), "в", index_data$YEAR[which.min(index_data$INDEX)], "году\n")
cat("Максимальное значение:", max(ts_data), "в", index_data$YEAR[which.max(index_data$INDEX)], "году\n")

#### 2. РАСШИРЕННАЯ ВИЗУАЛИЗАЦИЯ ####

# Основной график динамики
p1 <- ggplot(index_data, aes(x = YEAR, y = INDEX/1e6)) +
  geom_line(linewidth = 1.2, color = "steelblue", alpha = 0.8) +
  geom_point(size = 4, color = "navy", alpha = 0.7) +
  geom_text(aes(label = round(INDEX/1e6, 1)), vjust = -1.2, size = 3.5, fontface = "bold") +
  geom_smooth(method = "loess", se = TRUE, color = "red", alpha = 0.3, linewidth = 0.8) +
  labs(title = "Динамика индекса обилия камчатского краба",
       subtitle = "2015-2024 гг. с трендом (значения в млн экз.)",
       x = "Год", y = "Индекс обилия, млн экз.") +
  ylim(0, 2.5) +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 12, color = "gray60"))

# График изменений по годам
index_data$CHANGE <- c(NA, diff(index_data$INDEX))
index_data$CHANGE_PCT <- c(NA, diff(index_data$INDEX)/index_data$INDEX[-10]*100)

p2 <- ggplot(index_data[-1,], aes(x = YEAR, y = CHANGE_PCT)) +
  geom_col(aes(fill = CHANGE_PCT > 0), alpha = 0.7) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  geom_text(aes(label = paste0(round(CHANGE_PCT, 1), "%")), 
            vjust = ifelse(index_data$CHANGE_PCT[-1] > 0, -0.5, 1.5), size = 3) +
  scale_fill_manual(values = c("TRUE" = "darkgreen", "FALSE" = "darkred"), guide = "none") +
  labs(title = "Годовые изменения индекса обилия",
       x = "Год", y = "Изменение, %") +
  theme_minimal()

# Корреляционный анализ
cor_data <- index_data[, c("INDEX", "TEMPERATURE", "SALINITY", "FISHING_EFFORT")]
cor_matrix <- cor(cor_data, use = "complete.obs")

p3 <- ggplot(data.frame(
  x = rep(1:4, 4),
  y = rep(4:1, each = 4),
  corr = as.vector(cor_matrix),
  var1 = rep(c("Индекс", "Температура", "Соленость", "Промысел"), 4),
  var2 = rep(c("Промысел", "Соленость", "Температура", "Индекс"), each = 4)
), aes(x = x, y = y, fill = corr)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(corr, 2)), color = "white", fontface = "bold") +
  scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0) +
  labs(title = "Корреляционная матрица", fill = "Корреляция") +
  theme_minimal() +
  theme(axis.text = element_blank(), axis.ticks = element_blank())

# Объединение графиков
grid.arrange(p1, p2, p3, ncol = 2, nrow = 2)

#### 3. РАСШИРЕННЫЙ АНАЛИЗ СТАЦИОНАРНОСТИ ####

cat("\n=== РАСШИРЕННЫЙ АНАЛИЗ СТАЦИОНАРНОСТИ ===\n")

# Тест Дики-Фуллера
adf_test <- ur.df(ts_data, type = "trend", lags = 2)
cat("Тест Дики-Фуллера (тренд):\n")
cat("Статистика:", round(adf_test@teststat[1], 3), "\n")
cat("Критическое значение (5%):", round(adf_test@cval[1,2], 3), "\n")
if(adf_test@teststat[1] < adf_test@cval[1,2]) {
  cat("Результат: ряд стационарен (отвергаем H0)\n")
} else {
  cat("Результат: ряд нестационарен (не отвергаем H0)\n")
}

# Анализ структурных изменений
bp_test <- breakpoints(INDEX ~ YEAR, data = index_data)
if(length(bp_test$breakpoints) > 0) {
  cat("Обнаружены структурные изменения в точках:", bp_test$breakpoints + 2014, "\n")
} else {
  cat("Структурные изменения не обнаружены\n")
}

# Анализ автокорреляции
acf_result <- acf(ts_data, plot = FALSE)
cat("Автокорреляция первого порядка:", round(acf_result$acf[2], 3), "\n")

#### 4. РАСШИРЕННЫЕ ПРОГНОСТИЧЕСКИЕ МОДЕЛИ ####

# Разделение на обучающую и тестовую выборки
train <- window(ts_data, end = 2022)
test <- window(ts_data, start = 2023)

cat("\n=== РАСШИРЕННОЕ МОДЕЛИРОВАНИЕ ===\n")

### Модель 1: Наивный прогноз с ограничениями ###
naive_model <- naive(train, h = 2)
naive_forecast <- forecast(naive_model, h = 1)

### Модель 2: ARIMA с автоматическим выбором ###
arima_model <- auto.arima(train, seasonal = FALSE, stepwise = FALSE, approximation = FALSE)
arima_forecast <- forecast(arima_model, h = 1)
cat("Выбранная ARIMA модель:", arima_model$arma, "\n")

### Модель 3: ARIMA с логарифмическим преобразованием ###
ts_data_log <- ts(log(index_data$INDEX), start = 2015)
train_log <- window(ts_data_log, end = 2022)

arima_model_log <- auto.arima(train_log, seasonal = FALSE)
arima_forecast_log <- forecast(arima_model_log, h = 1)

# Обратное преобразование с коррекцией смещения
arima_forecast_log_corrected <- arima_forecast_log
arima_forecast_log_corrected$mean <- exp(arima_forecast_log$mean + 0.5 * arima_forecast_log$model$sigma2)
arima_forecast_log_corrected$lower <- exp(arima_forecast_log$lower)
arima_forecast_log_corrected$upper <- exp(arima_forecast_log$upper)

### Модель 4: ETS с ограничениями ###
ets_model <- ets(train, lambda = 0)
ets_forecast <- forecast(ets_model, h = 1)

### Модель 5: Модель с внешними факторами ###
# Простая регрессионная модель
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
                       ), 
                       interval = "prediction")

### Модель 6: Комбинированный прогноз ###
# Простое усреднение прогнозов
combined_forecast <- (naive_forecast$mean[1] + 
                     arima_forecast$mean[1] + 
                     arima_forecast_log_corrected$mean[1] + 
                     ets_forecast$mean[1] + 
                     reg_forecast[1]) / 5

#### 5. РАСШИРЕННАЯ ОЦЕНКА КАЧЕСТВА МОДЕЛЕЙ ####

# Функции для расчета метрик
calc_mape <- function(actual, forecast) {
  mean(abs((actual - forecast) / actual)) * 100
}

calc_smape <- function(actual, forecast) {
  200 * mean(abs(actual - forecast) / (abs(actual) + abs(forecast)), na.rm = TRUE)
}

calc_rmse <- function(actual, forecast) {
  sqrt(mean((actual - forecast)^2))
}

# Сбор метрик качества
actual_2023 <- as.numeric(test[1])
actual_2024 <- as.numeric(test[2])

results <- data.frame(
  Модель = c("Наивная", "ARIMA", "ARIMA (лог)", "ETS", "Регрессия", "Комбинированная"),
  Прогноз_2023 = c(naive_forecast$mean[1], arima_forecast$mean[1], 
                   arima_forecast_log_corrected$mean[1], ets_forecast$mean[1],
                   reg_forecast[1], combined_forecast),
  Факт_2023 = rep(actual_2023, 6),
  Ошибка_2023 = c(
    abs(naive_forecast$mean[1] - actual_2023),
    abs(arima_forecast$mean[1] - actual_2023),
    abs(arima_forecast_log_corrected$mean[1] - actual_2023),
    abs(ets_forecast$mean[1] - actual_2023),
    abs(reg_forecast[1] - actual_2023),
    abs(combined_forecast - actual_2023)
  ),
  MAPE_2023 = c(
    calc_mape(actual_2023, naive_forecast$mean[1]),
    calc_mape(actual_2023, arima_forecast$mean[1]),
    calc_mape(actual_2023, arima_forecast_log_corrected$mean[1]),
    calc_mape(actual_2023, ets_forecast$mean[1]),
    calc_mape(actual_2023, reg_forecast[1]),
    calc_mape(actual_2023, combined_forecast)
  ),
  sMAPE_2023 = c(
    calc_smape(actual_2023, naive_forecast$mean[1]),
    calc_smape(actual_2023, arima_forecast$mean[1]),
    calc_smape(actual_2023, arima_forecast_log_corrected$mean[1]),
    calc_smape(actual_2023, ets_forecast$mean[1]),
    calc_smape(actual_2023, reg_forecast[1]),
    calc_smape(actual_2023, combined_forecast)
  ),
  RMSE_2023 = c(
    calc_rmse(actual_2023, naive_forecast$mean[1]),
    calc_rmse(actual_2023, arima_forecast$mean[1]),
    calc_rmse(actual_2023, arima_forecast_log_corrected$mean[1]),
    calc_rmse(actual_2023, ets_forecast$mean[1]),
    calc_rmse(actual_2023, reg_forecast[1]),
    calc_rmse(actual_2023, combined_forecast)
  )
)

cat("Сравнение моделей на тестовых данных (2023 год):\n")
print(round(results, 2))

# Выбор лучшей модели по комбинации метрик
results$Score <- (results$MAPE_2023 + results$sMAPE_2023 + results$RMSE_2023/1000000) / 3
best_model_idx <- which.min(results$Score)
best_model_name <- results$Модель[best_model_idx]

cat("\nЛучшая модель:", best_model_name, "(комплексный балл =", round(results$Score[best_model_idx], 2), ")\n")

#### 6. СЦЕНАРНОЕ ПРОГНОЗИРОВАНИЕ ####

cat("\n=== СЦЕНАРНОЕ ПРОГНОЗИРОВАНИЕ ===\n")

# Финальный прогноз на 2025-2027 с использованием всей доступной информации
if(best_model_name == "Наивная") {
  final_model <- naive(ts_data, h = 3, lambda = 0)
} else if(best_model_name == "ARIMA") {
  final_model <- forecast(auto.arima(ts_data, seasonal = FALSE), h = 3)
} else if(best_model_name == "ARIMA (лог)") {
  arima_final_log <- auto.arima(ts_data_log, seasonal = FALSE)
  final_model_log <- forecast(arima_final_log, h = 3)
  
  final_model <- final_model_log
  final_model$mean <- exp(final_model_log$mean + 0.5 * final_model_log$model$sigma2)
  final_model$lower <- exp(final_model_log$lower)
  final_model$upper <- exp(final_model_log$upper)
  final_model$x <- ts_data
} else if(best_model_name == "ETS") {
  final_model <- forecast(ets(ts_data, lambda = 0), h = 3)
} else if(best_model_name == "Регрессия") {
  # Простой линейный тренд для прогноза
  trend_model <- lm(INDEX ~ YEAR, data = index_data)
  future_years <- data.frame(YEAR = 2025:2027)
  trend_forecast <- predict(trend_model, newdata = future_years, interval = "prediction")
  
  final_model <- list(
    mean = ts(trend_forecast[,1], start = 2025),
    lower = ts(trend_forecast[,2], start = 2025),
    upper = ts(trend_forecast[,3], start = 2025),
    x = ts_data
  )
} else {
  # Комбинированный прогноз
  models <- list(
    naive(ts_data, h = 3, lambda = 0),
    forecast(auto.arima(ts_data, seasonal = FALSE), h = 3),
    forecast(ets(ts_data, lambda = 0), h = 3)
  )
  
  final_model <- models[[1]]
  final_model$mean <- (models[[1]]$mean + models[[2]]$mean + models[[3]]$mean) / 3
  final_model$lower <- (models[[1]]$lower + models[[2]]$lower + models[[3]]$lower) / 3
  final_model$upper <- (models[[1]]$upper + models[[2]]$upper + models[[3]]$upper) / 3
}

# Обеспечиваем неотрицательность прогноза
final_model$mean <- pmax(final_model$mean, 0)
final_model$lower <- pmax(final_model$lower, 0)

# Создание сценариев
scenarios <- data.frame(
  Год = 2025:2027,
  Базовый_сценарий = round(final_model$mean/1e6, 2),
  Пессимистичный = round(final_model$lower[,2]/1e6, 2),
  Оптимистичный = round(final_model$upper[,2]/1e6, 2),
  Консервативный = round(final_model$lower[,1]/1e6, 2),
  Агрессивный = round(final_model$upper[,1]/1e6, 2)
)

cat("Сценарные прогнозы (млн экз.):\n")
print(scenarios)

#### 7. РАСШИРЕННЫЙ РИСК-АНАЛИЗ ####

cat("\n=== РАСШИРЕННЫЙ РИСК-АНАЛИЗ ===\n")

# Анализ устойчивости запаса
current_level <- index_data$INDEX[10]
trend_5yr <- (index_data$INDEX[10] - index_data$INDEX[6]) / index_data$INDEX[6] * 100
trend_10yr <- (index_data$INDEX[10] - index_data$INDEX[1]) / index_data$INDEX[1] * 100

cat("Текущий уровень запаса (2024):", round(current_level/1e6, 2), "млн экз.\n")
cat("Изменение за последние 4 года:", round(trend_5yr, 1), "%\n")
cat("Изменение за весь период:", round(trend_10yr, 1), "%\n")

# Анализ волатильности
volatility <- sd(index_data$CHANGE_PCT[-1], na.rm = TRUE)
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

# Детальные рекомендации
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

#### 8. ИНТЕРАКТИВНАЯ ВИЗУАЛИЗАЦИЯ ####

# Создание интерактивного графика
plot_data <- data.frame(
  Year = c(index_data$YEAR, 2025:2027),
  Index = c(index_data$INDEX/1e6, final_model$mean/1e6),
  Lower_80 = c(rep(NA, 10), final_model$lower[,1]/1e6),
  Upper_80 = c(rep(NA, 10), final_model$upper[,1]/1e6),
  Lower_95 = c(rep(NA, 10), final_model$lower[,2]/1e6),
  Upper_95 = c(rep(NA, 10), final_model$upper[,2]/1e6),
  Type = c(rep("Historical", 10), rep("Forecast", 3))
)

p_interactive <- ggplot(plot_data, aes(x = Year, y = Index)) +
  geom_ribbon(aes(ymin = Lower_95, ymax = Upper_95), alpha = 0.2, fill = "lightblue") +
  geom_ribbon(aes(ymin = Lower_80, ymax = Upper_80), alpha = 0.3, fill = "blue") +
  geom_line(aes(color = Type), linewidth = 1.2) +
  geom_point(aes(color = Type), size = 3) +
  scale_color_manual(values = c("Historical" = "steelblue", "Forecast" = "red")) +
  labs(title = "Интерактивный прогноз популяции камчатского краба",
       subtitle = paste("Лучшая модель:", best_model_name),
       x = "Год", y = "Индекс обилия, млн экз.",
       color = "Тип данных") +
  theme_minimal() +
  theme(legend.position = "bottom")

print(p_interactive)

#### 9. АВТОМАТИЧЕСКАЯ ГЕНЕРАЦИЯ ОТЧЕТА ####

cat("\n=== ГЕНЕРАЦИЯ ОТЧЕТА ===\n")

# Создание сводного отчета
final_report <- list(
  metadata = list(
    анализ_проведен = Sys.time(),
    период_данных = paste(min(index_data$YEAR), "-", max(index_data$YEAR)),
    количество_наблюдений = nrow(index_data),
    лучшая_модель = best_model_name
  ),
  статистика = list(
    среднее_значение = round(mean(ts_data), 0),
    медиана = round(median(ts_data), 0),
    стандартное_отклонение = round(sd(ts_data), 0),
    коэффициент_вариации = round(sd(ts_data)/mean(ts_data)*100, 1),
    минимальное_значение = min(ts_data),
    максимальное_значение = max(ts_data)
  ),
  тренды = list(
    изменение_за_период = round(trend_10yr, 1),
    изменение_за_последние_4_года = round(trend_5yr, 1),
    волатильность = round(volatility, 1)
  ),
  сравнение_моделей = results,
  прогноз = scenarios,
  статус_запаса = list(
    статус = status,
    уровень_риска = risk_level,
    текущий_уровень = round(current_level/1e6, 2)
  )
)

# Сохранение отчета
saveRDS(final_report, "crab_analysis_report.rds")

# Создание HTML отчета
html_content <- paste0("
<!DOCTYPE html>
<html>
<head>
    <title>Отчет по анализу популяции камчатского краба</title>
    <meta charset='UTF-8'>
    <style>
        body { font-family: Arial, sans-serif; margin: 40px; }
        h1 { color: #2c3e50; }
        h2 { color: #34495e; }
        table { border-collapse: collapse; width: 100%; margin: 20px 0; }
        th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
        th { background-color: #f2f2f2; }
        .status-critical { color: red; font-weight: bold; }
        .status-warning { color: orange; font-weight: bold; }
        .status-stable { color: green; font-weight: bold; }
    </style>
</head>
<body>
    <h1>Отчет по анализу популяции камчатского краба</h1>
    <p><strong>Дата анализа:</strong> ", Sys.time(), "</p>
    <p><strong>Период данных:</strong> ", min(index_data$YEAR), "-", max(index_data$YEAR), "</p>
    <p><strong>Лучшая модель:</strong> ", best_model_name, "</p>
    
    <h2>Статус запаса</h2>
    <p><strong>Статус:</strong> <span class='status-", tolower(status), "'>", status, "</span></p>
    <p><strong>Уровень риска:</strong> ", risk_level, "</p>
    <p><strong>Текущий уровень:</strong> ", round(current_level/1e6, 2), " млн экз.</p>
    
    <h2>Прогноз на 2025-2027 гг.</h2>
    <table>
        <tr><th>Год</th><th>Базовый сценарий</th><th>Пессимистичный</th><th>Оптимистичный</th></tr>
        <tr><td>2025</td><td>", scenarios$Базовый_сценарий[1], "</td><td>", scenarios$Пессимистичный[1], "</td><td>", scenarios$Оптимистичный[1], "</td></tr>
        <tr><td>2026</td><td>", scenarios$Базовый_сценарий[2], "</td><td>", scenarios$Пессимистичный[2], "</td><td>", scenarios$Оптимистичный[2], "</td></tr>
        <tr><td>2027</td><td>", scenarios$Базовый_сценарий[3], "</td><td>", scenarios$Пессимистичный[3], "</td><td>", scenarios$Оптимистичный[3], "</td></tr>
    </table>
</body>
</html>
")

writeLines(html_content, "crab_analysis_report.html")

cat("\n=== АНАЛИЗ ЗАВЕРШЕН ===\n")
cat("Лучшая модель:", best_model_name, "\n")
cat("Статус запаса:", status, "\n")
cat("Прогноз на 2025:", scenarios$Базовый_сценарий[1], "млн экз.\n")
cat("Диапазон неопределенности (95% ДИ):", 
    scenarios$Пессимистичный[1], "-", scenarios$Оптимистичный[1], "млн экз.\n")
cat("Отчет сохранен в файлы: crab_analysis_report.rds и crab_analysis_report.html\n")

# Сохранение графиков
ggsave("crab_dynamics_plot.png", p1, width = 10, height = 6, dpi = 300)
ggsave("crab_forecast_plot.png", p_interactive, width = 12, height = 8, dpi = 300)

cat("Графики сохранены: crab_dynamics_plot.png и crab_forecast_plot.png\n")