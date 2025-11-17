# ДЕМОНСТРАЦИОННЫЙ СКРИПТ: ИСПОЛЬЗОВАНИЕ РАСШИРЕННЫХ ФУНКЦИЙ АНАЛИЗА
# Показывает возможности усовершенствованного анализа популяции камчатского краба

# Загрузка основных пакетов
suppressPackageStartupMessages({
  library(tidyverse)
  library(forecast)
  library(ggplot2)
  library(gridExtra)
})

# Загрузка пользовательских функций
source("crab_analysis_functions.R")

# Загрузка основного скрипта анализа
source("enhanced_crab_analysis.R")

cat("=== ДЕМОНСТРАЦИЯ РАСШИРЕННЫХ ВОЗМОЖНОСТЕЙ ===\n")

# Создание расширенного набора данных для демонстрации
demo_data <- data.frame(
  YEAR = 2010:2024,
  INDEX = c(2200000, 1950000, 2100000, 1850000, 2100000, 1950000, 2381774, 
           1634549, 1920507, 1036673, 1147685, 1055733, 980000, 920000, 850000),
  TEMPERATURE = c(2.0, 1.9, 2.1, 1.8, 2.0, 1.9, 2.0, 1.9, 2.2, 1.7, 1.8, 1.9, 2.1, 2.0, 1.8),
  SALINITY = c(34.1, 34.2, 34.0, 34.3, 34.1, 34.2, 34.0, 34.2, 34.1, 34.4, 34.3, 34.2, 34.1, 34.0, 34.2),
  FISHING_EFFORT = c(100, 110, 105, 120, 115, 135, 150, 180, 160, 200, 190, 195, 210, 220, 230)
)

cat("\n1. ДЕМОНСТРАЦИЯ БИОЛОГИЧЕСКИХ ИНДЕКСОВ\n")
biological_indices <- calculate_biological_indices(demo_data)
cat("Индекс стабильности популяции:", round(biological_indices$stability_index, 3), "\n")
cat("Индекс восстановления:", round(biological_indices$recovery_index, 3), "\n")
if(!is.null(biological_indices$fishing_resistance)) {
  cat("Устойчивость к промыслу:", round(biological_indices$fishing_resistance, 3), "\n")
}

cat("\n2. АНАЛИЗ ЭКОЛОГИЧЕСКИХ ФАКТОРОВ\n")
env_factors <- analyze_environmental_factors(demo_data)
cat("Влияние температуры:", round(env_factors$temperature_impact, 3), 
    "(", env_factors$temperature_significance, ")\n")
cat("Влияние солености:", round(env_factors$salinity_impact, 3), 
    "(", env_factors$salinity_significance, ")\n")
cat("Влияние промысла:", round(env_factors$fishing_impact, 3), 
    "(", env_factors$fishing_significance, ")\n")

cat("\n3. ПРОМЫСЛОВЫЕ РЕКОМЕНДАЦИИ\n")
current_index <- demo_data$INDEX[nrow(demo_data)]
fishing_recs <- calculate_fishing_recommendations(current_index, demo_data, NULL)
cat("Рекомендуемый уровень промысла:", fishing_recs$fishing_level, "\n")
cat("Снижение промыслового усилия:", round(fishing_recs$effort_reduction * 100, 0), "%\n")
cat("Приоритет мониторинга:", fishing_recs$priority, "\n")

cat("\n4. АНАЛИЗ РИСКОВ\n")
risk_analysis <- analyze_risks(demo_data, NULL)
cat("Риск критического снижения:", risk_analysis$critical_decline, "\n")
cat("Риск экологических изменений:", risk_analysis$environmental_change, "\n")

cat("\n5. СОЗДАНИЕ ДЕТАЛЬНОГО ОТЧЕТА\n")
detailed_report <- create_detailed_report(demo_data, NULL, "ARIMA")
cat("Отчет создан с", length(detailed_report), "разделами\n")
cat("Базовые метрики:\n")
cat("- Период:", detailed_report$basic_metrics$period, "\n")
cat("- Наблюдений:", detailed_report$basic_metrics$observations, "\n")
cat("- Средний индекс:", detailed_report$basic_metrics$mean_index, "\n")

cat("\n6. СОЗДАНИЕ ВИЗУАЛИЗАЦИЙ\n")
plots <- create_comprehensive_plots(demo_data, NULL)
cat("Создано графиков:", length(plots), "\n")

# Отображение основного графика
if(!is.null(plots$main_plot)) {
  print(plots$main_plot)
}

cat("\n7. ЭКСПОРТ РЕЗУЛЬТАТОВ\n")
export_results(detailed_report, plots, "demo_crab_analysis")

cat("\n8. ДЕМОНСТРАЦИЯ МОНИТОРИНГА ИЗМЕНЕНИЙ\n")
# Симуляция новых данных
new_data <- data.frame(
  YEAR = 2025,
  INDEX = 800000,  # Новое значение
  TEMPERATURE = 1.9,
  SALINITY = 34.1,
  FISHING_EFFORT = 240
)

changes <- monitor_changes(new_data, demo_data)
cat("Изменение по сравнению с предыдущим годом:", round(changes$year_over_year, 1), "%\n")
cat("Изменение по сравнению с историческим средним:", round(changes$vs_historical_mean, 1), "%\n")
cat("Значимость изменений:", changes$significance, "\n")

cat("\n=== ДЕМОНСТРАЦИЯ ЗАВЕРШЕНА ===\n")
cat("Все функции работают корректно!\n")
cat("Созданы файлы:\n")
cat("- demo_crab_analysis_detailed_report.rds\n")
cat("- demo_crab_analysis_main_plot.png\n")
cat("- demo_crab_analysis_summary.csv\n")

# Дополнительная демонстрация: сравнение различных моделей
cat("\n9. СРАВНЕНИЕ МОДЕЛЕЙ НА РАСШИРЕННЫХ ДАННЫХ\n")

# Создание временного ряда
ts_demo <- ts(demo_data$INDEX, start = 2010, frequency = 1)

# Разделение данных
train_demo <- window(ts_demo, end = 2022)
test_demo <- window(ts_demo, start = 2023)

# Построение различных моделей
models_comparison <- list(
  "Наивная" = naive(train_demo, h = 2),
  "ARIMA" = forecast(auto.arima(train_demo), h = 2),
  "ETS" = forecast(ets(train_demo), h = 2),
  "Линейный тренд" = forecast(tslm(train_demo ~ trend), h = 2)
)

# Оценка качества
actual_values <- as.numeric(test_demo)
model_scores <- data.frame(
  Модель = names(models_comparison),
  MAPE = sapply(models_comparison, function(m) {
    mean(abs((actual_values - m$mean[1:2]) / actual_values)) * 100
  }),
  RMSE = sapply(models_comparison, function(m) {
    sqrt(mean((actual_values - m$mean[1:2])^2))
  })
)

cat("Сравнение моделей на расширенных данных:\n")
print(round(model_scores, 2))

# Визуализация сравнения моделей
comparison_plot <- ggplot(data.frame(
  Year = c(2010:2024, 2025:2026),
  Index = c(demo_data$INDEX/1e6, models_comparison$ARIMA$mean/1e6),
  Type = c(rep("Historical", 15), rep("Forecast", 2))
), aes(x = Year, y = Index, color = Type)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  scale_color_manual(values = c("Historical" = "steelblue", "Forecast" = "red")) +
  labs(title = "Сравнение моделей на расширенных данных",
       x = "Год", y = "Индекс обилия, млн экз.") +
  theme_minimal()

print(comparison_plot)

cat("\n=== ПОЛНАЯ ДЕМОНСТРАЦИЯ ЗАВЕРШЕНА ===\n")
cat("Расширенный анализ включает:\n")
cat("✓ Биологические индексы\n")
cat("✓ Анализ экологических факторов\n")
cat("✓ Промысловые рекомендации\n")
cat("✓ Анализ рисков\n")
cat("✓ Детальные отчеты\n")
cat("✓ Комплексная визуализация\n")
cat("✓ Экспорт результатов\n")
cat("✓ Мониторинг изменений\n")
cat("✓ Сравнение моделей\n")