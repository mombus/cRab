# ДОПОЛНИТЕЛЬНЫЕ ФУНКЦИИ ДЛЯ АНАЛИЗА ПОПУЛЯЦИИ КАМЧАТСКОГО КРАБА
# Библиотека специализированных функций

# Функция для расчета биологических индексов
calculate_biological_indices <- function(data) {
  indices <- list()
  
  # Индекс стабильности популяции
  indices$stability_index <- 1 - (sd(data$INDEX) / mean(data$INDEX))
  
  # Индекс восстановления (способность к восстановлению после снижения)
  declines <- which(diff(data$INDEX) < 0)
  if(length(declines) > 0) {
    recovery_rates <- numeric()
    for(i in declines) {
      if(i < length(data$INDEX) - 1) {
        recovery_rate <- (data$INDEX[i+2] - data$INDEX[i+1]) / data$INDEX[i+1]
        recovery_rates <- c(recovery_rates, recovery_rate)
      }
    }
    indices$recovery_index <- mean(recovery_rates, na.rm = TRUE)
  } else {
    indices$recovery_index <- 0
  }
  
  # Индекс устойчивости к промыслу
  if("FISHING_EFFORT" %in% names(data)) {
    fishing_correlation <- cor(data$INDEX, data$FISHING_EFFORT, use = "complete.obs")
    indices$fishing_resistance <- 1 - abs(fishing_correlation)
  }
  
  return(indices)
}

# Функция для анализа экологических факторов
analyze_environmental_factors <- function(data) {
  factors <- list()
  
  # Анализ влияния температуры
  if("TEMPERATURE" %in% names(data)) {
    temp_cor <- cor(data$INDEX, data$TEMPERATURE, use = "complete.obs")
    factors$temperature_impact <- temp_cor
    factors$temperature_significance <- ifelse(abs(temp_cor) > 0.5, "Высокое", 
                                             ifelse(abs(temp_cor) > 0.3, "Умеренное", "Низкое"))
  }
  
  # Анализ влияния солености
  if("SALINITY" %in% names(data)) {
    sal_cor <- cor(data$INDEX, data$SALINITY, use = "complete.obs")
    factors$salinity_impact <- sal_cor
    factors$salinity_significance <- ifelse(abs(sal_cor) > 0.5, "Высокое", 
                                          ifelse(abs(sal_cor) > 0.3, "Умеренное", "Низкое"))
  }
  
  # Анализ влияния промыслового усилия
  if("FISHING_EFFORT" %in% names(data)) {
    fish_cor <- cor(data$INDEX, data$FISHING_EFFORT, use = "complete.obs")
    factors$fishing_impact <- fish_cor
    factors$fishing_significance <- ifelse(abs(fish_cor) > 0.5, "Высокое", 
                                         ifelse(abs(fish_cor) > 0.3, "Умеренное", "Низкое"))
  }
  
  return(factors)
}

# Функция для расчета промысловых рекомендаций
calculate_fishing_recommendations <- function(current_index, historical_data, forecast_data) {
  recommendations <- list()
  
  # Расчет оптимального промыслового усилия
  historical_mean <- mean(historical_data$INDEX)
  current_ratio <- current_index / historical_mean
  
  if(current_ratio < 0.5) {
    recommendations$fishing_level <- "МИНИМАЛЬНЫЙ"
    recommendations$effort_reduction <- 0.8
    recommendations$priority <- "ВЫСОКИЙ"
  } else if(current_ratio < 0.7) {
    recommendations$fishing_level <- "ОГРАНИЧЕННЫЙ"
    recommendations$effort_reduction <- 0.5
    recommendations$priority <- "ВЫСОКИЙ"
  } else if(current_ratio < 0.9) {
    recommendations$fishing_level <- "УМЕРЕННЫЙ"
    recommendations$effort_reduction <- 0.25
    recommendations$priority <- "СРЕДНИЙ"
  } else {
    recommendations$fishing_level <- "СТАНДАРТНЫЙ"
    recommendations$effort_reduction <- 0
    recommendations$priority <- "НИЗКИЙ"
  }
  
  # Расчет рекомендуемого улова на основе прогноза
  if(!is.null(forecast_data)) {
    forecast_2025 <- forecast_data$mean[1]
    recommendations$recommended_catch_2025 <- forecast_2025 * 0.1  # 10% от прогнозируемого запаса
  }
  
  return(recommendations)
}

# Функция для анализа рисков
analyze_risks <- function(data, forecast_data) {
  risks <- list()
  
  # Риск критического снижения
  min_historical <- min(data$INDEX)
  current_level <- data$INDEX[nrow(data)]
  decline_risk <- (current_level - min_historical) / min_historical
  
  if(decline_risk < 0.2) {
    risks$critical_decline <- "ВЫСОКИЙ"
  } else if(decline_risk < 0.4) {
    risks$critical_decline <- "УМЕРЕННЫЙ"
  } else {
    risks$critical_decline <- "НИЗКИЙ"
  }
  
  # Риск неопределенности прогноза
  if(!is.null(forecast_data)) {
    forecast_range <- forecast_data$upper[1,2] - forecast_data$lower[1,2]
    forecast_cv <- forecast_range / forecast_data$mean[1]
    
    if(forecast_cv > 0.5) {
      risks$forecast_uncertainty <- "ВЫСОКИЙ"
    } else if(forecast_cv > 0.3) {
      risks$forecast_uncertainty <- "УМЕРЕННЫЙ"
    } else {
      risks$forecast_uncertainty <- "НИЗКИЙ"
    }
  }
  
  # Риск экологических изменений
  if("TEMPERATURE" %in% names(data)) {
    temp_trend <- cor(data$YEAR, data$TEMPERATURE, use = "complete.obs")
    if(abs(temp_trend) > 0.5) {
      risks$environmental_change <- "ВЫСОКИЙ"
    } else if(abs(temp_trend) > 0.3) {
      risks$environmental_change <- "УМЕРЕННЫЙ"
    } else {
      risks$environmental_change <- "НИЗКИЙ"
    }
  }
  
  return(risks)
}

# Функция для создания детального отчета
create_detailed_report <- function(data, forecast_data, best_model_name) {
  report <- list()
  
  # Базовые метрики
  report$basic_metrics <- list(
    period = paste(min(data$YEAR), "-", max(data$YEAR)),
    observations = nrow(data),
    mean_index = round(mean(data$INDEX), 0),
    median_index = round(median(data$INDEX), 0),
    cv = round(sd(data$INDEX)/mean(data$INDEX)*100, 1)
  )
  
  # Биологические индексы
  report$biological_indices <- calculate_biological_indices(data)
  
  # Экологические факторы
  report$environmental_factors <- analyze_environmental_factors(data)
  
  # Промысловые рекомендации
  report$fishing_recommendations <- calculate_fishing_recommendations(
    data$INDEX[nrow(data)], data, forecast_data
  )
  
  # Анализ рисков
  report$risk_analysis <- analyze_risks(data, forecast_data)
  
  # Прогноз
  if(!is.null(forecast_data)) {
    report$forecast <- list(
      best_model = best_model_name,
      forecast_2025 = round(forecast_data$mean[1]/1e6, 2),
      forecast_2026 = round(forecast_data$mean[2]/1e6, 2),
      uncertainty_2025 = round((forecast_data$upper[1,2] - forecast_data$lower[1,2])/1e6, 2)
    )
  }
  
  return(report)
}

# Функция для визуализации результатов
create_comprehensive_plots <- function(data, forecast_data) {
  plots <- list()
  
  # График динамики с прогнозом
  plot_data <- data.frame(
    Year = c(data$YEAR, 2025:2027),
    Index = c(data$INDEX/1e6, forecast_data$mean/1e6),
    Lower = c(rep(NA, nrow(data)), forecast_data$lower[,2]/1e6),
    Upper = c(rep(NA, nrow(data)), forecast_data$upper[,2]/1e6),
    Type = c(rep("Historical", nrow(data)), rep("Forecast", 3))
  )
  
  plots$main_plot <- ggplot(plot_data, aes(x = Year, y = Index)) +
    geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.2, fill = "lightblue") +
    geom_line(aes(color = Type), linewidth = 1.2) +
    geom_point(aes(color = Type), size = 3) +
    scale_color_manual(values = c("Historical" = "steelblue", "Forecast" = "red")) +
    labs(title = "Динамика и прогноз популяции камчатского краба",
         x = "Год", y = "Индекс обилия, млн экз.") +
    theme_minimal()
  
  # График экологических факторов
  if("TEMPERATURE" %in% names(data) && "FISHING_EFFORT" %in% names(data)) {
    env_data <- data.frame(
      Year = data$YEAR,
      Index = data$INDEX/1e6,
      Temperature = data$TEMPERATURE,
      Fishing = data$FISHING_EFFORT/100
    )
    
    plots$environmental_plot <- ggplot(env_data, aes(x = Year)) +
      geom_line(aes(y = Index, color = "Индекс обилия"), linewidth = 1) +
      geom_line(aes(y = Temperature, color = "Температура"), linewidth = 1) +
      geom_line(aes(y = Fishing, color = "Промысловое усилие"), linewidth = 1) +
      scale_y_continuous(
        name = "Индекс обилия, млн экз.",
        sec.axis = sec_axis(~., name = "Температура / Промысловое усилие")
      ) +
      scale_color_manual(values = c("Индекс обилия" = "steelblue", 
                                   "Температура" = "red", 
                                   "Промысловое усилие" = "green")) +
      labs(title = "Влияние экологических факторов",
           x = "Год") +
      theme_minimal() +
      theme(legend.position = "bottom")
  }
  
  return(plots)
}

# Функция для экспорта результатов в различные форматы
export_results <- function(report, plots, filename_prefix = "crab_analysis") {
  # Сохранение отчета в RDS
  saveRDS(report, paste0(filename_prefix, "_detailed_report.rds"))
  
  # Сохранение графиков
  if(!is.null(plots$main_plot)) {
    ggsave(paste0(filename_prefix, "_main_plot.png"), plots$main_plot, 
           width = 12, height = 8, dpi = 300)
  }
  
  if(!is.null(plots$environmental_plot)) {
    ggsave(paste0(filename_prefix, "_environmental_plot.png"), plots$environmental_plot, 
           width = 12, height = 8, dpi = 300)
  }
  
  # Создание CSV с основными результатами
  summary_data <- data.frame(
    Метрика = c("Период анализа", "Количество наблюдений", "Средний индекс", 
                "Коэффициент вариации", "Лучшая модель", "Прогноз 2025"),
    Значение = c(report$basic_metrics$period, 
                report$basic_metrics$observations,
                report$basic_metrics$mean_index,
                paste0(report$basic_metrics$cv, "%"),
                report$forecast$best_model,
                paste0(report$forecast$forecast_2025, " млн экз."))
  )
  
  write.csv(summary_data, paste0(filename_prefix, "_summary.csv"), 
            row.names = FALSE, fileEncoding = "UTF-8")
  
  cat("Результаты экспортированы:\n")
  cat("- ", filename_prefix, "_detailed_report.rds\n")
  cat("- ", filename_prefix, "_main_plot.png\n")
  if(!is.null(plots$environmental_plot)) {
    cat("- ", filename_prefix, "_environmental_plot.png\n")
  }
  cat("- ", filename_prefix, "_summary.csv\n")
}

# Функция для мониторинга изменений
monitor_changes <- function(new_data, historical_data) {
  changes <- list()
  
  # Сравнение с предыдущим годом
  if(nrow(new_data) > 0 && nrow(historical_data) > 0) {
    last_historical <- historical_data$INDEX[nrow(historical_data)]
    first_new <- new_data$INDEX[1]
    
    changes$year_over_year <- (first_new - last_historical) / last_historical * 100
    
    # Сравнение с историческим средним
    historical_mean <- mean(historical_data$INDEX)
    changes$vs_historical_mean <- (first_new - historical_mean) / historical_mean * 100
    
    # Оценка значимости изменений
    if(abs(changes$year_over_year) > 20) {
      changes$significance <- "КРИТИЧЕСКОЕ"
    } else if(abs(changes$year_over_year) > 10) {
      changes$significance <- "ЗНАЧИТЕЛЬНОЕ"
    } else if(abs(changes$year_over_year) > 5) {
      changes$significance <- "УМЕРЕННОЕ"
    } else {
      changes$significance <- "НЕЗНАЧИТЕЛЬНОЕ"
    }
  }
  
  return(changes)
}