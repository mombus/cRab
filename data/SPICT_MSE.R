# ===============================================================
#     СКРИПТ 3: MANAGEMENT STRATEGY EVALUATION (MSE)- ОЦЕНКА СТРАТЕГИИ УПРАВЛЕНИЯ
#     Сравнение трех ключевых стратегий управления
#     Курс: Оценка водных биоресурсов при недостатке данных в R
#     Автор: Баканёв С.В.
#     Дата создания: 28.08.2025
# ===============================================================

# ======================= ВВЕДЕНИЕ =============================
# В этом скрипте сравниваются три основные стратегии управления:
# 1. Fish at Fmsy - промысел на уровне оптимальной смертности
# 2. MSY hockey-stick rule - адаптивное правило с защитой запаса
# 3. ICES advice rule - предосторожный подход ICES
#
# MSE позволяет тестировать эти стратегии с учетом всех
# источников неопределенности в системе промысел-запас

# ------------------- 1. ПОДГОТОВКА СРЕДЫ --------------------

## 1.1 Очистка рабочей среды и настройка
rm(list = ls())
set.seed(123)  # Для воспроизводимости результатов

## 1.2 Загрузка необходимых библиотек
library(spict)       # Для работы с моделью SPiCT
library(tidyverse)   # Обработка данных
library(ggplot2)     # Визуализация
library(viridis)     # Цветовые схемы
library(patchwork)   # Компоновка графиков
library(gridExtra)   # Дополнительные возможности компоновки

## 1.3 Установка рабочей директории и загрузка модели
setwd("C:/SPICT")

# Загружаем подогнанную модель из первого скрипта
if (file.exists("spict_model_fit.rds")) {
  fit <- readRDS("spict_model_fit.rds")
  cat("\n✓ Модель SPiCT успешно загружена\n")
} else {
  stop("Файл модели не найден. Запустите первый скрипт для подгонки модели.")
}

cat("\n" , strrep("=", 60), "\n")
cat("      MANAGEMENT STRATEGY EVALUATION (MSE)\n")
cat("   Сравнение трех стратегий управления промыслом\n")
cat(strrep("=", 60), "\n")

# ------------------- 2. ИЗВЛЕЧЕНИЕ ПАРАМЕТРОВ ИЗ МОДЕЛИ --------------------

cat("\n========== ПАРАМЕТРЫ ОПЕРАЦИОННОЙ МОДЕЛИ ==========\n")

## 2.1 Извлечение оценок параметров популяционной динамики
# Эти параметры представляют "истинное" состояние в операционной модели
r_true <- get.par("logr", fit, exp = TRUE)[1]          # Внутренний темп роста
K_true <- get.par("logK", fit, exp = TRUE)[1]          # Ёмкость среды
Bmsy_true <- get.par("logBmsy", fit, exp = TRUE)[1]    # Биомасса MSY
Fmsy_true <- get.par("logFmsy", fit, exp = TRUE)[1]    # Промысловая смертность MSY
B_current <- get.par("logB", fit, exp = TRUE)[1]       # Текущая биомасса
F_current <- get.par("logF", fit, exp = TRUE)[1]       # Текущая пром. смертность

# Параметр улавливаемости для индекса CPUE
q_cpue <- get.par("logq", fit, exp = TRUE)[1]

## 2.2 Вывод ключевых параметров
cat("\nИзвлеченные параметры модели:\n")
cat(sprintf("r (темп роста): %.3f год⁻¹\n", r_true))
cat(sprintf("K (ёмкость средыь): %.1f тыс. т\n", K_true))
cat(sprintf("Bmsy: %.1f тыс. т\n", Bmsy_true))
cat(sprintf("Fmsy: %.3f год⁻¹\n", Fmsy_true))
cat(sprintf("Текущая биомасса: %.1f тыс. т\n", B_current))
cat(sprintf("Текущее B/Bmsy: %.2f\n", B_current/Bmsy_true))
cat(sprintf("Текущее F/Fmsy: %.2f\n", F_current/Fmsy_true))

# ------------------- 3. НАСТРОЙКИ СИМУЛЯЦИИ --------------------

cat("\n========== НАСТРОЙКИ MSE ==========\n")

## 3.1 Основные параметры симуляции
n_sim <- 500              # Количество симуляций (реализаций)
n_years <- 100             # Период прогнозирования (лет)
assessment_interval <- 2  # Интервал между оценками запаса (лет)

cat(sprintf("Количество симуляций: %d\n", n_sim))
cat(sprintf("Горизонт прогнозирования: %d лет\n", n_years))
cat(sprintf("Частота оценки запаса: каждые %d года\n", assessment_interval))

## 3.2 Параметры неопределенности
# Все источники неопределенности в системе

# Процессная ошибка - естественная изменчивость в динамике популяции
process_error_cv <- 0.15  

# Ошибка наблюдения - неточность в индексах биомассы
observation_error_cv <- 0.20  

# Ошибка оценки - неточность в оценке состояния запаса
assessment_bias_cv <- 0.15  

# Ошибка реализации - разница между рекомендованным и фактическим выловом
implementation_error_cv <- 0.10  

cat("\nПараметры неопределенности:\n")
cat(sprintf("CV процессной ошибки: %.0f%%\n", process_error_cv * 100))
cat(sprintf("CV ошибки наблюдения: %.0f%%\n", observation_error_cv * 100))
cat(sprintf("CV смещения в оценках: %.0f%%\n", assessment_bias_cv * 100))
cat(sprintf("CV ошибки реализации: %.0f%%\n", implementation_error_cv * 100))

# ------------------- 4. ОПЕРАЦИОННАЯ МОДЕЛЬ --------------------

cat("\n========== ОПРЕДЕЛЕНИЕ ОПЕРАЦИОННОЙ МОДЕЛИ ==========\n")

## 4.1 Функция истинной динамики популяции
# Модель Шефера с процессной стохастичностью
simulate_population_dynamics <- function(B_t, F_t, r, K, process_cv = 0.15) {
  
  # Расчет прибавочной продукции (модель Шефера)
  surplus_production <- r * B_t * (1 - B_t/K)
  
  # Расчет вылова
  catch <- F_t * B_t
  
  # Обновление биомассы с учетом процессной ошибки
  # Используем лог-нормальное распределение для мультипликативной ошибки
  process_error <- rlnorm(1, meanlog = -process_cv^2/2, sdlog = process_cv)
  
  B_next <- (B_t + surplus_production - catch) * process_error
  
  # Ограничения для реалистичности
  B_next <- max(B_next, 0.001 * K)  # Минимум 0.1% от K
  B_next <- min(B_next, 1.5 * K)    # Максимум 150% от K
  
  return(list(
    B_next = B_next,
    catch_realized = catch,
    surplus = surplus_production,
    process_multiplier = process_error
  ))
}

## 4.2 Функция генерации наблюдаемого индекса
# Имитирует процесс сбора данных с ошибками
generate_index_observation <- function(B_true, q, obs_cv = 0.20) {
  
  # Истинный индекс пропорционален биомассе
  true_index <- q * B_true
  
  # Добавление ошибки наблюдения
  obs_error <- rlnorm(1, meanlog = -obs_cv^2/2, sdlog = obs_cv)
  observed_index <- true_index * obs_error
  
  return(observed_index)
}

## 4.3 Функция оценки состояния запаса
# Упрощенная процедура оценки (в реальности здесь бы запускалась полная модель)
assess_stock_status <- function(index_history, catch_history, 
                               true_q, true_Bmsy, true_Fmsy,
                               assessment_cv = 0.15) {
  
  # Оценка текущей биомассы по последним наблюдениям индекса
  # Используем среднее за последние 3 года для сглаживания
  recent_indices <- tail(index_history[!is.na(index_history)], 3)
  mean_recent_index <- mean(recent_indices, na.rm = TRUE)
  
  # Оценка биомассы с учетом смещения
  assessment_bias <- rlnorm(1, meanlog = -assessment_cv^2/2, sdlog = assessment_cv)
  B_estimated <- (mean_recent_index / true_q) * assessment_bias
  
  # Оценка текущей F из последнего вылова
  recent_catch <- tail(catch_history[catch_history > 0], 1)
  F_estimated <- ifelse(length(recent_catch) > 0 && B_estimated > 0,
                        recent_catch / B_estimated,
                        0.1)  # Значение по умолчанию
  
  # Оценки референсных точек также содержат неопределенность
  Bmsy_estimated <- true_Bmsy * rlnorm(1, 0, assessment_cv/2)
  Fmsy_estimated <- true_Fmsy * rlnorm(1, 0, assessment_cv/2)
  
  return(list(
    B = B_estimated,
    F = F_estimated,
    Bmsy = Bmsy_estimated,
    Fmsy = Fmsy_estimated,
    B_Bmsy = B_estimated / Bmsy_estimated,
    F_Fmsy = F_estimated / Fmsy_estimated
  ))
}

# ------------------- 5. ПРАВИЛА УПРАВЛЕНИЯ (HCR) --------------------

cat("\n========== ОПРЕДЕЛЕНИЕ ПРАВИЛ УПРАВЛЕНИЯ ==========\n")

## 5.1 Стратегия 1: Fish at Fmsy
# Простейшее правило - всегда промысел на уровне Fmsy
HCR_Fmsy <- function(assessment, previous_TAC = NULL) {
  
  # Рекомендация по F
  F_advice <- assessment$Fmsy
  
  # Расчет TAC
  TAC <- F_advice * assessment$B
  
  # Обеспечиваем неотрицательность
  TAC <- max(TAC, 0)
  
  return(list(
    F_advice = F_advice,
    TAC = TAC,
    rule_name = "Fish at Fmsy",
    status = ifelse(assessment$B_Bmsy < 0.5, "Риск истощения", "Стандартный промысел")
  ))
}

cat("✓ Стратегия 1: Fish at Fmsy\n")
cat("  - Постоянный промысел на уровне Fmsy\n")
cat("  - Не учитывает состояние запаса\n")
cat("  - Простое в применении правило\n")

## 5.2 Стратегия 2: MSY Hockey-stick Rule
# Адаптивное правило с защитой при низкой биомассе
HCR_hockey_stick <- function(assessment, previous_TAC = NULL,
                            Blim_fraction = 0.5) {
  
  # Определение предельной биомассы
  Blim <- assessment$Bmsy * Blim_fraction
  
  # Применение правила hockey-stick
  if (assessment$B <= Blim) {
    # Линейное снижение F при B < Blim
    F_multiplier <- assessment$B / Blim
    F_advice <- assessment$Fmsy * F_multiplier
    status <- "Снижение промысла"
  } else {
    # Полный промысел при B >= Blim
    F_advice <- assessment$Fmsy
    status <- "Полный промысел"
  }
  
  # Расчет TAC
  TAC <- F_advice * assessment$B
  TAC <- max(TAC, 0)
  
  # Полное закрытие при критически низкой биомассе
  if (assessment$B_Bmsy < 0.2) {
    TAC <- 0
    F_advice <- 0
    status <- "Промысел закрыт"
  }
  
  return(list(
    F_advice = F_advice,
    TAC = TAC,
    rule_name = "MSY Hockey-stick",
    status = status,
    Blim = Blim
  ))
}

cat("\n✓ Стратегия 2: MSY Hockey-stick Rule\n")
cat("  - Адаптивное управление в зависимости от состояния\n")
cat("  - Линейное снижение F при B < 0.5*Bmsy\n")
cat("  - Полное закрытие при B < 0.2*Bmsy\n")

## 5.3 Стратегия 3: ICES Advice Rule
# Предосторожный подход с ограничениями межгодовых изменений
HCR_ICES <- function(assessment, previous_TAC = NULL,
                    Bpa_multiplier = 1.4,
                    Fpa_multiplier = 0.85,
                    max_TAC_change = 0.20) {
  
  # Расчет предосторожных ориентиров управления
  Bpa <- assessment$Bmsy / Bpa_multiplier   # Предосторожная биомасса
  Blim <- Bpa / 1.4                         # Предельная биомасса  
  Fpa <- assessment$Fmsy * Fpa_multiplier   # Предосторожная F
  
  # Определение F по правилу ICES
  if (assessment$B < Blim) {
    # Критическое состояние - закрытие промысла
    F_advice <- 0
    status <- "Критическое - промысел закрыт"
    
  } else if (assessment$B >= Blim && assessment$B < Bpa) {
    # Восстановление - пропорциональное снижение F
    F_multiplier <- (assessment$B - Blim) / (Bpa - Blim)
    F_advice <- Fpa * F_multiplier
    status <- "Восстановление запаса"
    
  } else {
    # Нормальное состояние - предосторожный промысел
    F_advice <- Fpa
    status <- "Устойчивый промысел"
  }
  
  # Расчет TAC
  TAC <- F_advice * assessment$B
  TAC <- max(TAC, 0)
  
  # Ограничение межгодовых изменений TAC (стабильность для промышленности)
  if (!is.null(previous_TAC) && previous_TAC > 0 && TAC > 0) {
    max_increase <- previous_TAC * (1 + max_TAC_change)
    max_decrease <- previous_TAC * (1 - max_TAC_change)
    TAC_constrained <- min(max(TAC, max_decrease), max_increase)
    
    if (TAC != TAC_constrained) {
      status <- paste(status, "(TAC ограничен)")
    }
    TAC <- TAC_constrained
  }
  
  return(list(
    F_advice = F_advice,
    TAC = TAC,
    rule_name = "ICES Advice Rule",
    status = status,
    Bpa = Bpa,
    Blim = Blim,
    Fpa = Fpa
  ))
}

cat("\n✓ Стратегия 3: ICES Advice Rule\n")
cat("  - Предосторожный подход (Fpa = 0.85*Fmsy)\n")
cat("  - Двухуровневая система пороговых значений\n")
cat("  - Ограничение межгодовых изменений TAC (±20%)\n")

# ------------------- 6. ФУНКЦИЯ MSE ДЛЯ ОДНОЙ СИМУЛЯЦИИ --------------------

## 6.1 Основная функция симуляции
run_mse_simulation <- function(sim_id, HCR_function, params, settings) {
  
  # Инициализация с небольшой вариацией начальных условий
  B_sim <- params$B_initial * rlnorm(1, 0, 0.05)
  
  # Вариация в "истинных" параметрах (представляет неопределенность в природе)
  r_sim <- params$r * rlnorm(1, 0, 0.05)
  K_sim <- params$K * rlnorm(1, 0, 0.05)
  
  # Массивы для хранения результатов
  n_years <- settings$n_years
  results <- data.frame(
    sim_id = sim_id,
    year = 1:n_years,
    B_true = numeric(n_years),
    B_estimated = numeric(n_years),
    F_advice = numeric(n_years),
    F_realized = numeric(n_years),
    catch = numeric(n_years),
    TAC = numeric(n_years),
    index_obs = numeric(n_years),
    B_Bmsy_true = numeric(n_years),
    F_Fmsy_true = numeric(n_years),
    HCR_status = character(n_years),
    stringsAsFactors = FALSE
  )
  
  # История наблюдений для оценки
  index_history <- numeric()
  catch_history <- numeric()
  previous_TAC <- NULL
  last_assessment <- NULL
  
  # Основной цикл симуляции по годам
  for (t in 1:n_years) {
    
    # 1. Генерация наблюдений (индекс биомассы)
    index_obs <- generate_index_observation(
      B_sim, 
      params$q, 
      settings$observation_error_cv
    )
    index_history <- c(index_history, index_obs)
    
    # 2. Оценка запаса (периодически или в первый год)
    if (t == 1 || (t - 1) %% settings$assessment_interval == 0) {
      
      last_assessment <- assess_stock_status(
        index_history,
        catch_history,
        params$q,
        params$Bmsy,
        params$Fmsy,
        settings$assessment_bias_cv
      )
      
      # 3. Применение правила управления (HCR)
      advice <- HCR_function(last_assessment, previous_TAC)
      F_advice <- advice$F_advice
      TAC <- advice$TAC
      previous_TAC <- TAC
      HCR_status <- advice$status
      
    } else {
      # Используем прошлогодние рекомендации
      F_advice <- results$F_advice[t-1]
      TAC <- results$TAC[t-1]
      HCR_status <- results$HCR_status[t-1]
    }
    
    # 4. Реализация промысла с ошибкой implementation
    implementation_error <- rlnorm(
      1, 
      meanlog = -settings$implementation_error_cv^2/2,
      sdlog = settings$implementation_error_cv
    )
    
    F_realized <- F_advice * implementation_error
    F_realized <- min(F_realized, 2.0)  # Ограничение максимальной F
    
    # 5. Расчет динамики популяции
    pop_update <- simulate_population_dynamics(
      B_sim, 
      F_realized,
      r_sim,
      K_sim,
      settings$process_error_cv
    )
    
    # 6. Сохранение результатов текущего года
    results$B_true[t] <- B_sim
    results$B_estimated[t] <- ifelse(!is.null(last_assessment), 
                                     last_assessment$B, NA)
    results$F_advice[t] <- F_advice
    results$F_realized[t] <- F_realized
    results$catch[t] <- pop_update$catch_realized
    results$TAC[t] <- TAC
    results$index_obs[t] <- index_obs
    results$B_Bmsy_true[t] <- B_sim / params$Bmsy
    results$F_Fmsy_true[t] <- F_realized / params$Fmsy
    results$HCR_status[t] <- HCR_status
    
    # 7. Обновление состояния для следующего года
    B_sim <- pop_update$B_next
    catch_history <- c(catch_history, pop_update$catch_realized)
  }
  
  return(results)
}

# ------------------- 7. ЗАПУСК MSE ДЛЯ ВСЕХ СТРАТЕГИЙ --------------------

cat("\n========== ЗАПУСК СИМУЛЯЦИЙ MSE ==========\n")

## 7.1 Подготовка параметров для симуляций
sim_params <- list(
  r = r_true,
  K = K_true,
  Bmsy = Bmsy_true,
  Fmsy = Fmsy_true,
  B_initial = B_current,
  q = q_cpue
)

sim_settings <- list(
  n_years = n_years,
  assessment_interval = assessment_interval,
  process_error_cv = process_error_cv,
  observation_error_cv = observation_error_cv,
  assessment_bias_cv = assessment_bias_cv,
  implementation_error_cv = implementation_error_cv
)

## 7.2 Список стратегий для сравнения
strategies <- list(
  "Fish at Fmsy" = HCR_Fmsy,
  "MSY Hockey-stick" = HCR_hockey_stick,
  "ICES Advice Rule" = HCR_ICES
)

## 7.3 Запуск симуляций для каждой стратегии
all_results <- list()

for (strategy_name in names(strategies)) {
  
  cat(sprintf("\nЗапуск стратегии: %s\n", strategy_name))
  cat("Прогресс: ")
  
  # Запуск n_sim симуляций для текущей стратегии
  strategy_results <- list()
  
  # Прогресс-индикатор
  pb <- txtProgressBar(min = 0, max = n_sim, style = 3, width = 50)
  
  for (i in 1:n_sim) {
    sim_result <- run_mse_simulation(
      sim_id = i,
      HCR_function = strategies[[strategy_name]],
      params = sim_params,
      settings = sim_settings
    )
    sim_result$strategy <- strategy_name
    strategy_results[[i]] <- sim_result
    
    setTxtProgressBar(pb, i)
  }
  
  close(pb)
  
  # Объединение результатов симуляций
  all_results[[strategy_name]] <- bind_rows(strategy_results)
  
  cat(sprintf("\n✓ Завершено %d симуляций для %s\n", n_sim, strategy_name))
}

## 7.4 Объединение всех результатов в один датафрейм
mse_results <- bind_rows(all_results)

cat("\n✓ Всего выполнено:", n_sim * length(strategies), "симуляций\n")

# ------------------- 8. РАСЧЕТ РИСК-МЕТРИК --------------------

cat("\n========== РАСЧЕТ РИСК-МЕТРИК ==========\n")

## 8.1 Функция для расчета комплексных метрик производительности
calculate_performance_metrics <- function(results) {
  
  metrics <- results %>%
    group_by(strategy) %>%
    summarise(
      
      # === БИОЛОГИЧЕСКИЕ МЕТРИКИ ===
      # Вероятность перелова (F > Fmsy)
      prob_overfishing = mean(F_Fmsy_true > 1, na.rm = TRUE),
      
      # Вероятность истощения запаса (B < 0.5*Bmsy)
      prob_overfished = mean(B_Bmsy_true < 0.5, na.rm = TRUE),
      
      # Вероятность коллапса (B < 0.2*Bmsy)
      prob_collapsed = mean(B_Bmsy_true < 0.2, na.rm = TRUE),
      
      # Вероятность нахождения в "зеленой зоне" Кобе
      prob_green_zone = mean(B_Bmsy_true > 1 & F_Fmsy_true < 1, na.rm = TRUE),
      
      # Средние показатели за последние 5 лет
      mean_B_Bmsy_final = mean(B_Bmsy_true[year > (max(year) - 5)], na.rm = TRUE),
      mean_F_Fmsy_final = mean(F_Fmsy_true[year > (max(year) - 5)], na.rm = TRUE),
      
      # === ЭКОНОМИЧЕСКИЕ МЕТРИКИ ===
      # Средний вылов за весь период
      mean_catch = mean(catch, na.rm = TRUE),
      
      # Суммарный вылов за период
      total_catch = sum(catch) / n_distinct(sim_id),
      
      # Стабильность вылова (обратная величина коэффициента вариации)
      catch_stability = {
        annual_catch <- aggregate(catch, by = list(year), mean)$x
        1 - (sd(annual_catch, na.rm = TRUE) / mean(annual_catch, na.rm = TRUE))
      },
      
      # Средняя межгодовая изменчивость вылова (AAV)
      catch_aav = {
        annual_catch <- aggregate(catch, by = list(year), mean)$x
        mean(abs(diff(annual_catch)) / annual_catch[-length(annual_catch)], na.rm = TRUE)
      },
      
      # === УПРАВЛЕНЧЕСКИЕ МЕТРИКИ ===
      # Вероятность закрытия промысла
      prob_closure = mean(TAC == 0, na.rm = TRUE),
      
      # Частота изменения рекомендаций
      advice_changes = {
        changes <- aggregate(F_advice, by = list(sim_id), 
                           function(x) sum(diff(x) != 0))$x
        mean(changes) / (max(year) - 1)
      },
      
      .groups = "drop"
    )
  
  return(metrics)
}

## 8.2 Расчет метрик для каждой стратегии
performance_metrics <- calculate_performance_metrics(mse_results)

## 8.3 Добавление комплексных оценок
performance_metrics <- performance_metrics %>%
  mutate(
    # Биологический индекс (0-1, где 1 - лучше)
    bio_score = 0.4 * (1 - prob_overfished) + 
                0.3 * (1 - prob_overfishing) + 
                0.3 * prob_green_zone,
    
    # Экономический индекс (0-1, где 1 - лучше)
    econ_score = 0.5 * (mean_catch / max(mean_catch)) + 
                 0.3 * catch_stability +
                 0.2 * (1 - catch_aav),
    
    # Общий индекс производительности
    overall_score = 0.6 * bio_score + 0.4 * econ_score
  ) %>%
  arrange(desc(overall_score))

## 8.4 Вывод таблицы метрик
cat("\n--- ТАБЛИЦА РИСК-МЕТРИК ---\n\n")

# Форматированный вывод ключевых метрик
for (strat in unique(performance_metrics$strategy)) {
  metrics_row <- performance_metrics[performance_metrics$strategy == strat, ]
  
  cat(sprintf("СТРАТЕГИЯ: %s\n", strat))
  cat(strrep("-", 40), "\n")
  cat(sprintf("Биологические риски:\n"))
  cat(sprintf("  P(перелов): %.1f%%\n", metrics_row$prob_overfishing * 100))
  cat(sprintf("  P(истощение): %.1f%%\n", metrics_row$prob_overfished * 100))
  cat(sprintf("  P(коллапс): %.1f%%\n", metrics_row$prob_collapsed * 100))
  cat(sprintf("  P(зеленая зона): %.1f%%\n", metrics_row$prob_green_zone * 100))
  cat(sprintf("Экономические показатели:\n"))
  cat(sprintf("  Средний вылов: %.1f тыс. т\n", metrics_row$mean_catch))
  cat(sprintf("  Стабильность: %.2f\n", metrics_row$catch_stability))
  cat(sprintf("  Межгодовая изменчивость: %.1f%%\n", metrics_row$catch_aav * 100))
  cat(sprintf("Итоговые оценки:\n"))
  cat(sprintf("  Биологический индекс: %.2f\n", metrics_row$bio_score))
  cat(sprintf("  Экономический индекс: %.2f\n", metrics_row$econ_score))
  cat(sprintf("  ОБЩИЙ ИНДЕКС: %.2f\n", metrics_row$overall_score))
  cat("\n")
}

# ------------------- 9. ВИЗУАЛИЗАЦИЯ РЕЗУЛЬТАТОВ --------------------

cat("========== СОЗДАНИЕ ГРАФИКОВ ==========\n")

## 9.1 Настройка темы для графиков
theme_mse <- theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 11),
    strip.text = element_text(size = 11, face = "bold"),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

## 9.2 Подготовка сводных данных для визуализации
summary_data <- mse_results %>%
  group_by(strategy, year) %>%
  summarise(
    # Квантили для биомассы
    B_median = median(B_true, na.rm = TRUE),
    B_q25 = quantile(B_true, 0.25, na.rm = TRUE),
    B_q75 = quantile(B_true, 0.75, na.rm = TRUE),
    B_q05 = quantile(B_true, 0.05, na.rm = TRUE),
    B_q95 = quantile(B_true, 0.95, na.rm = TRUE),
    
    # Квантили для B/Bmsy
    B_Bmsy_median = median(B_Bmsy_true, na.rm = TRUE),
    B_Bmsy_q25 = quantile(B_Bmsy_true, 0.25, na.rm = TRUE),
    B_Bmsy_q75 = quantile(B_Bmsy_true, 0.75, na.rm = TRUE),
    
    # Квантили для F/Fmsy
    F_Fmsy_median = median(F_Fmsy_true, na.rm = TRUE),
    F_Fmsy_q25 = quantile(F_Fmsy_true, 0.25, na.rm = TRUE),
    F_Fmsy_q75 = quantile(F_Fmsy_true, 0.75, na.rm = TRUE),
    
    # Квантили для вылова
    catch_median = median(catch, na.rm = TRUE),
    catch_q25 = quantile(catch, 0.25, na.rm = TRUE),
    catch_q75 = quantile(catch, 0.75, na.rm = TRUE),
    
    .groups = "drop"
  )

## 9.3 График 1: Траектории B/Bmsy
p1_depletion <- ggplot(summary_data, aes(x = year, color = strategy, fill = strategy)) +
  geom_ribbon(aes(ymin = B_Bmsy_q25, ymax = B_Bmsy_q75), alpha = 0.3, color = NA) +
  geom_line(aes(y = B_Bmsy_median), size = 1.5) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black", alpha = 0.5) +
  geom_hline(yintercept = 0.5, linetype = "dotted", color = "red", alpha = 0.5) +
  scale_color_manual(values = c("Fish at Fmsy" = "#E41A1C", 
                               "MSY Hockey-stick" = "#377EB8",
                               "ICES Advice Rule" = "#4DAF4A")) +
  scale_fill_manual(values = c("Fish at Fmsy" = "#E41A1C", 
                              "MSY Hockey-stick" = "#377EB8",
                              "ICES Advice Rule" = "#4DAF4A")) +
  labs(title = "Динамика относительной биомассы (B/Bmsy)",
       subtitle = "Медиана и межквартильный размах",
       x = "Год прогноза", 
       y = "B/Bmsy",
       color = "Стратегия",
       fill = "Стратегия") +
  theme_mse +
  annotate("text", x = max(summary_data$year), y = 1.02, 
           label = "Bmsy", hjust = 1, vjust = 0, size = 3) +
  annotate("text", x = max(summary_data$year), y = 0.52, 
           label = "0.5 Bmsy", hjust = 1, vjust = 0, size = 3, color = "red")


p1_depletion


## 9.4 График 2: Траектории F/Fmsy
p2_fishing <- ggplot(summary_data, aes(x = year, color = strategy, fill = strategy)) +
  geom_ribbon(aes(ymin = F_Fmsy_q25, ymax = F_Fmsy_q75), alpha = 0.3, color = NA) +
  geom_line(aes(y = F_Fmsy_median), size = 1.5) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black", alpha = 0.5) +
  scale_color_manual(values = c("Fish at Fmsy" = "#E41A1C", 
                               "MSY Hockey-stick" = "#377EB8",
                               "ICES Advice Rule" = "#4DAF4A")) +
  scale_fill_manual(values = c("Fish at Fmsy" = "#E41A1C", 
                              "MSY Hockey-stick" = "#377EB8",
                              "ICES Advice Rule" = "#4DAF4A")) +
  labs(title = "Динамика промысловой смертности (F/Fmsy)",
       subtitle = "Медиана и межквартильный размах",
       x = "Год прогноза", 
       y = "F/Fmsy",
       color = "Стратегия",
       fill = "Стратегия") +
  theme_mse +
  coord_cartesian(ylim = c(0, 2)) +
  annotate("text", x = max(summary_data$year), y = 1.02, 
           label = "Fmsy", hjust = 1, vjust = 0, size = 3)

p2_fishing

## 9.5 График 3: Динамика вылова
p3_catch <- ggplot(summary_data, aes(x = year, color = strategy, fill = strategy)) +
  geom_ribbon(aes(ymin = catch_q25, ymax = catch_q75), alpha = 0.3, color = NA) +
  geom_line(aes(y = catch_median), size = 1.5) +
  scale_color_manual(values = c("Fish at Fmsy" = "#E41A1C", 
                               "MSY Hockey-stick" = "#377EB8",
                               "ICES Advice Rule" = "#4DAF4A")) +
  scale_fill_manual(values = c("Fish at Fmsy" = "#E41A1C", 
                              "MSY Hockey-stick" = "#377EB8",
                              "ICES Advice Rule" = "#4DAF4A")) +
  labs(title = "Динамика вылова",
       subtitle = "Медиана и межквартильный размах",
       x = "Год прогноза", 
       y = "Вылов (тыс. т)",
       color = "Стратегия",
       fill = "Стратегия") +
  theme_mse

p3_catch 


## 9.6 График 4: Фазовая диаграмма Кобе для последнего года
kobe_data <- mse_results %>%
  filter(year == max(year)) %>%
  group_by(strategy) %>%
  sample_n(min(100, n()))  # Максимум 100 точек для читаемости

p4_kobe <- ggplot(kobe_data, aes(x = B_Bmsy_true, y = F_Fmsy_true)) +
  # Зоны Кобе
  annotate("rect", xmin = 0, xmax = 1, ymin = 1, ymax = 3,
           fill = "red", alpha = 0.15) +        # Красная зона
  annotate("rect", xmin = 1, xmax = 3, ymin = 1, ymax = 3,
           fill = "#FFA500", alpha = 0.15) +    # Оранжевая зона  
  annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1,
           fill = "#FFFF00", alpha = 0.15) +    # Желтая зона
  annotate("rect", xmin = 1, xmax = 3, ymin = 0, ymax = 1,
           fill = "green", alpha = 0.15) +      # Зеленая зона
  
  # Точки симуляций
  geom_point(aes(color = strategy), alpha = 0.5, size = 1.5) +
  
  # Референсные линии
  geom_vline(xintercept = 1, linetype = "solid", size = 0.8) +
  geom_hline(yintercept = 1, linetype = "solid", size = 0.8) +
  
  # Настройки
  facet_wrap(~ strategy, ncol = 3) +
  scale_color_manual(values = c("Fish at Fmsy" = "#E41A1C", 
                               "MSY Hockey-stick" = "#377EB8",
                               "ICES Advice Rule" = "#4DAF4A")) +
  labs(title = "Фазовая диаграмма Кобе (год 20)",
       subtitle = "Распределение конечных состояний",
       x = "B/Bmsy", 
       y = "F/Fmsy") +
  theme_mse +
  theme(legend.position = "none") +
  coord_cartesian(xlim = c(0, 2), ylim = c(0, 2))

p4_kobe

## 9.7 График 5: Сравнение риск-метрик (барплот)
metrics_long <- performance_metrics %>%
  select(strategy, prob_overfishing, prob_overfished, 
         prob_collapsed, prob_green_zone) %>%
  pivot_longer(cols = -strategy, names_to = "metric", values_to = "probability")

# Переименование метрик для графика
metric_labels <- c(
  prob_overfishing = "Перелов\n(F > Fmsy)",
  prob_overfished = "Истощение\n(B < 0.5 Bmsy)",
  prob_collapsed = "Коллапс\n(B < 0.2 Bmsy)",
  prob_green_zone = "Зеленая зона\n(B > Bmsy, F < Fmsy)"
)

p5_risks <- ggplot(metrics_long, aes(x = metric, y = probability, fill = strategy)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = sprintf("%.1f%%", probability * 100)),
            position = position_dodge(width = 0.8), 
            vjust = -0.5, size = 3) +
  scale_fill_manual(values = c("Fish at Fmsy" = "#E41A1C", 
                              "MSY Hockey-stick" = "#377EB8",
                              "ICES Advice Rule" = "#4DAF4A")) +
  scale_x_discrete(labels = metric_labels) +
  scale_y_continuous(labels = scales::percent, limits = c(0, max(metrics_long$probability) * 1.1)) +
  labs(title = "Сравнение вероятностей риска",
       subtitle = "По всему периоду симуляции",
       x = "", 
       y = "Вероятность",
       fill = "Стратегия") +
  theme_mse

p5_risks

## 9.8 График 6: Компромиссы (trade-offs)
p6_tradeoff <- ggplot(performance_metrics, 
                      aes(x = mean_catch, y = 1 - prob_overfished)) +
  geom_point(aes(color = strategy, size = catch_stability), alpha = 0.7) +
  geom_text(aes(label = strategy), vjust = -1.5, size = 3.5) +
  scale_color_manual(values = c("Fish at Fmsy" = "#E41A1C", 
                               "MSY Hockey-stick" = "#377EB8",
                               "ICES Advice Rule" = "#4DAF4A")) +
  scale_size_continuous(range = c(8, 15), name = "Стабильность\nвылова") +
  labs(title = "Анализ компромиссов между целями управления",
       subtitle = "Размер точки = стабильность вылова",
       x = "Средний вылов (тыс. т)", 
       y = "Вероятность избежать истощения",
       color = "Стратегия") +
  theme_mse +
  theme(legend.position = "right")

p6_tradeoff


# ------------------- 10. ДЕТАЛЬНЫЙ АНАЛИЗ --------------------

cat("\n========== ДЕТАЛЬНЫЙ АНАЛИЗ РЕЗУЛЬТАТОВ ==========\n")

## 10.1 Анализ времени восстановления
cat("\n--- Анализ восстановления запаса ---\n")

if (B_current/Bmsy_true < 1) {
  recovery_analysis <- mse_results %>%
    group_by(strategy, sim_id) %>%
    summarise(
      recovery_year = which(B_Bmsy_true > 1)[1],
      recovered = !is.na(recovery_year),
      .groups = "drop"
    ) %>%
    group_by(strategy) %>%
    summarise(
      prob_recovery = mean(recovered, na.rm = TRUE),
      median_recovery_time = median(recovery_year, na.rm = TRUE),
      q25_recovery = quantile(recovery_year, 0.25, na.rm = TRUE),
      q75_recovery = quantile(recovery_year, 0.75, na.rm = TRUE),
      .groups = "drop"
    )
  
  for (i in 1:nrow(recovery_analysis)) {
    cat(sprintf("%s:\n", recovery_analysis$strategy[i]))
    cat(sprintf("  Вероятность восстановления: %.1f%%\n", 
                recovery_analysis$prob_recovery[i] * 100))
    if (!is.na(recovery_analysis$median_recovery_time[i])) {
      cat(sprintf("  Медианное время восстановления: %.0f лет [%.0f-%.0f]\n",
                  recovery_analysis$median_recovery_time[i],
                  recovery_analysis$q25_recovery[i],
                  recovery_analysis$q75_recovery[i]))
    }
  }
} else {
  cat("Запас находится на уровне B/Bmsy >= 1. Восстановление не требуется.\n")
  cat(sprintf("Текущее состояние: B/Bmsy = %.2f\n", B_current / Bmsy_true))
}

## 10.2 Анализ частоты закрытия промысла
closure_analysis <- mse_results %>%
  group_by(strategy) %>%
  summarise(
    total_closures = sum(TAC == 0),
    closure_rate = mean(TAC == 0),
    avg_closure_duration = {
      rle_results <- rle(TAC == 0)
      if (any(rle_results$values)) {
        mean(rle_results$lengths[rle_results$values])
      } else {
        0
      }
    },
    .groups = "drop"
  )

cat("\n--- Анализ закрытия промысла ---\n")
for (i in 1:nrow(closure_analysis)) {
  cat(sprintf("%s:\n", closure_analysis$strategy[i]))
  cat(sprintf("  Частота закрытия: %.1f%%\n", closure_analysis$closure_rate[i] * 100))
  if (closure_analysis$avg_closure_duration[i] > 0) {
    cat(sprintf("  Средняя продолжительность закрытия: %.1f лет\n",
                closure_analysis$avg_closure_duration[i]))
  }
}

# ------------------- 11. СОХРАНЕНИЕ РЕЗУЛЬТАТОВ --------------------

cat("\n========== СОХРАНЕНИЕ РЕЗУЛЬТАТОВ ==========\n")


## 11.2 Сохранение таблиц
write.csv(performance_metrics, "MSE_performance_metrics.csv", row.names = FALSE)
write.csv(summary_data, "MSE_summary_by_year.csv", row.names = FALSE)
cat("✓ Таблицы метрик сохранены\n")

## 11.3 Сохранение полных результатов
saveRDS(mse_results, "MSE_full_results.rds")
cat("✓ Полные результаты симуляций сохранены: 'MSE_full_results.rds'\n")

## 11.4 Создание текстового отчета
sink("MSE_report.txt")
cat("ОТЧЕТ ПО MANAGEMENT STRATEGY EVALUATION\n")
cat(strrep("=", 60), "\n\n")
cat("Дата анализа:", format(Sys.Date(), "%d.%m.%Y"), "\n")
cat("Версия SPiCT:", as.character(packageVersion("spict")), "\n\n")

cat("ПАРАМЕТРЫ СИМУЛЯЦИИ:\n")
cat(sprintf("  Количество симуляций: %d\n", n_sim))
cat(sprintf("  Период прогноза: %d лет\n", n_years))
cat(sprintf("  Интервал оценки: каждые %d года\n", assessment_interval))
cat(sprintf("  CV процессной ошибки: %.0f%%\n", process_error_cv * 100))
cat(sprintf("  CV ошибки наблюдения: %.0f%%\n", observation_error_cv * 100))
cat(sprintf("  CV смещения оценок: %.0f%%\n", assessment_bias_cv * 100))
cat(sprintf("  CV ошибки реализации: %.0f%%\n\n", implementation_error_cv * 100))

cat("ИСХОДНОЕ СОСТОЯНИЕ ЗАПАСА:\n")
cat(sprintf("  B/Bmsy: %.2f\n", B_current/Bmsy_true))
cat(sprintf("  F/Fmsy: %.2f\n\n", F_current/Fmsy_true))

cat("РЕЗУЛЬТАТЫ ПО СТРАТЕГИЯМ:\n\n")
for (i in 1:nrow(performance_metrics)) {
  cat(sprintf("СТРАТЕГИЯ: %s\n", performance_metrics$strategy[i]))
  cat(strrep("-", 40), "\n")
  cat(sprintf("  Риск перелова: %.1f%%\n", performance_metrics$prob_overfishing[i] * 100))
  cat(sprintf("  Риск истощения: %.1f%%\n", performance_metrics$prob_overfished[i] * 100))
  cat(sprintf("  Вероятность устойчивого состояния: %.1f%%\n", 
              performance_metrics$prob_green_zone[i] * 100))
  cat(sprintf("  Средний вылов: %.1f тыс. т\n", performance_metrics$mean_catch[i]))
  cat(sprintf("  Стабильность вылова: %.2f\n", performance_metrics$catch_stability[i]))
  cat(sprintf("  Общий индекс производительности: %.3f\n\n", 
              performance_metrics$overall_score[i]))
}

cat("РЕКОМЕНДАЦИИ:\n")
best_strategy <- performance_metrics$strategy[1]
cat(sprintf("Оптимальная стратегия: %s\n", best_strategy))
cat("\nОбоснование:\n")
if (best_strategy == "ICES Advice Rule") {
  cat("- Предосторожный подход минимизирует биологические риски\n")
  cat("- Ограничение межгодовых изменений обеспечивает стабильность для промышленности\n")
  cat("- Баланс между сохранением запаса и экономической эффективностью\n")
} else if (best_strategy == "MSY Hockey-stick") {
  cat("- Адаптивное управление в зависимости от состояния запаса\n")
  cat("- Защита запаса при низкой биомассе\n")
  cat("- Возможность полного использования при хорошем состоянии\n")
} else {
  cat("- Простота применения\n")
  cat("- Максимизация долгосрочного вылова\n")
  cat("- Требует точной оценки Fmsy\n")
}

sink()
cat("✓ Текстовый отчет сохранен: 'MSE_report.txt'\n")

cat("\n" , strrep("=", 60), "\n")
cat("         MSE АНАЛИЗ УСПЕШНО ЗАВЕРШЕН\n")
cat(strrep("=", 60), "\n")


