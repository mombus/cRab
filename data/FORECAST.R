# =============================================================================
#  ПРАКТИКУМ: ПРОГНОСТИЧЕСКАЯ СИЛА ПРОДУКЦИОННОЙ МОДЕЛИ (SPiCT)
# =============================================================================
# 
#  Период данных: 2005–2024 
#  Оценка прогноза: 2022–2024 (горизонт h = 3 года)
# 
#  Цель: Сравнить «истинную» биомассу (из SPiCT) с прогнозами базовых методов:
#        Naive, Drift, Mean, Theta, ETS, ARIMA и SPiCT-dynamic
#
#  Автор: [Ваше имя]
#  Дата: [Текущая дата]
# =============================================================================

cat("\n" + strrep("=", 70))
cat("\n  ПРАКТИКУМ SPiCT: ОЦЕНКА ПРОГНОСТИЧЕСКОЙ СИЛЫ")
cat("\n  Период: 2005–2024 | Прогноз: 2022–2024 (h=3)")
cat("\n" + strrep("=", 70) + "\n")

# =============================================================================
#  БЛОК 0: ПОДГОТОВКА РАБОЧЕЙ СРЕДЫ И ПАКЕТОВ
# =============================================================================

cat("\n?? БЛОК 0: Подготовка рабочей среды...")

# Функция для установки и загрузки пакетов (если они не установлены)
install_and_load <- function(packages) {
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message("?? Устанавливаю пакет: ", pkg)
      install.packages(pkg, repos = "https://cloud.r-project.org")
    }
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
}

# Установка и загрузка необходимых пакетов
install_and_load(c("spict", "forecast", "ggplot2", "dplyr", "tidyr", "scales"))

# Настройки вывода чисел (отключаем научную нотацию)
options(scipen = 999, digits = 4)

cat(" ? Готово!\n")

# =============================================================================
#  БЛОК 1: ЗАГРУЗКА И ПОДГОТОВКА ИСХОДНЫХ ДАННЫХ
# =============================================================================

cat("\n?? БЛОК 1: Загрузка исходных данных...")

# Вектор лет наблюдений
years_all <- 2005:2024

# Улов (тыс. тонн) по годам
# 2005–2014: рост уловов, 2015–2024: снижение и стабилизация
catch_tonnes <- c(
  5, 7, 6, 10, 14, 25, 28, 30, 32, 35,  # 2005–2014
  25, 20, 15, 12, 10, 12, 10, 13, 11, 12 # 2015–2024
)

# Промысловый индекс (CPUE) - данные с середины года
cpue_index <- c(
  27.427120, 26.775958, 16.811997, 22.979653, 29.048568,
  29.996072, 16.476301, 17.174455, 10.537272, 14.590435,
  8.286352, 11.394168, 15.537878, 13.791166, 11.527548,
  15.336093, 12.154069, 15.568450, 16.221933, 13.421132
)

# Научный индекс (BESS) - данные с конца года (IV квартал)
# Первый год NA - съёмка не проводилась
bess_index <- c(
  NA, 16.270375, 20.691355, 15.141784, 18.594620,
  15.975548, 13.792012, 13.328805, 11.659744, 11.753855,
  9.309859, 7.104886, 7.963839, 9.161322, 10.271221,
  9.822960, 10.347376, 11.703610, 13.679876, 13.413696
)

# Проверка целостности данных (все векторы должны быть одинаковой длины)
stopifnot(
  length(years_all) == length(catch_tonnes),
  length(years_all) == length(cpue_index),
  length(years_all) == length(bess_index)
)

cat(" ? Данные загружены и проверены!\n")
cat("   • Годы наблюдений:", length(years_all), "лет (2005–2024)\n")
cat("   • Уловы:", length(catch_tonnes), "значений\n") 
cat("   • Индекс CPUE:", length(cpue_index), "значений\n")
cat("   • Индекс BESS:", sum(!is.na(bess_index)), "значений (1 пропуск)\n")

# =============================================================================
#  БЛОК 2: КАЛИБРОВКА ПОЛНОЙ МОДЕЛИ SPiCT НА ВСЕХ ДАННЫХ
# =============================================================================

cat("\n?? БЛОК 2: Калибровка полной модели SPiCT (2005–2024)...")

# Формируем входные данные для SPiCT
spict_input <- list(
  timeC = years_all,              # годы для улова (конец года)
  obsC  = catch_tonnes,           # значения улова
  timeI = list(
    years_all + 0.5,              # CPUE: середина года (июль)
    years_all + 0.75              # BESS: конец года (октябрь)
  ),
  obsI  = list(cpue_index, bess_index)
)

# Базовая проверка и подготовка входных данных
spict_inp <- check.inp(spict_input, verbose = FALSE)

# ?? НАСТРОЙКИ МОДЕЛИ ДЛЯ УЧЕБНОЙ ЗАДАЧИ:

# Фиксируем форму кривой производства (Шефер: n=2)
spict_inp$priors$logn <- c(log(2), 0.1, 1)  # априорное распределение
spict_inp$ini$logn    <- log(2)              # начальное значение
spict_inp$phases$logn <- -1                  # фиксируем параметр

# Информативный априор на емкость среды (K)
spict_inp$priors$logK <- c(log(150), 0.7, 1)

# Точность интегрирования (чем больше - тем точнее)
spict_inp$dteuler <- 1/16

# Калибровка модели (может занять несколько секунд)
fit_full <- fit.spict(spict_inp)

cat(" ? Модель откалибрована!\n")

# =============================================================================
#  БЛОК 3: ОПРЕДЕЛЕНИЕ «ИСТИННОЙ» БИОМАССЫ
# =============================================================================

cat("\n?? БЛОК 3: Определение эталонной ('истинной') биомассы...")

# Функция для извлечения биомассы на конец каждого года
extract_biomass_year_end <- function(fit, years, frac_end = 0.999) {
  # Получаем оценки биомассы в непрерывном времени
  gp <- get.par("logB", fit, exp = TRUE)
  est_values <- if (!is.null(dim(gp))) as.numeric(gp[, "est"]) else as.numeric(gp)
  time_grid  <- as.numeric(fit$inp$time)
  
  # Для каждого года находим ближайшую точку к концу года
  year_end_values <- vapply(years, function(y) {
    idx <- which(time_grid <= y + frac_end)
    if (length(idx) == 0) return(NA_real_)
    est_values[max(idx)]
  }, numeric(1))
  
  data.frame(year = years, B = year_end_values)
}

# Извлекаем "истинную" биомассу (эталон для сравнения)
biomass_truth <- extract_biomass_year_end(fit_full, years_all)

# Проверка качества данных
if (anyNA(biomass_truth$B)) {
  warning("?? В рассчитанной биомассе есть пропуски. Проверьте сходимость модели.")
} else {
  cat(" ? Биомасса рассчитана для всех", nrow(biomass_truth), "лет\n")
}

# =============================================================================
#  БЛОК 4: ПОДГОТОВКА ПРОГНОЗНЫХ МОДЕЛЕЙ (БАЗОВЫЕ МЕТОДЫ)
# =============================================================================

cat("\n?? БЛОК 4: Подготовка базовых прогнозных моделей...")

# Определяем периоды обучения и тестирования
horizon_years <- 3     # горизонт прогнозирования
years_train   <- years_all[years_all <= 2021]  # обучение: 2005-2021
years_test    <- 2022:2024                     # тест: 2022-2024

# Готовим обучающие данные (биомасса до 2021 года)
mask_train <- biomass_truth$year %in% years_train & is.finite(biomass_truth$B)
y_train_ts <- ts(
  biomass_truth$B[mask_train], 
  start = min(biomass_truth$year[mask_train]), 
  frequency = 1
)

# Истинные значения для тестового периода
true_test <- biomass_truth$B[biomass_truth$year %in% years_test]

cat("   • Обучение: 2005–2021 (", length(years_train), " лет)\n", sep = "")
cat("   • Тестирование: 2022–2024 (", length(years_test), " года)\n", sep = "")

# Калибровка моделей временных рядов
fit_ets   <- ets(y_train_ts)        # Экспоненциальное сглаживание
fit_arima <- auto.arima(y_train_ts) # Автоматический подбор ARIMA

# ?? СПИСОК ПРОГНОЗНЫХ МОДЕЛЕЙ:
fc_list <- list(
  Naive = naive(y_train_ts, h = horizon_years),    # Наивный прогноз
  Drift = rwf(y_train_ts, h = horizon_years, drift = TRUE),  # С дрейфом
  Mean  = meanf(y_train_ts, h = horizon_years),    # Среднее значение
  Theta = thetaf(y_train_ts, h = horizon_years),   # Тета-метод
  ETS   = forecast(fit_ets, h = horizon_years),    # Экспоненциальное сглаживание
  ARIMA = forecast(fit_arima, h = horizon_years)   # ARIMA
)

cat(" ? Все модели откалиброваны!\n")

# =============================================================================
#  БЛОК 5: SPiCT-DYNAMIC - ПРОГНОЗ С УЧЕТОМ ДИНАМИКИ ЗАПАСА
# =============================================================================

cat("\n?? БЛОК 5: Прогноз SPiCT-dynamic с учетом будущих уловов...")

# Функция для обрезки данных по году
trim_by_year <- function(time_vec, obs_vec, last_year) {
  keep <- time_vec <= last_year
  list(time = time_vec[keep], obs = obs_vec[keep])
}

# Создаем усеченный набор данных (только до 2021 года)
spict_inp_train <- list(
  timeC = spict_input$timeC,
  obsC  = spict_input$obsC,
  timeI = spict_input$timeI,
  obsI  = spict_input$obsI
)

# Обрезаем уловы строго до 2021 года
trim_c <- trim_by_year(spict_inp_train$timeC, spict_inp_train$obsC, 2021)
spict_inp_train$timeC <- trim_c$time
spict_inp_train$obsC  <- trim_c$obs

# Обрезаем индексы с небольшим запасом (до конца 2021)
for (k in seq_along(spict_inp_train$timeI)) {
  tpair <- trim_by_year(spict_inp_train$timeI[[k]], spict_inp_train$obsI[[k]], 2021.999)
  spict_inp_train$timeI[[k]] = tpair$time
  spict_inp_train$obsI[[k]]  = tpair$obs
}

# Применяем те же настройки модели
spict_inp_train <- check.inp(spict_inp_train, verbose = FALSE)
spict_inp_train$priors$logn <- c(log(2), 0.1, 1)
spict_inp_train$ini$logn    <- log(2)
spict_inp_train$phases$logn <- -1
spict_inp_train$priors$logK <- c(log(150), 0.7, 1)
spict_inp_train$dteuler     <- 1/16

# Калибровка модели на усеченных данных
fit_train <- try(fit.spict(spict_inp_train), silent = TRUE)

# Функция безопасного извлечения параметров
safe_get_par_exp <- function(fit, logname) {
  val <- NA_real_
  if (!is.null(fit$par.fixed) && logname %in% names(fit$par.fixed)) {
    val <- as.numeric(exp(fit$par.fixed[logname]))
  }
  if (!is.finite(val)) {
    gp <- try(get.par(logname, fit, exp = TRUE), silent = TRUE)
    if (!inherits(gp, "try-error")) {
      if (is.null(dim(gp))) {
        val <- if ("est" %in% names(gp)) as.numeric(gp["est"]) else as.numeric(gp[1])
      } else if ("est" %in% colnames(gp)) {
        val <- as.numeric(gp[, "est"][1])
      }
    }
  }
  val
}

# Прогноз по дискретной логистике с учетом фактических уловов
if (inherits(fit_train, "try-error")) {
  warning("? SPiCT (2005-2021) не сошёлся. SPiCT-dynamic не будет рассчитан.")
  spict_dynamic_df <- data.frame(year = years_test, point = NA_real_)
} else {
  # Извлекаем параметры модели
  r_hat <- safe_get_par_exp(fit_train, "logr")  # скорость роста
  K_hat <- safe_get_par_exp(fit_train, "logK")  # емкость среды
  B_2021 <- extract_biomass_year_end(fit_train, 2021)$B  # биомасса на старте
  
  # Резервный способ извлечения биомассы
  if (!is.finite(B_2021)) {
    gpB <- get.par("logB", fit_train, exp = TRUE)
    estB <- if (!is.null(dim(gpB))) as.numeric(gpB[, "est"]) else as.numeric(gpB)
    if (length(estB) > 0 && is.finite(tail(estB, 1))) {
      B_2021 <- tail(estB, 1)
    }
  }
  
  if (!is.finite(B_2021) || !is.finite(r_hat) || !is.finite(K_hat)) {
    warning("? Не удалось определить параметры для SPiCT-dynamic.")
    spict_dynamic_df <- data.frame(year = years_test, point = NA_real_)
  } else {
    # Фактические уловы в прогнозный период
    catch_future <- catch_tonnes[years_all %in% years_test]
    
    # ?? ДИСКРЕТНАЯ ЛОГИСТИЧЕСКАЯ МОДЕЛЬ:
    forward_logistic <- function(B0, r, K, Cvec) {
      out <- numeric(length(Cvec))
      Bcur <- B0
      for (i in seq_along(Cvec)) {
        # Основное уравнение: рост + добыча
        Bnext <- Bcur + r * Bcur * (1 - Bcur / K) - Cvec[i]
        Bnext <- max(Bnext, 0)  # биомасса не может быть отрицательной
        out[i] <- Bnext
        Bcur <- Bnext
      }
      out
    }
    
    # Расчет прогнозов
    spict_points <- forward_logistic(
      B0 = B_2021, 
      r = r_hat, 
      K = K_hat, 
      Cvec = catch_future
    )
    
    spict_dynamic_df <- data.frame(year = years_test, point = spict_points)
    cat(" ? SPiCT-dynamic рассчитан!\n")
    cat("   • Биомасса 2021:", round(B_2021, 1), "тыс. т\n")
    cat("   • Параметр r:", round(r_hat, 3), "\n")
    cat("   • Емкость K:", round(K_hat, 1), "тыс. т\n")
  }
}

# =============================================================================
#  БЛОК 6: ФОРМИРОВАНИЕ СВОДНОЙ ТАБЛИЦЫ ПРОГНОЗОВ
# =============================================================================

cat("\n?? БЛОК 6: Формирование сводной таблицы прогнозов...")

# Функция преобразования прогнозов в таблицу
forecast_to_df <- function(fc_obj, method_name, years_target) {
  # Извлекаем интервальные оценки (если есть)
  has_lower <- !is.null(fc_obj$lower)
  levs <- if (has_lower) colnames(fc_obj$lower) else character(0)
  lo80 <- if ("80%" %in% levs) fc_obj$lower[, "80%"] else NA_real_
  hi80 <- if ("80%" %in% levs) fc_obj$upper[, "80%"] else NA_real_
  lo95 <- if ("95%" %in% levs) fc_obj$lower[, "95%"] else NA_real_
  hi95 <- if ("95%" %in% levs) fc_obj$upper[, "95%"] else NA_real_

  data.frame(
    method = method_name,
    year   = years_target,
    point  = as.numeric(fc_obj$mean),  # точечный прогноз
    lo80   = as.numeric(lo80),         # нижняя граница 80% ДИ
    hi80   = as.numeric(hi80),         # верхняя граница 80% ДИ  
    lo95   = as.numeric(lo95),         # нижняя граница 95% ДИ
    hi95   = as.numeric(hi95),         # верхняя граница 95% ДИ
    stringsAsFactors = FALSE
  )
}

# Объединяем прогнозы всех базовых моделей
forecast_table <- dplyr::bind_rows(
  lapply(names(fc_list), function(nm) {
    forecast_to_df(fc_list[[nm]], nm, years_test)
  })
)

# Добавляем SPiCT-dynamic (без доверительных интервалов)
forecast_table <- dplyr::bind_rows(
  forecast_table,
  dplyr::mutate(
    spict_dynamic_df, 
    method = "SPiCT-dynamic",
    lo80 = NA_real_, hi80 = NA_real_, 
    lo95 = NA_real_, hi95 = NA_real_
  )
)

# Сортируем для удобства просмотра
forecast_table <- dplyr::arrange(forecast_table, method, year)

cat(" ? Таблица прогнозов сформирована!\n\n")
cat("?? ТАБЛИЦА ПРОГНОЗОВ (точечные оценки):\n")
cat(rep("-", 60), "\n")
print(forecast_table)

# =============================================================================
#  БЛОК 7: РАСЧЕТ МЕТРИК ТОЧНОСТИ И ПРОГНОСТИЧЕСКОЙ СИЛЫ
# =============================================================================

cat("\n?? БЛОК 7: Расчет метрик точности и прогностической силы...")

# Создаем таблицу с истинными значениями
truth_df <- data.frame(year = years_test, truth = true_test)

# ?? ОСНОВНЫЕ МЕТРИКИ ТОЧНОСТИ:
metrics_table <- forecast_table |>
  dplyr::select(method, year, point) |>
  dplyr::left_join(truth_df, by = "year") |>
  dplyr::group_by(method) |>
  dplyr::summarise(
    H    = dplyr::n(),                                   # число прогнозов
    MSE  = mean((point - truth)^2, na.rm = TRUE),        # средняя квадратичная ошибка
    RMSE = sqrt(MSE),                                    # корень из MSE
    MAE  = mean(abs(point - truth), na.rm = TRUE),       # средняя абсолютная ошибка
    .groups = "drop"
  )

# ?? ПРОГНОСТИЧЕСКАЯ СИЛА (относительно Naive-прогноза):
mse_naive <- metrics_table$MSE[metrics_table$method == "Naive"]
metrics_table <- dplyr::mutate(
  metrics_table,
  Skill_MSE_vs_Naive = 1 - MSE / mse_naive  # >0 лучше Naive, <0 хуже
)

cat(" ? Метрики рассчитаны!\n\n")
cat("?? МЕТРИКИ ТОЧНОСТИ И ПРОГНОСТИЧЕСКАЯ СИЛА (2022–2024):\n")
cat(rep("-", 70), "\n")
print(metrics_table)

# =============================================================================
#  БЛОК 8: ПРОЦЕНТНЫЕ ОТКЛОНЕНИЯ ОТ ИСТИНЫ
# =============================================================================

cat("\n?? БЛОК 8: Расчет процентных отклонений от истины...")

# Таблица процентных отклонений
percent_diff_table <- forecast_table |>
  dplyr::select(method, year, point) |>
  dplyr::left_join(truth_df, by = "year") |>
  dplyr::mutate(
    percent_diff     = 100 * (point - truth) / truth,  # относительное отклонение
    abs_percent_diff = abs(percent_diff)               # абсолютное отклонение
  ) |>
  dplyr::arrange(method, year)

cat(" ? Отклонения рассчитаны!\n\n")
cat("?? ПРОЦЕНТНЫЕ ОТКЛОНЕНИЯ ОТ ИСТИНЫ ПО ГОДАМ (2022–2024):\n")
cat(rep("-", 70), "\n")
print(percent_diff_table)

# =============================================================================
#  БЛОК 9: ВИЗУАЛИЗАЦИЯ РЕЗУЛЬТАТОВ
# =============================================================================

cat("\n?? БЛОК 9: Построение графиков и визуализаций...")

# 9.1 ОСНОВНОЙ ГРАФИК: ДИНАМИКА БИОМАССЫ И ПРОГНОЗЫ

cat("\n   ?? График 1: Динамика биомассы и прогнозы...")

# Точка последнего наблюдения для соединения линий прогнозов
last_train_point <- data.frame(
  year   = 2021,
  method = names(fc_list),
  point  = rep(as.numeric(tail(y_train_ts, 1)), length(fc_list))
)

# Данные для построения прогнозов
plot_forecasts <- forecast_table |>
  dplyr::select(year, method, point)

# Построение основного графика
p_biomass <- ggplot() +
  # Историческая динамика биомассы
  geom_line(
    data = dplyr::filter(biomass_truth, is.finite(B)),
    aes(x = year, y = B),
    color = "black", linewidth = 1.2
  ) +
  geom_point(
    data = dplyr::filter(biomass_truth, is.finite(B)),
    aes(x = year, y = B),
    color = "black", size = 2
  ) +
  # Линии прогнозов (пунктир)
  geom_line(
    data = dplyr::bind_rows(last_train_point, plot_forecasts),
    aes(x = year, y = point, color = method),
    linewidth = 0.9, linetype = "dashed", alpha = 0.7
  ) +
  # Точки прогнозов
  geom_point(
    data = plot_forecasts,
    aes(x = year, y = point, color = method),
    size = 2.5
  ) +
  # Вертикальная линия разделения обучения/прогноза
  geom_vline(xintercept = 2021.5, linetype = "dotted", color = "gray50") +
  scale_color_brewer(palette = "Set1") +
  labs(
    title = "Динамика биомассы и прогнозы на 2022–2024 годы",
    subtitle = "Сплошная линия: 'истинная' биомасса SPiCT | Пунктир: прогнозы",
    x = "Год", 
    y = "Биомасса (тыс. тонн)", 
    color = "Метод прогноза"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "gray40")
  )

cat(" ? Отрисован!\n")
print(p_biomass)

# 9.2 ГРАФИК ПРОГНОСТИЧЕСКОЙ СИЛЫ

cat("   ?? График 2: Прогностическая сила методов...")

p_skill <- metrics_table |>
  dplyr::mutate(
    method = factor(
      method, 
      levels = metrics_table$method[order(Skill_MSE_vs_Naive, decreasing = TRUE)]
    )
  ) |>
  ggplot(aes(x = method, y = Skill_MSE_vs_Naive, fill = method)) +
  geom_col(width = 0.7, alpha = 0.8) +
  geom_hline(yintercept = 0, linewidth = 0.6, linetype = "dashed", color = "red") +
  geom_text(
    aes(label = sprintf("%.1f%%", Skill_MSE_vs_Naive * 100)),
    vjust = -0.5, 
    size = 3.5,
    fontface = "bold"
  ) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1)) +
  scale_fill_brewer(palette = "Set1") +
  labs(
    title = "Прогностическая сила методов относительно Naive",
    subtitle = "Выше 0% = лучше Naive | Ниже 0% = хуже Naive",
    x = NULL, 
    y = "Прогностическая сила (1 - MSE_model / MSE_naive)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold")
  )

cat(" ? Отрисован!\n")
print(p_skill)

# 9.3 ТОЧЕЧНЫЙ ГРАФИК РЕЙТИНГА МЕТОДОВ

cat("   ?? График 3: Рейтинг методов по точности...")

# Ранжирование методов по средней абсолютной ошибке
method_ranking <- percent_diff_table %>%
  group_by(method) %>%
  summarise(
    mean_abs_diff = mean(abs_percent_diff),
    mean_diff = mean(percent_diff)
  ) %>%
  arrange(mean_abs_diff)

# Точечный график с рейтингом
p_dot_rank <- percent_diff_table %>%
  mutate(method = factor(method, levels = method_ranking$method)) %>%
  ggplot(aes(x = percent_diff, y = method, color = factor(year))) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(size = 4, alpha = 0.8) +
  geom_text(
    aes(label = sprintf("%+.1f%%", percent_diff)),
    position = position_nudge(y = 0.2),
    size = 3,
    color = "black"
  ) +
  scale_color_brewer(palette = "Set1", name = "Год") +
  scale_x_continuous(labels = function(x) paste0(x, "%"), limits = c(-25, 25)) +
  labs(
    title = "Рейтинг методов по точности прогнозов",
    subtitle = "Методы отсортированы по среднему абсолютному отклонению",
    x = "Процентное отклонение от истины",
    y = "Метод прогнозирования"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "top",
    panel.grid.major.y = element_line(color = "gray90"),
    plot.title = element_text(face = "bold")
  )

cat(" ? Отрисован!\n")
print(p_dot_rank)

# =============================================================================
#  ФИНАЛЬНЫЙ ВЫВОД РЕЗУЛЬТАТОВ
# =============================================================================

cat("\n")
cat(strrep("=", 70))
cat("\n?? АНАЛИЗ ЗАВЕРШЕН! РЕЗУЛЬТАТЫ ГОТОВЫ")
cat("\n")
cat(strrep("=", 70))

cat("\n\n?? ИТОГОВЫЕ РЕЗУЛЬТАТЫ:")
cat("\n")
cat(rep("-", 50), sep = "")

# Лучший метод по прогностической силе
best_method <- metrics_table |> 
  filter(Skill_MSE_vs_Naive == max(Skill_MSE_vs_Naive, na.rm = TRUE))

cat("\n?? ЛУЧШИЙ МЕТОД:", best_method$method)
cat(sprintf("\n   • Прогностическая сила: +%.1f%%", best_method$Skill_MSE_vs_Naive * 100))
cat(sprintf("\n   • RMSE: %.1f тыс. т", best_method$RMSE))

# Наихудший метод (исключая отрицательные значения)
worst_method <- metrics_table |> 
  filter(Skill_MSE_vs_Naive == min(Skill_MSE_vs_Naive[Skill_MSE_vs_Naive > -1], na.rm = TRUE))

cat("\n\n??  НАИХУДШИЙ МЕТОД:", worst_method$method)
cat(sprintf("\n   • Прогностическая сила: %.1f%%", worst_method$Skill_MSE_vs_Naive * 100))
cat(sprintf("\n   • RMSE: %.1f тыс. т", worst_method$RMSE))

# SPiCT-dynamic в сравнении
spict_perf <- metrics_table |> filter(method == "SPiCT-dynamic")
if (nrow(spict_perf) > 0 && !is.na(spict_perf$Skill_MSE_vs_Naive)) {
  cat("\n\n?? SPiCT-DYNAMIC:")
  cat(sprintf("\n   • Прогностическая сила: +%.1f%%", spict_perf$Skill_MSE_vs_Naive * 100))
  cat(sprintf("\n   • Ранг по точности: %d из %d", 
              which(method_ranking$method == "SPiCT-dynamic"),
              nrow(method_ranking)))
}

cat("\n\n?? ВИЗУАЛИЗАЦИИ:")
cat("\n   • График 1: Динамика биомассы и прогнозы")
cat("\n   • График 2: Прогностическая сила методов")  
cat("\n   • График 3: Рейтинг методов по точности")

cat("\n\n")
cat(strrep("=", 70))
cat("\n? СКРИПТ УСПЕШНО ВЫПОЛНЕН!")
cat("\n")
cat(strrep("=", 70))
cat("\n\n")

