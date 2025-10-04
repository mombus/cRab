# ПРАКТИЧЕСКОЕ ЗАНЯТИЕ: ОСНОВЫ ПРОГНОЗИРОВАНИЯ
# УЛУЧШЕННАЯ ВЕРСИЯ - УЧЕТ БИОЛОГИЧЕСКИХ ОГРАНИЧЕНИЙ, ДИАГНОСТИКА, CV, СЦЕНАРИИ

# Очистка рабочей среды (опционально)
rm(list = ls())

# Фиксируем воспроизводимость
set.seed(123)

# =========================
# 0. ЗАГРУЗКА/УСТАНОВКА ПАКЕТОВ
# =========================
required_pkgs <- c("tidyverse", "forecast", "tseries", "ggplot2", "jsonlite", "readr")
ensure_packages <- function(pkgs) {
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
      try({
        install.packages(p, repos = "https://cloud.r-project.org", quiet = TRUE)
      }, silent = TRUE)
    }
    suppressPackageStartupMessages(library(p, character.only = TRUE))
  }
}
ensure_packages(required_pkgs)

# =========================
# 1. ПОДГОТОВКА ДАННЫХ
# =========================
# Данные индекса обилия (пример)
index_data <- data.frame(
  YEAR = 2019:2024,
  INDEX = c(2381774, 1634549, 1920507, 1036673, 1147685, 1055733)
)

cat("=== АНАЛИЗ ДИНАМИКИ ЗАПАСА КАМЧАТСКОГО КРАБА ===\n")
cat("Период наблюдений:", min(index_data$YEAR), "-", max(index_data$YEAR), "\n")
cat("Изменение запаса за период:", 
    round((index_data$INDEX[6] - index_data$INDEX[1])/index_data$INDEX[1]*100, 1), "%\n")

# Временной ряд (годовая частота)
ts_data <- ts(index_data$INDEX, start = min(index_data$YEAR), frequency = 1)

# Директории для вывода
outputs_dir <- file.path(getwd(), "outputs")
if (!dir.exists(outputs_dir)) dir.create(outputs_dir, recursive = TRUE)

# Базовая визуализация
p_hist <- ggplot(index_data, aes(x = YEAR, y = INDEX/1e6)) +
  geom_line(linewidth = 1, color = "steelblue") +
  geom_point(size = 3, color = "navy") +
  geom_text(aes(label = round(INDEX/1e6, 1)), vjust = -1, size = 3.5) +
  labs(title = "Динамика индекса обилия камчатского краба",
       subtitle = paste0(min(index_data$YEAR), "-", max(index_data$YEAR), 
                         " гг. (значения в млн экз.)"),
       x = "Год", y = "Индекс обилия, млн экз.") +
  ylim(0, 2.5) +
  theme_minimal()

ggsave(filename = file.path(outputs_dir, "historical_index.png"), plot = p_hist, 
       width = 8, height = 4.5, dpi = 150)

# =========================
# 2. АНАЛИЗ СТАЦИОНАРНОСТИ
# =========================
cat("\n=== АНАЛИЗ СТАЦИОНАРНОСТИ РЯДА ===\n")
cat("Длина ряда:", length(ts_data), "наблюдений\n")
cat("Среднее значение:", round(mean(ts_data), 0), "\n")
cat("Стандартное отклонение:", round(sd(ts_data), 0), "\n")

if (mean(diff(ts_data)) < 0) {
  cat("Визуальный анализ: присутствует нисходящий тренд\n")
} else {
  cat("Визуальный анализ: присутствует восходящий тренд\n")
}

# ADF-тест (с осторожностью при малой выборке)
safe_adf <- function(x) {
  out <- try(tseries::adf.test(x), silent = TRUE)
  if (inherits(out, "try-error")) return(NULL) else return(out)
}
adf_res <- safe_adf(ts_data)
if (!is.null(adf_res)) {
  cat("ADF тест: p-value =", signif(adf_res$p.value, 3), "(малая выборка — интерпретация осторожна)\n")
}

# =========================
# 3. МОДЕЛИ С БИОЛОГИЧЕСКИМИ ОГРАНИЧЕНИЯМИ
# =========================
# Разделение на обучающую (до 2023) и тестовую (2024)
train <- window(ts_data, end = max(index_data$YEAR) - 1)
test  <- window(ts_data, start = max(index_data$YEAR))
actual_2024 <- as.numeric(test)

# Параметры биологических ограничений
K_default <- 1.5 * max(ts_data, na.rm = TRUE)  # Носительная способность (априорно)
K <- K_default

# Ограничение годового темпа прироста (опционально). Например, [-80%, +60%]
growth_cap_enable <- TRUE
g_bounds <- c(-0.8, 0.6)

clip_to_bounds <- function(x, lower = 0, upper = Inf) {
  x <- pmax(x, lower)
  x <- pmin(x, upper)
  return(x)
}

adjust_by_growth_cap <- function(path, last_actual, g_bounds) {
  # Поэтапно ограничиваем относительные изменения для каждого шага прогноза
  out <- numeric(length(path))
  prev <- last_actual
  for (i in seq_along(path)) {
    raw <- path[i]
    gr <- (raw - prev) / prev
    gr <- min(max(gr, g_bounds[1]), g_bounds[2])
    adj <- prev * (1 + gr)
    out[i] <- adj
    prev <- adj
  }
  out
}

# Логистические вспомогательные функции
logit <- function(p) log(p/(1 - p))
inv_logit <- function(z) 1/(1 + exp(-z))

# ========== Определение моделей ==========
# Модель 1: Наивная (Random Walk)
fit_naive <- forecast::naive(train, h = 2)

# Модель 2: RW с дрейфом
fit_rwd <- forecast::rwf(train, h = 2, drift = TRUE)

# Модель 3: ARIMA на лог-данных (с коррекцией смещения)
ts_data_log <- ts(log(index_data$INDEX), start = min(index_data$YEAR), frequency = 1)
train_log <- window(ts_data_log, end = max(index_data$YEAR) - 1)
fit_arima_log <- forecast::auto.arima(train_log, seasonal = FALSE)
fc_arima_log_1 <- forecast::forecast(fit_arima_log, h = 1)

# Модель 4: ETS на лог-шкале (lambda = 0)
fit_ets <- forecast::ets(train, lambda = 0)
fc_ets_1 <- forecast::forecast(fit_ets, h = 1)

# Модель 5: ARIMA на логит-преобразовании, ограниченном K (с насыщением)
# Преобразуем y -> z = logit(y/K). При прогнозе обратно: y = K * inv_logit(z)
train_K <- as.numeric(train)
K_train <- max(K, 1.05 * max(train_K))

# Защитим от значений 0/K и 1
eps <- 1e-6
p_train <- pmin(pmax(train_K / K_train, eps), 1 - eps)
z_train <- ts(logit(p_train), start = start(train), frequency = frequency(train))
fit_arima_logit <- forecast::auto.arima(z_train, seasonal = FALSE)
fc_arima_logit_1 <- forecast::forecast(fit_arima_logit, h = 1)

# ========== Прогнозы на 2024 (оценка на тесте) ==========
# Наивная
naive_fc_1 <- forecast::forecast(forecast::naive(train, h = 1), h = 1)

# RW with drift
rwd_fc_1 <- forecast::forecast(forecast::rwf(train, h = 1, drift = TRUE), h = 1)

# ARIMA (log) — обратное преобразование с коррекцией смещения
arima_log_fc_1_mean <- exp(fc_arima_log_1$mean + 0.5 * fit_arima_log$sigma2)
# 80/95 интервалы переводим экспонентой (приближение)
arima_log_fc_1_lower <- exp(fc_arima_log_1$lower)
arima_log_fc_1_upper <- exp(fc_arima_log_1$upper)

# ETS (lambda=0) уже возвращает на исходной шкале
ets_fc_1_mean <- as.numeric(fc_ets_1$mean)
ets_fc_1_lower <- as.matrix(fc_ets_1$lower)
ets_fc_1_upper <- as.matrix(fc_ets_1$upper)

# ARIMA (logit with K)
arima_logit_fc_1_mean <- K_train * inv_logit(fc_arima_logit_1$mean)
arima_logit_fc_1_lower <- K_train * inv_logit(fc_arima_logit_1$lower)
arima_logit_fc_1_upper <- K_train * inv_logit(fc_arima_logit_1$upper)

# Применяем биологические ограничения (0..K) и ограничение годового прироста
last_actual <- as.numeric(tail(train, 1))

apply_bio_constraints <- function(mean_vec, lower_mat, upper_mat, last_actual, K, growth_cap_enable, g_bounds) {
  # mean_vec: numeric length h
  # lower_mat/upper_mat: матрицы с колонками 80%, 95%
  h <- length(mean_vec)
  if (is.null(dim(lower_mat))) lower_mat <- matrix(lower_mat, nrow = h)
  if (is.null(dim(upper_mat))) upper_mat <- matrix(upper_mat, nrow = h)

  # Ограничения 0..K
  mean_vec <- clip_to_bounds(mean_vec, 0, K)
  lower_mat <- apply(lower_mat, 2, clip_to_bounds, lower = 0, upper = K)
  upper_mat <- apply(upper_mat, 2, clip_to_bounds, lower = 0, upper = K)

  if (growth_cap_enable) {
    mean_vec <- adjust_by_growth_cap(mean_vec, last_actual, g_bounds)
    for (j in seq_len(ncol(lower_mat))) {
      lower_mat[, j] <- adjust_by_growth_cap(lower_mat[, j], last_actual, g_bounds)
      upper_mat[, j] <- adjust_by_growth_cap(upper_mat[, j], last_actual, g_bounds)
    }
  }
  list(mean = mean_vec, lower = lower_mat, upper = upper_mat)
}

# Собираем прогнозы h=1 для оценки на тесте
collect_fc1 <- function(mean_val, lower_mat, upper_mat) {
  data.frame(
    mean = as.numeric(mean_val[1]),
    lower80 = as.numeric(lower_mat[1, 1]),
    upper80 = as.numeric(upper_mat[1, 1]),
    lower95 = as.numeric(lower_mat[1, 2]),
    upper95 = as.numeric(upper_mat[1, 2])
  )
}

# Приводим все прогнозы к общей форме с ограничениями
naive_c <- apply_bio_constraints(
  mean_vec = as.numeric(naive_fc_1$mean),
  lower_mat = as.matrix(naive_fc_1$lower),
  upper_mat = as.matrix(naive_fc_1$upper),
  last_actual = last_actual, K = K, 
  growth_cap_enable = growth_cap_enable, g_bounds = g_bounds
)

rwd_c <- apply_bio_constraints(
  mean_vec = as.numeric(rwd_fc_1$mean),
  lower_mat = as.matrix(rwd_fc_1$lower),
  upper_mat = as.matrix(rwd_fc_1$upper),
  last_actual = last_actual, K = K, 
  growth_cap_enable = growth_cap_enable, g_bounds = g_bounds
)

arima_log_c <- apply_bio_constraints(
  mean_vec = as.numeric(arima_log_fc_1_mean),
  lower_mat = as.matrix(arima_log_fc_1_lower),
  upper_mat = as.matrix(arima_log_fc_1_upper),
  last_actual = last_actual, K = K, 
  growth_cap_enable = growth_cap_enable, g_bounds = g_bounds
)

ets_c <- apply_bio_constraints(
  mean_vec = as.numeric(ets_fc_1_mean),
  lower_mat = ets_fc_1_lower,
  upper_mat = ets_fc_1_upper,
  last_actual = last_actual, K = K, 
  growth_cap_enable = growth_cap_enable, g_bounds = g_bounds
)

arima_logit_c <- apply_bio_constraints(
  mean_vec = as.numeric(arima_logit_fc_1_mean),
  lower_mat = as.matrix(arima_logit_fc_1_lower),
  upper_mat = as.matrix(arima_logit_fc_1_upper),
  last_actual = last_actual, K = K, 
  growth_cap_enable = growth_cap_enable, g_bounds = g_bounds
)

# Метрика sMAPE
calc_smape <- function(actual, forecast) {
  200 * mean(abs(actual - forecast) / (abs(actual) + abs(forecast)), na.rm = TRUE)
}

# Сравнение на тесте 2024
results <- tibble::tibble(
  Модель = c("Наивная", "RW с дрейфом", "ARIMA (лог)", "ETS (лог)", "ARIMA (логит, K)"),
  Прогноз_2024 = c(naive_c$mean[1], rwd_c$mean[1], arima_log_c$mean[1], ets_c$mean[1], arima_logit_c$mean[1]) %>% as.numeric(),
  Факт_2024 = rep(actual_2024, 5) %>% as.numeric()
) %>% 
  dplyr::mutate(
    Ошибка = abs(Прогноз_2024 - Факт_2024),
    sMAPE = c(
      calc_smape(actual_2024, naive_c$mean[1]),
      calc_smape(actual_2024, rwd_c$mean[1]),
      calc_smape(actual_2024, arima_log_c$mean[1]),
      calc_smape(actual_2024, ets_c$mean[1]),
      calc_smape(actual_2024, arima_logit_c$mean[1])
    )
  )

cat("\nСравнение моделей на тестовых данных (2024 год):\n")
print(results)
readr::write_csv(results, file.path(outputs_dir, "model_comparison_test2024.csv"))

# =========================
# 4. КРОСС-ВАЛИДАЦИЯ (tsCV, h=1)
# =========================
# Определяем функции-прогнозёры для tsCV (возвращают вектор длины h)
fc_fun_naive <- function(y, h) { forecast::forecast(forecast::naive(y, h = h), h = h)$mean }
fc_fun_rwd   <- function(y, h) { forecast::forecast(forecast::rwf(y, h = h, drift = TRUE), h = h)$mean }
fc_fun_ets   <- function(y, h) { forecast::forecast(forecast::ets(y, lambda = 0), h = h)$mean }
fc_fun_arima_log <- function(y, h) {
  yl <- log(y)
  fit <- forecast::auto.arima(yl, seasonal = FALSE)
  fc <- forecast::forecast(fit, h = h)
  as.numeric(exp(fc$mean + 0.5 * fit$sigma2))
}
fc_fun_arima_logit <- function(y, h) {
  K_loc <- max(1.5 * max(y), max(y) * 1.05)
  eps <- 1e-6
  p <- pmin(pmax(y / K_loc, eps), 1 - eps)
  z <- log(p/(1 - p))
  fit <- forecast::auto.arima(ts(z, start = start(ts(y)), frequency = frequency(ts(y))), seasonal = FALSE)
  fc <- forecast::forecast(fit, h = h)
  as.numeric(K_loc * inv_logit(fc$mean))
}

cv_smape <- function(y, f_fun, h = 1) {
  e <- forecast::tsCV(y, f_fun, h = h)
  # e содержит ошибки прогнозов: e[t] = y[t] - yhat[t]
  # Преобразуем в sMAPE, исключая NA
  idx <- which(!is.na(e[, h]))
  if (length(idx) == 0) return(NA_real_)
  y_act <- y[idx]
  y_hat <- y_act - e[idx, h]
  calc_smape(y_act, y_hat)
}

cv_table <- tibble::tibble(
  Модель = c("Наивная", "RW с дрейфом", "ARIMA (лог)", "ETS (лог)", "ARIMA (логит, K)"),
  sMAPE_CV_h1 = c(
    cv_smape(ts_data, fc_fun_naive, h = 1),
    cv_smape(ts_data, fc_fun_rwd,   h = 1),
    cv_smape(ts_data, fc_fun_arima_log, h = 1),
    cv_smape(ts_data, fc_fun_ets,   h = 1),
    cv_smape(ts_data, fc_fun_arima_logit, h = 1)
  )
)

cat("\nКросс-валидация (h=1), sMAPE:\n")
print(cv_table)
readr::write_csv(cv_table, file.path(outputs_dir, "cv_smape_h1.csv"))

# =========================
# 5. ВЫБОР ЛУЧШЕЙ МОДЕЛИ И ФИНАЛЬНЫЙ ПРОГНОЗ (2025-2026)
# =========================
# Выбираем по тестовой sMAPE; при равенстве — по CV
i_best <- which.min(results$sMAPE)
best_model_name <- results$Модель[i_best]
if (length(best_model_name) == 0 || is.na(best_model_name)) best_model_name <- "ARIMA (логит, K)"

cat("\nЛучшая модель:", best_model_name, "(sMAPE тест=", round(results$sMAPE[i_best], 2), "%)\n")

# Получаем финальный прогноз h=2 на полной выборке
get_final_fc <- function(model_name, ts_data, K, growth_cap_enable, g_bounds) {
  last_actual <- as.numeric(tail(ts_data, 1))
  if (model_name == "Наивная") {
    fit <- forecast::naive(ts_data, h = 2)
    mean_vec <- as.numeric(fit$mean)
    lower_mat <- as.matrix(fit$lower)
    upper_mat <- as.matrix(fit$upper)
  } else if (model_name == "RW с дрейфом") {
    fit <- forecast::rwf(ts_data, h = 2, drift = TRUE)
    mean_vec <- as.numeric(fit$mean)
    lower_mat <- as.matrix(fit$lower)
    upper_mat <- as.matrix(fit$upper)
  } else if (model_name == "ARIMA (лог)") {
    fit <- forecast::auto.arima(log(ts_data), seasonal = FALSE)
    fc  <- forecast::forecast(fit, h = 2)
    mean_vec <- as.numeric(exp(fc$mean + 0.5 * fit$sigma2))
    lower_mat <- as.matrix(exp(fc$lower))
    upper_mat <- as.matrix(exp(fc$upper))
  } else if (model_name == "ETS (лог)") {
    fit <- forecast::ets(ts_data, lambda = 0)
    fc  <- forecast::forecast(fit, h = 2)
    mean_vec <- as.numeric(fc$mean)
    lower_mat <- as.matrix(fc$lower)
    upper_mat <- as.matrix(fc$upper)
  } else { # "ARIMA (логит, K)"
    K_loc <- max(K, 1.05 * max(ts_data))
    eps <- 1e-6
    p <- pmin(pmax(as.numeric(ts_data) / K_loc, eps), 1 - eps)
    z <- ts(log(p/(1 - p)), start = start(ts_data), frequency = frequency(ts_data))
    fit <- forecast::auto.arima(z, seasonal = FALSE)
    fc  <- forecast::forecast(fit, h = 2)
    mean_vec <- as.numeric(K_loc * inv_logit(fc$mean))
    lower_mat <- as.matrix(K_loc * inv_logit(fc$lower))
    upper_mat <- as.matrix(K_loc * inv_logit(fc$upper))
  }

  cons <- apply_bio_constraints(mean_vec, lower_mat, upper_mat, last_actual, K, growth_cap_enable, g_bounds)
  list(mean = cons$mean, lower = cons$lower, upper = cons$upper)
}

final_fc <- get_final_fc(best_model_name, ts_data, K = K, 
                         growth_cap_enable = growth_cap_enable, g_bounds = g_bounds)

forecast_details <- tibble::tibble(
  Год = (max(index_data$YEAR) + 1):(max(index_data$YEAR) + 2),
  Прогноз = round(final_fc$mean/1e6, 2),
  Нижняя_80 = round(final_fc$lower[, 1]/1e6, 2),
  Верхняя_80 = round(final_fc$upper[, 1]/1e6, 2),
  Нижняя_95 = round(final_fc$lower[, 2]/1e6, 2),
  Верхняя_95 = round(final_fc$upper[, 2]/1e6, 2)
)

cat("\nПрогнозные значения (млн экз.):\n")
print(forecast_details)
readr::write_csv(forecast_details, file.path(outputs_dir, "final_forecast_2025_2026.csv"))

# Визуализация прогноза с доверительными интервалами
hist_df <- index_data %>% dplyr::select(YEAR, INDEX) %>% dplyr::mutate(phase = "История")
fc_df <- tibble::tibble(
  YEAR = forecast_details$Год,
  mean = final_fc$mean,
  lower80 = final_fc$lower[, 1], upper80 = final_fc$upper[, 1],
  lower95 = final_fc$lower[, 2], upper95 = final_fc$upper[, 2],
  phase = "Прогноз"
)

p_fc <- ggplot() +
  geom_ribbon(data = fc_df, aes(x = YEAR, ymin = lower95/1e6, ymax = upper95/1e6), fill = "#7fb3d5", alpha = 0.4) +
  geom_ribbon(data = fc_df, aes(x = YEAR, ymin = lower80/1e6, ymax = upper80/1e6), fill = "#2e86c1", alpha = 0.4) +
  geom_line(data = hist_df, aes(x = YEAR, y = INDEX/1e6, color = phase), linewidth = 1) +
  geom_point(data = hist_df, aes(x = YEAR, y = INDEX/1e6, color = phase), size = 2.5) +
  geom_line(data = fc_df, aes(x = YEAR, y = mean/1e6, color = phase), linewidth = 1, linetype = "dashed") +
  geom_point(data = fc_df, aes(x = YEAR, y = mean/1e6, color = phase), size = 2.5, shape = 17) +
  scale_color_manual(values = c("История" = "steelblue", "Прогноз" = "#b03a2e")) +
  labs(title = "Прогноз индекса обилия камчатского краба на 2025-2026 гг.",
       x = "Год", y = "Индекс обилия, млн экз.", color = "Этап") +
  theme_minimal()

ggsave(filename = file.path(outputs_dir, "forecast_plot.png"), plot = p_fc, width = 9, height = 5, dpi = 150)

# =========================
# 6. ДИАГНОСТИКА ОСТАТКОВ
# =========================
# Сохраним графики диагностики для ключевых моделей
safe_checkres <- function(fit, file) {
  try({
    grDevices::png(filename = file, width = 900, height = 700)
    forecast::checkresiduals(fit)
    grDevices::dev.off()
  }, silent = TRUE)
}

# Диагностика для моделей, которые имеют объект fit
safe_checkres(fit_naive, file.path(outputs_dir, "diag_naive.png"))
safe_checkres(forecast::rwf(train, h = 1, drift = TRUE), file.path(outputs_dir, "diag_rwd.png"))
safe_checkres(fit_ets, file.path(outputs_dir, "diag_ets.png"))
# Для ARIMA лог-шкалы
safe_checkres(fit_arima_log, file.path(outputs_dir, "diag_arima_log.png"))
# Для ARIMA логит
safe_checkres(fit_arima_logit, file.path(outputs_dir, "diag_arima_logit.png"))

# =========================
# 7. РИСК-АНАЛИЗ И СЦЕНАРИИ ЭКСПЛУАТАЦИИ
# =========================
cat("\n=== РИСК-АНАЛИЗ И РЕКОМЕНДАЦИИ ===\n")
current_level <- index_data$INDEX[nrow(index_data)]
trend_4yr <- (index_data$INDEX[6] - index_data$INDEX[2]) / index_data$INDEX[2] * 100

cat("Текущий уровень запаса (", max(index_data$YEAR), "): ", round(current_level/1e6, 2), " млн экз.\n", sep = "")
cat("Изменение за последние 4 года:", round(trend_4yr, 1), "%\n")

if (trend_4yr < -20) {
  status_text <- "КРИТИЧЕСКИЙ - значительное снижение запаса"
  recs <- c("Ужесточение промысловой нагрузки",
            "Усиление мониторинга",
            "Исследование причин снижения")
} else if (trend_4yr < -10) {
  status_text <- "ТРЕВОЖНЫЙ - умеренное снижение"
  recs <- c("Стабилизация промысла", "Мониторинг ключевых параметров")
} else {
  status_text <- "СТАБИЛЬНЫЙ"
  recs <- c("Поддержание текущего уровня промысла")
}

cat("СТАТУС:", status_text, "\n")
cat("РЕКОМЕНДАЦИИ:\n-", paste(recs, collapse = "\n- "), "\n")

# Вероятностный анализ на основе приближённого распределения прогноза
# Оценим sd из 95% ДИ: sd ≈ (U95 - L95) / (2*1.96)
calc_sd_from_ci <- function(l95, u95) {
  as.numeric((u95 - l95) / (2 * 1.96))
}

h <- 2
means <- final_fc$mean
sd_approx <- mapply(calc_sd_from_ci, final_fc$lower[, 2], final_fc$upper[, 2])
nsim <- 10000

sim_mat <- matrix(NA_real_, nrow = nsim, ncol = h)
for (i in 1:h) {
  sim_i <- rnorm(nsim, mean = means[i], sd = pmax(sd_approx[i], 1e-9))
  sim_i <- clip_to_bounds(sim_i, 0, K)
  if (growth_cap_enable) {
    base_prev <- if (i == 1) as.numeric(tail(ts_data, 1)) else sim_mat[, i - 1]
    # применяем покомпонентно ограничение роста относительно base_prev
    gr <- (sim_i - base_prev) / base_prev
    gr <- pmin(pmax(gr, g_bounds[1]), g_bounds[2])
    sim_i <- base_prev * (1 + gr)
  }
  sim_mat[, i] <- sim_i
}

# Пороговые значения
Bcrit_abs <- 1e6           # абсолютный критический уровень (пример)
Blim_relK <- 0.3 * K       # порог как доля K

risk_tbl <- tibble::tibble(
  Год = forecast_details$Год,
  P_below_abs = colMeans(sim_mat < Bcrit_abs),
  P_below_0.3K = colMeans(sim_mat < Blim_relK),
  E_mean = colMeans(sim_mat),
  E_median = apply(sim_mat, 2, median)
) %>% dplyr::mutate(
  E_mean_mln = round(E_mean/1e6, 2),
  E_median_mln = round(E_median/1e6, 2)
)

readr::write_csv(risk_tbl, file.path(outputs_dir, "risk_analysis.csv"))

cat("\nВероятностный анализ (2025):\n")
cat("- P(Index < ", round(Bcrit_abs/1e6, 2), " млн) = ", round(risk_tbl$P_below_abs[1], 3), "\n", sep = "")
cat("- P(Index < 0.3*K) = ", round(risk_tbl$P_below_0.3K[1], 3), "\n", sep = "")

# Сценарии эксплуатационной нагрузки (демонстрация): пост-эксплуатационный индекс = I * (1 - u)
exploit_rates <- c(0, 0.1, 0.2, 0.3)
scen_list <- list()
for (u in exploit_rates) {
  sim_post <- sim_mat * (1 - u)
  scen_list[[as.character(u)]] <- tibble::tibble(
    Год = forecast_details$Год,
    u = u,
    P_below_abs = colMeans(sim_post < Bcrit_abs),
    P_below_0.3K = colMeans(sim_post < Blim_relK),
    mean_post = colMeans(sim_post),
    median_post = apply(sim_post, 2, median)
  )
}
scenario_tbl <- dplyr::bind_rows(scen_list) %>%
  dplyr::mutate(mean_post_mln = round(mean_post/1e6, 2),
                median_post_mln = round(median_post/1e6, 2))

readr::write_csv(scenario_tbl, file.path(outputs_dir, "scenario_risk_analysis.csv"))

# Визуализация сценариев для 2025
scen2025 <- scenario_tbl %>% dplyr::filter(Год == (max(index_data$YEAR) + 1))
p_scen <- ggplot(scen2025, aes(x = factor(u), y = mean_post/1e6)) +
  geom_col(fill = "#3498db") +
  geom_errorbar(aes(ymin = (mean_post - 1.96*sd_approx[1])/1e6,
                    ymax = (mean_post + 1.96*sd_approx[1])/1e6), width = 0.2) +
  geom_hline(yintercept = Bcrit_abs/1e6, linetype = "dashed", color = "#b03a2e") +
  geom_hline(yintercept = Blim_relK/1e6, linetype = "dotted", color = "#7d3c98") +
  labs(title = "Сценарии эксплуатации: ожидаемый индекс (2025)",
       x = "Доля изъятия u", y = "Индекс, млн экз.") +
  theme_minimal()

ggsave(filename = file.path(outputs_dir, "scenario_bar_2025.png"), plot = p_scen, width = 7, height = 4.5, dpi = 150)

# =========================
# 8. ФОРМИРОВАНИЕ И ЭКСПОРТ ОТЧЁТА
# =========================
final_report <- list(
  исходные_данные = index_data,
  сравнение_моделей_тест2024 = results,
  cv_smape_h1 = cv_table,
  лучшая_модель = best_model_name,
  прогноз_2025_2026 = forecast_details,
  риск_анализ = risk_tbl,
  сценарии = scenario_tbl,
  параметры = list(K = K, growth_cap_enable = growth_cap_enable, g_bounds = g_bounds),
  рекомендации = list(статус = status_text, меры = recs)
)

# JSON, RDS, TXT-summary
jsonlite::write_json(final_report, path = file.path(outputs_dir, "final_report.json"), pretty = TRUE, auto_unbox = TRUE)
saveRDS(final_report, file = file.path(outputs_dir, "final_report.rds"))

summary_lines <- c(
  "=== АНАЛИЗ ЗАВЕРШЕН ===",
  paste0("Лучшая модель: ", best_model_name),
  paste0("Прогноз на ", forecast_details$Год[1], ": ", forecast_details$Прогноз[1], " млн экз."),
  paste0("80% ДИ: ", forecast_details$Нижняя_80[1], " - ", forecast_details$Верхняя_80[1], " млн экз."),
  paste0("P(Index < ", round(Bcrit_abs/1e6, 2), " млн) в ", forecast_details$Год[1], ": ", round(risk_tbl$P_below_abs[1], 3)),
  paste0("P(Index < 0.3*K) в ", forecast_details$Год[1], ": ", round(risk_tbl$P_below_0.3K[1], 3))
)

writeLines(summary_lines, con = file.path(outputs_dir, "summary.txt"))

# Сохраним sessionInfo
utils::capture.output(utils::sessionInfo(), file = file.path(outputs_dir, "sessionInfo.txt"))

# Дублируем ключевую информацию в консоль
cat("\n=== АНАЛИЗ ЗАВЕРШЕН ===\n")
cat("Лучшая модель:", best_model_name, "\n")
cat("Прогноз на ", forecast_details$Год[1], ": ", forecast_details$Прогноз[1], " млн экз.\n", sep = "")
cat("Диапазон неопределенности (80% ДИ): ", forecast_details$Нижняя_80[1], " - ", forecast_details$Верхняя_80[1], " млн экз.\n", sep = "")
