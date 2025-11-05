# ПРАКТИЧЕСКОЕ ЗАНЯТИЕ: ОСНОВЫ ПРОГНОЗИРОВАНИЯ
# УЛУЧШЕННАЯ ВЕРСИЯ — НОВЫЙ РЯД (BESS), РОБАСТНОСТЬ, БИО-ОГРАНИЧЕНИЯ, CV, СЦЕНАРИИ

# Очистка рабочей среды (опционально)
rm(list = ls())
set.seed(123)

# ===============
# 0. ПАКЕТЫ
# ===============
required_pkgs <- c("tidyverse", "forecast", "tseries", "ggplot2", "jsonlite", "readr")
ensure_packages <- function(pkgs) {
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
      try({ install.packages(p, repos = "https://cloud.r-project.org", quiet = TRUE) }, silent = TRUE)
    }
    suppressPackageStartupMessages(library(p, character.only = TRUE))
  }
}
ensure_packages(required_pkgs)

# ===============
# 1. ДАННЫЕ
# ===============
# 2.1 Вектор лет наблюдений
Year <- 2017:2025
# 2.4 Индекс BESS (научная съемка)
SURVEY <- c(29.4, 56.8, 44.3, 39.9, 42.0, 23.1, 26.1, 18.9, 21.7)

stopifnot(length(Year) == length(SURVEY))

index_data <- tibble::tibble(
  YEAR = Year,
  INDEX = SURVEY
)

cat("=== АНАЛИЗ ДИНАМИКИ ИНДЕКСА BESS ===\n")
cat("Период наблюдений:", min(index_data$YEAR), "-", max(index_data$YEAR), "\n")
cat("Изменение индекса за период:", 
    round((tail(index_data$INDEX, 1) - head(index_data$INDEX, 1)) / head(index_data$INDEX, 1) * 100, 1), "%\n")

# Временной ряд (годовая частота)
ts_data <- ts(index_data$INDEX, start = min(index_data$YEAR), frequency = 1)

# Директория вывода (относительно текущей рабочей директории, без setwd)
outputs_dir <- file.path(getwd(), "outputs")
if (!dir.exists(outputs_dir)) dir.create(outputs_dir, recursive = TRUE)

# Базовая визуализация (в условных единицах индекса)
p_hist <- ggplot(index_data, aes(x = YEAR, y = INDEX)) +
  geom_line(linewidth = 1, color = "steelblue") +
  geom_point(size = 3, color = "navy") +
  geom_text(aes(label = round(INDEX, 1)), vjust = -1, size = 3.5) +
  labs(title = "Динамика индекса BESS (научная съемка)",
       subtitle = paste0(min(index_data$YEAR), "-", max(index_data$YEAR)),
       x = "Год", y = "Индекс (усл. ед.)") +
  theme_minimal()

ggsave(filename = file.path(outputs_dir, "historical_index.png"), plot = p_hist, 
       width = 8, height = 4.5, dpi = 150)

# ===============
# 2. СТАЦИОНАРНОСТЬ
# ===============
cat("\n=== АНАЛИЗ СТАЦИОНАРНОСТИ РЯДА ===\n")
cat("Длина ряда:", length(ts_data), "наблюдений\n")
cat("Среднее значение:", round(mean(ts_data), 2), "\n")
cat("Стандартное отклонение:", round(sd(ts_data), 2), "\n")

if (mean(diff(ts_data)) < 0) {
  cat("Визуальный анализ: нисходящий тренд\n")
} else {
  cat("Визуальный анализ: восходящий тренд\n")
}

# ADF-тест (робастная печать p-value; NaN допустим)
safe_adf <- function(x) {
  out <- try(tseries::adf.test(x), silent = TRUE)
  if (inherits(out, "try-error")) return(NULL) else return(out)
}
adf_res <- safe_adf(ts_data)
if (!is.null(adf_res) && is.numeric(adf_res$p.value)) {
  cat("ADF тест: p-value =", ifelse(is.finite(adf_res$p.value), signif(adf_res$p.value, 3), NA), "(осторожная интерпретация)\n")
} else {
  cat("ADF тест: p-value недоступно (малая выборка), пропускаем вывод.\n")
}

# ===============
# 3. МОДЕЛИ С БИОЛОГИЧЕСКИМИ ОГРАНИЧЕНИЯМИ
# ===============
# Разделение на train (до 2024) и test (2025)
train <- window(ts_data, end = max(index_data$YEAR) - 1)
test  <- window(ts_data, start = max(index_data$YEAR))
actual_test <- as.numeric(test)

# Параметры ограничений
K <- 1.5 * max(ts_data, na.rm = TRUE)      # Носительная способность в шкале индекса
growth_cap_enable <- TRUE                  # Ограничение темпа роста/падения
# Например, в год не ниже -60% и не выше +50% относительно прошлого уровня
g_bounds <- c(-0.6, 0.5)

db_clip <- function(x, lower = 0, upper = Inf) {
  x <- pmax(x, lower)
  x <- pmin(x, upper)
  x
}

adjust_by_growth_cap <- function(path, last_actual, g_bounds) {
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

logit <- function(p) log(p/(1 - p))
inv_logit <- function(z) 1/(1 + exp(-z))

# Вспомогательная: гарантировать форму матрицы ДИ (h x 2)
ensure_ci_matrix <- function(ci_obj, h) {
  if (is.null(ci_obj)) {
    return(matrix(NA_real_, nrow = h, ncol = 2))
  }
  m <- as.matrix(ci_obj)
  if (nrow(m) == 0 && length(m) > 0) {
    m <- matrix(m, nrow = h, byrow = TRUE)
  }
  if (ncol(m) == 1) {
    m <- cbind(m[, 1], m[, 1])
  }
  if (nrow(m) != h) {
    # попытка привести к размеру h
    m <- matrix(m[seq_len(min(length(m), h*2))], nrow = h, byrow = TRUE)
  }
  colnames(m) <- c("80%", "95%")
  m
}

apply_bio_constraints <- function(mean_vec, lower_mat, upper_mat, last_actual, K, growth_cap_enable, g_bounds) {
  h <- length(mean_vec)
  lower_mat <- ensure_ci_matrix(lower_mat, h)
  upper_mat <- ensure_ci_matrix(upper_mat, h)

  # Границы 0..K
  mean_vec <- db_clip(mean_vec, 0, K)
  lower_mat <- apply(lower_mat, 2, db_clip, lower = 0, upper = K)
  upper_mat <- apply(upper_mat, 2, db_clip, lower = 0, upper = K)

  # Ограничение годового изменения
  if (growth_cap_enable) {
    mean_vec <- adjust_by_growth_cap(mean_vec, last_actual, g_bounds)
    for (j in seq_len(ncol(lower_mat))) {
      lower_mat[, j] <- adjust_by_growth_cap(lower_mat[, j], last_actual, g_bounds)
      upper_mat[, j] <- adjust_by_growth_cap(upper_mat[, j], last_actual, g_bounds)
    }
  }
  list(mean = mean_vec, lower = lower_mat, upper = upper_mat)
}

# ---- Обучаем модели на train ----
fit_naive <- forecast::naive(train, h = 1)
fit_rwd   <- forecast::rwf(train, drift = TRUE, h = 1)

# ARIMA на лог-данных (с коррекцией смещения)
train_log <- ts(log(as.numeric(train)), start = start(train), frequency = frequency(train))
fit_arima_log <- forecast::auto.arima(train_log, seasonal = FALSE)
fc_arima_log_1 <- forecast::forecast(fit_arima_log, h = 1)

# ETS на лог-шкале (lambda=0)
fit_ets <- forecast::ets(train, lambda = 0)
fc_ets_1 <- forecast::forecast(fit_ets, h = 1)

# ARIMA на логит-данных, масштабированных K
eps <- 1e-6
p_train <- pmin(pmax(as.numeric(train) / K, eps), 1 - eps)
z_train <- ts(logit(p_train), start = start(train), frequency = frequency(train))
fit_arima_logit <- forecast::auto.arima(z_train, seasonal = FALSE)
fc_arima_logit_1 <- forecast::forecast(fit_arima_logit, h = 1)

# ---- Прогнозы на тестовый 2025 (h=1) ----
# Все объекты формируем как mean (num), lower/upper (1x2)
h1 <- 1
last_actual <- as.numeric(tail(train, 1))

# Наивная
naive_mean  <- as.numeric(fit_naive$mean)
naive_lower <- ensure_ci_matrix(fit_naive$lower, h1)
naive_upper <- ensure_ci_matrix(fit_naive$upper, h1)
naive_c <- apply_bio_constraints(naive_mean, naive_lower, naive_upper, last_actual, K, growth_cap_enable, g_bounds)

# RW с дрейфом
rwd_mean  <- as.numeric(fit_rwd$mean)
rwd_lower <- ensure_ci_matrix(fit_rwd$lower, h1)
rwd_upper <- ensure_ci_matrix(fit_rwd$upper, h1)
rwd_c <- apply_bio_constraints(rwd_mean, rwd_lower, rwd_upper, last_actual, K, growth_cap_enable, g_bounds)

# ARIMA (лог) -> исходная шкала с коррекцией смещения
arima_log_mean  <- as.numeric(exp(fc_arima_log_1$mean + 0.5 * fit_arima_log$sigma2))
arima_log_lower <- ensure_ci_matrix(exp(fc_arima_log_1$lower), h1)
arima_log_upper <- ensure_ci_matrix(exp(fc_arima_log_1$upper), h1)
arima_log_c <- apply_bio_constraints(arima_log_mean, arima_log_lower, arima_log_upper, last_actual, K, growth_cap_enable, g_bounds)

# ETS (lambda=0 уже на исходной шкале)
ets_mean  <- as.numeric(fc_ets_1$mean)
ets_lower <- ensure_ci_matrix(fc_ets_1$lower, h1)
ets_upper <- ensure_ci_matrix(fc_ets_1$upper, h1)
ets_c <- apply_bio_constraints(ets_mean, ets_lower, ets_upper, last_actual, K, growth_cap_enable, g_bounds)

# ARIMA (логит, K)
arima_logit_mean  <- as.numeric(K * inv_logit(fc_arima_logit_1$mean))
arima_logit_lower <- ensure_ci_matrix(K * inv_logit(fc_arima_logit_1$lower), h1)
arima_logit_upper <- ensure_ci_matrix(K * inv_logit(fc_arima_logit_1$upper), h1)
arima_logit_c <- apply_bio_constraints(arima_logit_mean, arima_logit_lower, arima_logit_upper, last_actual, K, growth_cap_enable, g_bounds)

# Метрика sMAPE
calc_smape <- function(actual, forecast) {
  200 * mean(abs(actual - forecast) / (abs(actual) + abs(forecast)), na.rm = TRUE)
}

results <- tibble::tibble(
  Модель = c("Наивная", "RW с дрейфом", "ARIMA (лог)", "ETS (лог)", "ARIMA (логит, K)"),
  Прогноз_тест = c(naive_c$mean[1], rwd_c$mean[1], arima_log_c$mean[1], ets_c$mean[1], arima_logit_c$mean[1]) %>% as.numeric(),
  Факт = rep(actual_test, 5) %>% as.numeric()
) %>% dplyr::mutate(
  Ошибка = abs(Прогноз_тест - Факт),
  sMAPE = c(
    calc_smape(actual_test, naive_c$mean[1]),
    calc_smape(actual_test, rwd_c$mean[1]),
    calc_smape(actual_test, arima_log_c$mean[1]),
    calc_smape(actual_test, ets_c$mean[1]),
    calc_smape(actual_test, arima_logit_c$mean[1])
  )
)

cat("\nСравнение моделей на тесте (", max(index_data$YEAR), "):\n", sep = "")
print(results)
readr::write_csv(results, file.path(outputs_dir, "model_comparison_test.csv"))

# ===============
# 4. КРОСС-ВАЛИДАЦИЯ (tsCV, h=1)
# ===============
# Безопасные прогнозёры для tsCV (возвращают numeric длины h; при ошибке — NA)
safe_fc_wrapper <- function(fun) {
  force(fun)
  function(y, h) {
    out <- try(fun(y, h), silent = TRUE)
    if (inherits(out, "try-error") || is.null(out)) return(rep(NA_real_, h))
    as.numeric(out)
  }
}

fc_fun_naive <- safe_fc_wrapper(function(y, h) forecast::naive(y, h = h)$mean)
fc_fun_rwd   <- safe_fc_wrapper(function(y, h) forecast::rwf(y, h = h, drift = TRUE)$mean)
fc_fun_ets   <- safe_fc_wrapper(function(y, h) forecast::forecast(forecast::ets(y, lambda = 0), h = h)$mean)
fc_fun_arima_log <- safe_fc_wrapper(function(y, h) {
  yl <- log(as.numeric(y))
  fit <- forecast::auto.arima(yl, seasonal = FALSE)
  fc <- forecast::forecast(fit, h = h)
  exp(as.numeric(fc$mean) + 0.5 * fit$sigma2)
})
fc_fun_arima_logit <- safe_fc_wrapper(function(y, h) {
  K_loc <- max(1.5 * max(y, na.rm = TRUE), 1.05 * max(y, na.rm = TRUE))
  eps <- 1e-6
  p <- pmin(pmax(as.numeric(y) / K_loc, eps), 1 - eps)
  z <- ts(logit(p), start = start(ts(y)), frequency = frequency(ts(y)))
  fit <- forecast::auto.arima(z, seasonal = FALSE)
  fc <- forecast::forecast(fit, h = h)
  as.numeric(K_loc * inv_logit(fc$mean))
})

cv_smape <- function(y, f_fun, h = 1) {
  e <- try(forecast::tsCV(y, f_fun, h = h), silent = TRUE)
  if (inherits(e, "try-error")) return(NA_real_)
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

# ===============
# 5. ФИНАЛЬНЫЙ ПРОГНОЗ (h=2) НА ПОЛНОМ РЯДЕ => 2026–2027
# ===============
# Выбор по sMAPE теста; при равенстве — по CV
best_idx <- which.min(results$sMAPE)
best_model_name <- results$Модель[best_idx]

if (length(best_model_name) == 0 || is.na(best_model_name)) {
  # запасной вариант
  best_model_name <- "ARIMA (логит, K)"
}

cat("\nЛучшая модель:", best_model_name, "(sMAPE тест=", round(results$sMAPE[best_idx], 2), "%)\n")

get_final_fc <- function(model_name, ts_data, K, growth_cap_enable, g_bounds) {
  last_actual <- as.numeric(tail(ts_data, 1))
  h <- 2
  if (model_name == "Наивная") {
    fit <- forecast::naive(ts_data, h = h)
    mean_vec <- as.numeric(fit$mean)
    lower_mat <- ensure_ci_matrix(fit$lower, h)
    upper_mat <- ensure_ci_matrix(fit$upper, h)
  } else if (model_name == "RW с дрейфом") {
    fit <- forecast::rwf(ts_data, h = h, drift = TRUE)
    mean_vec <- as.numeric(fit$mean)
    lower_mat <- ensure_ci_matrix(fit$lower, h)
    upper_mat <- ensure_ci_matrix(fit$upper, h)
  } else if (model_name == "ARIMA (лог)") {
    fit <- forecast::auto.arima(log(as.numeric(ts_data)), seasonal = FALSE)
    fc  <- forecast::forecast(fit, h = h)
    mean_vec <- as.numeric(exp(fc$mean + 0.5 * fit$sigma2))
    lower_mat <- ensure_ci_matrix(exp(fc$lower), h)
    upper_mat <- ensure_ci_matrix(exp(fc$upper), h)
  } else if (model_name == "ETS (лог)") {
    fit <- forecast::ets(ts_data, lambda = 0)
    fc  <- forecast::forecast(fit, h = h)
    mean_vec <- as.numeric(fc$mean)
    lower_mat <- ensure_ci_matrix(fc$lower, h)
    upper_mat <- ensure_ci_matrix(fc$upper, h)
  } else { # "ARIMA (логит, K)"
    eps <- 1e-6
    p <- pmin(pmax(as.numeric(ts_data) / K, eps), 1 - eps)
    z <- ts(logit(p), start = start(ts_data), frequency = frequency(ts_data))
    fit <- forecast::auto.arima(z, seasonal = FALSE)
    fc  <- forecast::forecast(fit, h = h)
    mean_vec <- as.numeric(K * inv_logit(fc$mean))
    lower_mat <- ensure_ci_matrix(K * inv_logit(fc$lower), h)
    upper_mat <- ensure_ci_matrix(K * inv_logit(fc$upper), h)
  }
  apply_bio_constraints(mean_vec, lower_mat, upper_mat, last_actual, K, growth_cap_enable, g_bounds)
}

final_fc <- get_final_fc(best_model_name, ts_data, K = K, growth_cap_enable = growth_cap_enable, g_bounds = g_bounds)

forecast_years <- (max(index_data$YEAR) + 1):(max(index_data$YEAR) + 2)
forecast_details <- tibble::tibble(
  Год = forecast_years,
  Прогноз = round(final_fc$mean, 2),
  Нижняя_80 = round(final_fc$lower[, 1], 2),
  Верхняя_80 = round(final_fc$upper[, 1], 2),
  Нижняя_95 = round(final_fc$lower[, 2], 2),
  Верхняя_95 = round(final_fc$upper[, 2], 2)
)

cat("\nПрогнозные значения (индекс, усл. ед.):\n")
print(forecast_details)
readr::write_csv(forecast_details, file.path(outputs_dir, "final_forecast.csv"))

# Визуализация прогноза
hist_df <- index_data %>% dplyr::select(YEAR, INDEX) %>% dplyr::mutate(phase = "История")
fc_df <- tibble::tibble(
  YEAR = forecast_years,
  mean = final_fc$mean,
  lower80 = final_fc$lower[, 1], upper80 = final_fc$upper[, 1],
  lower95 = final_fc$lower[, 2], upper95 = final_fc$upper[, 2],
  phase = "Прогноз"
)

p_fc <- ggplot() +
  geom_ribbon(data = fc_df, aes(x = YEAR, ymin = lower95, ymax = upper95), fill = "#7fb3d5", alpha = 0.4) +
  geom_ribbon(data = fc_df, aes(x = YEAR, ymin = lower80, ymax = upper80), fill = "#2e86c1", alpha = 0.4) +
  geom_line(data = hist_df, aes(x = YEAR, y = INDEX, color = phase), linewidth = 1) +
  geom_point(data = hist_df, aes(x = YEAR, y = INDEX, color = phase), size = 2.5) +
  geom_line(data = fc_df, aes(x = YEAR, y = mean, color = phase), linewidth = 1, linetype = "dashed") +
  geom_point(data = fc_df, aes(x = YEAR, y = mean, color = phase), size = 2.5, shape = 17) +
  scale_color_manual(values = c("История" = "steelblue", "Прогноз" = "#b03a2e")) +
  labs(title = paste0("Прогноз индекса BESS на ", forecast_years[1], "-", forecast_years[2]),
       x = "Год", y = "Индекс (усл. ед.)", color = "Этап") +
  theme_minimal()

ggsave(filename = file.path(outputs_dir, "forecast_plot.png"), plot = p_fc, width = 9, height = 5, dpi = 150)

# ===============
# 6. ДИАГНОСТИКА ОСТАТКОВ
# ===============
safe_checkres <- function(fit, file) {
  try({
    grDevices::png(filename = file, width = 900, height = 700)
    forecast::checkresiduals(fit)
    grDevices::dev.off()
  }, silent = TRUE)
}

safe_checkres(fit_naive, file.path(outputs_dir, "diag_naive.png"))
safe_checkres(forecast::rwf(train, h = 1, drift = TRUE), file.path(outputs_dir, "diag_rwd.png"))
safe_checkres(fit_ets, file.path(outputs_dir, "diag_ets.png"))
safe_checkres(fit_arima_log, file.path(outputs_dir, "diag_arima_log.png"))
safe_checkres(fit_arima_logit, file.path(outputs_dir, "diag_arima_logit.png"))

# ===============
# 7. РИСК-АНАЛИЗ И СЦЕНАРИИ
# ===============
cat("\n=== РИСК-АНАЛИЗ И РЕКОМЕНДАЦИИ ===\n")
current_level <- tail(index_data$INDEX, 1)
trend_recent <- (tail(index_data$INDEX, 1) - index_data$INDEX[which(index_data$YEAR == (max(index_data$YEAR) - 4))]) /
                 index_data$INDEX[which(index_data$YEAR == (max(index_data$YEAR) - 4))] * 100

cat("Текущий уровень индекса (", max(index_data$YEAR), "): ", round(current_level, 2), "\n", sep = "")
cat("Изменение за последние 4 года:", round(trend_recent, 1), "%\n")

if (is.finite(trend_recent) && trend_recent < -20) {
  status_text <- "КРИТИЧЕСКИЙ — значительное снижение"
  recs <- c("Ужесточение промысловой нагрузки", "Усиление мониторинга", "Исследование причин снижения")
} else if (is.finite(trend_recent) && trend_recent < -10) {
  status_text <- "ТРЕВОЖНЫЙ — умеренное снижение"
  recs <- c("Стабилизация промысла", "Мониторинг ключевых параметров")
} else {
  status_text <- "СТАБИЛЬНЫЙ"
  recs <- c("Поддержание текущего уровня промысла")
}

cat("СТАТУС:", status_text, "\n")
cat("РЕКОМЕНДАЦИИ:\n- ", paste(recs, collapse = "\n- "), "\n", sep = "")

# Имитации на основе приближённых ДИ финпрогноза
calc_sd_from_ci <- function(l95, u95) as.numeric((u95 - l95) / (2 * 1.96))
means <- final_fc$mean
sd_approx <- mapply(calc_sd_from_ci, final_fc$lower[, 2], final_fc$upper[, 2])
nsim <- 10000
h <- length(means)

sim_mat <- matrix(NA_real_, nrow = nsim, ncol = h)
for (i in seq_len(h)) {
  sim_i <- rnorm(nsim, mean = means[i], sd = pmax(sd_approx[i], 1e-9))
  sim_i <- db_clip(sim_i, 0, K)
  if (growth_cap_enable) {
    base_prev <- if (i == 1) as.numeric(tail(ts_data, 1)) else sim_mat[, i - 1]
    gr <- (sim_i - base_prev) / base_prev
    gr <- pmin(pmax(gr, g_bounds[1]), g_bounds[2])
    sim_i <- base_prev * (1 + gr)
  }
  sim_mat[, i] <- sim_i
}

# Пороговые значения (в шкале индекса)
Bcrit_abs <- 0.2 * K
Blim_relK <- 0.3 * K

risk_tbl <- tibble::tibble(
  Год = forecast_years,
  P_below_abs = colMeans(sim_mat < Bcrit_abs),
  P_below_0.3K = colMeans(sim_mat < Blim_relK),
  E_mean = colMeans(sim_mat),
  E_median = apply(sim_mat, 2, median)
) %>% dplyr::mutate(
  Bcrit_abs = round(Bcrit_abs, 2),
  Blim_0.3K = round(Blim_relK, 2)
)

readr::write_csv(risk_tbl, file.path(outputs_dir, "risk_analysis.csv"))

cat("\nВероятностный анализ (", forecast_years[1], "):\n", sep = "")
cat("- P(Index < Bcrit_abs=", round(Bcrit_abs, 2), ") = ", round(risk_tbl$P_below_abs[1], 3), "\n", sep = "")
cat("- P(Index < 0.3*K=", round(Blim_relK, 2), ") = ", round(risk_tbl$P_below_0.3K[1], 3), "\n", sep = "")

# Сценарии эксплуатации: пост-эксплуатационный индекс = I * (1 - u)
exploit_rates <- c(0, 0.1, 0.2, 0.3)
scen_list <- list()
for (u in exploit_rates) {
  sim_post <- sim_mat * (1 - u)
  scen_list[[as.character(u)]] <- tibble::tibble(
    Год = forecast_years,
    u = u,
    P_below_abs = colMeans(sim_post < Bcrit_abs),
    P_below_0.3K = colMeans(sim_post < Blim_relK),
    mean_post = colMeans(sim_post),
    median_post = apply(sim_post, 2, median)
  )
}
scenario_tbl <- dplyr::bind_rows(scen_list) %>%
  dplyr::mutate(mean_post_round = round(mean_post, 2),
                median_post_round = round(median_post, 2))

readr::write_csv(scenario_tbl, file.path(outputs_dir, "scenario_risk_analysis.csv"))

# Визуализация сценариев для первого года прогноза
scen_year1 <- forecast_years[1]
scen_y1 <- scenario_tbl %>% dplyr::filter(Год == scen_year1)

p_scen <- ggplot(scen_y1, aes(x = factor(u), y = mean_post)) +
  geom_col(fill = "#3498db") +
  geom_errorbar(aes(ymin = mean_post - 1.96*sd_approx[1], ymax = mean_post + 1.96*sd_approx[1]), width = 0.2) +
  geom_hline(yintercept = Bcrit_abs, linetype = "dashed", color = "#b03a2e") +
  geom_hline(yintercept = Blim_relK, linetype = "dotted", color = "#7d3c98") +
  labs(title = paste0("Сценарии эксплуатации: ожидаемый индекс (", scen_year1, ")"),
       x = "Доля изъятия u", y = "Индекс (усл. ед.)") +
  theme_minimal()

ggsave(filename = file.path(outputs_dir, paste0("scenario_bar_", scen_year1, ".png")), plot = p_scen, width = 7, height = 4.5, dpi = 150)

# ===============
# 8. ОТЧЁТЫ/ЭКСПОРТЫ
# ===============
final_report <- list(
  исходные_данные = index_data,
  сравнение_моделей_на_тесте = results,
  cv_smape_h1 = cv_table,
  лучшая_модель = best_model_name,
  прогноз = forecast_details,
  риск_анализ = risk_tbl,
  сценарии = scenario_tbl,
  параметры = list(K = K, growth_cap_enable = growth_cap_enable, g_bounds = g_bounds),
  рекомендации = list(статус = status_text, меры = recs)
)

jsonlite::write_json(final_report, path = file.path(outputs_dir, "final_report.json"), pretty = TRUE, auto_unbox = TRUE)
saveRDS(final_report, file = file.path(outputs_dir, "final_report.rds"))

summary_lines <- c(
  "=== АНАЛИЗ ЗАВЕРШЕН ===",
  paste0("Лучшая модель: ", best_model_name),
  paste0("Прогноз на ", forecast_years[1], ": ", forecast_details$Прогноз[1]),
  paste0("80% ДИ: ", forecast_details$Нижняя_80[1], " - ", forecast_details$Верхняя_80[1]),
  paste0("P(Index < Bcrit_abs=", round(Bcrit_abs, 2), ") в ", forecast_years[1], ": ", round(risk_tbl$P_below_abs[1], 3)),
  paste0("P(Index < 0.3*K=", round(Blim_relK, 2), ") в ", forecast_years[1], ": ", round(risk_tbl$P_below_0.3K[1], 3))
)

writeLines(summary_lines, con = file.path(outputs_dir, "summary.txt"))
utils::capture.output(utils::sessionInfo(), file = file.path(outputs_dir, "sessionInfo.txt"))

# Консольный итог
cat("\n=== АНАЛИЗ ЗАВЕРШЕН ===\n")
cat("Лучшая модель:", best_model_name, "\n")
cat("Прогноз на ", forecast_years[1], ": ", forecast_details$Прогноз[1], "\n", sep = "")
cat("Диапазон неопределенности (80% ДИ): ", forecast_details$Нижняя_80[1], " - ", forecast_details$Верхняя_80[1], "\n", sep = "")
