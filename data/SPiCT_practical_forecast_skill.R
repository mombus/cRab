# ===============================================================
#  Практикум: Прогностическая сила продукционной модели (SPiCT)
#  Данные: 2005–2024; Оценка прогноза: 2022–2024
#  Авторский шаблон адаптирован под учебную задачу
# ===============================================================

# ----- Минимальная подготовка среды -----
options(scipen = 999, digits = 4)

need <- c("spict", "forecast", "ggplot2", "dplyr", "tidyr", "readr", "scales")
for (p in need) {
  if (!requireNamespace(p, quietly = TRUE)) {
    message("Устанавливаю пакет: ", p)
    install.packages(p, repos = "https://cloud.r-project.org")
  }
}

suppressPackageStartupMessages({
  library(spict)
  library(forecast)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(scales)
})

# ----- Папка для результатов -----
out_dir <- file.path(getwd(), "output")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

cat("\n=== Практикум SPiCT: прогнозная сила (2005–2024, h=3: 2022–2024) ===\n")

# ----- 1) Данные (пример из задания) -----
Year <- 2005:2024
Catch <- c(5, 7, 6, 10, 14, 25, 28, 30, 32, 35, 25, 20, 15, 12, 10, 12, 10, 13, 11, 12)
CPUEIndex <- c(27.427120, 26.775958, 16.811997, 22.979653, 29.048568, 
               29.996072, 16.476301, 17.174455, 10.537272, 14.590435,
               8.286352, 11.394168, 15.537878, 13.791166, 11.527548, 
               15.336093, 12.154069, 15.568450, 16.221933, 13.421132)
BESSIndex <- c(NA, 16.270375, 20.691355, 15.141784, 18.594620, 
               15.975548, 13.792012, 13.328805, 11.659744, 11.753855,
               9.309859, 7.104886, 7.963839, 9.161322, 10.271221, 
               9.822960, 10.347376, 11.703610, 13.679876, 13.413696)
stopifnot(length(Year)==length(Catch), length(Year)==length(CPUEIndex), length(Year)==length(BESSIndex))

# ----- 2) SPiCT: фит полных данных 2005–2024 и извлечение годовой биомассы (истина) -----
input_data <- list(
  timeC = Year,
  obsC  = Catch,
  timeI = list(Year + 0.5, Year + 0.75),
  obsI  = list(CPUEIndex, BESSIndex)
)

inp <- check.inp(input_data, verbose = FALSE)
# Минимальные настройки и умеренные priors
inp$priors$logn      <- c(log(2), 0.1, 1);  inp$ini$logn <- log(2); inp$phases$logn <- -1
inp$priors$logK      <- c(log(150), 0.7, 1)
# Чуть точнее интегрирование
inp$dteuler <- 1/16

fit_all <- try(fit.spict(inp), silent = TRUE)
if (inherits(fit_all, "try-error")) stop("SPiCT не сошелся на полном наборе данных")

# Функция: перевести сплошной ряд времени SPiCT в годовую биомассу (последняя точка года)
extract_annual_B <- function(fit) {
  gp <- get.par("logB", fit, exp = TRUE)
  if (is.null(dim(gp))) {
    stop("Не удалось извлечь временной ряд биомассы из SPiCT")
  }
  time <- fit$inp$time
  df <- data.frame(time = as.numeric(time), B = as.numeric(gp[,"est"]))
  df$year <- floor(df$time + 1e-9)
  annual <- df |>
    group_by(year) |>
    arrange(time, .by_group = TRUE) |>
    summarise(B = dplyr::last(B), .groups = "drop")
  annual <- annual |>
    filter(year %in% Year) |>
    arrange(year)
  return(annual)
}

annual_B_all <- extract_annual_B(fit_all)
truth <- left_join(data.frame(year = Year), annual_B_all, by = c("year" = "year"))
if (any(is.na(truth$B))) {
  warning("В годовом ряду биомассы есть NA; проверьте входные данные/сходимость.")
}

# ----- 3) Модели прогноза по ряду биомассы (обучение до 2021, прогноз 2022–2024) -----
H <- 3
train_years <- Year[Year <= 2021]
test_years  <- 2022:2024

y_train <- ts(truth$B[truth$year %in% train_years], start = min(train_years), frequency = 1)
true_test <- truth$B[truth$year %in% test_years]

fit_ets   <- ets(y_train)
fit_arima <- auto.arima(y_train)

fc_list <- list(
  Naive   = naive(y_train, h = H),
  Drift   = rwf(y_train, h = H, drift = TRUE),
  Mean    = meanf(y_train, h = H),
  Theta   = thetaf(y_train, h = H),
  ETS     = forecast(fit_ets, h = H),
  ARIMA   = forecast(fit_arima, h = H)
)

# ----- 4) SPiCT-динамика как отдельный сценарий прогноза -----
# Обучаем SPiCT на усеченных данных (до 2021), затем делаем детерминированный прогноз
# по дискретной логистике с вычетом фактического вылова 2022–2024 (hindcast under realized C)

trim_by_year <- function(v_time, v_obs, last_year) {
  keep <- v_time <= last_year
  list(time = v_time[keep], obs = v_obs[keep])
}

inp_train <- list(timeC = input_data$timeC, obsC = input_data$obsC, timeI = input_data$timeI, obsI = input_data$obsI)
# Ограничиваем уловы и индексы годом 2021
trimC <- trim_by_year(inp_train$timeC, inp_train$obsC, 2021)
inp_train$timeC <- trimC$time; inp_train$obsC <- trimC$obs
for (k in seq_along(inp_train$timeI)) {
  tr <- trim_by_year(inp_train$timeI[[k]], inp_train$obsI[[k]], 2021.999)
  inp_train$timeI[[k]] <- tr$time
  inp_train$obsI[[k]]  <- tr$obs
}

inp_train <- check.inp(inp_train, verbose = FALSE)
inp_train$priors$logn <- c(log(2), 0.1, 1); inp_train$ini$logn <- log(2); inp_train$phases$logn <- -1
inp_train$priors$logK <- c(log(150), 0.7, 1)
inp_train$dteuler <- 1/16

fit_train <- try(fit.spict(inp_train), silent = TRUE)
if (inherits(fit_train, "try-error")) {
  warning("SPiCT (train до 2021) не сошелся; сценарий SPiCT-динамика будет пропущен.")
  spict_pred <- data.frame(year = test_years, point = NA_real_)
} else {
  # Параметры динамики
  r_hat <- as.numeric(exp(fit_train$par.fixed["logr"]))
  K_hat <- as.numeric(exp(fit_train$par.fixed["logK"]))
  annual_B_train <- extract_annual_B(fit_train)
  B0 <- as.numeric(annual_B_train$B[annual_B_train$year == 2021])
  if (!is.finite(B0)) {
    warning("Не удалось определить B(2021) из SPiCT train; SPiCT-сценарий пропущен.")
    spict_pred <- data.frame(year = test_years, point = NA_real_)
  } else {
    C_future <- Catch[Year %in% test_years]
    logistic_forward <- function(B0, r, K, Cvec) {
      out <- numeric(length(Cvec)); Bc <- B0
      for (i in seq_along(Cvec)) {
        Bn <- Bc + r * Bc * (1 - Bc / K) - Cvec[i]
        Bn <- max(Bn, 0)
        out[i] <- Bn
        Bc <- Bn
      }
      out
    }
    spict_point <- logistic_forward(B0, r_hat, K_hat, C_future)
    spict_pred <- data.frame(year = test_years, point = spict_point)
  }
}

# ----- 5) Сбор прогнозов в единую таблицу -----
fc_to_df <- function(obj, name) {
  levs <- colnames(obj$lower)
  lo80 <- if ("80%" %in% levs) obj$lower[,"80%"] else NA_real_
  hi80 <- if ("80%" %in% levs) obj$upper[,"80%"] else NA_real_
  lo95 <- if ("95%" %in% levs) obj$lower[,"95%"] else NA_real_
  hi95 <- if ("95%" %in% levs) obj$upper[,"95%"] else NA_real_
  data.frame(
    year  = test_years,
    method= name,
    point = as.numeric(obj$mean),
    lo80  = as.numeric(lo80),
    hi80  = as.numeric(hi80),
    lo95  = as.numeric(lo95),
    hi95  = as.numeric(hi95),
    stringsAsFactors = FALSE
  )
}

fc_tbl <- bind_rows(lapply(names(fc_list), function(nm) fc_to_df(fc_list[[nm]], nm)))
fc_tbl <- bind_rows(
  fc_tbl,
  spict_pred |>
    mutate(method = "SPiCT-dynamic",
           lo80 = NA_real_, hi80 = NA_real_, lo95 = NA_real_, hi95 = NA_real_)
)

# Сохраним прогнозы
write_csv(fc_tbl, file.path(out_dir, "forecast_values.csv"))

# ----- 6) Оценка прогностической силы (относительно наивного) -----
true_df <- data.frame(year = test_years, truth = true_test)

metrics <- fc_tbl |>
  select(method, year, point) |>
  left_join(true_df, by = "year") |>
  group_by(method) |>
  summarise(
    H      = n(),
    MSE    = mean((point - truth)^2, na.rm = TRUE),
    RMSE   = sqrt(MSE),
    MAE    = mean(abs(point - truth), na.rm = TRUE),
    .groups = "drop"
  )

mse_naive <- metrics$MSE[metrics$method == "Naive"]
metrics <- metrics |>
  mutate(Skill_MSE_vs_Naive = 1 - MSE / mse_naive)

write_csv(metrics, file.path(out_dir, "skill_metrics.csv"))

cat("\nИтоги точности (2022–2024):\n")
print(metrics)

# ----- 7) Визуализация -----
# 7.1 Ряд биомассы и прогнозы
plot_df <- truth |>
  select(year, B) |>
  rename(truth = B) |>
  mutate(method = "Истина")

last_train <- data.frame(year = 2021, method = names(fc_list), point = rep(as.numeric(tail(y_train,1)), length(fc_list)))

plot_fc <- fc_tbl |>
  select(year, method, point)

p1 <- ggplot() +
  geom_line(data = truth, aes(x = year, y = B), color = "black", linewidth = 1) +
  geom_point(data = truth, aes(x = year, y = B), color = "black", size = 1.5) +
  geom_line(data = bind_rows(last_train, plot_fc), aes(x = year, y = point, color = method), linewidth = 0.9, linetype = "dashed") +
  geom_point(data = plot_fc, aes(x = year, y = point, color = method), size = 2) +
  scale_color_brewer(palette = "Set1") +
  labs(title = "Биомасса (SPiCT) и прогнозы на 2022–2024",
       x = "Год", y = "Биомасса (тыс. т)", color = "Сценарий") +
  theme_minimal(base_size = 12)

ggsave(filename = file.path(out_dir, "biomass_forecasts.png"), plot = p1, width = 10, height = 6, dpi = 300)

# 7.2 Прогностическая сила (относительно Naive, по MSE)
p2 <- metrics |>
  mutate(method = factor(method, levels = metrics$method[order(metrics$Skill_MSE_vs_Naive, decreasing = TRUE)])) |>
  ggplot(aes(x = method, y = Skill_MSE_vs_Naive, fill = method)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = 0, linewidth = 0.6, linetype = "dashed") +
  scale_y_continuous(labels = label_percent(accuracy = 1)) +
  labs(title = "Прогностическая сила vs Naive (MSE)", x = NULL, y = "1 - MSE_model / MSE_naive") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")

ggsave(filename = file.path(out_dir, "predictive_skill.png"), plot = p2, width = 9, height = 5, dpi = 300)

cat("\nФайлы сохранены в папке: ", out_dir, "\n",
    "- forecast_values.csv (точки и интервалы прогнозов)\n",
    "- skill_metrics.csv (метрики и прогностическая сила)\n",
    "- biomass_forecasts.png (график ряд + прогнозы)\n",
    "- predictive_skill.png (график прогностической силы)\n", sep = "")

cat("\n\nПримечание: сценарий 'SPiCT-dynamic' — детерминированный прогноз\n",
    "по дискретной логистике с параметрами r,K из SPiCT (train≤2021)\n",
    "и реализованными выловами 2022–2024 (hindcast).\n", sep = "")
