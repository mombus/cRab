# ==============================================================================
# Урезанная версия: только LM / GLM(Gamma) / GAM
# Без caret/train: стандартная оценка параметров lm/glm/gam, собственная time-slice CV,
# выбор лучшей модели, прогноз 2022–2024, эмпирические интервалы и график.
# ==============================================================================

# 0) Пакеты и окружение --------------------------------------------------------
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  readxl, tidyverse, mgcv, lmtest, car, ggplot2, corrplot
)

rm(list = ls())
set.seed(123)
setwd("C:/RECRUITMENT/")  # при необходимости измените путь


# 1) Подготовка данных ---------------------------------------------------------
# Загрузка, приведение типов, создание NAOspring, фиксируем набор признаков
DATA <- readxl::read_excel("RECRUITMENT.xlsx", sheet = "RECRUITMENT") %>%
  filter(YEAR > 1989 & YEAR < 2022) %>%
  mutate(
    across(starts_with("T"), as.numeric),
    across(starts_with("I"), as.numeric),
    across(starts_with("O"), as.numeric),
    across(where(is.character), ~na_if(., "NA"))
  )

if (all(c("NAO3","NAO4","NAO5") %in% names(DATA))) {
  DATA <- DATA %>%
    mutate(NAOspring = rowMeans(pick(NAO3, NAO4, NAO5), na.rm = TRUE)) %>%
    select(-NAO3, -NAO4, -NAO5)
}

needed <- c("codTSB", "T12", "I5", "NAOspring", "haddock68")
stopifnot(all(needed %in% names(DATA)))

model_data <- DATA %>%
  select(YEAR, all_of(needed), R3haddock) %>%
  drop_na() %>%
  arrange(YEAR)

write.csv(model_data, "selected_predictors_dataset.csv", row.names = FALSE)
glimpse(model_data)


# 2) Формулы моделей и вспомогательные функции --------------------------------
f_lm  <- as.formula("R3haddock ~ codTSB + T12 + I5 + NAOspring + haddock68")
f_gam <- as.formula("R3haddock ~ s(codTSB,bs='tp') + s(T12,bs='tp') + s(I5,bs='tp') + s(NAOspring,bs='tp') + s(haddock68,bs='tp')")

rmse <- function(a, p) sqrt(mean((a - p)^2, na.rm = TRUE))
mae  <- function(a, p) mean(abs(a - p), na.rm = TRUE)
r2   <- function(a, p) 1 - sum((a - p)^2, na.rm = TRUE) / sum((a - mean(a))^2, na.rm = TRUE)

safe_fit <- function(expr) {
  out <- try(eval(expr), silent = TRUE)
  if (inherits(out, "try-error")) NULL else out
}


# 3) Time-slice CV (expanding window, h=3) + хронологический тест -------------
md <- model_data
md_for_fit <- md %>% select(codTSB, T12, I5, NAOspring, haddock68, R3haddock)

n <- nrow(md_for_fit)
holdout_frac <- 0.2
n_test <- max(4, ceiling(n * holdout_frac))
train_ts <- head(md_for_fit, n - n_test)
test_ts  <- tail(md_for_fit, n_test)

n_train <- nrow(train_ts)
initial_frac <- 0.6
horizon      <- 3
initialWindow <- max(10, floor(initial_frac * n_train))
if (initialWindow + horizon > n_train) initialWindow <- n_train - horizon

# Аккумулируем метрики и остатки по срезам
cv_summ <- tibble(Model = character(), RMSE = double(), MAE = double())
resids_cv <- list(LM = c(), GLM = c(), GAM = c())

slice_id <- 0
for (i in seq(initialWindow, n_train - horizon)) {
  slice_id <- slice_id + 1
  idx_tr <- 1:i
  idx_te <- (i+1):(i+horizon)
  dtr <- train_ts[idx_tr, ]
  dte <- train_ts[idx_te, ]

  # LM
  lm_fit <- safe_fit(quote(lm(f_lm, data = dtr)))
  if (!is.null(lm_fit)) {
    pr <- try(predict(lm_fit, newdata = dte), silent = TRUE)
    if (!inherits(pr, "try-error")) {
      cv_summ <- add_row(cv_summ, Model = "LM", RMSE = rmse(dte$R3haddock, pr), MAE = mae(dte$R3haddock, pr))
      resids_cv$LM <- c(resids_cv$LM, dte$R3haddock - pr)
    }
  }

  # GLM (Gamma)
  glm_fit <- safe_fit(quote(glm(f_lm, data = dtr, family = Gamma(link = "log"))))
  if (!is.null(glm_fit)) {
    pr <- try(predict(glm_fit, newdata = dte, type = "response"), silent = TRUE)
    if (!inherits(pr, "try-error")) {
      cv_summ <- add_row(cv_summ, Model = "GLM", RMSE = rmse(dte$R3haddock, pr), MAE = mae(dte$R3haddock, pr))
      resids_cv$GLM <- c(resids_cv$GLM, dte$R3haddock - pr)
    }
  }

  # GAM
  gam_fit <- safe_fit(quote(mgcv::gam(f_gam, data = dtr, method = "REML", select = TRUE)))
  if (!is.null(gam_fit)) {
    pr <- try(predict(gam_fit, newdata = dte, type = "response"), silent = TRUE)
    if (!inherits(pr, "try-error")) {
      cv_summ <- add_row(cv_summ, Model = "GAM", RMSE = rmse(dte$R3haddock, pr), MAE = mae(dte$R3haddock, pr))
      resids_cv$GAM <- c(resids_cv$GAM, dte$R3haddock - pr)
    }
  }
}

# Средние метрики по моделям
cv_rank <- cv_summ %>% group_by(Model) %>% summarise(RMSE = mean(RMSE, na.rm = TRUE), MAE = mean(MAE, na.rm = TRUE), .groups = "drop") %>% arrange(RMSE, MAE)
print(cv_rank)

best_model_name <- cv_rank$Model[1]
cat(sprintf("\nЛучшая модель по time-slice CV: %s\n", best_model_name))

# Хронологический тест: обучаем на всём train_ts, прогнозируем на test_ts
fit_on <- function(model_name, data) {
  if (model_name == "LM") return(lm(f_lm, data = data))
  if (model_name == "GLM") return(glm(f_lm, data = data, family = Gamma(link = "log")))
  mgcv::gam(f_gam, data = data, method = "REML", select = TRUE)
}

predict_on <- function(fit, newdata, model_name) {
  if (model_name == "GLM") return(predict(fit, newdata = newdata, type = "response"))
  predict(fit, newdata = newdata)
}

fit_train <- fit_on(best_model_name, train_ts)
pred_te   <- predict_on(fit_train, test_ts, best_model_name)
test_metrics <- tibble(
  Model = best_model_name,
  RMSE  = rmse(test_ts$R3haddock, pred_te),
  MAE   = mae (test_ts$R3haddock, pred_te),
  R2    = r2  (test_ts$R3haddock, pred_te)
)
print(test_metrics)


# 4) Диагностика моделей (подгонка на всех данных до 2021) -------------------
full_fit_df <- md_for_fit

lm_full  <- lm(f_lm,  data = full_fit_df)
glm_full <- glm(f_lm, data = full_fit_df, family = Gamma(link = "log"))
gam_full <- mgcv::gam(f_gam, data = full_fit_df, method = "REML", select = TRUE)

cat("\n[LM] Сводка:\n"); print(summary(lm_full))
cat("\n[LM] VIF:\n"); print(car::vif(lm_full))
cat("\n[LM] Breusch–Pagan:\n"); print(lmtest::bptest(lm_full))
cat("\n[LM] Durbin–Watson:\n"); print(lmtest::dwtest(lm_full))

glm_resid <- residuals(glm_full, type = "pearson")
cat("\n[GLM-Gamma] Сводка:\n"); print(summary(glm_full))
cat(sprintf("[GLM-Gamma] Pearson dispersion: %.3f\n", sum(glm_resid^2, na.rm = TRUE) / glm_full$df.residual))

cat("\n[GAM] Сводка:\n"); print(summary(gam_full))
cat("\n[GAM] Concurvity:\n"); print(mgcv::concurvity(gam_full, full = TRUE))
mgcv::gam.check(gam_full)


# 5) Прогноз 2022–2024 и эмпирические интервалы ------------------------------
best_full <- switch(best_model_name,
  LM  = lm_full,
  GLM = glm_full,
  GAM = gam_full
)

# Остатки для PI: из CV выбранной модели, иначе из полного фита
resids <- if (length(resids_cv[[best_model_name]]) > 5) resids_cv[[best_model_name]] else residuals(best_full)

q025 <- as.numeric(quantile(resids, 0.025, na.rm = TRUE))
q250 <- as.numeric(quantile(resids, 0.250, na.rm = TRUE))
q750 <- as.numeric(quantile(resids, 0.750, na.rm = TRUE))
q975 <- as.numeric(quantile(resids, 0.975, na.rm = TRUE))

fc_start <- 2022
pred_cols <- c("codTSB","T12","I5","NAOspring","haddock68")
mu <- md %>% filter(YEAR > 1989 & YEAR < fc_start) %>% summarise(across(all_of(pred_cols), ~mean(.x, na.rm = TRUE))) %>% as.list()

build_future <- function(years, mu) {
  df <- tibble::tibble(YEAR = years)
  for (v in pred_cols) df[[v]] <- mu[[v]]
  df
}

future_years <- fc_start:2024
scenario_future <- build_future(future_years, mu)

predict_best <- function(fit, newdata, model_name) {
  if (model_name == "GLM") return(predict(fit, newdata = newdata, type = "response"))
  predict(fit, newdata = newdata)
}

pred_future <- predict_best(best_full, scenario_future, best_model_name)

forecast_tbl <- tibble::tibble(
  YEAR      = scenario_future$YEAR,
  Model     = best_model_name,
  pred_mean = as.numeric(pred_future),
  PI50_low  = pred_future + q250, PI50_high = pred_future + q750,
  PI95_low  = pred_future + q025, PI95_high = pred_future + q975
)
print(round(forecast_tbl, 0))


# 6) Визуализация 1990–2024 ---------------------------------------------------
pred_df <- bind_rows(
  md %>% select(YEAR, all_of(pred_cols)),
  scenario_future
) %>% distinct(YEAR, .keep_all = TRUE) %>% arrange(YEAR)

pred_df$Pred      <- as.numeric(predict_best(best_full, pred_df, best_model_name))
pred_df$PI50_low  <- pred_df$Pred + q250
pred_df$PI50_high <- pred_df$Pred + q750
pred_df$PI95_low  <- pred_df$Pred + q025
pred_df$PI95_high <- pred_df$Pred + q975

hist_df <- md %>% select(YEAR, R3haddock)

ggplot() +
  geom_ribbon(data = pred_df, aes(x = YEAR, ymin = PI95_low, ymax = PI95_high), fill = "grey80", alpha = 0.25) +
  geom_ribbon(data = pred_df, aes(x = YEAR, ymin = PI50_low, ymax = PI50_high), fill = "grey60", alpha = 0.35) +
  geom_line(data = subset(pred_df, YEAR < fc_start), aes(x = YEAR, y = Pred), color = "steelblue4", linewidth = 1) +
  geom_line(data = subset(pred_df, YEAR >= fc_start-1), aes(x = YEAR, y = Pred), color = "steelblue4", linewidth = 1, linetype = "dashed") +
  geom_point(data = hist_df, aes(x = YEAR, y = R3haddock), color = "black", size = 2, alpha = 0.9) +
  scale_x_continuous(expand = expansion(mult = c(0, 0))) +
  labs(title = paste0("Пополнение R3haddock: факт (1990–2021) и прогноз (2022–2024) — ", best_model_name),
       subtitle = "Прогноз — пунктир, интервалы — эмпирические из остатков",
       x = "Год", y = "R3haddock") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")

# ============================================================================
# Конец
# ============================================================================

