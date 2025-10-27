# ==============================================================================
# ПРАКТИЧЕСКОЕ ЗАНЯТИЕ (урезанная версия): LM / GLM / GAM
# Цель: подготовка данных, обучение только трёх моделей (LM, GLM[Gamma], GAM),
#       диагностика, выбор лучшей по time-slice CV и прогноз 2022–2024 + график.
# ==============================================================================


# ==============================================================================
# 1) ПОДГОТОВКА ДАННЫХ
# Создаём NAOspring, фиксируем финальный набор признаков, сохраняем CSV.
# ==============================================================================

if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  readxl, tidyverse, caret, corrplot, mgcv,
  lmtest, car, GGally, lattice
)

rm(list = ls())
set.seed(123)
setwd("C:/RECRUITMENT/")  # при необходимости замените на ваш путь

# 1.2 Загрузка исходных данных и приведение типов
DATA <- readxl::read_excel("RECRUITMENT.xlsx", sheet = "RECRUITMENT") %>%
  filter(YEAR > 1989 & YEAR < 2022) %>%
  mutate(
    across(starts_with("T"), as.numeric),
    across(starts_with("I"), as.numeric),
    across(starts_with("O"), as.numeric),
    across(where(is.character), ~na_if(., "NA"))
  )

# 1.3 Создаём NAOspring (если есть NAO3, NAO4, NAO5)
if (all(c("NAO3","NAO4","NAO5") %in% names(DATA))) {
  DATA <- DATA %>%
    mutate(NAOspring = rowMeans(pick(NAO3, NAO4, NAO5), na.rm = TRUE)) %>%
    select(-NAO3, -NAO4, -NAO5)
}

# 1.4 Финальный учебный набор предикторов (фиксируем)
needed <- c("codTSB", "T12", "I5", "NAOspring", "haddock68")
stopifnot(all(needed %in% names(DATA)))

model_data <- DATA %>%
  select(YEAR, all_of(needed), R3haddock) %>%
  drop_na()

write.csv(model_data, "selected_predictors_dataset.csv", row.names = FALSE)
glimpse(model_data)


# ==============================================================================
# 2) ОБУЧЕНИЕ ТОЛЬКО LM / GLM(Gamma) / GAM С БАЗОВЫМИ ДИАГНОСТИКАМИ
# ==============================================================================

model_data <- read.csv("selected_predictors_dataset.csv", header = TRUE, stringsAsFactors = FALSE)

# Разделим на train/test (holdout) для первичного контроля качества
df <- model_data %>% select(-YEAR)
set.seed(123)
train_idx <- caret::createDataPartition(df$R3haddock, p = 0.8, list = FALSE)
train <- df[train_idx, ]
test  <- df[-train_idx, ]

ctrl <- caret::trainControl(method = "cv", number = 5, savePredictions = "final")

# Кастомный GAM (mgcv) для caret
gam_spec <- list(
  type = "Regression", library = "mgcv", loop = NULL,
  parameters = data.frame(parameter = "none", class = "character", label = "none"),
  grid = function(x,y,len=NULL,search="grid") data.frame(none = NA),
  fit = function(x,y,...) {
    df <- x; df$R3haddock <- y
    mgcv::gam(
      R3haddock ~ s(codTSB,bs="tp") + s(T12,bs="tp") + s(I5,bs="tp") +
                  s(NAOspring,bs="tp") + s(haddock68,bs="tp"),
      data = df, method = "REML", select = TRUE, ...
    )
  },
  predict = function(modelFit, newdata, submodels = NULL) {
    predict(modelFit, newdata = newdata, type = "response")
  },
  prob = NULL, sort = function(x) x
)

# Обучаем три модели на train (5-fold CV)
lm_model  <- caret::train(R3haddock ~ ., data = train, method = "lm", trControl = ctrl)
glm_model <- caret::train(R3haddock ~ ., data = train, method = "glm",
                          family = Gamma(link = "log"), trControl = ctrl)
gam_model <- caret::train(x = train[, -which(names(train)=="R3haddock")],
                          y = train$R3haddock, method = gam_spec, trControl = ctrl)

# Метрики
rmse  <- function(a, p) sqrt(mean((a - p)^2, na.rm = TRUE))
mae   <- function(a, p) mean(abs(a - p), na.rm = TRUE)
r2    <- function(a, p) 1 - sum((a - p)^2, na.rm = TRUE) / sum((a - mean(a))^2, na.rm = TRUE)

# Оценка на holdout
y_test <- test$R3haddock
preds_test <- list(
  LM  = predict(lm_model,  test),
  GLM = predict(glm_model, test),
  GAM = predict(gam_model, test)
)
metrics_table <- do.call(rbind, lapply(names(preds_test), function(nm){
  data.frame(Model = nm,
             RMSE = rmse(y_test, preds_test[[nm]]),
             MAE  = mae (y_test, preds_test[[nm]]),
             R2   = r2  (y_test, preds_test[[nm]]),
             row.names = NULL)
})) %>% arrange(RMSE, MAE)
print(round(metrics_table, 3))


# ==============================================================================
# 3) ДИАГНОСТИКА МОДЕЛЕЙ (на полном обучении до 2021) 
# Примечание: диагностика проводится по finalModel внутри caret-объектов.
# ==============================================================================

# Подготовим полный датасет для обучения финальных моделей (<= 2021)
full_fit_df <- model_data %>% arrange(YEAR) %>% select(-YEAR)

# Переобучим модели на полном наборе (method = "none") для чистых остатков
lm_full  <- caret::train(R3haddock ~ ., data = full_fit_df, method = "lm",
                         trControl = caret::trainControl(method = "none"))
glm_full <- caret::train(R3haddock ~ ., data = full_fit_df, method = "glm",
                         family = Gamma(link = "log"),
                         trControl = caret::trainControl(method = "none"))
gam_full <- caret::train(x = full_fit_df[, -which(names(full_fit_df)=="R3haddock")],
                         y = full_fit_df$R3haddock, method = gam_spec,
                         trControl = caret::trainControl(method = "none"))

# --- Диагностика LM -----------------------------------------------------------
cat("\n[LM] Сводка модели:\n"); print(summary(lm_full$finalModel))
cat("\n[LM] VIF (мультиколлинеарность):\n"); print(car::vif(lm_full$finalModel))
cat("\n[LM] Тест Бройша–Пагана (гомоскедастичность):\n"); print(lmtest::bptest(lm_full$finalModel))
cat("\n[LM] Тест Дарбина–Уотсона (автокорреляция остатков):\n"); print(lmtest::dwtest(lm_full$finalModel))

op <- par(no.readonly = TRUE); par(mfrow = c(2,2))
plot(lm_full$finalModel)  # стандартные графики диагностики lm
par(op)

# --- Диагностика GLM (Gamma) -------------------------------------------------
cat("\n[GLM-Gamma] Сводка модели:\n"); print(summary(glm_full$finalModel))
glm_resid <- residuals(glm_full$finalModel, type = "pearson")
dispersion <- sum(glm_resid^2, na.rm = TRUE) / glm_full$finalModel$df.residual
cat(sprintf("[GLM-Gamma] Оценка дисперсионного параметра (Pearson): %.3f\n", dispersion))

par(mfrow = c(1,2))
plot(fitted(glm_full$finalModel), glm_resid,
     xlab = "Fitted values", ylab = "Pearson residuals",
     main = "GLM-Gamma: Residuals vs Fitted"); abline(h = 0, col = "red")
qqnorm(glm_resid, main = "GLM-Gamma: QQ-plot (Pearson residuals)"); qqline(glm_resid, col = "red")
par(op)

# --- Диагностика GAM ---------------------------------------------------------
cat("\n[GAM] Сводка модели:\n"); print(summary(gam_full$finalModel))
cat("\n[GAM] Concurvity (нелинейная коллинеарность):\n"); print(mgcv::concurvity(gam_full$finalModel, full = TRUE))
mgcv::gam.check(gam_full$finalModel)
par(mfrow = c(2,3))
plot(gam_full$finalModel, shade = TRUE, pages = 1, rug = TRUE, scale = 0)
par(op)


# ==============================================================================
# 4) ВЫБОР ЛУЧШЕЙ МОДЕЛИ ПО TIME-SLICE CV (ГОРИЗОНТ 3 ГОДА)
# ==============================================================================

md <- model_data %>% arrange(YEAR)
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

slices <- caret::createTimeSlices(1:n_train, initialWindow = initialWindow, horizon = horizon, fixedWindow = FALSE)
ctrl_ts <- caret::trainControl(method = "cv", index = slices$train, indexOut = slices$test, savePredictions = "final")

fit_ts <- function(method, form, data, ctrl, ...) {
  out <- try(caret::train(form, data = data, method = method, trControl = ctrl, ...), TRUE)
  if (inherits(out, "try-error")) NULL else out
}

lm_ts  <- fit_ts("lm",        R3haddock ~ ., train_ts, ctrl_ts)
glm_ts <- fit_ts("glm",       R3haddock ~ ., train_ts, ctrl_ts, family = Gamma(link = "log"))
gam_ts <- caret::train(x = train_ts[, -which(names(train_ts)=="R3haddock")],
                       y = train_ts$R3haddock, method = gam_spec, trControl = ctrl_ts)

rmse  <- function(a, p) sqrt(mean((a - p)^2, na.rm = TRUE))
mae   <- function(a, p) mean(abs(a - p), na.rm = TRUE)

cv_metrics <- function(m) {
  if (is.null(m$pred) || !"Resample" %in% names(m$pred)) return(c(RMSE=NA, MAE=NA))
  by_slice <- m$pred %>% dplyr::group_by(Resample) %>%
    dplyr::summarise(RMSE=rmse(obs,pred), MAE=mae(obs,pred), .groups="drop")
  c(RMSE = mean(by_slice$RMSE, na.rm = TRUE), MAE = mean(by_slice$MAE, na.rm = TRUE))
}

models_ts <- list(LM = lm_ts, GLM = glm_ts, GAM = gam_ts)
models_ts <- models_ts[!vapply(models_ts, is.null, logical(1))]

cv_rank <- do.call(rbind, lapply(models_ts, cv_metrics)) %>% as.data.frame()
cv_rank$Model <- rownames(cv_rank)
cv_rank <- cv_rank[is.finite(cv_rank$RMSE), ] %>% dplyr::relocate(Model) %>% dplyr::arrange(RMSE, MAE)
print(cv_rank)

# Выбираем лучшую модель по RMSE
best_model_name <- cv_rank$Model[1]
cat(sprintf("\nЛучшая модель по time-slice CV: %s\n", best_model_name))


# ==============================================================================
# 5) ПРОГНОЗ 2022–2024 И ВИЗУАЛИЗАЦИЯ
# ==============================================================================

# Обучаем выбранную модель на всех данных до 2021 (full_fit_df)
train_none <- caret::trainControl(method = "none")

if (best_model_name == "LM") {
  best_full <- caret::train(R3haddock ~ ., data = full_fit_df, method = "lm", trControl = train_none)
} else if (best_model_name == "GLM") {
  best_full <- caret::train(R3haddock ~ ., data = full_fit_df, method = "glm", family = Gamma(link = "log"), trControl = train_none)
} else {
  best_full <- caret::train(x = full_fit_df[, -which(names(full_fit_df)=="R3haddock")],
                            y = full_fit_df$R3haddock, method = gam_spec, trControl = train_none)
}

# Остатки для PI: из time-slice CV выбранной модели, если есть; иначе — по фитам full
get_residuals_for_pi <- function(best_name) {
  m <- models_ts[[best_name]]
  if (!is.null(m) && !is.null(m$pred)) {
    return(m$pred$obs - m$pred$pred)
  }
  full_pred <- predict(best_full, newdata = full_fit_df)
  return(full_fit_df$R3haddock - full_pred)
}

resids <- get_residuals_for_pi(best_model_name)
q025 <- as.numeric(quantile(resids, 0.025, na.rm = TRUE))
q250 <- as.numeric(quantile(resids, 0.250, na.rm = TRUE))
q750 <- as.numeric(quantile(resids, 0.750, na.rm = TRUE))
q975 <- as.numeric(quantile(resids, 0.975, na.rm = TRUE))

# Сценарии будущего (по умолчанию — средние до 2021)
fc_start <- 2022
pred_cols <- c("codTSB","T12","I5","NAOspring","haddock68")
mu <- model_data %>% dplyr::filter(YEAR > 1989 & YEAR < fc_start) %>%
  dplyr::summarise(across(all_of(pred_cols), ~mean(.x, na.rm = TRUE))) %>% as.list()

build_future <- function(years, mu) {
  df <- tibble::tibble(YEAR = years)
  for (v in pred_cols) df[[v]] <- mu[[v]]
  df
}

future_years <- fc_start:2024
scenario_future <- build_future(future_years, mu)

pred_future <- predict(best_full, newdata = scenario_future)
forecast_tbl <- tibble::tibble(
  YEAR      = scenario_future$YEAR,
  Model     = best_model_name,
  pred_mean = as.numeric(pred_future),
  PI50_low  = pred_future + q250, PI50_high = pred_future + q750,
  PI95_low  = pred_future + q025, PI95_high = pred_future + q975
)
print(round(forecast_tbl, 0))

# Непрерывный ряд 1990–2024 для графика (ленты ДИ + линия прогноза)
pred_df <- dplyr::bind_rows(
  model_data %>% dplyr::select(YEAR, all_of(pred_cols)),
  scenario_future
) %>% dplyr::distinct(YEAR, .keep_all = TRUE) %>% dplyr::arrange(YEAR)

pred_df$Pred      <- as.numeric(predict(best_full, newdata = pred_df))
pred_df$PI50_low  <- pred_df$Pred + q250
pred_df$PI50_high <- pred_df$Pred + q750
pred_df$PI95_low  <- pred_df$Pred + q025
pred_df$PI95_high <- pred_df$Pred + q975

hist_df <- model_data %>% dplyr::select(YEAR, R3haddock)

library(ggplot2)
ggplot() +
  geom_ribbon(data = pred_df, aes(x = YEAR, ymin = PI95_low, ymax = PI95_high),
              fill = "grey80", alpha = 0.25) +
  geom_ribbon(data = pred_df, aes(x = YEAR, ymin = PI50_low, ymax = PI50_high),
              fill = "grey60", alpha = 0.35) +
  geom_line(data = subset(pred_df, YEAR < fc_start), aes(x = YEAR, y = Pred),
            color = "steelblue4", linewidth = 1) +
  geom_line(data = subset(pred_df, YEAR >= fc_start-1), aes(x = YEAR, y = Pred),
            color = "steelblue4", linewidth = 1, linetype = "dashed") +
  geom_point(data = hist_df, aes(x = YEAR, y = R3haddock),
             color = "black", size = 2, alpha = 0.9) +
  scale_x_continuous(expand = expansion(mult = c(0, 0))) +
  labs(
    title = paste0("Пополнение R3haddock: факт (1990–2021) и прогноз (2022–2024) — ", best_model_name),
    subtitle = "Прогноз — пунктир, интервалы — эмпирические из остатков",
    x = "Год", y = "R3haddock"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")

# Конец

