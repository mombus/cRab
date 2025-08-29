# ========================================================================================================================
# ПРАКТИЧЕСКОЕ ЗАНЯТИЕ: Модели пространственного распределения видов (SDM)
# Курс: "Оценка водных биоресурсов в среде R (для начинающих)"
# Автор: Баканев С. В. Дата: 27.08.2025
# Структура:
# 1) Подготовка и визуализация данных
# 2) Отбор переменных и анализ важности признаков
# 3) Построение моделей и анализ результатов
# ========================================================================================================================

# Установка рабочей директории
setwd("C:/SDM")
rm(list = ls())

# Подключение необходимых библиотек
library(tidyverse)    # Обработка данных и визуализация
library(janitor)      # Очистка имен переменных
library(recipes)      # Предобработка данных
library(caret)        # Машинное обучение
library(car)          # VIF анализ
library(Boruta)       # Отбор признаков
library(glmnet)       # LASSO регрессия
library(randomForest) # Случайный лес
library(mgcv)         # GAM модели
library(terra)        # Пространственный анализ
library(scales)       # Форматирование графиков

# ========================================================================================================================
# 1. ЗАГРУЗКА И ПРЕДВАРИТЕЛЬНАЯ ОБРАБОТКА ДАННЫХ
# ========================================================================================================================

# Загрузка данных
DATA <- read.csv("SDM_all_pred_full_set.csv")

# Предварительная обработка: безопасные имена и удаление пропусков в целевой переменной
df0 <- DATA %>%
  janitor::clean_names() %>%
  filter(!is.na(occ)) %>%
  select(-x, -y)  # Удаление координат и ненужных переменных

# Удаление переменных с near-zero variance
preds0 <- df0 %>% select(-occ)
nzv_idx <- caret::nearZeroVar(preds0)
if (length(nzv_idx) > 0) preds0 <- preds0[, -nzv_idx, drop = FALSE]
df1 <- bind_cols(occ = df0$occ, preds0)

# Импутация пропусков и стандартизация данных
rec <- recipe(occ ~ ., data = df1) %>%
  step_impute_knn(all_predictors()) %>%
  step_normalize(all_predictors())

prep_rec <- prep(rec)
dat <- bake(prep_rec, new_data = NULL)

# ========================================================================================================================
# 2. ОТБОР ПЕРЕМЕННЫХ И АНАЛИЗ МУЛЬТИКОЛЛИНЕАРНОСТИ
# ========================================================================================================================

# Удаление высококоррелированных переменных (коэффициент > 0.8)
corr <- cor(dat %>% select(-occ), use = "pairwise.complete.obs")
highCorr <- caret::findCorrelation(corr, cutoff = 0.8, names = TRUE, exact = TRUE)
dat_cf <- dat %>% select(occ, any_of(setdiff(names(dat)[names(dat) != "occ"], highCorr)))

# Функция для фильтрации по VIF
vif_filter <- function(df, thresh = 5) {
  vars <- setdiff(names(df), "occ")
  repeat {
    fit <- lm(occ ~ ., data = df[, c("occ", vars)])
    v <- car::vif(fit)
    if (max(v) < thresh) break
    drop_var <- names(which.max(v))
    vars <- setdiff(vars, drop_var)
    if (length(vars) == 0) break
  }
  df[, c("occ", vars)]
}

dat_vif <- vif_filter(dat_cf, thresh = 5)

# Отбор признаков с помощью Boruta
set.seed(42)
bor <- Boruta(occ ~ ., data = dat_vif, maxRuns = 200, doTrace = 0)
bor_selected <- getSelectedAttributes(bor, withTentative = FALSE)

# Визуализация результатов Boruta
plot(bor, las = 2, cex.axis = 0.7)
print(bor)

# Отбор признаков с помощью LASSO
x <- as.matrix(dat_vif %>% select(-occ))
y <- dat_vif$occ
cv <- cv.glmnet(x, y, alpha = 1, family = "binomial", standardize = FALSE)
coef_1se <- as.matrix(coef(cv, s = "lambda.1se"))
lasso_selected <- rownames(coef_1se)[coef_1se[, 1] != 0]
lasso_selected <- setdiff(lasso_selected, "(Intercept)")

# Ранжирование переменных по важности в LASSO
coef_all <- as.matrix(coef(cv, s = "lambda.1se"))
coef_tbl <- tibble(var = rownames(coef_all), coef = as.numeric(coef_all[, 1])) %>%
  filter(var != "(Intercept)") %>%
  arrange(desc(abs(coef)))

# Формирование финального набора переменных (5-8 наиболее важных)
consensus <- intersect(bor_selected, lasso_selected)
fill_from_lasso <- setdiff(coef_tbl$var, consensus)
fill_from_boruta <- setdiff(bor_selected, c(consensus, fill_from_lasso))

final_vars <- unique(c(consensus, fill_from_lasso, fill_from_boruta))
final_vars <- head(final_vars, 8)  # Ограничение до 8 переменных

if (length(final_vars) < 5) {
  final_vars <- unique(c(final_vars, head(coef_tbl$var, 5)))[1:5]
}

final_df <- dat_vif %>% select(occ, all_of(final_vars))

# ========================================================================================================================
# 3. АНАЛИЗ ВАЖНОСТИ ПЕРЕМЕННЫХ
# ========================================================================================================================

# Важность переменных по LASSO
lasso_imp <- coef_tbl %>%
  mutate(abs_coef = abs(coef)) %>%
  filter(var != "(Intercept)", abs_coef > 0) %>%
  arrange(desc(abs_coef)) %>%
  slice_head(n = 20)

ggplot(lasso_imp, aes(x = reorder(var, abs_coef), y = abs_coef)) +
  geom_col(fill = "#3B82F6") +
  coord_flip() +
  labs(x = "Переменная", y = "|коэффициент| (lambda.1se)",
       title = "LASSO важность (топ-20 по |коэффициенту|)") +
  theme_minimal(base_size = 12)

# Важность переменных по Random Forest
rf_data <- dat_vif %>% dplyr::select(occ, dplyr::all_of(final_vars))
is_classif <- is.factor(rf_data$occ) || all(rf_data$occ %in% c(0, 1), na.rm = TRUE)
if (is_classif && !is.factor(rf_data$occ)) {
  rf_data$occ <- factor(rf_data$occ)
}

rf <- randomForest(occ ~ ., data = rf_data, importance = TRUE, na.action = na.omit)
imp_mat <- randomForest::importance(rf)
imp_df <- as.data.frame(imp_mat) %>%
  tibble::rownames_to_column("var") %>%
  pivot_longer(cols = -var, names_to = "metric", values_to = "importance")

ggplot(imp_df, aes(x = reorder(var, importance), y = importance, fill = metric)) +
  geom_col(show.legend = FALSE) +
  coord_flip() +
  facet_wrap(~ metric, scales = "free_y") +
  scale_fill_manual(values = c("#10B981", "#F59E0B", "#6366F1", "#EF4444")) +
  labs(x = "Переменная", y = "Важность", title = "Важность признаков по Random Forest") +
  theme_minimal(base_size = 12)

# ========================================================================================================================
# 4. ПОСТРОЕНИЕ И АНАЛИЗ GAM МОДЕЛЕЙ
# ========================================================================================================================

df <- final_df
if (is.factor(df$occ)) df$occ <- as.numeric(df$occ) - 1L
if (!all(df$occ %in% c(0, 1))) df$occ <- ifelse(df$occ > 0, 1L, 0L)

pred_vars <- setdiff(names(df), "occ")
k_basis <- 8

# Функция для построения унимодальных GAM моделей
fit_gam_uni <- function(var_name) {
  form <- as.formula(paste0("occ ~ s(", var_name, ", k=", k_basis, ")"))
  fams <- list(
    binomial(link = "identity"),
    binomial(link = "probit"),
    binomial(link = "logit")
  )
  model <- NULL
  fam_used <- NULL
  for (f in fams) {
    model_try <- try(
      gam(form, data = df, family = f, method = "REML", select = TRUE),
      silent = TRUE
    )
    if (!inherits(model_try, "try-error")) {
      model <- model_try
      fam_used <- model$family$link
      break
    }
  }
  if (is.null(model)) stop(paste("Не удалось обучить GAM для", var_name))

  x <- df[[var_name]]
  x_seq <- seq(quantile(x, 0.02, na.rm = TRUE),
               quantile(x, 0.98, na.rm = TRUE),
               length.out = 200)

  newd <- tibble(!!var_name := x_seq)
  pr <- predict(model, newdata = newd, type = "link", se.fit = TRUE)

  inv <- model$family$linkinv
  tibble(
    variable = var_name,
    x = x_seq,
    prob  = pmax(pmin(inv(pr$fit), 1), 0),
    lower = pmax(pmin(inv(pr$fit - 1.96 * pr$se.fit), 1), 0),
    upper = pmax(pmin(inv(pr$fit + 1.96 * pr$se.fit), 1), 0),
    link  = fam_used
  ) %>%
    mutate(model_link = paste0("binomial(", link, ")"))
}

# Обучение GAM моделей и сбор кривых влияния
pd_all <- map_dfr(pred_vars, fit_gam_uni)

# Подготовка данных для визуализации
rug_data <- df %>%
  filter(occ == 1) %>%
  pivot_longer(cols = all_of(pred_vars), names_to = "variable", values_to = "x")

# Первый график: базовые кривые GAM с риджинами
ggplot(pd_all, aes(x = x, y = prob)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#93C5FD", alpha = 0.35) +
  geom_line(color = "#1D4ED8", linewidth = 1) +
  geom_rug(data = rug_data, aes(x = x), sides = "b", alpha = 0.25, inherit.aes = FALSE) +
  facet_wrap(~ variable, scales = "free_x") +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Значение предиктора (стандартизовано)",
       y = "Вероятность присутствия",
       title = "Унимодельные GAM (binomial): влияние каждого предиктора",
       subtitle = "Ссылка: identity (fallback > probit > logit); ленты — 95% ДИ") +
  theme_minimal(base_size = 12)

# Второй график: с джиттером и биновой эмпирической вероятностью
obs_raw <- df %>%
  pivot_longer(cols = all_of(pred_vars), names_to = "variable", values_to = "x") %>%
  select(variable, x, occ)

obs_bin <- obs_raw %>%
  group_by(variable) %>%
  mutate(bin = cut_number(x, 20)) %>%
  group_by(variable, bin) %>%
  summarise(
    x_mid = mean(x, na.rm = TRUE),
    occ_mean = mean(occ),
    n = n(), .groups = "drop"
  )

ggplot(pd_all, aes(x = x, y = prob)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#93C5FD", alpha = 0.35) +
  geom_line(color = "#1D4ED8", linewidth = 1) +
  geom_jitter(
    data = obs_raw, inherit.aes = FALSE,
    aes(x = x, y = occ),
    width = 0, height = 0.04, alpha = 0.25, size = 0.9, color = "black"
  ) +
  geom_point(
    data = obs_bin, inherit.aes = FALSE,
    aes(x = x_mid, y = occ_mean),
    color = "#111827", size = 1.6
  ) +
  facet_wrap(~ variable, scales = "free_x") +
  scale_y_continuous(limits = c(0, 1), labels = percent_format(accuracy = 1)) +
  labs(
    x = "Значение предиктора (стандартизовано)",
    y = "Вероятность присутствия",
    title = "Унимодельные GAM (binomial): влияние предикторов",
    subtitle = "Серые точки — фактические (джиттер); Черные точки — биновая средняя встречаемость"
  ) +
  theme_minimal(base_size = 12)

# ========================================================================================================================
# 5. ПОСТРОЕНИЕ ТАБЛИЦЫ НА ОСНОВЕ ГРАФИКОВ
# ========================================================================================================================

# Биннинг (20 квантильных бинов) и эмпирическая вероятность
obs_raw <- df %>%
  pivot_longer(cols = all_of(pred_vars), names_to = "variable", values_to = "x") %>%
  select(variable, x, occ)

obs_bin <- obs_raw %>%
  group_by(variable) %>%
  mutate(bin = cut_number(x, 20)) %>%
  group_by(variable, bin, .add = FALSE) %>%
  summarise(
    bin_id = cur_group_id(),
    x_min = min(x, na.rm = TRUE),
    x_max = max(x, na.rm = TRUE),
    x_mid = median(x, na.rm = TRUE),
    n = n(),
    occ_mean = mean(occ, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(variable, bin_id)

# Обучение GAM моделей для таблицы
fit_gam_uni_model <- function(var_name) {
  form <- as.formula(paste0("occ ~ s(", var_name, ", k=", k_basis, ")"))
  fams <- list(
    binomial(link = "identity"),
    binomial(link = "probit"),
    binomial(link = "logit")
  )
  for (f in fams) {
    m_try <- try(gam(form, data = df, family = f, method = "REML", select = TRUE), silent = TRUE)
    if (!inherits(m_try, "try-error")) return(m_try)
  }
  stop(paste("Не удалось обучить GAM для", var_name))
}

models <- set_names(map(pred_vars, fit_gam_uni_model), pred_vars)

# Предсказания GAM в центрах бинов
table_20bins <- obs_bin %>%
  group_by(variable) %>%
  group_modify(function(.x, .y) {
    var <- .y$variable[[1]]
    mdl <- models[[var]]
    newd <- tibble(!!rlang::sym(var) := .x$x_mid)
    pr <- predict(mdl, newdata = newd, type = "link", se.fit = TRUE)
    inv <- mdl$family$linkinv
    mutate(
      .x,
      gam_prob  = pmin(pmax(inv(pr$fit), 0), 1),
      gam_lower = pmin(pmax(inv(pr$fit - 1.96 * pr$se.fit), 0), 1),
      gam_upper = pmin(pmax(inv(pr$fit + 1.96 * pr$se.fit), 0), 1),
      link      = mdl$family$link
    )
  }) %>%
  ungroup() %>%
  select(variable, bin_id, x_min, x_max, x_mid, n, occ_mean, gam_prob, gam_lower, gam_upper, link)

# Вывод таблицы
print(table_20bins, n = 5)

# ========================================================================================================================
# 6. ФОРМИРОВАНИЕ ФИНАЛЬНОЙ ТАБЛИЦЫ ДАННЫХ
# ========================================================================================================================

# Создание финальной таблицы с исходными данными
final_table <- DATA %>%
  janitor::clean_names() %>%
  select(x, y, occ, all_of(final_vars))

# Сохранение результатов
write.csv(final_table, "final_sdm_table_with_na.csv", row.names = FALSE)

# Вывод структуры финальной таблицы
str(final_table)