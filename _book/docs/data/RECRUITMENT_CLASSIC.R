# ==============================================================================
# ПРАКТИЧЕСКОЕ ЗАНЯТИЕ: АНАЛИЗ ФАКТОРОВ, ВЛИЯЮЩИХ НА ПОПОЛНЕНИЕ 
# (КЛАССИЧЕСКИЕ МОДЕЛИ ЗАПАС-ПОПОЛНЕНИЕ)
# Курс: "Оценка водных биоресурсов в среде R (для начинающих)"
# ==============================================================================

# Установка и подключение ТОЛЬКО необходимых библиотек
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  tidyverse,    # Манипуляции с данными и визуализация
  FSA,          # Начальные оценки для моделей запас-пополнение
  minpack.lm,   # Улучшенный алгоритм нелинейной регрессии (nlsLM)
  car,          # Проверка допущений моделей
  mgcv,         # Построение GAM-моделей
  investr,      # Доверительные интервалы для нелинейных моделей
   caret)       # Расчет RMSE

# Очистка среды и установка рабочей директории
rm(list = ls())
setwd("C:/RECRUITMENT/")

# 1. ЗАГРУЗКА ДАННЫХ -----------------------------------------------------------
model_data <- read.csv("selected_predictors_dataset.csv", 
                      header = TRUE, 
                      stringsAsFactors = FALSE)

# Проверка структуры данных
str(model_data)

# Проверка на пропущенные значения (должно быть 0 после предобработки)
sum(is.na(model_data))

# 2. ПОДГОТОВКА ДАННЫХ ДЛЯ МОДЕЛЕЙ ЗАПАС-ПОПОЛНЕНИЕ ----------------------------
rec_data <- data.frame(
  S = model_data$haddock68,  # Нерестовый запас
  R = model_data$R3haddock   # Пополнение
)

# 3. ПОДГОНКА МОДЕЛИ РИКЕРА ----------------------------------------------------
ricker_starts <- FSA::srStarts(R ~ S, data = rec_data, type = "Ricker")
ricker_model <- minpack.lm::nlsLM(
  R ~ a * S * exp(-b * S),
  data = rec_data,
  start = ricker_starts,
  lower = c(a = 0, b = 0)
)

# 4. ПОДГОНКА МОДЕЛИ БИВЕРТОНА-ХОЛТА -------------------------------------------
a_start <- mean(rec_data$R[rec_data$S < quantile(rec_data$S, 0.25)] / 
                rec_data$S[rec_data$S < quantile(rec_data$S, 0.25)], na.rm = TRUE)

if (is.na(a_start) || a_start <= 0) a_start <- 0.001

bh_model <- minpack.lm::nlsLM(
  R ~ (a * S) / (1 + b * S),
  data = rec_data,
  start = list(a = a_start, b = 0.0001),
  lower = c(a = 0.0001, b = 0.00001),
  control = nls.lm.control(maxiter = 200)
)

# 5. ОЦЕНКА КАЧЕСТВА МОДЕЛЕЙ ---------------------------------------------------

calculate_R2 <- function(model, data) {
  predicted <- predict(model, newdata = data)
  residuals <- data$R - predicted
  SSE <- sum(residuals^2)
  SST <- sum((data$R - mean(data$R))^2)
  R2 <- 1 - (SSE / SST)
  n <- nrow(data)
  p <- length(coef(model))
  adj_R2 <- 1 - ((n - 1) / (n - p - 1)) * (1 - R2)
  return(list(R2 = R2, adj_R2 = adj_R2))
}

calculate_pvalue <- function(model, data) {
  predicted <- predict(model, newdata = data)
  residuals <- data$R - predicted
  SSE <- sum(residuals^2)
  SST <- sum((data$R - mean(data$R))^2)
  SSR <- SST - SSE
  n <- nrow(data)
  p <- length(coef(model))
  F_stat <- (SSR / (p - 1)) / (SSE / (n - p))
  p_value <- pf(F_stat, df1 = p - 1, df2 = n - p, lower.tail = FALSE)
  return(p_value)
}

ricker_r2 <- calculate_R2(ricker_model, rec_data)
bh_r2 <- calculate_R2(bh_model, rec_data)

ricker_p <- calculate_pvalue(ricker_model, rec_data)
bh_p <- calculate_pvalue(bh_model, rec_data)

cat("AIC Рикера:", AIC(ricker_model), "\n")
cat("AIC Бивертона-Холта:", AIC(bh_model), "\n")

# 6. ВИЗУАЛИЗАЦИЯ РЕЗУЛЬТАТОВ --------------------------------------------------
new_data <- data.frame(S = seq(min(rec_data$S), max(rec_data$S), length.out = 100))
ricker_ci <- investr::predFit(ricker_model, newdata = new_data, interval = "confidence")
bh_ci <- investr::predFit(bh_model, newdata = new_data, interval = "confidence")

plot_data <- new_data %>%
  mutate(
    ricker_pred = predict(ricker_model, newdata = .),
    ricker_lwr = ricker_ci[, "lwr"],
    ricker_upr = ricker_ci[, "upr"],
    bh_pred = predict(bh_model, newdata = .),
    bh_lwr = bh_ci[, "lwr"],
    bh_upr = bh_ci[, "upr"]
  )

ggplot() +
  geom_point(data = rec_data, aes(x = S, y = R), 
             color = "darkgray", size = 3, alpha = 0.7) +
  geom_ribbon(data = plot_data, aes(x = S, ymin = ricker_lwr, ymax = ricker_upr), 
              fill = "red", alpha = 0.2) +
  geom_ribbon(data = plot_data, aes(x = S, ymin = bh_lwr, ymax = bh_upr), 
              fill = "blue", alpha = 0.2) +
  geom_line(data = plot_data, aes(x = S, y = ricker_pred), 
            color = "red", linewidth = 1.2) +
  geom_line(data = plot_data, aes(x = S, y = bh_pred), 
            color = "blue", linewidth = 1.2, linetype = "dashed") +
  labs(
    title = "Сравнение моделей запас-пополнение",
    subtitle = paste0(
      "Рикер: R² = ", round(ricker_r2$R2, 2), ", p = ", format.pval(ricker_p, digits = 3),
      " | Бивертон-Холт: R² = ", round(bh_r2$R2, 2), ", p = ", format.pval(bh_p, digits = 3)
    ),
    x = "Нерестовый запас",
    y = "Пополнение"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "none"
  )








# 7. СРАВНЕНИЕ С ДРУГИМИ ТИПАМИ МОДЕЛЕЙ ----------------------------------------

lm_model <- lm(R3haddock ~ ., data = model_data)
glm_model <- glm(R3haddock ~ ., family = Gamma(link = "log"), data = model_data)
gam_model <- mgcv::gam(R3haddock ~ s(haddock68) + s(codTSB) + s(T12) + s(I5) + s(NAOspring),
               data = model_data, method = "REML")

model_comparison <- data.frame(
  Model = c("Рикер", "Бивертон-Холт", "LM", "GLM", "GAM"),
  AIC = c(AIC(ricker_model), AIC(bh_model), AIC(lm_model), AIC(glm_model), AIC(gam_model)),
  R2 = c(
    ricker_r2$R2, 
    bh_r2$R2, 
    summary(lm_model)$r.squared,
    cor(model_data$R3haddock, predict(glm_model, type = "response"))^2,
    summary(gam_model)$r.sq
  )
)

print(model_comparison)

# 8. ИНТЕРПРЕТАЦИЯ РЕЗУЛЬТАТОВ -------------------------------------------------
cat("\nПараметры модели Рикера:")
cat("\na =", coef(ricker_model)[1], "- максимальная продукция потомства")
cat("\nb =", coef(ricker_model)[2], "- коэффициент плотностной зависимости")

cat("\n\nПараметры модели Бивертона-Холта:")
cat("\na =", coef(bh_model)[1], "- максимальное пополнение на особь")
cat("\nb =", coef(bh_model)[2], "- коэффициент внутривидовой конкуренции")

# ==============================================================================
# 7. СРАВНЕНИЕ С ДРУГИМИ ТИПАМИ МОДЕЛЕЙ ----------------------------------------

# Построение линейной модели LM
lm_model <- lm(R3haddock ~ ., data = model_data)

# Диагностика
par(mfrow = c(2, 2))
plot(lm_model)
vif(lm_model)  # Проверка мультиколлинеарности

# Интерпретация
summary(lm_model)

# Построение обобщенной линейной модели GLM
glm_model <- glm(R3haddock ~ ., 
                family = Gamma(link = "log"), 
                data = model_data)
summary(glm_model)

# Построение обобщенной аддитивной модели GАM
library(mgcv)
gam_model <- gam(R3haddock ~ 
                 s(codTSB) + 
                 s(T12) + 
                 s(I5) + 
                 s(NAOspring) + 
                 s(haddock68),
               data = model_data,
               method = "REML")
summary(gam_model)
plot(gam_model, pages = 1, residuals = TRUE)



# Таблица сравнения моделей
# Сравнение моделей
model_comparison <- data.frame(
  Model = c("Рикер", "Бивертон-Холт", "LM", "GLM", "GAM"),
  AIC = c(AIC(ricker_model), AIC(bh_model), AIC(lm_model), AIC(glm_model), AIC(gam_model)),
  R2 = c(ricker_r2$R2, bh_r2$R2, summary(lm_model)$r.squared, 
         cor(model_data$R3haddock, predict(glm_model))^2, 
         summary(gam_model)$r.sq),  # Используем summary(gam_model)$r.sq для R^2
  Adj_R2 = c(ricker_r2$adj_R2, bh_r2$adj_R2, summary(lm_model)$adj.r.squared, NA, 
             summary(gam_model)$r.sq),  # Используем summary(gam_model)$r.sq для Adjusted R^2
  RMSE = c(RMSE(predict(ricker_model), rec_data$R), 
           RMSE(predict(bh_model), rec_data$R),
           RMSE(predict(lm_model), model_data$R3haddock),
           RMSE(predict(glm_model, type = "response"), model_data$R3haddock),
           RMSE(predict(gam_model, type = "response"), model_data$R3haddock))
)

# Вывод таблицы
print(model_comparison)



# ==============================================================================
# ВИЗУАЛИЗАЦИЯ ВСЕХ МОДЕЛЕЙ НА ОДНОМ ГРАФИКЕ 
# ==============================================================================

# Фиксируем другие предикторы на их средних значениях (исключая haddock68)
mean_values <- model_data %>%
  select(-R3haddock, -haddock68) %>%
  summarise(across(everything(), ~ mean(.x, na.rm = TRUE)))

# Расширяем new_data, добавляя средние значения других предикторов
new_data_full <- new_data %>%
  bind_cols(mean_values[rep(1, nrow(new_data)), ]) %>%
  rename(haddock68 = S)  # Переименовываем S в haddock68 для совместимости

# Получаем предсказания для всех моделей
new_data_full <- new_data_full %>%
  mutate(
    # Предсказания для моделей запаса-пополнения
    ricker_pred = predict(ricker_model, newdata = data.frame(S = haddock68)),
    bh_pred = predict(bh_model, newdata = data.frame(S = haddock68)),
    
    # Предсказания для линейной модели (LM)
    lm_pred = predict(lm_model, newdata = .),
    
    # Предсказания для обобщенной линейной модели (GLM)
    glm_pred = predict(glm_model, newdata = ., type = "response"),
    
    # Предсказания для обобщенной аддитивной модели (GAM)
    gam_pred = predict(gam_model, newdata = ., type = "response")
  )

# Создаем длинный формат данных для ggplot
plot_data <- new_data_full %>%
  select(haddock68, ricker_pred, bh_pred, lm_pred, glm_pred, gam_pred) %>%
  pivot_longer(
    cols = -haddock68,
    names_to = "model",
    values_to = "prediction"
  ) %>%
  mutate(
    model = case_when(
      model == "ricker_pred" ~ "Рикер",
      model == "bh_pred" ~ "Бивертон-Холт",
      model == "lm_pred" ~ "Линейная (LM)",
      model == "glm_pred" ~ "Обобщенная линейная (GLM)",
      model == "gam_pred" ~ "Обобщенная аддитивная (GAM)",
      TRUE ~ model
    )
  )

# Создаем палитру цветов для моделей
model_colors <- c(
  "Рикер" = "#E41A1C",          # Красный
  "Бивертон-Холт" = "#377EB8",  # Синий
  "Линейная (LM)" = "#4DAF4A",  # Зеленый
  "Обобщенная линейная (GLM)" = "#984EA3", # Фиолетовый
  "Обобщенная аддитивная (GAM)" = "#FF7F00" # Оранжевый
)

# Создаем график
ggplot() +
  # Точки исходных данных
  geom_point(data = rec_data, aes(x = S, y = R), 
             color = "darkgray", size = 2.5, alpha = 0.7) +
  
  # Линии предсказаний моделей
  geom_line(data = plot_data, 
            aes(x = haddock68, y = prediction, color = model, linetype = model),
            linewidth = 1.2) +
  
  # Настройка цветов и типов линий
  scale_color_manual(values = model_colors) +
  scale_linetype_manual(values = c(
    "Рикер" = "solid",
    "Бивертон-Холт" = "dashed",
    "Линейная (LM)" = "dotdash",
    "Обобщенная линейная (GLM)" = "longdash",
    "Обобщенная аддитивная (GAM)" = "twodash"
  )) +
  
  # Подписи и темы
  labs(
    title = "Сравнение моделей зависимости пополнения от нерестового запаса",
    subtitle = "Фиксация других предикторов на средних значениях",
    x = "Нерестовый запас (тыс. тонн)",
    y = "Пополнение (млн особей)",
    color = "Модель",
    linetype = "Модель"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray30"),
    axis.title = element_text(size = 12),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "gray80", fill = NA, linewidth = 0.5)
  ) +
  guides(
    color = guide_legend(nrow = 2, byrow = TRUE),
    linetype = guide_legend(nrow = 2, byrow = TRUE)
  )

