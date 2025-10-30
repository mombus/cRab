# ПРАКТИЧЕСКОЕ ЗАНЯТИЕ: ГРАДИЕНТНЫЙ БУСТИНГ ДЛЯ ПРОГНОЗИРОВАНИЯ БИОМАССЫ РЫБ
# Прогноз биомассы леща в водохранилище по гидрологическим и кормовым показателям

# Установка и загрузка необходимых пакетов
if(!require(pacman)) install.packages("pacman")
pacman::p_load(tidyverse, xgboost, caret, ggplot2, patchwork, lubridate, corrplot, DiagrammeR)

# Создаем реалистичный датасет по биомассе леща в водохранилище
set.seed(123)
n_years <- 15
observations_per_year <- 12

fishery_data <- expand.grid(
  year = 1:n_years,
  month = 1:observations_per_year
) %>%
  mutate(
    observation_id = row_number(),
    # Температура воды (°C) - сезонные колебания
    temperature = 4 + 15 * sin(2 * pi * (month - 7) / 12) + rnorm(n(), 0, 2),
    # Содержание кислорода (мг/л) - зависит от температуры
    oxygen = 14 - 0.4 * temperature + rnorm(n(), 0, 0.5),
    # Биомасса зообентоса (г/м²)
    zoobenthos_biomass = 15 + 5 * sin(2 * pi * month / 12) + 0.3 * year + rnorm(n(), 0, 3),
    # Прозрачность воды (м)
    water_transparency = 2 + 0.5 * sin(2 * pi * month / 12) - 0.02 * year + rnorm(n(), 0, 0.3),
    # pH воды
    pH = 7.0 + 0.5 * sin(2 * pi * month / 12) + rnorm(n(), 0, 0.2),
    # Общая минерализация (мг/л)
    mineralization = 200 + 10 * sin(2 * pi * month / 12) + rnorm(n(), 0, 15)
  )

# Создаем сложные нелинейные зависимости для биомассы леща
fishery_data <- fishery_data %>%
  mutate(
    # Базовый тренд с учетом года
    base_trend = 100 + 2 * year - 0.1 * year^2,
    
    # Эффект температуры: оптимальный диапазон 18-22°C
    temp_effect = case_when(
      temperature < 10 ~ -30,
      temperature >= 10 & temperature < 18 ~ 20 * (temperature - 10)/8,
      temperature >= 18 & temperature <= 22 ~ 20,
      temperature > 22 & temperature <= 26 ~ 20 - 5 * (temperature - 22),
      temperature > 26 ~ -20
    ),
    
    # Эффект кислорода: критический порог при 4 мг/л
    oxygen_effect = ifelse(oxygen > 4, 15 * (1 - exp(-0.5 * (oxygen - 4))), -40),
    
    # Эффект кормовой базы с насыщением
    food_effect = 30 * (1 - exp(-0.1 * zoobenthos_biomass)),
    
    # Взаимодействие температура-кислород
    temp_oxygen_interaction = ifelse(temperature > 20 & oxygen < 6, -15, 0),
    
    # Сезонный эффект
    seasonal_effect = 10 * cos(2 * pi * (month - 5) / 12),
    
    # Шум
    noise = rnorm(n(), 0, 10)
  ) %>%
  mutate(
    # Итоговая биомасса леща (кг/га)
    bream_biomass = base_trend + temp_effect + oxygen_effect + food_effect + 
                    temp_oxygen_interaction + seasonal_effect + noise,
    bream_biomass = round(pmax(50, bream_biomass), 1)
  )

cat("Первые 6 строк данных:\n")
print(head(fishery_data))

cat("\nОписательная статистика биомассы леща:\n")
print(summary(fishery_data$bream_biomass))

# Визуализация исходных данных
p1 <- ggplot(fishery_data, aes(x = temperature, y = bream_biomass)) +
  geom_point(alpha = 0.6, color = "blue") +
  geom_smooth(method = "loess", color = "red", se = FALSE) +
  labs(title = "Зависимость биомассы от температуры",
       x = "Температура воды (°C)", y = "Биомасса леща (кг/га)") +
  theme_minimal()

p2 <- ggplot(fishery_data, aes(x = zoobenthos_biomass, y = bream_biomass)) +
  geom_point(alpha = 0.6, color = "darkgreen") +
  geom_smooth(method = "loess", color = "red", se = FALSE) +
  labs(title = "Зависимость биомассы от кормовой базы",
       x = "Биомасса зообентоса (г/м²)", y = "Биомасса леща (кг/га)") +
  theme_minimal()

p3 <- ggplot(fishery_data, aes(x = oxygen, y = bream_biomass)) +
  geom_point(alpha = 0.6, color = "purple") +
  geom_smooth(method = "loess", color = "red", se = FALSE) +
  labs(title = "Зависимость биомассы от кислорода",
       x = "Содержание кислорода (мг/л)", y = "Биомасса леща (кг/га)") +
  theme_minimal()

p4 <- fishery_data %>%
  group_by(month) %>%
  summarise(mean_biomass = mean(bream_biomass)) %>%
  ggplot(aes(x = month, y = mean_biomass)) +
  geom_line(color = "brown", linewidth = 1) +
  geom_point(color = "brown", size = 2) +
  labs(title = "Сезонная динамика биомассы",
       x = "Месяц", y = "Средняя биомасса (кг/га)") +
  theme_minimal()

# Объединяем графики
(p1 + p2) / (p3 + p4)

# Матрица корреляций
cor_matrix <- fishery_data %>%
  select(temperature, oxygen, zoobenthos_biomass, water_transparency, pH, mineralization, bream_biomass) %>%
  cor()

corrplot(cor_matrix, method = "color", type = "upper", 
         title = "Корреляции между переменными", mar = c(0,0,1,0))

# Подготовка данных для моделирования
features <- c("temperature", "oxygen", "zoobenthos_biomass", 
              "water_transparency", "pH", "mineralization", "month")

# Добавляем пропуски в данные
set.seed(123)
missing_indices <- sample(1:nrow(fishery_data), 40)
fishery_data$oxygen[missing_indices[1:20]] <- NA
fishery_data$zoobenthos_biomass[missing_indices[21:40]] <- NA

# Удаляем пропуски для обучения
model_data <- fishery_data %>%
  select(all_of(features), bream_biomass) %>%
  na.omit()

# Разделяем данные на обучающую и тестовую выборки
set.seed(123)
train_index <- createDataPartition(model_data$bream_biomass, p = 0.7, list = FALSE)
train_data <- model_data[train_index, ]
test_data <- model_data[-train_index, ]

# Подготовка матриц для XGBoost
xgb_train <- xgb.DMatrix(
  data = as.matrix(train_data[, features]),
  label = train_data$bream_biomass
)

xgb_test <- xgb.DMatrix(
  data = as.matrix(test_data[, features]),
  label = test_data$bream_biomass
)

# Параметры модели для регрессии
params <- list(
  objective = "reg:squarederror",
  eta = 0.1,
  max_depth = 6,
  gamma = 0.1,
  subsample = 0.8,
  colsample_bytree = 0.8,
  lambda = 1,
  alpha = 0.1
)

# Обучение модели
xgb_model <- xgb.train(
  params = params,
  data = xgb_train,
  nrounds = 500,
  watchlist = list(train = xgb_train, test = xgb_test),
  print_every_n = 50,
  early_stopping_rounds = 20
)

# Предсказания на тестовой выборке
predictions <- predict(xgb_model, xgb_test)

# Оценка качества модели
results <- data.frame(
  Actual = test_data$bream_biomass,
  Predicted = predictions
)

# Метрики качества
rmse <- sqrt(mean((results$Actual - results$Predicted)^2))
mae <- mean(abs(results$Actual - results$Predicted))
r_squared <- cor(results$Actual, results$Predicted)^2

cat("\n=== РЕЗУЛЬТАТЫ МОДЕЛИ ===\n")
cat("RMSE (среднеквадратичная ошибка):", round(rmse, 2), "\n")
cat("MAE (средняя абсолютная ошибка):", round(mae, 2), "\n")
cat("R² (коэффициент детерминации):", round(r_squared, 3), "\n")
cat("Количество деревьев в финальной модели:", xgb_model$niter, "\n")
cat("Лучшая итерация:", xgb_model$best_iteration, "\n")

# Визуализация важности признаков
importance_matrix <- xgb.importance(
  feature_names = features, 
  model = xgb_model
)

# Строим график важности признаков
importance_plot <- importance_matrix %>%
  ggplot(aes(x = reorder(Feature, Gain), y = Gain)) +
  geom_col(fill = "steelblue", alpha = 0.8) +
  coord_flip() +
  labs(title = "Важность признаков для прогноза биомассы леща",
       subtitle = "На основе вклада в уменьшение ошибки предсказания (Gain)",
       x = "Признаки",
       y = "Важность (Gain)") +
  theme_minimal()

print(importance_plot)

# График предсказаний vs фактические значения
p_pred <- ggplot(results, aes(x = Actual, y = Predicted)) +
  geom_point(alpha = 0.6, color = "blue") +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", linewidth = 1) +
  geom_smooth(method = "lm", color = "green", se = FALSE) +
  labs(title = "Предсказания vs Фактические значения",
       subtitle = paste("R² =", round(r_squared, 3)),
       x = "Фактическая биомасса (кг/га)",
       y = "Предсказанная биомасса (кг/га)") +
  theme_minimal()

# График остатков
p_resid <- ggplot(results, aes(x = Predicted, y = Actual - Predicted)) +
  geom_point(alpha = 0.6, color = "orange") +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed", linewidth = 1) +
  geom_smooth(method = "loess", color = "blue", se = FALSE) +
  labs(title = "Анализ остатков",
       x = "Предсказанная биомасса (кг/га)",
       y = "Остатки (Факт - Прогноз)") +
  theme_minimal()

# Объединяем графики предсказаний и остатков
p_pred + p_resid

# График динамики обучения
eval_log <- xgb_model$evaluation_log
p_learning <- ggplot(eval_log, aes(x = iter)) +
  geom_line(aes(y = train_rmse, color = "Обучающая"), linewidth = 1) +
  geom_line(aes(y = test_rmse, color = "Тестовая"), linewidth = 1) +
  labs(title = "Динамика обучения градиентного бустинга",
       subtitle = "Изменение RMSE с ростом числа деревьев",
       x = "Количество деревьев",
       y = "RMSE",
       color = "Выборка") +
  scale_color_manual(values = c("Обучающая" = "blue", "Тестовая" = "red")) +
  theme_minimal()

print(p_learning)

# Функция для прогнозирования биомассы с правильным форматированием вывода
predict_biomass <- function(temperature, oxygen, zoobenthos, transparency, ph, mineral, month) {
  new_data <- data.frame(
    temperature = temperature,
    oxygen = oxygen,
    zoobenthos_biomass = zoobenthos,
    water_transparency = transparency,
    pH = ph,
    mineralization = mineral,
    month = month
  )
  
  new_matrix <- xgb.DMatrix(as.matrix(new_data[, features]))
  prediction <- predict(xgb_model, new_matrix)
  
  # Определяем качество условий
  conditions_quality <- case_when(
    temperature >= 18 & temperature <= 22 & oxygen >= 6 & zoobenthos >= 20 ~ "ОПТИМАЛЬНЫЕ",
    temperature < 10 | temperature > 26 | oxygen < 4 | zoobenthos < 10 ~ "КРИТИЧЕСКИЕ",
    TRUE ~ "УДОВЛЕТВОРИТЕЛЬНЫЕ"
  )
  
  # Создаем разделители
  separator <- paste(rep("=", 50), collapse = "")
  line_separator <- paste(rep("-", 50), collapse = "")
  
  cat("\n", separator, "\n", sep = "")
  cat("ПРОГНОЗ БИОМАССЫ ЛЕЩА\n")
  cat(separator, "\n", sep = "")
  cat("Условия среды:", conditions_quality, "\n")
  cat("• Температура:", temperature, "°C\n")
  cat("• Кислород:", oxygen, "мг/л\n") 
  cat("• Зообентос:", zoobenthos, "г/м²\n")
  cat("• Прозрачность:", transparency, "м\n")
  cat("• pH:", ph, "\n")
  cat("• Минерализация:", mineral, "мг/л\n")
  cat("• Месяц:", month, "\n")
  cat(line_separator, "\n", sep = "")
  cat("ПРОГНОЗИРУЕМАЯ БИОМАССА ЛЕЩА:", round(prediction, 1), "кг/га\n")
  cat(separator, "\n\n", sep = "")
  
  return(prediction)
}

# Упрощенная функция для анализа чувствительности (без вывода)
quiet_predict_biomass <- function(temperature, oxygen, zoobenthos, transparency, ph, mineral, month) {
  new_data <- data.frame(
    temperature = temperature,
    oxygen = oxygen,
    zoobenthos_biomass = zoobenthos,
    water_transparency = transparency,
    pH = ph,
    mineralization = mineral,
    month = month
  )
  new_matrix <- xgb.DMatrix(as.matrix(new_data[, features]))
  predict(xgb_model, new_matrix)
}

# Примеры прогнозирования для разных сценариев
cat("\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
cat("ПРИМЕРЫ ПРОГНОЗИРОВАНИЯ БИОМАССЫ ЛЕЩА\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

predict_biomass(temperature = 20, oxygen = 8, zoobenthos = 25, 
                transparency = 2.5, ph = 7.5, mineral = 210, month = 7)

predict_biomass(temperature = 26, oxygen = 3.5, zoobenthos = 18, 
                transparency = 1.2, ph = 8.2, mineral = 250, month = 8)

predict_biomass(temperature = 2, oxygen = 12, zoobenthos = 12, 
                transparency = 3.0, ph = 7.0, mineral = 190, month = 1)

# Анализ чувствительности к температуре
cat("\nАНАЛИЗ ЧУВСТВИТЕЛЬНОСТИ К ТЕМПЕРАТУРЕ\n")
temp_range <- seq(2, 28, by = 2)

temp_sensitivity <- map_dbl(temp_range, ~ quiet_predict_biomass(
  temperature = ., oxygen = 8, zoobenthos = 20, 
  transparency = 2.0, ph = 7.5, mineral = 200, month = 6
))

sensitivity_df <- data.frame(
  temperature = temp_range,
  predicted_biomass = temp_sensitivity
)

# Визуализация чувствительности
ggplot(sensitivity_df, aes(x = temperature, y = predicted_biomass)) +
  geom_line(color = "red", linewidth = 1.5) +
  geom_point(color = "red", size = 3) +
  geom_vline(xintercept = c(18, 22), linetype = "dashed", color = "blue", alpha = 0.7) +
  annotate("text", x = 20, y = max(temp_sensitivity), 
           label = "Оптимальный\nдиапазон", color = "blue", size = 3) +
  labs(title = "Чувствительность биомассы леща к температуре",
       subtitle = "Прогноз при фиксированных прочих условиях",
       x = "Температура воды (°C)",
       y = "Прогнозируемая биомасса (кг/га)") +
  theme_minimal()

# Финальные выводы
cat("\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
cat("ИТОГИ МОДЕЛИРОВАНИЯ\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
cat("• Градиентный бустинг успешно уловил сложные нелинейные зависимости\n")
cat("• Модель демонстрирует хорошую точность прогнозирования\n") 
cat("• Наиболее важные признаки:", paste(head(importance_matrix$Feature, 3), collapse = ", "), "\n")
cat("• Модель готова для использования в практической гидробиологии\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

# Дополнительный анализ: сравнение с линейной регрессией
cat("\nСРАВНЕНИЕ С ЛИНЕЙНОЙ РЕГРЕССИЕЙ\n")
lm_model <- lm(bream_biomass ~ ., data = train_data)
lm_predictions <- predict(lm_model, test_data)
lm_rmse <- sqrt(mean((lm_predictions - test_data$bream_biomass)^2))
lm_r_squared <- cor(lm_predictions, test_data$bream_biomass)^2

cat("Линейная регрессия - RMSE:", round(lm_rmse, 2), "\n")
cat("Линейная регрессия - R²:", round(lm_r_squared, 3), "\n")
cat("XGBoost превосходит линейную модель на", round((lm_rmse - rmse)/lm_rmse * 100, 1), "% по RMSE\n")

# Сохранение модели для последующего использования
# xgb.save(xgb_model, "bream_biomass_model.xgb")
# cat("\nМодель сохранена в файл: bream_biomass_model.xgb\n")