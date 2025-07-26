# ЗАГРУЗКА БИБЛИОТЕК =================================================
library(neuralnet)   # Для построения нейронных сетей
library(ggplot2)     # Для продвинутой визуализации
library(dplyr)       # Для обработки данных
library(tidyr)       # Для преобразования данных
theme_set(theme_minimal())  # Устанавливаем минималистичную тему

# РАЗДЕЛ 1: ЛИНЕЙНАЯ РЕГРЕССИЯ =======================================
### Основы экологического моделирования
w <- c(85, 90, 85, 95, 95, 135, 165, 135, 140)
lt <- c(51, 51, 52, 54, 54, 59, 59, 60, 62)
snake_data <- data.frame(lt, w)

# Построение модели
lreg <- lm(w ~ lt, data = snake_data)

# Визуализация с ggplot2
ggplot(snake_data, aes(lt, w)) +
  geom_point(size = 3, alpha = 0.8, color = "#1E88E5") +
  geom_smooth(method = "lm", formula = y ~ x, 
              se = TRUE, color = "#D81B60", fill = "#F48FB1") +
  labs(title = "Линейная регрессия: зависимость массы от длины тела",
       subtitle = "Vipera berus (обыкновенная гадюка)",
       x = "Длина тела (см)", 
       y = "Масса (г)",
       caption = "Данные: Полевые исследования Карелия, 2023") +
  annotate("text", x = 58, y = 100, 
           label = paste0("w = ", round(coef(lreg)[1], 2), " + ", 
                          round(coef(lreg)[2], 2), "*lt\nR? = ", 
                          round(summary(lreg)$r.squared, 3)),
           color = "#004D40", size = 4)

# РАЗДЕЛ 2: МНОГОФАКТОРНЫЙ АНАЛИЗ ====================================
w <- c(40, 156, 105, 85, 80, 50, 75, 48, 75, 67)
lt <- c(44, 59, 49, 50, 54, 43, 49, 42, 47, 47)
lc <- c(70, 78, 66, 90, 83, 70, 62, 75, 40, 80)
multi_data <- data.frame(lt, lc, w)

# Множественная регрессия
multi_reg <- lm(w ~ lt + lc, data = multi_data)

# Визуализация 3D-проекции
library(plotly)  # Для интерактивной 3D визуализации
plot_ly(multi_data, x = ~lt, y = ~lc, z = ~w, type = "scatter3d", 
        mode = "markers", size = 5, 
        marker = list(color = ~w, colorscale = "Viridis")) %>%
  layout(title = "Множественная регрессия: масса от длины тела и хвоста",
         scene = list(xaxis = list(title = "Длина тела (см)"),
                      yaxis = list(title = "Длина хвоста (см)"),
                      zaxis = list(title = "Масса (г)")))

# РАЗДЕЛ 3: НЕЛИНЕЙНЫЕ ЗАВИСИМОСТИ ===================================
### Степенная зависимость
power_data <- data.frame(lt = lt, w = w)

# Подбор модели
power_model <- nls(w ~ a * lt^b, start = list(a = 0.001, b = 3), data = power_data)

# Создание данных для гладкой кривой
curve_data <- data.frame(lt = seq(min(lt), max(lt), length.out = 100))
curve_data$w <- predict(power_model, newdata = curve_data)

# Визуализация
ggplot(power_data, aes(lt, w)) +
  geom_point(size = 3, color = "#FF9800") +
  geom_line(data = curve_data, aes(x = lt, y = w), 
            color = "#5E35B1", size = 1.2) +
  labs(title = "Степенная зависимость массы от длины тела",
       x = "Длина тела (см)", 
       y = "Масса (г)") +
  annotate("text", x = 45, y = 150, 
           label = paste0("w = ", round(coef(power_model)[1], 5), " * lt^", 
                          round(coef(power_model)[2], 2)),
           color = "#5E35B1", size = 5)

# РАЗДЕЛ 4: ПОРОГОВЫЕ ЭФФЕКТЫ ========================================
### Логистическая регрессия для токсикологических данных
K <- c(100, 126, 158, 200, 251, 316, 398, 501, 631, 794, 1000)
p <- c(0, 0, 0, 0, 0, 0.5, 0.5, 1, 1, 1, 1)
tox_data <- data.frame(K, p)

# Построение модели
logit_model <- glm(p ~ K, family = binomial, data = tox_data)

# Подготовка данных для гладкой кривой
tox_curve <- data.frame(K = seq(100, 1000, length.out = 200))
tox_curve$p_pred <- predict(logit_model, newdata = tox_curve, type = "response")

# Визуализация
ggplot(tox_data, aes(K, p)) +
  geom_point(size = 4, color = "#E53935") +
  geom_line(data = tox_curve, aes(x = K, y = p_pred), 
            color = "#1E88E5", size = 1.5) +
  geom_vline(xintercept = 400, linetype = "dashed", color = "#43A047") +
  annotate("text", x = 420, y = 0.8, label = "Порог LC50", color = "#43A047") +
  labs(title = "Логистическая модель токсического воздействия",
       subtitle = "Смертность Daphnia magna при разных концентрациях лигнина",
       x = "Концентрация токсиканта (мг/л)", 
       y = "Доля погибших особей") +
  scale_y_continuous(labels = scales::percent)

# РАЗДЕЛ 5: НЕЙРОННЫЕ СЕТИ ===========================================
### Сравнение разных архитектур сетей

# Простейшая сеть (аналог линейной регрессии)
nn_simple <- neuralnet(p ~ K, data = tox_data, hidden = 0)

# Сеть с одним скрытым нейроном
nn_1hidden <- neuralnet(p ~ K, data = tox_data, hidden = 1)

# Сеть с двумя скрытыми нейронами
nn_2hidden <- neuralnet(p ~ K, data = tox_data, hidden = 2)

# Создание данных для сравнения
comparison_data <- tox_data %>%
  mutate(
    simple = predict(nn_simple, tox_data),
    hidden1 = predict(nn_1hidden, tox_data),
    hidden2 = predict(nn_2hidden, tox_data),
    logit = predict(logit_model, type = "response")
  ) %>%
  pivot_longer(cols = c(simple, hidden1, hidden2, logit), 
               names_to = "model", 
               values_to = "prediction")

# Визуализация сравнения моделей
ggplot(comparison_data, aes(K, prediction, color = model)) +
  geom_line(size = 1.2) +
  geom_point(aes(y = p), color = "black", size = 3) +
  scale_color_manual(
    name = "Модель",
    values = c("#D81B60", "#1E88E5", "#004D40", "#FFC107"),
    labels = c("1 скрытый нейрон", "2 скрытых нейрона", "Логистическая", "Линейная")
  ) +
  labs(title = "Сравнение моделей для токсикологических данных",
       subtitle = "Нейронные сети vs логистическая регрессия",
       x = "Концентрация токсиканта (мг/л)", 
       y = "Предсказанная смертность") +
  theme(legend.position = "bottom")

# РАЗДЕЛ 6: КЛАССИФИКАЦИЯ ============================================
### Определение пола гадюк по морфометрии
# Загрузка данных
vipers <- read.csv("vipkar.csv")

# Функция для обучения и оценки моделей
train_evaluate_model <- function(hidden_layers) {
  model <- neuralnet(
    ns ~ lc + lt + w,
    data = vipers,
    hidden = hidden_layers,
    linear.output = FALSE
  )
  
  predictions <- predict(model, vipers)
  accuracy <- mean(round(predictions) == vipers$ns)
  
  return(list(model = model, accuracy = accuracy))
}

# Сравнение разных архитектур
architectures <- list(
  "0 нейронов" = 0,
  "1 нейрон" = 1,
  "3 нейрона" = 3,
  "2 слоя (3,2)" = c(3, 2)
)

results <- lapply(architectures, train_evaluate_model)

# Создание данных для визуализации
accuracy_data <- data.frame(
  Architecture = names(architectures),
  Accuracy = sapply(results, function(x) x$accuracy)
)

# Визуализация точности моделей
ggplot(accuracy_data, aes(x = reorder(Architecture, -Accuracy), y = Accuracy, fill = Architecture)) +
  geom_col(alpha = 0.8) +
  geom_text(aes(label = scales::percent(Accuracy, accuracy = 1)), 
            vjust = -0.5, size = 5) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Точность определения пола гадюк",
       subtitle = "Сравнение разных архитектур нейронных сетей",
       x = "Архитектура сети", 
       y = "Точность классификации") +
  theme(axis.text.x = element_text(angle = 15, hjust = 1))

# Визуализация лучшей модели
best_model <- results[["2 слоя (3,2)"]]$model
plot(best_model, rep = "best")

# РАЗДЕЛ 7: ПРОСТРАНСТВЕННОЕ МОДЕЛИРОВАНИЕ ==========================
island_data <- read.csv("kihzsdat.csv")

# Обучение модели
set.seed(123)
model_island <- neuralnet(
  vb ~ fo + me + bo,
  data = island_data,
  hidden = c(5, 3),
  linear.output = TRUE
)

# Прогнозирование
predictions <- predict(model_island, island_data)
island_data$predicted <- round(predictions)

# Визуализация сравнения наблюдений и прогноза
ggplot(island_data, aes(x = vb, y = predicted)) +
  geom_point(size = 4, alpha = 0.7, color = "#1E88E5") +
  geom_abline(slope = 1, intercept = 0, color = "#D81B60", size = 1.2) +
  geom_smooth(method = "lm", se = FALSE, color = "#43A047", linetype = "dashed") +
  labs(title = "Прогнозирование численности гадюк на островах",
       subtitle = "Сравнение наблюдаемых и предсказанных значений",
       x = "Наблюдаемая численность (категория)", 
       y = "Предсказанная численность (категория)") +
  scale_x_continuous(breaks = 1:4) +
  scale_y_continuous(breaks = 0:4) +
  coord_fixed(ratio = 1)

# Визуализация важности переменных
imp_data <- data.frame(
  Variable = c("Леса", "Луга", "Болота"),
  Importance = c(mean(abs(model_island$weights[[1]][[1]][, 1])),
                mean(abs(model_island$weights[[1]][[1]][, 2])),
                mean(abs(model_island$weights[[1]][[1]][, 3]))
))

ggplot(imp_data, aes(x = reorder(Variable, -Importance), y = Importance, fill = Variable)) +
  geom_col(alpha = 0.8) +
  scale_fill_manual(values = c("#4CAF50", "#FF9800", "#2196F3")) +
  labs(title = "Важность факторов среды обитания",
       subtitle = "Для прогнозирования численности гадюк",
       x = "Фактор", 
       y = "Относительная важность") +
  theme(legend.position = "none")