# ПРАКТИЧЕСКОЕ ЗАНЯТИЕ: БАЙЕСОВСКИЙ ПОДХОД 
# Оценка параметров роста рыбы по уравнению фон Берталанффи при недостатке данных

# Установка и загрузка необходимых пакетов
required_packages <- c("ggplot2", "dplyr", "tidyr", "gridExtra", "MASS")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(MASS)

# -------------------- ЧАСТЬ 1: ГЕНЕРАЦИЯ ДАННЫХ --------------------

set.seed(175)

# Истинные параметры
true_Linf <- 60.0
true_K <- 0.27
true_sigma <- 6.0

# Уравнение фон Берталанффи
von_bertalanffy <- function(age, Linf, K) {
  Linf * (1 - exp(-K * age))
}

# Генерация данных
ages_full <- seq(1, 15, length.out = 100)
lengths_true <- von_bertalanffy(ages_full, true_Linf, true_K)

sample_ages <- c(1, 2, 4, 6, 10, 14)
sample_lengths <- von_bertalanffy(sample_ages, true_Linf, true_K) + 
  rnorm(length(sample_ages), 0, true_sigma)

full_data <- data.frame(Age = ages_full, Length = lengths_true)
observed_data <- data.frame(Age = sample_ages, Length = sample_lengths)

# -------------------- ЧАСТЬ 2: ВИЗУАЛИЗАЦИЯ ПРОБЛЕМЫ --------------------

plot_data <- ggplot() +
  geom_line(data = full_data, aes(x = Age, y = Length, color = "Истинная кривая"), 
            linewidth = 1.5, alpha = 0.8) +
  geom_point(data = observed_data, aes(x = Age, y = Length, color = "Наблюдения"), 
             size = 4, alpha = 0.8) +
  scale_color_manual(values = c("Истинная кривая" = "#3366CC", "Наблюдения" = "#FF4444")) +
  labs(title = "ПРОБЛЕМА НЕДОСТАТКА ДАННЫХ",
       subtitle = paste("Всего", length(sample_ages), "наблюдений для оценки кривой роста"),
       y = "Длина тела (см)", x = "Возраст (лет)",
       color = "") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 16),
        plot.subtitle = element_text(color = "gray40"))

print(plot_data)

# -------------------- ЧАСТЬ 3: АПРИОРНЫЕ РАСПРЕДЕЛЕНИЯ --------------------

# Знания эксперта-ихтиолога о похожих видах
expert_Linf_mean <- 55    # Ожидаемая предельная длина
expert_Linf_sd <- 10      # Неопределенность оценки
expert_K_mean <- 0.25      # Ожидаемая константа роста  
expert_K_sd <- 0.08       # Неопределенность оценки

cat("ЭКСПЕРТНЫЕ ЗНАНИЯ ДЛЯ ОБОИХ ПАРАМЕТРОВ:\n")
cat("• Linf ~ N(mean =", expert_Linf_mean, ", sd =", expert_Linf_sd, ")\n")
cat("• K ~ N(mean =", expert_K_mean, ", sd =", expert_K_sd, ")\n")
cat("• Истинные значения: Linf =", true_Linf, ", K =", true_K, "\n\n")

# Создаем сетку параметров
param_grid <- expand.grid(
  Linf = seq(40, 70, length.out = 150),
  K = seq(0.2, 0.4, length.out = 150)
)

# Априорные распределения (нормальные)
param_grid$prior <- dnorm(param_grid$Linf, expert_Linf_mean, expert_Linf_sd) *
                   dnorm(param_grid$K, expert_K_mean, expert_K_sd)
param_grid$prior <- param_grid$prior / sum(param_grid$prior)

# -------------------- ЧАСТЬ 4: ВИЗУАЛИЗАЦИЯ АПРИОРОВ --------------------

# Маргинальные распределения - Linf
prior_linf_plot <- ggplot(data.frame(x = seq(40, 70, length.out = 100)), aes(x = x)) +
  geom_area(aes(y = dnorm(x, expert_Linf_mean, expert_Linf_sd)), 
            fill = "#3366CC", alpha = 0.6) +
  geom_line(aes(y = dnorm(x, expert_Linf_mean, expert_Linf_sd)), 
            color = "#3366CC", linewidth = 1.2) +
  geom_vline(xintercept = true_Linf, linetype = "dashed", color = "#FF4444", linewidth = 1) +
  geom_vline(xintercept = expert_Linf_mean, linetype = "dashed", color = "#3366CC", linewidth = 1) +
  annotate("text", x = true_Linf, y = 0.03, label = "Истинное значение", 
           hjust = -0.1, color = "#FF4444", size = 3) +
  annotate("text", x = expert_Linf_mean, y = 0.025, label = "Мнение эксперта", 
           hjust = 1.1, color = "#3366CC", size = 3) +
  labs(title = "АПРИОРНОЕ РАСПРЕДЕЛЕНИЕ: Linf",
       subtitle = "Предельная длина тела",
       y = "Плотность вероятности", x = "Linf (см)") +
  theme_minimal()

# Маргинальные распределения - K (ИСПРАВЛЕНО: добавлена линия мнения эксперта)
prior_k_plot <- ggplot(data.frame(x = seq(0.2, 0.4, length.out = 100)), aes(x = x)) +
  geom_area(aes(y = dnorm(x, expert_K_mean, expert_K_sd)), 
            fill = "#3366CC", alpha = 0.6) +
  geom_line(aes(y = dnorm(x, expert_K_mean, expert_K_sd)), 
            color = "#3366CC", linewidth = 1.2) +
  geom_vline(xintercept = true_K, linetype = "dashed", color = "#FF4444", linewidth = 1) +
  geom_vline(xintercept = expert_K_mean, linetype = "dashed", color = "#3366CC", linewidth = 1) +  # ДОБАВЛЕНО
  annotate("text", x = true_K, y = 4.5, label = "Истинное значение", 
           hjust = -0.1, color = "#FF4444", size = 3) +
  annotate("text", x = expert_K_mean, y = 4.0, label = "Мнение эксперта", 
           hjust = 1.1, color = "#3366CC", size = 3) +  # ДОБАВЛЕНО
  labs(title = "АПРИОРНОЕ РАСПРЕДЕЛЕНИЕ: K",
       subtitle = "Константа роста",
       y = "Плотность вероятности", x = "K (1/год)") +
  theme_minimal()

# Совместное распределение
joint_prior_plot <- ggplot(param_grid, aes(x = Linf, y = K)) +
  geom_raster(aes(fill = prior), interpolate = TRUE) +
  geom_contour(aes(z = prior), color = "white", alpha = 0.5, linewidth = 0.3) +
  annotate("point", x = true_Linf, y = true_K, color = "#FF4444", size = 3, shape = 4) +
  annotate("point", x = expert_Linf_mean, y = expert_K_mean, color = "#3366CC", size = 3, shape = 16) +
  scale_fill_gradient(low = "#E6F2FF", high = "#003366", name = "Плотность") +
  labs(title = "СОВМЕСТНОЕ АПРИОРНОЕ РАСПРЕДЕЛЕНИЕ",
       subtitle = "Синий круг - мнение эксперта, Красный крест - истинные значения",
       y = "K (1/год)", x = "Linf (см)") +
  theme_minimal()

grid.arrange(prior_linf_plot, prior_k_plot, joint_prior_plot, 
             layout_matrix = rbind(c(1, 2), c(3, 3)),
             top = "БАЙЕСОВСКИЙ ПОДХОД: АПРИОРНЫЕ ЗНАНИЯ ЭКСПЕРТА ДЛЯ ОБОИХ ПАРАМЕТРОВ")

# -------------------- ЧАСТЬ 5: ФУНКЦИЯ ПРАВДОПОДОБИЯ --------------------

calculate_likelihood <- function(Linf, K, sigma) {
  predicted <- von_bertalanffy(sample_ages, Linf, K)
  sum(dnorm(sample_lengths, predicted, sigma, log = TRUE))
}

sigma_prior <- 2.5
param_grid$log_likelihood <- apply(param_grid, 1, function(row) {
  calculate_likelihood(row["Linf"], row["K"], sigma_prior)
})

# Преобразуем в likelihood
param_grid$likelihood <- exp(param_grid$log_likelihood - max(param_grid$log_likelihood))
param_grid$likelihood <- param_grid$likelihood / sum(param_grid$likelihood, na.rm = TRUE)

# Визуализация правдоподобия
likelihood_plot <- ggplot(param_grid, aes(x = Linf, y = K)) +
  geom_raster(aes(fill = likelihood), interpolate = TRUE) +
  geom_contour(aes(z = likelihood), color = "white", alpha = 0.6, linewidth = 0.4) +
  annotate("point", x = true_Linf, y = true_K, color = "#FF4444", size = 3, shape = 4) +
  annotate("point", x = expert_Linf_mean, y = expert_K_mean, color = "#3366CC", size = 3, shape = 16) +  # ДОБАВЛЕНО
  scale_fill_gradient(low = "#E6FFE6", high = "#006600", name = "Правдоподобие") +
  labs(title = "ФУНКЦИЯ ПРАВДОПОДОБИЯ",
       subtitle = "Вероятность наблюдать данные при заданных параметрах\nКрасный крест - истинные значения, Синий круг - мнение эксперта",
       y = "K (1/год)", x = "Linf (см)") +
  theme_minimal()

print(likelihood_plot)


# -------------------- ЧАСТЬ 6: АПОСТЕРИОРНОЕ РАСПРЕДЕЛЕНИЕ --------------------

# Удаляем NA значения
param_grid <- param_grid[complete.cases(param_grid), ]

# Вычисляем апостериорное распределение
param_grid$posterior <- param_grid$prior * param_grid$likelihood
param_grid$posterior <- param_grid$posterior / sum(param_grid$posterior, na.rm = TRUE)

# Визуализация апостериорного распределения
posterior_plot <- ggplot(param_grid, aes(x = Linf, y = K)) +
  geom_raster(aes(fill = posterior), interpolate = TRUE) +
  geom_contour(aes(z = posterior), color = "white", alpha = 0.6, linewidth = 0.4) +
  annotate("point", x = true_Linf, y = true_K, color = "#FF4444", size = 3, shape = 4) +
  scale_fill_gradient(low = "#FFE6E6", high = "#CC0000", name = "Плотность") +
  labs(title = "АПОСТЕРИОРНОЕ РАСПРЕДЕЛЕНИЕ",
       subtitle = "Обновленные знания после учета данных",
       y = "K (1/год)", x = "Linf (см)") +
  theme_minimal()

print(posterior_plot)

# -------------------- ЧАСТЬ 7: МАРГИНАЛЬНЫЕ РАСПРЕДЕЛЕНИЯ --------------------

# Маргинализация
marginal_linf <- param_grid %>%
  group_by(Linf) %>%
  summarise(prior = sum(prior),
            posterior = sum(posterior))

marginal_k <- param_grid %>%
  group_by(K) %>%
  summarise(prior = sum(prior),
            posterior = sum(posterior))

# Графики маргинальных распределений
marginal_linf_plot <- ggplot(marginal_linf) +
  geom_area(aes(x = Linf, y = prior, fill = "Априорное"), alpha = 0.6) +
  geom_line(aes(x = Linf, y = prior, color = "Априорное"), linewidth = 1.2) +
  geom_area(aes(x = Linf, y = posterior, fill = "Апостериорное"), alpha = 0.6) +
  geom_line(aes(x = Linf, y = posterior, color = "Апостериорное"), linewidth = 1.2) +
  geom_vline(xintercept = true_Linf, linetype = "dashed", color = "#FF4444", linewidth = 1) +
  scale_fill_manual(values = c("Априорное" = "#3366CC", "Апостериорное" = "#CC0000"), 
                    name = "Распределение") +
  scale_color_manual(values = c("Априорное" = "#3366CC", "Апостериорное" = "#CC0000"), 
                     name = "Распределение") +
  labs(title = "МАРГИНАЛЬНОЕ РАСПРЕДЕЛЕНИЕ: Linf",
       y = "Плотность вероятности", x = "Linf (см)") +
  theme_minimal() +
  theme(legend.position = "none")  # Скрываем легенду у первого графика

marginal_k_plot <- ggplot(marginal_k) +
  geom_area(aes(x = K, y = prior, fill = "Априорное"), alpha = 0.6) +
  geom_line(aes(x = K, y = prior, color = "Априорное"), linewidth = 1.2) +
  geom_area(aes(x = K, y = posterior, fill = "Апостериорное"), alpha = 0.6) +
  geom_line(aes(x = K, y = posterior, color = "Апостериорное"), linewidth = 1.2) +
  geom_vline(xintercept = true_K, linetype = "dashed", color = "#FF4444", linewidth = 1) +
  scale_fill_manual(values = c("Априорное" = "#3366CC", "Апостериорное" = "#CC0000"), 
                    name = "Распределение") +
  scale_color_manual(values = c("Априорное" = "#3366CC", "Апостериорное" = "#CC0000"), 
                     name = "Распределение") +
  labs(title = "МАРГИНАЛЬНОЕ РАСПРЕДЕЛЕНИЕ: K",
       y = "Плотность вероятности", x = "K (1/год)") +
  theme_minimal() +
  theme(legend.position = "bottom")  # Размещаем легенду внизу у второго графика

# Объединяем графики
grid.arrange(marginal_linf_plot, marginal_k_plot, ncol = 2)

# -------------------- ЧАСТЬ 8: РАСЧЕТ СТАТИСТИК --------------------

# Точечные оценки
posterior_mean_linf <- sum(param_grid$Linf * param_grid$posterior, na.rm = TRUE)
posterior_mean_k <- sum(param_grid$K * param_grid$posterior, na.rm = TRUE)

# Стандартные отклонения
posterior_sd_linf <- sqrt(sum(param_grid$posterior * (param_grid$Linf - posterior_mean_linf)^2, na.rm = TRUE))
posterior_sd_k <- sqrt(sum(param_grid$posterior * (param_grid$K - posterior_mean_k)^2, na.rm = TRUE))

# 95% доверительные интервалы
linf_samples <- sample(param_grid$Linf, size = 10000, replace = TRUE, prob = param_grid$posterior)
k_samples <- sample(param_grid$K, size = 10000, replace = TRUE, prob = param_grid$posterior)

linf_quantiles <- quantile(linf_samples, probs = c(0.025, 0.975))
k_quantiles <- quantile(k_samples, probs = c(0.025, 0.975))

# Вывод результатов
cat("РЕЗУЛЬТАТЫ БАЙЕСОВСКОГО АНАЛИЗА:\n")
cat("================================\n\n")

cat("ТОЧЕЧНЫЕ ОЦЕНКИ:\n")
cat("• Linf:", round(posterior_mean_linf, 1), "см (истинное:", true_Linf, "см)\n")
cat("• K:", round(posterior_mean_k, 3), "/год (истинное:", true_K, "/год)\n\n")

cat("НЕОПРЕДЕЛЕННОСТЬ ОЦЕНОК:\n")
cat("• SD(Linf):", round(posterior_sd_linf, 1), "см\n")
cat("• SD(K):", round(posterior_sd_k, 3), "/год\n\n")

cat("95% ДОВЕРИТЕЛЬНЫЕ ИНТЕРВАЛЫ:\n")
cat("• Linf: [", round(linf_quantiles[1], 1), ",", round(linf_quantiles[2], 1), "] см\n")
cat("• K: [", round(k_quantiles[1], 3), ",", round(k_quantiles[2], 3), "] /год\n\n")

# -------------------- ЧАСТЬ 9: ФИНАЛЬНАЯ ВИЗУАЛИЗАЦИЯ --------------------

# Генерируем кривые роста на основе апостериорного распределения
set.seed(456)
n_curves <- 150
sample_indices <- sample(1:nrow(param_grid), n_curves, prob = param_grid$posterior)

curve_data <- data.frame()
for (i in 1:n_curves) {
  idx <- sample_indices[i]
  curve <- data.frame(
    Age = ages_full,
    Length = von_bertalanffy(ages_full, param_grid$Linf[idx], param_grid$K[idx]),
    Curve = i
  )
  curve_data <- rbind(curve_data, curve)
}

final_plot <- ggplot() +
  geom_line(data = curve_data, aes(x = Age, y = Length, group = Curve), 
            color = "#FF6666", alpha = 0.08, linewidth = 0.3) +
  geom_line(data = full_data, aes(x = Age, y = Length, color = "Истинная кривая"), 
            linewidth = 2, alpha = 0.9) +
  geom_point(data = observed_data, aes(x = Age, y = Length, color = "Наблюдения"), 
             size = 3, alpha = 0.9) +
  stat_summary(data = curve_data, aes(x = Age, y = Length), 
               fun = mean, geom = "line", color = "#CC0000", linewidth = 1.5) +
  scale_color_manual(values = c("Истинная кривая" = "#3366CC", "Наблюдения" = "#00AA00")) +
  labs(title = "БАЙЕСОВСКАЯ ОЦЕНКА КРИВОЙ РОСТА",
       subtitle = paste("Красные линии - возможные кривые (n =", n_curves, ")\nТемно-красная - средняя апостериорная оценка"),
       y = "Длина тела (см)", x = "Возраст (лет)",
       color = "") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 16),
        plot.subtitle = element_text(color = "gray40"))

print(final_plot)

cat("\nЗАКЛЮЧЕНИЕ:\n")
cat("==========\n")
cat("Байесовский подход позволяет эффективно сочетать:\n")
cat("• Экспертные знания (априорные распределения)\n")
cat("• Ограниченные данные наблюдений (правдоподобие)\n")
cat("• Количественную оценку неопределенности\n")
cat("• Интерпретируемые результаты в терминах вероятностей\n\n")