# Очистка рабочей среды и установка рабочей директории
rm(list = ls())
setwd("C:/SNOW/")

# Загрузка необходимых библиотек
library(mclust)
library(ggplot2)
library(dplyr)
library(cluster)
library(teigen)
library(gridExtra)
library(ggpubr)

# Загрузка данных
DATA <- readxl::read_excel("SNOW.xlsx", sheet = "ECO")

# Подготовка данных для кластеризации
clustering_data <- DATA %>%
  dplyr::select(CW, chela) %>%
  na.omit()

# Проверка структуры данных
str(clustering_data)

# Стандартизация данных
scaled_data <- scale(clustering_data)

# =====================================================================
# ROBUST GMM (устойчивая к выбросам) - t-распределение смеси
# =====================================================================
cat("========== ROBUST GMM (t-distribution mixture) ==========\n")

# t-распределение смеси для робастности
set.seed(123)
robust_gmm <- teigen(clustering_data, Gs = 2, models = "UUUU", verbose = FALSE)

robust_clusters <- robust_gmm$classification

# Добавление кластеров к данным
clustering_data$Cluster <- factor(robust_clusters)

# =====================================================================
# ДИАГНОСТИКА МОДЕЛИ
# =====================================================================
cat("\n========== ДИАГНОСТИКА ROBUST GMM ==========\n\n")

# Параметры модели
cat("ПАРАМЕТРЫ МОДЕЛИ:\n")
cat("Логарифмическое правдоподобие:", robust_gmm$logl, "\n")
cat("BIC:", robust_gmm$bic, "\n")
cat("Число итераций:", robust_gmm$iter, "\n\n")

# Параметры компонент
cat("ПАРАМЕТРЫ КОМПОНЕНТ:\n")
for(k in 1:2) {
  cat("Кластер", k, ":\n")
  cat("  Размер:", sum(robust_clusters == k), 
      "(", round(100 * sum(robust_clusters == k) / length(robust_clusters), 1), "%)\n")
  
  # Получение параметров средних из scaled данных и преобразование обратно
  cw_mean <- robust_gmm$parameters$mean[1, k] * attr(scaled_data, "scaled:scale")[1] + attr(scaled_data, "scaled:center")[1]
  chela_mean <- robust_gmm$parameters$mean[2, k] * attr(scaled_data, "scaled:scale")[2] + attr(scaled_data, "scaled:center")[2]
  
  cat("  Среднее CW:", round(cw_mean, 2), "\n")
  cat("  Среднее chela:", round(chela_mean, 2), "\n")
  cat("  Степени свободы:", round(robust_gmm$parameters$df[k], 2), "\n\n")
}

# Анализ неопределенности классификации
if (!is.null(robust_gmm$fuzzy) && length(dim(robust_gmm$fuzzy)) > 0) {
  uncertainty <- 1 - apply(robust_gmm$fuzzy, 1, max)
  cat("АНАЛИЗ НЕОПРЕДЕЛЕННОСТИ:\n")
  cat("Средняя неопределенность:", round(mean(uncertainty), 4), "\n")
  cat("Максимальная неопределенность:", round(max(uncertainty), 4), "\n")
  cat("Точек с неопределенностью > 0.1:", sum(uncertainty > 0.1), "\n\n")
  
  # Добавление неопределенности к данным для визуализации
  clustering_data$uncertainty <- uncertainty
} else {
  cat("АНАЛИЗ НЕОПРЕДЕЛЕННОСТИ: Матрица нечеткости не доступна\n\n")
  clustering_data$uncertainty <- NA
}

# Силуэтный анализ
dist_matrix <- dist(scaled_data)
silhouette_score <- silhouette(robust_clusters, dist_matrix)
cat("СИЛУЭТНЫЙ АНАЛИЗ:\n")
cat("Средний силуэтный коэффициент:", round(mean(silhouette_score[, 3]), 4), "\n")

# =====================================================================
# ВИЗУАЛИЗАЦИЯ РЕЗУЛЬТАТОВ
# =====================================================================
cat("\n========== ВИЗУАЛИЗАЦИЯ РЕЗУЛЬТАТОВ ==========\n")

# Основной график кластеризации
p_main <- ggplot(clustering_data, aes(x = CW, y = chela, color = Cluster)) +
  geom_point(alpha = 0.7, size = 1.5) +
  #stat_ellipse(type = "t", level = 0.95, linewidth = 1) +
  scale_color_manual(values = c("blue", "red")) +
  labs(title = "Robust GMM кластеризация (t-распределение)",
       subtitle = paste0("BIC = ", round(robust_gmm$bic, 2),
                        ", Silhouette = ", round(mean(silhouette_score[, 3]), 4)),
       x = "CW", y = "Chela") +
  theme_minimal() +
  theme(legend.position = "bottom")

# График неопределенности
if (!all(is.na(clustering_data$uncertainty))) {
  p_uncertainty <- ggplot(clustering_data, aes(x = CW, y = chela, color = uncertainty)) +
    geom_point(alpha = 0.7, size = 1.5) +
    scale_color_gradient(low = "blue", high = "red", limits = c(0, 0.5)) +
    labs(title = "Неопределенность классификации",
         subtitle = paste0("Средняя = ", round(mean(clustering_data$uncertainty, na.rm = TRUE), 4)),
         x = "CW", y = "Chela",
         color = "Неопределенность") +
    theme_minimal() +
    theme(legend.position = "bottom")
}

# Распределение неопределенности
if (!all(is.na(clustering_data$uncertainty))) {
  p_uncertainty_dist <- ggplot(clustering_data, aes(x = uncertainty)) +
    geom_histogram(bins = 30, fill = "lightblue", color = "black") +
    labs(title = "Распределение неопределенности",
         x = "Неопределенность", y = "Частота") +
    theme_minimal()
}

# Распределение переменных по кластерам
p_cw_dist <- ggplot(clustering_data, aes(x = CW, fill = Cluster)) +
  geom_density(alpha = 0.5) +
  labs(title = "Распределение CW по кластерам", x = "CW", y = "Плотность") +
  theme_minimal() +
  theme(legend.position = "none")

p_chela_dist <- ggplot(clustering_data, aes(x = chela, fill = Cluster)) +
  geom_density(alpha = 0.5) +
  labs(title = "Распределение Chela по кластерам", x = "Chela", y = "Плотность") +
  theme_minimal() +
  theme(legend.position = "none")

# Силуэтный график
silhouette_df <- data.frame(
  cluster = as.factor(silhouette_score[, 1]),
  silhouette_width = silhouette_score[, 3],
  observation = 1:nrow(silhouette_score)
) %>%
  arrange(cluster, silhouette_width)

p_silhouette <- ggplot(silhouette_df, aes(x = observation, y = silhouette_width, fill = cluster)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("blue", "red")) +
  labs(title = "Силуэтный график",
       subtitle = paste0("Средняя ширина силуэта: ", round(mean(silhouette_score[, 3]), 3)),
       x = "Наблюдения", y = "Ширина силуэта") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

# Компоновка графиков
if (!all(is.na(clustering_data$uncertainty))) {
  top_plots <- grid.arrange(p_main, p_uncertainty, ncol = 2)
  middle_plots <- grid.arrange(p_uncertainty_dist, p_silhouette, ncol = 2)
  bottom_plots <- grid.arrange(p_cw_dist, p_chela_dist, ncol = 2)
  
  # Объединение всех графиков
  grid.arrange(top_plots, middle_plots, bottom_plots, nrow = 3, heights = c(2, 2, 1.5))
} else {
  top_plots <- grid.arrange(p_main, p_silhouette, ncol = 2)
  bottom_plots <- grid.arrange(p_cw_dist, p_chela_dist, ncol = 2)
  
  # Объединение всех графиков
  grid.arrange(top_plots, bottom_plots, nrow = 2, heights = c(2, 1.5))
}

# =====================================================================
# СТАТИСТИЧЕСКАЯ СВОДКА
# =====================================================================
cat("\n========== СТАТИСТИЧЕСКАЯ СВОДКА КЛАСТЕРОВ ==========\n\n")

final_stats <- clustering_data %>%
  group_by(Cluster) %>%
  summarise(
    n = n(),
    Proportion = round(n() / nrow(clustering_data) * 100, 1),
    Mean_CW = round(mean(CW), 2),
    SD_CW = round(sd(CW), 2),
    Mean_chela = round(mean(chela), 2),
    SD_chela = round(sd(chela), 2),
    .groups = 'drop'
  )

print(final_stats)

# Дополнительная статистика по неопределенности
if (!all(is.na(clustering_data$uncertainty))) {
  uncertainty_stats <- clustering_data %>%
    group_by(Cluster) %>%
    summarise(
      Mean_Uncertainty = round(mean(uncertainty), 4),
      Max_Uncertainty = round(max(uncertainty), 4),
      High_Uncertainty_Count = sum(uncertainty > 0.1),
      .groups = 'drop'
    )
  
  cat("\nСТАТИСТИКА НЕОПРЕДЕЛЕННОСТИ ПО КЛАСТЕРАМ:\n")
  print(uncertainty_stats)


###################ВЫВОД ИЛИ СОХРАНЕНИЕ КЛАССИФИЦИРУЕМЫХ КРАБОВ#################

str(clustering_data)
}