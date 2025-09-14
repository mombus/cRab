# ОПРЕДЕЛЕНИЕ ПАРАМЕТРОВ ЗРЕЛОСТИ  КРАБА-СТРИГУНА (L50 и L95)
# Анализ выполняется с использованием кластеризации и логистической регрессии


# =============================================================================
# 1. ПОДГОТОВКА СРЕДЫ И ЗАГРУЗКА БИБЛИОТЕК
# =============================================================================

# Очистка рабочей среды для исключения конфликтов с предыдущими данными
rm(list = ls())

# Установка рабочей директории (замените на путь к вашим данным)
setwd("C:/SNOW/")

# Загрузка необходимых библиотек
library(mclust)      # Для Gaussian Mixture Models
library(ggplot2)     # Для создания графиков
library(dplyr)       # Для манипуляций с данными
library(cluster)     # Для силуэтного анализа
library(teigen)      # Для t-распределенной смеси (робастная кластеризация)
library(gridExtra)   # Для компоновки графиков
library(ggpubr)      # Для публикационных графиков

# =============================================================================
# 2. ЗАГРУЗКА И ПОДГОТОВКА ДАННЫХ
# =============================================================================

# Загрузка данных из Excel-файла
DATA <- readxl::read_excel("SNOW.xlsx", sheet = "ECO")

# Проверка доступных колонок в данных
cat("Доступные колонки в данных:\n")
print(colnames(DATA))

# Фильтрация данных: используем только самцов или все данные, если пол не указан
if(!"Sex" %in% colnames(DATA)) {
  cat("Колонка 'Sex' не найдена. Используем все данные.\n")
  DATA_filtered <- DATA
} else {
  # Фильтрация только самцов (предполагается, что самки не являются объектом промысла)
  DATA_filtered <- DATA %>% filter(Sex == "M")
  cat("Отфильтрованы только самцы. Размер выборки:", nrow(DATA_filtered), "\n")
}

# Выбор переменных для анализа: ширина карапакса (CW) и размер клешни (chela)
clustering_data <- DATA_filtered %>%
  dplyr::select(CW, chela) %>%
  na.omit()  # Удаление пропущенных значений

# Проверка структуры данных
cat("Структура данных для кластеризации:\n")
str(clustering_data)
cat("Размерность данных:", dim(clustering_data), "\n")

# =============================================================================
# 3. КЛАСТЕРИЗАЦИЯ ДЛЯ РАЗДЕЛЕНИЯ НА ЗРЕЛЫХ И НЕЗРЕЛЫХ ОСОБЕЙ
# =============================================================================

# Стандартизация данных для улучшения сходимости алгоритма
scaled_data <- scale(clustering_data)

# ROBUST GMM с t-распределением (устойчив к выбросам)
cat("\n========== ВЫПОЛНЕНИЕ КЛАСТЕРИЗАЦИИ ==========\n")

# Установка seed для воспроизводимости результатов
set.seed(123)

# Кластеризация с использованием t-распределенной смеси
# Gs = 2: ожидаем 2 кластера (неполовозрелые и половозрелые)
# models = "UUUU": неограниченная ковариационная матрица
robust_gmm <- teigen(clustering_data, Gs = 2, models = "UUUU", verbose = FALSE)

# Получение результатов кластеризации
robust_clusters <- robust_gmm$classification

# Добавление информации о кластерах к данным
clustering_data$Cluster <- factor(robust_clusters)

# =============================================================================
# 4. ДИАГНОСТИКА КАЧЕСТВА КЛАСТЕРИЗАЦИИ
# =============================================================================

cat("\n========== ДИАГНОСТИКА КАЧЕСТВА КЛАСТЕРИЗАЦИИ ==========\n\n")

# Основные параметры модели
cat("ПАРАМЕТРЫ МОДЕЛИ:\n")
cat("Логарифмическое правдоподобие:", robust_gmm$logl, "\n")
cat("BIC (Байесовский информационный критерий):", robust_gmm$bic, "\n")
cat("Число итераций до сходимости:", robust_gmm$iter, "\n\n")

# Характеристики кластеров
cat("ПАРАМЕТРЫ КЛАСТЕРОВ:\n")
for(k in 1:2) {
  cat("Кластер", k, ":\n")
  cat("  Размер:", sum(robust_clusters == k), 
      "(", round(100 * sum(robust_clusters == k) / length(robust_clusters), 1), "%)\n")
  
  # Преобразование стандартизованных средних обратно в оригинальные единицы
  cw_mean <- robust_gmm$parameters$mean[1, k] * attr(scaled_data, "scaled:scale")[1] + 
    attr(scaled_data, "scaled:center")[1]
  chela_mean <- robust_gmm$parameters$mean[2, k] * attr(scaled_data, "scaled:scale")[2] + 
    attr(scaled_data, "scaled:center")[2]
  
  cat("  Средняя ширина карапакса (CW):", round(cw_mean, 2), "мм\n")
  cat("  Средний размер клешни (chela):", round(chela_mean, 2), "мм\n")
  cat("  Степени свободы t-распределения:", round(robust_gmm$parameters$df[k], 2), "\n\n")
}

# Анализ неопределенности классификации
if (!is.null(robust_gmm$fuzzy) && length(dim(robust_gmm$fuzzy)) > 0) {
  uncertainty <- 1 - apply(robust_gmm$fuzzy, 1, max)
  cat("АНАЛИЗ НЕОПРЕДЕЛЕННОСТИ КЛАССИФИКАЦИИ:\n")
  cat("Средняя неопределенность:", round(mean(uncertainty), 4), "\n")
  cat("Максимальная неопределенность:", round(max(uncertainty), 4), "\n")
  cat("Количество точек с неопределенностью > 0.1:", sum(uncertainty > 0.1), "\n\n")
  
  # Добавление неопределенности к данным для визуализации
  clustering_data$uncertainty <- uncertainty
} else {
  cat("АНАЛИЗ НЕОПРЕДЕЛЕННОСТИ: Матрица нечеткости не доступна\n\n")
  clustering_data$uncertainty <- NA
}

# Силуэтный анализ качества кластеризации
dist_matrix <- dist(scaled_data)
silhouette_score <- silhouette(robust_clusters, dist_matrix)
cat("СИЛУЭТНЫЙ АНАЛИЗ КАЧЕСТВА КЛАСТЕРИЗАЦИИ:\n")
cat("Средний силуэтный коэффициент:", round(mean(silhouette_score[, 3]), 4), "\n")
cat("Интерпретация: значения ближе к 1 указывают на лучшую кластеризацию\n")

# =============================================================================
# 5. ВИЗУАЛИЗАЦИЯ РЕЗУЛЬТАТОВ КЛАСТЕРИЗАЦИИ
# =============================================================================

cat("\n========== ВИЗУАЛИЗАЦИЯ РЕЗУЛЬТАТОВ ==========\n")

# Основной график кластеризации
p_main <- ggplot(clustering_data, aes(x = CW, y = chela, color = Cluster)) +
  geom_point(alpha = 0.7, size = 1.5) +
  scale_color_manual(values = c("blue", "red"), 
                     labels = c("Незрелые", "Зрелые")) +
  labs(title = "Кластеризация краба-стригуна по размеру карапакса и клешни",
       subtitle = paste0("BIC = ", round(robust_gmm$bic, 2),
                        ", Силуэт = ", round(mean(silhouette_score[, 3]), 4)),
       x = "Ширина карапакса (CW, мм)", 
       y = "Размер клешни (chela, мм)",
       color = "Кластер") +
  theme_minimal() +
  theme(legend.position = "bottom")

# График неопределенности классификации (если доступен)
if (!all(is.na(clustering_data$uncertainty))) {
  p_uncertainty <- ggplot(clustering_data, aes(x = CW, y = chela, color = uncertainty)) +
    geom_point(alpha = 0.7, size = 1.5) +
    scale_color_gradient(low = "blue", high = "red", limits = c(0, 0.5),
                         name = "Неопределенность") +
    labs(title = "Неопределенность классификации",
         subtitle = paste0("Средняя неопределенность = ", 
                          round(mean(clustering_data$uncertainty, na.rm = TRUE), 4)),
         x = "Ширина карапакса (CW, мм)", 
         y = "Размер клешни (chela, мм)") +
    theme_minimal() +
    theme(legend.position = "bottom")
}

# Распределение неопределенности
if (!all(is.na(clustering_data$uncertainty))) {
  p_uncertainty_dist <- ggplot(clustering_data, aes(x = uncertainty)) +
    geom_histogram(bins = 30, fill = "lightblue", color = "black") +
    labs(title = "Распределение неопределенности классификации",
         x = "Неопределенность", 
         y = "Частота") +
    theme_minimal()
}

# Распределение ширины карапакса по кластерам
p_cw_dist <- ggplot(clustering_data, aes(x = CW, fill = Cluster)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("blue", "red")) +
  labs(title = "Распределение ширины карапакса по кластерам", 
       x = "Ширина карапакса (CW, мм)", 
       y = "Плотность") +
  theme_minimal() +
  theme(legend.position = "none")

# Распределение размера клешни по кластерам
p_chela_dist <- ggplot(clustering_data, aes(x = chela, fill = Cluster)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("blue", "red")) +
  labs(title = "Распределение размера клешни по кластерам", 
       x = "Размер клешни (chela, мм)", 
       y = "Плотность") +
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
  scale_fill_manual(values = c("blue", "red"), 
                    labels = c("Незрелые", "Зрелые")) +
  labs(title = "Силуэтный анализ кластеризации",
       subtitle = paste0("Средняя ширина силуэта: ", round(mean(silhouette_score[, 3]), 3)),
       x = "Наблюдения", 
       y = "Ширина силуэта",
       fill = "Кластер") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

# Компоновка графиков
if (!all(is.na(clustering_data$uncertainty))) {
  top_plots <- grid.arrange(p_main, p_uncertainty, ncol = 2)
  middle_plots <- grid.arrange(p_uncertainty_dist, p_silhouette, ncol = 2)
  bottom_plots <- grid.arrange(p_cw_dist, p_chela_dist, ncol = 2)
  
  # Объединение всех графиков
  combined_plot <- grid.arrange(top_plots, middle_plots, bottom_plots, 
                               nrow = 3, heights = c(2, 2, 1.5))
} else {
  top_plots <- grid.arrange(p_main, p_silhouette, ncol = 2)
  bottom_plots <- grid.arrange(p_cw_dist, p_chela_dist, ncol = 2)
  
  # Объединение всех графиков
  combined_plot <- grid.arrange(top_plots, bottom_plots, nrow = 2, heights = c(2, 1.5))
}

# Сохранение комплексного графика
ggsave("clustering_analysis.png", combined_plot, width = 12, height = 10, dpi = 300)
cat("\nКомплексный график кластеризации сохранен как 'clustering_analysis.png'\n")

# =============================================================================
# 6. ОПРЕДЕЛЕНИЕ ПАРАМЕТРОВ ЗРЕЛОСТИ L50 и L95
# =============================================================================

cat("\n========== ОПРЕДЕЛЕНИЕ ПАРАМЕТРОВ ЗРЕЛОСТИ ==========\n")

# Определение зрелого кластера (предполагаем, что зрелые особи имеют большие размеры)
cluster_means <- clustering_data %>%
  group_by(Cluster) %>%
  summarise(mean_CW = mean(CW))

mature_cluster <- cluster_means$Cluster[which.max(cluster_means$mean_CW)]
cat("Зрелый кластер определен как кластер", mature_cluster, "\n")

# Создание бинарной переменной зрелости
clustering_data$Mature <- ifelse(clustering_data$Cluster == mature_cluster, 1, 0)

# Логистическая регрессия для определения зависимости зрелости от размера
model <- glm(Mature ~ CW, 
             data = clustering_data, 
             family = binomial)

# Вычисление параметров зрелости L50 и L95
beta0 <- coef(model)[1]  # Интерцепт
beta1 <- coef(model)[2]  # Коэффициент для CW

# L50: длина, при которой 50% особей зрелые
L50 <- -beta0/beta1

# L95: длина, при которой 95% особей зрелые
L95 <- (log(0.95/0.05) - beta0)/beta1

# Округление до целых мм
L50_mat_mm <- round(L50)
L95_mat_mm <- round(L95)

# Визуализация логистической кривой
cw_range <- seq(min(clustering_data$CW), max(clustering_data$CW), length.out = 100)
pred_probs <- predict(model, 
                     newdata = data.frame(CW = cw_range), 
                     type = "response")

p_maturity <- ggplot(clustering_data, aes(x = CW, y = Mature)) +
  geom_point(alpha = 0.3, position = position_jitter(height = 0.02)) +
  geom_line(data = data.frame(CW = cw_range, Probability = pred_probs),
            aes(y = Probability), 
            color = "blue",
            size = 1.5) +
  geom_vline(xintercept = L50, linetype = "dashed", color = "red") +
  geom_vline(xintercept = L95, linetype = "dashed", color = "green") +
  annotate("text", x = L50, y = 0.5, 
           label = paste("L50 =", L50_mat_mm, "мм"), 
           color = "red", hjust = -0.1, size = 4) +
  annotate("text", x = L95, y = 0.95, 
           label = paste("L95 =", L95_mat_mm, "мм"), 
           color = "green", hjust = -0.1, size = 4) +
  labs(title = "Логистическая кривая зрелости краба-стригуна",
       x = "Ширина карапакса (CW, мм)",
       y = "Вероятность зрелости") +
  theme_minimal(base_size = 14)

print(p_maturity)

# Сохранение графика логистической кривой
ggsave("maturity_curve.png", p_maturity, width = 10, height = 6, dpi = 300)
cat("График логистической кривой сохранен как 'maturity_curve.png'\n")

# =============================================================================
# 7. СТАТИСТИЧЕСКАЯ СВОДКА И СОХРАНЕНИЕ РЕЗУЛЬТАТОВ
# =============================================================================

cat("\n========== СТАТИСТИЧЕСКАЯ СВОДКА ==========\n\n")

# Сводка по кластерам
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

cat("СТАТИСТИКА КЛАСТЕРОВ:\n")
print(final_stats)

# Дополнительная статистика по неопределенности (если доступна)
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
}

# Вывод параметров зрелости
cat("\nПАРАМЕТРЫ ЗРЕЛОСТИ КРАБА-СТРИГУНА:\n")
cat("L50 (50% особей зрелые):", L50_mat_mm, "мм\n")
cat("L95 (95% особей зрелые):", L95_mat_mm, "мм\n")

# Сохранение результатов в файл
results <- data.frame(
  Parameter = c("L50", "L95"),
  Value_mm = c(L50_mat_mm, L95_mat_mm),
  Description = c("Длина, при которой 50% особей зрелые", 
                  "Длина, при которой 95% особей зрелые")
)

write.csv(results, "maturity_parameters.csv", row.names = FALSE)
cat("\nРезультаты сохранены в файл 'maturity_parameters.csv'\n")

# Сохранение полных данных с информацией о кластерах
output_data <- DATA_filtered
output_data$Cluster <- robust_clusters
output_data$Mature <- ifelse(output_data$Cluster == mature_cluster, 1, 0)

write.csv(output_data, "classified_crabs.csv", row.names = FALSE)
cat("Классифицированные данные сохранены в файл 'classified_crabs.csv'\n")

# =============================================================================
# 8. ЗАКЛЮЧЕНИЕ
# =============================================================================

cat("\n========== АНАЛИЗ ЗАВЕРШЕН ==========\n")
cat("Параметры зрелости успешно определены и сохранены.\n")
cat("Рекомендуется визуально проверить графики для оценки качества кластеризации.\n")
cat("L50 представляет собой важный параметр для управления промыслом.\n")