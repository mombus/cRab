# Загрузка необходимых библиотек
library(readxl)
library(mclust)
library(ggplot2)
library(dplyr)
library(factoextra)
library(dbscan)

# Загрузка данных
DATA <- readxl::read_excel("SNOW.xlsx", sheet = "ECO")

# Подготовка данных для кластеризации
# Используем только числовые переменные CW и chela, игнорируя дату
data_clust <- DATA %>% dplyr::select(CW, chela)

# Нормализация данных (важно для корректной работы алгоритмов кластеризации)
data_scaled <- scale(data_clust)
data_scaled <- as.data.frame(data_scaled)

# 1. K-means кластеризация (один из самых популярных алгоритмов кластеризации) [[1]]
set.seed(123)
kmeans_result <- kmeans(data_scaled, centers = 2, nstart = 25)
DATA$kmeans_cluster <- as.factor(kmeans_result$cluster)

# 2. Иерархическая кластеризация (один из наиболее известных методов кластеризации) [[6]]
dist_matrix <- dist(data_scaled)
hc <- hclust(dist_matrix, method = "ward.D2")
hc_clusters <- cutree(hc, k = 2)
DATA$hc_cluster <- as.factor(hc_clusters)

# 3. GMM (Gaussian Mixture Models) - обязательный алгоритм по заданию
# GMM является одним из распространенных алгоритмов кластеризации без учителя [[7]]
gmm <- Mclust(data_scaled, G = 2)
DATA$gmm_cluster <- as.factor(gmm$classification)



# Построение графиков для всех алгоритмов
par(mfrow = c(2, 2))

# K-means
plot(data_clust, col = DATA$kmeans_cluster, pch = 16, main = "K-means Clustering (2 Clusters)",
     xlab = "CW", ylab = "chela")
legend("topright", legend = levels(DATA$kmeans_cluster), col = 1:2, pch = 16)

# Иерархическая кластеризация
plot(data_clust, col = DATA$hc_cluster, pch = 16, main = "Hierarchical Clustering (2 Clusters)",
     xlab = "CW", ylab = "chela")
legend("topright", legend = levels(DATA$hc_cluster), col = 1:2, pch = 16)

# GMM
plot(data_clust, col = DATA$gmm_cluster, pch = 16, main = "GMM Clustering (2 Clusters)",
     xlab = "CW", ylab = "chela")
legend("topright", legend = levels(DATA$gmm_cluster), col = 1:2, pch = 16)


# Альтернативный вариант построения графиков с использованием ggplot2
library(gridExtra)

p1 <- ggplot(DATA, aes(x = CW, y = chela, color = kmeans_cluster)) +
  geom_point() +
  theme_minimal() +
  labs(title = "K-means Clustering", color = "Cluster") +
  theme(plot.title = element_text(hjust = 0.5))

p2 <- ggplot(DATA, aes(x = CW, y = chela, color = hc_cluster)) +
  geom_point() +
  theme_minimal() +
  labs(title = "Hierarchical Clustering", color = "Cluster") +
  theme(plot.title = element_text(hjust = 0.5))

p3 <- ggplot(DATA, aes(x = CW, y = chela, color = gmm_cluster)) +
  geom_point() +
  theme_minimal() +
  labs(title = "GMM Clustering", color = "Cluster") +
  theme(plot.title = element_text(hjust = 0.5))



# Отображение всех графиков в одном окне
grid.arrange(p1, p2, p3,  ncol = 1)