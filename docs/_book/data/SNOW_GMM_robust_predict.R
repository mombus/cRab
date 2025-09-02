# Функция для предсказания кластеров на новых данных
predict_teigen_cluster <- function(new_data, teigen_model, scaled_attrs) {
  # Проверка наличия необходимых столбцов
  if (!all(c("CW", "chela") %in% colnames(new_data))) {
    stop("Новые данные должны содержать столбцы 'CW' и 'chela'")
  }
  
  # Стандартизация новых данных с использованием параметров из обучающей выборки
  new_data_scaled <- scale(new_data, 
                           center = scaled_attrs$center, 
                           scale = scaled_attrs$scale)
  
  # Вычисление вероятностей принадлежности к кластерам
  # Для teigen нам нужно вручную вычислить вероятности, используя параметры модели
  probabilities <- matrix(0, nrow = nrow(new_data), ncol = teigen_model$G)
  
  for (k in 1:teigen_model$G) {
    # Вычисление плотности многомерного t-распределения для каждого кластера
    # Используем параметры модели: средние значения, матрицы ковариации и степени свободы
    cluster_mean <- teigen_model$parameters$mean[, k]
    cluster_cov <- teigen_model$parameters$sigma[, , k]
    cluster_df <- teigen_model$parameters$df[k]
    
    # Вычисление плотности с помощью многомерного t-распределения
    # (упрощенный подход, может потребовать дополнительной проверки)
    for (i in 1:nrow(new_data_scaled)) {
      # Вычисление квадратичной формы Махаланобиса
      diff <- new_data_scaled[i, ] - cluster_mean
      quad_form <- t(diff) %*% solve(cluster_cov) %*% diff
      
      # Вычисление плотности t-распределения
      p <- 2
      density <- gamma((cluster_df + p)/2) / (gamma(cluster_df/2) * (pi * cluster_df)^(p/2) * sqrt(det(cluster_cov))) *
        (1 + quad_form/cluster_df)^(-(cluster_df + p)/2)
      
      probabilities[i, k] <- density * teigen_model$parameters$pig[k]
    }
  }
  
  # Нормализация вероятностей
  probabilities <- probabilities / rowSums(probabilities)
  
  # Определение кластера с максимальной вероятностью
  predicted_clusters <- apply(probabilities, 1, which.max)
  
  # Добавление результатов к данным
  new_data$cluster <- predicted_clusters
  new_data$cluster_prob_1 <- probabilities[, 1]
  new_data$cluster_prob_2 <- probabilities[, 2]
  new_data$uncertainty <- 1 - apply(probabilities, 1, max)
  
  return(new_data)
}

# Пример использования функции предсказания
# Сохраняем параметры стандартизации из обучающей выборки
scaled_attrs <- list(
  center = attr(scaled_data, "scaled:center"),
  scale = attr(scaled_data, "scaled:scale")
)

# Создаем новые данные для предсказания
new_crabs <- data.frame(
  CW = c(39, 23, 110, 65, 120),
  chela = c(5, 4, 25, 10, 30)
)

# Предсказываем кластеры для новых данных
predicted_crabs <- predict_teigen_cluster(new_crabs, robust_gmm, scaled_attrs)
print(predicted_crabs)