# =============================================================================
# ПРАКТИЧЕСКОЕ ЗАНЯТИЕ: СЛУЧАЙНОЕ БЛУЖДАНИЕ В ГИДРОБИОЛОГИИ
# Упрощенная модель с постоянным выпуском частиц запаха
# =============================================================================

# Очистим рабочее пространство
rm(list = ls())

# Загружаем необходимые пакеты
library(ggplot2)
library(dplyr)

# =============================================================================
# ПАРАМЕТРЫ МОДЕЛИ
# =============================================================================

# Количество частиц, выпускаемых за ВЕСЬ период
total_particles <- 2400  # 100 частиц/час × 24 часа

# Время моделирования (часы)
n_steps <- 48  # 24 часа выпуска + 24 часа наблюдения

# Приманка работает 24 часа
release_duration <- 24

# Параметры случайного блуждания
sd_move <- 0.08  # "разброс" частиц

# Параметры течения
current_strength <- 0.05
current_angle <- 45

# Преобразуем угол течения
current_angle_rad <- current_angle * pi / 180
current_x <- current_strength * cos(current_angle_rad)
current_y <- current_strength * sin(current_angle_rad)

# =============================================================================
# ФУНКЦИЯ МОДЕЛИРОВАНИЯ (УПРОЩЕННАЯ)
# =============================================================================

simulate_scent <- function(total_particles, n_steps, release_duration, sd_move, current_x, current_y) {
  # Вычисляем сколько частиц выпускать на каждом шаге
  particles_per_step <- total_particles / release_duration
  
  # Массивы для координат
  x_pos <- matrix(NA, nrow = total_particles, ncol = n_steps)
  y_pos <- matrix(NA, nrow = total_particles, ncol = n_steps)
  
  # Время выпуска для каждой частицы
  release_time <- rep(1:release_duration, each = particles_per_step)
  
  cat("Моделируем распространение запаха...\n")
  
  # Основной цикл по времени
  for (step in 1:n_steps) {
    # Выпускаем новые частицы
    if (step <= release_duration) {
      new_indices <- which(release_time == step)
      x_pos[new_indices, step] <- 0
      y_pos[new_indices, step] <- 0
    }
    
    # Двигаем существующие частицы
    if (step > 1) {
      active <- which(!is.na(x_pos[, step-1]))
      if (length(active) > 0) {
        x_pos[active, step] <- x_pos[active, step-1] + 
          rnorm(length(active), 0, sd_move) + current_x
        y_pos[active, step] <- y_pos[active, step-1] + 
          rnorm(length(active), 0, sd_move) + current_y
      }
    }
  }
  
  return(list(x = x_pos, y = y_pos, release_time = release_time))
}

# =============================================================================
# ЗАПУСК МОДЕЛИРОВАНИЯ
# =============================================================================

cat("=== ЗАПУСК МОДЕЛИ ===\n")
results <- simulate_scent(total_particles, n_steps, release_duration, sd_move, current_x, current_y)

# =============================================================================
# РИСУНОК 1: 4 ФАСЕТКИ (6, 12, 18, 24 ЧАСА)
# =============================================================================

cat("Создаем первый рисунок с 4 фасетками...\n")

# Выбираем моменты времени во время работы приманки
time_points <- c(6, 12, 18, 24)
plot_data <- data.frame()

for (time_point in time_points) {
  # Находим активные частицы в этот момент
  active_indices <- which(!is.na(results$x[, time_point]))
  
  if (length(active_indices) > 0) {
    temp_data <- data.frame(
      x = results$x[active_indices, time_point],
      y = results$y[active_indices, time_point],
      time = paste(time_point, "часов")
    )
    plot_data <- rbind(plot_data, temp_data)
  }
}

# Вычисляем границы для компактного отображения
x_range <- range(plot_data$x, na.rm = TRUE)
y_range <- range(plot_data$y, na.rm = TRUE)

# Первый график: фасетки
p1 <- ggplot(plot_data, aes(x = x, y = y)) +
  geom_point(alpha = 0.4, size = 0.8, color = "blue") +
  
  # Приманка
  annotate("point", x = 0, y = 0, color = "red", size = 3, shape = 17) +
  annotate("text", x = 0, y = y_range[1] + 0.1, label = "ПРИМАНКА", 
           color = "red", size = 2.5, vjust = 1) +
  
  # Течение - располагаем ближе к частицам
  annotate("segment", 
           x = x_range[1] + 0.3, y = y_range[2] - 0.3,
           xend = x_range[1] + 0.3 + current_x * 8, 
           yend = y_range[2] - 0.3 + current_y * 8,
           arrow = arrow(length = unit(0.2, "cm")), 
           color = "darkgreen", size = 0.8) +
  annotate("text", 
           x = x_range[1] + 0.3, y = y_range[2] - 0.5,
           label = "ТЕЧЕНИЕ", color = "darkgreen", size = 2.5) +
  
  # Фасетки
  facet_wrap(~ time, ncol = 2) +
  
  labs(
    title = "Формирование шлейфа запаха от постоянно действующей приманки",
    subtitle = "Приманка работает 24 часа, постоянно выпуская новые частицы",
    x = "Координата X (км)",
    y = "Координата Y (км)"
  ) +
  coord_fixed(xlim = c(x_range[1], x_range[2]), 
              ylim = c(y_range[1], y_range[2])) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
    plot.subtitle = element_text(hjust = 0.5, size = 9),
    strip.background = element_rect(fill = "lightblue"),
    panel.spacing = unit(0.5, "lines")
  )

print(p1)

# =============================================================================
# РИСУНОК 2: ТРАЕКТОРИИ ЧАСТИЦ
# =============================================================================

cat("Создаем второй рисунок с траекториями...\n")

# Выбираем 20 случайных частиц
set.seed(123)
n_trajectories <- 20
sample_particles <- sample(1:total_particles, n_trajectories)

# Подготавливаем данные для траекторий
trajectory_data <- data.frame()
for (particle_id in sample_particles) {
  active_steps <- which(!is.na(results$x[particle_id, ]))
  
  if (length(active_steps) > 0) {
    temp_data <- data.frame(
      x = results$x[particle_id, active_steps],
      y = results$y[particle_id, active_steps],
      step = active_steps,
      particle = as.factor(particle_id)
    )
    trajectory_data <- rbind(trajectory_data, temp_data)
  }
}

# Вычисляем границы для траекторий
traj_x_range <- range(trajectory_data$x, na.rm = TRUE)
traj_y_range <- range(trajectory_data$y, na.rm = TRUE)

# Второй график: траектории
p2 <- ggplot(trajectory_data, aes(x = x, y = y, color = particle)) +
  geom_path(alpha = 0.7, size = 0.5, show.legend = FALSE) +
  geom_point(data = trajectory_data %>% filter(step == 1), 
             color = "red", size = 1.5, shape = 17) +
  geom_point(data = trajectory_data %>% filter(step == max(step)), 
             color = "blue", size = 1.2) +
  
  # Приманка
  annotate("point", x = 0, y = 0, color = "red", size = 3, shape = 17) +
  annotate("text", x = 0, y = traj_y_range[1] + 0.1, 
           label = "ПРИМАНКА", color = "red", size = 2.5, vjust = 1) +
  
  # Течение - располагаем ближе к частицам
  annotate("segment", 
           x = traj_x_range[1] + 0.3, y = traj_y_range[2] - 0.3,
           xend = traj_x_range[1] + 0.3 + current_x * 8, 
           yend = traj_y_range[2] - 0.3 + current_y * 8,
           arrow = arrow(length = unit(0.2, "cm")), 
           color = "darkgreen", size = 0.8) +
  annotate("text", 
           x = traj_x_range[1] + 0.3, y = traj_y_range[2] - 0.5,
           label = "ТЕЧЕНИЕ", color = "darkgreen", size = 2.5) +
  
  labs(
    title = "Траектории движения отдельных частиц запаха",
    subtitle = "Красные треугольники - начало движения, синие точки - конец",
    x = "Координата X (км)",
    y = "Координата Y (км)"
  ) +
  coord_fixed(xlim = c(traj_x_range[1], traj_x_range[2]), 
              ylim = c(traj_y_range[1], traj_y_range[2])) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
    plot.subtitle = element_text(hjust = 0.5, size = 9)
  )

print(p2)

# =============================================================================
# ОБЪЯСНЕНИЕ РЕЗУЛЬТАТОВ
# =============================================================================

cat("\n=== ЧТО МЫ ВИДИМ НА ГРАФИКАХ ===\n")
cat("РИСУНОК 1 (Фасетки):\n")
cat("• Шлейф запаха РАСТЕТ со временем - приманка постоянно работает\n")
cat("• Форма шлейфа определяется ТЕЧЕНИЕМ (зеленая стрелка)\n")
cat("• К 24 часам шлейф достигает максимальной плотности\n\n")

cat("РИСУНОК 2 (Траектории):\n")
cat("• Каждая частица движется СЛУЧАЙНО, но с ДРЕЙФОМ от течения\n")
cat("• Красные треугольники - начало движения у приманки\n")
cat("• Синие точки - конечные положения частиц\n")
cat("• Видно общее направление движения, заданное течением\n\n")

cat("КЛЮЧЕВОЙ ВЫВОД:\n")
cat("Даже СЛУЧАЙНОЕ движение частиц создает УПОРЯДОЧЕННУЮ структуру - шлейф,\n")
cat("когда есть постоянный источник (приманка) и направленное воздействие (течение).\n")