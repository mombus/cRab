# ============================================
# МОДЕЛЬ ДРЕЙФА ЛИЧИНОК В БАРЕНЦЕВОМ МОРЕ
# ============================================

rm(list=ls())

library(ncdf4)
library(terra)
library(ggplot2)
library(ggspatial)
library(sf)
library(dplyr)
library(tidyr)
library(viridis)
library(patchwork)

setwd("C:/CURRENT/")

# 1. ЗАГРУЗКА И ОБРАБОТКА ДАННЫХ
# ============================================

# Извлечение данных из NetCDF файлов
extract_current_data <- function() {
  
  # Открытие файлов
  nc_dir_mean <- nc_open("Current direction [mean].nc")
  nc_dir_range <- nc_open("Current direction [range].nc")
  nc_vel_mean <- nc_open("Current velocity [mean].nc")
  nc_vel_range <- nc_open("Current velocity [range].nc")
  
  # Извлечение переменных
  lon <- ncvar_get(nc_dir_mean, "longitude")
  lat <- ncvar_get(nc_dir_mean, "latitude")
  
  dir_mean <- ncvar_get(nc_dir_mean, "swd_mean")
  dir_range <- ncvar_get(nc_dir_range, "swd_range")
  vel_mean <- ncvar_get(nc_vel_mean, "sws_mean")
  vel_range <- ncvar_get(nc_vel_range, "sws_range")
  
  # Закрытие файлов
  nc_close(nc_dir_mean)
  nc_close(nc_dir_range)
  nc_close(nc_vel_mean)
  nc_close(nc_vel_range)
  
  return(list(
    lon = lon,
    lat = lat,
    dir_mean = dir_mean,
    dir_range = dir_range,
    vel_mean = vel_mean,
    vel_range = vel_range
  ))
}

# Загрузка данных
cat("Загрузка данных океанических течений...\n")
current_data <- extract_current_data()

# Проверка размерности данных
cat(sprintf("Размер сетки: %d x %d\n", length(current_data$lon), length(current_data$lat)))

# Определение области Баренцева моря
barents_extent <- c(
  lon_min = 40,
  lon_max = 46,
  lat_min = 67.5,
  lat_max = 69.5
)

# Исправленная функция обрезки данных для Баренцева моря
crop_to_barents <- function(data, extent) {
  lon_idx <- which(data$lon >= extent["lon_min"] & data$lon <= extent["lon_max"])
  lat_idx <- which(data$lat >= extent["lat_min"] & data$lat <= extent["lat_max"])
  
  # Проверка наличия данных в регионе
  if(length(lon_idx) == 0 || length(lat_idx) == 0) {
    stop("Нет данных для указанного региона!")
  }
  
  # Проверка размерности массивов
  cat(sprintf("Размерность dir_mean: %s\n", paste(dim(data$dir_mean), collapse=" x ")))
  cat(sprintf("Индексы долготы: %d:%d из %d\n", min(lon_idx), max(lon_idx), length(data$lon)))
  cat(sprintf("Индексы широты: %d:%d из %d\n", min(lat_idx), max(lat_idx), length(data$lat)))
  
  # Извлечение подмассивов - правильный порядок индексов
  # Если данные 3D (lon, lat, time), берем первый временной слой
  if(length(dim(data$dir_mean)) == 3) {
    dir_mean_subset <- data$dir_mean[lon_idx, lat_idx, 1]
    dir_range_subset <- data$dir_range[lon_idx, lat_idx, 1]
    vel_mean_subset <- data$vel_mean[lon_idx, lat_idx, 1]
    vel_range_subset <- data$vel_range[lon_idx, lat_idx, 1]
  } else {
    # Если данные 2D
    dir_mean_subset <- data$dir_mean[lon_idx, lat_idx]
    dir_range_subset <- data$dir_range[lon_idx, lat_idx]
    vel_mean_subset <- data$vel_mean[lon_idx, lat_idx]
    vel_range_subset <- data$vel_range[lon_idx, lat_idx]
  }
  
  return(list(
    lon = data$lon[lon_idx],
    lat = data$lat[lat_idx],
    dir_mean = dir_mean_subset,
    dir_range = dir_range_subset,
    vel_mean = vel_mean_subset,
    vel_range = vel_range_subset
  ))
}

cat("Обрезка данных для Баренцева моря...\n")
barents_data <- crop_to_barents(current_data, barents_extent)

# Проверка результата
cat(sprintf("Размер региональной сетки: %d x %d\n", 
            length(barents_data$lon), length(barents_data$lat)))

# Замена NA и отрицательных значений
barents_data$vel_mean[is.na(barents_data$vel_mean) | barents_data$vel_mean < 0] <- 0
barents_data$dir_mean[is.na(barents_data$dir_mean)] <- 0
barents_data$vel_range[is.na(barents_data$vel_range) | barents_data$vel_range < 0] <- 0.01
barents_data$dir_range[is.na(barents_data$dir_range) | barents_data$dir_range < 0] <- 10

# Проверка данных
cat("\nСтатистика данных течений:\n")
cat(sprintf("Скорость - мин: %.3f, макс: %.3f, среднее: %.3f м/с\n", 
            min(barents_data$vel_mean, na.rm=TRUE),
            max(barents_data$vel_mean, na.rm=TRUE),
            mean(barents_data$vel_mean, na.rm=TRUE)))
cat(sprintf("Направление - мин: %.1f, макс: %.1f град\n", 
            min(barents_data$dir_mean, na.rm=TRUE),
            max(barents_data$dir_mean, na.rm=TRUE)))

# 2. ФУНКЦИИ ДЛЯ МОДЕЛИ ДРЕЙФА
# ============================================

# Билинейная интерполяция для получения значений течений в произвольной точке
bilinear_interpolation <- function(lon_point, lat_point, data) {
  
  # Найти ближайшие индексы сетки
  lon_idx <- findInterval(lon_point, data$lon)
  lat_idx <- findInterval(lat_point, data$lat)
  
  # Проверка границ
  if(lon_idx < 1 || lon_idx >= length(data$lon) || 
     lat_idx < 1 || lat_idx >= length(data$lat)) {
    return(list(vel = 0.05, dir = runif(1, 0, 360), vel_range = 0.02, dir_range = 30))
  }
  
  # Безопасное извлечение значений с проверкой
  safe_get <- function(matrix, i, j) {
    if(i <= 0 || j <= 0 || i > nrow(matrix) || j > ncol(matrix)) {
      return(NA)
    }
    val <- matrix[i, j]
    if(is.na(val) || is.nan(val) || is.infinite(val)) {
      return(NA)
    }
    return(val)
  }
  
  # Веса для интерполяции
  dx <- (lon_point - data$lon[lon_idx]) / (data$lon[lon_idx + 1] - data$lon[lon_idx])
  dy <- (lat_point - data$lat[lat_idx]) / (data$lat[lat_idx + 1] - data$lat[lat_idx])
  
  dx <- max(0, min(1, dx))
  dy <- max(0, min(1, dy))
  
  # Билинейная интерполяция для каждой переменной
  interpolate_var <- function(var_matrix, default_val = 0) {
    v00 <- safe_get(var_matrix, lon_idx, lat_idx)
    v10 <- safe_get(var_matrix, lon_idx+1, lat_idx)
    v01 <- safe_get(var_matrix, lon_idx, lat_idx+1)
    v11 <- safe_get(var_matrix, lon_idx+1, lat_idx+1)
    
    # Если хотя бы одно значение доступно, используем среднее доступных
    valid_vals <- c(v00, v10, v01, v11)[!is.na(c(v00, v10, v01, v11))]
    
    if(length(valid_vals) == 0) {
      return(default_val)
    } else if(length(valid_vals) < 4) {
      # Если не все значения доступны, используем среднее доступных
      return(mean(valid_vals, na.rm = TRUE))
    }
    
    # Полная билинейная интерполяция
    val <- (1-dx)*(1-dy)*v00 + dx*(1-dy)*v10 + (1-dx)*dy*v01 + dx*dy*v11
    return(val)
  }
  
  vel <- interpolate_var(data$vel_mean, 0.05)
  dir <- interpolate_var(data$dir_mean, runif(1, 0, 360))
  vel_range <- interpolate_var(data$vel_range, 0.02)
  dir_range <- interpolate_var(data$dir_range, 30)
  
  # Убедимся, что значения в допустимых пределах
  vel <- max(0, min(vel, 2))  # Максимум 2 м/с
  vel_range <- max(0.01, min(vel_range, 1))
  dir_range <- max(1, min(dir_range, 180))
  dir <- dir %% 360  # Нормализация направления
  
  return(list(
    vel = vel,
    dir = dir,
    vel_range = vel_range,
    dir_range = dir_range
  ))
}

# Преобразование географических координат
haversine_distance <- function(lon1, lat1, lon2, lat2) {
  R <- 6371000  # Радиус Земли в метрах
  phi1 <- lat1 * pi/180
  phi2 <- lat2 * pi/180
  dphi <- (lat2 - lat1) * pi/180
  dlambda <- (lon2 - lon1) * pi/180
  
  a <- sin(dphi/2)^2 + cos(phi1) * cos(phi2) * sin(dlambda/2)^2
  c <- 2 * atan2(sqrt(a), sqrt(1-a))
  
  return(R * c)
}

# Расчет нового положения частицы
calculate_new_position <- function(lon, lat, velocity, direction, dt_hours, add_stochastic = TRUE) {
  
  # Перевод скорости в км/час
  vel_km_h <- velocity * 3.6
  
  # Расстояние за временной шаг
  distance_km <- vel_km_h * dt_hours
  
  # Направление в радианах
  # Преобразование из океанографического (откуда течет) в математическое
  dir_rad <- (90 - direction) * pi/180
  
  # Добавление стохастической составляющей
  if(add_stochastic) {
    # Случайное отклонение направления (±15°)
    dir_rad <- dir_rad + rnorm(1, 0, pi/12)
    # Случайное изменение скорости (±15%)
    distance_km <- distance_km * (1 + rnorm(1, 0, 0.15))
    distance_km <- max(0, distance_km)  # Не может быть отрицательным
  }
  
  # Расчет изменения координат
  R <- 6371  # Радиус Земли в км
  
  # Изменение широты
  dlat <- (distance_km / R) * cos(dir_rad) * 180/pi
  new_lat <- lat + dlat
  
  # Изменение долготы с учетом схождения меридианов
  if(abs(cos(lat * pi/180)) > 0.001) {  # Избегаем деления на ноль
    dlon <- (distance_km / (R * cos(lat * pi/180))) * sin(dir_rad) * 180/pi
  } else {
    dlon <- 0
  }
  new_lon <- lon + dlon
  
  return(c(new_lon, new_lat))
}

# 3. СИМУЛЯЦИЯ ДРЕЙФА ЛИЧИНОК
# ============================================

simulate_larval_drift <- function(start_points, 
                                 days = 90, 
                                 dt_hours = 6, 
                                 n_particles = 100,
                                 barents_data) {
  
  n_steps <- (days * 24) / dt_hours
  n_starts <- nrow(start_points)
  
  # Инициализация массива траекторий
  trajectories <- array(NA, dim = c(n_starts, n_particles, n_steps + 1, 2))
  
  # Для каждой стартовой точки
  for(s in 1:n_starts) {
    cat(sprintf("\nСимуляция для точки %d из %d (%s)\n", 
                s, n_starts, start_points$name[s]))
    
    pb <- txtProgressBar(min = 0, max = n_particles, style = 3)
    
    # Для каждой частицы
    for(p in 1:n_particles) {
      # Начальное положение
      trajectories[s, p, 1, ] <- as.numeric(start_points[s, c("lon", "lat")])
      
      # Симуляция движения
      for(t in 1:n_steps) {
        current_pos <- trajectories[s, p, t, ]
        
        # Получение параметров течения в текущей точке
        current_params <- bilinear_interpolation(
          current_pos[1], current_pos[2], barents_data
        )
        
        # Добавление случайной составляющей к параметрам
        vel_stoch <- abs(current_params$vel + 
                        rnorm(1, 0, max(current_params$vel_range/2, 0.01)))
        dir_stoch <- (current_params$dir + 
                     rnorm(1, 0, min(current_params$dir_range/2, 30))) %% 360
        
        # Расчет новой позиции
        new_pos <- calculate_new_position(
          current_pos[1], current_pos[2],
          vel_stoch, dir_stoch,
          dt_hours, add_stochastic = TRUE
        )
        
        # Мягкая проверка границ
        if(new_pos[1] < (barents_extent["lon_min"] - 5) || 
           new_pos[1] > (barents_extent["lon_max"] + 5) ||
           new_pos[2] < (barents_extent["lat_min"] - 2) || 
           new_pos[2] > (barents_extent["lat_max"] + 2)) {
          # Если вышли далеко за границы, уменьшаем шаг
          new_pos <- current_pos + 0.5 * (new_pos - current_pos)
        }
        
        trajectories[s, p, t + 1, ] <- new_pos
      }
      
      setTxtProgressBar(pb, p)
    }
    close(pb)
  }
  
  return(trajectories)
}

# 4. ЗАПУСК СИМУЛЯЦИИ
# ============================================

# Начальные точки выпуска личинок
start_points <- data.frame(
  lon = c(43.754, 43.525, 43.200, 43.066, 43.544, 44.002,44.442, 44.901),
  lat = c(68.251, 68.386, 68.512, 68.710,68.715,68.658,68.632,68.626),
  name = c("Точка 1 ", "Точка 2 ", "Точка 3 ", "Точка 4 ", "Точка 5 ", "Точка 6", "Точка 7 ", "Точка 8 ")
)


cat("\n=== НАЧИНАЕМ СИМУЛЯЦИЮ ДРЕЙФА ЛИЧИНОК ===\n")
cat(sprintf("Количество частиц на точку: %d\n", 100))
cat(sprintf("Период моделирования: %d дней\n", 90))
cat(sprintf("Временной шаг: %d часов\n", 6))

set.seed(42)  # Для воспроизводимости

trajectories <- simulate_larval_drift(
  start_points = start_points,
  days = 90,
  dt_hours = 6,
  n_particles = 100,
  barents_data = barents_data
)

cat("\nСимуляция завершена!\n")


# 5. ОБРАБОТКА РЕЗУЛЬТАТОВ
# ============================================

# Преобразование траекторий в data frame для визуализации
prepare_trajectory_data <- function(trajectories, start_points) {
  
  n_starts <- dim(trajectories)[1]
  n_particles <- dim(trajectories)[2]
  n_steps <- dim(trajectories)[3]
  dt_hours <- 6
  
  cat("Преобразование траекторий в формат data frame...\n")
  
  # Создание data frame со всеми траекториями
  traj_list <- list()
  
  for(s in 1:n_starts) {
    for(p in 1:n_particles) {
      # Извлечение траектории одной частицы
      particle_traj <- data.frame(
        start_point = start_points$name[s],
        particle = p,
        step = 1:n_steps,
        day = (0:(n_steps-1)) * dt_hours / 24,
        lon = trajectories[s, p, , 1],
        lat = trajectories[s, p, , 2]
      )
      
      # Удаление NA значений
      particle_traj <- particle_traj[!is.na(particle_traj$lon), ]
      
      if(nrow(particle_traj) > 0) {
        traj_list[[length(traj_list) + 1]] <- particle_traj
      }
    }
  }
  
  traj_df <- do.call(rbind, traj_list)
  return(traj_df)
}

traj_df <- prepare_trajectory_data(trajectories, start_points)
cat(sprintf("✓ Обработано %d точек траекторий\n\n", nrow(traj_df)))

# Расчет плотности распределения личинок
calculate_density <- function(traj_df, target_day, grid_size = 0.5) {
  
  # Фильтрация данных для конкретного дня
  day_data <- traj_df %>%
    filter(abs(day - target_day) < 0.1)
  
  if(nrow(day_data) == 0) {
    return(NULL)
  }
  
  # Создание сетки
  lon_range <- range(day_data$lon, na.rm = TRUE)
  lat_range <- range(day_data$lat, na.rm = TRUE)
  
  # Добавляем буфер к границам
  lon_breaks <- seq(lon_range[1] - grid_size, 
                    lon_range[2] + grid_size, 
                    by = grid_size)
  lat_breaks <- seq(lat_range[1] - grid_size, 
                    lat_range[2] + grid_size, 
                    by = grid_size)
  
  # Создание сетки для подсчета
  density_grid <- expand.grid(
    lon = lon_breaks[-length(lon_breaks)] + grid_size/2,
    lat = lat_breaks[-length(lat_breaks)] + grid_size/2
  )
  
  # Подсчет частиц в каждой ячейке сетки
  density_grid$count <- 0
  
  for(i in 1:nrow(density_grid)) {
    density_grid$count[i] <- sum(
      day_data$lon >= (density_grid$lon[i] - grid_size/2) &
      day_data$lon < (density_grid$lon[i] + grid_size/2) &
      day_data$lat >= (density_grid$lat[i] - grid_size/2) &
      day_data$lat < (density_grid$lat[i] + grid_size/2)
    )
  }
  
  # Расчет вероятности
  total_particles <- length(unique(paste(day_data$start_point, day_data$particle)))
  density_grid$probability <- density_grid$count / total_particles
  density_grid$day <- target_day
  
  return(density_grid)
}

# 6. СТАТИСТИЧЕСКИЙ АНАЛИЗ
# ============================================

cat("=== СТАТИСТИЧЕСКИЙ АНАЛИЗ ===\n")

# Анализ дисперсии и распространения
analyze_dispersion <- function(trajectories, start_points) {
  
  cat("Расчет статистических показателей дисперсии...\n")
  
  n_starts <- dim(trajectories)[1]
  n_particles <- dim(trajectories)[2]
  n_steps <- dim(trajectories)[3]
  dt_hours <- 6
  
  results <- data.frame()
  
  for(s in 1:n_starts) {
    # Анализ для каждого дня
    for(t in seq(1, n_steps, by = 4)) {  # Каждый день (4 шага по 6 часов)
      day <- (t - 1) * dt_hours / 24
      
      # Извлечение позиций всех частиц в данный момент времени
      positions <- trajectories[s, , t, ]
      valid_positions <- positions[!is.na(positions[,1]), , drop = FALSE]
      
      if(nrow(valid_positions) > 5) {  # Минимум 5 частиц для статистики
        # Расчет центроида
        centroid_lon <- mean(valid_positions[,1])
        centroid_lat <- mean(valid_positions[,2])
        
        # Расчет расстояний от начальной точки
        distances <- sapply(1:nrow(valid_positions), function(i) {
          haversine_distance(
            start_points$lon[s], start_points$lat[s],
            valid_positions[i,1], valid_positions[i,2]
          ) / 1000  # в км
        })
        
        # Расчет стандартного отклонения позиций
        std_lon <- sd(valid_positions[,1])
        std_lat <- sd(valid_positions[,2])
        
        # Расчет эллипса рассеивания (приблизительная площадь)
        mean_lat_rad <- centroid_lat * pi/180
        km_per_deg_lon <- 111.32 * cos(mean_lat_rad)
        km_per_deg_lat <- 111.32
        
        area_km2 <- pi * (std_lon * km_per_deg_lon) * (std_lat * km_per_deg_lat)
        
        results <- rbind(results, data.frame(
          start_point = start_points$name[s],
          day = day,
          n_particles = nrow(valid_positions),
          mean_distance_km = mean(distances),
          max_distance_km = max(distances),
          min_distance_km = min(distances),
          std_distance_km = sd(distances),
          centroid_lon = centroid_lon,
          centroid_lat = centroid_lat,
          dispersion_lon = std_lon,
          dispersion_lat = std_lat,
          area_km2 = area_km2
        ))
      }
    }
  }
  
  return(results)
}

dispersion_stats <- analyze_dispersion(trajectories, start_points)

# Вывод таблицы со статистикой для ключевых дней
cat("\n📊 Статистика распространения для ключевых дней:\n")
key_days_stats <- dispersion_stats %>%
  filter(day %in% c(1, 7, 30, 60, 90)) %>%
  select(start_point, day, mean_distance_km, max_distance_km, area_km2) %>%
  mutate(mean_distance_km = round(mean_distance_km, 1),
         max_distance_km = round(max_distance_km, 1),
         area_km2 = round(area_km2, 0))
print(key_days_stats)

# 7. ПОШАГОВАЯ ВИЗУАЛИЗАЦИЯ РЕЗУЛЬТАТОВ
# ============================================

cat("\n=== НАЧАЛО ПОШАГОВОЙ ВИЗУАЛИЗАЦИИ ===\n")
cat("Нажмите Enter после просмотра каждого графика...\n\n")

# Загрузка карты мира для фона
library(maps)
library(mapdata)

# Создание базовой карты Баренцева моря
create_base_map <- function() {
  
  # Получение данных карты
  world_map <- map_data("world")
  
  p <- ggplot() +
    geom_polygon(data = world_map, 
                 aes(x = long, y = lat, group = group),
                 fill = "lightgray", color = "darkgray", size = 0.2) +
    coord_fixed(xlim = c(barents_extent["lon_min"], barents_extent["lon_max"]),
                ylim = c(barents_extent["lat_min"], barents_extent["lat_max"]),
                ratio = 1/cos(mean(c(barents_extent["lat_min"], 
                                    barents_extent["lat_max"])) * pi/180)) +
    theme_minimal() +
    theme(
      panel.grid.major = element_line(color = "gray90", size = 0.5),
      panel.grid.minor = element_line(color = "gray95", size = 0.25),
      panel.background = element_rect(fill = "#E6F3FF"),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5),
      axis.title = element_text(size = 11),
      legend.position = "right"
    ) +
    labs(x = "Долгота (°в.д.)", y = "Широта (°с.ш.)")
  
  return(p)
}

# ========== ГРАФИК 1: Траектории День 1 ==========
cat("📍 График 1: Траектории дрейфа - День 1\n")

day_1_data <- traj_df %>%
  filter(day <= 1) %>%
  group_by(start_point, particle) %>%
  arrange(day)

end_points_1 <- day_1_data %>%
  group_by(start_point, particle) %>%
  slice_tail(n = 1)

mean_dist_1 <- dispersion_stats %>%
  filter(abs(day - 1) < 0.5) %>%
  summarise(dist = mean(mean_distance_km, na.rm = TRUE)) %>%
  pull(dist)

traj_plot_1 <- create_base_map() +
  geom_path(data = day_1_data,
            aes(x = lon, y = lat, 
                group = interaction(start_point, particle),
                color = start_point),
            alpha = 0.3, size = 0.4) +
  geom_point(data = end_points_1,
             aes(x = lon, y = lat, color = start_point),
             alpha = 0.5, size = 1) +
  geom_point(data = start_points,
             aes(x = lon, y = lat),
             color = "red", size = 2, shape = 16) +
  scale_color_viridis_d(name = "Источник", option = "turbo") +
  labs(title = "Траектории дрейфа личинок - День 1",
       subtitle = sprintf("Средняя дистанция: %.1f км", mean_dist_1))

print(traj_plot_1)
readline(prompt="Нажмите [Enter] для продолжения...")

# ========== ГРАФИК 2: Траектории День 7 ==========
cat("\n📍 График 2: Траектории дрейфа - День 7\n")

day_7_data <- traj_df %>%
  filter(day <= 7) %>%
  group_by(start_point, particle) %>%
  arrange(day)

end_points_7 <- day_7_data %>%
  group_by(start_point, particle) %>%
  slice_tail(n = 1)

mean_dist_7 <- dispersion_stats %>%
  filter(abs(day - 7) < 0.5) %>%
  summarise(dist = mean(mean_distance_km, na.rm = TRUE)) %>%
  pull(dist)

traj_plot_7 <- create_base_map() +
  geom_path(data = day_7_data,
            aes(x = lon, y = lat, 
                group = interaction(start_point, particle),
                color = start_point),
            alpha = 0.2, size = 0.3) +
  geom_point(data = end_points_7,
             aes(x = lon, y = lat, color = start_point),
             alpha = 0.4, size = 0.8) +
  geom_point(data = start_points,
             aes(x = lon, y = lat),
             color = "red", size = 2, shape = 16) +
  scale_color_viridis_d(name = "Источник", option = "turbo") +
  labs(title = "Траектории дрейфа личинок - День 7",
       subtitle = sprintf("Средняя дистанция: %.1f км", mean_dist_7))

print(traj_plot_7)
readline(prompt="Нажмите [Enter] для продолжения...")

# ========== ГРАФИК 3: Траектории День 30 ==========
cat("\n📍 График 3: Траектории дрейфа - День 30\n")

day_30_data <- traj_df %>%
  filter(day <= 30) %>%
  group_by(start_point, particle) %>%
  arrange(day)

end_points_30 <- day_30_data %>%
  group_by(start_point, particle) %>%
  slice_tail(n = 1)

mean_dist_30 <- dispersion_stats %>%
  filter(abs(day - 30) < 0.5) %>%
  summarise(dist = mean(mean_distance_km, na.rm = TRUE)) %>%
  pull(dist)

traj_plot_30 <- create_base_map() +
  geom_path(data = day_30_data,
            aes(x = lon, y = lat, 
                group = interaction(start_point, particle),
                color = start_point),
            alpha = 0.15, size = 0.3) +
  geom_point(data = end_points_30,
             aes(x = lon, y = lat, color = start_point),
             alpha = 0.3, size = 0.6) +
  geom_point(data = start_points,
             aes(x = lon, y = lat),
             color = "red", size = 2, shape = 16) +
  scale_color_viridis_d(name = "Источник", option = "turbo") +
  labs(title = "Траектории дрейфа личинок - День 30",
       subtitle = sprintf("Средняя дистанция: %.1f км", mean_dist_30))

print(traj_plot_30)
readline(prompt="Нажмите [Enter] для продолжения...")

# ========== ГРАФИК 4: Траектории День 90 ==========
cat("\n📍 График 4: Траектории дрейфа - День 90\n")

day_90_data <- traj_df %>%
  filter(day <= 90) %>%
  group_by(start_point, particle) %>%
  arrange(day)

end_points_90 <- day_90_data %>%
  group_by(start_point, particle) %>%
  slice_tail(n = 1)

mean_dist_90 <- dispersion_stats %>%
  filter(abs(day - 90) < 0.5) %>%
  summarise(dist = mean(mean_distance_km, na.rm = TRUE)) %>%
  pull(dist)

traj_plot_90 <- create_base_map() +
  geom_path(data = day_90_data,
            aes(x = lon, y = lat, 
                group = interaction(start_point, particle),
                color = start_point),
            alpha = 0.1, size = 0.2) +
  geom_point(data = end_points_90,
             aes(x = lon, y = lat, color = start_point),
             alpha = 0.3, size = 0.5) +
  geom_point(data = start_points,
             aes(x = lon, y = lat),
             color = "red", size = 2, shape = 16) +
  scale_color_viridis_d(name = "Источник", option = "turbo") +
  labs(title = "Траектории дрейфа личинок - День 90",
       subtitle = sprintf("Средняя дистанция: %.1f км", mean_dist_90))

print(traj_plot_90)
readline(prompt="Нажмите [Enter] для продолжения...")

# ========== ГРАФИК 5: Плотность День 7 ==========
cat("\n🗺️  График 5: Плотность распределения - День 7\n")

density_7 <- calculate_density(traj_df, 7, grid_size = 0.125)

if(!is.null(density_7)) {
  density_7_filtered <- density_7 %>% filter(count > 0)
  
  density_plot_7 <- create_base_map() +
    geom_tile(data = density_7_filtered,
              aes(x = lon, y = lat, fill = probability),
              alpha = 0.8) +
    geom_point(data = start_points,
               aes(x = lon, y = lat),
               color = "red", size = 2, shape = 16) +
    scale_fill_viridis(name = "Вероятность",
                       option = "plasma",
                       trans = "sqrt") +
    labs(title = "Плотность распределения личинок - День 7",
         subtitle = "Вероятность нахождения в ячейке 0.25° × 0.25°")
  
  print(density_plot_7)
  readline(prompt="Нажмите [Enter] для продолжения...")
}

# ========== ГРАФИК 6: Плотность День 30 ==========
cat("\n🗺️  График 6: Плотность распределения - День 30\n")

density_30 <- calculate_density(traj_df, 30, grid_size = 0.125)

if(!is.null(density_30)) {
  density_30_filtered <- density_30 %>% filter(count > 0)
  
  density_plot_30 <- create_base_map() +
    geom_tile(data = density_30_filtered,
              aes(x = lon, y = lat, fill = probability),
              alpha = 0.8) +
    stat_contour(data = density_30_filtered,
                 aes(x = lon, y = lat, z = probability),
                 color = "white", alpha = 0.5, size = 0.3,
                 breaks = c(0.001, 0.005, 0.01)) +
    geom_point(data = start_points,
               aes(x = lon, y = lat),
               color = "red", size = 2, shape = 16) +
    scale_fill_viridis(name = "Вероятность",
                       option = "plasma",
                       trans = "sqrt") +
    labs(title = "Плотность распределения личинок - День 30",
         subtitle = "Вероятность нахождения в ячейке 0.25° × 0.25°")
  
  print(density_plot_30)
  readline(prompt="Нажмите [Enter] для продолжения...")
}

# ========== ГРАФИК 7: Плотность День 90 ==========
cat("\n🗺️  График 7: Плотность распределения - День 90\n")

density_90 <- calculate_density(traj_df, 90, grid_size = 0.125)

if(!is.null(density_90)) {
  density_90_filtered <- density_90 %>% filter(count > 0)
  
  density_plot_90 <- create_base_map() +
    geom_tile(data = density_90_filtered,
              aes(x = lon, y = lat, fill = probability),
              alpha = 0.8) +
    stat_contour(data = density_90_filtered,
                 aes(x = lon, y = lat, z = probability),
                 color = "white", alpha = 0.5, size = 0.3,
                 breaks = c(0.0005, 0.001, 0.005, 0.01)) +
    geom_point(data = start_points,
               aes(x = lon, y = lat),
               color = "red", size = 2, shape = 16) +
    scale_fill_viridis(name = "Вероятность",
                       option = "plasma",
                       trans = "sqrt") +
    labs(title = "Плотность распределения личинок - День 90",
         subtitle = "Вероятность нахождения в ячейке 0.25° × 0.25°")
  
  print(density_plot_90)
  readline(prompt="Нажмите [Enter] для продолжения...")
}

# ========== ГРАФИК 8: Динамика расстояния ==========
cat("\n📈 График 8: Динамика средней дистанции дрейфа\n")

dist_plot <- ggplot(dispersion_stats, 
                    aes(x = day, y = mean_distance_km, color = start_point)) +
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = mean_distance_km - std_distance_km,
                  ymax = mean_distance_km + std_distance_km,
                  fill = start_point),
              alpha = 0.2) +
  geom_point(size = 2) +
  scale_color_viridis_d(name = "Начальная точка", option = "turbo") +
  scale_fill_viridis_d(name = "Начальная точка", option = "turbo") +
  labs(title = "Средняя дистанция дрейфа личинок",
       subtitle = "Полоса показывает ± стандартное отклонение",
       x = "Время (дни)",
       y = "Расстояние от начальной точки (км)") +
  theme_minimal() +
  theme(legend.position = "bottom")

print(dist_plot)
readline(prompt="Нажмите [Enter] для продолжения...")

# ========== ГРАФИК 9: Динамика площади ==========
cat("\n📈 График 9: Динамика площади распространения\n")

area_plot <- ggplot(dispersion_stats, 
                    aes(x = day, y = area_km2, color = start_point)) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  scale_color_viridis_d(name = "Начальная точка", option = "turbo") +
  scale_y_continuous(trans = "log10",
                    labels = scales::comma) +
  labs(title = "Площадь распространения личинок",
       subtitle = "Логарифмическая шкала",
       x = "Время (дни)",
       y = "Площадь распространения (км²)") +
  theme_minimal() +
  theme(legend.position = "bottom")

print(area_plot)
readline(prompt="Нажмите [Enter] для продолжения...")

# ========== ГРАФИК 10: Итоговое распределение ==========
cat("\n🎯 График 10: Итоговое распределение на день 90\n")

# Данные для финального дня
final_data <- traj_df %>%
  filter(abs(day - 90) < 0.1)

# Расчет статистики для каждой начальной точки
confidence_ellipses <- final_data %>%
  group_by(start_point) %>%
  summarise(
    mean_lon = mean(lon),
    mean_lat = mean(lat),
    sd_lon = sd(lon),
    sd_lat = sd(lat),
    n = n(),
    .groups = 'drop'
  )

cat("\n📊 Центроиды распределений на день 90:\n")
print(confidence_ellipses %>% 
      select(start_point, mean_lon, mean_lat, sd_lon, sd_lat) %>%
      mutate_if(is.numeric, round, 3))

# Создание эллипсов доверительных интервалов
create_ellipse <- function(center_lon, center_lat, sd_lon, sd_lat, n_points = 100) {
  angles <- seq(0, 2*pi, length.out = n_points)
  ellipse_x <- 2 * sd_lon * cos(angles)
  ellipse_y <- 2 * sd_lat * sin(angles)
  
  ellipse_points <- data.frame(
    lon = center_lon + ellipse_x,
    lat = center_lat + ellipse_y
  )
  return(ellipse_points)
}

ellipse_polygons <- confidence_ellipses %>%
  group_by(start_point) %>%
  do(create_ellipse(.$mean_lon, .$mean_lat, .$sd_lon, .$sd_lat)) %>%
  ungroup()

final_plot <- create_base_map() +
  {if(!is.null(density_90) && nrow(density_90_filtered) > 0) 
    geom_tile(data = density_90_filtered,
              aes(x = lon, y = lat, fill = probability),
              alpha = 0.6)} +
  geom_path(data = ellipse_polygons,
            aes(x = lon, y = lat, group = start_point, color = start_point),
            size = 1.5, linetype = "dashed") +
  geom_point(data = confidence_ellipses,
             aes(x = mean_lon, y = mean_lat, color = start_point),
             size = 4, shape = 19) +
  geom_point(data = start_points,
             aes(x = lon, y = lat),
             color = "red", size = 2, shape = 16) +
  scale_fill_viridis(name = "Вероятность",
                     option = "plasma",
                     trans = "sqrt") +
  scale_color_viridis_d(name = "Начальная точка", option = "turbo") +
  labs(title = "Итоговое распределение личинок - День 90",
       subtitle = "Пунктир: 95% доверительные интервалы; Точки: центроиды",
       x = "Долгота (°в.д.)",
       y = "Широта (°с.ш.)")

print(final_plot)
readline(prompt="Нажмите [Enter] для продолжения...")

# ========== АНАЛИЗ ЦЕЛЕВЫХ ОБЛАСТЕЙ ==========
cat("\n🎯 Анализ достижения целевых областей\n")

# Определение целевых областей
target_areas <- data.frame(
  lon = c(45, 35, 25, 50),
  lat = c(72, 75, 73, 70),
  name = c("Новоземельское мелководье", "Центральная впадина", 
           "Западный шельф", "Восточный район"),
  radius_km = c(50, 75, 50, 60)
)

cat("\n📍 Целевые области:\n")
print(target_areas)

# Расчет вероятностей достижения
arrival_prob <- data.frame()

for(i in 1:nrow(target_areas)) {
  target <- target_areas[i,]
  
  for(sp in unique(traj_df$start_point)) {
    point_data <- traj_df %>%
      filter(start_point == sp)
    
    arrivals <- point_data %>%
      group_by(particle) %>%
      mutate(
        dist_to_target = sqrt((lon - target$lon)^2 + (lat - target$lat)^2) * 111,
        arrived = dist_to_target <= target$radius_km
      ) %>%
      filter(arrived) %>%
      summarise(
        first_arrival_day = min(day),
        .groups = 'drop'
      )
    
    n_total <- length(unique(point_data$particle))
    n_arrived <- nrow(arrivals)
    
    arrival_prob <- rbind(arrival_prob, data.frame(
      start_point = sp,
      target_area = target$name,
      probability = n_arrived / n_total,
      mean_arrival_day = ifelse(n_arrived > 0, mean(arrivals$first_arrival_day), NA)
    ))
  }
}

cat("\n📊 Вероятности достижения целевых областей:\n")
arrival_prob_formatted <- arrival_prob %>%
  mutate(probability = paste0(round(probability * 100, 1), "%"),
         mean_arrival_day = round(mean_arrival_day, 1))
print(arrival_prob_formatted)

readline(prompt="Нажмите [Enter] для продолжения...")

# ========== ГРАФИК 11: Вероятности достижения ==========
cat("\n📊 График 11: Вероятности достижения целевых областей\n")

arrival_plot <- ggplot(arrival_prob, 
                       aes(x = target_area, y = probability, fill = start_point)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_viridis_d(name = "Начальная точка", option = "turbo") +
  scale_y_continuous(labels = scales::percent, limits = c(0, max(arrival_prob$probability) * 1.1)) +
  labs(title = "Вероятность достижения целевых областей",
       x = "Целевая область",
       y = "Вероятность") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")

print(arrival_plot)
readline(prompt="Нажмите [Enter] для продолжения...")

# ========== ИТОГОВАЯ СТАТИСТИКА ==========
cat("\n=== ИТОГОВАЯ СТАТИСТИКА МОДЕЛИРОВАНИЯ ===\n")

cat(sprintf("\n📊 Общие параметры:\n"))
cat(sprintf("  • Период моделирования: %d дней\n", max(traj_df$day)))
cat(sprintf("  • Количество начальных точек: %d\n", length(unique(traj_df$start_point))))
cat(sprintf("  • Количество частиц на точку: %d\n", length(unique(traj_df$particle))))
cat(sprintf("  • Общее количество частиц: %d\n", 
            length(unique(paste(traj_df$start_point, traj_df$particle)))))

# Финальная статистика
cat("\n📈 Статистика распространения на день 90:\n")
final_stats <- dispersion_stats %>%
  filter(day == max(day)) %>%
  select(start_point, mean_distance_km, max_distance_km, area_km2) %>%
  mutate(mean_distance_km = round(mean_distance_km, 1),
         max_distance_km = round(max_distance_km, 1),
         area_km2 = round(area_km2, 0))
print(final_stats)

# Средние скорости дрейфа
cat("\n🚢 Средние скорости дрейфа:\n")
mean_speeds <- dispersion_stats %>%
  group_by(start_point) %>%
  summarise(
    mean_speed_km_day = round(mean(mean_distance_km / day, na.rm = TRUE), 2),
    max_speed_km_day = round(max(mean_distance_km / day, na.rm = TRUE), 2),
    .groups = 'drop'
  )
print(mean_speeds)

# ========== СРАВНИТЕЛЬНЫЙ АНАЛИЗ ==========
cat("\n📊 Сравнительный анализ траекторий:\n")

# Анализ направления движения
direction_analysis <- traj_df %>%
  filter(day == 90) %>%
  group_by(start_point) %>%
  summarise(
    mean_final_lon = mean(lon),
    mean_final_lat = mean(lat),
    .groups = 'drop'
  ) %>%
  left_join(start_points %>% 
            select(name, start_lon = lon, start_lat = lat), 
            by = c("start_point" = "name")) %>%
  mutate(
    direction_deg = atan2(mean_final_lon - start_lon, 
                         mean_final_lat - start_lat) * 180 / pi,
    direction_deg = ifelse(direction_deg < 0, direction_deg + 360, direction_deg),
    direction_text = case_when(
      direction_deg >= 337.5 | direction_deg < 22.5 ~ "С",
      direction_deg >= 22.5 & direction_deg < 67.5 ~ "СВ",
      direction_deg >= 67.5 & direction_deg < 112.5 ~ "В",
      direction_deg >= 112.5 & direction_deg < 157.5 ~ "ЮВ",
      direction_deg >= 157.5 & direction_deg < 202.5 ~ "Ю",
      direction_deg >= 202.5 & direction_deg < 247.5 ~ "ЮЗ",
      direction_deg >= 247.5 & direction_deg < 292.5 ~ "З",
      direction_deg >= 292.5 & direction_deg < 337.5 ~ "СЗ"
    )
  )

cat("\nНаправления дрейфа:\n")
print(direction_analysis %>% 
      select(start_point, direction_deg, direction_text) %>%
      mutate(direction_deg = round(direction_deg, 1)))

readline(prompt="Нажмите [Enter] для продолжения...")

# ========== ГРАФИК 12: Сравнение траекторий ==========
cat("\n📊 График 12: Сравнение средних траекторий\n")

# Расчет средних траекторий
mean_trajectories <- traj_df %>%
  filter(day %in% seq(0, 90, by = 3)) %>%
  group_by(start_point, day) %>%
  summarise(
    mean_lon = mean(lon),
    mean_lat = mean(lat),
    sd_lon = sd(lon),
    sd_lat = sd(lat),
    .groups = 'drop'
  )

comparison_plot <- create_base_map() +
  # Средние траектории
  geom_path(data = mean_trajectories,
            aes(x = mean_lon, y = mean_lat, color = start_point),
            size = 2, alpha = 0.8) +
  # Полосы неопределенности
  geom_errorbar(data = mean_trajectories %>% filter(day %% 15 == 0),
                aes(x = mean_lon, 
                    ymin = mean_lat - sd_lat, 
                    ymax = mean_lat + sd_lat,
                    color = start_point),
                width = 0.3, alpha = 0.5) +
  geom_errorbarh(data = mean_trajectories %>% filter(day %% 15 == 0),
                 aes(y = mean_lat,
                     xmin = mean_lon - sd_lon,
                     xmax = mean_lon + sd_lon,
                     color = start_point),
                 height = 0.3, alpha = 0.5) +
  # Точки через каждые 15 дней
  geom_point(data = mean_trajectories %>% filter(day %% 15 == 0),
             aes(x = mean_lon, y = mean_lat, color = start_point),
             size = 3) +
  # Метки дней
  geom_text(data = mean_trajectories %>% filter(day %% 30 == 0),
            aes(x = mean_lon, y = mean_lat + 0.2, 
                label = paste0("День ", day)),
            size = 2.5) +
  # Начальные точки
  geom_point(data = start_points,
             aes(x = lon, y = lat),
             color = "red", size = 5, shape = 17) +
  scale_color_viridis_d(name = "Начальная точка", option = "turbo") +
  labs(title = "Сравнение средних траекторий дрейфа",
       subtitle = "Показаны средние пути с интервалами неопределенности",
       x = "Долгота (°в.д.)",
       y = "Широта (°с.ш.)")

print(comparison_plot)
readline(prompt="Нажмите [Enter] для продолжения...")

# ========== ГРАФИК 13: Розы дрейфа ==========
cat("\n🌹 График 13: Розы дрейфа для каждой начальной точки\n")

# Подготовка данных для роз дрейфа
drift_roses_data <- traj_df %>%
  filter(day == 90) %>%
  left_join(start_points %>% 
            select(name, start_lon = lon, start_lat = lat), 
            by = c("start_point" = "name")) %>%
  mutate(
    distance_km = sqrt((lon - start_lon)^2 + (lat - start_lat)^2) * 111,
    direction_deg = atan2(lon - start_lon, lat - start_lat) * 180 / pi,
    direction_deg = ifelse(direction_deg < 0, direction_deg + 360, direction_deg),
    direction_bin = cut(direction_deg, 
                       breaks = seq(0, 360, by = 30),
                       labels = FALSE,
                       include.lowest = TRUE)
  )

rose_plot <- ggplot(drift_roses_data, 
                    aes(x = direction_deg, fill = start_point)) +
  geom_histogram(binwidth = 30, boundary = 0, alpha = 0.7) +
  coord_polar(start = 0) +
  scale_x_continuous(breaks = seq(0, 330, by = 30),
                     labels = c("С", "30°", "60°", "В", "120°", "150°",
                               "Ю", "210°", "240°", "З", "300°", "330°")) +
  facet_wrap(~start_point, ncol = 3) +
  scale_fill_viridis_d(name = "Начальная точка", option = "turbo") +
  labs(title = "Розы дрейфа личинок",
       subtitle = "Распределение конечных направлений на день 90",
       x = "Направление",
       y = "Количество частиц") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 8))

print(rose_plot)
readline(prompt="Нажмите [Enter] для продолжения...")

# ========== АНАЛИЗ СВЯЗНОСТИ РАЙОНОВ ==========
cat("\n🔗 Анализ связности между районами\n")

# Создание сетки районов
grid_size <- 2  # градусы
lon_grid <- seq(barents_extent["lon_min"], barents_extent["lon_max"], by = grid_size)
lat_grid <- seq(barents_extent["lat_min"], barents_extent["lat_max"], by = grid_size)

# Определение начального и конечного района для каждой частицы
connectivity <- traj_df %>%
  filter(day %in% c(0, 90)) %>%
  select(start_point, particle, day, lon, lat) %>%
  pivot_wider(names_from = day, values_from = c(lon, lat)) %>%
  mutate(
    start_zone_lon = cut(lon_0, breaks = lon_grid, labels = FALSE, include.lowest = TRUE),
    start_zone_lat = cut(lat_0, breaks = lat_grid, labels = FALSE, include.lowest = TRUE),
    end_zone_lon = cut(lon_90, breaks = lon_grid, labels = FALSE, include.lowest = TRUE),
    end_zone_lat = cut(lat_90, breaks = lat_grid, labels = FALSE, include.lowest = TRUE),
    start_zone = paste0("Z", start_zone_lon, "-", start_zone_lat),
    end_zone = paste0("Z", end_zone_lon, "-", end_zone_lat)
  ) %>%
  filter(!is.na(end_zone_lon) & !is.na(end_zone_lat))

# Матрица связности
connectivity_matrix <- connectivity %>%
  count(start_point, start_zone, end_zone) %>%
  group_by(start_point, start_zone) %>%
  mutate(probability = n / sum(n)) %>%
  ungroup()

cat("\n📊 Топ-10 связей между зонами:\n")
top_connections <- connectivity_matrix %>%
  arrange(desc(n)) %>%
  head(10) %>%
  select(start_point, start_zone, end_zone, n, probability) %>%
  mutate(probability = paste0(round(probability * 100, 1), "%"))
print(top_connections)

readline(prompt="Нажмите [Enter] для продолжения...")

# ========== ВРЕМЕННОЙ АНАЛИЗ ==========
cat("\n⏱️ Временной анализ распространения\n")

# Анализ скорости распространения по периодам
temporal_analysis <- dispersion_stats %>%
  mutate(
    period = case_when(
      day <= 7 ~ "Неделя 1",
      day <= 30 ~ "Недели 2-4",
      day <= 60 ~ "Месяц 2",
      TRUE ~ "Месяц 3"
    )
  ) %>%
  group_by(start_point, period) %>%
  summarise(
    mean_speed = mean(mean_distance_km / day, na.rm = TRUE),
    mean_area_growth = mean(area_km2 / day, na.rm = TRUE),
    .groups = 'drop'
  )

cat("\n📊 Скорость распространения по периодам (км/день):\n")
temporal_table <- temporal_analysis %>%
  select(start_point, period, mean_speed) %>%
  mutate(mean_speed = round(mean_speed, 2)) %>%
  pivot_wider(names_from = period, values_from = mean_speed)
print(temporal_table)

# ========== СОХРАНЕНИЕ РЕЗУЛЬТАТОВ ==========
cat("\n💾 Сохранение результатов\n")

save_choice <- readline(prompt = "Сохранить все результаты в файлы? (y/n): ")

if(tolower(save_choice) == "y") {
  
  cat("\n📁 Сохранение в папку drift_results...\n")
  
  # Создание директории
  dir.create("drift_results", showWarnings = FALSE)
  
  # Список всех графиков для сохранения
  plots_to_save <- list(
    trajectories_day_1 = traj_plot_1,
    trajectories_day_7 = traj_plot_7,
    trajectories_day_30 = traj_plot_30,
    trajectories_day_90 = traj_plot_90,
    density_day_7 = density_plot_7,
    density_day_30 = density_plot_30,
    density_day_90 = density_plot_90,
    dispersion_distance = dist_plot,
    dispersion_area = area_plot,
    final_distribution = final_plot,
    arrival_probabilities = arrival_plot,
    mean_trajectories = comparison_plot,
    drift_roses = rose_plot
  )
  
  # Сохранение графиков
  for(name in names(plots_to_save)) {
    if(exists(deparse(substitute(plots_to_save[[name]])))) {
      ggsave(sprintf("drift_results/%s.png", name),
             plot = plots_to_save[[name]],
             width = 12, height = 10, dpi = 300)
      cat(sprintf("  ✓ Сохранен: %s.png\n", name))
    }
  }
  
  # Сохранение таблиц
  write.csv(dispersion_stats, "drift_results/dispersion_statistics.csv", row.names = FALSE)
  cat("  ✓ Сохранена статистика дисперсии\n")
  
  write.csv(arrival_prob, "drift_results/arrival_probabilities.csv", row.names = FALSE)
  cat("  ✓ Сохранены вероятности достижения\n")
  
  write.csv(connectivity_matrix, "drift_results/connectivity_matrix.csv", row.names = FALSE)
  cat("  ✓ Сохранена матрица связности\n")
  
  write.csv(final_stats, "drift_results/final_statistics.csv", row.names = FALSE)
  cat("  ✓ Сохранена финальная статистика\n")
  
  write.csv(mean_speeds, "drift_results/mean_speeds.csv", row.names = FALSE)
  cat("  ✓ Сохранены средние скорости\n")
  
  cat(sprintf("\n✅ Все результаты сохранены в: %s/drift_results\n", getwd()))
}

# ========== ЗАКЛЮЧЕНИЕ ==========
cat("\n" ,rep("=", 50), "\n", sep="")
cat("🎉 МОДЕЛИРОВАНИЕ ДРЕЙФА ЛИЧИНОК ЗАВЕРШЕНО!\n")
cat(rep("=", 50), "\n\n", sep="")

cat("📋 ОСНОВНЫЕ ВЫВОДЫ:\n\n")

# Вывод 1: Максимальная дистанция
max_dist_point <- final_stats %>%
  arrange(desc(mean_distance_km)) %>%
  slice(1)
cat(sprintf("1. Максимальная средняя дистанция дрейфа (%.1f км) наблюдается\n   для личинок из %s\n\n",
            max_dist_point$mean_distance_km, max_dist_point$start_point))

# Вывод 2: Максимальная площадь распространения
max_area_point <- final_stats %>%
  arrange(desc(area_km2)) %>%
  slice(1)
cat(sprintf("2. Максимальная площадь распространения (%.0f км²)\n   достигнута личинками из %s\n\n",
            max_area_point$area_km2, max_area_point$start_point))

# Вывод 3: Преобладающее направление
main_direction <- direction_analysis %>%
  count(direction_text) %>%
  arrange(desc(n)) %>%
  slice(1)
cat(sprintf("3. Преобладающее направление дрейфа: %s\n\n", main_direction$direction_text))

# Вывод 4: Наиболее вероятная целевая область
most_probable_target <- arrival_prob %>%
  arrange(desc(probability)) %>%
  slice(1)
cat(sprintf("4. Наиболее вероятно достижение области '%s'\n   личинками из %s (вероятность %.1f%%)\n\n",
            most_probable_target$target_area, 
            most_probable_target$start_point,
            most_probable_target$probability * 100))

cat("Спасибо за использование модели дрейфа личинок!\n")
cat(rep("=", 50), "\n", sep="")
