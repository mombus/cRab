# Очистка окружения и установка рабочей директории
rm(list = ls())
setwd("C:/KARTOGRAPH/")

# Загрузка необходимых библиотек
library(rnaturalearth)
library(tidyverse)
library(sf)
library(scatterpie)
library(ggspatial)  # ?? ДОБАВЛЕНА ЗАГРУЗКА ПАКЕТА ДЛЯ МАСШТАБНОЙ ЛИНЕЙКИ

####### ЗАГРУЗКА ДАННЫХ И ПОДГОТОВКА ################

# Чтение данных из нового файла
DATA <- readxl::read_excel("KARTOGRAPHICpie.xlsx", sheet = "pie")

# Получение границ России
russia <- ne_countries(scale = 10, country = "Russia")

# Установка границ отображаемой области (в градусах)
xmin <- 40
xmax <- 52
ymin <- 70.5
ymax <- 71

# Преобразование в проекцию EPSG:3575 (метры)
crs_wgs84 <- st_crs(4326)   # исходная система — градусы
crs_laea <- st_crs(3575)    # целевая — метры, для Арктики

# Преобразуем Россию
russia_proj <- st_transform(russia, crs = crs_laea)

# Преобразуем точки из DATA в sf и переводим в проекцию
DATA_sf <- st_as_sf(DATA, coords = c("X", "Y"), crs = crs_wgs84)
DATA_proj <- st_transform(DATA_sf, crs = crs_laea)

# Извлекаем координаты в метрах
DATA_coords <- st_coordinates(DATA_proj) %>% as_tibble() %>% setNames(c("X", "Y"))

# Рассчитываем общую численность самцов
DATA_attrs <- st_drop_geometry(DATA_sf) %>%
  mutate(
    total = `Самцы пром.` + `Самцы пререкруты I` + `Самцы пререкруты II` + `Самцы молодь`
  )

# Определяем диапазон радиусов в градусах (можно настроить)
min_radius_deg <- 0.1   # для точек с малым total
max_radius_deg <- 0.5   # для точек с максимальным total

# Масштабируем радиусы относительно total (линейно)
max_total <- max(DATA_attrs$total, na.rm = TRUE)
min_total <- min(DATA_attrs$total, na.rm = TRUE)

# Нормализуем total > радиус в градусах
DATA_attrs <- DATA_attrs %>%
  mutate(
    radius_deg = ifelse(total == 0, min_radius_deg,
                        scales::rescale(total, 
                                        to = c(min_radius_deg, max_radius_deg), 
                                        from = c(min_total, max_total)))
  )

# Переводим радиус из градусов в метры (на широте ~70°)
# Используем вспомогательный отрезок для точного перевода
radius_ref_deg <- 1.0  # эталонный градус для перевода
radius_ref_line <- st_sfc(
  st_linestring(rbind(c(40, 70), c(40 + radius_ref_deg, 70))),
  crs = crs_wgs84
) %>% st_transform(crs_laea)
meters_per_degree <- as.numeric(st_length(radius_ref_line))

# Переводим radius_deg > radius_m
DATA_attrs <- DATA_attrs %>%
  mutate(radius = radius_deg * meters_per_degree)

# Формируем финальный датафрейм для scatterpie
pie_data_proj <- DATA_coords %>%
  bind_cols(DATA_attrs)

# Преобразуем границы области в проекцию (для coord_sf)
bbox_wgs84 <- st_sfc(st_polygon(list(rbind(
  c(xmin, ymin),
  c(xmax, ymin),
  c(xmax, ymax),
  c(xmin, ymax),
  c(xmin, ymin)
))), crs = crs_wgs84)

bbox_proj <- st_transform(bbox_wgs84, crs_laea)
bbox_coords <- st_bbox(bbox_proj)

####### ПОСТРОЕНИЕ КАРТЫ ################

ggplot() +
  # Базовая карта России
  geom_sf(data = russia_proj, fill = "lightblue", color = "gray") + 
  # Ограничение области
  coord_sf(xlim = c(bbox_coords$xmin, bbox_coords$xmax),
           ylim = c(bbox_coords$ymin, bbox_coords$ymax)) +
  # ?????? ДОБАВЛЕНА МАСШТАБНАЯ ЛИНЕЙКА ВВЕРХУ СПРАВА
  annotation_scale(location = "tr", width_hint = 0.3) +
  # Добавление множества pie charts с переменным радиусом
  geom_scatterpie(aes(x = X, y = Y, r = radius), 
                  data = pie_data_proj,
                  cols = c("Самцы пром.", "Самцы пререкруты I", 
                           "Самцы пререкруты II", "Самцы молодь"),
                  color = "black",
                  alpha = 0.8) +
  # Цвета для секторов
  scale_fill_manual(
    values = c(
      "Самцы пром." = "#E74C3C",           # красный
      "Самцы пререкруты I" = "#3498DB",    # синий
      "Самцы пререкруты II" = "#2ECC71",   # зелёный
      "Самцы молодь" = "#F39C12"           # оранжевый
    ),
    name = "Типы самцов"
  ) +
  # Подписи и тема
  labs(title = "Распределение уловов самцов камчатского краба по категориям",
       x = "Долгота", 
       y = "Широта") +
  theme_minimal() +
  # Настройка легенды: внизу, в 2 строки
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 10),
    plot.title = element_text(hjust = 0.5),
    panel.border = element_rect(color = "black", linewidth = 1, fill = NA)) +
  guides(fill = guide_legend(nrow = 2))   # разбиваем легенду на 2 строки