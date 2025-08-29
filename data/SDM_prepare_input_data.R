# ========================================================================================================================
# ПОДГОТОВКА ДАННЫХ ДЛЯ МОДЕЛИРОВАНИЯ ПРОСТРАНСТВЕННОГО РАСПРЕДЕЛЕНИЯ ВИДОВ (SDM)
# 
# Скрипт выполняет:
# 1. Загрузку и фильтрацию данных наблюдений
# 2. Создание регулярной сетки для анализа
# 3. Агрегацию данных по ячейкам сетки
# 4. Исключение наземных территорий
# 5. Расчет расстояния до берега
# 6. Извлечение био-физических переменных из NetCDF файлов
# 7. Сохранение итогового набора данных
# 
# Курс: "Оценка водных биоресурсов в среде R (для начинающих)"
# Автор: Баканев С. В. 
# Дата: 27.08.2025
# ========================================================================================================================
# Очистка рабочей среды
rm(list = ls())

# Установка рабочей директории
setwd("C:/SDM")

# Загрузка необходимых библиотек
library(tidyverse)    # Обработка данных и визуализация
library(readxl)       # Чтение Excel-файлов
library(rnaturalearth) # Векторные карты мира
library(sf)           # Пространственный анализ
library(ggOceanMaps)  # Расчет дистанции до берега
library(terra)        # Работа с растровыми данными

# ---------------------------
# 2. ЗАГРУЗКА И ФИЛЬТРАЦИЯ ДАННЫХ
# ---------------------------
DATA <- read_excel("PECTEN.xlsx", sheet = "PECTEN")
str(DATA)
# Получение границ Европы
europe <- ne_countries(scale = 10, continent = "Europe")  # Загрузка векторных границ Европы (масштаб 1:10м)


# Установка границ отображаемой области (долгота/широта)
xmin=10  # Западная граница
xmax=45  # Восточная граница
ymin=66 # Южная граница
ymax=72 # Северная граница

# Фильтрация данных по заданным координатам
PECTEN <- DATA %>%
  filter(X >= xmin & X <= xmax & Y >= ymin & Y <= ymax) %>%
  select(X, Y, OCC)  # Выбор только нужных колонок

# Построение карты
ggplot() +
  # Базовая карта Европы
  geom_sf(data = europe, fill = "#E8E5D6") + 
  # Ограничение области отображения
  coord_sf(xlim = c(xmin, xmax), ylim = c(ymin, ymax)) +
  # Точки наблюдений с размером и цветом по переменной OCC
  geom_point(aes(x = X, y = Y, size = OCC, color = OCC),
             data = PECTEN, alpha = 0.6) +
  # Цветовая шкала (viridis, вариант H)
  scale_color_viridis_c(option = "H") +
  # Подписи осей
  labs(x = "Долгота", y = "Широта", 
       size = "Наличие вида", color = "Наличие вида",
       title = "Распределение вида")

# Просмотр отфильтрованной таблицы
str(PECTEN)

# ---------------------------
# 3. СОЗДАНИЕ СЕТКИ И АГРЕГАЦИЯ ДАННЫХ
# ---------------------------

# Если объект europe не в формате sf, приведем:
europe <- suppressWarnings(sf::st_as_sf(europe))

# Границы бинов
x_breaks <- seq(xmin, xmax, by = 0.1)
y_breaks <- seq(ymin, ymax, by = 0.1)

# 2.1. Присвоим наблюдениям индексы ячеек и посчитаем среднее OCC по ячейке
PECTEN_binned <- PECTEN %>%
  mutate(
    x_id = cut(X, breaks = x_breaks, include.lowest = TRUE, right = FALSE, labels = FALSE),
    y_id = cut(Y, breaks = y_breaks, include.lowest = TRUE, right = FALSE, labels = FALSE)
  ) %>%
  # точки, попавшие ровно в xmax/ymax, уйдут в NA — это нормально, т.к. верхняя граница полуоткрытая
  filter(!is.na(x_id), !is.na(y_id)) %>%
  group_by(x_id, y_id) %>%
  summarise(OCC = mean(OCC, na.rm = TRUE), .groups = "drop") %>%
  mutate(OCC = na_if(OCC, NaN))  # если все NA в ячейке > NA, а не NaN

# 2.2. Полная таблица ячеек (все комбинации), координаты центров
grid_cells <- tidyr::expand_grid(
  x_id = seq_along(head(x_breaks, -1)),
  y_id = seq_along(head(y_breaks, -1))
) %>%
  mutate(
    X_center = (x_breaks[x_id] + x_breaks[x_id + 1]) / 2,
    Y_center = (y_breaks[y_id] + y_breaks[y_id + 1]) / 2
  )

# 2.3. Приклеим средние OCC к полной решетке (пустые > NA)
grid_occ <- grid_cells %>%
  left_join(PECTEN_binned, by = c("x_id", "y_id"))

# 2.4. Оставим только ячейки океана (центры, не попадающие на сушу Европы)
centers_sf <- st_as_sf(grid_occ, coords = c("X_center", "Y_center"), crs = 4326, remove = FALSE)
on_land_mat <- st_intersects(centers_sf, europe, sparse = FALSE)
on_land <- apply(on_land_mat, 1, any)

grid_ocean <- grid_occ %>%
  filter(!on_land) %>%
  select(X_center, Y_center, OCC) %>%
  arrange(Y_center, X_center)

# Переименование столбцов
grid_ocean <- grid_ocean %>%
  rename(X = X_center,
         Y = Y_center,
         OCC = OCC)

# Вывод результата: таблица центров ячеек и средняя OCC (NA, если нет данных)

str(grid_ocean)

# Построение карты+ проверка сетки
ggplot() +
  # Базовая карта Европы
  geom_sf(data = europe, fill = "#E8E5D6") + 
  # Ограничение области отображения
  coord_sf(xlim = c(xmin, xmax), ylim = c(ymin, ymax)) +
  # Точки наблюдений с размером и цветом по переменной OCC
  geom_point(aes(x = X, y = Y, size = X, color = Y),
             data = grid_ocean, alpha = 0.6) +
  # Цветовая шкала (viridis, вариант H)
  scale_color_viridis_c(option = "H") +
  # Подписи осей
  labs(x = "Долгота", y = "Широта", 
       size = "Наличие вида", color = "Наличие вида",
       title = "Распределение сетки")



# Присоединение к базовой таблице (X, Y, OCC) переменной DIST (рсстояние до берега)
dt <- data.frame(lon = grid_ocean$X, lat = grid_ocean$Y)
dt <- dist2land(dt, verbose = FALSE)
qmap(dt, color = ldist) + scale_color_viridis_c()

DATA <- tibble(
  X = grid_ocean$X,
  Y = grid_ocean$Y,
  OCC = grid_ocean$OCC,
  DIST = dt$ldist
)
# Выводим первые несколько строк для проверки
str(DATA)
# Присоединение к базовой таблице (X, Y, OCC, DIST) переменных с сайта BioORACLE
#
TEMP <- rast("C:/SDM/BIO/Bathymetry.nc")
# визуализация переменной
plot(TEMP , xlim = c(xmin, xmax), ylim = c(ymin, ymax=72)) 
# Автоматическое присоединение переменных из всех nc файлов в папке BIO

# Получаем список всех nc файлов в папке BIO
nc_files <- list.files("C:/SDM/BIO", pattern = "\\.nc$", full.names = TRUE)

# Проверяем, что файлы найдены
if (length(nc_files) == 0) {
  stop("Не найдено nc файлов в папке C:/SDM/BIO")
}

# Создаем координатный датафрейм один раз
coord <- data.frame(x = DATA$X, y = DATA$Y)

# Проходим по всем nc файлам
for (file_path in nc_files) {
  # Извлекаем имя переменной из названия файла (убираем расширение .nc)
  var_name <- tools::file_path_sans_ext(basename(file_path))
  
  # Загружаем raster файл
  temp_rast <- rast(file_path)
  
  # Извлекаем значения для координат DATA
  extracted_values <- extract(temp_rast, coord, method = "bilinear")
  
  # Убираем столбец ID и оставляем только значения
  values <- extracted_values[, -1, drop = FALSE]
  
  # Переименовываем столбец в имя переменной
  colnames(values) <- var_name
  
  # Присоединяем к DATA
  DATA <- cbind(DATA, values)
  
  # Сообщение о прогрессе
  message("Добавлена переменная: ", var_name, " из файла: ", basename(file_path))
}

# Проверяем результат
str(DATA)
write.csv(DATA, "SDM_all_pred_full_set.csv", row.names = FALSE)