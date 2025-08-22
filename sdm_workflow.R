#!/usr/bin/env Rscript

# ========================================================================================================================
# ПРАКТИЧЕСКОЕ ЗАНЯТИЕ: Модели пространственного распределения видов (SDM)
# Курс: "Оценка водных биоресурсов в среде R (для начинающих)"
# Автор: Баканев С. В.    Дата: 21.08.2025
# Версия сценария: 1.0 (оптимизировано и снабжено подробными комментариями)
# ========================================================================================================================
# Структура:
# 1) Подготовка окружения, путей и пакетов
# 2) Загрузка и валидация данных
# 3) Подготовка данных для biomod2
# 4) Построение одиночных моделей (SDM) и их оценка
# 5) Проекции (текущие условия)
# 6) Модели ансамблей, их оценка и прогнозы
# 7) Экспорт результатов (оценки, важность переменных, предсказания)
# 8) Картографирование: фоновые данные, батиметрия, точки вероятности и растры
# 9) Диагностика, воспроизводимость и сохранение сессии
# ========================================================================================================================

suppressPackageStartupMessages({
  if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman", repos = "https://cloud.r-project.org")
  pacman::p_load(
    biomod2, sf, marmap, tidyverse, rnaturalearth, rnaturalearthdata,
    ggspatial, viridis, readr, dplyr, ggplot2, janitor, mgcv
  )
})

# -----------------------------------------------
# 1) Подготовка окружения, путей и пакетов
# -----------------------------------------------
set.seed(42)

# Рабочая директория: используем переменную окружения или создаём стандартную структуру в /workspace
sdm_root <- Sys.getenv("SDM_WORKDIR", unset = "/workspace/SDM")
if (!dir.exists(sdm_root)) dir.create(sdm_root, recursive = TRUE, showWarnings = FALSE)
setwd(sdm_root)

# Поддиректории проекта
paths <- list(
  data      = file.path(sdm_root, "data"),
  outputs   = file.path(sdm_root, "outputs"),
  figures   = file.path(sdm_root, "figures"),
  cache     = file.path(sdm_root, "cache")
)
for (pth in paths) if (!dir.exists(pth)) dir.create(pth, recursive = TRUE, showWarnings = FALSE)

# Имя файла с данными (ожидаются колонки: occ, x, y, далее предикторы)
input_csv <- file.path(paths$data, "final_sdm_table_with_na.csv")

# Вспомогательная функция безопасной записи CSV
safe_write_csv <- function(x, file, ...) {
  tryCatch({
    readr::write_csv(x, file, na = "NA")
    message("[write] ", normalizePath(file, mustWork = FALSE))
  }, error = function(e) {
    warning("Не удалось записать ", file, ": ", conditionMessage(e))
  })
}

# -----------------------------------------------
# 2) Загрузка и валидация данных
# -----------------------------------------------
if (!file.exists(input_csv)) {
  stop(paste0("Не найден входной файл: ", input_csv, "\nПоместите CSV с колонками occ, x, y и предикторами в: ", paths$data))
}

raw_data <- suppressMessages(readr::read_csv(input_csv, show_col_types = FALSE))
raw_data <- janitor::clean_names(raw_data)

required_cols <- c("occ", "x", "y")
missing_required <- setdiff(required_cols, names(raw_data))
if (length(missing_required) > 0) {
  stop(paste0("Отсутствуют обязательные колонки: ", paste(missing_required, collapse = ", ")))
}

# Базовые проверки
if (!is.numeric(raw_data$occ)) {
  # Попытаться привести к числу (0/1)
  raw_data <- raw_data %>% mutate(occ = as.numeric(occ))
}
if (!all(raw_data$occ %in% c(0,1) | is.na(raw_data$occ))) {
  stop("Колонка occ должна содержать только 0/1 или NA")
}

if (!is.numeric(raw_data$x) || !is.numeric(raw_data$y)) stop("Колонки x и y должны быть числовыми (долгота/широта)")

# Удалим строки без координат или ответа
data_valid <- raw_data %>% filter(!is.na(occ), !is.na(x), !is.na(y))

# Предикторы: все столбцы начиная с 4-го, если пользователь не задал иначе
if (ncol(data_valid) < 4) stop("В таблице недостаточно столбцов для предикторов (ожидается >= 4)")
predictor_cols <- names(data_valid)[4:ncol(data_valid)]

# Приведение предикторов к numeric, удаление некорректных
predictors <- data_valid %>% dplyr::select(all_of(predictor_cols))

non_numeric <- predictor_cols[!sapply(predictors, is.numeric)]
if (length(non_numeric) > 0) {
  warning("Ненумерические предикторы обнаружены и будут исключены: ", paste(non_numeric, collapse = ", "))
  predictor_cols <- setdiff(predictor_cols, non_numeric)
}

# Удалим строки с NA среди выбранных предикторов
model_df <- data_valid %>%
  dplyr::select(occ, x, y, all_of(predictor_cols)) %>%
  tidyr::drop_na()

if (nrow(model_df) < 50) warning("Мало наблюдений после очистки (", nrow(model_df), ") — результаты могут быть нестабильны")

# Лёгкая диагностика распределений и отсутствий
message("Размерность данных для модели: ", nrow(model_df), " x ", ncol(model_df))
message("Доля присутствий: ", round(mean(model_df$occ), 3))

# -----------------------------------------------
# 3) Подготовка данных для biomod2
# -----------------------------------------------
response_name <- "occ"
response_vec  <- as.numeric(model_df[[response_name]])
coords_df     <- model_df %>% dplyr::select(x, y)
expl_df       <- model_df %>% dplyr::select(all_of(predictor_cols))

# Форматирование данных
biomod_data <- BIOMOD_FormatingData(
  resp.var  = response_vec,
  expl.var  = expl_df,
  resp.xy   = coords_df,
  resp.name = response_name
)

# Сохраним предварительный обзор
try(plot(biomod_data))

# -----------------------------------------------
# 4) Построение одиночных моделей (SDM) и их оценка
# -----------------------------------------------
# Подбор стратегий кросс-валидации: spatial block если наблюдений достаточно, иначе random
cv_strategy <- if (nrow(model_df) >= 150) "block" else "random"
message("Стратегия кросс-валидации: ", cv_strategy)

# Набор алгоритмов: надёжные по умолчанию; MAXENT (java) исключён из-за зависимостей
candidate_algorithms <- c("GLM", "GAM", "GBM", "RF", "MAXNET", "XGBOOST")

# Параметры моделирования
modeling_id <- paste0("SDM_", format(Sys.time(), "%Y%m%d_%H%M%S"))

biomod_models <- BIOMOD_Modeling(
  bm.format    = biomod_data,
  modeling.id  = modeling_id,
  models       = candidate_algorithms,
  CV.strategy  = cv_strategy,
  CV.nb.rep    = 3,
  CV.perc      = 0.8,
  OPT.strategy = "bigboss",
  metric.eval  = c("TSS", "ROC", "KAPPA"),
  var.import   = 3,
  seed.val     = 42
)

# Оценки и важность переменных по одиночным моделям
single_eval <- get_evaluations(biomod_models)
single_imp  <- get_variables_importance(biomod_models)

# Визуализация метрик (сохраняем в файл)
try({
  p1 <- bm_PlotEvalMean(bm.out = biomod_models, dataset = "calibration")
  p2 <- bm_PlotEvalMean(bm.out = biomod_models, dataset = "validation")
  p3 <- bm_PlotEvalBoxplot(bm.out = biomod_models, group.by = c("algo", "run"))
  ggsave(filename = file.path(paths$figures, paste0(modeling_id, "_eval_calibration.png")), plot = p1, width = 8, height = 5, dpi = 200)
  ggsave(filename = file.path(paths$figures, paste0(modeling_id, "_eval_validation.png")),  plot = p2, width = 8, height = 5, dpi = 200)
  ggsave(filename = file.path(paths$figures, paste0(modeling_id, "_eval_boxplot.png")),     plot = p3, width = 9, height = 6, dpi = 200)
})

# -----------------------------------------------
# 5) Проекции (текущие условия)
# -----------------------------------------------
# Здесь проецируем на те же предикторы (табличные точки). Для растровой проекции требуется стек растров.
biomod_proj <- BIOMOD_Projection(
  bm.mod       = biomod_models,
  proj.name    = "Current",
  new.env      = expl_df,
  models.chosen = "all"
)

# -----------------------------------------------
# 6) Ансамблевые модели, их оценка и прогнозы
# -----------------------------------------------
biomod_em <- BIOMOD_EnsembleModeling(
  bm.mod                 = biomod_models,
  models.chosen          = "all",
  em.by                  = "all",
  em.algo                = c("EMmean", "EMwmean", "EMca"),
  metric.select          = c("TSS"),
  metric.select.thresh   = c(0.6),
  metric.eval            = c("TSS", "ROC"),
  var.import             = 3,
  seed.val               = 42
)

# Оценки и важность переменных по ансамблям
em_eval <- get_evaluations(biomod_em)
em_imp  <- get_variables_importance(biomod_em)

# Прогноз ансамблей на табличные точки
biomod_em_proj <- BIOMOD_EnsembleForecasting(
  bm.em         = biomod_em,
  bm.proj       = biomod_proj,
  models.chosen = "all",
  metric.binary = "all",
  metric.filter = "all"
)

# Альтернатива: построение одиночной проекции заново внутри ансамблевой функции
# biomod_em_proj <- BIOMOD_EnsembleForecasting(
#   bm.em         = biomod_em,
#   proj.name     = "CurrentEM",
#   new.env       = expl_df,
#   models.chosen = "all",
#   metric.binary = "all",
#   metric.filter = "all"
# )

# -----------------------------------------------
# 7) Экспорт результатов
# -----------------------------------------------
# Предсказания на точки (фильтруем EMmean как основную метрику)
preds_all <- get_predictions(biomod_em_proj)
if (!is.null(preds_all) && all(c("algo", "pred") %in% names(preds_all))) {
  preds_emmean <- preds_all %>% dplyr::filter(algo == "EMmean")
} else {
  warning("Структура предсказаний нестандартна — сохраняем как есть")
  preds_emmean <- preds_all
}

# Собираем таблицу для карты
map_data <- tibble::tibble(
  x = model_df$x,
  y = model_df$y,
  pred = if (!is.null(preds_emmean) && "pred" %in% names(preds_emmean)) preds_emmean$pred else NA_real_
)

# Сохраняем CSV с оценками и важностью
if (!is.null(single_eval)) safe_write_csv(single_eval, file.path(paths$outputs, paste0(modeling_id, "_single_eval.csv")))
if (!is.null(single_imp))  safe_write_csv(single_imp,  file.path(paths$outputs, paste0(modeling_id, "_single_imp.csv")))
if (!is.null(em_eval))     safe_write_csv(em_eval,    file.path(paths$outputs, paste0(modeling_id, "_ensemble_eval.csv")))
if (!is.null(em_imp))      safe_write_csv(em_imp,     file.path(paths$outputs, paste0(modeling_id, "_ensemble_imp.csv")))
if (!is.null(map_data))    safe_write_csv(map_data,   file.path(paths$outputs, paste0(modeling_id, "_map_points.csv")))

# -----------------------------------------------
# 8) Картографирование: мир, батиметрия, точки и растровая интерполяция
# -----------------------------------------------
# Границы области интереса (настройте под свои данные)
lon_min <- min(model_df$x, na.rm = TRUE) - 1
lon_max <- max(model_df$x, na.rm = TRUE) + 1
lat_min <- min(model_df$y, na.rm = TRUE) - 1
lat_max <- max(model_df$y, na.rm = TRUE) + 1

# Мировой контур
world <- tryCatch({
  rnaturalearth::ne_countries(scale = 50, returnclass = 'sf')
}, error = function(e) {
  warning("Не удалось получить слой стран: ", conditionMessage(e), "; карта будет без береговой линии")
  NULL
})

# Батиметрия NOAA с устойчивой обработкой ошибок
bathy_xyz <- NULL
try({
  bat <- marmap::getNOAA.bathy(x1 = lon_min, x2 = lon_max, y1 = lat_min, y2 = lat_max, resolution = 1)
  bathy_xyz <- marmap::as.xyz(bat)
  bathy_xyz$bin <- cut(bathy_xyz$z,
                       breaks = c(-3500,-500,-300,-200,-100,-50,-25,0,12.5,25,50,75,100,200,300,400,500,625,750,1000,1500),
                       include.lowest = TRUE)
}, silent = TRUE)

# Палитра для батиметрии
bathy_cols <- c("#577f92","#577f92","#5f8a9f","#79a0ac","#9ebac2","#c3d4d8","#d0dde1",
                "#387374","#4b8d85","#519777","#528c7e","#589684","#55a98f","#7dbaa2",
                "#a5cbb4","#bfd8b6","#d8e5b7","#d6df98","#d3d878","#d8ce4e","#83492e")

# Точки вероятности
p_points <- ggplot() +
  { if (!is.null(world)) geom_sf(data = world, linewidth = 0.2) } +
  coord_sf(xlim = c(lon_min, lon_max), ylim = c(lat_min, lat_max), expand = FALSE) +
  { if (!is.null(bathy_xyz)) geom_tile(data = bathy_xyz, aes(x = x, y = y, fill = bin), show.legend = FALSE, alpha = 0.8) } +
  { if (!is.null(bathy_xyz)) scale_fill_manual(name = "Глубина", values = bathy_cols, drop = FALSE) } +
  geom_point(data = map_data, aes(x = x, y = y, size = pred), color = "black", fill = "white", shape = 21, alpha = 0.8, stroke = 0.2) +
  ggspatial::annotation_scale(location = "tr", width_hint = 0.4) +
  scale_size(name = "Вероятность", range = c(1.2, 5)) +
  labs(title = "SDM: вероятности по точкам (EMmean)") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "right")

# Псевдо-растр по точкам (визуализация поверх карты, без интерполяции)
p_raster <- ggplot() +
  geom_raster(data = map_data, aes(x = x, y = y, fill = pred), interpolate = FALSE, alpha = 0.9) +
  scale_fill_viridis_c(option = "D", na.value = "transparent", name = "Вероятность") +
  { if (!is.null(world)) geom_sf(data = world, linewidth = 0.2, fill = NA) } +
  coord_sf(xlim = c(lon_min, lon_max), ylim = c(lat_min, lat_max), expand = FALSE) +
  ggspatial::annotation_scale(location = "tr", width_hint = 0.4) +
  labs(title = "SDM: растровое представление значений по точкам") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "right")

# Сохранение карт
# Сглаженная поверхность по GAM (тонкие пластинные сплайны)
p_smooth <- NULL
try({
  md_nn <- map_data %>% dplyr::filter(!is.na(pred))
  if (nrow(md_nn) >= 50) {
    k_basis <- min(200, max(20L, floor(nrow(md_nn) / 5)))
    res_x <- max((lon_max - lon_min) / 200, 0.05)
    res_y <- max((lat_max - lat_min) / 200, 0.05)
    grid <- expand.grid(
      x = seq(lon_min, lon_max, by = res_x),
      y = seq(lat_min, lat_max, by = res_y)
    )
    gam_fit <- mgcv::gam(pred ~ s(x, y, bs = "tp", k = k_basis), data = md_nn, family = gaussian())
    grid$pred_smooth <- as.numeric(predict(gam_fit, newdata = grid, type = "response"))
    grid$pred_smooth <- pmin(pmax(grid$pred_smooth, 0), 1)
    p_smooth <- ggplot() +
      geom_raster(data = grid, aes(x = x, y = y, fill = pred_smooth), interpolate = TRUE) +
      scale_fill_viridis_c(option = "D", name = "Вероятность") +
      { if (!is.null(world)) geom_sf(data = world, linewidth = 0.2, fill = NA) } +
      coord_sf(xlim = c(lon_min, lon_max), ylim = c(lat_min, lat_max), expand = FALSE) +
      ggspatial::annotation_scale(location = "tr", width_hint = 0.4) +
      labs(title = "SDM: сглаженная поверхность (GAM)") +
      theme_minimal(base_size = 11) +
      theme(legend.position = "right")
    ggsave(filename = file.path(paths$figures, paste0(modeling_id, "_map_smooth_gam.png")), plot = p_smooth, width = 8, height = 6, dpi = 220)
  }
}, silent = TRUE)
try(ggsave(filename = file.path(paths$figures, paste0(modeling_id, "_map_points.png")), plot = p_points, width = 8, height = 6, dpi = 220))
try(ggsave(filename = file.path(paths$figures, paste0(modeling_id, "_map_raster.png")), plot = p_raster, width = 8, height = 6, dpi = 220))

# -----------------------------------------------
# 9) Диагностика и воспроизводимость
# -----------------------------------------------
# Сводка сессии и сохранение объектов
sess_file <- file.path(paths$outputs, paste0(modeling_id, "_sessionInfo.txt"))
try({
  con <- file(sess_file, open = "wt", encoding = "UTF-8")
  writeLines(c(
    paste0("Дата/время: ", as.character(Sys.time())),
    paste0("Рабочая директория: ", getwd()),
    "\nSessionInfo:",
    capture.output(sessionInfo())
  ), con = con)
  close(con)
})

# Сохранение ключевых объектов RDS
try(saveRDS(biomod_models, file.path(paths$cache, paste0(modeling_id, "_biomod_models.rds"))))
try(saveRDS(biomod_em,     file.path(paths$cache, paste0(modeling_id, "_biomod_em.rds"))))
try(saveRDS(biomod_proj,   file.path(paths$cache, paste0(modeling_id, "_biomod_proj.rds"))))
try(saveRDS(biomod_em_proj,file.path(paths$cache, paste0(modeling_id, "_biomod_em_proj.rds"))))

message("Готово. Результаты сохранены в: ", normalizePath(paths$outputs, mustWork = FALSE))