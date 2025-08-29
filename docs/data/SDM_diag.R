# ========================================================================================================================
# ПРАКТИЧЕСКОЕ ЗАНЯТИЕ:  Модели пространственного распределения видов (SDM) с biomod2
# 
# Автор: Баканев С. В.  | Обновлено: Sys.Date()
# ========================================================================================================================

# Библиотеки ------------------------------------------------------------------------------------------------------------ #
suppressPackageStartupMessages({
  library(biomod2)
  library(sf)
  library(marmap)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(readr)
  library(pROC)
  library(precrec)
  library(ecospat)
  library(dismo)
  library(rnaturalearth)
  library(ggspatial)
  library(raster)
})

# Опции и воспроизводимость -------------------------------------------------------------------------------------------- #
options(stringsAsFactors = FALSE)
set.seed(42)

 setwd("C:/SDM")  # при необходимости

# Хелперы устойчивые к ошибкам ------------------------------------------------------------------------------------------ #
ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

scale_predictions_01 <- function(predictions_numeric) {
  mx <- suppressWarnings(max(predictions_numeric, na.rm = TRUE))
  if (is.finite(mx) && mx > 1.5) return(pmin(pmax(predictions_numeric / 1000, 0), 1))
  pmin(pmax(predictions_numeric, 0), 1)
}

calibration_table <- function(labels_binary, probs_01, num_bins = 10) {
  stopifnot(length(labels_binary) == length(probs_01))
  idx <- is.finite(probs_01) & is.finite(labels_binary)
  labels_binary <- as.integer(labels_binary[idx])
  probs_01 <- as.numeric(probs_01[idx])
  keep <- labels_binary %in% c(0, 1)
  labels_binary <- labels_binary[keep]
  probs_01 <- probs_01[keep]
  if (!length(probs_01)) return(list(table = tibble(), brier = NA_real_, ece = NA_real_))
  breaks <- seq(0, 1, length.out = num_bins + 1)
  bin_id <- cut(probs_01, breaks = breaks, include.lowest = TRUE, labels = FALSE)
  mids <- (breaks[-length(breaks)] + breaks[-1]) / 2
  tb <- tibble(bin_id = bin_id, prob = probs_01, label = labels_binary) %>%
    group_by(bin_id) %>%
    summarise(
      bin_mid = mids[unique(bin_id)],
      prob_mean = mean(prob, na.rm = TRUE),
      obs_rate = mean(label, na.rm = TRUE),
      n = dplyr::n(),
      .groups = "drop"
    ) %>% arrange(bin_id)
  brier <- mean((probs_01 - labels_binary)^2, na.rm = TRUE)
  ece <- sum((tb$n / sum(tb$n)) * abs(tb$obs_rate - tb$prob_mean))
  list(table = tb, brier = brier, ece = ece)
}

plot_calibration <- function(tbl, title = "Калибровка (reliability)") {
  ggplot(tbl, aes(x = prob_mean, y = obs_rate, size = n)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
    geom_point(color = "#2C7FB8", alpha = 0.85) +
    scale_size_continuous(name = "N") +
    coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(x = "Средняя предсказанная вероятность", y = "Наблюдаемая доля присутствий", title = title) +
    theme_minimal(base_size = 12)
}

roc_pr_metrics <- function(labels_binary, probs_01) {
  idx <- is.finite(probs_01) & is.finite(labels_binary)
  labels_binary <- as.integer(labels_binary[idx])
  probs_01 <- as.numeric(probs_01[idx])
  keep <- labels_binary %in% c(0, 1)
  labels_binary <- labels_binary[keep]
  probs_01 <- probs_01[keep]
  if (length(unique(labels_binary)) < 2) {
    return(list(roc = NULL, auc_roc = NA_real_, pr = NULL, auc_pr = NA_real_))
  }
  roc_obj <- tryCatch(pROC::roc(response = labels_binary, predictor = probs_01, quiet = TRUE), error = function(e) NULL)
  auc_roc <- if (!is.null(roc_obj)) as.numeric(pROC::auc(roc_obj)[1]) else NA_real_
  pr_obj <- tryCatch(precrec::evalmod(scores = probs_01, labels = labels_binary), error = function(e) NULL)
  auc_pr <- if (!is.null(pr_obj)) {
    dd <- precrec::auc(pr_obj)
    as.numeric(dd %>% dplyr::filter(curvetypes == "PRC") %>% dplyr::pull(aucs))
  } else NA_real_
  list(roc = roc_obj, auc_roc = auc_roc, pr = pr_obj, auc_pr = auc_pr)
}

plot_roc_curve <- function(roc_obj, title = "ROC кривая") {
  if (is.null(roc_obj)) return(ggplot() + labs(title = paste(title, "(недостаточно классов)")))
  pROC::ggroc(roc_obj, colour = "#1B9E77", size = 1) +
    geom_abline(slope = 1, intercept = 1, linetype = "dashed", color = "gray50") +
    coord_equal(xlim = c(1, 0), ylim = c(0, 1)) +
    labs(x = "1 - Specificity", y = "Sensitivity", title = title) +
    theme_minimal(base_size = 12)
}

plot_pr_curve <- function(pr_obj, title = "PR кривая") {
  if (is.null(pr_obj)) return(ggplot() + labs(title = paste(title, "(недостаточно классов)")))
  autoplot(pr_obj) + labs(title = title) + theme_minimal(base_size = 12)
}

optimal_threshold_tss <- function(labels_binary, probs_01, step = 0.01) {
  thresholds <- seq(0, 1, by = step)
  idx <- is.finite(probs_01) & is.finite(labels_binary)
  labels_binary <- as.integer(labels_binary[idx])
  probs_01 <- as.numeric(probs_01[idx])
  keep <- labels_binary %in% c(0, 1)
  labels_binary <- labels_binary[keep]
  probs_01 <- probs_01[keep]
  if (!length(probs_01)) {
    return(data.frame(threshold = NA_real_, TSS = NA_real_, Sensitivity = NA_real_, Specificity = NA_real_))
  }
  metrics <- lapply(thresholds, function(th) {
    pred_class <- as.integer(probs_01 >= th)
    tp <- sum(pred_class == 1 & labels_binary == 1)
    tn <- sum(pred_class == 0 & labels_binary == 0)
    fp <- sum(pred_class == 1 & labels_binary == 0)
    fn <- sum(pred_class == 0 & labels_binary == 1)
    tpr <- ifelse((tp + fn) > 0, tp / (tp + fn), NA_real_)
    tnr <- ifelse((tn + fp) > 0, tn / (tn + fp), NA_real_)
    c(threshold = th, TSS = (tpr + tnr - 1), Sensitivity = tpr, Specificity = tnr)
  })
  m <- do.call(rbind, metrics)
  m <- as.data.frame(m)
  m$threshold <- as.numeric(m$threshold)
  m$TSS <- as.numeric(m$TSS)
  m$Sensitivity <- as.numeric(m$Sensitivity)
  m$Specificity <- as.numeric(m$Specificity)
  best_row <- m[which.max(m$TSS), , drop = FALSE]
  best_row
}

boyce_index <- function(labels_binary, probs_01, num_class = 0, window_w = NULL) {
  idx <- is.finite(probs_01) & is.finite(labels_binary)
  labels_binary <- as.integer(labels_binary[idx])
  probs_01 <- as.numeric(probs_01[idx])
  pres <- probs_01[labels_binary == 1]
  back <- probs_01
  tryCatch({
    res <- ecospat::ecospat.boyce(fit = back, obs = pres, nclass = num_class, window.w = window_w)
    list(CBI = as.numeric(res$Spearman.cor), curve = res$F.ratio)
  }, error = function(e) list(CBI = NA_real_, curve = NULL))
}

uncertainty_per_point <- function(pred_long_df) {
  # Ожидаем столбцы: run, algo, points, pred
  if (!all(c("points", "pred") %in% names(pred_long_df))) {
    stop("Для неопределенности нужен длинный формат с колонками 'points' и 'pred'.")
  }
  run_levels <- unique(pred_long_df$run)
  run_pick <- if ("allRun" %in% run_levels) "allRun" else run_levels[1]
  df <- subset(pred_long_df, run == run_pick)
  # Агрегируем по точкам: SD и среднее по всем алгоритмам
  agg_sd <- aggregate(pred ~ points, data = df, FUN = function(x) sd(as.numeric(x), na.rm = TRUE))
  names(agg_sd)[2] <- "UNC_SD"
  agg_mean <- aggregate(pred ~ points, data = df, FUN = function(x) mean(as.numeric(x), na.rm = TRUE))
  names(agg_mean)[2] <- "MEAN_PRED"
  unc <- merge(agg_sd, agg_mean, by = "points", all = TRUE)
  unc$UNC_CV <- unc$UNC_SD / ifelse(unc$MEAN_PRED == 0, NA, unc$MEAN_PRED)
  unc
}

mess_scores <- function(reference_env_df, target_env_df) {
  tryCatch({
    as.numeric(dismo::mess(x = as.data.frame(target_env_df), v = as.data.frame(reference_env_df)))
  }, error = function(e) rep(NA_real_, nrow(target_env_df)))
}

# ДАННЫЕ: текущее состояние -------------------------------------------------------------------------------------------- #
DATA <- read.csv("final_sdm_table_with_na.csv")
str(DATA)
DataSpecies <- as.data.frame(DATA)
myRespName <- 'occ'
myResp <- as.numeric(DataSpecies[[myRespName]])
myRespXY <- DataSpecies[, c("x", "y")]
myExpl <- DataSpecies[, 4:12]

myBiomodData <- BIOMOD_FormatingData(
  resp.var = myResp,
  expl.var = myExpl,
  resp.xy = myRespXY,
  resp.name = myRespName
)
print(myBiomodData)
plot(myBiomodData)

# ОБУЧЕНИЕ ЕДИНИЧНЫХ МОДЕЛЕЙ ------------------------------------------------------------------------------------------- #
algos <- c("ANN", "CTA", "FDA", "GAM", "GBM", "GLM", "MAXENT", "MAXNET", "RF", "XGBOOST")
# Уберем MAXENT, если нет maxent.jar рядом с рабочей директорией
if (!file.exists(file.path(getwd(), "maxent.jar"))) {
  algos <- setdiff(algos, "MAXENT")
}

myBiomodModelOut <- BIOMOD_Modeling(
  bm.format = myBiomodData,
  modeling.id = 'AllModels',
  models = algos,
  CV.strategy = 'random',
  CV.nb.rep = 2,
  CV.perc = 0.8,
  OPT.strategy = 'bigboss',
  metric.eval = c('TSS','ROC'),
  var.import = 2,
  seed.val = 42
)
print(myBiomodModelOut)

# Оценки и важности
str(get_evaluations(myBiomodModelOut))
str(get_variables_importance(myBiomodModelOut))

bm_PlotEvalMean(bm.out = myBiomodModelOut, dataset = 'calibration')
bm_PlotEvalMean(bm.out = myBiomodModelOut, dataset = 'validation')
bm_PlotEvalBoxplot(bm.out = myBiomodModelOut, group.by = c('algo', 'run'))

# ПРОГНОЗ: текущее состояние -------------------------------------------------------------------------------------------- #
myBiomodProj <- BIOMOD_Projection(
  bm.mod = myBiomodModelOut,
  proj.name = 'Current',
  new.env = myExpl,
  models.chosen = 'all'
)

# Предсказания единичных моделей (long)
pred_current_single <- get_predictions(myBiomodProj)
print(head(pred_current_single))

# АНСАМБЛИРОВАНИЕ И ПРОГНОЗ АНСАМБЛЯ ----------------------------------------------------------------------------------- #
myBiomodEM <- BIOMOD_EnsembleModeling(
  bm.mod = myBiomodModelOut,
  models.chosen = 'all',
  em.by = 'all',
  em.algo = c('EMmean', 'EMca'),
  metric.select = c('TSS'),
  metric.select.thresh = c(0.4),
  metric.eval = c('TSS', 'ROC'),
  var.import = 3,
  seed.val = 42
)
str(get_evaluations(myBiomodEM))
str(get_variables_importance(myBiomodEM))

myBiomodEMProj <- BIOMOD_EnsembleForecasting(
  bm.em = myBiomodEM,
  bm.proj = myBiomodProj,
  models.chosen = 'all',
  metric.binary = 'all',
  metric.filter = 'all'
)

pred_current_em <- get_predictions(myBiomodEMProj)
pred_current_emmean <- dplyr::filter(pred_current_em, .data$algo == "EMmean")
print(head(pred_current_emmean))

# КАРТЫ: текущий период ------------------------------------------------------------------------------------------------ #
MAPDATA <- tibble(
  point_id = seq_len(nrow(DATA)),
  X = DATA$x,
  Y = DATA$y,
  PRED = pred_current_emmean$pred
)

world <- rnaturalearth::ne_countries(scale = 50, returnclass = 'sf')
xmin <- 10; xmax <- 45; ymin <- 66; ymax <- 72

bat <- tryCatch(getNOAA.bathy(xmin, xmax, ymin, ymax, resolution = 4), error = function(e) NULL)
bat_xyz <- if (!is.null(bat)) as.xyz(bat) else NULL

p_points <- ggplot() +
  geom_sf(data = world) +
  coord_sf(xlim = c(xmin, xmax), ylim = c(ymin, ymax)) +
  {if (!is.null(bat_xyz)) geom_tile(data = bat_xyz, aes(x = V1, y = V2, fill = V3), show.legend = FALSE)} +
  geom_point(data = MAPDATA, aes(x = X, y = Y, size = PRED), color = "black", fill = "white", shape = 21, alpha=0.8) +
  ggspatial::annotation_scale(location = "tr", width_hint = 0.5) +
  scale_size(name = "Вероятность", range = c(1, 5)) +
  labs(title = "Точки: интенсивность предсказания (Current)")

p_raster <- ggplot() +
  geom_raster(data = MAPDATA, aes(x = X, y = Y, fill = PRED), interpolate = FALSE) +
  scale_fill_viridis_c(option = "D", name = "PRED") +
  geom_sf(data = world, color = "gray30", fill = "#E8E5D6", lwd = 0.3) +
  coord_sf(xlim = c(xmin*1.2, xmax*0.96), ylim = c(ymin*1.02, ymax*0.99)) +
  labs(title = "Растер: предсказание EMmean (Current)")

print(p_points)
print(p_raster)

# ДИАГНОСТИКА НАДЕЖНОСТИ: текущее -------------------------------------------------------------------------------------- #
labels <- as.integer(myResp)
probs_current_01 <- scale_predictions_01(pred_current_emmean$pred)

calib <- calibration_table(labels, probs_current_01, num_bins = 10)
print(calib$table)
print(plot_calibration(calib$table) + labs(subtitle = sprintf("Brier = %.3f, ECE = %.3f", calib$brier, calib$ece)))

rp <- roc_pr_metrics(labels, probs_current_01)
print(plot_roc_curve(rp$roc) + labs(subtitle = sprintf("AUC = %.3f", rp$auc_roc)))
print(plot_pr_curve(rp$pr) + labs(subtitle = sprintf("AUPRC = %.3f", rp$auc_pr)))

opt_thr <- optimal_threshold_tss(labels, probs_current_01, step = 0.005)
print(opt_thr)

boy <- boyce_index(labels, probs_current_01)
print(tibble(metric = "Continuous Boyce Index", value = boy$CBI))

# НЕОПРЕДЕЛЕННОСТЬ МЕЖДУ АЛГОРИТМАМИ (SD, CV) -------------------------------------------------------------------------- #
unc <- uncertainty_per_point(pred_current_single)
MAPDATA_unc <- MAPDATA %>% left_join(unc %>% transmute(point_id = points, UNC_SD, UNC_CV), by = "point_id")

p_unc_sd <- ggplot(MAPDATA_unc, aes(x = X, y = Y, fill = UNC_SD)) +
  geom_raster() +
  scale_fill_viridis_c(option = "C", name = "SD") +
  coord_equal() +
  labs(title = "Неопределенность (SD) между алгоритмами — Current")

print(p_unc_sd)

# БУДУЩЕЕ: данные, прогноз и диагностика -------------------------------------------------------------------------------- #
DATA_F <- read.csv("future_sdm_table_with_na.csv")
str(DATA_F)
DataSpeciesF <- as.data.frame(DATA_F)
myRespF <- as.numeric(DataSpeciesF[[myRespName]])
myRespXYF <- DataSpeciesF[, c("x", "y")]
myExplP1 <- DataSpeciesF[, 4:12]

myBiomodProj1 <- BIOMOD_Projection(
  bm.mod = myBiomodModelOut,
  proj.name = 'Future',
  new.env = myExplP1,
  models.chosen = 'all'
)

myBiomodEMProj1 <- BIOMOD_EnsembleForecasting(
  bm.em = myBiomodEM,
  bm.proj = myBiomodProj1,
  models.chosen = 'all',
  metric.binary = 'all',
  metric.filter = 'all'
)

pred_future_em <- get_predictions(myBiomodEMProj1)
pred_future_emmean <- dplyr::filter(pred_future_em, .data$algo == "EMmean")

MAPDATA2 <- tibble(
  point_id = seq_len(nrow(DATA_F)),
  X = DATA_F$x,
  Y = DATA_F$y,
  PRED = pred_future_emmean$pred
)

# Карты: будущее и ? ---------------------------------------------------------------------------------------------------- #
p_points2 <- ggplot() +
  geom_sf(data = world) +
  coord_sf(xlim = c(xmin, xmax), ylim = c(ymin, ymax)) +
  {if (!is.null(bat_xyz)) geom_tile(data = bat_xyz, aes(x = V1, y = V2, fill = V3), show.legend = FALSE)} +
  geom_point(data = MAPDATA2, aes(x = X, y = Y, size = PRED), color = "black", fill = "white", shape = 21, alpha=0.8) +
  ggspatial::annotation_scale(location = "tr", width_hint = 0.5) +
  scale_size(name = "Вероятность", range = c(1, 5)) +
  labs(title = "Точки: EMmean (Future)")

p_raster2 <- ggplot() +
  geom_raster(data = MAPDATA2, aes(x = X, y = Y, fill = PRED)) +
  scale_fill_viridis_c(option = "D", name = "PRED") +
  geom_sf(data = world, color = "gray30", fill = "#E8E5D6", lwd = 0.3) +
  coord_sf(xlim = c(xmin*1.2, xmax*0.96), ylim = c(ymin*1.01, ymax*0.99)) +
  labs(title = "Растер: EMmean (Future)")

probs_future_01 <- scale_predictions_01(pred_future_emmean$pred)
probs_current_01 <- scale_predictions_01(pred_current_emmean$pred)
MAPDATA_delta <- tibble(X = MAPDATA$X, Y = MAPDATA$Y, delta = probs_future_01 - probs_current_01)

p_delta <- ggplot(MAPDATA_delta, aes(x = X, y = Y, fill = delta)) +
  geom_raster() +
  scale_fill_gradient2(low = "#D7301F", mid = "#FFFFBF", high = "#1A9850", midpoint = 0, name = "? Prob") +
  coord_equal() +
  labs(title = "? (Future ? Current) EMmean (scaled)")

print(p_points2)
print(p_raster2)
print(p_delta)

# Неопределенность (будущее) и MESS ------------------------------------------------------------------------------------ #
pred_future_single <- get_predictions(myBiomodProj1)
unc_f <- uncertainty_per_point(pred_future_single)
MAPDATA2_unc <- MAPDATA2 %>% left_join(unc_f %>% transmute(point_id = points, UNC_SD = UNC_SD, UNC_CV = UNC_CV), by = "point_id")

p_unc_sd_f <- ggplot(MAPDATA2_unc, aes(x = X, y = Y, fill = UNC_SD)) +
  geom_raster() +
  scale_fill_viridis_c(option = "C", name = "SD") +
  coord_equal() +
  labs(title = "Неопределенность (SD) между алгоритмами — Future")

common <- names(myExpl)
# убрать предикторы без размаха в референсе
ok <- sapply(myExpl[, common, drop = FALSE], function(z) length(unique(na.omit(z))) >= 2)
vars <- common[ok]

# матрица референса без NA-строк
ref_mat <- as.matrix(na.omit(myExpl[, vars, drop = FALSE]))

# собрать RasterStack для будущего из точек x,y и значений предикторов
layers <- lapply(vars, function(nm) {
  rasterFromXYZ(data.frame(x = DATA_F$x, y = DATA_F$y, z = myExplP1[, nm]))
})
x_stack <- stack(layers); names(x_stack) <- vars

# посчитать MESS как растер
mess_r <- dismo::mess(x_stack, ref_mat)

# извлечь значения MESS в порядке строк будущего набора
mess_vals <- raster::extract(mess_r, DATA_F[, c("x","y")])

MAPDATA2_MESS <- dplyr::mutate(MAPDATA2, MESS = mess_vals)

p_mess <- ggplot(MAPDATA2_MESS, aes(x = X, y = Y, fill = MESS)) +
  geom_raster() +
  scale_fill_gradient2(low = "#762A83", mid = "#F7F7F7", high = "#1B7837", midpoint = 0, name = "MESS") +
  coord_equal() +
  labs(title = "MESS (экстраполяционный риск) — Future vs Current")

print(p_unc_sd_f)
print(p_mess)

# Конец скрипта --------------------------------------------------------------------------------------------------------- #