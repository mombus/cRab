# =============================================================
# ОЦЕНКА ЗАПАСА С ПОМОЩЬЮ JABBA (Just Another Bayesian Biomass Assessment)
# Улучшенная, переносимая версия (без setwd, с управлением зависимостями)
# =============================================================

# ------------------------- 0. НАСТРОЙКА СЕССИИ ---------------------------

# Безопасные опции
options(stringsAsFactors = FALSE)

# Настроить пользовательскую библиотеку пакетов (без прав root)
user_lib <- Sys.getenv("R_LIBS_USER", unset = "/workspace/Rlibs")
if (!dir.exists(user_lib)) dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(unique(c(user_lib, .libPaths())))

# Управление воспроизводимостью
set.seed(12345)

# Функция установки/загрузки пакетов
ensure_packages <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
  if (length(missing)) {
    install.packages(missing, repos = "https://cloud.r-project.org")
  }
  invisible(vapply(pkgs, require, character.only = TRUE, quietly = TRUE, FUN.VALUE = logical(1)))
}

# Основные зависимости
ensure_packages(c("remotes", "reshape2", "ggplot2", "coda", "rjags", "R2jags"))

# Установка JABBA (если не установлена) из GitHub
if (!requireNamespace("JABBA", quietly = TRUE)) {
  msg <- NULL
  try({ remotes::install_github("jabbamodel/JABBA", upgrade = "never") }, silent = TRUE)
  if (!requireNamespace("JABBA", quietly = TRUE)) {
    try({ remotes::install_github("henning-winker/JABBA", upgrade = "never") }, silent = TRUE)
  }
  if (!requireNamespace("JABBA", quietly = TRUE)) {
    stop("Не удалось установить пакет JABBA из GitHub. Проверьте интернет/доступ к GitHub.")
  }
}

# Подгрузить JABBA
suppressPackageStartupMessages(library(JABBA))

# ------------------------- 1. КОНФИГУРАЦИЯ ---------------------------

# Имя оценки и базовая папка результатов
assessment <- Sys.getenv("JABBA_ASSESSMENT", unset = "NEW_JABBA")
base_output_dir <- Sys.getenv("JABBA_OUTPUT_DIR", unset = getwd())
run_tag <- format(Sys.time(), "%Y%m%d-%H%M%S")
output_dir <- file.path(base_output_dir, paste0(assessment, "_", run_tag))

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
retro_dir <- file.path(output_dir, "retro")
if (!dir.exists(retro_dir)) dir.create(retro_dir, recursive = TRUE, showWarnings = FALSE)
fig_dir <- file.path(output_dir, "figures")
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# Быстрый режим (для отладки в CI/демо)
FAST_RUN <- as.logical(Sys.getenv("FAST_RUN", unset = "FALSE"))

# ------------------------- 2. ДАННЫЕ ---------------------------

Year <- 2005:2024
Catch <- c(5, 7, 6, 10, 14, 25, 28, 30, 32, 35, 25, 20, 15, 12, 10, 12, 10, 13, 11, 12)
CPUE  <- c(27.4, 26.8, 16.8, 23.0, 29.0, 30.0, 16.5, 17.2, 10.5, 14.6, 8.3, 11.4, 15.5, 13.8, 11.5, 15.3, 12.2, 15.6, 16.2, 13.4)
BESS  <- c(NA, 16.3, 20.7, 15.1, 18.6, 16.0, 13.8, 13.3, 11.7, 11.8, 9.3, 7.1, 8.0, 9.2, 10.3, 9.8, 10.3, 11.7, 13.7, 13.4)

stopifnot(length(Year) == length(Catch), length(Year) == length(CPUE), length(Year) == length(BESS))

catch_data <- data.frame(year = Year, catch = Catch)
cpue_data  <- data.frame(year = Year, CPUE = CPUE, BESS = BESS)

# CV=0.2 -> SE ~ CV при лог-нормальном приближении; здесь оставляем как константу
se_data <- data.frame(
  year = Year,
  CPUE = ifelse(is.na(CPUE), NA, 0.2),
  BESS = ifelse(is.na(BESS), NA, 0.2)
)

# ------------------- 3. POSTERIOR SETUP И ВХОД В JABBA --------------------

# Параметры модели
model_type <- "Schaefer"
scenario <- "SPiCT_adapted"

# Настройки MCMC зависят от FAST_RUN
mcmc_cfg <- if (isTRUE(FAST_RUN)) list(ni = 5000, nb = 1000, nt = 5, nc = 2) else list(ni = 50000, nb = 10000, nt = 5, nc = 2)

jbinput <- JABBA::build_jabba(
  catch = catch_data,
  cpue = cpue_data,
  se = se_data,
  assessment = assessment,
  scenario = scenario,
  model.type = model_type,
  sigma.est = TRUE,
  r.prior = c(0.2, 0.5),
  K.prior = c(189.6, 0.795),
  psi.prior = c(0.75, 0.25),
  igamma = c(0.001, 0.001),
  verbose = FALSE
)

# ------------------- 4. ЗАПУСК МОДЕЛИ --------------------

fit <- JABBA::fit_jabba(
  jbinput,
  ni = mcmc_cfg$ni,
  nb = mcmc_cfg$nb,
  nt = mcmc_cfg$nt,
  nc = mcmc_cfg$nc
)

# Сохранить объект подбора
saveRDS(fit, file = file.path(output_dir, "fit.rds"))

# ------------------- 5. ДИАГНОСТИКА И ВИЗУАЛИЗАЦИЯ --------------------

# Основные графики (PNG)
png(file.path(fig_dir, "ensemble.png"), width = 1400, height = 900, res = 144); JABBA::jbplot_ensemble(fit); dev.off()
JABBA::jabba_plots(fit, output.dir = fig_dir)

png(file.path(fig_dir, "catch.png"), width = 1400, height = 900, res = 144); JABBA::jbplot_catch(fit); dev.off()
png(file.path(fig_dir, "cpuefits.png"), width = 1400, height = 900, res = 144); JABBA::jbplot_cpuefits(fit); dev.off()
png(file.path(fig_dir, "ppdist.png"), width = 1400, height = 900, res = 144); JABBA::jbplot_ppdist(fit); dev.off()
png(file.path(fig_dir, "residuals.png"), width = 1400, height = 900, res = 144); JABBA::jbplot_residuals(fit); dev.off()
png(file.path(fig_dir, "mcmc.png"), width = 1400, height = 900, res = 144); JABBA::jbplot_mcmc(fit); dev.off()
png(file.path(fig_dir, "traj_B.png"), width = 1400, height = 900, res = 144); JABBA::jbplot_trj(fit, type = "B"); dev.off()
png(file.path(fig_dir, "traj_F.png"), width = 1400, height = 900, res = 144); JABBA::jbplot_trj(fit, type = "F"); dev.off()
png(file.path(fig_dir, "traj_BBmsy.png"), width = 1400, height = 900, res = 144); JABBA::jbplot_trj(fit, type = "BBmsy"); dev.off()
png(file.path(fig_dir, "traj_FFmsy.png"), width = 1400, height = 900, res = 144); JABBA::jbplot_trj(fit, type = "FFmsy"); dev.off()
png(file.path(fig_dir, "spphase.png"), width = 1400, height = 900, res = 144); JABBA::jbplot_spphase(fit); dev.off()
png(file.path(fig_dir, "kobe.png"), width = 1400, height = 900, res = 144); JABBA::jbplot_kobe(fit); dev.off()

# Дополнительные диагностики
png(file.path(fig_dir, "runstest.png"), width = 1400, height = 900, res = 144); JABBA::jbplot_runstest(fit); dev.off()
png(file.path(fig_dir, "logfits.png"), width = 1400, height = 900, res = 144); JABBA::jbplot_logfits(fit); dev.off()
png(file.path(fig_dir, "procdev.png"), width = 1400, height = 900, res = 144); JABBA::jbplot_procdev(fit); dev.off()

# Сохранение временных рядов
utils::write.csv(fit$timeseries, file = file.path(output_dir, "results_timeseries.csv"), row.names = FALSE)

# ------------------- 6. РЕТРОСПЕКТИВА --------------------

peels_seq <- if (isTRUE(FAST_RUN)) 1:2 else 1:5
hc <- JABBA::hindcast_jabba(jbinput = jbinput, fit = fit, peels = peels_seq)

mohnsrho <- JABBA::jbplot_retro(hc, as.png = TRUE, output.dir = retro_dir, xlim = c(2007, 2022))
utils::write.csv(mohnsrho, file = file.path(retro_dir, "mohnsrho.csv"), row.names = FALSE)

mase <- JABBA::jbplot_hcxval(hc, as.png = TRUE, output.dir = retro_dir)
utils::write.csv(mase, file = file.path(retro_dir, "mase.csv"), row.names = FALSE)

# ------------------- 7. ПРОГНОЗ --------------------

fw1 <- JABBA::fw_jabba(
  fit,
  nyears = 10,
  imp.yr = 1,
  imp.values = seq(10, 16, 2),
  quant = "Catch",
  type = "abs",
  stochastic = TRUE
)

# Сохранить прогнозный объект
saveRDS(fw1, file = file.path(output_dir, "forecast_fw1.rds"))

# Визуализация ансамбля прогнозов
png(file.path(fig_dir, "forecast_ensemble.png"), width = 1400, height = 900, res = 144)
JABBA::jbpar(mfrow = c(3, 2))
JABBA::jbplot_ensemble(fw1)
dev.off()

# Доп. график для B/Bmsy
png(file.path(fig_dir, "forecast_BBmsy.png"), width = 1400, height = 900, res = 144)
JABBA::jbplot_ensemble(
  fw1,
  subplots = c(1),
  add = FALSE,
  xlim = c(2020, 2035),
  legend.loc = "topleft"
)
dev.off()

# ------------------- 8. ОБРАБОТКА ПРОГНОЗА --------------------

forecast_data <- subset(
  fw1,
  year %in% 2025:2034 &
    run %in% c("C10", "C12", "C14", "C16") &
    type == "prj"
)

median_B <- aggregate(B ~ year + run, data = forecast_data, FUN = median)
median_Catch <- aggregate(Catch ~ year + run, data = forecast_data, FUN = median)

b_table <- reshape2::dcast(median_B, year ~ run, value.var = "B")
catch_table <- reshape2::dcast(median_Catch, year ~ run, value.var = "Catch")

utils::write.csv(b_table, file.path(output_dir, "biomass_forecast.csv"), row.names = FALSE)
utils::write.csv(catch_table, file.path(output_dir, "catch_forecast.csv"), row.names = FALSE)

# Также выводим краткие результаты в stdout
cat("\nМедианная биомасса (B) по сценарию:\n"); print(b_table)
cat("\nМедианные уловы (Catch) по сценарию:\n"); print(catch_table)

# ------------------- 9. СВОДКА --------------------

cat(sprintf("\nГотово. Результаты сохранены в: %s\n", output_dir))
cat(sprintf("Графики: %s\n", fig_dir))
cat(sprintf("Ретроспектива: %s\n", retro_dir))