# ------------------------------------------------------------
# CPUE standardization by YEAR controlling for MONTH, CALL, REGION
# Methods kept: GLM, GAM, GAMM (no ML)
# Outputs: diagnostics, effect plots, standardized indices with CIs
# ------------------------------------------------------------

# Set working directory as requested
setwd("C:/GLM/")

suppressPackageStartupMessages({
	library(tidyverse)
	library(readxl)
	library(mgcv)
	library(gamm4)
	library(emmeans)
	library(broom)
	library(broom.mixed)
	library(performance)
	library(DHARMa)
})

# ---------------------------- Paths ----------------------------
DATA_PATH <- file.path(getwd(), "KARTOGRAPHIC.xlsx")
OUTPUT_DIR <- file.path(getwd(), "output")
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)
PLOT_DIR <- file.path(OUTPUT_DIR, "plots")
if (!dir.exists(PLOT_DIR)) dir.create(PLOT_DIR, recursive = TRUE)
TABLE_DIR <- file.path(OUTPUT_DIR, "tables")
if (!dir.exists(TABLE_DIR)) dir.create(TABLE_DIR, recursive = TRUE)

set.seed(42)

theme_set(ggplot2::theme_minimal(base_size = 12))

# ---------------------------- Helpers ----------------------------

save_plot <- function(plot_obj, filename, width = 9, height = 6) {
	ggplot2::ggsave(filename = file.path(PLOT_DIR, filename), plot = plot_obj, width = width, height = height, dpi = 300)
}

scale_to_index <- function(values, method = c("mean", "first")) {
	method <- match.arg(method)
	if (method == "mean") {
		return(as.numeric(values) / mean(as.numeric(values), na.rm = TRUE))
	}
	if (method == "first") {
		return(as.numeric(values) / as.numeric(values[1]))
	}
}

# Average predictions over empirical covariate distribution (for GAMM)
compute_standardized_index <- function(model, base_data, year_levels, predict_fun, response_transform = identity, n_boot = 200L, seed = 42L) {
	set.seed(seed)
	results <- vector("list", length(year_levels))
	for (i in seq_along(year_levels)) {
		year_value <- year_levels[i]
		newdata <- base_data
		newdata$YEAR <- factor(year_value, levels = levels(base_data$YEAR))
		preds <- suppressWarnings(predict_fun(model, newdata))
		mu <- mean(response_transform(preds), na.rm = TRUE)
		if (n_boot > 0) {
			boot_vals <- replicate(n_boot, {
				idx <- sample.int(n = nrow(base_data), size = nrow(base_data), replace = TRUE)
				bd <- newdata[idx, , drop = FALSE]
				p <- suppressWarnings(predict_fun(model, bd))
				mean(response_transform(p), na.rm = TRUE)
			})
			ci <- stats::quantile(boot_vals, probs = c(0.025, 0.975), na.rm = TRUE)
			res <- tibble::tibble(YEAR = year_value, value = mu, lcl = ci[[1]], ucl = ci[[2]])
		} else {
			res <- tibble::tibble(YEAR = year_value, value = mu, lcl = NA_real_, ucl = NA_real_)
		}
		results[[i]] <- res
	}
	dplyr::bind_rows(results)
}

# Emmeans-based standardized means on response scale (GLM/GAM)
emmeans_standardized_index <- function(model, variable = "YEAR") {
	out <- suppressWarnings(emmeans::emmeans(model, specs = stats::as.formula(paste0("~ ", variable)), type = "response"))
	df <- tibble::as_tibble(out) %>% dplyr::select(!!rlang::sym(variable), response = response, lower.CL, upper.CL)
	colnames(df) <- c("YEAR", "value", "lcl", "ucl")
	df
}

# ---------------------------- 1) Load & preprocess ----------------------------

if (!file.exists(DATA_PATH)) {
	stop(sprintf("Data file not found at %s", DATA_PATH))
}

DATA <- readxl::read_excel(DATA_PATH, sheet = "FISHERY") %>%
	tibble::as_tibble() %>%
	dplyr::filter(MONTH > 8 & MONTH < 12)

DATA <- DATA %>%
	dplyr::mutate(
		YEAR = as.factor(YEAR),
		MONTH = as.factor(MONTH),
		CALL = as.factor(CALL),
		REGION = as.factor(REGION),
		VESSELNUMBER = as.factor(VESSELNUMBER),
		CPUE = as.numeric(CPUE)
	) %>%
	dplyr::filter(!is.na(CPUE))

if (any(DATA$CPUE <= 0, na.rm = TRUE)) {
	min_pos <- min(DATA$CPUE[DATA$CPUE > 0], na.rm = TRUE)
	offset <- min_pos / 2
	DATA <- DATA %>% dplyr::mutate(CPUE_POS = dplyr::if_else(CPUE <= 0, CPUE + offset, CPUE))
} else {
	DATA <- DATA %>% dplyr::mutate(CPUE_POS = CPUE)
}

# EDA
p_eda <- DATA %>%
	ggplot2::ggplot(ggplot2::aes(x = YEAR, y = CPUE)) +
	ggplot2::geom_boxplot(outlier.alpha = 0.2) +
	ggplot2::labs(title = "CPUE by YEAR", x = "YEAR", y = "CPUE")
save_plot(p_eda, "eda_cpue_by_year.png")

# ---------------------------- 2) GLM (Gamma, log-link) ----------------------------

glm_gamma_fit <- stats::glm(
	formula = CPUE_POS ~ YEAR + MONTH + CALL + REGION,
	family = stats::Gamma(link = "log"),
	data = DATA
)

png(file.path(PLOT_DIR, "glm_gamma_residuals.png"), width = 1200, height = 800, res = 150)
par(mfrow = c(2, 2))
plot(glm_gamma_fit)
dev.off()

sim_glm <- suppressWarnings(DHARMa::simulateResiduals(fittedModel = glm_gamma_fit, n = 1000))
png(file.path(PLOT_DIR, "glm_gamma_DHARMa.png"), width = 1200, height = 400, res = 150)
par(mfrow = c(1, 3))
plot(sim_glm)
dev.off()

idx_glm_gamma <- emmeans_standardized_index(glm_gamma_fit) %>%
	dplyr::mutate(model = "GLM_Gamma",
		index_mean = scale_to_index(value, method = "mean"),
		index_first = scale_to_index(value, method = "first"))

# ---------------------------- 3) GAM (mgcv) ----------------------------

gam_fit <- mgcv::gam(
	formula = CPUE_POS ~ YEAR + MONTH + CALL + REGION,
	family = stats::Gamma(link = "log"),
	method = "REML",
	data = DATA
)

png(file.path(PLOT_DIR, "gam_diagnostics.png"), width = 1200, height = 800, res = 150)
par(mfrow = c(2, 2))
if (length(gam_fit$smooth) > 0) {
	plot(gam_fit, shade = TRUE, pages = 1)
}
mgcv::gam.check(gam_fit)
dev.off()

idx_gam <- emmeans_standardized_index(gam_fit) %>%
	dplyr::mutate(model = "GAM",
		index_mean = scale_to_index(value, method = "mean"),
		index_first = scale_to_index(value, method = "first"))

# ---------------------------- 4) GAMM (random vessel effect) ----------------------------

gamm_fit <- gamm4::gamm4(
	formula = CPUE_POS ~ YEAR + MONTH + CALL + REGION,
	random = ~(1 | VESSELNUMBER),
	family = stats::Gamma(link = "log"),
	data = DATA
)

gamm_mer <- gamm_fit$mer

sim_gamm <- suppressWarnings(DHARMa::simulateResiduals(fittedModel = gamm_mer, n = 1000))
png(file.path(PLOT_DIR, "gamm_DHARMa.png"), width = 1200, height = 400, res = 150)
par(mfrow = c(1, 3))
plot(sim_gamm)
dev.off()

predict_fun_gamm <- function(m, newdata) {
	predict(m, newdata = newdata, type = "response", re.form = NA)
}

idx_gamm <- compute_standardized_index(
	model = gamm_mer,
	base_data = DATA,
	year_levels = levels(DATA$YEAR),
	predict_fun = predict_fun_gamm,
	response_transform = identity,
	n_boot = 200,
	seed = 7
) %>% dplyr::mutate(
	model = "GAMM",
	index_mean = scale_to_index(value, method = "mean"),
	index_first = scale_to_index(value, method = "first")
)

# ---------------------------- 5) Combine & save indices ----------------------------

indices_all <- dplyr::bind_rows(idx_glm_gamma, idx_gam, idx_gamm) %>%
	dplyr::mutate(YEAR = factor(YEAR, levels = levels(DATA$YEAR)))

readr::write_csv(indices_all, file.path(TABLE_DIR, "standardized_indices_all_models.csv"))

p_abs <- indices_all %>%
	ggplot2::ggplot(ggplot2::aes(x = YEAR, y = value, color = model, group = model)) +
	ggplot2::geom_line(ggplot2::aes(group = model), position = ggplot2::position_dodge(width = 0.2)) +
	ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.2)) +
	ggplot2::geom_errorbar(ggplot2::aes(ymin = lcl, ymax = ucl), width = 0.1, alpha = 0.5, position = ggplot2::position_dodge(width = 0.2)) +
	ggplot2::labs(title = "Standardized CPUE by YEAR (absolute, response scale)", y = "Predicted CPUE", x = "YEAR") +
	ggplot2::theme(legend.position = "bottom")
save_plot(p_abs, "indices_absolute.png", width = 10, height = 6)

p_idx_mean <- indices_all %>%
	ggplot2::ggplot(ggplot2::aes(x = YEAR, y = index_mean, color = model, group = model)) +
	ggplot2::geom_line(ggplot2::aes(group = model), position = ggplot2::position_dodge(width = 0.2)) +
	ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.2)) +
	ggplot2::labs(title = "Standardized CPUE index (scaled to mean = 1)", y = "Index", x = "YEAR") +
	ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
	ggplot2::theme(legend.position = "bottom")
save_plot(p_idx_mean, "indices_scaled_mean1.png", width = 10, height = 6)

p_idx_first <- indices_all %>%
	ggplot2::ggplot(ggplot2::aes(x = YEAR, y = index_first, color = model, group = model)) +
	ggplot2::geom_line(ggplot2::aes(group = model), position = ggplot2::position_dodge(width = 0.2)) +
	ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.2)) +
	ggplot2::labs(title = "Standardized CPUE index (scaled to first year = 1)", y = "Index", x = "YEAR") +
	ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
	ggplot2::theme(legend.position = "bottom")
save_plot(p_idx_first, "indices_scaled_first1.png", width = 10, height = 6)

# ---------------------------- 6) Model summaries ----------------------------

readr::write_csv(broom::tidy(glm_gamma_fit), file.path(TABLE_DIR, "glm_gamma_tidy.csv"))
readr::write_csv(broom::glance(glm_gamma_fit), file.path(TABLE_DIR, "glm_gamma_glance.csv"))

suppressWarnings(try(readr::write_csv(broom.mixed::tidy(gam_fit), file.path(TABLE_DIR, "gam_tidy.csv")), silent = TRUE))
suppressWarnings(try(readr::write_csv(broom.mixed::glance(gam_fit), file.path(TABLE_DIR, "gam_glance.csv")), silent = TRUE))

capture.output(summary(gamm_mer), file = file.path(TABLE_DIR, "gamm_mer_summary.txt"))
suppressWarnings(try(capture.output(summary(gamm_fit$gam), file = file.path(TABLE_DIR, "gamm_gam_summary.txt")), silent = TRUE))

# ---------------------------- 7) Session info ----------------------------

writeLines(capture.output(sessionInfo()), con = file.path(OUTPUT_DIR, "sessionInfo.txt"))

message("Completed CPUE standardization (GLM/GAM/GAMM). Outputs: ", OUTPUT_DIR)