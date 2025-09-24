
#remotes::install_github("DeepWaterIMR/ggFishPlots")

library(ggFishPlots)

data(survey_ghl) # example data

head(survey_ghl)

plot_growth(survey_ghl, length = "length", age = "age")

htmlcat <- function(text){
  cat(gsub(pattern = "\n", replacement = "  \n", x = text))
}

htmlcat(plot_growth(survey_ghl)$text)

plot_growth(survey_ghl, split.by.sex = TRUE)$plot

plot_growth(survey_ghl, force.zero.group.length = 14, boxplot = FALSE)$plot

plot_maturity(survey_ghl, length = "length", maturity = "maturity")

plot_maturity(survey_ghl, bootstrap.n = 10)

plot_maturity(survey_ghl, split.by.sex = TRUE)$plot

plot_maturity(survey_ghl, length = "age", length.unit = "years",
              xlab = "Age", length.bin.width = 1, split.by.sex = TRUE)$plot


plot_maturity(survey_ghl, length = "age", length.unit = "years",
              xlab = "Age", length.bin.width = 1, 
              force.zero.group.length = 0,
              force.zero.group.n = c("F" = exp(11.363), "M" = exp(11.885)),
              split.by.sex = TRUE)$plot

plot_lw(survey_ghl, length = "length", weight = "weight")

plot_lw(survey_ghl, use.nls = TRUE)

plot_lw(survey_ghl, split.by.sex = TRUE, correct.units = TRUE)

plot_lw(survey_ghl, split.by.sex = TRUE, log.axes = TRUE)$plot

plot_lw(survey_ghl, outlier.percentile = 99.5, annotate.coefficients = TRUE)$plot







