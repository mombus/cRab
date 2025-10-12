# ------------------------------------------------------------
# Иерархическая байесовская модель (HBM) для оценки R0 лосося
# с использованием пакета nimble (ИСПРАВЛЕННАЯ ВЕРСИЯ)
# ------------------------------------------------------------

# 1. Установка и подключение пакетов
if (!require("nimble")) install.packages("nimble")
if (!require("tidyverse")) install.packages("tidyverse")
if (!require("bayesplot")) install.packages("bayesplot")
if (!require("coda")) install.packages("coda")

library(nimble)
library(tidyverse)
library(bayesplot)
library(coda)

# 2. Генерация синтетических данных
set.seed(2025)

n_rivers <- 10
n_years <- 15
years <- 2009:2023

# Истинные гиперпараметры
mu_logR0_true <- log(1000)     # ≈ среднее R0 = 1000
sigma_logR0_true <- 0.8
tau_obs_true <- 0.3

# Истинные R0 для каждой реки
logR0_true <- rnorm(n_rivers, mu_logR0_true, sigma_logR0_true)
R0_true <- exp(logR0_true)

# Наблюдения: лог-нормальный шум
obs_data <- tibble()
for (i in 1:n_rivers) {
  log_obs <- rnorm(n_years, mean = logR0_true[i], sd = tau_obs_true)
  obs_data <- bind_rows(obs_data,
                        tibble(
                          river_id = i,
                          year = years,
                          log_recruits = log_obs
                        ))
}

# 3. Визуализация данных
p1 <- obs_data %>%
  mutate(river = paste0("River_", river_id)) %>%
  ggplot(aes(x = year, y = exp(log_recruits), group = river)) +
  geom_line(alpha = 0.7) +
  geom_point(alpha = 0.7) +
  scale_y_log10() +
  labs(title = "Синтетические данные пополнения (log-шкала)",
       y = "Пополнение (особи)", x = "Год") +
  theme_minimal()

print(p1)

# 4. Определение модели в Nimble
hbm_code <- nimbleCode({
  # Гиперпараметры
  mu ~ dnorm(7, sd = 10)
  sigma ~ dunif(0.1, 3)
  tau_obs ~ dunif(0.01, 1)

  # Иерархия по рекам
  for (i in 1:N_rivers) {
    logR0[i] ~ dnorm(mu, sd = sigma)
  }

  # Наблюдения
  for (j in 1:N_obs) {
    log_recruits[j] ~ dnorm(logR0[river_id[j]], sd = tau_obs)
  }
})

# 5. Подготовка данных, констант и начальных значений
constants <- list(
  N_rivers = n_rivers,
  N_obs = nrow(obs_data)
)

data_list <- list(
  log_recruits = obs_data$log_recruits,
  river_id = obs_data$river_id
)

inits_list <- list(
  mu = rnorm(1, 7, 1),
  sigma = runif(1, 0.5, 1.5),
  tau_obs = runif(1, 0.1, 0.5),
  logR0 = rnorm(n_rivers, 7, 1)
)

params_to_monitor <- c("mu", "sigma", "tau_obs", "logR0")

# 6. Создание и компиляция модели
hbm_model <- nimbleModel(code = hbm_code,
                         constants = constants,
                         data = data_list,
                         inits = inits_list,
                         check = FALSE)

hbm_mcmc <- buildMCMC(hbm_model, monitors = params_to_monitor)

compiled_model <- compileNimble(hbm_model)
compiled_mcmc  <- compileNimble(hbm_mcmc, project = hbm_model)

# 7. Запуск MCMC
set.seed(123)
samples <- runMCMC(compiled_mcmc,
                   niter = 12000,
                   nburnin = 2000,
                   nchains = 3,
                   thin = 10,
                   samplesAsCodaMCMC = TRUE)

# 8. Диагностика сходимости
cat("Диагностика сходимости:\n")
print(gelman.diag(samples))

cat("\nЭффективные размеры выборки:\n")
print(effectiveSize(samples))

# Трассировки
p2 <- mcmc_trace(samples, pars = c("mu", "sigma", "tau_obs"))
print(p2)

# 9. Визуализация априоров и постериоров (ИСПРАВЛЕНО)
# Генерируем априорные распределения
prior_mu <- rnorm(5000, 7, 10)
prior_sigma <- runif(5000, 0.1, 3)
prior_tau <- runif(5000, 0.01, 1)

# Создаем данные для априорных распределений
prior_mu_dens <- density(prior_mu)
prior_sigma_dens <- density(prior_sigma)
prior_tau_dens <- density(prior_tau)

# График для mu
p3 <- mcmc_dens(samples, pars = "mu") +
  geom_vline(xintercept = mu_logR0_true, color = "red", linetype = "dashed", size = 1) +
  geom_line(data = tibble(x = prior_mu_dens$x, y = prior_mu_dens$y),
            aes(x = x, y = y), color = "gray50", linetype = "dotted", size = 1) +
  labs(title = "μ: априор vs постериор", 
       subtitle = "Красная линия — истинное значение",
       x = "log(R0)", y = "Плотность") +
  theme_minimal()

print(p3)

# График для sigma
p4 <- mcmc_dens(samples, pars = "sigma") +
  geom_vline(xintercept = sigma_logR0_true, color = "red", linetype = "dashed", size = 1) +
  geom_line(data = tibble(x = prior_sigma_dens$x, y = prior_sigma_dens$y),
            aes(x = x, y = y), color = "gray50", linetype = "dotted", size = 1) +
  labs(title = "σ: межпопуляционная SD", 
       subtitle = "Красная линия — истинное значение",
       x = "σ", y = "Плотность") +
  theme_minimal()

print(p4)

# График для tau_obs
p5 <- mcmc_dens(samples, pars = "tau_obs") +
  geom_vline(xintercept = tau_obs_true, color = "red", linetype = "dashed", size = 1) +
  geom_line(data = tibble(x = prior_tau_dens$x, y = prior_tau_dens$y),
            aes(x = x, y = y), color = "gray50", linetype = "dotted", size = 1) +
  labs(title = "τ: наблюдательная SD", 
       subtitle = "Красная линия — истинное значение",
       x = "τ", y = "Плотность") +
  theme_minimal()

print(p5)

# 10. Оценки R0 по рекам (ИСПРАВЛЕНО)
# Извлекаем постериорные выборки logR0 как числовую матрицу
logR0_samples <- as.matrix(samples[, paste0("logR0[", 1:n_rivers, "]")])
R0_est_median <- apply(exp(logR0_samples), 2, median)
R0_est_q025  <- apply(exp(logR0_samples), 2, quantile, 0.025)
R0_est_q975  <- apply(exp(logR0_samples), 2, quantile, 0.975)

river_summary <- tibble(
  river = paste0("River_", 1:n_rivers),
  R0_true = R0_true,
  R0_est = R0_est_median,
  lo = R0_est_q025,
  hi = R0_est_q975
)

p6 <- ggplot(river_summary, aes(x = reorder(river, R0_true))) +
  geom_point(aes(y = R0_true), color = "red", size = 2, shape = 4) +
  geom_pointrange(aes(y = R0_est, ymin = lo, ymax = hi), color = "steelblue") +
  scale_y_log10() +
  coord_flip() +
  labs(title = "Оценка R0: истинные vs. постериорные",
       y = "R0 (начальное пополнение)", x = "Река") +
  theme_minimal()

print(p6)

# 11. Дополнительная диагностика и сводка
cat("\n=== СВОДКА РЕЗУЛЬТАТОВ ===\n")

# Сводка по гиперпараметрам
cat("\nГиперпараметры:\n")
hyper_summary <- summary(samples[, c("mu", "sigma", "tau_obs")])
print(hyper_summary)

# Сравнение истинных и оцененных значений
cat("\nСравнение истинных и оцененных значений R0:\n")
comparison <- river_summary %>%
  mutate(
    error = abs(R0_est - R0_true),
    relative_error = error / R0_true * 100,
    coverage = (R0_true >= lo & R0_true <= hi)
  )

print(comparison)

cat("\nСтатистики точности:\n")
cat("Средняя абсолютная ошибка:", round(mean(comparison$error), 2), "\n")
cat("Средняя относительная ошибка:", round(mean(comparison$relative_error), 2), "%\n")
cat("Покрытие 95% ДИ:", round(mean(comparison$coverage) * 100, 1), "%\n")

# График точности оценок
p7 <- ggplot(comparison, aes(x = R0_true, y = R0_est)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_errorbar(aes(ymin = lo, ymax = hi), alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  scale_x_log10() +
  scale_y_log10() +
  labs(title = "Точность оценок R0",
       x = "Истинное R0", y = "Оцененное R0") +
  theme_minimal()

print(p7)

cat("\nАнализ завершен успешно!\n")