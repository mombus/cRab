# ============================================
# Улучшенная байесовская иерархическая модель для лосося (R + JAGS)
# - имитация данных
# - JAGS-модель (state-space, лог-AR1, Poisson наблюдения)
# - оценка вероятности выполнения CL и допустимого изъятия
# - расширенная визуализация и диагностика
# ВНИМАНИЕ: учебный упрощённый пример, не для реального менеджмента.
# ============================================

# 0) Пакеты ---------------------------------------------------------------
need <- c("rjags", "coda", "ggplot2", "dplyr", "tidyr", "gridExtra", 
          "viridis", "scales", "bayesplot", "MCMCvis")
inst <- need[!need %in% installed.packages()[,"Package"]]
if(length(inst)) install.packages(inst, repos = "https://cloud.r-project.org")

library(rjags)
library(coda)
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(viridis)
library(scales)
library(bayesplot)
library(MCMCvis)

set.seed(42)

# 1) Имитация данных ------------------------------------------------------
R <- 6    # число рек
Tt <- 12  # число лет

# Экологический ковариат (например, индекс климата), центрированный
env <- scale(arima.sim(model = list(ar=0.6), n = Tt), center = TRUE, scale = TRUE)[,1]

# Истинные параметры (для симуляции)
phi_true   <- 0.6          # автокорреляция (на лог-шкале)
beta_true  <- 0.25         # эффект среды
sigma_proc <- 0.30         # процессная дисперсия (на лог-шкале)
mu_alpha   <- log(1200)    # средний уровень по рекам
sigma_alpha<- 0.50         # межречная вариабельность (на лог-шкале)
alpha_r    <- rnorm(R, mu_alpha, sigma_alpha)

# Детектируемость (вероятность учёта) по рекам
q_r_true <- rbeta(R, shape1 = 20, shape2 = 5)  # вокруг ~0.8

# Генерация латентных спавнеров (S) и наблюдений (y)
logS <- matrix(NA_real_, R, Tt)
S    <- matrix(NA_real_, R, Tt)
y    <- matrix(NA_integer_, R, Tt)

for(r in 1:R){
  # t=1
  logS[r,1] <- rnorm(1, mean = alpha_r[r] + beta_true*env[1], sd = sigma_proc)
  S[r,1]    <- exp(logS[r,1])
  y[r,1]    <- rpois(1, lambda = q_r_true[r]*S[r,1])
  # t>1
  for(t in 2:Tt){
    mean_log <- alpha_r[r] + beta_true*env[t] + phi_true*(logS[r,t-1] - alpha_r[r])
    logS[r,t] <- rnorm(1, mean_log, sd = sigma_proc)
    S[r,t]    <- exp(logS[r,t])
    y[r,1]    <- rpois(1, lambda = q_r_true[r]*S[r,t])
  }
}

# Зададим Conservation Limits по рекам как 80% долгосрочного уровня
CL <- round(exp(alpha_r) * 0.8)

# 2) Подготовка данных для JAGS ------------------------------------------
data_jags <- list(
  R = R,
  Tt = Tt,
  y = y,
  env = as.numeric(env)
)

# Начальные значения для MCMC
inits_fun <- function(){
  list(
    mu_alpha = log(1000),
    sigma_alpha = runif(1, 0.2, 1.0),
    sigma_proc = runif(1, 0.1, 0.8),
    phi = runif(1, 0.2, 0.9),
    beta = rnorm(1, 0, 0.2),
    alpha = rnorm(R, log(1000), 0.5),
    q = pmin(pmax(rbeta(R, 8, 2), 0.05), 0.98),
    logS = log(pmax(y + 1, 5)) # грубая инициализация
  )
}

# 3) Исправленная JAGS-модель ---------------------------------------------
model_string <- "
model {
  # Гиперпараметры
  mu_alpha ~ dnorm(0, 1.0E-4)
  sigma_alpha ~ dunif(0, 5)
  tau_alpha <- pow(sigma_alpha, -2)

  phi ~ dunif(0, 0.99)         # AR(1) на лог-шкале
  beta ~ dnorm(0, 1.0E-4)      # эффект среды
  sigma_proc ~ dunif(0, 2)
  tau_proc <- pow(sigma_proc, -2)

  for(r in 1:R){
    alpha[r] ~ dnorm(mu_alpha, tau_alpha)
    q[r] ~ dbeta(20, 5)        # информативный приор на детектируемость (~0.8)

    # Состояние и наблюдения
    # t=1
    logS[r,1] ~ dnorm(alpha[r] + beta*env[1], tau_proc)
    S[r,1] <- exp(logS[r,1])
    y[r,1] ~ dpois(q[r] * S[r,1])

    for(t in 2:Tt){
      mean_logS[r,t] <- alpha[r] + beta*env[t] + phi * (logS[r,t-1] - alpha[r])
      logS[r,t] ~ dnorm(mean_logS[r,t], tau_proc)
      S[r,t] <- exp(logS[r,t])
      y[r,t] ~ dpois(q[r] * S[r,t])
    }
  }
}
"

# 4) Запуск MCMC ----------------------------------------------------------
cat("Компиляция JAGS...\n")
jm <- jags.model(
  textConnection(model_string),
  data = data_jags,
  inits = inits_fun,
  n.chains = 3,
  n.adapt = 1500
)

cat("Burn-in...\n")
update(jm, n.iter = 2000)  # burn-in

params <- c("mu_alpha", "sigma_alpha", "phi", "beta", "sigma_proc", "alpha", "q", "S")
mcmc_samp <- coda.samples(jm, variable.names = params, n.iter = 4000, thin = 2)

cat("Диагностика сходимости (Gelman-Rubin)...\n")
gelman_diag <- gelman.diag(mcmc_samp, multivariate = FALSE)
print(gelman_diag)

# 5) Постобработка: извлечь постериоры S ----------------------------------
# Комбинируем цепи в одну матрицу (iters x params)
post_mat <- as.matrix(mcmc_samp)

# Вспомогательная функция для выборки столбцов S[r,t]
S_colname <- function(r,t) paste0("S[", r, ",", t, "]")

# Вероятность выполнения CL по рекам и годам
P_CL <- matrix(NA_real_, R, Tt, dimnames = list(paste0("R",1:R), paste0("Y",1:Tt)))
for(r in 1:R){
  for(t in 1:Tt){
    col <- S_colname(r,t)
    s_vals <- post_mat[, col]
    P_CL[r,t] <- mean(s_vals >= CL[r])
  }
}

# Допустимое изъятие по правилу 75% (для последнего года)
Hstar75 <- numeric(R)
names(Hstar75) <- paste0("R",1:R)
S_last <- vector("list", R)
for(r in 1:R){
  s_vals <- post_mat[, S_colname(r, Tt)]
  surplus <- pmax(0, s_vals - CL[r])
  Hstar75[r] <- as.numeric(quantile(surplus, probs = 0.25))  # 25-й перцентиль => 75% шанс >= CL
  S_last[[r]] <- s_vals
}

# 6) Вывод результатов -----------------------------------------------------
cat("\nConservation Limits (CL) по рекам:\n")
print(setNames(CL, paste0("R",1:R)))

cat("\nВероятность выполнения CL по рекам в последнем году (T = ", Tt, "):\n", sep = "")
print(round(P_CL[,Tt], 3))

cat("\nОценка допустимого изъятия на реку (H* при 75% вероятности не нарушить CL), последний год:\n")
print(round(Hstar75))

# 7) Расширенная визуализация ---------------------------------------------

# 7.1) Траектории P(CL) по годам
df_p <- data.frame(
  river = rep(paste0("R",1:R), each = Tt),
  year  = rep(1:Tt, times = R),
  Pcl   = as.vector(P_CL)
)

p1 <- ggplot(df_p, aes(year, Pcl, color = river)) +
  geom_line(size = 1.2) + 
  geom_point(size = 2) +
  geom_hline(yintercept = 0.75, linetype = 2, color = "red", size = 1) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_color_viridis_d() +
  labs(title = "Вероятность выполнения CL по годам",
       x = "Год", y = "P(S >= CL)", color = "Река") +
  theme_minimal() +
  theme(legend.position = "bottom")

# 7.2) Распределения S в последнем году vs CL
df_s <- do.call(rbind, lapply(1:R, function(r){
  data.frame(river = paste0("R", r), S = S_last[[r]])
}))
df_cl <- data.frame(river = paste0("R",1:R), CL = CL)

p2 <- ggplot(df_s, aes(S, after_stat(density), fill = river, color = river)) +
  geom_histogram(bins = 40, alpha = 0.3, position = "identity") +
  geom_vline(data = df_cl, aes(xintercept = CL, color = river), 
             linetype = 2, size = 1.2) +
  scale_x_continuous(labels = scales::comma) +
  scale_fill_viridis_d() +
  scale_color_viridis_d() +
  labs(title = "Постериорные распределения S (последний год) и линии CL",
       x = "S (спавнеры)", y = "Плотность", fill = "Река", color = "Река") +
  theme_minimal() +
  theme(legend.position = "bottom")

# 7.3) Траектории наблюдений vs истинных значений
df_obs <- data.frame(
  river = rep(paste0("R",1:R), each = Tt),
  year = rep(1:Tt, times = R),
  observed = as.vector(y),
  true = as.vector(S),
  env = rep(env, times = R)
)

p3 <- ggplot(df_obs, aes(year)) +
  geom_line(aes(y = true, color = "Истинное S"), size = 1.2) +
  geom_point(aes(y = observed, color = "Наблюдения"), size = 2) +
  geom_line(aes(y = env*100 + 500, color = "Климат (масштаб)"), size = 1, alpha = 0.7) +
  facet_wrap(~river, scales = "free_y") +
  scale_color_manual(values = c("Истинное S" = "blue", "Наблюдения" = "red", 
                               "Климат (масштаб)" = "green")) +
  labs(title = "Траектории спавнеров: истинные vs наблюдения",
       x = "Год", y = "Количество", color = "Тип") +
  theme_minimal() +
  theme(legend.position = "bottom")

# 7.4) Диагностика MCMC
p4 <- mcmc_trace(mcmc_samp, pars = c("mu_alpha", "phi", "beta", "sigma_proc"))
p5 <- mcmc_dens(mcmc_samp, pars = c("mu_alpha", "phi", "beta", "sigma_proc"))

# 7.5) Матрица корреляций параметров
param_names <- c("mu_alpha", "phi", "beta", "sigma_proc", "sigma_alpha")
param_samples <- post_mat[, param_names]
cor_matrix <- cor(param_samples)

df_cor <- expand.grid(Var1 = param_names, Var2 = param_names)
df_cor$value <- as.vector(cor_matrix)

p6 <- ggplot(df_cor, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  geom_text(aes(label = round(value, 2)), color = "white", size = 3) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                       midpoint = 0, limits = c(-1, 1)) +
  labs(title = "Корреляционная матрица параметров",
       x = "", y = "", fill = "Корреляция") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 7.6) Управленческие рекомендации
df_management <- data.frame(
  river = paste0("R", 1:R),
  CL = CL,
  P_CL_last = P_CL[, Tt],
  Hstar75 = Hstar75,
  status = ifelse(P_CL[, Tt] > 0.75, "Безопасно", "Осторожно")
)

p7 <- ggplot(df_management, aes(CL, P_CL_last, size = Hstar75, color = status)) +
  geom_point(alpha = 0.7) +
  geom_hline(yintercept = 0.75, linetype = 2, color = "red") +
  geom_vline(xintercept = mean(CL), linetype = 2, color = "blue", alpha = 0.5) +
  scale_size_continuous(range = c(3, 10), name = "Допустимое\nизъятие") +
  scale_color_manual(values = c("Безопасно" = "green", "Осторожно" = "orange")) +
  labs(title = "Управленческие рекомендации по рекам",
       x = "Conservation Limit", y = "P(выполнения CL)", 
       color = "Статус") +
  theme_minimal() +
  theme(legend.position = "right")

# Объединяем все графики
grid.arrange(p1, p2, p3, p4, p5, p6, p7, ncol = 2)

# 8) Дополнительная диагностика и анализ ----------------------------------

# 8.1) Эффективные размеры выборки
cat("\nЭффективные размеры выборки:\n")
print(effectiveSize(mcmc_samp))

# 8.2) Автокорреляция
cat("\nАвтокорреляция (первые 10 лагов):\n")
print(autocorr(mcmc_samp, lags = c(1, 5, 10)))

# 8.3) Сводная статистика по параметрам
cat("\nСводная статистика по параметрам:\n")
print(summary(mcmc_samp))

# 8.4) Анализ чувствительности к приорам
cat("\nАнализ чувствительности: сравнение с истинными значениями\n")
true_params <- data.frame(
  parameter = c("mu_alpha", "phi", "beta", "sigma_proc", "sigma_alpha"),
  true_value = c(mu_alpha, phi_true, beta_true, sigma_proc, sigma_alpha)
)

post_summary <- summary(mcmc_samp)$statistics[true_params$parameter, c("Mean", "2.5%", "97.5%")]
comparison <- cbind(true_params, post_summary)
print(comparison)

# 9) Краткое резюме интерпретации -----------------------------------------
cat("
Интерпретация:
- P_CL[r,T] — вероятность, что в реке r в последний год спавнеров достаточно для выполнения CL.
- Hstar75[r] — оценка изъятия, которое можно взять (в сумме по реке) и всё ещё иметь ≥75% шанса выполнить CL
  (в учебной постановке это 'биологический запас под добычу' без явного разложения по смешанным промыслам).
- Модель иерархична: уровень по реке alpha[r] тянется к общему среднему mu_alpha; во времени лог-численность следует AR(1);
  наблюдения имеют детектируемость q[r]. Неопределённость пропагируется в решения.

Замечание:
- В реальных оценках ICES/NASCO добавляют уловы по районам, смертности выпуска, преднерестовые потери, возрастную структуру,
  генетическое разложение смешанных уловов и проводят MSE. Здесь — минимальный демонстрационный каркас.
")