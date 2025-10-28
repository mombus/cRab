# ------------------------------------------------------------
# Практическое занятие: Случайное блуждание vs AR(1)
# Цель: визуализировать, как 10 "мячиков" дрейфуют по течению
# ------------------------------------------------------------

# Параметры симуляции
n_balls  <- 10        # количество мячиков
n_steps  <- 50        # количество шагов (время)
sigma    <- 0.5       # стандартное отклонение шума
phi_ar1  <- 0.8       # коэффициент авторегрессии (|phi| < 1)

set.seed(2025)        # для воспроизводимости

# --- 1. Случайное блуждание (RW): x[t] = x[t-1] + ε[t] ---
rw <- matrix(0, nrow = n_steps, ncol = n_balls)
for (i in 1:n_balls) {
  for (t in 2:n_steps) {
    rw[t, i] <- rw[t - 1, i] + rnorm(1, mean = 0, sd = sigma)
  }
}

# --- 2. AR(1): x[t] = phi * x[t-1] + ε[t] ---
ar1 <- matrix(0, nrow = n_steps, ncol = n_balls)
for (i in 1:n_balls) {
  for (t in 2:n_steps) {
    ar1[t, i] <- phi_ar1 * ar1[t - 1, i] + rnorm(1, mean = 0, sd = sigma)
  }
}

# --- Визуализация ---
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))

# RW
matplot(rw, type = "l", lty = 1, col = rgb(0.8, 0.2, 0.2, 0.6),
        main = "Случайное блуждание (RW)", xlab = "Шаг", ylab = "Положение")
abline(h = 0, lty = 2, col = "gray50")

# AR(1)
matplot(ar1, type = "l", lty = 1, col = rgb(0.2, 0.4, 0.8, 0.6),
        main = "AR(1), φ = 0.8", xlab = "Шаг", ylab = "Положение")
abline(h = 0, lty = 2, col = "gray50")

# Сброс компоновки
par(mfrow = c(1, 1))