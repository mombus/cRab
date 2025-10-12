# ============================================
# Bayesian hierarchical salmon demo (NASCO-style)
# - simulate data
# - JAGS state-space model (log-AR1, Poisson obs)
# - posterior summaries: P(S>=CL) and H* (75% rule)
# - posterior predictive checks and plots
# NOTE: educational example, not for management use
# ============================================

suppressPackageStartupMessages({
  need <- c("rjags","coda","ggplot2")
  inst <- need[!need %in% installed.packages()[,"Package"]]
  if(length(inst)) install.packages(inst, repos = "https://cloud.r-project.org")
  library(rjags)
  library(coda)
  library(ggplot2)
})

# Optional: light dependencies if available, but not required
.has_scales <- requireNamespace("scales", quietly = TRUE)

args <- commandArgs(trailingOnly = TRUE)
fast_mode <- any(grepl("--fast|--quick|--demo", args))

set.seed(42)

# 1) Simulate data -----------------------------------------------------------
R <- 6    # rivers
Tt <- 12  # years

# centered environmental covariate
env <- scale(arima.sim(model = list(ar=0.6), n = Tt), center = TRUE, scale = TRUE)[,1]

# true parameters (for sim)
phi_true   <- 0.6
beta_true  <- 0.25
sigma_proc <- 0.30
mu_alpha   <- log(1200)
sigma_alpha<- 0.50
alpha_r    <- rnorm(R, mu_alpha, sigma_alpha)

# detectability by river
q_r_true <- rbeta(R, shape1 = 20, shape2 = 5)  # around ~0.8

# latent spawners and observations
logS <- matrix(NA_real_, R, Tt)
S    <- matrix(NA_real_, R, Tt)
y    <- matrix(NA_integer_, R, Tt)

for(r in 1:R){
  logS[r,1] <- rnorm(1, mean = alpha_r[r] + beta_true*env[1], sd = sigma_proc)
  S[r,1]    <- exp(logS[r,1])
  y[r,1]    <- rpois(1, lambda = q_r_true[r]*S[r,1])
  for(t in 2:Tt){
    mean_log <- alpha_r[r] + beta_true*env[t] + phi_true*(logS[r,t-1] - alpha_r[r])
    logS[r,t] <- rnorm(1, mean_log, sd = sigma_proc)
    S[r,t]    <- exp(logS[r,t])
    y[r,t]    <- rpois(1, lambda = q_r_true[r]*S[r,t])
  }
}

# Conservation Limits per river as 80% of long-term level
CL <- round(exp(alpha_r) * 0.8)

# 2) JAGS data ---------------------------------------------------------------
data_jags <- list(
  R = R,
  Tt = Tt,
  y = y,
  env = as.numeric(env)
)

# 3) JAGS model (ASCII only) ------------------------------------------------
model_string <- "
model {
  # Hyperparameters
  mu_alpha ~ dnorm(0, 1.0E-4)
  sigma_alpha ~ dunif(0, 5)
  tau_alpha <- pow(sigma_alpha, -2)

  phi ~ dunif(0, 0.99)        # AR(1) on log scale
  beta ~ dnorm(0, 1.0E-4)     # environmental effect
  sigma_proc ~ dunif(0, 2)
  tau_proc <- pow(sigma_proc, -2)

  for(r in 1:R){
    alpha[r] ~ dnorm(mu_alpha, tau_alpha)
    q[r] ~ dbeta(20, 5)       # informative prior on detectability (~0.8)

    # State and observations
    # t=1
    logS[r,1] ~ dnorm(alpha[r] + beta*env[1], tau_proc)
    S[r,1] <- exp(logS[r,1])
    y[r,1] ~ dpois(q[r] * S[r,1])
    y_rep[r,1] ~ dpois(q[r] * S[r,1]) # PPC

    for(t in 2:Tt){
      mean_logS[r,t] <- alpha[r] + beta*env[t] + phi * (logS[r,t-1] - alpha[r])
      logS[r,t] ~ dnorm(mean_logS[r,t], tau_proc)
      S[r,t] <- exp(logS[r,t])
      y[r,t] ~ dpois(q[r] * S[r,t])
      y_rep[r,t] ~ dpois(q[r] * S[r,t]) # PPC
    }
  }
}
"

# 4) Initial values ----------------------------------------------------------
inits_fun <- function(){
  list(
    mu_alpha = log(1000),
    sigma_alpha = runif(1, 0.2, 1.0),
    sigma_proc = runif(1, 0.1, 0.8),
    phi = runif(1, 0.2, 0.9),
    beta = rnorm(1, 0, 0.2),
    alpha = rnorm(R, log(1000), 0.5),
    q = pmin(pmax(rbeta(R, 8, 2), 0.05), 0.98),
    logS = log(pmax(y + 1, 5))
  )
}

# 5) Run MCMC ----------------------------------------------------------------
cat("Compiling JAGS...\n")
jm <- jags.model(
  textConnection(model_string),
  data = data_jags,
  inits = inits_fun,
  n.chains = 3,
  n.adapt = if (fast_mode) 800 else 1500
)

update(jm, n.iter = if (fast_mode) 800 else 2000)

params <- c("mu_alpha","sigma_alpha","phi","beta","sigma_proc","alpha","q","S","y_rep")
mcmc_samp <- coda.samples(jm, variable.names = params, n.iter = if (fast_mode) 1500 else 4000, thin = 2)

# 6) Post-processing ---------------------------------------------------------
cat("Convergence diagnostics (Gelman-Rubin)\n")
print(gelman.diag(mcmc_samp, multivariate = FALSE))

post_mat <- as.matrix(mcmc_samp)

S_colname <- function(r,t) paste0("S[", r, ",", t, "]")
y_colname <- function(r,t) paste0("y[", r, ",", t, "]")
yrep_colname <- function(r,t) paste0("y_rep[", r, ",", t, "]")

# P(S>=CL) by river and year
P_CL <- matrix(NA_real_, R, Tt, dimnames = list(paste0("R",1:R), paste0("Y",1:Tt)))
for(r in 1:R){
  for(t in 1:Tt){
    s_vals <- post_mat[, S_colname(r,t)]
    P_CL[r,t] <- mean(s_vals >= CL[r])
  }
}

# H* (75% rule) for last year
Hstar75 <- numeric(R)
names(Hstar75) <- paste0("R",1:R)
S_last <- vector("list", R)
for(r in 1:R){
  s_vals <- post_mat[, S_colname(r, Tt)]
  surplus <- pmax(0, s_vals - CL[r])
  Hstar75[r] <- as.numeric(quantile(surplus, probs = 0.25))
  S_last[[r]] <- s_vals
}

# 7) Plots -------------------------------------------------------------------
output_dir <- file.path("images", "NASCO_salmon")
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# P(CL) trajectories

df_p <- data.frame(
  river = rep(paste0("R",1:R), each = Tt),
  year  = rep(1:Tt, times = R),
  Pcl   = as.vector(P_CL)
)

p_pcl <- ggplot(df_p, aes(year, Pcl, color = river, group = river)) +
  geom_line() + geom_point() +
  geom_hline(yintercept = 0.75, linetype = 2) +
  scale_y_continuous(limits = c(0,1)) +
  labs(title = "Probability of meeting CL over time",
       x = "Year", y = "P(S >= CL)") +
  theme_minimal()

ggsave(filename = file.path(output_dir, "pcl_timeseries.png"), plot = p_pcl, width = 8, height = 5, dpi = 150)

# S distributions (last year) vs CL

df_s <- do.call(rbind, lapply(1:R, function(r){
  data.frame(river = paste0("R", r), S = S_last[[r]])
}))

df_cl <- data.frame(river = paste0("R",1:R), CL = CL)

p_s <- ggplot(df_s, aes(S, after_stat(density), fill = river, color = river)) +
  geom_histogram(bins = 40, alpha = 0.2, position = "identity") +
  geom_vline(data = df_cl, aes(xintercept = CL, color = river), linetype = 2) +
  {if (.has_scales) scale_x_continuous(labels = scales::comma) else scale_x_continuous()} +
  labs(title = "Posterior S (last year) vs CL",
       x = "S (spawners)", y = "Density") +
  theme_minimal()

ggsave(filename = file.path(output_dir, "s_last_hist_vs_cl.png"), plot = p_s, width = 8, height = 5, dpi = 150)

# Posterior predictive check: observed y vs mean(y_rep)
yrep_mean <- matrix(NA_real_, R, Tt)
for(r in 1:R){
  for(t in 1:Tt){
    yrep_mean[r,t] <- mean(post_mat[, yrep_colname(r,t)])
  }
}

df_ppc <- data.frame(
  river = rep(paste0("R",1:R), each = Tt),
  year  = rep(1:Tt, times = R),
  y_obs = as.vector(y),
  y_rep = as.vector(yrep_mean)
)

p_ppc <- ggplot(df_ppc, aes(y_obs, y_rep, color = river)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  geom_point() +
  labs(title = "Posterior predictive: mean(y_rep) vs observed y", x = "Observed", y = "Posterior mean y_rep") +
  theme_minimal()

ggsave(filename = file.path(output_dir, "ppc_scatter.png"), plot = p_ppc, width = 6, height = 5, dpi = 150)

# 8) Console output ----------------------------------------------------------
cat("\nConservation Limits (CL) by river:\n")
print(setNames(CL, paste0("R",1:R)))

cat("\nP(S>=CL) by river in last year (T=", Tt, "):\n", sep = "")
print(round(P_CL[,Tt], 3))

cat("\nHarvest allowance H* at 75% rule (last year):\n")
print(round(Hstar75))

cat("\nOutputs saved to:", output_dir, "\n")
