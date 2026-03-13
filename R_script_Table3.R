#############################################
# Challenges in addressing causal questions #
# using administrative data                 #
#                                           #
# REPRODUCING TABLE 3                       #
# Short term effect: S_t on Y_{t+1}         #
#                                           #
# Outputs (as in the paper):                #
#  - Rate Ratio (%) = 100 * RR              #
#  - Rate Difference (RD) on rate scale     #
# CIs: nonparametric bootstrap (B = 100)    #
#############################################
install.packages("tidyverse")
install.packages("haven")

library(tidyverse)
library(haven)

# setwed("")  #setting working directory


dat_long <- read_dta("data/Simulation_timevar_SEND_Study.dta")

rr_pct <- function(rr) 100 * rr

# Bootstrap function
boot_est <- function(data, B = 100, seed = 1212, est_fun) {
  set.seed(seed)
  n <- nrow(data)
  
  point <- tryCatch(est_fun(data), error = function(e) rep(NA_real_, 4))
  
  boot_mat <- replicate(B, {
    idx <- sample.int(n, n, replace = TRUE)
    tryCatch(est_fun(data[idx, , drop = FALSE]), error = function(e) rep(NA_real_, 4))
  }) %>% t()
  
  ci <- apply(boot_mat, 2, function(z) {
    if (all(is.na(z))) c(NA_real_, NA_real_)
    else as.numeric(quantile(z, c(0.025, 0.975), na.rm = TRUE))
  }) %>% t()
  
  list(point = point, ci = ci)
}

# exposure at time j, outcome at time j+1
make_pair <- function(j) {
  i <- j + 1
  
  dat_long %>%
    filter(t %in% c(j, i)) %>%
    arrange(id, t) %>%
    group_by(id) %>%
    summarise(
      male      = first(male),
      white     = first(white),
      idacibin  = first(idacibin),
      eyfsp_bin = first(eyfsp_bin),
      
      A  = SEN[t == j][1],
      H  = Hosp[t == j][1],
      Y2 = Y[t == i][1],
      D2 = den[t == i][1],
      
      Y0 = `Y0_`[t == i][1],
      Y1 = `Y1_`[t == i][1],
      .groups = "drop"
    ) %>%
    mutate(
      rate  = Y2 / D2,
      inter = as.integer(idacibin == 1 & eyfsp_bin == 0)
    ) %>%
    mutate(
      across(c(male, white, idacibin, eyfsp_bin, A, H, inter), ~ factor(.x, levels = c(0, 1)))
    )
}

# G-computation estimator
gcomp_fun <- function(d) {
  meanD     <- mean(d$D2)
  meanD_att <- mean(d$D2[d$A == 1])
  
  fit <- glm(
    Y2 ~ A * male * white * idacibin * eyfsp_bin * H,
    family = poisson(),
    data = d,
    offset = log(D2)
  )
  
  # ATE
  nd1 <- d; nd1$A <- factor(1, levels = c(0, 1))
  nd0 <- d; nd0$A <- factor(0, levels = c(0, 1))
  
  mu1 <- predict(fit, newdata = nd1, type = "response")
  mu0 <- predict(fit, newdata = nd0, type = "response")
  
  mu1_bar <- mean(mu1)
  mu0_bar <- mean(mu0)
  
  rr_ate <- mu1_bar / mu0_bar
  rd_ate <- (mu1_bar / meanD) - (mu0_bar / meanD)
  
  # ATT
  dt  <- d %>% filter(A == 1)
  nd1t <- dt; nd1t$A <- factor(1, levels = c(0, 1))
  nd0t <- dt; nd0t$A <- factor(0, levels = c(0, 1))
  
  mu1t <- predict(fit, newdata = nd1t, type = "response")
  mu0t <- predict(fit, newdata = nd0t, type = "response")
  
  mu1t_bar <- mean(mu1t)
  mu0t_bar <- mean(mu0t)
  
  rr_att <- mu1t_bar / mu0t_bar
  rd_att <- (mu1t_bar / meanD_att) - (mu0t_bar / meanD_att)
  
  c(ATE_RR = rr_pct(rr_ate), ATT_RR = rr_pct(rr_att),
    ATE_RD = rd_ate,        ATT_RD = rd_att)
}

# True values
true_vals <- function(effect_name) {
  if (effect_name == "S1 on Y2") return(c(21.5, 12.7, -0.046, -0.077))
  if (effect_name == "S2 on Y3") return(c(17.3, 14.4, -0.075, -0.111))
  if (effect_name == "S3 on Y4") return(c(14.3, 13.5, -0.089, -0.124))
  stop("Unknown effect_name: ", effect_name)
}

# ----- Building Table 3 -----
rows <- list()

for (j in 1:3) {
  d <- make_pair(j)
  eff_name <- paste0("S", j, " on Y", j + 1)
  
  tv <- true_vals(eff_name)
  
  gc <- boot_est(d, B = 100, seed = 1212, est_fun = gcomp_fun)
  gc_pt <- gc$point
  gc_ci <- gc$ci
  
  rows[[length(rows) + 1]] <- tibble(
    Effect = eff_name,
    Method = "True value",
    ATE_RR = sprintf("%.1f", tv[1]),
    ATE_RR_CI = "",
    ATT_RR = sprintf("%.1f", tv[2]),
    ATT_RR_CI = "",
    ATE_RD = sprintf("%.3f", tv[3]),
    ATE_RD_CI = "",
    ATT_RD = sprintf("%.3f", tv[4]),
    ATT_RD_CI = ""
  )
  
  rows[[length(rows) + 1]] <- tibble(
    Effect = eff_name,
    Method = "g-computation",
    ATE_RR = sprintf("%.1f", gc_pt[1]),
    ATE_RR_CI = sprintf("%.1f, %.1f", gc_ci[1,1], gc_ci[1,2]),
    ATT_RR = sprintf("%.1f", gc_pt[2]),
    ATT_RR_CI = sprintf("%.1f, %.1f", gc_ci[2,1], gc_ci[2,2]),
    ATE_RD = sprintf("%.3f", gc_pt[3]),
    ATE_RD_CI = sprintf("%.3f, %.3f", gc_ci[3,1], gc_ci[3,2]),
    ATT_RD = sprintf("%.3f", gc_pt[4]),
    ATT_RD_CI = sprintf("%.3f, %.3f", gc_ci[4,1], gc_ci[4,2])
  )
}

tab3_steps <- bind_rows(rows)

# ----- Average short term effect assuming constant effect -----
d_avg <- dat_long %>%
  arrange(id, t) %>%
  group_by(id) %>%
  mutate(
    A_lag = lag(SEN),
    H_lag = lag(Hosp)
  ) %>%
  ungroup() %>%
  filter(t > 1) %>%
  mutate(
    inter = as.integer(idacibin == 1 & eyfsp_bin == 0),
    D = den
  ) %>%
  mutate(
    across(c(male, white, idacibin, eyfsp_bin, A_lag, H_lag, inter), ~ factor(.x, levels = c(0, 1))),
    t = factor(t)
  )

gcomp_avg_fun <- function(d) {
  meanD     <- mean(d$D)
  meanD_att <- mean(d$D[d$A_lag == 1])
  
  fit <- glm(
    Y ~ A_lag * (male + white + idacibin + eyfsp_bin + H_lag + t),
    family = poisson(),
    data = d,
    offset = log(D)
  )
  
  nd1 <- d; nd1$A_lag <- factor(1, levels = c(0, 1))
  nd0 <- d; nd0$A_lag <- factor(0, levels = c(0, 1))
  
  mu1 <- predict(fit, newdata = nd1, type = "response")
  mu0 <- predict(fit, newdata = nd0, type = "response")
  
  mu1_bar <- mean(mu1)
  mu0_bar <- mean(mu0)
  
  rr_ate <- mu1_bar / mu0_bar
  rd_ate <- (mu1_bar / meanD) - (mu0_bar / meanD)
  
  dt <- d %>% filter(A_lag == 1)
  nd1t <- dt; nd1t$A_lag <- factor(1, levels = c(0, 1))
  nd0t <- dt; nd0t$A_lag <- factor(0, levels = c(0, 1))
  
  mu1t <- predict(fit, newdata = nd1t, type = "response")
  mu0t <- predict(fit, newdata = nd0t, type = "response")
  
  mu1t_bar <- mean(mu1t)
  mu0t_bar <- mean(mu0t)
  
  rr_att <- mu1t_bar / mu0t_bar
  rd_att <- (mu1t_bar / meanD_att) - (mu0t_bar / meanD_att)
  
  c(ATE_RR = rr_pct(rr_ate), ATT_RR = rr_pct(rr_att),
    ATE_RD = rd_ate,        ATT_RD = rd_att)
}

# cluster bootstrap by id (keeps within-id correlation across time)
boot_est_cluster <- function(data, id_var = "id", B = 100, seed = 1212, est_fun) {
  set.seed(seed)
  split_list <- split(data, data[[id_var]])
  n_ids <- length(split_list)
  
  point <- tryCatch(est_fun(data), error = function(e) rep(NA_real_, 4))
  
  boot_mat <- replicate(B, {
    samp <- sample.int(n_ids, n_ids, replace = TRUE)
    d_b <- bind_rows(split_list[samp])
    tryCatch(est_fun(d_b), error = function(e) rep(NA_real_, 4))
  }) %>% t()
  
  ci <- apply(boot_mat, 2, function(z) {
    if (all(is.na(z))) c(NA_real_, NA_real_)
    else as.numeric(quantile(z, c(0.025, 0.975), na.rm = TRUE))
  }) %>% t()
  
  list(point = point, ci = ci)
}

avg_res <- boot_est_cluster(d_avg, B = 100, seed = 1212, est_fun = gcomp_avg_fun)

avg_row <- tibble(
  Effect = "Average of S_t on Y_(t+1)",
  Method = "g-computation",
  ATE_RR = sprintf("%.1f", avg_res$point[1]),
  ATE_RR_CI = sprintf("%.1f, %.1f", avg_res$ci[1,1], avg_res$ci[1,2]),
  ATT_RR = sprintf("%.1f", avg_res$point[2]),
  ATT_RR_CI = sprintf("%.1f, %.1f", avg_res$ci[2,1], avg_res$ci[2,2]),
  ATE_RD = sprintf("%.3f", avg_res$point[3]),
  ATE_RD_CI = sprintf("%.3f, %.3f", avg_res$ci[3,1], avg_res$ci[3,2]),
  ATT_RD = sprintf("%.3f", avg_res$point[4]),
  ATT_RD_CI = sprintf("%.3f, %.3f", avg_res$ci[4,1], avg_res$ci[4,2])
)

tab3 <- bind_rows(tab3_steps, avg_row) %>%
  select(
    Effect, Method,
    ATE_RR, ATE_RR_CI, ATT_RR, ATT_RR_CI,
    ATE_RD, ATE_RD_CI, ATT_RD, ATT_RD_CI
  )

print(tab3, n = Inf)

# Optional: save for GitHub
# write_csv(tab3, "results/Table3_R.csv")
