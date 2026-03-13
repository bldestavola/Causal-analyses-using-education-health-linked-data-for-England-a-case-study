#############################################
# Challenges in addressing causal questions #
# using administrative data                 #
#                                           #
# REPRODUCING TABLE 2)                      #
# Long-term effect of SEND (Year 1) on      #
# cumulative unauthorised absences Y2–Y4    #
#                                           #
# Outputs: ATE and ATT on                   #
#  - Rate Ratio (%) = 100 * RR              #
#  - Rate Difference (RD) on rate scale     #
# CIs: nonparametric bootstrap (B = 1000)   #
#############################################
install.packages("tidyverse")
install.packages("haven")
install.packages("AER")   

library(tidyverse)
library(haven)
library(AER)

setwd("")   # set up working directory here

dat_long <- haven::read_dta("data/Simulation_timevar_SEND_Study.dta")

dat <- dat_long %>%
  pivot_wider(id_cols = c(id, male, white, idacibin, eyfsp_bin, region),
              names_from = t,
              values_from = c(Hosp, SEN, Y, den, Y0_, Y1_),
              names_glue = "{.value}{t}") %>%
  rename(
    A = SEN1,     # exposure at baseline
    H = Hosp1     # time-varying confounder at baseline
  ) %>%
  mutate(
    Y    = Y2 + Y3 + Y4,  #cumulative Y2+Y3+Y4
    D    = den2 + den3 + den4,  #cumulative denominator
    rate = Y / D,
    inter = as.integer(idacibin == 1 & eyfsp_bin == 0)) %>%
  mutate(across(c(male, white, idacibin, eyfsp_bin, H, A, inter), ~ factor(.x, levels = c(0, 1))),
         region = factor(region))

mean_D     <- mean(dat$D)
mean_D_att <- mean(dat$D[dat$A == 1])

dat <- dat %>%
  mutate(
    Y0 = Y0_2 + Y0_3 + Y0_4,
    Y1 = Y1_2 + Y1_3 + Y1_4
  )

# Computing True values
true_ate_rrp <- 100 * (mean(dat$Y1) / mean(dat$Y0))
true_ate_rd  <- (mean(dat$Y1) / mean_D) - (mean(dat$Y0) / mean_D)

true_att_rrp <- 100 * (mean(dat$Y1[dat$A == 1]) / mean(dat$Y0[dat$A == 1]))
true_att_rd  <- (mean(dat$Y1[dat$A == 1]) / mean_D_att) - (mean(dat$Y0[dat$A == 1]) / mean_D_att)

true_row <- tibble(
  Method = "True value",
  Model  = "",
  ATE_RR_95CI = sprintf("%.1f", true_ate_rrp),
  ATT_RR_95CI = sprintf("%.1f", true_att_rrp),
  ATE_RD_95CI = sprintf("%.3f", true_ate_rd),
  ATT_RD_95CI = sprintf("%.3f", true_att_rd)
)

## ---- bootstrap helper ----
boot_est <- function(data, B = 1000, seed = 1212, est_fun) {
  set.seed(seed)
  n <- nrow(data)
  
  # point estimate
  point <- tryCatch(est_fun(data), error = function(e) rep(NA_real_, 4))
  
  # bootstrap replicates
  boot_mat <- replicate(B, {
    idx <- sample.int(n, n, replace = TRUE)
    tryCatch(est_fun(data[idx, , drop = FALSE]), error = function(e) rep(NA_real_, 4))
  })
  boot_mat <- t(boot_mat) 
  
  # percentile CI
  ci <- apply(boot_mat, 2, function(z) {
    if (all(is.na(z))) c(NA_real_, NA_real_) else as.numeric(quantile(z, c(0.025, 0.975), na.rm = TRUE))
  })
  ci <- t(ci)
  
  list(point = point, ci = ci)
}

## ---- effect scales ----
# Rate Ratio (%) in the paper: positive numbers mean reduction
rr_pct <- function(rr) 100 * (rr)

## ---- 1) G-computation (Poisson on counts with offset log(D)) ----
f_g1 <- Y ~ A + male + white + idacibin + eyfsp_bin + H
f_g2 <- Y ~ male + white + idacibin + H + A * eyfsp_bin
f_g3 <- Y ~ A * male * white * idacibin * eyfsp_bin * H

gcomp_estimator <- function(formula) {
  function(d) {
    fit <- glm(formula, family = poisson(), data = d, offset = log(D))
    
    # ATE: average predicted rates under A=1 and A=0
    nd1 <- d; nd1$A <- factor(1, levels = c(0, 1))
    nd0 <- d; nd0$A <- factor(0, levels = c(0, 1))
    
    mu1 <- predict(fit, newdata = nd1, type = "response")
    mu0 <- predict(fit, newdata = nd0, type = "response")
    
    r1 <- mean(mu1 / d$D)
    r0 <- mean(mu0 / d$D)
    
    ate_rrp <- rr_pct(r1 / r0)
    ate_rd  <- r1 - r0
    
    # ATT: average predicted rates among treated under A=1 and A=0
    dt <- d %>% filter(A == 1)
    nd1t <- dt; nd1t$A <- factor(1, levels = c(0, 1))
    nd0t <- dt; nd0t$A <- factor(0, levels = c(0, 1))
    
    mu1t <- predict(fit, newdata = nd1t, type = "response")
    mu0t <- predict(fit, newdata = nd0t, type = "response")
    
    rt1 <- mean(mu1t / dt$D)
    rt0 <- mean(mu0t / dt$D)
    
    att_rrp <- rr_pct(rt1 / rt0)
    att_rd  <- rt1 - rt0
    
    c(ate_rrp, att_rrp, ate_rd, att_rd)
  }
}

g1 <- boot_est(dat, est_fun = gcomp_estimator(f_g1))
g2 <- boot_est(dat, est_fun = gcomp_estimator(f_g2))
g3 <- boot_est(dat, est_fun = gcomp_estimator(f_g3))


## ---- 2) IPW ----
ps_incorrect <- A ~ male + white + idacibin + eyfsp_bin + H
ps_correct   <- A ~ male + white + idacibin + eyfsp_bin + inter + H

ipw_estimator <- function(ps_formula) {
  function(d) {
    ps_fit <- glm(ps_formula, family = binomial(), data = d)
    ps <- as.numeric(predict(ps_fit, type = "response"))
    ps <- pmin(pmax(ps, 1e-6), 1 - 1e-6)
    
    A_num <- as.integer(as.character(d$A))
    y <- d$rate
    
    # ATE: E[Y1] and E[Y0]
    EY1 <- mean(A_num * y / ps)
    EY0 <- mean((1 - A_num) * y / (1 - ps))
    
    ate_rrp <- rr_pct(EY1 / EY0)
    ate_rd  <- EY1 - EY0
    
    # ATT: E[Y1|A=1] and E[Y0|A=1]
    n1 <- sum(A_num == 1)
    EY1_att <- mean(y[A_num == 1])
    EY0_att <- sum((A_num == 0) * y * ps / (1 - ps)) / n1
    
    att_rrp <- rr_pct(EY1_att / EY0_att)
    att_rd  <- EY1_att - EY0_att
    
    c(ate_rrp, att_rrp, ate_rd, att_rd)
  }
}

ipw1 <- boot_est(dat, est_fun = ipw_estimator(ps_incorrect))
ipw2 <- boot_est(dat, est_fun = ipw_estimator(ps_correct))


## ---- 3) AIPW (doubly robust; Poisson outcome model for rate) ----
y_incorrect <- rate ~ A + male + white + idacibin + eyfsp_bin + H
y_correct   <- rate ~ A * male * white * idacibin * eyfsp_bin * H

aipw_estimator <- function(y_formula, ps_formula) {
  function(d) {
    # outcome regression for rate
    y_fit <- glm(y_formula, family = poisson(), data = d)
    
    # propensity score model
    ps_fit <- glm(ps_formula, family = binomial(), data = d)
    ps <- as.numeric(predict(ps_fit, type = "response"))
    ps <- pmin(pmax(ps, 1e-6), 1 - 1e-6)
    
    A_num <- as.integer(as.character(d$A))
    y <- d$rate
    
    # predicted mu(a) for each unit
    nd1 <- d; nd1$A <- factor(1, levels = c(0, 1))
    nd0 <- d; nd0$A <- factor(0, levels = c(0, 1))
    mu1 <- as.numeric(predict(y_fit, newdata = nd1, type = "response"))
    mu0 <- as.numeric(predict(y_fit, newdata = nd0, type = "response"))
    
    # ATE AIPW
    EY1 <- mean(mu1 + (A_num / ps) * (y - mu1))
    EY0 <- mean(mu0 + ((1 - A_num) / (1 - ps)) * (y - mu0))
    
    ate_rrp <- rr_pct(EY1 / EY0)
    ate_rd  <- EY1 - EY0
    
    # ATT AIPW
    n1 <- sum(A_num == 1)
    EY1_att <- mean(y[A_num == 1])
    EY0_att <- (sum((A_num == 1) * mu0) + sum((A_num == 0) * ps / (1 - ps) * (y - mu0))) / n1
    
    att_rrp <- rr_pct(EY1_att / EY0_att)
    att_rd  <- EY1_att - EY0_att
    
    c(ate_rrp, att_rrp, ate_rd, att_rd)
  }
}

aipw_ii <- boot_est(dat, est_fun = aipw_estimator(y_incorrect, ps_incorrect)) # incorrect Y, incorrect PS
aipw_ic <- boot_est(dat, est_fun = aipw_estimator(y_incorrect, ps_correct))   # incorrect Y, correct PS
aipw_ci <- boot_est(dat, est_fun = aipw_estimator(y_correct,   ps_incorrect)) # correct Y, incorrect PS
aipw_cc <- boot_est(dat, est_fun = aipw_estimator(y_correct,   ps_correct))   # correct Y, correct PS

## ---- 4) IV (Two-stage least squares) on RD scale only ----
iv_res <- boot_est(dat, est_fun = function(d) {
  fit <- ivreg(rate ~ A + male + white + idacibin + eyfsp_bin + H |
                 region + male + white + idacibin + eyfsp_bin + H,
               data = d)
  bA <- coef(fit)[["A1"]]
  c(NA_real_, NA_real_, bA, bA)  
})

## ---- 5) Standard regression (conditional effect) ----
std_res <- boot_est(dat, est_fun = function(d) {
  fit <- glm(Y ~ A + male + white + idacibin + eyfsp_bin + H,
             family = poisson(), data = d, offset = log(D))
  b0 <- coef(fit)[["(Intercept)"]]
  bA <- coef(fit)[["A1"]]
  
  # conditional RR% and RD at reference covariate pattern 
  rr  <- exp(bA)
  rd  <- exp(b0 + bA) - exp(b0)
  
  c(rr_pct(rr), NA_real_, rd, NA_real_)
})


## ---- building Table 2-style output ----
make_row <- function(method, model, obj) {
  est <- obj$point
  ci  <- obj$ci
  
  tibble(
    Method = method,
    Model  = model,
    
    ATE_RR_pct = est[1],
    ATE_RR_95CI = ifelse(is.na(est[1]), NA_character_,
                         sprintf("%.1f (%.1f, %.1f)", est[1], ci[1,1], ci[1,2])),
    
    ATT_RR_pct = est[2],
    ATT_RR_95CI = ifelse(is.na(est[2]), NA_character_,
                         sprintf("%.1f (%.1f, %.1f)", est[2], ci[2,1], ci[2,2])),
    
    ATE_RD = est[3],
    ATE_RD_95CI = ifelse(is.na(est[3]), NA_character_,
                         sprintf("%.3f (%.3f, %.3f)", est[3], ci[3,1], ci[3,2])),
    
    ATT_RD = est[4],
    ATT_RD_95CI = ifelse(is.na(est[4]), NA_character_,
                         sprintf("%.3f (%.3f, %.3f)", est[4], ci[4,1], ci[4,2]))
  )
}

tab2 <- bind_rows(
  true_row,
  make_row("G-computation", "Incorrect Y-version 1", g1),
  make_row("G-computation", "Incorrect Y-version 2", g2),
  make_row("G-computation", "Correct Y-version 3",   g3),
  
  make_row("IPW", "Incorrect PS", ipw1),
  make_row("IPW", "Correct PS",   ipw2),
  
  make_row("AIPW", "Incorrect Y, incorrect PS", aipw_ii),
  make_row("AIPW", "Incorrect Y, correct PS",   aipw_ic),
  make_row("AIPW", "Correct Y, incorrect PS",   aipw_ci),
  make_row("AIPW", "Correct Y, correct PS",     aipw_cc),
  
  make_row("IV-based", "2SLS", iv_res),
  make_row("Standard regression", "Conditional effect", std_res)
) %>%
  select(Method, Model,
         ATE_RR_95CI, ATT_RR_95CI,
         ATE_RD_95CI, ATT_RD_95CI)

print(tab2, n = Inf)

# Optional: save tab2
# write_csv(tab2, "results/Table2_R.csv")
