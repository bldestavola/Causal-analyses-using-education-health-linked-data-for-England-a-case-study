#############################################
# Challenges in addressing causal questions #
# using administrative data                 #
#                                           #
# REPRODUCING TABLE 4                       #
# Joint effect of (S1, S2, S3) on (Y2+Y3+Y4)#
#                                           #
# Outputs (as in the paper):                #
#  - Risk Ratio (%) = 100 * RR              #
#  - Risk Difference (RD) on rate scale     #
# CIs: bootstrap (B = 100)                 #
#############################################

install.packages("tidyverse")
install.packages("haven")

library(tidyverse)
library(haven)

setwd("")  #setting working directory

dat_long0 <- read_dta("data/Simulation_timevar_SEND_Study.dta")

TRUE_RR_PCT <- 12.8
TRUE_RD <- -0.089

# -----------------------------
# bootstrap function by id 
# -----------------------------
boot_method <- function(split_list, B = 100, seed = 1012, est_fun) {
  set.seed(seed)
  n_ids <- length(split_list)
  
  d_point <- bind_rows(Map(function(df, k) mutate(df, idb = k), split_list, seq_along(split_list)))
  point <- est_fun(d_point)
  
  boot_mat <- replicate(B, {
    samp <- sample.int(n_ids, n_ids, replace = TRUE)
    d_b  <- bind_rows(Map(function(df, k) mutate(df, idb = k), split_list[samp], seq_along(samp)))
    est_fun(d_b)
  }) %>% t()
  
  ci <- apply(boot_mat, 2, function(z) quantile(z, c(0.025, 0.975), na.rm = TRUE)) %>% t()
  list(point = point, ci = ci)
}

# ----------------
# Preparing data
# ----------------
prep_data <- function(dlong) {
  d <- dlong %>%
    arrange(idb, t) %>%
    group_by(idb) %>%
    mutate(
      SEN_lag1  = lag(SEN, 1),
      SEN_lag2  = lag(SEN, 2),
      SEN_lag3  = lag(SEN, 3),
      Hosp_lag1 = lag(Hosp, 1),
      Hosp_lag2 = lag(Hosp, 2),
      Hosp_lag3 = lag(Hosp, 3),
      Hosp_base = first(Hosp),
      # cumulative outcome and denom to end of year 4 (t==4)
      cumY = cumsum(if_else(t > 1, as.numeric(Y), 0)),
      cumD = cumsum(if_else(t > 1, as.numeric(den), 0)),
      inter = as.integer(eyfsp_bin == 0 & SEN_lag1 == 1),
      inter_ps = as.integer(idacibin == 1 & eyfsp_bin == 0)
    ) %>%
    ungroup()
  
  d4 <- d %>%
    filter(t == 4) %>%
    transmute(
      idb,
      male, white, idacibin, eyfsp_bin,
      Hosp_base,
      SEN_lag1, SEN_lag2, SEN_lag3,
      Hosp_lag1, Hosp_lag2, Hosp_lag3,
      inter, inter_ps,
      cumY, cumD
    ) %>%
    mutate(
      across(c(male, white, idacibin, eyfsp_bin,
               Hosp_base, SEN_lag1, SEN_lag2, SEN_lag3,
               Hosp_lag1, Hosp_lag2, Hosp_lag3,
               inter, inter_ps), ~ factor(.x, levels = c(0, 1)))
    )
  
# for sustained g-comp 
  base <- d %>%
    group_by(idb) %>%
    summarise(
      male      = first(male),
      white     = first(white),
      idacibin  = first(idacibin),
      eyfsp_bin = first(eyfsp_bin),
      Hosp1 = Hosp[t == 1][1],
      den2  = den[t == 2][1],
      den3  = den[t == 3][1],
      den4  = den[t == 4][1],
      .groups = "drop"
    ) %>%
    mutate(
      cumD = as.numeric(den2) + as.numeric(den3) + as.numeric(den4),
      across(c(male, white, idacibin, eyfsp_bin, Hosp1), ~ factor(.x, levels = c(0, 1)))
    )
  
  list(long = d, d4 = d4, base = base)
}

# -----------------------------
# Method 1: Regression adjustment (Poisson on cumulative outcome at t==4)
# -----------------------------
est_reg_adj <- function(dlong) {
  p <- prep_data(dlong)
  d4 <- p$d4
  
  fit <- glm(
    cumY ~ male + white + idacibin + eyfsp_bin +
      Hosp_lag1 + Hosp_lag2 + Hosp_lag3 +
      SEN_lag1 + SEN_lag2 + SEN_lag3 + inter,
    family = poisson(),
    data = d4,
    offset = log(as.numeric(cumD))
  )
  
  nd0 <- d4 %>%
    mutate(
      SEN_lag1 = factor(0, levels = c(0, 1)),
      SEN_lag2 = factor(0, levels = c(0, 1)),
      SEN_lag3 = factor(0, levels = c(0, 1)),
      inter    = factor(0, levels = c(0, 1))
    )
  
  nd1 <- d4 %>%
    mutate(
      SEN_lag1 = factor(1, levels = c(0, 1)),
      SEN_lag2 = factor(1, levels = c(0, 1)),
      SEN_lag3 = factor(1, levels = c(0, 1)),
      inter    = factor(as.integer(as.numeric(as.character(eyfsp_bin)) == 0), levels = c(0, 1))
    )
  
  mu0 <- predict(fit, newdata = nd0, type = "response")
  mu1 <- predict(fit, newdata = nd1, type = "response")
  
  rate0 <- mean(mu0) / mean(as.numeric(d4$cumD))
  rate1 <- mean(mu1) / mean(as.numeric(d4$cumD))
  
  rr <- rate1 / rate0
  rd <- rate1 - rate0
  c(RR_pct = 100 * rr, RD = rd)
}

# -----------------------------
# Method 2 and 3: Sustained g-computation (version 1 simple, version 2 full)
# Parametric g-formula with iterated expectations over binary Hosp
# -----------------------------
est_gcomp_sustained <- function(dlong, version = 1) {
  p <- prep_data(dlong)
  d  <- p$long
  b  <- p$base
  
  dH <- d %>%
    filter(t %in% c(2, 3)) %>%
    mutate(
      # predictors used in models
      SEN_lag1  = factor(SEN_lag1, levels = c(0, 1)),
      Hosp_lag1 = factor(Hosp_lag1, levels = c(0, 1)),
      inter     = factor(as.integer(eyfsp_bin == 0 & SEN_lag1 == 1), levels = c(0, 1)),
      t         = factor(t),
      Hosp      = factor(Hosp, levels = c(0, 1)),
      male      = factor(male, levels = c(0, 1)),
      white     = factor(white, levels = c(0, 1)),
      idacibin  = factor(idacibin, levels = c(0, 1)),
      eyfsp_bin = factor(eyfsp_bin, levels = c(0, 1))
    )
  
  dY <- d %>%
    filter(t %in% c(2, 3, 4)) %>%
    mutate(
      SEN_lag1  = factor(SEN_lag1, levels = c(0, 1)),
      Hosp_lag1 = factor(Hosp_lag1, levels = c(0, 1)),
      inter     = factor(as.integer(eyfsp_bin == 0 & SEN_lag1 == 1), levels = c(0, 1)),
      t         = factor(t),
      male      = factor(male, levels = c(0, 1)),
      white     = factor(white, levels = c(0, 1)),
      idacibin  = factor(idacibin, levels = c(0, 1)),
      eyfsp_bin = factor(eyfsp_bin, levels = c(0, 1))
    )
  
  # Version 1: simple (incorrect)
  if (version == 1) {
    fitH <- glm(
      Hosp ~ SEN_lag1 + male + white + idacibin + eyfsp_bin + t,
      family = binomial(),
      data = dH
    )
    
    fitY <- glm(
      Y ~ SEN_lag1 + male + white + idacibin + eyfsp_bin + t,
      family = poisson(),
      data = dY,
      offset = log(as.numeric(den))
    )
  }
  
  # Version 2: fuller (more general)
  if (version == 2) {
    fitH <- glm(
      Hosp ~ SEN_lag1 + Hosp_lag1 + male + white + idacibin + eyfsp_bin + inter + t,
      family = binomial(),
      data = dH
    )
    
    fitY <- glm(
      Y ~ SEN_lag1 + Hosp_lag1 + male + white + idacibin + eyfsp_bin + inter + t,
      family = poisson(),
      data = dY,
      offset = log(as.numeric(den))
    )
  }
  
  regime_rate <- function(a1, a2, a3) {
    # t=2 pieces
    ndY2 <- b %>%
      transmute(
        male, white, idacibin, eyfsp_bin,
        SEN_lag1  = factor(a1, levels = c(0, 1)),
        Hosp_lag1 = factor(Hosp1, levels = c(0, 1)),
        inter     = factor(as.integer(as.numeric(as.character(eyfsp_bin)) == 0 & a1 == 1), levels = c(0, 1)),
        t         = factor(2),
        den       = as.numeric(den2),
        Y         = 0
      )
    muY2 <- predict(fitY, newdata = ndY2, type = "response")
    
    ndH2 <- b %>%
      transmute(
        male, white, idacibin, eyfsp_bin,
        SEN_lag1  = factor(a1, levels = c(0, 1)),
        Hosp_lag1 = factor(Hosp1, levels = c(0, 1)),
        inter     = factor(as.integer(as.numeric(as.character(eyfsp_bin)) == 0 & a1 == 1), levels = c(0, 1)),
        t         = factor(2),
        Hosp      = factor(0, levels = c(0, 1))
      )
    pH2 <- predict(fitH, newdata = ndH2, type = "response")
    
    # t=3 pieces (mixture over H2)
    ndY3_0 <- b %>%
      transmute(
        male, white, idacibin, eyfsp_bin,
        SEN_lag1  = factor(a2, levels = c(0, 1)),
        Hosp_lag1 = factor(0, levels = c(0, 1)),
        inter     = factor(as.integer(as.numeric(as.character(eyfsp_bin)) == 0 & a2 == 1), levels = c(0, 1)),
        t         = factor(3),
        den       = as.numeric(den3),
        Y         = 0
      )
    ndY3_1 <- ndY3_0 %>% mutate(Hosp_lag1 = factor(1, levels = c(0, 1)))
    
    muY3_0 <- predict(fitY, newdata = ndY3_0, type = "response")
    muY3_1 <- predict(fitY, newdata = ndY3_1, type = "response")
    EY3 <- (1 - pH2) * muY3_0 + pH2 * muY3_1
    
    ndH3_0 <- b %>%
      transmute(
        male, white, idacibin, eyfsp_bin,
        SEN_lag1  = factor(a2, levels = c(0, 1)),
        Hosp_lag1 = factor(0, levels = c(0, 1)),
        inter     = factor(as.integer(as.numeric(as.character(eyfsp_bin)) == 0 & a2 == 1), levels = c(0, 1)),
        t         = factor(3),
        Hosp      = factor(0, levels = c(0, 1))
      )
    ndH3_1 <- ndH3_0 %>% mutate(Hosp_lag1 = factor(1, levels = c(0, 1)))
    
    pH3_0 <- predict(fitH, newdata = ndH3_0, type = "response")
    pH3_1 <- predict(fitH, newdata = ndH3_1, type = "response")
    
    # t=4 mixture over H2 then H3|H2
    ndY4_0 <- b %>%
      transmute(
        male, white, idacibin, eyfsp_bin,
        SEN_lag1  = factor(a3, levels = c(0, 1)),
        Hosp_lag1 = factor(0, levels = c(0, 1)),
        inter     = factor(as.integer(as.numeric(as.character(eyfsp_bin)) == 0 & a3 == 1), levels = c(0, 1)),
        t         = factor(4),
        den       = as.numeric(den4),
        Y         = 0
      )
    ndY4_1 <- ndY4_0 %>% mutate(Hosp_lag1 = factor(1, levels = c(0, 1)))
    
    muY4_0 <- predict(fitY, newdata = ndY4_0, type = "response")
    muY4_1 <- predict(fitY, newdata = ndY4_1, type = "response")
    
    EY4 <- (1 - pH2) * ((1 - pH3_0) * muY4_0 + pH3_0 * muY4_1) +
      pH2 * ((1 - pH3_1) * muY4_0 + pH3_1 * muY4_1)
    
    cum_count <- muY2 + EY3 + EY4
    rate <- mean(cum_count) / mean(b$cumD)
    rate
  }
  
  r0 <- regime_rate(0, 0, 0)
  r1 <- regime_rate(1, 1, 1)
  
  rr <- r1 / r0
  rd <- r1 - r0
  c(RR_pct = 100 * rr, RD = rd)
}

est_gcomp_v1 <- function(dlong) est_gcomp_sustained(dlong, version = 1)
est_gcomp_v2 <- function(dlong) est_gcomp_sustained(dlong, version = 2)

# -----------------------------
# Method 4: IPW 
# -----------------------------
est_ipw <- function(dlong) {
  p <- prep_data(dlong)
  d  <- p$long
  d4 <- p$d4 %>%
    mutate(
      rate = as.numeric(cumY) / as.numeric(cumD),
      S1 = as.numeric(as.character(SEN_lag3)),  # SEN at t=1
      S2 = as.numeric(as.character(SEN_lag2)),  # SEN at t=2
      S3 = as.numeric(as.character(SEN_lag1))   # SEN at t=3
    )
  
  # exposure model data at t=1,2,3
  treat <- d %>%
    filter(t %in% 1:3) %>%
    arrange(idb, t) %>%
    group_by(idb) %>%
    mutate(
      Hosp_base = first(Hosp),
      Hosp_lag1 = lag(Hosp, 1),
      Hosp_lag1 = replace_na(Hosp_lag1, 0),
      inter_ps  = as.integer(idacibin == 1 & eyfsp_bin == 0),
      A         = SEN
    ) %>%
    ungroup() %>%
    mutate(
      t        = factor(t),
      A        = factor(A, levels = c(0, 1)),
      male     = factor(male, levels = c(0, 1)),
      white    = factor(white, levels = c(0, 1)),
      idacibin = factor(idacibin, levels = c(0, 1)),
      eyfsp_bin= factor(eyfsp_bin, levels = c(0, 1)),
      Hosp_base= factor(Hosp_base, levels = c(0, 1)),
      Hosp     = factor(Hosp, levels = c(0, 1)),
      Hosp_lag1= factor(Hosp_lag1, levels = c(0, 1)),
      inter_ps = factor(inter_ps, levels = c(0, 1))
    )
  
  # numerator: baseline only (+ time)
  fit_num <- glm(
    A ~ male + white + Hosp_base + idacibin + eyfsp_bin + inter_ps + t,
    family = binomial(),
    data = treat
  )
  
  # denominator: baseline + time-varying confounding 
  fit_den <- glm(
    A ~ male + white + Hosp_base + idacibin + eyfsp_bin + inter_ps + t + Hosp + Hosp_lag1,
    family = binomial(),
    data = treat
  )
  
  pnum <- as.numeric(predict(fit_num, type = "response"))
  pden <- as.numeric(predict(fit_den, type = "response"))
  pnum <- pmin(pmax(pnum, 1e-6), 1 - 1e-6)
  pden <- pmin(pmax(pden, 1e-6), 1 - 1e-6)
  
  A_num <- as.numeric(as.character(treat$A))
  
  cnum <- if_else(A_num == 1, pnum, 1 - pnum)
  cden <- if_else(A_num == 1, pden, 1 - pden)
  
  sw <- treat %>%
    mutate(cnum = cnum, cden = cden) %>%
    group_by(idb) %>%
    summarise(sw = prod(cnum) / prod(cden), .groups = "drop")
  
  d4w <- d4 %>%
    left_join(sw, by = "idb") %>%
    mutate(
      sw = as.numeric(sw),
      s000 = (S1 == 0 & S2 == 0 & S3 == 0),
      s111 = (S1 == 1 & S2 == 1 & S3 == 1)
    )
  
  # normalised IPW regime means
  r0 <- sum(d4w$sw[d4w$s000] * d4w$rate[d4w$s000]) / sum(d4w$sw[d4w$s000])
  r1 <- sum(d4w$sw[d4w$s111] * d4w$rate[d4w$s111]) / sum(d4w$sw[d4w$s111])
  
  rr <- r1 / r0
  rd <- r1 - r0
  
  c(RR_pct = 100 * rr, RD = rd)
}


split_list <- split(dat_long0, dat_long0$id)

reg_res  <- boot_method(split_list, B = 100, seed = 1012, est_fun = est_reg_adj)
g1_res   <- boot_method(split_list, B = 100, seed = 1012, est_fun = est_gcomp_v1)
g2_res   <- boot_method(split_list, B = 100, seed = 1012, est_fun = est_gcomp_v2)
ipw_res  <- boot_method(split_list, B = 100, seed = 3009, est_fun = est_ipw)

# -----------------------------
# Building Table 4 output
# -----------------------------
fmt_rr <- function(x) sprintf("%.1f", x)
fmt_rr_ci <- function(lo, hi) sprintf("%.1f, %.1f", lo, hi)
fmt_rd <- function(x) sprintf("%.3f", x)
fmt_rd_ci <- function(lo, hi) sprintf("%.3f, %.3f", lo, hi)

tab4 <- tibble(
  Method = c("True value",
             "Regression adjustment",
             "G-computation- version 1",
             "G-computation- version 2",
             "IPW"),
  `Risk Ratio (%) Estimate` = c(
    fmt_rr(TRUE_RR_PCT),
    fmt_rr(reg_res$point[1]),
    fmt_rr(g1_res$point[1]),
    fmt_rr(g2_res$point[1]),
    fmt_rr(ipw_res$point[1])
  ),
  `Risk Ratio (%) 95% CI` = c(
    "",
    fmt_rr_ci(reg_res$ci[1,1], reg_res$ci[1,2]),
    fmt_rr_ci(g1_res$ci[1,1],  g1_res$ci[1,2]),
    fmt_rr_ci(g2_res$ci[1,1],  g2_res$ci[1,2]),
    fmt_rr_ci(ipw_res$ci[1,1], ipw_res$ci[1,2])
  ),
  `Risk Difference Estimate` = c(
    fmt_rd(TRUE_RD),
    fmt_rd(reg_res$point[2]),
    fmt_rd(g1_res$point[2]),
    fmt_rd(g2_res$point[2]),
    fmt_rd(ipw_res$point[2])
  ),
  `Risk Difference 95% CI` = c(
    "",
    fmt_rd_ci(reg_res$ci[2,1], reg_res$ci[2,2]),
    fmt_rd_ci(g1_res$ci[2,1],  g1_res$ci[2,2]),
    fmt_rd_ci(g2_res$ci[2,1],  g2_res$ci[2,2]),
    fmt_rd_ci(ipw_res$ci[2,1], ipw_res$ci[2,2])
  )
)

print(tab4, n = Inf)

# Optional: save for GitHub
# write_csv(tab4, "results/Table4_R.csv")
