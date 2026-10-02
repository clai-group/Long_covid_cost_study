# =============================================================================
# Two-part GEE for longitudinal cost data, with g-computation and sensitivity
# analyses.
#
# Input `dat`: one row per person-period with columns
#   id        person identifier
#   period    time since index (0, 1, 2, ...)
#   exposed   1 = exposed group, 0 = comparison group
#   cost      cost in the period (0 if no utilization)
#   visits    number of encounters in the period
#   followup  last observable period for the person
#   plus the covariates listed in `covars` (categorical covariates as factors)
# =============================================================================
library(dplyr)
library(tidyr)
library(geepack)
library(broom)
library(MASS)

set.seed(123)
B       <- 400                       # parametric bootstrap replications
periods <- 0:15                      # periods to estimate
covars  <- c("age", "I(age^2)", "sex", "comorbidity", "baseline_cost")  # edit as needed
cov_cols <- c("age", "sex", "comorbidity", "baseline_cost")             # raw columns used above

make_rhs <- function(time_term, cv = covars)
  reformulate(c(paste("exposed *", time_term), cv))

# geeglm requires each person's rows to be contiguous
dat <- dat %>% arrange(id, period) %>% mutate(any_cost = as.integer(cost > 0))

# ---- Two-part model: P(any cost) x E(cost | cost > 0) ------------------------
fit_two_part <- function(data, rhs, corstr = "exchangeable") {
  list(
    p1 = geeglm(update(rhs, any_cost ~ .), data = data, id = id,
                family = binomial("logit"), corstr = corstr),
    p2 = geeglm(update(rhs, cost ~ .), data = filter(data, cost > 0), id = id,
                family = Gamma("log"), corstr = corstr)
  )
}

# ---- g-computation: average treatment effect on the treated (exposed = 1 vs 0)
#      per period, averaged over the exposed group's covariates. CIs from a
#      parametric bootstrap drawing both model parts (robust covariance).
gcomp <- function(fit, rhs, data, time_factor = FALSE) {
  cov_exp <- data %>% filter(exposed == 1) %>% distinct(id, .keep_all = TRUE) %>%
    dplyr::select(any_of(cov_cols))
  frame <- tidyr::crossing(cov_exp, period = periods)
  if (time_factor) frame$period_f <- factor(frame$period, levels = periods)
  b1 <- coef(fit$p1); b2 <- coef(fit$p2)
  X <- function(e) { frame$exposed <- e; model.matrix(rhs, frame)[, names(b1)] }
  X1 <- X(1); X0 <- X(0)
  ecost <- function(X, d1, d2) tapply(as.vector(plogis(X %*% d1) * exp(X %*% d2)), frame$period, mean)
  effect <- function(d1, d2) ecost(X1, d1, d2) - ecost(X0, d1, d2)
  boot <- replicate(B, effect(mvrnorm(1, b1, vcov(fit$p1)), mvrnorm(1, b2, vcov(fit$p2))))
  list(per_period = tibble(period = periods, effect = effect(b1, b2),
                           lo = apply(boot, 1, quantile, .025),
                           hi = apply(boot, 1, quantile, .975)),
       boot = boot)
}

# Cumulative effect over selected periods (point estimate = bootstrap mean)
cumulative <- function(g, which_periods) {
  i <- match(which_periods, periods)
  tot <- colSums(g$boot[i, , drop = FALSE])
  c(estimate = sum(rowMeans(g$boot)[i]),
    lo = quantile(tot, .025, names = FALSE), hi = quantile(tot, .975, names = FALSE))
}

# ---- Primary model -----------------------------------------------------------
rhs     <- make_rhs("period")
primary <- fit_two_part(dat, rhs)
tidy(primary$p1, exponentiate = TRUE, conf.int = TRUE)   # odds ratios
tidy(primary$p2, exponentiate = TRUE, conf.int = TRUE)   # cost ratios
g_primary <- gcomp(primary, rhs, dat)
cumulative(g_primary, periods)

# ---- Sensitivity: categorical time --------------------------------------------
dat_cat <- dat %>% mutate(period_f = factor(period, levels = periods))
rhs_cat <- make_rhs("period_f")
g_cat   <- gcomp(fit_two_part(dat_cat, rhs_cat), rhs_cat, dat_cat, time_factor = TRUE)
cumulative(g_cat, periods)

# ---- Sensitivity: drop a covariate (e.g., a possible mediator) ---------------
rhs_drop <- make_rhs("period", setdiff(covars, "comorbidity"))
g_drop   <- gcomp(fit_two_part(dat, rhs_drop), rhs_drop, dat)
cumulative(g_drop, periods)

# ---- Sensitivity: persons observable through the last period -----------------
dat_full <- dat %>% filter(followup >= max(periods)) %>% droplevels()
g_full   <- gcomp(fit_two_part(dat_full, rhs), rhs, dat_full)
cumulative(g_full, periods)

# ---- Sensitivity: working correlation structure ------------------------------
for (cs in c("ar1", "independence")) {
  fit <- fit_two_part(dat, rhs, corstr = cs)
  print(bind_rows(tidy(fit$p1, exponentiate = TRUE, conf.int = TRUE),
                  tidy(fit$p2, exponentiate = TRUE, conf.int = TRUE)) %>%
          filter(term == "exposed:period") %>% mutate(corstr = cs))
}

# ---- Decomposition: visit frequency and cost per visit -----------------------
visit_model <- geeglm(update(rhs, visits ~ .), data = dat, id = id,
                      family = poisson("log"), corstr = "exchangeable")
cpv_model   <- geeglm(update(rhs, cost_per_visit ~ .),
                      data = dat %>% filter(visits > 0, cost > 0) %>%
                        mutate(cost_per_visit = cost / visits),
                      id = id, family = Gamma("log"), corstr = "exchangeable")
tidy(visit_model, exponentiate = TRUE, conf.int = TRUE)   # rate ratios
tidy(cpv_model,   exponentiate = TRUE, conf.int = TRUE)   # cost-per-visit ratios



