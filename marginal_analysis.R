# =============================================================================
# Pre/post marginal analysis: piecewise two-part GEE with group-specific pre-
# and post-index slopes and a group-specific level change at index;
# g-computation over the exposed group; normalization to the mean pre-index
# difference.
#
# Input `dat`: one row per person-period with columns
#   id, period (negative = before index, 0 = index period), exposed (1/0),
#   cost, and the covariates listed in `covars` (categorical as factors)
# =============================================================================
library(dplyr)
library(tidyr)
library(geepack)
library(broom)
library(MASS)

B       <- 400
periods <- -8:15                    # pre- and post-index periods
covars  <- c("age", "I(age^2)", "sex", "comorbidity", "baseline_cost")  # edit as needed
cov_cols <- c("age", "sex", "comorbidity", "baseline_cost")

rhs <- reformulate(c("exposed + post + exposed:post + t_pre + t_post",
                     "exposed:t_pre + exposed:t_post", covars))

add_time <- function(d) d %>% mutate(post   = as.integer(period >= 0),
                                     t_pre  = pmin(period, -1) + 1,
                                     t_post = pmax(period, 0))

run_marginal <- function(dat) {
  set.seed(123)
  d <- dat %>% add_time() %>% mutate(any_cost = as.integer(cost > 0)) %>%
    arrange(id, period)                                   # contiguous clusters
  stopifnot(!anyDuplicated(rle(as.character(d$id))$values))

  p1 <- geeglm(update(rhs, any_cost ~ .), data = d, id = id,
               family = binomial("logit"), corstr = "exchangeable")
  p2 <- geeglm(update(rhs, cost ~ .), data = filter(d, cost > 0), id = id,
               family = Gamma("log"), corstr = "exchangeable")

  # Difference in pre-index slopes between groups
  pretrend <- bind_rows(tidy(p1, exponentiate = TRUE, conf.int = TRUE),
                        tidy(p2, exponentiate = TRUE, conf.int = TRUE)) %>%
    filter(term == "exposed:t_pre")

  # g-computation over the exposed group's covariates
  cov_exp <- d %>% filter(exposed == 1) %>% distinct(id, .keep_all = TRUE) %>%
    dplyr::select(any_of(cov_cols))
  frame <- tidyr::crossing(cov_exp, period = periods) %>% add_time()
  b1 <- coef(p1); b2 <- coef(p2)
  X <- function(e) { frame$exposed <- e; model.matrix(rhs, frame)[, names(b1)] }
  X1 <- X(1); X0 <- X(0)
  ecost <- function(X, e1, e2) tapply(as.vector(plogis(X %*% e1) * exp(X %*% e2)), frame$period, mean)
  gap <- function(e1, e2) ecost(X1, e1, e2) - ecost(X0, e1, e2)

  diff <- gap(b1, b2)
  boot <- replicate(B, gap(mvrnorm(1, b1, vcov(p1)), mvrnorm(1, b2, vcov(p2))))
  pre  <- periods < 0
  base <- mean(diff[pre]); base_boot <- colMeans(boot[pre, ])

  cumulative <- function(which_periods, normalize = TRUE) {
    i <- periods %in% which_periods
    tot <- colSums(boot[i, , drop = FALSE]) - normalize * sum(i) * base_boot
    c(estimate = sum(diff[i] - normalize * base), quantile(tot, c(.025, .975)))
  }

  list(
    per_period = tibble(period = periods, diff = diff,
                        lo = apply(boot, 1, quantile, .025),
                        hi = apply(boot, 1, quantile, .975)),
    normalized = tibble(period = periods[!pre], diff = diff[!pre] - base,
                        lo = apply(sweep(boot[!pre, ], 2, base_boot), 1, quantile, .025),
                        hi = apply(sweep(boot[!pre, ], 2, base_boot), 1, quantile, .975)),
    cumulative = rbind(normalized_all      = cumulative(0:max(periods)),
                       normalized_post_idx = cumulative(1:max(periods)),
                       unnormalized_all    = cumulative(0:max(periods), FALSE)),
    pretrend = pretrend
  )
}

# result <- run_marginal(dat)
