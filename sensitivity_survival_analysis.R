# =============================================================================
# Mortality-adjusted expected cost (unadjusted group means)
# Expected cost in period t = mean survivor cost x S(t)
#                           + mean cost in the period of death x P(death in t)
# Cumulative difference between groups with a Monte Carlo CI.
#
# Inputs:
#   persons      one row per person: id, group, time (years to death or
#                censoring), event (1 = died)
#   alive_panel  person-periods before death: id, group, period, cost
#   death_costs  one row per decedent: id, group, death_cost (cost in the
#                period of death)
# Period length is `len` years (0.25 = quarterly).
# =============================================================================
library(dplyr)
library(tidyr)
library(survival)

set.seed(123)
periods <- 0:15
len     <- 0.25
R       <- 2000

# Kaplan-Meier survival at period boundaries
km  <- summary(survfit(Surv(time, event) ~ group, data = persons),
               times = (0:(max(periods) + 1)) * len)
S_t <- tibble(group = sub("^group=", "", km$strata), period = round(km$time / len),
              S = km$surv, S_lo = km$lower, S_hi = km$upper)

p_death <- S_t %>% group_by(group) %>% arrange(period) %>%
  mutate(p    = S - lead(S),
         p_lo = pmax(0, S_lo - lead(S_hi)),
         p_hi = pmax(0, S_hi - lead(S_lo))) %>%
  ungroup() %>% filter(period %in% periods)

# Survivor component
alive <- alive_panel %>% group_by(group, period) %>%
  summarise(m = mean(cost), se = sd(cost) / sqrt(n()), .groups = "drop") %>%
  left_join(S_t, by = c("group", "period")) %>%
  mutate(alive = m * S, alive_lo = pmax(0, m - 1.96 * se) * S_lo,
         alive_hi = (m + 1.96 * se) * S_hi)

# Death-period component
dc <- death_costs %>% group_by(group) %>%
  summarise(e = mean(death_cost), se = sd(death_cost) / sqrt(n()), .groups = "drop")
death <- p_death %>% left_join(dc, by = "group") %>%
  mutate(dth = p * e, dth_lo = p_lo * pmax(0, e - 1.96 * se), dth_hi = p_hi * (e + 1.96 * se))

expected <- alive %>%
  left_join(dplyr::select(death, group, period, dth, dth_lo, dth_hi), by = c("group", "period")) %>%
  mutate(expected = alive + dth)

# Cumulative difference (first group level minus second) with Monte Carlo CI.
# Components are treated as independent across periods, so the CI is approximate.
groups <- sort(unique(expected$group))
draw <- function(g) {
  d <- filter(expected, group == g)
  replicate(R, sum(pmax(0, rnorm(nrow(d), d$alive, (d$alive_hi - d$alive_lo) / 3.92)) +
                     pmax(0, rnorm(nrow(d), d$dth,   (d$dth_hi   - d$dth_lo)   / 3.92))))
}
totals <- expected %>% group_by(group) %>% summarise(total = sum(expected))
mc <- draw(groups[1]) - draw(groups[2])
est <- totals$total[totals$group == groups[1]] - totals$total[totals$group == groups[2]]
c(estimate = est, quantile(mc, c(.025, .975)))




