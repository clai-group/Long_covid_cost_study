library(dplyr)
library(tidyr)
library(survival)

# Inputs: quarterly_panel_censored, quarterly_panel_partitioned, patient_covariates

# KM survival by PASC status
km_fit <- survfit(Surv(followup_time, event) ~ status, data = patient_covariates)

S_by_group <- summary(km_fit, times = (0:20) * 0.25) %>%
  {tibble(time = .$time, S = .$surv, S_lo = .$lower, S_hi = .$upper,
          status = sub("^status=", "", .$strata))} %>%
  mutate(quarter_since_index = round(time / 0.25)) %>%
  select(status, quarter_since_index, S, S_lo, S_hi) %>%
  complete(status, quarter_since_index = 0:19,
           fill = list(S = 0, S_lo = 0, S_hi = 0))

# Per-quarter death probability: S(t) - S(t+1)
death_prob <- S_by_group %>%
  bind_rows(
    S_by_group %>% group_by(status) %>%
      summarise(quarter_since_index = 20, S = 0, S_lo = 0, S_hi = 0, .groups = "drop")
  ) %>%
  arrange(status, quarter_since_index) %>%
  group_by(status) %>%
  mutate(
    death_prob_q     = S - lead(S, default = 0),
    death_prob_q_lo  = pmax(0, S_lo - lead(S_hi, default = 0)),
    death_prob_q_hi  = pmax(0, S_hi - lead(S_lo, default = 0))
  ) %>%
  ungroup() %>%
  filter(quarter_since_index <= 19) %>%
  select(status, quarter_since_index, death_prob_q, death_prob_q_lo, death_prob_q_hi)

# Component 1: Mean survivor costs × S(t)
alive_component <- quarterly_panel_censored %>%
  group_by(status, quarter_since_index) %>%
  summarise(
    mean_cost      = mean(quarterly_cost, na.rm = TRUE),
    se_cost        = sd(quarterly_cost, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(
    mean_cost_lo = pmax(0, mean_cost - 1.96 * se_cost),
    mean_cost_hi = mean_cost + 1.96 * se_cost
  ) %>%
  left_join(S_by_group, by = c("status", "quarter_since_index")) %>%
  mutate(
    alive_comp     = mean_cost    * S,
    alive_comp_lo  = pmax(0, mean_cost_lo * S_lo),
    alive_comp_hi  = mean_cost_hi * S_hi
  ) %>%
  select(status, quarter_since_index, alive_comp, alive_comp_lo, alive_comp_hi)

# Component 2: Mean EOL costs × Pr(death in quarter t)
# EOL = total costs in final 2 quarters before death
eol_mean <- quarterly_panel_partitioned %>%
  filter(cost_category == "end_of_life", died) %>%
  group_by(patient_num) %>%
  summarise(eol_cost = sum(quarterly_cost, na.rm = TRUE), .groups = "drop") %>%
  left_join(patient_covariates %>% select(patient_num, status), by = "patient_num") %>%
  group_by(status) %>%
  summarise(
    mean_eol     = mean(eol_cost),
    se_eol       = sd(eol_cost) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(
    mean_eol_lo  = pmax(0, mean_eol - 1.96 * se_eol),
    mean_eol_hi  = mean_eol + 1.96 * se_eol
  )

eol_component <- death_prob %>%
  left_join(eol_mean, by = "status") %>%
  mutate(
    eol_comp     = death_prob_q    * mean_eol,
    eol_comp_lo  = pmax(0, death_prob_q_lo * mean_eol_lo),
    eol_comp_hi  = death_prob_q_hi * mean_eol_hi
  ) %>%
  select(status, quarter_since_index, eol_comp, eol_comp_lo, eol_comp_hi)

# Population expected cost = alive component + EOL component
pop_cost <- alive_component %>%
  left_join(eol_component, by = c("status", "quarter_since_index")) %>%
  mutate(
    exp_cost     = alive_comp    + eol_comp,
    exp_cost_lo  = alive_comp_lo + eol_comp_lo,
    exp_cost_hi  = alive_comp_hi + eol_comp_hi
  )