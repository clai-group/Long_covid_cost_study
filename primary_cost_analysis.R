library(dplyr)
library(geepack)
library(broom)

# Input: analysis_data (quarterly panel, analysis period)
#        pre_pandemic_data (quarterly panel, 2018-2019)

### PRIMARY: Two-part GEE 
part1_model <- geeglm(
  has_any_cost ~ LC * quarter_since_index +
    age_centered + I(age_centered^2) + sex + race + ethnicity +
    hospitalization + charlson_centered + log_baseline_cost,
  data = analysis_data, id = patient_num,
  family = binomial(link = "logit"), corstr = "exchangeable"
)

part2_model <- geeglm(
  quarterly_cost ~ LC * quarter_since_index +
    age_centered + I(age_centered^2) + sex + race + ethnicity +
    hospitalization + charlson_centered + log_baseline_cost,
  data = analysis_data %>% filter(quarterly_cost > 0), id = patient_num,
  family = Gamma(link = "log"), corstr = "exchangeable"
)

### SENSITIVITY: Parallel trends (pre-pandemic) 
parallel_trends_model <- geeglm(
  quarterly_cost ~ LC * quarter_num +
    age_centered + I(age_centered^2) + sex + race + ethnicity +
    hospitalization + charlson_centered + log_baseline_cost,
  data = pre_pandemic_data %>% filter(quarterly_cost > 0), id = patient_num,
  family = Gamma(link = "log"), corstr = "exchangeable"
)

### DECOMPOSITION
visit_model <- geeglm(
  visits ~ LC * quarter_since_index +
    age_centered + I(age_centered^2) + sex + race + ethnicity +
    hospitalization + charlson_centered + log_baseline_cost,
  data = analysis_data, id = patient_num,
  family = poisson(link = "log"), corstr = "exchangeable"
)

cost_per_visit_model <- geeglm(
  cost_per_visit ~ LC * quarter_since_index +
    age_centered + I(age_centered^2) + sex + race + ethnicity +
    hospitalization + charlson_centered + log_baseline_cost,
  data = analysis_data %>% filter(visits > 0, quarterly_cost > 0), id = patient_num,
  family = Gamma(link = "log"), corstr = "exchangeable"
)