# Long COVID healthcare cost study: modeling code

Modeling scripts for "Long COVID and Healthcare Costs and Utilization Among Adults
with Prior COVID-19 Infection." Data preparation is not included; each script lists
the analysis-ready inputs it expects.

- `primary_cost_analysis.R`: two-part GEE (logistic for any cost; gamma-log for cost
  given any cost), g-computation of the effect in the exposed group with
  parametric bootstrap CIs, sensitivity analyses (categorical time, alternative
  covariate set, restriction to complete follow-up, alternative working
  correlations), and visit-frequency / cost-per-visit models.
- `marginal_analysis.R`: piecewise two-part GEE over pre- and post-index
  periods with group-specific slopes and level change, g-computation, and
  normalization to the mean pre-index difference.
- `sensitivity_survival_analysis.R`: expected cost combining Kaplan-Meier survival,
  survivor costs, and costs in the period of death, with a Monte Carlo CI.

Requires R packages dplyr, tidyr, geepack, broom, MASS, survival.