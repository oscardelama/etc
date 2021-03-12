data {
  int<lower=0> cases_vac; // Casos positivos en grupo de vacunados
  int<lower=0> cases_pla; // Casos positivos en grupo placebo
  real<lower=0> surv_time_vac; // Tiempo bajo vigilancia en grupo de vacunados
  real<lower=0> surv_time_pla; // Tiempo bajo vigilancia en grupo placebo
  real<lower=0> ve_prior; // Prior de Eficacia de vacuna
}
transformed data {
  real<lower=0> t_v;
  real<lower=0> t_p;
  real<lower=0> r_t;
  real<lower=0> alpha_prior;

  t_v = surv_time_vac * 1000;
  t_p = surv_time_pla * 1000;
  r_t = t_v / t_p;

  { real theta_prior;
    theta_prior = (1 - ve_prior) / (2 - ve_prior);
    alpha_prior = theta_prior / (1 - theta_prior);
  }
}
parameters {
  real<lower=0> ri_v; // Ratio de incidencia en vacunados
  real<lower=0> ri_p; // Ratio de incidencia en placebo
}
model {
  real ve;
  real theta;

  cases_vac ~ poisson(ri_v * t_v);
  cases_pla ~ poisson(ri_p * t_p);

  if (alpha_prior != 1) {
    ve = 1 - ri_v/ri_p;
    theta = r_t*(1-ve) / (1 + r_t*(1-ve));
    theta ~ beta(alpha_prior + ri_v * t_v, 1 + ri_p * t_p);
  }
}
generated quantities {
  real VE;
  VE = 1 - ri_v/ri_p;
}
