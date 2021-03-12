## Herramientas de uso general

`%cat%` <- function(a, b) paste0(a, b)
## -----  ---

## Calcula muestras de la eficacia usando modelo bayesiano
## El modelo estima los casos con distribución de Poisson y
## prior para theta usando la distribución Beta con alpha
## estimado a partir de un hyperparámetro `ve_prior` para la eficacia,
## y beta = 1.
##
## Cuando ve_prior==0, que equivale a una distribución uniforme,
## equivale a no tener ninguna distribución anterior, pues los
## calculos de ésta se omiten
poisson_bayes <- function(key, data.src, ve_prior = 0.3) {
  library(rstan)
  vacc_data <- data.src(key)
  vacc_data[["ve_prior"]] <- ve_prior
  # browser()
  poisson_bayes_raw(vacc_data)
}

poisson_bayes_raw <- function(vacc_data, ret_samp = F) {
  library(rstan)
  options(mc.cores = parallel::detectCores())
  rstan_options(auto_write = TRUE)

  # Caso Sputnik
  if (is.na(vacc_data[["surv_time_vac"]]) & is.na(vacc_data[["surv_time_pla"]])) {
    vacc_data[["surv_time_vac"]] <- vacc_data[["size_vac"]]
    vacc_data[["surv_time_pla"]] <- vacc_data[["size_pla"]]
  }

  model <- stan(
    file = "poisson-ve.stan",
    # model_code = poisson_stan_model,
    data = vacc_data,
    seed = 12345,
    iter = 2500,
    refresh = 0,
    control = list(max_treedepth = 15, adapt_delta = 0.99)
  )

  samples <- as.data.frame(model)

  if (ret_samp == TRUE) {
    return(samples)
  } else {
    VE_mean <- mean(samples$VE)
    VE_median <- median(samples$VE)
    prob_ve30plus <- mean(samples$VE > 0.3)

    library(HDInterval)
    cred_int <- hdi(samples$VE, credMass = 0.95) %>% unname()

    c(
      CI_2.5 = cred_int[1], mean = VE_mean,
      CI_97.5 = cred_int[2], `Pr_30+` = prob_ve30plus, median = VE_median
    )
  }
}

vaccine_pharma <- c("Pfizer", "Moderna", "AstraZeneca", "Sputnik", "Sino1", "Sino2")

vaccine_all_data <- function() {
  tibble(
    pharma = vaccine_pharma,
    cases_vac = c(8, 11, 30, 16, 142, 188),
    cases_pla = c(162, 185, 101, 62, 211, 211),
    size_vac = c(18198, 14134, 5807, 14964, 2457, 2460),
    size_pla = c(18325, 14073, 5829, 4902, 2458, 2458),
    surv_time_vac = c(2.214, 3.274, 0.68, NA, 0.5477, 0.5433),
    surv_time_pla = c(2.222, 3.333, 0.677, NA, 0.5396, 0.5396),
    pharma_ve = c(95.0, 94.1, 70.4, 91.6, NA, NA),
    pharma_95_li = c(90.3, 89.3, 54.8, 85.6, NA, NA),
    pharma_95_ls = c(97.6, 96.8, 80.6, 95.2, NA, NA),
  )
}

sino_outcome <- c("Desat93", "TAC", "Hospital", "Fallec")
peru_outcome <- c("Sino1:" %cat% sino_outcome, "Sino2:" %cat% sino_outcome)

peru_sino_outcome_data <- function() {
  tibble(
    outcome = "Sino1:" %cat% sino_outcome,
    cases_vac = c(6, 34, 1, 0),
    cases_pla = c(16, 93, 11, 2),
    surv_time_vac = NA,
    surv_time_pla = NA,
    size_vac = c(142, 142, 142, 142),
    size_pla = c(211, 211, 211, 211)
  ) %>%
    bind_rows(
      tibble(
        outcome = "Sino2:" %cat% sino_outcome,
        cases_vac = c(17, 79, 9, 1),
        cases_pla = c(16, 93, 11, 2),
        surv_time_vac = NA,
        surv_time_pla = NA,
        size_vac = c(188, 188, 188, 188),
        size_pla = c(211, 211, 211, 211)
      )
    )
}

peru_outc_data <- function(key) {
  res <- peru_sino_outcome_data() %>%
    filter(outcome == key)

  as.list(res)
}

vaccine_data <- function(key) {
  res <- vaccine_all_data() %>% filter(pharma == key)
  as.list(res)
}

pois_raw_approx <- function(key, data.src) {
  vacc_data <- data.src(key)
  with(vacc_data, {
    VE <- (cases_vac / surv_time_vac) / (cases_pla / surv_time_pla)
    delta <- 1.96 * sqrt(1 / cases_vac + 1 / cases_pla)

    tibble_row(
      CI_2.5 = 1 - VE * exp(delta),
      mean = 1 - VE, CI_97.5 = 1 - VE * exp(-delta)
    )
  })
}

get_estim_all <- function(fun, keys, data.src, ...) {
  local_fun <- function(x) {
    fun(x, data.src, ...) %>% pcnt()
  }
  res <- sapply(keys, local_fun) %>%
    t() %>%
    as.data.frame()
  res <- add_column(res, key = keys, .before = 1)
  rownames(res) <- NULL
  res
}

pcnt <- function(x, digits = 2) {
  # round(x * 100, digits)
  round100 <- function(e) {
    if (is.numeric(e)) {
      round(e * 100, digits)
    } else {
      e
    }
  }
  sapply(x, round100)
}

theta_from_ve <- function(VE, r_t = 1) {
  return(r_t * (1 - VE) / (1 + r_t * (1 - VE)))
}

ve_from_theta <- function(theta, r_t = 1) {
  # calculate VE given case rate (theta) and surveillance time ratio
  return(1 + theta / (r_t * (theta - 1)))
}

## La misma función que binomial_gibbs, pero con otra interface
## Basada en https://boyangzhao.github.io/posts/vaccine_efficacy_bayesian
binomial_gibbs_raw <- function(vacc_data) {
  alpha_from_ve_prior <- function(ve_prior) {
    theta <- (1 - ve_prior) / (2 - ve_prior)
    theta / (1 - theta)
  }

  with(vacc_data, {

    # Cuando no hay surveillance time => usar group size
    if (is.na(surv_time_vac) & is.na(surv_time_pla)) {
      surv_time_vac <- size_vac
      surv_time_pla <- size_pla
    }

    alpha <- alpha_from_ve_prior(ve_prior) + cases_vac
    beta <- 1 + cases_pla
    irr_v <- cases_vac / surv_time_vac
    irr_p <- cases_pla / surv_time_pla
    r <- surv_time_vac / surv_time_pla

    VE <- 1 - irr_v / irr_p

    # confidence interval
    theta_ci_lower <- qbeta(0.025, alpha, beta)
    theta_ci_higher <- qbeta(0.975, alpha, beta)

    VE_ci_lower <- ve_from_theta(theta_ci_higher, r)
    VE_ci_upper <- ve_from_theta(theta_ci_lower, r)

    # Pr(VE>30%)
    prob_ve30plus <- pbeta(theta_from_ve(0.3, r), alpha, beta)


    c(CI_2.5 = VE_ci_lower, mean = VE, CI_97.5 = VE_ci_upper, `Pr_30+` = prob_ve30plus)
  })
}

## La misma función que binomial_gibbs_raw, pero con otra interface
binomial_gibbs <- function(key, data.src, ve_prior = 0.3) {
  vacc_data <- data.src(key)
  vacc_data[["ve_prior"]] <- ve_prior
  binomial_gibbs_raw(vacc_data)
}

pharma_ve_results <- function() {
  tibble()
}
