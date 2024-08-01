library(tidyverse)
library(survival)
library(kableExtra)

# density and cdf of Weibull distribution left truncated at 28 days
dtwei28 <- function(x, shape, scale, log = FALSE) {
  dweibull(x, shape, scale, log)/
    pweibull(28, shape, scale, lower=FALSE)
}
ptwei28 <- function(q, shape, scale , log = FALSE) {
  pweibull(q, shape, scale, log)/
    pweibull(28, shape, scale, lower=FALSE)
}

# Reparametrize from R [d/p/q]weibull() functions to log hazard scale
phi <- function(v) {names(v) <- NULL; a=v[1]; b=v[2]; log(a) - a*log(b)}
del_phi <- function(v) {names(v) <- NULL; a=v[1]; b=v[2]; c(1/a - log(b), -a/b)}

# Summary statistics from maximum likelihood estimation by fitdistrplus::fitdistcens()
# calculate log hazard scale parameter and standard error by delta method
extract_weibull_stats <- function(distfit) {
  stopifnot("fitdistcens" %in% class(distfit))

  eta_hat <- phi(distfit$est)
  del_phi_hat <- del_phi(distfit$est)

  var_eta_hat <- (1/distfit$n) * t(del_phi_hat) %*% distfit$vcov %*% del_phi_hat
  se_eta_hat <- sqrt(var_eta_hat[1,1])

  tibble(param = c(names(distfit$estimate), "loghaz"),
         est = c(distfit$estimate, eta_hat),
         sd = c(distfit$sd, se_eta_hat),
         n = distfit$n)
}

# fit the Weibull state transition model for sojourn times in a given line of therapy
estimate_1line_2events <- function(surv_df, k){
  stopifnot(all(surv_df$os_entrydate < surv_df$startdate))

  lot_df <- surv_df %>%
    filter(linenumber == k & sojourn_time > 0) %>%
    mutate(event_factor = fct_drop(event_factor)) %>%
    mutate(event_death = as.numeric(event_factor == "death"))

  e1_fitdf <- lot_df %>%
    mutate(left = sojourn_time,
           right = case_when(event_death == 1 ~ sojourn_time)) %>%
    dplyr::select(left, right) %>%
    as.data.frame()

  e2_fitdf <- lot_df %>%
    mutate(left = sojourn_time,
           right = case_when(event_death == 0 ~ sojourn_time)) %>%
    dplyr::select(left, right) %>%
    filter(left > 28) %>%
    as.data.frame()

  e1 <- fitdistrplus::fitdistcens(e1_fitdf, dist = "weibull")
  e2 <- tryCatch(fitdistrplus::fitdistcens(e2_fitdf, dist = "twei28",
                                           start = list(shape=1 , scale=600)),
                 error = function(e) NA)

  lapply(list(e1, e2)[which(!is.na(list(e1, e2)))], extract_weibull_stats) %>%
    bind_rows %>%
    mutate(event = rep(c("death", "newline")[which(!is.na(list(e1, e2)))], each = 3) )
}

# wrapper function per disease, creates analysis data and iterates over line numbers with sufficient n
weibull_lines <- function(disease, verbose = T, maxK = 12){
  lot_survival <- cgdb_prep(disease, build_quarter = "q122")

  # find maximum line number where both events and censoring are observed
  tmp_df <- lot_survival %>%
    filter(os_entrydate < startdate) %>%
    count(linenumber, event_fct2) %>%
    spread(event_fct2, n) %>%
    mutate(all_events = death > 5 & newline > 5) %>%
    slice_min(all_events)

  c_max_lines <- min(min(tmp_df$linenumber) - 1, maxK)

  weibull_params <- lapply(1:c_max_lines, function(k) {
    estimate_1line_2events(lot_survival %>% filter(os_entrydate < startdate), k)}) %>%
    bind_rows(.id = "linenumber") %>%
    mutate(disease = disease)
}

plot_weibull_lines <- function(params_df) {
  maxK <- max(params_df$linenumber)
  n_big_line <- max(params_df$n)
  disease <- unique(params_df$disease)

  params_df %>% filter(n >= 10) %>%
    mutate(linenumber = factor(linenumber, levels = 1:maxK)) %>%
    ggplot(aes(x = linenumber, y = est, color = event)) + geom_point() +
    facet_wrap(vars(param), scales = "free",
               labeller = as_labeller(c(loghaz = "log hazard",
                                        scale = "Weibull scale",
                                        shape = "Weibull shape"))) +
    labs(title = paste0("Multi-state transitions in ", disease),
         subtitle = paste0("n = ", n_big_line))
}

#----------------
n_iter = 500

load(file = here::here("results/breast_params16.Rda"))
J <- length(unique(sim_params$linenumber))

#' generate random numbers with truncated Weibull distribution
#' @param n number of observations.
#' @param shape shape parameter
#' @param scale scale parameter, defaulting to 1.
#' @param lower_bound vector of left truncation limits.
#' @returns length-n vector of iid T | T > lower_bound where T ~ Weibull(shape, scale)
rtruncweibull <- function(n, shape, scale, lower_bound = 0) {
  p_thresh <- pweibull(lower_bound, shape, scale)
  u_trunc <- runif(n = n, min = p_thresh, max = 1)
  qweibull(p = u_trunc, shape = shape, scale = scale)
}

#' Transform between alternative parameterizations of the Weibull distribution
#' @param a shape parameter, on the log hazard scale
#' @param theta log-hazard scale parameter
#' @returns scale parameter `b` as defined for R functions `[dpqr]weibull()`.
loghaz_rweibull_scale <- function(a, theta) {
  exp((log(a)-theta)/a)
}

#' Simulate event times following possibly-truncated Weibull distribution with covariate effects on the log-hazard scale
#' @param n number of observations to simulate
#' @param shape numeric > 0, shape parameter of Weibull distribution.
#' This determines the baseline hazard as shape as h_0(t) = t^(shape-1)
#' @param loghaz_scale function of covariates defining effects on the log-hazard scale.
#' Overall hazard is thus a function of loghaz_effects and the baseline hazard (via \code{shape}).
#' @param X data frame containing all covariates used in the function loghaz_scale()
#' @param entry a non-negative number or length-n vector thereof, giving lower bound for truncated Weibull variates
#' @returns a numeric vector of length n containing simulated event times
rweibull_cond <- function(n, shape, loghaz_scale, entry = 0, X) {
  stopifnot(n==nrow(X))
  stopifnot(all(names(formals(loghaz_scale)) %in% names(X)))

  b <- loghaz_rweibull_scale(
    a = shape,
    theta = do.call(loghaz_scale,
                    X[,which(names(X) %in% names(formals(loghaz_scale))), drop = F]))

  rtruncweibull(n = n, shape = shape, scale = b, lower_bound = entry)
}

#' create a survival analysis estimation procedure object
#' @param coxph_formula a formula object to be used by survival::coxph()
#' @param index_rule a function to define index dates
#' @param iptw_model NULL for an unweighted model, or a function to fit propensity scores
#' @param iptw_trunc numeric, maximum value for inverse propensity weights, default is 10.
#' @returns a function of class `estiproc` that takes one cgdb-like data frame
#' and produces a fitted `coxph` object.
#'
estiproc <- function(coxph_formula, index_rule = NULL,
                     iptw_model = NULL, iptw_trunc = 10) {
  tmp <- as.list(environment())
  validate_estiproc_args(tmp)

  proc <- function(df) {
    if (!is.null(index_rule)) {
      df <- index_rule(df)
    }
    if (!is.null(iptw_model)) {
      ps_mod <- iptw_model(df)
      stopifnot(is(ps_mod, "glm"))

      pi_hat <- predict(ps_mod, type = "response")
      stopifnot(length(pi_hat) == nrow(df))
      stopifnot(is.numeric(iptw_trunc) & iptw_trunc > 0)

      df$iptw_hat <- pmin(iptw_trunc,
                          as.numeric(df$A1 == 1)/pi_hat + as.numeric(df$A1 == 0)/(1-pi_hat))

    }
    if (length(unique(df$patientid)) < nrow(df)) {
      if (is.null(iptw_model)) {
        survival::coxph(coxph_formula, cluster = patientid, data = df)
      } else {
        survival::coxph(coxph_formula, cluster = patientid, weights = iptw_hat, data = df)
      }
    } else {
      if (is.null(iptw_model)) {
        survival::coxph(coxph_formula, data = df)
      } else {
        survival::coxph(coxph_formula, weights = iptw_hat, data = df)
      }
    }
  }

  class(proc) <- c("estiproc", class(proc))
  proc
}

#' check that estimation procedure elements are basically the right format
validate_estiproc_args <- function(proc) {
  stopifnot(is(proc$coxph_formula, "formula"))
  stopifnot(is.null(proc$index_rule) | is(proc$index_rule, "function"))
  stopifnot(is.null(proc$iptw_model) | is(proc$iptw_model, "function"))
}

## return point estimates and confidence intervals from a fitted Cox model
extract_A1_stats <- function(coxmod, beta) {
  stopifnot(is(coxmod, "coxph"))

  beta_ci <- confint(coxmod)["A1",]
  cover <- beta > beta_ci[1] & beta < beta_ci[2]
  names(cover) <- "cover95"

  c(coef(coxmod)["A1"], cover)
}

#' evaluate a 2-layer-deep list of estimation procedures,
estimate_A1_stats <- function(esti_ls, cohort_df, sim_beta) {
  fits_ls <- esti_ls %>%
    map(function(ls) {
      map(ls, function(f) f(cohort_df))
    })

  map(fits_ls, function(ls) {
    map(ls, function(x) extract_A1_stats(x, beta = sim_beta))
  }) %>%
    map(unlist) %>%
    unlist
}

## return point estimates and confidence intervals from a fitted Cox model
extract_int_stats <- function(coxmod, beta) {
  stopifnot(is(coxmod, "coxph"))

  beta_ci <- confint(coxmod)["A1:Z",]
  cover <- beta > beta_ci[1] & beta < beta_ci[2]
  names(cover) <- "cover95"

  c(coef(coxmod)["A1:Z"], cover)
}

#' evaluate a 2-layer-deep list of estimation procedures,
estimate_int_stats <- function(esti_ls, cohort_df, sim_beta) {
  fits_ls <- esti_ls %>%
    map(function(ls) {
      map(ls, function(f) f(cohort_df))
    })

  map(fits_ls, function(ls) {
    map(ls, function(x) extract_int_stats(x, beta = sim_beta))
  }) %>%
    map(unlist) %>%
    unlist
}

#' check format of elements of baseline dgp
validate_baseline_dgp <- function(dgp) {
  stopifnot(is.numeric(dgp$n))
}

#' create object with baseline_dgp class
baseline_dgp <- function(n, sigma = NULL, p_responder = NULL,
                         frailty_rng = NULL, responder_rng = NULL) {
  tmp <- as.list(environment())
  validate_baseline_dgp(tmp)
  attr(tmp, "class") <- c("baseline_dgp", class(tmp))
  tmp
}

#' create a simulated data frame with specified dgp
baseline_sim <- function(dgp) {
  stopifnot("baseline_dgp" %in% class(dgp))

  df <- data.frame(patientid = 1:dgp$n)
  if (!is.null(dgp$frailty_rng)) {
    df$frailty <- do.call(dgp$frailty_rng, dgp[
      which(names(dgp) %in% names(formals(dgp$frailty_rng))), drop = F])
  }
  if (!is.null(dgp$responder_rng)) {
    df$Z <- do.call(dgp$responder_rng, dgp[
      which(names(dgp) %in% names(formals(dgp$responder_rng))), drop = F])
  }
  df
}

#' check the format of elements of candidate multistate dgp object
validate_multistate_dgp <- function(dgp) {
  stopifnot(is.data.frame(dgp$X))
  stopifnot(all(names(formals(dgp$Y1_rng)) %in%
                  c("n", "A1", "everA1", names(dgp$X), names(dgp$line_params))))
  stopifnot(is.null(dgp$Y2_rng) |
              all(names(formals(dgp$Y2_rng)) %in%
                    c("n", "A1", "everA1", names(dgp$X), names(dgp$line_params))))
}

#' creates S3 object of multistate DGP class
multistate_dgp <- function(X, line_params = NULL,
                           A1_rng = NULL, A2_rng = NULL,
                           Y1_rng = NULL, Y2_rng = NULL) {
  tmp <- as.list(environment())
  validate_multistate_dgp(tmp)
  attr(tmp, "class") <- c("multistate_dgp", class(tmp))
  tmp
}

#' create a simulated data frame with specified dgp
multistate_sim <- function(dgp, verbose = F) {
  stopifnot("multistate_dgp" %in% class(dgp))

  df <- dgp$X
  if (!is.null(dgp$line_params)) {
    if ("linenumber" %in% names(df)) {
      df <- full_join(df, dgp$line_params)
    } else df <- full_join(df, dgp$line_params, by = character())
  }
  n <- nrow(df)

  if (!is.null(dgp$A1_rng) & !is.null(dgp$A2_rng)) {
    df$A1 <- do.call(dgp$A1_rng,
                     c(n, df[,which(names(df) %in% names(formals(dgp$A1_rng))), drop = F]))
    df$A2 <- do.call(dgp$A2_rng,
                     c(n, df[,which(names(df) %in% names(formals(dgp$A2_rng))), drop = F]))
    df <- df %>%
      group_by(patientid) %>%
      arrange(patientid, linenumber) %>%
      mutate(everA1 = as.numeric(cumsum(A1) > 0),
             everA2 = as.numeric(cumsum(A2) > 0))
  }

  df$Y1 <- do.call(dgp$Y1_rng,
                   c(n, df[,which(names(df) %in% names(formals(dgp$Y1_rng))), drop = F]))
  df$Y2 <- do.call(dgp$Y1_rng,
                   c(n, df[,which(names(df) %in% names(formals(dgp$Y1_rng))), drop = F]))

  if (verbose) print(names(df))
  cols_to_remove <- names(df)[grepl("alpha|eta|tx|Y|frailty", names(df))]

  df <- df %>% group_by(patientid) %>%
    mutate(Tj = pmin(Y1, Y2),
           delta1 = as.numeric(Tj == Y1),
           delta2 = as.numeric(Tj == Y2)) %>%
    mutate(dead = cumsum(replace_na(lag(delta1), 0))) %>%
    filter(dead == 0) %>%
    arrange(desc(linenumber), .by_group = TRUE) %>%
    mutate(OS = cumsum(Tj)) %>%
    #mutate(OS = rev(cumsum(rev(Tj)))) %>%
    dplyr::select(!all_of(cols_to_remove), -dead) %>%
    ungroup() %>%
    mutate(sojourn_time = Tj,
           event_factor = factor(delta1, labels = c("death", "newline")))

  dgp$X <- df
  class(dgp) <- c("multistate_sim", class(dgp))
  dgp
}

#' create a multistate dgp with persistent treatment effects from a data frame of parameters
persistent_dgp <- function(bl_dgp, params_df) {
  multistate_dgp(X = baseline_sim(bl_dgp),
                 line_params = params_df, # danger zone, the value of beta_A1 is set by default argument
                 A1_rng = function(n, tx1) rbinom(n = n, prob = tx1, size = 1),
                 A2_rng = function(n, tx2) rbinom(n = n, prob = tx2, size = 1),
                 Y1_rng = function(n, alpha_death, eta_death, beta_A1, frailty, everA1) {
                   rweibull_cond(n = n, shape = alpha_death,
                                 loghaz_scale = function(eta_death, beta_A1, frailty, everA1) {
                                   eta_death + frailty + beta_A1*everA1},
                                 X = data.frame(alpha_death, eta_death, beta_A1, frailty, everA1)) },
                 Y2_rng = function(n, alpha_newline, eta_newline, beta_A1, frailty, everA1) {
                   rweibull_cond(n = n, shape = alpha_newline,
                                 loghaz_scale = function(eta_newline, beta_A1, frailty, everA1) {
                                   eta_newline + frailty + beta_A1*everA1},
                                 X = data.frame(alpha_newline, eta_newline, beta_A1, frailty, everA1)) })
}

#' create a multistate dgp with heterogeneous treatment effects from a data frame of parameters
heterogeneous_dgp <- function(bl_dgp, params_df) {
  multistate_dgp(X = baseline_sim(bl_dgp), line_params = params_df,
                 A1_rng = function(n, tx1) rbinom(n = n, prob = tx1, size = 1),
                 A2_rng = function(n, tx2) rbinom(n = n, prob = tx2, size = 1),
                 Y1_rng = function(n, alpha_death, eta_death, beta_A1, frailty, everA1, Z) {
                   rweibull_cond(n = n, shape = alpha_death,
                                 loghaz_scale = function(eta_death, beta_A1, frailty, everA1, Z) {
                                   eta_death + frailty + Z*beta_A1*everA1},
                                 X = data.frame(alpha_death, eta_death, beta_A1, frailty, everA1, Z)) },
                 Y2_rng = function(n, alpha_newline, eta_newline, beta_A1, frailty, everA1, Z) {
                   rweibull_cond(n = n, shape = alpha_newline,
                                 loghaz_scale = function(eta_newline, beta_A1, frailty, everA1, Z) {
                                   eta_newline + frailty + Z*beta_A1*everA1},
                                 X = data.frame(alpha_newline, eta_newline, beta_A1, frailty, everA1, Z)) })
}

#' a list of functions, each taking a parameters data frame and returning a multistate dgp
dgp_fct_list <- list(persistent_dgp, heterogeneous_dgp)
names(dgp_fct_list) <- c("persistent", "heterogeneous")

tx_trend <- data.frame(linenumber = 1:J,
                       tx1 = seq(0, .1, length = J),
                       tx2 = rep(0.2, length = J))

sim_params <- sim_params %>%
  mutate(param = case_when(param == "shape" ~ "alpha",
                           param == "scale" ~ "eta")) %>%
  pivot_wider(names_from = c(param, event), values_from = est)


index_rules_list <- list(
  all = function(df) filter(df, A1 + A2 > 0),
  first = function(df) {
    df %>%
      group_by(patientid) %>%
      arrange(linenumber) %>%
      filter(A1 + A2 > 0) %>%
      slice_min(linenumber) %>%
      ungroup()
  },
  random = function(df) {
    treated <- df %>%
      filter(A1 == 1) %>%
      group_by(patientid) %>%
      slice_sample(n=1) %>%
      ungroup()
    comparator <- df %>%
      filter(!(patientid %in% treated$patientid) & A2 == 1) %>%
      group_by(patientid) %>%
      slice_sample(n=1) %>%
      ungroup()
    rbind(treated, comparator)
  })


## --------------------------------------------------------------
## Estimation procedures for main effects

# formulas for survival::coxph()
models_list1 <- c(m0 = "Surv(OS) ~ A1",
                  m1 = "Surv(OS) ~ A1 + linenumber",
                  m2 = "Surv(OS) ~ A1 + strata(linenumber)")

# unweighted estimation models
esti_list1a <- map(index_rules_list, function(rule) {
  map(models_list1, function(model) {
    estiproc(as.formula(model), index_rule = rule)
  })
})

# IPTW estimation models
esti_list1b <- map(index_rules_list, function(rule) {
  inner_list <- map(models_list1, function(model) {
    estiproc(as.formula(model), index_rule = rule,
             iptw_model = function(df) glm(A1 ~ linenumber, family = binomial, data = df))
  })
  names(inner_list) <- paste(names(models_list1), "iptw", sep = ".")
  inner_list
})



## --------------------------------------------------------------

#' tibble data frame enumerating simulation scenarios
sims_ctrl <- tibble(expand.grid(dgp = "persistent",
                                tx_trend = "linear",
                                beta = -2:2, #-2:2
                                n = 100*5^{1:3},
                                seed = 1:n_iter,
                                stats1 = list(NULL),
                                stats2 = list(NULL)))

for (i in 1:nrow(sims_ctrl)) {
  print(i); print(Sys.time())
  set.seed(sims_ctrl$seed[i])

  params_df <- sim_params %>%
    mutate(beta_A1 = sims_ctrl$beta[i]) %>%
    inner_join(txtrend_ls[[sims_ctrl$tx_trend[i]]])
  sim_bl <- baseline_dgp(n = sims_ctrl$n[i],
                         sigma = 0.4,
                         p_responder = 1,
                         frailty_rng = function(n, sigma) {
                           rnorm(n = n, mean = 0, sd = sigma)},
                         responder_rng = function(n, p_responder) {
                           rbinom(n = n, prob = p_responder, size = 1)
                         })

  sim_dgp <- dgp_fct_list[[sims_ctrl$dgp[i]]](sim_bl, params_df)
  sim_df <- multistate_sim(sim_dgp)$X

  sims_ctrl$stats1[[i]] <- estimate_A1_stats(
    esti_list1a, sim_df, sims_ctrl$beta[i])

  sims_ctrl$stats2[[i]] <- estimate_A1_stats(
    esti_list1b, sim_df, sims_ctrl$beta[i])
}

## --------------------------------------------------------------
## Estimation procedures for interaction effects

# formulas for survival::coxph()
models_list2 <- c(m0 = "Surv(OS) ~ A1*Z",
                  m1 = "Surv(OS) ~ A1*Z + linenumber",
                  m2 = "Surv(OS) ~ A1*Z + strata(linenumber)")

# unweighted estimation models
esti_list2a <- map(index_rules_list, function(rule) {
  map(models_list2, function(model) {
    estiproc(as.formula(model), index_rule = rule)
  })
})

# IPTW estimation models
esti_list2b <- map(index_rules_list, function(rule) {
  inner_list <- map(models_list2, function(model) {
    estiproc(as.formula(model), index_rule = rule,
             iptw_model = function(df) glm(A1 ~ linenumber,
                                           family = binomial, data = df))
  })
  names(inner_list) <- paste(names(models_list2), "iptw", sep = ".")
  inner_list
})

#' tibble data frame enumerating simulation scenarios
sims_ctrl <- tibble(expand.grid(dgp = "heterogeneous",
                                tx_trend = "linear",
                                n = 100*5^{1:3},
                                beta = -2:2, #-2:2
                                p_Z = c(0.1, 0.5),
                                seed = 100+1:n_iter,
                                stats2a = list(NULL),
                                stats2b = list(NULL)))

for (i in 1:nrow(sims_ctrl)) {
  print(i); print(Sys.time())
  set.seed(sims_ctrl$seed[i])

  params_df <- sim_params %>%
    mutate(beta_A1 = sims_ctrl$beta[i]) %>%
    inner_join(txtrend_ls[[sims_ctrl$tx_trend[i]]])

  sim_bl <- baseline_dgp(n = sims_ctrl$n[i],
                         sigma = 0.4,
                         p_responder = sims_ctrl$p_Z[i],
                         frailty_rng = function(n, sigma) {
                           rnorm(n = n, mean = 0, sd = sigma)},
                         responder_rng = function(n, p_responder) {
                           rbinom(n = n, prob = p_responder, size = 1)
                         })

  sim_dgp <- dgp_fct_list[[sims_ctrl$dgp[i]]](sim_bl, params_df)
  sim_df <- multistate_sim(sim_dgp)$X

  sims_ctrl$stats2a[[i]] <- estimate_int_stats(
    esti_list2a, sim_df, sims_ctrl$beta[i])

  sims_ctrl$stats2b[[i]] <- estimate_int_stats(
    esti_list2b, sim_df, sims_ctrl$beta[i])
}
