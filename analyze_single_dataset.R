library(plyr)
library(dplyr)
library(ggplot2)
library(gtools)
library(msm)
library(gam)
library(nnls)
library(ROCR)
library(cvAUC)
library(SuperLearner)


# helper function to estimate the surrogate index
estimate_g0 <- function(df, outcome, covariates, surrogate, learner = "glm"){
  
  # Ensure the dataset contains a weights column
  if (!"weights" %in% colnames(df)) {
    stop("The dataset must contain a 'weights' variable.")
  }
  
  # Ensure surrogate variable exists in the dataset
  if (!(surrogate %in% colnames(df))) {
    stop("The specified surrogate variable is not in the dataset.")
  }
  
  # Filter dataset for non-missing surrogate values (S_measured == 1)
  df_nonmissing <- df %>% filter(S_measured == 1)
  
  # Define formula dynamically, including the surrogate variable
  formula_str <- as.formula(paste(outcome, "~", paste(c(covariates, surrogate), collapse = " + ")))
  
  # Fit g0 model
  if (learner == "glm"){
    g0_mod <- glm(formula_str,
                  data = df_nonmissing,
                  family = "binomial",
                  weights = df_nonmissing$weights)
  } 
  else {
    g0_mod <- SuperLearner(Y = df_nonmissing[[outcome]], 
                           X = df_nonmissing %>% select(all_of(c(covariates, surrogate))), 
                           family = binomial(),
                           SL.library = c("SL.mean", "SL.glm", "SL.gam"), 
                           obsWeights = df_nonmissing$weights)
  } 
  
  return(g0_mod)
}


# helper function to estimate the outer expectation 
estimate_E_hat <- function(g0_mod, df_p3, treatment, covariates, learner = "glm") {
  
  # Predict Y using the g0 model
  if (learner == "glm") {
    df_p3$Y_pred <- predict(g0_mod, df_p3, type = "response")
  }
  else{
    df_p3$Y_pred <- predict(g0_mod, df_p3, type = "response")$pred[,1]
  } 
  
  # Filter based on treatment and surrogate measurement
  df_treated <- df_p3 %>% filter(get(treatment) == 1 & S_measured == 1)
  df_placebo <- df_p3 %>% filter(get(treatment) == 0 & S_measured == 1)
  
  if (learner == "glm"){
    formula_str <- as.formula(paste("Y_pred ~", paste(covariates, collapse = " + ")))
    
    E_hat_A0 <- glm(formula_str, data = df_placebo, weights = df_placebo$weights)
    E_hat_A1 <- glm(formula_str, data = df_treated, weights = df_treated$weights)
  } 
  else {
    # Use SuperLearner for flexible estimation
    E_hat_A0 <- SuperLearner(
      Y = df_placebo$Y_pred, 
      X = df_placebo %>% select(all_of(covariates)), 
      family = gaussian(),
      SL.library = c("SL.mean", "SL.glm", "SL.gam"), 
      obsWeights = df_placebo$weights
    )
    
    E_hat_A1 <- SuperLearner(
      Y = df_treated$Y_pred, 
      X = df_treated %>% select(all_of(covariates)), 
      family = gaussian(),
      SL.library = c("SL.mean", "SL.glm", "SL.gam"), 
      obsWeights = df_treated$weights
    )
  } 
  
  return(list(E_hat_A0 = E_hat_A0, E_hat_A1 = E_hat_A1))
}


# function to compute plugin estimates
compute_plug_in_estimates <- function(df_obs, df_p3, treatment, outcome, covariates, 
                                      surrogate, learner, ct_bias, uc_bias) {
  
  # Estimate g0 function
  g0_mod <- estimate_g0(df_obs, outcome, covariates, surrogate, learner)
  
  # Estimate E_hat functions
  E_hats <- estimate_E_hat(g0_mod, df_p3, treatment, covariates, learner)
  
  # Predict E0 and E1
  df_p3$E0_pred <- predict(E_hats$E_hat_A0, df_p3)
  df_p3$E1_pred <- predict(E_hats$E_hat_A1, df_p3)
  
  # Compute counterfactual means
  Y_counter_1 <- mean(df_p3$E1_pred) + ct_bias - uc_bias
  Y_counter_0 <- mean(df_p3$E0_pred) - uc_bias
  VE_estimate <- 1 - (Y_counter_1 / Y_counter_0)
  
  return(list(Y_1 = Y_counter_1, 
              Y_0 = Y_counter_0, 
              VE = VE_estimate, 
              log_1_minus_VE = log(1 - VE_estimate)))
}


# function to obtain bootstrap standard errors for the plugin estimator
bootstrap_se <- function(df_obs, df_p3, treatment, outcome, covariates, surrogate, 
                         learner, ct_bias, uc_bias, bootstrap_reps) {
  
  Y_1_boot <- c()
  Y_0_boot <- c()
  log_1_minus_VE_boot <- c()
  
  for (i in 1:bootstrap_reps) {
    
    # resample with replacement
    df_obs_boot <- df_obs[sample(nrow(df_obs), replace = TRUE), ]
    df_p3_boot <- df_p3[sample(nrow(df_p3), replace = TRUE), ]
    
    # compute VE estimates using the new function
    boot_results <- compute_plug_in_estimates(df_obs_boot, df_p3_boot, treatment, 
                                              outcome, covariates, surrogate, learner, 
                                              ct_bias, uc_bias)
    
    # store bootstrap estimates
    Y_1_boot[i] <- boot_results$Y_1
    Y_0_boot[i] <- boot_results$Y_0
    log_1_minus_VE_boot[i] <- boot_results$log_1_minus_VE
  }
  
  # compute standard errors
  Y_1_se <- sd(Y_1_boot)
  Y_0_se <- sd(Y_0_boot)
  log_1_minus_VE_se <- sd(log_1_minus_VE_boot)
  
  # store covariance matrix for tipping point analysis
  cov_mat <- matrix(cov(data.frame(x = Y_0_boot, y =Y_1_boot)), nrow = 2)
  
  return(list(Y_1_se = Y_1_se, 
              Y_0_se = Y_0_se, 
              log_1_minus_VE_se = log_1_minus_VE_se, 
              cov_mat = cov_mat))
}


# helper function to estimate sampling probabilities 
estimate_P_S_measured_ZAX <- function(df_obs, df_p3, treatment, covariates) {
  
  # Construct formula strings dynamically
  formula_Z0 <- as.formula(paste("S_measured ~", paste(c(covariates, treatment), collapse = " + ")))
  formula_Z1 <- as.formula(paste("S_measured ~", paste(covariates, collapse = " + ")))
  
  # Fit logistic regression models
  log_model_Z0 <- glm(formula_Z0, data = df_p3, family = binomial)
  log_model_Z1 <- glm(formula_Z1, data = df_obs, family = binomial)
  
  return(list(
    log_model_Z0 = log_model_Z0, 
    log_model_Z1 = log_model_Z1
  ))
}

# helper function to get additional nuisance estimates for the one-step estimator 
estimate_P_Z_A_XS <- function(df_obs, df_p3, treatment, outcome, covariates, 
                              surrogate, learner = "glm") {
  
  sampling_mod <- estimate_P_S_measured_ZAX(df_obs, df_p3, treatment, covariates)
  
  log_model_Z0 <- sampling_mod[["log_model_Z0"]]
  log_model_Z1 <- sampling_mod[["log_model_Z1"]]
  
  df_obs$fitted_probs <- predict(log_model_Z1, df_obs, type = "response")
  df_obs$weights <- 1 / df_obs$fitted_probs 
  
  df_p3$fitted_probs <- predict(log_model_Z0, df_p3, type = "response")
  df_p3$weights <- 1 / df_p3$fitted_probs 
  
  df_combined <- bind_rows(df_obs, df_p3)
  
  df_nonmissing <- df_combined %>% 
    filter(S_measured == 1) %>%
    mutate(Z0A1 = ifelse(Z == 0 & get(treatment) == 1, 1, 0), 
           Z0A0 = ifelse(Z == 0 & get(treatment) == 0, 1, 0))
  
  model_covariates <- c(covariates, surrogate)
  
  if (learner == "glm") {
    log_model_Z0A1 <- glm(
      Z0A1 ~ ., 
      data = df_nonmissing %>% select(Z0A1, all_of(model_covariates)), 
      family = "binomial",
      weights = df_nonmissing$weights
    )
    
    log_model_Z0A0 <- glm(
      Z0A0 ~ ., 
      data = df_nonmissing %>% select(Z0A0, all_of(model_covariates)), 
      family = "binomial",
      weights = df_nonmissing$weights
    )
    
    log_model_Z1 <- glm(
      Z ~ ., 
      data = df_nonmissing %>% select(Z, all_of(model_covariates)), 
      family = "binomial",
      weights = df_nonmissing$weights
    )
  }
  else {
    log_model_Z0A1 <- SuperLearner(
      Y = df_nonmissing$Z0A1, 
      X = df_nonmissing %>% select(all_of(model_covariates)), 
      family = binomial(),
      SL.library = c("SL.mean", "SL.glm", "SL.gam"), 
      obsWeights = df_nonmissing$weights
    )
    
    log_model_Z0A0 <- SuperLearner(
      Y = df_nonmissing$Z0A0, 
      X = df_nonmissing %>% select(all_of(model_covariates)), 
      family = binomial(),
      SL.library = c("SL.mean", "SL.glm", "SL.gam"), 
      obsWeights = df_nonmissing$weights
    )
    
    log_model_Z1 <- SuperLearner(
      Y = df_nonmissing$Z, 
      X = df_nonmissing %>% select(all_of(model_covariates)), 
      family = binomial(),
      SL.library = c("SL.mean", "SL.glm", "SL.gam"), 
      obsWeights = df_nonmissing$weights
    )
    
  } 
  
  return(list(
    log_model_Z0A1 = log_model_Z0A1, 
    log_model_Z0A0 = log_model_Z0A0,
    log_model_Z1 = log_model_Z1
  ))
}



# function to run the one-step estimator 
compute_onestep_estimates <- function(df_obs, df_p3, treatment, outcome, covariates, surrogate, 
                                      learner = "glm", ct_bias = 0, uc_bias = 0, alpha = 0.05) {
  
  # estimate P(Z, A | X, S)
  log_mod_P_ZA_XS <- estimate_P_Z_A_XS(df_obs, df_p3, treatment, outcome, covariates, surrogate, learner)
  
  # predict conditional probabilities based on learner type
  if (learner == "glm") {
    df_obs$P_Z0A1_pred <- predict(log_mod_P_ZA_XS[["log_model_Z0A1"]], df_obs, type = "response")
    df_obs$P_Z0A0_pred <- predict(log_mod_P_ZA_XS[["log_model_Z0A0"]], df_obs, type = "response")
    df_obs$P_Z1A0_pred <- predict(log_mod_P_ZA_XS[["log_model_Z1"]], df_obs, type = "response")
  } else {
    df_obs$P_Z0A1_pred <- predict(log_mod_P_ZA_XS[["log_model_Z0A1"]], df_obs, type = "response")$pred[,1]
    df_obs$P_Z0A0_pred <- predict(log_mod_P_ZA_XS[["log_model_Z0A0"]], df_obs, type = "response")$pred[,1]
    df_obs$P_Z1A0_pred <- predict(log_mod_P_ZA_XS[["log_model_Z1"]], df_obs, type = "response")$pred[,1]
  }
  
  # estimate g0 and E_hat models
  g0_mod <- estimate_g0(df_obs, outcome, covariates, surrogate, learner)
  E_hats <- estimate_E_hat(g0_mod, df_p3, treatment, covariates, learner)
  
  # predict outcomes based on learner type
  if (learner == "glm") {
    df_p3$g_pred <- predict(g0_mod, df_p3, type = "response")
    df_obs$g_pred <- predict(g0_mod, df_obs, type = "response")
    
    df_p3$E0_pred <- predict(E_hats$E_hat_A0, df_p3, type = "response")
    df_p3$E1_pred <- predict(E_hats$E_hat_A1, df_p3, type = "response")
  } else {
    df_p3$g_pred <- predict(g0_mod, df_p3, type = "response")$pred[,1]
    df_obs$g_pred <- predict(g0_mod, df_obs, type = "response")$pred[,1]
    
    df_p3$E0_pred <- predict(E_hats$E_hat_A0, df_p3, type = "response")$pred[,1]
    df_p3$E1_pred <- predict(E_hats$E_hat_A1, df_p3, type = "response")$pred[,1]
  }
  
  # counterfactual means
  Y_counter_1 <- mean(df_p3$E1_pred)
  Y_counter_0 <- mean(df_p3$E0_pred)
  
  # estimate P(Z = 0) and P(A = a | X, Z = 0)
  est_P_Z_0 <- nrow(df_p3) / (nrow(df_obs) + nrow(df_p3))
  est_P_A <- mean(df_p3[[treatment]])
  
  # compute influence function corrections
  df_obs <- df_obs %>% 
    mutate(
      phi0 = (1 / est_P_Z_0) * (1 / est_P_A) * (P_Z0A0_pred / P_Z1A0_pred) * (Y - g_pred),
      phi1 = (1 / est_P_Z_0) * (1 / est_P_A) * (P_Z0A1_pred / P_Z1A0_pred) * (Y - g_pred)
    )
  
  df_p3 <- df_p3 %>% 
    mutate(
      phi0 = ((1 / est_P_Z_0) * ((1 - get(treatment)) / (1 - est_P_A)) * (g_pred - E0_pred)) + 
        ((1 / est_P_Z_0) * (E0_pred - Y_counter_0)),
      phi1 = ((1 / est_P_Z_0) * (get(treatment) / est_P_A) * (g_pred - E1_pred)) + 
        ((1 / est_P_Z_0) * (E1_pred - Y_counter_1))
    )
  
  # compute conditional expectations of phi0 and phi1
  df_obs_meas <- df_obs %>% filter(S_measured == 1)
  df_p3_meas <- df_p3 %>% filter(S_measured == 1)
  
  E_phi0_obs_mod <- glm(phi0 ~ ., data = df_obs_meas %>% select(phi0, all_of(covariates), outcome))
  E_phi0_p3_mod <- glm(phi0 ~ ., data = df_p3_meas %>% select(phi0, all_of(covariates), treatment))
  
  E_phi1_obs_mod <- glm(phi1 ~ ., data = df_obs_meas %>% select(phi1, all_of(covariates), outcome))
  E_phi1_p3_mod <- glm(phi1 ~ ., data = df_p3_meas %>% select(phi1, all_of(covariates), treatment))
  
  # predict E[phi0] and E[phi1]
  df_obs$E_phi0 <- predict(E_phi0_obs_mod, df_obs, type = "response")
  df_obs$E_phi1 <- predict(E_phi1_obs_mod, df_obs, type = "response")
  
  df_p3$E_phi0 <- predict(E_phi0_p3_mod, df_p3, type = "response")
  df_p3$E_phi1 <- predict(E_phi1_p3_mod, df_p3, type = "response")
  
  # combine datasets
  df_combined <- bind_rows(df_obs, df_p3)
  
  # compute efficient influence function adjustments
  df_combined <- df_combined %>% 
    mutate(
      psi0 = (S_measured * weights) * phi0 + (1 - (S_measured * weights)) * E_phi0, 
      psi1 = (S_measured * weights) * phi1 + (1 - (S_measured * weights)) * E_phi1
    )
  
  # compute final estimates
  Y_0_os_estimate <- Y_counter_0 + mean(df_combined$psi0) - uc_bias
  Y_1_os_estimate <- Y_counter_1 + mean(df_combined$psi1) + ct_bias - uc_bias
  
  Y_0_os_se <- sqrt(var(df_combined$psi0)) / sqrt(nrow(df_combined))
  Y_1_os_se <- sqrt(var(df_combined$psi1)) / sqrt(nrow(df_combined))
  
  log_1_minus_VE_os_estimate <- log(Y_1_os_estimate / Y_0_os_estimate)
  
  cov_mat <- matrix(cov(df_combined[,c("psi0", "psi1")]) / nrow(df_combined), nrow = 2)
  
  log_1_minus_VE_os_se <- deltamethod(
    ~ log(x2 / x1), 
    mean = c(Y_0_os_estimate, Y_1_os_estimate), 
    cov = cov_mat
  ) 
  
  z_mult <- qnorm(1 - alpha/2)
  
  # Store results in a data frame
  results_df <- data.frame(
    ct_bias = ct_bias,
    uc_bias = uc_bias,
    Y_1 = Y_1_os_estimate,
    Y_0 = Y_0_os_estimate,
    VE = 1 - (Y_1_os_estimate / Y_0_os_estimate),
    log_1_minus_VE = log_1_minus_VE_os_estimate,
    Y_1_se = Y_1_os_se,
    Y_0_se = Y_0_os_se,
    log_1_minus_VE_se = log_1_minus_VE_os_se,
    
    # Wald 95% CI for Y_1
    Y_1_lower = Y_1_os_estimate - z_mult * Y_1_os_se,
    Y_1_upper = Y_1_os_estimate + z_mult * Y_1_os_se,
    
    # Wald 95% CI for Y_0
    Y_0_lower = Y_0_os_estimate - z_mult * Y_0_os_se,
    Y_0_upper = Y_0_os_estimate + z_mult * Y_0_os_se,
    
    # Transform log-scale confidence intervals back to VE scale
    VE_lower = 1 - exp(log_1_minus_VE_os_estimate + z_mult * log_1_minus_VE_os_se),
    VE_upper = 1 - exp(log_1_minus_VE_os_estimate - z_mult * log_1_minus_VE_os_se)
  )
  
  list(
    results = results_df, 
    cov_mat = cov_mat
  )
}


# helper function for tipping point analysis
solve_for_shift <- function(shift, Y_0, Y_1, cov_mat, lower_ve_thresh, 
                            ct = T, alpha = 0.05) {
  
  z_mult <- qnorm(1 - alpha/2)  
  fixed_val <- log(1 - lower_ve_thresh)
  
  if (ct){
    
    log_1_minus_VE_os_se <- deltamethod(
      ~ log(x2 / x1), 
      mean = c(Y_0, Y_1 + shift), 
      cov = cov_mat 
    )
    
    lhs <- log((Y_1 + shift) / Y_0) + (z_mult * log_1_minus_VE_os_se)
    
    lhs - fixed_val
  }
  else{
    log_1_minus_VE_os_se <- deltamethod(
      ~ log(x2 / x1), 
      mean = c(Y_0 - shift, Y_1 - shift), 
      cov = cov_mat 
    )
    
    lhs <- log((Y_1 - shift) / (Y_0 - shift)) + (z_mult * log_1_minus_VE_os_se)
    
    lhs - fixed_val
  }
}


# function to conduct tipping point analysis 
find_tipping_point_val <- function(Y_0, Y_1, cov_mat, tipping_point = c(NA, NA), 
                                   alpha = 0.05){ 
  
  lower_ve_thresh <- tipping_point[1]
  point_ve_thresh <- tipping_point[2]
  
  if (!is.na(lower_ve_thresh)){
    result_ct <- uniroot(solve_for_shift, 
                         interval = c(0, 1), 
                         Y_0 = Y_0,
                         Y_1 = Y_1,
                         cov_mat = cov_mat,
                         lower_ve_thresh = lower_ve_thresh, 
                         ct = T, 
                         alpha = alpha)
    
    root_ct_lower <- result_ct$root
    
    result_uc <- uniroot(solve_for_shift, 
                         interval = c(-1, 0), 
                         Y_0 = Y_0,
                         Y_1 = Y_1,
                         cov_mat = cov_mat,
                         lower_ve_thresh = lower_ve_thresh, 
                         ct = F, 
                         alpha = alpha)
    
    root_uc_lower <- result_uc$root  
  }
  else{
    root_ct_lower <- Inf
    root_uc_lower <- -Inf 
  }
  if (!is.na(point_ve_thresh)){
    root_ct_point <- ((1 - point_ve_thresh) * Y_0) - Y_1
    root_uc_point <- (Y_1 + (point_ve_thresh * Y_0) - Y_0) / point_ve_thresh
  }
  else{
    root_ct_point <- Inf
    root_uc_point <- -Inf 
  }
  
  root_ct <- min(root_ct_lower, root_ct_point)
  root_uc <- max(root_uc_lower, root_uc_point)
  
  list(ct_thresh = root_ct, 
       uc_thresh = root_uc)
}

# main function to run the analysis for a single dataset 
estimate_ve_surrogate <- function(df_obs, 
                                  df_p3, 
                                  treatment, 
                                  outcome,
                                  covariates, 
                                  surrogate, 
                                  estimate_weights = T, 
                                  weight = NULL, 
                                  estimator = "plugin", 
                                  learner = "glm", 
                                  ct_bias_values = c(0), 
                                  uc_bias_values = c(0),
                                  alpha = 0.05, 
                                  success_criteria = c(NA, NA),
                                  bootstrap_reps = 100){
  
  if (!all(df_obs[[treatment]] %in% c(0, 1))) {
    stop("Treatment variable must be binary (0/1).")
  }
  
  if (!all(df_obs[[outcome]] %in% c(0, 1))) {
    stop("Outcome variable must be binary (0/1).")
  }
  
  if (!estimator %in% c("plugin", "onestep")){
    stop("The estimator must be specified as 'plugin' or 'onestep'.")
  }
  
  if (!0 %in% ct_bias_values) {
    ct_bias_values <- c(0, ct_bias_values)
  }
  
  if (!0 %in% uc_bias_values) {
    uc_bias_values <- c(0, uc_bias_values)
  }
  
  df_obs$S_measured <- ifelse(is.na(df_obs[[surrogate]]), 0, 1)
  df_p3$S_measured <- ifelse(is.na(df_p3[[surrogate]]), 0, 1)
  
  covariates_formula <- paste(covariates, collapse = " + ")
  
  if (estimate_weights) {
    
    # formula for observational dataset 
    sampling_formula_obs <- as.formula(paste("S_measured ~", 
                                             paste(c(covariates, outcome), 
                                                   collapse = " + ")))
    
    # formula for P3 dataset 
    sampling_formula_p3 <- as.formula(paste("S_measured ~", 
                                            paste(covariates, 
                                                  collapse = " + ")))
    
    # fit logistic regression model to estimate sampling probabilities for df_obs
    sampling_mod_obs <- glm(sampling_formula_obs, 
                            family = "binomial", 
                            data = df_obs)
    
    # predict probabilities and compute weights for df_obs
    df_obs$fitted_probs <- predict(sampling_mod_obs, df_obs, type = "response")
    df_obs$weights <- 1 / df_obs$fitted_probs
    
    # fit logistic regression model to estimate sampling probabilities for df_p3
    sampling_mod_p3 <- glm(sampling_formula_p3, family = "binomial", data = df_p3)
    
    # predict probabilities and compute weights for df_p3
    df_p3$fitted_probs <- predict(sampling_mod_p3, df_p3, type = "response")
    df_p3$weights <- 1 / df_p3$fitted_probs
    
  } else {
    
    if (missing(weight)) stop("If estimate_weights = FALSE, you must provide a weight variable.")
    
    # assign user-provided weights
    df_p3$weights <- df_p3[[weight]]
    df_obs$weights <- df_obs[[weight]]
  }
  
  bias_combinations <- expand.grid(ct_bias = ct_bias_values, uc_bias = uc_bias_values)
  
  results_df <- data.frame()
  
  for (i in 1:nrow(bias_combinations)) {
    
    ct_bias <- bias_combinations$ct_bias[i]
    uc_bias <- bias_combinations$uc_bias[i]
    
    if (estimator == "plugin") {
      if (learner != "glm") {
        warning("Plugin estimator requires learner = 'glm'. Overriding to 'glm'.")
        learner <- "glm"
      }
      
      plug_in_results <- compute_plug_in_estimates(df_obs, df_p3, treatment, outcome, 
                                                   covariates, surrogate, learner, 
                                                   ct_bias, uc_bias)
      
      bootstrap_results <- bootstrap_se(df_obs, df_p3, treatment, outcome, 
                                        covariates, surrogate, learner, 
                                        ct_bias, uc_bias, bootstrap_reps)
      
      z_mult <- qnorm(1 - alpha/2)  
      
      temp_result_df <- data.frame(
        ct_bias = ct_bias,
        uc_bias = uc_bias,
        Y_1 = plug_in_results$Y_1,
        Y_0 = plug_in_results$Y_0,
        VE = plug_in_results$VE,
        log_1_minus_VE = plug_in_results$log_1_minus_VE,
        Y_1_se = bootstrap_results$Y_1_se,
        Y_0_se = bootstrap_results$Y_0_se,
        log_1_minus_VE_se = bootstrap_results$log_1_minus_VE_se, 
        
        Y_1_lower = plug_in_results$Y_1 - z_mult * bootstrap_results$Y_1_se,
        Y_1_upper = plug_in_results$Y_1 + z_mult * bootstrap_results$Y_1_se,
        
        Y_0_lower = plug_in_results$Y_0 - z_mult * bootstrap_results$Y_0_se,
        Y_0_upper = plug_in_results$Y_0 + z_mult * bootstrap_results$Y_0_se,
        
        VE_lower = 1 - exp(plug_in_results$log_1_minus_VE + z_mult * bootstrap_results$log_1_minus_VE_se),
        VE_upper = 1 - exp(plug_in_results$log_1_minus_VE - z_mult * bootstrap_results$log_1_minus_VE_se)
      )
      
      cov_mat <- bootstrap_results[["cov_mat"]]
      results_df <- rbind(results_df, temp_result_df)
    } else{ 
      temp_results <- compute_onestep_estimates(df_obs, df_p3, 
                                                treatment, outcome, covariates, surrogate, 
                                                learner = learner, 
                                                ct_bias = ct_bias, 
                                                uc_bias = uc_bias,
                                                alpha = alpha)
      
      cov_mat <- temp_results[["cov_mat"]]
      temp_result_df <- temp_results[["results"]]
      results_df <- rbind(results_df, temp_result_df)
    }
    
    if (ct_bias == 0 & uc_bias == 0 ) {
      
      if (is.na(success_criteria[1]) & is.na(success_criteria[2])){
        cat("Tipping point analysis not performed - no success criteria specified.")
      }
      else if (!is.na(success_criteria[1]) & temp_result_df$VE_lower < success_criteria[1]){
        cat("Tipping point analysis not performed - success criteria not met with bias functions specified as 0.")
      }
      else if (!is.na(success_criteria[2]) & temp_result_df$VE < success_criteria[2]){
        cat("Tipping point analysis not performed - success criteria not met with bias functions specified as 0.")
      }
      else{
        cat("\n-----------------------------\n")
        cat(" Tipping Point Analysis\n")
        if (!is.na(success_criteria[1])){
          cat(sprintf("Lower VE CI = %s Threshold\n", success_criteria[1]))
        }
        if (!is.na(success_criteria[2])){
          cat(sprintf("Point Estimate = %s Threshold\n", success_criteria[2]))
        }
        
        tipping_points <- find_tipping_point_val(Y_0 = temp_result_df$Y_0, 
                                                 Y_1 = temp_result_df$Y_1, 
                                                 cov_mat = cov_mat, 
                                                 tipping_point = success_criteria, 
                                                 alpha = alpha)
        
        cat("-----------------------------\n")
        cat(" CT Bias: ", tipping_points$ct_thresh, "\n")
        cat(" UC Bias: ", tipping_points$uc_thresh, "\n")
      }
    }
  }
  
  # Return the dataframe with results
  return(results_df)
}

