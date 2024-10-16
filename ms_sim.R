library(plyr)
library(dplyr)
library(ggplot2)
library(gtools)
library(rootSolve)
library(geex)
library(msm)
library(haldensify)
library(gam)
library(nnls)
library(ROCR)
library(cvAUC)
library(SuperLearner)
library(parallel)


###############################################################################
# Data generating functions 
###############################################################################

#' Generate Y probabilities from X1, X2, X3, S, using logistic function
#'
#' P(Y = 1|X1, X2, X3, S) = scale * logistic(beta0 + beta1*S + beta2*S^2 + beta3*S^3 +
#' beta4*X1 + beta5*X2 + beta6*X3)
#'
#' @param X1 Vector of X1's
#' @param X2 Vector of X2's
#' @param X3 Vector of X3's
#' @param S Vector of S's
#' @param betas Vector of 7 beta values 
#' @param scale Scaling factor for the Y probabilities  
#'
#' @return Vector of Y probabilities
gen_Y_probs <- function(X1, X2, X3, S, 
                        betas = c(-5, -6, 6, 1, 0.69, -0.03, 0), 
                        scale = 1){
  
  beta_0 <- betas[1]
  beta_1 <- betas[2]
  beta_2 <- betas[3]
  beta_3 <- betas[4]
  beta_4 <- betas[5]
  beta_5 <- betas[6]
  beta_6 <- betas[7]
  
  scale * plogis((S * beta_1) + ((S^2) * beta_2) + ((S^3) * beta_3) + 
                   (X1 * beta_4) + (X2 * beta_5) + (X3 * beta_6) + beta_0)
}



#' Generate observational Z = 1 dataset 
#'
#' @param seed Random seed 
#' @param n_obs Number of observations in the observational (Z = 1) dataset 
#' @param S_A0_mean_var Vector of 2 values to represent the mean and SD of S for 
#' Z=1, A=0
#' @param g0_betas Vector of 7 beta values used to generate P(Y=1) probabilities
#' @param scale Scaling factor for the generated Y probabilities  
#'
#' @return Dataframe with simulated observational (Z = 1) dataset 
ms_sim_data_epi <- function(seed = 100, 
                            n_obs = 39000, 
                            S_A0_mean_var = c(-1.45, 0.15),
                            g0_betas = c(-17.1, -8.2, 0, 0, 0.69, -0.03, 0), 
                            scale = 1){
  
  if (!is.na(seed)){
    set.seed(seed)
  }
  
  Z <- 1
  A <- 0
  
  X1 <- rbinom(n = n_obs, size = 1, prob = 0.05)
  X2 <- runif(n = n_obs, min = 18, max = 40)
  X3 <- rnorm(n_obs)
  
  S <- rnorm(n_obs, mean =  S_A0_mean_var[1], sd = S_A0_mean_var[2])
  
  probs <- gen_Y_probs(X1, X2, X3, S, g0_betas, scale)
  
  Y <- rbinom(n = n_obs, size = 1, prob = probs)
  
  S_measured <- ifelse(Y == 1, 1, NA)
  S_measured[is.na(S_measured)] <- permute(c(rep(1, sum(Y == 1) * 5),
                                             rep(0, sum(is.na(S_measured)) - sum(Y == 1) * 5)))
  
  data.frame(
    S = S, 
    S_measured = S_measured, 
    X1 = X1, 
    X2 = X2, 
    X3 = X3, 
    Y = Y, 
    A = A, 
    Z = Z, 
    Y_prob = probs
  )
}


#' Generate phase 3 Z = 0 dataset 
#'
#' @param seed Random seed 
#' @param n_obs Number of observations in the Phase 3 (Z = 0) dataset
#' @param num_sampled Number sampled for S in each arm 
#' @param S_A0_mean_var Vector of 2 values to represent the mean and SD of S for
#' Z=0, A=0
#' @param S_A1_mean_var Vector of 2 values to represent the mean and SD of S for 
#' Z=0, A=1
#' @param g0_A0_betas Vector of 7 beta values used to generate P(Y=1|A=0) 
#' probabilities
#' @param g0_A1_shift True u_CT bias in data generation 
#' @param scale Scaling factor for the generated Y probabilities
#'
#' @return Dataframe with simulated Phase 3 (Z = 0) dataset
ms_sim_data_phase3 <- function(seed = 100, 
                               n_obs = 6200, 
                               num_sampled = 250,
                               S_A0_mean_var = c(-1.45, 0.15),
                               S_A1_mean_var =  c(-1.296, 0.2), 
                               g0_A0_betas = c(-17.1, -8.2, 0, 0, 0.69, -0.03, 0), 
                               g0_A1_shift = 0, 
                               scale = 1){
  if (!is.na(seed)){
    set.seed(seed)
  }
  
  Z <- 0
  A <- rep(c(0,1), each = n_obs/2)
  
  X1 <- rbinom(n = n_obs, size = 1, prob = 0.05)
  X2 <- runif(n = n_obs, min = 18, max = 40)
  X3 <- rnorm(n_obs)
  
  S_0 <- rnorm(n_obs / 2, mean = S_A0_mean_var[1], sd = S_A0_mean_var[2])
  S_1 <- rnorm(n_obs / 2, mean = S_A1_mean_var[1], sd = S_A1_mean_var[2])
  
  probs_A0 <- gen_Y_probs(X1[1:(n_obs/2)], 
                          X2[1:(n_obs/2)], 
                          X3[1:(n_obs/2)], 
                          S_0, 
                          g0_A0_betas, 
                          scale)
  
  probs_A1 <- gen_Y_probs(X1[((n_obs/2)+1):n_obs], 
                          X2[((n_obs/2)+1):n_obs], 
                          X3[((n_obs/2)+1):n_obs], 
                          S_1, 
                          g0_A0_betas, 
                          scale) + g0_A1_shift
  
  S <- c(S_0, S_1)
  probs <- c(probs_A0, probs_A1)
  
  S_measured <- c(rep(1, num_sampled), rep(0, (n_obs/2) - num_sampled), 
                  rep(1, num_sampled), rep(0, (n_obs/2) - num_sampled))
  
  Y <- rbinom(n = n_obs, size = 1, prob = probs)
  
  data.frame(
    S = S, 
    X1 = X1, 
    X2 = X2, 
    X3 = X3, 
    Y = Y, 
    A = A, 
    Z = Z, 
    S_measured = S_measured
  )
}


###############################################################################
# Nuisance function estimation 
###############################################################################

#' Estimate sampling probabilities 
#'
#' @param df Dataframe 
#' @param formula_str Formula to use in glm
#'
#' @return Logistic regression glm model object for sampling probabilities
get_sampling_mod <- function(df, formula_str = "S_measured ~ X1 + X2 + X3 + Y"){
  logit_mod <- glm(as.formula(formula_str), 
                   family = "binomial",
                   data = df)
}


#' Estimate g0 function from observational data 
#'
#' @param df Dataframe 
#' @param flexible T if we want to fit g0 using SuperLearner, F if we want
#' to use logistic regression
#' @param fit_sampling_probs T if we want to fit sampling probabilities, F if 
#' we want to use the true sampling probabilities
#'
#' @return Linear regression glm model object for g0 function
estimate_g0_from_obs <- function(df, flexible = F, fit_sampling_probs = F){
  
  if (fit_sampling_probs){
    sampling_mod <- get_sampling_mod(df)
    df$fitted_probs <- predict(sampling_mod, df, type = "response")
    df$weights <- 1 / df$fitted_probs 
  }
  else{
    df <- df %>% 
      mutate(fitted_probs = ifelse(Y == 1, 1, (5*sum(df$Y == 1)) / sum(df$Y == 0)), 
             weights = 1 / fitted_probs)
  }
  
  df_nonmissing <- df %>% filter(S_measured == 1)
  
  if (flexible){
    g0_mod <- SuperLearner(Y = df_nonmissing$Y, 
                           X = df_nonmissing %>% select(X1, X2, X3, S), 
                           family = binomial(),
                           SL.library = c("SL.mean", "SL.glm", "SL.gam"), 
                           obsWeights = df_nonmissing$weights)
  }
  else{
    g0_mod <- glm(Y ~ X1 + X2 + X3 + S,
                  data = df_nonmissing,
                  family = "binomial",
                  weights = df_nonmissing$weights)
  }
  g0_mod
}

#' Function to estimate the outer expectation
#'
#' @param g0_mod glm model object 
#' @param df_p3 Phase 3 dataframe 
#' @param formula_str Formula to use in glm  
#' @param flexible T if we want to fit outer expectation using SuperLearner, F 
#' if we want to use linear regression
#'
#' @return A list of 2 glm model objects, (1) for E_hat_A0 and (2) for E_hat_A1
estimate_E_hat <- function(g0_mod, df_p3, formula_str = "Y_pred ~ X1 + X2 + X3", 
                           flexible = F){
  
  if (flexible){
    df_p3$Y_pred <- predict(g0_mod, df_p3, type = "response")$pred[,1]
  }
  else{
    df_p3$Y_pred <- predict(g0_mod, df_p3, type = "response")
  }
  
  df_treated <- df_p3 %>% filter(A == 1 & S_measured == 1)
  df_placebo <- df_p3 %>% filter(A == 0 & S_measured == 1)
  
  if (flexible){
    E_hat_A0 <- SuperLearner(Y = df_placebo$Y_pred, 
                             X = df_placebo %>% select(X1, X2, X3), 
                             family = gaussian(),
                             SL.library = c("SL.mean", "SL.glm", "SL.gam"), 
                             obsWeights = df_placebo$weights)
    
    E_hat_A1 <- SuperLearner(Y = df_treated$Y_pred, 
                             X = df_treated %>% select(X1, X2, X3), 
                             family = gaussian(),
                             SL.library = c("SL.mean", "SL.glm", "SL.gam"), 
                             obsWeights = df_treated$weights)
  }
  else{
    E_hat_A0 <- glm(as.formula(formula_str),
                    data = df_placebo, 
                    weights = df_placebo$weights)
    
    E_hat_A1 <- glm(as.formula(formula_str),
                    data = df_treated, 
                    weights = df_treated$weights)
  }
  
  list(
    E_hat_A0 = E_hat_A0, 
    E_hat_A1 = E_hat_A1
  )
}


#' Estimate P(Z=0)
#'
#' @param df_obs Dataframe for observational study
#' @param df_p3 Dataframe for phase 3 study
#'
#' @return Estimate of P(Z=0)
estimate_P_Z <- function(df_obs, df_p3){
  nrow(df_p3) / (nrow(df_obs) + nrow(df_p3))
}


#' Estimate P(A=a | X, Z = 0)
#'
#' @param df_obs Dataframe for observational study
#' @param df_p3 Dataframe for phase 3 study
#'
#' @return Assumed probability of P(A = 1 | X, Z = 0) = 0.5
estimate_P_A <- function(df_obs, df_p3){
  0.5
}


#' Estimate P(S_measured = 1| Z, A, X)
#'
#' @param df_obs Dataframe for observational study
#' @param df_p3 Dataframe for phase 3 study
#'
#' @return Estimate for P(S_measured = 1 | Z, A, X)
estimate_P_S_measured_ZAX <- function(df_obs, df_p3){
  
  log_model_Z0 <- glm(S_measured ~ X1 + X2 + X3 + A, 
                      data = df_p3, 
                      family = binomial)
  
  log_model_Z1 <- glm(S_measured ~ X1 + X2 + X3, 
                      data = df_obs, 
                      family = binomial)
  
  list(
    log_model_Z0 = log_model_Z0, 
    log_model_Z1 = log_model_Z1
  )
  
}


#' Estimate P(Z = z, A = a | X, S)
#'
#' @param df_obs Dataframe for observational study
#' @param df_p3 Dataframe for phase 3 study
#' @param flexible T if we want to fit using SuperLearner, F if we want
#' to use logistic regression
#'
#' @return Estimates for P(Z = z, A = a | X, S)
estimate_P_Z_A_XS <- function(df_obs, df_p3, flexible = F){
  
  sampling_mod <- estimate_P_S_measured_ZAX(df_obs, df_p3)
  
  log_model_Z0 <- sampling_mod[["log_model_Z0"]]
  log_model_Z1 <- sampling_mod[["log_model_Z1"]]
  
  df_obs$fitted_probs <- predict(log_model_Z1, df_obs, type = "response")
  df_obs$weights <- 1 / df_obs$fitted_probs 
  
  df_p3$fitted_probs <- predict(log_model_Z0, df_p3, type = "response")
  df_p3$weights <- 1 / df_p3$fitted_probs 
  
  df_combined <- df_obs %>% rbind.fill(df_p3)
  
  df_nonmissing <- df_combined %>% 
    filter(S_measured == 1) %>%
    mutate(Z0A1 = ifelse(Z == 0 & A == 1, 1, 0), 
           Z0A0 = ifelse(Z == 0 & A == 0, 1, 0))
  
  
  if (flexible){
    log_model_Z0A1 <- SuperLearner(Y = df_nonmissing$Z0A1, 
                                   X = df_nonmissing %>% select(X1, X2, X3, S), 
                                   family = binomial(),
                                   SL.library = c("SL.mean", "SL.glm", "SL.gam"), 
                                   obsWeights = df_nonmissing$weights)
    
    log_model_Z0A0 <- SuperLearner(Y = df_nonmissing$Z0A0, 
                                   X = df_nonmissing %>% select(X1, X2, X3, S), 
                                   family = binomial(),
                                   SL.library = c("SL.mean", "SL.glm", "SL.gam"), 
                                   obsWeights = df_nonmissing$weights)
    
    log_model_Z1 <- SuperLearner(Y = df_nonmissing$Z, 
                                 X = df_nonmissing %>% select(X1, X2, X3, S), 
                                 family = binomial(),
                                 SL.library = c("SL.mean", "SL.glm", "SL.gam"), 
                                 obsWeights = df_nonmissing$weights)
    
  }
  else{
    log_model_Z0A1 <- glm(Z0A1 ~ S + X1 + X2 + X3, 
                          data = df_nonmissing, 
                          family = "binomial",
                          weights = df_nonmissing$weights)
    
    log_model_Z0A0 <- glm(Z0A0 ~ S + X1 + X2 + X3, 
                          data = df_nonmissing, 
                          family = "binomial",
                          weights = df_nonmissing$weights)
    
    log_model_Z1 <- glm(Z ~ S + X1 + X2 + X3, 
                        data = df_nonmissing, 
                        family = "binomial",
                        weights = df_nonmissing$weights)
  }
  
  list(
    log_model_Z0A1 = log_model_Z0A1, 
    log_model_Z0A0 = log_model_Z0A0,
    log_model_Z1 = log_model_Z1
  )
}


###############################################################################
# NEW ESTIMATION/INFERENCE FUNCTIONS 
###############################################################################

#' Plug-in estimator 
#'
#' @param df_obs Dataframe for observational study
#' @param df_p3 Dataframe for phase 3 study
#' @param ct_bias CT bias used for estimation 
#' @param uc_bias UC bias used for estimation 
#' @param fit_sampling_probs T if you want to fit the sampling probabilities, F 
#' if you want to use known sampling probabilities
#'
#' @return List of estimates for E(Y(1)|Z=0), E(Y(0)|Z=0), and VE
estimate_plug_in <- function(df_obs, df_p3, ct_bias = 0, uc_bias = 0, 
                             fit_sampling_probs = F){ 
  
  g0_mod_est <- estimate_g0_from_obs(df_obs, fit_sampling_probs)
  
  if (fit_sampling_probs){
    # estimate sampling probability from Phase 3 dataset 
    pi_mod_est <- get_sampling_mod(df_p3, 
                                   formula_str = "S_measured ~ X1 + X2 + X3 + Y")
    
    # add weight columns (inverse of sampling probability)
    df_p3$fitted_probs <- predict(pi_mod_est, 
                                  newdata = df_p3, 
                                  type = "response")
    
    df_p3$weights <- 1 / df_p3$fitted_probs 
  }
  else{
    num_sampled <- nrow(df_p3 %>% filter(S_measured == 1))
    df_p3$weights <- nrow(df_p3) / num_sampled
  }
  
  E_hats <- estimate_E_hat(g0_mod_est, df_p3)
  
  df_p3$E0_pred <- predict(E_hats[["E_hat_A0"]], df_p3)
  df_p3$E1_pred <- predict(E_hats[["E_hat_A1"]], df_p3)
  
  # estimate counterfactual E[Y(1)]
  Y_counter_1 <- mean(df_p3$E1_pred) + ct_bias
  
  # estimate counterfactual E[Y(0)]
  Y_counter_0 <- mean(df_p3$E0_pred) - uc_bias
  
  # estimate vaccine efficacy 
  VE_est <- 1 - (Y_counter_1 / Y_counter_0)
  
  list(
    Y_1 = Y_counter_1, 
    Y_0 = Y_counter_0, 
    VE = VE_est, 
    log_1_minus_VE = log(1 - VE_est)
  )
}


#' Bootstrap estimator
#'
#' @param df_obs Dataframe for observational study
#' @param df_p3 Dataframe for phase 3 study
#' @param boots Number of bootstrap runs 
#' @param ct_bias CT bias used for estimation 
#'
#' @return
bootstrap_estimator <- function(df_obs, df_p3, boots = 500, ct_bias = 0, old = F){
  
  if (!old){
    estimate_fun <- estimate_plug_in
  }
  else{
    estimate_fun <- estimate_plug_in_old
  }
  
  orig_res <- estimate_fun(df_obs, df_p3, ct_bias)
  
  Y_1s <- c() 
  Y_0s <- c() 
  VEs <- c()
  log_1_minus_VEs <- c() 
  
  for (i in 1:boots){
    # for observational study, sample separately cases and controls
    df_cases <- df_obs %>% filter(Y == 1)
    df_controls <- df_obs %>% filter(Y == 0)
    
    idx_obs_cases <- sample(1:nrow(df_cases), replace = T)
    idx_obs_controls <- sample(1:nrow(df_controls), replace = T)
    
    boot_df_obs <- rbind(df_cases[idx_obs_cases,], 
                         df_controls[idx_obs_controls,])
    
    # sample phase 3 study
    idx_p3 <- sample(1:nrow(df_p3), replace = T)
    boot_df_p3 <- df_p3[idx_p3,]
    
    boot_res <- estimate_fun(boot_df_obs, 
                                 boot_df_p3, 
                                 ct_bias)
    
    Y_1s <- c(Y_1s, boot_res$Y_1)
    Y_0s <- c(Y_0s, boot_res$Y_0)
    VEs <- c(VEs, boot_res$VE)
    log_1_minus_VEs <- c(log_1_minus_VEs, boot_res$log_1_minus_VE)
  }
  
  list(
    Y_0_bs_estimate = orig_res$Y_0,
    Y_1_bs_estimate = orig_res$Y_1,
    VE_bs_estimate = orig_res$VE,
    log_1_minus_VE_bs_estimate = orig_res$log_1_minus_VE,
    
    Y_0_bs_se = sd(Y_0s),
    Y_1_bs_se = sd(Y_1s),
    VE_bs_se = sd(VEs),
    log_1_minus_VE_bs_se = sd(log_1_minus_VEs), 
    
    Y_0_emp_bs_lower_ci = as.numeric(quantile(Y_0s, 0.025)), 
    Y_0_emp_bs_upper_ci = as.numeric(quantile(Y_0s, 0.975)), 
    Y_1_emp_bs_lower_ci = as.numeric(quantile(Y_1s, 0.025)), 
    Y_1_emp_bs_upper_ci = as.numeric(quantile(Y_1s, 0.975)), 
    VE_emp_bs_lower_ci = as.numeric(quantile(VEs, 0.025)), 
    VE_emp_bs_upper_ci = as.numeric(quantile(VEs, 0.975))
  )
}


#' Helper function 1 for sandwich variance  
#'
#' @param X1 Vector of X1's
#' @param X2 Vector of X2's
#' @param X3 Vector of X3's
#' @param S Vector of S's
#' @param beta0 
#' @param beta1 
#' @param beta2 
#' @param beta3 
#' @param beta4 
#'
#' @return Logistic transformation of X1, X2, X3, and S with beta coefficients 
pi1 <- function(X1, X2, X3, S, beta0, beta1, beta2, beta3, beta4){
  exp(beta0 + (beta1 * X1) + (beta2 * X2) + (beta3 * X3) + (beta4 * S)) / 
    (1 + exp(beta0 + (beta1 * X1) + (beta2 * X2) + (beta3 * X3) + (beta4 * S)))
}


#' Helper function 2 for sandwich variance  
#'
#' @param X1 vector of X1's
#' @param X2 vector of X2's
#' @param X3 vector of X3's
#' @param beta0 
#' @param beta1 
#' @param beta2 
#' @param beta3 
#'
#' @return Linear transformation of X1, X2, and X3 with beta coefficients 
pi2 <- function(X1, X2, X3, beta0, beta1, beta2, beta3){
  beta0 + (beta1 * X1) + (beta2 * X2) + (beta3 * X3)
}


#' Stacked estimating functions for Stefanski and Boos M-estimating equations 
#'
#' @param data Combined observational and phase 3 dataframe 
#' @param sampling_prob_Z1_Y0 True sampling probability for Z = 1, Y = 0 group
#' @param sampling_prob_Z0 True sampling probability for Z = 0 group 
#' @param ct_bias CT bias used for estimation
#'
#' @return Estimating equation function
sb_estfun <- function(data, sampling_prob_Z1_Y0, sampling_prob_Z0, ct_bias = 0){
  X1 <- data$X1
  X2 <- data$X2
  X3 <- data$X3
  Z <- data$Z 
  A <- data$A
  S <- data$S
  Y <- data$Y
  S_measured <- data$S_measured
  
  sampling_probs <- ifelse(
    Z == 0, sampling_prob_Z0, 
    ifelse(
      Z == 1 & Y == 1, 1, 
      sampling_prob_Z1_Y0
    )
  )
  
  # thetas[1-5]: parameters for g 
  # thetas[6-9]: parameters for E hat (A = 0)
  # thetas[10-13]: parameters for E hat (A = 1)
  # thetas[14-16]: E[Y(0)], E[Y(1)] and VE
  
  function(theta){
    c(
      # estimating equations for g parameters 
      Z * S_measured * (Y - pi1(X1, X2, X3, S, 
                                theta[1], theta[2], theta[3], 
                                theta[4], theta[5])) / 
        sampling_probs, 
      
      X1 * Z * S_measured * (Y - pi1(X1, X2, X3, S, 
                                     theta[1], theta[2], theta[3], 
                                     theta[4], theta[5])) / 
        sampling_probs, 
      
      X2 * Z * S_measured * (Y - pi1(X1, X2, X3, S, 
                                     theta[1], theta[2], theta[3], 
                                     theta[4], theta[5])) / 
        sampling_probs, 
      
      X3 * Z * S_measured * (Y - pi1(X1, X2, X3, S, 
                                     theta[1], theta[2], theta[3], 
                                     theta[4], theta[5])) / 
        sampling_probs, 
      
      S * Z * S_measured * (Y - pi1(X1, X2, X3, S, 
                                    theta[1], theta[2], theta[3], 
                                    theta[4], theta[5])) / 
        sampling_probs,
      
      # estimating equations for E hat parameters (A = 0)
      (1-Z) * (1-A) * S_measured * (pi1(X1, X2, X3, S, 
                                theta[1], theta[2], theta[3], 
                                theta[4], theta[5]) - 
                                  
                                  pi2(X1, X2, X3, theta[6], theta[7], 
                                      theta[8], theta[9])) / 
        sampling_probs, 
      
      X1 * (1-Z) * (1-A) * S_measured * (pi1(X1, X2, X3, S, 
                                             theta[1], theta[2], theta[3], 
                                             theta[4], theta[5]) - 
                                           
                                           pi2(X1, X2, X3, theta[6], theta[7], 
                                               theta[8], theta[9])) / 
        sampling_probs, 
      
      X2 * (1-Z) * (1-A) * S_measured * (pi1(X1, X2, X3, S, 
                                             theta[1], theta[2], theta[3], 
                                             theta[4], theta[5]) - 
                                           
                                           pi2(X1, X2, X3, theta[6], theta[7], 
                                               theta[8], theta[9])) / 
        sampling_probs, 
      
      X3 * (1-Z) * (1-A) * S_measured * (pi1(X1, X2, X3, S, 
                                             theta[1], theta[2], theta[3], 
                                             theta[4], theta[5]) - 
                                           
                                           pi2(X1, X2, X3, theta[6], theta[7], 
                                               theta[8], theta[9])) / 
        sampling_probs, 
      
      # estimating equations for E hat parameters (A = 1)
      (1-Z) * A * S_measured * (pi1(X1, X2, X3, S, 
                                        theta[1], theta[2], theta[3], 
                                        theta[4], theta[5]) - 
                                      
                                  pi2(X1, X2, X3, theta[10], theta[11], 
                                      theta[12], theta[13])) / 
        sampling_probs, 
      
      X1 * (1-Z) * A * S_measured * (pi1(X1, X2, X3, S, 
                                             theta[1], theta[2], theta[3], 
                                             theta[4], theta[5]) - 
                                           
                                       pi2(X1, X2, X3, theta[10], theta[11], 
                                           theta[12], theta[13])) / 
        sampling_probs, 
      
      X2 * (1-Z) * A * S_measured * (pi1(X1, X2, X3, S, 
                                             theta[1], theta[2], theta[3], 
                                             theta[4], theta[5]) - 
                                           
                                       pi2(X1, X2, X3, theta[10], theta[11], 
                                           theta[12], theta[13])) / 
        sampling_probs, 
      
      X3 * (1-Z) * A * S_measured * (pi1(X1, X2, X3, S, 
                                             theta[1], theta[2], theta[3], 
                                             theta[4], theta[5]) - 
                                           
                                       pi2(X1, X2, X3, theta[10], theta[11], 
                                           theta[12], theta[13])) / 
        sampling_probs,
      
      # estimating equations for Y(0), Y(1), VE
      (1 - Z) * (pi2(X1, X2, X3, theta[6], theta[7], theta[8], theta[9]) - 
                   theta[14]),

      (1 - Z) * (pi2(X1, X2, X3, theta[10], theta[11], theta[12], theta[13]) + 
                   ct_bias - theta[15])
    )
  }
}


#' Helper function to get true sampling probabilities 
#'
#' @param df_obs Dataframe for observational study
#' @param df_p3 Dataframe for phase 3 study
#'
#' @return List of true sampling probabilites for different groups 
get_true_sampling_probs <- function(df_obs, df_p3){
  
  sampling_prob_Z1_Y0 <- (sum(df_obs$Y==1) * 5) /
    sum(df_obs$Y == 0)

  sampling_prob_Z0 <- nrow(df_p3 %>% filter(S_measured == 1)) / nrow(df_p3) 
  
  sampling_prob_Z1_Y1 <- 1
  
  list(
    sampling_prob_Z1_Y0 = sampling_prob_Z1_Y0, 
    sampling_prob_Z0 = sampling_prob_Z0, 
    sampling_prob_Z1_Y1 = sampling_prob_Z1_Y1
  )
  
}


#' Function to run sandwich estimator 
#'
#' @param df_obs Data for observational study
#' @param df_p3 Data for phase 3 study
#' @param ct_bias CT bias used for estimation 
#'
#' @return List of estimates and standard errors for E[Y(0)|Z=0], E[Y(1)|Z=0], 
#' and VE
sandwich_estimator <- function(df_obs, df_p3, ct_bias = 0){
  
  df_combined <- rbind.fill(
    df_obs %>% mutate(Z = 1) %>% filter(S_measured == 1), 
    df_p3 %>% mutate(Z = 0) %>% filter(S_measured == 1)
  ) 
  
  sampling_probs <- get_true_sampling_probs(df_obs, df_p3)
  
  sampling_prob_Z1_Y0 <- sampling_probs[["sampling_prob_Z1_Y0"]]
  
  sampling_prob_Z0 <- sampling_probs[["sampling_prob_Z0"]]
  
  results <- m_estimate(
    estFUN = sb_estfun, 
    data   = df_combined, 
    outer_args = list(ct_bias = ct_bias, 
                      sampling_prob_Z1_Y0 = sampling_prob_Z1_Y0, 
                      sampling_prob_Z0 = sampling_prob_Z0),
    root_control = setup_root_control(start = c(rep(0,13), 0.1, 0.1)), 
    deriv_control = setup_deriv_control(method = "simple")
  )
  
  log_1_minus_VE_estimate = log(coef(results)[15] / coef(results)[14])
  
  log_1_minus_VE_se <- deltamethod(
    ~ log(x2 / x1), mean = coef(results)[14:15], cov = vcov(results)[14:15, 14:15]
  ) 
  
  list(
    Y_0_sw_estimate = coef(results)[14],
    Y_1_sw_estimate = coef(results)[15],
    log_1_minus_VE_sw_estimate = log_1_minus_VE_estimate,
    
    Y_0_sw_se = sqrt(vcov(results)[14,14]),
    Y_1_sw_se = sqrt(vcov(results)[15,15]),
    log_1_minus_VE_sw_se = log_1_minus_VE_se, 
    
    VE_sw_estimate = 1 - exp(log_1_minus_VE_estimate), 
    VE_sw_lower_ci = 1 - exp(log_1_minus_VE_estimate + (1.96 * log_1_minus_VE_se)), 
    VE_sw_upper_ci = 1 - exp(log_1_minus_VE_estimate - (1.96 * log_1_minus_VE_se))
  )
}

###############################################################################
# ONE-STEP ESTIMATOR  
###############################################################################

#' Function to run one-step estimator
#'
#' @param df_obs Dataframe for observational study
#' @param df_p3 Dataframe for phase 3 study
#' @param flexible T if we want to fit g0 using SuperLearner, F if we want
#' to use logistic regression
#' @param ct_bias CT bias used for estimation 
#' @param uc_bias UC bias used for estimation 
#'
#' @return One-step estimator results
one_step_estimator <- function(df_obs, df_p3, flexible = T, ct_bias = 0, 
                               uc_bias = 0){
  
  # get true sampling probabilities 
  sampling_probs <- get_true_sampling_probs(df_obs, df_p3)
  
  sampling_prob_Z1_Y0 <- sampling_probs[["sampling_prob_Z1_Y0"]]
  sampling_prob_Z0 <- sampling_probs[["sampling_prob_Z0"]]
  sampling_prob_Z1_Y1 <- sampling_probs[["sampling_prob_Z1_Y1"]]
  
  df_obs <- df_obs %>% mutate(sampling_prob = ifelse(Y == 1, sampling_prob_Z1_Y1, 
                                                     sampling_prob_Z1_Y0))
  
  df_p3 <- df_p3 %>% mutate(sampling_prob = sampling_prob_Z0)
  
  df_p3$weights <- 1 / df_p3$sampling_prob
  
  # fit g0
  g0_mod_est <- estimate_g0_from_obs(df_obs, flexible = flexible)
  
  # fit outer expectation E
  E_hats <- estimate_E_hat(g0_mod_est, df_p3, flexible = flexible)
  
  # get E0 and E1 predictions 
  if (flexible){
    df_p3$g_pred <- predict(g0_mod_est, df_p3, type = "response")$pred[,1]
    df_obs$g_pred <- predict(g0_mod_est, df_obs, type = "response")$pred[,1]
    
    df_p3$E0_pred <- predict(E_hats[["E_hat_A0"]], df_p3, type = "response")$pred[,1]
    df_p3$E1_pred <- predict(E_hats[["E_hat_A1"]], df_p3,  type = "response")$pred[,1]
  }
  else{
    df_p3$g_pred <- predict(g0_mod_est, df_p3, type = "response")
    df_obs$g_pred <- predict(g0_mod_est, df_obs, type = "response")
    
    df_p3$E0_pred <- predict(E_hats[["E_hat_A0"]], df_p3, type = "response")
    df_p3$E1_pred <- predict(E_hats[["E_hat_A1"]], df_p3,  type = "response")
  }
  
  # estimate counterfactual E[Y(1)]: plug-in estimate
  Y_counter_1 <- mean(df_p3$E1_pred) 
  
  # estimate counterfactual E[Y(0)]: plug-in estimate
  Y_counter_0 <- mean(df_p3$E0_pred)
  
  # estimate P(Z = 0)
  est_P_Z_0 <- estimate_P_Z(df_obs, df_p3)
  
  # estimate P(A = a | X, Z = 0)
  est_P_A <- estimate_P_A(df_obs, df_p3)
  
  # estimate P(Z=z,A=a | X, S)
  log_mod_P_ZA_XS <- estimate_P_Z_A_XS(df_obs, df_p3, flexible = flexible)
  
  # get P_Z0A1_preds 
  if (flexible){
    df_obs$P_Z0A1_pred <- predict(log_mod_P_ZA_XS[["log_model_Z0A1"]], 
                                  df_obs,
                                  type = "response")$pred[,1]
    
    df_obs$P_Z0A0_pred <- predict(log_mod_P_ZA_XS[["log_model_Z0A0"]],
                                  df_obs, 
                                  type = "response")$pred[,1]
    
    df_obs$P_Z1A0_pred <- predict(log_mod_P_ZA_XS[["log_model_Z1"]],
                                  df_obs, 
                                  type = "response")$pred[,1]
  }
  else{
    df_obs$P_Z0A1_pred <- predict(log_mod_P_ZA_XS[["log_model_Z0A1"]], 
                                  df_obs,
                                  type = "response")
    
    df_obs$P_Z0A0_pred <- predict(log_mod_P_ZA_XS[["log_model_Z0A0"]],
                                  df_obs, 
                                  type = "response")
    
    df_obs$P_Z1A0_pred <- predict(log_mod_P_ZA_XS[["log_model_Z1"]],
                                  df_obs, 
                                  type = "response")
  }
  # get P_Z0A1_preds
  df_obs <- df_obs %>% 
    
    mutate(
      phi0 = (1 / est_P_Z_0) * (1 / est_P_A) * (P_Z0A0_pred / P_Z1A0_pred) * (Y - g_pred),
      phi1 = (1 / est_P_Z_0) * (1 / est_P_A) * (P_Z0A1_pred / P_Z1A0_pred) * (Y - g_pred) 
    )
  
  
  df_p3 <- df_p3 %>% 
    
    mutate(
      phi0 = ((1 / est_P_Z_0) * ((1-A) / (1 - est_P_A)) * (g_pred - E0_pred)) + 
        ((1 / est_P_Z_0) * (E0_pred - Y_counter_0)),
      phi1 =  ((1 / est_P_Z_0) * (A / est_P_A) * (g_pred - E1_pred)) + 
        ((1 / est_P_Z_0) * (E1_pred - Y_counter_1))
    )
  
  # calculate phi_0 and phi_1
  df_obs_meas <- df_obs %>% filter(S_measured == 1)
  df_p3_meas <- df_p3 %>% filter(S_measured == 1)
  
  E_phi0_obs_mod <- glm(phi0 ~ X1 + X2 + X3 + Y, data = df_obs_meas)
  E_phi0_p3_mod <- glm(phi0 ~ X1 + X2 + X3 + A, data = df_p3_meas)
  
  E_phi1_obs_mod <- glm(phi1 ~ X1 + X2 + X3 + Y, data = df_obs_meas)
  E_phi1_p3_mod <- glm(phi1 ~ X1 + X2 + X3 + A, data = df_p3_meas)
  
  df_obs$E_phi0 <- predict(E_phi0_obs_mod, df_obs,  type = "response")
  df_obs$E_phi1 <- predict(E_phi1_obs_mod, df_obs,  type = "response")
  
  df_p3$E_phi0 <- predict(E_phi0_p3_mod, df_p3, type = "response")
  df_p3$E_phi1 <- predict(E_phi1_p3_mod, df_p3,  type = "response")
  
  df_combined <- df_obs %>% rbind.fill(df_p3)
  
  df_combined <- df_combined %>% 
    mutate(psi0 = (S_measured / sampling_prob) * phi0 +
             (1 - (S_measured / sampling_prob)) * E_phi0, 
           psi1 = (S_measured / sampling_prob) * phi1 +
             (1 - (S_measured / sampling_prob)) * E_phi1)
  
  Y_0_os_estimate <- Y_counter_0 + mean(df_combined$psi0) - uc_bias
  Y_1_os_estimate <- Y_counter_1 + mean(df_combined$psi1) + ct_bias - uc_bias
  
  Y_0_os_se <- sqrt(var(df_combined$psi0)) / sqrt(nrow(df_combined))
  Y_1_os_se <- sqrt(var(df_combined$psi1)) / sqrt(nrow(df_combined))
  
  log_1_minus_VE_os_estimate = log(Y_1_os_estimate / Y_0_os_estimate)
  
  log_1_minus_VE_os_se <- deltamethod(
    ~ log(x2 / x1), 
    mean = c(Y_0_os_estimate, Y_1_os_estimate), 
    cov = matrix(cov(df_combined[,c("psi0", "psi1")]) / nrow(df_combined), 
                 nrow = 2)
  ) 
  
  list(
    Y_0_os_estimate = Y_0_os_estimate,
    Y_1_os_estimate = Y_1_os_estimate,
    
    Y_0_os_se = Y_0_os_se,
    Y_1_os_se = Y_1_os_se, 
    log_1_minus_VE_os_estimate = log_1_minus_VE_os_estimate, 
    log_1_minus_VE_os_se = log_1_minus_VE_os_se
  )
  
  
}



###############################################################################
# RUN SINGLE SIMULATION   
###############################################################################

#' Function to run single simulation
#'
#' @param seed Random seed
#' @param g0_betas Vector of 7 beta values used to generate P(Y=1) probabilities
#' @param S_A1_mean_var Vector of 2 values for representing the mean and SD for S 
#' in the A = 1 arm 
#' @param S_A0_mean_var Vector of 2 values for representing the mean and SD for S 
#' in the A = 1 arm 
#' @param n_obs_obs Number of observations in the observational study
#' @param n_obs_p3 Number of observations in the phase 3 study
#' @param num_sampled Number sampled for S in the phase 3 study
#' @param boots Number of bootstrap runs 
#' @param run_type "bs" for bootstrap, "sw" for sandwich, "one-step" for one-step, 
#' "all" for all 
#' @param true_ct_bias True CT bias in data generation
#' @param est_ct_bias CT bias used for estimation
#' @param scale Scale factor when generating P(Y=1) probabilities 
#'
#' @return List of estimates and standard errors for E[Y(0)|Z=0], E[Y(1)|Z=0], 
#' and VE from randomly generated datasets
run_single_sim <- function(seed = NA,
                           g0_betas = c(-17.1, -8.2, 0, 0, 0.69, -0.03, 0),
                           S_A0_mean_var = c(-1.5, 0.3), 
                           S_A1_mean_var = c(-1.45, 0.15),
                           n_obs_obs = 39000,
                           n_obs_p3 = 6200, 
                           num_sampled = 250,
                           boots = 500, 
                           run_type = "bs", 
                           true_ct_bias = 0, 
                           est_ct_bias = 0, 
                           scale = 1){
  
  df_obs <- ms_sim_data_epi(seed = seed, 
                            g0_betas = g0_betas,
                            S_A0_mean_var = S_A0_mean_var, 
                            n_obs = n_obs_obs,
                            scale = scale)
  
  df_p3 <- ms_sim_data_phase3(seed = seed,
                              S_A0_mean_var = S_A0_mean_var, 
                              S_A1_mean_var = S_A1_mean_var, 
                              n_obs = n_obs_p3,
                              num_sampled = num_sampled,
                              g0_A1_shift = true_ct_bias, 
                              g0_A0_betas = g0_betas,
                              scale = scale)
  
  run_estimator <- function(df_obs, df_p3){
    if (run_type == "bs"){
      bootstrap_estimator(df_obs, df_p3, boots = boots, ct_bias = est_ct_bias)
    }
    else if  (run_type == "sandwich"){
      sandwich_estimator(df_obs, df_p3, ct_bias = est_ct_bias)
    }
    else if (run_type == "one-step"){
      one_step_estimator(df_obs, df_p3, flexible = T, 
                         ct_bias = est_ct_bias)
    }
    else{
      bs_res <- bootstrap_estimator(df_obs, df_p3, boots = boots, 
                                    ct_bias = est_ct_bias)
      
      sandwich_res <- sandwich_estimator(df_obs, df_p3, 
                                         ct_bias = est_ct_bias)
      
      one_step_res <- one_step_estimator(df_obs, df_p3, flexible = T, 
                                         ct_bias = est_ct_bias)
      
      c(bs_res, sandwich_res, one_step_res)
    }
  }
  
  res <- tryCatch(
    {
      run_estimator(df_obs, df_p3)
    },
    error = function(cond) {
      message(conditionMessage(cond))
      message("Rerunning with different seed...")
      seed_new <- seed + 1000000
      
      print(seed_new)
      df_obs <- ms_sim_data_epi(seed = seed_new, 
                                g0_betas = g0_betas,
                                S_A0_mean_var = S_A0_mean_var, 
                                n_obs = n_obs_obs,
                                scale = scale)
      
      df_p3 <- ms_sim_data_phase3(seed = seed_new,
                                  S_A0_mean_var = S_A0_mean_var, 
                                  S_A1_mean_var = S_A1_mean_var, 
                                  n_obs = n_obs_p3,
                                  num_sampled = num_sampled,
                                  g0_A1_shift = true_ct_bias, 
                                  g0_A0_betas = g0_betas,
                                  scale = scale)
      run_estimator(df_obs, df_p3)
    }
  )
  
  res
}



###############################################################################
# SOME SAMPLE CODE TO TRY
###############################################################################

# create simulated datasets 
df_obs <- ms_sim_data_epi(seed = 2)
df_p3 <- ms_sim_data_phase3(seed = 2)

# estimate plug-in 
estimate_plug_in(df_obs, df_p3)

# bootstrap estimator
bootstrap_estimator(df_obs, df_p3, boots = 100)

# sandwich estimator
sandwich_estimator(df_obs, df_p3)

# one-step estimator 
one_step_estimator(df_obs, df_p3)

# runs a full simulation for true VE = 0 case
run_single_sim(seed = 1, 
               S_A1_mean_var = c(-1.45, 0.15),
               run_type = "all")

# runs a full simulation for true VE = 0.5 case
run_single_sim(seed = 1, 
               S_A1_mean_var = c(-1.29, 0.2),
               run_type = "all")


