# Simulation code for "Surrogate Endpoint Based Provisional Approval Causal Roadmap"

This repository contains code for "A Surrogate Endpoint Based Provisional Approval Causal Roadmap" by 
Gilbert et al. (under review). 

Simulation code is in the `ms_sim.R` file. Most of the code consists of helper functions. The last section "Some Sample
Code to Try" has a few lines of code to generate simulated datasets, run the estimators, obtain bootstrap and sandwich 
standard errors, and run a single simulation.

The pre-generated simulated datasets include `sim_obs_data.csv` (observational dataset) and `sim_p3_data.csv` (phase 3 dataset). 

The code to analyze a single dataset is in `analyze_single_dataset.R`. More details are included in the vignette in G.3 in the paper's
supplementary materials. Some starter code to analyze the single dataset is shown below: 

```
source("analyze_single_dataset.R")

df_obs <- read.csv("sim_obs_data.csv")
df_p3 <- read.csv("sim_p3_data.csv")

estimate_ve_surrogate(df_obs, 
                      df_p3,
                      treatment = "A", 
                      outcome = "Y", 
                      covariates = c("X1", "X2", "X3"), 
                      surrogate = "S", 
                      estimate_weights = T, 
                      weight = NULL, 
                      estimator = "onestep", 
                      learner = "glm", 
                      ct_bias_values = c(0,0.001), 
                      uc_bias_values = c(0,-0.001), 
                      success_criteria = c(0.3, NA))

```

If you would like to get in touch about this code, please email James Peng at jpspeng@uw.edu. 
