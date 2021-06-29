setwd("~/Documents/Oxford Yr 2/THESIS_micro_synth/Mariel boatlift shock/Analysis Data")
library(purrr)
library(Metrics)
library(dplyr)
library(CVXR)

pscm <- function(y, X, lambda = 1, nn_1 = 23){
  wei <- Variable(rows = ncol(X))
  nn_0 <-  ncol(X)
  objective <- Minimize(sum_squares(as.matrix(y) - X %*%  wei  ) + lambda * 
                          t(sum_entries(( as.matrix(y) %*% as.matrix(t(rep(1, nn_0))) - X )^2, axis = 2)) %*%  wei )
  constraints <- list(wei  >= rep(0, nn_0), t(as.matrix(rep(1, nn_0))) %*% wei == 1)
  problem <- Problem(objective, constraints)
  solution  <- psolve(problem, num_iter = 150000)
  if(solution$status == "solver_error"){
    print(paste0("upsy fail of pscm"))
    return(matrix(1/ncol(X), ncol = 1, nrow = ncol(X)))
  }
  return(solution$getValue(wei))
}
run_pscm <- function(num_donors, num_treated, outcome_variable_matrix,
                     donors_pre_treatment, treated_pre_treatment,
                     donors_post_treatment, treated_post_treatment,
                     covariates_continous, untreated_treated_post_treatment,
                     true_ate ){
  
  # 0. Preparate the data
  set.seed(1)
  stack_X_pre_train <- t(cbind(  donors_pre_treatment[, c(1:(ncol(donors_pre_treatment)-1 ))], 
                                 covariates_continous[outcome_variable_matrix$treatment_indicator < 0,c(1:2) ]))
  stack_Y_pre_train <- t(cbind(   treated_pre_treatment[, c(1:(ncol(donors_pre_treatment)-1 ))], 
                                  covariates_continous[outcome_variable_matrix$treatment_indicator > 0,c(1:2) ]))
  stack_X_pre_test <- t(cbind(  donors_pre_treatment[, (ncol(donors_pre_treatment) )], 
                                covariates_continous[outcome_variable_matrix$treatment_indicator < 0,c(1:2) ]))
  stack_Y_pre_test<- t(cbind(   treated_pre_treatment[, (ncol(donors_pre_treatment) )], 
                                covariates_continous[outcome_variable_matrix$treatment_indicator > 0,c(1:2) ]))
  
  # 1. Create object to store the results from CV
  pscm_weights <- matrix(0, ncol = num_treated, nrow = num_donors)
  rmse_pscm<- matrix(0, ncol = 2, nrow = 3)
  rmse_pscm[,1] <- c(0.4, 0.8, 6.4)
  
  # 2. Perform Cross-Validation
  for(j in 1:3){
    for(i in 1:num_treated){
      pscm_weights[,i] <- pscm(y = as.matrix(stack_Y_pre_train[,i]), 
                               X = as.matrix(stack_X_pre_train), lambda = rmse_pscm[j,1],
                               nn_1 = ncol(stack_Y_pre_train) )
    }
    # Predict values
    rmse_pscm[j, 2] <- rmse(c(stack_X_pre_test %*% pscm_weights), 
                            c(stack_Y_pre_test) )
  }
  
  # 3. Get the optimal lambda
  rmse_pscm <- as.data.frame(rmse_pscm)
  optimal_lambda <- (rmse_pscm %>% arrange(V2) %>% slice(1))[1]
  print(paste0("optimal lambda is ", optimal_lambda))
  
  # 4. Restimate with optimal lambda
  pscm_weights <- matrix(0, ncol = num_treated, nrow = num_donors)
  for(i in 1:num_treated){
    pscm_weights[,i] <- pscm( 
      y = stack_Y_pre_train[,i], X = as.matrix(stack_X_pre_train), 
      lambda = optimal_lambda)
  }
  
  # 5. Calculate the ATE-s 
  est_ind_ate <- treated_post_treatment - t(t(donors_post_treatment) %*% pscm_weights)
  est_ate <- mean(est_ind_ate  )
  
  # 6. Calculate RMSE
  pscm_rmse <- rmse(untreated_treated_post_treatment, 
                    t(t(donors_post_treatment) %*% pscm_weights) )
  mse <- (true_ate - est_ate)^2 + var(est_ind_ate)
  NewList <-list("pscm_ate"=est_ate, "pscm_mse"= mse, "pscm_rmse"=pscm_rmse)
}
################################################################################


pscm_calculate_performance_metrics <- function(simulated_data){
  
  # Derive objects
  y_mat <- simulated_data[["treated_outcome_variable_mat"]]
  x1_mat <- simulated_data[["continous_predictor_mat"]]
  x2_mat <- simulated_data[["x2_dummy_mat"]]
  covariates_single <- simulated_data[["covariates_single"]]
  true_ate <- simulated_data[["true_ate"]]
  
  
  
  # True value in the matrix
  untreated_y_mat <- simulated_data[["untreated_outcome_variable_mat"]]
  untreated_y_treatment_group <- as.data.frame(cbind(untreated_y_mat, y_mat$treatment_indicator))
  untreated_y_treatment_gr_vec <- untreated_y_treatment_group[
    untreated_y_treatment_group[, ncol(untreated_y_treatment_group)]>0, 
    ncol(untreated_y_treatment_group)-1]
  
  # Note we use the coding from the old unsuccessful first attempt at a simulation
  n0 <- length(y_mat$treatment_indicator[y_mat$treatment_indicator < 0.5] )
  n1 <- length(y_mat$treatment_indicator[y_mat$treatment_indicator > 0.5] )
  y_mat$treatment_indicator[y_mat$treatment_indicator == 0] <- -c(1:n0)
  y_mat$treatment_indicator[y_mat$treatment_indicator == 1] <- c(1:n1)
  time_periods <- ncol(y_mat)
  
  # Form the matrices with outcome values for control and treated group
  X_pre <- y_mat[y_mat$treatment_indicator < 0, c(1:(ncol(y_mat)-2))]
  X_post <-  y_mat[y_mat$treatment_indicator < 0, (ncol(y_mat)-1)]
  Y_pre <- y_mat[y_mat$treatment_indicator > 0, c(1:(ncol(y_mat)-2))]
  Y_post <-  y_mat[y_mat$treatment_indicator >0, (ncol(y_mat)-1)]
  
  # Matrix with continous covariates
  # Matrix with predictors
  H_covariates  <- cbind(covariates_single , y_mat$treatment_indicator)
  colnames(H_covariates)[ncol(H_covariates)] <- "treatment_indicators"
  H_covariates_continous <- cbind(covariates_single[,1], y_mat$treatment_indicator)
  
  # PSCM continous
  # pscm_results_cont <- run_pscm(num_donors = n0, num_treated = n1,
  #                               outcome_variable_matrix = y_mat,
  #                               donors_pre_treatment = X_pre, treated_pre_treatment = Y_pre,
  #                               donors_post_treatment = X_post, treated_post_treatment = Y_post,
  #                               untreated_treated_post_treatment = untreated_y_treatment_gr_vec,
  #                               covariates_continous = as.data.frame(H_covariates_continous),
  #                               true_ate = true_ate)
  
  # PSCM all
  pscm_results_all <- run_pscm(num_donors = n0, num_treated = n1,
                               outcome_variable_matrix = y_mat,
                               donors_pre_treatment = X_pre, treated_pre_treatment = Y_pre,
                               donors_post_treatment = X_post, treated_post_treatment = Y_post,
                               untreated_treated_post_treatment = untreated_y_treatment_gr_vec,
                               covariates_continous = as.data.frame(H_covariates), 
                               true_ate = true_ate)
  
  # Print results
  perf_measures <- c(
    pscm_results_all[["pscm_ate"]], pscm_results_all[["pscm_mse"]], pscm_results_all[["pscm_rmse"]])
  
  return(perf_measures)
}

# Treatment
ate <- 1 
sd_ate <- 0.3

# Mean
mu_mean <- 1
lambda_mean <- 2

################################################################################




################################################################################
# Matrix to save the results
# Paramaters we will change are T, N, 
pscm_tunning_results <- matrix(0, nrow = 8, ncol = 4)
pscm_tunning_results[,1] <- c("baseline", "N = 150", "N = 200",  
                              "T = 4",  "T = 7", 
                              "F = 3", "F = 4",
                              "treatment_at_random")

colnames(pscm_tunning_results) <- c("change", 
                                    "pscm_all_ate", "pscm_all_mse_ate", "pscm_all_rmse")

# Matrixt to use in the tunning 
tunning_params <- matrix(0 , nrow = 8, ncol = 7 )
# Note we change number of columns, so the beta's decrease
colnames(tunning_params) <- c( "iteration_detail",  "N", "T",
                               "num_categ", "num_factors", "prob_treated",
                               "treat_not_at_random")
tunning_params[,1] <- c("baseline", "N = 150", "N = 200",  
                        "T = 4",  "T = 7", 
                        "F = 3", "F = 4",
                        "treatment_at_random")

tunning_params[1, c(2:7)] <- c(100, 5, 5, 2, 0.15,  1 )
tunning_params[2, c(2:7)] <- c(150, 5, 5, 2, 0.15,  1 ) 
tunning_params[3, c(2:7)] <- c(200, 5, 5, 2, 0.15,  1 )
tunning_params[4, c(2:7)] <- c(100, 4, 5, 2, 0.15,  1 )
tunning_params[5, c(2:7)] <- c(100, 7, 5, 2, 0.15,  1 )
tunning_params[6, c(2:7)] <- c(100, 5, 5, 3, 0.15,  1 )
tunning_params[7, c(2:7)] <- c(100, 5, 5, 4, 0.15,  1 )
tunning_params[8, c(2:7)] <- c(100, 5, 5, 2, 0.15,  0 )

# tunning_params[12, c(2:9)] <- c(200, 4, 4, 2, 0.1, 3, 1, 1 )
# tunning_params[13, c(2:9)] <- c(200, 4, 4, 2, 0.1, 4, 1, 1 )
# tunning_params[14, c(2:9)] <- c(200, 4, 4, 2, 0.1, 2, 2, 1 )
# tunning_params[15, c(2:9)] <- c(200, 4, 4, 2, 0.1, 2, 3, 1 )

num_runs <- 250
for( j in 5:8){
  print("????????????????????????????????????????????????????????????????")
  print("????????????????????????????????????????????????????????????????")
  print("????????????????????????????????????????????????????????????????")
  print(paste0("????  now working on simulation where we change ", tunning_params[j,1], "  ???"))
  print("????????????????????????????????????????????????????????????????")
  print("????????????????????????????????????????????????????????????????")
  print("????????????????????????????????????????????????????????????????")
  rm(results_simulation)
  results_simulation <- matrix(0, ncol = 4, nrow = num_runs)
  colnames(results_simulation) <- c("change", 
                                    "pscm_all_ate", "pscm_all_mse_ate", "pscm_all_rmse")
  results_simulation[,1] <- c(1:num_runs)
  
  for(i in 1:num_runs){
    if( (i == 200) | (i == 400) | (i == 600) | (i == 800) ){
      print(paste0("we've reached succesfully ", i, "th iteration for this simulation"))
    }
    set.seed(i)
    uniform_variable <-  matrix(runif(as.numeric(tunning_params[j, 4])+1, 
                                      min=0, max = 1),
                                ncol = 1, nrow = as.numeric(tunning_params[j, 4]) + 1)
    set.seed(i+1000)
    seeds <- sample( 20000, 9)
    simulated_data <- simulate_data_int_random_effects(
      time_periods = as.numeric(tunning_params[j, 3]), 
      num_obs = as.numeric(tunning_params[j, 2]), 
      num_factors = as.numeric(tunning_params[j, 5]), 
      prob_tr = as.numeric(tunning_params[j, 6]),
      lambda_mean = lambda_mean ,  
      mu_mean = mu_mean,
      seeds = seeds, num_categories = as.numeric(tunning_params[j, 4]),
      beta_coef_outcome_variable = uniform_variable, 
      ate = ate, sd_ate = sd_ate,
      treat_not_at_random = as.numeric(tunning_params[j, 7])
    )
    performance <- t(as.matrix(pscm_calculate_performance_metrics(simulated_data)))
    results_simulation[i,c(2:4)] <- performance
  }
  pscm_tunning_results[j, c(2:4)] <- round(colMeans( results_simulation[, c(2:4)]), 2)
  
  assign(  paste("results_simulation_", j, sep = ""), results_simulation )
}
