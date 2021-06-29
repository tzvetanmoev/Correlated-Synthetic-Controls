# ########################   1. Choice variables  ################################
# Dimensions of matrix. 
# Note originally I wanted to follow the dimensions for Sturm and Redding but
# with treatment being non-randomly assigned this is not possible

# Factors tunning parameter

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
tunning_results <- matrix(0, nrow = 8, ncol = 10)
tunning_results[,1] <- c("baseline", "N = 150", "N = 200",  
                         "T = 4",  "T = 7", 
                        # "num_categ = 6", "num_categ = 7",
                        "FF = 3", "FF = 4",
                        # "prob_treated = 0.2", "prob_treated = 0.3", 
                        "treatment_at_random")

colnames(tunning_results) <- c("change", 
                               "rwscm_ate", "rwscm_mse_ate", "rwscm_rmse",
                               "fdid_ate", "fdid_mse_ate", "fdid_rmse",
                               "idid_ate", "idid_mse_ate", "idid_rmse")

# Matrixt to use in the tunning 
tunning_params <- matrix(0 , nrow = 8, ncol = 7 )
# Note we change number of columns, so the beta's decrease
colnames(tunning_params) <- c( "iteration_detail",  "N", "T",
                                "num_categ", "num_factors", "prob_treated",
                                "treat_not_at_random")
tunning_params[,1] <- c("baseline", "N = 150", "N = 200",  
                        "T = 4",  "T = 7", 
                        # "num_categ = 6", "num_categ = 7",
                        "FF = 3", "FF = 4",
                        # "prob_treated = 0.2", "prob_treated = 0.3", 
                        "treatment_at_random")

tunning_params[1, c(2:7)] <- c(100, 5, 5, 2, 0.15,  1 )
tunning_params[2, c(2:7)] <- c(150, 5, 5, 2, 0.15,  1 )
tunning_params[3, c(2:7)] <- c(200, 5, 5, 2, 0.15,  1 )
tunning_params[4, c(2:7)] <- c(100, 4, 5, 2, 0.15,  1 )
tunning_params[5, c(2:7)] <- c(100, 7, 5, 2, 0.15,  1 )
# tunning_params[6, c(2:7)] <- c(200, 5, 6, 2, 0.1,  1 )
# tunning_params[7, c(2:7)] <- c(200, 5, 7, 2, 0.1,  1 )
tunning_params[6, c(2:7)] <- c(100, 5, 5, 3, 0.15,  1 )
tunning_params[7, c(2:7)] <- c(100, 5, 5, 4, 0.15, 1 )
# tunning_params[10, c(2:7)] <- c(200, 5, 5, 2, 0.2,  1 )
# tunning_params[11, c(2:7)] <- c(200, 5, 5, 2, 0.3,  1 )
tunning_params[8, c(2:7)] <- c(100, 4, 5, 2, 0.15,  0 )

# tunning_params[12, c(2:9)] <- c(200, 4, 4, 2, 0.1, 3, 1, 1 )
# tunning_params[13, c(2:9)] <- c(200, 4, 4, 2, 0.1, 4, 1, 1 )
# tunning_params[14, c(2:9)] <- c(200, 4, 4, 2, 0.1, 2, 2, 1 )
# tunning_params[15, c(2:9)] <- c(200, 4, 4, 2, 0.1, 2, 3, 1 )

num_runs <- 2
for( j in 1:2){
  print("????????????????????????????????????????????????????????????????")
  print("????????????????????????????????????????????????????????????????")
  print("????????????????????????????????????????????????????????????????")
  print(paste0("????  now working on simulation where we change ", tunning_params[j,1], "  ???"))
  print("????????????????????????????????????????????????????????????????")
  print("????????????????????????????????????????????????????????????????")
  print("????????????????????????????????????????????????????????????????")
  rm(results_simulation)
  results_simulation <- matrix(0, ncol = 10, nrow = num_runs)
  colnames(results_simulation) <- c("change", 
                                    "rwscm_ate", "rwscm_mse_ate", "rwscm_rmse",
                                    "fdid_ate", "fdid_mse_ate", "fdid_rmse",
                                    "idid_ate", "idid_mse_ate", "idid_rmse")
  results_simulation[,1] <- c(1:num_runs)
  
  for(i in c(1:num_runs)){
    if( (i == 100) | (i == 200) | (i == 400)  | (i == 600)  | (i == 800)    ){
      print(paste0("we've reached succesfully ", i, "th iteration for this simulation"))
    }
    set.seed(i)
    uniform_variable <-  matrix(runif(as.numeric(tunning_params[j, 4])+1, 
                                        min=0, max = 1),
                                  ncol = 1, nrow = as.numeric(tunning_params[j, 4]) + 1)
    set.seed(i+500)
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
    performance <- t(as.matrix(calculate_performance_metrics(simulated_data)))
    results_simulation[i,c(2:10)] <- performance
  }
  tunning_results[j, c(2:10)] <- round(colMeans( results_simulation[, c(2:10)]), 2)
  
  assign(  paste("results_simulation_", j, sep = ""), results_simulation )
}


#     1) Clearly, including discrete predictor in treatment makes no difference
#     2) Performance clearly improves as we increase N 
#     3) Common shocks has a dramatic effect on the coefficients
#     4) Increasing number of categories improves performance for all estimators,
#        apart from the continuous rwscm (not sure why for pscm?)
#     5) Imposing $\forall i,j: \ \beta^{3i} =\beta^{3i}$  )
#     6) Decreasing \beta_1 leads to rwscm-s performing better than in baseline
#     7) Changing betas have ambiguous effects with nothing very clear. However,
#        moving beta towards 0 clearly makes all of them perform better
#     8) Lastly, note that MSE of tau is smaller for rwscm than PSCM





