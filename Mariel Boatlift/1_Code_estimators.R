########################   3. Prepare Data  #####################################
# Load packages
library(CVXR)
library(dplyr)
library(Metrics)

# Matrix with continous covariates
merge_all_continous_predictors <- function(
  continous_predictor, binary_predictor, treatment_indicators, time_periods){
  H_covariates_continous_2 <- rowSums(binary_predictor)/time_periods
  H_covariates_continous_2 <- ifelse(H_covariates_continous_2 > 0.49, 1, 0)
  H_covariates_continous <- as.matrix(rbind(rowSums(continous_predictor)/time_periods, H_covariates_continous_2))
  return( cbind(t(H_covariates_continous), treatment_indicators) )
}

# Matrix with discrete covariates
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
################################################################################





#####################   4. SIMULATE PSCM  ######################################
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
                     covariates_continous, untreated_treated_post_treatment){
  
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
  print(paste0("Doing Cross-Validation for PSCM takes more than half as much time as the rest of the code to run..."))
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
  est_ate <- mean(treated_post_treatment - t(t(donors_post_treatment) %*% pscm_weights) )
  est_sd_ate <- sd(treated_post_treatment - t(t(donors_post_treatment) %*% pscm_weights) )
  
  
  # 6. Calculate RMSE
  pscm_rmse <- rmse(untreated_treated_post_treatment, 
                    t(t(donors_post_treatment) %*% pscm_weights) )
  NewList <-list("est_ate"=est_ate, "est_sd_ate"=est_sd_ate, "pscm_rmse"=pscm_rmse)
}
################################################################################








################## 3. RWSCM WITH  CONTINOUS COV ################################
run_rwscm <- function(num_donors, num_treated, 
                      outcomes_matrix, covariates, untreated_treated_post_treatment){
  # Create matrices
  X_pre <- outcomes_matrix[outcomes_matrix$treatment_indicator < 0, c(1:(ncol(outcomes_matrix) - 2))]
  X_post <-  as.matrix(outcomes_matrix[outcomes_matrix$treatment_indicator < 0, (ncol(outcomes_matrix)-1) ])
  Y_pre <- outcomes_matrix[outcomes_matrix$treatment_indicator > 0, c(1:(ncol(outcomes_matrix) - 2))]
  Y_post <- as.matrix( outcomes_matrix[outcomes_matrix$treatment_indicator >0, (ncol(outcomes_matrix) -1)])
  Y_pre_long <- c(t(as.matrix(Y_pre)))
  
  # Outcome matrix for control group
  l_N <- matrix(1, ncol = 1, nrow = num_treated )
  l_K <- matrix(1, ncol = 1, nrow = num_donors)
  X_pre_full <- do.call(rbind,replicate(num_treated, t(X_pre),  simplify=FALSE))
  X_post_full <- do.call(rbind,replicate(num_treated, t(X_post),  simplify=FALSE))
  
  # Set up the (Q x nnum_treated) matrix of covariates
  QQ_cont <- ncol(covariates ) - 1
  covariates <- as.data.frame(covariates)
  H_covariates_treatment_sample <- t(covariates[covariates$treatment_indicators > 0, -ncol(covariates)])
  
  # Set Up the Constraints
  wei <- Variable(rows = num_donors)
  alphas <- Variable(rows = num_donors, cols = QQ_cont)
  
  # Intercepts
  O_NK <- matrix(0, ncol = num_treated, nrow = num_donors)
  l_KQ <- matrix(0, ncol = 1, nrow = num_donors * QQ_cont)
  
  # PROBLEM IS THE WEIRD CONSTRAINT
  objective <- Minimize(sum_squares(
    as.matrix(Y_pre_long) -
      X_pre_full %*% wei -
      vec( t(X_pre) %*% alphas %*% H_covariates_treatment_sample ) ) )
  constraints_int <- list(
    wei %*% t(l_N) + alphas %*%  H_covariates_treatment_sample >= O_NK,
    t(l_K) %*% (wei %*% t(l_N) + alphas %*% H_covariates_treatment_sample) == t(l_N) )
  print(paste0("Random Weights Synthetic Control is running"))
  problem <- Problem(objective, constraints_int)
  solution  <- psolve(problem, verbose = F, num_iter = 500000)
  # Next step is to calculate the weights for each observation
  if(solution$status == "solver_error"){
    print("upsy small fail equal weights sad times")
    equal_weights <- matrix(1/num_donors, ncol = num_treated, nrow = num_donors)
    est_ind_ate_rwscm_cont <- Y_post - t(t(X_post) %*% equal_weights) 
    ate <- mean(Y_post - t(t(X_post) %*% equal_weights) )
    sd_ate <- sd(Y_post - t(t(X_post) %*% equal_weights) )
    rmse <- rmse(untreated_treated_post_treatment , t(t(X_post) %*% equal_weights)  )
    return(list("est_ate" = ate, "sd_ate" = sd_ate, "rwscm_rmse" = rmse ))
    print(paste("Estimate ATE in RWSCM IS", ate))
  }
  else{
    
    alpha_QK <- solution$getValue(alphas)
    cons_weights <- solution$getValue(wei)
    
    # Calculate each weight
    ind_spec <- as.matrix(alpha_QK  %*% as.matrix(H_covariates_treatment_sample))
    ind_inv <- do.call(cbind, replicate(num_treated, cons_weights, simplify=FALSE))
    random_weights <- ind_inv + ind_spec
    
    # Output results
    est_ind_ate_rwscm_cont <- Y_post - t(t(X_post) %*% random_weights) 
    ate <- mean(Y_post - t(t(X_post) %*% random_weights) )
    sd_ate <- sd(Y_post - t(t(X_post) %*% random_weights) )
    rmse <- rmse(untreated_treated_post_treatment , t(t(X_post) %*% random_weights)  )
    return(list("est_ate" = ate, "sd_ate" = sd_ate, "rwscm_rmse" = rmse ))
  }
}

################################################################################








