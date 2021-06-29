########################   3. Prepare Data  #####################################
# Load packages
library(CVXR)
library(dplyr)
library(Metrics)
library(purrr)

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







################## 3. RWSCM WITH  CONTINOUS COV ################################
run_rwscm <- function(num_donors, num_treated, 
                      outcomes_matrix, covariates, untreated_treated_post_treatment,
                      true_ate, pre_treatment_length ){
  # Create matrices
  X_pre <- outcomes_matrix[outcomes_matrix$treatment_indicator < 0, 
                           c(1:(ncol(outcomes_matrix) - 2))]
  X_post <-  as.matrix(outcomes_matrix[outcomes_matrix$treatment_indicator < 0,
                                       (ncol(outcomes_matrix)-1) ])
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
  mu <- Variable(rows = num_treated)
  
  # Intercepts
  O_NK <- matrix(0, ncol = num_treated, nrow = num_donors)
  l_KQ <- matrix(0, ncol = 1, nrow = num_donors * QQ_cont)
  O1_n1Tn1 <-  matrix(0, ncol = num_treated, nrow = num_treated* (pre_treatment_length-1) )
  pre_treatment_length <- pre_treatment_length-1
  for( i in 1: (ncol(O1_n1Tn1)-1) ){
    O1_n1Tn1[ c((i*pre_treatment_length+1):(i*pre_treatment_length+pre_treatment_length)) , i+1] <- 1
  }
  
  # PROBLEM IS THE WEIRD CONSTRAINT
  objective <- Minimize(sum_squares(
    as.matrix(Y_pre_long) -
      O1_n1Tn1 %*% mu -
      X_pre_full %*% wei -
      vec( t(X_pre) %*% alphas %*% H_covariates_treatment_sample ) ) )
  constraints_int <- list(
    wei %*% t(l_N) + alphas %*%  H_covariates_treatment_sample >= O_NK,
    t(l_K) %*% (wei %*% t(l_N) + alphas %*% H_covariates_treatment_sample) == t(l_N),
    mu[1] ==0)
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
    mse_ate <- var(est_ind_ate_rwscm_cont) +  (ate - true_ate)^2
    return(list("est_ate" = ate, "mse_ate" = mse_ate, "rwscm_rmse" = rmse ))
  }
  else{
    alpha_QK <- solution$getValue(alphas)
    print(colSums(alpha_QK))
    cons_weights <- solution$getValue(wei)
    print(colSums(cons_weights))
    intercepts <- solution$getValue(mu)
    if(is.null(dim(alpha_QK)) & is.null(dim( intercepts) ) & is.null(dim(cons_weights))){
      print("upsy small fail equal weights sad times")
      equal_weights <- matrix(1/num_donors, ncol = num_treated, nrow = num_donors)
      est_ind_ate_rwscm_cont <- Y_post - t(t(X_post) %*% equal_weights) 
      ate <- mean(Y_post - t(t(X_post) %*% equal_weights) )
      sd_ate <- sd(Y_post - t(t(X_post) %*% equal_weights) )
      rmse <- rmse(untreated_treated_post_treatment , t(t(X_post) %*% equal_weights)  )
      mse_ate <- var(est_ind_ate_rwscm_cont) +  (ate - true_ate)^2
      return(list("est_ate" = ate, "mse_ate" = mse_ate, "rwscm_rmse" = rmse ))
    }
    
    # Manage intercepts
    ones_matrix <- matrix(1, ncol = pre_treatment_length, nrow = num_treated)
    all_intercepts <- do.call(cbind, replicate(pre_treatment_length, intercepts, simplify=FALSE))
    
    # Calculate each weight
    print("Head of alpha_QK")
    print(head(alpha_QK))
    print("Head of H_covariates_treatment_sample")
    print(head(H_covariates_treatment_sample))
    ind_spec <- as.matrix(alpha_QK  %*% as.matrix(H_covariates_treatment_sample))
    ind_inv <- do.call(cbind, replicate(num_treated, cons_weights, simplify=FALSE))
    print("Individual specific head")
    print(head(ind_spec))
    random_weights <- ind_inv + ind_spec
    print("Random Weights head")
    print(head(random_weights))
    # Output results
    synthetic_controls <-  t(t(X_post) %*% random_weights) + intercepts
    est_ind_ate <- Y_post -synthetic_controls
    ate <- mean(est_ind_ate )
    sd_ate <- sd(est_ind_ate )
    rmse <- rmse(untreated_treated_post_treatment,  synthetic_controls  )
    mse_ate <- var(est_ind_ate) +  (ate - true_ate)^2
    return(list("est_ate" = ate, "mse_ate" = mse_ate, "rwscm_rmse" = rmse ))
  }
}

################## 3. RWSCM WITH  CONTINOUS COV ################################
run_did <- function(num_donors, num_treated, outcomes_matrix, 
                    untreated_treated_post_treatment, pre_treatment_length,
                    true_ate){
  # Pre-treatment length 
  TT <- ncol(outcomes_matrix)-1
  post_treatment_period <- TT - pre_treatment_length
  N <- num_treated + num_donors
  
  # Create matrices
  outcomes_matrix <- as.data.frame(outcomes_matrix)
  outcomes_matrix <- outcomes_matrix %>% arrange(treatment_indicator)
  # The idea is that we put in the end the treated units, so the matrix has the 
  # same form as in the slides
  all_pre <- outcomes_matrix[, c((1:(ncol(outcomes_matrix) - 1)) )]
  X_post <-  as.matrix(outcomes_matrix[outcomes_matrix$treatment_indicator < 0, (ncol(outcomes_matrix)-1) ])
  Y_post <- as.matrix( outcomes_matrix[outcomes_matrix$treatment_indicator > 0, (ncol(outcomes_matrix) -1)])
  Y_long <- c(t(as.matrix(all_pre)))
  
  # Set Up the Constraints
  alpha <- Variable(1)
  gammas <- Variable(rows = N)
  deltas <- Variable(rows = TT)
  tau <- Variable(1)
  
  # Intercepts
  O_NK <- matrix(0, ncol = num_treated, nrow = num_donors)
  
  # Dummies for individual effect
  O1_NTxN <-  matrix(0, ncol =  N, nrow =  N*  TT)
  for( i in 1: (ncol(O1_NTxN)-1) ){
    O1_NTxN[ c((i*TT+1):(i*TT+TT)) , i+1] <- 1
  }
  
  # Dummies for time effect
  O1_NTxT <-  matrix(0, ncol = TT, nrow = N*TT)
  for( i in 1: (ncol(O1_NTxT)-1) ){
    O1_NTxT[seq(from = 1+i, to = nrow(O1_NTxT), by =  TT), i+1] <- 1
  }
  
  # Treatment indicator
  D <- matrix(0, ncol = 1, nrow = length(Y_long))
  D[seq(from = length(Y_long) - TT * num_treated, to = length(Y_long), by =  TT)] <- 1
  # PROBLEM IS THE WEIRD CONSTRAINT
  objective <- Minimize(sum_squares(
    as.matrix(Y_long) -  alpha - O1_NTxN %*% gammas -
      O1_NTxT %*% deltas -
      D %*% tau) ) 
  constraints_int <- list( gammas[1] == 0, deltas[1] == 0 )
  problem <- Problem(objective, constraints_int)
  solution  <- psolve(problem, verbose = F, num_iter = 500000)
  # Next step is to calculate the weights for each observation
  if(solution$status == "solver_error"){
    print("upsy small fail of covergence for did sad times")
    equal_weights <- matrix(1/num_donors, ncol = num_treated, nrow = num_donors)
    est_ind_ate_rwscm_cont <- Y_post - t(t(X_post) %*% equal_weights) 
    ate <- mean(Y_post - t(t(X_post) %*% equal_weights) )
    sd_ate <- sd(Y_post - t(t(X_post) %*% equal_weights) )
    rmse <- rmse(untreated_treated_post_treatment , t(t(X_post) %*% equal_weights)  )
    mse_ate <- var(est_ind_ate_rwscm) + (mean(est_ind_ate_rwscm) - true_ate)^2
    return(list("est_ate" = ate, "sd_ate" = sd_ate, "did_rmse" = rmse, "mse_ate" =  mse_ate ))
  }
  else{
    gammas_treated <- solution$getValue(gammas)[c( (N - num_treated+1) :N)]
    predicted <- gammas_treated + solution$getValue(deltas)[TT] + solution$getValue(alpha)
    
    # Output results
    est_ind_ate_rwscm <- Y_post - predicted 
    mse_ate <- var(est_ind_ate_rwscm) + (mean(est_ind_ate_rwscm) - true_ate)^2
    ate <- mean(est_ind_ate_rwscm)
    sd_ate <- sd(est_ind_ate_rwscm)
    rmse <- rmse(untreated_treated_post_treatment, predicted) 
    return(list("est_ate" = ate, "sd_ate" = sd_ate, 
                "did_rmse" = rmse, "mse_ate" =  mse_ate))
  }
}
################################################################################







##################### 5. Combine into a single function ########################
########################   3. Prepare Data  #####################################

##################### 5. Combine into a single function ########################
########################   3. Prepare Data  #####################################
calculate_performance_metrics <- function(simulated_data){
  # Derive objects
  y_mat <- simulated_data[["treated_outcome_variable_mat"]]
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
  
  # Matrix with predictors
  H_covariates  <- cbind(covariates_single , y_mat$treatment_indicator)
  colnames(H_covariates)[ncol(H_covariates)] <- "treatment_indicators"
  # Continous predictors RWSCM
  rwscm_res <-  run_rwscm(num_donors = n0, num_treated = n1, 
                          outcomes_matrix = y_mat, covariates = H_covariates,
                          untreated_treated_post_treatment = untreated_y_treatment_gr_vec,
                          true_ate = true_ate, pre_treatment_length= time_periods-1)
  
  # Feasible Difference in Difference
  fdid_res <- run_did(num_donors = n0, num_treated= n1, outcomes_matrix= y_mat, 
                      untreated_treated_post_treatment= untreated_y_treatment_gr_vec, 
                      pre_treatment_length = time_periods-1, true_ate = true_ate)
  
  # Let us 'correct' the outcome amtrix
  tilde_mu_i <- simulated_data[["factor_loadings"]] - simulated_data[["mu_mean"]]
  tilde_lambda_t <- simulated_data[["common_factors"]] - simulated_data[["lambda_mean"]]
  # Calculate correction
  correction <- tilde_mu_i %*% t(tilde_lambda_t)
  corrected_y_mat <- y_mat[,c(1:(ncol(y_mat) -1 ))] - correction
  corrected_y_mat <- cbind(corrected_y_mat, y_mat$treatment_indicator)
  corrected_y_mat <- as.data.frame(corrected_y_mat)
  colnames(corrected_y_mat)[ncol(corrected_y_mat)] <- "treatment_indicator"
  # Calculate correction of true outcome variable
  corrected_untreated_y_treatment_gr_vec <-  untreated_y_treatment_group[,c(1:(ncol( untreated_y_treatment_group ) -1 ))] - correction
  colnames(untreated_y_treatment_group)[ncol(untreated_y_treatment_group)] <- "treatment_indicator"
  corrected_untreated_y_treatment_gr_vec <- cbind(corrected_untreated_y_treatment_gr_vec,  untreated_y_treatment_group$treatment_indicator)
  corrected_untreated_y_treatment_gr_vec <- as.data.frame(corrected_untreated_y_treatment_gr_vec)
  corrected_untreated_y_treatment_gr_vec <- corrected_untreated_y_treatment_gr_vec[
    corrected_untreated_y_treatment_gr_vec[, ncol(corrected_untreated_y_treatment_gr_vec)]>0, 
    ncol(corrected_untreated_y_treatment_gr_vec)-1] 
  
  # Feasible Difference in Difference
  idid_res <- run_did(num_donors = n0, num_treated= n1, outcomes_matrix= corrected_y_mat, 
                      untreated_treated_post_treatment= corrected_untreated_y_treatment_gr_vec, 
                      pre_treatment_length = time_periods-1, true_ate = true_ate)
  
  # Print results
  perf_measures <- c(
    rwscm_res[["est_ate"]],  rwscm_res[["mse_ate"]],  rwscm_res[["rwscm_rmse"]], 
    fdid_res[["est_ate"]],  fdid_res[["mse_ate"]],  fdid_res[["did_rmse"]],
    idid_res[["est_ate"]],  idid_res[["mse_ate"]],  idid_res[["did_rmse"]])
  
  return(perf_measures)
}
#################################################################################


