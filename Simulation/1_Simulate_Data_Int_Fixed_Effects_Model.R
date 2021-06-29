# The DGP for Interactive Random Effects is:
# y^{D=0}_it = \beta' x_i + \mu'_i \lambda_t + u_it
# 1) \mu'_i - stochastic coefficient s.t. E[\mu'_i] \neq 0; E[\mu'_i x_i] \neq 0
#       \mu'_i = \gamma x_i + v_i
#       v_i = N(v,1)
# 2) $\lambda_t$ - fixed parameter; $ \frac{\sum_t f_t}{T} \neq 0 $
#       \lambda_t = N(f,1)
# 3) Strict Exogeneity: \forall t: \ E[u_{it}x_i] = 0
#        u_it = N(0,1)
# x1_i = U[{1,2,3,4}]
# x2_i = N(0,1)


# Functions for treating individuals
treatment_indicator_not_random <- function(prob_tr,  factor_loadings, gamma, covariates_single,  num_obs, seed = 3){
  
  factor_loadings <- rowMeans(factor_loadings)
  e_i <- rnorm(nrow(gamma), mean = 0, sd = 0.1 )
  gamma <- as.matrix(rowMeans(gamma) +  e_i)
  set.seed(seed*4)
  u_i <- rnorm(nrow(covariates_single), sd = 1)
  pi_i <- (exp(factor_loadings +  covariates_single %*% gamma +   u_i) * (exp(factor_loadings + covariates_single %*% gamma +   u_i) + 1) )^(-1) 
  pi_i  <- pi_i + (prob_tr - mean(pi_i))
  set.seed(seed*3)
  random <- sample(length(pi_i), 1)
  pi_i[random] <- 1
  pi_i <- ifelse(pi_i == max(pi_i), 1, pi_i )
  set.seed(seed)
  D_i <- rbernoulli(n=num_obs, p =  pi_i )
  return(D_i)
}
# 2.8. Generate treatment effect
treat_some_observations <- function(matrix_outcomes_all_obs, treatment_indicator,
                                    ate, sd_ate,
                                    seed = 3){
  set.seed(seed)
  y_mat_temp <- cbind(matrix_outcomes_all_obs, treatment_indicator)
  e_it <- rnorm(sum(treatment_indicator), mean = 0, sd = sd_ate)
  y_mat_temp <- as.matrix(y_mat_temp)
  y_mat_temp[ y_mat_temp[, ncol(y_mat_temp)] == 1, ncol(y_mat_temp)-1]  <- 
  y_mat_temp[ y_mat_temp[, ncol(y_mat_temp)] == 1,  ncol(y_mat_temp)-1] + ate + e_it
  return( as.data.frame(y_mat_temp))
}

simulate_data_int_random_effects <-function(
  time_periods, num_obs, num_factors, 
  lambda_mean, mu_mean, prob_tr, 
  seeds, num_categories,
  beta_coef_outcome_variable,  
  ate, sd_ate,
  treat_not_at_random = 1
){
  # Exogenous shocks u_it
  set.seed(seeds[1])
  u_it <- rnorm(num_obs*time_periods, mean = 0, sd = 1)
  
  # Common Factors lambda_t 
  set.seed(seeds[2])
  lambda_t <- matrix( 0, ncol = num_factors, nrow = time_periods )
  for(i in 1:num_factors){
    set.seed(seeds[2]*2)
    exp_lambda  <- sample(c((-3): 3), 1, replace =  T)
    lambda_t[,i] <- rnorm( time_periods, mean = exp_lambda + lambda_mean)
  }

  # Continous Covariate 
  # set.seed(seeds[3])
  # xi1_vec <- rnorm(num_obs, mean = 0, sd =  1)
  
  # Discrete Covariate 
  set.seed(seeds[4])
  xi2_vec <- sample(num_categories, size = num_obs, replace = T)
  xi2_vec <- ifelse(xi2_vec > num_categories, num_categories, xi2_vec)
  xi2_vec <- ifelse(xi2_vec < 1, 1, xi2_vec)
  xi2_vec <- as.data.frame(xi2_vec)
  xi2_dummy_mat_single <- model.matrix(~ as.factor(xi2_vec) + 0, xi2_vec)
  xi2_vec <- as.matrix(xi2_vec)
  xi2_dummy_mat_single <- as.matrix(xi2_dummy_mat_single)
  # covariates_single <- cbind(xi1_vec,  xi2_dummy_mat_single)
  covariates_single <- cbind( xi2_dummy_mat_single)
  covariates_single <- as.matrix(covariates_single)
  
  # Generate factor 
  mu_i <- matrix(0, ncol = num_factors, nrow = num_obs)
  v_i <- matrix(0, ncol = num_factors, nrow = num_obs )
  set.seed(seeds[5])
  gamma <- matrix( rnorm(ncol(covariates_single) * num_factors ), ncol = num_factors, nrow = ncol(covariates_single) )
  for(i in 1:num_factors){
    set.seed(seeds[6] * i)
    expect_mu <- mu_mean + sample(c((-3): 3), 1, replace =  T)
    set.seed(seeds[9] * i)
    v_i[,i] <-  rnorm( num_obs, mean = expect_mu  )
    mu_i[,i] <-   covariates_single %*% gamma[, i] + v_i[, i]
  }


  # Add panel dimension to covariates
  # xi1_mat <- do.call(cbind, replicate(time_periods, xi1_vec, simplify=FALSE))
  # xi1_vec <- c(t(xi1_mat))
  xi2_mat <- do.call(cbind, replicate(time_periods, xi2_vec, simplify=FALSE))
  xi2_vec <- c(t(xi2_mat))
  xi2_vec <- as.data.frame(xi2_vec)
  xi2_dummy_mat <- model.matrix(~ as.factor(xi2_vec) + 0, xi2_vec)
  xi2_vec <- as.matrix(xi2_vec)
  xi2_dummy_mat <- as.matrix(xi2_dummy_mat)
  # covariates <- as.matrix(cbind(xi1_vec, xi2_dummy_mat))
  covariates <- as.matrix(cbind(xi2_dummy_mat))
  xi2_dummy_mat <- as.data.frame(xi2_dummy_mat)
  
  # Generate outcome variable 
  factor_structure <- c( t(mu_i %*% t(lambda_t) ) ) 
  beta_coef_outcome_variable <- beta_coef_outcome_variable[-1]
  yit_vec <-  covariates %*% beta_coef_outcome_variable + factor_structure + u_it
  yit_mat <- t(matrix(yit_vec, nrow = time_periods, ncol =num_obs))
  
  # Treat some observations not-at-random
  if(treat_not_at_random == 1){
    D_i <- treatment_indicator_not_random( 
      prob_tr = prob_tr, factor_loadings = mu_i, num_obs = num_obs, seed = seeds[7],
      gamma = gamma, covariates_single = covariates_single)
  }
  else{
    set.seed(seeds[7])
    D_i <- rbernoulli(num_obs, p = prob_tr)
  }
  
  # Generate treatment matrix with
  yit_mat_temp <- treat_some_observations( 
    matrix_outcomes_all_obs = yit_mat, treatment_indicator = as.numeric(D_i),
    ate = ate, sd_ate = sd_ate, seed = seeds[8])

  newList <- list(
   #  "continous_predictor_mat" = xi1_mat,  
    "discrete_predictor_mat" = xi2_mat, 
    "treated_outcome_variable_mat" = yit_mat_temp, 
    "untreated_outcome_variable_mat" = yit_mat, 
    "interactive_random_effects" = t(mu_i %*% t(lambda_t) ),
    "common_factors" = lambda_t, "factor_loadings" = mu_i,
    "treatment_matrix" = D_i, "x2_dummy_mat" = xi2_dummy_mat_single,
    "covariates_single" = covariates_single, "true_ate" = ate, 
    "lambda_mean" = lambda_mean, "mu_mean" = mu_mean)
}

