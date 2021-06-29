# install.packages("tensor")
setwd("~/Documents/Oxford Yr 2/THESIS_micro_synth/Mariel boatlift shock/Analysis Data")
library(dplyr)
library(Matrix)
library(CVXR)
library(Metrics)

# The purpose of this file is to:
#     1) Prepare the data on labour supply and divide it into: 
#         i) training period (1972-1978);
#         ii) testing period (1979-1980);
#         iii) post-treatment period (1981-1984)
#     2) Run CSC on the data for the training period
#     3) Compute synthetic controls for the testing period and get rmse
#     4) Get ATT for the post-treatment period
#

################## 1. Prepare the Data  #######################################
##################  1.1. TREATMENT GROUP #######################################
ls_treated_outcome <- readRDS("treatment_group_outcome_ls_ABADIE.Rds")
# The annoying missing observation
ls_treated_outcome[4,36] <- 1950

# Create key variables
nn_1 <- ncol(ls_treated_outcome)-1
train_per_end <- 7
test_per_start <- train_per_end+1
test_per_end <- train_per_end+2
post_tr_per_start <- train_per_end+3
TT <- nrow(ls_treated_outcome)

# Treatment Group: Create testing and training data
ls_treated_outcome_train <- ls_treated_outcome[c(1:train_per_end),]
ls_treated_outcome_train <- as.matrix(c(as.matrix(ls_treated_outcome_train[,-1])))
ls_treated_outcome_test <- ls_treated_outcome[c(test_per_start:test_per_end),]
ls_treated_outcome_test <- as.matrix(c(as.matrix(ls_treated_outcome_test[,-1])))
ls_treated_outcome_post_tr <- ls_treated_outcome[c(post_tr_per_start:TT),]
ls_treated_outcome_post_tr <- as.matrix(c(as.matrix(ls_treated_outcome_post_tr[,-1])))
################################################################################



##################  1.2. DONOR POOL ###########################################
# Donors Group: Outcome Variables
ls_control_outcome_single <- readRDS("control_group_outcome_ls.Rds")
# Create the full set of interactions
ls_control_outcome_single <- ls_control_outcome_single[,-1]
ls_control_outcome_single_train <- ls_control_outcome_single[c(1:train_per_end), ]
ls_control_outcome_single_test <- ls_control_outcome_single[c(test_per_start :test_per_end), ]
ls_control_outcome_post_tr <- ls_control_outcome_single[c(post_tr_per_start :TT), ]

# Create supplementary variables
nn_0 <- ncol(ls_control_outcome_single )
l_N <- matrix(1, ncol = 1, nrow = nn_1 )
l_K <- matrix(1, ncol = 1, nrow = nn_0)

# Add a donor pool for every observation
ls_control_outcome_full<- do.call(rbind,
                                  replicate(nn_1, ls_control_outcome_single, 
                                                  simplify=FALSE))
ls_control_outcome_full_train <- do.call(rbind, 
                                  replicate(nn_1, ls_control_outcome_single[1:train_per_end, ], 
                                                          simplify=FALSE))
ls_control_outcome_full_test <- do.call(rbind, 
                                    replicate(nn_1, ls_control_outcome_single[test_per_start:test_per_end,], 
                                                         simplify=FALSE))
ls_control_outcome_full_post_tr <- do.call(rbind, 
                                        replicate(nn_1, ls_control_outcome_single[post_tr_per_start :TT,], 
                                                  simplify=FALSE))
################################################################################



##################  1.3. COVARIATES ###########################################
# We define the time-invariant matrix and note that we only want continous
# predictors for the moment
covariates_treated <- readRDS("treatment_group_covariares_ABADIE.Rds")


# Create dummies for head_occupation
covariates_treated <- as.data.frame(t(covariates_treated ))
occupations_dictionary <- cbind(as.character(unique(covariates_treated$head_occupation)), c(1:9))
industry_dictionary <- cbind(as.character(unique(covariates_treated$head_industry)), c(1:10))

# Create Dummy for Occupation
j<- 1
for(t in unique(covariates_treated$head_occupation)) {
  covariates_treated[paste("occupation",j,sep="_")] <- ifelse(covariates_treated$head_occupation==t,1,0)
  j <- j+1
}

# Create Dummy for Industry
j<- 1
for(t in unique(covariates_treated$head_industry)) {
  covariates_treated[paste("industry",j,sep="_")] <- ifelse(covariates_treated$head_industry==t,1,0)
  j <- j+1
}
covariates_treated <- as.data.frame(t(covariates_treated ))

# Create a single data frame 
covariates_treated_continous_single <- covariates_treated[
  c(  "married",  "race_not_white",  "share_ill", "weight", "education", "age"
      #"occupation_1", "occupation_2", "occupation_3", "occupation_4", "occupation_5", "occupation_6", "occupation_7", "occupation_8",  "industry_2", 
      # "industry_3", "industry_4", "industry_5", "industry_6", "industry_7", "industry_8", "industry_9", "industry_10"
      ),]

# turn character variables into numeric once
covariates_treated_continous_single <- round(
  covariates_treated_continous_single  %>%
    mutate_all(funs(as.numeric(as.character(.) ))),2)
QQ_cont <- nrow(covariates_treated_continous_single )

# I think no neeed for further coding here
################################################################################






####################### 2. ESTIMATION  #########################################
####################### 2.1. Solve the problem #################################
application_run_rwscm<- function(num_donors, num_treated, 
                                 outcomes_matrix, covariates, pre_treatment_length){
  TT <- ncol(outcomes_matrix)-1
  post_treatment_period <- TT - pre_treatment_length
  
  # Create matrices
  X_pre <- outcomes_matrix[outcomes_matrix$treatment_indicator < 1, 
                           c(1:(ncol(outcomes_matrix) - 1 - post_treatment_period ))]
  Y_pre <- outcomes_matrix[outcomes_matrix$treatment_indicator > 0, 
                           c(1:(ncol(outcomes_matrix) - 1 - post_treatment_period ))]
  Y_pre_long <- as.matrix(c(t(as.matrix(Y_pre))))
  
  # Outcome matrix for control group
  l_N <- matrix(1, ncol = 1, nrow = num_treated )
  l_K <- matrix(1, ncol = 1, nrow = num_donors)
  X_pre_full <- do.call(rbind,replicate(num_treated, t(X_pre),  simplify=FALSE))
  
  # Set up the (Q x nnum_treated) matrix of covariates
  QQ_cont <- ncol(covariates ) 
  H_covariates_treatment_sample <- as.data.frame(t(covariates))
  
  # Set Up the Constraints
  wei <- Variable(rows = num_donors)
  alphas <- Variable(rows = num_donors, cols = QQ_cont)
  mu <- Variable(rows = num_treated)
  
  # Intercepts
  O_NK <- matrix(0, ncol = num_treated, nrow = num_donors)
  l_KQ <- matrix(0, ncol = 1, nrow = num_donors * QQ_cont)
  
  O1_n1Tn1 <-  matrix(0, ncol = num_treated, nrow = num_treated* pre_treatment_length)
  for( i in 0: (ncol(O1_n1Tn1)-1) ){
    O1_n1Tn1[ c((i*pre_treatment_length+1):(i*pre_treatment_length+pre_treatment_length)) , i+1] <- 1
  }
  
  # PROBLEM IS THE WEIRD CONSTRAINT
  print(dim(Y_pre_long))
  print(dim(O1_n1Tn1 %*% mu ))
  objective <- Minimize(sum_squares(
    as.matrix(Y_pre_long) -
      O1_n1Tn1 %*% mu -
      X_pre_full %*% wei -
      vec( t(X_pre) %*% alphas %*% H_covariates_treatment_sample ) ) )
  constraints_int <- list(
    wei %*% t(l_N) + alphas %*%  H_covariates_treatment_sample >= O_NK,
    t(l_K) %*% (wei %*% t(l_N) + alphas %*% H_covariates_treatment_sample) == t(l_N) )
  print(paste0("Random Weights Synthetic Control is running"))
  problem <- Problem(objective, constraints_int)
  solution  <- psolve(problem, verbose = T, num_iter = 500000)
  # Next step is to calculate the weights for each observation
  if(solution$status == "solver_error"){
    print("upsy small fail equal weights sad times")
    equal_weights <- matrix(1/num_donors, ncol = num_treated, nrow = num_donors)
    return(equal_weights)
  }
  else{
    alpha_QK <- solution$getValue(alphas)
    cons_weights <- solution$getValue(wei)
    intercepts <- solution$getValue(mu)
    
    # Calculate each weight
    ind_spec <- as.matrix(alpha_QK  %*% as.matrix(H_covariates_treatment_sample))
    ind_inv <- do.call(cbind, replicate(num_treated, cons_weights, simplify=FALSE))
    random_weights <- ind_inv + ind_spec
    newList <- list("random_weights" = random_weights, "intercepts" = intercepts)
  }
}
























################################################################################
# Prepare 
wei <- Variable(rows = nn_0)
alphas <- Variable(rows = nn_0, cols = QQ_cont)
O_n1n0 <- matrix(0, ncol = nn_1, nrow = nn_0)
l_KQ <- matrix(0, ncol = 1, nrow = nn_0 * QQ_cont)

# Form the objective function
objective <- Minimize(sum_squares(
  as.matrix(ls_treated_outcome_train) - 
    as.matrix(ls_control_outcome_full_train) %*% wei - 
    vec( ls_control_outcome_single_train %*% alphas %*% covariates_treated_continous_single) ) )

# Form the constraints 
constraints_int <- list(
  wei %*% t(l_N) + alphas %*% covariates_treated_continous_single >= O_n1n0 , 
  t(l_K) %*% (wei %*% t(l_N) + alphas %*% covariates_treated_continous_single) == t(l_N) )
# Solve the problem 
problem <- Problem(objective, constraints_int)
solution  <- psolve(problem, verbose = T, num_iter = 100000)

# Get the reduced-form parameters
cons_weights <- as.matrix(solution$getValue(wei))
alpha_QK <- solution$getValue(alphas)
################################################################################


######################### 2.2. Calculate weights ###############################
# Calculate each weight
random_weights <- matrix(0, nrow = 45, ncol = 1502)
ind_spec <- as.matrix(alpha_QK  %*% as.matrix(covariates_treated_continous_single))
ind_inv <- do.call(cbind, replicate(nn_1, cons_weights, simplify=FALSE))
random_weights <- ind_inv+ ind_spec
################################################################################


######################### 2.3. Make prediction #################################
calculate_rmse <-  function(treated_test, control_test , weights ){
  treated_units <- ncol(treated_test)
  testing_period <- nrow(treated_test)
  predicted <- as.matrix(control_test)  %*% as.matrix(weights)
  results <- matrix(ncol = treated_units + 1, nrow = testing_period + 1)
  results[c(1:testing_period), c(1:treated_units)] <- (as.matrix(predicted) -
                                                         as.matrix(treated_test))^2
  results[testing_period+1, c(1:treated_units)] <-
    sqrt(((results[1, c(1:treated_units)] +
             results[2, c(1:treated_units)] + results[3, c(1:treated_units)]))/3 )
   results[c(1:testing_period), treated_units+1] <- sqrt((results[c(1:testing_period), c(1:treated_units)] %*%
                                                           matrix(1, ncol = 1, nrow = treated_units))/treated_units)
  results[testing_period+1, treated_units+1] <- rmse(c(predicted),
                                                     c(as.matrix(treated_test)))
  return(results)
}
a <- as.matrix(ls_control_outcome_single_test)  %*% as.matrix(random_weights)

dim(ls_treated_outcome[c(6:8), -1])
dim(ls_control_outcome_single_test)
ls_treated_outcome_test <- ls_treated_outcome[c(6:8), -1]
rmse_rwscm_cont <- calculate_rmse(treated_test = ls_treated_outcome[c(6:8), -1], 
               control_test = ls_control_outcome_single_test , weights = random_weights)


################################################################################









######################### 3. DISCRETE CASE  ####################################

####################### 3.1. Prepare Data ######################################
# Focus on the case of only discrete predictors
a <- as.numeric(t(covariates_treated["share_ill",]))
covariates_treated_discrete_single <- covariates_treated[
  c(  "married",  "race_not_white",  "share_ill",
      "occupation_1", "occupation_2", "occupation_3", "occupation_4", 
      "occupation_5", "occupation_6", "occupation_7", "occupation_8", 
      "industry_2", "industry_3", "industry_4", "industry_5", "industry_6", 
      "industry_7", "industry_8", "industry_9", "industry_10"
  ),]
# turn character variables into numeric once
covariates_treated_discrete_single <- round(
  covariates_treated_discrete_single %>%
    mutate_all(funs(as.numeric(as.character(.) ))),2)
treated_names <- colnames(covariates_treated_discrete_single)
covariates_treated_discrete_single <- t(covariates_treated_discrete_single) 
covariates_treated_discrete_single[,3] <- ifelse(covariates_treated_discrete_single[,3] > summary(a )[5], 1, 0 )
covariates_treated_discrete_single <- t(covariates_treated_discrete_single) 
QQ_disc <- nrow(covariates_treated_discrete_single)

################################################################################


####################### 3.2. Solve the problem #################################
# Prepare 
wei <- Variable(rows = nn_0)
alphas <- Variable(rows = nn_0, cols = QQ_disc)
O_n1n0 <- matrix(0, ncol = nn_1, nrow = nn_0)
l_KQ <- matrix(0, ncol = 1, nrow = nn_0 * QQ_disc)

# Form the objective function
objective_disc <- Minimize(sum_squares(
  as.matrix(ls_treated_outcome_train) - 
    as.matrix(ls_control_outcome_full_train) %*% wei - 
    vec( ls_control_outcome_single_train %*% alphas %*% covariates_treated_discrete_single) ) )

# Form the constraints 
constraints_int <- list(
  wei %*% t(l_N) + alphas %*% covariates_treated_discrete_single >= O_n1n0 , 
  t(l_K) %*% (wei %*% t(l_N) + alphas %*% covariates_treated_discrete_single) == t(l_N) )
# Solve the problem 
problem <- Problem(objective_disc, constraints_int)
solution  <- psolve(problem, verbose = T, num_iter = 100000)

# Get the reduced-form parameters
cons_weights <- as.matrix(solution$getValue(wei))
alpha_QK <- solution$getValue(alphas)
################################################################################


######################### 3.3. Calculate weights ###############################
# Calculate each weight
random_weights <- matrix(0, nrow = 45, ncol = 1502)
ind_spec <- as.matrix(alpha_QK  %*% as.matrix(covariates_treated_discrete_single))
ind_inv <- do.call(cbind, replicate(nn_1, cons_weights, simplify=FALSE))
random_weights <- ind_inv+ ind_spec
################################################################################


######################### 3.4. Make prediction #################################
a <- as.matrix(ls_control_outcome_single_test)  %*% as.matrix(random_weights)
ls_treated_outcome_test <- ls_treated_outcome[c(6:8), -1]
rmse_rwscm_discr <- calculate_rmse(treated_test = ls_treated_outcome[c(6:8), -1], 
                             control_test = ls_control_outcome_single_test , weights = random_weights)
################################################################################




