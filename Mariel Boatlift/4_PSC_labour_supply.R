setwd("~/Documents/Oxford Yr 2/THESIS_micro_synth/Mariel boatlift shock/Analysis Data")
# Load packages 
library(CVXR)
# install.packages("Metrics")
library(Metrics)
library(dplyr)

# The purpose of this file is to:
#     1) Prepare the data on wage and divide it into: 
#         i) training+testing period (1972-1980);
#         ii) post-treatment period (1981-1984);
#     2) Run CSC on the data for the training period;
#     3) Compute synthetic controls for the testing period and get rmse;

#################################################################################



#####################   1.  PREPARE DATA  ######################################
##################  1.1. TREATMENT GROUP #######################################
# Outcomes matrix FOR TRAINING PERIOD
ls_treated_outcomes <- readRDS("treatment_group_outcome_ls_ABADIE.Rds")
# The annoying missing observation
ls_treated_outcomes[4,36] <- 1950
ls_treated_outcomes <- ls_treated_outcomes[-c(1:6),]
true_years_ls <- cbind(c(1:nrow(ls_treated_outcomes)) , ls_treated_outcomes[,1] - 1)
ls_treated_outcomes <- ls_treated_outcomes[,-1]
# Create key variables
nn_1 <- ncol(ls_treated_outcomes)
ls_train_per_end <- 1
ls_test_per_start <- ls_train_per_end+1
ls_test_per_end <- ls_train_per_end+2
ls_post_tr_per_start <- ls_train_per_end+3
TT <- nrow(ls_treated_outcomes)

# Treatment Group: Create testing and training data
ls_treated_outcomes_pre_tr <- ls_treated_outcomes[c(1:ls_test_per_end),]
# ls_treated_outcomes_train <- ls_treated_outcomes[c(1:train_per_end),]
# ls_treated_outcomes_test <- ls_treated_outcomes[c(test_per_start:test_per_end),]
ls_treated_outcomes_post_tr <- ls_treated_outcomes[c(ls_post_tr_per_start:TT),]

##################  1.2. CONTROL GROUP #######################################
# Control group
ls_control_outcomes_single <- readRDS("control_group_outcome_ls.Rds")
ls_control_outcomes_single <- ls_control_outcomes_single[-c(1:6),-1]

# Create the full set of interactions
nn_0 <- ncol(ls_control_outcomes_single )
ls_control_outcomes_pre_tr <- ls_control_outcomes_single[c(1:ls_test_per_end), ]
ls_control_outcomes_post_tr <- ls_control_outcomes_single[c(ls_post_tr_per_start :TT), ]

# Combine into outcome matrix
outcomes_matrix_ls <- t(as.matrix(cbind(ls_control_outcomes_pre_tr, ls_treated_outcomes_pre_tr )))
outcomes_matrix_ls <- cbind(outcomes_matrix_ls, c(rep(0, nn_0), rep(1, nn_1)) )
################################################################################



################# 1.3. COVARIATES   ############################################
# We define the time-invariant matrix and note that we only want continous
# predictors for the moment
treatment_covariates <- readRDS("treatment_group_covariares_matrix.Rds")
treatment_covariates <- treatment_covariates[c(
  "father_occupation", "head_occupation", "head_industry", 
  "weight", "education", "age", "race_not_white", "married", "share_ill"),]
control_covariates <- readRDS("control_group_covariates.Rds")
control_covariates<- control_covariates[c(
  "father_occupation", "head_occupation", "head_industry", 
  "weight", "education", "age", "race_not_white", "married", "share_ill"),]
all_covariates<-as.data.frame(t(cbind(control_covariates, treatment_covariates )))
all_covariates$treatment_indicator <- c(rep(0, nn_0), rep(1, nn_1)) 

# Contruct covariates matrix
occupations_industry_per_person <- all_covariates[,c(2,3)]
occupations_dictionary <- as.data.frame(cbind(as.character(unique(all_covariates$head_occupation)), 
                                              c(1:length(levels(all_covariates$head_occupation))) ))
colnames(occupations_dictionary) <- c("occupation", "code_occ") 
occupations_dictionary$low_skilled_occ <- 1
occupations_dictionary$low_skilled_occ[c(1,3)] <- 0
all_covariates$head_industry <- as.character(all_covariates$head_industry) 
all_covariates$head_industry <- ifelse(is.na(all_covariates$head_industry), "Manufacturing" , all_covariates$head_industry)
all_covariates$head_industry <- as.factor(all_covariates$head_industry )
industry_dictionary <- as.data.frame(cbind(as.character(unique(all_covariates$head_industry)), 
                                           c(1:length(levels(all_covariates$head_industry)) ) ))
colnames(industry_dictionary) <- c("industry", "code_ind") 
industry_dictionary$low_skilled_ind <- 1
industry_dictionary$low_skilled_ind[c(1,4)] <- 0

# Match to original data
occupations_industry_per_person  <- occupations_industry_per_person  %>% 
  left_join(occupations_dictionary, by = c("head_occupation" = "occupation")) %>%
  left_join(industry_dictionary, by = c("head_industry" = "industry"))
occupations_industry_per_person$low_skilled <- ifelse(
  (occupations_industry_per_person$low_skilled_occ == 1) & (occupations_industry_per_person$low_skilled_ind == 1), 1,0 )

# Create Dummy for Occupation
j<- 1
for(t in unique(all_covariates$head_occupation)) {
  all_covariates[paste("occupation",j,sep="_")] <- ifelse(all_covariates$head_occupation==t,1,0)
  j <- j+1
}

# Create Dummy for Industry
j<- 1
for(t in unique(all_covariates$head_industry)) {
  all_covariates[paste("industry",j,sep="_")] <- ifelse(all_covariates$head_industry==t,1,0)
  j <- j+1
}
# Education
all_covariates$education <- as.numeric(as.character(all_covariates$education))
table(all_covariates$education)
# Note for one person we DON'T know education, as it is coded as 99
all_covariates$some_college <- ifelse(all_covariates$education > 12, 1, 0)

# Ill
all_covariates$share_ill <- as.numeric(as.character(all_covariates$share_ill))
summary(all_covariates$share_ill)
all_covariates$was_often_ill <- ifelse(all_covariates$share_ill > 0.75, 1, 0)

# Create a single data frame 
all_covariates_single <- all_covariates[, c(  "some_college", "married", "race_not_white", "was_often_ill",
                                              "occupation_1", "occupation_2", "occupation_3", "occupation_4", 
                                              "occupation_5", "occupation_6", "occupation_7", "occupation_8",  
                                              "industry_2",  "industry_3", "industry_4", "industry_5", 
                                              "industry_6", "industry_7", "industry_8", "industry_9", "industry_10",
                                              "treatment_indicator"    )]

# turn character variables into numeric once
all_covariates_single <- round(
  all_covariates_single  %>%
    mutate_all(funs(as.numeric(as.character(.) ))),2)
################################################################################






#####################   2. ESTIMATION  #########################################
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
                     covariates_continous, train_period_end){
  
  # 0. Preparate the data
  all_T <- ncol(outcome_variable_matrix) - 1
  QQ <- ncol(covariates_continous)
  outcome_variable_matrix <- as.data.frame(outcome_variable_matrix )
  colnames(outcome_variable_matrix)[ncol(outcome_variable_matrix )] <- "treatment_indicator"

  # Create useful matrices
  donors_pre_treatment <- outcome_variable_matrix[
    outcome_variable_matrix$treatment_indicator == 0,c(1:all_T)]
  treated_pre_treatment <- outcome_variable_matrix[
    outcome_variable_matrix$treatment_indicator == 1, c(1:all_T)]

  stack_X_pre_train <-cbind(  donors_pre_treatment[, c(1:( train_period_end ))], 
                                 covariates_continous[covariates_continous$treatment_indicator == 0, c(1:QQ) ])
  stack_Y_pre_train <- cbind(   treated_pre_treatment[, c(1:( train_period_end ))], 
                                  covariates_continous[covariates_continous$treatment_indicator == 1,c(1:QQ) ])
  stack_X_pre_test <- cbind(  donors_pre_treatment[, c( (train_period_end+1):all_T )], 
                                covariates_continous[covariates_continous$treatment_indicator == 0,c(1:QQ) ])
  stack_Y_pre_test<- cbind(   treated_pre_treatment[, ((train_period_end+1):all_T  )], 
                                covariates_continous[covariates_continous$treatment_indicator == 1,c(1:QQ) ])
  
  ####### TRRANPOSE!!!
  stack_X_pre_train <- t(stack_X_pre_train)
  stack_Y_pre_train <- t(stack_Y_pre_train )
  stack_X_pre_test <- t( stack_X_pre_test)
  stack_Y_pre_test <- t(stack_Y_pre_test)
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
  return( pscm_weights )
}

penalised_weights_ls <- run_pscm( num_donors = nn_0, num_treated = nn_1, 
                         outcome_variable_matrix = outcomes_matrix_ls ,
                         covariates_continous =  all_covariates_single,
                         train_period_end = ls_train_per_end)

pen_synthetic_controls_ls <-  as.matrix(ls_control_outcomes_single) %*% penalised_weights_ls
################################################################################





############## 3. SAVE RESULTS    ##############################################
rmse_training_period_ls <- read.csv("ls_rmse_training_per.csv")
rmse_training_period_ls <- rmse_training_period_ls[,-1]
rmse_training_period_ls[1,3] <- round(rmse(pen_synthetic_controls_ls[ls_test_per_start,], 
                                           as.numeric(ls_treated_outcomes[ls_test_per_start,]) ),2 )
rmse_training_period_ls[2,3] <- round(rmse(pen_synthetic_controls_ls[ls_test_per_end,], 
                                           as.numeric(ls_treated_outcomes[ls_test_per_end,]) ),2 )
rmse_training_period_ls[3,3] <- round(rmse(pen_synthetic_controls_ls[c(ls_test_per_start:ls_test_per_end),], 
                                           as.matrix(ls_treated_outcomes[c(ls_test_per_start:ls_test_per_end),]) ),2 )

