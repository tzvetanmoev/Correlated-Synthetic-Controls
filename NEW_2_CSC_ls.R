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



##################  1.3. COVARIATES ###########################################
# We define the time-invariant matrix and note that we only want continous
# predictors for the moment
treatment_covariates <- readRDS("treatment_group_covariares_matrix.Rds")
# Contruct covariates matrix
treatment_covariates<- as.data.frame(t(treatment_covariates))
occupations_industry_per_person <- treatment_covariates[,c(3,4)]
occupations_dictionary <- as.data.frame(cbind(as.character(unique(treatment_covariates$head_occupation)), c(1:9)))
colnames(occupations_dictionary) <- c("occupation", "code_occ") 
occupations_dictionary$low_skilled_occ <- 1
occupations_dictionary$low_skilled_occ[c(1,3)] <- 0
industry_dictionary <- as.data.frame(cbind(as.character(unique(treatment_covariates$head_industry)), c(1:10)))
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
for(t in unique(treatment_covariates$head_occupation)) {
  treatment_covariates[paste("occupation",j,sep="_")] <- ifelse(treatment_covariates$head_occupation==t,1,0)
  j <- j+1
}

# Create Dummy for Industry
j<- 1
for(t in unique(treatment_covariates$head_industry)) {
  treatment_covariates[paste("industry",j,sep="_")] <- ifelse(treatment_covariates$head_industry==t,1,0)
  j <- j+1
}
# Education
treatment_covariates$education <- as.numeric(as.character(treatment_covariates$education))
table(treatment_covariates$education)
# Note for one person we DON'T know education, as it is coded as 99
treatment_covariates$some_college <- ifelse(treatment_covariates$education > 12, 1, 0)

# Ill
treatment_covariates$share_ill <- as.numeric(as.character(treatment_covariates$share_ill))
summary(treatment_covariates$share_ill)
treatment_covariates$was_often_ill <- ifelse(treatment_covariates$share_ill > 0.75, 1, 0)

# Create a single data frame 
treatment_covariates_single <- treatment_covariates[,
                                                    c(   "some_college", "married", "race_not_white", "was_often_ill",
                                                         "occupation_1", "occupation_2", "occupation_3", "occupation_4", "occupation_5", "occupation_6", "occupation_7", "occupation_8",  "industry_2", 
                                                         "industry_3", "industry_4", "industry_5", "industry_6", "industry_7", "industry_8", "industry_9", "industry_10"
                                                    )]
save_cov <- treatment_covariates_single
# turn character variables into numeric once
treatment_covariates_single <- round(
  treatment_covariates_single  %>%
    mutate_all(funs(as.numeric(as.character(.) ))),2)

# Last necessary manipulations 
outcomes_matrix_ls <- as.data.frame(outcomes_matrix_ls)
colnames(outcomes_matrix_ls)[ncol(outcomes_matrix_ls)] <- "treatment_indicator"
# I think no neeed for further coding here
################################################################################






####################### 2. ESTIMATION  #########################################
# We create a function to sovle the problem for us....
application_run_rwscm<- function(num_donors, num_treated, 
                                 outcomes_matrix, covariates, pre_treatment_length){
  all_T <- ncol(outcomes_matrix)-1
  post_treatment_period <-  all_T - pre_treatment_length
  
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

# Calculate CSC (or its old name Random Weights Synthetic Control Method)
rwscm_res_ls <- application_run_rwscm(num_donors = nn_0, 
                      num_treated = nn_1, 
                      outcomes_matrix = outcomes_matrix_ls, 
                      covariates = treatment_covariates_single, 
                      pre_treatment_length = ls_train_per_end )
random_weights <- rwscm_res_ls[['random_weights']]
intercepts <- rwscm_res_ls[['intercepts']]
ones_matrix <- matrix(1, ncol = TT, nrow = sum(outcomes_matrix_ls$treatment_indicator) )
all_intercepts <- do.call(cbind, replicate(TT, intercepts, simplify=FALSE))

# Calculate the full set of synthetic controls
synthetic_controls_ls <-  t(t(as.matrix(ls_control_outcomes_single) %*% random_weights) +  all_intercepts)

################################################################################






####################### 3. GET USEFUL STATS ####################################
########################### 3. 1. RMSE    ######################################

# Get rmse
rmse_training_period_ls <- matrix(0, ncol = 3, nrow = 3)
colnames(rmse_training_period_ls) <- c("Year", "CSC", "PSC")
rmse_training_period_ls[,1] <- c("1978", "1979", "Total")

#  Save RMSE results
rmse_training_period_ls[1,2] <- round(rmse(synthetic_controls_ls[ls_test_per_start,], 
                                    as.numeric(ls_treated_outcomes[ls_test_per_start,]) ),2 )
rmse_training_period_ls[2,2] <- round(rmse(synthetic_controls_ls[ls_test_per_end,], 
                                    as.numeric(ls_treated_outcomes[ls_test_per_end,]) ),2 )
rmse_training_period_ls[3,2] <- round(rmse(synthetic_controls_ls[c(ls_test_per_start:ls_test_per_end),], 
                                    as.matrix(ls_treated_outcomes[c(ls_test_per_start:ls_test_per_end),]) ),2 )
write.csv(rmse_training_period_ls, "ls_rmse_training_per.csv")
###############################################################################


########################### 3. 2. ATT    ######################################

# Calculated ATT-s
csc_att_ls <-  matrix(0, ncol = 5, nrow = nn_1)
colnames(csc_att_ls) <- c("mr_num", "y80", "y81", "y82", "y83")

# Calculate treatment effects
csc_att_ls[,1] <- t(colnames(synthetic_controls_ls))
csc_att_ls[,c(2:5)] <-round(t(as.matrix(ls_treated_outcomes_post_tr) - synthetic_controls_ls[c((ls_test_per_end+1):TT), ]),2)
csc_att_ls <- as.data.frame(csc_att_ls)

# Merge data on education
education_data <- as.data.frame(cbind(rownames(treatment_covariates_single), treatment_covariates_single$some_college))
csc_att_ls <- csc_att_ls %>% left_join(education_data, by = c("mr_num" = "V1"))

# Save results
write.csv(csc_att_ls, "ls_csc_att_1.csv")
###############################################################################














