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
#     4) Get ATT for the post-treatment period;

# Change definition of time periods
wage_train_per_end <- 4
wage_test_per_start <- wage_train_per_end+1
wage_test_per_end <- wage_train_per_end+2
wage_post_tr_per_start <- wage_train_per_end+3


# We load the datasets
wage_control_outcomes <- readRDS("control_group_outcome_wage.Rds")
wage_control_outcomes <- wage_control_outcomes[-c(1:2),-1]
TT <- nrow(wage_control_outcomes)
true_years_wage <- cbind(c(1:TT), wage_control_outcomes[,1])
wage_treated_outcomes <- readRDS("treatment_group_outcome_wage_matrix.Rds")
wage_treated_outcomes <- wage_treated_outcomes[-c(1:2),-1]


# Construct outcomes_matrix_wage
outcomes_matrix_wage <- as.data.frame(rbind(t(wage_control_outcomes) , t(wage_treated_outcomes)))
# Construct outcomes_matrix_wage
log_wage_control_outcomes <- log(wage_control_outcomes +1)
log_wage_treated_outcomes <- log(wage_treated_outcomes+1)
outcomes_matrix_wage$treatment_indicator <- 1
outcomes_matrix_wage$treatment_indicator[1:ncol(log_wage_control_outcomes)] <- 0
# Create a function for calculating the weights 
log_wage_treated_outcomes_post_tr <- log_wage_treated_outcomes [c(wage_post_tr_per_start:TT), ]

# Turn into weights
log_outcomes_matrix_wage<- log(outcomes_matrix_wage[,c(1:wage_test_per_end)]+1)
log_outcomes_matrix_wage <- cbind(log_outcomes_matrix_wage,
                                  outcomes_matrix_wage[,ncol(outcomes_matrix_wage)] )
log_outcomes_matrix_wage <- as.data.frame(log_outcomes_matrix_wage)
colnames(log_outcomes_matrix_wage)[ncol(log_outcomes_matrix_wage )] <- "treatment_indicator"
# Note: we shall use treatment_covariates_single from the previous coding
rwscm_res_wages <- application_run_rwscm(num_donors = nn_0, 
                                          num_treated = nn_1 , 
                                          outcomes_matrix = log_outcomes_matrix_wage, 
                                          covariates = treatment_covariates_single, 
                                          pre_treatment_length = wage_train_per_end)

random_weights_wages <- rwscm_res_wages[['random_weights']]
intercepts_wages  <- rwscm_res_wages[['intercepts']]
ones_matrix_wages <- matrix(1, ncol = TT, nrow = sum(log_outcomes_matrix_wage$treatment_indicator) )
all_intercepts_wages <- do.call(cbind, replicate(TT, intercepts_wages, simplify=FALSE))

# Calculate the full set of synthetic controls
synthetic_controls_wages <-  t(t(as.matrix(log_wage_control_outcomes ) %*% random_weights_wages) +  all_intercepts_wages)

################################################################################





########################### 3. USEFUL STATS  ###################################
########################### 3. 1. RMSE    ######################################

# Get rmse
rmse_training_period_wages <- matrix(0, ncol = 3, nrow = 3)
colnames(rmse_training_period_wages) <- c("Year", "CSC", "PSC")
rmse_training_period_wages[,1] <- c("1978", "1979", "Total")

#  Save RMSE results
rmse_training_period_wages[1,2] <- round(rmse(synthetic_controls_wages[wage_test_per_start,], 
                                           as.numeric(log_wage_treated_outcomes[wage_test_per_start,]) ),2 )
rmse_training_period_wages[2,2] <- round(rmse(synthetic_controls_wages[wage_test_per_end,], 
                                           as.numeric(log_wage_treated_outcomes[wage_test_per_end,]) ),2 )
rmse_training_period_wages[3,2] <- round(rmse(synthetic_controls_wages[c(wage_test_per_start:wage_test_per_end),], 
                                           as.matrix(log_wage_treated_outcomes[c(wage_test_per_start:wage_test_per_end),]) ),2 )
write.csv(rmse_training_period_wages, "wage_rmse_training_per.csv")
################################################################################


########################### 3. 2. ATT    ######################################
# Calculated ATT-s
csc_att_wage <-  matrix(0, ncol = 6, nrow = nn_1)
colnames(csc_att_wage) <- c("mr_num", "y80", "y81", "y82", "y83", "y84")

# Calculate treatment effects
csc_att_wage[,1] <- t(colnames(synthetic_controls_wages))
csc_att_wage[,c(2:6)] <-round(t(as.matrix(log_wage_treated_outcomes_post_tr) - synthetic_controls_wages[c((wage_test_per_end+1):TT), ]) ,2)
csc_att_wage <- as.data.frame(csc_att_wage)

# Merge data on education
education_data <- as.data.frame(cbind(rownames(treatment_covariates_single), treatment_covariates_single$some_college))
csc_att_wage <- csc_att_wage %>% left_join(education_data, by = c("mr_num" = "V1"))

# Save results
write.csv(csc_att_wage, "wage_csc_att_4.csv")
################################################################################





