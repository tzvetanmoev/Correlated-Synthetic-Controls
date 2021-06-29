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
# Change definition of time periods
wage_train_per_end <- 2
wage_test_per_start <- wage_train_per_end+1
wage_test_per_end <- wage_train_per_end+2
wage_post_tr_per_start <- wage_train_per_end+3


# We load the datasets
wage_control_outcomes <- readRDS("control_group_outcome_wage.Rds")
wage_control_outcomes <- wage_control_outcomes[-c(1:4),-1]
TT <- nrow(wage_control_outcomes)
true_years_wage <- cbind(c(1:TT), wage_control_outcomes[,1])
wage_treated_outcomes <- readRDS("treatment_group_outcome_wage_matrix.Rds")
wage_treated_outcomes <- wage_treated_outcomes[-c(1:4),-1]


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
################################################################################




#####################   2. ESTIMATION  #########################################
penalised_weights_wages <- run_pscm( num_donors = nn_0, num_treated = nn_1, 
                               outcome_variable_matrix = log_outcomes_matrix_wage,
                               covariates_continous =  all_covariates_single,
                               train_period_end = wage_train_per_end)

pen_synthetic_controls_wages <-  as.matrix(log_wage_control_outcomes) %*% penalised_weights_wages
colnames(pen_synthetic_controls_wages) <- colnames(log_wage_treated_outcomes)
#### SAVE THE RESULTS
rmse_training_period_wages <- read.csv( "wage_rmse_training_per.csv")
rmse_training_period_wages <- rmse_training_period_wages[,-1]
rmse_training_period_wages[1,3] <- round(rmse(pen_synthetic_controls_wages[wage_test_per_start,], 
                                           as.numeric(log_wage_treated_outcomes[wage_test_per_start,]) ),2 )
rmse_training_period_wages[2,3] <- round(rmse(pen_synthetic_controls_wages [wage_test_per_end,], 
                                           as.numeric(log_wage_treated_outcomes[wage_test_per_end,]) ),2 )
rmse_training_period_wages[3,3] <- round(rmse(pen_synthetic_controls_wages [c(wage_test_per_start:wage_test_per_end),], 
                                           as.matrix(log_wage_treated_outcomes[c(wage_test_per_start:wage_test_per_end),]) ),2 )

write.csv(rmse_training_period_wages, "rmse_training_period_wages_1.csv" )
write.csv(rmse_training_period_ls, "rmse_training_period_ls_1.csv" )
