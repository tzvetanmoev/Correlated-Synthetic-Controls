# Our experiment will estimate synthetic controls for the period 1972 to 1976
# It will then estimate the RMSE for the periods 1977-1979
# The Mariol Boatlift occurred in 1980, so possible treatment effects should
# not be influencing the time series
library(tidyr)
library(janitor)
# We begin by restricting our dataset to observations we will actually need
# Note we need to do a bit of imputation, as some values are zero and some NA
treatment_group <- readRDS("treatment_group.Rds")
treatment_group <- treatment_group[-c(784:nrow(treatment_group )),]
# 4. Dropping of observations
# Let us first see how many observations we will have under a balanced panel
# with the 8 periods
treatment_group <- treatment_group %>% 
  group_by(pid) %>% 
  mutate(min_year = min(year), max_year = max(year) ) %>%
  filter(min_year < 1973 & max_year > 1979 ) 
  
# So, we have a lot of missing records for two heads: 22002 and 501001
# We will drop them and impute for the three guys who have just one missing row
a <- as.data.frame(table(treatment_group$year, treatment_group$pid))
# Dropping 
treatment_group <- treatment_group %>% 
  filter(pid!="22002" & pid!="501001" & pid!="1658001" & 
           pid!= "2890170" & pid!= "1660001" & pid!= "1900001" & pid!= "2885001" & 
           pid!= "965021")

# 5. Imputation for the four treated observations
# Add empty rows for the missing observations
all_combos <- crossing(treatment_group$pid, treatment_group$year)
colnames(all_combos) <- c("pid", "year")
treatment_group <- all_combos %>% 
  left_join(treatment_group, by = c("pid" = "pid", "year" = "year"))
treatment_group$new_match <- ifelse(is.na(treatment_group$current_state), 1, 0)
table(treatment_group$new_match)

# Add info on time invariant covariates
# Can be done with a loop...
# First Dude
treatment_group[
  (treatment_group$new_match==1) & 
    (treatment_group$pid == "2266001"), c(3, 6:16)] <- 
  treatment_group[
    lag(treatment_group$new_match==1) & 
      (treatment_group$pid == "2266001"), c(3, 6:16)]
# Second Dude
treatment_group[
  (treatment_group$new_match==1) & 
    (treatment_group$pid == "1988001"), c(3:16)] <- 
  treatment_group[
    lead(treatment_group$new_match==1) & 
      (treatment_group$pid == "1988001"), c(3:16)]
# RERUN FOR the same dude
treatment_group[
  (treatment_group$new_match==1) & 
    (treatment_group$pid == "1988001"), c(3:16)] <- 
  treatment_group[
    lead(treatment_group$new_match==1) & 
      (treatment_group$pid == "1988001"), c(3:16)]
# Second Dude
treatment_group[
  (treatment_group$new_match==1) & 
    (treatment_group$pid == "2594001"), c(3:19)] <- 
  treatment_group[
    lead(treatment_group$new_match==1) & 
      (treatment_group$pid == "2594001"), c(3:19)]
# RERUN FOR THE SAME DUDE
treatment_group[
  (treatment_group$new_match==1) & 
    (treatment_group$pid == "2594001"), c(3:19)] <- 
  treatment_group[
    lead(treatment_group$new_match==1) & 
      (treatment_group$pid == "2594001"), c(3:19)]
# Third dude impute
treatment_group[
  (treatment_group$new_match==1) & 
    (treatment_group$pid == "2277003"), c(3:16)] <- 
  treatment_group[
    lead(treatment_group$new_match==1) & 
      (treatment_group$pid == "2277003"), c(3:16)]
# Fourth dude impute
treatment_group[
  (treatment_group$new_match==1) & 
    (treatment_group$pid == "22001"), c(3:16)] <- 
  treatment_group[
    lead(treatment_group$new_match==1) & 
      (treatment_group$pid == "22001"), c(3:16)]

# Impute missing hours
treatment_group$head_obs_avr_earnings <- ifelse(
  is.na(treatment_group$head_obs_avr_earnings) & treatment_group$pid==lead(treatment_group$pid) &
  treatment_group$pid==lag(treatment_group$pid), 
  0.5*lag(treatment_group$head_obs_avr_earnings) + 0.5*lead(treatment_group$head_obs_avr_earnings),
  treatment_group$head_obs_avr_earnings)
sum(is.na(treatment_group$head_obs_avr_earnings))
# Impute missing wages
treatment_group$worked_hours <- ifelse(
  is.na(treatment_group$worked_hours) & treatment_group$pid==lead(treatment_group$pid) &
    treatment_group$pid==lag(treatment_group$pid), 
  0.5*lag(treatment_group$worked_hours) + 0.5*lead(treatment_group$worked_hours),
  treatment_group$worked_hours) 
treatment_group$worked_hours <- ifelse(
  is.na(treatment_group$worked_hours) & treatment_group$pid==lag(treatment_group$pid) & treatment_group$pid!=lead(treatment_group$pid), 
  lag(treatment_group$worked_hours),
  treatment_group$worked_hours) 
sum(is.na(treatment_group$worked_hours))

# Remove unnecessary variables
treatment_group <- as.data.frame(treatment_group)
treatment_group <- treatment_group[,-c(17:19)]


# Step 1: Good - no missing values for hours BUT note we have 0 values
treatment_group$miss_worked_hours <- ifelse(
  is.na(treatment_group$worked_hours) | (treatment_group$worked_hours==0), 1, 0 )
table(treatment_group$miss_worked_hours )

# Step 2: Wage
colnames(treatment_group)[5] <- "wage_pr"
treatment_group$miss_wage_pr <- ifelse(
  is.na(treatment_group$wage_pr) , 1, 0 )
table(treatment_group$miss_wage_pr )

# Step 3: Impute wages
# 3.1. When we have surrounding wages, we can do:
treatment_group$wage_pr <- ifelse(treatment_group$wage_pr==0, 
                                  NA, treatment_group$wage_pr)
treatment_group$lag_lead_wage_pr <-  (lag(treatment_group$wage_pr) + 
                                      lead(treatment_group$wage_pr))/2
treatment_group$wage_pr[
  (treatment_group$miss_wage_pr==1) & 
  (treatment_group$pid==lead(treatment_group$pid) &
  (treatment_group$pid!="1064003") & (treatment_group$pid!="965021") )] <- 
treatment_group$lag_lead_wage_pr[
  (treatment_group$miss_wage_pr==1) & 
  (treatment_group$pid==lead(treatment_group$pid))&
  (treatment_group$pid!="1064003") & (treatment_group$pid!="965021")] 
# 23 left

# 3.2. For the very first group we can, do
treatment_group$wage_pr[
  (treatment_group$miss_wage_pr==1) & (treatment_group$pid=="1064003") ]<-7.705


# 3.3. For the last observation in each group, we take the previous one
treatment_group$wage_pr[
  is.na(treatment_group$wage_pr) & 
   (treatment_group$pid!=lead(treatment_group$pid))] <- 
  treatment_group$wage_pr[
    lead(is.na(treatment_group$wage_pr)) & 
    (treatment_group$pid!=lead(treatment_group$pid, 2))]
# NOTE: this code does not really work I think

# 3.4. For the rest we simply take the mean of the other observations per group
# This is the most straightforward way to do imputation in this case
treatment_group <- treatment_group %>% 
  group_by(pid) %>% 
  mutate(avr_wage_pr = mean(wage_pr, na.rm = TRUE))
treatment_group$wage_pr[ 
  (is.na(treatment_group$wage_pr)) | (treatment_group$wage_pr==0) ]<- 
treatment_group$avr_wage_pr[ 
  (is.na(treatment_group$wage_pr)) | (treatment_group$wage_pr==0)]
treatment_group$miss_wage_pr <- ifelse(
  is.na(treatment_group$wage_pr) , 1, 0 )
table(treatment_group$miss_wage_pr )
# Ok, we are fucking done here

saveRDS(treatment_group, "treatment_group_PARTIALLY_POOLED.Rds")

# Now the above works well for my estimator but not so much for Abadie et al or 
# Abadie and L'Hour. For this reason, we need to further transform the dataset

# 1. Covariates for treatmetn group for Abadie
# Note this step is necessary for the later analysis... I think that it should
# have been better done at an earlier stage  but whatever
# Recode variables as binary 
table(treatment_group$race)
treatment_group$race_not_white <- ifelse(treatment_group$race=="white", 0, 1)
# control_group <- control_group %>% select(-race)
table(treatment_group$race_not_white)

# Recode religion similarly
table(treatment_group$religion)
# Actually, we can just drop it, as less than 5% are religious in our sample...
# That's a bit weird tbh

# Recode variables as binary 
table(treatment_group$marr_status)
treatment_group$married <- ifelse(treatment_group$marr_status=="married", 1, 0)
table(treatment_group$married)

# Remove bad observations
treatment_group  <- subset(treatment_group,  
                  select = -c(miss_worked_hours, miss_wage_pr, lag_lead_wage_pr, state_grew_up))

# Transpose
wide_treatment_covariates <- as.data.frame(
  treatment_group %>% 
  group_by(pid) %>%
  slice(1))
wide_treatment_covariates <- as.data.frame(t(subset(wide_treatment_covariates,  
                        select = -c(current_state, worked_hours, wage_pr))))
wide_treatment_covariates <- row_to_names(wide_treatment_covariates, 
                                        row_number = 1, remove_row = TRUE)

# Well, let us ignore the grow_up state for the moment, as we are not using it
# either way....
saveRDS(wide_treatment_covariates, "treatment_group_covariares_matrix.Rds")

# 2. Treatment Group for Abadie
# 2.1. Wages
treatment_outcomes_wage <- pivot_wider(treatment_group, 
                                     id_cols = c("pid", "year"),
                                     names_from = c("pid" ),
                                     values_from = c("wage_pr"))
treatment_outcomes_wage<- treatment_outcomes_wage  %>% arrange(year) 
treatment_outcomes_wage <- as.data.frame(treatment_outcomes_wage)
saveRDS(treatment_outcomes_wage, "treatment_group_outcome_wage_matrix.Rds")

# 2.2. Labour Supply
treatment_outcomes_ls <- pivot_wider(treatment_group, 
                                       id_cols = c("pid", "year"),
                                       names_from = c("pid" ),
                                       values_from = c("worked_hours"))
treatment_outcomes_ls<- treatment_outcomes_ls  %>% arrange(year) 
treatment_outcomes_ls <- as.data.frame(treatment_outcomes_ls)
saveRDS(treatment_outcomes_ls, "treatment_group_outcome_ls_matrix.Rds")








