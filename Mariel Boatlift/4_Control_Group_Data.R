library(tidyr)
library(janitor)
################################################################################
# Control Group - not so izi pizi
d <- readRDS("full_data_for_application.Rds")

control_group <- d[d$current_state!="AAA_Florida",]
control_group <- control_group[ -c(28221:nrow(control_group)), ]

# Now given how big the donor pool is we will only get observations which
# have a nearly complete set of observation
control_group <- control_group %>% 
  group_by(pid) %>% 
  mutate(min_year = min(year), max_year = max(year) ) %>%
  filter(min_year < 1973 & max_year > 1983 )  


# We will drop them and impute for the three guys who have just one missing row
a <- as.data.frame(table(control_group$year, control_group$pid))
missing_years <- a %>% filter(Freq == 0) %>%  group_by(Var2) %>% slice(1) %>% select(Var2) 
# So, we have a lot of missing records for several heads: 
# 1043005; 13004; 1574001; 1579003; 1701002; 2198003; 2898003
control_group <- subset(control_group, !(pid %in% c(as.character(missing_years$Var2 ) )))


# No we proceed with the imputation exericise
# Step 1: Good - no missing values for hours BUT note we have 0 values
control_group$miss_worked_hours <- ifelse(
  is.na(control_group$worked_hours) | (control_group$worked_hours==0), 1, 0 )
table(control_group$miss_worked_hours )

# Step 2: Wage
colnames(control_group)[5] <- "wage_pr"
control_group$miss_wage_pr <- ifelse(
  is.na(control_group$wage_pr) | (control_group$wage_pr==0), 1, 0 )
table(control_group$miss_wage_pr )
# Quite a few missing values. Let's see if we have some observations that are
# particularly bad
# AND Drop observations with more than 2 missing values
missing_outcomes <- control_group %>% group_by(pid) %>% 
  summarise(missing = sum(miss_wage_pr )) %>% arrange(missing) %>%
  filter(missing > 2) %>% select(pid)
control_group = subset(control_group, !(pid %in% c(missing_outcomes$pid )))

# Step 3: Impute wages
# 3.1. When we have surrounding wages, we can do:
control_group$wage_pr <- ifelse(control_group$wage_pr==0, 
                                  NA, control_group$wage_pr)
control_group$lag_lead_wage_pr <-  (lag(control_group$wage_pr) + 
                                        lead(control_group$wage_pr))/2
control_group <- control_group %>% group_by(pid, year)
control_group$wage_pr[
  (control_group$miss_wage_pr==1) & 
    (control_group$pid==lead(control_group$pid) &
       (control_group$pid!="1173001") & (control_group$pid!="351001") )] <- 
  control_group$lag_lead_wage_pr[
    (control_group$miss_wage_pr==1) & 
      (control_group$pid==lead(control_group$pid))&
      (control_group$pid!="1173001") & (control_group$pid!="351001")] 
# 229 left

# 3.3. For the last observation in each group, we take the previous one
control_group$wage_pr[
  is.na(control_group$wage_pr) & 
    (control_group$pid!=lead(control_group$pid))] <- 
  control_group$wage_pr[
    lead(is.na(control_group$wage_pr)) & 
      (control_group$pid!=lead(control_group$pid, 2))]
# For the first observation we do it the other way around
control_group$wage_pr[
  is.na(control_group$wage_pr) & 
    (control_group$pid!=lag(control_group$pid))] <- 
  control_group$wage_pr[
    lag(is.na(control_group$wage_pr)) & 
      (control_group$pid!=lag(control_group$pid, 2))]

# 3.4. For the rest we simply take the mean of the other observations per group
# This is the most straightforward way to do imputation in this case
control_group<- control_group %>% 
  group_by(pid) %>% 
  mutate(avr_wage_pr = mean(wage_pr, na.rm = TRUE))
control_group$wage_pr[ 
  (is.na(control_group$wage_pr)) | (control_group$wage_pr==0) ]<- 
  control_group$avr_wage_pr[ 
    (is.na(control_group$wage_pr)) | (control_group$wage_pr==0)]

# 3.5. We can now impute hours for the cases where earnings are not missing
# by simply taking the average across groups 
control_group <- control_group %>% 
  group_by(pid) %>% filter(worked_hours>5) %>% 
  summarise(avr_hours = mean(worked_hours, na.rm = TRUE)) %>% 
  right_join(control_group, by = c("pid" = "pid")) 
control_group$worked_hours[(control_group$worked_hours==0) & 
      (!is.na(control_group$worked_hours)) & (!is.na(control_group$wage_pr))] <- 
control_group$avr_hours[(control_group$worked_hours==0) & 
       (!is.na(control_group$worked_hours))& (!is.na(control_group$wage_pr))]
# THIS Imputation step is definitely not ideal, however.
# Remove unnecessary variables
control_group <- control_group[, -c(2, 18:ncol(control_group))]

# Now turn this into wide format
# Wage
control_group <- control_group %>% arrange(pid, year)
control_outcomes_wage <- pivot_wider(control_group, 
                                       id_cols = c("pid", "year"),
                                       names_from = c("pid" ),
                                       values_from = c("wage_pr"))
control_outcomes_wage<- control_outcomes_wage  %>% arrange(year) 
control_outcomes_wage$na_count <- apply(control_outcomes_wage, 1, function(x) sum(is.na(x)))
print(control_outcomes_wage$na_count)
control_outcomes_wage <- control_outcomes_wage[, -ncol(control_outcomes_wage)]
control_outcomes_wage <- as.data.frame(control_outcomes_wage)
saveRDS(control_outcomes_wage, "control_group_outcome_wage.Rds")

# Labour Supply
control_group <- control_group %>% arrange(pid, year)
control_outcomes_ls <- pivot_wider(control_group, 
                                     id_cols = c("pid", "year"),
                                     names_from = c("pid" ),
                                     values_from = c("worked_hours"))
control_outcomes_ls<- control_outcomes_ls  %>% arrange(year) 
control_outcomes_ls$na_count <- apply(control_outcomes_ls, 1, function(x) sum(is.na(x)))
print(control_outcomes_ls$na_count)
control_outcomes_ls <- control_outcomes_ls[, -ncol(control_outcomes_ls)]
control_outcomes_ls <- as.data.frame(control_outcomes_ls)
saveRDS(control_outcomes_ls, "control_group_outcome_ls.Rds")
################################################################################






################################################################################
# Get time invariant covariates
control_group$state_grew_up[is.na(control_group$state_grew_up)] <- 
  control_group$current_state[is.na(control_group$state_grew_up)]
#control_group <- control_group %>% 
# mutate(has_industry = is.na(head_industry)) %>%
#  filter(has_industry == FALSE)
# Industry will not enter either way into Abadie's SCM, so we better keep it 
# to avoid some difficulties at a later stage

# Recode variables as binary 
table(control_group$race)
control_group$race_not_white <- ifelse(control_group$race=="white", 0, 1)
# control_group <- control_group %>% select(-race)
table(control_group$race_not_white)

# Recode religion similarly
table(control_group$religion)
# Actually, we can just drop it, as less than 5% are religious in our sample...
# That's a bit weird tbh

# Recode variables as binary 
table(control_group$marr_status)
control_group$married <- ifelse(control_group$marr_status=="married", 1, 0)
# control_group <- control_group %>% select(-race)
table(control_group$married)

# 
control_group <- control_group %>% select(-race, -religion, -marr_status)
control_group <- control_group %>% arrange(pid, year)
wide_control_covariates <- as.data.frame(t(control_group %>% 
                                             group_by(pid) %>%
                                             slice(1) %>%
                                             select( -worked_hours, -wage_pr)))
wide_control_covariates <- row_to_names(wide_control_covariates, 
                                        row_number = 1, remove_row = TRUE)
wide_control_covariates$na_count <- apply(wide_control_covariates, 1, function(x) sum(is.na(x)))
print(wide_control_covariates$na_count )
wide_control_covariates <- wide_control_covariates %>% select(-na_count)

# File 4: Time Invariant Convariates Controls
saveRDS(wide_control_covariates, file = "control_group_covariates.Rds")

