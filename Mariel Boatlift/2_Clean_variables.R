library(dplyr)
setwd("~/Documents/Oxford Yr 2/THESIS_micro_synth/Mariel boatlift shock/Analysis Data")
# START FROM HERE
d <- read.csv("panel_data.csv")
d <- subset(d , select = -X)
# Ok, now we have saved our results

# Note that ideally we want to work with outcome variable 
d <- d %>% arrange(pid, year)

# Matching of Current State
state_dictionary <- read.csv("state_dictionary.csv")
state_dictionary$number <- as.numeric(state_dictionary$number)
d <- d %>% left_join(state_dictionary, by = c("current_state" = "number")) %>% 
        mutate(current_state = state_name) 
d <- subset(d , select = -state_name)
# Matching of grew_up state
d <- d %>% left_join(state_dictionary, by = c("state_grew_up" = "number")) %>%
  mutate(state_grew_up = state_name) 
d <- subset(d , select = -state_name)


# Combine outcome variables and check for missing observations
# Sort dataset in a good way
d <- d %>% dplyr::select(pid, year, current_state, 
                   work_hours, head_hours_week,
                   head_avg_earnings_77_weird, head_reg_wage,
                   wife_avr_earnings_weird, wife_reg_wage, 
                   everything())

# Let us now remove some bad observations
# NOTE: worth exploring further the correlation in NAs between work_hours
# and head work_hours
d$pid <- as.character(d$pid)
missing_obs_hours <- d %>% filter(year <1985) %>%  
                           group_by(pid) %>% 
                           summarise( total_obs = length(work_hours), 
                                      total_zeros = sum(work_hours == 0 ),
                                      share_zeros = total_zeros/total_obs) %>%
                           select(pid, share_zeros)
 
# Leave only observations with enough data on working hours
d <- d %>% left_join(missing_obs_hours, by = c("pid" = "pid")) %>%
              filter(share_zeros < 0.35)
# Check if it's family or individual

# Get a rough estimate of number of treated units
treated_group <-  d %>% filter(current_state == "Florida")  %>%
                    select(pid, year, current_state, work_hours, 
                           head_avg_earnings_77_weird, head_occupation,
                           head_industry) %>% 
                    group_by(pid) %>%
                    summarise(total = length(pid))
# Not great, not terrible: around 40 treated observations pass all the 
# necessary requirments

# Sort out average wage issue - the idea is that we don't have an observation
# for 1977. I think this is because of the formula by which it is calculated
# In the ideal world, we will work with regular wage but we are imposing this
# Restriction, given that we have way too many missing values of regular wage
# However, it could be one potential extension
d <- d %>% group_by(pid) %>% 
           mutate(lead_avr_earnings = 
                  dplyr::lead(head_avg_earnings_77_weird, n =1, default = NA))
d$head_obs_avr_earnings <- 0
d$head_obs_avr_earnings[d$year<1977] <- d$head_avg_earnings_77_weird[d$year<1977]
d$head_obs_avr_earnings[d$year>1976] <- d$lead_avr_earnings[d$year>1976]
d <- d %>% filter(year < 1985)
d <- d %>% dplyr::select(-share_zeros, -lead_avr_earnings, -head_avg_earnings_77_weird,
                  -int_date, -release_num, -ID1968, -sequence, -pernum, 
                  -interview, -unemployed_hours, -seq_num)
d$head_obs_avr_earnings[is.na(d$head_obs_avr_earnings) & (d$head_reg_wage!=0)]<- 
  d$head_reg_wage[is.na(d$head_obs_avr_earnings) & (d$head_reg_wage!=0)] 

# Note that we have made the decision to ignore wives for the present moment
# So, let us drop all the wives-related variables
d <- d %>% select(-rel_to_head)
d  <- d %>% select( -wife_avr_earnings_weird, -wife_reg_wage, -wife_reg_wage,
                    -wife_ill,  -wife_industry, -wife_occupation, 
                    -wife_hours_week, -relation.head, -marr_pair_indic)

# Note that we impose a furthe restriction: we focus on outcome variables
# for total yearly hours worked AND average earnings instead of regular wages
# This simplifies stuff a lot
d <- d %>% select(-head_hours_week, -head_reg_wage) %>%
  select(pid, year, current_state, 
         work_hours, head_obs_avr_earnings, 
         everything())

# Remove all observations with less than 6 observations (could have been done
# In the initial stage of defining the panel)
d <- d %>% group_by(pid) %>% mutate(total_obs = length(pid)) %>% 
              filter(total_obs > 6)

# Next step is to clean each predictor in turn. 
# This would not be so easy for occupation and industry

# 1. Illness
d$ill <- 0
d$ill <- ifelse(d$head_ill > 0 , 1, 0)

# 2. Marriage Status
# Calculate married, singled, or smth else
d$marr_status[d$marr_status>1] <- 2
d$marr_status[d$marr_status==1] <- "married"
d$marr_status[d$marr_status==2] <- "single"
table(d$marr_status)

# 3. Race
table(d$race)
d$race[d$race>2] <- 3
d$race[d$race==1] <- "white"
d$race[d$race==2] <- "black"
d$race[d$race==3] <- "other"

# 4. religion
table(d$religion)
d$religion[(d$religion>0) & (d$religion<8)] <- "Protestant"
d$religion[(d$religion == 8)] <- "Catholic"
d$religion[(d$religion == 0) | (d$religion==9)] <- "Atheist/Other"

# 5. Industry
d$head_industry_disc <- cut(d$head_industry, 
                            breaks = c(0, 28, 57, 77, 398, 479, 698, 718, 
                                       759, 798, 809, 897, 937) , 
                            labels = c("Agriculture, Forestry, and Fisheries",
                                       "Mining", "Construction", "Manufacturing",
                                       "Transportation, Communications, and Other Public Utilities",
                                       "Wholesale and Retail Trade", 
                                       "Finance, Insurance, and Real Estate",
                                       "Business and Repair Services", 
                                       "Personal Services",
                                       "Entertainment and Recreation Services",
                                       "Professional and Related Services",
                                       "Public Administration" ),
                            right = T) 
table(d$head_industry)
table(d$head_industry_disc)
# All good with industry


# 6. Occupation
d$head_occupation_disc <- cut(d$head_occupation, 
                            breaks = c(1, 195, 245, 285, 395, 600, 695, 715,
                                       785, 802, 824, 965, 984) , 
                            labels = c("Professional, Technical, and Kindred Workers",
                                       "Managers and Administrators, Except Farm",
                                       "Sales Workers",
                                       "Clerical and Kindred Workers", 
                                       "Craftsmen and Kindred Workers",
                                       "Operatives, Except Transport", 
                                       "Transport Equipment Operatives",
                                       "Laborers, Except Farm",
                                       "Farmers and Farm Managers",
                                       "Farm Laborers and Farm Foremen",
                                       "Service Workers, Except Private Household",
                                       "Private Household Workers"),
                            right = T) 
table(d$head_occupation_disc)
table(d$head_occupation)

# 7. Father Occupation
d$father_occupation_disc <- cut(d$occupation_father, 
                              breaks = c(1,2,3,4,5,6,7,8,9,10) , 
                              labels = c("Professional, Technical, and Kindred Workers",
                                         "Managers, officials and proprietors",
                                         "Self-employed businessman",
                                         "Clerical and sales workers", 
                                         "Craftsmen and Kindred Workers",
                                         "Operatives and kindred workers", 
                                         "Laborers and service workers, farm laborers",
                                         "Farmers and farm managers",
                                         "Miscellaneous (armed services, protective workers)"),
                              right = FALSE) 
table(d$father_occupation_disc)
table(d$occupation_father)







