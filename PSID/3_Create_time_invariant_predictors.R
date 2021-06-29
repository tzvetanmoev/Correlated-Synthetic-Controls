# Prepare datast for analysis by transforming time-varying covariates into 
# time-invariant once
# Order dataset to be ready for analysis
d <- d %>% select(-head_industry, -head_occupation, 
                  -occupation_father, -head_ill)
d <- d %>% select(-total_obs)
d <- d %>% arrange(current_state, pid, year, work_hours, head_obs_avr_earnings,
                   race, ill)
# Note that living in Florida is like an indicator for being treated


# 1. Marr_status
d2 <- d
d <- d %>% group_by(pid) %>% count(marr_status) %>% 
  arrange(pid, -n, desc(marr_status) ) %>% 
  top_n(n=1) %>% 
  distinct(pid, .keep_all = T) %>% 
  rename(tinv_marr_status = marr_status) %>% 
  select(-n) %>% 
  right_join(d, by = "pid") %>% 
  select(-marr_status)


# Religion
d <- d %>% group_by(pid) %>% 
  count(religion) %>% 
  arrange(pid, -n, desc(religion) ) %>% 
  top_n(n=1) %>% 
  distinct(pid, .keep_all = T) %>% 
  select(-n) %>% 
  rename(tinv_religion = religion) %>% 
  right_join(d, by = "pid") %>% 
  select(-religion)

# Current state
d <- d %>% group_by(pid) %>% 
  count(current_state) %>% 
  arrange(pid, -n, desc(current_state) ) %>% 
  top_n(n=1) %>% 
  distinct(pid, .keep_all = T) %>% 
  select(-n) %>% 
  rename(tinv_current_state = current_state) %>% 
  right_join(d, by = "pid") %>% 
  select(-current_state) %>% 
  arrange(pid, year)

# State grew-up in - Constant over time, so no need to take action

# Kids - this is a bit trickier, as the variables used vary across waves
# THis should be looked at in more detail
# For this reason, I drop it
d  <- d %>% select(-kids, -int_num)

# Age - just take the average
d <- d %>% group_by(pid) %>% 
  summarise(mean_age = mean(age)) %>% 
  rename(tinv_age = mean_age) %>% 
  right_join(d, by = "pid") %>% 
  select(-age) %>% 
  arrange(pid, year)

# Years of Education - we need to take the mode
d <- d %>% filter(year > 1974) %>% 
  group_by(pid) %>% 
  count(educ) %>% 
  arrange(pid, -n, desc(educ) ) %>% 
  top_n(n=1) %>% 
  distinct(pid, .keep_all = T) %>% 
  select(-n) %>% 
  rename(tinv_education = educ) %>% 
  right_join(d, by = "pid") %>% 
  select(-educ) %>% 
  arrange(pid, year)

# Weight Take the average
d <- d %>% group_by(pid) %>% 
  summarise(mean_weight = mean(weight)) %>% 
  rename(tinv_weight = mean_weight) %>% 
  right_join(d, by = "pid") %>% 
  select(-weight) %>% 
  arrange(pid, year)

# Share of ill
d <- d %>% group_by(pid) %>% 
  mutate(share_ill = length(pid[ill>0])/length(pid) ) %>% 
  select(-ill)

# Head's Industry - we need to take the mode
d <- d %>% 
  group_by(pid) %>% 
  count(head_industry_disc) %>% 
  arrange(pid, -n, desc(head_industry_disc) ) %>% 
  drop_na(head_industry_disc) %>%
  top_n(n=1) %>% 
  distinct(pid, .keep_all = T) %>% 
  select(-n) %>% 
  rename(tinv_head_industry = head_industry_disc) %>% 
  right_join(d, by = "pid") %>% 
  select(-head_industry_disc) %>% 
  arrange(pid, year)

# Head's Occupation
d <- d %>% 
  group_by(pid) %>% 
  count(head_occupation_disc) %>% 
  arrange(pid, -n, desc(head_occupation_disc) ) %>% 
  drop_na(head_occupation_disc) %>%
  top_n(n=1) %>% 
  distinct(pid, .keep_all = T) %>% 
  select(-n) %>% 
  rename(tinv_head_occupation = head_occupation_disc) %>% 
  right_join(d, by = "pid") %>% 
  select(-head_occupation_disc) %>% 
  arrange(pid, year)


# Head's Father's Occupation
d <- d %>% 
  group_by(pid) %>% 
  count(father_occupation_disc) %>% 
  arrange(pid, -n, desc(father_occupation_disc) ) %>% 
  drop_na(father_occupation_disc) %>%
  top_n(n=1) %>% 
  distinct(pid, .keep_all = T) %>% 
  select(-n) %>% 
  rename(tinv_father_occupation_disc = father_occupation_disc) %>% 
  right_join(d, by = "pid") %>% 
  select(-father_occupation_disc) %>% 
  arrange(pid, year)

# Rearrange the columns in the right order
d2 <- d
# d <- d2
d <- d %>% select(pid, year, tinv_current_state, 
                  work_hours, head_obs_avr_earnings, everything())
d$tinv_current_state <- as.character(d$tinv_current_state)
d$tinv_current_state[d$tinv_current_state == "Florida"] <- "AAA_Florida"
table(d$tinv_current_state)
d <- d %>% arrange(tinv_current_state, pid, year)
# Note we change the name of the FLorida, so that it is the first on when 
# sorting out alphabetically. Obviously, not the prettiest solution...
colnames(d) <- c("pid", "year", "current_state", "worked_hours", 
                 "head_obs_avr_earnings", "father_occupation", 
                 "head_occupation", "head_industry", "weight", "education" ,           
                 "age", "religion", "marr_status", "race", "state_grew_up", 
                 "share_ill")

# Lastly, save time-invariant predictors from treatment group
treatment_group <- d[d$current_state=="AAA_Florida",]
saveRDS(treatment_group, file="treatment_group.Rds")

saveRDS(d, file="full_data_for_application.Rds")

