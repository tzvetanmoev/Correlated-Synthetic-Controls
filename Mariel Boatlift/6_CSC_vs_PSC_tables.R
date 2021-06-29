full_results <- matrix(0, ncol = 6, nrow = 12)
colnames(full_results) <- c("Number T", "Year",  "CSC", "PSC", "CSC", "PSC")
full_results[,1] <- c(4,4,4,  3,3,3,  2,2,2,  1,1,1)

# T = 4
# note the space after csv
ls_4 <- as.matrix(read.csv("rmse_training_period_ls_4.csv "))
full_results[c(1:3),c(2:4)] <- ls_4[,-1]
wages_4 <- as.matrix(read.csv("rmse_training_period_wages_4.csv "))
full_results[c(1:3),c(5:6)] <- wages_4 [,-c(1,2)]

# T = 3
ls_3 <- as.matrix(read.csv("rmse_training_period_ls_3.csv"))
full_results[c(4:6),c(2:4)] <- ls_3[,-1]
wages_3 <- as.matrix(read.csv("rmse_training_period_wages_3.csv"))
full_results[c(4:6),c(5:6)] <- wages_3[,-c(1,2)]

# T = 2
ls_2 <- as.matrix(read.csv("rmse_training_period_ls_2.csv"))
full_results[c(7:9),c(2:4)] <- ls_2[,-1]
wages_2 <- as.matrix(read.csv("rmse_training_period_wages_2.csv"))
full_results[c(7:9),c(5:6)] <- wages_2[,-c(1,2)]


# T = 1
ls_1 <- as.matrix(read.csv("rmse_training_period_ls_1.csv"))
full_results[c(10:12),c(2:4)] <- ls_1[,-1]
wages_1 <- as.matrix(read.csv("rmse_training_period_wages_1.csv"))
full_results[c(10:12),c(5:6)] <- wages_1[,-c(1,2)]

stargazer(full_results, summary = F, digits = 2)


