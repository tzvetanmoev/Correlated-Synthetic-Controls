### CREATES EXTRA TABLES FOR APPENDIX A.9

################################################################################
colSd <- function (x, na.rm=FALSE) apply(X=x, MARGIN=2, FUN=sd)
std <- function(x) sd(x)/sqrt(length(x))
colSe <- function (x, na.rm=FALSE) apply(X=x, MARGIN=2, FUN=std)
library(ggplot2)
library(ggpubr)
# install.packages("ggpubr")

#################### 1. WAGES   ################################################
wages1 <- read.csv("wage_csc_att.csv")
wages1 <- wages1[,-1]
colnames(wages1)[ncol(wages1)] <- "some_college"
wages1 <- wages1[-c(5, 11),]

# Clculate for low-skilled workers
wages_low_skilled_workers <- matrix(0, ncol = 3, nrow = ncol(wages1) - 2)
colnames(wages_low_skilled_workers) <- c("point", "low_ci", "high_ci")
low_skill <- wages1[wages1$some_college==0 , c(2:6)]
wages_low_skilled_workers[,1]<- colMeans( low_skill)
wages_low_skilled_workers[,2]<- -colSe(low_skill)*2.58 + wages_low_skilled_workers[,1]
wages_low_skilled_workers[,3] <- colSe(low_skill)*2.58 + wages_low_skilled_workers[,1]
wages_low_skilled_workers <- as.data.frame(wages_low_skilled_workers)
wages_low_skilled_workers$years <- c(1980, 1981, 1982, 1983, 1984)

# Clculate for low-skilled workers
wages_high_skilled_workers <- matrix(0, ncol = 3, nrow = ncol(wages1) - 2)
colnames(wages_high_skilled_workers) <- c("point", "low_ci", "high_ci")
high_skill <- wages1[wages1$some_college==1 , c(2:6)]
wages_high_skilled_workers[,1]<- colMeans( high_skill)
wages_high_skilled_workers[,2]<- -colSe(high_skill)*2.58 + wages_high_skilled_workers[,1]
wages_high_skilled_workers[,3] <- colSe(high_skill)*2.58 + wages_high_skilled_workers[,1]
wages_high_skilled_workers <- as.data.frame(wages_high_skilled_workers)
wages_high_skilled_workers$years <- c(1980, 1981, 1982, 1983, 1984)
# install.packages("ggplot2")
# a) Wage, Low-Skilled
a <- ggplot(wages_low_skilled_workers, 
            aes(years, point)) +
  geom_point(color = "#00AFBB") +
  geom_line(color = "#00AFBB") +
  geom_errorbar(aes(ymin = low_ci, ymax = high_ci), 
                width = 0.3, 
                color = "#00AFBB") +
  labs(x = "",
       y = "ATT") +
  geom_hline(yintercept = 0, colour = "darksalmon", linetype = "dashed") +
  ylim(-1.5, 0.2)+
  ggtitle("a) Wage, Low-Skilled") +
  theme(plot.title = element_text(hjust = 0.5))

# b) Wage, High-Skilled
b <- ggplot(wages_high_skilled_workers, aes(years, point)) +
  geom_point(color = "navy") +
  geom_line(color = "navy") +
  geom_errorbar(aes(ymin = low_ci, ymax = high_ci), 
                width = 0.3,
                color = "navy") +
  labs(x = "",
       y = "ATT") +
  geom_hline(yintercept = 0, colour = "darksalmon", linetype = "dashed") +
  ylim(-1.5, 0.4) +
  ggtitle("b) Wage, High-Skilled") +
  theme(plot.title = element_text(hjust = 0.5))
#################################################################################



#################################################################################
# Labour Supply Analysis
ls1 <- read.csv("ls_csc_att.csv")
ls1 <- ls1[,-1]
colnames(ls1)[ncol(ls1)] <- "some_college"
ls1 <- ls1[-c(5, 11),]

# Clculate for low-skilled workers
ls_low_skilled_workers <- matrix(0, ncol = 3, nrow = ncol(ls1) - 2)
colnames(ls_low_skilled_workers) <- c("point", "low_ci", "high_ci")
ls_low_skill <- ls1[ls1$some_college==0 , c(2:5)]
ls_low_skilled_workers[,1]<- colMeans( ls_low_skill)
ls_low_skilled_workers[,2]<- -colSe(ls_low_skill)*2.58 + ls_low_skilled_workers[,1]
ls_low_skilled_workers[,3] <- colSe(ls_low_skill)*2.58 + ls_low_skilled_workers[,1]
ls_low_skilled_workers <- as.data.frame(ls_low_skilled_workers)
ls_low_skilled_workers$years <- c(1980, 1981, 1982, 1983)

# Clculate for low-skilled workers
ls_high_skilled_workers <- matrix(0, ncol = 3, nrow = ncol(ls1) - 2)
colnames(ls_high_skilled_workers) <- c("point", "low_ci", "high_ci")
ls_high_skill <- ls1[ls1$some_college==1 , c(2:5)]
ls_high_skilled_workers[,1]<- colMeans( ls_high_skill)
ls_high_skilled_workers[,2]<- -colSe(ls_high_skill)*2.58 + ls_high_skilled_workers[,1]
ls_high_skilled_workers[,3] <- colSe(ls_high_skill)*2.69 + ls_high_skilled_workers[,1]
ls_high_skilled_workers <- as.data.frame(ls_high_skilled_workers)
ls_high_skilled_workers$years <- c(1980, 1981, 1982, 1983)

# c) LS, Low-Skilled
c<-ggplot(ls_low_skilled_workers, aes(years, point)) +
  geom_point(color = "seagreen3") +
  geom_line(color = "seagreen3") +
  geom_errorbar(aes(ymin = low_ci, ymax = high_ci), 
                width = 0.3,
                color = "seagreen3") +
  labs(x = "",
       y = "ATT" ) +
  geom_hline(yintercept = 0, colour = "darksalmon", linetype = "dashed") +
  ylim(-1500, 600) +
  ggtitle("c) LS, Low-Skilled") +
  theme(plot.title = element_text(hjust = 0.5))

# d) LS, High-Skilled
d <- ggplot(ls_high_skilled_workers, aes(years, point)) +
  geom_point(color = "darkgreen") +
  geom_line(color = "darkgreen") +
  geom_errorbar(aes(ymin = low_ci, ymax = high_ci), 
                width = 0.3,
                color  = "darkgreen") +
  labs(x = "",
       y = "ATT") +
  geom_hline(yintercept = 0, colour = "darksalmon", linetype = "dashed") +
  ylim(-1500, 600) +
  ggtitle("d) LS, High-Skilled") +
  theme(plot.title = element_text(hjust = 0.5))

figure <- ggarrange(a, b, c, d,
                    ncol = 2, nrow = 2)
figure