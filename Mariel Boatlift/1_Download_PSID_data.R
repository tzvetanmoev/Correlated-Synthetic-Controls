# Note that this follows largely the procedure outline in the README file from
# psidR's github directory: https://github.com/floswald/psidR
# It is extremely conveneint... Let's pray it works tho


library(data.table)
# install.packages("psidR")
library(psidR)
#  library(RCurl)
# install.packages("RCurl")

# install.packages("futile.logger", dependencies=TRUE, repos='http://cran.rstudio.com/')
# library(futile.logger)
library(foreign)
# library(SAScii)
# library(openxlsx)
# library(stats)
# library(utils)
# library(lodown)
# library(plyr)
library(dplyr)
library(collapse)
library(tidyr)
# install.packages("janitor")
library(janitor)

# install.packages("collapse")

# Create variables
setwd("~/Documents/Oxford Yr 2/THESIS_micro_synth/Mariel boatlift shock/Analysis Data")
r <- system.file(package="psidR")
f <- fread(file.path(r,"pls_work","famvars.txt"))
i <- fread(file.path(r,"pls_work","indvars.txt"))

# Note alternative strategy: get variables from table
# wf <- read.xlsx("http://psidonline.isr.umich.edu/help/xyr/psid.xlsx")
# can be worth a try later on

# Get data in good format
new_i = dcast(i[,list(year,name,variable)],year~name, value.var = "variable")

# Remove some annoying duplicates
f <- f[(f$name!="int_num") & 
         (f$label!="# OF KIDS OF HD") & 
         (f$variable!="V4765"),]
new_f = dcast(f[,list(year,name,variable)],year~name, value.var = "variable")

new_f[2,16]<- "V3233"
d = build.panel(fam.vars=new_f, ind.vars=new_i,
                heads.only = TRUE, sample="SRC",design="all")
# Note we are getting only the HEADS (as heads.only = TRUE) and we are working
# with unbalanced panel data (as design="all")
# As a result, we have rel_head = 1 for all observation!

# Alternative sampling
# d2 = build.panel(fam.vars=new_f,
#                  ind.vars=new_i,
#                  heads.only = FALSE, sample="SRC",design=7)

# Most restrictive sample would be something like
# d3 = build.panel(fam.vars=new_f, ind.vars=new_i,
#                  heads.only = TRUE, sample="SRC",design=7)


# Save the data
write.csv(d,
          file="panel_data_2.csv")
save(d, file="pls_work.Rds")

##################################################################################