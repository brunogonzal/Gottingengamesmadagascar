# Code by Katie Devenish. ktd19ycv@bangor.ac.uk

# This code is used in the publication entitled "Madagascar's biggest mine is set to achieve No Net Loss of Forest"
# Authors: Katie Devenish, Sebastien Desbureaux, Simon Willcock and Julia PG Jones. 2021.

# Here, we explore which of the modelling choices tested in the robustness checks most influence the regression estimates. 

# Load libraries #

library("stargazer")


# --------------------- 1) Data Construction for Matching -----------------------------------------#

setwd("c:/Users/desbures/Documents/Sebastien/11_Autres/Katie/Ambatovy/Output_data/")

# We load all csv from the different loop
temp = list.files(pattern="*.csv")
myfiles = lapply(temp, read.csv)

# We include the name of the file to keep track of the offset
for (i in 1:12) {
  myfiles[[i]]$csv <- temp[i]
}

# We create one data frame
d <- do.call(rbind, myfiles)

# We create a variable for the offfset
d$offset <- substring(d$csv,5,7)
  d$offset[d$offset=="CFA"] <- "CFAM"
  d$offset[d$offset=="CZ_"] <- "CZ"
  d$offset[d$offset=="Ank"] <- "ANK"  
d$offset <- as.factor(d$offset)
d$csv <- NULL
    
## We indentify valid models as in the main text
ANK$less_10pc_unmatched  <- ifelse(ANK$unmatched<2862*0.1,T,F)

d$less_10pc_unmatched <- ifelse(d$offset=="ANK" & d$unmatched<2862*0.1,T,
                                ifelse(d$offset=="CFAM" & d$unmatched<2626*0.1,T,
                                       ifelse(d$offset=="CZ" & d$unmatched<1340*0.1,T,
                                              ifelse(d$offset=="TTF" & d$unmatched<1170*0.1,T,F))))


# A simple reg
reg1 <- lm(coef ~ . , data=d[,-c(1,3:7,24)])
  summary(reg1)

reg2 <- lm(coef ~ . , data=d[d$less_10pc_unmatched==T,-c(1,3:7,24)])
  summary(reg2)
  
reg3 <- lm(coef ~ . , data=d[d$mean_std_diff<=0.25,-c(1,3:7,24)])
  summary(reg3)
  
reg4 <- lm(coef ~ . , data=d[d$max_std_diff<=0.25,-c(1,3:7,24)])
  summary(reg4)

reg5 <- lm(coef ~ . , data=d[d$parralel_trend==T,-c(1,3:7,24)])
  summary(reg5)  
    
reg6 <- lm(coef ~ . , data=d[d$offset=="ANK",-c(1,3:7,23,24)])
  summary(reg6)
  
reg7 <- lm(coef ~ . , data=d[d$offset=="CFAM",-c(1,3:7,23,24)])
  summary(reg7)
  
reg8 <- lm(coef ~ . , data=d[d$offset=="CZ",-c(1,3:7,23,24)])
  summary(reg8)
  
reg9 <- lm(coef ~ . , data=d[d$offset=="TTF",-c(1,3:7,23,24)])
  summary(reg9)


stargazer(reg1,reg2,reg3,reg4,reg5,reg6,reg7,reg8,reg9,
          title = "Effect moddeling choices on estimated impact",
          out="mod_choices.doc",
          type = "html")
