# Orginal code by Katie Devenish# and Sebastien Desbureaux. 
# ktd19ycv@bangor.ac.uk

#|# Code modified by Yannic Damm in April 2025, contact: yannic.damm@ilr.uni-bonn.de

# This code is used in the publication entitled "Madagascar's biggest mine is set to achieve No Net Loss of Forest"
# Authors: Katie Devenish, Sebastien Desbureaux, Simon Willcock and Julia Jones. 2021.

# Load libraries #

library("foreign")
library("dplyr")
library("xlsx")
library("ggplot2")
library("broom")
library("MatchIt")
library("ggrepel")
library("calibrate")
library("plm")
library("gridExtra")
library("tidyr")
library("gplots")
library("lmtest")

#|# Add fixest library for regressions with robust standard errors
library("fixest")

# Code Structure

# 1) Data Construction for Matching
# 2) Matching
# 3) Extract Matched Pairs from matching output
# 4) Data construction for DiD regression using matched pairs
# 5) DiD Regression 
# 6) Quantifying Avoided Deforestation
# 7) Fixed Effects Panel Regression


# --------------------- 1) Data Construction for Matching -----------------------------------------#

# Function to read in input data for each offset, remove unnecessary columns and add column to indicate treatment status
# and offset of origin. 

tidy_data <-function(path, name, number, label){
  name = read.dbf(path)
  name = subset(name, select= -c(1, 16))
  name$treated = number              # Number = 1 for pixels from an offset and 0 for control pixels. 
  name$offset = label
  return(name)
}

TTF <- tidy_data("Sample_TTF3.dbf", TTF, 1, "TTF")           
ANK <- tidy_data("Sample_ANK3.dbf", ANK, 1, "ANK")
CFAM <- tidy_data("Sample_CFAM3.dbf", CFAM, 1, "CFAM")
CZ <- tidy_data("Sample_CZ3.dbf", CZ, 1, "CZ")

# TTF = Torotorofotsy
# ANK = Ankerana
# CFAM = Corridor Forestier Analamay-Mantadia
# CZ = Conservation Zone


# Load Control separately because columns are different # 

Control = read.dbf("Final_control.dbf")
Control$treated = 0
Control$offset = "Cont"
Control <- subset(Control, select = -c(4))        # Remove unwanted fire variable
Control <- Control[,c(2,3,1,4:16)]                # Re-order columns to match offset dataframes

# Rename columns #

cols <- c("X","Y","Tree_loss", "Pop_density", "Dist_sett", "Slope", "Elevation", "Aspect", "Annual_Rain",
          "Dist_track", "Dist_road", "Dist_river", "Dist_edge", "Dist_defor", "treated", "offset")

names(Control) <- paste0(cols)
names(ANK) <- paste0(cols)
names(CZ) <- paste0(cols)
names(CFAM) <- paste0(cols)
names(TTF) <- paste0(cols)

# Merge each offset with the control dataset 
# These offset + control datasets will be the input for matching #

ANKCONT <- rbind(ANK, Control)
CZCONT <- rbind(CZ, Control)
CFAMCONT <- rbind(CFAM, Control)
TTFCONT <- rbind(TTF, Control)


# Data cleaning

# - remove erroneous -9999 values from missing values in the raster data layers
# - create calipers to drop extreme observations that will never be matched
# - caliper based on the distribution of values for the 5 essential 
#   covariates within the treated sample.

Seb_dataclean <- function(data, x){
  
  # Replace all -9999 with NA
  data <- data %>% mutate(across(everything(), ~ifelse(. == -9999, NA, .)))
  
  # We keep observations without na
  data <- data %>% drop_na()
  
  # We determine the caliper
  caliper_dist_def<-x*sd(data$Dist_defor[data$treated==1])
  caliper_slope<-x*sd(data$Slope[data$treated==1])
  caliper_elev<-x*sd(data$Elevation[data$treated==1])
  caliper_edge<-x*sd(data$Dist_edge[data$treated==1])
  caliper_road <- x*sd(data$Dist_road[data$treated==1])
  
  # We keep control observations that are on the same support
  data <- data[data$Dist_defor < max(data$Dist_defor[data$treated==1])+caliper_dist_def,]  
  data <- data[data$Dist_defor > min(data$Dist_defor[data$treated==1])- caliper_dist_def,]
  
  data <- data[data$Slope < max(data$Slope[data$treated==1])+caliper_slope,]      
  data <- data[data$Slope > min(data$Slope[data$treated==1])- caliper_slope,]
  
  data <- data[data$Elevation < max(data$Elevation[data$treated==1])+caliper_elev,]  
  data <- data[data$Elevation > min(data$Elevation[data$treated==1])-caliper_elev,]  
  
  data <- data[data$Dist_edge < max(data$Dist_edge[data$treated==1])+caliper_edge,]
  data <- data[data$Dist_edge > min(data$Dist_edge[data$treated==1])-caliper_edge,]
  
  data <- data[data$Dist_road < max(data$Dist_road[data$treated==1])+caliper_road,]
  data <- data[data$Dist_road > min(data$Dist_road[data$treated==1])-caliper_road,]
  
  
}

ANKCONT <- Seb_dataclean(data = ANKCONT, x= 1)
ANKCONT$ID <- seq(nrow(ANKCONT))              # Add column for observation ID
rownames(ANKCONT) <- ANKCONT$ID

CZCONT <- Seb_dataclean(CZCONT, 1)
CZCONT$ID <- seq(nrow(CZCONT))
rownames(CZCONT) <- CZCONT$ID

CFAMCONT <- Seb_dataclean(CFAMCONT, 1)
CFAMCONT$ID <- seq(nrow(CFAMCONT))
rownames(CFAMCONT) <- CFAMCONT$ID 

TTFCONT <- Seb_dataclean(TTFCONT, 1)
TTFCONT$ID <- seq(nrow(TTFCONT))
rownames(TTFCONT) <- TTFCONT$ID


#----------------------------2) The Matching -----------------------------------------#

# The matching algorithm and selected arguments Ho et al (2007). 
# Matching run separately for each offset.


variables <- c("Slope", "Elevation", "Dist_road", "Dist_edge", "Dist_defor")

# Set caliper for each covariate

cal <- rep(1, length(variables))        
names(cal) <- variables

# The main matching specification. 1:1 nearest-neighbour matching without replacement using Mahalanobis distance measure
# and a caliper of 1sd. The other 4 specifications tested for balance can be found in the script entitled "Choosing_matching_spec"


matching <- function(offset){
  output = matchit(treated ~ Slope + Elevation + Dist_road + Dist_edge + Dist_defor,
                   data= offset, method= "nearest", distance = "mahalanobis", replace = FALSE, caliper=cal)
  return(output)
}          

m.out.ANK <- matching(offset = ANKCONT)
m.out.CFAM <- matching(offset = CFAMCONT)
m.out.CZ <- matching(offset = CZCONT)
m.out.TTF <- matching(offset = TTFCONT)


# Summarise results #

List1 <- list(m.out.ANK, m.out.CFAM, m.out.CZ, m.out.TTF)

sum1_matching <- function(matched_data){
  sum1 = summary(matched_data, standardize = TRUE)      # Standardise = TRUE to get standardised mean difference in output
  return(sum1)
}

Output_names <- c("m.out.ANK", "m.out.CFAM", "m.out.CZ", "m.out.TTF")

allOutputs <- lapply(List1, sum1_matching)              # Apply the summary function to all matched datasets in List1

names(allOutputs) <- paste0(Output_names)  


# Repeat but with standardise = FALSE to get the mean EQQ values. This is a second measure to assess covariate balance. #

sum2_matching <- function(matched_data){
  sum2 = summary(matched_data, standardize = FALSE)
  return(sum2)
}

allOutputs2 <- lapply(List1, sum2_matching)
names(allOutputs2) <- paste0(Output_names)

# Function to extract useful summary values and combine into one dataset # 

cov_balance <- function(x, y){
  test <- rbind(allOutputs[[x]]$sum.all, allOutputs[[x]]$sum.matched)   # Combine the non-matched and matched values into one dataset #
  test <- subset(test, select= c(1,2,3))    # select means treated, means control and standardised mean difference columns #
  colnames(test)[3] <- c(y)
  test <- data.frame(test)
  test$sample <- 0
  test$sample[1:(nrow(test)/2)] = "Non-matched"           
  test$sample[((nrow(test)/2)+1):nrow(test)] = "matched"
  test$covariates <- 0
  test[1:(nrow(test)/2), 5] <- variables        # Add column with the variables. 
  test[((nrow(test)/2)+1):nrow(test), 5] <- variables
  test2 <- rbind(allOutputs2[[x]]$sum.all, allOutputs2[[x]]$sum.matched)  # Combine the non-matched and matched values from allOutputs 2 (this has the mean EQQ values) #
  test2 <- subset(test2, select= c(5))                # Extract mean EQQ column
  test3 <- cbind(test, test2) # Add mean EQQ column to the test. This produces a dataset of means treated, means control, standardised mean difference, mean EQQ values for each covariate for the non-matched and matched data. This can be used to assess covariate balance post-matching.
}

cov_balance_ANK <- cov_balance(x=1, "SMD_1_cal")
cov_balance_CFAM = cov_balance(x=2, "SMD_1_cal")
cov_balance_CZ = cov_balance(x=3, "SMD_1_cal")
cov_balance_TTF = cov_balance(x=4, "SMD_1_cal")


#---------------------------3) Extract matched pairs from matching output ------------------------------------#

# Extract matched pairs from input data


m.data.ANK <- match.data(m.out.ANK, group = "all", distance = "distance", weights = "weights", subclass = "subclass",
                         data = ANKCONT, drop.unmatched = TRUE)

m.data.CZ <- match.data(m.out.CZ, group = "all", distance = "distance", weights = "weights", subclass = "subclass",
                        data = CZCONT, drop.unmatched = TRUE)

m.data.CFAM <- match.data(m.out.CFAM, group = "all", distance = "distance", weights = "weights", subclass = "subclass",
                          data = CFAMCONT, drop.unmatched = TRUE)

m.data.TTF <- match.data(m.out.TTF, group = "all", distance = "distance", weights = "weights", subclass = "subclass",
                         data = TTFCONT, drop.unmatched = TRUE)


# ---------------------------4) Data Construction for DiD Regression using Matched Pairs -------------------------#


# Aggregate pixels into treated (offset) and control samples and tabulate observations within each sample by tree loss year. 
# This gives the count of pixels within each sample deforested each year. 

annual_defor_ANK <- data.frame(table(factor(m.data.ANK$Tree_loss, levels = 1:19), m.data.ANK$offset))

# levels = 1:19 removes observations with 0 value for Tree Loss year (which were not deforested over the study period).
# This is because we are interested in comparing deforestation outcomes between offsets and the matched controls. 

annual_defor_CZ <- data.frame(table(factor(m.data.CZ$Tree_loss, levels = 1:19), m.data.CZ$offset))
annual_defor_CZ$Var2 <- factor(annual_defor_CZ$Var2, levels = c("CZ", "Cont"))    # For the plots the order of the factors needs to match the other offsets

annual_defor_CFAM <- data.frame(table(factor(m.data.CFAM$Tree_loss, levels = 1:19), m.data.CFAM$offset))

annual_defor_TTF <- data.frame(table(factor(m.data.TTF$Tree_loss, levels = 1:19), m.data.TTF$offset))
annual_defor_TTF$Var2 <- factor(annual_defor_TTF$Var2, levels = c("TTF", "Cont"))

label <- c("Year", "Sample", "Annual_Deforestation")
names(annual_defor_ANK) <- paste0(label)
names(annual_defor_CFAM) <- paste0(label)
names(annual_defor_CZ) <- paste0(label)
names(annual_defor_TTF) <- paste0(label)

# But offsets are different sizes - need to calculate annual deforestation as a percentage of total pixels in the sample to plot

# Calculating Percentage Annual Deforestation with 1:1 matching.

annual_defor_ANK$Perc_Annual_Defor <- (annual_defor_ANK$Annual_Deforestation/(nrow(m.data.ANK)/2))*100          # In the matched dataset the no. of control pixels = No. of treatment pixels because I matched 1:1.  
annual_defor_CFAM$Perc_Annual_Defor <- (annual_defor_CFAM$Annual_Deforestation/(nrow(m.data.CFAM)/2))*100       # Therefore, nrow(m.data.CFAM)/2 is the total number of treatment and the total number of control pixels 
annual_defor_CZ$Perc_Annual_Defor <- (annual_defor_CZ$Annual_Deforestation/(nrow(m.data.CZ)/2))*100             # in the matched dataset. 
annual_defor_TTF$Perc_Annual_Defor <- (annual_defor_TTF$Annual_Deforestation/(nrow(m.data.TTF)/2))*100



#-------------------------------5) Difference in Differences Regression -----------------------------------------#

                                 # a) Data Construction 
                                  

Data_construction_DiD <- function(offset, y){    # y corresponds to the year of protection - so Time = 0 before protection and 1 after protection. 
  offset$Year <- as.numeric(rep(1:19,2))        # Have to make Year numeric for >= to work
  offset$TimeF <- factor(ifelse(offset$Year >= y, 1,0))
  offset$TreatedF <- factor(ifelse(offset$Sample != "Cont", 1,0))
  return(offset)
}

annual_defor_ANK <- Data_construction_DiD(annual_defor_ANK, 11)
annual_defor_CFAM <- Data_construction_DiD(annual_defor_CFAM, 13)
annual_defor_CZ <- Data_construction_DiD(annual_defor_CZ, 9)
annual_defor_TTF <- Data_construction_DiD(annual_defor_TTF, 14)


                              # b) Outcome variable transformation

# log(y+1) transformation of outcome variable required because non-normal 
# properties of count data violate assumptions of 
# homoscedascity of linear models.

annual_defor_ANK$log_annual_defor <- log(annual_defor_ANK$Annual_Deforestation + 1)
annual_defor_CZ$log_annual_defor <- log(annual_defor_CZ$Annual_Deforestation + 1)
annual_defor_CFAM$log_annual_defor <- log(annual_defor_CFAM$Annual_Deforestation +1)
annual_defor_TTF$log_annual_defor <- log(annual_defor_TTF$Annual_Deforestation +1)

#|# b.2) Outcome variable transformation

#|# asinh (inverse hyperbolic sine, ihs) transformation of outcome variable

annual_defor_ANK$ihs_annual_defor <- asinh(annual_defor_ANK$Annual_Deforestation)
annual_defor_CZ$ihs_annual_defor <- asinh(annual_defor_CZ$Annual_Deforestation)
annual_defor_CFAM$ihs_annual_defor <- asinh(annual_defor_CFAM$Annual_Deforestation)
annual_defor_TTF$ihs_annual_defor <- asinh(annual_defor_TTF$Annual_Deforestation)

#|# b.3) Outcome variable transformation

#|# log transformation of outcome variable, converting zero deforestation to NAs, This examines the intensive deforestation margin

annual_defor_ANK$log.annual_defor <- log(annual_defor_ANK$Annual_Deforestation)
annual_defor_CZ$log.annual_defor <- log(annual_defor_CZ$Annual_Deforestation)
annual_defor_CFAM$log.annual_defor <- log(annual_defor_CFAM$Annual_Deforestation)
annual_defor_TTF$log.annual_defor <- log(annual_defor_TTF$Annual_Deforestation)


                              # c) Test for parallel trends

# Parallel trends in outcomes between treated and control samples in the years before the intervention
# is a key assumption of difference-in-differences regressions.

# Use only data from the years before the offsets were protected #

ANK_data_before <- annual_defor_ANK[(annual_defor_ANK$Year <11),]
CFAM_data_before <- annual_defor_CFAM[(annual_defor_CFAM$Year <13),]
CZ_data_before <- annual_defor_CZ[(annual_defor_CZ$Year <9),]
TTF_data_before <- annual_defor_TTF[(annual_defor_TTF$Year <14),]


# ANK #

ANKa <- lm(log_annual_defor ~ Year*TreatedF, data= ANK_data_before)
summary(ANKa)       # If the interaction between Year and TreatedF is not significant, there is no 
                    # significant difference in the relationship between Year and the log-transformed
                    # count of deforestation between treated and control samples --> parallel trends assumption holds.

# CFAM #

CFAMa <- lm(log_annual_defor ~ Year*TreatedF, data= CFAM_data_before)
summary(CFAMa)

# There is a significant difference in the trend in deforestation over time between treated and control samples 
# No parallel trends. 


# CZ #

CZa <- lm(log_annual_defor ~ Year*TreatedF, data= CZ_data_before)
summary(CZa)


# TTF #

TTFa <- lm(log_annual_defor ~ Year*TreatedF, data= TTF_data_before)
summary(TTFa)


# All showed parallel trends except CFAM which cannot be used in individual DiD regressions



                              # d) DiD Regression

# Formula = y ~ treatment + time + (treatment x time)

# Interaction between treated and time is the coefficent of interest. This represents the effect
# of an observation being in an offset, after protection on the log-transformed count of deforestation.
# If this is significant and negative it means protection significantly reduced deforestation within the offset, 
# relative to the counterfactual.


# ANK # 

#|# log(y+1)
modelANK <- lm(log_annual_defor ~ TreatedF*TimeF, data= annual_defor_ANK)
summary(modelANK) 

#|# log(y+1) with robust SEs
modelANK_robust_se <- feols(log_annual_defor ~ TreatedF*TimeF, data = annual_defor_ANK, vcov = "hetero")
summary(modelANK_robust_se) 

#|# ihs with robust SEs
modelANK_ihs <- feols(ihs_annual_defor ~ TreatedF*TimeF, data = annual_defor_ANK, vcov = "hetero")
summary(modelANK_ihs) 

#|# poisson regression with robust SEs
modelANK_poisson<- fepois(Annual_Deforestation ~ TreatedF*TimeF, data = annual_defor_ANK, vcov = "hetero")
summary(modelANK_poisson) 

#|# log(y) with robust SEs
modelANK_log.y <- feols(log.annual_defor ~ TreatedF*TimeF, data = annual_defor_ANK, vcov = "hetero")
summary(modelANK_log.y) 


# Back-transform the estimate to get the treatment effect - 
# the  percentage difference in average annual deforestation between
# the offset and the estimated counterfactual following protection. 
# The estimated counterfactual is the average annual deforestation in the matched control 
# sample after the intervention, adjusted to 
# to account for pre-intervention differences between the two samples. 

exp(coef(modelANK)[4])-1
exp(confint(modelANK)[4,])-1

# CZ # 

#|# log(y+1)
modelCZ <- lm(log_annual_defor ~ TreatedF*TimeF, data= annual_defor_CZ)
summary(modelCZ) 

#|# log(y+1) with robust SEs
modelCZ_robust_se <- feols(log_annual_defor ~ TreatedF*TimeF, data = annual_defor_CZ, vcov = "hetero")
summary(modelCZ_robust_se) 

#|# ihs with robust SEs
modelCZ_ihs <- feols(ihs_annual_defor ~ TreatedF*TimeF, data = annual_defor_CZ, vcov = "hetero")
summary(modelCZ_ihs) 

#|# poisson regression with robust SEs
modelCZ_poisson<- fepois(Annual_Deforestation ~ TreatedF*TimeF, data = annual_defor_CZ, vcov = "hetero")
summary(modelCZ_poisson) 

#|# log(y) with robust SEs
modelCZ_log.y <- feols(log.annual_defor ~ TreatedF*TimeF, data = annual_defor_CZ, vcov = "hetero")
summary(modelCZ_log.y) 


exp(coef(modelCZ)[4])-1
exp(confint(modelCZ)[4,])-1

# TTF #

#|# log(y+1)
modelTTF <- lm(log_annual_defor ~ TreatedF*TimeF, data= annual_defor_TTF)
summary(modelTTF) 

#|# log(y+1) with robust SEs
modelTTF_robust_se <- feols(log_annual_defor ~ TreatedF*TimeF, data = annual_defor_TTF, vcov = "hetero")
summary(modelTTF_robust_se) 

#|# ihs with robust SEs
modelTTF_ihs <- feols(ihs_annual_defor ~ TreatedF*TimeF, data = annual_defor_TTF, vcov = "hetero")
summary(modelTTF_ihs) 

#|# poisson regression with robust SEs
modelTTF_poisson<- fepois(Annual_Deforestation ~ TreatedF*TimeF, data = annual_defor_TTF, vcov = "hetero")
summary(modelTTF_poisson) 

#|# log(y) with robust SEs
modelTTF_log.y <- feols(log.annual_defor ~ TreatedF*TimeF, data = annual_defor_TTF, vcov = "hetero")
summary(modelTTF_log.y) 


exp(coef(modelTTF)[4])-1
exp(confint(modelTTF)[4,])-1

# Results show a significant reduction in average annual deforestation of 96% (- 89 to - 98%) in Ankerana
# and 66% (-27 to -84%) in the COnservation Zone following protection. 

# In Torotorofotsy, protection had no significant effect on deforestation. 

DiD_res <- rbind(t(tidy(modelANK)), t(tidy(modelCZ)), t(tidy(modelTTF)))



#               e) Extract back-transformed coefficients and confidence intervals from DiD models


tidy_data5e <- function(model, name){
  b <- (exp(coef(model)[4])-1)*100              # Back-transform coefficient of interaction term (exp(coefficient)-1)*100 to convert to % difference
  c <- (exp(confint(model)[4,c(1,2)])-1)*100    # Back-transform confidence intervals
  c <- data.frame(t(c))
  b <- data.frame(b)
  b <- cbind(b,c)
  names(b) <- c("ATT", "Lower_CI", "Upper_CI")
  rownames(b) <- name
  return(b)
}

ATT_ANK <- tidy_data5e(modelANK, "ANK")
ATT_CZ <- tidy_data5e(modelCZ, "CZ")
ATT_TTF <- tidy_data5e(modelTTF, "TTF")

ATT_all <- rbind(ATT_ANK, ATT_CZ, ATT_TTF)



# -----------------------------6) Quantify Avoided Deforestation -----------------------------------------------------------------# 

# Read in datasets containing forest cover and annual forest loss values for whole offset area. 
# This is so we can use the estimated average treatment effect (ATT) to convert the actual deforestation observed
# within the offsets following protection to counterfactual levels. The difference between these two values is the amount
# of deforestation which has been avoided through protection. 


ANK_dat <- read.dbf("ANK_var.dbf")      
CFAM_dat <- read.dbf("CFAM_var.dbf")
CZ_dat <- read.dbf("CZ_var.dbf")
CZ_dat$VALUE_5 <- 0                                # Fill in missing value. No tree loss in CZ in 2005
CZ_dat <- CZ_dat[,c(1:7,34,8:33)]                 # Re-order columns to match other datasets
TTF_dat <- read.dbf("Torotorofotsy_var.dbf")

#|# Extract initial forest areas
initial_forest_areas <- list(
  "initial_forest_area_ANK" = ANK_dat$AREA_FOR_H
  ,"initial_forest_area_CFAM" = CFAM_dat$AREA_FOR_H
  ,"initial_forest_area_CZ" = CZ_dat$AREA_FOR_H
  ,"initial_forest_area_TFF" = TTF_dat$AREA_FOR_H
  ,"total_initial_forest_area" = ANK_dat$AREA_FOR_H + CFAM_dat$AREA_FOR_H + CZ_dat$AREA_FOR_H +TTF_dat$AREA_FOR_H
)


            # a) Calculate observed, counterfactual and avoided deforestation (plus Upper and Lower CIs)


# Calculate for each offset that showed parallel trends on which the site-based DiD regression was run.
# Calculate for years following protection. 

tidy_data7 <- function(data, model, y){                     
  data <- data[ ,4:22]                                         # Extract only Forest Loss per Year columns, excluding Year 0.
  data <- data.frame(t(data))                                  # Transpose
  names(data) <- "Tree_loss"
  data$Tree_loss <- data$Tree_loss/10000                      # Convert tree loss in m2 to hectares
  data$Year <- 1:19
  data <- data[data$Year>= y,]                                # Remove years before protection of the offset
  rownames(data) <- 1:nrow(data)
  
  if (length(coef(model)) > 1){
    b      <- coef(model)[4]
    b_upr  <- confint(model)[4,2]
    b_lwr  <- confint(model)[4,1]
  } else if (length(coef(model)) == 1){
    b      <- coef(model)[[1]]
    b_upr  <- confint(model)[[1,2]]
    b_lwr  <- confint(model)[,1]
  } 
  
  data$counterfactual_defor <- (1/(exp(b)))*data$Tree_loss       # Multiply annual tree loss in hectares. eg. In Ankerana, observed deforestation in the offset after protection was 95.8% lower than the counterfactual.
  data$avoided_defor <- data$counterfactual_defor - data$Tree_loss            # Observed defor was 4.14% of the counterfactual. To scale up to get the amount of deforestation which would have occurred under the 
  
  data$counterfactual_upr <- (1/(exp(b_upr)))*data$Tree_loss    # the counterfactual scenario need to do 1/0.0414 and multiply by observed deforestation. 
  data$avoided_defor_upr <- data$counterfactual_upr - data$Tree_loss          # This is the same as doing data$Tree_loss/exp(estimate)
  
  data$counterfactual_lwr <- (1/(exp(b_lwr)))*data$Tree_loss
  data$avoided_defor_lwr <- data$counterfactual_lwr - data$Tree_loss
  
  return(data)
}


ANK_forest <- tidy_data7(ANK_dat, modelANK, 11)            
CZ_forest <- tidy_data7(CZ_dat, modelCZ, 9)
TTF_forest <- tidy_data7(TTF_dat, modelTTF, 14)

ANK_forest[10,] <- 0                                                    # Create new row for sums
ANK_forest[10, c(1,3:8)] <- apply(ANK_forest[, c(1,3:8)], 2, sum)       # Calculate total observed, counterfactual and avoided deforestation following protection.

CZ_forest[12,] <- 0
CZ_forest[12, c(1,3:8)] <- apply(CZ_forest[, c(1,3:8)], 2, sum)

TTF_forest[7,] <- 0
TTF_forest[7, c(1,3:8)] <- apply(TTF_forest[, c(1,3:8)], 2, sum)


# b) Combine above datasets into one for plotting
b <- rbind(ANK_forest[10,], CZ_forest[12,], TTF_forest[7,])
b$Offset <- c("ANK", "CZ", "TTF")
b <- b[, c(9, 1, 3, 5, 7, 4, 6, 8)]                                   # Re-arrange columns
b[, 6:8] <- b[, 6:8]*-1                                               # Multiply by -1 to turn avoided deforestation columns into 'impact on deforestation' which is -ve if the offset has reduced deforestation and positive if it has increased it.
names(b) <- c("Offset", "Observed_defor", "Counterfactual_Defor", "Counterfactual_Defor_Upper", "Counterfactual_Defor_Lower", "Impact_defor", "Impact_defor_Upper", "Impact_defor_Lower")
Impact_defor <- b
# This represents how much forest would have been lost in the years following protection if the offsets had not been protected (the counterfactual scenario)  


#|# Also for the robust se version

ANK_forest_robust_se <- tidy_data7(ANK_dat, modelANK_robust_se, 11)            
CZ_forest_robust_se <- tidy_data7(CZ_dat, modelCZ_robust_se, 9)
TTF_forest_robust_se <- tidy_data7(TTF_dat, modelTTF_robust_se, 14)


ANK_forest_robust_se[10,] <- 0                                                    # Create new row for sums
ANK_forest_robust_se[10, c(1,3:8)] <- apply(ANK_forest_robust_se[, c(1,3:8)], 2, sum)       # Calculate total observed, counterfactual and avoided deforestation following protection.

CZ_forest_robust_se[12,] <- 0
CZ_forest_robust_se[12, c(1,3:8)] <- apply(CZ_forest_robust_se[, c(1,3:8)], 2, sum)

TTF_forest_robust_se[7,] <- 0
TTF_forest_robust_se[7, c(1,3:8)] <- apply(TTF_forest_robust_se[, c(1,3:8)], 2, sum)

#|# Also for the robust se version
b <- rbind(ANK_forest_robust_se[10,], CZ_forest_robust_se[12,], TTF_forest_robust_se[7,])
b$Offset <- c("ANK", "CZ", "TTF")
b <- b[, c(9, 1, 3, 5, 7, 4, 6, 8)]                                   # Re-arrange columns
b[, 6:8] <- b[, 6:8]*-1                                               # Multiply by -1 to turn avoided deforestation columns into 'impact on deforestation' which is -ve if the offset has reduced deforestation and positive if it has increased it.
names(b) <- c("Offset", "Observed_defor", "Counterfactual_Defor", "Counterfactual_Defor_Upper", "Counterfactual_Defor_Lower", "Impact_defor", "Impact_defor_Upper", "Impact_defor_Lower")
Impact_defor_robust_se <- b


#|# Also for the log(y) version

ANK_forest_log.y <- tidy_data7(ANK_dat, modelANK_log.y, 11)            
CZ_forest_log.y <- tidy_data7(CZ_dat, modelCZ_log.y, 9)
TTF_forest_log.y <- tidy_data7(TTF_dat, modelTTF_log.y, 14)


ANK_forest_log.y[10,] <- 0                                                    # Create new row for sums
ANK_forest_log.y[10, c(1,3:8)] <- apply(ANK_forest_log.y[, c(1,3:8)], 2, sum)       # Calculate total observed, counterfactual and avoided deforestation following protection.

CZ_forest_log.y[12,] <- 0
CZ_forest_log.y[12, c(1,3:8)] <- apply(CZ_forest_log.y[, c(1,3:8)], 2, sum)

TTF_forest_log.y[7,] <- 0
TTF_forest_log.y[7, c(1,3:8)] <- apply(TTF_forest_log.y[, c(1,3:8)], 2, sum)

#|# Also for the robust se version
b <- rbind(ANK_forest_log.y[10,], CZ_forest_log.y[12,], TTF_forest_log.y[7,])
b$Offset <- c("ANK", "CZ", "TTF")
b <- b[, c(9, 1, 3, 5, 7, 4, 6, 8)]                                   # Re-arrange columns
b[, 6:8] <- b[, 6:8]*-1                                               # Multiply by -1 to turn avoided deforestation columns into 'impact on deforestation' which is -ve if the offset has reduced deforestation and positive if it has increased it.
names(b) <- c("Offset", "Observed_defor", "Counterfactual_Defor", "Counterfactual_Defor_Upper", "Counterfactual_Defor_Lower", "Impact_defor", "Impact_defor_Upper", "Impact_defor_Lower")
Impact_defor_log.y <- b




#|# a.2) The same for Inverse Hyperbolic Sine

tidy_data7_ihs <- function(data, model, y) {                     
  data <- data[ ,4:22]                                       # Extract only Forest Loss per Year columns
  data <- data.frame(t(data))                                # Transpose
  names(data) <- "Tree_loss"
  data$Tree_loss <- data$Tree_loss/10000                     # Convert tree loss from m2 to hectares
  data$Year <- 1:19
  data <- data[data$Year >= y,]                              # Remove years before protection of the offset
  rownames(data) <- 1:nrow(data)
  
  # Calculate counterfactual deforestation using the inverse hyperbolic sine transformation
  if (length(coef(model)) > 1){
    b      <- coef(model)[4]
    b_upr  <- confint(model)[4,2]
    b_lwr  <- confint(model)[4,1]
  } else if (length(coef(model)) == 1){
    b      <- coef(model)
    b_upr  <- confint(model)[,2]
    b_lwr  <- confint(model)[,1]
  } 
  
  
  # Undo the asinh transformation for the observed tree loss, then subtract the treatment effect
  data$counterfactual_defor <- sinh(asinh(data$Tree_loss) - b)
  data$avoided_defor        <- data$counterfactual_defor - data$Tree_loss
  
  data$counterfactual_upr   <- sinh(asinh(data$Tree_loss) - b_upr)
  data$avoided_defor_upr    <- data$counterfactual_upr - data$Tree_loss
  
  data$counterfactual_lwr   <- sinh(asinh(data$Tree_loss) - b_lwr)
  data$avoided_defor_lwr    <- data$counterfactual_lwr - data$Tree_loss
  
  return(data)
}

ANK_forest_ihs <- tidy_data7_ihs(ANK_dat, modelANK_ihs, 11)            
CZ_forest_ihs <- tidy_data7_ihs(CZ_dat, modelCZ_ihs, 9)
TTF_forest_ihs <- tidy_data7_ihs(TTF_dat, modelTTF_ihs, 14)

ANK_forest_ihs[10,] <- 0                                                    # Create new row for sums
ANK_forest_ihs[10, c(1,3:8)] <- apply(ANK_forest_ihs[, c(1,3:8)], 2, sum)       # Calculate total observed, counterfactual and avoided deforestation following protection.

CZ_forest_ihs[12,] <- 0
CZ_forest_ihs[12, c(1,3:8)] <- apply(CZ_forest_ihs[, c(1,3:8)], 2, sum)

TTF_forest_ihs[7,] <- 0
TTF_forest_ihs[7, c(1,3:8)] <- apply(TTF_forest_ihs[, c(1,3:8)], 2, sum)

#|# b.2) Combine above datasets into one for plotting                                                                        

b <- rbind(ANK_forest_ihs[10,], CZ_forest_ihs[12,], TTF_forest_ihs[7,])
b$Offset <- c("ANK", "CZ", "TTF")
b <- b[, c(9, 1, 3, 5, 7, 4, 6, 8)]                                   # Re-arrange columns
b[, 6:8] <- b[, 6:8]*-1                                               # Multiply by -1 to turn avoided deforestation columns into 'impact on deforestation' which is -ve if the offset has reduced deforestation and positive if it has increased it.
names(b) <- c("Offset", "Observed_defor", "Counterfactual_Defor", "Counterfactual_Defor_Upper", "Counterfactual_Defor_Lower", "Impact_defor", "Impact_defor_Upper", "Impact_defor_Lower")
Impact_defor_ihs <- b



#|# a.4) The same for Poisson regression
#|
tidy_data7_poisson <- function(data, model, y){                     
  data <- data[, 4:22]                                         # Extract only Annual Deforestation columns, excluding Year 0.
  data <- data.frame(t(data))                                  # Transpose
  names(data) <- "Tree_loss"
  data$Tree_loss <- data$Tree_loss / 10000                     # Convert tree loss from m2 to hectares
  data$Year <- 1:19
  data <- data[data$Year >= y,]                                # Remove years before protection (offset)
  rownames(data) <- 1:nrow(data)
  
  # Treatment effect from the Poisson regression (coefficient of interest)
  if (length(coef(model)) > 1){
    b      <- coef(model)[4]
    b_upr  <- confint(model)[4,2]
    b_lwr  <- confint(model)[4,1]
  } else if (length(coef(model)) == 1){
    b      <- coef(model)
    b_upr  <- confint(model)[,2]
    b_lwr  <- confint(model)[,1]
  }
  
  # Calculate counterfactual deforestation on the absolute (hectare) scale
  data$counterfactual_defor <- data$Tree_loss / exp(b)
  data$avoided_defor        <- data$counterfactual_defor - data$Tree_loss
  
  # Calculate upper and lower confidence bounds
  data$counterfactual_upr <- data$Tree_loss / exp(b_upr)
  data$avoided_defor_upr  <- data$counterfactual_upr - data$Tree_loss
  
  data$counterfactual_lwr <- data$Tree_loss / exp(b_lwr)
  data$avoided_defor_lwr  <- data$counterfactual_lwr - data$Tree_loss
  
  return(data)
}

ANK_forest_poisson <- tidy_data7_poisson(ANK_dat, modelANK_poisson, 11)            
CZ_forest_poisson <- tidy_data7_poisson(CZ_dat, modelCZ_poisson, 9)
TTF_forest_poisson <- tidy_data7_poisson(TTF_dat, modelTTF_poisson, 14)

ANK_forest_poisson[10,] <- 0                                                    # Create new row for sums
ANK_forest_poisson[10, c(1,3:8)] <- apply(ANK_forest_poisson[, c(1,3:8)], 2, sum)       # Calculate total observed, counterfactual and avoided deforestation following protection.

CZ_forest_poisson[12,] <- 0
CZ_forest_poisson[12, c(1,3:8)] <- apply(CZ_forest_poisson[, c(1,3:8)], 2, sum)

TTF_forest_poisson[7,] <- 0
TTF_forest_poisson[7, c(1,3:8)] <- apply(TTF_forest_poisson[, c(1,3:8)], 2, sum)

#|# b.4) Combine above datasets into one for plotting                                                                        

b <- rbind(ANK_forest_poisson[10,], CZ_forest_poisson[12,], TTF_forest_poisson[7,])
b$Offset <- c("ANK", "CZ", "TTF")
b <- b[, c(9, 1, 3, 5, 7, 4, 6, 8)]                                   # Re-arrange columns
b[, 6:8] <- b[, 6:8]*-1                                               # Multiply by -1 to turn avoided deforestation columns into 'impact on deforestation' which is -ve if the offset has reduced deforestation and positive if it has increased it.
names(b) <- c("Offset", "Observed_defor", "Counterfactual_Defor", "Counterfactual_Defor_Upper", "Counterfactual_Defor_Lower", "Impact_defor", "Impact_defor_Upper", "Impact_defor_Lower")
Impact_defor_poisson <- b



#--------------------------------7) Fixed Effects Panel regression ----------------------------------------------#

# Second outcome regression. This allows us to estimate the effect of protection across the entire offset portfolio,
# controlling for site and year fixed effects. This helps to control for any unobserved bias.   


                                      # a) Data Construction 

ANK_FE_dat <- annual_defor_ANK
CFAM_FE_dat <- annual_defor_CFAM
CZ_FE_dat <- annual_defor_CZ
TTF_FE_dat <- annual_defor_TTF

# Pool data. 

levels(ANK_FE_dat$Sample)[levels(ANK_FE_dat$Sample) == "Cont"] <- "ANK1"
levels(CFAM_FE_dat$Sample)[levels(CFAM_FE_dat$Sample) == "Cont"] <- "CFAM1"
levels(CZ_FE_dat$Sample)[levels(CZ_FE_dat$Sample) == "Cont"] <- "CZ1"
levels(TTF_FE_dat$Sample)[levels(TTF_FE_dat$Sample) == "Cont"] <- "TTF1"

# Change the name of the four matched control samples (Cont) to ensure that when the data are pooled, each control sample is considered a separate site
# In the pooled dataset we have 152 observations, 1 observation per sample (n= 8, 4 offsets and 4 control), per year (n=19)

FE_dat <- rbind(ANK_FE_dat, CFAM_FE_dat, CZ_FE_dat, TTF_FE_dat)
FE_dat$Year <- factor(FE_dat$Year)


# Make one predictor indicating treated status. Tr = 1 for observations from an offset after protection and 0 for
# observations from an offset before protection and the control sample. Cannot use Treated and Time and the interaction because these
# are collinear with the fixed effects. 


FE_dat$Tr <- ifelse(FE_dat$TreatedF== "1" & FE_dat$TimeF== "1",1,0) 


# Check for differences between groups and over time. 

plotmeans(log_annual_defor ~ Sample, data = FE_dat)

plotmeans(log_annual_defor ~ Year, data = FE_dat)

  
    
                                # b) Fixed Effects Panel Regression 

#|# Plot histogram to see how many 0s are in the panel data
hist(FE_dat$Annual_Deforestation,
     breaks = seq(min(FE_dat$Annual_Deforestation),
                  max(FE_dat$Annual_Deforestation) + 1,
                  by = 1),
     # main = "Histogram of Annual Deforestation",
     xlab = "Annual Deforestation",
     col = "lightblue")

#|# log(y+1)
FE_all <- plm(log_annual_defor ~ Tr, #
              data= FE_dat, index = c("Sample", "Year"), model= "within", effect = "twoways")
summary(FE_all)

#|# log(y+1) with robust SEs
FE_all_robust_se <- feols(log_annual_defor ~ Tr | Year + Sample, data = FE_dat, vcov = "hetero")
summary(FE_all_robust_se)

#|# inverse hyperbolic sine with robust SEs
FE_all_ihs <- feols(ihs_annual_defor ~ Tr | Sample + Year, data = FE_dat, vcov = "hetero")
summary(FE_all_ihs)


#|# poisson regressionwith robust SEs
FE_all_poisson<- fepois(Annual_Deforestation ~ Tr | Sample + Year, data = FE_dat, vcov = "hetero")
summary(FE_all_poisson)

#|# log(y) with robust SEs
FE_all_log.y <- feols(log.annual_defor ~ Tr | Year + Sample, data = FE_dat, vcov = "hetero")
summary(FE_all_log.y)


# Tr is the coefficient of interest. A significant negative coefficient indicates a significant reduction in the log-transformed
# count of deforestation following protection of the four biodiversity offsets.

# Back-transform the estimates.

(exp(coef(FE_all))-1)*100              # Results show that protection reduced average annual deforestation by 58% (-37 to -72%) across
(exp(confint(FE_all))-1)*100           # the entire offset portfolio 


                                # c) Tests 

# Compare to simple ols regression to test whether the fixed effects are needed. 

ols2 <- lm(log_annual_defor ~ Tr, data = FE_dat)
summary(ols2)

pFtest(FE_all, ols2)

# p<0.05 so there is significant heterogeneity between groups and over time - the fixed effects are needed.


#                               d) Try also using random effects

# m1<-lmer(log_annual_defor ~ Tr  + (1|Sample)+(1|Year),
#          data= FE_dat)
# summary(m1)
# 
# (exp(-0.7476)-1)*100
# (exp(confint(m1))-1)*100

# Results are very similar, and within the confidence intervals to the results from the FE panel regression. 
# They indicate a 53% (-27 to -69%) reduction in deforestation following protection. 



#                                e) Calculate avoided deforestation across entire offset portfolio

# Using the estimate of treatment effect from the fixed effects panel regression. 

# Read in datasets containing forest cover and annual forest loss values for whole offset area. 
# Remove unwanted columns and convert annual tree loss in m2 to hectares.

tidy_data8 <- function(data){                     
  data <- data[ ,4:22]                                         # Extract only Forest Loss per Year columns, excluding Year 0.
  data <- data.frame(t(data))                                  # Transpose
  names(data) <- "Tree_loss"
  data$Tree_loss <- data$Tree_loss/10000                      # Convert tree loss in m2 to hectares
  data$Year <- 1:19
  rownames(data) <- 1:nrow(data)
  return(data)
}

ANK_dat <- tidy_data8(ANK_dat)
CZ_dat <- tidy_data8(CZ_dat)
CFAM_dat <- tidy_data8(CFAM_dat)
TTF_dat <- tidy_data8(TTF_dat)


# Merge into one dataset

offsets_defor <- cbind(ANK_dat[,c(2,1)], CZ_dat[,1], CFAM_dat[,1], TTF_dat[,1])   
names(offsets_defor) <-  c("Year", "ANK", "CZ", "CFAM", "TTF")               # ANK, CZ, CFAM, TTF columns show the total amount deforestation in the offset each year (in hectares)

# Set annual deforestation before protection of the offset to 0 because we're only interested in avoided deforestation after protection.

offsets_defor$ANK[offsets_defor$Year<11] <- 0                                
offsets_defor$CZ[offsets_defor$Year<9] <- 0                                   
offsets_defor$CFAM[offsets_defor$Year<13] <- 0
offsets_defor$TTF[offsets_defor$Year<14] <- 0

# Remove Years before the first offset was protected

offsets_defor <- offsets_defor[offsets_defor$Year >= 9,]                    
rownames(offsets_defor) <- 1:nrow(offsets_defor)

# Calculate the total deforestation across all protected offsets each year.

offsets_defor$Sum_defor <- apply(offsets_defor[,2:5], 1, sum)   

#|# Create objects for different outcome specifications
offsets_defor_robust_se <- offsets_defor
offsets_defor_ihs <- offsets_defor
offsets_defor_poisson <- offsets_defor
offsets_defor_log.y <- offsets_defor



# Use the estimated treatment effect from the FE panel regression to convert the annual observed deforestation
# to the counterfactual. The difference between these two values is the amount of avoided deforestation. 

offsets_defor$counterfactual_defor <- (1/exp(coef(FE_all)))*offsets_defor$Sum_defor
offsets_defor$avoided_defor <- offsets_defor$counterfactual_defor - offsets_defor$Sum_defor

offsets_defor$counterfactual_defor_upr <- (1/exp(confint(FE_all)[1,2]))*offsets_defor$Sum_defor
offsets_defor$avoided_defor_upr <- offsets_defor$counterfactual_defor_upr - offsets_defor$Sum_defor

offsets_defor$counterfactual_defor_lwr <- (1/exp(confint(FE_all)[1,1]))*offsets_defor$Sum_defor
offsets_defor$avoided_defor_lwr <- offsets_defor$counterfactual_defor_lwr - offsets_defor$Sum_defor

# Create new dataset containing total observed, counterfactual and avoided deforestation across all four offsets following protection. 
# (i.e. the column totals from offsets_defor). This can then be joined to the individual results from the site-based
# difference-in-difference regressions to create a single dataframe for plotting.

sum_overall_defor <- data.frame(matrix(nrow = 1, ncol =8))
names(sum_overall_defor) <- names(Impact_defor)
sum_overall_defor$Offset <- "All"
sum_overall_defor[1,2:8] <- apply(offsets_defor[, c(6,7, 9, 11, 8, 10, 12)], 2, sum)
sum_overall_defor[1, c(6,7,8)] <- sum_overall_defor[1, c(6,7,8)]*-1                   # Multiply by -1 to turn avoided deforestation columns into 'impact on deforestation' 
# which is -ve if the offset has reduced deforestation and positive if it has increased it.


Impact_defor <- rbind(Impact_defor, sum_overall_defor)



#|# The same for log(y+1) with robust SEs
offsets_defor_robust_se$counterfactual_defor <- (1/exp(coef(FE_all_robust_se)))*offsets_defor_robust_se$Sum_defor
offsets_defor_robust_se$avoided_defor <- offsets_defor_robust_se$counterfactual_defor - offsets_defor_robust_se$Sum_defor

offsets_defor_robust_se$counterfactual_defor_upr <- (1/exp(confint(FE_all_robust_se)[1,2]))*offsets_defor_robust_se$Sum_defor
offsets_defor_robust_se$avoided_defor_upr <- offsets_defor_robust_se$counterfactual_defor_upr - offsets_defor_robust_se$Sum_defor

offsets_defor_robust_se$counterfactual_defor_lwr <- (1/exp(confint(FE_all_robust_se)[1,1]))*offsets_defor_robust_se$Sum_defor
offsets_defor_robust_se$avoided_defor_lwr <- offsets_defor_robust_se$counterfactual_defor_lwr - offsets_defor_robust_se$Sum_defor

sum_overall_defor <- data.frame(matrix(nrow = 1, ncol =8))
names(sum_overall_defor) <- names(Impact_defor)
sum_overall_defor$Offset <- "All"
sum_overall_defor[1,2:8] <- apply(offsets_defor_robust_se[, c(6,7, 9, 11, 8, 10, 12)], 2, sum)
sum_overall_defor[1, c(6,7,8)] <- sum_overall_defor[1, c(6,7,8)]*-1                   # Multiply by -1 to turn avoided deforestation columns into 'impact on deforestation' 
# which is -ve if the offset has reduced deforestation and positive if it has increased it.


Impact_defor_robust_se <- rbind(Impact_defor_robust_se, sum_overall_defor)


#|# The same for log(y) with robust SEs
offsets_defor_log.y$counterfactual_defor <- (1/exp(coef(FE_all_log.y)))*offsets_defor_log.y$Sum_defor
offsets_defor_log.y$avoided_defor <- offsets_defor_log.y$counterfactual_defor - offsets_defor_log.y$Sum_defor

offsets_defor_log.y$counterfactual_defor_upr <- (1/exp(confint(FE_all_log.y)[1,2]))*offsets_defor_log.y$Sum_defor
offsets_defor_log.y$avoided_defor_upr <- offsets_defor_log.y$counterfactual_defor_upr - offsets_defor_log.y$Sum_defor

offsets_defor_log.y$counterfactual_defor_lwr <- (1/exp(confint(FE_all_log.y)[1,1]))*offsets_defor_log.y$Sum_defor
offsets_defor_log.y$avoided_defor_lwr <- offsets_defor_log.y$counterfactual_defor_lwr - offsets_defor_log.y$Sum_defor

sum_overall_defor <- data.frame(matrix(nrow = 1, ncol =8))
names(sum_overall_defor) <- names(Impact_defor)
sum_overall_defor$Offset <- "All"
sum_overall_defor[1,2:8] <- apply(offsets_defor_log.y[, c(6,7, 9, 11, 8, 10, 12)], 2, sum)
sum_overall_defor[1, c(6,7,8)] <- sum_overall_defor[1, c(6,7,8)]*-1                   # Multiply by -1 to turn avoided deforestation columns into 'impact on deforestation' 
# which is -ve if the offset has reduced deforestation and positive if it has increased it.


Impact_defor_log.y <- rbind(Impact_defor_log.y, sum_overall_defor)



#|# The same for IHS with robust SEs 
#|
offsets_defor_ihs$counterfactual_defor <- sinh( asinh(offsets_defor_ihs$Sum_defor) - coef(FE_all_ihs) )
offsets_defor_ihs$avoided_defor <- offsets_defor_ihs$counterfactual_defor - offsets_defor_ihs$Sum_defor

offsets_defor_ihs$counterfactual_defor_upr <- sinh( asinh(offsets_defor_ihs$Sum_defor) - confint(FE_all_ihs)[1,2] )
offsets_defor_ihs$avoided_defor_upr <- offsets_defor_ihs$counterfactual_defor_upr - offsets_defor_ihs$Sum_defor

offsets_defor_ihs$counterfactual_defor_lwr <- sinh( asinh(offsets_defor_ihs$Sum_defor) - confint(FE_all_ihs)[1,1] )
offsets_defor_ihs$avoided_defor_lwr <- offsets_defor_ihs$counterfactual_defor_lwr - offsets_defor_ihs$Sum_defor



sum_overall_defor <- data.frame(matrix(nrow = 1, ncol =8))
names(sum_overall_defor) <- names(Impact_defor)
sum_overall_defor$Offset <- "All"
sum_overall_defor[1,2:8] <- apply(offsets_defor_ihs[, c(6,7, 9, 11, 8, 10, 12)], 2, sum)
sum_overall_defor[1, c(6,7,8)] <- sum_overall_defor[1, c(6,7,8)]*-1                   # Multiply by -1 to turn avoided deforestation columns into 'impact on deforestation' 
# which is -ve if the offset has reduced deforestation and positive if it has increased it.


Impact_defor_ihs <- rbind(Impact_defor_ihs, sum_overall_defor)#




#|# The same for Poisson with robust SEs 
offsets_defor_poisson$counterfactual_defor <- offsets_defor_poisson$Sum_defor / exp(coef(FE_all_poisson))
offsets_defor_poisson$avoided_defor <- offsets_defor_poisson$counterfactual_defor - offsets_defor_poisson$Sum_defor

offsets_defor_poisson$counterfactual_defor_upr <- offsets_defor_poisson$Sum_defor / exp(confint(FE_all_poisson)[1,2])
offsets_defor_poisson$avoided_defor_upr <- offsets_defor_poisson$counterfactual_defor_upr - offsets_defor_poisson$Sum_defor

offsets_defor_poisson$counterfactual_defor_lwr <- offsets_defor_poisson$Sum_defor / exp(confint(FE_all_poisson)[1,1])
offsets_defor_poisson$avoided_defor_lwr <- offsets_defor_poisson$counterfactual_defor_lwr - offsets_defor_poisson$Sum_defor



sum_overall_defor <- data.frame(matrix(nrow = 1, ncol =8))
names(sum_overall_defor) <- names(Impact_defor)
sum_overall_defor$Offset <- "All"
sum_overall_defor[1,2:8] <- apply(offsets_defor_poisson[, c(6,7, 9, 11, 8, 10, 12)], 2, sum)
sum_overall_defor[1, c(6,7,8)] <- sum_overall_defor[1, c(6,7,8)]*-1                   # Multiply by -1 to turn avoided deforestation columns into 'impact on deforestation' 
# which is -ve if the offset has reduced deforestation and positive if it has increased it.


Impact_defor_poisson <- rbind(Impact_defor_poisson , sum_overall_defor)#




                              # Convert treatment effect from FE estimate to Cohen's D for comparison with Borner et al

# Cohen's d effect size is (mean of treated group - mean of control group)/ standard deviation of pooled data
# However, Borner et al use the standard deviation of the control group so we will use that instead. 


test <- offsets_defor

# Don't actually need to do deforestation as a % of area but will keep because could be useful

test$counterfactual_defor_perc <- 0
test <- test[ ,c(1:7,13,8:12)]

test$avoid_defor_perc <- 0
test <- test[ ,c(1:9,14,10:13)]

test$Area <- 0
test$Area[c(1,2)] <- 3787                   # Area of CZ
test$Area[c(3,4)] <- 3787 + 6904            # Area of CZ + ANK
test$Area[5] <- (3787 + 6904 + 9423)        # Area of CZ + ANK + CFAM
test$Area[6:11] <- (3787 + 6904 + 9423 + 8626)     # Total area of offsets

test$counterfactual_defor_perc <- (test$counterfactual_defor/test$Area)*100

test$avoid_defor_perc <- (test$avoided_defor/test$Area)*100

# Calculation of Cohen's d:
# Treated sample is the observed annual deforestation across the entire offset portfolio following protection
# Control sample is the counterfactual annual deforestation for the entire offset portfolio following protection

# We are using this instead of deforestation in the matched control sample because because my estimates of treatment effect (the % difference and the hectares of avoided deforestation) 
# weren’t derived from the matched control sample itself but the estimated counterfactual (mean annual defor in the matched control adjusted to account for the pre-intervention differences between groups – the difference-in-differences). 


(mean(test$Sum_defor) - mean(test$counterfactual_defor))/sd(test$counterfactual_defor)

# Cohen's d = - 0.51. 

# Test using the standard deviation of the pooled data:

x <- matrix(nrow = 22, ncol = 1)
x <- c(test$Sum_defor, test$counterfactual_defor)

sd(x)

(mean(test$Sum_defor) - mean(test$counterfactual_defor))/sd(x)

# Repeat for results from each individual offset

ANK_forest2 <- ANK_forest[1:9,]


(mean(ANK_forest2$Tree_loss) - mean(ANK_forest2$counterfactual_defor))/sd(ANK_forest2$counterfactual_defor)

# Cohen's D for Ankerana = -1.03

CZ_forest2 <- CZ_forest[1:11,]

(mean(CZ_forest2$Tree_loss) - mean(CZ_forest2$counterfactual_defor))/sd(CZ_forest2$counterfactual_defor)


# Cohen's D for Conservation Zone is -0.63
 
TTF_forest2 <- TTF_forest[1:6,]

(mean(TTF_forest2$Tree_loss) - mean(TTF_forest2$counterfactual_defor))/sd(TTF_forest2$counterfactual_defor)

# Cohen's D for Torotorofotsy is 1.29



################################################################################################################

# ------------------------- 5) put all results of fixed effect and of pooled reg.

################################################################################################################

#|# get package here: https://osf.io/nxwvd
#|# install.packages("i4results_0.1.0.tar.gz", repos = NULL, type = "source")
# library(i4results) #load libraries
library(openxlsx)


#|# Flexible version of i4results_flexible allowing t-value, t value, z value and Pr(>|t|), Pr(>|z|), and fixest models
i4results_flexible <- function(original_model, robustness_models, out = NULL, append = FALSE) {
  if (!is.null(out) && !requireNamespace("openxlsx", quietly = TRUE)) {
    stop("Package 'openxlsx' needed for Excel import/export. Please install it.")
  }
  if (length(robustness_models) > 5) {
    stop("You can supply up to 5 robustness models only.")
  }
  
  # Extract original model info
  o_cmdline <- paste(deparse(original_model$call), collapse = " ")
  o_n       <- stats::nobs(original_model)
  
  # Summary and coefficients extraction, including fixest
  is_fixest_o <- inherits(original_model, "fixest")
  if (is_fixest_o) {
    summary_o <- summary(original_model)
    coefs_o   <- as.data.frame(summary_o$coeftable)
  } else {
    summary_o <- summary(original_model)
    coefs_o   <- as.data.frame(summary_o$coefficients)
  }
  
  # Detect t/z and p-value column names
  t_o_col <- grep("t[- ]?value|z value", names(coefs_o), ignore.case = TRUE, value = TRUE)[1]
  p_o_col <- grep("Pr\\(>\\|t\\|\\)|Pr\\(>\\|z\\|\\)", names(coefs_o), value = TRUE)[1]
  
  # Remove intercept
  coefs_o <- coefs_o[!rownames(coefs_o) %in% "(Intercept)", , drop = FALSE]
  
  # Confidence intervals
  ci_o <- as.data.frame(confint(original_model, level = 0.95))
  colnames(ci_o) <- c("lower", "upper")
  ci_o <- ci_o[!rownames(ci_o) %in% "(Intercept)", , drop = FALSE]
  
  paramlist_all <- rownames(coefs_o)
  results_list  <- list()
  
  # Loop over robustness models
  for (rep_name in names(robustness_models)) {
    rep_model <- robustness_models[[rep_name]]
    r_cmdline <- paste(deparse(rep_model$call), collapse = " ")
    r_n       <- stats::nobs(rep_model)
    
    # Summary and coefficients extraction for robustness, including fixest
    is_fixest_r <- inherits(rep_model, "fixest")
    if (is_fixest_r) {
      summary_r <- summary(rep_model)
      coefs_r   <- as.data.frame(summary_r$coeftable)
    } else {
      summary_r <- summary(rep_model)
      coefs_r   <- as.data.frame(summary_r$coefficients)
    }
    
    # Detect t/z and p-value column names for robustness
    t_r_col <- grep("t[- ]?value|z value", names(coefs_r), ignore.case = TRUE, value = TRUE)[1]
    p_r_col <- grep("Pr\\(>\\|t\\|\\)|Pr\\(>\\|z\\|\\)", names(coefs_r), value = TRUE)[1]
    
    # Remove intercept for robustness
    coefs_r <- coefs_r[!rownames(coefs_r) %in% "(Intercept)", , drop = FALSE]
    
    # Confidence intervals for robustness
    ci_r <- as.data.frame(confint(rep_model, level = 0.95))
    colnames(ci_r) <- c("lower", "upper")
    ci_r <- ci_r[!rownames(ci_r) %in% "(Intercept)", , drop = FALSE]
    
    # Combine per parameter
    for (p in paramlist_all) {
      o_coeff    <- coefs_o[p, "Estimate"]
      o_std_err  <- coefs_o[p, "Std. Error"]
      o_t        <- coefs_o[p, t_o_col]
      o_p_val    <- coefs_o[p, p_o_col]
      o_ci_lower <- ci_o[p, "lower"]
      o_ci_upper <- ci_o[p, "upper"]
      
      if (p %in% rownames(coefs_r)) {
        r_coeff    <- coefs_r[p, "Estimate"]
        r_std_err  <- coefs_r[p, "Std. Error"]
        r_t        <- coefs_r[p, t_r_col]
        r_p_val    <- coefs_r[p, p_r_col]
        r_ci_lower <- ci_r[p, "lower"]
        r_ci_upper <- ci_r[p, "upper"]
      } else {
        r_coeff <- r_std_err <- r_t <- r_p_val <- r_ci_lower <- r_ci_upper <- NA
      }
      
      record <- data.frame(
        paramname  = p,
        study      = rep_name,
        o_cmdline  = substr(o_cmdline, 1, 244),
        r_cmdline  = substr(r_cmdline, 1, 244),
        o_n        = o_n,
        r_n        = r_n,
        o_coeff    = round(o_coeff, 3),
        o_std_err  = round(o_std_err, 3),
        o_t        = round(o_t, 3),
        o_p_val    = round(o_p_val, 3),
        o_ci_lower = round(o_ci_lower, 3),
        o_ci_upper = round(o_ci_upper, 3),
        r_coeff    = round(r_coeff, 3),
        r_std_err  = round(r_std_err, 3),
        r_t        = round(r_t, 3),
        r_p_val    = round(r_p_val, 3),
        r_ci_lower = round(r_ci_lower, 3),
        r_ci_upper = round(r_ci_upper, 3),
        stringsAsFactors = FALSE
      )
      
      results_list[[length(results_list) + 1]] <- record
    }
  }
  
  results_df <- do.call(rbind, results_list)
  
  # Output or export
  if (!is.null(out)) {
    if (append && file.exists(out)) {
      oldres   <- openxlsx::read.xlsx(out)
      combined <- rbind(oldres, results_df)
      openxlsx::write.xlsx(combined, file = out, asTable = FALSE)
      message("Appended results to Excel file: ", out)
    } else {
      openxlsx::write.xlsx(results_df, file = out, asTable = FALSE)
      if (append) message("File did not exist, so created a new Excel file: ", out)
      else       message("Exported results to Excel file: ", out)
    }
  } else {
    message("Data with original & robustness comparisons returned as a data frame.")
  }
  
  invisible(results_df)
}


#|#------------- a. 3 offsets:
robustness_ANK <- 
  list(
    "log(y+1), robust SEs" = modelANK_robust_se
    ,"ihs(y), robust SEs" = modelANK_ihs
    ,"Poisson, robust SEs" = modelANK_poisson
    ,"log(y), robust SEs" = modelANK_log.y
) # create a list of robustness checks

robustness_CZ <- 
  list(
    "log(y+1), robust SEs" = modelCZ_robust_se
    ,"ihs(y), robust SEs" = modelCZ_ihs
    ,"Poisson, robust SEs" = modelCZ_poisson
    ,"log(y), robust SEs" = modelCZ_log.y
  ) # create a list of robustness checks

robustness_TTF <- 
  list(
    "log(y+1), robust SEs" = modelTTF_robust_se
    ,"ihs(y), robust SEs" = modelTTF_ihs
    ,"Poisson, robust SEs" = modelTTF_poisson
    ,"log(y), robust SEs" = modelTTF_log.y
  ) # create a list of robustness checks


#|# export
i4results_flexible(modelANK, robustness_ANK, out = "4-Outcome_transformation/ANKresults.xlsx") # put 1. the original model, 2. the list of robusteness checks (up to 5), and where to save it
i4results_flexible(modelCZ, robustness_CZ, out = "4-Outcome_transformation/CZresults.xlsx")
i4results_flexible(modelTTF, robustness_TTF, out = "4-Outcome_transformation/TTFresults.xlsx")

#------------- b. pooled FE:
#check i4results function; # not working for plm() regression!
# View(i4results) # modify the "t value" to "t-value"


#|# export
robustness_PooledFE <- 
  list(
    "log(y+1), robust SEs" = FE_all_robust_se
    ,"ihs(y), robust SEs" = FE_all_ihs
    ,"Poisson, robust SEs" = FE_all_poisson
    ,"log(y), robust SEs" = FE_all_log.y
  ) # create a list of robustness checks


i4results_flexible(FE_all, robustness_PooledFE, out = "4-Outcome_transformation/Pooled_FE.xlsx")



