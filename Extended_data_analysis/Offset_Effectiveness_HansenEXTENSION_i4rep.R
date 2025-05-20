# This code is Solène's extension of Hansen dataset to 2023 and "changes" to the original script.
#|# <- comments by Solène are marked like this


# Original code is used in the publication entitled "Madagascar's biggest mine is set to achieve No Net Loss of Forest"
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
library("readr") #|# added
library("lme4") #|# added

# Code Structure

#|# A) ORIGINAL 
# 1) Data Construction for Matching
# 2) Matching
# 3) Extract Matched Pairs from matching output
# 4) Data construction for DiD regression using matched pairs
# 5) DiD Regression 
# 6) Fixed Effects Panel Regression
#|#########################################################
#|# B) EXTENSION
# 7) Data Construction for Matching
# 8) Matching
# 9) Extract Matched Pairs from matching output
# 10) Data construction for DiD regression using matched pairs
# 11) DiD Regression 
# 12) Fixed Effects Panel Regression
#|########################################################
#|# 13) EXPORTING

setwd('D:/Solene/replication_games/Biodiversity_offset_effectiveness-main') #|# added



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

TTF <- tidy_data("Input_data/Final_covariates/Sample_TTF3.dbf", TTF, 1, "TTF")           
ANK <- tidy_data("Input_data/Final_covariates/Sample_ANK3.dbf", ANK, 1, "ANK")
CFAM <- tidy_data("Input_data/Final_covariates/Sample_CFAM3.dbf", CFAM, 1, "CFAM")
CZ <- tidy_data("Input_data/Final_covariates/Sample_CZ3.dbf", CZ, 1, "CZ")

# TTF = Torotorofotsy
# ANK = Ankerana
# CFAM = Corridor Forestier Analamay-Mantadia
# CZ = Conservation Zone


# Load Control separately because columns are different # 

Control = read.dbf("Input_data/Final_covariates/Final_control.dbf")
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
  data <- data %>% mutate(across(everything(), ~ifelse(. == -9999, NA, .))) # changed from data <- data %>% na_if(-9999)
  
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

# log(y+1) transformation of outcome variable required because non-normal properties of count data violate assumptions of 
# homoscedascity of linear models.

annual_defor_ANK$log_annual_defor <- log(annual_defor_ANK$Annual_Deforestation + 1)
annual_defor_CZ$log_annual_defor <- log(annual_defor_CZ$Annual_Deforestation + 1)
annual_defor_CFAM$log_annual_defor <- log(annual_defor_CFAM$Annual_Deforestation +1)
annual_defor_TTF$log_annual_defor <- log(annual_defor_TTF$Annual_Deforestation +1)


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

modelANK <- lm(log_annual_defor ~ TreatedF*TimeF, data= annual_defor_ANK)
summary(modelANK) 

# Back-transform the estimate to get the treatment effect - the  percentage difference in average annual deforestation between
# the offset and the estimated counterfactual following protection. 
# The estimated counterfactual is the average annual deforestation in the matched control sample after the intervention, adjusted to 
# to account for pre-intervention differences between the two samples. 

exp(coef(modelANK)[4])-1
exp(confint(modelANK)[4,])-1

# CZ # 

modelCZ <- lm(log_annual_defor ~ TreatedF*TimeF, data= annual_defor_CZ)
summary(modelCZ)

exp(coef(modelCZ)[4])-1
exp(confint(modelCZ)[4,])-1

# TTF #

modelTTF <- lm(log_annual_defor ~ TreatedF*TimeF, data= annual_defor_TTF)
summary(modelTTF)

exp(coef(modelTTF)[4])-1
exp(confint(modelTTF)[4,])-1

# Results show a significant reduction in average annual deforestation of 96% (- 89 to - 98%) in Ankerana
# and 66% (-27 to -84%) in the COnservation Zone following protection. 

# In Torotorofotsy, protection had no significant effect on deforestation. 

DiD_res <- rbind(t(tidy(modelANK)), t(tidy(modelCZ)), t(tidy(modelTTF)))



#                               e) Extract back-transformed coefficients and confidence intervals from DiD models


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



#--------------------------------6) Fixed Effects Panel regression ----------------------------------------------#

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

FE_all <- plm(log_annual_defor ~ Tr, 
              data= FE_dat, index = c("Sample", "Year"), model= "within", effect = "twoways")
summary(FE_all)


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

m1<-lmer(log_annual_defor ~ Tr  + (1|Sample)+(1|Year),
         data= FE_dat)
summary(m1)

(exp(-0.7476)-1)*100
(exp(confint(m1))-1)*100

# Results are very similar, and within the confidence intervals to the results from the FE panel regression. 
# They indicate a 53% (-27 to -69%) reduction in deforestation following protection. 







################################################################################################################

#---------------------------B) EXTENSION: MODELS USING HANSEN et al. 2013 UNTIL 2019 ------------------------------------# #|# same section, see 4b for extension

################################################################################################################

# --------------------- 7) Data Construction for Matching: same for both outcome -----------------------------------------

# Function to read in input data for each offset, remove unnecessary columns and add column to indicate treatment status
# and offset of origin. #|# changed from original!

df <- read_csv("devenish_2022_allsites_controls_jrc_def_deg_gfw_def.csv") #|# added

# Rename columns #
df <- df %>% dplyr::select(ANNUAL_RAI, ASPECT_MAD, DEM_MADA_R, DIST_DEFOR, DIST_EDGE, 
                           DIST_RIVER, DIST_ROUTE, DIST_SETT, DIST_TRACK, GFW_defory, 
                           JRC_defory, JRC_degrad, POPDENS_MA, SLOPE_MADA, 
                           SOURCE, SOURCE_ID, TREELOSS19, UNIQUE_ID)

# Renaming columns to match the order of 'cols'
colnames(df) <- c("Annual_Rain", "Aspect", "Elevation", "Dist_defor", "Dist_edge", 
                  "Dist_river", "Dist_road", "Dist_sett", "Dist_track", "GFW_defory", 
                  "JRC_defory", "JRC_degrad", "Pop_density", "Slope", 
                  "Source", "Source_ID", "Tree_loss", "Unique_ID")

# original names
#cols <- c("Tree_loss", "Pop_density", "Dist_sett", "Slope", "Elevation", "Aspect", "Annual_Rain",
# "Dist_track", "Dist_road", "Dist_river", "Dist_edge", "Dist_defor", "treated", "offset")

# clean Hansen dataset:
df <- df %>% dplyr::mutate(
  GFW_defory = ifelse(GFW_defory == -1, NA, GFW_defory-2000)
)


tidy_data <-function(name, number, label){
  name <- df %>% filter(Source == name)
  #|# name = subset(name, select= -c(1, 16))
  name$treated = number              # Number = 1 for pixels from an offset and 0 for control pixels. 
  name$offset = label
  return(name)
}

TTF <- tidy_data('TTF', 1, "TTF")           
ANK <- tidy_data("ANK", 1, "ANK")
CFAM <- tidy_data("CFAM", 1, "CFAM")
CZ <- tidy_data("CZ", 1, "CZ")

# TTF = Torotorofotsy
# ANK = Ankerana
# CFAM = Corridor Forestier Analamay-Mantadia
# CZ = Conservation Zone


# Load Control separately because columns are different # 

Control <- df %>% filter(Source == "CONTROL")
Control$treated = 0
Control$offset = "Cont"
#Control <- subset(Control, select = -c(4))        # Remove unwanted fire variable #|# unnecessary in our case
#Control <- Control[,c(2,3,1,4:16)]                # Re-order columns to match offset dataframes #|# unnecessary in our case


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
  data <- data %>% mutate(across(Annual_Rain:Dist_track, ~ifelse(. == -9999, NA, .))) #|# changed from data <- data %>% na_if(-9999)
  
  # We keep observations without na
  data <- data %>% drop_na(Annual_Rain:Dist_track)
  
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


#----------------------------8) The Matching -----------------------------------------#

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


#---------------------------9) Extract matched pairs from matching output ------------------------------------#

# Extract matched pairs from input data


m.data.ANK <- match.data(m.out.ANK, group = "all", distance = "distance", weights = "weights", subclass = "subclass",
                         data = ANKCONT, drop.unmatched = TRUE)

m.data.CZ <- match.data(m.out.CZ, group = "all", distance = "distance", weights = "weights", subclass = "subclass",
                        data = CZCONT, drop.unmatched = TRUE)

m.data.CFAM <- match.data(m.out.CFAM, group = "all", distance = "distance", weights = "weights", subclass = "subclass",
                          data = CFAMCONT, drop.unmatched = TRUE)

m.data.TTF <- match.data(m.out.TTF, group = "all", distance = "distance", weights = "weights", subclass = "subclass",
                         data = TTFCONT, drop.unmatched = TRUE)



#---------------------------10) Data construction for DiD regression using matched pairs ------------------------------


# Aggregate pixels into treated (offset) and control samples and tabulate observations within each sample by tree loss year. 
# This gives the count of pixels within each sample deforested each year. 

annual_ext_defor_ANK_ext <- data.frame(table(factor(m.data.ANK$GFW_defory, levels = 1:23), m.data.ANK$offset))

# levels = 1:23 removes observations with 0 value for Tree Loss year (which were not deforested over the study period).
# This is because we are interested in comparing deforestation outcomes between offsets and the matched controls. 

annual_ext_defor_CZ_ext <- data.frame(table(factor(m.data.CZ$GFW_defory, levels = 1:23), m.data.CZ$offset))
annual_ext_defor_CZ_ext$Var2 <- factor(annual_ext_defor_CZ_ext$Var2, levels = c("CZ", "Cont"))    # For the plots the order of the factors needs to match the other offsets

annual_ext_defor_CFAM_ext <- data.frame(table(factor(m.data.CFAM$GFW_defory, levels = 1:23), m.data.CFAM$offset))

annual_ext_defor_TTF_ext <- data.frame(table(factor(m.data.TTF$GFW_defory, levels = 1:23), m.data.TTF$offset))
annual_ext_defor_TTF_ext$Var2 <- factor(annual_ext_defor_TTF_ext$Var2, levels = c("TTF", "Cont"))

label <- c("Year", "Sample", "annual_ext_Deforestation")
names(annual_ext_defor_ANK_ext) <- paste0(label)
names(annual_ext_defor_CFAM_ext) <- paste0(label)
names(annual_ext_defor_CZ_ext) <- paste0(label)
names(annual_ext_defor_TTF_ext) <- paste0(label)

# But offsets are different sizes - need to calculate annual_ext deforestation as a percentage of total pixels in the sample to plot

# Calculating Percentage annual_ext Deforestation with 1:1 matching.

annual_ext_defor_ANK_ext$Perc_annual_ext_Defor <- (annual_ext_defor_ANK_ext$annual_ext_Deforestation/(nrow(m.data.ANK)/2))*100          # In the matched dataset the no. of control pixels = No. of treatment pixels because I matched 1:1.  
annual_ext_defor_CFAM_ext$Perc_annual_ext_Defor <- (annual_ext_defor_CFAM_ext$annual_ext_Deforestation/(nrow(m.data.CFAM)/2))*100       # Therefore, nrow(m.data.CFAM_ext)/2 is the total number of treatment and the total number of control pixels 
annual_ext_defor_CZ_ext$Perc_annual_ext_Defor <- (annual_ext_defor_CZ_ext$annual_ext_Deforestation/(nrow(m.data.CZ)/2))*100             # in the matched dataset. 
annual_ext_defor_TTF_ext$Perc_annual_ext_Defor <- (annual_ext_defor_TTF_ext$annual_ext_Deforestation/(nrow(m.data.TTF)/2))*100



#------------------------------11) DiD Regression  -----------------------------------------#

# a) Data Construction 


Data_construction_DiD <- function(offset, y){    # y corresponds to the year of protection - so Time = 0 before protection and 1 after protection. 
  offset$Year <- as.numeric(rep(1:23,2))        # Have to make Year numeric for >= to work
  offset$TimeF <- factor(ifelse(offset$Year >= y, 1,0))
  offset$TreatedF <- factor(ifelse(offset$Sample != "Cont", 1,0))
  return(offset)
}

annual_ext_defor_ANK_ext <- Data_construction_DiD(annual_ext_defor_ANK_ext, 11)
annual_ext_defor_CFAM_ext <- Data_construction_DiD(annual_ext_defor_CFAM_ext, 13)
annual_ext_defor_CZ_ext <- Data_construction_DiD(annual_ext_defor_CZ_ext, 9)
annual_ext_defor_TTF_ext <- Data_construction_DiD(annual_ext_defor_TTF_ext, 14)


# b) Outcome variable transformation

# log(y+1) transformation of outcome variable required because non-normal properties of count data violate assumptions of 
# homoscedascity of linear models.

annual_ext_defor_ANK_ext$log_annual_ext_defor <- log(annual_ext_defor_ANK_ext$annual_ext_Deforestation + 1)
annual_ext_defor_CZ_ext$log_annual_ext_defor <- log(annual_ext_defor_CZ_ext$annual_ext_Deforestation + 1)
annual_ext_defor_CFAM_ext$log_annual_ext_defor <- log(annual_ext_defor_CFAM_ext$annual_ext_Deforestation +1)
annual_ext_defor_TTF_ext$log_annual_ext_defor <- log(annual_ext_defor_TTF_ext$annual_ext_Deforestation +1)


# c) Test for parallel trends

# Parallel trends in outcomes between treated and control samples in the years before the intervention
# is a key assumption of difference-in-differences regressions.

# Use only data from the years before the offsets were protected #

ANK_ext_data_before <- annual_ext_defor_ANK_ext[(annual_ext_defor_ANK_ext$Year <11),]
CFAM_ext_data_before <- annual_ext_defor_CFAM_ext[(annual_ext_defor_CFAM_ext$Year <13),]
CZ_ext_data_before <- annual_ext_defor_CZ_ext[(annual_ext_defor_CZ_ext$Year <9),]
TTF_ext_data_before <- annual_ext_defor_TTF_ext[(annual_ext_defor_TTF_ext$Year <14),]


# ANK_ext #

ANK_exta <- lm(log_annual_ext_defor ~ Year*TreatedF, data= ANK_ext_data_before)
summary(ANK_exta)       # If the interaction between Year and TreatedF is not significant, there is no 
# significant difference in the relationship between Year and the log-transformed
# count of deforestation between treated and control samples --> parallel trends assumption holds.

#|# NEW: ANK_ext still holds

# CFAM_ext #

CFAM_exta <- lm(log_annual_ext_defor ~ Year*TreatedF, data= CFAM_ext_data_before)
summary(CFAM_exta)

# There is a significant difference in the trend in deforestation over time between treated and control samples 
# No parallel trends. 

#|# NEW: CFAM_ext significant, no parallel trend 


# CZ_ext #

CZ_exta <- lm(log_annual_ext_defor ~ Year*TreatedF, data= CZ_ext_data_before)
summary(CZ_exta)

#|# NEW: no parallel trend

# TTF_ext #

TTF_exta <- lm(log_annual_ext_defor ~ Year*TreatedF, data= TTF_ext_data_before)
summary(TTF_exta)

#|# NEW: no parallel trend

# All showed parallel trends except CFAM_ext which cannot be used in individual DiD regressions
#|# same for the replication


# d) DiD Regression

# Formula = y ~ treatment + time + (treatment x time)

# Interaction between treated and time is the coefficent of interest. This represents the effect
# of an observation being in an offset, after protection on the log-transformed count of deforestation.
# If this is significant and negative it means protection significantly reduced deforestation within the offset, 
# relative to the counterfactual.


# ANK_ext # 

modelANK_ext <- lm(log_annual_ext_defor ~ TreatedF*TimeF, data= annual_ext_defor_ANK_ext)
summary(modelANK_ext) 

# Back-transform the estimate to get the treatment effect - the  percentage difference in average annual_ext deforestation between
# the offset and the estimated counterfactual following protection. 
# The estimated counterfactual is the average annual_ext deforestation in the matched control sample after the intervention, adjusted to 
# to account for pre-intervention differences between the two samples. 

exp(coef(modelANK_ext)[4])-1
exp(confint(modelANK_ext)[4,])-1

#|# NEW: -97.53414% to -83.88969% of reduction


# CZ_ext # 

modelCZ_ext <- lm(log_annual_ext_defor ~ TreatedF*TimeF, data= annual_ext_defor_CZ_ext)
summary(modelCZ_ext)

exp(coef(modelCZ_ext)[4])-1
exp(confint(modelCZ_ext)[4,])-1

#|# NEW CZ_ext: -81.16364% to -21.52678% of reduction


# TTF_ext #

modelTTF_ext <- lm(log_annual_ext_defor ~ TreatedF*TimeF, data= annual_ext_defor_TTF_ext)
summary(modelTTF_ext)

exp(coef(modelTTF_ext)[4])-1
exp(confint(modelTTF_ext)[4,])-1


# Results show a significant reduction in average annual_ext deforestation of 96% (- 89 to - 98%) in ANK_exterana
# and 66% (-27 to -84%) in the COnservation Zone following protection. 

# In Torotorofotsy, protection had no significant effect on deforestation
#|# NEW: -53.6359% to  298.0839%, 35.85582% => no effect



DiD_res <- rbind(t(tidy(modelANK_ext)), t(tidy(modelCZ_ext)), t(tidy(modelTTF_ext)))



#                               e) Extract back-transformed coefficients and confidence intervals from DiD models


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

ATT_ANK_ext <- tidy_data5e(modelANK_ext, "ANK_ext")
ATT_CZ_ext <- tidy_data5e(modelCZ_ext, "CZ_ext")
ATT_TTF_ext <- tidy_data5e(modelTTF_ext, "TTF_ext")

ATT_ANK_ext$offset <- 'ANK_ext'
ATT_CZ_ext$offset <- 'CZ_ext'
ATT_TTF_ext$offset <- 'TTF_ext'


ATT_all <- rbind(ATT_ANK_ext, ATT_CZ_ext, ATT_TTF_ext)

#------------------------------------------12) Fixed Effects Panel regression ----------------------------------------------#

# Second outcome regression. This allows us to estimate the effect of protection across the entire offset portfolio,
# controlling for site and year fixed effects. This helps to control for any unobserved bias.   


# a) Data Construction 

ANK_ext_FE_dat <- annual_ext_defor_ANK_ext
CFAM_ext_FE_dat <- annual_ext_defor_CFAM_ext
CZ_ext_FE_dat <- annual_ext_defor_CZ_ext
TTF_ext_FE_dat <- annual_ext_defor_TTF_ext

# Pool data. 

levels(ANK_ext_FE_dat$Sample)[levels(ANK_ext_FE_dat$Sample) == "Cont"] <- "ANK_ext1"
levels(CFAM_ext_FE_dat$Sample)[levels(CFAM_ext_FE_dat$Sample) == "Cont"] <- "CFAM_ext1"
levels(CZ_ext_FE_dat$Sample)[levels(CZ_ext_FE_dat$Sample) == "Cont"] <- "CZ_ext1"
levels(TTF_ext_FE_dat$Sample)[levels(TTF_ext_FE_dat$Sample) == "Cont"] <- "TTF_ext1"

# Change the name of the four matched control samples (Cont) to ensure that when the data are pooled, each control sample is considered a separate site
# In the pooled dataset we have 152 observations, 1 observation per sample (n= 8, 4 offsets and 4 control), per year (n=23)

FE_dat <- rbind(ANK_ext_FE_dat, CFAM_ext_FE_dat, CZ_ext_FE_dat, TTF_ext_FE_dat)
FE_dat$Year <- factor(FE_dat$Year)


# Make one predictor indicating treated status. Tr = 1 for observations from an offset after protection and 0 for
# observations from an offset before protection and the control sample. Cannot use Treated and Time and the interaction because these
# are collinear with the fixed effects. 


FE_dat$Tr <- ifelse(FE_dat$TreatedF== "1" & FE_dat$TimeF== "1",1,0) 


# Check for differences between groups and over time. 

plotmeans(log_annual_ext_defor ~ Sample, data = FE_dat)

plotmeans(log_annual_ext_defor ~ Year, data = FE_dat)



# b) Fixed Effects Panel Regression 

FE_all_ext <- plm(log_annual_ext_defor ~ Tr, 
                  data= FE_dat, index = c("Sample", "Year"), model= "within", effect = "twoways")
FE_all_ext_summary <- summary(FE_all_ext)

# Tr is the coefficient of interest. A significant negative coefficient indicates a significant reduction in the log-transformed
# count of deforestation following protection of the four biodiversity offsets.

# Back-transform the estimates.

coef_FE <- (exp(coef(FE_all))-1)*100              # Results show that protection reduced average annual_ext deforestation by 58% (-37 to -72%) across
continft_FE <- (exp(confint(FE_all))-1)*100           # the entire offset portfolio 

# c) Tests 

# Compare to simple ols regression to test whether the fixed effects are needed. 

ols2 <- lm(log_annual_ext_defor ~ Tr, data = FE_dat)
summary(ols2)

pFtest(FE_all, ols2)

# p<0.05 so there is significant heterogeneity between groups and over time - the fixed effects are needed.


#                               d) Try also using random effects

m1<-lmer(log_annual_ext_defor ~ Tr  + (1|Sample)+(1|Year),
         data= FE_dat)
summary(m1)

(exp(-0.7476)-1)*100
(exp(confint(m1))-1)*100

# Results are very similar, and within the confidence intervals to the results from the FE panel regression. 
# They indicate a 53% (-27 to -69%) reduction in deforestation following protection. 



################################################################################################################

# ------------------------- 13) Exporting all results of fixed effect and of pooled regressions                   #|# exporting extension results

################################################################################################################

#|# get package here: https://osf.io/nxwvd
#|# install.packages("i4results_0.1.0.tar.gz", repos = NULL, type = "source")
library(i4results) #load libraries
library(openxlsx)

#|#------------- a. 3 offsets:
robustness_ANK <- list(hansen_extended = modelANK_ext) #|# create a list of robustness checks, in my case Hansen et al. 2013 extended to 2023
robustness_CZ <- list(hansen_extended = modelCZ_ext)
robustness_TTF <- list(hansen_extended = modelTTF_ext)

#|# export
i4results(modelANK, robustness_ANK, out = "ANKresults_HansenExtended.xlsx") #|# put 1. the original model, 2. the list of robusteness checks (up to 5), and where to save it
i4results(modelCZ, robustness_CZ, out = "CZresults_HansenExtended.xlsx")
i4results(modelTTF, robustness_TTF, out = "TTFresults_HansenExtended.xlsx")


#|#------------- b. pooled FE:
#|# check i4results function; #|# not working for plm() regression!
#|# View(i4results) #|# modify the "t value" to "t-value"

#|# new function:
i4results_for_plm <- function(original_model, robustness_models, out = NULL, append = FALSE){
  if (!is.null(out) && !requireNamespace("openxlsx", quietly = TRUE)) {
    stop("Package 'openxlsx' needed for Excel import/export. Please install it.")
  }
  if (length(robustness_models) > 5) {
    stop("You can supply up to 5 robustness models only.")
  }
  o_cmdline <- paste(deparse(original_model$call), collapse = " ")
  o_n <- stats::nobs(original_model)
  summary_o <- summary(original_model)
  coefs_o <- as.data.frame(summary_o$coefficients)
  coefs_o <- coefs_o[!(rownames(coefs_o) %in% "(Intercept)"), 
                     , drop = FALSE]
  ci_o <- confint(original_model, level = 0.95)
  ci_o <- as.data.frame(ci_o)
  colnames(ci_o) <- c("lower", "upper")
  ci_o <- ci_o[!(rownames(ci_o) %in% "(Intercept)"), , drop = FALSE]
  paramlist_all <- rownames(coefs_o)
  results_list <- list()
  for (rep_name in names(robustness_models)) {
    rep_model <- robustness_models[[rep_name]]
    r_cmdline <- paste(deparse(rep_model$call), collapse = " ")
    r_n <- stats::nobs(rep_model)
    summary_r <- summary(rep_model)
    coefs_r <- as.data.frame(summary_r$coefficients)
    coefs_r <- coefs_r[!(rownames(coefs_r) %in% "(Intercept)"), 
                       , drop = FALSE]
    ci_r <- confint(rep_model, level = 0.95)
    ci_r <- as.data.frame(ci_r)
    colnames(ci_r) <- c("lower", "upper")
    ci_r <- ci_r[!(rownames(ci_r) %in% "(Intercept)"), , 
                 drop = FALSE]
    for (p in paramlist_all) {
      o_coeff <- coefs_o[p, "Estimate"]
      o_std_err <- coefs_o[p, "Std. Error"]
      o_t <- coefs_o[p, "t-value"]
      o_p_val <- coefs_o[p, "Pr(>|t|)"]
      o_ci_lower <- ci_o[p, "lower"]
      o_ci_upper <- ci_o[p, "upper"]
      if (p %in% rownames(coefs_r)) {
        r_coeff <- coefs_r[p, "Estimate"]
        r_std_err <- coefs_r[p, "Std. Error"]
        r_t <- coefs_r[p, "t-value"]
        r_p_val <- coefs_r[p, "Pr(>|t|)"]
        r_ci_lower <- ci_r[p, "lower"]
        r_ci_upper <- ci_r[p, "upper"]
      }
      else {
        r_coeff <- NA
        r_std_err <- NA
        r_t <- NA
        r_p_val <- NA
        r_ci_lower <- NA
        r_ci_upper <- NA
      }
      record <- data.frame(paramname = p, study = rep_name, 
                           o_cmdline = substring(o_cmdline, 1, 244), r_cmdline = substring(r_cmdline, 
                                                                                           1, 244), o_n = o_n, r_n = r_n, o_coeff = round(o_coeff, 
                                                                                                                                          3), o_std_err = round(o_std_err, 3), o_t = round(o_t, 
                                                                                                                                                                                           3), o_p_val = round(o_p_val, 3), o_ci_lower = round(o_ci_lower, 
                                                                                                                                                                                                                                               3), o_ci_upper = round(o_ci_upper, 3), r_coeff = round(r_coeff, 
                                                                                                                                                                                                                                                                                                      3), r_std_err = round(r_std_err, 3), r_t = round(r_t, 
                                                                                                                                                                                                                                                                                                                                                       3), r_p_val = round(r_p_val, 3), r_ci_lower = round(r_ci_lower, 
                                                                                                                                                                                                                                                                                                                                                                                                           3), r_ci_upper = round(r_ci_upper, 3), stringsAsFactors = FALSE)
      results_list[[length(results_list) + 1]] <- record
    }
  }
  results_df <- do.call(rbind, results_list)
  if (!is.null(out)) {
    if (append && file.exists(out)) {
      oldres <- openxlsx::read.xlsx(out)
      combined <- rbind(oldres, results_df)
      openxlsx::write.xlsx(combined, file = out, asTable = FALSE)
      message(paste("Appended results to Excel file:", 
                    out))
    }
    else {
      openxlsx::write.xlsx(results_df, file = out, asTable = FALSE)
      if (append) {
        message(paste("File did not exist, so created a new Excel file:", 
                      out))
      }
      else {
        message(paste("Exported results to Excel file:", 
                      out))
      }
    }
  }
  else {
    message("Data with original & robustness comparisons returned as a data frame.")
  }
  invisible(results_df)
} # redo the function for plm() changing t value to t-value


#|# export
robustness_PooledFE <- list(hansen_extended = FE_all_ext)
plm_rob <- i4results_for_plm(original_model = FE_all, robustness_models = robustness_PooledFE, out = "Pooled_FE_HansenExtended.xlsx")











