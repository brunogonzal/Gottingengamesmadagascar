# Code by Katie Devenish. ktd19ycv@bangor.ac.uk

# This code is used in the publication entitled "Madagascar's biggest mine is set to achieve No Net Loss of Forest"
# Authors: Katie Devenish, Sebastien Desbureaux, Simon Willcock and Julia Jones. 2021.

# 2- We demonstratate that adding a time trend in the DID regressions doesnt change the results

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
library("stargazer")

# Code Structure

# 1) Old code in bracket
# 2) Did with the additions


# --------------------- 1) Data Construction for Matching -----------------------------------------#

setwd("c:/Users/desbures/Documents/Sebastien/11_Autres/Katie/Ambatovy/")


# Function to read in input data for each offset, remove unnecessary columns and add column to indicate treatment status
# and offset of origin. 

{
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
  data <- data %>% na_if(-9999)
  
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

}

#-------------------------------2) Difference in Differences with time trends -----------------------------------------#

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

modelANKtt <- lm(log_annual_defor ~ TreatedF*TimeF + Year, data= annual_defor_ANK)
summary(modelANKtt)


# CZ # 

modelCZ <- lm(log_annual_defor ~ TreatedF*TimeF, data= annual_defor_CZ)
summary(modelCZ)

modelCZtt <- lm(log_annual_defor ~ TreatedF*TimeF + Year, data= annual_defor_CZ)
summary(modelCZtt)

# TTF #

modelTTF <- lm(log_annual_defor ~ TreatedF*TimeF, data= annual_defor_TTF)
summary(modelTTF)

modelTTFtt <- lm(log_annual_defor ~ TreatedF*TimeF + Year, data= annual_defor_TTF)
summary(modelTTFtt)


stargazer(modelANK,modelANKtt,modelCZ,modelCZtt,modelTTF,modelTTFtt,
          title = "DiD with time trends",
          out="did_tt.doc",
          type = "html")
