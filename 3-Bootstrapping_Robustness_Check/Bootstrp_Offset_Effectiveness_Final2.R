
# This code is Niklas robustness check or "change" to the original script, mainly calculating the proper estimate using bootstrapping 
#|# <- comments by Niklas are marked like this


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

#|######## CALCULATE PROPENSITY SCORES BY HAND

propensity_scores_ANK <- glm(treated ~ Slope + Elevation + Dist_road + Dist_edge + Dist_defor, data = ANKCONT, family = binomial)
propensity_scores_CFAM <- glm(treated ~ Slope + Elevation + Dist_road + Dist_edge + Dist_defor, data = CFAMCONT, family = binomial)
propensity_scores_CZ <- glm(treated ~ Slope + Elevation + Dist_road + Dist_edge + Dist_defor, data = CZCONT, family = binomial)
propensity_scores_TTF <- glm(treated ~ Slope + Elevation + Dist_road + Dist_edge + Dist_defor, data = TTFCONT, family = binomial)


#|######## Change matching function to incorporate manually calculated scores

matching <- function(offset, distance){
  output = matchit(treated ~ Slope + Elevation + Dist_road + Dist_edge + Dist_defor,
                   data= offset, method= "nearest",distance = distance, replace = FALSE, caliper=cal)
  return(output)
}          


#|# Function to re-sample propensity scores

resample_log_model <- function(log_model){ # We simply resample the log-model randomly to produce a new log_model that we use in the meta-analysis
  
  # Sample new coefficients from standard normal distribution
  temp_1 <- as_tibble(log_model$coefficients) |> # make into dataframe
    add_column(sqrt(diag(vcov(log_model)))) |>
    rename(coefficients = value) |>
    rename(std_err = "sqrt(diag(vcov(log_model)))") |>
    rowwise() |>
    mutate(resampled_coefficients = coefficients + rnorm(1,mean=0, sd = std_err))  # resample coefficient estimates

  # Return a new model object with tweaked coefficients
  log_model$coefficients <- temp_1$resampled_coefficients
  
  # Re-calculate fitted values by hand
  
  log_odds <- model.matrix(log_model) %*% log_model$coefficients #matrix multiplication of data and newly calculated coefficients - returns log odds
  
  log_model$fitted.values <- as.vector(exp(log_odds) / (1 + exp(log_odds))) # We calculate probilities by inverting the logistic function
  
  return(log_model)
}


#|# initialize empty results dataframe 

DiD_bootstrap_results <- tribble(
  ~ANK, ~CZ, ~TTF, ~ANK_std_err, ~CZ_std_err, ~TTF_std_err
)


#|##|# BOOTSTRAP BEGINS HERE 

set.seed(142857)

for (boot in c(1:1000)){

  #|# REDUCE THERMAL LOAD ON MY NOT-SO-POWERFUL LAPTOP, just delete when running proper hardware
  Sys.sleep(0.1)

#|# Actually re-sample propensity scores

propensity_scores_ANK <- resample_log_model(propensity_scores_ANK)
propensity_scores_CFAM <- resample_log_model(propensity_scores_CFAM)
propensity_scores_CZ <- resample_log_model(propensity_scores_CZ)
propensity_scores_TTF <- resample_log_model(propensity_scores_TTF)



#|#|### Re-sample matches with new propensity scores

m.out.ANK <- matching(offset = ANKCONT, distance= propensity_scores_ANK$fitted.values)
m.out.CFAM <- matching(offset = CFAMCONT, distance= propensity_scores_CFAM$fitted.values) 
m.out.CZ <- matching(offset = CZCONT, distance = propensity_scores_CZ$fitted.values)
m.out.TTF <- matching(offset = TTFCONT, distance = propensity_scores_TTF$fitted.values)


#|# We kick "Summarise results" because we are interested in outcomes only #


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




# d) DiD Regression

# Formula = y ~ treatment + time + (treatment x time)

# Interaction between treated and time is the coefficent of interest. This represents the effect
# of an observation being in an offset, after protection on the log-transformed count of deforestation.
# If this is significant and negative it means protection significantly reduced deforestation within the offset, 
# relative to the counterfactual.


# ANK # 

modelANK <- lm(log_annual_defor ~ TreatedF*TimeF, data= annual_defor_ANK)
#|#summary(modelANK) 

# Back-transform the estimate to get the treatment effect - the  percentage difference in average annual deforestation between
# the offset and the estimated counterfactual following protection. 
# The estimated counterfactual is the average annual deforestation in the matched control sample after the intervention, adjusted to 
# to account for pre-intervention differences between the two samples. 

#|#exp(coef(modelANK)[4])-1
#|#exp(confint(modelANK)[4,])-1

# CZ # 

modelCZ <- lm(log_annual_defor ~ TreatedF*TimeF, data= annual_defor_CZ)
#|#summary(modelCZ)

#|#exp(coef(modelCZ)[4])-1
#|#exp(confint(modelCZ)[4,])-1

# TTF #

modelTTF <- lm(log_annual_defor ~ TreatedF*TimeF, data= annual_defor_TTF)
#|#summary(modelTTF)

#|#exp(coef(modelTTF)[4])-1
#|#exp(confint(modelTTF)[4,])-1

# Results show a significant reduction in average annual deforestation of 96% (- 89 to - 98%) in Ankerana
# and 66% (-27 to -84%) in the COnservation Zone following protection. 

# In Torotorofotsy, protection had no significant effect on deforestation. 

#|#DiD_res <- rbind(t(tidy(modelANK)), t(tidy(modelCZ)), t(tidy(modelTTF)))



#                               e) Extract back-transformed coefficients and confidence intervals from DiD models


#|#tidy_data5e <- function(model, name){
#|#  b <- (exp(coef(model)[4])-1)*100              # Back-transform coefficient of interaction term (exp(coefficient)-1)*100 to convert to % difference
#|#  c <- (exp(confint(model)[4,c(1,2)])-1)*100    # Back-transform confidence intervals
#|#  c <- data.frame(t(c))
#|#  b <- data.frame(b)
#|#  b <- cbind(b,c)
#|#  names(b) <- c("ATT", "Lower_CI", "Upper_CI")
#|#  rownames(b) <- name
#|#  return(b)
#|#}

#|#ATT_ANK <- tidy_data5e(modelANK, "ANK")
#|#ATT_CZ <- tidy_data5e(modelCZ, "CZ")
#|#ATT_TTF <- tidy_data5e(modelTTF, "TTF")

#|#ATT_all <- rbind(ATT_ANK, ATT_CZ, ATT_TTF)


#|# Save the computed estimates

DiD_bootstrap_results <- DiD_bootstrap_results |>
  add_row(ANK=coef(modelANK)[4],#save coefficient estimates
            CZ= coef(modelCZ)[4],
            TTF=coef(modelTTF)[4],
            ANK_std_err = tidy(modelANK)$std.error[4],#save standard errors, we resample later from all estimates 
            CZ_std_err = tidy(modelCZ)$std.error[4],
            TTF_std_err = tidy(modelTTF)$std.error[4])



#|# We kick steps 6 and ongoing from the script and focus on the DiD estimate

#|# BOOTSTRAP ENDS HERE
}


#|#### Calculate final estimates including uncertainty in matching 


set.seed(428571)


ANK <- DiD_bootstrap_results |>
  rowwise() |>
  mutate(ANK_distribution = list(rnorm(1000, mean =ANK, sd = ANK_std_err))) |> # Take 1000 samples from each estimate 
  unnest(ANK_distribution) |> # equals 1000*1000=10^6 point estimates of the true effect
  dplyr::select(ANK_distribution) # kick the other variables

CZ <- DiD_bootstrap_results |>
  rowwise() |>
  mutate(CZ_distribution = list(rnorm(1000, mean = CZ, sd = CZ_std_err))) |> # Take 1000 samples from each estimate 
  unnest(CZ_distribution) |> # equals 1000*1000=10^6 point estimates of the true effect
  dplyr::select(CZ_distribution) # kick the other variables

TTF <- DiD_bootstrap_results |>
  rowwise() |>
  mutate(TTF_distribution = list(rnorm(1000, mean =TTF, sd = TTF_std_err))) |> # Take 1000 samples from each estimate 
  unnest(TTF_distribution) |> # equals 1000*1000=10^6 point estimates of the true effect
  dplyr::select(TTF_distribution) # kick the other variables


#|# OPTIONAL: Make plots (uncomment with cmd + shift +c)

# Plot each result 
# 
# ggplot(TTF, aes(x = TTF_distribution)) +
#   geom_histogram(aes(y = ..density..), bins = 50, fill = "lightblue", color = "black") +
#   geom_density(color = "red3", linewidth = 1) +
#   ggtitle("Mixture of Normal Distributions") +
#   xlab("Coefficient estimate")

# All distributions look pretty normal so the end results are 

mean(ANK$ANK_distribution)
sd(ANK$ANK_distribution)

mean(CZ$CZ_distribution)
sd(CZ$CZ_distribution)

mean(TTF$TTF_distribution)
sd(TTF$TTF_distribution)


