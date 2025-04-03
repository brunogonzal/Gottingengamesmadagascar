# Code by Katie Devenish. ktd19ycv@bangor.ac.uk

# This code is used in the publication entitled "Madagascar's biggest mine is set to achieve No Net Loss of Forest"
# Authors: Katie Devenish, Sebastien Desbureaux, Simon Willcock and Julia Jones. 2021.

# 1- We demonstratate that population (a candidate for essential cov) is colinear with the rest of the covariates 

# Load libraries #

library("foreign")
library("dplyr")
library("xlsx")
library("stargazer")

# Code Structure

# 1) Data Construction for Matching
# 2) Population
# 3) Extract Matched Pairs from matching output
# 4) Data construction for DiD regression using matched pairs
# 5) DiD Regression 
# 6) Quantifying Avoided Deforestation
# 7) Fixed Effects Panel Regression


# --------------------- 1) Data Construction for Matching -----------------------------------------#

setwd("c:/Users/desbures/Documents/Sebastien/11_Autres/Katie/Ambatovy/")


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


dropna <- function(data){    ## Just like in the main code
  
  # Replace all -9999 with NA
  data <- data %>% na_if(-9999)
  
  # We keep observations without na
  data <- data %>% drop_na()
  
}

ANK <- dropna(data = ANK)
CZ <- dropna(CZ)
CF <- dropna(CFAM)
TTF <- dropna(TTF)
Control <- dropna(Control)

#----------------------------2) Examining Population -----------------------------------------#

## We regress population density on the other covariates to see if they correlate.
## To ease the reading of the results, we multiply pop density by 1000 (equivalent to transforming other covariates in meters)


## Pooled dataset
d <- rbind(ANK,CFAM,CZ,TTF,Control)
d$Pop_density <- d$Pop_density*1000

poolr <- lm(Pop_density ~ Dist_road + Dist_edge + Dist_defor + Slope + Elevation
           , data = d)

## Each offset and potential controls

## Ank
d <- rbind(ANK,Control)
d$Pop_density <- d$Pop_density*1000

ankr <- lm(Pop_density ~ Dist_road + Dist_edge + Dist_defor + Slope + Elevation
           , data = d)

## CFAM
d <- rbind(CFAM,Control)
d$Pop_density <- d$Pop_density*1000

cfamr <- lm(Pop_density ~ Dist_road + Dist_edge + Dist_defor + Slope + Elevation
                  , data = d)

## CZ
d <- rbind(CZ,Control)
d$Pop_density <- d$Pop_density*1000

czr <- lm(Pop_density ~ Dist_road + Dist_edge + Dist_defor + Slope + Elevation
                  , data = d)

## TTF
d <- rbind(TTF,Control)
d$Pop_density <- d$Pop_density*1000

ttfr <- lm(Pop_density ~ Dist_road + Dist_edge + Dist_defor + Slope + Elevation
                  , data = d)

## A summary table
stargazer(poolr, ankr, cfamr, czr, ttfr,
          title = "Correlation between Population and essential covariates",
          out="pop_results.doc",
          type = "html")

