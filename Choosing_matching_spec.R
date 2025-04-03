# Code by Katie Devenish. ktd19ycv@bangor.ac.uk

# This code is used in the publication entitled "Madagascar's biggest mine is set to achieve No Net Loss of Forest"
# Authors: Katie Devenish, Sebastien Desbureaux, Simon Willcock and Julia Jones. 2021. 



                        ### Choosing the main matching specification ###

# Following the recommendations of Schleicher et al (2019) we tested several matching
# specifications and chose the one which maximised the number of treated units matched 
# and the covariate balance as the main specification. The main specification was chosen
# before any outcome regressions were run.


# All specifications used nearest-neighbour matching and Mahalanobis distance.
# Unless otherwise specified matching was without replacement and with a ratio of one to one.


# 1) 1sd calipers:

cal <- rep(1, length(variables))        
names(cal) <- variables

matching <- function(offset){
  output = matchit(treated ~ Slope + Elevation + Dist_road + Dist_edge + Dist_defor,
                   data= offset, method= "nearest", distance = "mahalanobis", replace = FALSE, caliper=cal)
  return(output)
}    

m.out.ANK <- matching(offset = ANKCONT)      
m.out.CFAM <- matching(offset = CFAMCONT)     
m.out.CZ <- matching(offset = CZCONT)
m.out.TTF <- matching(offset = TTFCONT)


# 2) 0.5sd calipers:

cal <- rep(0.5, length(variables))        
names(cal) <- variables


m.out.ANK2 <- matching(offset = ANKCONT)      
m.out.CFAM2 <- matching(offset = CFAMCONT)     
m.out.CZ2 <- matching(offset = CZCONT)
m.out.TTF2 <- matching(offset = TTFCONT)


# 3) 0.25sd calipers:

cal <- rep(0.25, length(variables))        
names(cal) <- variables

m.out.ANK3 <- matching(offset = ANKCONT)     
m.out.CFAM3 <- matching(offset = CFAMCONT)     
m.out.CZ3 <- matching(offset = CZCONT)
m.out.TTF3 <- matching(offset = TTFCONT)


# 4) 0.5sd calipers + matching with replacement

cal <- rep(0.5, length(variables))        
names(cal) <- variables

matching <- function(offset){
  output = matchit(treated ~ Slope + Elevation + Dist_road + Dist_edge + Dist_defor,
                   data= offset, method= "nearest", distance = "mahalanobis", replace = TRUE, caliper=cal)
  return(output)
} 

m.out.ANK4 <- matching(offset = ANKCONT)     
m.out.CFAM4 <- matching(offset = CFAMCONT)     
m.out.CZ4 <- matching(offset = CZCONT)
m.out.TTF4 <- matching(offset = TTFCONT)


# 5) 1sd calipers + matching 1:5 nearest-neighbours. 

cal <- rep(1, length(variables))        
names(cal) <- variables

matching <- function(offset){
  output = matchit(treated ~ Slope + Elevation + Dist_road + Dist_edge + Dist_defor,
                   data= offset, method= "nearest", distance = "mahalanobis", replace = FALSE, caliper=cal, ratio = 5)
  return(output)
} 

m.out.ANK5 <- matching(offset = ANKCONT)     
m.out.CFAM5 <- matching(offset = CFAMCONT)     
m.out.CZ5 <- matching(offset = CZCONT)
m.out.TTF5 <- matching(offset = TTFCONT)


                            # Summarise results # 

# Group matching output from each test iteration into a List.
# Run the following section for each List changing List where marked with #####


List1 <- list(m.out.ANK, m.out.CFAM, m.out.CZ, m.out.TTF)
List2 <- list(m.out.ANK2, m.out.CFAM2, m.out.CZ2, m.out.TTF2)
List3 <- list(m.out.ANK3, m.out.CFAM3, m.out.CZ3, m.out.TTF3)
List4 <- list(m.out.ANK4, m.out.CFAM4, m.out.CZ4, m.out.TTF4)
List5 <- list(m.out.ANK5, m.out.CFAM5, m.out.CZ5, m.out.TTF5)

sum1_matching <- function(matched_data){
  sum1 = summary(matched_data, standardize = TRUE)    # Standardise = TRUE to get standardised mean difference in output
  return(sum1)
}

Output_names <- c("m.out.ANK", "m.out.CFAM", "m.out.CZ", "m.out.TTF")

allOutputs <- lapply(List5, sum1_matching)      #####  Apply the summary to all matched datasets in the chosen list. 

names(allOutputs) <- paste0(Output_names)  

# Repeat but with standardise = FALSE to get the mean EQQ values. This is a second measure to assess covariate balance. #

sum2_matching <- function(matched_data){
  sum2 = summary(matched_data, standardize = FALSE)
  return(sum2)
}

allOutputs2 <- lapply(List5, sum2_matching)     #####
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

cov_balance_ANK5 <- cov_balance(x=1, "SMD_1cal_Ratio5")         # Rename output according to iteration (1,2,3,4, or 5)
cov_balance_CFAM5 = cov_balance(x=2, "SMD_1cal_Ratio5")
cov_balance_CZ5 = cov_balance(x=3, "SMD_1cal_Ratio5")
cov_balance_TTF5 = cov_balance(x=4, "SMD_1cal_Ratio5")


# Combine results from different matching specifications.#

ANK_cov_balances <- cbind(cov_balance_ANK, cov_balance_ANK2, cov_balance_ANK3, cov_balance_ANK4, cov_balance_ANK5)
ANK_cov_balances <- subset(ANK_cov_balances, select = c(4,5,3,9,15,21,27))         # Select columns with standardised mean difference from each matching iteration. This is used to assess match quality.

ANK_matches <- c("NA", "NA", m.out.ANK$nn[5,2], m.out.ANK2$nn[5,2], m.out.ANK3$nn[5,2], m.out.ANK4$nn[5,2], m.out.ANK5$nn[5,2])   # Add info on number of matched treated units.
ANK_cov_balances <- rbind(ANK_cov_balances, ANK_matches)
rownames(ANK_cov_balances)[11] <- "Num_unmatched_treated"

CZ_cov_balances <- cbind(cov_balance_CZ, cov_balance_CZ2, cov_balance_CZ3, cov_balance_CZ4, cov_balance_CZ5)
CZ_cov_balances <- subset(CZ_cov_balances, select = c(4,5,3,9,15,21,27))

CZ_matches <- c("NA", "NA", m.out.CZ$nn[5,2], m.out.CZ2$nn[5,2], m.out.CZ3$nn[5,2], m.out.CZ4$nn[5,2], m.out.CZ5$nn[5,2])
CZ_cov_balances <- rbind(CZ_cov_balances, CZ_matches)
rownames(CZ_cov_balances)[11] <- "Num_unmatched_treated"

CFAM_cov_balances <- cbind(cov_balance_CFAM, cov_balance_CFAM2, cov_balance_CFAM3, cov_balance_CFAM4, cov_balance_CFAM5)
CFAM_cov_balances <- subset(CFAM_cov_balances, select = c(4,5,3,9,15,21,27))

CFAM_matches <- c("NA", "NA", m.out.CFAM$nn[5,2], m.out.CFAM2$nn[5,2], m.out.CFAM3$nn[5,2], m.out.CFAM4$nn[5,2], m.out.CFAM5$nn[5,2])
CFAM_cov_balances <- rbind(CFAM_cov_balances, CFAM_matches)
rownames(CFAM_cov_balances)[11] <- "Num_unmatched_treated"

TTF_cov_balances <- cbind(cov_balance_TTF, cov_balance_TTF2, cov_balance_TTF3, cov_balance_TTF4, cov_balance_TTF5)
TTF_cov_balances <- subset(TTF_cov_balances, select = c(4,5,3,9,15,21,27))

TTF_matches <- c("NA", "NA", m.out.TTF$nn[5,2], m.out.TTF2$nn[5,2], m.out.TTF3$nn[5,2], m.out.TTF4$nn[5,2], m.out.TTF5$nn[5,2])
TTF_cov_balances <- rbind(TTF_cov_balances, TTF_matches)
rownames(TTF_cov_balances)[11] <- "Num_unmatched_treated"


