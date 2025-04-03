#### HARMONIZING WITH THE COV LOOPS FOR AGGREGATE AND REPLACE




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
library("foreach")

# Code Structure

# 1) Data Construction for Matching
# 2) Matching
# 3) Plot Matching Results - Covariate Balance pre and post matching
# 4) Extract Matched Pairs from matching output
# 5) Data construction for DiD regression using matched pairs
# 6) DiD Regression 
# 7) Quantifying Avoided Deforestation
# 8) Fixed Effects Panel Regression


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

# Merge each offset with the control dataset and create a column with row ID.
# These offset + control datasets will the input for matching #

ANKCONT <- rbind(ANK, Control)
ANKCONT$ID <- seq(nrow(ANKCONT))
CZCONT <- rbind(CZ, Control)
CZCONT$ID <- seq(nrow(CZCONT))
CFAMCONT <- rbind(CFAM, Control)
CFAMCONT$ID <- seq(nrow(CFAMCONT))
TTFCONT <- rbind(TTF, Control)
TTFCONT$ID <- seq(nrow(TTFCONT))





#----------------------------2) Robustness: addition of covariates -----------------------------------------#


# 1 - Load the functions in the separete script

# 2 -  We evaluate the function
the_cov_loop(ANKCONT,11)
write.csv(synthesis_loop_cov,"YOURPATH/loopAnk_cov.csv")


the_cov_loop(CFAMCONT,13)
write.csv(synthesis_loop_cov,"YOURPATH/loopCFAM_cov.csv")

the_cov_loop(CZCONT,9)
write.csv(synthesis_loop_cov,"YOURPATH/loopCZ_cov.csv")

the_cov_loop(TTFCONT,14)
write.csv(synthesis_loop_cov,"YOURPATH/loopTTF_cov.csv")




#----------------------------3) Robustness: different parameters -----------------------------------------#

memory.limit(100000)
the_parameters_loop(ANKCONT,11,
                    caliper_value = c(0.25,0.5,1), 
                    replacement = c(T,F), 
                    distance = c("mahalanobis","glm","randomforest"),
                    ratio=c(1,5,10))
 write.csv(loop_parameters,"YOURPATH/loopANK_parameters.csv")

the_parameters_loop(CFAMCONT,13,
                    caliper_value = c(0.25,0.5,1), 
                    replacement = c(T,F), 
                    distance = c("mahalanobis","glm","randomforest"),
                    ratio=c(1,5,10))
write.csv(loop_parameters,"YOURPATH/loopCFAM_parameters.csv")

the_parameters_loop(CZCONT,9,
                    caliper_value = c(0.25,0.5,1), 
                    replacement = c(T,F), 
                    distance = c("mahalanobis","glm","randomforest"),
                    ratio=c(1,5,10))
 write.csv(loop_parameters,"YOURPATH/loopCZ_parameters.csv")

the_parameters_loop(TTFCONT,14,
                    replacement = c(T,F), 
                    distance = c("mahalanobis","glm","randomforest"),
                    caliper_value = c(0.25,0.5,1), 
                    ratio=c(1,5,10))
write.csv(loop_parameters,"YOURPATH/loopTTF_parameters.csv")

 
 
 
 
 
 
 
#----------------------------4) Robustness: different parameters and different covariates -----------------------------------------#

the_cov_parameters_loop(ANKCONT,11,
                    caliper_value = c(0.25,0.5,1), 
                    replacement = c(T,F), 
                    distance = c("mahalanobis","glm","randomforest"),
                    ratio=c(1,5,10))
write.csv(synthesis_loop_cov_param,"YOURPATH/loopANK_cov_parameters.csv")

the_cov_parameters_loop(CFAMCONT,13,
                    caliper_value = c(0.25,0.5,1), 
                    replacement = c(T,F), 
                    distance = c("mahalanobis","glm","randomforest"),
                    ratio=c(1,5,10))
write.csv(synthesis_loop_cov_param,"YOURPATH/loopCFAM_cov_parameters.csv")

the_cov_parameters_loop(CZCONT,9,
                    caliper_value = c(0.25,0.5,1), 
                    replacement = c(T,F), 
                    distance = c("mahalanobis","glm","randomforest"),
                    ratio=c(1,5,10))
write.csv(synthesis_loop_cov_param,"YOURPATH/loopCZ_cov_parameters.csv")

the_cov_parameters_loop(TTFCONT,14,
                    replacement = c(T,F), 
                    distance = c("mahalanobis","glm","randomforest"),
                    caliper_value = c(0.25,0.5,1), 
                    ratio=c(1,5,10))
write.csv(synthesis_loop_cov_param,"YOURPATH/loopTTF_cov_parameters.csv")

 


#----------------------------5) Robustness: dropping years -----------------------------------------#
the_year_loop(ANKCONT,11)
did_res_ank <- did_res 
 
the_year_loop(CZCONT,9)
did_res_cz <- did_res 

the_year_loop(TTFCONT,14)
did_res_ttf <- did_res 


 
 
#-------------------------- Graphs  


labels <- list("Balance:" = c("<10% obs. unmatched",">90% obs. unmatched", "Mean bal. achieved", "Balanced for all cov."),
               "Parallel trend:" = c("Achieved"), 
               "Additional covariates:" = c("Population density", "Dist. nearest settlement", "Annual rainfall", "Dist. nearest track", "Dist. nearest river"),
               "Model parameters:" = c("0.25SD caliper","0.5SD caliper", "1SD caliper", "With replacement"),
               "Distance:" = c("Mahalanobis", "Standard PSM", "Random forest PSM"),
               "Nb. nearest neighbors:" = c("1", "5", "10"))


## ANK
ANK_cov  <- read.csv("YOURPATH/loopANK_cov.csv")
ANK_param  <- read.csv("YOURPATH/loopANK_parameters.csv")
ANK_covparam  <- read.csv("YOURPATH/loopANK_cov_parameters.csv")
ANK <- rbind(ANK_cov[,-1],ANK_param[,-1],ANK_covparam[,-1])   
ANK$less_10pc_unmatched  <- ifelse(ANK$unmatched<2862*0.1,T,F)
ANK$more_90pc_unmatched  <- ifelse(ANK$unmatched>2862*0.9,T,F)
ANK$bal__mean_achieved  <- ifelse(ANK$mean_std_diff<0.25,T,F)
ANK$bal_achieved_all_cov  <- ifelse(ANK$max_std_diff<0.25,T,F)
ANK  <- ANK[,-c(3:5)]
ANK  <- ANK[,c(1,2,19:22,3:18)]

rows_invalid_models = which( 
                         ANK$more_90pc_unmatched==T | 
                         ANK$bal__mean_achieved==F |
                         ANK$parralel_trend==F)  

row_main_model = which(ANK$cal1==T & ANK$replacement==F & ANK$maha==T & ANK$ratio1==T & ANK$Pop_density==F & ANK$Dist_sett==F & ANK$Annual_Rain==F &ANK$Dist_track==F & ANK$Dist_river==F)  



tiff(file="YOURPATH/ANK_rob4.tiff",
     height = 550, width = 950, pointsize=0.3)

schart(ANK, labels, 
       order = "increasing",
       highlight = rows_invalid_models,
       highlight2 = row_main_model,
       col.est = c("#E69F00","#fcd67e","black"),
       bg.dot=c("grey60","grey95","grey95","grey60"),
       ylab = "Raw treatment effect",
       fonts=c(2,3))
#text(x=row_main_model,y=-1.5,"*",pos=3,cex=1.5)

dev.off()


## CFAM
CFAM_cov  <- read.csv("YOURPATH/loopCFAM_cov.csv")
CFAM_param  <- read.csv("YOURPATH/loopCFAM_parameters.csv")
CFAM_covparam  <- read.csv("YOURPATH/loopCFAM_cov_parameters.csv")
CFAM <- rbind(CFAM_cov[,-c(1)],CFAM_param[,-c(1)],CFAM_covparam[,-c(1)])   
CFAM$less_10pc_unmatched  <- ifelse(CFAM$unmatched<2626*0.1,T,F)
CFAM$more_90pc_unmatched  <- ifelse(CFAM$unmatched>2626*0.9,T,F)
CFAM$bal__mean_achieved  <- ifelse(CFAM$mean_std_diff<0.25,T,F)
CFAM$bal_achieved_all_cov  <- ifelse(CFAM$max_std_diff<0.25,T,F)
CFAM  <- CFAM[,-c(3:5)]
CFAM  <- CFAM[,c(1,2,19:22,3:18)]

rows_invalid_models = which( 
                            CFAM$more_90pc_unmatched==T | 
                            CFAM$bal__mean_achieved==F |
                            CFAM$parralel_trend==F)  

row_main_model = which(CFAM$cal1==T & CFAM$replacement==F & CFAM$maha==T & CFAM$ratio1==T & CFAM$Pop_density==F & CFAM$Dist_sett==F & CFAM$Annual_Rain==F &CFAM$Dist_track==F & CFAM$Dist_river==F)  




tiff(file="YOURPATH/CFAM_rob4.tiff",
     height = 600, width = 950, pointsize=0.3)

schart(CFAM, labels,
       order = "increasing",
       highlight = rows_invalid_models,
       highlight2 = row_main_model,
       col.est = c("chartreuse3","#b7ff6e","black"),
       bg.dot=c("grey60","grey95","grey95","grey60"),
       ylab = "Raw treatment effect",
       fonts=c(2,3))
dev.off()


## CZ
CZ_cov  <- read.csv("YOURPATH/loopCZ_cov.csv")
CZ_param  <- read.csv("YOURPATH/loopCZ_parameters.csv")
CZ_covparam  <- read.csv("YOURPATH/loopCZ_cov_parameters.csv")
  CZ <- rbind(CZ_cov[,-c(1)],CZ_param[,-c(1)],CZ_covparam[,-c(1)])   
  CZ$less_10pc_unmatched  <- ifelse(CZ$unmatched<1340*0.1,T,F)
  CZ$more_90pc_unmatched  <- ifelse(CZ$unmatched>1340*0.9,T,F)
  CZ$bal__mean_achieved  <- ifelse(CZ$mean_std_diff<0.25,T,F)
  CZ$bal_achieved_all_cov  <- ifelse(CZ$max_std_diff<0.25,T,F)
  CZ  <- CZ[,-c(3:5)]
  CZ  <- CZ[,c(1,2,19:22,3:18)]
  
  rows_invalid_models = which( 
                              CZ$more_90pc_unmatched==T | 
                              CZ$bal__mean_achieved==F |
                              CZ$parralel_trend==F) 
  
  row_main_model = which(CZ$cal1==T & CZ$replacement==F & CZ$maha==T & CZ$ratio1==T & CZ$Pop_density==F & CZ$Dist_sett==F & CZ$Annual_Rain==F &CZ$Dist_track==F & CZ$Dist_river==F)  
  

tiff(file="YOURPATH/CZ_rob4.tiff",
     height = 550, width = 950, pointsize=0.3)

schart(CZ, labels, 
       order = "increasing",
       highlight = rows_invalid_models,
       highlight2 = row_main_model,
       col.est = c("#F0E442","#fffab3","black"),
       bg.dot=c("grey60","grey95","grey95","grey60"),
       ylab = "Raw treatment effect",
       fonts=c(2,3))
dev.off()


 
## TTF
TTF_cov  <- read.csv("YOURPATH/loopTTF_cov.csv")
TTF_param  <- read.csv("YOURPATH/loopTTF_parameters.csv")
TTF_covparam  <- read.csv("YOURPATH/loopTTF_cov_parameters.csv")
TTF <- rbind(TTF_cov[,-c(1)],TTF_param[,-c(1)],TTF_covparam[,-c(1)])   
TTF$less_10pc_unmatched  <- ifelse(TTF$unmatched<1170*0.1,T,F)
TTF$more_90pc_unmatched  <- ifelse(TTF$unmatched>1170*0.9,T,F)
TTF$bal__mean_achieved  <- ifelse(TTF$mean_std_diff<0.25,T,F)
TTF$bal_achieved_all_cov  <- ifelse(TTF$max_std_diff<0.25,T,F)
TTF  <- TTF[,-c(3:5)]
TTF  <- TTF[,c(1,2,19:22,3:18)]

rows_invalid_models = which( 
                          TTF$more_90pc_unmatched==T | 
                          TTF$bal__mean_achieved==F  |
                          TTF$parralel_trend==F)  

row_main_model = which(TTF$cal1==T & TTF$replacement==F & TTF$maha==T & TTF$ratio1==T & TTF$Pop_density==F & TTF$Dist_sett==F & TTF$Annual_Rain==F &TTF$Dist_track==F & TTF$Dist_river==F)  


tiff(file="YOURPATH/TTF_rob4.tiff",
     height = 600, width = 950, pointsize=0.3)

schart(TTF, labels, 
       order = "increasing",
       highlight = rows_invalid_models,
       highlight2 = row_main_model,
       col.est = c("#56B4E9","#bfe8ff","black"),
       bg.dot=c("grey60","grey95","grey95","grey60"),
       ylab = "Raw treatment effect",
       fonts=c(2,3))
dev.off() 





# --------- Graphs for years

### Ank
plot <- ggplot(did_res_ank, aes(x=Year, y=coef)) +
  geom_errorbar(aes(ymin = coef + 1.96*se, ymax = coef - 1.96*se, col=Drop)) +
  scale_color_manual(values=c("#E69F00","#fcd67e")) +
  geom_point() +
  labs(title = "Ankerana", x = "Year dropped", y = "Coefficient") +
  ylim(-4.5,0) +
  geom_hline(yintercept = 0, linetype=2) +
  guides(color=F) +
  #scale_x_continuous(breaks = 2000:2019) +
  theme_classic()
plot

ggsave(plot, 
       filename = "ANK_years.png")


### CZ
plot <- ggplot(did_res_cz, aes(x=Year, y=coef)) +
  geom_errorbar(aes(ymin = coef + 1.96*se, ymax = coef - 1.96*se, col=Drop)) +
  scale_color_manual(values=c("#F0E442","#fffab3")) +
  geom_point() +
  labs(title = "Conservation Zone", x = "Year dropped", y = "Coefficient") +
  ylim(-2.5,0) +
  geom_hline(yintercept = 0, linetype=2) +
  guides(color=F) +
  #scale_x_continuous(breaks = 2000:2019) +
  theme_classic()


ggsave(plot, 
       filename = "CZ_years.png")


### TTF
plot <- ggplot(did_res_ttf, aes(x=Year, y=coef)) +
  geom_errorbar(aes(ymin = coef + 1.96*se, ymax = coef - 1.96*se, col=Drop)) +
  scale_color_manual(values=c("#56B4E9","#bfe8ff")) +
  geom_point() +
  labs(title = "Torotorofotsy", x = "Year dropped", y = "Coefficient") +
  geom_hline(yintercept = 0, linetype=2) +
  guides(color=F) +
  #scale_x_continuous(breaks = 2000:2019) +
  theme_classic()


ggsave(plot, 
       filename = "TTF_years.png")


