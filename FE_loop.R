### FE regressions

library(dplyr)
library(foreach)

## 1 - We load the data
path_files = "YOUR PATH TO THE CREATED FE FILES FROM THE LOOPS/"
dbs <- list.files(path = paste0(path_files,"TTF"))


FEloop <- foreach (i = seq_along(dbs), .combine=rbind, .packages=c("plm")) %do% {
  iteration = dbs[i]  
  ANK <- read.csv(paste0(path_files,"ANK/",iteration))
  CFAM <- read.csv(paste0(path_files,"CFAM/",iteration))
  CZ <- read.csv(paste0(path_files,"CZ/",iteration))
  TTF <- read.csv(paste0(path_files,"TTF/",iteration))
  
ANK$Sample <- as.factor(ANK$Sample)
levels(ANK$Sample) <- c("ANK1", "ANK")

CFAM$Sample <- as.factor(CFAM$Sample)
levels(CFAM$Sample) <- c("CFAM1", "CFAM")

CZ$Sample <- as.factor(CZ$Sample)
levels(CZ$Sample) <- c("CZ1", "CZ")

TTF$Sample <- as.factor(TTF$Sample)
levels(TTF$Sample) <- c("TTF1", "TTF")


  FE_dat <- rbind(ANK, CFAM, CZ, TTF)
  
                # Test for parallel trends
  
  FE_dat_before <- FE_dat[FE_dat$TimeF == "before",]
  
  all_dat_PT <- lm(log(Annual_Deforestation+1) ~ Year*TreatedF, data = FE_dat_before)
  summary(all_dat_PT)
  
  
  # FE Panel Regression
  
  
  #FE_all1 <- plm(log(Annual_Deforestation+1) ~ TreatedF + TimeF + TreatedF:TimeF, 
  #              data= FE_dat, index = c("Sample", "Year"), model= "within", effect = "twoways")
  #summary(FE_all1)
  
  FE_dat$Tr <- ifelse(FE_dat$TreatedF=="treatment" & FE_dat$TimeF=="after",1,0) 
  FE_all <- plm(log(Annual_Deforestation+1) ~ Tr, 
                data= FE_dat, index = c("Sample", "Year"), model= "within", effect = "twoways")
  summary(FE_all)
  

  
  ## We join the info on covariates and methods
  
  
  
  # We export
  m1i <-   data.frame(coef=summary(FE_all)$coefficients[1,1],      # The coefficient of the estimation
                      se=summary(FE_all)$coefficients[1,2])        # The standard error of the coefficient
  }

## We join the info on covariates and methods
cov_inf  <- read.csv("YOURPATH/loopTTF_cov.csv")
param  <- read.csv("YOURPATH/loopTTF_parameters.csv")
covparam  <- read.csv("YOURPATH/loopTTF_cov_parameters.csv")
specinf <- rbind(cov_inf[,-c(1:7)],param[,-c(1:7)],covparam[,-c(1:7)])   

FEloop  <- cbind(FEloop, specinf)




labels <- list( "Additional covariates:" = c("Population density", "Dist. nearest settlement", "Annual rainfall", "Dist. nearest track", "Dist. nearest river"),
               "Model parameters:" = c("0.25SD caliper","0.5SD caliper", "1SD caliper", "With replacement"),
               "Distance:" = c("Mahalanobis", "Standard PSM", "Random forest PSM"),
               "Nb. nearest neighbors:" = c("1", "5", "10"))

tiff(file="YOURPATH/FE_rob.tiff",
     height = 600, width = 1000, pointsize=0.3)

schart(FEloop, labels,
       order = "increasing",
       highlight = 77,
       col.est = c("#999999","black"),
       bg.dot=c("grey60","grey95","grey95","grey60"),
       ylab = "Raw treatment effect",
       fonts=c(2,3))
dev.off()