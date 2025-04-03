## Jan 16 2021
## This R-script contains the functions created for robustness checks. Load them first



the_cov_loop <- function(offset, ycreation) { #},caliper_value,replacement){
  # A matrix will all possible covariates, and distinguish between essential and additional
  
  
  # We store in a string the name of the offset
  name_offset = offset$offset[1]
  
  #We pre-clean the data to speed up the computations by manually getting rid off observations that will be too far when using calipers
  offset_l <- offset %>% na_if(-9999)
  
  # We keep observations without na
  offset_l <- offset_l %>% drop_na()
  
  
  # We determine the caliprs for each variables
  caliper_value = 1
  caliper_slope<-caliper_value*sd(offset_l$Slope[offset_l$treated==1])
  caliper_elev<-caliper_value*sd(offset_l$Elevation[offset_l$treated==1])
  caliper_dist_road<-caliper_value*sd(offset_l$Dist_road[offset_l$treated==1])
  caliper_dist_edge<-caliper_value*sd(offset_l$Dist_edge[offset_l$treated==1])
  caliper_dist_def<-caliper_value*sd(offset_l$Dist_defor[offset_l$treated==1])
  
  # We keep control observations that are on the same support
  offset_l <- offset_l[offset_l$Slope > min(offset_l$Slope[offset_l$treated==1])-caliper_slope,]
  offset_l <- offset_l[offset_l$Slope < max(offset_l$Slope[offset_l$treated==1])+caliper_slope,]
  
  offset_l <- offset_l[offset_l$Elevation > min(offset_l$Elevation[offset_l$treated==1])-caliper_elev,] 
  offset_l <- offset_l[offset_l$Elevation < max(offset_l$Elevation[offset_l$treated==1])+caliper_elev,]  
  
  offset_l <- offset_l[offset_l$Dist_road > min(offset_l$Dist_road[offset_l$treated==1])-caliper_dist_road,]
  offset_l <- offset_l[offset_l$Dist_road < max(offset_l$Dist_road[offset_l$treated==1])+caliper_dist_road,] 
  
  offset_l <- offset_l[offset_l$Dist_edge > min(offset_l$Dist_edge[offset_l$treated==1])-caliper_dist_edge,] 
  offset_l <- offset_l[offset_l$Dist_edge < max(offset_l$Dist_edge[offset_l$treated==1])+caliper_dist_edge,] 
  
  offset_l <- offset_l[offset_l$Dist_defor > min(offset_l$Dist_defor[offset_l$treated==1])-caliper_dist_def,] 
  offset_l <- offset_l[offset_l$Dist_defor < max(offset_l$Dist_defor[offset_l$treated==1])+caliper_dist_def,]
  
  
  # We only keep the ppol of covaiates, treatment index and outcome
  db<-offset_l[,c(4:7,9:14)]
  treated <- offset_l$treated
  treeloss <- offset_l$Tree_loss
  
  
  
  # We define two set of covariates
  covariates_min <- names(db[,c(3,4,7,9,10)])  # A set of minimal covariates that should be included in any matching
  covariates_add <- names(db[,c(1,2,5,6,8)])    # A set of standard covariates that can be included as well


  
  # We determine how many different models can be constructed based on the inclusion of additional covariates
  id <- unlist(
    lapply(1:length(covariates_add),
           function(i)combn(1:length(covariates_add),i,simplify=FALSE)
    )
    ,recursive=FALSE)
  
  # We create all combinations of required + optional covariates
  mydbs_full <- lapply(id,function(i) c(covariates_min, covariates_add[i]))
  
  
  # We loop over all possible models based on any combination of 
  
  maha <- lapply(mydbs_full, function(cov) {
    
    idx <- match(cov, names(db))   # We determines the columns of the covariates
    idx <- sort(idx)               # We need to arrange them in increasing order
    NewDF <- db[,c(idx)]           # We create a database that contains only the covariates needed in this iteration of the loop
    
    formulae  <- as.formula(paste("treated ~ ", paste0(names(NewDF),collapse ="+")))  # We create a formula to be inserted in the matching
    
    
    
    # Calipers
    cal=rep(caliper_value,length(cov)) # In the version 4.0 of matchIT, we now need to put calipers for each variables included in the model
    names(cal) <- cov         # We name the vestor to be inserted in the matching function
    
    
    # Applying the matching function to 
    eval_code <- TRUE
    m.out <- tryCatch(
      matchit(formulae, data= cbind(treated,NewDF,treeloss), 
              method= "nearest", 
              distance = "mahalanobis", 
              replace = F, 
              caliper = cal), error = function(e) {eval_code <<- FALSE})
    
    
    if(eval_code) {
      # Summarise results #
      match_output <-summary(m.out, standardize = TRUE)
      
      # Extract useful summary values and combine into one dataset # 
      test <- as.data.frame(match_output$sum.matched[,3])
      colnames(test) <- "Stand_Mean_Diff"
      mean_std_diff <- mean(abs(test$Stand_Mean_Diff))        # Mean value of the std diff
      max_std_diff <- max(abs(test$Stand_Mean_Diff))          # Max value of the std diff
      unmatched <- match_output$nn[5,2]
      
      cov_balance <- cbind(unmatched, mean_std_diff,max_std_diff)
      
      
      # Create db
      m.data <- match.data(m.out)
      
      
      # Aggregate pixels into treatment (offset_l_l) and control groups. Tabulate the number of pixels by tree loss year within each group (count of number of pixels within each group deforested each year)
      annual_defor <- m.data %>%  # As weights are created with the replace or with PSM, I change the tabulate approach of base R by a dplyr count
        count(as.factor(treated), as.factor(treeloss), wt=weights, .drop = F)
      
      label <- c("Sample", "Year", "Annual_Deforestation")
      names(annual_defor) <- paste0(label)

      
      # Have to make Year numeric for >= to work
      annual_defor$Year <- as.numeric(as.character(annual_defor$Year))
      
      
      # We delete the 0 values (aka pixels that remain forested at the end)
      annual_defor <- annual_defor[annual_defor$Year!=0,]
      
      # Construct data with dummy variables for time and treated to use in DiD regression
      annual_defor$Time <- ifelse(annual_defor$Year >= ycreation, 1,0)
      annual_defor$TimeF <-factor(annual_defor$Time, levels = c(0,1), labels = c("before","after"))
      annual_defor$TreatedF <- factor(annual_defor$Sample, levels = c(0,1), labels = c("control","treatment"))
            
      # Export the data for the FE regression later on
      save_path_name <- paste0("c:/Users/Virunga/Desktop/Data_for_FE/",
                               name_offset,
                               "/",
                               paste(substr(cov,1,1), collapse = ""), paste(substr(cov,6,6), collapse = ""),
                               ".csv")       
      write.csv(annual_defor, save_path_name)
      
      
      # c) Test for parallel trends
      
      # Use only data from the before period #
      
      data_before <- annual_defor[(annual_defor$Year < ycreation),]
      
      Partest <- lm(log(Annual_Deforestation+1) ~ Year*TreatedF, data= data_before)
      summary(Partest)        
      parralel_trend <- ifelse(summary(Partest)$coefficients[4,4]> 0.05,T,F)    # TBC If the p-val of the interaction term is larger than 0.1, then // trend #
      
      
      
      
      # d) DiD Regression
      
      # Formula = y ~ treatment + time + (treatment x time)
      
      modeldid <- lm(log(Annual_Deforestation+1) ~ TreatedF*TimeF, data= annual_defor)
      summary(modeldid)     
      
      
      

      
      
      #---------------  We export all the information we want
      
      # A list of covariates that were included
      var <- as.data.frame(t(covariates_add %in% names(NewDF)))
      names(var)  <- covariates_add  
      
      # We export
      m1i <- data.frame(coef=summary(modeldid)$coefficients[4,1],      # The coefficient of the estimation
                          se=summary(modeldid)$coefficients[4,2],        # The standard error of the coefficient
                          cov_balance,                                # The covariate balance
                          parralel_trend,                                 # Parralel trend info
                          #                      Gen=as.logical(F), NN1=as.logical(F), ebal=as.logical(F), psm=as.logical(T), lm=as.logical(F),
                          var,                                             # The extra covariates included
                          cal025=as.logical(F), cal05=as.logical(F),cal1=as.logical(T),
                          replacement=as.logical(F),
                          maha=as.logical(T), glm=as.logical(F), rf=as.logical(F),
                          #NN=ifelse(meth=="nearest",as.logical(T),as.logical(F)), genetic=ifelse(meth=="genetic", as.logical(T), as.logical(F)),
                          ratio1=as.logical(T), ratio5=as.logical(F), ratio10=as.logical(F))
      
      
      
    }
  }
)
    
  
  synthesis_loop_cov <<- data.frame(do.call("rbind",maha))
}


the_parameters_loop <- function(offset, ycreation, caliper_value, replacement, distance, ratio) { #},method,){
  
  loop_parameters <<- foreach(cali=caliper_value, .combine=rbind, .packages=c("MatchIt","dplyr")) %do% {
    foreach(rep=replacement, .combine=rbind, .packages=c("MatchIt","dplyr")) %do% {
      foreach(dist=distance, .combine=rbind, .packages=c("MatchIt","dplyr")) %do% {
#       foreach(meth=method, .combine=rbind, .packages=c("MatchIt","dplyr")) %do% {    
          foreach(n=ratio, .combine=rbind, .packages=c("MatchIt","dplyr")) %do% {  
      
          
          # We store in a string the name of the offset
          name_offset = offset$offset[1]
            
              
          # Replace all -9999 from AKNCONT with NA
          offset_l <- offset %>% na_if(-9999)
          
          # We keep observations without na
          offset_l <- offset_l %>% drop_na()
          
          
          
          
          
          # We determine the caliprs for each variables
          caliper_slope<-cali*sd(offset_l$Slope[offset_l$treated==1])
          caliper_elev<-cali*sd(offset_l$Elevation[offset_l$treated==1])
          caliper_dist_road<-cali*sd(offset_l$Dist_road[offset_l$treated==1])
          caliper_dist_edge<-cali*sd(offset_l$Dist_edge[offset_l$treated==1])
          caliper_dist_def<-cali*sd(offset_l$Dist_defor[offset_l$treated==1])
          
          # We keep control observations that are on the same support
          offset_l <- offset_l[offset_l$Slope > min(offset_l$Slope[offset_l$treated==1])-caliper_slope,]
          offset_l <- offset_l[offset_l$Slope < max(offset_l$Slope[offset_l$treated==1])+caliper_slope,]
          
          offset_l <- offset_l[offset_l$Elevation > min(offset_l$Elevation[offset_l$treated==1])-caliper_elev,] 
          offset_l <- offset_l[offset_l$Elevation < max(offset_l$Elevation[offset_l$treated==1])+caliper_elev,]  
          
          offset_l <- offset_l[offset_l$Dist_road > min(offset_l$Dist_road[offset_l$treated==1])-caliper_dist_road,]
          offset_l <- offset_l[offset_l$Dist_road < max(offset_l$Dist_road[offset_l$treated==1])+caliper_dist_road,] 
          
          offset_l <- offset_l[offset_l$Dist_edge > min(offset_l$Dist_edge[offset_l$treated==1])-caliper_dist_edge,] 
          offset_l <- offset_l[offset_l$Dist_edge < max(offset_l$Dist_edge[offset_l$treated==1])+caliper_dist_edge,] 
          
          offset_l <- offset_l[offset_l$Dist_defor > min(offset_l$Dist_defor[offset_l$treated==1])-caliper_dist_def,] 
          offset_l <- offset_l[offset_l$Dist_defor < max(offset_l$Dist_defor[offset_l$treated==1])+caliper_dist_def,]
          
          
          # We only keep the pool of esential covaiates, treatment index and outcome
          treated <- offset_l$treated
          treeloss <- offset_l$Tree_loss
          offset_l<-offset_l[,c(6,7,11,13,14)]
          
          # We create a formula to be inserted in the matching  
          formulae  <- as.formula(paste("treated ~ ", paste0(names(offset_l),collapse ="+")))
          
          
          #maha <- lapply(mydbs_full, function(cov) {
          
          # Calipers
          cal=rep(cali,ncol(offset_l)) # In the version 4.0 of matchIT, we now need to put calipers for each variables included in the model
          names(cal) <- names(offset_l)         # We name the vestor to be inserted in the matching function
          
          
          # Applying the matching function to 
          eval_code <- TRUE
          m.out <- #tryCatch(
            matchit(formulae, data= cbind(treated,offset_l,treeloss), 
                    method= "nearest", 
                    distance = dist, 
                    replace = rep, 
                    caliper = cal, 
                    ratio=n)#, error = function(e) { eval_code <<- FALSE})
          
          
          if(eval_code) {
            # Summarise results #
            match_output <-summary(m.out, standardize = TRUE)
            
            # Extract useful summary values and combine into one dataset # 
            test <- as.data.frame(match_output$sum.matched[,3])
            colnames(test) <- "Stand_Mean_Diff"
            mean_std_diff <- mean(abs(test$Stand_Mean_Diff))        # Mean value of the std diff
            max_std_diff <- max(abs(test$Stand_Mean_Diff))          # Max value of the std diff
            unmatched <- match_output$nn[5,2]
            
            cov_balance <- cbind(unmatched, mean_std_diff,max_std_diff)
            
            
            # Create db
            m.data <- match.data(m.out)
            m.data$weights[m.data$treated==0]<-m.data$weights[m.data$treated==0]/5
            
            # Aggregate pixels into treatment (offset_l) and control groups. Tabulate the number of pixels by tree loss year within each group (count of number of pixels within each group deforested each year)
            annual_defor <- m.data %>%  # As weights are created with the replace or with PSM, I change the tabulate approach of base R by a dplyr count
              count(as.factor(treated), as.factor(treeloss), wt=weights, .drop = F)
            
            label <- c("Sample", "Year", "Annual_Deforestation")
            names(annual_defor) <- paste0(label)
            
            # Have to make Year numeric for >= to work
            annual_defor$Year <- as.numeric(as.character(annual_defor$Year))
            
            
            # We delete the 0 values (aka pixels that remain forested at the end)
            annual_defor <- annual_defor[annual_defor$Year!=0,]
            
            
            
            # Construct data with dummy variables for time and treated to use in DiD regression
            
            annual_defor$Time <- ifelse(annual_defor$Year >= ycreation, 1,0)
            annual_defor$TimeF <-factor(annual_defor$Time, levels = c(0,1), labels = c("before","after"))
            annual_defor$TreatedF <- factor(annual_defor$Sample, levels = c(0,1), labels = c("control","treatment"))
            

            # Export the data for the FE regression later on
            save_path_name <- paste0("c:/Users/Virunga/Desktop/Data_for_FE/",
                                     name_offset,"/",
                         cali,rep,dist,n,".csv")       
            write.csv(annual_defor, save_path_name)
            

                        # c) Test for parallel trends
            
            # Use only data from the before period #
            
            data_before <- annual_defor[(annual_defor$Year <ycreation),]
            
            Partest <- lm(log(Annual_Deforestation+1) ~ Year*TreatedF, data= data_before)
            summary(Partest)        
            parralel_trend <- ifelse(summary(Partest)$coefficients[4,4]> 0.05,T,F)    # TBC If the p-val of the interaction term is larger than 0.1, then // trend #
            
            
            
            
            # d) DiD Regression
            
            # Formula = y ~ treatment + time + (treatment x time)
            
            modeldid <- lm(log(Annual_Deforestation+1) ~ TreatedF*TimeF, data= annual_defor)
            summary(modeldid)     
            
            
          
            
            #---------------  We export all the information we want
            
            # A list of covariates that were included
            var <- as.data.frame(t(rep(FALSE,5)))
            names(var)  <- c("Pop_density", "Dist_sett", "Annual_Rain", "Dist_track", "Dist_river")  
            
            # We export
            m1i <-   data.frame(coef=summary(modeldid)$coefficients[4,1],      # The coefficient of the estimation
                                se=summary(modeldid)$coefficients[4,2],        # The standard error of the coefficient
                                cov_balance,                                # The covariate balance
                               parralel_trend,                                 # Parralel trend info
                               #                      Gen=as.logical(F), NN1=as.logical(F), ebal=as.logical(F), psm=as.logical(T), lm=as.logical(F),
                               var,
                               cal025=ifelse(cali==0.25,as.logical(T),as.logical(F)), cal05=ifelse(cali==0.5,as.logical(T),as.logical(F)),cal1=ifelse(cali==1,as.logical(T),as.logical(F)),
                               replacement=ifelse(rep==T,as.logical(T),as.logical(F)),
                               maha=ifelse(dist=="mahalanobis",as.logical(T),as.logical(F)), glm=ifelse(dist=="glm",as.logical(T),as.logical(F)), rf=ifelse(dist=="randomforest",as.logical(T),as.logical(F)),
#                               NN=ifelse(meth=="nearest",as.logical(T),as.logical(F)), genetic=ifelse(meth=="genetic", as.logical(T), as.logical(F)),
                               ratio1=ifelse(n==1,as.logical(T),as.logical(F)),ratio5=ifelse(n==5,as.logical(T),as.logical(F)),ratio10=ifelse(n==10,as.logical(T),as.logical(F))
            )
            
            
            
          }
        }}}}}





the_cov_parameters_loop <- function(offset, ycreation, distance, replacement, caliper_value, ratio) { #},caliper_value,replacement){
  # A matrix will all possible covariates, and distinguish between essential and additional
  
  # We store in a string the name of the offset
  name_offset = offset$offset[1]
  
  
  # Replace all -9999 from AKNCONT with NA
  offset_l <- offset %>% na_if(-9999)
  
  # We keep observations without na
  offset_l <- offset_l %>% drop_na()
  
  
  # We define two set of covariates
  covariates_min <- names(offset_l[,c(6,7,11,13,14)])  # A set of minimal covariates that should be included in any matching
  covariates_add <- names(offset_l[,c(4,5,9,10,12)])    # A set of standard covariates that can be included as well
  
  
  # We determine how many different models can be constructed based on the inclusion of additional covariates
  id <- unlist(
    lapply(1:length(covariates_add),
           function(i)combn(1:length(covariates_add),i,simplify=FALSE)
    )
    ,recursive=FALSE)
  
  # We create all combinations of required + optional covariates
  
  set.seed(12345)
  mydbs_full <- lapply(id,function(i) c(sample(distance,1),
                                        sample(replacement,1),
                                        sample(ratio,1),
                                        sample(caliper_value,1),
                                        covariates_min, 
                                        covariates_add[i]))
  

    

  # We loop over all possible models based on any combination of 
  
  maha <- lapply(mydbs_full, function(cov) {
  
    
    # We determine the caliprs for each variables
    cali=as.numeric(cov[4])
    caliper_slope<-cali*sd(offset_l$Slope[offset_l$treated==1])
    caliper_elev<-cali*sd(offset_l$Elevation[offset_l$treated==1])
    caliper_dist_road<-cali*sd(offset_l$Dist_road[offset_l$treated==1])
    caliper_dist_edge<-cali*sd(offset_l$Dist_edge[offset_l$treated==1])
    caliper_dist_def<-cali*sd(offset_l$Dist_defor[offset_l$treated==1])
    
    # We keep control observations that are on the same support
    offset_l <- offset_l[offset_l$Slope > min(offset_l$Slope[offset_l$treated==1])-caliper_slope,]
    offset_l <- offset_l[offset_l$Slope < max(offset_l$Slope[offset_l$treated==1])+caliper_slope,]
    
    offset_l <- offset_l[offset_l$Elevation > min(offset_l$Elevation[offset_l$treated==1])-caliper_elev,] 
    offset_l <- offset_l[offset_l$Elevation < max(offset_l$Elevation[offset_l$treated==1])+caliper_elev,]  
    
    offset_l <- offset_l[offset_l$Dist_road > min(offset_l$Dist_road[offset_l$treated==1])-caliper_dist_road,]
    offset_l <- offset_l[offset_l$Dist_road < max(offset_l$Dist_road[offset_l$treated==1])+caliper_dist_road,] 
    
    offset_l <- offset_l[offset_l$Dist_edge > min(offset_l$Dist_edge[offset_l$treated==1])-caliper_dist_edge,] 
    offset_l <- offset_l[offset_l$Dist_edge < max(offset_l$Dist_edge[offset_l$treated==1])+caliper_dist_edge,] 
    
    offset_l <- offset_l[offset_l$Dist_defor > min(offset_l$Dist_defor[offset_l$treated==1])-caliper_dist_def,] 
    offset_l <- offset_l[offset_l$Dist_defor < max(offset_l$Dist_defor[offset_l$treated==1])+caliper_dist_def,]
    
    
    # We only keep the ppol of covaiates, treatment index and outcome
    db<-offset_l[,c(4:7,9:14)]
    treated <- offset_l$treated
    treeloss <- offset_l$Tree_loss
      
    idx <- match(cov, names(db))   # We determines the columns of the covariates
    idx <- sort(idx)               # We need to arrange them in increasing order
    NewDF <- db[,c(idx)]           # We create a database that contains only the covariates needed in this iteration of the loop
    
    formulae  <- as.formula(paste("treated ~ ", paste0(names(NewDF),collapse ="+")))  # We create a formula to be inserted in the matching
    
      
    # Calipers
    cal=rep(cali,length(cov)-4) # In the version 4.0 of matchIT, we now need to put calipers for each variables included in the model
    names(cal) <- cov[-c(1:4)]         # We name the vestor to be inserted in the matching function
    
    
    # Applying the matching function to 
    eval_code <- TRUE
    m.out <- tryCatch(
      matchit(formulae, data= cbind(treated,NewDF,treeloss), 
              method= "nearest", 
              distance = as.character(cov[1]), 
              replace = as.logical(cov[2]),
              ratio=as.numeric(cov[3]),
              caliper = cal), error = function(e) { eval_code <<- FALSE})
    
    
    if(eval_code) {
      # Summarise results #
      match_output <-summary(m.out, standardize = TRUE)
      
      # Extract useful summary values and combine into one dataset # 
      test <- as.data.frame(match_output$sum.matched[,3])
      colnames(test) <- "Stand_Mean_Diff"
      mean_std_diff <- mean(abs(test$Stand_Mean_Diff))        # Mean value of the std diff
      max_std_diff <- max(abs(test$Stand_Mean_Diff))          # Max value of the std diff
      unmatched <- match_output$nn[5,2]
      
      cov_balance <- cbind(unmatched, mean_std_diff,max_std_diff)
      
      
      # Create db
      m.data <- match.data(m.out)
      m.data$weights[m.data$treated==0]<-m.data$weights[m.data$treated==0]/5
      
      # Aggregate pixels into treatment (offset_l_l) and control groups. Tabulate the number of pixels by tree loss year within each group (count of number of pixels within each group deforested each year)
      annual_defor <- m.data %>%  # As weights are created with the replace or with PSM, I change the tabulate approach of base R by a dplyr count
        count(as.factor(treated), as.factor(treeloss), wt=weights, .drop = F)

      
      label <- c("Sample", "Year", "Annual_Deforestation")
      names(annual_defor) <- paste0(label)

      # Have to make Year numeric for >= to work
      annual_defor$Year <- as.numeric(as.character(annual_defor$Year))
      
      
      # We delete the 0 values (aka pixels that remain forested at the end)
      annual_defor <- annual_defor[annual_defor$Year!=0,]
      
    
      # Construct data with dummy variables for time and treated to use in DiD regression
      annual_defor$Time <- ifelse(annual_defor$Year >= ycreation, 1,0)
      annual_defor$TimeF <-factor(annual_defor$Time, levels = c(0,1), labels = c("before","after"))
      annual_defor$TreatedF <- factor(annual_defor$Sample, levels = c(0,1), labels = c("control","treatment"))

      # Export the data for the FE regression later on
      save_path_name <- paste0("c:/Users/Virunga/Desktop/Data_for_FE/",
                               name_offset,
                               "/",
                               paste(substr(cov,1,1), collapse = ""), paste(substr(cov,6,6), collapse = ""),
                               ".csv")
      write.csv(annual_defor, save_path_name)
      


      if(nrow(annual_defor)>15) {
        
        
      # c) Test for parallel trends
      
      # Use only data from the before period #
      data_before <- annual_defor[(annual_defor$Year <ycreation),]
      
      Partest <- lm(log(Annual_Deforestation +1) ~ Year*TreatedF, data= data_before)
      summary(Partest)        

      parralel_trend <- #tryCatch(
        ifelse(summary(Partest)$coefficients[4,4]> 0.05,T,F)#,
#        error = function(e) {Dont_skip <<- FALSE}) # TBC If the p-val of the interaction term is larger than 0.1, then // trend #
      
      
      
      
      # d) DiD Regression
      
      # Formula = y ~ treatment + time + (treatment x time)
      
 
      modeldid <- lm(log(Annual_Deforestation +1) ~ TreatedF*TimeF, data= annual_defor)
      summary(modeldid)     
      
      
      

      #---------------  We export all the information we want
      
      # A list of covariates that were included
      var <- as.data.frame(t(covariates_add %in% names(NewDF)))
      names(var)  <- covariates_add  
      
      # We export
      m1i <-   data.frame(coef=summary(modeldid)$coefficients[4,1],      # The coefficient of the estimation
                          se=summary(modeldid)$coefficients[4,2],        # The standard error of the coefficient
                          cov_balance,                                # The covariate balance
                          parralel_trend,                                 # Parralel trend info
                          #                      Gen=as.logical(F), NN1=as.logical(F), ebal=as.logical(F), psm=as.logical(T), lm=as.logical(F),
                          var,                                             # The extra covariates included
                          cal025=ifelse(cali==0.25,as.logical(T),as.logical(F)), cal05=ifelse(cali==0.5,as.logical(T),as.logical(F)),cal1=ifelse(cali==1,as.logical(T),as.logical(F)),
                          replacement=ifelse(as.logical(cov[2])==T,as.logical(T),as.logical(F)),
                          maha=ifelse(cov[1]=="mahalanobis",as.logical(T),as.logical(F)), glm=ifelse(cov[1]=="glm",as.logical(T),as.logical(F)), rf=ifelse(cov[1]=="randomforest",as.logical(T),as.logical(F)),
                          #                               NN=ifelse(meth=="nearest",as.logical(T),as.logical(F)), genetic=ifelse(meth=="genetic", as.logical(T), as.logical(F)),
                          ratio1=ifelse(as.numeric(cov[3])==1,as.logical(T),as.logical(F)),ratio5=ifelse(as.numeric(cov[3])==5,as.logical(T),as.logical(F)),ratio10=ifelse(as.numeric(cov[3])==10,as.logical(T),as.logical(F))
      )
      }             
      
    }
    
  })
  
  synthesis_loop_cov_param <<- data.frame(do.call("rbind",maha))
  }  



the_year_loop <- function(offset, ycreation) {
  
  # We store in a string the name of the offset
  name_offset = offset$offset[1]
  
  #We pre-clean the data to speed up the computations by manually getting rid off observations that will be too far when using calipers
  offset_l <- offset %>% na_if(-9999)
  
  # We keep observations without na
  offset_l <- offset_l %>% drop_na()
  
  
  # We determine the caliprs for each variables
  caliper_value = 1
  caliper_slope<-caliper_value*sd(offset_l$Slope[offset_l$treated==1])
  caliper_elev<-caliper_value*sd(offset_l$Elevation[offset_l$treated==1])
  caliper_dist_road<-caliper_value*sd(offset_l$Dist_road[offset_l$treated==1])
  caliper_dist_edge<-caliper_value*sd(offset_l$Dist_edge[offset_l$treated==1])
  caliper_dist_def<-caliper_value*sd(offset_l$Dist_defor[offset_l$treated==1])
  
  # We keep control observations that are on the same support
  offset_l <- offset_l[offset_l$Slope > min(offset_l$Slope[offset_l$treated==1])-caliper_slope,]
  offset_l <- offset_l[offset_l$Slope < max(offset_l$Slope[offset_l$treated==1])+caliper_slope,]
  
  offset_l <- offset_l[offset_l$Elevation > min(offset_l$Elevation[offset_l$treated==1])-caliper_elev,] 
  offset_l <- offset_l[offset_l$Elevation < max(offset_l$Elevation[offset_l$treated==1])+caliper_elev,]  
  
  offset_l <- offset_l[offset_l$Dist_road > min(offset_l$Dist_road[offset_l$treated==1])-caliper_dist_road,]
  offset_l <- offset_l[offset_l$Dist_road < max(offset_l$Dist_road[offset_l$treated==1])+caliper_dist_road,] 
  
  offset_l <- offset_l[offset_l$Dist_edge > min(offset_l$Dist_edge[offset_l$treated==1])-caliper_dist_edge,] 
  offset_l <- offset_l[offset_l$Dist_edge < max(offset_l$Dist_edge[offset_l$treated==1])+caliper_dist_edge,] 
  
  offset_l <- offset_l[offset_l$Dist_defor > min(offset_l$Dist_defor[offset_l$treated==1])-caliper_dist_def,] 
  offset_l <- offset_l[offset_l$Dist_defor < max(offset_l$Dist_defor[offset_l$treated==1])+caliper_dist_def,]
  
  
  # We only keep the ppol of covaiates, treatment index and outcome
  db<-offset_l[,c(4:7,9:14)]
  treated <- offset_l$treated
  treeloss <- offset_l$Tree_loss
  
  
  
  # We define the set of covariates
  covariates <- names(db[,c(3,4,7,9,10)])  # A set of minimal covariates that should be included in any matching

  
  
  # We keep the good covariates
    
    idx <- match(covariates, names(db))   # We determines the columns of the covariates
    idx <- sort(idx)               # We need to arrange them in increasing order
    NewDF <- db[,c(idx)]           # We create a database that contains only the covariates needed in this iteration of the loop
    
    formulae  <- as.formula(paste("treated ~ ", paste0(names(NewDF),collapse ="+")))  # We create a formula to be inserted in the matching
    
    
    
    # Calipers
    cal=rep(caliper_value,length(covariates)) # In the version 4.0 of matchIT, we now need to put calipers for each variables included in the model
    names(cal) <- covariates         # We name the vestor to be inserted in the matching function
    
    
    # Applying the matching function to 
    m.out <-  matchit(formulae, data= cbind(treated,NewDF,treeloss), 
              method= "nearest", 
              distance = "mahalanobis", 
              replace = F, 
              caliper = cal)
    
    
      # Summarise results #
      match_output <-summary(m.out, standardize = TRUE)
      
      # Create db
      m.data <- match.data(m.out)
      
      
      # Aggregate pixels into treatment (offset_l_l) and control groups. Tabulate the number of pixels by tree loss year within each group (count of number of pixels within each group deforested each year)
      annual_defor <- m.data %>%  # As weights are created with the replace or with PSM, I change the tabulate approach of base R by a dplyr count
        count(as.factor(treated), as.factor(treeloss), wt=weights, .drop = F)
      
      label <- c("Sample", "Year", "Annual_Deforestation")
      names(annual_defor) <- paste0(label)
      
      
      # Have to make Year numeric for >= to work
      annual_defor$Year <- as.numeric(as.character(annual_defor$Year))
      
      
      # We delete the 0 values (aka pixels that remain forested at the end)
      annual_defor <- annual_defor[annual_defor$Year!=0,]
      
      # Construct data with dummy variables for time and treated to use in DiD regression
      annual_defor$Time <- ifelse(annual_defor$Year >= ycreation, 1,0)
      annual_defor$TimeF <-factor(annual_defor$Time, levels = c(0,1), labels = c("before","after"))
      annual_defor$TreatedF <- factor(annual_defor$Sample, levels = c(0,1), labels = c("control","treatment"))
      
      
      
      # d) Loop over DiD Regression
      

did_res <- data.frame(coef=numeric(),      # The coefficient of the estimation
                  se=numeric(),
                  Year = character())
      
#for (y in c(1:19)) {

annual_defor$Year <- annual_defor$Year+2000

for (y in 2001:2019) {
annual_defor2 <- annual_defor[annual_defor$Year != y,]
i = y-2000  
      modeldid <- lm(log(Annual_Deforestation+1) ~ TreatedF*TimeF, data= annual_defor2)
      did_res[i,1] <- summary(modeldid)$coefficients[4,1]      # The coefficient of the estimation
      did_res[i,2] <- summary(modeldid)$coefficients[4,2]
      did_res[i,3] <- as.character(y)
  }



# We add the regression without any drop
modeldid <- lm(log(Annual_Deforestation+1) ~ TreatedF*TimeF, data= annual_defor)
summary(modeldid)     
did_res[20,1] <- summary(modeldid)$coefficients[4,1]      # The coefficient of the estimation
did_res[20,2] <- summary(modeldid)$coefficients[4,2]
did_res[20,3] <- "None"

# We create a variable of whether there is a drop or not for the graphs
did_res$Drop <- ifelse(did_res$Year=="None",F,T)

did_res <<- did_res
}
