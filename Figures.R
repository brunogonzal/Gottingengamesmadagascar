

#                                 1) Covariate Balance plot

cov_balance_plot <- function(input, colour1, name){
  result <- ggplot(input, aes(x= factor(sample, level= c("Non-matched", "matched")), y=SMD_1_cal, group= covariates))+
    geom_line(color= colour1)+
    geom_text_repel(data=subset(input, sample=="Non-matched"), label= variables, hjust=1.3, direction="y", size = 3.5, colour = "grey15", segment.color = colour1)+
    geom_point(size=2.5, colour = colour1) +
    annotate("rect", xmin = 0, xmax = Inf, ymin = -0.25, ymax = 0.25,alpha = .2)+
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40")+
    scale_x_discrete(name="Sample")+
    scale_y_continuous(name="Standardised Difference in Means", breaks = round(seq(min(input[,3]), max(input[,3]) , 0.2), digits = 1))+
    #scale_y_continuous(name="Standardised Difference in Means", breaks = seq(-1.4, 0.8, 0.2), labels = scales::comma)+
    theme_classic()+
    theme(legend.title = element_blank(), axis.title.x = element_text(size=10, face="bold", vjust= 0.5), axis.title.y= element_text(size=10, face="bold", vjust = 3),
          axis.text.x = element_text(face="bold", color = "grey15"), axis.text.y=element_text(face="bold", color = "grey15"))+
    ggtitle(name)
  return(result)
  }


input <- cov_balance_TTF2
colour1 <- "#56B4E9"
name <- "Torotorofotsy"


cov_balance_plot(cov_balance_ANK2, "#E69F00", "Ankerana")
cov_balance_plot(cov_balance_CFAM2, "chartreuse3", "CFAM")
cov_balance_plot(cov_balance_CZ2,"#F0E442","Conservation Zone")
cov_balance_plot(cov_balance_TTF2, "#56B4E9", "Torotorofotsy")



# ----------------------------------------------------------------------------------------------------------------------------------------------------------------#


#                                  2) Percentage Annual Deforestation over whole study period 

defor_plot2 <- function(a,c,d){
  ggplot(a, aes(x=Year, y=Perc_Annual_Defor, group= Sample, colour=Sample))+
    geom_point(size= 3.5, shape = 18)+
    #geom_line(size=0.75)+
    geom_vline(xintercept = c, linetype="dashed", color = "brown", size= 0.75)+
    scale_color_manual(values=c("#CC3300","#FFCC33" ))+
    scale_y_continuous(name="Annual Deforestation Rate")+
    scale_x_continuous(name= "Year", breaks = seq(1,19,2), labels = c("2001", "2003", "2005", "2007", "2009", "2011", "2013", "2015", "2017", "2019"))+
    theme_classic()+
    theme(legend.title = element_blank(), axis.title.x = element_text(size=10, face="bold", vjust=0.5), axis.title.y= element_text(size=10, face="bold", vjust = 3),
          axis.text.x = element_text(face="bold", color = "grey15"), axis.text.y=element_text(face="bold", color = "grey15"))+
    ggtitle(d)
}


defor_plot2(annual_defor_ANK, 11, "Ankerana")
defor_plot2(annual_defor_CFAM, 13, "CFAM")
defor_plot2(annual_defor_CZ, 9, "Conservation Zone")
defor_plot2(annual_defor_TTF, 14, "Torotorofotsy")

ggplot(annual_defor_ANK, aes(x=Year, y=Perc_Annual_Defor, group= Sample, colour=Sample))+
  geom_point(size= 3.5, shape = 18)+
  geom_vline(xintercept = 11, linetype="dashed", color = "grey25", size= 0.75)+
  scale_color_manual(values=c("#E69F00","grey50"))+
  scale_y_continuous(name="Annual Deforestation Rate (%)", breaks = seq(0, 4, 1))+
  scale_x_continuous(name= "Year", breaks = seq(1,19,2), labels = c("2001", "2003", "2005", "2007", "2009", "2011", "2013", "2015", "2017", "2019"))+
  theme_classic()+
  theme(legend.title = element_blank(), axis.title.x = element_text(size=10, face="bold", vjust=0.5), axis.title.y= element_text(size=10, face="bold", vjust = 3),
        axis.text.x = element_text(face="bold", color = "grey15"), axis.text.y=element_text(face="bold", color = "grey15"))+
  ggtitle("Ankerana")+
  ggsave("Graphics/ANK_defor_rate_final.jpeg")

ggplot(annual_defor_CFAM, aes(x=Year, y=Perc_Annual_Defor, group= Sample, colour=Sample))+
  geom_point(size= 3.5, shape = 18)+
  geom_vline(xintercept = 13, linetype="dashed", color = "grey25", size= 0.75)+
  scale_color_manual(values=c("chartreuse3", "grey50"))+
  scale_y_continuous(name="Annual Deforestation Rate (%)", breaks = seq(0, 7, 1))+
  scale_x_continuous(name= "Year", breaks = seq(1,19,2), labels = c("2001", "2003", "2005", "2007", "2009", "2011", "2013", "2015", "2017", "2019"))+
  theme_classic()+
  theme(legend.title = element_blank(), axis.title.x = element_text(size=10, face="bold", vjust=0.5), axis.title.y= element_text(size=10, face="bold", vjust = 3),
        axis.text.x = element_text(face="bold", color = "grey15"), axis.text.y=element_text(face="bold", color = "grey15"))+
  ggtitle("CFAM")+
  ggsave("Graphics/CFAM_defor_rate_final.jpeg")

ggplot(annual_defor_CZ, aes(x=Year, y=Perc_Annual_Defor, group= Sample, colour=Sample))+
  geom_point(size= 3.5, shape = 18)+
  geom_vline(xintercept = 9, linetype="dashed", color = "grey25", size= 0.75)+
  scale_color_manual(values=c("#F0E442", "grey50"))+
  scale_y_continuous(name="Annual Deforestation Rate (%)", breaks = seq(0, 5, 1))+
  scale_x_continuous(name= "Year", breaks = seq(1,19,2), labels = c("2001", "2003", "2005", "2007", "2009", "2011", "2013", "2015", "2017", "2019"))+
  theme_classic()+
  theme(legend.title = element_blank(), axis.title.x = element_text(size=10, face="bold", vjust=0.5), axis.title.y= element_text(size=10, face="bold", vjust = 3),
        axis.text.x = element_text(face="bold", color = "grey15"), axis.text.y=element_text(face="bold", color = "grey15"))+
  ggtitle("CZ")+
  ggsave("Graphics/CZ_defor_rate_final.jpeg")

ggplot(annual_defor_TTF, aes(x=Year, y=Perc_Annual_Defor, group= Sample, colour=Sample))+
  geom_point(size= 3.5, shape = 18)+
  geom_vline(xintercept = 14, linetype="dashed", color = "grey25", size= 0.75)+
  scale_color_manual(values=c("#56B4E9", "grey50"))+
  scale_y_continuous(name="Annual Deforestation Rate (%)", breaks = seq(0, 12.5, 1))+
  scale_x_continuous(name= "Year", breaks = seq(1,19,2), labels = c("2001", "2003", "2005", "2007", "2009", "2011", "2013", "2015", "2017", "2019"))+
  theme_classic()+
  theme(legend.title = element_blank(), axis.title.x = element_text(size=10, face="bold", vjust=0.5), axis.title.y= element_text(size=10, face="bold", vjust = 3),
        axis.text.x = element_text(face="bold", color = "grey15"), axis.text.y=element_text(face="bold", color = "grey15"))+
  ggtitle("Torotorofotsy")+
  ggsave("Graphics/TTF_defor_rate_final.jpeg")



# ----------------------------------------------------------------------------------------------------------------------------------------------------------------#


#                                 3) Parallel trends     

# Plot annual defor in before period in offset and matched control samples
# Add the significant predictors from the parallel trends regression

# -----ANK --------#
# No significant predictors - just plot points:
                                        
ANKa_plot <- ggplot(ANK_data_before, aes(x=Year, y=log_annual_defor, group = Sample, color= Sample))+
  geom_point(size= 3.5, shape = 18)+
  scale_color_manual(values=c("#E69F00","grey50"))+
  scale_y_continuous(name = "log( deforestation+1 )", breaks = seq(0,4,0.5))+ 
  scale_x_continuous(name= "Year", breaks = seq(1,13,2), labels = c("2001", "2003", "2005", "2007", "2009", "2011", "2013"))+
  expand_limits(x= 13) +
  theme_classic()+
  theme(legend.title = element_blank(), axis.title.x = element_text(size=10, face="bold", vjust=0.5), axis.title.y= element_text(size=10, face="bold", vjust = 3),
        axis.text.x = element_text(face="bold", color = "grey15"), axis.text.y=element_text(face="bold", color = "grey15"))+
  ggtitle("Ankerana")+
  ggsave("ANK_pt_fit_FINAL.jpeg")


# -----TTF ---------#

summary(TTFa)        # Only Year is significant predictor. 
#No significant difference in the relationship between Year and log annual defor between control and treated samples.

# Plot the significant line from the parallel trends regression lm(log_annual_defor ~ Year*TreatedF, data= TTF_data_before)
# Add the fitted values from the parallel trends regression to TTF_data_before dataframe as this supplies the points to be plotted.

TTF_data_before$Year_fit <- 0

# Only interested in the fitted values for the Year line for the control group (as there is no sig. difference between treated and control.
# This is the first 13 values

TTF_data_before$Year_fit <- rep(TTFa$fitted.values[1:13])

# Add columns for Upper and Lower Confidence intervals

TTF_data_before$Upper <- 0
TTF_data_before$Lower <- 0

# Use predict function to get the confidence intervals around the fitted values. 

TTF_new_dat <- data.frame(Year = c(1:13), TreatedF = 0)                # Tell predict we want to predict values for Years 1-13 for the control group (as this is the significant trend).
TTF_new_dat$TreatedF <- factor(TTF_new_dat$TreatedF, levels = c(0,1))  # Make TreatedF a factor

TTF_cont_pred <- data.frame(predict(TTFa, newdata = TTF_new_dat, interval = 'confidence'))  

# Add the predicted confidence intervals around the fitted values to the dataframe
# Because there is no significant difference in the slope of the relationship between Year and Deforestation between the control and treated samples
# I want the values for the treated sample to be the same. (NB. this is why I didn't use the whole set of fitted values from the model because they
# are slightly different for the treated sample, however this difference is not significant)

TTF_data_before[, c(9,10)] <- rep(TTF_cont_pred[, c(2,3)])

names(TTF_data_before)[c(9,10)] <- c("Lower", "Upper")   # Re-order column names

TTFa_plot <- ggplot(TTF_data_before, aes(x=Year, y=log_annual_defor, group = Sample, color= Sample))+
  geom_point(size=3.5, shape = 18)+
  scale_color_manual(values=c("TTF" = "#56B4E9","Cont" = "grey50"))+
  geom_line(aes(x = Year, y = Year_fit), colour = "navy", size = 1)+                    # Plot fitted values Year and Upper and Lower CIs
  geom_line(aes(x = Year, y = Lower), linetype = "dashed", colour = "navy", size = 1)+
  geom_line(aes(x = Year, y = Upper), linetype = "dashed", colour = "navy", size = 1)+
  scale_y_continuous(name = "log( deforestation+1 )", breaks = seq(0,4.5,0.5))+ 
  scale_x_continuous(name= "Year", breaks = seq(1,13,2), labels = c("2001", "2003", "2005", "2007","2009", "2011", "2013"))+
  theme_classic()+
  theme(legend.position = "none", axis.title.x = element_text(size=10, face="bold", vjust=0.5), axis.title.y= element_text(size=10, face="bold", vjust = 3),
        axis.text.x = element_text(face="bold", color = "grey15"), axis.text.y=element_text(face="bold", color = "grey15"))+
  ggtitle("Torotorofotsy")


# ----- CZ -------#

# For the Conservation Zone only the Treated term is significant. This means there is no significant trend 
# in deforestation over time at the 95% confidence interval in either the control or the treated group. However, on average
# deforestation is lower in the treated sample than the control sample (because the coefficent of TreatedF1 is negative)

# Add the confidence intervals of the Treated estimate to the results.

CZa_res <- tidy(summary(CZa))
CZa_res <- cbind(CZa_res, confint(CZa))
names(CZa_res)[c(6,7)] <- c("Lower", "Upper")

# Intercept = the average deforestation in the pre-intervention period for the control sample 
# Upper and Lower confidence intervals of this estimate provided by confint

# Intercept + Treated coefficient gives the average deforestation for the treated sample (1.42 +- 1.599)
# Intercept + Upper Confidence interval of the estimated Treated coefficient gives the upper confidence interval of average deforestation in the treated sample.
# Intercept + Lower confidence interval of the estimated Treated coefficient gives the lower confidence interval of the average deforestation in the treated sample

# When plotted the two confidence intervals overlap but because the upper confidence interval of the Treated Sample is below the mean estimate of the Control Sample
# solid grey line, the difference is still significant.# 


CZa_plot <- ggplot(CZ_data_before, aes(x=Year, y=log_annual_defor, group = Sample, color= Sample))+
  geom_point(size=3.5, shape = 18)+
  geom_segment(x= 1, xend = 8, y = CZa_res[1,2], yend = CZa_res[1,2],  color = "grey50", size = 1)+     # Intercept
  geom_segment(x= 1, xend = 8, y =  CZa_res[1,2] + CZa_res[3,2], yend = CZa_res[1,2] + CZa_res[3,2], color = "#F0E442", size= 1)+  # Intercept - estimate of TreatedF1.
  
  geom_segment(x= 1, xend = 8, y = CZa_res[1,6], yend = CZa_res[1,6],  color = "grey50", size = 1, linetype = "dashed")+     # Lower CI for control 
  geom_segment(x= 1, xend = 8, y =  CZa_res[1,7], yend = CZa_res[1,7], color = "grey50", size= 1, linetype = "dashed")+  # Upper CI for control
  
  geom_segment(x= 1, xend = 8, y =  CZa_res[1,2] + CZa_res[3,6], yend = CZa_res[1,2] + CZa_res[3,6], color = "#F0E442", size= 1, linetype = "dashed")+  # Lower CI for treated (Intercept + lower CI of treated estimate) 
  geom_segment(x= 1, xend = 8, y =  CZa_res[1,2] + CZa_res[3,7], yend = CZa_res[1,2] + CZa_res[3,7], color = "#F0E442", size= 1, linetype = "dashed")+  # Upper CI for treated (Intercept + upper CI of treated estimate)
  
  scale_color_manual(values=c("#F0E442","grey50" ))+
  scale_y_continuous(name = "log( deforestation+1 )", breaks = seq(-0.5,3.5,0.5))+ 
  scale_x_continuous(name= "Year", breaks = seq(1,13,2), labels = c("2001", "2003", "2005", "2007","2009", "2011", "2013"))+
  expand_limits(x=13, y = -1.5)+
  theme_classic()+
  theme(legend.title = element_blank(), axis.title.x = element_text(size=10, face="bold", vjust=0.5), axis.title.y= element_text(size=10, face="bold", vjust = 3),
        axis.text.x = element_text(face="bold", color = "grey15"), axis.text.y=element_text(face="bold", color = "grey15"))+
  ggtitle("Conservation Zone")


# ------ CFAM --------#

summary(CFAMa)

# Year and the interaction between Year and Treated status are significant. This means that there is a significant trend in deforestation over time in the control sample (Year)
# and the treated sample, but the slope of this relationship differs significantly between the two samples. Need to plot both the time trend for each sample

# Extract fitted values from regression and append to data frame used for plotting

CFAM_data_before$Year_fit <- CFAMa$fitted.values
CFAM_data_before$Upper <- 0
CFAM_data_before$Lower <- 0

# Use predict to get the confidence intervals around the fitted values.
# Predict deforestation each year for the control (0) and treated sample (1) using the regression and get the upper and lower confidence
# intervals around these estimates.

new_dat_cont <- data.frame(Year = c(1:12), TreatedF = 0)
new_dat_cont$TreatedF <- factor(new_dat_cont$TreatedF, levels = c(0,1))

new_dat_treat <- data.frame(Year = c(1:12), TreatedF = 1)
new_dat_treat$TreatedF <- factor(new_dat_treat$TreatedF, levels = c(0,1))

CFAM_cont_pred <- data.frame(predict(CFAMa, newdata= new_dat_cont, interval = 'confidence'))
CFAM_treat_pred <- data.frame(predict(CFAMa, newdata= new_dat_treat, interval = 'confidence'))

# Append fitted values plus upper and lower CIs generated using predict to the dataframe used for plotting. 

CFAM_data_before[c(1:12), c(9,10)] <- CFAM_treat_pred[,c(2,3)]     # Treated values
CFAM_data_before[c(13:24), c(9,10)] <- CFAM_cont_pred[,c(2,3)]     # Control values

CFAMc_plot <- ggplot(CFAM_data_before, aes(x=Year, y=log_annual_defor, group = Sample, color= Sample))+
  geom_point(size=3.5, shape = 18)+
  geom_line(aes(x = Year, y = Year_fit), size =1)+
  geom_line(aes(x= Year, y = Lower), linetype = "dashed", size =1)+
  geom_line(aes(x= Year, y = Upper), linetype = "dashed", size = 1)+
  scale_color_manual(values=c("chartreuse3","grey50"))+
  scale_y_continuous(name = "log( deforestation+1 )", breaks = seq(0,5,0.5))+ 
  scale_x_continuous(name= "Year", breaks = seq(1,13,2), labels = c("2001", "2003", "2005", "2007","2009", "2011", "2013"))+
  expand_limits(x=13)+
  theme_classic()+
  theme(legend.title = element_blank(), axis.title.x = element_text(size=10, face="bold", vjust=0.5), axis.title.y= element_text(size=10, face="bold", vjust = 3),
        axis.text.x = element_text(face="bold", color = "grey15"), axis.text.y=element_text(face="bold", color = "grey15"))+
  ggtitle("CFAM")



# ----------------------------------------------------------------------------------------------------------------------------------------------------------------#


#                                 4) Plot ATT from individual DiD and overall fixed effects


# Add fixed effects results to ATT_all

ATT_all[4,] <- 0
rownames(ATT_all)[4] <- "Overall"
ATT_all[4,1] <- (exp(coef(FE_all))-1)*100
ATT_all[4, 2:3] <- (exp(confint(FE_all))-1)*100
ATT_all$Sample <- factor(c("ANK", "CZ", "TTF", "Overall"), levels = c("ANK", "CZ", "TTF", "Overall"))


cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442")
#               grey      orange      sky blue  green       yellow

#                                               "chartreuse3" <- this is a much better green for CFAM

# Re-order for bar plot


ggplot(ATT_all, aes(x = Sample, y= ATT, fill = Sample))+
  geom_col(colour = "black", width = c(0.525, 0.27, 0.2, 1.5))+
  geom_errorbar(aes(ymin = Lower_CI, ymax = Upper_CI), width = c(0.4, 0.15, 0.1, 0.9))+
  scale_fill_manual(values = c("#E69F00", "#F0E442", "#56B4E9", NA))+
  scale_x_discrete(name = "Offset")+
  coord_cartesian(ylim=c(-150, 150))+
  scale_y_continuous(name = "Treatment Effect (% difference in deforestation)", breaks = seq(-150, 150, 50))+
  geom_hline(yintercept = 0, linetype = "dashed")+
  theme_classic()+
  theme(legend.position = "none", axis.title.x = element_text(size=11, vjust=1), axis.title.y= element_text(size=11, vjust = 1),
        axis.text.x = element_text(size = 11, color = "grey15"), axis.text.y=element_text(size= 11, color = "grey15"))

1.5* 0.13


# ----------------------------------------------------------------------------------------------------------------------------------------------------------------#


#                               5) Plot observed vs. counterfactual deforestation

rm(a)

a <- data.frame(t(Impact_defor[, c(2,3,6)]))
a[1:12, 1] <- data.frame(c(a[, "X10"], a[, "X12"], a[, "X7"], a[, "X1"]))
a$X <- rep(c("Observed", "Counterfactual", "Avoided Deforestation"),4)
a <- a[, c(1,5)]
a$Offset <- rep(c("Ankerana", "Conservation Zone", "Torotorofotsy", "Overall"), each=3)
a <- a[, c(3,2,1)]
names(a) <- c("Offset", "X", "Total_Defor_Ha")
a$Upper_CI <- 0
a$Lower_CI <- 0
a$Upper_CI[c(2,5,8,11)] <- Impact_defor$Counterfactual_Defor_Upper
a$Lower_CI[c(2,5,8,11)] <- Impact_defor$Counterfactual_Defor_Lower
a$Upper_CI[c(3,6,9, 12)] <- Impact_defor$Impact_defor_Upper
a$Lower_CI[c(3,6,9, 12)] <- Impact_defor$Impact_defor_Lower
a$Offset <- factor(a$Offset, levels = c("Ankerana", "Conservation Zone", "Torotorofotsy", "Overall"))

a$X <- factor(a$X, levels = c("Observed", "Counterfactual", "Avoided Deforestation"))

line <- data.frame(Y = -2064, Offset = "Overall")
line$Offset <- factor(line$Offset, levels = c(NA, NA, NA, "Overall"))

ggplot(a, aes(x = X, y= Total_Defor_Ha, fill = Offset))+
  geom_bar(stat = "identity", width = 0.9, colour = "black")+
  geom_errorbar(aes(ymin = Lower_CI, ymax = Upper_CI), width = 0.4)+
  geom_hline(yintercept = 0, linetype = "dashed", colour = "red", size = 0.71)+
  geom_hline(data = line, aes(yintercept = Y), linetype = "dashed", colour = "chartreuse 4", size = 0.75)+
  facet_wrap(~ Offset)+
  scale_fill_manual(values = c("#E69F00", "#F0E442", "#56B4E9","#999999"))+
  scale_x_discrete(labels= c("Observed", "Counterfactual", "Avoided\nDeforestation"))+
  scale_y_continuous(name = "Total deforestation after protection (ha)", limits = c(-5500, 5500), breaks = seq(-5000, 5000, 2500))+
  theme_classic()+
  theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y= element_text(size=11, vjust = 2),
        axis.text.x = element_text(size = 11, color = "grey15", angle = 90, hjust = 0.5), axis.text.y=element_text(size= 11, color = "grey15"),
        )
  
line <- data.frame(Y = -2064, Offset = "Overall")
line$Offset <- factor(line$Offset, levels = c(NA, NA, NA, "Overall"))


# ------------------------------------------------------------------------------------------------------------------------------------#


                                        # 6) Comparing our Cohen's normalised effect size to results from Borner et al


cohens_comp <- read.csv("C:/Users/ktd19ycv/OneDrive - Bangor University/Documents/Industrial Mining/Ambatovy offset analysis/Paper prep/Revisions/cohens_d_comparison.csv")

cohens_comp <- cohens_comp[c(137:140,1:136),]

cohens_comp$ID <- 1:nrow(cohens_comp)

cohens_comp$Intervention_group[cohens_comp$Intervention_group == "Land Titling and Reform"] <- "LTR"

cohens_comp$Intervention_group <- factor(cohens_comp$Intervention_group)

cohens_comp$Intervention_group[1:4] <- "Ambatovy's offsets"

cohens_comp <- cohens_comp[-4,]

cohens_comp$Colour <- "None"
cohens_comp$Colour[1:3] <- c("Grey", "Orange", "Yellow")
cohens_comp$Colour <- factor(cohens_comp$Colour, levels = c("Grey", "Orange", "Yellow", "None"))

ggplot(cohens_comp, aes(x = Intervention_group, y = Cohens, group = Colour, colour = Colour))+
  geom_hline(yintercept = 0, linetype = "dashed", colour = "red", size = 0.6)+
  geom_point(size = 2.8, shape = 18)+
  scale_colour_manual(values = col_pal)+
  theme_classic()+
  scale_y_continuous(name = "Cohen's d effect size")+
  scale_x_discrete(name = "Type of intervention")+
  theme(legend.position = "none", axis.title.x = element_text(size=10, face="bold", vjust = 3), axis.title.y= element_text(size=10, face="bold", vjust = 3), 
        axis.text.x = element_text(angle = 45, hjust = 0.95, vjust = 0.95, color = "grey15"), axis.text.y=element_text(color = "grey15"))

col_pal <- c("#999999", "#FF9900", "gold", "grey15")


