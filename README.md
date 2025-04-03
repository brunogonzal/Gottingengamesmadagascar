# Biodiversity_offset_effectiveness
This repository contains code and input data used in the paper under review 'Madagascar's biggest mine is on track to achieve No Net Loss of Forest' by Katie Devenish, Sebastien Desbureaux, Simon Willcock and Julia PG Jones. 

# Input data

- Final_control.dbf
- Sample_TTF3.dbf
- Sample_ANK3.dbf
- Sample_CFAM3.dbf
- Sample_CZ3.dbf

These are the input data for the matching and outcome regressions. Each observation contains covariate and outcome variable values at sample pixel locations (denoted in X,Y co-ordinates) from the four offsets and the control area. The outcome variable is Tree Loss Year from the Global Forest Change dataset of Hansen et al (2013) and refers to the year in which a pixel was deforested between 2001 and 2019. If the pixel was not deforested during this period the value is set as 0.

- ANK_var.dbf
- CFAM_var.dbf
- CZ_var.dbf
- Torotorofotsy_var.dbf

These database files contain forest cover and annual tree loss data for the whole area of each biodiversity offset. These data were used to quantify the avoided deforestation following protection. 

# Scripts

1) Choosing_matching_spec - details the alternative matching specifications tested for covariate balance before choosing the main specification. 
2) Offset_Effectiveness_Final2 - code for the matching and outcome regressions using the main matching specification. This produces the main results presented in the paper.
3) Figures - code used to produce the figures presented in the paper
4) Robustness_functions_v4_withyearloop -
5) Offset_effectiveness_slim_lopp_func4 - 
6) FE_loop - 4), 5) and 6) are the loops for the robustness checks
9) Adding_time_trend - tests the effect of including a time trend in the difference-in-differences regression. 
10) Pop_density - shows the correlation of Population Density with the 5 essential covariates. 
11) Impact_modelling_choices - regression to explore which modelling choices have the greatest influence on the estimated treatment effect using the 496 estimates from the robustness checks.
12) 2B_AOB_spec_graph_modifiedhighlight - This is the code for the specification graphs developed by Ariel Ortiz-Bobea used to plot the results from our robustness checks (Ortiz-Bobea, A., Ault, T. R., Carrillo, C. M., Chambers, R. G. & Lobell, D. B. Anthropogenic climate change has slowed global agricultural productivity growth. doi:10.1038/s41558-021-01000-1.)


