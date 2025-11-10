* install the repframe package
net install repframe, from("https://raw.githubusercontent.com/guntherbensch/repframe/main") replace

* load the data and use the first row as column names
* this Excelfile comes from the R-code for the specification curve, in the same folder
import excel "D:/Solene/Figs_Dashboard_SpecCurve_replication/formatted_results.xlsx", clear firstrow

* create the original path indicator (origpath)
gen origpath = (specification_n == 0)  // If specification_n is 0, it's the original analysis

* Replace the value of outcome to 0 when extension = 0
replace outcome = 0 if extension == 0
replace outcome = 1 if extension == 1


*** duplicate origpath - we make here 2 "outcomes" but those are just robustness vs. extension
preserve

keep if origpath == 1

replace outcome = 1
replace extension = 1

tempfile temp_rows
save `temp_rows', replace

restore

* append the saved rows back to the dataset (duplicate)
append using `temp_rows'


***** separate robustness vs extension to have 2 "outcomes"
* Define the label for outcome variable
label define outcome_lbl 0 "Effects on deforestation (robustness checks)" ///
                         1 "Effects on deforestation (dataset extensions)"


* Apply the label to the outcome variable
label values outcome outcome_lbl

* Rename columns for better compatibility with repframe
rename beta b
rename se se
rename pval p




* Create significance level for the original analysis (siglevel_orig)
gen siglevel_orig = cond(pval_orig < 0.05, 5, 10)  // Assuming pval_orig is the original p-value

* Create significance level for the robustness analysis (siglevel)
gen siglevel = cond(p < 0.05, 5, 10)  // Assuming p is the p-value for the robustness analysis


* create graph and save
repframe outcome, beta(b) se(se) pval(p) origpath(origpath) siglevel_orig(5) siglevel(5) shortref(extension) 
graph export "D:/Solene/Figs_Dashboard_SpecCurve_replication/plots/dashboard_extension_robustness.png", replace









