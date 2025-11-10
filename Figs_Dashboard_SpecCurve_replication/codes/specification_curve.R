# Replication Games
#
# Deforestation paper:
# Devenish, K., Desbureaux, S., Willcock, S. et al. 
# On track to achieve no net loss of forest at Madagascar’s biggest mine. 
# Nat Sustain 5, 498–508 (2022). https://doi.org/10.1038/s41893-022-00850-7

# Specification curve after replication

##################################################################################
# Set up
##################################################################################

setwd("D:/Solene/Figs_Dashboard_SpecCurve_replication")

if (!require("r2especcurve")) remotes::install_github('NickCH-K/r2especcurve') # install package required for spec curve

if (!require("pacman")) install.packages("pacman")

pacman::p_load(
  tidyverse,
  dplyr,
  readr,
  ggplot2,
  r2especcurve,
  readxl
)


##################################################################################
# Pooled Effects - Spec Curve
##################################################################################

pooled <- read_xlsx("D:/Solene/Figs_Dashboard_SpecCurve_replication/Table1_pooled_TWFE_regression_RobChecks.xlsx", skip = 1)

colnames(pooled) <- c("Specification", "What", "n", "beta", "se", "t_value", "pval", "CI_lower", "CI_upper")

pooled <- pooled %>%
  mutate(extension = ifelse(grepl("Extended", Specification), 1, 0)) %>%
  mutate(specification_n = row_number()-1)

pooled$outcome <- 1
pooled$n_orig <- 152
pooled$beta_orig <- -0.878
pooled$se_orig <- 0.216
pooled$pval_orig <- 0.000

# Select and print the 'specification_n' and 'Specification' columns
print(pooled %>% select(specification_n, Specification))

# Alternatively, view the selected columns
head(pooled %>% select(specification_n, Specification))



# Example: calling spec_curve function
result <- spec_curve(
  pooled[pooled$outcome == 1,], 
  decision_cols = c('specification_n', 'extension'),
  decision_labels = c('Specification', 'Extension to 2023 (1)'),
  orig_spec_values = c(0, 0),
  return_separate = TRUE
)

# Now you have two plots: result[[1]] (p1) and result[[2]] (p2)
p1 <- result[[1]]
p2 <- result[[2]]

# Modify text properties for the p1 plot
p1 <- p1 +
  ggplot2::theme(
  ) +
  ggplot2::labs(
    x = "",
    y = "Reduction\nin deforestation",   # Change title to "Reduction in deforestation"
    caption = '90% and 95% confidence intervals shown.',
    fill = NULL
  )

p1

# Modify text properties for the p2 plot (labels)
p2 <- p2 +
  ggplot2::theme_void() +
  ggplot2::theme(
    text = ggplot2::element_text(face = "bold"),    # All text bold and size 16
    axis.text.y = ggplot2::element_text(face = "bold"),  # Y-axis label text
    #axis.text.x = ggplot2::element_text(face = "bold")   # X-axis label text
  )

p2

# Combine both plots
p1/p2 + patchwork::plot_layout(widths = c(1, 1), heights = c(4, 1))


ggsave("D:/Solene/Figs_Dashboard_SpecCurve_replication/plots/specification_curve.png") #, width=16, height=14, dpi=300)

writexl::write_xlsx(pooled, "D:/Solene/Figs_Dashboard_SpecCurve_replication/formatted_results.xlsx")