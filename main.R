# nolint: object_name_linter
#installation of necessary packages
source("utils.R")
source("visualization.R")
library(tidyverse)
library(ggplot2)



#import data
raw_df <- read.csv(PATH)


#Data Visualization
# For continuous variable (tarv_1)
plots1 <- plot_continuous_distribution(data, "your_var_name", "strata_column")
plots1$density  # View density plot
plots1$boxplot  # View boxplot
plots1$stratified  # View stratified boxplot

# For proportion variable (tarv_2)
plots2 <- plot_proportion_distribution(data, "your_var_name", threshold = 100, "strata_column")
plots2$overall  # View overall proportions
plots2$stratified  # View stratified proportions



#Sampling Simulation
# Run simulation
sim_results <- compare_sampling_methods(
  data = target_poplation_df,
  tarv_1 = "your_continuous_var",
  tarv_2 = "your_binary_var",
  thres = your_threshold,
  n = your_sample_size,
  strata_var = "your_strata_var",
  n_simulations = 1000
)

# Plot results for continuous variable
plots_tarv1 <- plot_simulation_results(sim_results, "tarv_1")
plots_tarv1$density_plot
plots_tarv1$coverage_plot

# Plot results for proportion variable
plots_tarv2 <- plot_simulation_results(sim_results, "tarv_2")
plots_tarv2$density_plot
plots_tarv2$coverage_plot
