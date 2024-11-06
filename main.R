# nolint: object_name_linter
#installation of necessary packages
source("utils.R")
library(tidyverse)
library(ggplot2)
library(ipumsr)

####################################################################################
#——————————————————————————————Research Description————————————————————————————————#
####################################################################################





####################################################################################
#———————————————————————————————————Data Import————————————————————————————————————#
####################################################################################

#----------DO NOT try to import the .csv file directly. It will not work.----------#
#----------You need to import the data using the ipumsr package----------#
#----------Also, DO NOT change the file name of <*.dat>----------#

ddi <- read_ipums_ddi("./data/usa_222016_01.xml")
rawdf_222016 <- read_ipums_micro(ddi)

#split data into 3 years: 2016, 2020, 2022

rawdf_2016 <- rawdf_222016[rawdf_222016$YEAR == 2016,]
rawdf_222016 <- rawdf_222016[rawdf_222016$YEAR != 2016,]
rawdf_2020 <- rawdf_222016[rawdf_222016$YEAR == 2020,]
rawdf_222016 <- rawdf_222016[rawdf_222016$YEAR != 2020,]
rawdf_2022 <- rawdf_222016[rawdf_222016$YEAR == 2022,]
rawdf_222016 <- rawdf_222016[rawdf_222016$YEAR != 2022,]

#check if there's any NA in each df
colSums(is.na(rawdf_2016))
colSums(is.na(rawdf_2020))
colSums(is.na(rawdf_2022))
#No NA in any df!!!


#Select only relevant columns that are relevant to our study:
#i.e. keep what we care about
rawdf_2016 <- rawdf_2016 %>% select(YEAR, SERIAL, CBSERIAL, HHWT, STATEFIP, PERNUM）

                                    
#find plausible auxiliary variable for Ratio Estimation 




                                    

####################################################################################
#————————————————————————————————Data Visualization————————————————————————————————#
####################################################################################
# For continuous variable (tarv_1)
plots1 <- plot_continuous_distribution(data, "your_var_name", "strata_column")
plots1$density  # View density plot
plots1$boxplot  # View boxplot
plots1$stratified  # View stratified boxplot

# For proportion variable (tarv_2)
plots2 <- plot_proportion_distribution(data, "your_var_name", threshold = 100, "strata_column")
plots2$overall  # View overall proportions
plots2$stratified  # View stratified proportions



####################################################################################
#——————————————————————————————Sampling Simulation—————————————————————————————————#
####################################################################################
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
