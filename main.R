# nolint: object_name_linter
#installation of necessary packages
source("utils.R")
library(tidyverse)
library(ggplot2)
library(ipumsr)

####################################################################################
#——————————————————————————————Research Description————————————————————————————————#
####################################################################################



# We have 6 kinds of income types, and we should choose the one(or some) that can 
# most reflect the well-being of the average citizen in the US.
# The 6 income types are:
#   1. INCTOT: total personal income
#         INCTOT reports each respondent's total pre-tax personal income or losses from all sources for the previous year
#   2. INCWAGE: wage and salary income
#         INCWAGE reports each respondent's total pre-tax wage and salary income, that is, money received as an employee for work or services performed
#   3. INCWELFR: welfare(public assistance) income
#         INCWELFR reports how much pre-tax income (if any) the respondent received during the previous year from various public assistance programs
#   4. INCRETIR: retirement income
#         INCRETIR reports how much pre-tax retirement, survivor, and disability pension income, other than Social Security, the respondent received during the previous year
#   5. INCSUPP: supplementary security income
#         INCSUPP reports how much pre-tax income (if any) the respondent received from Supplemental Security Income (SSI) during the previous year
#   6. INCEARN: total personal earned income
#         INCEARN reports income earned from wages or a person's own business or farm for the previous year


# 在相关性分析部分，发现INCEARN与每周工作小时数的相关性比较高，适合应用ratio estimation
# 故初步考虑最关心INCEARN这一变量，THRESHOLD也要围绕他进行查找&说明！！！！！！


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
rm(rawdf_222016)

#check if there's any NA in each df
colSums(is.na(rawdf_2016))
colSums(is.na(rawdf_2020))
colSums(is.na(rawdf_2022))
#No NA in any df!!!


#Select only relevant columns that are relevant to our study:
#i.e. keep what we care about
select_relevant_columns <- function(df) {
  df %>%
    select(-c(YEAR, SAMPLE, CLUSTER, STRATA, GQ, WKSWORK1, RACED))
}
df_2016 <- select_relevant_columns(rawdf_2016)
df_2020 <- select_relevant_columns(rawdf_2020)
df_2022 <- select_relevant_columns(rawdf_2022)


                                    
#find plausible auxiliary variable for Ratio Estimation 
explore_correlations <- function(df, title = "Finding Auxiliary Variable", save_plot = FALSE) {
  library(corrplot)
  
  # Get numeric columns
  numeric_cols <- sapply(df, is.numeric)
  numeric_data <- df[, numeric_cols]
  
  # Calculate correlation matrix
  cor_matrix <- cor(numeric_data, use = "complete.obs")
  
  # Create correlation plot
  if (save_plot) {
    # Save to file
    png(paste0(title, "_correlation_plot.png"), 
        width = 1200, height = 1200, res = 100)
    corrplot(cor_matrix, 
             method = "color", 
             type = "upper", 
             tl.col = "black",
             tl.srt = 45,
             addCoef.col = "black",
             number.cex = 0.7,
             title = paste("Correlation Plot -", title))
    dev.off()
  } else {
    # Display in R interface
    corrplot(cor_matrix, 
             method = "color", 
             type = "upper", 
             tl.col = "black",
             tl.srt = 45,
             addCoef.col = "black",
             number.cex = 0.7,
             title = paste("Correlation Plot -", title))
  }
  
  # Return correlation matrix
  return(cor_matrix)
}

# To see whether there is a single variable that is highly correlated with a person's earned income
explore_correlations(df_2016[,-c(1,2,3,4,5,6,8,12:17)])
explore_correlations(df_2020[,-c(1,2,3,4,5,6,8,12:17)])
explore_correlations(df_2022[,-c(1,2,3,4,5,6,8,12:17)])
# THen can find that "UHRSWORK" is highly correlated with "INCEARN" in all 3 years
# So we can use "UHRSWORK" as the auxiliary variable for Ratio Estimation

                                    
# Now, what we're going to do is:
# 1. Find the threshold for INCEARN, say its value is income.thres, a pre-set constant
#    We have 2 variables of interst:
#       1st: average personal earned income in the US, that is, average INCEARN (colname: INCEARN_CPIU_2020)
#       2nd: the proportion of people whose INCEARN is above income.thres 
#    Then we want to sample a certain proportion of the population to estimate these 2 variables
#    Suppose our total sample size is n, also predetermined
# 2. Find the true value of these 2 variables using our full dataset for results evaluation
# 3. Apply 6 estimation methods under 2 sampling methods:
#       Simple Random Sampling:
#            vanilla estimation
#            ratio estimation
#               (take UHRSWORK as auxiliary variable)
#       Stratified Sampling with proportional allocation: 
#          strata differentiated by:
#            - STATE (colname: STATEICP)
#            - SEX (colname: SEX)
#            - AGE GROUP (colname: AGE)
#              (need to add a new column to the df to categorize age into groups)
#              (Age groups: 0-18, 19-30, 31-50, 51-70, 71+)
#            - RACE (colname:RACE)
#              (integer 1 to 9 each representing a different race:
#               1: White, 
#               2: Black/African American, 
#               3: American Indian or Alaska Native, 
#               4: Chinese, 
#               5: Japanese, 
#               6: Other Asian or Pacific Islander, 
#               7: Other race, nec, 
#               8: Two major races, 
#               9: Three or more major races
#               )
# 4. Simulate step 2 for <1000> times
# 5. Compare:
#       estimate
#       SE
#       CI 
#    of the 6 estimation methods under 2 sampling methods with
#    corresponding true population parameter
# 6. Plot the results in a vertical-layout 6 row, 1 col grid
#    with each graph in a row representing a different estimation method




####################################################################################
#—————————————————————————————Analysis Parameters——————————————————————————————————#
####################################################################################

# Set random seed for reproducibility
set.seed(2024)

# Pre-determined parameters
income.thres <- 50000  # threshold for INCEARN
sample_size <- 10000   # total sample size for each simulation

####################################################################################
#————————————————————————————Data Preprocessing———————————————————————————————————#
####################################################################################

# Add age groups to the datasets
add_age_groups <- function(df) {
  df %>%
    mutate(AGE_GROUP = case_when(
      AGE <= 18 ~ "0-18",
      AGE <= 30 ~ "19-30",
      AGE <= 50 ~ "31-50",
      AGE <= 70 ~ "51-70",
      TRUE ~ "71+"
    ))
}

df_2016 <- add_age_groups(df_2016)
df_2020 <- add_age_groups(df_2020)
df_2022 <- add_age_groups(df_2022)

####################################################################################
#—————————————————————————True Population Parameters————————————————————————————————#
####################################################################################

# Calculate true parameters for each year
get_true_params <- function(df) {
  list(
    mean_income = mean(df$INCEARN_CPIU_2010),
    prop_above_thres = mean(df$INCEARN_CPIU_2010 > income.thres)
  )
}

true_params <- list(
  "2016" = get_true_params(df_2016),
  "2020" = get_true_params(df_2020),
  "2022" = get_true_params(df_2022)
)

####################################################################################
#————————————————————————————Sampling Simulation——————————————————————————————————#
####################################################################################

# Run simulations for each stratification variable and year
strata_vars <- c("STATEICP", "SEX", "AGE_GROUP", "RACE")
years <- c("2016", "2020", "2022")
n_sims <- 1000

# Initialize results storage
results <- list()

for(year in years) {
  df <- get(paste0("df_", year))
  results[[year]] <- list()
  
  for(strata_var in strata_vars) {
    # Run simulation
    sim_results <- compare_sampling_methods(
      data = df,
      tarv_1 = "INCEARN_CPIU_2010",
      tarv_2 = "INCEARN_CPIU_2010",
      thres = income.thres,
      n = sample_size,
      strata_var = strata_var,
      n_simulations = n_sims
    )
    
    results[[year]][[strata_var]] <- sim_results
  }
}

####################################################################################
#————————————————————————————Results Visualization————————————————————————————————#
####################################################################################

# Create visualization for each year
plot_year_results <- function(year_results, year) {
  # Combine results from all stratification methods
  plots <- list()
  
  for(strata_var in names(year_results)) {
    # Plot continuous variable results
    p1 <- plot_simulation_results(year_results[[strata_var]], "tarv_1")
    plots[[paste0(strata_var, "_continuous")]] <- p1
    
    # Plot proportion results
    p2 <- plot_simulation_results(year_results[[strata_var]], "tarv_2")
    plots[[paste0(strata_var, "_proportion")]] <- p2
  }
  
  return(plots)
}

# Generate plots for each year
plots <- lapply(names(results), function(year) {
  plot_year_results(results[[year]], year)
})

# Save results for report
saveRDS(list(
  parameters = list(
    threshold = income.thres,
    sample_size = sample_size,
    n_simulations = n_sims
  ),
  true_parameters = true_params,
  simulation_results = results,
  plots = plots
), "sampling_analysis_results.rds")
