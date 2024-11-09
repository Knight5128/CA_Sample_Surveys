# nolint: object_name_linter
#installation of necessary packages
source("utils.R")
library(tidyverse)
library(ggplot2)
library(ipumsr)


####################################################################################
#—————————————————————————————Notes before submitting——————————————————————————————#
####################################################################################
# 将 所有中文注释 & 非指导、引导性的文字 翻译为英文
# 删除这个模块！！！！！！！！！！！！！


################ 未完待续:
# 确认income.thres以及n后，进行一次大型模拟，保存好big_result！！！！！！！
# 画出第一张图
# 修改第二张表格图中的Method名
# 看如何将这两张图直接加入到报告中






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
#——————————————————Data Import & Data Processing - 1st Stage———————————————————————#
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


####################################################################################
#———————————————Find plausible auxiliary variable for Ratio Estimation—————————————#
####################################################################################


explore_correlations <- function(df, title = "Finding Auxiliary Variable", save_plot = FALSE) {
  library(corrplot)
  
  # 如果需要保存图片且results目录不存在，则创建该目录
  if (save_plot && !dir.exists("./results")) {
    dir.create("./results")
  }
  
  # 获取数值型列
  numeric_cols <- sapply(df, is.numeric)
  numeric_data <- df[, numeric_cols]
  
  # 计算相关矩阵
  cor_matrix <- cor(numeric_data, use = "complete.obs")
  
  # 创建相关图
  if (save_plot) {
    # 保存到文件
    png(paste0("./results/", title, "_correlation_plot.png"), 
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
    # 在R界面显示
    corrplot(cor_matrix, 
             method = "color", 
             type = "upper", 
             tl.col = "black",
             tl.srt = 45,
             addCoef.col = "black",
             number.cex = 0.7,
             title = paste("Correlation Plot -", title))
  }
  
  # 返回相关矩阵
  return(cor_matrix)
}

# To see whether there is a single variable that is highly correlated with a person's earned income
explore_correlations(df_2016[,-c(1,2,3,4,5,6,8,12:17,20,23,25,27,29)], "2016 Data", save_plot = TRUE)
explore_correlations(df_2020[,-c(1,2,3,4,5,6,8,12:17,20,23,25,27,29)], "2020 Data", save_plot = TRUE)
explore_correlations(df_2022[,-c(1,2,3,4,5,6,8,12:17,20,23,25,27,29)], "2022 Data", save_plot = TRUE)
# THen can find that "UHRSWORK" is highly correlated with "INCEARN" in all 3 years
# So we can use "UHRSWORK" as the auxiliary variable for Ratio Estimation

                           
         
####################################################################################
#————————————————————————————————Estimation Planning———————————————————————————————#
####################################################################################

# Next, what we're going to do is:
# 1. Find the threshold for INCEARN, say its value is income.thres, a pre-set constant
#    Determine a uniform total sample size in consideration of real-problem cost constraints & other factors possible
# 2. We have 2 variables of interst:
#       1st: average personal earned income in the US, that is, average INCEARN (colname: INCEARN_CPIU_2020)
#       2nd: the proportion of people whose INCEARN is above income.thres 
#    
#    Then we want to sample a certain proportion of the population to estimate these 2 variables
#    Find the true value of these 2 variables using our full dataset for results evaluation
# 3. Apply these estimation methods on 2 different target variables under 11 estimation methods:
#       Target Variable 1: mean INCEARN
#       Methods used:
#         Simple Random Sampling:
#            vanilla estimation
#            ratio estimation
#            (take UHRSWORK as auxiliary variable)
#         Stratified Sampling with proportional allocation: 
#           use 4 strata differention variables separately:
#             - STATE 
#               (colname: STATEICP)
#             - SEX 
#               (colname: SEX)
#             - AGE GROUP 
#               (colname: AGE)
#               (need to add a new column to the df to categorize age into groups)
#               (Age groups: 0-18, 19-30, 31-50, 51-70, 71+)
#             - RACE 
#               (colname:RACE)
#               (integer 1 to 9 each representing a different race:
#                1: White, 
#                2: Black/African American, 
#                3: American Indian or Alaska Native, 
#                4: Chinese, 
#                5: Japanese, 
#                6: Other Asian or Pacific Islander, 
#                7: Other race, nec, 
#                8: Two major races, 
#                9: Three or more major races
#                )
#       Target Variable 2: proportion of people whose INCEARN is above income.thres
#       Methods used:
#         Simple Random Sampling:
#            vanilla estimation
#         Stratified Sampling with proportional allocation: 
#          use 4 strata differention variables separately:
#            - STATE 
#              (colname: STATEICP)
#            - SEX 
#              (colname: SEX)
#            - AGE GROUP 
#              (colname: AGE)
#              (need to add a new column to the df to categorize age into groups)
#              (Age groups: 0-18, 19-30, 31-50, 51-70, 71+)
#            - RACE 
#              (colname:RACE)
# 4. Each time step (3) is taken, there should be 11 different estimates and their corresponding SEs and CIs(along with other results) in total, 
#    and they should be stored for further comparison
# 4. Simulate step 2 for <> times
# 5. Compare the average value of the following results for each of the 11 estimation scheme:
#       estimate
#       SE
#       CI 
#    with the corresponding true population parameter
# 6. Plot the results in a vertical-layout 6 row, 1 col grid
#    with each graph in a row representing a different estimation method



####################################################################################
#———————————————————————————Determining Sample Size————————————————————————————————#
####################################################################################














####################################################################################
#—————————————————————————————Data Processing - 2nd Stage——————————————————————————#
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
####################################################################################
#—————————————————————————————————Main Analysis————————————————————————————————————#
####################################################################################
####################################################################################

# All the functions necessary for the main analysis are defined in utils.R
# Referring to utils.R for more details

####################################################################################
#—————————————————————————————Define Global Params—————————————————————————————————#
####################################################################################
# 1. Set random seed for reproducibility
set.seed(2024)

# 2. Pre-determined parameters
income.thres <- 50000  # threshold for INCEARN
sample_size <- 10000   # total sample size for each simulation

# 3. Set stratification variables and which years' data to analyse 
strata_vars <- c("STATEICP", "SEX", "AGE_GROUP", "RACE")
years <- c("2016", "2020", "2022")

# 4. number of simulations
n_sims <- 1000


####################################################################################
#—————————————————————————————Sampling Simulation——————————————————————————————————#
####################################################################################


# Initialize results storage
big_results <- list()

# Run simulations for each year
for(year in years) {
  df <- get(paste0("df_", year))
  big_results[[year]] <- compare_sampling_methods(
    data = df,
    tarv_1 = "INCEARN_CPIU_2010",
    tarv_2 = "INCEARN_CPIU_2010",
    thres = income.thres,
    n = sample_size,
    strata_vars = c("STATEICP", "SEX", "AGE_GROUP", "RACE"),
    n_simulations = 3
  )
}


####################################################################################
#—————————————————————————————Result Visualization—————————————————————————————————#
####################################################################################


# Visualize estimation results across years
yearly_plots <- plot_yearly_performance(big_results)

#————————————————————————————————————————————————————————————————————————#
#————————————————use this to export overall result plot——————————————————#
#————————————————————————————————————————————————————————————————————————#
yearly_plots


# Save plots (unnecessary)
ggsave("mean_estimation_yearly_comparison.png", 
       yearly_plots$mean_plots, 
       width = 15, 
       height = 20)
ggsave("proportion_estimation_yearly_comparison.png", 
       yearly_plots$proportion_plots, 
       width = 15, 
       height = 20)


# Display results as tables (unnecessary)
yearly_comparisons <- compare_yearly_performance(big_results)

# Save mean estimation results as HTML (unnecessary)
mean_table <- knitr::kable(yearly_comparisons$mean_estimates,
                          format = "html",
                          caption = "Mean Estimation Performance Across Years")
writeLines(mean_table, "./results/mean_estimation_table.html")

# Save proportion estimation results as HTML (unnecessary)
prop_table <- knitr::kable(yearly_comparisons$proportion_estimates, 
                          format = "html",
                          caption = "Proportion Estimation Performance Across Years")
writeLines(prop_table, "./results/proportion_estimation_table.html")

# Save results as CSV (unnecessary)
write.csv(yearly_comparisons$mean_estimates, 
          "./results/mean_estimation_results.csv", 
          row.names = FALSE)
write.csv(yearly_comparisons$proportion_estimates,
          "./results/proportion_estimation_results.csv",
          row.names = FALSE)


####################################################################################
#————————————————————————————Performance Comparison————————————————————————————————#
####################################################################################

# Here, we will look at the boarder picture of the performance metrics of 
# different estimation methods throughout the whole simulation process in all 3 years
# Including:
#   1. Average Bias
#      (Avg. difference between the estimate and the true population parameter)
#   2. Average SE
#      (Avg. standard errors of the estimates)
#   2. RMSE
#      (Rooted Mean Squared Error of estimate compared to the true population parameter)
#   3. Coverage Rate
#      (Proportion of CIs that contain the true population parameter)
#   5. **Relative Efficiency**
#      (Ratio of the variance of the RMSE to the RMSE of the SRS estimator)

perf_table = create_performance_table(big_results)
perf_table

# And we can see that... (see main report for further explanation)