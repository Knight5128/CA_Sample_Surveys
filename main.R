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
# 3. Apply these estimation methods on 2 different target variables under 2 sampling methods:
#       Target Variable 1: mean INCEARN
#       Methods used:
#         Simple Random Sampling:
#            vanilla estimation
#            ratio estimation
#            (take UHRSWORK as auxiliary variable)
#       Stratified Sampling with proportional allocation: 
#          use 4 strata differention variables separately:
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
#       Target Variable 2: proportion of people whose INCEARN is above income.thres
#       Methods used:
#         Simple Random Sampling:
#            vanilla estimation
#         Stratified Sampling with proportional allocation: 
#          use 4 strata differention variables separately:
#            - STATE (colname: STATEICP)
#            - SEX (colname: SEX)
#            - AGE GROUP (colname: AGE)
#              (need to add a new column to the df to categorize age into groups)
#              (Age groups: 0-18, 19-30, 31-50, 51-70, 71+)
#            - RACE (colname:RACE)
# 4. Each time step (3) is taken, there should be 11 different estimates and their corresponding SEs and CIs(along with other results) in total, 
#    and they should be stored for further comparison
# 4. Simulate step 2 for <1000> times
# 5. Compare the average value of the following results for each of the 11 estimation scheme:
#       estimate
#       SE
#       CI 
#    with the corresponding true population parameter
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
#———————————————————————Sampling Simulation Preperation————————————————————————————#
####################################################################################


compare_sampling_methods <- function(data, tarv_1, tarv_2, thres, n, strata_vars, n_simulations = 1000) {
  # 加载progress包
  if (!require("progress")) {
    install.packages("progress")
    library(progress)
  }
  
  # 计算总体参数真值
  true_mean <- mean(data[[tarv_1]])
  true_prop <- mean(data[[tarv_2]] > thres)
  
  # 为每个分层变量初始化结果存储
  results <- list()
  for(strata_var in strata_vars) {
    results[[strata_var]] <- list(
      # 目标变量1: 均值估计结果
      mean_estimates = list(
        srs_vanilla = data.frame(
          estimate = numeric(n_simulations),
          se = numeric(n_simulations),
          ci_lower = numeric(n_simulations),
          ci_upper = numeric(n_simulations)
        ),
        srs_ratio = data.frame(
          estimate = numeric(n_simulations),
          se = numeric(n_simulations),
          ci_lower = numeric(n_simulations),
          ci_upper = numeric(n_simulations)
        ),
        str_prop = data.frame(
          estimate = numeric(n_simulations),
          se = numeric(n_simulations),
          ci_lower = numeric(n_simulations),
          ci_upper = numeric(n_simulations)
        )
      ),
      # 目标变量2: 比例估计结果
      prop_estimates = list(
        srs_vanilla = data.frame(
          estimate = numeric(n_simulations),
          se = numeric(n_simulations),
          ci_lower = numeric(n_simulations),
          ci_upper = numeric(n_simulations)
        ),
        str_prop = data.frame(
          estimate = numeric(n_simulations),
          se = numeric(n_simulations),
          ci_lower = numeric(n_simulations),
          ci_upper = numeric(n_simulations)
        )
      )
    )
  }
  
  # 创建进度条
  pb <- progress_bar$new(
    format = paste0("Year ", year, " - Simulation progress [:bar] :percent eta: :eta"),
    total = n_simulations,
    clear = FALSE,
    width = 80
  )
  
  # 运行模拟
  for(i in 1:n_simulations) {
    # SRS方法只需要运行一次（与分层变量无关）
    srs_vanilla_mean <- srs_sampling_est_mean_vanilla(data, tarv_1, n)
    srs_ratio_mean <- srs_sampling_est_mean_ratio(data, tarv_1, "UHRSWORK", n)
    srs_vanilla_prop <- srs_sampling_est_prop_vanilla(data, tarv_2, thres, n)
    
    # 对每个分层变量进行分层抽样
    for(strata_var in strata_vars) {
      # 分层抽样估计
      str_prop_mean <- str_prop_sampling_est_mean_vanilla(data, tarv_1, n, strata_var)
      str_prop_prop <- str_prop_sampling_est_prop_vanilla(data, tarv_2, thres, n, strata_var)
      
      # 存储均值估计结果
      results[[strata_var]]$mean_estimates$srs_vanilla[i,] <- c(
        srs_vanilla_mean$estimate,
        srs_vanilla_mean$se,
        srs_vanilla_mean$ci
      )
      
      results[[strata_var]]$mean_estimates$srs_ratio[i,] <- c(
        srs_ratio_mean$estimate,
        srs_ratio_mean$se,
        srs_ratio_mean$ci
      )
      
      results[[strata_var]]$mean_estimates$str_prop[i,] <- c(
        str_prop_mean$estimate,
        str_prop_mean$se,
        str_prop_mean$ci
      )
      
      # 存储比例估计结果
      results[[strata_var]]$prop_estimates$srs_vanilla[i,] <- c(
        srs_vanilla_prop$estimate,
        srs_vanilla_prop$se,
        srs_vanilla_prop$ci
      )
      
      results[[strata_var]]$prop_estimates$str_prop[i,] <- c(
        str_prop_prop$estimate,
        str_prop_prop$se,
        str_prop_prop$ci
      )
    }
    
    # 更新进度条
    pb$tick()
  }
  
  # 计算每个分层变量的性能指标
  performance <- list()
  for(strata_var in strata_vars) {
    performance[[strata_var]] <- list(
      mean_estimates = list(),
      prop_estimates = list()
    )
    
    # 均值估计的性能指标
    for(method in c("srs_vanilla", "srs_ratio", "str_prop")) {
      performance[[strata_var]]$mean_estimates[[method]] <- list(
        estimates = results[[strata_var]]$mean_estimates[[method]]$estimate,
        ci_lower = results[[strata_var]]$mean_estimates[[method]]$ci_lower,
        ci_upper = results[[strata_var]]$mean_estimates[[method]]$ci_upper,
        bias = mean(results[[strata_var]]$mean_estimates[[method]]$estimate - true_mean),
        rmse = sqrt(mean((results[[strata_var]]$mean_estimates[[method]]$estimate - true_mean)^2)),
        coverage = mean(results[[strata_var]]$mean_estimates[[method]]$ci_lower <= true_mean & 
                       results[[strata_var]]$mean_estimates[[method]]$ci_upper >= true_mean),
        avg_ci_width = mean(results[[strata_var]]$mean_estimates[[method]]$ci_upper - 
                           results[[strata_var]]$mean_estimates[[method]]$ci_lower)
      )
    }
    
    # 比例估计的性能指标
    for(method in c("srs_vanilla", "str_prop")) {
      performance[[strata_var]]$prop_estimates[[method]] <- list(
        estimates = results[[strata_var]]$prop_estimates[[method]]$estimate,
        ci_lower = results[[strata_var]]$prop_estimates[[method]]$ci_lower,
        ci_upper = results[[strata_var]]$prop_estimates[[method]]$ci_upper,
        bias = mean(results[[strata_var]]$prop_estimates[[method]]$estimate - true_prop),
        rmse = sqrt(mean((results[[strata_var]]$prop_estimates[[method]]$estimate - true_prop)^2)),
        coverage = mean(results[[strata_var]]$prop_estimates[[method]]$ci_lower <= true_prop & 
                       results[[strata_var]]$prop_estimates[[method]]$ci_upper >= true_prop),
        avg_ci_width = mean(results[[strata_var]]$prop_estimates[[method]]$ci_upper - 
                           results[[strata_var]]$prop_estimates[[method]]$ci_lower)
      )
    }
  }
  
  # 返回结果
  return(list(
    true_values = list(
      mean = true_mean,
      proportion = true_prop
    ),
    simulation_results = results,
    performance = performance
  ))
}


# 构建分层变量性能比较表格
compare_strata_performance <- function(performance_results) {
  # 初始化结果数据框
  mean_comparisons <- data.frame()
  prop_comparisons <- data.frame()
  
  # 遍历每个分层变量
  for(strata_var in names(performance_results)) {
    # 均值估计比较
    mean_perf <- data.frame(
      strata_var = strata_var,
      method = c("SRS-Vanilla", "SRS-Ratio", "Stratified"),
      # 使用均值而不是原始向量
      estimate = c(
        mean(performance_results[[strata_var]]$mean_estimates$srs_vanilla$estimates),
        mean(performance_results[[strata_var]]$mean_estimates$srs_ratio$estimates),
        mean(performance_results[[strata_var]]$mean_estimates$str_prop$estimates)
      ),
      ci_lower = c(
        mean(performance_results[[strata_var]]$mean_estimates$srs_vanilla$ci_lower),
        mean(performance_results[[strata_var]]$mean_estimates$srs_ratio$ci_lower),
        mean(performance_results[[strata_var]]$mean_estimates$str_prop$ci_lower)
      ),
      ci_upper = c(
        mean(performance_results[[strata_var]]$mean_estimates$srs_vanilla$ci_upper),
        mean(performance_results[[strata_var]]$mean_estimates$srs_ratio$ci_upper),
        mean(performance_results[[strata_var]]$mean_estimates$str_prop$ci_upper)
      ),
      bias = c(
        performance_results[[strata_var]]$mean_estimates$srs_vanilla$bias,
        performance_results[[strata_var]]$mean_estimates$srs_ratio$bias,
        performance_results[[strata_var]]$mean_estimates$str_prop$bias
      ),
      rmse = c(
        performance_results[[strata_var]]$mean_estimates$srs_vanilla$rmse,
        performance_results[[strata_var]]$mean_estimates$srs_ratio$rmse,
        performance_results[[strata_var]]$mean_estimates$str_prop$rmse
      ),
      coverage = c(
        performance_results[[strata_var]]$mean_estimates$srs_vanilla$coverage,
        performance_results[[strata_var]]$mean_estimates$srs_ratio$coverage,
        performance_results[[strata_var]]$mean_estimates$str_prop$coverage
      ),
      ci_width = c(
        performance_results[[strata_var]]$mean_estimates$srs_vanilla$avg_ci_width,
        performance_results[[strata_var]]$mean_estimates$srs_ratio$avg_ci_width,
        performance_results[[strata_var]]$mean_estimates$str_prop$avg_ci_width
      )
    )
    mean_comparisons <- rbind(mean_comparisons, mean_perf)
    
    # 比例估计比较
    prop_perf <- data.frame(
      strata_var = strata_var,
      method = c("SRS-Vanilla", "Stratified"),
      # 使用均值而不是原始向量
      estimate = c(
        mean(performance_results[[strata_var]]$prop_estimates$srs_vanilla$estimates),
        mean(performance_results[[strata_var]]$prop_estimates$str_prop$estimates)
      ),
      ci_lower = c(
        mean(performance_results[[strata_var]]$prop_estimates$srs_vanilla$ci_lower),
        mean(performance_results[[strata_var]]$prop_estimates$str_prop$ci_lower)
      ),
      ci_upper = c(
        mean(performance_results[[strata_var]]$prop_estimates$srs_vanilla$ci_upper),
        mean(performance_results[[strata_var]]$prop_estimates$str_prop$ci_upper)
      ),
      bias = c(
        performance_results[[strata_var]]$prop_estimates$srs_vanilla$bias,
        performance_results[[strata_var]]$prop_estimates$str_prop$bias
      ),
      rmse = c(
        performance_results[[strata_var]]$prop_estimates$srs_vanilla$rmse,
        performance_results[[strata_var]]$prop_estimates$str_prop$rmse
      ),
      coverage = c(
        performance_results[[strata_var]]$prop_estimates$srs_vanilla$coverage,
        performance_results[[strata_var]]$prop_estimates$str_prop$coverage
      ),
      ci_width = c(
        performance_results[[strata_var]]$prop_estimates$srs_vanilla$avg_ci_width,
        performance_results[[strata_var]]$prop_estimates$str_prop$avg_ci_width
      )
    )
    prop_comparisons <- rbind(prop_comparisons, prop_perf)
  }
  
  return(list(
    mean_estimates = mean_comparisons,
    proportion_estimates = prop_comparisons
  ))
}


# 首先创建一个函数来整合多年的结果
compare_yearly_performance <- function(yearly_results) {
  # 初始化结果数据框
  all_years_mean <- data.frame()
  all_years_prop <- data.frame()
  
  # 遍历每一年
  for(year in names(yearly_results)) {
    # 获取当年的性能比较结果
    year_comparisons <- compare_strata_performance(yearly_results[[year]]$performance)
    
    # 添加年份列
    year_comparisons$mean_estimates$year <- year
    year_comparisons$proportion_estimates$year <- year
    
    # 合并结果
    all_years_mean <- rbind(all_years_mean, year_comparisons$mean_estimates)
    all_years_prop <- rbind(all_years_prop, year_comparisons$proportion_estimates)
  }
  
  return(list(
    mean_estimates = all_years_mean,
    proportion_estimates = all_years_prop
  ))
}

# 创建多年比较的可视化函数
plot_yearly_performance <- function(yearly_results) {
  library(ggplot2)
  library(tidyr)
  library(dplyr)
  library(gridExtra)
  
  # 获取多年比较数据
  comparisons <- compare_yearly_performance(yearly_results)
  print(names(comparisons$mean_estimates))
  # 处理数据为绘图格式
  mean_plot_data <- comparisons$mean_estimates %>%
    select(year, strata_var, method, estimate, ci_lower, ci_upper) %>%
    mutate(
      year = factor(year),
      estimate_type = "Mean Income"
    )
  
  prop_plot_data <- comparisons$proportion_estimates %>%
    select(year, strata_var, method, estimate, ci_lower, ci_upper) %>%
    mutate(
      year = factor(year),
      estimate_type = "Proportion Above Threshold"
    )
  
  # 合并数据
  plot_data <- rbind(mean_plot_data, prop_plot_data)
  
  # 为每个estimate_type获取真实值
  true_values <- sapply(names(yearly_results), function(year) {
    c(yearly_results[[year]]$true_values$mean,
      yearly_results[[year]]$true_values$proportion)
  })
  
  true_value_df <- data.frame(
    year = rep(names(yearly_results), each = 2),
    estimate_type = rep(c("Mean Income", "Proportion Above Threshold"), length(yearly_results)),
    true_value = c(true_values)
  )
  
  # 创建可视化
  ggplot(plot_data, aes(x = method, y = estimate, color = strata_var)) +
    # 添加置信区间
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), 
                 width = 0.2, 
                 position = position_dodge(width = 0.8)) +
    # 添加点估计
    geom_point(position = position_dodge(width = 0.8), size = 2) +
    # 添加真实值参考线
    geom_hline(data = true_value_df, 
               aes(yintercept = true_value),
               linetype = "dashed",
               color = "red") +
    # 分面
    facet_grid(estimate_type ~ year, scales = "free_y") +
    # 主题设置
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 8),
      strip.text = element_text(size = 12)
    ) +
    # 标签
    labs(
      x = "Estimation Method",
      y = "Estimate with 95% CI",
      color = "Stratification Variable",
      title = "Comparison of Estimation Methods Across Years",
      subtitle = "Red dashed lines indicate true population values"
    )
}


####################################################################################
####################################################################################
#—————————————————————————————Sampling Simulation——————————————————————————————————#
####################################################################################
####################################################################################


# Run simulations for each stratification variable and year
strata_vars <- c("STATEICP", "SEX", "AGE_GROUP", "RACE")
years <- c("2016", "2020", "2022")
n_sims <- 1000

# Initialize results storage
big_results <- list()


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


##########################################以下需要修改
# 2. 创建比较可视化
yearly_plots <- plot_yearly_performance(big_results)
yearly_plots
# 3. 保存结果
# 保存图表
ggsave("mean_estimation_yearly_comparison.png", 
       yearly_plots$mean_plots, 
       width = 15, 
       height = 20)
ggsave("proportion_estimation_yearly_comparison.png", 
       yearly_plots$proportion_plots, 
       width = 15, 
       height = 20)

# 如果需要导出数据表格
yearly_comparisons <- compare_yearly_performance(big_results)
write.csv(yearly_comparisons$mean_estimates, 
          "mean_estimation_yearly_comparison.csv")
write.csv(yearly_comparisons$proportion_estimates, 
          "proportion_estimation_yearly_comparison.csv")

####################################################################################
#————————————————————————————Results Visualization————————————————————————————————#
####################################################################################

# 可视化比较不同分层指标的表现
plot_strata_performance <- function(performance_results) {
  library(ggplot2)
  library(tidyr)
  library(dplyr)
  
  # 获取比较表格
  comparisons <- compare_strata_performance(performance_results)
  
  # 均值估计的可视化
  mean_plot_data <- comparisons$mean_estimates %>%
    gather(metric, value, bias:ci_width)
  
  mean_plots <- ggplot(mean_plot_data, aes(x = strata_var, y = value, fill = method)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap(~metric, scales = "free_y", ncol = 2) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Performance Comparison for Mean Estimation",
         x = "Stratification Variable",
         y = "Value") +
    scale_fill_brewer(palette = "Set2")
  
  # 比例估计的可视化
  prop_plot_data <- comparisons$proportion_estimates %>%
    gather(metric, value, bias:ci_width)
  
  prop_plots <- ggplot(prop_plot_data, aes(x = strata_var, y = value, fill = method)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap(~metric, scales = "free_y", ncol = 2) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Performance Comparison for Proportion Estimation",
         x = "Stratification Variable",
         y = "Value") +
    scale_fill_brewer(palette = "Set2")
  
  return(list(
    mean_plots = mean_plots,
    proportion_plots = prop_plots
  ))
}
