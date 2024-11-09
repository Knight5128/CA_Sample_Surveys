library(dplyr)
library(stats)
library(ggplot2)


####################################################################################
#———————————————————————————————File Description———————————————————————————————————#
####################################################################################

# This is a file containing necessary code snippets & sampling functions & simulation functions for the project
# The functions are designed to be used in the main.R file under the same directory
# The file is structured as follows:
# 1. Snippet functions
# 2. Sampling functions
# 3. Simulation functions
# 4. Visualization functions
# 5. Other sampling functions

####################################################################################
#——————————————————————————————————————Basic———————————————————————————————————————#
####################################################################################


############################################################################# Function 1.1
# Calculate proportion of a variable above a threshold
calculate_proportion <- function(x, threshold) {
  mean(x > threshold)
}

############################################################################# Function 1.2
# Calculate variance of a binary variable
calculate_proportion_var <- function(x, threshold) {
  p <- mean(x > threshold)
  var <- p * (1 - p)
  return(var)
}

####################################################################################
#—————————————————————————————————————Sampling—————————————————————————————————————#
####################################################################################


############################################################################# Function 2.1
# Simple Random Sampling —— estimating mean —— vanilla estimatior
srs_sampling_est_mean_vanilla <- function(data, tarv, n) {
  # Input validation
  if (!tarv %in% names(data)) {
    stop("Target variable not found in dataset")
  }
  
  N <- nrow(data)
  if (n >= N) {
    stop("Sample size must be less than population size")
  }
  
  # Take simple random sample
  sample.ids <- sample(N, size = n, replace = FALSE)
  sample_data <- data[sample.ids, ]
  
  # Calculate estimate
  est <- mean(sample_data[[tarv]])
  
  # Calculate standard error
  var.smp <- var(sample_data[[tarv]])
  se <- sqrt((1 - n/N) * var.smp/n)
  
  # Create confidence interval
  CI <- c(est - 1.96 * se, 
          est + 1.96 * se)
  
  return(list(
    estimate = est,
    se = se,
    ci = CI
  ))
}


############################################################################# Function 2.2
# Simple Random Sampling —— estimating proportion —— vanilla estimator
srs_sampling_est_prop_vanilla <- function(data, tarv, thres, n) {
  # Input validation
  if (!tarv %in% names(data)) {
    stop("Target variable not found in dataset")
  }
  
  N <- nrow(data)
  if (n >= N) {
    stop("Sample size must be less than population size")
  }
  
  # Take simple random sample
  sample.ids <- sample(N, size = n, replace = FALSE)
  sample_data <- data[sample.ids, ]
  
  # Calculate proportion estimate
  est <- calculate_proportion(sample_data[[tarv]], thres)
  
  # Calculate standard error
  var.smp <- est * (1 - est)  # variance for binary variable
  se <- sqrt((1 - n/N) * var.smp/n)
  
  # Create confidence interval
  CI <- c(est - 1.96 * se,
          est + 1.96 * se)
  
  return(list(
    estimate = est,
    se = se,
    ci = CI
  ))
}


############################################################################# Function 2.3
# Simple Random Sampling —— estimating mean —— ratio estimator
srs_sampling_est_mean_ratio <- function(data, tarv_1, aux_var, n) {
  # Input validation
  if (!all(c(tarv_1, aux_var) %in% names(data))) {
    stop("Target or auxiliary variables not found in dataset")
  }
  
  N <- nrow(data)
  if (n >= N) {
    stop("Sample size must be less than population size")
  }
  
  # Take simple random sample
  sample.ids <- sample(N, size = n, replace = FALSE)
  sample_data <- data[sample.ids, ]
  
  # Calculate ratio estimate
  r_hat <- sum(sample_data[[tarv_1]]) / sum(sample_data[[aux_var]])
  X_bar <- mean(data[[aux_var]])  # Population mean of auxiliary variable
  tarv_1.est <- r_hat * X_bar
  
  # Calculate variance of ratio estimate
  y <- sample_data[[tarv_1]]
  x <- sample_data[[aux_var]]
  d <- y - r_hat * x
  var_d <- var(d)
  
  # Standard error
  tarv_1.se <- sqrt((1 - n/N) * var_d/n)
  
  # Create confidence interval
  tarv_1.CI <- c(tarv_1.est - 1.96 * tarv_1.se,
                 tarv_1.est + 1.96 * tarv_1.se)
  
  return(list(
    estimate = tarv_1.est,
    se = tarv_1.se,
    ci = tarv_1.CI,
    ratio = r_hat
  ))
}


############################################################################# Function 2.4
# Stratified Sampling with Proportional Allocation —— estimating mean —— vanilla estimator
str_prop_sampling_est_mean_vanilla <- function(data, tarv, n, strata_var) {
  # Input validation
  if (!all(c(tarv, strata_var) %in% names(data))) {
    stop("Required variables not found in dataset")
  }
  
  # Calculate stratum sizes
  N <- nrow(data)
  N_hs <- table(data[[strata_var]])
  strata_names <- names(N_hs)
  
  # Calculate proportional allocation
  n_hs <- round((N_hs/N) * n)
  
  # Initialize results storage
  tarv.bar.hs <- numeric(length(strata_names))
  tarv.var.smp.hs <- numeric(length(strata_names))
  
  # Sample from each stratum
  STR.sample <- data.frame()
  for (i in seq_along(strata_names)) {
    stratum_data <- if(strata_var == "AGE_GROUP") {
      data[data[[strata_var]] == strata_names[i], ]
    } else {
      data[data[[strata_var]] == as.integer(strata_names[i]), ]
    }
    
    # 从当前层中随机抽样
    if (nrow(stratum_data) > 0) {  # 确保该层有数据
      sample_indices <- sample.int(nrow(stratum_data), size = n_hs[i], replace = FALSE)
      stratum_sample <- stratum_data[sample_indices, ]
      STR.sample <- rbind(STR.sample, stratum_sample)
    }
  }

  # Calculate estimate for each stratum
  tarv.bar.hs <- tapply(STR.sample[[tarv]], STR.sample[[strata_var]], mean)
  
  # Calculate sample variance for each stratum
  tarv.var.smp.hs <- tapply(STR.sample[[tarv]], STR.sample[[strata_var]], var)
  
  # Calculate overall estimate
  tarv.str.est <- sum((N_hs/N) * tarv.bar.hs)
  
  # Calculate standard error
  tarv.se <- sqrt(sum((N_hs/N)^2 * (1 - n_hs/N_hs) * tarv.var.smp.hs/n_hs))
  
  # Create confidence interval
  tarv.CI <- c(tarv.str.est - 1.96 * tarv.se,
               tarv.str.est + 1.96 * tarv.se)
  
  return(list(
    estimate = tarv.str.est,
    se = tarv.se,
    ci = tarv.CI,
    strata_estimates = tarv.bar.hs
  ))
}


############################################################################# Function 2.5
# Stratified Sampling with Proportional Allocation —— estimating proportion —— vanilla estimator
str_prop_sampling_est_prop_vanilla <- function(data, tarv, thres, n, strata_var) {
  # Input validation
  if (!all(c(tarv, strata_var) %in% names(data))) {
    stop("Required variables not found in dataset")
  }
  
  # Calculate stratum sizes
  N <- nrow(data)
  N_hs <- table(data[[strata_var]])
  strata_names <- names(N_hs)
  
  # Calculate proportional allocation
  n_hs <- round((N_hs/N) * n)
  
  # Initialize results storage
  tarv.bar.hs <- numeric(length(strata_names))
  tarv.var.smp.hs <- numeric(length(strata_names))
  
  # Sample from each stratum
  STR.sample <- data.frame()
  for (i in seq_along(strata_names)) {
    stratum_data <- if(strata_var == "AGE_GROUP") {
      data[data[[strata_var]] == strata_names[i], ]
    } else {
      data[data[[strata_var]] == as.integer(strata_names[i]), ]
    }
    
    # 从当前层中随机抽样
    if (nrow(stratum_data) > 0) {  # 确保该层有数据
      sample_indices <- sample.int(nrow(stratum_data), size = n_hs[i], replace = FALSE)
      stratum_sample <- stratum_data[sample_indices, ]
      STR.sample <- rbind(STR.sample, stratum_sample)
    }
  }


  # Calculate proportion estimate for each stratum
  tarv.bar.hs <- tapply(STR.sample[[tarv]], STR.sample[[strata_var]], 
                        calculate_proportion, threshold = thres)
  
  # Calculate sample variance for each stratum
  tarv.var.smp.hs <- tapply(STR.sample[[tarv]], STR.sample[[strata_var]], 
                           calculate_proportion_var, threshold = thres)
  
  # Calculate overall estimate
  tarv.str.est <- sum((N_hs/N) * tarv.bar.hs)
  
  # Calculate standard error
  tarv.se <- sqrt(sum((N_hs/N)^2 * (1 - n_hs/N_hs) * tarv.var.smp.hs/n_hs))
  
  # Create confidence interval
  tarv.CI <- c(tarv.str.est - 1.96 * tarv.se,
               tarv.str.est + 1.96 * tarv.se)
  
  return(list(
    estimate = tarv.str.est,
    se = tarv.se,
    ci = tarv.CI,
    strata_estimates = tarv.bar.hs
  ))
}


####################################################################################
#————————————————————————————————————Simulation————————————————————————————————————#
####################################################################################


############################################################################# Function 3.1
# Simulating different sampling & estimation methods
# + evaluating performance across different sampling methods & strata variables used
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


############################################################################# Function 3.2
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


############################################################################# Function 3.3
# 首先创建一个函数来整合多年的结果
compare_yearly_performance <- function(yearly_results) {
  # 初始化结果数据框
  all_years_mean <- data.frame()
  all_years_prop <- data.frame()
  
  # 遍历每一年
  for(year in names(yearly_results)) {
    # 获取当年的性能比较结果
    year_comparisons <- compare_strata_performance(yearly_results[[year]]$performance)
    
    # 对于SRS-Vanilla和SRS-Ratio方法,只保留第一个分层变量的结果
    year_comparisons$mean_estimates <- year_comparisons$mean_estimates %>%
      filter((method != "SRS-Vanilla" & method != "SRS-Ratio") | strata_var == "STATEICP")
    
    year_comparisons$proportion_estimates <- year_comparisons$proportion_estimates %>%
      filter((method != "SRS-Vanilla" & method != "SRS-Ratio") | strata_var == "STATEICP")
    
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




####################################################################################
#————————————————————————————————————Visualization—————————————————————————————————#
####################################################################################


############################################################################# Function 4.1
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
      axis.text.x = element_text(size = 10),
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
      title = "Comparison of Different Estimation Methods Across 2016, 2020, 2022",
      subtitle = "Red dashed lines indicate true population values"
    )
}



####################################################################################
#————————————————————————————————Other Sampling————————————————————————————————————#
####################################################################################


############################################################################# Function 5.1
# Stratified Sampling with Equal Allocation
str_equal_sampling <- function(data, tarv_1, tarv_2, thres, n, strata_var) {
  # Input validation
  if (!all(c(tarv_1, tarv_2, strata_var) %in% names(data))) {
    stop("Required variables not found in dataset")
  }
  
  # Calculate stratum sizes
  N <- nrow(data)
  N_hs <- table(data[[strata_var]])
  strata_names <- names(N_hs)
  
  # divide n equally among strata
  H <- length(strata_names)  # number of strata
  n_hs <- rep(floor(n/H), H)  # equal allocation
  # Add remaining samples to maintain total n
  remainder <- n - sum(n_hs)
  if (remainder > 0) {
    n_hs[1:remainder] <- n_hs[1:remainder] + 1
  }
  
  # Initialize results storage
  tarv_1.bar.hs <- numeric(length(strata_names))
  tarv_2.bar.hs <- numeric(length(strata_names))
  tarv_1.var.smp.hs <- numeric(length(strata_names))
  tarv_2.var.smp.hs <- numeric(length(strata_names))
  
  # Sample from each stratum
  STR.sample <- data.frame()
  
  for (i in seq_along(strata_names)) {
    row.ids <- which(data[[strata_var]] == strata_names[i])
    sample.ids <- sample(row.ids, size = n_hs[i], replace = FALSE)
    STR.sample <- rbind(STR.sample, data[sample.ids,])
  }

  # Calculate estimate for each stratum
  tarv_1.bar.hs <- tapply(STR.sample[[tarv_1]], STR.sample[[strata_var]], mean)
  tarv_2.bar.hs <- tapply(STR.sample[[tarv_2]], STR.sample[[strata_var]], calculate_proportion, threshold = thres)

  # Calculate sample variance for each stratum
  tarv_1.var.smp.hs <- tapply(STR.sample[[tarv_1]], STR.sample[[strata_var]], var)
  tarv_2.var.smp.hs <- tapply(STR.sample[[tarv_2]], STR.sample[[strata_var]], calculate_proportion_var, threshold = thres)
  
  # Calculate overall estimates
  tarv_1.str.est <- sum((N_hs/N) * tarv_1.bar.hs)
  tarv_2.str.est <- sum((N_hs/N) * tarv_2.bar.hs)
  
  # Calculate standard errors
  tarv_1.se <- sqrt(sum((N_hs/N)^2 * (1 - n_hs/N_hs) * tarv_1.var.smp.hs/n_hs))
  tarv_2.se <- sqrt(sum((N_hs/N)^2 * (1 - n_hs/N_hs) * tarv_2.var.smp.hs/n_hs))
  
  # Create confidence intervals
  tarv_1.CI <- c(tarv_1.str.est - 1.96 * tarv_1.se,
               tarv_1.str.est + 1.96 * tarv_1.se)
  tarv_2.CI <- c(tarv_2.str.est - 1.96 * tarv_2.se,
               tarv_2.str.est + 1.96 * tarv_2.se)
  
  return(list(
    tarv_1 = list(
      estimate = tarv_1.str.est,
      se = tarv_1.se,
      ci = tarv_1.CI,
      strata_estimates = tarv_1.bar.hs
    ),
    tarv_2 = list(
      estimate = tarv_2.str.est,
      se = tarv_2.se,
      ci = tarv_2.CI,
      strata_estimates = tarv_2.bar.hs
    )
  ))
}


############################################################################# Function 5.2
# Stratified Sampling with Optimal Allocation
str_opt_sampling <- function(data, tarv_1, tarv_2, thres, n, strata_var, 
                           strata_s2_guess = NULL, strata_cost = NULL) {
  # Input validation
  if (!all(c(tarv_1, tarv_2, strata_var) %in% names(data))) {
    stop("Required variables not found in dataset")
  }
  
  # Calculate stratum sizes and standard deviations
  N <- nrow(data)
  strata_names <- unique(data[[strata_var]])
  N_h <- tapply(data[[tarv_1]], data[[strata_var]], length)
  


  # Use provided variance guesses or calculate from data
  if (is.null(strata_s2_guess)) {
    S_h <- tapply(data[[tarv_1]], data[[strata_var]], sd)
  } else {
    S_h <- sqrt(strata_s2_guess)
  }
  
  # Default cost is 1 for all strata if not provided
  if (is.null(strata_cost)) {
    strata_cost <- rep(1, length(N_h))
  }
  


  # Calculate optimal allocation considering N_h, S_h, and cost
  opt_allocation_factor <- N_h * S_h / sqrt(strata_cost)
  n_hs <- round(n * opt_allocation_factor / sum(opt_allocation_factor))
  
  # Initialize results storage
  tarv_1.bar.hs <- numeric(length(strata_names))
  tarv_2.bar.hs <- numeric(length(strata_names))
  tarv_1.var.smp.hs <- numeric(length(strata_names))
  tarv_2.var.smp.hs <- numeric(length(strata_names))
  
  # Sample from each stratum
  STR.sample <- data.frame()
  
  for (i in seq_along(strata_names)) {
    row.ids <- data[data[[strata_var]] == strata_names[i], ]
    sample.ids <- sample(nrow(row.ids), size = n_hs[i], replace = FALSE)
    stratum_sample <- row.ids[sample.ids, ]
    STR.sample <- rbind(STR.sample, stratum_sample)
  }
  

  # Calculate estimate for each stratum
  tarv_1.bar.hs <- tapply(STR.sample[[tarv_1]], STR.sample[[strata_var]], mean)
  tarv_2.bar.hs <- tapply(STR.sample[[tarv_2]], STR.sample[[strata_var]], calculate_proportion, threshold = thres)

  # Calculate sample variance for each stratum
  tarv_1.var.smp.hs <- tapply(STR.sample[[tarv_1]], STR.sample[[strata_var]], var)
  tarv_2.var.smp.hs <- tapply(STR.sample[[tarv_2]], STR.sample[[strata_var]], calculate_proportion_var, threshold = thres)

  # Calculate overall estimates
  tarv_1.est <- sum((N_h/N) * tarv_1.bar.hs)
  tarv_2.est <- sum((N_h/N) * tarv_2.bar.hs)
  
  # Calculate standard errors
  tarv_1.se <- sqrt(sum((N_h/N)^2 * (1 - n_hs/N_h) * tarv_1.var.smp.hs/n_hs))
  tarv_2.se <- sqrt(sum((N_h/N)^2 * (1 - n_hs/N_h) * tarv_2.var.smp.hs/n_hs))
  
  # Create confidence intervals
  tarv_1.CI <- c(tarv_1.est - 1.96 * tarv_1.se,
               tarv_1.est + 1.96 * tarv_1.se)
  tarv_2.CI <- c(tarv_2.est - 1.96 * tarv_2.se,
               tarv_2.est + 1.96 * tarv_2.se)
  
  return(list(
    tarv_1 = list(
      estimate = tarv_1.est,
      se = tarv_1.se,
      ci = tarv_1.CI,
      strata_estimates = tarv_1.bar.hs,
      strata_ns = n_hs
    ),
    tarv_2 = list(
      estimate = tarv_2.est,
      se = tarv_2.se,
      ci = tarv_2.CI,
      strata_estimates = tarv_2.bar.hs,
      strata_ns = n_hs
    )
  ))
}







