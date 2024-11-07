library(dplyr)
library(stats)
library(ggplot2)

####################################################################################
#——————————————————————————————————————Basic——————————————————————————————————————#
####################################################################################

calculate_proportion <- function(x, threshold) {
  mean(x > threshold)
}

calculate_proportion_var <- function(x, threshold) {
  p <- mean(x > threshold)
  var <- p * (1 - p)
  return(var)
}

####################################################################################
#—————————————————————————————————————Sampling—————————————————————————————————————#
####################################################################################


##############################################################################
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

##############################################################################
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


##############################################################################
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
  tarv_1.se <- sqrt((1 - n/N) * var_d/(n * X_bar^2))
  
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


##############################################################################
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
    row.ids <- which(data[[strata_var]] == strata_names[i])
    sample.ids <- sample(row.ids, size = n_hs[i], replace = FALSE)
    STR.sample <- rbind(STR.sample, data[sample.ids,])
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


##############################################################################
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
    row.ids <- which(data[[strata_var]] == strata_names[i])
    sample.ids <- sample(row.ids, size = n_hs[i], replace = FALSE)
    STR.sample <- rbind(STR.sample, data[sample.ids,])
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


##############################################################################
##############################################################################
# Simulation function for comparing **mean** estimation methods
compare_mean_estimation_methods <- function(data, tarv, aux_var, n, strata_var, 
                                         n_simulations = 1000) {
  # Calculate true population parameter
  true_mean <- mean(data[[tarv]])
  
  # Initialize storage for results
  results <- list(
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
  )
  
  # Run simulations
  for(i in 1:n_simulations) {
    # Simple Random Sampling - Vanilla
    srs_vanilla <- srs_sampling_est_mean_vanilla(data, tarv, n)
    
    # Simple Random Sampling - Ratio
    srs_ratio <- srs_sampling_est_mean_ratio(data, tarv, aux_var, n)
    
    # Stratified Sampling - Proportional
    str_prop <- str_prop_sampling_est_mean_vanilla(data, tarv, n, strata_var)
    
    # Store results
    results$srs_vanilla[i,] <- c(
      srs_vanilla$estimate,
      srs_vanilla$se,
      srs_vanilla$ci
    )
    
    results$srs_ratio[i,] <- c(
      srs_ratio$estimate,
      srs_ratio$se,
      srs_ratio$ci
    )
    
    results$str_prop[i,] <- c(
      str_prop$estimate,
      str_prop$se,
      str_prop$ci
    )
  }
  
  # Calculate performance metrics
  performance <- list()
  methods <- c("srs_vanilla", "srs_ratio", "str_prop")
  
  for(method in methods) {
    performance[[method]] <- list(
      # Bias
      bias = mean(results[[method]]$estimate - true_mean),
      # Root Mean Square Error
      rmse = sqrt(mean((results[[method]]$estimate - true_mean)^2)),
      # Coverage probability
      coverage = mean(results[[method]]$ci_lower <= true_mean & 
                     results[[method]]$ci_upper >= true_mean),
      # Average width of confidence intervals
      ci_width = mean(results[[method]]$ci_upper - results[[method]]$ci_lower)
    )
  }
  
  # Add true value and simulation results to output
  results$true_value <- true_mean
  results$performance <- performance
  
  # Create summary visualizations
  plots <- plot_mean_simulation_results(results)
  results$plots <- plots
  
  return(results)
}


##############################################################################
##############################################################################
# Simulation function for comparing **proportion** estimation methods
compare_prop_estimation_methods <- function(data, tarv, thres, n, strata_var, 
                                         n_simulations = 1000) {
  # Calculate true population parameter
  true_prop <- mean(data[[tarv]] > thres)
  
  # Initialize storage for results
  results <- list(
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
  
  # Run simulations
  for(i in 1:n_simulations) {
    # Simple Random Sampling - Vanilla
    srs_vanilla <- srs_sampling_est_prop_vanilla(data, tarv, thres, n)
    
    # Stratified Sampling - Proportional
    str_prop <- str_prop_sampling_est_prop_vanilla(data, tarv, thres, n, strata_var)
    
    # Store results
    results$srs_vanilla[i,] <- c(
      srs_vanilla$estimate,
      srs_vanilla$se,
      srs_vanilla$ci
    )
    
    results$str_prop[i,] <- c(
      str_prop$estimate,
      str_prop$se,
      str_prop$ci
    )
  }
  
  # Calculate performance metrics
  performance <- list()
  methods <- c("srs_vanilla", "str_prop")
  
  for(method in methods) {
    performance[[method]] <- list(
      # Bias
      bias = mean(results[[method]]$estimate - true_prop),
      # Root Mean Square Error
      rmse = sqrt(mean((results[[method]]$estimate - true_prop)^2)),
      # Coverage probability
      coverage = mean(results[[method]]$ci_lower <= true_prop & 
                     results[[method]]$ci_upper >= true_prop),
      # Average width of confidence intervals
      ci_width = mean(results[[method]]$ci_upper - results[[method]]$ci_lower)
    )
  }
  
  # Add true value and simulation results to output
  results$true_value <- true_prop
  results$performance <- performance
  
  # Create summary visualizations
  plots <- plot_prop_simulation_results(results)
  results$plots <- plots
  
  return(results)
}





####################################################################################
#————————————————————————————————————Visualization—————————————————————————————————#
####################################################################################


##############################################################################
##############################################################################
# Helper function to plot **mean** simulation results
plot_mean_simulation_results <- function(results) {
  # Prepare data for plotting
  methods <- c("SRS-Vanilla", "SRS-Ratio", "Stratified-Proportional")
  
  # Combine all estimates into one dataframe for density plot
  estimates_df <- data.frame(
    estimate = c(results$srs_vanilla$estimate,
                results$srs_ratio$estimate,
                results$str_prop$estimate),
    method = rep(methods, each = nrow(results$srs_vanilla))
  )
  
  # Create density plot
  p1 <- ggplot(estimates_df, aes(x = estimate, fill = method)) +
    geom_density(alpha = 0.5) +
    geom_vline(xintercept = results$true_value, linetype = "dashed", color = "red") +
    labs(title = "Distribution of Estimates",
         x = "Estimate",
         y = "Density") +
    theme_minimal()
  
  # Prepare performance metrics for bar plots
  performance_df <- data.frame(
    method = methods,
    bias = c(results$performance$srs_vanilla$bias,
             results$performance$srs_ratio$bias,
             results$performance$str_prop$bias),
    rmse = c(results$performance$srs_vanilla$rmse,
             results$performance$srs_ratio$rmse,
             results$performance$str_prop$rmse),
    coverage = c(results$performance$srs_vanilla$coverage,
                results$performance$srs_ratio$coverage,
                results$performance$str_prop$coverage),
    ci_width = c(results$performance$srs_vanilla$ci_width,
                 results$performance$srs_ratio$ci_width,
                 results$performance$str_prop$ci_width)
  )
  
  # Create RMSE plot
  p2 <- ggplot(performance_df, aes(x = method, y = rmse, fill = method)) +
    geom_bar(stat = "identity") +
    labs(title = "Root Mean Square Error",
         x = "Method",
         y = "RMSE") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Create coverage probability plot
  p3 <- ggplot(performance_df, aes(x = method, y = coverage, fill = method)) +
    geom_bar(stat = "identity") +
    geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
    labs(title = "95% CI Coverage Probability",
         x = "Method",
         y = "Coverage Probability") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(list(
    density_plot = p1,
    rmse_plot = p2,
    coverage_plot = p3
  ))
}


##############################################################################
############################################################################
# Helper function to plot **proportion** simulation results
plot_prop_simulation_results <- function(results) {
  # Prepare data for plotting
  methods <- c("SRS-Vanilla", "Stratified-Proportional")
  
  # Combine all estimates into one dataframe for density plot
  estimates_df <- data.frame(
    estimate = c(results$srs_vanilla$estimate,
                results$str_prop$estimate),
    method = rep(methods, each = nrow(results$srs_vanilla))
  )
  
  # Create density plot
  p1 <- ggplot(estimates_df, aes(x = estimate, fill = method)) +
    geom_density(alpha = 0.5) +
    geom_vline(xintercept = results$true_value, linetype = "dashed", color = "red") +
    labs(title = "Distribution of Proportion Estimates",
         x = "Estimate",
         y = "Density") +
    theme_minimal()
  
  # Prepare performance metrics for bar plots
  performance_df <- data.frame(
    method = methods,
    bias = c(results$performance$srs_vanilla$bias,
             results$performance$str_prop$bias),
    rmse = c(results$performance$srs_vanilla$rmse,
             results$performance$str_prop$rmse),
    coverage = c(results$performance$srs_vanilla$coverage,
                results$performance$str_prop$coverage),
    ci_width = c(results$performance$srs_vanilla$ci_width,
                 results$performance$str_prop$ci_width)
  )
  
  # Create RMSE plot
  p2 <- ggplot(performance_df, aes(x = method, y = rmse, fill = method)) +
    geom_bar(stat = "identity") +
    labs(title = "Root Mean Square Error",
         x = "Method",
         y = "RMSE") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Create coverage probability plot
  p3 <- ggplot(performance_df, aes(x = method, y = coverage, fill = method)) +
    geom_bar(stat = "identity") +
    geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
    labs(title = "95% CI Coverage Probability",
         x = "Method",
         y = "Coverage Probability") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Create CI width plot
  p4 <- ggplot(performance_df, aes(x = method, y = ci_width, fill = method)) +
    geom_bar(stat = "identity") +
    labs(title = "Average Confidence Interval Width",
         x = "Method",
         y = "CI Width") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(list(
    density_plot = p1,
    rmse_plot = p2,
    coverage_plot = p3,
    ci_width_plot = p4
  ))
}




####################################################################################
#————————————————————————————————Other Sampling————————————————————————————————————#
####################################################################################

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


##############################################################################
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



}



