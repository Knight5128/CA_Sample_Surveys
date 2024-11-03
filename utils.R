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
# Simple Random Sampling
srs_sampling <- function(data, tarv_1, tarv_2, thres, n) {
  # Input validation
  if (!all(c(tarv_1, tarv_2) %in% names(data))) {
    stop("Target variables not found in dataset")
  }
  
  N <- nrow(data)
  if (n >= N) {
    stop("Sample size must be less than population size")
  }
  
  # Take simple random sample
  sample.ids <- sample(N, size = n, replace = FALSE)
  sample_data <- data[sample.ids, ]
  
  # Calculate estimates
  tarv_1.est <- mean(sample_data[[tarv_1]])
  tarv_2.est <- calculate_proportion(sample_data[[tarv_2]], thres)
  
  # Calculate standard errors
  tarv_1.var.smp <- var(sample_data[[tarv_1]])
  tarv_1.se <- sqrt((1 - n/N) * tarv_1.var.smp/n)
  
  tarv_2.var.smp <- tarv_2.est * (1 - tarv_2.est)
  tarv_2.se <- sqrt((1 - n/N) * tarv_2.var.smp/n)
  
  # Create confidence intervals
  tarv_1.CI <- c(tarv_1.est - 1.96 * tarv_1.se, 
               tarv_1.est + 1.96 * tarv_1.se)
  tarv_2.CI <- c(tarv_2.est - 1.96 * tarv_2.se,
               tarv_2.est + 1.96 * tarv_2.se)
  
  return(list(
    tarv_1 = list(
      estimate = tarv_1.est,
      se = tarv_1.se,
      ci = tarv_1.CI
    ),
    tarv_2 = list(
      estimate = tarv_2.est,
      se = tarv_2.se,
      ci = tarv_2.CI
    )
  ))
}


##############################################################################
# Stratified Sampling with Proportional Allocation
str_prop_sampling <- function(data, tarv_1, tarv_2, thres, n, strata_var) {
  # Input validation
  if (!all(c(tarv_1, tarv_2, strata_var) %in% names(data))) {
    stop("Required variables not found in dataset")
  }
  
  # Calculate stratum sizes
  N <- nrow(data)
  N_hs <- table(data[[strata_var]])
  strata_names <- names(N_hs)
  
  # Calculate proportional allocation
  n_hs <- round((N_hs/N) * n)
  
  # Initialize results storage
  tarv_1.bar.hs <- numeric(length(strata_names))
  tarv_2.bar.hs <- numeric(length(strata_names))
  tarv_1.var.smp.hs <- numeric(length(strata_names))
  tarv_2.var.smp.hs <- numeric(length(strata_names))
  
  # Sample from each stratum
  STR.sample <- data.frame()
  
  for (i in seq_along(strata_names)) {
    row.ids <- data[data[[strata_var]] == strata_names[i], ]
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


####################################################################################
#————————————————————————————————————Simulation————————————————————————————————————#
####################################################################################


##############################################################################
compare_sampling_methods <- function(data, tarv_1, tarv_2, thres, n, strata_var, 
                                   n_simulations = 1000, 
                                   strata_s2_guess = NULL, 
                                   strata_cost = NULL) {
  # Calculate true population parameters
  true_mean <- mean(data[[tarv_1]])
  true_prop <- mean(data[[tarv_2]] > thres)
  
  # Initialize storage for results
  results <- list(
    srs = list(
      tarv_1 = data.frame(estimate = numeric(n_simulations), 
                         se = numeric(n_simulations),
                         ci_lower = numeric(n_simulations),
                         ci_upper = numeric(n_simulations)),
      tarv_2 = data.frame(estimate = numeric(n_simulations),
                         se = numeric(n_simulations),
                         ci_lower = numeric(n_simulations),
                         ci_upper = numeric(n_simulations))
    ),
    str_prop = list(
      tarv_1 = data.frame(estimate = numeric(n_simulations),
                         se = numeric(n_simulations),
                         ci_lower = numeric(n_simulations),
                         ci_upper = numeric(n_simulations)),
      tarv_2 = data.frame(estimate = numeric(n_simulations),
                         se = numeric(n_simulations),
                         ci_lower = numeric(n_simulations),
                         ci_upper = numeric(n_simulations))
    ),
    str_equal = list(
      tarv_1 = data.frame(estimate = numeric(n_simulations),
                         se = numeric(n_simulations),
                         ci_lower = numeric(n_simulations),
                         ci_upper = numeric(n_simulations)),
      tarv_2 = data.frame(estimate = numeric(n_simulations),
                         se = numeric(n_simulations),
                         ci_lower = numeric(n_simulations),
                         ci_upper = numeric(n_simulations))
    ),
    str_opt = list(
      tarv_1 = data.frame(estimate = numeric(n_simulations),
                         se = numeric(n_simulations),
                         ci_lower = numeric(n_simulations),
                         ci_upper = numeric(n_simulations)),
      tarv_2 = data.frame(estimate = numeric(n_simulations),
                         se = numeric(n_simulations),
                         ci_lower = numeric(n_simulations),
                         ci_upper = numeric(n_simulations))
    )
  )
  
  # Run simulations
  for(i in 1:n_simulations) {
    # Simple Random Sampling
    srs_result <- srs_sampling(data, tarv_1, tarv_2, thres, n)
    
    # Stratified Sampling - Proportional
    str_prop_result <- str_prop_sampling(data, tarv_1, tarv_2, thres, n, strata_var)
    
    # Stratified Sampling - Equal
    str_equal_result <- str_equal_sampling(data, tarv_1, tarv_2, thres, n, strata_var)
    
    # Stratified Sampling - Optimal
    str_opt_result <- str_opt_sampling(data, tarv_1, tarv_2, thres, n, strata_var,
                                     strata_s2_guess, strata_cost)
    
    # Store results for tarv_1
    results$srs$tarv_1[i,] <- c(srs_result$tarv_1$estimate, 
                               srs_result$tarv_1$se,
                               srs_result$tarv_1$ci)
    results$str_prop$tarv_1[i,] <- c(str_prop_result$tarv_1$estimate,
                                    str_prop_result$tarv_1$se,
                                    str_prop_result$tarv_1$ci)
    results$str_equal$tarv_1[i,] <- c(str_equal_result$tarv_1$estimate,
                                     str_equal_result$tarv_1$se,
                                     str_equal_result$tarv_1$ci)
    results$str_opt$tarv_1[i,] <- c(str_opt_result$tarv_1$estimate,
                                   str_opt_result$tarv_1$se,
                                   str_opt_result$tarv_1$ci)
    
    # Store results for tarv_2
    results$srs$tarv_2[i,] <- c(srs_result$tarv_2$estimate,
                               srs_result$tarv_2$se,
                               srs_result$tarv_2$ci)
    results$str_prop$tarv_2[i,] <- c(str_prop_result$tarv_2$estimate,
                                    str_prop_result$tarv_2$se,
                                    str_prop_result$tarv_2$ci)
    results$str_equal$tarv_2[i,] <- c(str_equal_result$tarv_2$estimate,
                                     str_equal_result$tarv_2$se,
                                     str_equal_result$tarv_2$ci)
    results$str_opt$tarv_2[i,] <- c(str_opt_result$tarv_2$estimate,
                                   str_opt_result$tarv_2$se,
                                   str_opt_result$tarv_2$ci)
  }
  
  # Add true values to results
  results$true_values <- list(
    tarv_1 = true_mean,
    tarv_2 = true_prop
  )
  
  return(results)
}




####################################################################################
#————————————————————————————————————Visualization—————————————————————————————————#
####################################################################################


##############################################################################
plot_continuous_distribution <- function(data, tarv_1, strata_var = NULL) {
  # Basic density plot
  p1 <- ggplot(data, aes(x = .data[[tarv_1]])) +
    geom_density(fill = "lightblue", alpha = 0.5) +
    geom_rug() +
    labs(
      title = paste("Distribution of", tarv_1),
      x = tarv_1,
      y = "Density"
    ) +
    theme_minimal()
  
  # Add boxplot for better outlier visualization
  p2 <- ggplot(data, aes(x = "", y = .data[[tarv_1]])) +
    geom_boxplot(fill = "lightblue", alpha = 0.5) +
    labs(
      title = "Boxplot",
      x = "",
      y = tarv_1
    ) +
    theme_minimal()
  
  # If strata variable is provided, add stratified visualization
  if (!is.null(strata_var)) {
    p3 <- ggplot(data, aes(x = .data[[strata_var]], y = .data[[tarv_1]])) +
      geom_boxplot(fill = "lightblue", alpha = 0.5) +
      labs(
        title = paste("Distribution by", strata_var),
        x = strata_var,
        y = tarv_1
      ) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    return(list(density = p1, boxplot = p2, stratified = p3))
  }
  
  return(list(density = p1, boxplot = p2))
}


##############################################################################
plot_binary_distribution <- function(data, tarv_2, thres, strata_var = NULL) {
  # Create binary variable
  data$binary_var <- ifelse(data[[tarv_2]] > thres, "Above", "Below")
  
  # Overall proportion plot
  p1 <- ggplot(data, aes(x = binary_var, fill = binary_var)) +
    geom_bar(aes(y = after_stat(prop), group = 1)) +
    scale_y_continuous(labels = scales::percent) +
    labs(
      title = paste("Proportion Distribution for", tarv_2),
      subtitle = paste("Threshold =", thres),
      x = "",
      y = "Percentage",
      fill = paste(tarv_2, "Value")
    ) +
    theme_minimal()
  
  # If strata variable is provided, add stratified visualization
  if (!is.null(strata_var)) {
    p2 <- ggplot(data, aes(x = .data[[strata_var]], fill = binary_var)) +
      geom_bar(position = "fill") +
      scale_y_continuous(labels = scales::percent) +
      labs(
        title = paste("Proportion Distribution by", strata_var),
        x = strata_var,
        y = "Percentage",
        fill = paste(tarv_2, "Value")
      ) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    return(list(overall = p1, stratified = p2))
  }
  
  return(list(overall = p1))
}


##############################################################################
plot_simulation_results <- function(sim_results, variable = "tarv_1") {
  # Prepare data for plotting
  methods <- c("SRS", "Stratified-Proportional", "Stratified-Equal", "Stratified-Optimal")
  true_value <- sim_results$true_values[[variable]]
  
  # Combine all estimates into one dataframe
  estimates_df <- data.frame(
    estimate = c(sim_results$srs[[variable]]$estimate,
                sim_results$str_prop[[variable]]$estimate,
                sim_results$str_equal[[variable]]$estimate,
                sim_results$str_opt[[variable]]$estimate),
    method = rep(methods, each = nrow(sim_results$srs[[variable]]))
  )
  
  # Create density plot
  p1 <- ggplot(estimates_df, aes(x = estimate, fill = method)) +
    geom_density(alpha = 0.5) +
    geom_vline(xintercept = true_value, linetype = "dashed", color = "red") +
    labs(title = paste("Distribution of Estimates -", variable),
         x = "Estimate",
         y = "Density") +
    theme_minimal()
  
  # Calculate coverage probabilities
  coverage <- data.frame(method = methods)
  coverage$coverage <- c(
    mean(sim_results$srs[[variable]]$ci_lower <= true_value & 
         sim_results$srs[[variable]]$ci_upper >= true_value),
    mean(sim_results$str_prop[[variable]]$ci_lower <= true_value & 
         sim_results$str_prop[[variable]]$ci_upper >= true_value),
    mean(sim_results$str_equal[[variable]]$ci_lower <= true_value & 
         sim_results$str_equal[[variable]]$ci_upper >= true_value),
    mean(sim_results$str_opt[[variable]]$ci_lower <= true_value & 
         sim_results$str_opt[[variable]]$ci_upper >= true_value)
  )
  
  # Create coverage plot
  p2 <- ggplot(coverage, aes(x = method, y = coverage, fill = method)) +
    geom_bar(stat = "identity") +
    geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
    labs(title = "95% CI Coverage Probability",
         x = "Method",
         y = "Coverage Probability") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(list(density_plot = p1, coverage_plot = p2))
}



