library(dplyr)
library(stats)
library(ggplot2)

Simple_Random_Sampling <- function(df, n, seed = NULL) {
    # Set seed for reproducibility, if provided
    if (!is.null(seed)) {
        set.seed(seed)
    }
    
    # Perform simple random sampling
    SRS.smp <- df %>%
        slice_sample(n = n)
    
    return(SRS.smp)
}

Stratified_Sampling_Proportional_Allocation <- function(
    df,
    strata_col,
    n,
    seed = NULL
) {
    # Create stratified random sample
    if (!is.null(seed)) {
        set.seed(seed)
    }
    
    # Calculate proportional sizes for each stratum
    strata_props <- df %>%
        count(!!sym(strata_col)) %>%
        mutate(prop = n / sum(n),
               size = round(prop * n))
    
    # Adjust for rounding errors
    total <- sum(strata_props$size)
    if (total != n) {
        diff <- n - total
        strata_props$size[1] <- strata_props$size[1] + diff
    }
    
    # Sample from each stratum
    samples <- lapply(1:nrow(strata_props), function(i) {
        stratum_name <- strata_props[[strata_col]][i]
        size <- strata_props$size[i]
        
        df %>%
            filter(!!sym(strata_col) == stratum_name) %>%
            slice_sample(n = size)
    })
    
    # Combine all strata
    bind_rows(samples)
}

Stratified_Sampling_Optimal_Allocation <- function(
    df,
    strata_col,
    n,
    seed = NULL
) {

}