# Load required libraries
library(gtools)          # For generating permutations
library(PlackettLuce)    # For fitting Plackett-Luce models
library(parallel)        # For parallel computing
library(ggplot2)         # For visualization

# Generate all possible rankings (permutations) for a given number of items
generate_permutations <- function(n_items)
{
  return (gtools::permutations(n = n_items, r = n_items))
}

# Calculate sampling weights for rankings based on bias factors
sampling_weights <- function(permutations_set, 
                             biased_items = c(1), # Positions of items to be sampled with bias (default: first item) 
                             bias_factors = c(1)  # By default, no weighted sampling 
                            ) 
{
  
  if (length(biased_items) != length(bias_factors)) {
    stop("The number of items to be sampled with bias does not match the number of weights. 
         Please provide a weight for each item you want to bias.")
  }
  
  # Initialize weights (uniform sampling by default)
  weights <- rep(1, nrow(permutations_set))
  
  # Adjust weights of items based on their specified bias factors
  # Example: If an item is twice as likely to appear first, its weight is increased
  for (i in seq_along(biased_items)) {
    weights[permutations_set[, 1] == biased_items[i]] <- bias_factors[i]
  }
  
  # Normalize weights to sum to 1
  weights <- weights / sum(weights)
  
  return(weights)
}

# Function to sample rankings based on provided weights
sample_permutations <- function(permutations_set, weights, N) {
  
  # Weighted sampling of N row indices 
  # (row indexes indicate which permutation to sample)
  sampled_indices <- sample(
    nrow(permutations_set), 
    size = N, 
    replace = TRUE, 
    prob = weights
  )
  
  # Extract the sampled rows (i.e., the actual rankings/permutations)
  sampled_permutations <- permutations_set[sampled_indices, ]
  
  # Convert to data frame and add frequency count
  result <- as.data.frame(sampled_permutations)
  names(result) <- paste0("Position", 1:ncol(result))
  result$Frequency <- 1
  
  # Aggregate to get unique permutations with their frequencies
  result <- aggregate(Frequency ~ ., data = result, FUN = sum)
  
  return(result)
}

# Function to calculate statistical power for the Plackett-Luce model
power_plackett_luce <- function(
    min_effect_size,
    n_items,
    result = NULL,
    bonferroni = TRUE
) 
{
  if min_effect_size < 0 {
    stop("min_effect_size must be positive or zero. We are testing for those
         effect sizes which are strictly greater than min_effect_size or 
         strictly lower than -min_effect_size")
  }
  
  if (is.null(result)) {
    stop("You must provide a distribution over rankings/permutations.")
  }
  
  # Ranking object for PlackettLuce
  R <- PlackettLuce::as.rankings(
    result[,-(n_items + 1)], 
    input = "orderings", 
    items = attr(result, "items")
  )
  
  # Fit PlackettLuce model
  mod <- tryCatch({
    PlackettLuce::PlackettLuce(R, weights = result$Frequency)
  }, error = function(e) {
    stop("Error in fitting the Plackett-Luce model: ", e$message)
  })
  
  mod_summary <- summary(mod)
  
  # Extract effect sizes and standard errors (excluding the reference category)
  # For simplicity, we assume the references category is "1", the first item
  effect_sizes <- mod_summary$coefficients[-1, "Estimate"]
  standard_errors <- mod_summary$coefficients[-1, "Std. Error"]
  
  # Filter point estimates and SE of effect sizes based on min_effect_size
  filtered_indices <- which((effect_sizes > min_effect_size) | (effect_sizes < -min_effect_size))
  filtered_effect_sizes <- effect_sizes[filtered_indices]
  filtered_standard_errors <- standard_errors[filtered_indices]
  
  # Calculate Bonferroni corrected intervals for filtered effect sizes
  if (bonferroni) {
    no_comparisons <- n_items - 1
    critical_value <- qnorm(1 - (0.05 / no_comparisons) / 2)
  } else {
    critical_value <- qnorm(0.975)
  }
  
  # CI of effect sizes  
  lower_bound <- filtered_effect_sizes - critical_value * filtered_standard_errors
  upper_bound <- filtered_effect_sizes + critical_value * filtered_standard_errors
  
  # Check if MEOF falls within confidence intervals
  meof_in_ci <- ((min_effect_size >= lower_bound) & (min_effect_size <= upper_bound)) |
    ((-min_effect_size >= lower_bound) & (-min_effect_size <= upper_bound))
  
  # Report effect sizes whose CI do nt include the MEOF 
  return(sum(meof_in_ci == FALSE))
}


# Set up a cluster with the number of cores you want to use
num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)

# Export necessary functions and variables to the cluster

permutations_set <- generate_permutations(n_items = 5)

clusterExport(cl, c("power_plackett_luce", 
                    "sample_permutations", 
                    "generate_permutations",
                    "sampling_weights", 
                    "permutations_set"
)
)

# Run the power analysis in parallel

n_sims <- 1000
power_analysis <- parLapply(
  cl, 
  1:n_sims, 
  function(x) power_plackett_luce(
    min_effect_size = 0.03, # If zero, checks only for stat. significance
    n_items = 5, 
    result = sample_permutations(permutations_set, 
                                 weights = sampling_weights(
                                   permutations_set,
                                   biased_items = c(1, 3, 5), # Biased items
                                   bias_factors = c(2, 1.5, 0.5) # Corresponding biases
                                 ), 
                                 N = 16250
    ), 
    bonferroni = TRUE
  )
)

# Stop the cluster when done
stopCluster(cl)

# Aggregate and visualize the results
formatted_results <- as.data.frame(table(unlist(power_analysis)) / n_sims)
names(formatted_results) <- c("EffectsNotInCI", "Frequency")

power_substantial_effects <- sum(formatted_results$Frequency[as.numeric(as.character(formatted_results$EffectsNotInCI)) >= 1])

# Create a bar plot for the results
ggplot(formatted_results, aes(x = EffectsNotInCI, y = Frequency)) +
  geom_bar(
    stat = "identity",
    aes(fill = ifelse(as.numeric(as.character(EffectsNotInCI)) >= 1, "red", "blue")),
    width = 0.3
  ) +
  scale_fill_identity() +
  labs(
    title = "Reference = 2 times more likely to be ranked first",
    x = "# Substantial Effect Sizes",
    y = "Frequency"
  ) +
  theme_minimal() +
  annotate(
    "text",
    x = 1, # Position annotation at the last x-axis value
    y = max(formatted_results$Frequency) * 1.1,         
    label = paste("Power:", power_substantial_effects),
    color = "red",
    hjust = 1
  )
