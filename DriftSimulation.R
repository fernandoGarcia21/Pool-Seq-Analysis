setwd('C:/Users/dgarciac/Documents/PhD Project/Third Chapter/PoolSeq-Analyses/')


######################################################################
# Function that estimates the size of the population at generation p_T
# according to the logistic growth in the neutral model of the  
# demographic inference
######################################################################
population_logistic_growth <- function(p_T, p_set_parameters){
  n_T = 0
  tmp_divisor = 1+((p_set_parameters[['K']] - p_set_parameters[['N0']])/p_set_parameters[['N0']])*exp(-p_set_parameters[['r']]*p_T)
  n_T = p_set_parameters$K / tmp_divisor
  return (n_T)
}



#Fixed values
N0 = 400 #Initial effective population size
K = 2000 #Karrying capacity
growth_rate = 0.075 #r 

total_generations = 60 #REDWOOD
N_g = N0 #The size of the population, it will increase over each generation


list_parameters = list()
list_parameters[['K']] = K
list_parameters[['N0']] = N0
list_parameters[['r']] = growth_rate

array_drift = c()
array_drift = append(array_drift, 1/N_g)
for(g in seq(1:total_generations)){
  
  #The population size at generation g with rate r
  N_g = population_logistic_growth(g, list_parameters)
  cat('Drift at population size  ', N_g, 'is', 1/N_g, '\n')
  array_drift = append(array_drift, 1/N_g)
}

cat('Total drift ', sum(array_drift), '\n')
















# 1. Define the logistic growth function (as you already have it)
population_logistic_growth <- function(p_T, p_set_parameters){
  n_T = 0
  # Ensure N0 is not zero to prevent division by zero in the denominator
  if (p_set_parameters[['N0']] == 0) {
    stop("Initial population size (N0) cannot be zero.")
  }
  tmp_divisor = 1 + ((p_set_parameters[['K']] - p_set_parameters[['N0']]) / p_set_parameters[['N0']]) * exp(-p_set_parameters[['r']] * p_T)
  n_T = p_set_parameters$K / tmp_divisor
  return (n_T)
}

# 2. Create a new function to calculate total drift for given K and r
calculate_total_drift <- function(N0, K, r, total_generations) {
  list_parameters = list(
    'K' = K,
    'N0' = N0,
    'r' = r
  )
  
  # Initialize array for drift values (1/Ne for each generation)
  array_drift = c()
  
  # Genetic drift is inversely proportional to effective population size (Ne)
  # Assuming N_g (population size at generation g) represents Ne for drift calculation
  # Drift at generation 0 (N0)
  array_drift = append(array_drift, 1 / N0)
  
  # Calculate population size and drift for subsequent generations
  for(g in 1:total_generations){ # Use 1:total_generations for sequence
    # The population size at generation g with rate r
    N_g = population_logistic_growth(g, list_parameters)
    
    # In some drift calculations, if N_g falls below 1, it's problematic.
    # We'll assume N_g represents effective population size here, and it should be >= 1.
    # If N_g can theoretically be less than 1 (e.g., due to rounding or model limits),
    # you might want to handle it (e.g., set 1/N_g to Inf or a very large number).
    # For this model, N_g typically grows towards K, so it should stay positive.
    if (N_g < 1) {
      warning(paste("Population size (N_g) fell below 1 at generation", g, "for K=", K, "r=", r, ". Setting drift to Inf for this generation."))
      current_drift = Inf
    } else {
      current_drift = 1 / N_g
    }
    
    # cat('Drift at population size ', N_g, 'is', current_drift, '\n') # Optional: uncomment for per-generation details
    array_drift = append(array_drift, current_drift)
  }
  
  total_drift = sum(array_drift)
  return(total_drift)
}


#Population 
pop = 'RED'
generations = list('RED' = 60,'OAK' = 50.8)
actual_drift = list('RED' = 0.15,'OAK' = 0.13)
# 3. Define the ranges for K and r you want to test
N0 = 200 # Initial effective population size (fixed as per your original code)
total_generations = generations[[pop]] # Number of generations to simulate

# Define combinations for K and r
# You can make these vectors as long as you need
k_values = c(500, 800, 1000, 2000, 5000)
r_values = c(0.015, 0.015, 0.02, 0.025, 0.05, 0.075, 0.1, 0.15)

# 4. Store results
# Create a data frame to store the results
results_df = data.frame(K = numeric(),
                        r = numeric(),
                        Total_Drift = numeric())

# 5. Loop through combinations of K and r
for (k_val in k_values) {
  for (r_val in r_values) {
    cat(paste0("Calculating drift for K = ", k_val, ", r = ", r_val, "...\n"))
    drift_result = calculate_total_drift(N0, k_val, r_val, total_generations)
    
    # Add the results to the data frame
    results_df = rbind(results_df, data.frame(K = k_val, r = r_val, Total_Drift = drift_result))
  }
}

# 6. View the results
print("Drift Estimation Results:")
print(results_df)

# Optional: You can sort the results or analyze them further
# For example, to see which combination resulted in the lowest drift:
# results_df[order(results_df$Total_Drift), ]


 library(ggplot2)

 ggplot(results_df, aes(x = r, y = Total_Drift, color = as.factor(K))) +
   geom_hline(aes(yintercept = actual_drift[[pop]]), color = "darkblue", size = 0.5, linetype='dashed')+
   geom_line() +
   geom_point() +
   ylim(0.04,0.25)+
   scale_x_continuous(
     limits = c(0.015, 0.16),                 # Set the exact range of the x-axis
     breaks = seq(0.015, 0.16, by = 0.015)        # Set tick marks at every integer from 0 to 10
   ) +
   labs(title = paste("Total Drift Across Different K and r Values in" , pop),
        x = "Growth Rate (r)",
        y = "Total Drift (Sum of 1/Ne)",
        color = "Carrying Capacity (K)") +
   theme_minimal()
 
 