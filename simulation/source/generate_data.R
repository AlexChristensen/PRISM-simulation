#%%%%%%%%%%%%%%%%%%%%%%%#
#### Data Generation ####
#%%%%%%%%%%%%%%%%%%%%%%%#

# Load packages
library(latentFactoR); library(mvtnorm); library(compiler)

# Generate data
generate_data <- function(total_variables, condition_structure, variable_names, N, Sigma, total_sequence, skew)
{

  # Initialize zero standard deviation
  zero_sd <- TRUE

  # Generate until there is variation
  while(any(zero_sd)){

    # Check for skew
    if(skew == "extreme"){

      # Set skewness
      skewness <- 11

      # Generate skew
      while(any(skewness > 10)){
        skewness <- qgamma(latentFactoR:::runif_xoshiro(total_variables), shape = 1.70, rate = 1)
      }

    }else{

      # Generate skew
      skewness <- latentFactoR:::runif_xoshiro(total_variables, min = 1, max = 2)

    }

    # Generate data and obtain results
    data_condition <- lapply(seq_along(condition_structure), function(k){

      # Obtain Cholesky decomposition
      cholesky <- chol(condition_structure[[k]]$population_correlation)

      # Loop over sample sizes
      data_list <- lapply(N, function(n){

        # Generate data
        data <- rmvnorm(n, sigma = Sigma) %*% cholesky

        # Add skew
        for(l in total_sequence){
          data[,l] <- categorize(data[,l], categories = 5, skew_value = skewness[l])
        }

        # Add variable names
        dimnames(data)[[2]] <- variable_names

        # Return data
        return(data)

      })

      # Set names
      names(data_list) <- N

      # Return data list
      return(data_list)

    })

    # Check that data does not have any zero variance
    zero_sd <- unlist(lapply(data_condition, function(condition_set){
      EGAnet:::lvapply(condition_set, function(data){
        any(apply(data, 2, sd) == 0)
      })
    }))

  }

  # Set names
  names(data_condition) <- names(condition_structure)

  # Return data condition
  return(data_condition)

}

# Compile function
generate_data <- cmpfun(generate_data)