#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#### Hierarchical | Simulate Conditions ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Provide password
system("sudo sh -c 'echo 3 > /proc/sys/vm/drop_caches'")

# Make sure environment has single threaded applications
Sys.setenv(BLIS_NUM_THREADS = "1")
Sys.setenv(AOCL_NUM_THREADS = "1")
Sys.setenv(OMP_NUM_THREADS = "1")
RhpcBLASctl::blas_set_num_threads(1)
RhpcBLASctl::omp_set_num_threads(1)

# Load packages
library(latentFactoR); library(parallel)

# Load sources
source("./source/generate_condition.R")
source("./source/methods.R")
source("./source/results.R")

# Create directory
dir.create("./conditions")

# Conditions
conditions <- expand.grid(
  LowV = c(5, 10), # number of variables per lower order factor
  LowF = c(3, 5), # number of lower order factors per higher order factor
  LowL = c(0.50, 0.70), # lower loading sizes
  HighF = c(1, 2, 4, 6), # number of higher order factors
  HighL = c(0.50, 0.70), # higher order loadings
  HighC = c(0.10, 0.30, 0.50) # higher order correlations
)

# Skip redundant conditions for unidimensional
conditions <- conditions[
  !((conditions$HighF == 1) & (conditions$HighC %in% c(0.30, 0.50))),
]

# Number of conditions
n_conditions <- nrow(conditions)
n_digits <- latentFactoR:::digits(n_conditions) - 1

# Initialize results
results <- vector("list", n_conditions)

# Loop over results
for(i in seq_len(n_conditions)){

  # Obtain condition
  condition <- conditions[i,]

  # Set condition number
  condition_number <- latentFactoR:::format_integer(i, n_digits)

  # Generate replicate results
  replicate <- mclapply(
    seq_len(100), function(j){

      # Regenerate with error
      data_condition <- list()
      class(data_condition) <- "try-error"

      # Loop until good condition
      while(is(data_condition, "try-error")){
        data_condition <- try(generate_condition(condition), silent = TRUE)
      }

      # Combine results
      return(data_condition)

    }, mc.cores = 50
  )

  # Save result
  save(replicate, file = paste0("./conditions/condition_", condition_number, ".RData"))

  # Remove replicate
  rm(replicate)
  gc(verbose = FALSE)
  if((i %% 24) == 0){
    system("sudo sh -c 'echo 3 > /proc/sys/vm/drop_caches'")
  }

  # Print update
  cat(paste0("\r Run ", condition_number, " of ", n_conditions, " completed."))

}

3