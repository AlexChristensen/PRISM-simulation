#%%%%%%%%%%%%%%%%%%%%%%%#
#### Data Generation ####
#%%%%%%%%%%%%%%%%%%%%%%%#

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
source("./source/generate_data.R")

# Create directory
dir.create("./data/")
dir.create("./data/extreme")
dir.create("./data/moderate")

# Set structure names
structure_names <- c("no_overlap", "no_overlap_error", "overlap", "overlap_error")

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

# Set sample sizes
N <- c(500, 1000, 2500, 5000, 10000)

# Number of conditions
n_conditions <- nrow(conditions)
n_digits <- latentFactoR:::digits(n_conditions) - 1

# Set files
FILES <- list.files("./conditions/")

# Loop over conditions
for(i in seq_along(FILES)){

  # Obtain condition
  condition <- conditions[i,]

  # Set lower factors
  lower_factors <- condition$LowF * condition$HighF

  # Set total variables
  total_variables <- condition$LowV * lower_factors
  total_sequence <- seq_len(total_variables)
  variable_names <- paste0("V", total_sequence)
  Sigma <- diag(total_variables)

  # Set condition number
  condition <- cbind.data.frame(
    number = latentFactoR:::format_integer(i, n_digits), condition
  )

  # Load replicates
  load(paste0("./conditions/", FILES[i]))

  # Generate replicate data
  invisible(mclapply(
    seq_along(replicate), function(j){

      # Set condition structure
      condition_structure <- replicate[[j]]

      # Remove bad generation check
      condition_structure <- condition_structure[structure_names]

      # Generate extreme skew data
      extreme_data <- generate_data(
        total_variables, condition_structure, variable_names, N,
        Sigma, total_sequence, skew = "extreme"
      )

      # Save data
      save(
        extreme_data,
        file = paste0(
          "./data/extreme/condition-", condition$number,
          "_replicate-", EGAnet:::format_integer(j, 2),
          ".RData"
        )
      )

      # Generate moderate skew data
      moderate_data <- generate_data(
        total_variables, condition_structure, variable_names, N,
        Sigma, total_sequence, skew = "moderate"
      )

      # Save data
      save(
        moderate_data,
        file = paste0(
          "./data/moderate/condition-", condition$number,
          "_replicate-", EGAnet:::format_integer(j, 2),
          ".RData"
        )
      )

      # Nothing to return; data is saved to disk
      invisible(NULL)

    }, mc.cores = 50, mc.preschedule = FALSE
  ))

  # Free memory
  gc(verbose = FALSE)
  if((i %% 10) == 0){
    system("sudo sh -c 'echo 3 > /proc/sys/vm/drop_caches'")
  }

  # Print update
  cat(paste0("\r Run ", condition$number, " of ", n_conditions, " completed."))

}
