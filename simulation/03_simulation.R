#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#### Hierarchical | Results ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Provide password
system("sudo sh -c 'echo 3 > /proc/sys/vm/drop_caches'")

# Make sure environment has single threaded applications
Sys.setenv(BLIS_NUM_THREADS = "1")
Sys.setenv(AOCL_NUM_THREADS = "1")
Sys.setenv(OMP_NUM_THREADS = "1")
RhpcBLASctl::blas_set_num_threads(1)
RhpcBLASctl::omp_set_num_threads(1)

# Load packages
library(latentFactoR); library(parallel); library(tidyverse)

# Load sources
source("./source/methods.R")
source("./source/results.R")
source("./source/obtain_results.R")

# Create directory
dir.create("./timing/")
dir.create("./results/")

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

# Loop over sample sizes first so every condition (and condition structure)
# finishes at the current N before moving on to the next N
for(n in N[-c(1:2)]){

  # Loop over conditions
  for(i in seq_along(FILES)){

    # Obtain condition
    condition <- conditions[i,]

    # Set lower factors
    lower_factors <- condition$LowF * condition$HighF

    # Set condition number
    condition <- cbind.data.frame(
      number = latentFactoR:::format_integer(i, n_digits), condition
    )

    # Create condition list to attach
    condition_df <- do.call(rbind.data.frame, lapply(seq_len(8), function(x){condition}))

    # Load replicates
    load(paste0("./conditions/", FILES[i]))

    # Generate replicate results
    result <- mclapply(
      seq_along(replicate), function(j){

        # Set condition structure
        condition_structure <- replicate[[j]]

        # Remove bad generation check
        condition_structure <- condition_structure[structure_names]

        # Load pre-generated extreme skew data
        load(paste0(
          "./data/extreme/condition-", condition$number,
          "_replicate-", EGAnet:::format_integer(j, 2),
          ".RData"
        ))

        # Load pre-generated moderate skew data
        load(paste0(
          "./data/moderate/condition-", condition$number,
          "_replicate-", EGAnet:::format_integer(j, 2),
          ".RData"
        ))

        # Obtain extreme results
        extreme_results <- obtain_results(extreme_data, condition_structure, condition_df, n)
        extreme_results$SKEW <- "extreme"
        extreme_results$LOWER_CORRECT <- extreme_results$lower_dimensions == lower_factors
        extreme_results$LOWER_MBE <- extreme_results$lower_dimensions - lower_factors
        extreme_results$HIGHER_CORRECT <- extreme_results$higher_dimensions == condition$HighF
        extreme_results$HIGHER_MBE <- extreme_results$higher_dimensions - condition$HighF

        # Remove row names
        row.names(extreme_results) <- NULL

        # Obtain moderate results
        moderate_results <- obtain_results(moderate_data, condition_structure, condition_df, n)
        moderate_results$SKEW <- "moderate"
        moderate_results$LOWER_CORRECT <- moderate_results$lower_dimensions == lower_factors
        moderate_results$LOWER_MBE <- moderate_results$lower_dimensions - lower_factors
        moderate_results$HIGHER_CORRECT <- moderate_results$higher_dimensions == condition$HighF
        moderate_results$HIGHER_MBE <- moderate_results$higher_dimensions - condition$HighF

        # Remove row names
        row.names(moderate_results) <- NULL

        # Combine for final result
        final_result <- rbind.data.frame(extreme_results, moderate_results)
        attr(final_result, "timing") <- rbind(
          attributes(extreme_results)$timing, attributes(moderate_results)$timing
        )

        # Combine results
        return(final_result)

      }, mc.cores = 50, mc.preschedule = FALSE
    )

    # Obtain timing
    timing <- do.call(rbind, lapply(result, function(x){attributes(x)$timing}))

    # Save timing
    save(timing, file = paste0("./timing/N-", n, "_condition-", condition$number, ".RData"))

    # Combine result
    result <- do.call(rbind.data.frame, result)

    # Save result
    save(result, file = paste0("./results/N-", n, "_condition-", condition$number, ".RData"))

    # Remove result
    rm(result)
    gc(verbose = FALSE)
    if((i %% 10) == 0){
      system("sudo sh -c 'echo 3 > /proc/sys/vm/drop_caches'")
    }

    # Print update
    cat(paste0("\r N = ", n, ": Run ", condition$number, " of ", n_conditions, " completed."))

  }

}

# Load results
FILES <- list.files("./results/")
results <- do.call(rbind.data.frame, lapply(FILES, function(file){
  load(paste0("./results/", file))
  return(result)
}))

# Save score results
save(results, file = "./results.RData")

# Load results
FILES <- list.files("./timing/")
timing <- do.call(rbind.data.frame, lapply(FILES, function(file){
  load(paste0("./timing/", file))
  return(timing)
}))

# Save score results
save(timing, file = "./timing.RData")