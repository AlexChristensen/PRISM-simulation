#%%%%%%%%%%%%%%%%%%%%%%#
#### Obtain Results ####
#%%%%%%%%%%%%%%%%%%%%%%#

# Load packages
library(compiler)

# Obtain results
obtain_results <- function(data_condition, condition_structure, condition_df, n)
{

  # Generate data and obtain results
  condition_results <- lapply(seq_along(data_condition), function(k){

    # Set data for the requested sample size
    data <- data_condition[[k]][[as.character(n)]]

    # Add data to structure
    condition_structure[[k]]$data <- data

    # Obtain results
    result_time <- system.time(
      results <- try(EGAnet:::silent_call(get_results(condition_structure[[k]])), silent = TRUE), gcFirst = FALSE
    )

    # Check for error
    if(is(results, "try-error")){

      if(!file.exists("~/Desktop/bad_data.RData")){
        bad_data <- condition_structure[[k]]
        save(bad_data, file = "~/Desktop/bad_data.RData")
      }

    }

    # Set up results data frame
    results_df <- cbind.data.frame(
      condition_df, N = n, Error = "Rpop" %in% names(condition_structure[[k]]),
      Overlap = condition_structure[[k]]$overlap, results
    )

    # Add results time to timing
    timing <- attributes(results)$timing

    # Attach timing
    attr(results_df, "timing") <- c(timing, result_time = result_time[["elapsed"]] - sum(timing))

    # Return results
    return(results_df)

  })

  # Combine for full results
  full_results <- do.call(rbind.data.frame, condition_results)

  # Attach timing
  attr(full_results, "timing") <- do.call(
    rbind, lapply(condition_results, function(x){attributes(x)$timing})
  )

  # Return results
  return(full_results)

}

# Compile function
obtain_results <- cmpfun(obtain_results)