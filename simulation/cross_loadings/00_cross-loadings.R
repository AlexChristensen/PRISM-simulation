#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#### Hierarchical | Main Results ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Provide password
system("sudo sh -c 'echo 3 > /proc/sys/vm/drop_caches'")

# Load packages
library(latentFactoR); library(EGAnet); library(parallel); library(tidyverse)

# Conditions
conditions <- expand.grid(
  NVAR = c(4, 8, 12),
  NFAC = c(2, 4, 6),
  LOAD = c(0.50, 0.60, 0.70),
  CROSS = c(0.05, 0.10, 0.15)
)

# Number of conditions
n_conditions <- nrow(conditions)
n_digits <- latentFactoR:::digits(n_conditions) - 1

# Initialize full results
full_results <- vector("list", n_conditions)

# Loop over conditions
for(i in seq_len(n_conditions)){

  # Obtain condition
  condition <- conditions[i,]

  # Set condition number (for labeling/output)
  condition_number <- latentFactoR:::format_integer(i, n_digits)

  # Set membership
  membership <- rep(seq_len(condition$NFAC), each = condition$NVAR)
  simple_structure <- outer(membership, seq_len(condition$NFAC), "==")

  # Obtain results
  results <- mclapply(seq_len(100), function(j){

    # Set error
    simulated <- try(stop("error"), silent = TRUE)

    # Generate data
    while(is(simulated, "try-error")){
      simulated <- try(simulate_factors(
        factors = condition$NFAC, variables = condition$NVAR,
        loadings = condition$LOAD, cross_loadings = condition$CROSS,
        correlations_range = c(0.10, 0.70),
        sample_size = 1000
      ), silent = TRUE)
    }

    # Estimate network
    ega <- EGA(simulated$data, plot.EGA = FALSE)

    # Estimate loadings
    loadings <- EGAnet:::silent_call(net.loads(ega$network, membership)$std[names(ega$wc),])

    # Simple loadings
    simple_loadings <- loadings * simple_structure

    # Compute correlations
    simple_corr <- cov2cor(crossprod(simple_loadings, ega$correlation) %*% simple_loadings)
    full_corr <- cov2cor(crossprod(loadings, ega$correlation) %*% loadings)

    # Compute metrics
    lower_triangle <- lower.tri(simulated$parameters$factor_correlations)
    true_corr <- simulated$parameters$factor_correlations[lower_triangle]

    # Return result
    return(
      data.frame(
        rbind(condition, condition),
        METHOD = c("SIMPLE", "FULL"),
        MAE = c(
          mean(abs(simple_corr[lower_triangle] - true_corr)),
          mean(abs(full_corr[lower_triangle] - true_corr))
        ),
        MBE = c(
          mean(simple_corr[lower_triangle] - true_corr),
          mean(full_corr[lower_triangle] - true_corr)
        )
      )
    )

  }, mc.cores = 10)

  # Combine result
  full_results[[i]] <- do.call(rbind.data.frame, results)

  # Print update
  cat(paste0("\r Condition ", condition_number, " completed."))

}

# Bind all results
final_results <- do.call(rbind.data.frame, full_results)

# Check general results
final_results %>%
  group_by(METHOD, CROSS) %>%
  summarize(
    MAE = mean(MAE, na.rm = TRUE),
    MBE = mean(MBE, na.rm = TRUE)
  ) %>% as.data.frame() %>% print(digits = 3)

# METHOD CROSS    MAE     MBE
# 1   FULL  0.05 0.0857  0.0735
# 2   FULL  0.10 0.1198  0.0900
# 3   FULL  0.15 0.1440  0.0689
# 4 SIMPLE  0.05 0.0853 -0.0781
# 5 SIMPLE  0.10 0.0914 -0.0679
# 6 SIMPLE  0.15 0.1109 -0.0697

# Check results
test <- final_results %>%
  group_by(METHOD, NVAR, NFAC, LOAD, CROSS) %>%
  summarize(
    MAE = mean(MAE, na.rm = TRUE),
    MBE = mean(MBE, na.rm = TRUE)
  )

# Compute t-test
t.test(
  test$MAE[test$METHOD == "SIMPLE"],
  test$MAE[test$METHOD == "FULL"],
  paired = TRUE, var.equal = var.test(
    test$MAE[test$METHOD == "SIMPLE"],
    test$MAE[test$METHOD == "FULL"]
  )$p.value < 0.05
)

# Paired t-test
#
# data:  test$MAE[test$METHOD == "SIMPLE"] and test$MAE[test$METHOD == "FULL"]
# t = -6.9812, df = 80, p-value = 0.0000000007712
# alternative hypothesis: true mean difference is not equal to 0
# 95 percent confidence interval:
#   -0.02645962 -0.01472069
# sample estimates:
#   mean difference
# -0.02059015

# Compute Cohen's d
abs(EGAnet:::d(
  test$MAE[test$METHOD == "SIMPLE"],
  test$MAE[test$METHOD == "FULL"],
  paired = TRUE
))

# 0.775685