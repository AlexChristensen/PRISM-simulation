#%%%%%%%%%%%%%%%%%%%%%%%%%%#
#### Generate Condition ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Load packages
library(latentFactoR); library(fungible); library(psych); library(compiler)

# Check condition
check_condition <- function(simulated)
{

  # Higher-order implied lower-order factor correlations
  lower_R <- simulated$parameters$higher$loadings %*% tcrossprod(
    simulated$parameters$higher$correlations,
    simulated$parameters$higher$loadings
  )

  # Check higher-order communalities
  check_higher_communalities <- any(diag(lower_R) > 0.90)

  # Set diagonal to zero
  diag(lower_R) <- 1

  # Compute population correlation
  R <- simulated$parameters$lower$loadings %*% tcrossprod(
    lower_R, simulated$parameters$lower$loadings
  )

  # Check communalities
  check_communalities <- check_higher_communalities | any(diag(R) > 0.90)

  # Reset diagonal
  diag(R) <- 1

  # Check eigenvalues (disturbance, factor-level, and variable-level matrices)
  check_eigenvalues <- any(latentFactoR:::matrix_eigenvalues(lower_R) <= 0) |
    any(latentFactoR:::matrix_eigenvalues(R) <= 0)

  # Return output
  return(check_higher_communalities | check_communalities | check_eigenvalues)

}

# Compile function
check_condition <- cmpfun(check_condition)

# CFI
cfi <- function(efa, S, nobs = 8.3e09)
{

  # Get parameters
  p <- nrow(efa$loadings)
  k <- ncol(efa$loadings)

  # Null model (already computed by psych::fa with n.obs = nobs)
  chisq_null <- efa$null.chisq
  df_null <- efa$null.dof
  null_value <- max(chisq_null - df_null, 0)

  # Model-implied correlation matrix
  model_R <- efa$loadings %*% tcrossprod(efa$Phi, efa$loadings)
  diag(model_R) <- 1

  # ML discrepancy function
  ML <- log(det(model_R)) + sum(diag(S %*% solve(model_R))) - log(det(S)) - p

  # Chi-square: scale by N (Bartlett-corrected, matching psych::fa's own
  # convention for null.chisq), NOT by degrees of freedom
  chisq <- ML * ((nobs - 1) - (2 * p + 5) / 6 - (2 * k) / 3)
  df <- efa$dof

  # Return CFI
  return((null_value - max(chisq - df, 0)) / null_value)

}

# Compile function
cfi <- cmpfun(cfi)

# Not recoverable
not_recoverable <- function(simulated)
{

  # Estimate EFA on population correlation matrix
  EFA <- EGAnet:::silent_call(efa(
    data = simulated$population_correlation,
    nfactors = simulated$parameters$lower$factors,
    cor = "pearson"
  ))

  # Set up true loadings
  true_loadings <- simulated$parameters$lower$loadings

  # Align with true loadings
  aligned <- EGAnet:::faAlign_fungible(
    F1 = true_loadings,
    F2 = EFA$loadings,
    Phi2 = EFA$Phi
  )

  # Put back into EFA
  EFA$loadings <- aligned$F2
  EFA$Phi <- aligned$Phi2

  # Set true membership
  true_membership <- abs(true_loadings) >= 0.30

  # Obtain whether it's recoverable
  recoverable <- omega_index(
    obtain_memberships(aligned$F2, true_membership, simulated$overlap),
    true_membership
  ) == 1

  # Return recoverable
  return(c(recoverable = recoverable, cfi = cfi(EFA, simulated$population_correlation)))

}

# Compile function
not_recoverable <- cmpfun(not_recoverable)

# Generate condition
generate_condition <- function(condition)
{

  # Total variables
  total_variables <- condition$LowV * condition$LowF * condition$HighF

  # Set lower factors
  lower_factors <- condition$LowF * condition$HighF

  # Set factors
  n_lower  <- condition$LowF
  n_higher <- condition$HighF

  # Number of overlapping cross-loadings to add
  lower_add <- condition$LowV * 0.20

  # Factor sequences
  lower_sequence <- seq_len(lower_factors)

  # Loop over communalities
  while(TRUE){

    # Generate data
    simulated_base <- simulate_hierarchical_factors(
      lower_factors = lower_factors, variables = condition$LowV,
      lower_loadings = condition$LowL, lower_cross_loadings = 0.00,
      higher_factors = condition$HighF, higher_loadings = condition$HighL,
      higher_cross_loadings = 0.00, higher_correlations = condition$HighC,
      off_disturbances = 0.00, sample_size = 1000,
      variable_categories = 5, skew_range = c(0, 2)
    )

    # Set cross-loadings to zero
    lower_loadings <- simulated_base$parameters$lower$loadings
    lower_loadings[abs(lower_loadings) < 0.40] <- 0
    higher_loadings <- simulated_base$parameters$higher$loadings
    higher_loadings[abs(higher_loadings) < 0.40] <- 0

    # Re-generate for no error and overlap condition
    no_overlap <- simulate_hierarchical_factors(
      lower_factors = simulated_base$parameters$lower$factors,
      variables = simulated_base$parameters$lower$variables,
      lower_loadings = lower_loadings, lower_cross_loadings = 0.00,
      higher_factors = simulated_base$parameters$higher$factors,
      higher_loadings = higher_loadings, higher_cross_loadings = 0.00,
      higher_correlations = simulated_base$parameters$higher$correlations,
      off_disturbances = 0, sample_size = 1000,
      variable_categories = 5, skew = simulated_base$parameters$skew
    )

    # Add overlap
    no_overlap$overlap <- FALSE

    # Check condition
    if(check_condition(no_overlap)){
      next
    }

    # Set loadings
    lower_loadings <- no_overlap$parameters$lower$loadings

    # Add overlapping
    for(k in lower_sequence){

      # Set row index
      row_index <- seq_len(condition$LowV) + condition$LowV * (k - 1)
      row_index <- row_index[order(lower_loadings[row_index, k])[seq_len(lower_add)]]

      # k's own higher-order group, and k's position within that group
      h_source   <- ceiling(k / n_lower)
      pos_source <- ((k - 1) %% n_lower) + 1

      # Cyclic successor position within k's own group only
      pos_target <- (pos_source %% n_lower) + 1
      target <- ((h_source - 1) * n_lower) + pos_target

      # Every item flagged for cross-loading on factor k gets the same
      # single deterministic target
      column_index <- rep(target, lower_add)

      # Loop over loadings to add
      for(l in seq_len(lower_add)){

        # Assign cross-loading magnitude (still uniform on 0.30-0.50)
        lower_loadings[row_index[l], column_index[l]] <-
          EGAnet:::runif_xoshiro(1, min = 0.30, max = 0.50)

      }

    }

    # Re-generate for no error but overlap condition
    overlap <- simulate_hierarchical_factors(
      lower_factors = no_overlap$parameters$lower$factors,
      variables = no_overlap$parameters$lower$variables,
      lower_loadings = lower_loadings, lower_cross_loadings = 0.00,
      higher_factors = no_overlap$parameters$higher$factors,
      higher_loadings = higher_loadings, higher_cross_loadings = 0.00,
      higher_correlations = no_overlap$parameters$higher$correlations,
      off_disturbances = 0, sample_size = 1000,
      variable_categories = 5, skew = no_overlap$parameters$skew
    )

    # Add overlap
    overlap$overlap <- TRUE

    # Check condition
    if(check_condition(overlap)){
      next
    }

    # All clear
    break

  }

  # Construct data frame
  result <- as.data.frame(rbind(
    not_recoverable(no_overlap),
    not_recoverable(overlap)
  ))

  # Add overlap
  result$overlap <- c(FALSE, TRUE)

  # Return
  return(result)

}

# Compile function
generate_condition <- cmpfun(generate_condition)
