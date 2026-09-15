#%%%%%%%%%%%%%%%%%%%%%%%%%%#
#### Generate Condition ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Load packages
library(latentFactoR); library(fungible); library(psych); library(compiler)

# Check condition
check_condition <- function(lower_loadings, higher_loadings, higher_correlations)
{

  # Higher-order implied lower-order factor correlations
  lower_R <- higher_loadings %*% tcrossprod(higher_correlations, higher_loadings)

  # Check higher-order communalities
  check_higher_communalities <- any(diag(lower_R) > 0.90)

  # Set diagonal to zero
  diag(lower_R) <- 1

  # Compute population correlation
  R <- lower_loadings %*% tcrossprod(lower_R, lower_loadings)

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

  # Calculate implied R
  implied_R <- EFA$loadings %*% tcrossprod(EFA$Phi, EFA$loadings)
  diag(implied_R) <- 1

  # Set true membership
  true_membership <- abs(true_loadings) >= 0.40

  # Obtain whether it's recoverable
  recoverable <- omega_index(
    obtain_memberships(aligned$F2, true_membership, simulated$overlap),
    true_membership
  ) == 1

  # Return recoverable
  return(c(recoverable = recoverable, max_abs_res = max(abs(implied_R - simulated$population_correlation))))

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

  # CFI sequence
  cfi_sequence <- seq(0.95, 0.99, 0.01)

  # Loop over communalities
  while(TRUE){

    # Generate data
    simulated_base <- simulate_hierarchical_factors(
      lower_factors = lower_factors, variables = condition$LowV,
      lower_loadings = condition$LowL, lower_cross_loadings = 0.00,
      higher_factors = condition$HighF, higher_loadings = condition$HighL,
      higher_cross_loadings = 0.00, higher_correlations = condition$HighC,
      off_disturbances = 0.00, sample_size = 100,
      variable_categories = 5, skew_range = c(0, 2)
    )

    # Set cross-loadings to zero
    lower_loadings <- simulated_base$parameters$lower$loadings
    lower_loadings[abs(lower_loadings) < 0.40] <- 0
    higher_loadings <- simulated_base$parameters$higher$loadings
    higher_loadings[abs(higher_loadings) < 0.40] <- 0

    # Check condition
    if(check_condition(lower_loadings, higher_loadings, simulated_base$parameters$higher$correlations)){
      next
    }

    # Re-generate for no error and overlap condition
    no_overlap <- simulate_hierarchical_factors(
      lower_factors = simulated_base$parameters$lower$factors,
      variables = simulated_base$parameters$lower$variables,
      lower_loadings = lower_loadings, lower_cross_loadings = 0.00,
      higher_factors = simulated_base$parameters$higher$factors,
      higher_loadings = higher_loadings, higher_cross_loadings = 0.00,
      higher_correlations = simulated_base$parameters$higher$correlations,
      off_disturbances = 0, sample_size = 100,
      variable_categories = 5, skew = simulated_base$parameters$skew
    )

    # Add overlap
    no_overlap$overlap <- FALSE

    # Generate with population error
    no_overlap_error <- simFA(
      Model = list(
        NFac = no_overlap$parameters$lower$factors,
        NItemPerFac = no_overlap$parameters$lower$variables,
        Model = "oblique"
      ),
      Loadings = list(FacPattern = no_overlap$parameters$lower$loadings),
      Phi = list(
        PhiType = "user",
        UserPhi = no_overlap$parameters$lower$correlations
      )
    )

    # Generate error objects
    no_overlap_errors <- lapply(cfi_sequence, function(target_cfi){

      # Generate error
      error <- EGAnet:::silent_call(noisemaker(
        mod = no_overlap_error, method = "TKL",
        target_cfi = target_cfi,
        tkl_ctrl = list(optim_type = "optim", factr = 1e10)
      ))

      # Set parameters
      no_overlap_error$population_correlation <- error$Sigma
      no_overlap_error$CFI <- error$cfi
      no_overlap_error$parameters <- no_overlap$parameters
      no_overlap_error$overlap <- FALSE

      # Return object
      return(no_overlap_error)

    })

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

        # Assign cross-loading magnitude
        lower_loadings[row_index[l], c(k, column_index[l])] <- sort(
          EGAnet:::runif_xoshiro(2, min = 0.40, max = 0.50),
          decreasing = TRUE
        )

      }

    }

    # Check condition
    if(check_condition(lower_loadings, higher_loadings, no_overlap$parameters$higher$correlations)){
      next
    }

    # Re-generate for no error but overlap condition
    overlap <- simulate_hierarchical_factors(
      lower_factors = no_overlap$parameters$lower$factors,
      variables = no_overlap$parameters$lower$variables,
      lower_loadings = lower_loadings, lower_cross_loadings = 0.00,
      higher_factors = no_overlap$parameters$higher$factors,
      higher_loadings = higher_loadings, higher_cross_loadings = 0.00,
      higher_correlations = no_overlap$parameters$higher$correlations,
      off_disturbances = 0, sample_size = 100,
      variable_categories = 5, skew = no_overlap$parameters$skew
    )

    # Generate with population error
    overlap_error <- simFA(
      Model = list(
        NFac = overlap$parameters$lower$factors,
        NItemPerFac = overlap$parameters$lower$variables,
        Model = "oblique"
      ),
      Loadings = list(FacPattern = overlap$parameters$lower$loadings),
      Phi = list(
        PhiType = "user",
        UserPhi = overlap$parameters$lower$correlations
      )
    )

    # Generate error objects
    overlap_errors <- lapply(cfi_sequence, function(target_cfi){

      # Generate error
      error <- EGAnet:::silent_call(noisemaker(
        mod = overlap_error, method = "TKL",
        target_cfi = target_cfi,
        tkl_ctrl = list(optim_type = "optim", factr = 1e10)
      ))

      # Set parameters
      overlap_error$population_correlation <- error$Sigma
      overlap_error$CFI <- error$cfi
      overlap_error$parameters <- overlap$parameters
      overlap_error$overlap <- TRUE

      # Return object
      return(overlap_error)

    })

    # All clear
    break

  }

  # Construct data frame
  no_errors <- as.data.frame(do.call(rbind, lapply(no_overlap_errors, not_recoverable)))
  no_errors$target_cfi <- cfi_sequence
  no_errors$actual_cfi <- EGAnet:::nvapply(no_overlap_errors, function(x){x$CFI})
  no_errors$overlap <- FALSE
  o_errors <- as.data.frame(do.call(rbind, lapply(overlap_errors, not_recoverable)))
  o_errors$target_cfi <- cfi_sequence
  o_errors$actual_cfi <- EGAnet:::nvapply(overlap_errors, function(x){x$CFI})
  o_errors$overlap <- TRUE

  # Return data
  return(rbind.data.frame(no_errors, o_errors))

}

# Compile function
generate_condition <- cmpfun(generate_condition)
