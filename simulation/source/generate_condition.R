#%%%%%%%%%%%%%%%%%%%%%%%%%%#
#### Generate Condition ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Load packages
library(latentFactoR); library(fungible); library(bifactor); library(compiler)

# Check condition
check_condition <- function(higher_loadings, higher_R, lower_loadings)
{

  # Higher-order implied lower-order factor correlations
  lower_R <- higher_loadings %*% tcrossprod(higher_R, higher_loadings)

  # Check higher-order
  if(any(diag(lower_R) > 0.90)){

    # Commonalities
    res <- TRUE
    attr(res, "cause") <- "Higher commonalities > 0.90"
    return(res)

  }

  # Set diagonal to zero
  diag(lower_R) <- 1

  # Check positive definite
  if(any(latentFactoR:::matrix_eigenvalues(lower_R) <= 0)){

    # Not positive definite
    res <- TRUE
    attr(res, "cause") <- "Higher not positive definite"
    return(res)

  }

  # Compute population correlation
  R <- lower_loadings %*% tcrossprod(lower_R, lower_loadings)

  # Check overall
  if(any(diag(R) > 0.90)){

    # Commonalities
    res <- TRUE
    attr(res, "cause") <- "Overall commonalities > 0.90"
    return(res)

  }

  # Reset diagonal
  diag(R) <- 1

  # Check positive definite
  if(any(latentFactoR:::matrix_eigenvalues(R) <= 0)){

    # Not positive definite
    res <- TRUE
    attr(res, "cause") <- "Overall not positive definite"
    return(res)

  }

  # Return good check
  return(FALSE)

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
  aligned <- faAlign(
    F1 = true_loadings,
    F2 = EFA$loadings,
    Phi2 = EFA$Phi
  )

  # Put back into EFA
  EFA$loadings <- aligned$F2
  EFA$Phi <- aligned$Phi2

  # Set true membership
  true_membership <- abs(true_loadings) >= 0.40

  # Check for residual
  if(max(abs(EFA$implied_R - simulated$population_correlation)) > 0.10){

    # Residual
    res <- TRUE
    attr(res, "cause") <- "Maximum residual greater than 0.10"
    return(res)

  }

  # Obtain whether it's recoverable
  recoverable <- omega_index(
    obtain_memberships(aligned$F2, true_membership, simulated$overlap),
    true_membership
  ) == 1

  # Check for recoverable
  if(!recoverable){

    # Residual
    res <- TRUE
    attr(res, "cause") <- "Omega index did not equal 1"
    return(res)

  }

  # Compute CFI
  CFI <- cfi(EFA, simulated$population_correlation)

  # Check for CFI
  if(abs(CFI - 0.97) > 1){

    # Residual
    res <- TRUE
    attr(res, "cause") <- "CFI outside of tolerance"
    return(res)

  }

  # Return good check
  return(FALSE)

}

# Compile function
not_recoverable <- cmpfun(not_recoverable)

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

  # Initialize count
  count <- 0

  # Initialize bad generation
  bad_gen <- character()
  factors <- c(
    "Higher commonalities > 0.90", "Higher not positive definite",
    "Overall commonalities > 0.90", "Overall not positive definite",
    "Maximum residual greater than 0.10",
    "Omega index did not equal 1", "CFI outside of tolerance"
  )

  # Loop over communalities
  while(TRUE){

    # Update count
    count <- count + 1

    # Generate data
    simulated_base <- simulate_hierarchical_factors(
      lower_factors = lower_factors, variables = condition$LowV,
      lower_loadings = condition$LowL, lower_cross_loadings = 0.00,
      higher_factors = condition$HighF, higher_loadings = condition$HighL,
      higher_cross_loadings = 0.00, higher_correlations = condition$HighC,
      off_disturbances = 0.00, sample_size = 100
      # sample size, variable categories, and skew are added in scripts:
      # 03_scores_results.R
      # 04_main_results.R
    )

    # No error ----

    # Set cross-loadings to zero
    lower_loadings <- simulated_base$parameters$lower$loadings
    lower_loadings[abs(lower_loadings) < 0.40] <- 0
    higher_loadings <- simulated_base$parameters$higher$loadings
    higher_loadings[abs(higher_loadings) < 0.40] <- 0
    higher_R <- simulated_base$parameters$higher$correlations

    # Check condition
    condition_check <- check_condition(higher_loadings, higher_R, lower_loadings)
    if(condition_check){
      bad_gen <- c(bad_gen, attr(condition_check, "cause"))
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
      off_disturbances = 0, sample_size = 100
    )

    # Add overlap
    no_overlap$overlap <- FALSE

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
    condition_check <- check_condition(higher_loadings, higher_R, lower_loadings)
    if(condition_check){
      bad_gen <- c(bad_gen, attr(condition_check, "cause"))
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
      off_disturbances = 0, sample_size = 100
    )

    # Add overlap
    overlap$overlap <- TRUE

    # Error ----

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

    # Generate correlation matrix
    no_overlap_error$population_correlation <- noisemaker(
      mod = no_overlap_error, method = "TKL",
      target_cfi = 0.97,
      tkl_ctrl = list(optim_type = "optim", factr = 1e10)
      # target_cfi determine from 00_simulate_error.R
    )$Sigma

    # Add objects to mirror necessary output
    no_overlap_error$overlap <- FALSE
    no_overlap_error$parameters <- no_overlap$parameters

    # Check recoverability
    recoverable_check <- not_recoverable(no_overlap_error)
    if(recoverable_check){
      bad_gen <- c(bad_gen, attr(recoverable_check, "cause"))
      next
    }

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

    # Generate correlation matrix
    overlap_error$population_correlation <- noisemaker(
      mod = overlap_error, method = "TKL",
      target_cfi = 0.97,
      tkl_ctrl = list(optim_type = "optim", factr = 1e10)
      # target_cfi determine from 00_simulate_error.R
    )$Sigma

    # Add objects to mirror necessary output
    overlap_error$overlap <- TRUE
    overlap_error$parameters <- overlap$parameters

    # Check recoverability
    recoverable_check <- not_recoverable(overlap_error)
    if(recoverable_check){
      bad_gen <- c(bad_gen, attr(recoverable_check, "cause"))
      next
    }

    # All clear
    break

  }

  # Set data to NULL for more compact storage
  overlap$data <- no_overlap$data <- NULL

  # Return data
  return(
    list(
      no_overlap = no_overlap, no_overlap_error = no_overlap_error,
      overlap = overlap, overlap_error = overlap_error,
      bad_generation = c(iterations = count, table(factor(bad_gen, levels = factors)))
    )
  )

}

# Compile function
generate_condition <- cmpfun(generate_condition)