#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#### Dimensionality Results ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Metrics
# Dimensions -- accuracy, mean bias error, ARI (with overlap)
# Loading recovery -- congruence (cosine)
# Factor correlations -- mean absolute error

# Load packages
library(igraph); library(compiler)

# Obtain memberships
obtain_memberships <- function(loadings, true_membership, overlap)
{

  # Obtain memberships
  factor_sequence <- seq_len(ncol(loadings))
  primary <- outer(max.col(abs(loadings)), factor_sequence, "==")

  # Check for overlap
  if(overlap){
    secondary <- outer(max.col(abs(loadings * (!primary))), factor_sequence, "==")
    secondary[rowSums(true_membership) == 1,] <- FALSE
    return(primary | secondary)
  }else{
    return(primary)
  }

}

# Compile
obtain_memberships <- cmpfun(obtain_memberships)

# Obtain omega index
# Collins, L. M., & Dent, C. W. (1988).
# Omega: A general formulation  of the rand index of cluster recovery suitable for non-disjoint solutions.
# Multivariate Behavioral Research, 23(2), 231–242.
omega_index <- function(membership, true_membership)
{

  # Check for NA
  if(all(is.na(membership))){
    return(NA)
  }

  # Co-membership counts
  comembership <- tcrossprod(membership)
  true_comembership <- tcrossprod(true_membership)

  # Obtain lower triangle
  lower_triangle <- lower.tri(comembership)
  lower_comembership <- comembership[lower_triangle]
  true_lower_comembership <- true_comembership[lower_triangle]

  # Guard against single dimensions
  if(length(unique(true_lower_comembership)) == 1){
    return(NA)
  }

  # Total comparisons
  total <- length(lower_comembership)

  # Maximum overlapping communities (used for bins below)
  max_overlap <- max(lower_comembership, true_lower_comembership, na.rm = TRUE) + 1

  # omega_u: fraction of pairs whose co-occurrence count matches exactly
  omega_u <- mean(lower_comembership == true_lower_comembership, na.rm = TRUE)

  # omega_e: expected match fraction under independent marginal distributions
  freq <- tabulate(lower_comembership + 1, nbins = max_overlap) / total
  true_freq <- tabulate(true_lower_comembership + 1, nbins = max_overlap) / total
  omega_e <- sum(freq * true_freq)

  # Return omega index
  return((omega_u - omega_e) / (1 - omega_e))

}

# Compile
omega_index <- cmpfun(omega_index)

# Return dimensions
return_dimensions <- function(dimensions)
{

  ifelse(
    all(is.na(dimensions)), 0,
    EGAnet:::unique_length(dimensions)
  )

}

# Compile
return_dimensions <- cmpfun(return_dimensions)

# Obtain higher memberships
higher_memberships <- function(lower_memberships, higher_loadings, higher_true)
{

  # Check for missing (e.g., zero dimensions)
  if(all(is.na(higher_loadings))){
    return(NA)
  }else{
    return(obtain_memberships(lower_memberships %*% higher_loadings, higher_true, FALSE))
  }

}

# Compile
higher_memberships <- cmpfun(higher_memberships)

# Get single result
get_single_result <- function(method, methods, lower_true, higher_true, simulated)
{

  # Check for loadings
  if(all(is.na(methods$loadings$lower[[method]]))){

    # Return dimensions
    return(c(
        lower_dimensions = methods$dimensions$lower[[method]],
        lower_ari = NA,
        higher_dimensions = NA,
        higher_ari = NA
    ))

  }

  # Check results
  if(
    (!is.na(methods$dimensions$lower[[method]])) & # not missing method
    (methods$dimensions$lower[[method]] > 0) # at least one dimension
  ){

    # Set memberships
    lower_memberships <- obtain_memberships(
      methods$loadings$lower[[method]], lower_true, simulated$overlap
    )

    # Obtain ARI
    lower_dimensions <- return_dimensions(max.col(lower_memberships))
    lower_ari <- omega_index(lower_memberships, lower_true)

    # Check higher order dimensions
    if(
      (!is.na(methods$dimensions$higher[[method]])) & # not missing method
      (methods$dimensions$higher[[method]] > 0) # at least one dimension
    ){

      # Set memberships
      higher_membership <- higher_memberships(
        lower_memberships, methods$loadings$higher[[method]], higher_true
      )


      # Check for missing
      if(all(is.na(higher_membership))){

        # Return missing
        higher_dimensions <- NA
        higher_ari <- NA

      }else{

        # Obtain ARI
        higher_dimensions <- return_dimensions(max.col(higher_membership))
        higher_ari <- omega_index(higher_membership, higher_true)

      }

    }else{

      # Set results
      higher_dimensions <- methods$dimensions$higher[[method]]
      higher_ari <- NA

    }

  }else{

    # Set results
    lower_dimensions <- methods$dimensions$lower[[method]]
    lower_ari <- NA
    higher_dimensions <- methods$dimensions$higher[[method]]
    higher_ari <- NA

  }

  # Return results
  return(c(
    lower_dimensions = lower_dimensions,
    lower_ari = lower_ari,
    higher_dimensions = higher_dimensions,
    higher_ari = higher_ari
  ))

}

# Compile function
get_single_result <- cmpfun(get_single_result)

# Obtain results
get_results <- function(simulated)
{

  # Obtain methods
  methods <- obtain_methods(simulated)

  # Set true loadings
  lower_true <- abs(simulated$parameters$lower$loadings) >= 0.30
  higher_true <- abs(lower_true %*% simulated$parameters$higher$loadings) >= 0.30

  # Set methods
  METHODS <- names(methods$loadings$lower)

  # Loop over results
  results <- as.data.frame(do.call(
    rbind, lapply(
      METHODS, get_single_result,
      methods = methods, lower_true = lower_true,
      higher_true = higher_true, simulated = simulated
    )
  ))

  ## Add method
  results$METHOD <- METHODS
  attr(results, "timing") <- methods$timing

  # Return result
  return(results)

}

# Compile
get_results <- cmpfun(get_results)
