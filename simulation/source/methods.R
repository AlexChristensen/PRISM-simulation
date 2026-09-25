#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#### Dimensionality Methods ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Load packages
library(EGAnet); library(latentFactoR); library(psych)
library(GPArotation); library(compiler); library(Matrix)

# Set polychoric
polychoric <- function(data)
{

  # Obtain polychoric matrix
  R <- polychoric.matrix(data)

  # Check for smoothing
  if(any(eigen(R, symmetric = TRUE, only.values = TRUE)$values < 0)){
    R <- as.matrix(nearPD(R, corr = TRUE, keepDiag = TRUE, ensureSymmetry = TRUE)$mat)
  }

  # Return correlation matrix
  return(R)

}

# Compile function
polychoric <- cmpfun(polychoric)

# SMC
SMC <- function(R)
{

  # Obtain SmC
  smc <- 1 - 1 / diag(solve(R))

  # Set bound indices
  greater <- smc > 1
  less <- smc < 0

  # Check bounds
  if(any(greater)){
    smc[greater] <- 1
  }
  if(any(less)){
    smc[less] <- 0
  }

  # Return SMC
  return(smc)

}

# Compile function
SMC <- cmpfun(SMC)

# Parallel Analysis
parallel_analysis <- function(data, corr, n.iter = 20, quant = 0.95)
{

  # Set up correlation function
  corr_FUN <- switch(
    corr,
    "pearson" = cor,
    "polychoric" = polychoric
  )

  # Sub-routine for correct counting of dimensions
  count_dimensions <- function(greater_than, variables)
  {

    # Initialize count
    count <- 0

    # Loop over to count consecutive dimensions
    for(i in seq_len(variables)){

      # Check for element TRUE
      if(greater_than[i]){
        count <- count + 1
      }else{
        break
      }

    }

    # Return dimensions
    return(count)

  }

  # Set dimensions
  dimensions <- dim(data)
  variables <- dimensions[2]

  # Compute polychoric correlations
  R_empirical <- corr_FUN(data)

  # Obtain PCA empirical eigenvalues
  PCA_empirical <- eigen(R_empirical, symmetric = TRUE, only.values = TRUE)$values

  # Compute SMC
  diag(R_empirical) <- psych::smc(R_empirical)

  # Obtain PAF empirical eigenvalues
  PAF_empirical <- eigen(R_empirical, symmetric = TRUE, only.values = TRUE)$values

  # Collect empirical resamples
  resampled <- lapply(seq_len(n.iter), function(i){

    # Shuffle data
    shuffled <- apply(data, 2, EGAnet:::shuffle)

    # Compute polychoric correlations
    R_shuffled <- corr_FUN(shuffled)

    # Obtain PCA shuffled eigenvalues
    PCA_shuffled <- eigen(R_shuffled, symmetric = TRUE, only.values = TRUE)$values

    # Compute SMC
    diag(R_shuffled) <- psych::smc(R_shuffled)

    # Return eigenvalues
    return(list(
      pca = PCA_shuffled,
      paf = eigen(R_shuffled, symmetric = TRUE, only.values = TRUE)$values
    ))

  })

  # Obtain means
  PCA_mean <- rowMeans(do.call(cbind, lapply(resampled, function(x){x$pca})), na.rm = TRUE)
  PAF_mean <- rowMeans(do.call(cbind, lapply(resampled, function(x){x$paf})), na.rm = TRUE)

  # # Obtain quantiles
  # PCA_quantile <- apply(
  #   do.call(cbind, lapply(resampled, function(x){x$pca})), 1, quantile, probs = quant
  # )
  # PAF_quantile <- apply(
  #   do.call(cbind, lapply(resampled, function(x){x$paf})), 1, quantile, probs = quant
  # )

  # Return results
  return(c(
    pca = count_dimensions(PCA_empirical >= PCA_mean, variables),
    paf = count_dimensions(PAF_empirical >= PAF_mean, variables)
  ))

}

# Compile
parallel_analysis <- cmpfun(parallel_analysis)

# MAP and VSS
#
# Custom replacement for `psych::vss(fm = "pc")`. psych re-derives a full
# eigendecomposition of R (via `principal()`) once per candidate k, and
# separately re-derives another eigendecomposition for MAP -- (max_k + 1)
# decompositions of the same matrix per call. It also computes the "very
# simple structure" residual fit for every complexity level c = 1..k at
# every k, even though only c = 1 (VSS1) and c = 2 (VSS2) are ever used
# downstream, which is unnecessary O(max_k^2) work.
#
# This version computes R's eigendecomposition once and reuses it for both
# MAP and the VSS1/VSS2 sweep, and only ever builds the c = 1 and c = 2
# simplified-loading fits. Loadings are varimax-rotated twice in a row to
# match psych's actual pipeline: `vss()` calls `principal()` (which
# varimax-rotates internally, since that is its default) and then
# varimax-rotates the already-rotated result again itself. A single
# rotation reproduces psych's picks in 11/12 cross-checked real-data
# cases; the double rotation reproduces psych's cfit values to machine
# precision (~2e-16) and its picks exactly in all cases checked.
map_vss <- function(R, N, max_k)
{

  # Number of variables
  nvar <- ncol(R)

  # Single eigendecomposition of R, reused for MAP and VSS
  eigen_R <- eigen(R, symmetric = TRUE)
  components <- eigen_R$vectors %*% diag(sqrt(pmax(eigen_R$values, 0)))

  # Velicer's MAP
  map_values <- numeric(max_k)
  for(k in seq_len(max_k)){

    # Partial out first k principal components
    partialed <- R - tcrossprod(components[, seq_len(k), drop = FALSE])
    scaled <- diag(1 / sqrt(diag(partialed)))
    partial_R <- scaled %*% partialed %*% scaled
    diag(partial_R) <- 0

    # Average squared partial correlation
    map_values[k] <- sum(partial_R * partial_R) / (nvar * (nvar - 1))

  }

  # Very Simple Structure (complexity 1 and 2 only)
  total_original <- sum(R * R)
  cfit1 <- numeric(max_k)
  cfit2 <- numeric(max_k)

  for(k in seq_len(max_k)){

    # Loadings for k components, rotated twice when k > 1 to match
    # psych's principal() + vss() double-rotation pipeline
    loadings <- components[, seq_len(k), drop = FALSE]
    if(k > 1){
      loadings <- unclass(varimax(loadings)$loadings)
      loadings <- unclass(varimax(loadings)$loadings)
    }

    # Complexity 1: keep only the largest |loading| per variable
    top1 <- max.col(abs(loadings), ties.method = "first")
    simple1 <- matrix(0, nvar, k)
    simple1[cbind(seq_len(nvar), top1)] <- loadings[cbind(seq_len(nvar), top1)]
    resid1 <- R - tcrossprod(simple1)
    cfit1[k] <- 1 - sum(resid1 * resid1) / total_original

    # Complexity 2: keep only the two largest |loadings| per variable
    if(k >= 2){
      top2 <- t(apply(abs(loadings), 1, order, decreasing = TRUE))[, 1:2, drop = FALSE]
      simple2 <- matrix(0, nvar, k)
      index2 <- cbind(rep(seq_len(nvar), 2), as.vector(top2))
      simple2[index2] <- loadings[index2]
      resid2 <- R - tcrossprod(simple2)
      cfit2[k] <- 1 - sum(resid2 * resid2) / total_original
    }else{
      cfit2[k] <- cfit1[k]
    }

  }

  # Return dimensions
  return(
    c(
      vss1 = which.max(cfit1),
      vss2 = which.max(cfit2),
      map = which.min(map_values)
    )
  )

}

# Compile
map_vss <- cmpfun(map_vss)

# EFA
efa <- function(data, nfactors, corr, n.obs = NULL, randomStarts = 100L, ...)
{

  # Data dimensions
  dimensions <- dim(data)

  # Check for correlations
  if(dimensions[1] != dimensions[2]){
    data <- switch(
      corr,
      "pearson" = cor(data),
      "polychoric" = polychoric(data)
    )
    n.obs <- dimensions[1]
  }else if(is.null(n.obs)){
    n.obs <- 8.3e09
  }

  # Estimate EFA on population correlation matrix
  EFA <- fa(
    r = data, n.obs = n.obs, nfactors = nfactors,
    rotate = "none", fm = EGAnet:::swiftelse(corr == "pearson", "ml", "uls")
  )

  # Check for whether rotation is necessary
  if(nfactors != 1){

    # Rotate loadings
    rotated <- GPArotation::geominQ(
      EFA$loadings, randomStarts = randomStarts, eps = switch(
        as.character(nfactors),
        "2" = 0.0001,
        "3" = 0.001,
        0.01
      )
    )

    # Update output
    EFA$loadings <- rotated$loadings
    EFA$Phi <- rotated$Phi

    # Obtain model-implied correlations
    EFA$implied_R <- EFA$loadings %*% tcrossprod(EFA$Phi, EFA$loadings)
    diag(EFA$implied_R) <- 1

  }else{

    # Obtain model-implied correlations
    EFA$Phi <- matrix(1)
    EFA$implied_R <- tcrossprod(EFA$loadings)
    diag(EFA$implied_R) <- 1

  }

  # Return
  return(EFA)

}

# Compile
efa <- cmpfun(efa)

# Custom bass-ackwards maximum loading
max_loading <- function(R, N, min_k, max_k)
{

  # Initialize count
  count <- min_k - 1

  # Cache of the previous iteration's fit, so the degenerate-solution
  # branch below can reuse it instead of re-fitting an identical model
  previous_EFA <- NULL
  previous_count <- NA

  # Loop over until criterion is met
  while(TRUE){

    # Increase count
    count <- count + 1

    # Estimate EFA on correlation matrix
    EFA <- fa(
      r = R, n.obs = N, nfactors = count,
      rotate = EGAnet:::swiftelse(count == 1, "none", "varimax"),
      fm = "minres"
    )

    # Check loadings
    estimated <- EGAnet:::unique_length(
      apply(abs(EFA$loadings), 1, which.max)
    )

    # Keep valid EFA
    if(estimated == count){
      valid_EFA <- EFA
    }

    # Check for break condition
    if(estimated < count){

      # The degenerate count-factor solution collapses to `estimated`
      # factors, which matches the previous iteration's fit whenever
      # estimated == count - 1 -- reuse it instead of refitting
      if(!is.null(previous_EFA) && estimated == previous_count){

        valid_EFA <- previous_EFA

      }else{

        # Obtain solution where valida
        valid_EFA <- fa(
          r = R, n.obs = N, nfactors = estimated,
          rotate = EGAnet:::swiftelse(count == 1, "none", "varimax"),
          fm = "minres"
        )

      }

      # Break
      break

    }

    # Cache this fit for potential reuse next iteration
    previous_EFA <- EFA
    previous_count <- count

    # Don't let it run
    if(count == max_k){
      break
    }

  }

  # Return EFA
  return(valid_EFA)

}

# Compile
max_loading <- cmpfun(max_loading)

# Custom EFA
efa_custom <- function(data, lower_order, method = c("paf", "pca"))
{

  # Check for zero lower-order dimensions
  if(lower_order %in% c(0, dim(data)[2])){

    # Send back missing
    lower_order_efa <- list(loadings = NA)
    higher_order <- 0
    higher_order_efa <- list(loadings = NA)

  }else{

    # Estimate EFA
    lower_order_efa <- efa(data = data, nfactors = lower_order, corr = "polychoric")

    # Compute scores
    lower_scores <- scale(data) %*% solve(
      lower_order_efa$implied_R, lower_order_efa$loadings %*% lower_order_efa$Phi
    )

    # Check for single dimension
    if(lower_order == 1){

      # Send back missing
      higher_order <- 0
      higher_order_efa <- list(loadings = NA)

    }else{

      # Obtain higher order
      higher_order <- parallel_analysis(lower_scores, corr = "pearson")[[method]]

      # Estimate EFA
      if((higher_order == 0) || (dim(lower_scores)[2] == higher_order)){

        # Set back missing
        higher_order_efa <- list(loadings = NA)

      }else{

        # Estimate EFA
        higher_order_efa <- try(
          efa(data = lower_scores, nfactors = higher_order, corr = "pearson"),
          silent = TRUE
        )

        # Check for optimization issue
        if(is(higher_order_efa, "try-error")){
          higher_order_efa <- list(loadings = NA)
        }

      }

    }

  }

  # Collect output
  return(
    list(
      dimensions = list(
        lower = lower_order,
        higher = higher_order
      ),
      loadings = list(
        lower = lower_order_efa$loadings,
        higher = higher_order_efa$loadings
      )
    )
  )

}

# Compile
efa_custom <- cmpfun(efa_custom)

# Compute network loadings
network_loadings <- function(network, wc, variables)
{

  # Obtain network loadings
  loadings <- EGAnet:::silent_call(net.loads(A = network, wc = wc)$std[variables,, drop = FALSE])

  # Return loadings
  return(loadings[,which(dimnames(loadings)[[2]] != "NA"), drop = FALSE])


}

# Compile
network_loadings <- cmpfun(network_loadings)

# Compute formative scores
simple_loadings <- function(loadings, memberships, dimensions)
{

  # Obtain simple structure
  simple_structure <- outer(memberships, seq_len(dimensions), "==")

  # Check for missing in simple structure
  simple_NA <- is.na(simple_structure)
  if(any(simple_NA)){simple_structure[simple_NA] <- FALSE}

  # Obtain simple loadings
  simple_loadings <- loadings * simple_structure
  simple_loadings[is.na(simple_loadings)] <- 0

  # Send simple loadings
  return(simple_loadings)

}

# Compile
simple_loadings <- cmpfun(simple_loadings)

# Compute formative scores
composite_correlations <- function(R, loadings)
{
  return(cov2cor(crossprod(loadings, R) %*% loadings))
}

# Compile
composite_correlations <- cmpfun(composite_correlations)

# Estimate lower order Louvain
lower_order_louvain <- function(network, R, N, variables)
{

  ## Lower order Louvain
  ega_lower_louvain <- community.consensus(network, order = "lower")
  ega_lower_louvain_dimensions <- EGAnet:::unique_length(ega_lower_louvain)


  # Check for no dimension
  if(ega_lower_louvain_dimensions == 0){

    # Return empty results
    ega_lower_louvain_loadings <- NA
    ega_lower_louvain_higher <- NA

  }else{

    # Estimate loadings
    ega_lower_louvain_loadings <- network_loadings(network, ega_lower_louvain, variables)

    # Check for single dimension
    if(ega_lower_louvain_dimensions == 1){

      # Send back missing
      ega_lower_louvain_higher <- NA

    }else{

      # Compute for higher order
      ega_lower_louvain_simple <- simple_loadings(
        ega_lower_louvain_loadings, ega_lower_louvain, ega_lower_louvain_dimensions
      )
      ega_lower_louvain_R <- composite_correlations(R, ega_lower_louvain_simple)
      ega_lower_louvain_ega <- EGA(ega_lower_louvain_R, n = N, plot.EGA = FALSE)

      ## Check for unidimensional
      if(ega_lower_louvain_ega$n.dim == 0){

        # Send back missing
        ega_lower_louvain_higher <- NA

      }else if(ega_lower_louvain_ega$n.dim == 1){

        ega_lower_louvain_higher <- network_loadings(
          ega_lower_louvain_ega$network,
          ega_lower_louvain_ega$wc,
          colnames(ega_lower_louvain_ega$network)
        )

      }else{

        ega_lower_louvain_higher <- network_loadings(
          ega_lower_louvain_ega$network,
          community.consensus(ega_lower_louvain_ega$network, order = "higher"),
          colnames(ega_lower_louvain_ega$network)
        )

      }

    }

  }

  # Return results
  return(
    list(
      loadings = list(
        lower = ega_lower_louvain_loadings,
        higher = ega_lower_louvain_higher
      ),
      dimensions = list(
        lower = ega_lower_louvain_dimensions,
        higher = EGAnet:::swiftelse(
          all(is.na(ega_lower_louvain_higher)), 0, dim(ega_lower_louvain_higher)[2]
        )
      )
    )
  )

}

# Compile
lower_order_louvain <- cmpfun(lower_order_louvain)

# Estimate PRISM
prism <- function(network, R, N, variables)
{

  ## PRISM
  ega_prism <- community.prism(network)
  ega_prism_dimensions <- EGAnet:::unique_length(ega_prism)


  # Check for no dimension
  if(ega_prism_dimensions == 0){

    # Return empty results
    ega_prism_loadings <- NA
    ega_prism_higher <- NA

  }else{

    # Estimate loadings
    ega_prism_loadings <- network_loadings(network, ega_prism, variables)

    # Check for single dimension
    if(ega_prism_dimensions == 1){

      # Send back missing
      ega_prism_higher <- NA

    }else{

      # Compute for higher order
      ega_prism_simple <- simple_loadings(
        ega_prism_loadings, ega_prism, ega_prism_dimensions
      )
      ega_prism_R <- composite_correlations(R, ega_prism_simple)
      ega_prism_ega <- EGA(ega_prism_R, n = N, plot.EGA = FALSE)

      ## Check for unidimensional
      if(ega_prism_ega$n.dim == 0){

        # Send back missing
        ega_prism_higher <- NA

      }else if(ega_prism_ega$n.dim == 1){

        ega_prism_higher <- network_loadings(
          ega_prism_ega$network,
          ega_prism_ega$wc,
          colnames(ega_prism_ega$network)
        )

      }else{

        ega_prism_higher <- network_loadings(
          ega_prism_ega$network,
          community.prism(ega_prism_ega$network),
          colnames(ega_prism_ega$network)
        )

      }

    }

  }

  # Return results
  return(
    list(
      loadings = list(
        lower = ega_prism_loadings,
        higher = ega_prism_higher
      ),
      dimensions = list(
        lower = ega_prism_dimensions,
        higher = EGAnet:::swiftelse(
          all(is.na(ega_prism_higher)), 0, dim(ega_prism_higher)[2]
        )
      )
    )
  )

}

# Compile
prism <- cmpfun(prism)

# Obtain methods
obtain_methods <- function(simulated)
{

  # Get variable names
  variables <- colnames(simulated$data)

  # Get sample size
  N <- dim(simulated$data)[1]

  # Create timing list
  time_list <- list()

  # Estimate correlations
  time_list$R_time <- system.time(
    R <- polychoric(simulated$data), gcFirst = FALSE
  )

  # Hierarchical EGA
  time_list$network_time <- system.time(
    network <- EGA(R, n = N, plot.EGA = FALSE)$network,
    gcFirst = FALSE
  )

  ## Lower order Louvain
  time_list$lower_louvain_time <- system.time(
    lower_louvain <- lower_order_louvain(network, R, N, variables),
    gcFirst = FALSE
  )

  ## PRISM
  time_list$prism_time <- system.time(
    prism_result <- prism(network, R, N, variables),
    gcFirst = FALSE
  )

  # Estimate lower order for parallel analysis
  time_list$pa_time <- system.time(
    pa_lower <- parallel_analysis(simulated$data, corr = "polychoric"),
    gcFirst = FALSE
  )

  # Estimate higher order for parallel analysis with PCA
  time_list$pca_time <- system.time(
    pca_higher <- efa_custom(simulated$data, pa_lower[["pca"]], "pca"),
    gcFirst = FALSE
  )

  # Estimate higher order for parallel analysis with PAF
  time_list$paf_time <- system.time(
    paf_higher <- efa_custom(simulated$data, pa_lower[["paf"]], "paf"),
    gcFirst = FALSE
  )

  # Set maximum dimensions
  k_max <- 2 * simulated$parameters$lower$factors

  # Estimate lower order for MAP and VSS
  time_list$map_vss_time <- system.time(
    map_vss_lower <- map_vss(R, N, k_max),
    gcFirst = FALSE
  )

  # Estimate loadings
  time_list$map_time <- system.time(
    map_loadings <- efa(R, map_vss_lower[["map"]], corr = "polychoric", n.obs = N)$loadings,
    gcFirst = FALSE
  )
  time_list$vss1_time <- system.time(
    vss1_loadings <- efa(R, map_vss_lower[["vss1"]], corr = "polychoric", n.obs = N)$loadings,
    gcFirst = FALSE
  )
  time_list$vss2_time <- system.time(
    vss2_loadings <- efa(R, map_vss_lower[["vss2"]], corr = "polychoric", n.obs = N)$loadings,
    gcFirst = FALSE
  )

  # Determine minimum and maximum to search between
  k_range <- range(c(pa_lower, map_vss_lower))

  # Estimate lower order for maximum loading
  time_list$max_time <- system.time(
    max_lower <- max_loading(
      R, N,
      min_k = max(1, k_range[1] - 2),
      max_k = min(k_range[2] + 2, length(variables))
    ),
    gcFirst = FALSE
  )

  # Return results
  return(
    list(
      loadings = list(
        lower = list(
          lower_louvain = lower_louvain$loadings$lower,
          prism = prism_result$loadings$lower,
          pca = pca_higher$loadings$lower,
          paf = paf_higher$loadings$lower,
          map = map_loadings,
          vss1 = vss1_loadings,
          vss2 = vss2_loadings,
          max = max_lower$loadings
        ),
        higher = list(
          lower_louvain = lower_louvain$loadings$higher,
          prism = prism_result$loadings$higher,
          pca = pca_higher$loadings$higher,
          paf = paf_higher$loadings$higher,
          map = NA,
          vss1 = NA,
          vss2 = NA,
          max = NA
        )
      ),
      dimensions = list(
        lower = list(
          lower_louvain = lower_louvain$dimensions$lower,
          prism = prism_result$dimensions$lower,
          pca = pa_lower[["pca"]],
          paf = pa_lower[["paf"]],
          map = map_vss_lower[["map"]],
          vss1 = map_vss_lower[["vss1"]],
          vss2 = map_vss_lower[["vss2"]],
          max = dim(max_lower$loadings)[2]
        ),
        higher = list(
          lower_louvain = lower_louvain$dimensions$higher,
          prism = prism_result$dimensions$higher,
          pca = pca_higher$dimensions$higher,
          paf = paf_higher$dimensions$higher,
          map = NA,
          vss1 = NA,
          vss2 = NA,
          max = NA
        )
      ),
      timing = EGAnet:::nvapply(time_list, function(x){x[["elapsed"]]})
    )
  )

}

# Compile
obtain_methods <- cmpfun(obtain_methods)
