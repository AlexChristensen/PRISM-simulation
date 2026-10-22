# Load packages
library(psych)

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



# Set directory and files
DIR <- "./data/extreme/"
FILES <- list.files(DIR)

# Loop over files
skew_df <- do.call(
  rbind.data.frame, lapply(FILES, function(file){

    # Print file
    cat("\r", file)

    # Load condition
    load(paste0(DIR, file))

    # Return data frame
    return(
      data.frame(
        condition_number = gsub("_.*", "", gsub("condition-", "", file)),
        replicate_number = gsub(".RData*", "", gsub(".*replicate-", "", file)),
        values = skew(extreme_data$no_overlap$`10000`)
      )
    )

  })
)

# Save skew
save(skew_df, file = "skew_df.RData")
