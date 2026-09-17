#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#### Hierarchical | Timing ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Load packages
library(tidyverse)

# Set directory and files
DIR <- "./timing/"
FILES <- list.files(DIR)

# Loop over conditions
timing_df <- do.call(
  rbind.data.frame, lapply(paste0("condition-", EGAnet:::format_integer(seq_len(160), 2)), function(number){

    # Select files
    select <- FILES[grep(number, FILES)]

    # Load files
    timing_condition <- as.data.frame(do.call(
      rbind, lapply(select, function(file){
        load(paste0(DIR, file))
        timing <- as.data.frame(timing)
        timing$N <- gsub("_.*", "", gsub("N-", "", select))
        return(timing)
      })
    ))

    # Add condition
    timing_condition$condition <- gsub("condition-", "", number)

    # Return timing condition
    return(timing_condition)

  })
)

# Save timing
save(timing_df, file = "timing_df.RData")

# Make longer
timing_long <- timing_df %>%
  pivot_longer(
    cols = colnames(timing_df)[-c(14:15)],
    names_to = "method",
    values_to = "seconds"
  )

# Breakdown by condition
timing_condition <- timing_long %>%
  group_by(method, condition) %>%
  filter(method == "prism_time") %>%
  summarize(average_seconds = mean(seconds, na.rm = TRUE)) %>%
  as.data.frame()

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

# Combine conditions with timing
combined <- cbind.data.frame(conditions, timing_condition)

# Determine effect
summary(lm(log(average_seconds) ~ log(LowF * HighF), data = combined))
plot(
  x = log(combined$LowF * combined$HighF),
  y = log(combined$average_seconds)
)
