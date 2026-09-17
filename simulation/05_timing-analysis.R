#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#### Hierarchical | Timing ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

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





