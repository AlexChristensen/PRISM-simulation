#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#### Hierarchical Simulation | Tables and Figures ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Create directory
dir.create("./figures")

# Load packages
library(ggplot2); library(ggpubr); library(tidyverse)

# Load data
load("./results.RData")

# Check max loading correct depth to search
## Set range to search
search_methods <- results[!(results$METHOD %in% c("prism", "lower_louvain", "max")),]
search_length <- dim(search_methods)[1]
start_sequence <- seq(1, search_length, 5)
end_sequence <- seq(5, search_length, 5)

## Set number of lower order factors
lower_order <- search_methods$LowF * search_methods$HighF

## Update lower order
lower_order <- lower_order[start_sequence]

## Initialize low end and high end
high_end <- low_end <- numeric(length(start_sequence))

# Total variables
total_variables <- search_methods$LowV * search_methods$LowF * search_methods$HighF

## Update total variables
total_variables <- total_variables[start_sequence]

## Populate low end and high end
for(i in seq_along(start_sequence)){

  ## Set values
  low_end[i] <- max(
    1, min(search_methods$lower_dimensions[start_sequence[i]:end_sequence[i]]) - 2
  )
  high_end[i] <- min(
    total_variables[i], max(search_methods$lower_dimensions[start_sequence[i]:end_sequence[i]]) + 2
  )

}

# Check for truth in range
max_range <- (lower_order <= high_end) & (lower_order >= low_end)
mean(max_range)

#%%%%%%%%%%%%%%%%%%%%#
#### Figure Theme ----
#%%%%%%%%%%%%%%%%%%%%#

# Set colors
COLORS <- c(
  "prism" = "#f35b04",
  "lower_louvain" = "#fca311",
  "max" = "#b56576",
  "map" = "#e56b6f",
  "pca" = "#7cb518",
  "paf" = "#5c8001",
  "vss1" = "#526a98",
  "vss2" = "#4caec2"
)

# Set labels
LABELS <- c(
  "prism" = "PRISM",
  "lower_louvain" = "Lower Order Louvain",
  "max" = "Max Loading",
  "map" = "Maximum A Posteriori",
  "paf" = "Parallel Analysis (PAF)",
  "pca" = "Parallel Analysis (PCA)",
  "vss1" = "Very Simple Structure 1",
  "vss2" = "Very Simple Structure 2"
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#### Table 1: Overall First Order Results ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Break down by method and correlation
full_results <- results %>%
  group_by(METHOD, Overlap, Error, SKEW) %>%
  summarize(
    lower_correct_mean = mean(LOWER_CORRECT, na.rm = TRUE),
    lower_ari_mean = mean(lower_ari, na.rm = TRUE),
    lower_mbe_mean = mean(LOWER_MBE, na.rm = TRUE),
    lower_mae_mean = mean(abs(lower_dimensions - (LowF * HighF)))
  ) %>% as.data.frame() %>% print(digits = 3)

# Results for underfactoring and overfactoring
results %>%
  group_by(METHOD) %>%
  summarize(
    prop_under = mean(LOWER_MBE < 0),
    prop_over = mean(LOWER_MBE > 0),
    prop_correct = mean(LOWER_MBE == 0)
  ) %>% as.data.frame()

# Check on max range
max_results <- results[max_range,] %>%
  group_by(METHOD, Overlap, Error, SKEW) %>%
  summarize(
    lower_correct_mean = mean(LOWER_CORRECT, na.rm = TRUE),
    lower_ari_mean = mean(lower_ari, na.rm = TRUE),
    lower_mbe_mean = mean(LOWER_MBE, na.rm = TRUE)
  ) %>% as.data.frame()

# Compare differences
full_results[full_results$METHOD == "max",-c(1:4)] - max_results[max_results$METHOD == "max",-c(1:4)]

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#### Figure 3: Effect of Factors and Sample Size ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Set up summary
condition_summary <- results %>%
  group_by(METHOD, LowF, HighF, N) %>%
  summarize(
    lower_correct_mean = mean(LOWER_CORRECT, na.rm = TRUE),
    lower_ari_mean = mean(lower_ari, na.rm = TRUE),
    higher_correct_mean = mean(HIGHER_CORRECT, na.rm = TRUE),
    higher_ari_mean = mean(higher_ari, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = c(lower_correct_mean, lower_ari_mean, higher_correct_mean, higher_ari_mean),
    names_to = c("Order", "Metric"),
    names_pattern = "(lower|higher)_(correct|ari)_mean",
    values_to = "value"
  ) %>% as.data.frame()

# Set factors
factor_cols <- seq_len(6)
condition_summary[,factor_cols] <- do.call(cbind.data.frame, lapply(condition_summary[,factor_cols], factor))
condition_summary$Metric <- factor(condition_summary$Metric, levels = c("correct", "ari"))
condition_summary$Order <- factor(condition_summary$Order, levels = c("lower", "higher"))

# Restrict to first order only for Figure 1
condition_summary_lower <- condition_summary %>%
  filter(Order == "lower")

# Factor the methods
condition_summary_lower$METHOD <- factor(
  condition_summary_lower$METHOD,
  levels = c(
    "prism", "lower_louvain",
    "max", "map", "paf", "pca",
    "vss1", "vss2"
  )
)

# Figure 3 (First Order Accuracy and Omega Index only)
figure3 <- ggplot(
  data = condition_summary_lower,
  aes(x = N, y = value, group = METHOD)
) +
  facet_grid(
    rows = vars(Metric), cols = vars(HighF, LowF), switch = "y",
    labeller = labeller(
      Metric = as_labeller(c("correct" = "Accuracy", "ari" = "omega"),
                           default = label_parsed),
      LowF = c(
        "3" = "First Order = 3",
        "5" = "First Order = 5"
      ),
      HighF = c(
        "1" = "Second Order = 1",
        "2" = "Second Order = 2",
        "4" = "Second Order = 4",
        "6" = "Second Order = 6"
      )
    )
  ) +
  geom_hline(yintercept = seq(0, 1, 0.25), linewidth = 0.3, color = "lightgrey") +
  geom_line(
    aes(color = METHOD), position = position_dodge(0.9),
    linewidth = 0.5, alpha = 0.5
  ) +
  geom_point(
    aes(fill = METHOD), position = position_dodge(0.9),
    size = 2, shape = 21, stroke = 0.25, color = "white"
  ) +
  scale_color_manual(name = "Method", labels = LABELS, values = COLORS) +
  scale_fill_manual(name = "Method", labels = LABELS, values = COLORS) +
  scale_y_continuous(
    limits = c(0.00, 1.05),
    breaks = seq(0.00, 1.00, 0.25),
    labels = EGAnet:::format_decimal(seq(0.00, 1.00, 0.25), 2),
    expand = c(0, 0),
    position = "right"
  ) +
  labs(x = "Sample Size", title = "First Order Accuracy") +
  theme(
    panel.background = element_blank(),
    panel.spacing.y = unit(0.5, "cm"),
    plot.title = element_text(family = "ubuntu", size = 12, face = "bold", hjust = 0.5),
    axis.line = element_line(linewidth = 0.3, color = "black"),
    axis.line.y = element_blank(),
    axis.ticks = element_line(linewidth = 0.3),
    axis.ticks.y = element_blank(),
    axis.title = element_text(family = "ubuntu", size = 10),
    axis.title.y = element_blank(),
    axis.text = element_text(family = "ubuntu", size = 8),
    axis.text.x = element_text(family = "ubuntu", angle = 45, hjust = 1),
    strip.background = element_rect(color = "black", fill = "white", linewidth = 0.3),
    strip.text = element_text(family = "ubuntu", size = 8),
    legend.title = element_text(family = "ubuntu", size = 10, hjust = 0.5),
    legend.text = element_text(family = "ubuntu", size = 8),
    legend.position = "bottom",
    legend.title.position = "bottom"
  ); figure3

# Save plot
ggsave(
  figure3, filename = "./figures/figure3.pdf",
  height = 6, width = 10, bg = "white", device = cairo_pdf
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#### Figure 4: Effect of Factors and Sample Size ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Restrict to second order only for Figure 5
condition_summary_higher <- condition_summary %>%
  filter(Order == "higher")

# Set methods
condition_summary_higher <- condition_summary_higher %>%
  filter(METHOD %in% c("lower_louvain", "prism", "paf", "pca"))

# Factor the methods
condition_summary_higher$METHOD <- factor(
  condition_summary_higher$METHOD,
  levels = c("prism", "lower_louvain", "paf", "pca")
)

# Figure 4 (Second Order Accuracy and Omega Index only)
figure4 <- ggplot(
  data = condition_summary_higher,
  aes(x = N, y = value, group = METHOD)
) +
  facet_grid(
    rows = vars(Metric), cols = vars(HighF, LowF), switch = "y",
    labeller = labeller(
      Metric = as_labeller(c("correct" = "Accuracy", "ari" = "omega"),
                           default = label_parsed),
      LowF = c(
        "3" = "First Order = 3",
        "5" = "First Order = 5"
      ),
      HighF = c(
        "1" = "Second Order = 1",
        "2" = "Second Order = 2",
        "4" = "Second Order = 4",
        "6" = "Second Order = 6"
      )
    )
  ) +
  geom_hline(yintercept = seq(0, 1, 0.25), linewidth = 0.3, color = "lightgrey") +
  geom_line(
    aes(color = METHOD), position = position_dodge(0.9),
    linewidth = 0.5, alpha = 0.5
  ) +
  geom_point(
    aes(fill = METHOD), position = position_dodge(0.9),
    size = 2, shape = 21, stroke = 0.25, color = "white"
  ) +
  scale_color_manual(name = "Method", labels = LABELS, values = COLORS) +
  scale_fill_manual(name = "Method", labels = LABELS, values = COLORS) +
  scale_y_continuous(
    limits = c(0.00, 1.05),
    breaks = seq(0.00, 1.00, 0.25),
    labels = EGAnet:::format_decimal(seq(0.00, 1.00, 0.25), 2),
    expand = c(0, 0),
    position = "right"
  ) +
  labs(x = "Sample Size", title = "Second Order Accuracy") +
  theme(
    panel.background = element_blank(),
    panel.spacing.y = unit(0.5, "cm"),
    plot.title = element_text(family = "ubuntu", size = 12, face = "bold", hjust = 0.5),
    axis.line = element_line(linewidth = 0.3, color = "black"),
    axis.line.y = element_blank(),
    axis.ticks = element_line(linewidth = 0.3),
    axis.ticks.y = element_blank(),
    axis.title = element_text(family = "ubuntu", size = 10),
    axis.title.y = element_blank(),
    axis.text = element_text(family = "ubuntu", size = 8),
    axis.text.x = element_text(family = "ubuntu", angle = 45, hjust = 1),
    strip.background = element_rect(color = "black", fill = "white", linewidth = 0.3),
    strip.text = element_text(family = "ubuntu", size = 8),
    legend.title = element_text(family = "ubuntu", size = 10, hjust = 0.5),
    legend.text = element_text(family = "ubuntu", size = 8),
    legend.position = "bottom",
    legend.title.position = "bottom"
  ); figure4

# Save plot
ggsave(
  figure4, filename = "./figures/figure4.pdf",
  height = 6, width = 10, bg = "white", device = cairo_pdf
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#### Figure 5: Worst and Best Accuracy by Condition ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Set factor labels
FACTOR_LABELS <- c(
  LowV = "Variables per first order factor", LowF = "First order factors\nper second order factor",
  LowL = "First order loading size", HighF = "Number of second order factors",
  HighL = "Second order loading size", HighC = "Second order factor correlation",
  N = "Sample size", SKEW = "Skew", Overlap = "Overlapping items", Error = "Population error"
)

# For each design factor, get the range (worst level to best level) per method
compute_range <- function(factor_var, data, outcome_var) {
  data %>%
    mutate(.f = .data[[factor_var]]) %>%
    group_by(METHOD, .f) %>%
    summarize(m = mean(.data[[outcome_var]], na.rm = TRUE), .groups = "drop") %>%
    group_by(METHOD) %>%
    summarize(low = min(m), high = max(m), .groups = "drop") %>%
    mutate(Factor = factor_var)
}

# Compute every design factor's effect against both outcomes, so first and
# second order share the same set of conditions on the y-axis
all_factor_vars <- c("LowV", "LowF", "LowL", "HighL", "HighF", "HighC", "N", "SKEW", "Overlap", "Error")

lower_effects <- bind_rows(lapply(
  all_factor_vars, compute_range, data = results, outcome_var = "LOWER_CORRECT"
)) %>% mutate(Order = "lower")
higher_effects <- bind_rows(lapply(
  all_factor_vars, compute_range, data = results, outcome_var = "HIGHER_CORRECT"
)) %>% mutate(Order = "higher")

effects <- bind_rows(lower_effects, higher_effects)

# Order factors by their average impact across methods and orders
factor_order <- c(
  "Error", "Overlap", "HighC", "HighF", "HighL", "LowF", "LowL", "LowV", "SKEW", "N"
)

effects$Factor <- factor(effects$Factor, levels = factor_order)
effects$Order <- factor(effects$Order, levels = c("lower", "higher"))

# Offset methods vertically so dumbbells don't overlap
METHOD_OFFSET <- setNames(
  seq(0.42, -0.42, length.out = 8),
  c("prism", "lower_louvain", "max", "map", "paf", "pca", "vss1", "vss2")
)
effects$y <- as.numeric(effects$Factor) + METHOD_OFFSET[as.character(effects$METHOD)]

y_breaks <- seq_len(nlevels(effects$Factor))
y_labels <- FACTOR_LABELS[levels(effects$Factor)]

# Shade every other row (no Order column, so it repeats across both facets)
row_bands <- data.frame(y_pos = y_breaks) %>% filter(y_pos %% 2 == 0)

# Set methods factor
effects$METHOD <- factor(
  effects$METHOD,
  levels = c("prism", "lower_louvain", "max", "map", "paf", "pca", "vss1", "vss2")
)

# Figure 5
figure5 <- ggplot(effects) +
  facet_grid(
    cols = vars(Order),
    labeller = labeller(Order = c("lower" = "First Order Accuracy", "higher" = "Second Order Accuracy"))
  ) +
  geom_rect(
    data = row_bands,
    aes(ymin = y_pos - 0.5, ymax = y_pos + 0.5, xmin = -Inf, xmax = 1),
    inherit.aes = FALSE, fill = "grey93"
  ) +
  geom_vline(xintercept = seq(0.00, 1, 0.10), linewidth = 0.3, color = "lightgrey") +
  geom_segment(
    aes(x = low, xend = high, y = y, yend = y, color = METHOD),
    linewidth = 0.75, alpha = 0.5, lineend = "round"
  ) +
  geom_point(aes(x = low, y = y, color = METHOD, shape = "Worst"), fill = "white", size = 2.5) +
  geom_point(aes(x = high, y = y, fill = METHOD, shape = "Best"), color = "white", size = 2.5, stroke = 0.25) +
  scale_color_manual(
    name = "Method", labels = LABELS, values = COLORS,
    guide = guide_legend(order = 1, override.aes = list(shape = NA, alpha = 1))
  ) +
  scale_fill_manual(
    name = "Method", labels = LABELS, values = COLORS,
    guide = guide_legend(order = 1, override.aes = list(shape = NA))
  ) +
  scale_shape_manual(
    name = "Performance", breaks = c("Worst", "Best"),
    values = c("Worst" = 21, "Best" = 21),
    guide = guide_legend(
      order = 2,
      override.aes = list(fill = c("white", "black"), color = c("black", "white"))
    )
  ) +
  scale_y_continuous(
    limits = c(0.25, 10.75),
    breaks = y_breaks, labels = y_labels, expand = c(0,0)
  ) +
  scale_x_continuous(
    limits = c(0.00, 1.025), breaks = seq(0.00, 1.00, 0.10),
    labels = EGAnet:::format_decimal(seq(0.00, 1.00, 0.10), 2),
    expand = c(0, 0)
  ) +
  theme(
    panel.background = element_blank(),
    panel.spacing.x = unit(0.75, "cm"),
    axis.line = element_line(linewidth = 0.3, color = "black"),
    axis.ticks = element_line(linewidth = 0.3),
    axis.title = element_text(family = "ubuntu", size = 14),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text = element_text(family = "ubuntu", size = 12),
    strip.background = element_rect(color = "black", fill = "white", linewidth = 0.3),
    strip.text = element_text(family = "ubuntu", size = 12),
    legend.title = element_text(family = "ubuntu", size = 12, hjust = 0.5),
    legend.text = element_text(family = "ubuntu", size = 10),
    legend.position = "bottom",
    legend.title.position = "left",
    legend.box = "vertical",
    legend.spacing.y = unit(0, "cm"),
    legend.margin = margin(t = 0, b = 0)
  ); figure5

# Save plot
ggsave(
  figure5, filename = "./figures/figure5.pdf",
  height = 8, width = 12, bg = "white", device = cairo_pdf
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#### Figure 6: Most Challenging Conditions ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Evaluate what's happening with PRISM
results %>%
  filter(METHOD == "prism", N == 10000) %>%
  group_by(METHOD, Overlap, LowV, LowF, LowL, HighF, HighL, HighC) %>%
  summarize(
    lower_acc = mean(LOWER_CORRECT),
    lower_mbe = mean(LOWER_MBE)
  ) %>%
  filter(lower_acc < 0.80) %>%
  as.data.frame()


# Flag the four conditions that, together, collapse first-order recovery:
# cross-loadings present (O), few indicators per factor (V), weak first-order
# loadings (F), strong second-order loadings (S)
# NOTE: flags are added to the full results (all N) rather than
# overwriting it with the N = 10000 subset, since the figure below needs
# both the pooled-across-N and the N = 10000-only versions.
results$O <- results$Overlap == TRUE
results$V <- results$LowV == 5
results$F <- results$LowL == 0.50
results$S <- results$HighL == 0.70

adverse_factors <- c("O", "V", "F", "S")

# All 16 disjoint cells: the empty set ("None adverse") through all four
# conditions adverse together (OVFS). Each row belongs to exactly one
# cell, so nothing here is pooled across incompatible conditions.
adverse_subsets <- c(
  list(character(0)),
  unlist(lapply(4:1, function(k) combn(adverse_factors, k, simplify = FALSE)), recursive = FALSE)
)

# Spell out the single-condition (and "None") labels; combos stay as letter codes
CONDITION_LABELS <- c(
  None = "All Other\nConditions",
  O = "O", # Overlap (O) = Present",
  V = "V", # "First Order Variables (V) = 5",
  F = "F", # "First Order Loadings (F) = 0.50",
  S = "S" # "Second Order Loadings (S) = 0.70"
)

# A "pure" cell: the named factors are adverse AND every other factor is at
# its favorable level -- mutually exclusive across all 16 subsets
pure_mask <- function(data, subset_vars) {
  other_vars <- setdiff(adverse_factors, subset_vars)
  mask_list <- c(
    lapply(subset_vars, function(v) data[[v]]),
    lapply(other_vars, function(v) !data[[v]])
  )
  Reduce(`&`, mask_list)
}

compute_pure <- function(subset_vars, data, outcome_var) {
  mask <- pure_mask(data, subset_vars)
  data[mask, ] %>%
    group_by(METHOD) %>%
    summarize(accuracy = mean(.data[[outcome_var]], na.rm = TRUE), .groups = "drop") %>%
    mutate(
      subset_label = if (length(subset_vars) == 0) "None" else paste(subset_vars, collapse = ""),
      subset_size  = length(subset_vars)
    )
}

# Row 1: pooled across every sample size; Row 2: N = 10,000 only
other_results <- results[results$N != 10000, ]
results_10k <- results[results$N == 10000, ]

# Restrict to the four focal methods for this figure
STORM_METHODS <- c("prism", "lower_louvain", "paf", "pca")
other_results_storm <- other_results %>% filter(METHOD %in% STORM_METHODS)
results_10k_storm <- results_10k %>% filter(METHOD %in% STORM_METHODS)

# Cross sample-size panel with order (first/second), so all four
# combinations of panel x order are computed
storm_data <- bind_rows(
  bind_rows(lapply(adverse_subsets, compute_pure, data = other_results_storm, outcome_var = "LOWER_CORRECT")) %>%
    mutate(panel = "All Other Sample Sizes", Order = "lower"),
  bind_rows(lapply(adverse_subsets, compute_pure, data = results_10k_storm, outcome_var = "LOWER_CORRECT")) %>%
    mutate(panel = "Sample Size = 10,000", Order = "lower"),
  bind_rows(lapply(adverse_subsets, compute_pure, data = other_results_storm, outcome_var = "HIGHER_CORRECT")) %>%
    mutate(panel = "All Other Sample Sizes", Order = "higher"),
  bind_rows(lapply(adverse_subsets, compute_pure, data = results_10k_storm, outcome_var = "HIGHER_CORRECT")) %>%
    mutate(panel = "Sample Size = 10,000", Order = "higher")
) %>%
  mutate(
    panel = factor(panel, levels = c("All Other Sample Sizes", "Sample Size = 10,000")),
    Order = factor(Order, levels = c("lower", "higher"))
  )

# Order subsets by size band (0 -> 4), and within each band by PRISM
# Louvain's (EGA Louvain) first-order accuracy at N = 10,000; all four
# panels share this order so they stay directly comparable
FOCAL_METHOD <- "prism"

# Fixed condition-combination order (replaces the focal-accuracy sort)
subset_label_order <- c(
  "None", "O", "V", "F", "S", "OV", "OF", "OS",
  "VF", "VS", "FS", "OVF", "OVS", "OFS", "VFS", "OVFS"
)

subset_order <- data.frame(
  subset_label = subset_label_order,
  subset_size  = ifelse(subset_label_order == "None", 0L, nchar(subset_label_order)),
  y_base       = seq_along(subset_label_order)
)

y_lookup <- setNames(subset_order$y_base, subset_order$subset_label)
y_labels_ordered <- ifelse(
  subset_order$subset_size <= 1,
  CONDITION_LABELS[subset_order$subset_label],
  subset_order$subset_label
)

storm_data$y_base <- y_lookup[storm_data$subset_label]

# Offset the four focal methods symmetrically around each category's
# base y-position; kept as a separate object so Figure 4's layout is untouched
METHOD_OFFSET_FIG3 <- setNames(
  seq(-0.2, 0.2, length.out = 4),
  STORM_METHODS
)
storm_data$y <- storm_data$y_base + METHOD_OFFSET_FIG3[as.character(storm_data$METHOD)]
storm_data <- storm_data %>% arrange(panel, Order, METHOD, y_base)

# Set method factors
storm_data$METHOD <- factor(
  storm_data$METHOD,
  levels = c("prism", "lower_louvain", "paf", "pca")
)

# Shade the two middle size-bands (1, 3 adverse conditions); no panel/Order
# column here, so it repeats identically across all four facets
subset_bounds <- storm_data %>%
  distinct(subset_label, subset_size, y_base) %>%
  group_by(subset_size) %>%
  summarize(ymin = min(y_base) - 0.5, ymax = max(y_base) + 0.5, .groups = "drop")
shaded_bounds <- subset_bounds %>% filter(subset_size %in% c(1, 3))

# One order's storm plot (the two sample-size panels, stacked); shared
# between the first- and second-order halves of Figure 6 so they read
# identically before being stacked via ggarrange
build_storm_plot <- function(order_label, x_lab, show_x_axis = TRUE) {
  ggplot() +
    facet_grid(rows = vars(panel)) +
    geom_rect(
      data = shaded_bounds,
      aes(ymin = ymin, ymax = ymax, xmin = -Inf, xmax = 1),
      inherit.aes = FALSE, fill = "grey93"
    ) +
    geom_vline(xintercept = seq(0, 1, 0.25), linewidth = 0.3, color = "grey60") +
    geom_path(
      data = storm_data %>% filter(Order == order_label),
      aes(x = accuracy, y = y, color = METHOD, group = METHOD),
      linewidth = 0.8, alpha = 0.50
    ) +
    geom_point(
      data = storm_data %>% filter(Order == order_label),
      aes(x = accuracy, y = y, fill = METHOD),
      shape = 21, color = "white", size = 3
    ) +
    scale_color_manual(
      name = "Method", labels = LABELS, values = COLORS,
      guide = guide_legend(override.aes = list(linetype = 1))
    ) +
    scale_fill_manual(
      name = "Method", labels = LABELS, values = COLORS,
      guide = guide_legend(override.aes = list(shape = 21))
    ) +
    scale_y_continuous(breaks = seq_along(y_labels_ordered), labels = y_labels_ordered, expand = expansion(add = 0.7)) +
    scale_x_continuous(
      limits = c(0.00, 1.05), breaks = seq(0.00, 1.00, 0.25),
      labels = EGAnet:::format_decimal(seq(0.00, 1.00, 0.25), 2),
      expand = c(0, 0)
    ) +
    labs(x = x_lab, y = NULL) +
    coord_flip() +
    theme(
      panel.background = element_blank(),
      panel.spacing.y = unit(0.75, "cm"),
      axis.line = element_line(linewidth = 0.3, color = "black"),
      axis.ticks = element_line(linewidth = 0.3),
      axis.title = element_text(family = "ubuntu", size = 12),
      axis.text = element_text(family = "ubuntu", size = 10),
      # axis.text.x = element_text(angle = 45, hjust = 1),
      strip.background = element_rect(color = "black", fill = "white", linewidth = 0.3),
      strip.text = element_text(family = "ubuntu", size = 10),
      legend.title = element_text(family = "ubuntu", size = 12, hjust = 0.5),
      legend.text = element_text(family = "ubuntu", size = 10),
      legend.position = "bottom",
      legend.title.position = "bottom",
      legend.key.width = unit(1, "cm")
    ) +
    if (!show_x_axis) {
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.title.x = element_blank()
      )
    }
}

# The condition-combination axis is identical across both halves of the
# stack, so it is only shown once, on the bottom (second-order) plot
figure6_lower <- build_storm_plot("lower", "First Order Accuracy", show_x_axis = FALSE)
figure6_higher <- build_storm_plot("higher", "Second Order Accuracy")

# Stack the first- and second-order plots into one 4-panel figure with a
# single shared legend
figure6 <- ggarrange(
  figure6_lower, figure6_higher,
  ncol = 1, nrow = 2, heights = c(0.95, 1.075),
  common.legend = TRUE, legend = "bottom"
)
figure6 <- annotate_figure(
  figure6,
  top = text_grob(
    "Overlapping Items | Variables per First Order Factor = 5 | First Order Loadings = 0.50 | Second Order Loadings = 0.70",
    size = 12, hjust = 0.50, color = "grey60"
  )
)
figure6 <- annotate_figure(
  figure6,
  top = text_grob(
    "Most Challenging Condition Combinations",
    face = "bold", size = 14, hjust = 1.08
  )
); figure6

# Save plot
ggsave(
  figure6, filename = "./figures/figure6.pdf",
  height = 9, width = 10, bg = "white", device = cairo_pdf
)
