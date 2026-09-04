## =============================================================================
## LENGTH-BASED INDICATOR ANALYSIS
## Species: Rhizoprionodon porosus ("tuba") and Scomberomorus brasiliensis ("serra")
##
## NOTE ON SPECIES NAMING:
##   The user's original request referred to "Scomberomorus cavalla" as the
##   species indicated by "serra". However, the uploaded dataset's `species`
##   column contains "Scomberomorus brasiliensis", not "S. cavalla". In Brazil,
##   "serra" is the common name conventionally used for S. brasiliensis, while
##   S. cavalla is usually called "cavala". This script therefore treats the
##   "serra" group as Scomberomorus brasiliensis, matching what is actually
##   present in the data (confirmed with the user).
##
## Length metric used: fork length (fl)
##
## Categories assigned to each individual (four exhaustive, non-overlapping
## bins built from the three thresholds L50, Lopt, and MS = Lopt*1.1):
##   1. "Below_L50"        : fl <  L50
##   2. "L50_to_Lopt"       : L50 <= fl <  Lopt
##   3. "At_Lopt"           : Lopt <= fl < MS      (optimum-length range up to
##                                                    the mega-spawner cutoff)
##   4. "MS_megaspawner"    : fl >= MS             (Lopt + 10%; Froese 2004
##                                                    "mega-spawner" cutoff)
##
## Sensitivity analysis: thresholds (L50, Lopt, MS) are each rescaled by
## -20%, -10%, 0% (base), +10%, +20% and the whole categorization is redone,
## so you can see how the % of individuals in each bin shifts as thresholds
## move.
##
## Bootstrap resampling: for each species, the length-frequency data are
## resampled with replacement (n_boot iterations) to build 95% confidence
## intervals around (a) the proportion of individuals in each length category
## and (b) the length-frequency histogram/density itself.
## =============================================================================

## ---- 0. SETUP --------------------------------------------------------------

required_packages = c("dplyr", "tidyr", "readr", "purrr", "tibble", "ggplot2", "scales")
new_packages = required_packages[!(required_packages %in% installed.packages()[, "Package"])]
if (length(new_packages) > 0) install.packages(new_packages, repos = "https://cloud.r-project.org")

library(dplyr)
library(tidyr)
library(readr)
library(purrr)
library(tibble)
library(ggplot2)
library(scales)

set.seed(123)  # reproducibility for bootstrap resampling

## ---- USER-EDITABLE PATHS ----------------------------------------------------
input_csv   = "length_data_fixed.csv"     # path to the uploaded data file
output_dir  = "outputs"                   # folder for all script outputs
dir.create(output_dir, showWarnings = FALSE)

n_boot = 1000                              # number of bootstrap resamples
ci_level = 0.95                            # confidence level for bootstrap CIs

## ---- 1. PARAMETERS & THRESHOLDS --------------------------------------------

params = tibble::tibble(
  species_common = c("tuba", "serra"),
  species_sci    = c("Rhizoprionodon porosus", "Scomberomorus brasiliensis"),
  l50            = c(54.4, 42.1),
  m              = c(0.765, 0.36),
  k              = c(0.36, 0.15),
  loo            = c(106.82, 109.18)
) %>%
  mutate(
    lopt = loo * (3 / (3 + (m / k))),
    ms   = lopt + (lopt * 0.1)
  )

print(params)

## ---- 2. LOAD & CLEAN DATA ---------------------------------------------------

raw_data = read_csv(input_csv, show_col_types = FALSE)

# Keep only the two species of interest, using fork length (fl)
data_clean = raw_data %>%
  filter(species %in% params$species_sci) %>%
  select(species, fl) %>%
  mutate(fl = as.numeric(fl)) %>%
  left_join(
    params %>% select(species_sci, species_common),
    by = c("species" = "species_sci")
  )

n_missing_fl = sum(is.na(data_clean$fl))
if (n_missing_fl > 0) {
  message(sprintf(
    "Dropping %d rows with missing/non-numeric fork length (fl) values.",
    n_missing_fl
  ))
}

data_clean = data_clean %>% filter(!is.na(fl), fl > 0)

cat("\nSample sizes after cleaning:\n")
print(table(data_clean$species_common))

## ---- 3. CATEGORIZATION FUNCTION --------------------------------------------

# Assigns each length value to one of four categories given threshold vectors
# (l50, lopt, ms). Vectorized with case_when so it works even when different
# rows carry different (species-specific) thresholds.
categorize_lengths = function(fl, l50, lopt, ms) {
  out = case_when(
    fl <  l50               ~ "Below_L50",
    fl >= l50  & fl < lopt  ~ "L50_to_Lopt",
    fl >= lopt & fl < ms    ~ "At_Lopt",
    fl >= ms                ~ "MS_megaspawner",
    TRUE                    ~ NA_character_
  )
  factor(out, levels = c("Below_L50", "L50_to_Lopt", "At_Lopt", "MS_megaspawner"))
}

## ---- 4. BASE-CASE CATEGORIZATION -------------------------------------------

data_categorized = data_clean %>%
  left_join(params, by = "species_common") %>%
  mutate(category = categorize_lengths(fl, l50, lopt, ms))

base_summary = data_categorized %>%
  count(species_common, category, .drop = FALSE) %>%
  group_by(species_common) %>%
  mutate(
    total = sum(n),
    proportion = n / total,
    percent = percent(proportion, accuracy = 0.1)
  ) %>%
  ungroup()

cat("\n--- Base-case length-category summary ---\n")
print(base_summary)

write_csv(base_summary, file.path(output_dir, "base_category_summary.csv"))

## ---- 5. SENSITIVITY ANALYSIS ON THRESHOLDS ---------------------------------

# Rescale ALL three thresholds (L50, Lopt, MS) together by each factor
# and recompute the categorization. This mirrors "what if the thresholds
# were X% higher/lower" as requested.
sensitivity_factors = c(0.80, 0.90, 1.00, 1.10, 1.20)
sensitivity_labels  = c("-20%", "-10%", "Base", "+10%", "+20%")

run_sensitivity = function(factor_value, factor_label) {
  params_adj = params %>%
    mutate(
      l50_adj  = l50  * factor_value,
      lopt_adj = lopt * factor_value,
      ms_adj   = ms   * factor_value
    )
  
  data_clean %>%
    left_join(params_adj, by = "species_common") %>%
    mutate(category = categorize_lengths(fl, l50_adj, lopt_adj, ms_adj)) %>%
    count(species_common, category, .drop = FALSE) %>%
    group_by(species_common) %>%
    mutate(
      total = sum(n),
      proportion = n / total
    ) %>%
    ungroup() %>%
    mutate(
      threshold_factor = factor_value,
      threshold_label  = factor_label
    )
}

sensitivity_results = map2_dfr(
  sensitivity_factors, sensitivity_labels, run_sensitivity
) %>%
  mutate(threshold_label = factor(threshold_label, levels = sensitivity_labels))

cat("\n--- Sensitivity analysis (thresholds scaled -20% to +20%) ---\n")
print(sensitivity_results %>% select(species_common, threshold_label, category, proportion))

write_csv(sensitivity_results, file.path(output_dir, "sensitivity_analysis.csv"))

# Plot: how the proportion in each category shifts as thresholds move
sensitivity_plot = ggplot(
  sensitivity_results,
  aes(x = threshold_label, y = proportion, fill = category, group = category)
) +
  geom_col(position = "stack") +
  facet_wrap(~species_common) +
  scale_y_continuous(labels = percent_format()) +
  labs(
    title = "Sensitivity of length-category composition to threshold scaling",
    x = "Threshold adjustment (L50, Lopt, MS scaled together)",
    y = "Proportion of individuals",
    fill = "Category"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(file.path(output_dir, "sensitivity_plot.png"), sensitivity_plot,
       width = 9, height = 5, dpi = 300)

## ---- 6. BOOTSTRAP RESAMPLING FOR CONFIDENCE INTERVALS ----------------------

## 6a. Bootstrap CI for the PROPORTION of individuals in each category
##     (uses the ORIGINAL, un-scaled thresholds)

bootstrap_category_ci = function(species_name, n_boot = 1000, ci_level = 0.95) {
  sp_params = params %>% filter(species_common == species_name)
  fl_values = data_clean %>%
    filter(species_common == species_name) %>%
    pull(fl)
  
  n_obs = length(fl_values)
  categories = c("Below_L50", "L50_to_Lopt", "At_Lopt", "MS_megaspawner")
  
  boot_props = matrix(NA_real_, nrow = n_boot, ncol = length(categories),
                       dimnames = list(NULL, categories))
  
  for (i in seq_len(n_boot)) {
    resampled = sample(fl_values, size = n_obs, replace = TRUE)
    cats = categorize_lengths(resampled, sp_params$l50, sp_params$lopt, sp_params$ms)
    tab = table(factor(cats, levels = categories))
    boot_props[i, ] = as.numeric(tab) / n_obs
  }
  
  alpha = 1 - ci_level
  ci_summary = tibble(
    species_common = species_name,
    category = categories,
    point_estimate = colMeans(boot_props),
    ci_lower = apply(boot_props, 2, quantile, probs = alpha / 2),
    ci_upper = apply(boot_props, 2, quantile, probs = 1 - alpha / 2)
  )
  
  ci_summary
}

bootstrap_ci_results = map_dfr(
  params$species_common,
  ~bootstrap_category_ci(.x, n_boot = n_boot, ci_level = ci_level)
)

cat(sprintf("\n--- Bootstrap %.0f%% CIs for category proportions (n_boot = %d) ---\n",
            ci_level * 100, n_boot))
print(bootstrap_ci_results)

write_csv(bootstrap_ci_results, file.path(output_dir, "bootstrap_category_ci.csv"))

# Plot: proportion estimates with bootstrap CI error bars
bootstrap_ci_plot = ggplot(
  bootstrap_ci_results,
  aes(x = category, y = point_estimate, fill = category)
) +
  geom_col() +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2) +
  facet_wrap(~species_common) +
  scale_y_continuous(labels = percent_format()) +
  labs(
    title = sprintf("Category proportions with %.0f%% bootstrap CIs (n_boot = %d)",
                    ci_level * 100, n_boot),
    x = "Length category", y = "Proportion of individuals"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(file.path(output_dir, "bootstrap_category_ci_plot.png"), bootstrap_ci_plot,
       width = 9, height = 5, dpi = 300)

## 6b. Bootstrap CI band around the length-frequency HISTOGRAM/DENSITY

bootstrap_length_frequency = function(species_name, n_boot = 1000,
                                       ci_level = 0.95, bin_width = 2) {
  fl_values = data_clean %>%
    filter(species_common == species_name) %>%
    pull(fl)
  
  n_obs = length(fl_values)
  breaks = seq(floor(min(fl_values) / bin_width) * bin_width,
                ceiling(max(fl_values) / bin_width) * bin_width,
                by = bin_width)
  mids = head(breaks, -1) + bin_width / 2
  
  boot_densities = matrix(NA_real_, nrow = n_boot, ncol = length(mids))
  
  for (i in seq_len(n_boot)) {
    resampled = sample(fl_values, size = n_obs, replace = TRUE)
    h = hist(resampled, breaks = breaks, plot = FALSE)
    boot_densities[i, ] = h$counts / n_obs
  }
  
  alpha = 1 - ci_level
  tibble(
    species_common = species_name,
    fl_mid = mids,
    freq_mean = colMeans(boot_densities),
    freq_lower = apply(boot_densities, 2, quantile, probs = alpha / 2),
    freq_upper = apply(boot_densities, 2, quantile, probs = 1 - alpha / 2)
  )
}

length_freq_ci = map_dfr(
  params$species_common,
  ~bootstrap_length_frequency(.x, n_boot = n_boot, ci_level = ci_level, bin_width = 2)
)

write_csv(length_freq_ci, file.path(output_dir, "bootstrap_length_frequency_ci.csv"))

# Plot: length-frequency histogram with bootstrap CI ribbon and threshold lines
threshold_lines = params %>%
  select(species_common, l50, lopt, ms) %>%
  pivot_longer(cols = c(l50, lopt, ms), names_to = "threshold", values_to = "value")

length_freq_plot = ggplot(length_freq_ci, aes(x = fl_mid, y = freq_mean)) +
  geom_ribbon(aes(ymin = freq_lower, ymax = freq_upper), fill = "steelblue", alpha = 0.3) +
  geom_line(color = "steelblue", linewidth = 0.8) +
  geom_vline(
    data = threshold_lines,
    aes(xintercept = value, linetype = threshold),
    color = "firebrick"
  ) +
  facet_wrap(~species_common, scales = "free_x") +
  labs(
    title = sprintf("Length-frequency distribution with %.0f%% bootstrap CI band", ci_level * 100),
    subtitle = "Vertical lines: L50, Lopt, MS (mega-spawner) thresholds",
    x = "Fork length (cm)", y = "Relative frequency",
    linetype = "Threshold"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(file.path(output_dir, "length_frequency_ci_plot.png"), length_freq_plot,
       width = 10, height = 5, dpi = 300)

## ---- 7. INDIVIDUAL-LEVEL OUTPUT (each fish tagged with its category) -------

write_csv(
  data_categorized %>% select(species, species_common, fl, l50, lopt, ms, category),
  file.path(output_dir, "individual_length_categories.csv")
)

## ---- 8. WRAP-UP -------------------------------------------------------------

cat("\n=============================================================\n")
cat("Analysis complete. Outputs written to:", normalizePath(output_dir), "\n")
cat("  - base_category_summary.csv          : base-case category proportions\n")
cat("  - sensitivity_analysis.csv           : proportions under +/-10%/20% thresholds\n")
cat("  - sensitivity_plot.png               : stacked-bar sensitivity plot\n")
cat("  - bootstrap_category_ci.csv          : bootstrap 95% CIs for proportions\n")
cat("  - bootstrap_category_ci_plot.png     : bar plot with bootstrap CI error bars\n")
cat("  - bootstrap_length_frequency_ci.csv  : bootstrap CI band for length-freq curve\n")
cat("  - length_frequency_ci_plot.png       : length-frequency plot with CI ribbon\n")
cat("  - individual_length_categories.csv   : per-individual category assignment\n")
cat("=============================================================\n")