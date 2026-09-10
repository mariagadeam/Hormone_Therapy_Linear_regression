# Main Project: UrinaryMicrobiota (KIJURI)
# Author: María Gadea Martínez
# Sub-project: Impact of HT (hormone therapy) in the urinary microbiome of 
# peri- or menopausal women
# Analysis of continuous duration effect using natural cubic splines
# Initial date: 2026-05-23
# Last update: 2026-05-25


# ---------------------------------------------------------------------------
## 1. Load libraries
library(rio)
library(dplyr)
library(ggplot2)
library(splines)
library(scales)
library(openxlsx)
library(tidyr)

# ---------------------------------------------------------------------------
## 2. Set working directory and load data (if already run once)
setwd("C:/Users/riade/OneDrive - Uppsala universitet/Thesis project TFM")
load("Rdata/splines.RData")

taxonomies <- rio::import("Data/taxonomy_students_females.csv")
# Covariates for cases and controls
cases_covariates  <- rio::import("Data/Cases_HRT.csv",    dec = ",")
controls_covariates <- rio::import("Data/Controls_HRT.csv", dec = ",")
match_covariates  <- bind_rows(cases_covariates, controls_covariates)
rownames(match_covariates) <- match_covariates[, 1]
match_covariates  <- match_covariates[, -1]

# ---------------------------------------------------------------------------
## 3. Helper: prepare abundance data (collapse to genus/family if needed)
prepare_data <- function(filepath, level, type) {
  abundance_raw <- rio::import(filepath, dec = ",")
  rownames(abundance_raw) <- abundance_raw[, 1]
  abundance_out <- abundance_raw[, -1]
  if ((level == "genus" || level == "family") && type == "RA") {
    lookup <- setNames(taxonomies[[level]], taxonomies$MGS)
    colnames(abundance_out) <- ifelse(
      colnames(abundance_out) %in% names(lookup),
      lookup[colnames(abundance_out)],
      colnames(abundance_out)
    )
    abundance_out <- as.data.frame(
      rowsum(t(abundance_out), group = sub("\\..*$", "", colnames(abundance_out)))
    )
    abundance_out <- as.data.frame(t(abundance_out))
  }
  return(abundance_out)
}
# ---------------------------------------------------------------------------
## 4. Helper function to build combined dataset with covariates and duration
build_combined <- function(taxa) {
  
  common <- intersect(rownames(match_covariates), rownames(taxa))
  combined <- cbind(match_covariates[common, ], taxa[common, ])
  combined$case_type[is.na(combined$case_type)] <- "control"
  combined$continuous_days <- as.numeric(combined$continuous_days)
  # Remove cases that are missing the duration predictor
  combined <- combined[!(combined$case_type != "control" & is.na(combined$continuous_days)), ]
  # Controls get 0 days
  combined$continuous_days[combined$case_type == "control"] <- 0
  return(combined)
}
# ---------------------------------------------------------------------------

## 5. Function of spline models 

run_spline_all_taxa <- function(taxa, level, type) {
  combined_data <- build_combined(taxa)
  covariates <- c("ageatvisitone", "DMcat", "bmi", "derived_smoke_status",
                  "physact", "alcohol", "placebirth", "education")
  # Custom knots (on cases)
  cases_days  <- combined_data$continuous_days[combined_data$continuous_days > 0]
  custom_knots <- quantile(cases_days, probs = c(0.33, 0.67), na.rm = TRUE)
  
  # 5.1 Model formulas
  
  # Adjusted: spline + batch + covariates
  fixed_part_adj   <- paste(
    c("ns(continuous_days, knots = custom_knots)", "BATCH", covariates),
    collapse = " + "
  )
  # Basic: spline + batch only
  fixed_part_basic <- "ns(continuous_days, knots = custom_knots) + BATCH"
  
  # 5.2 Reference for prediction (median of each covariate)
  
  ref_data_base <- combined_data[1, , drop = FALSE]
  for (col in c("BATCH", covariates)) {
    if (is.numeric(combined_data[[col]])) {
      ref_data_base[[col]] <- median(combined_data[[col]], na.rm = TRUE)
    } else {
      mode_val <- names(sort(table(combined_data[[col]]), decreasing = TRUE))[1]
      ref_data_base[[col]] <- combined_data[[col]][
        which(as.character(combined_data[[col]]) == mode_val)[1]
      ]
    }
  }
  # 5.3 Prediction grid: log-spaced from 0 to max days
  pred_days <- unique(
    c(0, exp(seq(log(1), log(max(combined_data$continuous_days) + 1), length.out = 100)) - 1)
  )
  # Containers
  all_preds  <- data.frame()
  stats_rows <- list()
  taxa_cols <- colnames(taxa)
  
  # 5.4 Loop for each taxon
  for (taxon in taxa_cols) {
    if (taxon %in% c("N", "Other")) next
    # Adjusted model
    formula_adj <- as.formula(paste0("`", taxon, "` ~ ", fixed_part_adj))
    model_adj   <- tryCatch(
      lm(formula_adj, data = combined_data),
      error = function(e) NULL
    )
    # Basic model
    formula_basic <- as.formula(paste0("`", taxon, "` ~ ", fixed_part_basic))
    model_basic   <- tryCatch(
      lm(formula_basic, data = combined_data),
      error = function(e) NULL
    )
    # Extract spline p-value via F-test 
    extract_spline_pval <- function(model_full, has_covariates) {
      
      if (is.null(model_full)) return(NA_real_)
      # Build reduced formula (remove the ns() term)
      if (has_covariates) {
        reduced_formula <- as.formula(
          paste0("`", taxon, "` ~ BATCH + ",
                 paste(covariates, collapse = " + "))
        )
      } else {
        reduced_formula <- as.formula(paste0("`", taxon, "` ~ BATCH"))
      }
      model_reduced <- tryCatch(
        lm(reduced_formula, data = combined_data),
        error = function(e) NULL
      )
      if (is.null(model_reduced)) return(NA_real_)
      a <- tryCatch(
        anova(model_reduced, model_full),
        error = function(e) NULL
      )
      if (is.null(a) || nrow(a) < 2) return(NA_real_)
      return(a[2, "Pr(>F)"])
    }
    pval_adj   <- extract_spline_pval(model_adj,   has_covariates = TRUE)
    pval_basic <- extract_spline_pval(model_basic, has_covariates = FALSE)
    
    # Record stats
    stats_rows[[length(stats_rows) + 1]] <- data.frame(
      taxon_id   = taxon,
      level      = level,
      type       = type,
      adjustment = "adjusted",
      pvalue     = pval_adj,
      stringsAsFactors = FALSE
    )
    stats_rows[[length(stats_rows) + 1]] <- data.frame(
      taxon_id   = taxon,
      level      = level,
      type       = type,
      adjustment = "basic",
      pvalue     = pval_basic,
      stringsAsFactors = FALSE
    )
    # Generate predictions
    for (mod_info in list(
      list(model = model_adj,   adj_label = "adjusted"),
      list(model = model_basic, adj_label = "basic")
    )) {
      if (is.null(mod_info$model)) next
      ref_data <- ref_data_base[rep(1, length(pred_days)), ]
      ref_data$continuous_days <- pred_days
      predicted <- tryCatch(
        predict(mod_info$model, newdata = ref_data, se.fit = TRUE),
        error = function(e) NULL
      )
      if (is.null(predicted)) next
      tmp <- ref_data[, c("continuous_days"), drop = FALSE]
      tmp$predicted  <- predicted$fit
      tmp$conf.low   <- predicted$fit - 1.96 * predicted$se.fit
      tmp$conf.high  <- predicted$fit + 1.96 * predicted$se.fit
      tmp$taxon_id   <- taxon
      tmp$adjustment <- mod_info$adj_label
      all_preds <- rbind(all_preds, tmp)
    }
  }
  
  # Combine stats into one data frame
  stats_df <- bind_rows(stats_rows)
  
  # 5.5 Resolve display labels
  if (level == "species") {
    label_map <- setNames(taxonomies$species, taxonomies$MGS)
    stats_df$display_label <- ifelse(
      stats_df$taxon_id %in% names(label_map),
      label_map[stats_df$taxon_id],
      stats_df$taxon_id
    )
    # Filter out "Other" genus
    other_mgs <- taxonomies$MGS[taxonomies$genus == "Other"]
    stats_df <- stats_df[!stats_df$taxon_id %in% other_mgs, ]
    all_preds$plot_label <- ifelse(
      all_preds$taxon_id %in% names(label_map),
      label_map[all_preds$taxon_id],
      all_preds$taxon_id
    )
    # Filter predictions too
    all_preds <- all_preds[!all_preds$taxon_id %in% other_mgs, ]
  } else {
    stats_df$display_label <- stats_df$taxon_id
    all_preds$plot_label   <- all_preds$taxon_id
  }
  return(list(
    stats   = stats_df,
    preds   = all_preds,
    knots   = custom_knots,
    unit    = ifelse(type == "CLR", "CLR Value", "Relative Abundance")
  ))
}
# ---------------------------------------------------------------------------

## 6. FDR correction 
apply_fdr <- function(all_stats_df) {
  # Separate by adjustment type then apply FDR within each group
  all_stats_df <- all_stats_df %>%
    group_by(adjustment) %>%
    mutate(qvalue = p.adjust(pvalue, method = "fdr")) %>%
    ungroup()
  return(all_stats_df)
}
# ---------------------------------------------------------------------------

## 7. Plotting and save SVG 

make_spline_plot <- function(preds, sig_taxa, level, type, adjustment, unit_str, output_prefix) {
  preds_sub <- preds %>% filter(adjustment == !!adjustment)
  # 7.1 Plot all taxa
  p_all <- ggplot(preds_sub, aes(x = continuous_days, y = predicted)) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high), fill = "firebrick", alpha = 0.2) +
    geom_line(color = "firebrick", linewidth = 0.8) +
    facet_wrap(~plot_label, scales = "free_y") +
    scale_x_continuous(
      trans  = "pseudo_log",
      breaks = c(0, 7, 30, 90, 365, 1460, 3650),
      labels = c("Control", "1wk", "1mo", "3mo", "1yr", "4yrs", "10yrs")
    ) +
    theme_bw(base_size = 9) +
    theme(
      axis.text.x     = element_text(angle = 45, hjust = 1),
      strip.text       = element_text(size = 7),
      panel.spacing    = unit(0.3, "lines")
    ) +
    labs(
      title    = paste(tools::toTitleCase(level), type,
                       paste0("(", tools::toTitleCase(adjustment), ")"),
                       "- Natural Spline Predictions (All Taxa)"),
      subtitle = "Non-linear time effect while controlling for clinical covariates",
      x        = "Duration (Days / Log-Scaled)",
      y        = paste("Predicted", unit_str)
    )
  
  n_taxa <- length(unique(preds_sub$plot_label))
  n_cols <- 4
  n_rows <- ceiling(n_taxa / n_cols)
  w      <- max(10, n_cols * 3)
  h      <- max(6,  n_rows * 2.5)
  ggsave(
    paste0(output_prefix, "_ALL_", adjustment, ".svg"),
    p_all, width = w, height = h, device = "svg"
  )
  message(sprintf("  [SVG] All taxa (%s, %s, %s): %d taxa -> %s",
                  level, type, adjustment, n_taxa,
                  paste0(output_prefix, "_ALL_", adjustment, ".svg")))
  
  # 7.2 Plot only significant taxa (q < 0.05)
  if (length(sig_taxa) > 0) {
    preds_sig <- preds_sub %>% filter(plot_label %in% sig_taxa)
    if (nrow(preds_sig) > 0) {
      n_sig  <- length(unique(preds_sig$plot_label))
      n_cols_s <- min(2, n_sig)
      n_rows_s <- ceiling(n_sig / n_cols_s)
      w_s    <- max(6, n_cols_s * 3)
      h_s    <- max(4, n_rows_s * 2.5)
      p_sig <- ggplot(preds_sig, aes(x = continuous_days, y = predicted)) +
        geom_ribbon(aes(ymin = conf.low, ymax = conf.high), fill = "steelblue", alpha = 0.25) +
        geom_line(color = "steelblue", linewidth = 1) +
        facet_wrap(~plot_label, scales = "free_y", ncol = 2) +
        scale_x_continuous(
          trans  = "pseudo_log",
          breaks = c(0, 7, 30, 90, 365, 1460, 3650),
          labels = c("Control", "1wk", "1mo", "3mo", "1yr", "4yrs", "10yrs")
        ) +
        theme_bw(base_size = 10) +
        theme(
          axis.text.x  = element_text(angle = 45, hjust = 1),
          strip.text   = element_text(size = 9, face = "italic")
        ) +
        labs(
          title    = paste(tools::toTitleCase(level), type,
                           paste0("(", tools::toTitleCase(adjustment), ")"),
                           "- Significant Spline Taxa (q < 0.05)"),
          subtitle = "Only taxa with FDR-corrected q-value < 0.05 shown",
          x        = "Duration (Days / Log-Scaled)",
          y        = paste("Predicted", unit_str)
        )
      ggsave(
        paste0(output_prefix, "_SIG_", adjustment, ".svg"),
        p_sig, width = w_s, height = h_s, device = "svg"
      )
      message(sprintf("  [SVG] Significant taxa (%s, %s, %s): %d taxa -> %s",
                      level, type, adjustment, n_sig,
                      paste0(output_prefix, "_SIG_", adjustment, ".svg")))
    }
  } else {
    message(sprintf("  [INFO] No significant taxa (q < 0.05) for %s %s %s", level, type, adjustment))
  }
}
# ---------------------------------------------------------------------------
## 8. Main execution loop

# 8.1 Define configurations for all 6 analyses
datasets <- list(
  list(file = "Data/Cohort_data/Female891_CLR_species_theses.csv", level = "species", type = "CLR"),
  list(file = "Data/Cohort_data/fem_species_relabundance.csv",     level = "species", type = "RA"),
  list(file = "Data/Cohort_data/Female891_GENUS_CLR_theses.csv",               level = "genus",   type = "CLR"),
  list(file = "Data/Cohort_data/fem_genus_relabundance.csv",       level = "genus",   type = "RA"),
  list(file = "Data/Cohort_data/Female891_FAM_CLR_theses.csv",              level = "family",  type = "CLR"),
  list(file = "Data/Cohort_data/fem_family_relabundance.csv",      level = "family",  type = "RA")
)

# Collect  stats 
all_stats_list <- list()
all_preds_list <- list()
analysis_meta  <- list()

# 8.2 Run spline models for all datasets 
for (i in seq_along(datasets)) {
  analysis <- datasets[[i]]
  message(sprintf("\n--- Dataset %d/%d: %s %s ---", i, length(datasets), analysis$level, analysis$type))
  taxa_data <- prepare_data(analysis$file, analysis$level, analysis$type)
  message(sprintf("  Loaded %d taxa x %d samples", ncol(taxa_data), nrow(taxa_data)))
  res <- run_spline_all_taxa(taxa_data, analysis$level, analysis$type)
  all_stats_list[[i]] <- res$stats
  all_preds_list[[i]] <- res$preds
  analysis_meta[[i]]  <- list(
    level  = analysis$level,
    type   = analysis$type,
    unit   = res$unit,
    prefix = paste0("Images/spline_", analysis$level, "_", tolower(analysis$type))
  )
}
# 8.3 FDR correction 
all_stats_df <- bind_rows(all_stats_list)
all_stats_df <- apply_fdr(all_stats_df)

# Significant taxa lookup (q < 0.05)
sig_lookup <- all_stats_df %>%
  filter(qvalue < 0.05) %>%
  select(taxon_id, level, type, adjustment, display_label)

# 8.4 Generate plots 
for (i in seq_along(datasets)) {
  meta  <- analysis_meta[[i]]
  preds <- all_preds_list[[i]]
  for (adj in c("adjusted", "basic")) {
    # Get significant display labels for this combination
    sig_labels <- sig_lookup %>%
      filter(level == meta$level, type == meta$type, adjustment == adj) %>%
      pull(display_label) %>%
      unique()
    make_spline_plot(
      preds        = preds,
      sig_taxa     = sig_labels,
      level        = meta$level,
      type         = meta$type,
      adjustment   = adj,
      unit_str     = meta$unit,
      output_prefix = meta$prefix
    )
  }
}

# ---------------------------------------------------------------------------

## 9. Save workspace

save.image(file = "Rdata/splines.RData")