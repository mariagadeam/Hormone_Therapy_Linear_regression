# Main Project: UrinaryMicrobiota (KIJURI)
# Author: María Gadea Martínez
# Sub-project: Impact of HT (hormone therapy) in the urinary microbiome of
# peri- or menopausal women
# Analysis of duration of therapy (recent users, continuous users,& unexposed)
# Initial date: 2026-03-24
# Last update: 2026-05-21

# ---------------------------------------------------------------------------
## 1. Load libraries

install.packages("svglite")
library(svglite)
library(rio)
library(dplyr)
library(ggplot2)
library(emmeans)
library(tibble)


# ---------------------------------------------------------------------------
## 2. Set working directory

setwd("C:/Users/riade/OneDrive - Uppsala universitet/Thesis project TFM")
load("Rdata/duration_B_lm_rdata.RData")

# ---------------------------------------------------------------------------
## 3. Load data

# Combined covariates file (cases + controls)
cases_covariates <- rio::import("Data/Cases_HRT.csv", dec = ",")
controls_covariates <- rio::import("Data/Controls_HRT.csv", dec = ",")
all_covariates <- bind_rows(cases_covariates, controls_covariates)
rownames(all_covariates) <- all_covariates[, 1]
all_covariates <- all_covariates[, -1]

# Taxonomy table
taxonomies <- rio::import("Data/taxonomy_students_females.csv")

# ---------------------------------------------------------------------------
## 4. Function to prepare abundance (CLR or relative abundance) dataset

prepare_data <- function(filepath, level, type) {
  
  #  4.1 Load specific dataset
  abundance_raw <- rio::import(filepath, dec = ",")
  rownames(abundance_raw) <- abundance_raw[, 1]
  abundance_out <- abundance_raw[, -1]

  #  4.2 Collapse dataset if it's family or genus level
  if ((level == "genus" || level == "family") && type == "Relab") {
    lookup <- setNames(taxonomies[[level]], taxonomies$MGS)

    colnames(abundance_out) <- ifelse(
      colnames(abundance_out) %in% names(lookup),
      lookup[colnames(abundance_out)],
      colnames(abundance_out)
    )

    abundance_out <- as.data.frame(
      rowsum(
        t(abundance_out),
        group = sub("\\..*$", "", colnames(abundance_out))
      )
    )
    abundance_out <- as.data.frame(t(abundance_out))
  }

  return(abundance_out)
}

# ---------------------------------------------------------------------------
## 5. Main analysis function (plain linear regression, no subclass term)

run_duration_model <- function(taxa, level, type, output_prefix) {

  # 5.1 Merge covariates, duration data, and abundance

  common <- intersect(rownames(all_covariates), rownames(taxa))

  combined_data <- cbind(
    all_covariates[common, , drop = FALSE],
    taxa[common, , drop = FALSE]
  )


  # 5.2 Specify exposure groups

  combined_data$case_type[is.na(combined_data$case_type)] <- "control"
  combined_data$case_type <- factor(
    combined_data$case_type,
    levels = c("control", "recent", "continuous")
  )

  # 5.3 Covariate lists

  covariates <- c("ageatvisitone", "DMcat", "bmi", "derived_smoke_status",
                  "physact", "alcohol", "placebirth", "education")

  # 5.4 Results containers

  results <- list(
    adj   = list(recent = data.frame(), cont = data.frame(), duration = data.frame()),
    basic = list(recent = data.frame(), cont = data.frame(), duration = data.frame())
  )

  # 5.5 Function to extract coefficients and continuous-vs-recent contrast

  extract_stats <- function(model_fit, taxon, nuis_vars) {
    coef_table <- summary(model_fit)$coefficients

    get_coef <- function(term) {
      res <- coef_table[term, ]
      data.frame(
        taxon_id = taxon,
        estimate = res["Estimate"],
        se       = res["Std. Error"],
        p        = res["Pr(>|t|)"]
      )
    }

    emm          <- emmeans(model_fit, ~ case_type, nuisance = nuis_vars)
    contrast_res <- pairs(emm, reverse = TRUE)

    cont_vs_recent <- as.data.frame(contrast_res) %>%
      filter(grepl("continuous.*recent|recent.*continuous", contrast))

    results_duration <- data.frame()
    if (nrow(cont_vs_recent) > 0) {
      est_val <- cont_vs_recent$estimate
      if (grepl("recent.*continuous", cont_vs_recent$contrast[1])) {
        est_val <- -est_val
      }
      results_duration <- data.frame(
        taxon_id = taxon,
        estimate = est_val,
        se       = cont_vs_recent$SE,
        p        = cont_vs_recent$p.value
      )
    }

    list(
      recent   = get_coef("case_typerecent"),
      cont     = get_coef("case_typecontinuous"),
      duration = results_duration
    )
  }

  # 5.6 Loop through taxon

  taxa_cols <- colnames(taxa)

  for (taxon in taxa_cols) {
    if (taxon %in% c("N", "Other") && level != "species") next

    fixed_part_adj   <- paste(c("case_type", "BATCH", covariates), collapse = " + ")
    fixed_part_basic <- "case_type + BATCH"

    # Adjusted model: case_type + batch + covariates
    formula_adj   <- as.formula(paste0("`", taxon, "` ~ ", fixed_part_adj))
    model_adj     <- lm(formula_adj,   data = combined_data)
    nuis_adj      <- c("BATCH", covariates)

    # Basic model: case_type + batch only
    formula_basic <- as.formula(paste0("`", taxon, "` ~ ", fixed_part_basic))
    model_basic   <- lm(formula_basic, data = combined_data)
    nuis_basic    <- c("BATCH")

    stats_adj   <- extract_stats(model_adj,   taxon, nuis_adj)
    stats_basic <- extract_stats(model_basic, taxon, nuis_basic)

    for (comp in c("recent", "cont", "duration")) {
      results$adj[[comp]]   <- rbind(results$adj[[comp]],   stats_adj[[comp]])
      results$basic[[comp]] <- rbind(results$basic[[comp]], stats_basic[[comp]])
    }
  }

  # 5.7 FDR correction

  res_list <- list(
    adj_recent   = results$adj$recent,
    adj_cont     = results$adj$cont,
    adj_dur      = results$adj$duration,
    basic_recent = results$basic$recent,
    basic_cont   = results$basic$cont,
    basic_dur    = results$basic$duration
  )
  res_list <- lapply(res_list, function(df) {
    df$p_adj <- p.adjust(df$p, method = "fdr")
    df
  })

  results_recent         <- res_list$adj_recent
  results_cont           <- res_list$adj_cont
  results_duration       <- res_list$adj_dur
  results_recent_basic   <- res_list$basic_recent
  results_cont_basic     <- res_list$basic_cont
  results_duration_basic <- res_list$basic_dur

  # 5.8 Plotting

  # Prepare the data for plotting and significance labels
  prepare_plot_data <- function(regression_results) {
    results_plot <- regression_results %>%
      mutate(
        lower = estimate - 1.96 * se,
        upper = estimate + 1.96 * se,
        significant = case_when(
          p_adj < 0.05                    ~ "Significant (q < 0.05)",
          p < 0.05 & p_adj >= 0.05        ~ "Non-significant after multiple testing adjustment (p < 0.05, q >= 0.05)",
          TRUE                            ~ "Non-significant (p > 0.05)"
        )
      )

    if (level == "species") {
      results_plot <- results_plot %>%
        left_join(taxonomies, by = c("taxon_id" = "MGS")) %>%
        filter(genus != "Other") %>%
        mutate(plot_label = species)
    } else {
      results_plot <- results_plot %>% mutate(plot_label = taxon_id)
    }
    results_plot
  }

  # Plot using ggplot, add significance colors
  plot_forest <- function(plot_data, title, xlabel) {
    sig_colors <- c(
      "Non-significant (p > 0.05)"                                              = "red",
      "Non-significant after multiple testing adjustment (p < 0.05, q >= 0.05)" = "green",
      "Significant (q < 0.05)"                                                  = "blue"
    )
    ggplot(plot_data, aes(x = estimate, y = plot_label)) +
      geom_point(aes(color = significant), size = 3) +
      geom_errorbar(aes(xmin = lower, xmax = upper), orientation = "y", height = 0.2) +
      geom_vline(xintercept = 0, linetype = "dashed") +
      scale_color_manual(name = "Significance", values = sig_colors) +
      labs(x = xlabel, y = tools::toTitleCase(level), title = title) +
      theme(legend.position = "none")
  }

  unit_str <- ifelse(type == "CLR", "CLR", "Relative Abundance")

  model_results <- list(
    list(results = list(recent = results_recent,       cont = results_cont,       dur = results_duration),
         label = "Adjusted", suffix = "adj"),
    list(results = list(recent = results_recent_basic, cont = results_cont_basic, dur = results_duration_basic),
         label = "Basic",    suffix = "basic")
  )

  # Three plots from each analysis: recents vs control, continuous vs control,
  # recents vs continuous (duration effect)
  
  for (mod_res in model_results) {
    comparisons <- list(
      list(res = mod_res$results$recent, suffix = "recent_vs_control", title = "Control vs Recent"),
      list(res = mod_res$results$cont,   suffix = "cont_vs_control",   title = "Control vs Continuous"),
      list(res = mod_res$results$dur,    suffix = "duration_effect",   title = "Recent vs Continuous")
    )
    for (comp in comparisons) {
      plot_df <- prepare_plot_data(comp$res)
      p <- plot_forest(
        plot_df,
        paste(tools::toTitleCase(level), type, comp$title,
              paste0("(Linear Regression ", mod_res$label, ")")),
        paste("Difference in", unit_str, "(", comp$title, ")")
      )
      ggsave(paste0(output_prefix, "_", mod_res$suffix, "_", comp$suffix, ".svg"),
             p, width = 12, height = 7)
    }
  }

  list(
    results_recent         = results_recent,
    results_cont           = results_cont,
    results_duration       = results_duration,
    results_recent_basic   = results_recent_basic,
    results_cont_basic     = results_cont_basic,
    results_duration_basic = results_duration_basic
  )
}

# ---------------------------------------------------------------------------
## 6. Execute analyses for all datasets

 # Define configurations for all 6 analyses
dataset <- list(
  list(file = "Data/Cohort_data/Female891_CLR_species_theses.csv",  level = "species", type = "CLR"),
  list(file = "Data/Cohort_data/fem_species_relabundance.csv",       level = "species", type = "Relab"),
  list(file = "Data/Cohort_data/Female891_GENUS_CLR_theses.csv",                  level = "genus",   type = "CLR"),
  list(file = "Data/Cohort_data/fem_genus_relabundance.csv",         level = "genus",   type = "Relab"),
  list(file = "Data/Cohort_data/Female891_FAM_CLR_theses.csv",                 level = "family",  type = "CLR"),
  list(file = "Data/Cohort_data/fem_family_relabundance.csv",        level = "family",  type = "Relab")
)

all_results <- list()

 #  Execution block, iterate through the configurations and run each analysis
for (analysis in dataset) {
  taxa_data <- prepare_data(analysis$file, analysis$level, analysis$type)

  prefix <- paste0("forest_", analysis$level, "_", tolower(analysis$type), "_lm")
  res    <- run_duration_model(taxa_data, analysis$level, analysis$type, prefix)
  all_results[[paste(analysis$level, analysis$type, sep = "_")]] <- res
}

# ---------------------------------------------------------------------------
## 7. Save workspace

save.image(file = "Rdata/duration_B_lm_rdata.RData")

