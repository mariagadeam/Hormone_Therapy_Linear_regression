# Main Project: UrinaryMicrobiota (KIJURI)
# Author: María Gadea Martínez
# Sub-project: Impact of HT (hormone therapy) in the urinary microbiome of 
# peri- or menopausal women
# Analysis of duration of therapy adjusted with fixed and mixed effects 
# for matching at three taxonomic ranks and with CLR and relab data
# Initial date: 2026-03-24
# Last update: 2026-03-25


# ---------------------------------------------------------------------------
## 1. Load the libraries

library(rio)
library(dplyr)
library(ggplot2)
library(lme4)
library(lmerTest)
library(emmeans)
library(tibble)
library(pbkrtest)

# ---------------------------------------------------------------------------
## 2. Set working directory and load common data (set sample number to row name)

setwd("C:/Users/riade/OneDrive - Uppsala universitet/Thesis project TFM")
load("Rdata/duration_all_analyses_rdata.RData")

# Duration time 

cases_duration <- rio::import("Data/FEMALE_Cases_Details.csv", dec = ",")
rownames(cases_duration) <- cases_duration[, 1]
cases_duration <- cases_duration[, -1]

# Taxonomies

taxonomies <- rio::import("Data/taxonomy_students_females.csv")

# Metadata / covariates

match_covariates <- rio::import("Data/FEMALE_Covariates_Cases_Control.csv",
                                dec = ",")
rownames(match_covariates) <- match_covariates[, 1]
match_covariates <- match_covariates[, -1]


# --------------------------------------------------------------------------
## 4. Function to prepare abundance (CLR or relative abundance) dataset
# If genus or family level, it needs to be collapsed

prepare_data <- function(filepath, level, type) {
  # Load dataset
  abundance_raw <- rio::import(filepath, dec = ",")
  rownames(abundance_raw) <- abundance_raw[, 1]
  abundance_out <- abundance_raw[, -1]
  
  # Collapse if genus or family
  if (level == "genus" || level == "family") {
    lookup <- setNames(taxonomies[[level]], taxonomies$MGS)
    
    # Replace column names with corresponding taxonomy names where possible
    colnames(abundance_out) <- ifelse(
      colnames(abundance_out) %in% names(lookup),
      lookup[colnames(abundance_out)],
      colnames(abundance_out)
    )
    
    # Aggregate by identical column names (same genera / family)
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

# --------------------------------------------------------------------------
# 5. Main Analysis Function

run_duration_model <- function(taxa, level, type, model_type, output_prefix) {
  
  # Merge metadata and duration data
  common <- intersect(rownames(match_covariates), rownames(taxa))
  combined_data <- cbind(match_covariates[common, ], taxa[common, ])
  
  cases_filtered <- cases_duration[rownames(combined_data), ]
  combined_data <- cbind(combined_data, cases_filtered)
  
  # Change empty exposure duration variable with 'control' and specify it as 
  # the baseline / reference variable for regressions (by defining the order)
  combined_data$case_type[is.na(combined_data$case_type)] <- "control"
  combined_data$case_type <- factor(
    combined_data$case_type,
    levels = c("control", "recent", "continuous")
  )
  
  # Subclass (pair identifier) in this case is formatted as e.g 2_1, 
  # so transformation to factor is not needed. But in the case of the subclass
  # being a continuous number, this is needed in order to treat each pair as 
  # a distinct categorical dummy variable
  if (model_type == "fixed") {
    combined_data$subclass <- as.factor(combined_data$subclass)
  }
  
  # Define formula components
  covariates <- c("ageatvisitone", "DMcat", "bmi", "derived_smoke_status", 
                  "physact", "alcohol", "placebirth", "education")
  
  # Results containers
  results <- list(
    adj = list(recent = data.frame(), cont = data.frame(), duration = data.frame()),
    basic = list(recent = data.frame(), cont = data.frame(), duration = data.frame())
  )
  
  # Helper function to extract coefficients and contrasts
  extract_stats <- function(model_fit, taxon, nuis_vars) {
    coef_table <- summary(model_fit)$coefficients
    
    # Simple comparisons from coefficients
    get_coef <- function(term) {
        res <- coef_table[term, ]
        return(data.frame(taxon_id = taxon, estimate = res["Estimate"], 
                          se = res["Std. Error"], p = res["Pr(>|t|)"]))
    }
    
    # Continuous vs Recent (emmeans)
    # emm <- emmeans(model_fit, ~ case_type, weights = "proportional")
    # Passed nuis_vars to the nuisance argument so emmeans skips building a massive grid
    emm <- emmeans(model_fit, ~ case_type, 
                   nuisance = nuis_vars)
    # we can run emmeans either by increasing the range limit or by using
    # nuisance variables
    # emm <- emmeans(model_fit, ~ case_type, rg.limit = 4000000)
    contrast_res <- pairs(emm, reverse = TRUE)
    
    # For the directionality of the analysis (emmeans can switch it based on
    # alphabetical order)
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
        se = cont_vs_recent$SE,
        p = cont_vs_recent$p.value
      )
    }

    list(recent = get_coef("case_typerecent"), 
         cont = get_coef("case_typecontinuous"), 
         duration = results_duration)
  }

  taxa_cols <- colnames(taxa)
  
  for (taxon in taxa_cols) {
    # Skip "Other" or "N" groups to save time and because they are not plotted
    if (taxon %in% c("N", "Other") && level != "species") next 
    
    # Define models
    fixed_part_adj <- paste(c("case_type", "BATCH", covariates), collapse = " + ")
    fixed_part_basic <- "case_type + BATCH"
    
    # Mixed model
    if (model_type == "mixed") {
      # Adjusted model (covariates)
      formula_adj <- as.formula(paste0("`", taxon, "` ~ ", fixed_part_adj, " + (1 | subclass)"))
      model_adj <- lmer(formula_adj, data = combined_data)
      # Basic model (no covariates)
      formula_basic <- as.formula(paste0("`", taxon, "` ~ ", fixed_part_basic, " + (1 | subclass)"))
      model_basic <- lmer(formula_basic, data = combined_data)
      # Specify nuisance variables for emmeans to avoid grid explosion
      nuis_adj <- c("BATCH", covariates)
      nuis_basic <- "BATCH"
    } else {
    # Fixed model
      #Adjusted model (covariates)
      formula_adj <- as.formula(paste0("`", taxon, "` ~ ", fixed_part_adj, " + subclass"))
      model_adj <- lm(formula_adj, data = combined_data)
      # Basic model (no covariates)
      formula_basic <- as.formula(paste0("`", taxon, "` ~ ", fixed_part_basic))
      model_basic <- lm(formula_basic, data = combined_data)
      # Specify nuisance variables for emmeans to avoid grid explosion
      nuis_adj <- c("BATCH", covariates, "subclass")
      nuis_basic <- "BATCH"
    }
    
    # Extract the stats from each model ( for the adjusted and the basic one)
    stats_adj <- extract_stats(model_adj, taxon, nuis_adj)
    stats_basic <- extract_stats(model_basic, taxon, nuis_basic)
    
    # Bind results
    for (comp in c("recent", "cont", "duration")) {
      results$adj[[comp]] <- rbind(results$adj[[comp]], stats_adj[[comp]])
      results$basic[[comp]] <- rbind(results$basic[[comp]], stats_basic[[comp]])
    }
  }
  
  # Specify results variable per each analysis
  results_recent <- results$adj$recent
  results_cont <- results$adj$cont
  results_duration <- results$adj$duration
  results_recent_basic <- results$basic$recent
  results_cont_basic <- results$basic$cont
  results_duration_basic <- results$basic$duration
  
  # Multiple testing correction for all results
  res_list <- list(
    adj_recent = results_recent, adj_cont = results_cont, adj_dur = results_duration,
    basic_recent = results_recent_basic, basic_cont = results_cont_basic, basic_dur = results_duration_basic
  )
  res_list <- lapply(res_list, function(adjustment) {
    adjustment$p_adj <- p.adjust(adjustment$p, method = "fdr")
    return(adjustment)
  })
  

  # Update variables with p adjusted from list
  results_recent <- res_list$adj_recent
  results_cont <- res_list$adj_cont
  results_duration <- res_list$adj_dur
  results_recent_basic <- res_list$basic_recent
  results_cont_basic <- res_list$basic_cont
  results_duration_basic <- res_list$basic_dur
  
  # Preparation and plotting
  prepare_plot_data <- function(regression_results) {
    results_plot <- regression_results %>%
      mutate(
        lower = estimate - 1.96 * se,
        upper = estimate + 1.96 * se,
        significant = case_when(
          p_adj < 0.05 ~ "Significant (q < 0.05)",
          p < 0.05 & p_adj >= 0.05 ~ "p < 0.05 but not FDR",
          TRUE ~ "Non-significant"
        )
      )
      
    # Map species names for species level (genus/family already use final names)
    if (level == "species") {
      results_plot <- results_plot %>%
        left_join(taxonomies, by = c("taxon_id" = "MGS")) %>%
        filter(genus != "Other") %>%
        mutate(plot_label = species)
    } else {
      results_plot <- results_plot %>% mutate(plot_label = taxon_id)
    }
    return(results_plot)
  }
  
  
  # Forest plot function
  plot_forest <- function(plot_data, title, xlabel) {
    ggplot(plot_data, aes(x = estimate, y = plot_label)) +
      geom_point(aes(color = significant)) +
      geom_errorbar(aes(xmin = lower, xmax = upper), orientation = "y", height = 0.2) +
      geom_vline(xintercept = 0, linetype = "dashed") +
      labs(x = xlabel, y = tools::toTitleCase(level), title = title)
  }

  # Use the prepare plot data and forest plot function for each comparison and save
  unit_str <- ifelse(type == "CLR", "CLR", "Relative Abundance")
  
  model_results <- list(
    list(results = list(recent = results_recent, cont = results_cont, dur = results_duration), 
         label = "Adjusted", suffix = "adj"),
    list(results = list(recent = results_recent_basic, cont = results_cont_basic, dur = results_duration_basic), 
         label = "Basic", suffix = "basic")
  )

  for (mod_res in model_results) {
    comparisons <- list(
      list(res = mod_res$results$recent, suffix = "recent_vs_control", title = "Control vs Recent"),
      list(res = mod_res$results$cont, suffix = "cont_vs_control", title = "Control vs Continuous"),
      list(res = mod_res$results$dur, suffix = "duration_effect", title = "Recent vs Continuous")
    )

    for (comp in comparisons) {
      plot_df <- prepare_plot_data(comp$res)
      p <- plot_forest(plot_df, 
                         paste(tools::toTitleCase(level), type, comp$title, 
                                paste0("(", tools::toTitleCase(model_type), " ", mod_res$label, " Effects)")),
                         paste("Difference in", unit_str, "(", comp$title, ")"))
      ggsave(paste0(output_prefix, "_", mod_res$suffix, "_", comp$suffix, ".png"), p, width = 13, height = 7)
    }
  }
  
  return(list(
    results_recent = results_recent,
    results_cont = results_cont,
    results_duration = results_duration,
    results_recent_basic = results_recent_basic,
    results_cont_basic = results_cont_basic,
    results_duration_basic = results_duration_basic
  ))
}

# ----------------------------------------------------------------------------
# 6. Execute the analyses (fixed and mixed models) per each data set

# Specify all datasets and characteristics
dataset <- list(
  list(file = "Data/Female_CLR_species_HT_theses.csv", 
       level = "species", type = "CLR"),
  list(file = "Data/fem_species_relabundance.csv",
       level = "species", type = "Relab"),
  list(file = "Data/fem_genus_clr.csv", 
       level = "genus", type = "CLR"),
  list(file = "Data/fem_genus_relabundance.csv",
       level = "genus", type = "Relab"),
  list(file = "Data/fem_family_clr.csv", 
       level = "family", type = "CLR"),
  list(file = "Data/fem_family_relabundance.csv", 
       level = "family", type = "Relab")
)

# Create empty results list
all_results <- list()


for (analysis in dataset) {
  # Load and prepare dataset
  taxa_data <- prepare_data(analysis$file, analysis$level, analysis$type)
  
  # Run Fixed model and save in results list
  prefix_fixed <- paste0("forest_", analysis$level, "_", 
                         tolower(analysis$type), "_fixed")
  res_fixed <- run_duration_model(taxa_data, analysis$level, analysis$type,
                                  "fixed", prefix_fixed)
  all_results[[paste(analysis$level, analysis$type, "fixed",
                     sep = "_")]] <- res_fixed
  
  # Run Mixed model and save in results list
  prefix_mixed <- paste0("forest_", analysis$level, "_",
                         tolower(analysis$type), "_mixed")
  res_mixed <- run_duration_model(taxa_data, analysis$level, analysis$type,
                                  "mixed", prefix_mixed)
  all_results[[paste(analysis$level, analysis$type, "mixed", 
                     sep = "_")]] <- res_mixed
}


# Save workspace

save.image(file = "Rdata/duration_all_analyses_rdata.RData")
