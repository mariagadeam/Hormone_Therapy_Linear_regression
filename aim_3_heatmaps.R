# Main Project: UrinaryMicrobiota (KIJURI)
# Author: María Gadea Martínez
# Sub-project: Impact of HT (hormone therapy) in the urinary microbiome of 
# peri- or menopausal women
# Analysis of continuous duration effect using defined time bins 
# Initial date: 2026-04-16
# Last update: 2026-05-25

# ---------------------------------------------------------------------------
## 1. Load the necessary libraries
library(rio)
library(dplyr)
library(ggplot2)
library(tibble)
library(svglite)
library(openxlsx)

# ---------------------------------------------------------------------------
## 2. Set working directory and load common data (set sample number to row name)

setwd("C:/Users/riade/OneDrive - Uppsala universitet/Thesis project TFM")
load("Rdata/duration_binned_analyses.RData")

# Taxonomy assignments
taxonomies <- rio::import("Data/taxonomy_students_females.csv")

# Covariates for cases and controls
cases_covariates <- rio::import("Data/Cases_HRT.csv", dec = ",")
controls_covariates <- rio::import("Data/Controls_HRT.csv", dec = ",")
match_covariates <- bind_rows(cases_covariates, controls_covariates)
rownames(match_covariates) <- match_covariates[, 1]
match_covariates <- match_covariates[, -1]


# ---------------------------------------------------------------------------
## 3. Plots and counts tables to analyze distribution of cases per duration (months)

# 3.1 Define unified bin labels matching the rest of the analyses
dist_bin_labels <- c("1-2 months", "3-4 months", "5-6 months", "7-8 months",
                     "9-10 months", "11-12 months", "2-3 years", "> 3 years")

# 3.2 Assign each case to its bin using the same break points used in the regression
cases_covariates$duration_bin <- cut(
  cases_covariates$continuous_days,
  breaks = c(0, 60, 120, 180, 240, 300, 365, 1095, Inf),
  labels = dist_bin_labels,
  right  = TRUE,
  include.lowest = TRUE
)

# 3.3 Bar plot of case counts per duration bin
dist_plot <- ggplot(cases_covariates, aes(x = duration_bin)) +
  geom_bar(fill = "steelblue") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11)) +
  labs(
    title = "Distribution of HT Cases by Duration of Therapy",
    x     = "Duration of Hormone Therapy",
    y     = "Number of Cases"
  )
dist_plot

# 3.4 Counts table
dist_counts <- table(cases_covariates$duration_bin)
print(dist_counts)


# --------------------------------------------------------------------------
## 4. Function to prepare abundance (CLR or relative abundance) dataset
# If genus or family level, it needs to be collapsed

prepare_data <- function(filepath, level, type) {
  
  # 4.1 Load abundance data
  abundance_raw <- rio::import(filepath, dec = ",")
  rownames(abundance_raw) <- abundance_raw[, 1]
  abundance_out <- abundance_raw[, -1]
  
  # 4.2 Collapse if genus or family
  if ((level == "genus" || level == "family") && type == "Relab") {
    lookup <- setNames(taxonomies[[level]], taxonomies$MGS)
    # Replace column names with corresponding taxonomy names where possible
    colnames(abundance_out) <- ifelse(
      colnames(abundance_out) %in% names(lookup),
      lookup[colnames(abundance_out)], 
      colnames(abundance_out))
    # Aggregate abundances by identical column names (same genera / family)
    abundance_out <- as.data.frame(
      rowsum(t(abundance_out),
             group = sub("\\..*$", "", colnames(abundance_out))))
    
    abundance_out <- as.data.frame(t(abundance_out))
  }
  
  return(abundance_out)
}

# --------------------------------------------------------------------------
# 5. Main Linear Regression Analysis Function

run_binned_model <- function(taxa, level, type, output_prefix) {
  
  # 5.1 Prepare data, time bins and covariates
  
   # Merge metadata and duration data
  common <- intersect(rownames(match_covariates), rownames(taxa))
  combined_data <- cbind(match_covariates[common, ], taxa[common, ])

   # Label samples without case type as "control" and give value 0
  combined_data$case_type[is.na(combined_data$case_type)] <- "control"
  combined_data$continuous_days[combined_data$case_type == "control"] <- 0
  
   # Convert continuous_days to a numeric type
  combined_data$continuous_days <- as.numeric(combined_data$continuous_days)
  
   # Categorical time bins — same break points and labels used throughout the analyses
  bin_levels <- c("Control", "1-2 months", "3-4 months",
                  "5-6 months", "7-8 months", "9-10 months", "11-12 months",
                  "2-3 years", "> 3 years")
  
  combined_data$duration_bin <- cut(
    combined_data$continuous_days,
    breaks = c(-Inf, 0.1, 60, 120, 180, 240, 300, 365, 1095, Inf),
    labels = bin_levels
  )
  
   # Convert the bin variable to a factor with the desired levels
  combined_data$duration_bin <- factor(combined_data$duration_bin, levels = bin_levels)
   # Ensure controls are labelled correctly
  combined_data$duration_bin[combined_data$case_type == "control"] <- "Control"
  
  # Define covariates for regression adjustment
  covariates <- c("ageatvisitone", "DMcat", "bmi", "derived_smoke_status",
                  "physact", "alcohol", "placebirth", "education")
  
  # Initialize result data frames
  results <- list(adj = data.frame(), basic = data.frame())
  
  # 5.2 Helper function to extract per-bin coefficients from a fitted model
  extract_stats <- function(model_fit, taxon) {
    coef_table <- summary(model_fit)$coefficients
    res_list <- lapply(bin_levels[-1], function(bin) {
      coef_name <- paste0("duration_bin", bin)
      if (coef_name %in% rownames(coef_table)) {
        res <- coef_table[coef_name, ]
        return(data.frame(taxon_id = taxon, bin = bin,
                          estimate = res["Estimate"],
                          se       = res["Std. Error"],
                          p        = res["Pr(>|t|)"]))
      }
      return(data.frame(taxon_id = taxon, bin = bin,
                        estimate = NA, se = NA, p = NA))
    })
    res_df <- do.call(rbind, res_list)
    return(res_df)
  }
  
  taxa_cols <- colnames(taxa)
  
  # 5.3 Loop through each taxon 
  for (taxon in taxa_cols) {
    if (taxon %in% c("N", "Other")) next
    
    # Adjusted model formula (with covariates, plain lm — no subclass)
    fixed_part_adj   <- paste(c("duration_bin", "BATCH", covariates), collapse = " + ")
    formula_adj      <- as.formula(paste0("`", taxon, "` ~ ", fixed_part_adj))
    model_adj        <- lm(formula_adj, data = combined_data)
    
    # Basic model formula (batch only)
    formula_basic    <- as.formula(paste0("`", taxon, "` ~ duration_bin + BATCH"))
    model_basic      <- lm(formula_basic, data = combined_data)
    
    results$adj   <- rbind(results$adj,   extract_stats(model_adj,   taxon))
    results$basic <- rbind(results$basic, extract_stats(model_basic, taxon))
  }
  
  # 5.4 FDR correction 
  res_list <- lapply(results, function(adjustment) {
    if (nrow(adjustment) > 0) {
      adjustment$p_adj <- p.adjust(adjustment$p, method = "fdr")
      adjustment <- adjustment[!is.na(adjustment$estimate), ]
    }
    return(adjustment)
  })
  
  # 5.5 Collect significant results (q < 0.05) 

  sig_rows <- do.call(rbind, lapply(names(res_list), function(adj_name) {
    df <- res_list[[adj_name]]
    df_sig <- df[!is.na(df$p_adj) & df$p_adj < 0.05, ]
    if (nrow(df_sig) == 0) return(NULL)
    df_sig$level      <- level
    df_sig$type       <- type
    df_sig$adjustment <- adj_name
    df_sig
  }))
  
  # 5.6 Heatmap plotting function
  plot_heatmap <- function(plot_data, title, output_name) {
    plot_data <- plot_data %>% mutate(
      significant = case_when(
        p_adj < 0.05                  ~ "**",
        p < 0.05 & p_adj >= 0.05     ~ "*",
        TRUE                          ~ ""
      ),
      bin = factor(bin, levels = bin_levels[-1])
    )
    
    if (level == "species") {
      plot_data <- plot_data %>%
        left_join(taxonomies, by = c("taxon_id" = "MGS")) %>%
        filter(genus != "Other") %>%
        mutate(plot_label = species)
    } else {
      plot_data <- plot_data %>% mutate(plot_label = taxon_id)
    }
    
    heatmap <- ggplot(plot_data, aes(x = bin, y = plot_label, fill = estimate)) +
      geom_tile(color = "grey80") +
      scale_fill_gradient2(low = "royalblue", mid = "white", high = "firebrick",
                           midpoint = 0, name = "Effect\n(vs Control)") +
      geom_text(aes(label = significant), vjust = 0.7, size = 5) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
        axis.text.y = element_text(size = 9),
        panel.grid  = element_blank()
      ) +
      labs(
        title   = title,
        x       = "Duration of Hormone Therapy",
        y       = tools::toTitleCase(level),
        caption = "* p < 0.05   ** q < 0.05\nPlot displays ALL analyzed taxa."
      )
    
    ggsave(paste0(output_name, ".svg"), heatmap,
           width = 14,
           height = max(6, length(unique(plot_data$plot_label)) * 0.25),
           limitsize = FALSE)
  }
  
  plot_heatmap(res_list$adj,
               paste(tools::toTitleCase(level), type, "Time Bin Effects (Adjusted)"),
               paste0(output_prefix, "_adj"))
  plot_heatmap(res_list$basic,
               paste(tools::toTitleCase(level), type, "Time Bin Effects (Basic)"),
               paste0(output_prefix, "_basic"))
  
  # Return both the full results (for storage) and the pre-filtered sig rows
  return(list(results = res_list, sig_rows = sig_rows))
}

# ----------------------------------------------------------------------------
## 6. Execute the analysis (plain linear regression) per each dataset

# 6.1 Specify all datasets and characteristics
dataset <- list(
  list(file = "Data/Cohort_data/Female891_CLR_species_theses.csv", level = "species", type = "CLR"),
  list(file = "Data/Cohort_data/fem_species_relabundance.csv",      level = "species", type = "Relab"),
  list(file = "Data/Cohort_data/Female891_GENUS_CLR_theses.csv",                 level = "genus",   type = "CLR"),
  list(file = "Data/Cohort_data/fem_genus_relabundance.csv",        level = "genus",   type = "Relab"),
  list(file = "Data/Cohort_data/Female891_FAM_CLR_theses.csv",                level = "family",  type = "CLR"),
  list(file = "Data/Cohort_data/fem_family_relabundance.csv",       level = "family",  type = "Relab")
)

# 6.2 Create empty results lists
all_results  <- list()   # full results (for storage / heatmaps)
all_sig_rows <- list()   # significant rows (q < 0.05) collected per analysis

# 6.3 Loop across each analytical condition
for (analysis in dataset) {
  
  # Load and prepare dataset
  taxa_data <- prepare_data(analysis$file, analysis$level, analysis$type)
  
  # Run plain linear regression and save results
  prefix    <- paste0("heatmap_binned_", analysis$level, "_", tolower(analysis$type))
  key       <- paste(analysis$level, analysis$type, sep = "_")
  out       <- run_binned_model(taxa_data, analysis$level, analysis$type, prefix)
  
  all_results[[key]]  <- out$results    # list(adj = df, basic = df)
  all_sig_rows[[key]] <- out$sig_rows   # data.frame of q < 0.05 rows
}

# 6.4 Save workspace

save.image(file = "Rdata/duration_binned_analyses.RData")


