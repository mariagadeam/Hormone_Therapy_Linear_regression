# Main Project: UrinaryMicrobiota (KIJURI)
# Author: María Gadea Martínez
# Sub-project: Impact of HT (hormone therapy) in the urinary microbiome of 
# peri- or menopausal women
# Linear regression for all taxonomic levels and data types
# Initial date: 2026-03-31
# Last update: 2026-05-23

# Load the libraries 
library(svglite)
library(rio)
library(dplyr)
library(ggplot2)

# Set working directory and load Rdata (if previously run)
setwd("C:/Users/riade/OneDrive - Uppsala universitet/Thesis project TFM")
load("Rdata/lm_matched_case_control.RData")

# ---------------------------------------------------------------------------
## 1. Shared data (taxonomy and metadata)
taxonomies <- rio::import("Data/taxonomy_students_females.csv")
match_covariates <- rio::import("Data/FEMALE_Covariates_Cases_Control.csv", dec = ",")
rownames(match_covariates) <- match_covariates[, 1]
match_covariates <- match_covariates[, -1]

# ---------------------------------------------------------------------------
## 2. Function to run the regression analysis for a given taxonomic level and data type
run_regression_analysis <- function(data_file, data_type, tax_level, output_rdata, plot_widths, plot_heights) {
  
  # 2.1 Load specific dataset
  dataset <- rio::import(data_file, dec = ",")
  rownames(dataset) <- dataset[, 1]
  dataset <- dataset[, -1]
  
  # 2.2 Collapse dataset if it's family or genus level
  if (tax_level %in% c("family", "genus")) {
    lookup <- setNames(taxonomies[[tax_level]], taxonomies$MGS)
    colnames(dataset) <- ifelse(
      colnames(dataset) %in% names(lookup),
      lookup[colnames(dataset)],
      colnames(dataset)
    )
    dataset <- as.data.frame(
      rowsum(
        t(dataset),
        group = sub("\\..*$", "", colnames(dataset))
      )
    )
    dataset <- t(dataset)
  }
  
  # 2.3 Merge relative abundance data and metadata
  common <- intersect(rownames(match_covariates), rownames(dataset))
  data <- cbind(match_covariates[common, ], dataset[common, ])
  
  id_col <- paste0(tax_level, "_id")
  features <- colnames(dataset)
  
  # 2.4 Helper function to run the model
  run_model <- function(formula_str) {
    res <- data.frame(
      id_temp = character(), estimate = numeric(),
      se = numeric(), p = numeric(), stringsAsFactors = FALSE
    )
    
    for (feat in features) {
      current_formula <- as.formula(paste0("`", feat, "`", formula_str))
      fit <- lm(current_formula, data = data)
      coef_obj <- summary(fit)$coefficients["case_bin", ]
      
      res <- rbind(res, data.frame(
        id_temp = feat, estimate = coef_obj["Estimate"],
        se = coef_obj["Std. Error"], p = coef_obj["Pr(>|t|)"]
      ))
    }
    
    colnames(res)[1] <- id_col
    res$p_adj <- p.adjust(res$p, method = "fdr")
    return(res)
  }
  
  # 2.5 Run Basic and Adjusted model
  basic_results <- run_model(" ~ case_bin + BATCH")

  adjusted_formula <- " ~ case_bin + BATCH + ageatvisitone + DMcat + bmi + derived_smoke_status + physact + alcohol + placebirth + education"
  adjusted_results <- run_model(adjusted_formula)
  
  # 2.6 Formatting data for forest plots
  format_plot_data <- function(results) {
    plot_data <- results |>
      left_join(taxonomies, by = setNames("MGS", id_col)) |>
      mutate(
        lower = estimate - 1.96 * se,
        upper = estimate + 1.96 * se,
        significant = case_when(
          p_adj < 0.05 ~ "Significant (q < 0.05)",
          p < 0.05 & p_adj >= 0.05 ~ "Non-significant after multiple testing adjustment (p < 0.05 and q >= 0.05)",
          TRUE ~ "Non-significant (p > 0.05)"
        )
      )
    
    # Specific filtering & variable settings
    if (tax_level == "family") {
      plot_data <- dplyr::filter(plot_data, !get(id_col) %in% c("N", "Other"))
      plot_data$plot_y_var <- plot_data[[id_col]]
    } else if (tax_level == "genus") {
      plot_data <- dplyr::filter(plot_data, !get(id_col) %in% c("Other"))
      plot_data <- plot_data |> mutate(
        plot_y_var = ifelse(get(id_col) == "Other", "Other", get(id_col))
      )
    } else if (tax_level == "species") {
      plot_data <- dplyr::filter(plot_data, genus != "Other")
      plot_data <- plot_data |> mutate(
        species = ifelse(get(id_col) == "Other", "Other", species),
        plot_y_var = species
      )
    }
    
    return(plot_data)
  }
  
  basic_plot_data <- format_plot_data(basic_results)
  adjusted_plot_data <- format_plot_data(adjusted_results)
  
  # 2.7 Generate Plots function 
  
    # Specify data trasnformation in x label
  x_label <- ifelse(data_type == "clr", 
                    "Difference in CLR ( no HT vs HT)", 
                    "Difference in RA ( no HT vs HT)")
  
   # Specify taxonomic rank in title
  title_suffix <- if (tax_level == "family") "Family" else if (tax_level == "genus") "Genus" else "Species"
  
   # Function for generating plot with significance specified
  generate_plot <- function(p_data, p_title, p_filename, width, height) {
    sig_colors <- c(
      "Non-significant (p > 0.05)" = "red",
      "Non-significant after multiple testing adjustment (p < 0.05 and q >= 0.05)" = "green",
      "Significant (q < 0.05)" = "blue"
    )
    p <- ggplot(p_data, aes(x = estimate, y = plot_y_var)) +
      geom_point(aes(color = significant), size = 3) +
      geom_errorbar(aes(xmin = lower, xmax = upper), orientation = "y", width = 0.2) +
      geom_vline(xintercept = 0, linetype = "dashed") +
      scale_color_manual(name = "Significance", values = sig_colors) +
      labs(x = x_label, y = title_suffix, title = p_title) +
      theme(legend.position = "none") 
    
    # Using print to view when running interactively
    print(p)
    ggsave(
      filename = p_filename,
      plot = p,
      width = 11,
      height = 7,
      device = "svg"
    )
    return(p)
  }
  
  # Name for image file
  file_basic <- paste0(data_type, "_", tax_level, "_basic.svg")
  file_adjusted <- paste0(data_type, "_", tax_level, "_adjusted.svg")
  
  # 2.8 Generate basic and adjusted models plots
  basic_plot <- generate_plot(
    basic_plot_data, 
    paste0(title_suffix, " associated with HRT (basic model)"), 
    file_basic, 
    plot_widths["basic_w"], plot_heights["basic_h"]
  )
  
  adjusted_plot <- generate_plot(
    adjusted_plot_data, 
    paste0(title_suffix, " associated with HRT (adjusted)"), 
    file_adjusted, 
    plot_widths["adj_w"], plot_heights["adj_h"]
  )
  
  # 2.9 Format adjusted results for significance output
  formatted_results <- adjusted_plot_data %>%
    filter(p_adj < 0.05) %>%
    mutate(
      output = paste0(
        plot_y_var, " (β = ", round(estimate, 3), ", 95% CI: ", round(lower, 3), " to ", round(upper, 3), ")"
      )
    ) %>%
    pull(output)
  
  # 2.10 Save output environment for this particular regression to reproduce functionality
  save(data, basic_results, adjusted_results, basic_plot_data, adjusted_plot_data, 
       basic_plot, adjusted_plot, formatted_results, file = output_rdata)
}

## ---------------------------------------------------------------------------
## 3. Execute analyses for all datasets

 # 3.1 Define configurations for all 6 analyses
analyses <- list(
  # CLR model
  list(
    file = "Data/fem_family_clr.csv", type = "clr", level = "family", 
    rdata = "Rdata/family_clr_rdata.RData", 
    w = c(basic_w = 9, adj_w = 9), h = c(basic_h = 6, adj_h = 6)
  ),
  list(
    file = "Data/fem_genus_clr.csv", type = "clr", level = "genus", 
    rdata = "Rdata/genus_clr_rdata.RData", 
    w = c(basic_w = 13, adj_w = 11), h = c(basic_h = 6, adj_h = 6)
  ),
  list(
    file = "Data/Female_CLR_species_HT_theses.csv", type = "clr", level = "species", 
    rdata = "Rdata/species_clr_rdata.RData", 
    w = c(basic_w = 13, adj_w = 13), h = c(basic_h = 9, adj_h = 9)
  ),
  
  # Relative abundance model 
  list(
    file = "Data/fem_family_relabundance.csv", type = "relab", level = "family", 
    rdata = "Rdata/family_relab_rdata.RData", 
    w = c(basic_w = 9, adj_w = 9), h = c(basic_h = 6, adj_h = 6)
  ),
  list(
    file = "Data/fem_genus_relabundance.csv", type = "relab", level = "genus", 
    rdata = "Rdata/genus_relab_rdata.RData", 
    w = c(basic_w = 13, adj_w = 11), h = c(basic_h = 6, adj_h = 6)
  ),
  list(
    file = "Data/fem_species_relabundance.csv", type = "relab", level = "species", 
    rdata = "Rdata/species_relab_rdata.RData", 
    w = c(basic_w = 13, adj_w = 13), h = c(basic_h = 7, adj_h = 7)
  )
)

 # 3.2 Execution block, iterate through the configurations and run each analysis
for (analysis in analyses) {
  run_regression_analysis(
    data_file = analysis$file,
    data_type = analysis$type,
    tax_level = analysis$level,
    output_rdata = analysis$rdata,
    plot_widths = analysis$w,
    plot_heights = analysis$h
  )
}


# ---------------------------------------------------------------------------

## 4. Save workspace

save.image(file = "Rdata/lm_matched_case_control.RData")