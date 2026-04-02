# Main Project: UrinaryMicrobiota (KIJURI)
# Author: María Gadea Martínez
# Sub-project: Impact of HT (hormone therapy) in the urinary microbiome of 
# peri- or menopausal women
# Linear regression for all taxonomic levels and data types
# Initial date: 2026-03-31

# Load the libraries 

library(rio)
library(dplyr)
library(ggplot2)

setwd("C:/Users/riade/OneDrive - Uppsala universitet/Thesis project TFM")

####################################################################################################
# 1. Shared data (taxonomy and metadata)
taxonomies <- rio::import("Data/taxonomy_students_females.csv")
match_covariates <- rio::import("Data/FEMALE_Covariates_Cases_Control.csv", dec = ",")
rownames(match_covariates) <- match_covariates[, 1]
match_covariates <- match_covariates[, -1]

####################################################################################################
# 2. Function to run the regression analysis for a given taxonomic level and data type
run_regression_analysis <- function(data_file, data_type, tax_level, output_rdata, plot_widths, plot_heights) {
  
  # Load specific dataset
  dataset <- rio::import(data_file, dec = ",")
  rownames(dataset) <- dataset[, 1]
  dataset <- dataset[, -1]
  
  # Collapse dataset if it's family or genus level
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
  
  # Merge relative abundance data and metadata
  common <- intersect(rownames(match_covariates), rownames(dataset))
  data <- cbind(match_covariates[common, ], dataset[common, ])
  
  id_col <- paste0(tax_level, "_id")
  features <- colnames(dataset)
  
  # Helper function to run the model
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
  
  # Run Basic model
  basic_results <- run_model(" ~ case_bin + BATCH")
  
  # Run Adjusted model
  adjusted_formula <- " ~ case_bin + BATCH + ageatvisitone + DMcat + bmi + derived_smoke_status + physact + alcohol + placebirth + education"
  adjusted_results <- run_model(adjusted_formula)
  
  # Formatting data for forest plots
  format_plot_data <- function(results) {
    plot_data <- results |>
      left_join(taxonomies, by = setNames("MGS", id_col)) |>
      mutate(
        lower = estimate - 1.96 * se,
        upper = estimate + 1.96 * se,
        significant = case_when(
          p_adj < 0.05 ~ "Significant (q < 0.05)",
          p < 0.05 & p_adj >= 0.05 ~ "Non-significant after multiple testing adjustment (p < 0.05 and q > 0.05)",
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
  
  # Generate Plots
  x_label <- ifelse(data_type == "clr", 
                    "Difference in CLR ( no HT vs HT)", 
                    "Difference in relative abundance ( no HT vs HT)")
  
  title_suffix <- if (tax_level == "family") "Family" else if (tax_level == "genus") "Genus" else "Species"
  
  generate_plot <- function(p_data, p_title, p_filename, width, height) {
    p <- ggplot(p_data, aes(x = estimate, y = plot_y_var)) +
      geom_point(aes(color = significant)) +
      geom_errorbar(aes(xmin = lower, xmax = upper), orientation = "y", width = 0.2) +
      geom_vline(xintercept = 0, linetype = "dashed") +
      labs(x = x_label, y = title_suffix, title = p_title)
    
    # Using print to view when running interactively
    print(p)
    ggsave(filename = p_filename, plot = p, width = width, height = height)
    return(p)
  }
  
  file_basic <- paste0(data_type, "_", tax_level, "_basic.png")
  file_adjusted <- paste0(data_type, "_", tax_level, "_adjusted.png")
  
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
  
  # Format adjusted results for significance output
  formatted_results <- adjusted_plot_data %>%
    filter(p_adj < 0.05) %>%
    mutate(
      output = paste0(
        plot_y_var, " (β = ", round(estimate, 3), ", 95% CI: ", round(lower, 3), " to ", round(upper, 3), ")"
      )
    ) %>%
    pull(output)
  
  # Save output environment for this particular regression to reproduce functionality
  save(data, basic_results, adjusted_results, basic_plot_data, adjusted_plot_data, 
       basic_plot, adjusted_plot, formatted_results, file = output_rdata)
}

####################################################################################################
# 3. Define configurations for all 6 analyses
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

####################################################################################################
# 4. Execution block: Iterate through the configurations and run each analysis
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
