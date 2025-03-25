# =====================================================================
# Flexible FST-to-Survey Pipeline Script (Modernized)
#
# This pipeline does the following:
# 1. Reads one or more flagged FST files from a specified directory,
#    using purrr for fast, safe file reading.
# 2. Drops diagnosis columns (matching regex patterns) with optional
#    additional patterns provided by the user.
# 3. Applies an optional subsetting filter (e.g., selecting only HIV+ cases).
# 4. Creates a survey design object (via srvyr) that supports one- or
#    two-stage (clustered) designs.
#
# The final survey design is equivalent to:
#
#   as_survey_design(ids = ~HOSPID_final + HOSP_NIS,
#                    weights = ~DISCWT,
#                    strata = ~NIS_STRATUM,
#                    nest = TRUE,
#                    multicore = getOption("survey.multicore"))
#
# All defaults are overridable.
# =====================================================================

suppressPackageStartupMessages({
  library(tidyverse)  # For data manipulation and purrr
  library(fst)        # For fast reading of fst files
  library(srvyr)      # For survey design (tidy interface)
  library(rlang)      # For non-standard evaluation (tidy eval)
})

# -----------------------------------------------------------------------------
# Function: load_and_clean_fst
#
# Reads all .fst files in a directory, combines them, drops diagnosis
# columns (default & additional patterns), and optionally subsets the data.
#
# Parameters:
#   fst_dir         : Directory containing the .fst files.
#   drop_patterns   : Default regex patterns for diagnosis columns.
#   additional_drop : Optional additional regex patterns to drop.
#   subset_expr     : Optional unquoted filtering expression.
#   verbose         : Logical; if TRUE, prints progress messages.
#
# Returns:
#   A tibble containing the combined and cleaned data.
# -----------------------------------------------------------------------------
load_and_clean_fst <- function(fst_dir,
                               drop_patterns   = c("^DX", "^I10_DX", "^PR", "^I10_PR"),
                               additional_drop = NULL,
                               subset_expr     = NULL,
                               verbose         = TRUE) {
  # List all .fst files in the directory
  fst_files <- list.files(fst_dir, pattern = "\\.fst$", full.names = TRUE)
  if (length(fst_files) == 0) {
    stop("No FST files found in directory: ", fst_dir)
  }
  if (verbose) message("Found ", length(fst_files), " FST file(s) in ", fst_dir)

  # Read files using purrr::map() with safe error handling
  safe_read <- safely(~ as_tibble(fst::read_fst(.x)))
  df_list <- fst_files %>%
    map(function(file) {
      if (verbose) message("Reading file: ", file)
      result <- safe_read(file)
      if (!is.null(result$error)) {
        warning("Failed to read ", file, ": ", result$error$message)
        return(NULL)
      }
      result$result
    }) %>%
    compact()  # Remove any NULL entries

  if (length(df_list) == 0) {
    stop("No valid data loaded from FST files.")
  }
  combined_df <- bind_rows(df_list)
  if (verbose) {
    message("Combined data: ", nrow(combined_df), " rows; ", ncol(combined_df), " columns.")
  }

  # Create combined drop patterns
  all_drop_patterns <- c(drop_patterns, additional_drop)
  drop_cols <- names(combined_df)[
    map_lgl(names(combined_df), ~ any(str_detect(.x, all_drop_patterns)))
  ]
  if (verbose) {
    if (length(drop_cols) > 0) {
      message("Dropping columns: ", paste(drop_cols, collapse = ", "))
    } else {
      message("No diagnosis columns found to drop.")
    }
  }
  cleaned_df <- combined_df %>% select(-all_of(drop_cols))

  # Apply subset filter if provided
  if (!is.null(subset_expr)) {
    sub_expr <- enquo(subset_expr)
    cleaned_df <- cleaned_df %>% filter(!!sub_expr)
    if (verbose) {
      message("Data subset applied; new row count: ", nrow(cleaned_df))
    }
  }

  cleaned_df
}

# -----------------------------------------------------------------------------
# Function: create_survey_object
#
# Creates a survey design object using srvyr with flexible design options.
#
# Parameters:
#   data         : Cleaned data as a tibble.
#   id_var       : Primary clustering variable.
#   weight_var   : Weight variable.
#   strata_var   : Stratification variable.
#   cluster_var  : Secondary clustering variable (for two-stage designs).
#   nest         : Logical; passed to as_survey_design().
#   multicore    : Logical; use multicore (if available).
#   verbose      : Logical; print progress messages.
#
# Returns:
#   A survey design object (tbl_svy).
# -----------------------------------------------------------------------------
create_survey_object <- function(data,
                                 id_var      = "HOSPID_final",
                                 weight_var  = "DISCWT",
                                 strata_var  = "NIS_STRATUM",
                                 cluster_var = "HOSP_NIS",
                                 nest        = TRUE,
                                 multicore   = getOption("survey.multicore", FALSE),
                                 verbose     = TRUE) {
  # Ensure required variables exist
  required_vars <- c(id_var, weight_var, strata_var)
  if (!is.null(cluster_var) && cluster_var != "") {
    required_vars <- c(required_vars, cluster_var)
  }
  missing_vars <- setdiff(required_vars, names(data))
  if (length(missing_vars) > 0) {
    stop("Missing required survey design variables: ", paste(missing_vars, collapse = ", "))
  }

  # Build cluster ids formula
  ids_formula <- if (!is.null(cluster_var) && cluster_var != "") {
    as.formula(paste("~", id_var, "+", cluster_var))
  } else {
    as.formula(paste("~", id_var))
  }
  if (verbose) message("Using clustering: ", deparse(ids_formula))

  # Build formulas for weights and strata
  weight_formula <- as.formula(paste("~", weight_var))
  strata_formula <- as.formula(paste("~", strata_var))
  if (verbose) {
    message("Creating survey design with weights: ", weight_var,
            " and strata: ", strata_var)
  }

  # Create survey design object with srvyr
  design_obj <- data %>%
    as_survey_design(
      ids       = ids_formula,
      weights   = weight_formula,
      strata    = strata_formula,
      nest      = nest,
      multicore = multicore
    )

  design_obj
}

# -----------------------------------------------------------------------------
# Function: process_fst_to_survey
#
# Main wrapper: loads, cleans FST files and creates a survey design object.
#
# Parameters:
#   fst_dir         : Directory with FST files.
#   drop_patterns   : Default patterns to drop diagnosis columns.
#   additional_drop : Additional patterns to drop.
#   subset_expr     : Optional subsetting filter.
#   id_var, weight_var, strata_var, cluster_var : Survey design variable names.
#   nest, multicore : Options for survey design.
#   verbose         : Logical; print progress messages.
#
# Returns:
#   A srvyr survey design object.
# -----------------------------------------------------------------------------
process_fst_to_survey <- function(fst_dir,
                                  drop_patterns   = c("^DX", "^I10_DX", "^PR", "^I10_PR"),
                                  additional_drop = NULL,
                                  subset_expr     = NULL,
                                  id_var          = "HOSPID_final",
                                  weight_var      = "DISCWT",
                                  strata_var      = "NIS_STRATUM",
                                  cluster_var     = "HOSP_NIS",
                                  nest            = TRUE,
                                  multicore       = getOption("survey.multicore", FALSE),
                                  verbose         = TRUE) {
  cleaned_data <- load_and_clean_fst(
    fst_dir         = fst_dir,
    drop_patterns   = drop_patterns,
    additional_drop = additional_drop,
    subset_expr     = subset_expr,
    verbose         = verbose
  )
  if (verbose) message("Data cleaning complete. Creating survey design object...")

  survey_design <- create_survey_object(
    data        = cleaned_data,
    id_var      = id_var,
    weight_var  = weight_var,
    strata_var  = strata_var,
    cluster_var = cluster_var,
    nest        = nest,
    multicore   = multicore,
    verbose     = verbose
  )

  if (verbose) message("Survey design object successfully created.")
  survey_design
}

# Example usage:
# my_survey <- process_fst_to_survey("path/to/fst_files",
#                                    subset_expr = HIV_flag == 1)
