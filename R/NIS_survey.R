suppressPackageStartupMessages({
  library(tidyverse)  # Data manipulation and purrr for iteration
  library(fst)        # Fast reading of fst files with column selection
  library(srvyr)      # For creating tidy survey design objects
  library(rlang)      # For tidy evaluation
})

# -----------------------------------------------------------------------------
# Function: load_and_clean_fst
#
# Reads each .fst file in a directory, selecting only the columns that the user
# wants to keep. By default, if keep_cols is NULL, the function reads all columns
# and then drops those matching heavy column patterns (i.e. DX or PR columns).
#
# The function applies:
#   - an optional subset filter (subset_expr), a single logical expression, and
#   - a complete-cases filter on the specified columns (complete_cols)
#
# Parameters:
#   fst_dir       : Directory containing the .fst files.
#   keep_cols     : Character vector of column names to import. If NULL,
#                   then all columns are imported and later heavy columns are dropped.
#                   (Defaults drop columns matching DX/PR patterns.)
#   drop_patterns : Default regex patterns for columns to drop when keep_cols is NULL.
#                   Defaults to c("^DX", "^I10_DX", "^PR", "^I10_PR").
#   additional_drop:
#                   Optional additional regex patterns to drop.
#   subset_expr   : Optional unquoted logical filtering expression (e.g., HIV == 1 & !is.na(AGE_grp)).
#   complete_cols : Character vector of column names that must have no missing values.
#                   Rows with NA in any of these columns are dropped.
#   verbose       : Logical; if TRUE, prints progress messages.
#
# Returns:
#   A tibble containing the combined and cleaned data from all FST files.
# -----------------------------------------------------------------------------
load_and_clean_fst <- function(fst_dir,
                               keep_cols      = NULL,
                               drop_patterns  = c("^DX", "^I10_DX", "^PR", "^I10_PR"),
                               additional_drop = NULL,
                               subset_expr    = NULL,
                               complete_cols  = c("FEMALE", "AGE_grp", "RACE", "DIED"),
                               verbose        = TRUE) {
  # Define the mandatory survey design columns
  mandatory_cols <- c("HOSPID_final", "DISCWT", "NIS_STRATUM")

  # If user supplied a vector of columns, ensure the mandatory columns are included.
  if (!is.null(keep_cols)) {
    missing_mandatory <- setdiff(mandatory_cols, keep_cols)
    if (length(missing_mandatory) > 0) {
      message("Mandatory column(s) ", paste(missing_mandatory, collapse = ", "),
              " not found in the provided 'keep_cols'; adding them automatically.")
      keep_cols <- union(keep_cols, mandatory_cols)
    }
  }

  # List all .fst files in the directory
  fst_files <- list.files(fst_dir, pattern = "\\.fst$", full.names = TRUE)
  if (length(fst_files) == 0) {
    stop("No FST files found in directory: ", fst_dir)
  }
  if (verbose) message("Found ", length(fst_files), " FST file(s) in ", fst_dir)

  # Initialize list to store cleaned data frames
  cleaned_list <- vector("list", length(fst_files))

  for (i in seq_along(fst_files)) {
    file <- fst_files[i]
    if (verbose) message("Processing file: ", file)

    # Attempt to read the FST file.
    # If keep_cols is provided, only those columns are imported.
    df <- tryCatch({
      if (!is.null(keep_cols)) {
        as_tibble(fst::read_fst(file, columns = keep_cols))
      } else {
        as_tibble(fst::read_fst(file))
      }
    }, error = function(e) {
      message("FST reading failure for file: ", file,
              ". The file may be corrupt, missing expected structure, or not a valid fst file. Detailed error: ", e$message)
      return(NULL)
    })

    if (is.null(df)) next  # Skip to next file on failure

    # Check that all mandatory columns are present in the data frame.
    missing_in_df <- setdiff(mandatory_cols, names(df))
    if (length(missing_in_df) > 0) {
      warning("File ", file, " is missing mandatory column(s): ",
              paste(missing_in_df, collapse = ", "),
              ". This file will be skipped.")
      next
    }

    # If keep_cols was not specified, drop heavy columns matching drop_patterns (plus any additional_drop)
    if (is.null(keep_cols)) {
      all_drop_patterns <- c(drop_patterns, additional_drop)
      drop_cols <- names(df)[
        map_lgl(names(df), ~ any(str_detect(.x, all_drop_patterns)))
      ]
      if (verbose && length(drop_cols) > 0) {
        message("Dropping columns: ", paste(drop_cols, collapse = ", "))
      }
      df <- df %>% select(-all_of(drop_cols))
    }

    # Apply the user-specified subsetting filter if provided.
    if (!is.null(subset_expr)) {
      sub_expr <- enquo(subset_expr)
      df <- df %>% filter(!!sub_expr)
      if (verbose) message("After subsetting, rows remaining: ", nrow(df))
    }

    # Enforce complete cases for specified columns (only for those that exist in df)
    existing_complete <- intersect(complete_cols, names(df))
    if (length(existing_complete) > 0) {
      before <- nrow(df)
      df <- df %>% drop_na(all_of(existing_complete))
      if (verbose) message("Dropped ", before - nrow(df), " row(s) with missing values in: ",
                           paste(existing_complete, collapse = ", "))
    }

    # Optionally extract year from the filename if it contains a 4-digit year
    year_extracted <- str_extract(basename(file), "\\d{4}")
    if (!is.na(year_extracted)) {
      df <- df %>% mutate(year = as.integer(year_extracted))
    }

    cleaned_list[[i]] <- df

    # Clean up memory
    rm(df); gc()
  }

  # Remove any NULLs from failed reads and combine data
  cleaned_list <- compact(cleaned_list)
  if (length(cleaned_list) == 0) {
    stop("No valid data loaded from FST files (all files missing mandatory columns or failed to load).")
  }

  combined_df <- bind_rows(cleaned_list)
  if (verbose) {
    message("Combined data: ", nrow(combined_df), " rows; ", ncol(combined_df), " columns.")
  }
  combined_df
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
#   multicore    : Logical; if TRUE, uses multicore processing.
#                  (Default is FALSE to avoid excessive core usage.)
#   verbose      : Logical; prints progress messages.
#
# Returns:
#   A survey design object (tbl_svy).
# -----------------------------------------------------------------------------
create_survey_object <- function(data,
                                 id_var      = "HOSPID_final",
                                 weight_var  = "DISCWT",
                                 strata_var  = "NIS_STRATUM",
                                 cluster_var = NULL,
                                 nest        = TRUE,
                                 multicore   = FALSE,
                                 verbose     = TRUE) {
  if (multicore) {
    cores <- parallel::detectCores(logical = TRUE)
    message("Multicore processing enabled. Detected ", cores, " logical cores.")
    if (cores > 8) {
      message("Using a high number of cores (", cores, ") may overwhelm your system. Consider reducing multicore usage.")
    }
  }

  # Ensure required survey design variables exist
  required_vars <- c(id_var, weight_var, strata_var)
  if (!is.null(cluster_var) && cluster_var != "") {
    required_vars <- c(required_vars, cluster_var)
  }
  missing_vars <- setdiff(required_vars, names(data))
  if (length(missing_vars) > 0) {
    stop("Missing required survey design variables: ", paste(missing_vars, collapse = ", "))
  }

  # Build clustering formula
  ids_formula <- if (!is.null(cluster_var) && cluster_var != "") {
    as.formula(paste("~", id_var, "+", cluster_var))
  } else {
    as.formula(paste("~", id_var))
  }
  if (verbose) message("Using clustering: ", deparse(ids_formula))

  # Build weight and strata formulas
  weight_formula <- as.formula(paste("~", weight_var))
  strata_formula <- as.formula(paste("~", strata_var))
  if (verbose) {
    message("Creating survey design with weights: ", weight_var,
            " and strata: ", strata_var)
  }

  # Create and return the survey design object
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
# Main wrapper: Processes each FST file in a directory individually, cleans them
# (by selecting desired columns, applying a subset filter, and enforcing complete cases),
# combines all years of data, creates a survey design object, and saves the final
# survey object to disk in the same directory.
#
# Parameters:
#   fst_dir         : Directory containing FST files.
#   keep_cols       : Character vector of column names to import. If NULL, all columns
#                     are imported and then heavy DX/PR columns are dropped.
#   drop_patterns   : Default regex patterns for columns to drop if keep_cols is NULL.
#                     Defaults to c("^DX", "^I10_DX", "^PR", "^I10_PR").
#   additional_drop : Additional regex patterns for columns to drop.
#   subset_expr     : Single unquoted filtering expression (e.g., HIV == 1 & !is.na(AGE_grp)).
#   complete_cols   : Columns that must have no missing data.
#   id_var, weight_var, strata_var, cluster_var :
#                     Survey design variable names.
#   nest            : Logical; passed to as_survey_design().
#   multicore       : Logical; if TRUE, uses multicore processing (default FALSE).
#   verbose         : Logical; prints progress messages.
#
# Returns:
#   A survey design object which is also saved as "final_survey_object.rds" in fst_dir.
# -----------------------------------------------------------------------------
process_fst_to_survey <- function(fst_dir,
                                  keep_cols      = NULL,
                                  drop_patterns  = c("^DX", "^I10_DX", "^PR", "^I10_PR"),
                                  additional_drop = NULL,
                                  subset_expr    = NULL,
                                  complete_cols  = c("FEMALE", "AGE_grp", "RACE", "DIED"),
                                  id_var         = "HOSPID_final",
                                  weight_var     = "DISCWT",
                                  cluster_var    = NULL,
                                  strata_var     = "NIS_STRATUM",
                                  nest           = TRUE,
                                  multicore      = FALSE,
                                  verbose        = TRUE) {
  if (verbose) message("Starting processing of FST files in directory: ", fst_dir)

  # Process each file and combine cleaned data
  combined_df <- load_and_clean_fst(
    fst_dir        = fst_dir,
    keep_cols      = keep_cols,
    drop_patterns  = drop_patterns,
    additional_drop = additional_drop,
    subset_expr    = subset_expr,
    complete_cols  = complete_cols,
    verbose        = verbose
  )

  if (verbose) message("Data cleaning complete. Creating survey design object...")

  survey_design <- create_survey_object(
    data        = combined_df,
    id_var      = id_var,
    weight_var  = weight_var,
    strata_var  = strata_var,
    cluster_var = cluster_var,
    nest        = nest,
    multicore   = multicore,
    verbose     = verbose
  )

  # Save the survey design object as an RDS file in the FST directory.
  output_file <- file.path(fst_dir, "final_survey_object.rds")
  saveRDS(survey_design, output_file)
  if (verbose) message("Final survey design object saved to: ", output_file)

  survey_design
}
