suppressPackageStartupMessages({
  library(tidyverse)  # Data manipulation, purrr, etc.
  library(fst)        # Fast reading of fst files with column selection
  library(srvyr)      # Tidy survey design object
  library(rlang)      # For tidy evaluation
})

# -----------------------------------------------------------------------------
# Function: load_and_clean_fst
#
# Reads each .fst file in a directory (one per year), selecting only the columns
# the user wants. By default the DX columns (and similar) are dropped unless
# explicitly included via include_cols.
#
# The function then applies:
#   - an optional subset filter (subset_expr)
#   - a complete-cases filter on specified columns (complete_cols)
#
# Parameters:
#   fst_dir       : Directory containing the .fst files.
#   include_cols  : Character vector of column names to import. If NULL,
#                   then all columns are imported, and later heavy columns are dropped.
#                   (Default behavior excludes DX columns.)
#   drop_patterns : Default regex patterns to drop (if include_cols is NULL).
#                   Defaults to diagnosis columns.
#   additional_drop: Optional additional regex patterns to drop.
#   subset_expr   : Optional unquoted filtering expression (e.g., PWH == 1).
#   complete_cols : Character vector of column names that must have no missing values.
#                   Rows with NA in any of these columns are dropped.
#   verbose       : Logical; if TRUE, prints progress messages.
#
# Returns:
#   A tibble containing the combined and cleaned data from all files.
# -----------------------------------------------------------------------------
load_and_clean_fst <- function(fst_dir,
                               include_cols  = NULL,
                               drop_patterns = c("^DX", "^I10_DX", "^PR", "^I10_PR"),
                               additional_drop = NULL,
                               subset_expr   = NULL,
                               complete_cols = c("FEMALE", "AGE_grp", "RACE", "DIED"),
                               verbose       = TRUE) {
  # List all .fst files in the directory
  fst_files <- list.files(fst_dir, pattern = "\\.fst$", full.names = TRUE)
  if (length(fst_files) == 0) {
    stop("No FST files found in directory: ", fst_dir)
  }
  if (verbose) message("Found ", length(fst_files), " FST file(s) in ", fst_dir)

  # Prepare a list to hold cleaned data from each file
  cleaned_list <- vector("list", length(fst_files))

  # Determine which columns to read:
  # If include_cols is provided, we read only those columns.
  # Otherwise, we read all columns then drop heavy ones matching drop_patterns.
  for (i in seq_along(fst_files)) {
    file <- fst_files[i]
    if (verbose) message("Processing file: ", file)

    # Read the file: if include_cols is specified, use that to limit columns.
    # This avoids loading unwanted columns into memory.
    df <- tryCatch({
      if (!is.null(include_cols)) {
        as_tibble(fst::read_fst(file, columns = include_cols))
      } else {
        as_tibble(fst::read_fst(file))
      }
    }, error = function(e) {
      warning("Failed to read file: ", file, " - ", e$message)
      return(NULL)
    })

    if (is.null(df)) next

    # If include_cols is NULL, then drop columns matching heavy drop patterns.
    if (is.null(include_cols)) {
      all_drop_patterns <- c(drop_patterns, additional_drop)
      drop_cols <- names(df)[
        map_lgl(names(df), ~ any(str_detect(.x, all_drop_patterns)))
      ]
      if (verbose && length(drop_cols) > 0) {
        message("Dropping columns: ", paste(drop_cols, collapse = ", "))
      }
      df <- df %>% select(-all_of(drop_cols))
    }

    # Apply subset filter if provided
    if (!is.null(subset_expr)) {
      sub_expr <- enquo(subset_expr)
      df <- df %>% filter(!!sub_expr)
      if (verbose) message("After subsetting, rows: ", nrow(df))
    }

    # Drop rows with missing values in the required complete_cols (if they exist)
    existing_complete <- intersect(complete_cols, names(df))
    if (length(existing_complete) > 0) {
      before <- nrow(df)
      df <- df %>% drop_na(all_of(existing_complete))
      if (verbose) message("Dropped ", before - nrow(df), " rows missing required columns: ",
                           paste(existing_complete, collapse = ", "))
    }

    # Optionally, you could add a column indicating the source (e.g., year)
    # For instance, if the file name contains the year, extract it:
    year_extracted <- str_extract(basename(file), "\\d{4}")
    if (!is.na(year_extracted)) {
      df <- df %>% mutate(year = as.integer(year_extracted))
    }

    cleaned_list[[i]] <- df

    # Clean up memory explicitly if desired
    rm(df); gc()
  }

  # Remove any NULLs from failed reads
  cleaned_list <- compact(cleaned_list)
  if (length(cleaned_list) == 0) {
    stop("No valid data loaded from FST files.")
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
                                 multicore   = FALSE,
                                 verbose     = TRUE) {
  # Warn if multicore is requested
  if (multicore) {
    cores <- parallel::detectCores(logical = TRUE)
    message("Multicore processing enabled. Detected ", cores, " logical cores.")
    # Optionally, limit the number of cores if too many (for safety)
    if (cores > 8) {
      warning("Using a high number of cores (", cores, ") may overwhelm your system. Consider reducing multicore usage.")
    }
  }

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
# Main wrapper: processes each FST file in a directory individually,
# cleans them (by selecting desired columns, subsetting, and enforcing complete cases),
# combines all years of data, creates a survey design object, and saves the final
# survey object to disk in the same directory.
#
# Parameters:
#   fst_dir         : Directory with FST files.
#   include_cols    : Character vector of column names to import.
#                     If NULL, all columns are imported and heavy DX columns dropped.
#   drop_patterns   : Default regex patterns for diagnosis columns to drop
#                     if include_cols is not provided.
#   additional_drop : Additional regex patterns to drop.
#   subset_expr     : Optional subsetting filter (e.g., PWH == 1).
#   complete_cols   : Columns that must have no missing data.
#   id_var, weight_var, strata_var, cluster_var : Survey design variable names.
#   nest            : Logical; passed to as_survey_design().
#   multicore       : Logical; if TRUE, uses multicore processing.
#                     (Default is FALSE for safety.)
#   verbose         : Logical; print progress messages.
#
# Returns:
#   A survey design object.
# -----------------------------------------------------------------------------
process_fst_to_survey <- function(fst_dir,
                                  include_cols    = NULL,
                                  drop_patterns   = c("^DX", "^I10_DX", "^PR", "^I10_PR"),
                                  additional_drop = NULL,
                                  subset_expr     = NULL,
                                  complete_cols   = c("FEMALE", "AGE_grp", "RACE", "DIED"),
                                  id_var          = "HOSPID_final",
                                  weight_var      = "DISCWT",
                                  strata_var      = "NIS_STRATUM",
                                  cluster_var     = "HOSP_NIS",
                                  nest            = TRUE,
                                  multicore       = FALSE,
                                  verbose         = TRUE) {
  if (verbose) message("Starting processing of FST files in directory: ", fst_dir)

  # Process each file one by one and combine cleaned data
  combined_df <- load_and_clean_fst(
    fst_dir       = fst_dir,
    include_cols  = include_cols,
    drop_patterns = drop_patterns,
    additional_drop = additional_drop,
    subset_expr   = subset_expr,
    complete_cols = complete_cols,
    verbose       = verbose
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

  # Save the final survey object as an RDS file in the same directory as the FST files.
  output_file <- file.path(fst_dir, "final_survey_object.rds")
  saveRDS(survey_design, output_file)
  if (verbose) message("Final survey design object saved to: ", output_file)

  survey_design
}

# -----------------------------------------------------------------------------
# Example usage:
#
# Suppose you have a directory "data/fst_files" containing one FST file per year.
# You want to include only columns "HOSPID_final", "DISCWT", "NIS_STRATUM", "HOSP_NIS",
# plus additional variables "FEMALE", "AGE_grp", "RACE", "DIED", "PWH", and maybe others.
#
# You also want to subset the data to only include PWH == 1, and enforce that
# FEMALE, AGE_grp, RACE, and DIED have no missing values.
#
# Finally, you decide not to use multicore processing by default.
#
# Run:
#
# my_survey <- process_fst_to_survey(
#   fst_dir         = "data/fst_files",
#   include_cols    = c("HOSPID_final", "DISCWT", "NIS_STRATUM", "HOSP_NIS",
#                       "FEMALE", "AGE_grp", "RACE", "DIED", "PWH", "OtherVar1", "OtherVar2"),
#   subset_expr     = PWH == 1,
#   complete_cols   = c("FEMALE", "AGE_grp", "RACE", "DIED"),
#   id_var          = "HOSPID_final",
#   weight_var      = "DISCWT",
#   strata_var      = "NIS_STRATUM",
#   cluster_var     = "HOSP_NIS",
#   nest            = TRUE,
#   multicore       = FALSE,   # Default is FALSE to avoid excessive resource usage
#   verbose         = TRUE
# )
#
# The resulting survey design object will be saved as "final_survey_object.rds"
# in the directory "data/fst_files".
# -----------------------------------------------------------------------------
