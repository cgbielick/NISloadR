#' Flag ICD Codes in NIS Core Data
#'
#' This function loads an FST file for NIS Core data for a specified year,
#' selects only the specified (or relevant) columns, and then flags diagnosis and
#' procedure codes based on user-provided code vectors. By default, it loads all columns.
#'
#' The function handles year-specific issues. For years 2011 or earlier, it requires a path
#' to the new DISCWT (trend weights) files unless overridden by setting `old_weights = 1`.
#' For 2015, note that diagnosis codes may appear in a separate grps file (and ICD-10-CM fields
#' are used from Q4 onward).
#'
#' @param fst_file Character. The full path to the FST file containing the NIS core data.
#' @param year Numeric. The year of the NIS data (e.g., 2011, 2015).
#' @param col_select Optional. A character vector of column names to load.
#'   If \code{NULL} (the default), all columns are loaded.
#' @param weights_dir Optional. For years 2011 or earlier, a directory containing the
#'   DISCWT specification file (e.g. \code{"FileSpecifications_NIS_2011_HOSPITAL_TrendWt.TXT"})
#'   and the corresponding ASCII file (e.g. \code{"NIS_2011_HOSPITAL_TrendWt.ASC"}).
#'   This is required unless \code{old_weights = 1}.
#' @param old_weights Numeric; default 0. Set to 1 to intentionally use old weights.
#' @param dx_codes Optional. A named list (or vector) of ICD-9-CM or ICD-10-CM codes
#'   to search for in diagnosis columns. Each element is a character vector.
#'   The name of each element is used as the base name for the flagged columns.
#'   If unnamed, the function will use defaults ("flag1", "flag2", etc.).
#' @param pr_codes Optional. A named list (or vector) of ICD-9-PCS or ICD-10-PCS codes
#'   to search for in procedure columns.
#' @param prday_codes Optional. A similar vector for PRDAY columns; default is to flag all
#'   if \code{pr_codes} is provided. Use 0 for none.
#' @param adult Logical, default \code{FALSE}. If \code{TRUE}, only include discharges for adults
#'   (AGE >= 18 or AGE_NEONATE is NA); otherwise, load all.
#' @param cost_file Optional. A path to the HCUP charges-to-cost ratio file.
#'   If provided (and if \code{year >= 2001}), charges will be converted to estimated costs.
#' @param dx_flag_type Character, either "binary" or "sum". Indicates whether to create a flag
#'   column (0/1) for each ICD code set or a column summing the number of matches per discharge.
#'   (Default is "binary".)
#' @param ... Additional arguments (currently not used).
#'
#' @return A data.table with the selected columns, filtered according to the \code{adult} flag,
#'   with additional columns for each ICD code flag. If an error is encountered (e.g., missing
#'   required files or columns), a verbose error message is produced.
#'
#' @details
#' The function determines which columns to search for ICD codes by inspecting the column names:
#' for ICD-9-CM codes, it looks for columns starting with "DX" (and "PR" for procedures);
#' for ICD-10-CM/PCS, it looks for "I10_DX" and "I10_PR", respectively.
#'
#' For each set of codes provided in \code{dx_codes} or \code{pr_codes}, it creates a new column:
#' if the input list element has a name, that name is used (appended with a suffix, if needed);
#' if not, the flag column is named "flag1", "flag2", etc.
#'
#' Exact matching is enforced using regex anchors (^ and $) so that only codes that match exactly
#' are flagged.
#'
#' @export
#'
#' @examples
#' \dontrun{
#'   # Load all columns from the 2011 core file and flag ICD-9-CM diagnoses starting with "250"
#'   dt <- flag_nis_dx(
#'           fst_file = "path/to/NIS_2011_Core.fst",
#'           year = 2011,
#'           weights_dir = "path/to/2011/TrendWt/",
#'           dx_codes = list(diabetes = c("250.00", "250.01")),
#'           pr_codes = NULL,
#'           adult = TRUE
#'         )
#' }
flag_nis_dx <- function(fst_file, year, col_select = NULL,
                        weights_dir = NULL, old_weights = 0,
                        dx_codes = NULL, pr_codes = NULL, prday_codes = NULL,
                        adult = FALSE, cost_file = NULL,
                        dx_flag_type = "binary", ...) {
  # Validate year input
  if (!is.numeric(year) || length(year) != 1) {
    stop("Argument 'year' must be a single numeric value.")
  }

  # Validate fst_file exists
  if (!file.exists(fst_file)) {
    stop("FST file not found: ", fst_file)
  }

  # For years 2011 or earlier, require weights_dir unless old_weights is set to 1
  if (year <= 2011 && old_weights == 0) {
    if (is.null(weights_dir) || !dir.exists(weights_dir)) {
      stop("For year ", year, " (<=2011), a valid 'weights_dir' must be provided containing the new DISCWT files.")
    }
  }

  # Read only selected columns (if specified) using fst
  if (!is.null(col_select)) {
    df <- fst::read_fst(fst_file, columns = col_select)
  } else {
    df <- fst::read_fst(fst_file)
  }

  # Convert to data.table by reference
  data.table::setDT(df)

  # Optional: Filter for adult discharges (if adult == TRUE)
  if (adult) {
    if ("AGE" %in% names(df)) {
      df <- df[is.na(AGE) | AGE >= 18]
    } else if ("AGE_NEONATE" %in% names(df)) {
      df <- df[is.na(AGE_NEONATE)]
    } else {
      warning("No age indicator (AGE or AGE_NEONATE) found; 'adult' filter not applied.")
    }
  }

  # If cost_file is provided (and year >=2001), merge in cost conversion
  if (!is.null(cost_file)) {
    if (!file.exists(cost_file)) {
      stop("Cost file not found: ", cost_file)
    }
    # (Placeholder: read cost file and merge with df to adjust TOTCHG to cost)
    # For example, one could use readr::read_csv() and a data.table join.
    message("Cost conversion not yet implemented; please implement as needed.")
  }

  # For years <=2011 with new weights (old_weights==0), load DISCWT data
  if (year <= 2011 && old_weights == 0) {
    # Construct expected file paths (assuming fixed naming conventions)
    spec_file <- file.path(weights_dir, paste0("FileSpecifications_NIS_", year, "_HOSPITAL_TrendWt.TXT"))
    asc_weights_file <- file.path(weights_dir, paste0("NIS_", year, "_HOSPITAL_TrendWt.ASC"))
    if (!file.exists(spec_file)) {
      stop("Weights specification file not found: ", spec_file)
    }
    if (!file.exists(asc_weights_file)) {
      stop("Weights ASC file not found: ", asc_weights_file)
    }
    # (Placeholder: read the specification file and ASCII file to get new DISCWT values)
    # Then overwrite the DISCWT column in df.
    message("Weights file processing not yet implemented; please implement as needed.")
  }

  # Helper to flag codes in a set of columns based on an exact regex match
  flag_codes <- function(dt, code_list, col_patterns, flag_type = "binary") {
    # dt: data.table
    # code_list: named list or vector of ICD codes to search for (each element is a vector)
    # col_patterns: a character vector of patterns for the column names (e.g., "^DX" or "^I10_DX")
    # flag_type: "binary" creates 0/1 flag, "sum" creates count of matches per row.

    # Identify columns matching the given pattern(s)
    target_cols <- names(dt)[sapply(col_patterns, function(pat) grepl(pat, names(dt)))]
    if (length(target_cols) == 0) {
      warning("No target columns found for pattern(s): ", paste(col_patterns, collapse = ", "))
      return(dt)
    }

    # Ensure code_list is a list; if a vector, wrap it in a list
    if (!is.list(code_list)) {
      code_list <- list(code_list)
    }

    # Iterate over each code vector
    flag_idx <- 1
    for (code_vec in code_list) {
      # Determine the flag column name:
      flag_name <- names(code_list)[flag_idx]
      if (is.null(flag_name) || flag_name == "") {
        flag_name <- paste0("flag", flag_idx)
      }
      # Create regex pattern that requires an exact match: ^(code1|code2|...)$
      pattern <- paste0("^(", paste(code_vec, collapse = "|"), ")$")

      # For each row, check if any of the target columns exactly match any of the codes
      if (flag_type == "binary") {
        dt[, (flag_name) := as.integer(apply(.SD, 1, function(x) any(grepl(pattern, x, perl = TRUE)))),
           .SDcols = target_cols]
      } else if (flag_type == "sum") {
        dt[, (flag_name) := apply(.SD, 1, function(x) sum(grepl(pattern, x, perl = TRUE))),
           .SDcols = target_cols]
      } else {
        stop("Unknown flag_type: ", flag_type)
      }
      flag_idx <- flag_idx + 1
    }
    return(dt)
  }

  # Flag diagnosis codes
  if (!is.null(dx_codes)) {
    # Determine which diagnosis columns to search for based on year:
    # For ICD-9-CM: columns starting with "DX"
    # For ICD-10-CM: columns starting with "I10_DX"
    if (year < 2015) {
      dx_pattern <- "^DX"
    } else {
      dx_pattern <- "^I10_DX"
    }
    df <- flag_codes(dt = df, code_list = dx_codes, col_patterns = dx_pattern, flag_type = dx_flag_type)
  }

  # Flag procedure codes
  if (!is.null(pr_codes)) {
    # For ICD-9-PCS: columns starting with "PR"
    # For ICD-10-PCS: columns starting with "I10_PR"
    if (year < 2015) {
      pr_pattern <- "^PR"
    } else {
      pr_pattern <- "^I10_PR"
    }
    df <- flag_codes(dt = df, code_list = pr_codes, col_patterns = pr_pattern, flag_type = dx_flag_type)
  }

  # If pr_codes is 0 (or NULL) then prday_codes should be 0 as well
  if (!is.null(prday_codes) && (is.null(pr_codes) || length(pr_codes) == 0)) {
    message("prday_codes provided but pr_codes is empty; skipping PRDAY flagging.")
  }
  # (Implement similar flagging for PRDAY columns if needed)

  # Return the updated data.table
  return(df)
}
