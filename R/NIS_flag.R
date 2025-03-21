#' Flag ICD Codes in NIS Core Data
#'
#' This function processes National Inpatient Sample (NIS) FST files for a series of years,
#' applies user-specified flagging for ICD-9 and ICD-10 diagnosis (and optionally procedure) codes,
#' and writes out processed files. The function supports both binary (0/1) and summed flags.
#'
#' @param fst_dir Character. Directory containing the pre-processed FST files.
#' @param projectname Character. Name of the project; used to build relative output directories.
#' @param years Numeric vector. Years of data to process (e.g., 2011:2021).
#' @param dx_codes_9 Optional. A vector or list of ICD-9 diagnosis codes to flag.
#' @param dx_codes_10 Optional. A vector or list of ICD-10 diagnosis codes to flag.
#' @param pr_codes_9 Optional. A vector or list of ICD-9 procedure codes to flag.
#' @param pr_codes_10 Optional. A vector or list of ICD-10 procedure codes to flag.
#' @param flag_type Character. Either "binary" (default) or "sum". Determines whether the flag is 0/1 or a count.
#' @param old_weights Numeric. Default 0. Set to 1 to use old weights even for years <= 2011.
#' @param weights_dir Optional. Directory containing DISCWT files for years <= 2011 (if old_weights == 0).
#' @param discwt_data Optional. A tibble with columns YEAR, HOSPID_final, and DISCWT for years <= 2011.
#' @param adult Logical. Default TRUE. If TRUE, only processes adult discharges (AGE > 17).
#' @param dry_run Logical. If TRUE, processing is simulated and no output files are written.
#' @param num_cores Numeric. Number of cores for parallel processing. Defaults to max(1, parallel::detectCores() - 1).
#' @param import_cols Character vector. Names of columns to import from each FST file.
#'        Default is a set of key columns and diagnosis/procedure columns.
#' @param rename_mapping List. A named list where each element is a vector of possible raw column names
#'        that should be standardized. (Note: FEMALE is no longer remapped to SEX.)
#'
#' @return Invisibly returns a list with processing summaries per year.
#'
#' @examples
#' \dontrun{
#'   # Assuming you have loaded discwt_data from the appropriate DISCWT file
#'   results <- flag_nis(
#'     fst_dir = "~/NIS_LARGE/Data/FST",
#'     projectname = "NIS_HIV",
#'     years = 2011:2021,
#'     dx_codes_9 = c("042", "07953"),
#'     dx_codes_10 = c("B20", "B9735"),
#'     pr_codes_9 = c("88.72", "39.61"),
#'     pr_codes_10 = c("0DTJ0ZZ", "0DTJ4ZZ"),
#'     flag_type = "binary",
#'     old_weights = 0,
#'     weights_dir = "~/NIS_LARGE/Data/Weights",
#'     discwt_data = discwt_data,  # pre-loaded DISCWT tibble
#'     adult = TRUE,
#'     dry_run = FALSE,
#'     num_cores = 4,
#'     import_cols = c("AGE", "HOSP_NIS", "KEY_NIS", "FEMALE", "HOSPID", "KEY",
#'                     "PL_NCHS", "NIS_STRATUM", "RACE", "YEAR", "ZIPINC_QRTL", "DIED",
#'                     "DISCWT", "PAY1",
#'                     paste0("DX", 1:25), paste0("I10_DX", 1:25),
#'                     paste0("PR", 1:25), paste0("I10_PR", 1:25)),
#'     rename_mapping = list(
#'       HOSPID_final = c("HOSP_NIS", "HOSPID"),
#'       KEY_final = c("KEY_NIS", "KEY"),
#'       PL_NCHS_final = c("PL_NCHS", "PL_NCHS2006")
#'     )
#'   )
#' }
#'
#' @export
flag_nis <- function(fst_dir,
                     projectname,
                     years,
                     dx_codes_9 = NULL,
                     dx_codes_10 = NULL,
                     pr_codes_9 = NULL,
                     pr_codes_10 = NULL,
                     flag_type = "binary",
                     old_weights = 0,
                     weights_dir = NULL,
                     discwt_data = NULL,
                     adult = TRUE,
                     dry_run = FALSE,
                     num_cores = max(1, parallel::detectCores() - 1),
                     import_cols = c("AGE", "HOSP_NIS", "KEY_NIS", "FEMALE", "HOSPID", "KEY",
                                     "PL_NCHS", "NIS_STRATUM", "RACE", "YEAR", "ZIPINC_QRTL", "DIED",
                                     "DISCWT", "PAY1",
                                     paste0("DX", 1:25), paste0("I10_DX", 1:25),
                                     paste0("PR", 1:25), paste0("I10_PR", 1:25)),
                     rename_mapping = list(
                       HOSPID_final = c("HOSP_NIS", "HOSPID"),
                       KEY_final = c("KEY_NIS", "KEY"),
                       PL_NCHS_final = c("PL_NCHS", "PL_NCHS2006")
                     )
) {

  requireNamespace("fst", quietly = TRUE)
  requireNamespace("data.table", quietly = TRUE)
  requireNamespace("dplyr", quietly = TRUE)
  requireNamespace("foreach", quietly = TRUE)
  requireNamespace("doParallel", quietly = TRUE)

  output_dir <- file.path("Projects", projectname, "output")
  if (!dir.exists(output_dir)) {
    message("Creating output directory: ", output_dir)
    dir.create(output_dir, recursive = TRUE)
  }

  # Standardize columns using provided mapping
  standardize_columns <- function(df, mapping) {
    for (std_name in names(mapping)) {
      raw_names <- mapping[[std_name]]
      present_names <- intersect(raw_names, names(df))
      if (length(present_names) > 0 && !(std_name %in% names(df))) {
        names(df)[names(df) == present_names[1]] <- std_name
      }
    }
    return(df)
  }

  # Safe read FST: read only available columns
  safe_read_fst <- function(file, import_cols) {
    df_header <- fst::read_fst(file, from = 1, to = 1)
    available_cols <- colnames(df_header)
    cols_to_read <- intersect(import_cols, available_cols)
    fst::read_fst(file, columns = cols_to_read)
  }

  # Process covariates (e.g., age groups)
  process_covariates <- function(df, current_year) {
    data.table::setDT(df)
    df[, YEAR_ALL := "all"]
    message("Processing covariates for year ", current_year, ".")

    df[, AGE_grp := data.table::fifelse(AGE <= 24, 1L,
                                        data.table::fifelse(AGE >= 25 & AGE <= 34, 2L,
                                                            data.table::fifelse(AGE >= 35 & AGE <= 44, 3L,
                                                                                data.table::fifelse(AGE >= 45 & AGE <= 54, 4L, 5L))))]
    df[, AGE_grp := factor(AGE_grp, levels = 1:5,
                           labels = c("18-24", "25-34", "35-44", "45-54", "55+"))]
    df[, AGE_grp := relevel(AGE_grp, ref = "45-54")]

    # Note: FEMALE column remains unchanged.
    message("Finished processing covariates for year ", current_year, ".")
    return(df)
  }

  # Flag codes using exact match with %in%
  flag_codes <- function(df_chunk, target_cols, codes_list, flag_type) {
    for (flag in names(codes_list)) {
      if (!(flag %in% names(df_chunk)))
        df_chunk[[flag]] <- 0L
    }
    valid_cols <- intersect(target_cols, names(df_chunk))
    if (length(valid_cols) > 0) {
      for (flag in names(codes_list)) {
        code_vector <- codes_list[[flag]]
        if (flag_type == "binary") {
          match_found <- rowSums(sapply(valid_cols, function(col) {
            !is.na(df_chunk[[col]]) & (df_chunk[[col]] %in% code_vector)
          })) > 0
          df_chunk[[flag]] <- as.integer(match_found)
        } else if (flag_type == "sum") {
          df_chunk[[flag]] <- rowSums(sapply(valid_cols, function(col) {
            as.integer(!is.na(df_chunk[[col]]) & (df_chunk[[col]] %in% code_vector))
          }))
        } else {
          stop("Unknown flag_type: ", flag_type)
        }
      }
    }
    return(df_chunk)
  }

  # Read data for a given year with filtering and DISCWT merging for years <= 2011
  read_year_data <- function(current_year, fst_dir, import_cols, rename_mapping, discwt_data, adult) {
    message("Loading data for year ", current_year, ".")
    file_name <- paste0("nis_core_", current_year, ".fst")
    full_file <- file.path(fst_dir, file_name)
    if (!file.exists(full_file)) {
      stop("FST file not found for year ", current_year, ": ", full_file)
    }
    df <- safe_read_fst(full_file, import_cols)
    df <- standardize_columns(df, rename_mapping)
    if (adult && "AGE" %in% names(df)) {
      df <- dplyr::filter(df, AGE > 17)
    }
    # For years <= 2011, enforce DISCWT data presence unless overridden by old_weights = 1
    if (current_year <= 2011 && old_weights == 0) {
      if (is.null(discwt_data) || !(current_year %in% unique(discwt_data$YEAR))) {
        stop("For year ", current_year, ", DISCWT data for altered weights must be provided.")
      }
      message("Merging DISCWT data for year ", current_year, ".")
      df$HOSPID_final <- as.numeric(df$HOSP_NIS)
      df <- dplyr::left_join(df, discwt_data, by = c("YEAR", "HOSPID_final"))
      df$DISCWT <- as.numeric(df$DISCWT)
    }
    message("Data for year ", current_year, " loaded with ", nrow(df), " rows.")
    return(df)
  }

  # Process a single year
  process_year <- function(current_year, fst_dir, output_dir,
                           import_cols, rename_mapping, discwt_data, adult,
                           dx_codes_9, dx_codes_10, pr_codes_9, pr_codes_10,
                           flag_type, num_cores, dry_run) {
    overall_start <- Sys.time()
    message("--------------------------------------------------")
    message("Processing year: ", current_year)

    tryCatch({
      full_data <- read_year_data(current_year, fst_dir, import_cols, rename_mapping, discwt_data, adult)

      full_data <- process_covariates(full_data, current_year)

      # Identify diagnosis and procedure columns exactly
      dx9_cols <- names(full_data)[grepl("^DX[0-9]+$", names(full_data))]
      dx10_cols <- names(full_data)[grepl("^I10_DX[0-9]+$", names(full_data))]
      pr9_cols <- names(full_data)[grepl("^PR[0-9]+$", names(full_data))]
      pr10_cols <- names(full_data)[grepl("^I10_PR[0-9]+$", names(full_data))]

      # Build code lists from user input
      dx_list <- list()
      if (!is.null(dx_codes_9)) dx_list[["dx9"]] <- dx_codes_9
      if (!is.null(dx_codes_10)) dx_list[["dx10"]] <- dx_codes_10
      pr_list <- list()
      if (!is.null(pr_codes_9)) pr_list[["pr9"]] <- pr_codes_9
      if (!is.null(pr_codes_10)) pr_list[["pr10"]] <- pr_codes_10

      total_rows <- nrow(full_data)
      chunk_size <- ceiling(total_rows / num_cores)
      df_chunks <- split(full_data, ceiling(seq_len(total_rows) / chunk_size))

      cl <- parallel::makeCluster(num_cores)
      doParallel::registerDoParallel(cl)

      flagged_chunks <- foreach::foreach(chunk = df_chunks,
                                         .packages = c("data.table", "dplyr"),
                                         .export = c("flag_codes")) %dopar% {
                                           if (length(dx_list) > 0) {
                                             if ("dx9" %in% names(dx_list)) {
                                               chunk <- flag_codes(chunk, dx9_cols, list(dx9 = dx_list[["dx9"]]), flag_type)
                                             }
                                             if ("dx10" %in% names(dx_list)) {
                                               chunk <- flag_codes(chunk, dx10_cols, list(dx10 = dx_list[["dx10"]]), flag_type)
                                             }
                                           }
                                           if (length(pr_list) > 0) {
                                             if ("pr9" %in% names(pr_list)) {
                                               chunk <- flag_codes(chunk, pr9_cols, list(pr9 = pr_list[["pr9"]]), flag_type)
                                             }
                                             if ("pr10" %in% names(pr_list)) {
                                               chunk <- flag_codes(chunk, pr10_cols, list(pr10 = pr_list[["pr10"]]), flag_type)
                                             }
                                           }
                                           chunk
                                         }
      parallel::stopCluster(cl)

      processed_data <- dplyr::bind_rows(flagged_chunks)

      # Placeholder: Future functionality (cost conversion, survey object conversion, etc.)

      output_file <- file.path(output_dir, paste0("core_", current_year, "_processed.fst"))
      if (dry_run) {
        message("[DRY RUN] Processed data for year ", current_year,
                " would be written to ", output_file)
      } else {
        message("Writing processed data for year ", current_year, " to ", output_file)
        fst::write_fst(processed_data, output_file, compress = 100)
      }
      message("Year ", current_year, " processing complete.")
    }, error = function(e) {
      message("Error processing year ", current_year, ": ", e$message)
    })

    elapsed <- Sys.time() - overall_start
    message("Year ", current_year, " elapsed time: ", round(as.numeric(elapsed, units = "secs"), 2), " seconds")
  }

  # Loop over years
  for (yr in years) {
    process_year(yr, fst_dir, output_dir, import_cols, rename_mapping, discwt_data,
                 adult, dx_codes_9, dx_codes_10, pr_codes_9, pr_codes_10,
                 flag_type, num_cores, dry_run)
  }

  message("All processing complete. Parallel clusters stopped (if any).")
  invisible(NULL)
}
