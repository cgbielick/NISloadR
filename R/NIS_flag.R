# =============================================================================
# Flag ICD Codes in NIS Core Data
#
# This function processes National Inpatient Sample (NIS) FST files for a series
# of years, applies user-specified flagging for ICD-9/ICD-10 (and optionally
# procedure) codes, merges DISCWT weights for years ≤ 2011, and writes out
# processed files. It preserves all covariate processing (including nuanced handling
# for 2015, REGION, URBAN, AGE groups, SEX, RACE, and PAY1 recoding) from the prior
# version. Users must supply the directory for DISCWT files (the new ASC files)
# using a naming convention flexible for each year ≤ 2011.
#
# =============================================================================
flag_nis <- function(fst_dir,
                     projectname,
                     years,
                     dx_codes_9 = NULL,
                     dx_codes_10 = NULL,
                     pr_codes_9 = NULL,
                     pr_codes_10 = NULL,
                     flag_type = "binary",         # "binary" (0/1) or "sum" (counts)
                     old_weights = 0,              # Set to 1 to use old weights even for years ≤ 2011
                     weights_dir = NULL,           # Directory for DISCWT ASC files
                     discwt_data = NULL,           # Optional pre-loaded tibble for DISCWT merging
                     adult = TRUE,                 # Only process records with AGE > 17 if TRUE
                     dry_run = FALSE,              # If TRUE, do not write output files
                     num_cores = max(1, parallel::detectCores() - 1),
                     icd_version = "both",         # "9", "10", or "both"
                     import_cols = c("AGE", "HOSP_NIS", "KEY_NIS", "FEMALE",
                                     "HOSPID", "KEY", "PL_NCHS", "NIS_STRATUM",
                                     "RACE", "YEAR", "ZIPINC_QRTL", "DIED",
                                     "DISCWT", "PAY1",
                                     paste0("DX", 1:25), paste0("I10_DX", 1:25),
                                     paste0("PR", 1:25), paste0("I10_PR", 1:25)),
                     rename_mapping = list(
                       HOSPID_final = c("HOSP_NIS", "HOSPID"),
                       KEY_final = c("KEY_NIS", "KEY"),
                       PL_NCHS_final = c("PL_NCHS", "PL_NCHS2006")
                     )
) {
  # Load required libraries quietly
  suppressPackageStartupMessages({
    library(fst)
    library(data.table)
    library(dplyr)
    library(foreach)
    library(doParallel)
    library(readr)
  })

  # Validate icd_version argument
  if(!icd_version %in% c("9", "10", "both")) {
    stop("Error: icd_version must be '9', '10', or 'both'.")
  }

  # Validate dx_codes_9 and dx_codes_10 inputs
  if (!is.null(dx_codes_9)) {
    if (!is.list(dx_codes_9)) {
      stop("Error: dx_codes_9 must be a list of character vectors (e.g., list(HIV = c(...), PJP = c(...))).")
    }
    for (flag in names(dx_codes_9)) {
      if (!is.character(dx_codes_9[[flag]])) {
        stop("Error: Element '", flag, "' in dx_codes_9 is not a character vector.")
      }
      if (length(dx_codes_9[[flag]]) == 0) {
        stop("Error: Element '", flag, "' in dx_codes_9 is empty.")
      }
    }
  }

  if (!is.null(dx_codes_10)) {
    if (!is.list(dx_codes_10)) {
      stop("Error: dx_codes_10 must be a list of character vectors (e.g., list(HIV = c(...), PJP = c(...))).")
    }
    for (flag in names(dx_codes_10)) {
      if (!is.character(dx_codes_10[[flag]])) {
        stop("Error: Element '", flag, "' in dx_codes_10 is not a character vector.")
      }
      if (length(dx_codes_10[[flag]]) == 0) {
        stop("Error: Element '", flag, "' in dx_codes_10 is empty.")
      }
    }
  }

  # Create output directory if needed
  output_dir <- file.path("Projects", projectname, "output")
  if (!dir.exists(output_dir)) {
    message("Creating output directory: ", output_dir)
    dir.create(output_dir, recursive = TRUE)
  }

  # ---------------------------------------------------------------------------
  # Helper Function: Standardize Column Names Based on Mapping
  # ---------------------------------------------------------------------------
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

  # ---------------------------------------------------------------------------
  # Helper Function: Safe FST Read (read only available columns)
  # ---------------------------------------------------------------------------
  safe_read_fst <- function(file, import_cols) {
    df_header <- fst::read_fst(file, from = 1, to = 1)
    available_cols <- colnames(df_header)
    cols_to_read <- intersect(import_cols, available_cols)
    fst::read_fst(file, columns = cols_to_read)
  }

  # ---------------------------------------------------------------------------
  # Helper Function: Process Covariates (includes AGE, SEX, RACE, REGION, URBAN, etc.)
  # ---------------------------------------------------------------------------
  process_covariates <- function(df, current_year) {
    setDT(df)
    df[, YEAR_ALL := "all"]
    message("Processing covariates for year ", current_year, ".")

    # AGE Grouping
    df[, AGE_grp := fifelse(AGE <= 24, 1L,
                            fifelse(AGE >= 25 & AGE <= 34, 2L,
                                    fifelse(AGE >= 35 & AGE <= 44, 3L,
                                            fifelse(AGE >= 45 & AGE <= 54, 4L, 5L))))]
    df[, AGE_grp := factor(AGE_grp, levels = 1:5,
                           labels = c("18-24", "25-34", "35-44", "45-54", "55+"))]
    df[, AGE_grp := relevel(AGE_grp, ref = "45-54")]

    # SEX Processing
    if ("SEX" %in% names(df)) {
      df[, SEX := factor(SEX, levels = c(0, 1), labels = c("Male", "Female"))]
      df[, SEX := relevel(SEX, ref = "Male")]
    }

    # RACE Processing
    if ("RACE" %in% names(df)) {
      df[, RACE := factor(RACE,
                          levels = c(1, 2, 3, 4, 5, 6),
                          labels = c("White", "Black", "Hispanic",
                                     "Asian American or Pacific Islander", "Other", "Other"))]
      df[, RACE := relevel(RACE, ref = "White")]
    }

    # REGION & URBAN processing via NIS_STRATUM
    if ("NIS_STRATUM" %in% names(df)) {
      df[, ns := as.character(NIS_STRATUM)]
      df[, is_2011minus := (YEAR <= 2011)]
      df[, is_2012plus := (YEAR >= 2012)]
      df[, REGION_num := fcase(
        is_2011minus & startsWith(ns, "1"), 1L,
        is_2011minus & startsWith(ns, "2"), 2L,
        is_2011minus & startsWith(ns, "3"), 3L,
        is_2011minus & startsWith(ns, "4"), 4L,
        is_2012plus & startsWith(ns, "1"), 1L,
        is_2012plus & startsWith(ns, "2"), 1L,
        is_2012plus & startsWith(ns, "3"), 2L,
        is_2012plus & startsWith(ns, "4"), 2L,
        is_2012plus & startsWith(ns, "5"), 3L,
        is_2012plus & startsWith(ns, "6"), 3L,
        is_2012plus & startsWith(ns, "7"), 3L,
        is_2012plus & startsWith(ns, "8"), 4L,
        is_2012plus & startsWith(ns, "9"), 4L,
        default = NA_integer_
      )]
      df[, REGION := factor(REGION_num, levels = c(1L, 2L, 3L, 4L),
                            labels = c("Northeast", "Midwest", "South", "West"))]
      df[, REGION := relevel(REGION, ref = "Northeast")]

      df[, URBAN_val := fcase(
        nchar(ns) == 4 & substring(ns, 3, 3) == "1", 0L,
        nchar(ns) == 4 & substring(ns, 3, 3) %in% c("2", "3"), 1L,
        default = NA_integer_
      )]
      df[, URBAN := factor(URBAN_val, levels = c(0, 1),
                           labels = c("Rural", "Urban"))]
      df[, URBAN := relevel(URBAN, ref = "Urban")]
    }

    # DIED and PAY1 Processing
    if ("DIED" %in% names(df)) df[, DIED := as.numeric(DIED)]
    if ("PAY1" %in% names(df)) {
      df[, PAY1_recoded := fifelse(PAY1 == 6, 5L, PAY1)]
      df[, PAY1 := factor(PAY1_recoded,
                          levels = c(1, 2, 3, 4, 5),
                          labels = c("Medicare", "Medicaid", "Private Insurance",
                                     "Self-pay", "Other"))]
      df[, PAY1 := relevel(PAY1, ref = "Medicaid")]
      df[, PAY1_recoded := NULL]
    }

    # ZIPINC_QRTL Processing
    if ("ZIPINC_QRTL" %in% names(df)) {
      df[, ZIPINC_QRTL := factor(ZIPINC_QRTL,
                                 levels = c(1, 2, 3, 4),
                                 labels = c("0-25th Percentile",
                                            "26th-50th Percentile",
                                            "51st to 75th Percentile",
                                            "76th to 100th Percentile"))]
      df[, ZIPINC_QRTL := relevel(ZIPINC_QRTL, ref = "0-25th Percentile")]
    }

    # PL_NCHS_final Processing
    if ("PL_NCHS_final" %in% names(df)) {
      df[, PL_NCHS_final := fifelse(PL_NCHS_final %in% c(1, 2, 3, 4), 1L,
                                    fifelse(PL_NCHS_final %in% c(5, 6), 2L, PL_NCHS_final))]
    }

    df[, c("URBAN_val", "ns", "is_2011", "is_2012_2018", "REGION_num") := NULL]
    message("Finished processing covariates for year ", current_year, ".")
    return(df)
  }

  # ---------------------------------------------------------------------------
  # Helper Function: Flag Codes Using Exact %in% Matching
  # ---------------------------------------------------------------------------
  flag_codes <- function(df_chunk, target_cols, codes_list, flag_type, current_year = NA) {
    if (!is.list(codes_list)) {
      stop("Error in flag_codes: codes_list must be a list.")
    }
    for (flag in names(codes_list)) {
      codes <- codes_list[[flag]]
      if (!is.character(codes) || length(codes) == 0) {
        stop("Error in flag_codes: Code vector for flag '", flag, "' must be a non-empty character vector.")
      }
    }

    message("Starting flagging for year ", current_year, " with ", length(codes_list), " flag(s).")
    valid_cols <- intersect(target_cols, names(df_chunk))
    if (length(valid_cols) > 0) {
      for (flag in names(codes_list)) {
        code_vector <- codes_list[[flag]]
        if (flag_type == "binary") {
          match_found <- Reduce("|", lapply(valid_cols, function(col) {
            df_chunk[[col]] %in% code_vector
          }))
          df_chunk[[flag]] <- as.integer(match_found)
        } else if (flag_type == "sum") {
          df_chunk[[flag]] <- rowSums(sapply(valid_cols, function(col) {
            as.integer(df_chunk[[col]] %in% code_vector)
          }))
        } else {
          stop("Unknown flag_type: ", flag_type)
        }
      }
    }
    return(df_chunk)
  }


  # ---------------------------------------------------------------------------
  # Helper Function: Read Data for a Given Year (with special handling for 2015
  # and DISCWT merging for years ≤ 2011)
  # ---------------------------------------------------------------------------
  # Read data for a given year with filtering and DISCWT merging for years <= 2011

  read_year_data <- function(current_year, fst_dir, import_cols, rename_mapping,
                             old_weights, weights_dir, discwt_data, adult) {
    message("Loading data for year ", current_year, ".")
    if (current_year == 2015) {
      file_q1q3 <- file.path(fst_dir, "nis_core_2015q1q3.fst")
      file_q4   <- file.path(fst_dir, "nis_core_2015q4.fst")
      if (!file.exists(file_q1q3)) stop("File not found: ", file_q1q3)
      if (!file.exists(file_q4)) stop("File not found: ", file_q4)
      data_q1q3 <- safe_read_fst(file_q1q3, import_cols)
      data_q4   <- safe_read_fst(file_q4, import_cols)
      data_q1q3 <- standardize_columns(data_q1q3, rename_mapping)
      data_q4   <- standardize_columns(data_q4, rename_mapping)
      if (adult && "AGE" %in% names(data_q1q3)) {
        data_q1q3 <- dplyr::filter(data_q1q3, AGE > 17)
        data_q4   <- dplyr::filter(data_q4, AGE > 17)
      }
      df <- bind_rows(data_q1q3, data_q4)
    } else {
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
    }

    if (current_year <= 2011 && old_weights == 0) {
      if (is.null(discwt_data)) {
        if (is.null(weights_dir)) {
          stop("weights_dir must be provided for years ≤ 2011 when old_weights == 0")
        }
        weights_file <- file.path(weights_dir, paste0("NIS_", current_year, "_HOSPITAL_TrendWt.ASC"))
        if (!file.exists(weights_file)) {
          stop("DISCWT file for year ", current_year, " not found in weights_dir: ", weights_file)
        }
        discwt_data <- readr::read_csv(weights_file, col_names = FALSE, trim_ws = FALSE) %>%
          mutate(
            YEAR = as.numeric(substr(X1, 1, 4)),
            HOSPID_final = as.numeric(substr(X1, 5, 9)),
            DISCWT = as.numeric(substr(X1, 11, 19))
          ) %>%
          select(-X1)
      }
      if (!("HOSPID_final" %in% names(df))) {
        stop("HOSPID_final column not found after standardization.")
      }
      df$HOSPID_final <- as.numeric(df$HOSPID_final)
      message("Merging DISCWT data for year ", current_year, ".")
      df <- left_join((df %>% select(-DISCWT)), discwt_data, by = c("YEAR", "HOSPID_final"))
      if (!("DISCWT" %in% names(df)) || length(df$DISCWT) == 0) {
        message("Warning: DISCWT column not found after join. Creating DISCWT column with NA values.")
        df$DISCWT <- NA_real_
      } else {
        df$DISCWT <- as.numeric(df$DISCWT)
      }
    }

    message("Data for year ", current_year, " loaded with ", nrow(df), " rows.")
    return(df)
  }

  # ---------------------------------------------------------------------------
  # Helper Function: Process a Single Year (with parallel flagging)
  # ---------------------------------------------------------------------------
  process_year <- function(current_year, fst_dir, output_dir, import_cols, rename_mapping,
                           discwt_data, adult, dx_codes_9, dx_codes_10, pr_codes_9, pr_codes_10,
                           flag_type, old_weights, weights_dir, num_cores, dry_run, icd_version) {
    overall_start <- Sys.time()
    message("--------------------------------------------------")
    message("Processing year: ", current_year)

    tryCatch({
      full_data <- read_year_data(current_year, fst_dir, import_cols, rename_mapping,
                                  old_weights, weights_dir, discwt_data, adult)

      full_data <- process_covariates(full_data, current_year)

      dx9_cols <- names(full_data)[grepl("^DX[0-9]+$", names(full_data))]
      dx10_cols <- names(full_data)[grepl("^I10_DX[0-9]+$", names(full_data))]
      pr9_cols <- names(full_data)[grepl("^PR[0-9]+$", names(full_data))]
      pr10_cols <- names(full_data)[grepl("^I10_PR[0-9]+$", names(full_data))]

      total_rows <- nrow(full_data)
      chunk_size <- ceiling(total_rows / num_cores)
      df_chunks <- split(full_data, ceiling(seq_len(total_rows) / chunk_size))

      message("Starting ICD code flagging for year ", current_year, ".")

      cl <- makeCluster(num_cores)
      registerDoParallel(cl)

      flagged_chunks <- foreach(chunk = df_chunks,
                                .packages = c("data.table", "dplyr"),
                                .export = c("flag_codes")) %dopar% {
                                  # Flag using ICD-9 codes if requested or if combining both
                                  if (icd_version %in% c("9", "both")) {
                                    if (!is.null(dx_codes_9)) {
                                      for (flag in names(dx_codes_9)) {
                                        colname <- if (icd_version == "9") flag else paste0(flag, "_9")
                                        chunk <- flag_codes(chunk, dx9_cols,
                                                            setNames(list(dx_codes_9[[flag]]), colname),
                                                            flag_type, current_year)
                                      }
                                    }
                                  }
                                  # Flag using ICD-10 codes if requested or if combining both
                                  if (icd_version %in% c("10", "both")) {
                                    if (!is.null(dx_codes_10)) {
                                      for (flag in names(dx_codes_10)) {
                                        colname <- if (icd_version == "10") flag else paste0(flag, "_10")
                                        chunk <- flag_codes(chunk, dx10_cols,
                                                            setNames(list(dx_codes_10[[flag]]), colname),
                                                            flag_type, current_year)
                                      }
                                    }
                                  }

                                  # When combining both, merge the temporary flags into one final column.
                                  if (icd_version == "both") {
                                    for (flag in union(names(dx_codes_9), names(dx_codes_10))) {
                                      col9 <- paste0(flag, "_9")
                                      col10 <- paste0(flag, "_10")
                                      if (col9 %in% names(chunk) && col10 %in% names(chunk)) {
                                        chunk[[flag]] <- as.integer(chunk[[col9]] == 1 | chunk[[col10]] == 1)
                                        chunk[[col9]] <- NULL
                                        chunk[[col10]] <- NULL
                                      } else if (col9 %in% names(chunk)) {
                                        chunk[[flag]] <- chunk[[col9]]
                                        chunk[[col9]] <- NULL
                                      } else if (col10 %in% names(chunk)) {
                                        chunk[[flag]] <- chunk[[col10]]
                                        chunk[[col10]] <- NULL
                                      }
                                    }
                                  }

                                  chunk
                                }

      stopCluster(cl)

      processed_data <- bind_rows(flagged_chunks)

      old_width <- getOption("width")
      options(width = 200)
      message("Preview for year ", current_year, ":")
      print(head(processed_data, n = 3))
      options(width = old_width)

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

  # ---------------------------------------------------------------------------
  # Main Loop: Process Each Year
  # ---------------------------------------------------------------------------
  for (yr in years) {
    process_year(yr, fst_dir, output_dir, import_cols, rename_mapping, discwt_data, adult,
                 dx_codes_9, dx_codes_10, pr_codes_9, pr_codes_10,
                 flag_type, old_weights, weights_dir, num_cores, dry_run, icd_version)
  }

  message("All processing complete. Parallel clusters stopped (if any).")
  invisible(NULL)
}
