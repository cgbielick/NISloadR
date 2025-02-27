#' Load NIS Hospital Data from .Do and ASC Files
#'
#' This function reads a Stata .Do file and its associated fixed-width ASCII (.ASC) file for the
#' NIS Hospital dataset. It extracts the variable schema and missing-value definitions from the .Do file,
#' reads in the ASCII file accordingly, renames the columns as specified in the schema,
#' and writes the loaded data as an fst file for fast reloading.
#'
#' @param do_path Character. The file path to the Stata .Do file.
#' @param asc_path Character. The file path to the fixed-width ASCII (.ASC) file.
#' @param output_dir Character. The working directory where the output fst file will be saved.
#'
#' @return A list with three elements:
#'   \item{schema}{A data.table containing the variable schema with columns \code{type}, \code{varname}, \code{start}, \code{end}, and \code{width}.}
#'   \item{missing_values}{A character vector of missing value definitions.}
#'   \item{data}{A data.table with the loaded NIS Hospital data (columns renamed per the schema).}
#'
#' @export
#'
#' @examples
#' \dontrun{
#'   nis_hosp <- load_nis_hospital(
#'                  do_path = "path/to/Stataload_NIS_2021_Hospital.Do",
#'                  asc_path = "path/to/NIS_2021_Hospital.ASC",
#'                  output_dir = "your/output/directory")
#' }
load_nis_hospital <- function(do_path, asc_path, output_dir) {
  # Load required packages
  requireNamespace("data.table", quietly = TRUE)
  requireNamespace("readr", quietly = TRUE)
  requireNamespace("fst", quietly = TRUE)
  requireNamespace("stringr", quietly = TRUE)

  # Set working directory
  setwd(output_dir)

  # --- Load the .Do file ---
  do_lines <- readLines(do_path)

  # --- Extract the year from line 2 ---
  year_match <- stringr::str_extract(do_lines[2], "\\d{4}")
  if (is.na(year_match)) stop("Year not found on line 2 of the .Do file.")
  year <- year_match

  # --- Part 1: Extract the variable schema from the .Do file ---
  read_section_start <- grep("^\\*\\*\\* Read data elements from the ASCII file \\*\\*\\*", do_lines)
  if (length(read_section_start) == 0) stop("Read data elements section not found.")
  using_line <- grep("^\\s*using", do_lines)
  if (length(using_line) == 0) stop("End of infix section not found.")
  infix_lines <- do_lines[(read_section_start + 1):(using_line - 1)]
  infix_lines <- trimws(infix_lines)
  infix_lines <- sub("///.*$", "", infix_lines)
  infix_lines[1] <- sub("^infix\\s+", "", infix_lines[1])

  pattern <- "^(int|byte|double|long|str)\\s+(\\S+)\\s+(\\d+)-\\s*(\\d+)"
  matches_list <- stringr::str_match_all(infix_lines, pattern)
  matches <- do.call(rbind, matches_list)
  if (nrow(matches) == 0) stop("No variable definitions matched the expected pattern.")

  schema_dt <- data.table::data.table(
    type    = matches[, 2],
    varname = matches[, 3],
    start   = as.integer(matches[, 4]),
    end     = as.integer(matches[, 5])
  )
  schema_dt[, width := end - start + 1]
  schema_dt[, type := data.table::fifelse(type == "int", "i",
                                          data.table::fifelse(type %in% c("byte", "double", "long"), "d", "c"))]

  schema_name <- paste0("nis_hosp_", year, "_schema")
  assign(schema_name, schema_dt, envir = parent.frame())

  # --- Part 2: Extract missing values ---
  missing_section_start <- grep("^\\*\\*\\* Convert special values to missing values \\*\\*\\*", do_lines)
  if (length(missing_section_start) == 0) stop("Missing values section not found.")
  missing_lines <- do_lines[(missing_section_start + 1):length(do_lines)]
  recode_lines <- missing_lines[grepl("^recode", missing_lines)]

  extracted_vals <- unlist(lapply(recode_lines, function(line) {
    m <- stringr::str_match(line, "\\(([^=]+)=")
    if (!is.na(m[1, 2])) {
      vals <- unlist(strsplit(m[1, 2], "\\s+"))
      vals <- vals[vals != ""]
      return(vals)
    } else {
      return(character(0))
    }
  }))

  missing_vals <- as.character(unique(extracted_vals))
  missing_name <- paste0("nis_hosp_", year, "_missing")
  assign(missing_name, missing_vals, envir = parent.frame())

  # --- Part 3: Read in the ASC file using the extracted schema and missing values ---
  nis_schema <- get(schema_name, envir = parent.frame())
  missing_vector <- get(missing_name, envir = parent.frame())

  asc_data <- readr::read_fwf(
    file = asc_path,
    col_positions = readr::fwf_widths(nis_schema$width),
    col_types = paste0(nis_schema$type, collapse = ""),
    trim_ws = TRUE,
    na = missing_vector
  )

  data.table::setDT(asc_data)
  data.table::setnames(asc_data, names(asc_data), nis_schema$varname)

  data_name <- paste0("nis_hosp_", year)
  assign(data_name, asc_data, envir = parent.frame())

  fst_filename <- paste0(data_name, ".fst")
  fst::write_fst(asc_data, fst_filename, compress = 100)

  return(list(
    schema = schema_dt,
    missing_values = missing_vector,
    data = asc_data
  ))
}
