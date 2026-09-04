required_columns <- c("CHR", "POS", "SNP", "PIP")
max_upload_bytes <- 100 * 1024^2

validate_finemap_data <- function(data) {
  if (!is.data.frame(data)) {
    stop("Fine-mapping data must be a data frame.", call. = FALSE)
  }
  missing_columns <- setdiff(required_columns, names(data))
  if (length(missing_columns) > 0) {
    stop(
      paste(
        "Missing required column(s):",
        paste(missing_columns, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  if (nrow(data) == 0) {
    stop("The fine-mapping data contains no rows.", call. = FALSE)
  }
  if (anyNA(data$CHR) || any(!as.character(data$CHR) %in% c("22", "chr22"))) {
    stop("CHR must contain chromosome 22 only (22 or chr22).", call. = FALSE)
  }
  pos_numeric <- as.numeric(data$POS)
  if (anyNA(data$POS) || any(!is.finite(pos_numeric)) || any(pos_numeric < 1)) {
    stop(
      "POS must contain positive numeric base-pair positions.",
      call. = FALSE
    )
  }
  if (anyNA(data$SNP) || any(!nzchar(trimws(as.character(data$SNP))))) {
    stop(
      "SNP must contain a non-empty variant identifier for every row.",
      call. = FALSE
    )
  }
  pip_numeric <- as.numeric(data$PIP)
  if (anyNA(data$PIP) || any(!is.finite(pip_numeric)) ||
        any(pip_numeric < 0 | pip_numeric > 1)) {
    stop("PIP must contain numeric values between 0 and 1.", call. = FALSE)
  }
  invisible(TRUE)
}

normalize_finemap_data <- function(data) {
  validate_finemap_data(data)
  data.frame(
    CHR = 22L,
    POS = as.numeric(data$POS),
    SNP = as.character(data$SNP),
    PIP = as.numeric(data$PIP),
    stringsAsFactors = FALSE
  )
}

read_finemap_csv <- function(path) {
  if (!file.exists(path)) {
    stop("The uploaded file could not be found.", call. = FALSE)
  }
  data <- tryCatch(
    read.csv(path, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(error) {
      stop("Could not read CSV: ", error$message, call. = FALSE)
    }
  )
  normalize_finemap_data(data)
}

filter_finemap_data <- function(data, position_range) {
  if (length(position_range) != 2 || any(!is.finite(position_range))) {
    stop("Position range must contain two finite values.", call. = FALSE)
  }
  in_range <- data$POS >= position_range[1] & data$POS <= position_range[2]
  data[in_range, , drop = FALSE]
}

top_finemap_variants <- function(data, top_n) {
  ordered <- data[order(-data$PIP, data$POS), , drop = FALSE]
  ordered[seq_len(min(top_n, nrow(ordered))), , drop = FALSE]
}
