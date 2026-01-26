##################################################
## R scripts for miRNet
## Description: Apache Arrow infrastructure for zero-copy data exchange
## Part of the Rserve/qs to Apache Arrow migration
## Author: XiaLab
###################################################

#' Shadow save function - saves to both qs and Arrow format
#' @description Saves an R object to qs format and additionally creates
#' an Arrow/Feather file for zero-copy reading from Java.
#' ALWAYS preserves rownames as the first column (row_names_id).
#' @param obj The R object to save (matrix or data.frame)
#' @param file The qs file path (will create .arrow sibling)
#' @export
shadow_save <- function(obj, file) {
  # Always save to qs for R compatibility
  qs::qsave(obj, file)

  # Additionally save to Arrow for Java zero-copy access
  arrow_path <- sub("\\.qs$", ".arrow", file)

  tryCatch({
    if (is.matrix(obj) || is.data.frame(obj)) {
      df <- as.data.frame(obj)

      # CRITICAL: Preserve rownames as first column
      rn <- rownames(obj)
      if (!is.null(rn)) {
        df <- cbind(row_names_id = as.character(rn), df)
      }

      # Write uncompressed for fastest memory-mapped access
      arrow::write_feather(df, arrow_path, compression = "uncompressed")
    }
  }, error = function(e) {
    warning("Arrow save failed (qs succeeded): ", e$message)
  })
}

#' Shadow read function - reads from Arrow if available, falls back to qs
#' @description Reads an R object, preferring Arrow format for performance.
#' Automatically restores rownames from row_names_id column.
#' @param file The qs file path (will try .arrow sibling first)
#' @return The R object
#' @export
shadow_read <- function(file) {
  arrow_path <- sub("\\.qs$", ".arrow", file)

  if (file.exists(arrow_path)) {
    tryCatch({
      df <- arrow::read_feather(arrow_path)

      # Restore rownames if present
      if ("row_names_id" %in% colnames(df)) {
        rn <- df$row_names_id
        df <- df[, colnames(df) != "row_names_id", drop = FALSE]
        rownames(df) <- rn
      }

      return(as.data.frame(df))
    }, error = function(e) {
      warning("Arrow read failed, falling back to qs: ", e$message)
    })
  }

  # Fallback to qs
  return(qs::qread(file))
}

#' Save matrix to Arrow format
#' @description Saves a numeric matrix to Arrow/Feather format for
#' zero-copy access from Java via memory-mapped files.
#' @param mat The matrix to save
#' @param file_path The output file path (.arrow extension recommended)
#' @export
save_matrix_arrow <- function(mat, file_path) {
  if (!is.matrix(mat)) {
    stop("Input must be a matrix")
  }

  df <- as.data.frame(mat)

  # Preserve rownames
  rn <- rownames(mat)
  if (!is.null(rn)) {
    df <- cbind(row_names_id = as.character(rn), df)
  }

  arrow::write_feather(df, file_path, compression = "uncompressed")
}

#' Read matrix from Arrow format
#' @description Reads a matrix from Arrow/Feather format.
#' @param file_path The input file path
#' @return The matrix with restored rownames
#' @export
read_matrix_arrow <- function(file_path) {
  df <- arrow::read_feather(file_path)

  # Restore rownames if present
  rn <- NULL
  if ("row_names_id" %in% colnames(df)) {
    rn <- df$row_names_id
    df <- df[, colnames(df) != "row_names_id", drop = FALSE]
  }

  mat <- as.matrix(df)

  if (!is.null(rn)) {
    rownames(mat) <- rn
  }

  return(mat)
}
