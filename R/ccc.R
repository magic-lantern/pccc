#' Complex Chronic Conditions (CCC)
#'
#' Generate CCC and CCC subcategory flags and the number of categories.
#'
#' It is recommended that you view the codes defining the CCC via
#' \code{\link{get_codes}} and make sure that the ICD codes in your data set are
#' formatted in the same way.  The ICD codes used for CCC are character strings
#' must be formatted as follows:
#' \itemize{
#' \item *Do not* use decimal points or other separators
#' \item ICD 9 codes: Codes less than 10 should be left padded with 2 zeros. Codes
#' less than 100 should be left padded with 1 zero.
#' }
#'
#' See See `vignette("pccc-overview")` for more details.
#'
#' @references
#' See \code{\link{pccc-package}} for published paper on the topic of identifying
#' Complex Chronix Conditions
#'
#' @param data a \code{data.frame} containing a patient id and all the ICD-9-CM
#' or ICD-10-CM codes.  The \code{data.frame} passed to the function should be
#' in wide format.
#' @param id bare name of the column containing the patient id
#' @param dx_cols,pc_cols column names with the diagnostic codes and procedure
#' codes respectively.  These argument are passed to \code{\link[dplyr]{select}}.
#' @param icdv ICD version 9 or 10
#'
#' @seealso \code{\link{get_codes}} to view the ICD codes used to define the
#' CCC.  \code{\link[dplyr]{select}} for more examples and details on how to
#' identify and select the diagnostic and procedure code columns.
#'
#' @return A \code{data.frame} with a column for the subject id and integer (0
#' or 1) columns for each each of the categories.
#'
#' @example examples/ccc.R
#'
#' @export
ccc <- function(data, id, dx_cols = NULL, pc_cols = NULL, icdv) {
  UseMethod("ccc")
}

#' @method ccc data.frame
#' @export
ccc.data.frame <- function(data, id, dx_cols, pc_cols, icdv) {

  if (missing(dx_cols) & missing(pc_cols)) {
    stop("dx_cols and pc_cols are both missing.  At least one must not be.",
         call. = FALSE)
  }

  if (!missing(dx_cols)) {
    dxmat <- sapply(dplyr::select(data, !!dplyr::enquo(dx_cols)), as.character)
    if(! is.matrix(dxmat)) {
      dxmat <- as.matrix(dxmat)
    }
  } else {
    dxmat <- matrix("", nrow = nrow(data))
  }

  if (!missing(pc_cols)) {
    pcmat <- sapply(dplyr::select(data, !!dplyr::enquo(pc_cols)), as.character)
    if(! is.matrix(pcmat)) {
      pcmat <- as.matrix(pcmat)
    }
  } else {
    pcmat <- matrix("", nrow = nrow(data))
  }

  ids <- dplyr::select(data, !!dplyr::enquo(id))

  dplyr::bind_cols(ids, ccc_mat_r(dxmat, pcmat, icdv))
  #dplyr::bind_cols(ids, ccc_mat_rcpp(dxmat, pcmat, icdv))
}

ccc_mat_r <- function(dx, pc, version = 9L) {
  codes <- get_code_set(version)

  df <- data.frame(neuromusc = integer(),
             cvd = integer(),
             respiratory = integer(),
             renal = integer(),
             gi = integer(),
             hemato_immu = integer(),
             metabolic = integer(),
             congeni_genetic = integer(),
             malignancy = integer(),
             neonatal = integer(),
             tech_dep = integer(),
             transplant = integer(),
             ccc_flag = integer())
  df[1:nrow(dx),] <- c(0L)

  for(i in 1:nrow(dx)) {
    dx_row <- dx[i,]
    pc_row <- pc[i,]

    if(find_match(dx_row, pc_row, codes$dx_neuromusc, codes$pc_neuromusc))
      df$neuromusc[i] <- 1L
    else if(find_match(dx_row, pc_row, codes$dx_cvd, codes$pc_cvd))
      df$cvd[i] <-  1L
    else if(find_match(dx_row, pc_row, codes$dx_respiratory, codes$pc_respiratory))
      df$respiratory[i] <- 1L
    else if(find_match(dx_row, pc_row, codes$dx_renal, codes$pc_renal))
      df$renal[i] <- 1L
    else if(find_match(dx_row, pc_row, codes$dx_gi, codes$pc_gi))
      df$gi[i] <- 1L
    else if (find_match(dx_row, pc_row, codes$dx_hemato_immu, codes$pc_hemato_immu))
      df$hemato_immu[i] <- 1L
    else if (find_match(dx_row, pc_row, codes$dx_metabolic, codes$pc_metabolic))
      df$metabolic[i] <- 1L
    else if (find_match(dx_row, NA, codes$dx_congeni_genetic))
      df$congeni_genetic[i] <- 1L
    else if (find_match(dx_row, pc_row, codes$dx_malignancy, codes$pc_malignancy))
      df$malignancy[i] <- 1L
    else if (find_match(dx_row, NA, codes$dx_neonatal))
      df$neonatal[i] <- 1L

    df$tech_dep[i] <- find_match(dx_row, pc_row, codes$dx_tech_dep, codes$pc_tech_dep)
    df$transplant[i] <- find_match(dx_row, pc_row, codes$dx_transplant, codes$pc_transplant)

    if (sum(df[i,]))
      df$ccc_flag[i] <- 1L
  }
  df
}

find_match <- function(dx,
                       pc,
                       dx_codes,
                       pc_codes = NULL){

  for (c in dx_codes) {
    if(any(stringi::stri_startswith_fixed(dx, c),na.rm = TRUE))
      return(1L)
  }

  for (c in pc_codes) {
    if(any(stringi::stri_startswith_fixed(pc, c),na.rm = TRUE))
      return(1L)
  }
  # return 0 if no match
  0L
}

