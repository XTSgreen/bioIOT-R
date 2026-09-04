# bioIOT GSM sample id helper
#
# 从 CEL 文件名提取 GSM 样本 ID，如 'GSM1523727_INT....CEL.gz' -> 'GSM1523727'
# （与 13_bulk_cohort_gex.R 的 gsm_id 同源，向量化；非匹配返回 NA）。

#' Extract GSM ids from CEL filenames
#'
#' @param filenames Character vector of filenames such as
#'   \code{"GSM1523727_INT_A.CEL.gz"}.
#'
#' @return Character vector of \code{"GSM"} ids; \code{NA} where no match.
#'
#' @examples
#' gsm_id(c("GSM1523727_INT_A.CEL.gz", "sample_no_id.CEL", NA))
#' @export
gsm_id <- function(filenames) {
  filenames <- as.character(filenames)
  out <- rep(NA_character_, length(filenames))
  ok <- !is.na(filenames)
  if (any(ok)) {
    x <- filenames[ok]
    pos <- regexpr("^GSM[0-9]+", x)
    hit <- which(pos > 0)
    if (length(hit)) {
      mlen <- attr(pos, "match.length")[hit]
      out[which(ok)[hit]] <- substr(x[hit], pos[hit], pos[hit] + mlen - 1)
    }
  }
  out
}
