# bioIOT RAW tar inspection helpers
#
# 检查 GEO *_RAW.tar 内是否含 CEL 数据文件 / 平台文件（bgx.gz），
# 不实际解压（与 13_bulk_cohort_gex.R 同源）。

#' Inspect a GEO RAW tar without extracting
#'
#' @param tar_path Path to a \code{*_RAW.tar} archive.
#'
#' @return \code{has_cel_file}: TRUE if the tar contains \code{*.CEL}
#'   (optionally \code{.gz}) members. \code{find_platform_file}: the first
#'   \code{*.bgx.gz} platform filename, or \code{NA_character_}.
#'
#' @examples
#' \dontrun{
#' has_cel_file("GSE62254_RAW.tar")
#' find_platform_file("GSE26901_RAW.tar")
#' }
#' @name tar_utils
NULL

#' @rdname tar_utils
#' @export
has_cel_file <- function(tar_path) {
  members <- untar(tar_path, list = TRUE)
  any(grepl("\\.cel(\\.gz)?$", members, ignore.case = TRUE))
}

#' @rdname tar_utils
#' @export
find_platform_file <- function(tar_path) {
  members <- untar(tar_path, list = TRUE)
  bgx <- members[grepl("\\.bgx\\.gz$", members, ignore.case = TRUE)]
  if (length(bgx) == 0) return(NA_character_)
  bgx[1]
}
