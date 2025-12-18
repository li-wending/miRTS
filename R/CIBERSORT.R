#' A small function to download the CIBERSORT function
#' @description
#' The CIBERSORT.R is licensed but free of charge for non-commercial use only.
#' Following registration, the **CIBERSORT.R** file is available via the [cibersort website](https://cibersortx.stanford.edu/).
#' @return the cibersortx website.
#'
#' @examples
#' CIBERSORT_download()
#' @export

CIBERSORT_download <- function(){
  print("Please download 'CIBERSORT.R' from https://cibersortx.stanford.edu/ , navigate 'Menu -> CS Archive -> CS Download -> Download CIBERSORT source code', and execute the file to obtain the 'CIBERSORT()' functionality!")
  utils::browseURL("https://cibersortx.stanford.edu/")
}
