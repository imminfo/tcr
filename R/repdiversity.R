#' General function for the repertoire diversity estimation.
#' 
#' @description
#' 
#' @param .data Input clonesets or a list(s) of clonesets.
#' @param .method Which method to use for the diversity estimation. See "Details" for methods.
#' @param .seq Which column to use for overlap: 
#' 
#' @seealso  \link{diversity}, \link{entropy}
#' 
repDiversity <- function (...,
                          .method = c("chao1", "gini.simp", "inv.simp", "gini", "div", "entropy"), 
                          .quant = c("read.count", "bc.count", "read.prop", "bc.prop"), 
                          .q = 5,
                          .do.norm = NA,
                          .laplace = 0
                          ) {
  
}