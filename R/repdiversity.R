#' General function for the repertoire diversity estimation.
#' 
#' @description
#' General interface to all diversity functions.
#' 
#' @param .data Cloneset or a list of clonesets.
#' @param .method Which method to use for the diversity estimation. See "Details" for methods.
#' @param .quant Which column to use for the quantity of clonotypes: "read.count" for the "Read.count" column, 
#' "bc.count" for the "Barcode.count" column, "read.prop" for the "Read.proportion" column, "bc.prop" for 
#' the "Barcode.proportion" column.
#' 
#' @seealso \link{diversity}, \link{entropy}
#' 
#' @examples
#' \dontrun{
#' data(twb)
#' twb.div <- repDiversity(twb, "chao1", "read.count")
#' ggplot(aes(), data = twb.div) + 
#' vis.group.boxplot()
#' }
repDiversity <- function (.data,
                          .method = c("chao1", "gini.simp", "inv.simp", "gini", "div", "entropy"), 
                          .quant = c("read.count", "bc.count", "read.prop", "bc.prop"), 
                          .q = 5,
                          .norm = F,
                          .do.norm = NA,
                          .laplace = 0,
                          .verbose = T) {
  
  quant <- .column.choice(.quant, .verbose)
  
  fun <- switch(.method[1], 
                chao1 = function (x, ...) chao1(x),
                gini.simp = gini.simpson,
                inv.simp = inverse.simpson,
                gini = gini,
                div = function (x, ...) diversity(x, .q = .q, ...),
                entropy = function (x, ...) entropy(x, .norm = .norm, ...),
                { .verbose.msg("You have specified an invalid method identifier. Choosed method: chao1\n", .verbose); chao1 })
  
  if (has.class(.data, 'data.frame')) { .data <- list(Data = .data) }
  
  sapply(.data, function (x) fun(x[[quant]], .do.norm = .do.norm, .laplace = .laplace))
}