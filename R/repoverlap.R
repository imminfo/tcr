#' General function for the repertoire overlap evaluation.
#' 
#' @description
#' 
#' @param ... Input clonesets or a list(s) of clonesets.
#' @param .method Which method to use for the overlap evaluation. See "Details" for methods.
#' @param .seq Which column to use for overlap: 
#' 
#' @seealso  \link{intersect}, \link{similarity}
#' 
repOverlap <- function (...,
                        .method = c("exact", "hamm", "lev", "jaccard", "morisita", "cosine", "tversky", "overlap", "horn"), 
                        .seq = c("nuc", "aa"), 
                        .quant = c("read.count", "bc.count", "read.prop", "bc.prop"),
                        .vgene = F,
                        .norm = T,
                        .a = .5,
                        .b = .5,
                        .do.unique = T,
                        .do.norm = T,
                        .laplace = 0,
                        .verbose = T
                        ) {
  
  .fun <- NULL
  
  switch(.method[1], 
         exact = { .fun <- function (x) intersect(x, .norm = .norm, .verbose = .verbose) })
  
}