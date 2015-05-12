#' General function for the repertoire overlap evaluation.
#' 
#' @description
#' General interface to all cloneset overlap functions.
#' 
#' @param .data List of clonesets.
#' @param .method Which method to use for the overlap evaluation. See "Details" for methods.
#' @param .seq Which clonotype sequences to use for the overlap: "nuc" for "CDR3.nucleotide.sequence", "aa" for 
#' "CDR3.amino.acid.sequence".
#' @param .quant Which column to use for the quantity of clonotypes: "read.count" for the "Read.count" column, 
#' "bc.count" for the "Barcode.count" column, "read.prop" for the "Read.proportion" column, "bc.prop" for 
#' the "Barcode.proportion" column. Ignored for all methods excluding "morisita".
#' @param .vgene
#' @param .norm
#' @param .a
#' @param .b
#' @param .do.unique
#' @param .verbose If T than output the data processing progress bar.
#' 
#' @details
#' Overlap methods:
#' 
#' 
#' 
#' @seealso  \link{intersectClonesets}, \link{similarity}, \link{repDiversity}
#' 
repOverlap <- function (.data,
                        .method = c("exact", "hamm", "lev", "jaccard", "morisita", "tversky", "overlap", "horn"), 
                        .seq = c("nuc", "aa"), 
                        .quant = c("read.count", "bc.count", "read.prop", "bc.prop"),
                        .vgene = F,
                        .norm = T,
                        .a = .5,
                        .b = .5,
                        .do.unique = T,
                        .verbose = T
                        ) {
  
  quant <- .column.choice(.quant, .verbose)
  
  if (!has.class(.data, 'list')) { cat("Error! Input data MUST be a list of clonesets!\n"); return(NA); }
  
  if (length(.data) < 2) { cat("Error! Input data MUST be a list of clonesets with minimum length equal to 2!\n"); return(NA); }
  
  .data <- .fix.listnames(.data)
  
  .fun <- NULL
  
  switch(.method[1], 
         exact = { .fun <- function (x) intersect(x, .norm = .norm, .verbose = .verbose) },
         { .verbose.msg("You have specified an invalid method identifier. Choosed method: normalised number of shared nuc clonotypes. \n", .verbose); { 1 + 1 } })
  
}