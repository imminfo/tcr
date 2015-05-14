########## Intersections among sets of sequences ##########


#' Intersection between sets of sequences or any elements.
#' 
#' @aliases intersectClonesets intersectCount intersectLogic intersectIndices
#' 
#' @description
#' Functions for the intersection of data frames with TCR / Ig data. 
#' See the \code{repOverlap} function for a general interface to all overlap analysis functions.
#' 
#' \code{intersectClonesets} - returns number of similar elements in the given two clonesets / data frames or matrix
#' with counts of similar elements among each pair of objects in the given list.
#' 
#' \code{intersectCount} - similar to \code{tcR::intersectClonesets}, but with fewer parameters and only for two objects.
#' 
#' \code{intersectIndices} - returns matrix M with two columns, where element with index M[i, 1] in the first
#' given object is similar to an element with index M[i, 2] in the second given object.
#' 
#' \code{intersectLogic} - returns logic vector with TRUE values in positions, where element in the first given data frame
#' is found in the second given data frame.
#' 
#' @usage
#' intersectClonesets(.alpha = NULL, .beta = NULL, .type = "n0e", .head = -1, .norm = F,
#'           .verbose = F)
#' 
#' intersectCount(.alpha, .beta, .method = c('exact', 'hamm', 'lev'), .col = NULL)
#' 
#' intersectIndices(.alpha, .beta, .method = c('exact', 'hamm', 'lev'), .col = NULL)
#' 
#' intersectLogic(.alpha, .beta, .method = c('exact', 'hamm', 'lev'), .col = NULL)
#' 
#' @param .alpha Either first vector or data.frame or list with data.frames.
#' @param .beta Second vector or data.frame or type of intersection procedure (see the \code{.type} parameter) if \code{.alpha} is a list.
#' @param .type Types of intersection procedure if \code{.alpha} and \code{.beta} is data frames. String with 3 characters (see 'Details' for more information).
#' @param .head Parameter for the \code{head} function, applied before intersecting.
#' @param .method Method to use for intersecting string elements: 'exact' for exact matching, 'hamm' for matching strings which have <= 1 hamming distance,
#' 'lev' for matching strings which have <= 1 levenshtein (edit) distance between them.
#' @param .col Which columns use for fetching values to intersect. First supplied column matched with \code{.method}, others as exact values.
#' @param .norm If TRUE than normalise result by product of length or nrows of the given data.
#' @param .verbose if T then produce output of processing the data.
#' 
#' @details
#' Parameter \code{.type} of the \code{intersectClonesets} function is a string of length 3
#' [0an][0vja][ehl], where:
#' \enumerate{
#'  \item First character defines which elements intersect ("a" for elements from the column "CDR3.amino.acid.sequence", 
#'  "n" for elements from the column "CDR3.nucleotide.sequence", other characters - intersect elements as specified);
#'  \item Second character defines which columns additionaly script should use
#' ('0' for cross with no additional columns, 'v' for cross using the "V.gene" column, 
#' 'j' for cross using "J.gene" column, 'a' for cross using both "V.gene" and "J.gene" columns);
#'  \item Third character defines a method of search for similar sequences is use:
#'  "e" stands for the exact match of sequnces, "h" for match elements which have the Hamming distance between them
#'  equal to or less than 1, "l" for match elements which have the Levenshtein distance between tham equal to or less than 1.
#' }
#' 
#' @seealso  \link{repOverlap}, \link{vis.heatmap}, \link{vis.group.boxplot}
#' 
#' @return
#' \code{intersectClonesets} returns (normalised) number of similar elements or matrix with numbers of elements.
#' 
#' \code{intersectCount} returns number of similar elements.
#' 
#' \code{intersectIndices} returns 2-row matrix with the first column stands for an index of an element in the given \code{x}, and the second column stands for an index of an element of \code{y} which is similar to a relative element in \code{x}; 
#' 
#' \code{intersectLogic} returns logical vector of \code{length(x)} or \code{nrow(x)}, where TRUE at position \code{i} means that element with index {i} has been found in the \code{y}
#' 
#' @examples
#' \dontrun{
#' data(twb)
#' # Equivalent to intersectClonesets(twb[[1]]$CDR3.nucleotide.sequence,
#' #                         twb[[2]]$CDR3.nucleotide.sequence)
#' # or intersectCount(twb[[1]]$CDR3.nucleotide.sequence,
#' #                    twb[[2]]$CDR3.nucleotide.sequence)
#' # First "n" stands for a "CDR3.nucleotide.sequence" column, "e" for exact match.
#' twb.12.n0e <- intersectClonesets(twb[[1]], twb[[2]], 'n0e')
#' stopifnot(twb.12.n0e == 46)
#' # First "a" stands for "CDR3.amino.acid.sequence" column.
#' # Second "v" means that intersect should also use the "V.gene" column.
#' intersectClonesets(twb[[1]], twb[[2]], 'ave')
#' # Works also on lists, performs all possible pairwise intersections.
#' intersectClonesets(twb, 'ave')
#' # Plot results.
#' vis.heatmap(intersectClonesets(twb, 'ave'), .title = 'twb - (ave)-intersection', .labs = '')
#' # Get elements which are in both twb[[1]] and twb[[2]].
#' # Elements are tuples of CDR3 nucleotide sequence and corresponding V-segment
#' imm.1.2 <- intersectLogic(twb[[1]], twb[[2]],
#'                            .col = c('CDR3.amino.acid.sequence', 'V.gene'))  
#' head(twb[[1]][imm.1.2, c('CDR3.amino.acid.sequence', 'V.gene')])
#' }
intersectClonesets <- function (.alpha = NULL, .beta = NULL, .type = 'n0e', .head = -1, .norm = F, .verbose = F) {
  if (class(.alpha) == 'list') {
    if (class(.beta) == 'character') {
      .type <- .beta
    }
    apply.symm(.alpha, intersectClonesets, .head = .head, .type = .type, .norm = .norm, .verbose = .verbose)
  } else {
    if (.head != -1) {
      .alpha <- head(.alpha, .head)
      .beta <- head(.beta, .head)
    }
    
    cols <- NULL
    if (class(.alpha) == 'data.frame') {
      if (substr(.type, 1, 1) == 'a' || substr(.type, 1, 1) == 'n') {
        if (substr(.type, 1, 1) == 'a') {
          cols <- 'CDR3.amino.acid.sequence'
        } else {
          cols <- 'CDR3.nucleotide.sequence'
        }
        
        if (substr(.type, 2, 2) == 'v') {
          cols <- c(cols, 'V.gene')
        } else if (substr(.type, 2, 2) == 'j') {
          cols <- c(cols, 'J.gene')
        } else if (substr(.type, 2, 2) == 'a') {
          cols <- c(cols, 'V.gene', 'J.gene')
        } else if (substr(.type, 2, 2) != '0') {
          cat("Second character in .type:", .type, 'is unknown!\n')
        }
      }
    }
    
    method <- switch(substr(.type, 3, 3),
                  e = 'exact',
                  h = 'hamm',
                  l = 'lev')
    res <- intersectCount(.alpha, .beta, method, cols)
    if (.norm) {
      if (is.null(dim(.alpha))) {
        res <- res / (as.numeric(length(.alpha)) * length(.beta))
      } else {
        res <- res / (as.numeric(nrow(.alpha)) * nrow(.beta))
      }
    }
    res
  }
}

intersectCount <- function (.alpha, .beta, .method = c('exact', 'hamm', 'lev'), .col = NULL) {
  if (.method[1] == 'exact' && (is.null(.col) || length(.col) == 1)) {
    if (is.null(.col)) {
      length(base::intersect(.alpha, .beta))
    } else {
      length(base::intersect(.alpha[,.col], .beta[,.col]))
    }
  } else {
    if (.method[1] == 'exact') {
      nrow(dplyr::intersect(.alpha[, .col], .beta[, .col]))
    } else {
      res <- intersectIndices(unique(.alpha), unique(.beta), .method, .col)
      if (!is.na(res[1,1])) {
        nrow(res)
      } else {
        0
      }
    }
  }
}

intersectIndices <- function (.alpha, .beta, .method = c('exact', 'hamm', 'lev'), .col = NULL) {    
  if (!is.null(.col)) {
    .alpha <- as.data.frame(.alpha, stringsAsFactors = F)
    .beta <- as.data.frame(.beta, stringsAsFactors = F)
    if (length(.col) == 1) {
      .alpha <- .alpha[, .col]
      .beta <- .beta[, .col]
    } else {
      if (length(.col) == 2) {
        alpha.cols <- as.matrix(.alpha[, .col[2]])
        beta.cols <- as.matrix(.beta[, .col[2]])
      } else {
        alpha.cols <- .alpha[, .col[2:length(.col)]]
        beta.cols <- .beta[, .col[2:length(.col)]]
      }
      .alpha <- .alpha[, .col[1]]
      .beta <- .beta[, .col[1]]
    }
  }
  
  res <- find.similar.sequences(.alpha, .beta, .method, 1, F)
  
  if (is.na(res[1,1])) {
    return(res)
  }
  
  if (!is.null(.col) && length(.col) > 1) {
    for (i in 1:ncol(alpha.cols)) {
      res <- res[alpha.cols[res[,1], i] == beta.cols[res[,2], i], ]
    }
  }
  
  if(!is.matrix(res)) {
    res <- rbind(res)
  }
  
  if (dim(res)[1] == 0) {
    res <- matrix(c(NA, NA), 1, 2)
  }
  
  res
}

intersectLogic <- function (.alpha, .beta, .method = c('exact', 'hamm', 'lev'), .col = NULL) {
  ind <- intersectIndices(.alpha, .beta, .method, .col)
  if (is.null(.col)) {
    logic <- rep(FALSE, length(.alpha))
  } else {
    logic <- rep(FALSE, nrow(.alpha))
  }
  logic[ind[,1]] <- TRUE
  logic
}


#' Compute convergence characteristics of repertoires.
#' 
#' @description 
#' Get a number of rows with similar aminoacid sequence but different nucleotide sequence.
#' 
#' @param .alpha Either data frame with columns \code{.col.nuc} and \code{.col.aa} or list with such data frames.
#' @param .beta Either data frame or none.
#' @param .col.nuc Name of the column with nucleotide sequences.
#' @param .col.aa Name of the columnw ith aminoacid sequences.
#' @return If \code{.alpha} is data frame, than integer vector of length 2 with . If \code{.alpha} is a list
#' than matrix M with M[i,j] = convergence.index(.alpha[[i]], .alpha[[j]]).
convergence.index <- function (.alpha, .beta, .col.nuc = 'CDR3.nucleotide.sequence', .col.aa = 'CDR3.amino.acid.sequence') {
  if (has.class(.alpha, 'list')) {
    return(apply.asymm(.alpha, function (x,y) convergence.index(x, y)[1]))
  }
  
  a.nuc.logic <- .alpha[, .col.nuc] %in% .beta[, .col.nuc]
  b.nuc.logic <- .beta[, .col.nuc] %in% .alpha[, .col.nuc]
  a.aa.logic <- .alpha[, .col.aa] %in% .beta[, .col.aa]
  b.aa.logic <- .beta[, .col.aa] %in% .alpha[, .col.aa]
  c(Amino.acid.not.nucleotide.count.1 = length(unique(.alpha[a.aa.logic & !(a.nuc.logic), .col.aa])),
    Amino.acid.not.nucleotide.count.2 = length(unique(.beta[b.aa.logic & !(b.nuc.logic), .col.aa])))
}