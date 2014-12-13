########## Support functions for managing the data ##########


#' Print the given message if second parameter is a TRUE.
#' 
#' @param .message Character vector standing for a message.
#' @param .verbose If T then print the given mesasge.
#' @return Nothing.
.verbose.msg <- function (.message, .verbose = T) {
  if (.verbose) cat(.message)
}


#' Shuffling data frames.
#' 
#' @aliases permutedf unpermutedf
#' 
#' @description
#' Shuffle the given data.frame and order it by the Read.count column or un-shuffle
#' a data frame and return it to the initial order.
#' 
#' @usage
#' permutedf(.data)
#' 
#' unpermutedf(.data)
#' 
#' @param .data MiTCR data.frame or list of such data frames.
#' 
#' @return Shuffled data.frame or un-shuffled data frame if \code{.data} is a data frame, else list of such data frames.
permutedf <- function (.data) {
  if (has.class(.data, 'list')) {
    return(lapply(.data, permutedf))
  }
  shuffle<-.data[sample(nrow(.data)),]
  shuffle[order(shuffle$Read.count, decreasing=T),]
}

unpermutedf <- function (.data) {
  if (has.class(.data, 'list')) {
    return(lapply(.data, unpermutedf))
  }
  .data[do.call(order, .data),]
}


#' Get a random subset from a data.frame.
#' 
#' @description
#' Sample rows of the given data frame with replacement.
#' 
#' @param .data Data.frame or a list with data.frames
#' @param .n Sample size if integer. If in bounds [0;1] than percent of rows to extract. "1" is a percent, not one row!
#' @param .replace If T than choose with replacement, else without.
#' 
#' @return Data.frame of nrow .n or a list with such data.frames.
sample.clones <- function (.data, .n, .replace = T) { 
  if (has.class(.data, 'list')) { return(lapply(.data, sample.clones, .n = .n)) }
  .data[sample(1:nrow(.data), if (.n > 1) .n else round(nrow(.data) * .n), replace = .replace), ]
}


#' Check if a given object has a given class.
#' 
#' @param .data Object.
#' @param .class String naming a class.
#' 
#' @return Logical.
has.class <- function (.data, .class) { .class %in% class(.data) }


#' Copy the up-triangle matrix values to low-triangle.
#' 
#' @param mat Given up-triangle matrix.
#' 
#' @return Full matrix.
matrixdiagcopy<-function(mat){
  for (i in 1:ncol(mat))
    for (j in 1:nrow(mat))
      if (i<j)
      {
        mat[j,i]<-mat[i,j]
      }
  mat
}


#' Simple functions for manipulating progress bars.
#' 
#' @aliases add.pb set.pb
#' 
#' @description
#' Set the progress bar with the given length (from zero) or add value to the created progress bar.
#' 
#' @usage
#' set.pb(.max)
#' 
#' add.pb(.pb, .value = 1)
#' 
#' @param .max Length of the progress bar.
#' @param .pb Progress bar object.
#' @param .value Value to add to the progress bar.
#' 
#' @return Progress bar (for set.pb) or length-one numeric vector giving the previous value (for add.pb).
set.pb <- function (.max) {
  txtProgressBar(min=0, max=.max, style=3)
}

add.pb <- function (.pb, .value = 1) {
  setTxtProgressBar(.pb, .pb$getVal() + .value)
}


#' Get a sample from matrix with probabilities.
#' 
#' @description
#' Get a sample from matrix or data frame with pair-wise probabilities.
#' 
#' @param .table Numeric matrix or data frame with probabilities and columns and rows names.
#' @param .count Number of sample to fetch.
#' 
#' @return Character matrix with nrow == .count and 2 columns. row[1] in row.names(.table), row[2] in colnames(.table).
sample2D <- function (.table, .count = 1) {
  if (has.class(.table, 'data.frame')) {
    .table <- as.matrix(.table)
  }
  melted.table <- melt(.table)
  as.matrix(melted.table[sample(1:nrow(melted.table), .count, T, melted.table[,3]),c(1,2)])
}


#' Apply function to every pair of data frames from a list.
#' 
#' @aliases apply.symm apply.asymm
#' 
#' @description
#' Apply the given function to every pair in the given datalist. Function either
#' symmetrical (i.e. fun(x,y) == fun(y,x)) or assymmetrical (i.e. fun(x,y) != fun(y,x)).
#' 
#' @usage 
#' apply.symm(.datalist, .fun, ..., .diag = NA, .verbose = T)
#' 
#' apply.asymm(.datalist, .fun, ..., .diag = NA, .verbose = T)
#' 
#' @param .datalist List with some data.frames.
#' @param .fun Function to apply, which return basic class value.
#' @param ... Arguments passsed to .fun.
#' @param .diag Either NA for NA or something else != NULL for .fun(x,x).
#' @param .verbose If T than output a progress bar.
#' 
#' @return Matrix with values M[i,j] = fun(datalist[i], datalist[j])
#' 
#' @examples
#' \dontrun{
#' apply.symm(immdata, intersect, .type = 'a0e') # equivalent to intersect(immdata, 'a0e')
#' }
apply.symm <- function (.datalist, .fun, ..., .diag = NA, .verbose = T) {
  res <- matrix(0, length(.datalist), length(.datalist))
  if (.verbose) pb <- set.pb(length(.datalist)^2 / 2 + length(.datalist)/2)
  for (i in 1:length(.datalist)) 
    for (j in i:length(.datalist)) {
      if (i == j && is.na(.diag)) { res[i,j] <- NA }
      else { res[i,j] <- .fun(.datalist[[i]], .datalist[[j]], ...) }
      if (.verbose) add.pb(pb)
    }
  if (.verbose) close(pb)
  row.names(res) <- names(.datalist)
  colnames(res) <- names(.datalist)
  matrixdiagcopy(res)
}

apply.asymm <- function (.datalist, .fun, ..., .diag = NA, .verbose = T) {
  res <- matrix(0, length(.datalist), length(.datalist))
  if (.verbose) pb <- set.pb(length(.datalist)^2)
  for (i in 1:length(.datalist)) 
    for (j in 1:length(.datalist)) {
      if (i == j && is.na(.diag)) { res[i,j] <- NA }
      else { res[i,j] <- .fun(.datalist[[i]], .datalist[[j]], ...) }
      if (.verbose) add.pb(pb)
    }
  if (.verbose) close(pb)
  row.names(res) <- names(.datalist)
  colnames(res) <- names(.datalist)
  res
}


#' Check for adequaty of distrubution.
#' 
#' @description
#' Check if the given .data is a distribution and normalise it if necessary with optional laplace correction.
#' 
#' @param .data Numeric vector of values.
#' @param .do.norm One of the three values - NA, T or F. If NA than check for distrubution (sum(.data) == 1)
#' and normalise if needed with the given laplace correction value. If T than do normalisation and laplace
#' correction. If F than don't do normalisaton and laplace correction.
#' @param .laplace Value for laplace correction.
#' @param .na.val Replace all NAs with this value.
#' 
#' @return Numeric vector.
check.distribution <- function (.data, .do.norm = NA, .laplace = 1, .na.val = 0) {
  if (is.na(.do.norm)) {
    .data[is.na(.data)] <- .na.val
    if (sum(.data) != 1) {
      .data <- (.data + .laplace) / sum(.data + .laplace)
    }
  } else if (.do.norm) {
    .data[is.na(.data)] <- .na.val
    .data <- (.data + .laplace) / sum(.data + .laplace)
  }
  .data
}


#' Get samples from a repertoire slice-by-slice or top-by-top and apply function to them.
#' 
#' @aliases top.fun slice.fun
#' 
#' @description
#' Functions for getting samples from data frames either by consequently applying
#' head functions (\code{top.fun}) or by getting equal number of rows in the moving window (\code{slice.fun})
#' and applying specified function to this samples.
#' 
#' @usage
#' top.fun(.data, .n, .fun, ..., .simplify = T)
#' 
#' slice.fun(.data, .size, .n, .fun, ..., .simplify = T)
#' 
#' @param .data Data.frame, matrix, vector or any enumerated type or a list of this types.
#' @param .n Vector of values passed to head function for top.fun or the number of slices for slice.fun.
#' @param .size Size of the slice for sampling for slice.fun.
#' @param .fun Funtions to apply to every sample subset. First input argument is a data.frame, others are passed as \code{...}.
#' @param ... Additional parameters passed to the .fun.
#' @param .simplify If T than try to simplify result to a vector or to a matrix if .data is a list.
#' 
#' @return List of length length(.n) for top.fun or .n for slice.fun.
#' 
#' @examples
#' \dontrun{
#' # Get entropy of V-usage for the first 1000, 2000, 3000, ... clones.
#' res <- top.fun(immdata[[1]], 1000, entropy.seg)
#' # Get entropy of V-usage for the interval of clones with indices [1,1000], [1001,2000], ...
#' res <- top.fun(immdata[[1]], 1000, entropy.seg)
#' }
top.fun <- function (.data, .n, .fun, ..., .simplify = T) {
  if (has.class(.data, 'list')) { 
    res <- lapply(.data, function (d) top.fun(d, .n, .fun, ..., .simplify = .simplify))
    if (.simplify) {
      res <- as.matrix(data.frame(res))
    }
    return(res)
  }
  
  res <- lapply(.n, function (nval) .fun(head(.data, nval), ...))
  names(res) <- .n
  if (.simplify) { res <- do.call(c, res) }
  res
}

slice.fun <- function(.data, .size, .n, .fun, ..., .simplify = T) {
  if (has.class(.data, 'list')) {
    res <- lapply(.data, function(d) { slice.fun(d, .size, .n, .fun, ..., .simplify = .simplify) })
    if (.simplify) {
      res <- as.matrix(data.frame(res))
    }
    return(res)
  }
  
  res <- lapply(1:.n, function (i) {
    i <- i - 1
    down <- i * .size + 1
    up <- (i + 1) * .size
    .fun(.data[down:up, ], ...)
  })
  names(res) <- sapply(1:.n, function (i) {
    paste0((i - 1) * .size + 1, ':', i * .size)
  })
  if (.simplify) { res <- do.call(c, res) }
  res
}


#' Internal function. Add legend to a grid of plots and remove legend from all plots of a grid.
#' 
#' @description
#' Given a list of ggplot2 plots, remove legend from each of them and return 
#' grid of such plots plus legend from the first vis. Suitable for plots
#' with similar legends.
.add.legend <- function (.vis.list, .vis.nrow = 2, .legend.ncol = 1) {
  leg <- gtable_filter(ggplot_gtable(ggplot_build(.vis.list[[1]] + guides(fill=guide_legend(ncol=.legend.ncol)))), "guide-box")
  grid.arrange(do.call(arrangeGrob, c(.vis.list, nrow = .vis.nrow)), leg, widths=unit.c(unit(1, "npc") - leg$width, leg$width), nrow = 1, main ='Top crosses')
}