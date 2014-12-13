#' Distribution evaluation.
#' 
#' @aliases inverse.simpson diversity gini
#' 
#' @description
#' Function for evaluating the diversity of species or objects in the given distribution.
#' 
#' - True diversity, or the effective number of types, refers to the number
#' of equally-abundant types needed for the average proportional abundance
#' of the types to equal that observed in the dataset of interest 
#' where all types may not be equally abundant.
#' 
#' - Inverse Simpson index is the effective number of types that is obtained when
#' the weighted arithmetic mean is used to quantify average 
#' proportional abundance of types in the dataset of interest.
#' 
#' - The Gini coefficient measures the inequality among values
#' of a frequency distribution (for example levels of income). A Gini coefficient of zero
#' expresses perfect equality, where all values are the same (for example, where everyone
#' has the same income). A Gini coefficient of one (or 100 percents ) expresses maximal inequality
#' among values (for example where only one person has all the income).
#' 
#' Functions will check if .data if a distribution of random variable (sum == 1) or not.
#' To force normalisation and / or to prevent this, set .do.norm to TRUE (do normalisation)
#' or FALSE (don't do normalisation).
#' 
#' @usage
#' inverse.simpson(.data, .do.norm = NA, .laplace = 0)
#' 
#' diversity(.data, .q = 5, .do.norm = NA, .laplace = 0)
#' 
#' gini(.data, .do.norm = NA, .laplace = 0)
#' 
#' @param .data Vector of values.
#' @param .q q-parameter for the Diversity index.
#' @param .do.norm One of the three values - NA, T or F. If NA than check for distrubution (sum(.data) == 1)
#' and normalise if needed with the given laplace correction value. If T than do normalisation and laplace
#' correction. If F than don't do normalisaton and laplace correction.
#' @param .laplace Value for Laplace correction which will be added to every value in the .data.
#' 
#' @return Numeric vector of length 1 with value.
#' 
#' @seealso \link{entropy}, \link{similarity}
inverse.simpson <- function (.data, .do.norm = NA, .laplace = 0) {
  .data <- check.distribution(.data, .do.norm, .laplace)
  1 / sum(.data ^ 2)
}

diversity <- function (.data, .q = 5, .do.norm = NA, .laplace = 0) {
  .data <- check.distribution(.data, .do.norm, .laplace)
  1 / (sum(.data ^ .q) ^ (1 / (.q - 1)))
}

gini <- function (.data, .do.norm = NA, .laplace = 0) {
  .data <- sort(check.distribution(.data, .do.norm, .laplace))
  n <- length(.data)
  1 / n * (n + 1 - 2 * sum((n + 1 - 1:n) * .data) / sum(.data))
}


#' Diversity evaluation using rarefaction.
#' 
#' @description
#' Sequentially resample the given data with growing sample size the given data and compute mean number of unique clones.
#' For more details on the procedure see "Details".
#' 
#' @param .data Data frame or a list with data frames.
#' @param .n How many samples choose for each step.
#' @param .step Step's size.
#' @param .quantile Numeric vector of length 2 with quantiles for confidence intervals.
#' @param .col Column's name from which choose frequency of each clone.
#' @param .verbose If T than print progress bar.
#' 
#' @return
#' Data frame with first column for sizes, second columns for the first quantile,
#' third column for the mean, fourth columns for the second quantile, fifth columns
#' for the name of subject.
#' 
#' @details
#' This subroutine is designed for diversity evaluation of repertoires. On each step it computes a
#' mean unique clones from sample of fixed size using bootstrapping. Unique clones for each sample from bootstrap computed
#' as a number of non-zero elements in a vector from multinomial distribution with input vector of probabilities from the \code{.col} column
#' using function \code{rmultinom} with parameters n = .n, size = i * .step, prob = .data[, .col] (i is an index of current iteration)
#' and choosing for lower and upper bound \code{quantile} bounds of the computed distribution of unique clones.
#' 
#' @seealso \link{vis.rarefaction} \link{rmultinom}
#' 
#' @examples
#' \dontrun{
#' rarefaction(immdata, .col = "Read.count")
#' }
rarefaction <- function (.data, .n = 10, .step = 30000, .quantile = c(.025, .975), .col = 'Barcode.count', .verbose = T) {
  if (has.class(.data, 'data.frame')) {
    .data <- list(Data = .data)
  }
  
  if (.verbose) {
    pb <- set.pb(sum(sapply(1:length(.data), function (i) {
      bc.vec <- .data[[i]][, .col]
      bc.sum <- sum(.data[[i]][, .col])
      sizes <- seq(.step, bc.sum, .step)
      if (sizes[length(sizes)] != bc.sum) {
        sizes <- c(sizes, bc.sum)
      }
      length(sizes)
    } )))
  }
  
  muc.list <- lapply(1:length(.data), function (i) {
    bc.vec <- .data[[i]][, .col]
    bc.sum <- sum(.data[[i]][, .col])
    sizes <- seq(.step, bc.sum, .step)
    if (sizes[length(sizes)] != bc.sum) {
      sizes <- c(sizes, bc.sum)
    }
    muc.res <- t(sapply(sizes, function (sz) {
      muc <- apply(rmultinom(.n, sz, bc.vec / bc.sum), 2, function (col) sum(col > 0))
      res <- c(sz, quantile(muc, .quantile[1]), mean(muc), quantile(muc, .quantile[2]))
      names(res) <- c('Size', paste0('Q', .quantile[1]), 'Mean', paste0('Q', .quantile[2]))
      if (.verbose) add.pb(pb)
      res
    }))
    data.frame(muc.res, People = names(.data)[i], stringsAsFactors = F)
  })
  if (.verbose) close(pb)
  
  do.call(rbind, muc.list)
}