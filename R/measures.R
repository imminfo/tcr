########## Evaluation of distributions, vectors and sets ##########


#' Information measures.
#' 
#' @aliases entropy js.div kl.div
#' 
#' @description
#' Functions for information measures of and between distributions of values:
#' 
#' - The Shannon entropy quantifies the uncertainty (entropy or degree of surprise)
#' associated with this prediction.
#' 
#' - Kullback-Leibler divergence (information gain, information divergence, 
#' relative entropy, KLIC) is a non-symmetric measure of the difference between
#' two probability distributions P and Q (measure of information lost when Q is used to
#' approximate P).
#' 
#' - Jensen-Shannon divergence is a symmetric version of KLIC. Square root of this
#' is a metric often referred to as Jensen-Shannon distance.
#' 
#' Functions will check if \code{.data} if a distribution of random variable (sum == 1) or not.
#' To force normalisation and / or to prevent this, set \code{.do.norm} to TRUE (do normalisation)
#' or FALSE (don't do normalisation). For \code{js.div} and \code{kl.div} vectors of values must have
#' equal length.
#' 
#' @usage
#' entropy(.data, .norm = F, .do.norm = NA, .laplace = 1e-12)
#' 
#' kl.div(.alpha, .beta, .do.norm = NA, .laplace = 0)
#' 
#' js.div(.alpha, .beta, .do.norm = NA, .laplace = 0, .norm.entropy = F)
#' 
#' @param .data,.alpha,.beta Vector of values.
#' @param .norm If T than compute normalised entropy (H / Hmax).
#' @param .do.norm One of the three values - NA, T or F. If NA than check for distrubution \code{(sum(.data) == 1)}.
#' and normalise if needed with the given laplace correction value. If T than do normalisation and laplace
#' correction. If F than don't do normalisaton and laplace correction.
#' @param .laplace Value for Laplace correction which will be added to every value in the .data.
#' @param .norm.entropy If T than normalise JS-divergence by entropy.
#' 
#' @return Shannon entropy, Jensen-Shannon divergence or Kullback-Leibler divergence values.
#' 
#' @seealso \link{similarity}, \link{diversity}
entropy <- function (.data, .norm = F, .do.norm = NA, .laplace = 1e-12) {
  .data <- check.distribution(.data, .do.norm, .laplace)
  res <- - sum(.data * log2(.data))
  if (.norm) {
    res / log2(length(.data))
  } else {
    res
  }
}

kl.div <- function (.alpha, .beta, .do.norm = NA, .laplace = 0) {
  .alpha <- check.distribution(.alpha, .do.norm, .laplace)
  .beta <- check.distribution(.beta, .do.norm, .laplace)
  sum(log2(.alpha / .beta) * .alpha)
}

js.div <- function (.alpha, .beta, .do.norm = NA, .laplace = 0, .norm.entropy = F) {
  .alpha <- check.distribution(.alpha, .do.norm, .laplace)
  .beta <- check.distribution(.beta, .do.norm, .laplace)
  nrm = if (.norm.entropy) 0.5 * (entropy(.alpha, F) + entropy(.beta, F)) else 1
  M <- (.alpha + .beta) / 2
  0.5 * (kl.div(.alpha, M, F) + kl.div(.beta, M, F)) / nrm
}


#' Log-likelihood.
#' 
#' @description
#' Compute the log-likelihood of the given distribution or vector of counts.
#' 
#' @param .data Vector for distribution or counts.
#' @param .base Logarightm's base for the loglikelihood.
#' @param .do.norm Parameter to the \code{check.distribution} function.
#' @param .laplace Laplace correction, Parameter to the \code{check.distribution} function.
#' 
#' @return Loglikelihood value.
loglikelihood <- function (.data, .base = 2, .do.norm = NA, .laplace = 0.000000000001) {
  .data <- check.distribution(.data, .do.norm, .laplace)
  sum(log(.data, .base))
}


#' Set and vector similarity measures.
#' 
#' @aliases similarity cosine.similarity tversky.index overlap.coef morisitas.index jaccard.index
#' 
#' @description
#' Functions for computing similarity between two vectors or sets. See "Details" for exact formulas.
#' 
#' - Cosine similarity is a measure of similarity between two vectors of an inner product space that measures the cosine of the angle between them.
#' 
#' - Tversky index is an asymmetric similarity measure on sets that compares a variant to a prototype.
#' 
#' - Overlap cofficient is a similarity measure related to the Jaccard index that measures the overlap between two sets, and is defined as the size of the intersection divided by the smaller of the size of the two sets.
#' 
#' - Jaccard index is a statistic used for comparing the similarity and diversity of sample sets.
#' 
#' - Morisita's overlap index is a statistical measure of dispersion of individuals in a population. It is used to compare overlap among samples (Morisita 1959). This formula is based on the assumption that increasing the size of the samples will increase the diversity because it will include different habitats (i.e. different faunas).
#' 
#' @usage
#' cosine.similarity(.alpha, .beta, .do.norm = NA, .laplace = 0)
#' 
#' tversky.index(x, y, .a = 0.5, .b = 0.5)
#' 
#' overlap.coef(.alpha, .beta)
#' 
#' jaccard.index(.alpha, .beta, .intersection.number = NA)
#' 
#' morisitas.index(.alpha, .beta, .do.unique = T)
#' 
#' @param .alpha,.beta,x,y Vector of numeric values for cosine similarity, vector of any values
#' (like characters) for tversky.index and overlap.coef, matrix or data.frame with 2 columns for morisitas.index,
#' either sets or number of elements in sets for Jaccard index (see "Details" section).
#' @param .a,.b Alpha and beta parameters for Tversky Index. Default values gives the Jaccard index measure.
#' @param .do.norm One of the three values - NA, T or F. If NA than check for distrubution (sum(.data) == 1)
#' and normalise if needed with the given laplace correction value. If T than do normalisation and laplace
#' correction. If F than don't do normalisaton and laplace correction.
#' @param .laplace Value for Laplace correction.
#' @param .do.unique If T than call unique on the first columns of the given data.frame or matrix.
#' @param .intersection.number Number of intersected elements between two sets. See "Details" for more information.
#' @details
#' For \code{morisitas.index} input data are matrices or data.frames with two columns: first column is
#' elements (species or individuals), second is a number of elements (species or individuals) in a population.
#' 
#' For \code{jaccard.index} there are two ways for computing the index. ???
#' 
#' Formulas:
#' 
#' Cosine similarity: \code{cos(a, b) = a * b / (||a|| * ||b||)}
#' 
#' Tversky index: \code{S(X, Y) = |X and Y| / (|X and Y| + a*|X - Y| + b*|Y - X|)}
#' 
#' Overlap coefficient: \code{overlap(X, Y) = |X and Y| / min(|X|, |Y|)}
#' 
#' Jaccard index: \code{J(A, B) = |A and B| / |A U B|}
#' 
#' Formual for Morisita's overlap index is quite complicated and can't be easily shown here, so just look at webpage: http://en.wikipedia.org/wiki/Morisita%27s_overlap_index
#' 
#' 
#' @return Value of similarity between the given sets or vectors.
#' 
#' @seealso \link{intersect}, \link{entropy}, \link{diversity}
#' 
#' @examples
#' \dontrun{
#' jaccard.index(1:10, 2:20)
#' a <- length(unique(immdata[[1]][, c('CDR3.amino.acid.sequence', 'V.segments')]))
#' b <- length(unique(immdata[[2]][, c('CDR3.amino.acid.sequence', 'V.segments')]))
#' jaccard.index(a, b, intersect(immdata[[1]], immdata[[2]], 'ave'))
#' }
cosine.similarity <- function (.alpha, .beta, .do.norm = NA, .laplace = 0) {
  .alpha <- check.distribution(.alpha, .do.norm, .laplace)
  .beta <- check.distribution(.beta, .do.norm, .laplace)
  sum(.alpha * .beta) / (sum(.alpha ^ 2) * sum(.beta ^ 2))
}

tversky.index <- function (x, y, .a = 0.5, .b = 0.5) {
  XiY <- length(intersect(x, y))
  XiY / (XiY + .a * length(setdiff(x, y)) + .b * length(setdiff(y, x)))
}

overlap.coef <- function (.alpha, .beta) {
  length(intersect(.alpha, .beta)) / min(length(.alpha), length(.beta))
}

jaccard.index <- function (.alpha, .beta, .intersection.number = NA) {
  if (is.na(.intersection.number)) {
    length(intersect(.alpha, .beta)) / (length(unique(.alpha)) + length(unique(.beta)) - length(intersect(.alpha, .beta)))
  } else {
    .intersection.number / (.alpha + .beta - .intersection.number)
  }
}

morisitas.index <- function (.alpha, .beta, .do.unique = T) {
  colnames(.alpha) <- c('Species', 'Count')
  colnames(.beta) <- c('Species', 'Count')
  if (.do.unique) {
    .alpha <- .alpha[!duplicated(.alpha[,1]),]
    .beta <- .beta[!duplicated(.beta[,1]),]
  }
  .alpha[,2] <- as.numeric(.alpha[,2])
  .beta[,2] <- as.numeric(.beta[,2])
  merged <- merge(.alpha, .beta, by = 'Species', all = T)
  merged[is.na(merged)] <- 0
  sum.alpha <- sum(.alpha[,2])
  sum.beta <- sum(.beta[,2])
  2 * sum(merged[,2] * merged[,3] / sum.alpha) / sum.beta / (
    (sum((.alpha[,2] / sum.alpha)^2) + sum((.beta[,2] / sum.beta)^2)))
}