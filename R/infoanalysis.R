#' Normalised log assymetry.
#' 
#' @description
#' Compute the value of the normalised log assymetry measure for the given data.frames
#' of the counts of shared clones.
#'
#' @param .alpha First mitcr data.frame or a list with data.frames.
#' @param .beta Second mitcr data.frame or NULL if \code{.alpha} is a list.
#' @param .by Which column use to merge. See "Details".
#' 
#' @details
#' Merge two data frames by the given column and compute
#' value Sum(Log((Percentage for shared clone S from alpha) / (Percentage for shared clone S from beta))) / (# of shared clones).
#' 
#' @return Value of the normalised log assymetry measure for the given data.frames.
assymetry<-function(.alpha, .beta = NULL, .by = 'CDR3.nucleotide.sequence'){
  if (class(.alpha) == 'list') {
    return(apply.asymm(.alpha, assymetry, .by = .by))
  }
  m<-merge(.alpha, .beta, by=.by)
  sum(log(m$Percentage.x/m$Percentage.y))/nrow(m)
}


#' Repertoires' analysis using information measures applied to V- and J- segment frequencies.
#' 
#' @aliases entropy.seg js.div.seg
#' 
#' @description
#' Information approach to repertoire analysis. Function \code{entropy.seg} applies Shannon entropy to V-usage and hence measures variability of V-usage.
#' Function \code{js.div.seg} applied Jensen-Shannon divergence to V-usage of two or more data frames and hence measures distance among this V-usages.
#' 
#' @usage
#' entropy.seg(.data, .frame = c('all', 'in', 'out'),
#'             .alphabet = if (.VJ) "beta" else 'TRBV',
#'             .meat = F, .other = T, .VJ = F, .laplace = 1)
#' 
#' js.div.seg(.data, .data2 = NULL, .frame = c('all', 'in', 'out'), .norm.entropy = T,
#'            .alphabet = if (.VJ) "beta" else 'TRBV', .meat = F, .other = T, .VJ = F,
#'            .verbose = T, .laplace = 1)
#' 
#' @param .data Mitcr data.frame or a list with mitcr data.frames.
#' @param .data2 NULL if .data is a list, or a second mitcr data.frame.
#' @param .frame Character vector of length 1 specified which *-frames should be used:
#' only in-frame ('in'), out-of-frame ('out') or all sequences ('all').
#' @param .norm.entropy If T than divide result by mean entropy of 2 segments' frequencies. 
#' @param .alphabet Parameter to \code{freq.segments()} and \code{freq.segments.2D()} functions.
#' @param .meat Parameter to \code{freq.segments()} and \code{freq.segments.2D()} functions.
#' @param .other Parameter to \code{freq.segments()} and \code{freq.segments.2D()} functions.
#' @param .VJ If F than apply \code{freq.segments} function, else apply \code{freq.segments.2D} function.
#' @param .verbose If T than print progress of function executing.
#' @param .laplace Parameter passed to \code{freq.segments}.
#' 
#' @return For \code{entropy.seg} - numeric integer with entropy value(s). For \code{js.div.seg} - integer of vector one if \code{.data} and \code{.data2} are provided;
#' esle matrix length(.data) X length(.data) if \code{.data} is a list.
#' 
#' @seealso \link{vis.heatmap}, \link{vis.group.boxplot}, \link{freq.segments}
entropy.seg <- function (.data, .frame = c('all', 'in', 'out'),
                         .alphabet = if (.VJ) "beta" else 'TRBV',
                         .meat = F, .other = T, .VJ = F, .laplace = 1) {
  if (.VJ) .fun <- freq.segments.2D
  else     .fun <- freq.segments
  
  if (class(.data) == 'list') {
    return(sapply(.data, entropy.seg, .frame = .frame, .alphabet = .alphabet, .meat = .meat, .other = .other, .VJ = .VJ))
  }
  
  .data <- get.frames(.data, .frame)
  entropy(as.matrix(.fun(.data, .alphabet = .alphabet, .meat = .meat, .other = .other)[,-1]))
}

js.div.seg <- function (.data, .data2 = NULL, .frame = c('all', 'in', 'out'), .norm.entropy = T, .alphabet = if (.VJ) "beta" else 'TRBV', .meat = F, .other = T, .VJ = F, .verbose = T, .laplace = 1) {
  if (.VJ) .fun <- freq.segments.2D
  else     .fun <- freq.segments
  
  if (class(.data) == 'list') {
    return(apply.symm(.data, js.div.seg, .frame = .frame, .norm.entropy = .norm.entropy, .alphabet = .alphabet, .meat = .meat, .other = .other, .VJ = .VJ, .verbose= .verbose, .laplace = .laplace))
  }
  
  .data <- get.frames(.data, .frame)
  
  freq.alpha <- as.matrix(.fun(.data, .alphabet = .alphabet, .meat = .meat, .other = .other, .laplace = .laplace)[,-1])
  freq.beta <- as.matrix(.fun(.data2, .alphabet = .alphabet, .meat = .meat, .other = .other, .laplace = .laplace)[,-1])
  nrm = if (.norm.entropy) 0.5 * (entropy(freq.alpha) + entropy(freq.beta)) else 1
  js.div(freq.alpha, freq.beta) / nrm
}