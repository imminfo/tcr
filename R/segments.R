########## Statistics and analysis of V- and J-segments usage ##########


if (getRversion() >= "2.15.1") {
  utils::globalVariables(c('PC1', 'PC2'))
}


#' V- and J-segments frequency.
#' 
#' @aliases freq.segments freq.segments.2D freq.Va freq.Vb freq.Ja freq.Jb
#' 
#' @description Get frequencies or counts of segments (V / J - usage).
#' 
#' @usage
#' freq.segments(.data, .alphabet = 'TRBV', .count = F, .meat = F, .other = T,
#'               .laplace = 1, .column = NULL, .sum.col = "Read.count")
#' 
#' freq.segments.2D(.data, .alphabet = 'beta', .count = F, .meat = F,
#'                  .laplace = 1, .columns = NULL, .sum.col = "Read.count", ...)
#' 
#' freq.Va(.data, .count = F, .meat = F, .other = T, .laplace = 1, .sum.col = "Read.count")
#' 
#' freq.Vb(.data, .count = F, .meat = F, .other = T, .laplace = 1, .sum.col = "Read.count")
#' 
#' freq.Ja(.data, .count = F, .meat = F, .other = T, .laplace = 1, .sum.col = "Read.count")
#' 
#' freq.Jb(.data, .count = F, .meat = F, .other = T, .laplace = 1, .sum.col = "Read.count")
#' 
#' @param .data Mitcr data.frame or list with data.frames.
#' @param .alphabet Vector of elements in the alphabet for freq.segments, one of the strings 'TRBV' (for using V_BETA_ALPHABET variable, that user should load before calling functions (same for other strings)), 'TRAV', 'TRBJ', 'TRAJ' for V- and J-segments alphabets for freq.segments
#' or one of the 'alpha' or 'beta' for freq.segments.2D or a list of length 2 with alphabets strings for freq.segments.2D.
#' @param .count Should we return count or percentage?
#' @param .meat Compute statistics using counts of elements (e.g., Read.count) or not.
#' @param .other Should elements not in the given alphabet be shown in the result matrix.
#' @param .laplace Value for the Laplace correction.
#' @param .column,.columns Column's name with elements from the given alphabet for the freq.segments or a character vector of length two for freq.segments.2D.
#' @param .sum.col Which column use to count frequencies if \code{.meat} = T. Default: 'Read.count'.
#' @param ... Don't use it, for internal purpose.
#' 
#' @return Data.frame with columns Segments and their frequencies in the column Freq.
#' 
#' @seealso \code{\link{vis.V.usage}} \code{\link{vis.J.usage}} \code{\link{pca.segments}}
#' 
#' @examples
#' \dontrun{
#' # Load your data
#' data(twb)
#' # Load human alphabets
#' data(human.alphabets)
#' # compute V-segments frequencies
#' seg <- freq.segments(twb)
#' # plot V-segments frequencies as a grid
#' vis.grid.stats(seg)
#' # plot V-segments frequencies from the data
#' vis.V.usage(twb)
#' }
freq.segments <- function (.data, .alphabet='TRBV', .count=F, .meat=F, .other=T, .laplace=1, .column = NULL, .sum.col = "Read.count") {
  if (class(.data) == 'list') {
    res <- freq.segments(.data[[1]], .alphabet=.alphabet, .count=.count, .meat = .meat, .other = .other, .laplace=.laplace, .column = .column)
    names(res)[2] <- names(.data)[1]
    for (i in 2:length(.data)) {
      res<-merge(res, freq.segments(.data[[i]],.alphabet=.alphabet,.count=.count, .meat = .meat, .other = .other, .column = .column, .sum.col= .sum.col, .laplace = .laplace), by = 'Segment', all = T)
      names(res)[length(names(res))] <- names(.data)[i]
    }
    res$Segment <- as.character(res$Segment)
    return(res)
  }
  
  if (.alphabet == "TRBV")      { .column <- 'V.segments'; alphabet <- V_BETA_ALPHABET }
  else if (.alphabet == "TRAV") { .column <- 'V.segments'; alphabet <- V_ALPHA_ALPHABET }
  else if (.alphabet == "TRBJ") { .column <- 'J.segments'; alphabet <- J_BETA_ALPHABET }
  else if (.alphabet == "TRAJ") { .column <- 'J.segments'; alphabet <- J_ALPHA_ALPHABET }
  else                          { alphabet <- .alphabet }
  seg <- .data[[.column]]
  
  read.count <- rep.int(1, nrow(.data))
  if (.meat) { read.count <- .data[, .sum.col] }
  
  counts.l <- tapply(read.count, seg, sum)
  freqs <- counts.l[alphabet]
  freqs[is.na(freqs)] <- 0
  res <- data.frame(Segment = alphabet, Freq = freqs, stringsAsFactors = F, row.names=NULL)
  if (.other) {
    res <- rbind(res, Other = data.frame(Segment = 'Other', Freq = sum(unlist(lapply(counts.l[is.na(match(names(counts.l), alphabet))], sum)), na.rm=T), stringsAsFactors = F))
  }
  res$Freq <- res$Freq + .laplace
  if (!.count) { res$Freq <- res$Freq / sum(res$Freq) }
  res[order(res$Segment),]
}

freq.Va <- function (.data, .count = F, .meat = F, .other = T, .laplace = 1, .sum.col = "Read.count") { freq.segments(.data, .alphabet='TRAV', .count, .meat, .other, .laplace, .sum.col = .sum.col) }
freq.Vb <- function (.data, .count = F, .meat = F, .other = T, .laplace = 1, .sum.col = "Read.count") { freq.segments(.data, .alphabet='TRBV', .count, .meat, .other, .laplace, .sum.col = .sum.col) }
freq.Ja <- function (.data, .count = F, .meat = F, .other = T, .laplace = 1, .sum.col = "Read.count") { freq.segments(.data, .alphabet='TRAJ', .count, .meat, .other, .laplace, .sum.col = .sum.col) }
freq.Jb <- function (.data, .count = F, .meat = F, .other = T, .laplace = 1, .sum.col = "Read.count") { freq.segments(.data, .alphabet='TRBJ', .count, .meat, .other, .laplace, .sum.col = .sum.col) }

freq.segments.2D <- function (.data, .alphabet = "beta", .count = F, .meat = F, .laplace = 1, .columns = NULL, .sum.col = "Read.count", ...) {
  if (has.class(.data, 'list')) {
    return(lapply(.data, freq.segments.2D, .alphabet = .alphabet, .count = .count, .meat = .meat, .laplace = .laplace, .columns = .columns, .sum.col = .sum.col, ...))
  }
  
  if (.alphabet=="beta") {
    .columns <- c('V.segments', 'J.segments')
    alphabetV <- V_BETA_ALPHABET
    alphabetJ <- J_BETA_ALPHABET
  } else if (.alphabet=="alpha") {
    .columns <- c('V.segments', 'J.segments')
    alphabetV <- V_ALPHA_ALPHABET
    alphabetJ <- J_ALPHA_ALPHABET
  } else {
    alphabetV <- .alphabet[[1]]
    alphabetJ <- .alphabet[[2]]
  }
  segV <- c(as.character(.data[, .columns[1]]))
  segJ <- c(as.character(.data[, .columns[2]]))
  counts <- c(rep.int(1, length(segV)))
  if (.meat) counts <- c(.data[, .sum.col])
  
  logic <- !is.na(match(segJ, alphabetJ)) & !is.na(match(segV, alphabetV))
  segV <- segV[logic]
  segJ <- segJ[logic]
  counts <- counts[logic]

  res <- matrix(.laplace, nrow=length(alphabetV), ncol=length(alphabetJ))
  row.names(res) <- alphabetV
  colnames(res) <- alphabetJ
  
  for (i in 1:length(segV)) {
    res[segV[i], segJ[i]] <- res[segV[i], segJ[i]] + counts[i]
  }
  
  if (!.count) { res <- prop.table(res) }
  res <- data.frame(alphabetV, res, stringsAsFactors=F)
  names(res) <- c("Segment", alphabetJ)
  row.names(res) <- NULL
  res
}


#' Perform PCA on segments frequency data.
#' 
#' @aliases pca.segments pca.segments.2D
#' 
#' @description
#' Perform PCA on segments frequency data for V- and J-segments and either return pca object or plot the results.
#' 
#' @usage
#' pca.segments(.data, .cast.freq.seg = T, ..., .do.plot = T)
#' 
#' pca.segments.2D(.data, .cast.freq.seg = T, ..., .do.plot = T)
#' 
#' @param .data Either data.frame or a list of data.frame or a result obtained from freq.segments or freq.segments.2D functions.
#' @param .cast.freq.seg If T than case freq.segments or freq.segments.2D to the supplied data.
#' @param ... Further arguments passed to prcomp.
#' @param .do.plot If T than plot a graphic, else return a pca object.
#' 
#' @return If .do.plot is T than ggplot object; else pca object.
#' 
#' @examples
#' \dontrun{
#' # Load the twins data.
#' data(twb)
#' # Plot a plot of results of PCA on V-segments usage.
#' pca.segments(twb, T, scale. = T)
#' }
pca.segments <- function(.data, .cast.freq.seg = T, ..., .do.plot = T){
  if (.cast.freq.seg) { .data <- freq.segments(.data)[,-1] }
  pca.res <- prcomp(t(as.matrix(.data)), ...)
  if (.do.plot) {
    pca.res <- data.frame(PC1 = pca.res$x[,1], PC2 = pca.res$x[,2], Subject = names(.data))
    ggplot() + geom_point(aes(x = PC1, y = PC2, colour = Subject), size = 3, data = pca.res) +
      geom_text(aes(x = PC1, y = PC2, label = Subject), data = pca.res, hjust=.5, vjust=-.3) +
      theme_linedraw() + guides(size=F) + ggtitle("VJ-usage: Principal Components Analysis") + .colourblind.discrete(length(pca.res$Subject), T)
  } else {
    pca.res
  }
}

pca.segments.2D <- function(.data, .cast.freq.seg = T, ..., .do.plot = T){
  if (.cast.freq.seg) { .data <- lapply(freq.segments.2D(.data), function (x) as.vector(as.matrix(x[,-1]))) }
  pca.res <- prcomp(do.call(rbind, .data), ...)
  if (.do.plot) {
    pca.res <- data.frame(PC1 = pca.res$x[,1], PC2 = pca.res$x[,2], Subject = names(.data))
    ggplot() + geom_point(aes(x = PC1, y = PC2, colour = Subject), size = 3, data = pca.res) +
      geom_text(aes(x = PC1, y = PC2, label = Subject), data = pca.res, hjust=.5, vjust=-.3) +
      theme_linedraw() + guides(size=F) + ggtitle("VJ-usage: Principal Components Analysis") + .colourblind.discrete(length(pca.res$Subject), T)
  } else {
    pca.res
  }
}