########## Statistics and analysis of Variable and Joining genes usage ##########


if (getRversion() >= "2.15.1") {
  utils::globalVariables(c('PC1', 'PC2'))
}


#' V- and J-segments frequency.
#' 
#' @aliases freq.segments freq.segments.2D freq.Va freq.Vb freq.Ja freq.Jb
#' 
#' @description Get frequencies or counts of gene segments ("V / J - usage").
#' 
#' @usage
#' freq.segments(.data, .alphabet = 'TRBV', .count = F, .meat = F, .other = F,
#'               .laplace = 0, .column = NULL, .sum.col = "Read.count")
#' 
#' freq.segments.2D(.data, .alphabet = 'beta', .count = F, .meat = F,
#'                  .laplace = 0, .columns = NULL, .sum.col = "Read.count", ...)
#' 
#' freq.Va(.data, .count = F, .meat = F, .other = F, .laplace = 0, .sum.col = "Read.count")
#' 
#' freq.Vb(.data, .count = F, .meat = F, .other = F, .laplace = 0, .sum.col = "Read.count")
#' 
#' freq.Ja(.data, .count = F, .meat = F, .other = F, .laplace = 0, .sum.col = "Read.count")
#' 
#' freq.Jb(.data, .count = F, .meat = F, .other = F, .laplace = 0, .sum.col = "Read.count")
#' 
#' @param .data Cloneset data frame or a list with clonesets.
#' @param .alphabet Vector of elements in the alphabet for freq.segments, one of the strings 'TRBV' (for using HUMAN_TRBV_MITCR variable, that user should load before calling functions (same for other strings)), 'TRAV', 'TRBJ', 'TRAJ' for V- and J-segments alphabets for freq.segments
#' or one of the 'alpha' or 'beta' for freq.segments.2D or a list of length 2 with alphabets strings for freq.segments.2D.
#' @param .count Should we return count or percentage?
#' @param .meat Compute statistics using counts of elements (e.g., Read.count) or not.
#' @param .other Should elements not in the given alphabet be shown in the result matrix.
#' @param .laplace Value for the Laplace correction.
#' @param .column,.columns Column's name with elements from the given alphabet for the freq.segments or a character vector of length two for freq.segments.2D.
#' @param .quant Which column to use for the quantity of clonotypes: NA for computing only number of genes without using clonotype counts, 
#' "read.count" for the "Read.count" column, "bc.count" for the "Barcode.count" column, "read.prop" for the "Read.proportion" column,
#' "bc.prop" for the "Barcode.proportion" column.
#' @param ... Don't use it, for internal purpose.
#' 
#' @return 
#' If \code{.data} is a cloneset and \code{.alphabet} is NOT a list than return a data frame with first column "Gene" with genes and second with counts / proportions.
#' 
#' If \code{.data} is a list with clonesets and \code{.alphabet} is NOT a list than return a data frame with first column "Gene" 
#' with genes and other columns with counts / proportions for each cloneset in the input list.
#' 
#' If \code{.data} is a cloneset and \code{.alphabet} IS a list than return a matrix with gene segments for the first gene in \code{.alphabet}
#' and column names for the second gene in \code{.alphabet}. See "Examples".
#' 
#' If \code{.data} is a list with clonesets and \code{.alphabet} IS a list than return a list with matrices like in the previous case.
#' 
#' @seealso \code{\link{genealphabets}}, \code{\link{vis.gene.usage}}, \code{\link{pca.segments}}
#' 
#' @examples
#' \dontrun{
#' # Load your data
#' data(twb)
#' # compute V-segments frequencies of human TCR beta.
#' seg <- freq.segments(twb, "TRBV")
#' # equivalent to the previos one
#' seg <- freq.segments(twb, HUMAN_TRBV, .column = "V.gene")
#' # plot V-segments frequencies as a grid
#' vis.grid.stats(seg)
#' # plot V-segments frequencies from the data
#' vis.V.usage(twb)
#' 
#' # Compute V-J joint usage.
#' geneUsage(twb, list(HUMAN_TRBV, HUMAN_TRBJ))
#' }
geneUsage <- function (.data, .genes = HUMAN_TRBV, .quant = c(NA, "read.count", "bc.count", "read.prop", "bc.prop"), .other = F, .laplace = 0) {
  
  .process.df <- function (.df, .quant, .cols) {
    cast.fun <- dcast
    if (length(.cols) == 2) { cast.fun <- acast }
    
    count.fun <- "n()"
    if (!is.na(.quant)) { count.fun <- paste0("sum(", .quant, ")", collapse = "", sep = "")}
    
    if (length(.cols) == 1) { .cols <- c(.cols, '.')}
    
    cast.fun(summarise_(grouped_df(.df, lapply(.cols, as.name)), Freq = count.fun), as.formula(paste0(.cols[1], "~", .cols[2])), value.var = 'Freq')
  }
  
  
  quant <- NA
  if (!is.na(.quant)) { quant <- .column.choice(.quant, .verbose) }
  
  if (has.class(.data, 'data.frame')) { .data <- list(Sample = .data) }
  
  if (has.class(.alphabet, 'list')) {    
    genecols <- c(paste0(substr(.alphabet[[1]], 3, 3), ".gene"), paste0(substr(.alphabet[[2]], 3, 3), ".gene"))
  } else {
    genecol1 <- paste0(substr(.alphabet[[1]], 3, 3), ".gene")
  }
  
  tbls <- lapply(.data, .process.df, .quant = quant, .cols = genecols)
  
  if (.other) {
    for (i in 1:length(.data)) {
      
    }
  }
}


freq.segments <- function (.data, .alphabet='TRBV', .count=F, .meat=F, .other=F, .laplace=0, .column = NULL, .sum.col = "Read.count") {
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
  
  if (.alphabet == "TRBV")      { .column <- 'V.gene'; alphabet <- HUMAN_TRBV_MITCR }
  else if (.alphabet == "TRAV") { .column <- 'V.gene'; alphabet <- HUMAN_TRAV }
  else if (.alphabet == "TRBJ") { .column <- 'J.gene'; alphabet <- HUMAN_TRBJ }
  else if (.alphabet == "TRAJ") { .column <- 'J.gene'; alphabet <- HUMAN_TRAJ }
  else                          { alphabet <- .alphabet }
  seg <- .data[[.column]]
  
  read.count <- rep.int(1, nrow(.data))
  if (.meat) { read.count <- .data[, .sum.col] }
  
  counts.l <- tapply(read.count, seg, sum)
  freqs <- as.numeric(counts.l[alphabet])
  freqs[is.na(freqs)] <- 0
  res <- data.frame(Segment = alphabet, Freq = freqs, stringsAsFactors = F, row.names=NULL)
  if (.other) {
    res <- rbind(res, Other = data.frame(Segment = 'Other', Freq = sum(unlist(lapply(counts.l[is.na(match(names(counts.l), alphabet))], sum)), na.rm=T), stringsAsFactors = F))
  }
  res$Freq <- res$Freq + .laplace
  if (!.count) { res$Freq <- res$Freq / sum(res$Freq) }
  res[order(res$Segment),]
}

freq.Va <- function (.data, .count = F, .meat = F, .other = F, .laplace = 0, .sum.col = "Read.count") { freq.segments(.data, .alphabet='TRAV', .count, .meat, .other, .laplace, .sum.col = .sum.col) }
freq.Vb <- function (.data, .count = F, .meat = F, .other = F, .laplace = 0, .sum.col = "Read.count") { freq.segments(.data, .alphabet='TRBV', .count, .meat, .other, .laplace, .sum.col = .sum.col) }
freq.Ja <- function (.data, .count = F, .meat = F, .other = F, .laplace = 0, .sum.col = "Read.count") { freq.segments(.data, .alphabet='TRAJ', .count, .meat, .other, .laplace, .sum.col = .sum.col) }
freq.Jb <- function (.data, .count = F, .meat = F, .other = F, .laplace = 0, .sum.col = "Read.count") { freq.segments(.data, .alphabet='TRBJ', .count, .meat, .other, .laplace, .sum.col = .sum.col) }

freq.segments.2D <- function (.data, .alphabet = "beta", .count = F, .meat = F, .laplace = 0, .columns = NULL, .sum.col = "Read.count", ...) {
  if (has.class(.data, 'list')) {
    return(lapply(.data, freq.segments.2D, .alphabet = .alphabet, .count = .count, .meat = .meat, .laplace = .laplace, .columns = .columns, .sum.col = .sum.col, ...))
  }
  
  if (.alphabet=="beta") {
    .columns <- c('V.gene', 'J.gene')
    alphabetV <- HUMAN_TRBV_MITCR
    alphabetJ <- HUMAN_TRBJ
  } else if (.alphabet=="alpha") {
    .columns <- c('V.gene', 'J.gene')
    alphabetV <- HUMAN_TRAV
    alphabetJ <- HUMAN_TRAJ
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