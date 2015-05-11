########## Statistics and analysis of Variable and Joining genes usage ##########


if (getRversion() >= "2.15.1") {
  utils::globalVariables(c('PC1', 'PC2'))
}


#' Variable and Joining gene usage.
#' 
#' @aliases geneUsage
#' 
#' @description Get frequencies or counts of gene segments ("V / J - usage").
#' 
#' @usage
#' freq.segments(.data, .genes = 'TRBV', .count = F, .meat = F, .other = F,
#'               .laplace = 0, .column = NULL, .sum.col = "Read.count")
#' 
#' freq.segments.2D(.data, .genes = 'beta', .count = F, .meat = F,
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
#' @param .genes Vector of elements in the alphabet for freq.segments, one of the strings 'TRBV' (for using HUMAN_TRBV_MITCR variable, that user should load before calling functions (same for other strings)), 'TRAV', 'TRBJ', 'TRAJ' for V- and J-segments alphabets for freq.segments
#' or one of the 'alpha' or 'beta' for freq.segments.2D or a list of length 2 with alphabets strings for freq.segments.2D.
#' @param .quant Which column to use for the quantity of clonotypes: NA for computing only number of genes without using clonotype counts, 
#' "read.count" for the "Read.count" column, "bc.count" for the "Barcode.count" column, "read.prop" for the "Read.proportion" column,
#' "bc.prop" for the "Barcode.proportion" column.
#' @param .norm If T then return proportions of resulting counting of genes.
#' 
#' @return 
#' If \code{.data} is a cloneset and \code{.genes} is NOT a list than return a data frame with first column "Gene" with genes and second with counts / proportions.
#' 
#' If \code{.data} is a list with clonesets and \code{.genes} is NOT a list than return a data frame with first column "Gene" 
#' with genes and other columns with counts / proportions for each cloneset in the input list.
#' 
#' If \code{.data} is a cloneset and \code{.genes} IS a list than return a matrix with gene segments for the first gene in \code{.genes}
#' and column names for the second gene in \code{.genes}. See "Examples".
#' 
#' If \code{.data} is a list with clonesets and \code{.genes} IS a list than return a list with matrices like in the previous case.
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
#' # plot V-segments frequencies as a heatmap
#' vis.heatmap(seg, .labs = c("Sample", "V gene"))
#' # plot V-segments frequencies from the data
#' vis.V.usage(twb)
#' 
#' # Compute V-J joint usage.
#' geneUsage(twb, list(HUMAN_TRBV, HUMAN_TRBJ))
#' }
geneUsage <- function (.data, .genes = HUMAN_TRBV, .quant = c(NA, "read.count", "bc.count", "read.prop", "bc.prop"), .norm = F) {
  
  .process.df <- function (.df, .quant, .cols) {
    cast.fun <- dcast
    if (length(.cols) == 2) { cast.fun <- acast; len <- 2 }
    
    for (i in 1:length(.cols)) {
      .df[ which(!(.df[[.cols[i]]] %in% .genes[[i]])), .cols[i] ] <- "Ambiguous"
    }
    
    count.fun <- "n()"
    if (!is.na(.quant)) { count.fun <- paste0("sum(", .quant, ")", collapse = "", sep = "")}
    
    if (length(.cols) == 1) { .cols <- c(.cols, '.'); len <- 1}
    
    cast.fun(summarise_(grouped_df(select_(.df, .dots = as.list(na.exclude(c(.quant, .cols[1:len])))), lapply(.cols[1:len], as.name)), Freq = count.fun), as.formula(paste0(.cols[1], " ~ ", .cols[2])), value.var = 'Freq')
  }
  
  
  quant <- NA
  if (!is.na(.quant[1])) { quant <- .column.choice(.quant, .verbose) }
  
  if (has.class(.data, 'data.frame')) { .data <- list(Sample = .data) }
  
  if (has.class(.genes, 'list')) {    
    genecols <- c(paste0(substr(.genes[[1]][1], 4, 4), ".gene"), paste0(substr(.genes[[2]][1], 4, 4), ".gene"))
  } else {
    genecols <- paste0(substr(.genes[1], 4, 4), ".gene")
    .genes <- list(.genes)
  }
  
  tbls <- lapply(.data, .process.df, .quant = quant, .cols = genecols)
  
  if (length(.genes) == 2) {
    if (.norm) {
      tbls <- lapply(tbls, function (x) x / sum(x))
    }
    if (length(.data) == 1) { return(tbls[[1]]) }
    else { return(tbls) }
  }
  
  res <- tbls[[1]]
  colnames(res) <- c("Gene", names(.data)[1])
  
  if (length(.data) > 1) {
    for (i in 2:length(.data)) {
      colnames(tbls[[i]]) <- c("Gene", names(.data)[i])
      res <- merge(res, tbls[[i]], by = "Gene", all = T)
    }
  }
  
  if (.norm) {
    res[,-1] <- apply(res[,-1], 2, function (col) col / sum(col))
  }
  
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
#' @param .cast.freq.seg if T then case freq.segments or freq.segments.2D to the supplied data.
#' @param ... Further arguments passed to prcomp.
#' @param .do.plot if T then plot a graphic, else return a pca object.
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