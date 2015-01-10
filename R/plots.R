########## Various plotting functions ##########


if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("Segment", 'Size', 'Freq', 'Subject', 'V.segments', 'J.segments', '..count..', 'Time.point', 'Percentage', 'Sequence',
                           'Lower', 'Upper', 'Lengths', 'Read.count', 'Var', 'Value', 'Group', 'variable', 'name', 'value', 'Kmers',
                           'Count', 'People', 'First', 'Second', 'Var1', 'Q0.025', 'Q0.975', 'Mean', 'Type'))
}


# red - yellow - green
.ryg.gradient <- function (.min = NA, .max = NA) {
  cs <- c('#66FF00', '#FFFF66', '#FF6633')
  if (!is.na(.min)) {
    scale_fill_gradientn(limits = c(.min, .max), colours = cs, na.value = 'grey60')
  } else {
    scale_fill_gradientn(colours = cs, na.value = 'grey60')
  }
}

# white/orange/yellow - green - blue
# colourblind - friendly
# for fill
.colourblind.gradient <- function (.min = NA, .max = NA) {
  #   cs <- c("#FFFFD9", "#41B6C4", "#225EA8")
#   cs <- c("#FFFFBB", "#41B6C4", "#225EA8")
  cs <- c("#FFBB00", "#41B6C4", "#225EA8")
  if (!is.na(.min)) {
    scale_fill_gradientn(limits = c(.min, .max), colours = cs, na.value = 'grey60')
  } else {
    scale_fill_gradientn(colours = cs, na.value = 'grey60')
  }
}

# white/orange/yellow - green - blue
# colourblind - friendly
# for fill and colour
.colourblind.discrete <- function (.n, .colour = F) {
  #   cs <- c("#FFFFD9", "#41B6C4", "#225EA8")
  #   cs <- c("#FFFFBB", "#41B6C4", "#225EA8")
  cs <- c("#FFBB00", "#41B6C4", "#225EA8")
  if (.colour) {
    scale_colour_manual(values = colorRampPalette(cs)(.n))
  } else {
    scale_fill_manual(values = colorRampPalette(cs)(.n))
  }
}

# light blues - dark blues
# for fill
.blues.gradient <- function (.min = NA, .max = NA) {
  cs <- c("#F7FBFF", "#9ECAE1", "#2171B5")
  if (!is.na(.min)) {
    scale_fill_gradientn(limits = c(.min, .max), colours = cs, na.value = 'grey60')
  } else {
    scale_fill_gradientn(colours = cs, na.value = 'grey60')
  }
}


#' Plot a histogram of lengths.
#' 
#' @description
#' Plot a histogram of distribution of lengths of CDR3 nucleotide sequences. On y-axis are sum of read counts for each length.
#' 
#' @param .data Data frame with columns 'CDR3.nucleotide.sequence' and 'Read.count' or list with such data frames.
#' @param .ncol If .data is a list, than number of columns in a grid of histograms for each data frame in \code{.data}. Else not used.
#' @param .name Title for this plot.
#' 
#' @details
#' If \code{.data} is a data frame, than one histogram will be plotted. Is \code{.data} is a list, than grid of histograms
#' will be plotted.
#' 
#' @return ggplot object.
#' 
#' @examples
#' \dontrun{
#' load('immdata.rda')
#' # Plot one histogram with main title.
#' vis.count.len(immdata[[1]], 'Main title here')
#' # Plot a grid of histograms with 2 columns.
#' vis.count.len(immdata, 2)
#' }
vis.count.len <- function (.data, .ncol = 3, .name = "") {
  if (has.class(.data, 'list')) {
    return(do.call(grid.arrange, c(lapply(1:length(.data), function (i) vis.count.len(.data[[i]], .name = names(.data)[i])), ncol = .ncol)))
  }
  tmp <- aggregate(Read.count ~ nchar(CDR3.nucleotide.sequence), .data, sum)
  names(tmp) <- c('Lengths', 'Read.count')
  ggplot(tmp) + aes(x = Lengths, y = Read.count, fill = Read.count) +
    geom_histogram(stat = 'identity', colour = 'black') +
    .colourblind.gradient(min(tmp$Read.count), max(tmp$Read.count)) +
    ggtitle(.name) + theme_linedraw()
}


#' Plot a histogram of counts.
#' 
#' @description
#' Plot a histogram of distribution of counts of CDR3 nucleotide sequences. On y-axis are number of counts.
#' 
#' @param .data Data frame with columns 'CDR3.nucleotide.sequence' and 'Read.count' or list with such data frames.
#' @param .ncol If .data is a list, than number of columns in a grid of histograms for each data frame in \code{.data}. Else not used.
#' @param .name Title for this plot.
#' 
#' @details
#' If \code{.data} is a data frame, than one histogram will be plotted. Is \code{.data} is a list, than grid of histograms
#' will be plotted.
#' 
#' @return ggplot object.
#' 
#' @examples
#' \dontrun{
#' load('immdata.rda')
#' # Plot one histogram with main title.
#' vis.number.count(immdata[[1]], 'Main title here')
#' # Plot a grid of histograms with 2 columns.
#' vis.number.count(immdata, 2)
#' }
vis.number.count <- function (.data, .ncol = 3, .name = '') {
#   cat('Limits for x-axis set to (0,50). Transform y-axis to sqrt(y).\n')
  
  if (has.class(.data, 'list')) {
    return(do.call(grid.arrange, c(lapply(1:length(.data), function (i) vis.number.count(.data[[i]], .name = names(.data)[i])), ncol = .ncol)))
  }
  
  ggplot(.data, aes(x = Read.count)) + 
    xlim(min(.data$Read.count), 300) +
    ylab('Number of clones') +
    geom_histogram(aes(fill = ..count..), binwidth = 1, colour = 'black') +
    coord_trans(xtrans = 'log10') + scale_y_log10() +
    ggtitle(.name) + 
    .colourblind.gradient() +
    theme_linedraw()
}


#' Heatmap.
#' 
#' @aliases vis.heatmap
#' 
#' @description
#' Plot a heatmap from a matrix or a data.frame
#'
#' @param .data Either a matrix with colnames and rownames specifyed or a data.frame with the first column of
#' strings for row names and other columns stands for values.
#' @param .title Main title of the plot.
#' @param .labs Labs names. Character vector of length 1 (for naming both axis with same name) or 2 (first elements stands for x-axis).
#' @param .legend Title for the legend.
#' @param .na.value Replace NAs with this values.
#' @param .text If T than print \code{.data} values at tiles.
#' 
#' @return ggplot object.
#' 
#' @examples
#' \dontrun{
#' # Load your data.
#' load('immdata.rda')
#' # Perform intersection by amino acid sequences with V-segments.
#' imm.av <- intersect(immdata, 'ave')
#' # Plot a heatmap.
#' vis.heatmap(imm.av, .title = 'Immdata - (ave)-intersection')
#' }
vis.heatmap <- function (.data, .title = "Number of shared clonotypes", .labs = 'Subject', .legend = 'Shared clonotypes', .na.value = NA, .text = T) {
  if (has.class(.data, 'data.frame')) {
    names <- .data[,1]
    .data <- as.matrix(.data[,-1])
    row.names(.data) <- names
  }

  if (is.null(colnames(.data))) {
    colnames(.data) <- paste0('C', 1:ncol(.data))
  }
  
  .data[is.na(.data)] <- .na.value
  
  tmp <- as.data.frame(.data)
  tmp$name <- row.names(.data)

  m <- melt(tmp, id.var = c('name'))
  m[,1] <- factor(m[,1], levels = rev(rownames(.data)))
  m[,2] <- factor(m[,2], levels = colnames(.data))
  
  p <- ggplot(m, aes(x = variable, y = name, fill = value))
  p <- p + geom_tile(aes(fill = value))
  if (.text) {
    p <- p + geom_text(aes(fill = value, label = value), size = 4)
  }
#   p <- p + geom_text(aes(fill = value, label = value))
#   p <- p + .ryg.gradient(min(m$value), max(m$value))
#   p <- p + .colourblind.gradient(min(m$value), max(m$value))
  p <- p + .blues.gradient(min(m$value), max(m$value))
  p + ggtitle(.title) + 
    guides(fill = guide_legend(title=.legend)) +
    xlab(.labs) + ylab(.labs) + coord_equal() +
    theme(axis.text.x  = element_text(angle=90)) + theme_linedraw() +
    scale_x_discrete(expand=c(0,0)) + scale_y_discrete(expand=c(0,0))
}


#' Boxplot for groups of observations.
#' 
#' @aliases vis.group.boxplot
#' 
#' @description
#' Plot boxplots for each group.
#'
#' @param .data Either a matrix with colnames and rownames specifyed or a data frame with the first column of
#' strings for row names and other columns stands for values.
#' @param .groups Named list with character vector for names of elements for each group. If NA than each
#' member is in the individual group.
#' @param .title Main title of the plot.
#' @param .labs Labs names. Character vector of length 1 (for naming both axis with same name) or 2 (first elements stands for x-axis).
#' @param .rotate.x If T than rotate x-axis.
#' @param ... Parameters passed to \code{melt}, applied to \code{.data} before plotting in \code{vis.group.boxplot}.
#' 
#' @return ggplot object.
#' 
#' @examples
#' \dontrun{
#' names(immdata)  # "A1" "A2" "B1" "B2" "C1" "C2"
#' # Plot a boxplot for V-usage for each plot
#' # three boxplots for each group.
#' vis.group.boxplot(freq.Vb(immdata),
#'    list(A = c('A1', 'A2'), B = c('B1', 'B2'), C = c('C1', 'C2')),
#'    c('V segments', 'Frequency')) 
#' }
vis.group.boxplot <- function (.data, .groups = list(A = c('A1', 'A2'), D = c('D1', 'D2'), C = c('C1', 'C2')), .labs = c('V segments', 'Frequency'), .title = '', .rotate.x = T, ...) {
  .data <- melt(.data, ...)
  
  colnames(.data) <- c('Var', 'Subject', 'Value')
  .data$Group <- as.character(.data$Subject)
  if (!is.na(.groups)[1]) {
    for (i in 1:length(.groups)) {
      for (name in .groups[[i]]) {
        .data$Group[.data$Subject == name] <- names(.groups)[i]
      }
    }
  }
  
  p <- ggplot() + geom_boxplot(aes(x = Var, y = Value, fill = Group), data = .data, colour = 'black')
  if (length(.labs) >= 2) {
    p <- p + xlab(.labs[1]) + ylab(.labs[2])
  }
  p <- p + ggtitle(.title) + theme_linedraw()
  if (.rotate.x) {
    p <- p + theme(axis.text.x  = element_text(angle=90))
  }
  p + .colourblind.discrete(length(unique(.data$Group)))
}


#' Histogram of segments usage.
#' 
#' @aliases vis.V.usage vis.J.usage
#' 
#' @description
#' Plot a histogram or a grid of histograms of V- / J-usage.
#' 
#' @param .data Mitcr data frame or a list with mitcr data frames.
#' 
#' @param .cast.freq If T than cast \code{freq.Vb} (for \code{vis.V.usage}) or \code{freq.Jb} (for \code{vis.J.usage}) on \code{.data} before plotting.
#' @param .main Main title of the plot.
#' @param .ncol Number of columns in a grid of histograms if \code{.data} is a list and \code{.dodge} is F.
#' @param .coord.flip If T than flip coordinates.
#' @param .dodge If \code{.data} is a list, than if this is T plot V-usage for all data frames to the one histogram.
#' @param ... Parameter passed to \code{freq.segments}. By default the function compute V-usage or J-usage for beta chains
#' w/o using read counts and w/ "Other" segments.
#' 
#' @return ggplot object.
#' 
#' @examples
#' \dontrun{
#' # Load your data.
#' load('immdata.rda')
#' # Compute V-usage statistics.
#' imm1.vs <- freq.Vb(immdata[[1]])
#' # Two eqivalent calls for plotting the V-usage for all data frames on the one plot:
#' vis.V.usage(immdata, .cast.freq = T, .main = 'Immdata V-usage [1]', .dodge = T)
#' vis.V.usage(imm1.vs, .cast.freq = F, .main = 'Immdata V-usage [2]', .dodge = T)
#' # Plot a grid of histograms - one histogram for V-usage for each data frame in .data.
#' vis.V.usage(immdata, .cast.freq = T, .main = 'Immdata V-usage [3]', .dodge = F, .other = F)
#' }
vis.V.usage <- function (.data, .cast.freq = T, .main = 'V-usage', .ncol = 3, .coord.flip = F, .dodge = F, ...) {
  if (has.class(.data, 'list')) {
    if (.dodge) {
      res <- melt(freq.Vb(.data))
      colnames(res) <- c('Segment', 'Subject', 'Freq')
      p <- ggplot(res, aes(x = Segment, y = Freq, fill = Subject)) + geom_bar(stat = 'identity', position = position_dodge(), colour = 'black') +
        theme_linedraw() + theme(axis.text.x  = element_text(angle=90)) + scale_fill_brewer(palette = 'YlGnBu')
      return(p)
    } else {
      return(grid.arrange(do.call(arrangeGrob, c(lapply(1:length(.data), function (i) vis.V.usage(.data[[i]], .cast.freq, names(.data)[i], 0, .coord.flip, ...)), ncol = .ncol)), 
                          main = .main))
    }
  }
  
  if (.cast.freq) { 
    if (! '.alphabet' %in% names(list(...))) {
      .data <- freq.segments(.data, 'TRBV', ...)
    } else {
      .data <- freq.segments(.data, ...)
    }
  }
  
  # If result from freq.segmens functions.
  if (names(.data)[1] == 'Segment') {
    p <- ggplot(.data, aes(x = Segment, y = Freq, fill = Freq)) + geom_bar(stat = 'identity', colour = 'black')
  }
  # If mitcr data.frame.
  else {
    p <- ggplot(.data, aes(x = V.segments)) + geom_histogram(aes(fill = ..count..), colour = 'black')
  }
  
  if (.coord.flip) { p <- p + coord_flip() }
  
  p + theme_linedraw() + theme(axis.text.x  = element_text(angle=90)) + ggtitle(.main) +
#     .colourblind.gradient()
   .colourblind.gradient()
#     scale_fill_gradientn(colours = brewer.pal(11, 'Spectral')[c(1,6,11)])
}

vis.J.usage <- function (.data, .cast.freq = T, .main = 'J-usage', .ncol = 3, .coord.flip = F, .dodge = F, ...) {
  if (has.class(.data, 'list')) {
    if (.dodge) {
      res <- melt(freq.Jb(.data))
      colnames(res) <- c('Segment', 'Subject', 'Freq')
      p <- ggplot(res, aes(x = Segment, y = Freq, fill = Subject)) + geom_bar(stat = 'identity', position = position_dodge(), colour = 'black') +
        theme_linedraw() + 
        theme(axis.text.x  = element_text(angle=90)) +
        scale_fill_brewer(palette = 'YlGnBu')
      return(p)
    } else {
      return(grid.arrange(do.call(arrangeGrob, c(lapply(1:length(.data), function (i) vis.J.usage(.data[[i]], .cast.freq, names(.data)[i], 0, .coord.flip, ...)), ncol = .ncol)), 
                          main =  .main, ...))
    }
  }
  
  if (.cast.freq) { 
    if (! '.alphabet' %in% names(list(...))) {
      .data <- freq.segments(.data, 'TRBJ', ...)
    } else {
      .data <- freq.segments(.data, ...)
    }
  }
  
  # If result from freq.segmens functions.
  if (names(.data)[1] == 'Segment') {
    p <- ggplot(.data, aes(x = Segment, y = Freq, fill = Freq)) + geom_bar(stat = 'identity', colour = 'black')
  }
  # If mitcr data.frame.
  else {
    p <- ggplot(.data, aes(x = J.segments)) + geom_histogram(aes(fill = ..count..), , colour = 'black')
  }
  
  if (.coord.flip) { p <- p + coord_flip() }
  
  p + theme_linedraw() + theme(axis.text.x  = element_text(angle=90)) + ggtitle(.main) +
    #     .colourblind.gradient()
    .colourblind.gradient()
  #     scale_fill_gradientn(colours = brewer.pal(11, 'Spectral')[c(1,6,11)])
}


#' PCA result visualisation
#' 
#' @description
#' Plot the given pca results with colour divided by the given groups.
#' 
#' @param .data Result from prcomp() function or a data frame with two columns 'First' and 'Second'
#' stands for the first PC and the second PC.
#' @param .groups List with names for groups and indices of the group members. If NA than each
#' member is in the individual group.
#' 
#' @return ggplot object.
vis.pca <- function (.data, .groups = NA) {
  if (has.class(.data, 'data.frame')) {
    dnames <- row.names(.data)
    .data <- data.frame(First = .data[,1], Second = .data[,2], Subject = row.names(.data),
                        Group = rep('group0', times = length(.data[,2])), stringsAsFactors=F)
  } else {
    dnames <- row.names(.data$x)
    .data <- data.frame(First = .data$x[,1], Second = .data$x[,2], Subject = row.names(.data$x),
                        Group = rep('group0', times = length(.data$x[,2])), stringsAsFactors=F)
  }
  
  if (is.na(.groups)) { 
    .groups <- lapply(1:nrow(.data), function (i) i)
    names(.groups) <- dnames
  }
  for (i in 1:length(.groups)) {
    for (j in 1:length(.groups[[i]])) {
      .data$Group[.groups[[i]][j]] <- names(.groups)[i]
    }
  }
  
  ggplot(data = .data) + 
    geom_point(aes(x = First, y = Second, colour = Group)) + 
    geom_text(aes(x = First, y = Second, label = Subject), .data, hjust=0, vjust=0) +
    theme_linedraw() +
    .colourblind.discrete(length(.groups), T)
}


#' Radar-like / spider-like plots.
#' 
#' @description
#' Plot a grid of radar(-like) plots for visualising a distance among objects.
#' 
#' @param .data Square data frame or matrix with row names and col names stands for objects and values for distances.
#' @param .ncol Number of columns in the grid.
#' @param .expand Interger vector of length 2, for \code{scale_y_continous(expand = .expand)} function.
#' 
#' @seealso \code{intersect} \code{js.div}
#' 
#' @examples
#' \dontrun{
#' load('immdata.rda')
#' # Compute Jensen-Shannon divergence among V-usage of repertoires.
#' imm.js <- js.div.seg(immdata, .verbose = F)
#' # Plot it.
#' vis.radarlike(imm.js)
#' }
vis.radarlike <- function (.data, .ncol = 3, .expand = c(.25, 0)) {
  step = ncol(.data)
  data.names <- colnames(.data)
  .data <- as.data.frame(melt(.data))
  .data[is.na(.data[,3]),3] <- 0
  ps <- lapply(seq(1, nrow(.data), step), function (l) {
    ggplot(.data[l:(l+step-1),], aes(x = Var1, y = value, fill = Var1)) + 
      geom_bar(colour = 'black', stat = 'identity') + 
      coord_polar() + 
      ggtitle(names(.data)[l]) +
      scale_y_continuous(expand = .expand) + 
      guides(fill = guide_legend(title="Subject")) +
      theme_linedraw() + xlab('') + ylab('')
    })
  for (i in 1:length(data.names)) {
    ps[[i]] <- ps[[i]] + ggtitle(data.names[i]) + .colourblind.discrete(length(data.names))
  }
  
  do.call(grid.arrange, c(ps, ncol = .ncol))
}


#' Visualisation of top clones proportions.
#' 
#' @description
#' Visualisation of proportion of the top clones.
#' 
#' @param .data Data frame with clones.
#' @param .head Integer vector of clones for the \code{.head} parameter for the \code{top.proportion} function.
#' @param .col Parameter \code{.col} for the \code{top.proportion} function.
#' 
#' @seealso \code{top.proportion}
#' 
#' @examples
#' \dontrun{
#' vis.top.proportions(immdata)
#' }
vis.top.proportions <- function (.data, .head = c(10, 100, 1000, 10000, 30000, 100000, 300000, 1000000), .col = "Read.count") {
  if (has.class(.data, 'data.frame')) {
    .data <- list(Data = .data)
  }
  
  res <- sapply(.head, function (h) top.proportion(.data, h, .col))
  tmp <- res
  if (is.null(dim(tmp))) {
    tmp <- t(as.matrix(tmp))
    res <- t(as.matrix(res))
  }
  for (i in 2:ncol(res)) {
    tmp[,i] <- res[,i] - res[,i-1]
  }
  res <- tmp
  colnames(res) <- paste0('[', c(1, .head[-length(.head)] + 1), ':', .head, ')')
  res <- as.data.frame(res)
  res$People <- factor(row.names(res), levels = row.names(res))
  res <- melt(res)
  #   res$variable <- factor(as.character(res$variable), labels = paste0('[', c(1, .head[-length(.head)] + 1), ':', .head, ')'), ordered = T)
  ggplot() + geom_bar(aes(x = People, y = value, fill = variable),data = res, stat = 'identity', position = 'stack', colour = 'black')+ 
    theme_linedraw()  + 
    theme(axis.text.x  = element_text(angle=90)) +
    ylab("Clonal percentage") + 
    xlab("Subject") + 
    ggtitle("Summary percentage of the top N clones")  + 
    guides(fill = guide_legend("Top N clones")) + .colourblind.discrete(length(.head))
}


#' Rarefaction statistics visualisation.
#' 
#' @description
#' Plot a line with mean unique clones.
#' 
#' @param .muc.res Output from the \code{muc} function.
#' @param .groups List with names for groups and names of the group members. If NULL than each
#' member is in the individual group.
#' @param .log If T than log-scale the y axis.
#' 
#' @seealso \link{rarefaction}
#' 
#' @examples
#' \dontrun{
#' data(twb)
#' names(twb)  # "Subj.A" "Subj.B" "Subj.C" "Subj.D"
#' twb.rar <- rarefaction(twb, .col = "Read.count")
#' vis.rarefaction(twb.rar, list(A = c("Subj.A", "Subj.B"), B = c("Subj.C", "Subj.D")))
#' }
vis.rarefaction <- function (.muc.res, .groups = NULL, .log = F) {
  .muc.res$Group <- .muc.res$People
  
  if (!is.null(.groups)) { 
    for (i in 1:length(.groups)) {
      for (j in 1:length(.groups[[i]])) {
        .muc.res$Group[ .muc.res$People == .groups[[i]][j] ] <- names(.groups)[i]
      }
    }
  }
  
  .muc.res$Type <- factor(.muc.res$Type, levels = c('interpolation', 'extrapolation'), ordered = T)
  
  p <- ggplot() + 
#     geom_point(aes(x = Size, y = Mean, colour = Group), data = .muc.res, size = 2) + 
    geom_line(aes(x = Size, y = Mean, colour = Group, Group = People, linetype = Type), data = .muc.res) + 
#     geom_errorbar(aes(x = Size, y = Mean, ymin = Q0.025, ymax = Q0.975, colour = Group), data = .muc.res) +
    xlab('Sample size') + ylab('Clones') + ggtitle("Rarefaction analysis") +
    theme_linedraw() + .colourblind.discrete(length(unique(.muc.res$Group)), T)
  
  for (subj in unique(.muc.res$People)) {
    tmp <- tail(.muc.res[.muc.res$People == subj, ], 1)
    p <- p + geom_text(aes(x = Size, y = Mean, label = People), data = tmp, hjust=1, vjust=1)
  }
  
  if (.log) {
    p <- p + scale_x_log10()
  }
  p
}


#' Plot of the most frequent kmers.
#' 
#' @description
#' Plot a distribution (bar plot) of the most frequent kmers in a data.
#' 
#' @param .kmers Data frame with two columns "Kmers" and "Count" or a list with such data frames. See Examples.
#' @param .head Number of the most frequent kmers to choose for plotting from each data frame.
#' @param .position Character vector of length 1. Position of bars for each kmers. Value for the \code{ggplot2} argument \code{position}.
#' 
#' @seealso \code{get.kmers}
#' 
#' @examples
#' \dontrun{
#' # Load necessary data and package.
#' library(gridExtra)
#' load('immdata.rda')
#' # Get 5-mers.
#' imm.km <- get.kmers(immdata)
#' # Plots for kmer proportions in each data frame in immdata.
#' p1 <- vis.kmer.histogran(imm.km, .position = 'stack')
#' p2 <- vis.kmer.histogran(imm.km, .position = 'fill')
#' grid.arrange(p1, p2)
#' }
vis.kmer.histogram <- function (.kmers, .head = 100, .position = c('stack', 'dodge', 'fill')) {
  kmers.df <- data.frame(Kmers = '')
  for (i in 2:ncol(.kmers)) {
    kmers.df <- merge(head(.kmers[order(.kmers[, i], decreasing = T), c(1,i)], .head), kmers.df, all = T)
  }
  kmers.df[is.na(kmers.df)] <- 0
  kmers.df <- melt(kmers.df[-1,])
  names(kmers.df) <- c('Kmers', 'People', 'Count')
  p <- ggplot() + geom_bar(aes(x = Kmers, y = Count, fill = People), data = kmers.df, stat = 'identity', position = .position[1]) + theme_linedraw()
  if (.position[1] == 'stack' || .position[1] == 'dodge') {
    p <- p + ylab('Count') + theme(axis.text.x  = element_text(angle=90))
  } else {
    p <- p + ylab('Proportions') + theme(axis.text.x  = element_text(angle=90))
  }
  p + scale_y_continuous(expand = c(0, 0)) + .colourblind.discrete(length(unique(kmers.df)))
}


#' Visualise clonal dynamics among time points.
#' 
#' @description
#' Visualise clonal dynamics (i.e., changes in frequency or count) with error bars of given
#' clones among time points.
#' 
#' @param .changed Result from the \code{find.clonotypes} function, i.e. data frame with first
#' columns with sequences (nucleotide or amino acid) and other columns are columns with frequency / count
#' for each time point for each clone.
#' @param .lower Similar to .changed but values are lower bound for clonal count / frequency.
#' @param .upper Similar to .changed but values are upper bound for clonal count / frequency.
#' @param .log If T than log-scale y-axis.
#' 
#' @return ggplot object.
vis.clonal.dynamics <- function (.changed, .lower, .upper, .log = T) {
  .changed <- melt(.changed, id.vars = names(.changed)[1])
  .lower <- melt(.lower, id.vars = names(.changed)[1])
  .upper <- melt(.upper, id.vars = names(.changed)[1])
  names(.changed) <- c('Sequence', 'Time.point', 'Percentage')
  d <- cbind(.changed, Lower = .lower[,3], Upper = .upper[,3])
  p <- ggplot() + geom_line(aes(x = Time.point, y = Percentage, colour = Sequence, group = Sequence), data = d) +
    geom_errorbar(aes(x = Time.point, y = Percentage, colour = Sequence, ymin = Lower, ymax = Upper), data = d, width = .25) +
    theme_linedraw() + theme(axis.text.x  = element_text(angle=90)) +
    .colourblind.discrete(length(unique(.changed$Sequence)), .colour = T)
  if (.log) {
    p <- p + scale_y_log10()
  }
  p
}