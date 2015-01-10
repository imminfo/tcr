########## Shared TCR repertoire managing and analysis ##########


#' Shared TCR repertoire managing and analysis
#' 
#' @aliases shared.repertoire shared.matrix
#' 
#' @description
#' Generate a repertoire of shared sequences - sequences presented in more than one subject. If sequence is appeared more than once in the one
#' repertoire, than only the first appeared one will be choosed for a shared repertoire.
#' 
#' \code{shared.repertoire} - make a shared repertoire of sequences from the given list of data frames.
#' 
#' \code{shared.matrix} - leave columns, which related to the count of sequences in people, and return them as a matrix.
#' I.e., this functions will remove such columns as 'CDR3.amino.acid.sequence', 'V.segments', 'People'.
#' 
#' @usage
#' shared.repertoire(.datalist, .type = 'avc', .min.ppl = 1, .head = -1,
#'                   .clear = T, .verbose = T, .by.col = '', .sum.col = '',
#'                   .max.ppl = length(.datalist))
#' 
#' shared.matrix(.shared.rep)
#' 
#' @param .datalist List with data frames.
#' @param .type String of length 3 denotes how to create a shared repertoire. See "Details" for
#' more information. If supplied, than parameters \code{.by.col} and \code{.sum.col} will be ignored. If not supplied, than columns
#' in \code{.by.col} and \code{.sum.col} will be used.
#' @param .min.ppl At least how many people must have a sequence to leave this sequence in the shared repertoire.
#' @param .head Parameter for the \code{head} function, applied to all data frames before clearing.
#' @param .clear If T than remove all sequences which have symbols "~" or "*" (i.e., out-of-frame sequences for amino acid sequences).
#' @param .verbose If T than output progress.
#' @param .by.col Character vector with names of columns with sequences and their parameters (like segment) for using for creating a shared repertoire.
#' @param .sum.col Character vector of length 1 with names of the column with count, percentage or any other numeric chaaracteristic of sequences for using for creating a shared repertoire.
#' @param .max.ppl At most how many people must have a sequence to leave this sequence in the shared repertoire.
#' @param .shared.rep Shared repertoire.
#' 
#' @details
#' Parameter \code{.type} is a string of length 3, where:
#' \enumerate{
#'  \item First character stands either for the letter 'a' for taking the "CDR3.amino.acid.sequence" column or
#' for the letter 'n' for taking the "CDR3.nucleotide.sequence" column.
#'  \item Second character stands whether or not take the V.segments column. Possible values are '0' (zero) stands
#' for taking no additional columns, 'v' stands for taking the "V.segments" column.
#'  \item Third character stands for name of the column to choose as numeric characteristic of sequences. Possible values are
#' "c" for the "Read.count" column, "p" for the "Percentage" column, "r" for the "Rank" column or "i" for the "Index" column.
#' If "Rank" or "Index" isn't in the given repertoire, than it will be created using \code{set.rank} function using default "Read.count" column.
#' }
#' 
#' @return
#' Data.table for \code{shared.repertoire}, matrix for \code{shared.matrix}.
#' 
#' @seealso \link{shared.representation}, \link{set.rank}
#' 
#' @examples
#' \dontrun{
#' # Set "Rank" column in data by "Read.count" column.
#' # This is doing automatically in shared.repertoire() function
#' # if the "Rank" column hasn't been found.
#' immdata <- set.rank(immdata)
#' # Generate shared repertoire using "CDR3.amino.acid.sequence" and
#' # "V.segments" columns and with rank.
#' imm.shared.av <- shared.repertoire(immdata, 'avr')
#' }
shared.repertoire <- function (.datalist, .type = 'avc', .min.ppl = 1, .head = -1, .clear = T,
                               .verbose = T, .by.col = '', .sum.col = '', .max.ppl = length(.datalist)) {
  
  .process.df <- function (.data, .bc, .sc) {
    if (.head != -1) {
      .data <- head(.data, .head)
    }
    
    if (.clear) {
      .data <- .data[grep('[*, ~]', .data[, .bc[1]], invert = T), ] 
    }
    
    # If .data$Rank or .data$Index is NULL, than generate this columns
    # using the "Read.count" column.
    if (is.null(.data[[.sc]])) {
      if (.sc == 'Rank') {
        .data <- set.rank(.data)
      } else {
        .data <- set.index(.data)
      }
    }
    
#     minidata <- as.data.table(.data[.bc])
    minidata <- as.data.table(.data[.bc])
    minidata$value <- .data[, .sc]
    res <- as.data.table(dplyr::summarise(grouped_df(minidata, lapply(.bc, as.name)), value = value[1]))
    class(res) <- c('data.table', 'data.frame')
    setnames(res, c(.bc, .sc))
    res
  }
  
  
  if (nchar(.by.col[1]) == 0) {
    if (substr(.type, 1, 1) == 'a') { .by.col <- 'CDR3.amino.acid.sequence' }
    else { .by.col <- 'CDR3.nucleotide.sequence' }
    
    if (substr(.type, 2, 2) == 'v') { .by.col <- c(.by.col, 'V.segments') }
  }
  
  if (nchar(.sum.col) == 0) {
    # read count
    if (substr(.type, 3, 3) == 'c') {
      .sum.col <- 'Read.count'
    } else if (substr(.type, 3, 3) == 'p') {
      .sum.col <- 'Percentage'
    } else if (substr(.type, 3, 3) == 'r') {
      .sum.col <- 'Rank'
    } else if (substr(.type, 3, 3) == 'i') {
      .sum.col <- 'Index'
    } else {
      # As a default option.
      .sum.col <- 'Read.count'
    }
  }
  
  if (.verbose) { cat("Aggregating sequences...\n"); pb <- set.pb(length(.datalist)) }
  
  l <- list()
  for (i in 1:length(.datalist)) {
    l[[i]] <- .process.df(.datalist[[i]], .by.col, .sum.col)
    if (.verbose) add.pb(pb)
  }
  names(l) <- names(.datalist)
  if (.verbose) close(pb)
  
#   l <- lapply(.datalist, .process.df, .bc = .by.col, .sc = .sum.col)
  res <- l[[1]]
  setnames(res, c(.by.col, names(.datalist)[1]))
  
  if (.verbose) { cat("Merging data tables...\n"); pb <- set.pb(length(.datalist) - 1) }
  
  if (length(l) > 1) {
    for (i in 2:length(l)) {
      res <- merge(res, l[[i]], by = .by.col, all=T, allow.cartesian=T)
      setnames(res, c(colnames(res)[-ncol(res)], names(.datalist)[i]))
      if (.verbose) add.pb(pb)
    }
  }
  
  if (.verbose) { close(pb)}
  
  res$People <- rowSums(!is.na(as.matrix(res[, (ncol(res) - length(l) + 1) : ncol(res), with = F])))
  res <- res[People >= .min.ppl & People <= .max.ppl][order(People, decreasing=T)][, c(1:(ncol(res) - length(l) - 1), ncol(res), (ncol(res) - length(l)) : (ncol(res) - 1)), with = F]
  setattr(res, 'by.col', .by.col)
  setattr(res, 'sum.col', .sum.col)
  res
}

shared.matrix <- function (.shared.rep) {
  as.matrix(.shared.rep[, -(1:(match('People', colnames(.shared.rep)))), with=F])
}


#' Shared repertoire analysis.
#' 
#' @aliases cosine.sharing shared.representation shared.clones.count shared.summary
#' 
#' @description
#' Functions for computing statistics and analysis of shared repertoire of sequences.
#' 
#' \code{cosine.sharing} - apply the cosine similarity measure to the vectors of sequences' counts or indices.
#' 
#' \code{shared.representation} - for every repertoire in the shared repetoire get a number of sequences in this repertoire which are in the other repertoires.
#' Row names of the input matrix is the number of people.
#' 
#' \code{shared.clones.count} - get the number of shared clones for every number of people.
#' 
#' \code{shared.summary} - get a matrix with counts of pairwise shared sequences (like a result from \code{cross} function, applied to a list of data frames).
#' 
#' @usage
#' cosine.sharing(.shared.rep, .log = T)
#' 
#' shared.representation(.shared.rep)
#' 
#' shared.clones.count(.shared.rep)
#' 
#' shared.summary(.shared.rep, .min.ppl = min(.shared.rep$People),
#'                .max.ppl = max(.shared.rep$People))
#' 
#' @param .shared.rep Shared repertoire, obtained from the function \code{shared.repertoire}.
#' @param ... Parameters passed to the \code{prcomp} function.
#' @param .log If T than apply log to the after adding laplace correction equal to one.
#' @param .min.ppl Filter: get sequences with # people >= .min.ppl.
#' @param .max.ppl Filter: get sequences with # people <= .max.ppl.
#' 
#' @return Plot or PCA resulr for the \code{shared.seq.pca} function or a matrix with cosine similarity values for the \code{cosine.sharing} function.
#' 
#' @seealso \link{shared.repertoire}
#' 
#' @examples
#' \dontrun{
#' # Load the twb data.
#' data(twb)
#' # Create shared repertoire on the twins data using CDR3 amino acid sequences with CDR1-2.
#' twb.shared <- shared.repertoire(twb, 'av', .verbose = T)
#' sh.repr <- shared.representation(twb.shared)
#' sh.repr
#' # Get proportion of represented shared sequences.
#' apply(sh.repr, 2, function (col) col / col[1])
#' }
cosine.sharing <- function (.shared.rep, .log = T) {
  x <- shared.matrix(.shared.rep)
  d <- lapply(seq_len(ncol(x)), function(i) x[,i])
  d <- lapply(d, function (l) {
    l <- unlist(l)
    l[is.na(l)] <- 0
    l <- l + 1
    if (.log) { l <- log(l) }
    l
  })
  apply.symm(d, cosine.similarity, .do.norm = F, .verbose = F)
}

shared.representation <- function (.shared.rep) {
  mat <- shared.matrix(.shared.rep)
  ppl <- as.integer(.shared.rep$People)
  apply(mat, 2, function (col) {
    table(c(1:ncol(mat), ppl[!is.na(col) & (col > 0)])) - 1
  })
}

shared.clones.count <- function (.shared.rep) {
  table(as.integer(.shared.rep$People))
}

shared.summary <- function (.shared.rep, .min.ppl = min(.shared.rep$People), .max.ppl = max(.shared.rep$People)) {
  x <- apply(shared.matrix(.shared.rep), 2, function (col) which(!is.na(col) & col != 0))
  if (class(x) != 'list') {
    x <- lapply(seq_len(ncol(x)), function(i) x[,i])
  }
  apply.symm(x, function (a, b) length(intersect(x = a, y = b)))
}