########## Data processing functions ##########


#' Parse input table files with immune receptor repertoire data.
#'
#' @description
#' General parser for cloneset table files. 
#'
#' @param .filepath Path to the input file with cloneset data.
#' 
#' @return Data frame with immune receptor repertoire data. See \link{\code{parse.file}} for more details.
#' 
#' @seealso \link{\code{parse.file}}
#' 
#' @examples
#' \dontrun{
#' # Parse file in "~/mitcr/immdata1.txt" as a MiTCR file.
#' immdata1 <- parse.file("~/mitcr/immdata1.txt", 'mitcr')
#' }
parse.cloneset <- function (.filename,
                            .nuc.seq,
                            .aa.seq,
                            .count,
                            .proportion,
                            .reads,
                            .events,
                            .vgenes,
                            .jgenes,
                            .dgenes,
                            .vend,
                            .jstart,
                            .dalignments,
                            .vd.insertions,
                            .dj.insertions,
                            .total.insertions,
                            .skip = 0,
                            .sep = '\t') {
  
  .dalignments1 <- .dalignments[1]
  .dalignments2 <- .dalignments[2]
  
  table.colnames <- read.table(.filename, sep = .sep, skip = .skip, nrows = 1, stringsAsFactors = F, strip.white = T)[1,]
  
  swlist <- list('character', 'character',
                 'integer', 'numeric',
                 'integer', 'integer',
                 'character', 'character', 'character',
                 'integer', 'integer', 'integer', 'integer',
                 'integer', 'integer', 'integer')
  names(swlist) <- c(.nuc.seq, .aa.seq,
                     .count, .proportion,
                     .reads, .events,
                     .vgenes, .jgenes, .dgenes,
                     .vend, .jstart, .dalignments,
                     .vd.insertions, .dj.insertions, .total.insertions)
  swlist <- c(swlist, 'NULL')
  
  col.classes <- unlist(sapply(table.colnames, function (x) {
    do.call(switch, c(x, swlist))
  }, USE.NAMES = F))
  
  suppressWarnings(df <- read.table(file = .filename, header = T, colClasses = col.classes, sep = .sep, skip = .skip, strip.white = T))
  
  if(is.na(.events)) {
    .events <- "Event.count"
    df$Events <- -1
  }
  
  if (is.na(.aa.seq)) {
    df$CDR3.amino.acid.sequence <- bunch.translate(df$CDR3.nucleotide.sequence)
    .aa.seq <- 'CDR3 amino acid sequence'
  }
  
  if (is.na(.proportion)) {
    df$Proportion <- df$Count / sum(df$Count)
    .proportion <- 'Proportion'
  }
  
  if (!(.total.insertions %in% table.colnames)) {
    if (.vd.insertions %in% table.colnames && .dj.insertions %in% table.colnames) {
      df$Total.insertions <- df$VD.insertions + df$DJ.insertions
    } else {
      df$Total.insertions <- -1
    }
  }
  
  if (!(.vd.insertions %in% table.colnames)) { df$VD.insertions <- -1 }
  
  if (!(.dj.insertions %in% table.colnames)) { df$DJ.insertions <- -1 }
  
  df <- df[, make.names(c(.count, .proportion, .nuc.seq, .aa.seq,
                          .vgenes, .jgenes, .dgenes,
                          .vend, .jstart, .dalignments,
                          .vd.insertions, .dj.insertions, .total.insertions,
                          .reads, .events))]
  
  colnames(df) <- c('Count', 'Proportion', 'CDR3.nucleotide.sequence', 'CDR3.amino.acid.sequence',
                    'V.segments', 'J.segments', 'D.segments',
                    'V.end', 'J.start', 'D5.end', 'D3.end',
                    'VD.insertions', 'DJ.insertions', 'Total.insertions',
                    'Read.count', 'Event.count')
  
  df
}


#' Parse input table files with immune receptor repertoire data.
#'
#' @aliases parse.folder parse.file.list parse.file parse.mitcr parse.mitcrbc parse.migec
#'
#' @description
#' Load the MITCR TCR data from the file with the given filename
#' to a data frame. For general parser see \code{\link{parse.cloneset}}.
#'
#' @param .filepath Path to the input file with cloneset data.
#' @param .filenames Vector or list with paths to files with cloneset data.
#' @param .folderpath Path to the folder with text cloneset files.
#' @param .format String specifing input format of files. Parsers for MiTCR output and MiGEC output are available.
#' @param ... Parameters passed to \code{parse.cloneset}.
#' 
#' @return Data frame with immune receptor repertoire data and following columns:
#' 
#' @seealso \link{\code{parse.cloneset}}
#' 
#' @examples
#' \dontrun{
#' # Parse file in "~/mitcr/immdata1.txt" as a MiTCR file.
#' immdata1 <- parse.file("~/mitcr/immdata1.txt", 'mitcr')
#' # Parse files "~/data/immdata1.txt" and "~/data/immdat2.txt" as MiGEC files.
#' immdata12 <- parse.file.list(c("~/data/immdata1.txt",
#'                              "~/data/immdata2.txt"), 'migec')
#' # Parse all files in "~/data/" as MiGEC files.
#' immdata <- parse.folder("~/data/", 'migec')
#' }
parse.folder <- function (.folderpath, .format = c('mitcr', 'mitcrbc', 'migec'), ...) {
  parse.file.list(list.files(.folderpath, full.names = T), .format)
}

parse.file.list <- function (.filenames, .format = c('mitcr', 'mitcrbc', 'migec'), .namelist = NA) {
  # Remove full paths and extension from the given string.
  .remove.ext <- function (.str) {
    gsub(pattern = '.*/|[.].*$', replacement = '', x = .str)
    #     .str
  }
  
  .filenames <- as.list(.filenames)
  
  datalist <- lapply(X = 1:length(.filenames), FUN = function (i) parse.file(.filenames[[i]], .format) )
  if (is.na(.namelist)) {
    namelist <- lapply(X = .filenames, FUN = .remove.ext)
    names(datalist) <- unlist(namelist)
  }
  datalist
}

parse.file <- function(.filename, .format = c('mitcr', 'mitcrbc', 'migec'), ...) {
  
  parse.fun <- switch(.format[1], 
                      mitcr = parse.mitcr,
                      mitcrbc = parse.mitcrbc,
                      migec = parse.migec,
                      parse.cloneset)
  
  parse.fun(.filename, ...)
}

parse.mitcr <- function (.filename) {
  
  filename <- .filename
  nuc.seq <- 'CDR3 nucleotide sequence'
  aa.seq <- 'CDR3 amino acid sequence'
  count <- 'Read count'
  Proportion <- 'Percentage'
  reads <- 'Read count'
  events <- NA
  vgenes <- 'V segments'
  jgenes <- 'J segments'
  dgenes <- 'D segments'
  vend <- 'Last V nucleotide position'
  jstart <- 'First J nucleotide position'
  dalignments <- c('First D nucleotide position', 'Last D nucleotide position')
  vd.insertions <- 'VD insertions'
  dj.insertions <- 'DJ insertions'
  total.insertions <- 'Total insertions'
  .skip = 1
  .sep = '\t'
    
  parse.cloneset(.filename = filename, 
                 .nuc.seq = nuc.seq,
                 .aa.seq = aa.seq,
                 .count = count,
                 .proportion = Proportion,
                 .reads = reads,
                 .events = events,
                 .vgenes = vgenes,
                 .jgenes = jgenes,
                 .dgenes = dgenes,
                 .vend = vend,
                 .jstart = jstart,
                 .dalignments = dalignments,
                 .vd.insertions = vd.insertions,
                 .dj.insertions = dj.insertions,
                 .total.insertions = total.insertions,
                 .skip = .skip,
                 .sep = .sep)
}

parse.mitcrbc <- function (.filename) {
  
  filename <- .filename
  nuc.seq <- 'CDR3 nucleotide sequence'
  aa.seq <- 'CDR3 amino acid sequence'
  count <- 'NNNs'
  Proportion <- NA
  reads <- 'Count'
  events <- 'NNNs'
  vgenes <- 'V segments'
  jgenes <- 'J segments'
  dgenes <- 'D segments'
  vend <- 'Last V nucleotide position'
  jstart <- 'First J nucleotide position'
  dalignments <- c('First D nucleotide position', 'Last D nucleotide position')
  vd.insertions <- 'VD insertions'
  dj.insertions <- 'DJ insertions'
  total.insertions <- 'Total insertions'
  .skip = 0
  .sep = '\t'
  
  parse.cloneset(.filename = filename, 
                 .nuc.seq = nuc.seq,
                 .aa.seq = aa.seq,
                 .count = count,
                 .proportion = Proportion,
                 .reads = reads,
                 .events = events,
                 .vgenes = vgenes,
                 .jgenes = jgenes,
                 .dgenes = dgenes,
                 .vend = vend,
                 .jstart = jstart,
                 .dalignments = dalignments,
                 .vd.insertions = vd.insertions,
                 .dj.insertions = dj.insertions,
                 .total.insertions = total.insertions,
                 .skip = .skip,
                 .sep = .sep)
}

parse.migec <- function (.filename) {
  
  filename <- .filename
  nuc.seq <- 'CDR3 nucleotide sequence'
  aa.seq <- 'CDR3 amino acid sequence'
  count <- 'Count'
  Proportion <- NA
  reads <- 'Good reads'
  events <- 'Good events'
  vgenes <- 'V segments'
  jgenes <- 'J segments'
  dgenes <- 'D segments'
  vend <- 'Last V nucleotide position'
  jstart <- 'First J nucleotide position'
  dalignments <- c('First D nucleotide position', 'Last D nucleotide position')
  vd.insertions <- 'VD insertions'
  dj.insertions <- 'DJ insertions'
  total.insertions <- 'Total insertions'
  .skip = 0
  .sep = '\t'
  
  parse.cloneset(.filename = filename, 
                 .nuc.seq = nuc.seq,
                 .aa.seq = aa.seq,
                 .count = count,
                 .proportion = Proportion,
                 .reads = reads,
                 .events = events,
                 .vgenes = vgenes,
                 .jgenes = jgenes,
                 .dgenes = dgenes,
                 .vend = vend,
                 .jstart = jstart,
                 .dalignments = dalignments,
                 .vd.insertions = vd.insertions,
                 .dj.insertions = dj.insertions,
                 .total.insertions = total.insertions,
                 .skip = .skip,
                 .sep = .sep)
}