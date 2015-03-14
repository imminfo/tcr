########## Data processing functions ##########


#' Parse input table files with immune receptor repertoire data.
#'
#' @description
#' Load the MITCR TCR data from the file with the given filename
#' to a data frame.
#'
#' @param .filepath Path to the input file with TCR data.
#' @param .barcode.flag If T than load MiTCR data.frames with new column Barcode.count
#' and without column Min.quality.
#' @param .i,.all Don't supply it, they are just for the output.
#' @return Data.frame with TCR data.
#' 
#' @seealso \link{parse.file.list}, \link{parse.folder}
#' 
#' @examples
#' \dontrun{
#' # Parse file in "~/data/immdata1.txt".
#' immdata1 <- parse.file("~/data/immdata1.txt")
#' }
NULL


#' Parse files from the given vector or list with filepaths
#' and return list with data.frames.
#' 
#' @description
#' Given the vector or list with filepaths, parse each file to a data frame
#' and combine them in a list. List items have names similar to names in the
#' given vector if filenames.
#' 
#' @param .filenames Vector or list with paths to files with TCR data.
#' @param .barcode.flag If T than load MiTCR data.frames with new column Barcode.count
#' and without column Min.quality.
#' @return List with data frame for the each name in the given filepaths.
#' 
#' @seealso \link{parse.file}, \link{parse.folder}
#' 
#' @examples
#' \dontrun{
#' # Parse files "~/data/immdata1.txt" and "~/data/immdat2.txt".
#' immdata12 <- parse.file.list(c("~/data/immdata1.txt", "~/data/immdata2.txt"))
#' }
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


#' Parse all files to dataframes from the given path to folder.
#' 
#' @description
#' Given the path to a folder with files with TCR data, parse all of them
#' to a list and return it. List items have names similar to names in the
#' given vector if filenames.
#' 
#' @param .folderpath Path to the folder with files with TCR data.
#' @param .barcode.flag If T than load MiTCR data.frames with new column Barcode.count
#' and without column Min.quality.
#' @return List with data frame for the each file in the given folder.
#' 
#' @seealso \link{parse.file}, \link{parse.file.list}
#' 
#' @examples
#' \dontrun{
#' # Parse all files in "~/data/".
#' immdata <- parse.folder("~/data/")
#' }
parse.folder <- function (.folderpath, .format = c('mitcr', 'mitcrbc', 'migec')) {
  parse.file.list(list.files(.folderpath, full.names = T), .format)
}


parse.file <- function(.filename, .format = c('mitcr', 'mitcrbc', 'migec')) {
  
  parse.fun <- switch(.format[1], 
                      mitcr = parse.mitcr,
                      mitcrbc = parse.mitcrbc,
                      migec = parse.migec)
  
  parse.fun(.filename)
}

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
#   print(rbind(table.colnames, col.classes))
  
  suppressWarnings(df <- read.table(file = .filename, header = T, colClasses = col.classes, sep = .sep, skip = .skip, strip.white = T))

  if(is.na(.events)) {
    .events <- "Events"
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
                    'Reads', 'Events')
  
  df
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