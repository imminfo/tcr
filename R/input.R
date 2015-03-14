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


parse.file <- (.filename, .format = c('mitcr', 'mitcrbc', 'migec')) {
  
  parse.fun <- switch(.format[1], 
                      mitcr = parse.file.mitcr,
                      mitcrbc = parse.file.mitcrbc,
                      migec = parse.file.migec)
  
  parse.fun(.filename)
}

parse.file.table <- (.filename,
                     .nuc.seq,
                     .aa.seq,
                     .count,
                     .percentage,
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
                     .skip = 0) {
  
  .replace.spaces <- function (.s) {
    gsub(' ', '.', .s, fixed = T, useBytes = T)
  }
  
  .nuc.seq <- .replace.spaces(.nuc.seq)
  .aa.seq <- .replace.spaces(.aa.seq)
  # ???
  # ???
  # ???
  if (is.na(.aa.seq)) {
    df$CDR3.amino.acid.sequence <- bunch.translate(df$CDR3.nucleotide.sequence)
  }
  if (is.na(.percentage)) {
    df$Percentage <- df$Count / sum(df$Count)
  }
  
  
  
  df <- df[c(.count, .percentage, .nuc.seq, .aa.seq,
             .vgenes, .jgenes, .dgenes,
             .vend, .jstart, .dalignments,
             .vd.insertions, .dj.insertions, .total.insertions,
             .reads, .events), ]
  
  colnames(.df) <-  c('Count', 'Percentage', 'CDR3.nucleotide.sequence', 'CDR3.amino.acid.sequence',
                      'V.segments', 'J.segments', 'D.segments',
                      'V.end', 'J.start', 'D5.end', 'D3.end',
                      'VD.insertions', 'DJ.insertions', 'Total.insertions',
                      'Reads', 'Events')
  
  df
}

parse.file.mitcr <- function (.filename) {
  
  nuc.seq <- 'CDR3 nucleotide sequence'
  aa.seq <- 'CDR3 amino acid sequence'
  count <- 'Read count'
  percentage <- 'Percentage'
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
    
  parse.file.table(.filename = filename, 
                   .nuc.seq = nuc.seq,
                   .aa.seq = aa.seq,
                   .count = count,
                   .percentage = percentage,
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
                   .skip = .skip)
}

parse.file.mitcrbc <- function (.filename) {
  
  nuc.seq <- 'CDR3 nucleotide sequence'
  aa.seq <- 'CDR3 amino acid sequence'
  count <- 'Barcode.count'
  percentage <- NA
  reads <- 'Read count'
  events <- 'Barcode.count'
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
  
  parse.file.table(.filename = filename, 
                   .nuc.seq = nuc.seq,
                   .aa.seq = aa.seq,
                   .count = count,
                   .percentage = percentage,
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
                   .skip = .skip)
}

parse.file.migec <- function (.filename) {
  
  nuc.seq <- 'CDR3 nucleotide sequence'
  aa.seq <- 'CDR3 amino acid sequence'
  count <- 'Count'
  percentage <- NA
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
  
  parse.file.table(.filename = filename, 
                   .nuc.seq = nuc.seq,
                   .aa.seq = aa.seq,
                   .count = count,
                   .percentage = percentage,
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
                   .skip = .skip)
}