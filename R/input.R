########## Data processing functions ##########


.onLoad <- function (libname, pkgname) {
  human.alphabets <- NULL
  utils::data(human.alphabets, package = pkgname)
}



#' Parse input file with the given filename to a data frame.
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
parse.file <- function (.filepath, .barcode.flag = F, .i = 1, .all = 1) {
  # Level without header and some columns
  LEVEL_1_NAMES <- c('Read.count', 'Percentage', 'CDR3.nucleotide.sequence',
                     'CDR3.amino.acid.sequence', 'V.segments', 'J.segments', 'D.segments')
  LEVEL_1_CLASSES <- c('integer', 'numeric', 'character',
                       'character', 'character', 'character', 'character')
  
  # Level without header and some columns
  LEVEL_2_NAMES <- c('Read.count', 'Percentage', 'CDR3.nucleotide.sequence',
                     'CDR3.amino.acid.sequence', 'V.segments', 'J.segments', 'D.segments',
                     'Last.V.nucleotide.position', 'First.D.nucleotide.position',
                     'Last.D.nucleotide.position', 'First.J.nucleotide.position',
                     'VD.insertions', 'DJ.insertions', 'Total.insertions')
  LEVEL_2_CLASSES <- c('integer', 'numeric', 'character',
                       'character', 'character', 'character', 'character',
                       'integer', 'integer',
                       'integer', 'integer',
                       'integer', 'integer', 'integer')
  
  # Level with all columns.
  if (.barcode.flag) {
    LEVEL_3_NAMES <- c('Barcode.count', 'Read.count', 'Percentage', 'CDR3.nucleotide.sequence',
                       'CDR3.nucleotide.quality',
                       'CDR3.amino.acid.sequence', 'V.alleles', 'V.segments', 
                       'J.alleles', 'J.segments', 'D.alleles', 'D.segments',
                       'Last.V.nucleotide.position', 'First.D.nucleotide.position',
                       'Last.D.nucleotide.position', 'First.J.nucleotide.position',
                       'VD.insertions', 'DJ.insertions', 'Total.insertions')
    LEVEL_3_CLASSES <- c('integer', 'integer', 'numeric', 'character',
                         'character', 
                         'character', 'character', 'character', 
                         'character', 'character', 'character', 'character',
                         'integer', 'integer',
                         'integer', 'integer',
                         'integer', 'integer', 'integer')
  } else {
    LEVEL_3_NAMES <- c('Read.count', 'Percentage', 'CDR3.nucleotide.sequence',
                       'CDR3.nucleotide.quality', 'Min.quality', 
                       'CDR3.amino.acid.sequence', 'V.alleles', 'V.segments', 
                       'J.alleles', 'J.segments', 'D.alleles', 'D.segments',
                       'Last.V.nucleotide.position', 'First.D.nucleotide.position',
                       'Last.D.nucleotide.position', 'First.J.nucleotide.position',
                       'VD.insertions', 'DJ.insertions', 'Total.insertions')
    LEVEL_3_CLASSES <- c('integer', 'numeric', 'character',
                         'character', 'integer', 
                         'character', 'character', 'character', 
                         'character', 'character', 'character', 'character',
                         'integer', 'integer',
                         'integer', 'integer',
                         'integer', 'integer', 'integer')
  }
  
  # Vector with all levels
  LEVEL_NAME_LIST <- list(LEVEL_1_NAMES, LEVEL_2_NAMES, LEVEL_3_NAMES)
  LEVEL_CLASS_LIST <- list(LEVEL_1_CLASSES, LEVEL_2_CLASSES, LEVEL_3_CLASSES)
  # ===========================================
  
  
  # Check if the given string is a MiTCR header.
  .header.check <- function (.str) {
    length(unlist(strsplit(.str, '\t'))) < length(LEVEL_1_NAMES)
  }
  
  cat('Parsing file', .i, '/', .all, .filepath, '...\t')
  
  input.file <- file(.filepath, 'r')
  
  first.line <- readLines(input.file, n = 1)[1]
  header.flag <- F
  if (.header.check(first.line)) {
    first.line <- readLines(input.file, n = 1)[1]
    header.flag <- T
  }
  
  level <- -1
  words <- unlist(strsplit(first.line, '\t'))
  if (length(words) == length(LEVEL_1_NAMES) || (length(words) - 1) == length(LEVEL_1_NAMES)) {
    level <- 1
  } else {
    if (length(words) == length(LEVEL_2_NAMES) || (length(words) - 1) == length(LEVEL_2_NAMES)) {
      level <- 2
    } else {
      level <- 3
    }
  }
  close(input.file)
  
  input.data <- NA
  if (header.flag) {
    input.data <- read.csv(file = .filepath,
                           header = T,
                           sep = '\t',
                           skip = 1,
                           colClasses = LEVEL_CLASS_LIST[[level]])
  } else {
    input.data <- read.csv(file = .filepath,
                           header = T,
                           sep = '\t',
                           colClasses = LEVEL_CLASS_LIST[[level]])
  }
  names(input.data) <- LEVEL_NAME_LIST[[level]]
  input.data$Rank <- rank(1 / input.data$Read.count, ties.method = 'average')
  
  cat('Done. Data.frame with', nrow(input.data), 'rows.\n')
  
  input.data
}


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
parse.file.list <- function (.filenames, .barcode.flag = F) {
  # Remove full paths and extension from the given string.
  .remove.ext <- function (.str) {
#     gsub(pattern = '.*/|[.].*$', replacement = '', x = .str)
    .str
  }
  
  .filenames <- as.list(.filenames)
  
  datalist <- lapply(X = 1:length(.filenames), FUN = function (i) parse.file(.filenames[[i]], .barcode.flag, i, length(.filenames)) )
  namelist <- lapply(X = .filenames, FUN = .remove.ext)
  names(datalist) <- unlist(namelist)
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
parse.folder <- function (.folderpath, .barcode.flag = F) {
  flist <- list.files(.folderpath, full.names = T)
  
  parse.file.list(flist, .barcode.flag)
}