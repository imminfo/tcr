# general function for loading data frames - either from text files
# or from MonetDB
repLoad <- function (.path, .format = c("monetdb", "mitcr", "migec"), .user = "default", .password = "default") {
  res <- list()
  
  for (i in 1:length(.path)) {
    if (dir.exists(.path[i])) {
      if (.format[1] == "monetdb") {
        # ???
        # ???
        # ???
      } else {
        res <- c(res, parse.folder(.path[i], .format[1]))
      }
    } else if (file.exists(.path[i])) {
      if (.format[1] == "monetdb") {
        conn <- dbConnect(MonetDB.R(), host="localhost", dbname="demo", user="monetdb", password="monetdb")
        monetdb_conn <- src_monetdb("demo")
        res <- c(res, tbl(monetdb_conn, "some_table"))
      } else {
        res <- c(res, parse.file(.path[1], .format[1]))
      }
    } else {
      cat('Can\'t find folder or file:\t"', .path[i], '"', sep = '', end = '\n')
    }
  }
  
  res
}


#' Save tcR data frames to disk as text files or as MonetDB database.
#' 
#' @description
#' Save repertoire files to either text files or gzipped text files or as MonetDB database.
#' You can read them later by \code{repLoad} function with either \code{.format = "tcr"} or with
#' \code{.format = "monetdb"}.
#' 
#' @param .data Either tcR data frame or a list of tcR data frames.
#' @param .format
#' @param .names Names of output files. By default it's an empty string so names will be taken from names of the input list.
#' @param .folder Path to the folder with output files.
#' 
#' @seealso \link{repLoad}
repSave <- function (.data, .format = c("txt", "gz", "monetdb"), .names = "", .folder = "./", .user = "default", .password = "default") {
  if (has.class(.data, 'data.frame')) { .data <- list(Sample = .data) }
  
  .folder <- paste0(.folder, "/")
  
  postfix <- ".txt"
  filefun <- function (...) file(...)
  if (.compress) { 
    postfix <- ".txt.gz"
    filefun <- function (...) gzfile(...)
  }
  
  if (.names[1] == "") {
    .names = paste0(.folder, names(.data), postfix)
  } else {
    if (length(.data) != length(.names)) {
      cat("Number of input data frames isn't equal to number of names\n")
      return(NULL)
    } else {
      .names = paste0(.folder, .names, postfix)
    }
  }
  
  for (i in 1:length(.data)) {
    cat("Writing", .names[i], "file...\t")
    fc <- filefun(description = .names[i], open = "w")
    write.table(.data[[i]], fc, quote = F, row.names = F, sep = '\t')
    close(fc)
    cat("Done.\n")
  }
}