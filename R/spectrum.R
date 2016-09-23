########## Spectratyping ##########


if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("Length", "Val"))
}


#' Spectratype
#' 
#' @description Plot a spectratype plot - a histogram of read counts / umi counts by CDR3 length.
#' 
#' @param .data tcR data frame.
#' @param .quant Either "read.count" or "umi.count" for choosing the corresponding columns, or "id" to compute avoid using counts.
#' @param .gene Either NA for not using genes, "V" or "J" for corresponding genes.
#' @param .plot If T than plot the spectratype plot, otherwise return a table with data for lengths and counts.
#' @param .main Main title.
#' @param .legend Legend title.
#' @param .labs Character vector of length 2 for x-lab and y-lab.
#' 
#' @examples
#' \dontrun{
#' data(twb)
#' tmp = twb[[1]]
#' spectratype(tmp)
#' spectratype(tmp, .quant = "id", .plot = T, .gene = 'V')
#' spectratype(tmp, .quant = "read.count", .plot = F)
#' }
spectratype <- function (.data, .quant = c("read.count", "umi.count", "id"), .gene = "V", .plot = T, 
                           .main = "Spectratype", .legend = "Gene segment", .labs = c("CDR3 length", NA)) {
  .data$Length = nchar(.data$CDR3.nucleotide.sequence)
  if (.quant[1] == "id") {
    .data$Count = 1
  } else {
    .data$Count = .data[[.column.choice(.quant[1])]]
  }
  
  if (is.na(.gene)) {
    df = summarise(group_by(do.call(select_, list(.data, "Length", "Count")), Length), Val = sum(Count))
    p = ggplot() + geom_bar(aes(x = Length, y = Val), data = df, stat = "identity")
  } else {
    .gene = switch(tolower(substr(.gene, 1, 1)), v = "V.gene", j = "J.gene")
    df = do.call(select_, list(.data, "Length", "Count", .gene))
    df = group_by_(df, "Length", .gene)
    df = summarise(df, Val = sum(Count))
    df$Gene = df[[.gene]]
    df[[.gene]] = NULL
    df = df[order(df$Val, decreasing = T),]
    
    dup = which(cumsum(!duplicated(df$Gene)) == 12)[1]
    if (length(dup)) {
      uniq = unique(df$Gene[1:(dup - 1)])
      df$Gene[!(df$Gene %in% uniq)] = "Z"  # <- dirty hack to avoid factors
    }
    
    p = ggplot() + geom_bar(aes(x = Length, y = Val, fill = Gene), data = df, stat = "identity") + 
      scale_fill_manual(name = .legend, breaks = c(sort(uniq), "Z"),
                        labels=c(sort(uniq), "Other"), 
                        values = c("#9E0142", "#D53E4F", "#F46D43", "#FDAE61", "#FEE08B", "#FFFFBF", "#E6F598", "#ABDDA4", "#66C2A5", "#3288BD", "#5E4FA2", "grey75")) + 
      ggtitle(.main)
  }
  
  if (.plot) {
    if (!is.na(.labs[2])) {
      p = p + ylab(.labs[2])
    } else {
      p = p + ylab(switch(.quant[1], id = "Number of clonotypes", 
                          read.count = "Sum of reads", 
                          umi.count = "Sum of UMIs"))
    }
    p = p + xlab(.labs[1])
    p + theme_linedraw()
  } else {
    as.data.frame(df)
  }
}