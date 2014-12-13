#' Tables with genetic code.
#' 
#' @docType data
#' 
#' @aliases AA_TABLE AA_TABLE_REVERSED
#' 
#' @description
#' Tables with genetic code.
#' 
#' @format
#' AA_TABLE:
#' 
#' \code{Class 'table'  Named chr [1:65] "K" "N" "K" "N" ...
#' ..- attr(*, "names")= chr [1:65] "AAA" "AAC" "AAG" "AAT" ...}
#' 
#' AA_TABLE_REVERSED:
#' 
#' \code{List of 22
#' $ *: chr [1:3] "TAA" "TAG" "TGA"
#' $ A: chr [1:4] "GCA" "GCC" "GCG" "GCT"
#' $ C: chr [1:2] "TGC" "TGT"
#' $ D: chr [1:2] "GAC" "GAT"
#' ...
#' }
#' 
#' @examples
#' \dontrun{
#' AA_TABLE['ATG']  #  => "M"
#' AA_TABLE_REVERSED['K']  #  => list(K = c("AAA", "AAG"))
#' }
AA_TABLE <- table(c('TCA', 'TCG', 'TCC', 'TCT', 'TTT', 'TTC', 'TTA', 'TTG', 'TAT', 'TAC', 'TAA', 'TAG', 'TGT',
                    'TGC', 'TGA', 'TGG', 'CTA', 'CTG', 'CTC', 'CTT', 'CCA', 'CCG', 'CCC', 'CCT', 'CAT', 'CAC',
                    'CAA', 'CAG', 'CGA', 'CGG', 'CGC', 'CGT', 'ATT', 'ATC', 'ATA', 'ATG', 'ACA', 'ACG', 'ACC',
                    'ACT', 'AAT', 'AAC', 'AAA', 'AAG', 'AGT', 'AGC', 'AGA', 'AGG', 'GTA', 'GTG', 'GTC', 'GTT',
                    'GCA', 'GCG', 'GCC', 'GCT', 'GAT', 'GAC', 'GAA', 'GAG', 'GGA', 'GGG', 'GGC', 'GGT', 'NNN'))
AA_TABLE[c('TCA', 'TCG', 'TCC', 'TCT', 'TTT', 'TTC', 'TTA', 'TTG', 'TAT', 'TAC', 'TAA', 'TAG', 'TGT',
           'TGC', 'TGA', 'TGG', 'CTA', 'CTG', 'CTC', 'CTT', 'CCA', 'CCG', 'CCC', 'CCT', 'CAT', 'CAC',
           'CAA', 'CAG', 'CGA', 'CGG', 'CGC', 'CGT', 'ATT', 'ATC', 'ATA', 'ATG', 'ACA', 'ACG', 'ACC',
           'ACT', 'AAT', 'AAC', 'AAA', 'AAG', 'AGT', 'AGC', 'AGA', 'AGG', 'GTA', 'GTG', 'GTC', 'GTT',
           'GCA', 'GCG', 'GCC', 'GCT', 'GAT', 'GAC', 'GAA', 'GAG', 'GGA', 'GGG', 'GGC', 'GGT', 'NNN')] <- c('S', 'S', 'S', 'S', 'F', 'F', 'L', 'L', 'Y', 'Y', '*', '*', 'C', 'C', '*',
                                                                                                            'W', 'L', 'L', 'L', 'L', 'P', 'P', 'P', 'P', 'H', 'H', 'Q', 'Q', 'R', 'R',
                                                                                                            'R', 'R', 'I', 'I', 'I', 'M', 'T', 'T', 'T', 'T', 'N', 'N', 'K', 'K', 'S',
                                                                                                            'S', 'R', 'R', 'V', 'V', 'V', 'V', 'A', 'A', 'A', 'A', 'D', 'D', 'E', 'E',
                                                                                                            'G', 'G', 'G', 'G', '~')
AA_TABLE_REVERSED <- sapply(unique(AA_TABLE), function (aa) { names(AA_TABLE)[AA_TABLE == aa] })
AA_TABLE_REVERSED <- AA_TABLE_REVERSED[order(names(AA_TABLE_REVERSED))]


#' Alphabets of V-J segments.
#' 
#' @docType data
#' 
#' @name segments.alphabets
#' 
#' @aliases V_BETA_ALPHABET J_BETA_ALPHABET V_ALPHA_ALPHABET J_ALPHA_ALPHABET
#' 
#' @description
#' Character vector with names for segments. With \code{tcR} we provided
#' two dataset of segments for humans ("./data/human.alphabets.rda") and mouses ("./data/mouse.alphabets.rda").
#' 
#' @format
#' Each \code{*_*_ALPHABET} is a character vector.
#' 
#' @examples
#' \dontrun{
#' data(human.alphabets)
#' V_BETA_ALPHABET[1]  #  => "TRBV10-1"
#' }
NULL


#' Segment data.
#' 
#' @docType data
#' 
#' @aliases segments
#' 
#' @name segments.list
#' 
#' @description
#' \code{segments} is a list with 5 data frames with data of human alpha-beta chain segments.
#' Elements names as "TRAV", "TRAJ", "TRBV", "TRVJ", "TRVD". Each data frame consists of 5 columns:
#' 
#' - V.allelles / J.allelles / D.allelles - character column with names of V/D/J-segments.
#' 
#' - CDR3.position - position in the full nucleotide segment sequence where CDR3 starts.
#' 
#' - Full.nucleotide.sequence - character column with segment CDR1-2-3 sequence.
#' 
#' - Nucleotide.sequence - character column with segment CDR3 sequences.
#' 
#' - Nucleotide.sequence.P - character column with segment CDR3 sequences with P-insertions.
#' 
#' @format
#' \code{segments} is a list with data frames.
#' 
#' @examples
#' \dontrun{
#' data(segments)
#' segments$Nucleotide.sequence[segments$TRBV[,1] == "TRBV10-1"]
#' }
NULL


#' List with assembling probabilities of beta chain TCRs.
#' 
#' @docType data
#' 
#' @name beta.prob
#' 
#' @aliases beta.prob
#' 
#' @description
#' \code{beta.prob} is a list with probabilities of TCR assembling taken from
#' \code{Murugan et al. Statistical inference of the generation probability
#' of T-cell receptors from sequence repertoires}. It's a list with the following elements:
#' 
#' - P.V - matrix with 1 column and row names stands for V-beta segments. Each element is
#' a probability of choosing corresponding V-beta segment.
#' 
#' - P.del.D1 - matrix 17x17 with probabilities of choosing D5-D3 deletions for TRBD1.
#' 
#' - P.del.D1 - matrix 17x17 with probabilities of choosing D5-D3 deletions for TRBD2.
#' 
#' - P.ins.len - matrix with first columns stands for number of insertions, second and third columns filled
#' with probability values of choosing corresponding number of insertions in VD- and DJ-junctions
#' correspondingly.
#' 
#' - P.ins.nucl - data frame with probability of choosing a nucleotide in the insertion on junctions with known
#' previous nucleotide. First column with names of nucleotides, 2-5 columns are probabilities of choosing
#' adjacent nucleotide in VD-junction, 6-9 columns are probabilities of choosing adjacent nucleotide in DJ-junction.
#' 
#' - P.del.J - matrix with the first column "P.del.V" with number of deletions, and other columns with
#' names for V-segments and with probabilities of corresponding deletions.
#' 
#' - P.del.J - matrix with the first column "P.del.J" with number of deletions, and other columns with
#' names for J-segments and with probabilities of corresponding deletions.
#' 
#' - P.J.D - matrix with two columns ("TRBD1" and "TRBD2") and 13 rows with row names stands for 
#' J-beta segments. Each element is a mutual probability of choosing J-D segments.
#' 
#' @format
#' \code{beta.prob} is a list of matrices and data frames.
#' 
#' @examples
#' \dontrun{
#' # Generate 10 kmers with adjacent nucleotide probability.
#' generate.kmers.prob(rep.int(10, 10), .probs=beta.prob$P.ins.nucl[,c(1, 2:5)])
#' }
NULL


#' Twins alpha-beta chain data
#' 
#' @docType data
#' 
#' @name twinsdata
#' 
#' @aliases twa twb
#' 
#' @description
#' \code{twa.rda}, \code{twb.rda} - data frames with downsampled to the 10000 most 
#' abundant clonesets and 4 samples data of twins data (alpha and beta chains).
#' Link: http://labcfg.ibch.ru/tcr.html
#' 
#' @format
#' \code{twa} and \code{twb} are lists of 4 data frames with 10000 row in each.
#' 
#' @examples
#' \dontrun{
#' data(twa)
#' data(twb)
#' }
NULL