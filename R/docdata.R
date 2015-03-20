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


#' Alphabets of TCR and Ig gene segments.
#' 
#' @docType data
#' 
#' @name segments.alphabets
#' 
#' @aliases genealphabets HUMAN_TRAV_ALPHABET HUMAN_TRAJ_ALPHABET HUMAN_TRBV_ALPHABET HUMAN_TRBD_ALPHABET HUMAN_TRBJ_ALPHABET HUMAN_TRBV_ALPHABET_MITCR HUMAN_TRGV_ALPHABET HUMAN_TRGJ_ALPHABET HUMAN_TRDV_ALPHABET HUMAN_TRDD_ALPHABET HUMAN_TRDJ_ALPHABET
#' 
#' @usage
#' HUMAN_TRAV_ALPHABET
#' 
#' HUMAN_TRAJ_ALPHABET
#' 
#' HUMAN_TRBV_ALPHABET
#' 
#' HUMAN_TRBD_ALPHABET
#' 
#' HUMAN_TRBJ_ALPHABET
#' 
#' HUMAN_TRBV_ALPHABET_MITCR
#' 
#' HUMAN_TRGV_ALPHABET
#' 
#' HUMAN_TRGJ_ALPHABET
#' 
#' HUMAN_TRDV_ALPHABET
#' 
#' HUMAN_TRDD_ALPHABET
#' 
#' HUMAN_TRDJ_ALPHABET
#' 
#' @description
#' Character vector with names for segments. With \code{tcR} we provided alphabets for all alpha, beta,
#' gamma and delta chains gene segments.
#' 
#' @format
#' Each \code{<SPECIES>_<GENES>_ALPHABET} is a character vector. <SPECIES> is an identifier of species, <GENES> is concatenated three
#' identifiers of cell type ("TR**" for TCR, "IG**" for Ig), chain (e.g., "**A*" for alpha chains) and gene segment ("***V" for V(ariable) gene segment, 
#' "***J" for J(oining) gene segment, "***D" for D(iversity) gene segment).
#' 
#' @examples
#' \dontrun{
#' HUMAN_TRBV_ALPHABET[1]  #  => "TRBV10-1"
#' }
HUMAN_TRAV_ALPHABET <- c('TRAV1-1', 'TRAV1-2', 'TRAV10', 'TRAV11', 'TRAV12-1', 'TRAV12-2', 'TRAV12-3', 'TRAV13-1',
                         'TRAV13-2', 'TRAV14/DV4', 'TRAV16', 'TRAV17', 'TRAV18', 'TRAV19', 'TRAV2', 'TRAV20',
                         'TRAV21', 'TRAV22', 'TRAV23/DV6', 'TRAV24', 'TRAV25', 'TRAV26-1', 'TRAV26-2', 'TRAV27',
                         'TRAV29/DV5', 'TRAV3', 'TRAV30', 'TRAV34', 'TRAV35', 'TRAV36/DV7', 'TRAV38-1', 'TRAV38-2/DV8',
                         'TRAV39', 'TRAV4', 'TRAV40', 'TRAV41', 'TRAV5', 'TRAV6', 'TRAV7', 'TRAV8-1',
                         'TRAV8-2', 'TRAV8-3', 'TRAV8-4', 'TRAV8-6', 'TRAV8-7', 'TRAV9-1', 'TRAV9-2')
HUMAN_TRAJ_ALPHABET <- c('TRAJ10', 'TRAJ11', 'TRAJ12', 'TRAJ13', 'TRAJ14', 'TRAJ15', 'TRAJ16', 'TRAJ17',
                         'TRAJ18', 'TRAJ20', 'TRAJ21', 'TRAJ22', 'TRAJ23', 'TRAJ24', 'TRAJ26', 'TRAJ27',
                         'TRAJ28', 'TRAJ29', 'TRAJ3', 'TRAJ30', 'TRAJ31', 'TRAJ32', 'TRAJ33', 'TRAJ34',
                         'TRAJ36', 'TRAJ37', 'TRAJ38', 'TRAJ39', 'TRAJ4', 'TRAJ40', 'TRAJ41', 'TRAJ42',
                         'TRAJ43', 'TRAJ44', 'TRAJ45', 'TRAJ46', 'TRAJ47', 'TRAJ48', 'TRAJ49', 'TRAJ5',
                         'TRAJ50', 'TRAJ52', 'TRAJ53', 'TRAJ54', 'TRAJ56', 'TRAJ57', 'TRAJ6', 'TRAJ7',
                         'TRAJ8', 'TRAJ9')

HUMAN_TRBV_ALPHABET <- c('TRBV10-1', 'TRBV10-2', 'TRBV10-3', 'TRBV11-1', 'TRBV11-2', 'TRBV11-3', 'TRBV12-4', 'TRBV12-3',
                         'TRBV12-5', 'TRBV13', 'TRBV14', 'TRBV15', 'TRBV16', 'TRBV18', 'TRBV19', 'TRBV2',
                         'TRBV20-1', 'TRBV21-1', 'TRBV23-1', 'TRBV24-1', 'TRBV25-1', 'TRBV27', 'TRBV28', 'TRBV29-1',
                         'TRBV3-1', 'TRBV30', 'TRBV4-1', 'TRBV4-2', 'TRBV4-3', 'TRBV5-1', 'TRBV5-4', 'TRBV5-5',
                         'TRBV5-6', 'TRBV5-8', 'TRBV6-1', 'TRBV6-3', 'TRBV6-2', 'TRBV6-4', 'TRBV6-5', 'TRBV6-6',
                         'TRBV6-7', 'TRBV7-1', 'TRBV7-2', 'TRBV7-3', 'TRBV7-4', 'TRBV7-6', 'TRBV7-7', 'TRBV7-8',
                         'TRBV7-9', 'TRBV9')
HUMAN_TRBD_ALPHABET <- c('TRBD1', 'TRBD2')
HUMAN_TRBJ_ALPHABET <- c('TRBJ1-1', 'TRBJ1-2', 'TRBJ1-3', 'TRBJ1-4', 'TRBJ1-5', 'TRBJ1-6', 'TRBJ2-1', 'TRBJ2-2',
                         'TRBJ2-3', 'TRBJ2-4', 'TRBJ2-5', 'TRBJ2-6', 'TRBJ2-7')

HUMAN_TRBV_ALPHABET_MITCR <- c('TRBV10-1', 'TRBV10-2', 'TRBV10-3', 'TRBV11-1', 'TRBV11-2', 'TRBV11-3', 'TRBV12-4, TRBV12-3',
                               'TRBV12-5', 'TRBV13', 'TRBV14', 'TRBV15', 'TRBV16', 'TRBV18', 'TRBV19', 'TRBV2',
                               'TRBV20-1', 'TRBV21-1', 'TRBV23-1', 'TRBV24-1', 'TRBV25-1', 'TRBV27', 'TRBV28', 'TRBV29-1',
                               'TRBV3-1', 'TRBV30', 'TRBV4-1', 'TRBV4-2', 'TRBV4-3', 'TRBV5-1', 'TRBV5-4', 'TRBV5-5',
                               'TRBV5-6', 'TRBV5-8', 'TRBV6-1', 'TRBV6-3, TRBV6-2', 'TRBV6-4', 'TRBV6-5', 'TRBV6-6', 'TRBV6-7',
                               'TRBV7-1', 'TRBV7-2', 'TRBV7-3', 'TRBV7-4', 'TRBV7-6', 'TRBV7-7', 'TRBV7-8', 'TRBV7-9', 'TRBV9')

HUMAN_TRGV_ALPHABET <- c()
HUMAN_TRGJ_ALPHABET <- c()

HUMAN_TRDV_ALPHABET <- c()
HUMAN_TRDD_ALPHABET <- c()
HUMAN_TRDJ_ALPHABET <- c()


#' Segment data.
#' 
#' @docType data
#' 
#' @aliases genesegments
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
#' \code{genesegments} is a list with data frames.
#' 
#' @examples
#' \dontrun{
#' data(genesegments)
#' genesegments$Nucleotide.sequence[segments$TRBV[,1] == "TRBV10-1"]
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