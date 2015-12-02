.set.attr <- function (.data, .attr, .value) {
  attr(.data, .attr) <- .value
  .data
}


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
#' @aliases genealphabets HUMAN_TRAV HUMAN_TRAJ HUMAN_TRBV HUMAN_TRBD HUMAN_TRBJ HUMAN_TRBV_MITCR HUMAN_TRGV HUMAN_TRGJ HUMAN_TRDV HUMAN_TRDD HUMAN_TRDJ HUMAN_IGHV HUMAN_IGHD HUMAN_IGHJ HUMAN_IGKV HUMAN_IGKJ HUMAN_IGLJ HUMAN_IGLV HUMAN_TRBV_ALS HUMAN_TRBV_FAM HUMAN_TRBV_GEN MACMUL_TRBJ MACMUL_TRBV MOUSE_TRBJ MOUSE_TRBV
#' 
#' @usage
#' HUMAN_TRAV
#' 
#' HUMAN_TRAJ
#' 
#' HUMAN_TRBV
#' 
#' HUMAN_TRBD
#' 
#' HUMAN_TRBJ
#' 
#' HUMAN_TRBV_MITCR
#' 
#' HUMAN_TRBV_ALS
#' 
#' HUMAN_TRGV
#' 
#' HUMAN_TRGJ
#' 
#' HUMAN_TRDV
#' 
#' HUMAN_TRDD
#' 
#' HUMAN_TRDJ
#' 
#' MOUSE_TRBV
#' 
#' MOUSE_TRBJ
#' 
#' MACMUL_TRBV
#' 
#' MACMUL_TRBJ
#' 
#' HUMAN_IGHV
#' 
#' HUMAN_IGHD
#' 
#' HUMAN_IGHJ
#' 
#' HUMAN_IGLV
#' 
#' HUMAN_IGLJ
#' 
#' @description
#' Character vector with names for segments. With \code{tcR} we provided alphabets for all alpha, beta,
#' gamma and delta chains gene segments.
#' 
#' @format
#' Each \code{<SPECIES>_<GENES>} is a character vector. <SPECIES> is an identifier of species, <GENES> is concatenated three
#' identifiers of cell type ("TR**" for TCR, "IG**" for Ig), chain (e.g., "**A*" for alpha chains) and gene segment ("***V" for V(ariable) gene segment, 
#' "***J" for J(oining) gene segment, "***D" for D(iversity) gene segment).
#' 
#' @examples
#' \dontrun{
#' HUMAN_TRBV[1]  #  => "TRBV10-1"
#' }
HUMAN_TRAV <- c('TRAV1-1', 'TRAV1-2', 'TRAV10', 'TRAV11', 'TRAV12-1', 'TRAV12-2', 'TRAV12-3', 'TRAV13-1',
                'TRAV13-2', 'TRAV14/DV4', 'TRAV16', 'TRAV17', 'TRAV18', 'TRAV19', 'TRAV2', 'TRAV20',
                'TRAV21', 'TRAV22', 'TRAV23/DV6', 'TRAV24', 'TRAV25', 'TRAV26-1', 'TRAV26-2', 'TRAV27',
                'TRAV29/DV5', 'TRAV3', 'TRAV30', 'TRAV34', 'TRAV35', 'TRAV36/DV7', 'TRAV38-1', 'TRAV38-2/DV8',
                'TRAV39', 'TRAV4', 'TRAV40', 'TRAV41', 'TRAV5', 'TRAV6', 'TRAV7', 'TRAV8-1',
                'TRAV8-2', 'TRAV8-3', 'TRAV8-4', 'TRAV8-6', 'TRAV8-7', 'TRAV9-1', 'TRAV9-2')
HUMAN_TRAJ <- c('TRAJ10', 'TRAJ11', 'TRAJ12', 'TRAJ13', 'TRAJ14', 'TRAJ15', 'TRAJ16', 'TRAJ17',
                'TRAJ18', 'TRAJ20', 'TRAJ21', 'TRAJ22', 'TRAJ23', 'TRAJ24', 'TRAJ26', 'TRAJ27',
                'TRAJ28', 'TRAJ29', 'TRAJ3', 'TRAJ30', 'TRAJ31', 'TRAJ32', 'TRAJ33', 'TRAJ34',
                'TRAJ36', 'TRAJ37', 'TRAJ38', 'TRAJ39', 'TRAJ4', 'TRAJ40', 'TRAJ41', 'TRAJ42',
                'TRAJ43', 'TRAJ44', 'TRAJ45', 'TRAJ46', 'TRAJ47', 'TRAJ48', 'TRAJ49', 'TRAJ5',
                'TRAJ50', 'TRAJ52', 'TRAJ53', 'TRAJ54', 'TRAJ56', 'TRAJ57', 'TRAJ6', 'TRAJ7',
                'TRAJ8', 'TRAJ9')

HUMAN_TRBV <- c('TRBV10-1', 'TRBV10-2', 'TRBV10-3', 'TRBV11-1', 'TRBV11-2', 'TRBV11-3', 'TRBV12-4', 'TRBV12-3',
                'TRBV12-5', 'TRBV13', 'TRBV14', 'TRBV15', 'TRBV16', 'TRBV18', 'TRBV19', 'TRBV2',
                'TRBV20-1', 'TRBV21-1', 'TRBV23-1', 'TRBV24-1', 'TRBV25-1', 'TRBV27', 'TRBV28', 'TRBV29-1',
                'TRBV3-1', 'TRBV30', 'TRBV4-1', 'TRBV4-2', 'TRBV4-3', 'TRBV5-1', 'TRBV5-4', 'TRBV5-5',
                'TRBV5-6', 'TRBV5-8', 'TRBV6-1', 'TRBV6-3', 'TRBV6-2', 'TRBV6-4', 'TRBV6-5', 'TRBV6-6',
                'TRBV6-7', 'TRBV7-1', 'TRBV7-2', 'TRBV7-3', 'TRBV7-4', 'TRBV7-6', 'TRBV7-7', 'TRBV7-8',
                'TRBV7-9', 'TRBV9')
HUMAN_TRBV <- .set.attr(HUMAN_TRBV, "column", "V.gene")
HUMAN_TRBD <- c('TRBD1', 'TRBD2')
HUMAN_TRBJ <- c('TRBJ1-1', 'TRBJ1-2', 'TRBJ1-3', 'TRBJ1-4', 'TRBJ1-5', 'TRBJ1-6', 'TRBJ2-1', 'TRBJ2-2',
                'TRBJ2-3', 'TRBJ2-4', 'TRBJ2-5', 'TRBJ2-6', 'TRBJ2-7')

HUMAN_TRBV_MITCR <- c('TRBV10-1', 'TRBV10-2', 'TRBV10-3', 'TRBV11-1', 'TRBV11-2', 'TRBV11-3', 'TRBV12-4, TRBV12-3',
                      'TRBV12-5', 'TRBV13', 'TRBV14', 'TRBV15', 'TRBV16', 'TRBV18', 'TRBV19', 'TRBV2',
                      'TRBV20-1', 'TRBV21-1', 'TRBV23-1', 'TRBV24-1', 'TRBV25-1', 'TRBV27', 'TRBV28', 'TRBV29-1',
                      'TRBV3-1', 'TRBV30', 'TRBV4-1', 'TRBV4-2', 'TRBV4-3', 'TRBV5-1', 'TRBV5-4', 'TRBV5-5',
                      'TRBV5-6', 'TRBV5-8', 'TRBV6-1', 'TRBV6-3, TRBV6-2', 'TRBV6-4', 'TRBV6-5', 'TRBV6-6', 'TRBV6-7',
                      'TRBV7-1', 'TRBV7-2', 'TRBV7-3', 'TRBV7-4', 'TRBV7-6', 'TRBV7-7', 'TRBV7-8', 'TRBV7-9', 'TRBV9')

HUMAN_TRBV_ALS <- c('TRBV10-1', 'TRBV10-2', 'TRBV10-3', 'TRBV11-1', 'TRBV11-2', 'TRBV11-3', 'TRBV12-4', 'TRBV12-3',
                    'TRBV12-5', 'TRBV13', 'TRBV14', 'TRBV15', 'TRBV16', 'TRBV18', 'TRBV19', 'TRBV2*01',
                    'TRBV20-1', 'TRBV21-1', 'TRBV23-1', 'TRBV24-1', 'TRBV25-1', 'TRBV27', 'TRBV28', 'TRBV29-1',
                    'TRBV3-1*01', 'TRBV30', 'TRBV4-1', 'TRBV4-2', 'TRBV4-3', 'TRBV5-1', 'TRBV5-4', 'TRBV5-5',
                    'TRBV5-6', 'TRBV5-8', 'TRBV6-1', 'TRBV6-3', 'TRBV6-2', 'TRBV6-4', 'TRBV6-5', 'TRBV6-6',
                    'TRBV6-7', 'TRBV7-1', 'TRBV7-2', 'TRBV7-3', 'TRBV7-4', 'TRBV7-6', 'TRBV7-7', 'TRBV7-8',
                    'TRBV7-9', 'TRBV9')
HUMAN_TRBV_ALS <- .set.attr(HUMAN_TRBV_ALS, "column", "V.allele")

HUMAN_TRBV_FAM <- c()
HUMAN_TRBV_GEN <- c()
HUMAN_TRBV_ALS <- c()


HUMAN_TRGV <- c()
HUMAN_TRGJ <- c()

HUMAN_TRDV <- c()
HUMAN_TRDD <- c()
HUMAN_TRDJ <- c()

HUMAN_IGHV <- c('IGHV1-18', 'IGHV1-2', 'IGHV1-24', 'IGHV1-3', 'IGHV1-45', 'IGHV1-46', 'IGHV1-58', 'IGHV1-69', 
                'IGHV1-69-2', 'IGHV1-69D', 'IGHV1-8', 'IGHV1-NL1', 'IGHV1/OR15-1', 'IGHV1/OR15-2', 'IGHV1/OR15-3', 
                'IGHV1/OR15-4', 'IGHV1/OR15-5', 'IGHV1/OR15-9', 'IGHV1/OR21-1', 'IGHV2-10', 'IGHV2-26', 'IGHV2-5', 
                'IGHV2-70', 'IGHV2/OR16-5', 'IGHV3-11', 'IGHV3-13', 'IGHV3-15', 'IGHV3-16', 'IGHV3-19', 'IGHV3-20',
                'IGHV3-21', 'IGHV3-22', 'IGHV3-23', 'IGHV3-23D', 'IGHV3-25', 'IGHV3-30', 'IGHV3-30-3', 'IGHV3-33', 
                'IGHV3-35', 'IGHV3-38', 'IGHV3-38-3', 'IGHV3-43', 'IGHV3-43D', 'IGHV3-47', 'IGHV3-48', 'IGHV3-49', 
                'IGHV3-52', 'IGHV3-53', 'IGHV3-54', 'IGHV3-62', 'IGHV3-63', 'IGHV3-64', 'IGHV3-66', 'IGHV3-69-1', 
                'IGHV3-7', 'IGHV3-71', 'IGHV3-72', 'IGHV3-73', 'IGHV3-74', 'IGHV3-9', 'IGHV3-NL1', 'IGHV3/OR15-7', 
                'IGHV3/OR16-10', 'IGHV3/OR16-12', 'IGHV3/OR16-13', 'IGHV3/OR16-14', 'IGHV3/OR16-15', 'IGHV3/OR16-16', 
                'IGHV3/OR16-8', 'IGHV3/OR16-9', 'IGHV4-28', 'IGHV4-30-2', 'IGHV4-30-4', 'IGHV4-31', 'IGHV4-34', 'IGHV4-38-2', 
                'IGHV4-39', 'IGHV4-4', 'IGHV4-55', 'IGHV4-59', 'IGHV4-61', 'IGHV4/OR15-8', 'IGHV5-10-1', 'IGHV5-51', 
                'IGHV5-78', 'IGHV6-1', 'IGHV7-34-1', 'IGHV7-4-1', 'IGHV7-81')
HUMAN_IGHD <- c('IGHD1-1', 'IGHD1-14', 'IGHD1-20', 'IGHD1-26', 'IGHD1-7', 'IGHD1/OR15-1a', 'IGHD1/OR15-1b', 'IGHD2-15', 
                'IGHD2-2', 'IGHD2-21', 'IGHD2-8', 'IGHD2/OR15-2a', 'IGHD2/OR15-2b', 'IGHD3-10', 'IGHD3-16', 'IGHD3-22', 
                'IGHD3-3', 'IGHD3-9', 'IGHD3/OR15-3a', 'IGHD3/OR15-3b', 'IGHD4-11', 'IGHD4-17', 'IGHD4-23', 'IGHD4-4', 
                'IGHD4/OR15-4a', 'IGHD4/OR15-4b', 'IGHD5-12', 'IGHD5-18', 'IGHD5-24', 'IGHD5-5', 'IGHD5/OR15-5a', 
                'IGHD5/OR15-5b', 'IGHD6-13', 'IGHD6-19', 'IGHD6-25', 'IGHD6-6', 'IGHD7-27')
HUMAN_IGHJ <- c('IGHJ1', 'IGHJ2', 'IGHJ3', 'IGHJ4', 'IGHJ5', 'IGHJ6')

HUMAN_IGKV <- c('IGKV1-12', 'IGKV1-13', 'IGKV1-16', 'IGKV1-17', 'IGKV1-27', 'IGKV1-33', 'IGKV1-39', 'IGKV1-5', 'IGKV1-6', 
                'IGKV1-8', 'IGKV1-9', 'IGKV1-NL1', 'IGKV1/OR-2', 'IGKV1/OR-3', 'IGKV1/OR-4', 'IGKV1/OR1-1', 'IGKV1/OR10-1', 
                'IGKV1/OR15-118', 'IGKV1/OR2-0', 'IGKV1/OR2-1', 'IGKV1/OR2-108', 'IGKV1/OR2-11', 'IGKV1/OR2-118', 'IGKV1/OR2-2',
                'IGKV1/OR2-3', 'IGKV1/OR2-9', 'IGKV1/OR22-5', 'IGKV1/OR9-1', 'IGKV1/OR9-2', 'IGKV1/ORY-1', 'IGKV1D-12', 'IGKV1D-13', 
                'IGKV1D-16', 'IGKV1D-17', 'IGKV1D-33', 'IGKV1D-39', 'IGKV1D-42', 'IGKV1D-43', 'IGKV1D-8', 'IGKV2-18', 'IGKV2-24', 'IGKV2-28', 
                'IGKV2-30', 'IGKV2-4', 'IGKV2-40', 'IGKV2/OR2-7D', 'IGKV2/OR22-4', 'IGKV2D-18', 'IGKV2D-24', 'IGKV2D-26', 'IGKV2D-28', 
                'IGKV2D-29', 'IGKV2D-30', 'IGKV2D-40', 'IGKV3-11', 'IGKV3-15', 'IGKV3-20', 'IGKV3-7', 'IGKV3-NL1', 'IGKV3-NL2', 'IGKV3-NL3', 
                'IGKV3-NL4', 'IGKV3-NL5', 'IGKV3/OR2-268', 'IGKV3D-11', 'IGKV3D-15', 'IGKV3D-20', 'IGKV3D-7', 'IGKV4-1', 'IGKV5-2', 'IGKV6-21', 
                'IGKV6D-21', 'IGKV6D-41', 'IGKV7-3')
HUMAN_IGKJ <- c('IGKJ1', 'IGKJ2', 'IGKJ3', 'IGKJ4', 'IGKJ5')

HUMAN_IGLV <- c('IGLV1-36', 'IGLV1-40', 'IGLV1-41', 'IGLV1-44', 'IGLV1-47', 'IGLV1-50', 'IGLV1-51', 'IGLV10-54', 'IGLV11-55', 'IGLV2-11', 'IGLV2-14', 
                'IGLV2-18', 'IGLV2-23', 'IGLV2-33', 'IGLV2-34', 'IGLV2-5', 'IGLV2-8', 'IGLV2-NL1', 'IGLV3-1', 'IGLV3-10', 'IGLV3-12', 'IGLV3-13', 'IGLV3-16', 
                'IGLV3-19', 'IGLV3-21', 'IGLV3-22', 'IGLV3-25', 'IGLV3-27', 'IGLV3-31', 'IGLV3-9', 'IGLV4-3', 'IGLV4-60', 'IGLV4-69', 'IGLV5-37', 'IGLV5-39', 
                'IGLV5-45', 'IGLV5-48', 'IGLV5-52', 'IGLV6-57', 'IGLV7-43', 'IGLV7-46', 'IGLV8-61', 'IGLV9-49')
HUMAN_IGLJ <- c('IGLJ1', 'IGLJ2', 'IGLJ3', 'IGLJ4', 'IGLJ5', 'IGLJ6', 'IGLJ7')

MOUSE_TRBV <- c('TRBV1', 'TRBV12-1', 'TRBV12-2', 'TRBV13-1', 'TRBV13-2', 'TRBV13-3', 
                'TRBV14', 'TRBV15', 'TRBV16', 'TRBV17', 'TRBV19', 'TRBV2', 'TRBV20',
                'TRBV23', 'TRBV24', 'TRBV26', 'TRBV29', 'TRBV3', 'TRBV30', 'TRBV31', 
                'TRBV4', 'TRBV5')
MOUSE_TRBJ <- c('TRBJ1-1', 'TRBJ1-2', 'TRBJ1-3', 'TRBJ1-4', 'TRBJ1-5', 'TRBJ2-1', 'TRBJ2-2',
                'TRBJ2-3', 'TRBJ2-4', 'TRBJ2-5', 'TRBJ2-7')

MACMUL_TRBV <- c('TRBD1', 'TRBD2', 'TRBJ1-1', 'TRBJ1-2', 'TRBJ1-3', 'TRBJ1-4', 'TRBJ1-5', 
                 'TRBJ1-6', 'TRBJ2-1', 'TRBJ2-2', 'TRBJ2-2P', 'TRBJ2-3', 'TRBJ2-4', 'TRBJ2-5',
                 'TRBJ2-6', 'TRBJ2-7', 'TRBV1-1', 'TRBV10-1', 'TRBV10-2', 'TRBV10-3', 
                 'TRBV11-1', 'TRBV11-2', 'TRBV11-3', 'TRBV12-1', 'TRBV12-2', 'TRBV12-3', 
                 'TRBV12-4', 'TRBV13', 'TRBV14', 'TRBV15', 'TRBV16', 'TRBV18', 'TRBV19', 
                 'TRBV2-1', 'TRBV2-2', 'TRBV2-3', 'TRBV20-1', 'TRBV21-1', 'TRBV22-1', 
                 'TRBV23-1', 'TRBV24-1', 'TRBV25-1', 'TRBV27', 'TRBV28', 'TRBV29-1', 
                 'TRBV3-1', 'TRBV3-2', 'TRBV3-3', 'TRBV3-4', 'TRBV30', 'TRBV4-1', 'TRBV4-2', 
                 'TRBV4-3', 'TRBV5-1', 'TRBV5-10', 'TRBV5-2', 'TRBV5-3', 'TRBV5-4', 
                 'TRBV5-5', 'TRBV5-6', 'TRBV5-7', 'TRBV5-8', 'TRBV5-9', 'TRBV6-1', 'TRBV6-2',
                 'TRBV6-3', 'TRBV6-4', 'TRBV6-5', 'TRBV6-6', 'TRBV6-7', 'TRBV7-10', 'TRBV7-2',
                 'TRBV7-3', 'TRBV7-4', 'TRBV7-5', 'TRBV7-6', 'TRBV7-7', 'TRBV7-9', 'TRBV9')

MACMUL_TRBJ <- c('TRBJ1-1', 'TRBJ1-2', 'TRBJ1-3', 'TRBJ1-4', 'TRBJ1-5', 'TRBJ1-6', 'TRBJ2-1', 'TRBJ2-2',
                 'TRBJ2-3', 'TRBJ2-4', 'TRBJ2-5', 'TRBJ2-6', 'TRBJ2-7')


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