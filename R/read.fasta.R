#' Read fasta formate and convert it to data frame
#'
#' @usage fasta.read(virus.fasta,host.fasta)
#'
#' @param virus.fasta  directory path to the virus fasta file.
#' @param host.fasta  directory path to the host fasta file.
#'
#' @return A list with two data frames.
#'
#' @note The list with two data.frames; the first one for virus DNA sequences and the second one for the host.
#'
#' @importFrom  Biostrings readDNAStringSet
#'
#' @examples
#' \dontrun{
#' fasta <- fasta.read("virus.fasta", "host.fasta")
#' fasta.v <- fasta[[1]]
#' fasta.h <- fasta[[2]]
#' }
#' @export
#'
#' @author Ali Mostafa Anwar \email{ali.mo.anwar@std.agr.cu.edu.eg} and Mohmed Soudy \email{MohmedSoudy2009@gmail.com}

fasta.read <- function(virus.fasta, host.fasta) {
  virus.fasta <- readDNAStringSet(virus.fasta)
  seq_name <- names(virus.fasta)
  sequence <- paste(virus.fasta)
  df.virus <- data.frame(seq_name, sequence)

  host.fasta <- readDNAStringSet(host.fasta)
  seq_name <- names(host.fasta)
  sequence <- paste(host.fasta)
  df.host <- data.frame(seq_name, sequence)

  return(list(df.virus, df.host))
}
