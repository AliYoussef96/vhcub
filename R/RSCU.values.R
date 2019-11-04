#' Relative Synonymous Codon Usage (RSCU)
#'
#' Measure the Relative Synonymous Codon Usage (RSCU) of DNA sequence.
#'
#' For more information about ENc \href{https://academic.oup.com/nar/article-abstract/14/13/5125/1143812?redirectedFrom=fulltext}{Sharp et al., 1986}.
#'
#' @usage RSCU.values(df.fasta)
#'
#' @param df.fasta  a data frame with seq_name and its DNA sequence.
#'
#' @return A data.frame containing the computed RSCU values for each codon for each DNA sequences within df.fasta.
#'
#' @import seqinr
#'
#' @examples
#' \dontrun{
#' # read DNA from fasta file
#' fasta <- fasta.read("virus.fasta", "host.fasta")
#' fasta.v <- fasta[[1]]
#' fasta.h <- fasta[[2]]
#' # Calculate RSCU
#' RSCU.H <- RSCU.values(fasta.h)
#' RSCU.V <- RSCU.values(fasta.v)
#' }
#' @export
#'
#' @author Ali Mostafa Anwar \email{ali.mo.anwar@std.agr.cu.edu.eg} and Mohmed Soudy \email{MohmedSoudy2009@gmail.com}
#'

RSCU.values <- function(df.fasta) {
  codons <- uco(s2c("aaa"), index = "rscu", NA.rscu = 0)
  codons <- as.data.frame(codons)
  codons.names <- as.vector(row.names(codons))

  rscu.df <- data.frame(row.names = codons.names)
  length <- 1:length(df.fasta$seq_name)
  for (i_seq in length) {
    sequence <- as.character(df.fasta$sequence[[i_seq]])
    seq_name <- as.character(df.fasta$seq_name[[i_seq]])
    rscu <- uco(s2c(sequence),
      index = "rscu",
      as.data.frame = FALSE, NA.rscu = 0
    )

    rscu <- as.data.frame(rscu)
    colnames(rscu) <- seq_name
    rscu.df <- cbind(rscu.df, rscu)
  }
  rscu.df <- as.data.frame(t(rscu.df))
  return(rscu.df)
}
