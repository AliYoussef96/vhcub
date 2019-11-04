#' Effective Number of Codons (ENc).
#'
#' Measure the Effective Number of Codons (ENc) of DNA sequence. Using its modified version (Novembre, 2002).
#'
#' For more information about ENc \href{https://academic.oup.com/mbe/article/19/8/1390/997706}{Novembre, 2002}.
#'
#' @usage ENc.values(df.fasta,genetic.code = "1",threshold=0)
#'
#' @param df.fasta  a data frame with seq_name and its DNA sequence.
#' @param genetic.code  a single string that uniquely identifies a genetic code to use.
#' @param threshold optional numeric, specifying sequence length, in codons, used for filtering.
#'
#' @return A data.frame containing the computed ENc values for each DNA sequences within df.fasta.
#'
#' @import coRdon
#' @importFrom  Biostrings DNAStringSet
#'
#' @examples
#' \dontrun{
#' # read DNA from fasta file
#' fasta <- fasta.read("virus.fasta", "host.fasta")
#' fasta.v <- fasta[[1]]
#' fasta.h <- fasta[[2]]
#' # Calculate ENc
#' enc.df <- ENc.values(fasta.v)
#' }
#' @export
#'
#' @author Ali Mostafa Anwar \email{ali.mo.anwar@std.agr.cu.edu.eg} and Mohmed Soudy \email{MohmedSoudy2009@gmail.com}
#'

ENc.values <- function(df.fasta, genetic.code = "1", threshold = 0) {
  length <- 1:length(df.fasta$seq_name)
  df.enc.all <- data.frame()
  for (i_seq in length) {
    sequence <- as.character(df.fasta$sequence[[i_seq]])
    sequence <- str_sub(sequence, start = 1, end = (nchar(sequence) - nchar(sequence) %% 3))
    seq_name <- df.fasta$seq_name[[i_seq]]

    dna <- DNAStringSet(c(sequence, "NNN"))
    cT <- codonTable(dna)
    ENc <- ENC(cT,
      id_or_name2 = genetic.code,
      alt.init = TRUE, stop.rm = TRUE,
      filtering = "none", len.threshold = threshold
    )[[1]]

    df.enc <- NULL
    df.enc <- data.frame(gene.name = seq_name, ENc = ENc)
    df.enc.all <- rbind(df.enc.all, df.enc)
  }
  return(df.enc.all)
}
