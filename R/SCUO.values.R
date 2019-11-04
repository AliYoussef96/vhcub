#' Synonymous codon usage eorderliness (SCUO)
#'
#' Measure the Synonymous Codon Usage Eorderliness (SCUO) of DNA sequence (Wan et al., 2004).
#'
#' For more information about ENc \href{https://bmcevolbiol.biomedcentral.com/articles/10.1186/1471-2148-4-19}{Wan et al., 2004}.
#'
#' @usage SCUO.values(df.fasta,genetic.code = "1",threshold=0)
#'
#' @param df.fasta  a data frame with seq_name and its DNA sequence.
#' @param genetic.code  a single string that uniquely identifies a genetic code to use.
#' @param threshold optional numeric, specifying sequence length, in codons, used for filtering.
#'
#' @return A data.frame containing the computed SCUO values for each DNA sequences within df.fasta.
#'
#' @import coRdon
#' @import seqinr
#' @importFrom  Biostrings DNAStringSet
#'
#' @examples
#' \dontrun{
#' # read DNA from fasta file
#' fasta <- fasta.read("virus.fasta", "host.fasta")
#' fasta.v <- fasta[[1]]
#' fasta.h <- fasta[[2]]
#' # Calculate SCUO
#' SCUO.df <- SCUO.values(fasta.v)
#' }
#' @export
#'
#' @author Ali Mostafa Anwar \email{ali.mo.anwar@std.agr.cu.edu.eg} and Mohmed Soudy \email{MohmedSoudy2009@gmail.com}

SCUO.values <- function(df.fasta, genetic.code = "1", threshold = 0) {
  length <- 1:length(df.fasta$seq_name)
  df.SCUO.all <- data.frame()
  for (i_seq in length) {
    sequence <- as.character(df.fasta$sequence[[i_seq]])
    sequence <- str_sub(sequence, start = 1, end = (nchar(sequence) - nchar(sequence) %% 3))
    seq_name <- df.fasta$seq_name[[i_seq]]

    dna <- DNAStringSet(c(sequence, "NNN"))
    cT <- codonTable(dna)
    SCUO <- SCUO(cT,
      id_or_name2 = genetic.code,
      alt.init = TRUE, stop.rm = TRUE, filtering = "none",
      len.threshold = threshold
    )[[1]]

    df.SCUO <- NULL
    df.SCUO <- data.frame(gene.name = seq_name, SCUO = SCUO)
    df.SCUO.all <- rbind(df.SCUO.all, df.SCUO)
  }
  return(df.SCUO.all)
}
