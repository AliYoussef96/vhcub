#' GC content
#'
#' Calculates overall GC content as well as GC at first, second, and third codon positions.
#'
#' @usage GC.content(df.virus)
#'
#' @param df.virus  data frame with seq_name and its DNA sequence.
#'
#' @return A data.frame with overall GC content as well as GC at first, second, and third codon positions of all DNA sequence from df.virus.
#'
#' @import seqinr
#'
#' @examples
#' \dontrun{
#' # read DNA from fasta file
#' fasta <- fasta.read("virus.fasta", "host.fasta")
#' fasta.v <- fasta[[1]]
#' fasta.h <- fasta[[2]]
#' # Calculate GC content
#' gc.df <- GC.content(fasta.v)
#' }
#' @export
#'
#' @author Ali Mostafa Anwar \email{ali.mo.anwar@std.agr.cu.edu.eg} and Mohmed Soudy \email{MohmedSoudy2009@gmail.com}

GC.content <- function(df.virus) {
  df.all.GC <- data.frame()
  length <- 1:length(df.virus$seq_name)
  for (i_seq in length) {

    sequence <- as.character(df.virus$sequence[[i_seq]])
    seq_name <- df.virus$seq_name[[i_seq]]
    gc <- GC(s2c(sequence))
    gc1 <- GCpos(s2c(sequence), "1")
    gc2 <- GCpos(s2c(sequence), "2")
    gc3 <- GCpos(s2c(sequence), "3")
    df.gc <- data.frame(gene.name = seq_name, GC = gc, GC1 = gc1, GC2 = gc2, GC3 = gc3)
    df.all.GC <- rbind(df.all.GC,df.gc)

  }
  return(df.all.GC)
}
