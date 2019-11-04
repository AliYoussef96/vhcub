#' Codon Adaptation  Index  (CAI)
#'
#' Measure the Codon Adaptation  Index  (CAI)  Sharp  and  Li  (1987), of DNA sequence.
#'
#' For more information about CAI \href{https://academic.oup.com/nar/article-abstract/15/3/1281/1166844?redirectedFrom=fulltext}{Sharp  and  Li, 1987}.
#'
#' @usage CAI.values(df.virus, ENc.set.host,
#'                        df.host,genetic.code = "1",set.len = 5, threshold = 0)
#'
#' @param df.virus  a data frame with seq_name and its virus DNA sequence.
#' @param ENc.set.host  a data frame with ENc values of a host.
#' @param df.host a data frame with seq_name and its host DNA sequence.
#' @param genetic.code  a single string that uniquely identifies a genetic code to use.
#' @param set.len  a number represents a percent that will be used as reference genes from the total host genes.
#' @param threshold optional numeric, specifying sequence length, in codons, used for filtering.
#'
#' @return A data.frame containing the computed CAI values for each DNA sequences within df.fasta.
#'
#' @import coRdon
#' @import stringr
#' @importFrom  Biostrings DNAStringSet
#'
#' @examples
#' \dontrun{
#' # read DNA from fasta file
#' fasta <- fasta.read("virus.fasta", "host.fasta")
#' fasta.v <- fasta[[1]]
#' fasta.h <- fasta[[2]]
#' # Calculate CAI
#' enc.df.h <- ENc.values(fasta.h)
#' cai.df <- CAI.values(fasta.v, enc.df.h, fasta.h)
#' }
#' @export
#' @author Ali Mostafa Anwar \email{ali.mo.anwar@std.agr.cu.edu.eg} and Mohmed Soudy \email{MohmedSoudy2009@gmail.com}
#'

CAI.values <- function(df.virus, ENc.set.host, df.host, genetic.code = "1",
                       set.len = 5, threshold = 0) {

  # get the refrence set from ENc values of host
  newENc <- ENc.set.host[order(ENc.set.host$ENc), ]
  set.len <- length(newENc$gene.name) * (set.len / 100)
  gene.set <- newENc$gene.name[1:set.len]
  gene.set <- df.host[df.host$seq_name %in% gene.set, ]
  dna.set <- as.vector(gene.set$sequence)
  # this function will make the sequence len&&3 = 0
  firstframe <- function(sequence) {
    sequence <- str_sub(sequence, start = 1, end = (nchar(sequence) - nchar(sequence) %% 3))
    return(sequence)
  }

  dna.set <- lapply(dna.set, function(x) firstframe(x))
  dna.set <- unlist(dna.set, use.names = FALSE)

  # calc. codontable for dna set (ref gene set)
  dna.set <- DNAStringSet(dna.set)
  cT.set <- codonTable(dna.set)

  # calc CAI for virus dna

  length <- 1:length(df.virus$seq_name)
  df.cai.all <- data.frame()
  for (i_seq in length) {
    sequence <- as.character(df.virus$sequence[[i_seq]])
    sequence <- str_sub(sequence, start = 1, end = (nchar(sequence) - nchar(sequence) %% 3))
    seq_name <- df.virus$seq_name[[i_seq]]

    dna <- DNAStringSet(c(sequence, "NNN"))
    cT <- codonTable(dna)

    cai <- CAI(cT,
      subsets = list(cT.set), ribosomal = FALSE,
      id_or_name2 = genetic.code, alt.init = TRUE,
      stop.rm = TRUE, filtering = "none",
      len.threshold = threshold
    )[[1]]

    df.cai <- NULL
    df.cai <- data.frame(gene.name = seq_name, CAI = cai)
    df.cai.all <- rbind(df.cai.all, df.cai)
  }
  return(df.cai.all)
}
