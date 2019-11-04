#' Statistical dinucleotide over- and underrepresentation (syncodon model).
#'
#' A measure of statistical dinucleotide over- and underrepresentation; by allows for random sequence generation by shuffling (with/without replacement) of synonymous codons.
#'
#' For more information \href{https://www.rdocumentation.org/packages/seqinr/versions/3.6-1/topics/dinucleotides}{seqinr}.
#'
#' @usage dinuc.syncodon(df.virus,permutations=500,exact_numbers = FALSE)
#'
#' @param df.virus  data frame with seq_name and its DNA sequence.
#' @param permutations  the number of permutations for the z-score computation.
#' @param exact_numbers if TRUE exact analytical calculation will be used.
#'
#' @return A data.frame containing the computed statistic for each dinucleotide in all DNA sequences within df.virus.
#'
#' @import seqinr
#'
#' @examples
#' \dontrun{
#' # read DNA from fasta file
#' fasta <- fasta.read("virus.fasta", "host.fasta")
#' fasta.v <- fasta[[1]]
#' fasta.h <- fasta[[2]]
#' # Calculate zscore using (syncodon model)
#' syncodon <- dinuc.syncodon(fasta.v, permutations = 500)
#' }
#' @export
#'
#' @author Ali Mostafa Anwar \email{ali.mo.anwar@std.agr.cu.edu.eg} and Mohmed Soudy \email{MohmedSoudy2009@gmail.com}


dinuc.syncodon <- function(df.virus, permutations = 500, exact_numbers = FALSE) {
  dinuc.syncodonall <- data.frame()
  length <- 1:length(df.virus$seq_name)
  for (i.seq in length) {
    sequence <- s2c(tolower(as.character(df.virus$sequence[[i.seq]])))
    sequence.name <- as.character(df.virus$seq_name[[i.seq]])

    syncodon <- zscore(sequence,
      simulations = permutations,
      modele = "syncodon", exact = exact_numbers
    )

    syncodon <- as.data.frame(syncodon)
    aa <- syncodon$Freq[syncodon$Var1 == "aa"]
    ac <- syncodon$Freq[syncodon$Var1 == "ac"]
    ag <- syncodon$Freq[syncodon$Var1 == "ag"]
    at <- syncodon$Freq[syncodon$Var1 == "at"]

    tt <- syncodon$Freq[syncodon$Var1 == "tt"]
    ta <- syncodon$Freq[syncodon$Var1 == "ta"]
    tc <- syncodon$Freq[syncodon$Var1 == "tc"]
    tg <- syncodon$Freq[syncodon$Var1 == "tg"]

    gg <- syncodon$Freq[syncodon$Var1 == "gg"]
    ga <- syncodon$Freq[syncodon$Var1 == "ga"]
    gt <- syncodon$Freq[syncodon$Var1 == "gt"]
    gc <- syncodon$Freq[syncodon$Var1 == "gc"]

    cc <- syncodon$Freq[syncodon$Var1 == "cc"]
    ca <- syncodon$Freq[syncodon$Var1 == "ca"]
    ct <- syncodon$Freq[syncodon$Var1 == "ct"]
    cg <- syncodon$Freq[syncodon$Var1 == "cg"]

    dinuc.syncodon <- NULL
    dinuc.syncodon <- data.frame(
      gene.name = sequence.name,
      aa = aa,
      ac = ac,
      ag = ag,
      at = at,
      tt = tt,
      ta = ta,
      tc = ta,
      tg = tg,
      gg = gg,
      ga = ga,
      gt = gt,
      gc = gc,
      cc = cc,
      ca = ca,
      ct = ct,
      cg = cg
    )
    dinuc.syncodonall <- rbind(dinuc.syncodonall, dinuc.syncodon)
  }
  return(dinuc.syncodonall)
}
