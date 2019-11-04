#' Statistical dinucleotide over- and underrepresentation (codon model).
#'
#' A measure of statistical dinucleotide over- and underrepresentation; by allows for random sequence generation by shuffling (with/without replacement) of codons.
#'
#' For more information \href{https://www.rdocumentation.org/packages/seqinr/versions/3.6-1/topics/dinucleotides}{seqinr}.
#'
#' @usage dinuc.codon(df.virus,permutations=500,exact_numbers = FALSE)
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
#' # Calculate zscore using (codon model)
#' codon <- dinuc.codon(fasta.v, permutations = 500)
#' }
#' @export
#'
#' @author Ali Mostafa Anwar \email{ali.mo.anwar@std.agr.cu.edu.eg} and Mohmed Soudy \email{MohmedSoudy2009@gmail.com}

dinuc.codon <- function(df.virus, permutations = 500, exact_numbers = FALSE) {
  dinuc.codonall <- data.frame()
  length <- 1:length(df.virus$seq_name)
  for (i.seq in length) {
    sequence <- s2c(tolower(as.character(df.virus$sequence[[i.seq]])))
    sequence.name <- as.character(df.virus$seq_name[[i.seq]])

    codon <- zscore(sequence, simulations = permutations, modele = "codon", exact = exact_numbers)

    codon <- as.data.frame(codon)
    aa <- codon$Freq[codon$Var1 == "aa"]
    ac <- codon$Freq[codon$Var1 == "ac"]
    ag <- codon$Freq[codon$Var1 == "ag"]
    at <- codon$Freq[codon$Var1 == "at"]

    tt <- codon$Freq[codon$Var1 == "tt"]
    ta <- codon$Freq[codon$Var1 == "ta"]
    tc <- codon$Freq[codon$Var1 == "tc"]
    tg <- codon$Freq[codon$Var1 == "tg"]

    gg <- codon$Freq[codon$Var1 == "gg"]
    ga <- codon$Freq[codon$Var1 == "ga"]
    gt <- codon$Freq[codon$Var1 == "gt"]
    gc <- codon$Freq[codon$Var1 == "gc"]

    cc <- codon$Freq[codon$Var1 == "cc"]
    ca <- codon$Freq[codon$Var1 == "ca"]
    ct <- codon$Freq[codon$Var1 == "ct"]
    cg <- codon$Freq[codon$Var1 == "cg"]

    dinuc.codon <- NULL
    dinuc.codon <- data.frame(
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

    dinuc.codonall <- rbind(dinuc.codonall, dinuc.codon)
  }
  return(dinuc.codonall)
}
