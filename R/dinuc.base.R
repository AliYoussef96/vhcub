#' Statistical dinucleotide over- and underrepresentation (base model).
#'
#' A measure of statistical dinucleotide over- and underrepresentation; by allows for random sequence generation by shuffling (with/without replacement) of all bases in the sequence.
#'
#' For more information \href{https://www.rdocumentation.org/packages/seqinr/versions/3.6-1/topics/dinucleotides}{seqinr}.
#'
#' @usage dinuc.base(df.virus,permutations=500,exact_numbers = FALSE)
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
#' # Calculate zscore using (base model)
#' base <- dinuc.base(fasta.v, permutations = 500)
#' }
#' @export
#'
#' @author Ali Mostafa Anwar \email{ali.mo.anwar@std.agr.cu.edu.eg} and Mohmed Soudy \email{MohmedSoudy2009@gmail.com}


dinuc.base <- function(df.virus, permutations = 500, exact_numbers = FALSE) {
  dinuc.baseall <- data.frame()
  length <- 1:length(df.virus$seq_name)
  for (i.seq in length) {
    sequence <- s2c(tolower(as.character(df.virus$sequence[[i.seq]])))
    sequence.name <- as.character(df.virus$seq_name[[i.seq]])

    base <- zscore(sequence, simulations = permutations, modele = "base", exact = exact_numbers)


    base <- as.data.frame(base)
    aa <- base$Freq[base$Var1 == "aa"]
    ac <- base$Freq[base$Var1 == "ac"]
    ag <- base$Freq[base$Var1 == "ag"]
    at <- base$Freq[base$Var1 == "at"]

    tt <- base$Freq[base$Var1 == "tt"]
    ta <- base$Freq[base$Var1 == "ta"]
    tc <- base$Freq[base$Var1 == "tc"]
    tg <- base$Freq[base$Var1 == "tg"]

    gg <- base$Freq[base$Var1 == "gg"]
    ga <- base$Freq[base$Var1 == "ga"]
    gt <- base$Freq[base$Var1 == "gt"]
    gc <- base$Freq[base$Var1 == "gc"]

    cc <- base$Freq[base$Var1 == "cc"]
    ca <- base$Freq[base$Var1 == "ca"]
    ct <- base$Freq[base$Var1 == "ct"]
    cg <- base$Freq[base$Var1 == "cg"]

    dinuc.base <- NULL
    dinuc.base <- data.frame(
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
    dinuc.baseall <- rbind(dinuc.baseall, dinuc.base)
  }
  return(dinuc.baseall)
}
