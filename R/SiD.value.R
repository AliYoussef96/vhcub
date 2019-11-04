#' Similarity Index (SiD)
#'
#' Measure the Similarity Index (SiD) between a virus and its host codon usage.
#'
#' For more information about SiD \href{https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0077239}{Zhou et al., 2013}.
#'
#'
#' @usage SiD.value(rscu.host,rscu.virus)
#'
#' @param rscu.host  a data frame with RSCU a host codon values.
#' @param rscu.virus  a data frame with RSCU a virus codon values.
#'
#' @return A numeric represent a SiD value.
#'
#' @examples
#' \dontrun{
#' # read DNA from fasta file
#' fasta <- fasta.read("virus.fasta", "host.fasta")
#' fasta.v <- fasta[[1]]
#' fasta.h <- fasta[[2]]
#' # Calculate SiD
#' RSCU.H <- RSCU.values(fasta.h)
#' RSCU.V <- RSCU.values(fasta.v)
#' SiD <- SiD.value(RSCU.H, RSCU.V)
#' }
#' @export
#'
#' @author Ali Mostafa Anwar \email{ali.mo.anwar@std.agr.cu.edu.eg} and Mohmed Soudy \email{MohmedSoudy2009@gmail.com}
#'

SiD.value <- function(rscu.host, rscu.virus) {
  df.rscu.host <- data.frame()
  length <- 1:length(rscu.host)
  for (i_mean in length) {
    df.rscu <- NULL
    means <- mean(rscu.host[[i_mean]], na.rm = TRUE)
    codon <- colnames(rscu.host)
    codon <- codon[i_mean]
    df.rscu <- data.frame(codon = codon, rscu.host = means)
    df.rscu.host <- rbind(df.rscu.host, df.rscu)
  }

  df.rscu.virus <- data.frame()
  length <- 1:length(rscu.virus)
  for (i_mean in length) {
    df.rscu <- NULL
    means <- mean(rscu.virus[[i_mean]], na.rm = TRUE)
    codon <- colnames(rscu.virus)
    codon <- codon[i_mean]
    df.rscu <- data.frame(codon = codon, rscu.virus = means)
    df.rscu.virus <- rbind(df.rscu.virus, df.rscu)
  }

  rscu.df.all <- merge(df.rscu.host, df.rscu.virus, by = "codon")
  rscu.df.all$rscu.all <- rscu.df.all$rscu.host * rscu.df.all$rscu.virus

  up <- sum(rscu.df.all$rscu.all)
  down <- sqrt((sum(rscu.df.all$rscu.host)^2) * (sum(rscu.df.all$rscu.virus)^2))
  R.a.b <- up / down

  D.a.b <- (1 - R.a.b) / 2

  return(D.a.b)
}
