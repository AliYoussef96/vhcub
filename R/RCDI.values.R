#' Relative Codon Deoptimization Index (RCDI)
#'
#' Measure the Relative Codon Deoptimization Index (RCDI) of DNA sequence.
#'
#' For more information about RCDI \href{https://bmcresnotes.biomedcentral.com/articles/10.1186/1756-0500-3-87}{Puigbò et al., 2010}
#'
#' @usage RCDI.values(fasta.virus, fasta.host, enc.host, set.len= 5)
#'
#' @param fasta.virus  a data frame with virus seq_name and its DNA sequence.
#' @param fasta.host  a data frame with host seq_name and its DNA sequence.
#' @param enc.host   a data frame of a hosts' ENc values.
#' @param set.len  a number represents a percent that will be used as reference genes from the total host genes.
#'
#' @return A data.frame containing the computed ENc values for each DNA sequences within df.fasta.
#'
#' @importFrom  Biostrings DNAStringSet
#'
#' @examples
#' \dontrun{
#' # read DNA from fasta file
#' fasta <- fasta.read("virus.fasta", "host.fasta")
#' fasta.v <- fasta[[1]]
#' fasta.h <- fasta[[2]]
#' # Calculate RCDI
#' enc.df.h <- ENc.values(fasta.h)
#' rcdi.df <- RCDI.values(fasta.v, fasta.h, enc.df.h)
#' }
#' @export
#'
#' @author Ali Mostafa Anwar \email{ali.mo.anwar@std.agr.cu.edu.eg} and Mohmed Soudy \email{MohmedSoudy2009@gmail.com}
#'
#'

RCDI.values <- function(fasta.virus, fasta.host, enc.host, set.len = 5) {
  newENc <- enc.host[order(enc.host$ENc), ]
  set.len <- length(newENc$gene.name) * (set.len / 100)
  gene.set <- newENc$gene.name[1:set.len]
  gene.set <- fasta.host[fasta.host$seq_name %in% gene.set, ]

  rscu.virus <- RSCU.values(fasta.virus)

  rscu.ref <- RSCU.values(gene.set)
  df.rscu.ref <- data.frame()
  length <- 1:length(rscu.ref)
  for (i_mean in length) {
    df.rscu <- NULL
    means <- mean(rscu.ref[[i_mean]], na.rm = TRUE)
    codon <- colnames(rscu.ref)
    codon <- codon[i_mean]
    df.rscu <- data.frame(codon = codon, rscu.ref = means)
    df.rscu.ref <- rbind(df.rscu.ref, df.rscu)
  }
  codon.name <- df.rscu.ref$codon
  df.rscu.ref <- as.data.frame(t(df.rscu.ref))
  colnames(df.rscu.ref) <- codon.name
  df.rscu.ref <- df.rscu.ref[-c(1), ]


  RCDI.df <- data.frame()
  length <- 1:length(fasta.virus$seq_name)
  for (i_seq in length) {
    sequence <- as.character(fasta.virus$sequence[[i_seq]])
    firstframe <- function(sequence) {
      sequence <- str_sub(sequence, start = 1, end = (nchar(sequence) - nchar(sequence) %% 3))
      return(sequence)
    }
    sequence <- firstframe(sequence)
    seq_name <- as.character(fasta.virus$seq_name[[i_seq]])

    rscu <- uco(s2c(sequence),
      index = "rscu",
      as.data.frame = FALSE, NA.rscu = 0
    )
    rscu <- as.data.frame(t(rscu))
    rownames(rscu) <- "rscu"


    dna <- DNAStringSet(c(sequence, "NNN"))
    count <- codonTable(dna)
    count <- as.data.frame(count@counts[1, ])
    count <- as.data.frame(t(count))
    colnames(count) <- tolower(colnames(count))
    rownames(count) <- "count"

    CiFa <- rbind(rscu, count)
    CiFa.CiFh <- rbind(CiFa, df.rscu.ref)
    CiFa.CiFh <- as.data.frame(t(CiFa.CiFh))
    N <- as.numeric(floor(nchar(sequence) / 3))


    CiFa.CiFh$rscu <- as.numeric(CiFa.CiFh$rscu)
    CiFa.CiFh$rscu.ref <- as.numeric(CiFa.CiFh$rscu.ref)
    CiFa.CiFh$count <- as.numeric(CiFa.CiFh$count)

    CiFa.CiFh$RCDI <- ((CiFa.CiFh$rscu / CiFa.CiFh$rscu.ref) * CiFa.CiFh$count) / N
    RCDI <- sum(CiFa.CiFh$RCDI)

    df <- NULL
    df <- data.frame(gene.name = seq_name, RCDI = RCDI)
    RCDI.df <- rbind(RCDI.df, df)
  }
  return(RCDI.df)
}
