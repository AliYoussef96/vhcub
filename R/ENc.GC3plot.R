#' ENc-GC3 scatterplot.
#'
#' Make an ENc-GC3 scatterplot. Where the y-axis represents the ENc values and the x-axis represents the GC3 content.
#' The red fitting line shows the expected ENc values when codon usage bias affected solely by GC3.
#'
#' For more information about ENc-GC3 plot \href{https://www.tandfonline.com/doi/full/10.1038/emi.2016.106}{Butt et al., 2016}.
#'
#' @usage ENc.GC3plot(enc.df, gc.df)
#'
#' @param enc.df  a data frame with ENc values.
#' @param gc.df  a data frame with GC3 values.
#'
#' @return A ggplot object.
#'
#' @import ggplot2
#'
#' @examples
#' \dontrun{
#' # read DNA from fasta file
#' fasta <- fasta.read("virus.fasta", "host.fasta")
#' fasta.v <- fasta[[1]]
#' fasta.h <- fasta[[2]]
#'
#' enc.df <- ENc.values(fasta.v)
#' gc.df <- GC.content(fasta.v)
#' ENc.GC3plot(enc.df, gc.df)
#' }
#' @export
#'
#' @author Ali Mostafa Anwar \email{ali.mo.anwar@std.agr.cu.edu.eg} and Mohmed Soudy \email{MohmedSoudy2009@gmail.com}
#'

ENc.GC3plot <- function(enc.df, gc.df) {
  x <- NULL
  eq <- function(x) {
    2 + x + (29 / (x^2 + (1 - x)^2))
  }
  plot <- ggplot() + geom_point(data = enc.df, aes(x = gc.df$GC3, y = enc.df$ENc)) +
    stat_function(fun = eq, geom = "line", color = "red", size = 1, data = data.frame(x = c(seq(0, 1, 0.001))), aes(x)) +
    theme_classic() + xlab("GC3") + ylab("ENc")

  return(plot)
}
