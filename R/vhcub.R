#' vhcub: A package to analysis the co-adaptation of codon usage between a virus and its host.
#'
#'
#' vhcub can calculate various codon usage bias measurements as; effective number of codons (ENc),
#' codon adaptation index (CAI), relative  codon deoptimization index (RCDI), similarity index (SiD),
#' synonymous codon usage eorderliness (SCUO) and, relative synonymous codon usage (RSCU).
#' Also, it provides a statistical dinucleotide over- and underrepresentation with three different models.
#' Implement several methods for visualization of codon usage as ENc.GC3plot and PR2.plot.
#'
#'
#' @section vhcub functions:
#'
#' fasta.read: read fasta format files and convert it to data.frame.
#'
#' GC.content: calculates overall GC content as well as GC at first, second, and third codon positions.
#'
#' RSCU.values: measure the Relative Synonymous Codon Usage (RSCU) of DNA sequence.
#'
#' SCUO.values: measure the Synonymous Codon Usage Eorderliness (SCUO) of DNA sequence.
#'
#' RCDI.values: measure the Relative Codon Deoptimization Index (RCDI) of DNA sequence.
#'
#' CAI.values: measure the Codon Adaptation  Index  (CAI)  Sharp  and  Li  (1987), of DNA sequence.
#'
#' ENc.values: measure the Effective Number of Codons (ENc) of DNA sequence. Using its modified version.
#'
#' dinuc.syncodon: measure of statistical dinucleotide over- and underrepresentation; by allows for random sequence generation by shuffling (with/without replacement) of synonymous codons.
#'
#' dinuc.codon: measure of statistical dinucleotide over- and underrepresentation; by allows for random sequence generation by shuffling (with/without replacement) of codons.
#'
#' dinuc.base: measure of statistical dinucleotide over- and underrepresentation; by allows for random sequence generation by shuffling (with/without replacement) of all bases in the sequence.
#'
#' ENc.GC3plot: make an ENc-GC3 scatterplot. Where the y-axis represents the ENc values and the x-axis represents the GC3 content. The red fitting line shows the expected ENc values when codon usage bias affected solely by GC3.
#'
#' PR2.plot: make a Parity rule 2 (PR2) plot, where the AT-bias [A3/(A3 +T3)] at the third codon position of the four-codon amino acids of entire genes is the ordinate and the GC-bias [G3/(G3 +C3)] is the abscissa. The center of the plot, where both coordinates are 0.5, is where A = U and G = C (PR2), with no bias between the influence of the mutation and selection rates.
#'
#' @examples
#' \dontrun{
#' # read DNA from fasta files
#' fasta <- fasta.read("virus.fasta", "host.fasta")
#' fasta.v <- fasta[[1]]
#' fasta.h <- fasta[[2]]
#' # calculate GC content
#' gc.df <- GC.content(fasta.v)
#' # measure of statistical dinucleotide over- and underrepresentation
#' syncodon <- dinuc.syncodon(fasta.v,permutations=10)
#' base <- dinuc.base(fasta.v,permutations=10)
#' codon <- dinuc.codon(fasta.v,permutations=10)
#' # calculate ENc
#' enc.df <- ENc.values(fasta.v)
#' enc.df.h <- ENc.values(fasta.h)
#' # calculate SCUO and CAI
#' SCUO.df <- SCUO.values(fasta.v)
#' cai.df <- CAI.values(fasta.v,enc.df.h, fasta.h)
#' # calculate RSCU
#' RSCU.H <- RSCU.values(fasta.h)
#' RSCU.V <- RSCU.values(fasta.v)
#' # calculate SiD
#' SiD <- SiD.value(RSCU.H,RSCU.V)
#' # calculate RCDI
#' rcdi.df <- RCDI.values(fasta.v,fasta.h, enc.df.h)
#' # plot ENc.GC3plot
#' ENc.GC3plot(enc.df,gc.df)
#' # plot PR2.plot
#' PR2.plot(fasta.v)
#' }
#'
#' @author Ali Mostafa Anwar \email{ali.mo.anwar@std.agr.cu.edu.eg} and Mohmed Soudy \email{MohmedSoudy2009@gmail.com}
#'
#' @docType package
#' @name vhcub
#'
NULL
