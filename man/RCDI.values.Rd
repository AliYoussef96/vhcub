% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RCDI.values.R
\name{RCDI.values}
\alias{RCDI.values}
\title{Relative Codon Deoptimization Index (RCDI)}
\usage{
RCDI.values(fasta.virus, fasta.host, enc.host, set.len= 5)
}
\arguments{
\item{fasta.virus}{a data frame with virus seq_name and its DNA sequence.}

\item{fasta.host}{a data frame with host seq_name and its DNA sequence.}

\item{enc.host}{a data frame of a hosts' ENc values.}

\item{set.len}{a number represents a percent that will be used as reference genes from the total host genes.}
}
\value{
A data.frame containing the computed ENc values for each DNA sequences within df.fasta.
}
\description{
Measure the Relative Codon Deoptimization Index (RCDI) of DNA sequence.
}
\details{
For more information about RCDI \href{https://bmcresnotes.biomedcentral.com/articles/10.1186/1756-0500-3-87}{Puigbò et al., 2010}
}
\examples{
\dontrun{
# read DNA from fasta file
fasta <- fasta.read("virus.fasta", "host.fasta")
fasta.v <- fasta[[1]]
fasta.h <- fasta[[2]]
# Calculate RCDI
enc.df.h <- ENc.values(fasta.h)
rcdi.df <- RCDI.values(fasta.v, fasta.h, enc.df.h)
}
}
\author{
Ali Mostafa Anwar \email{ali.mo.anwar@std.agr.cu.edu.eg} and Mohmed Soudy \email{MohmedSoudy2009@gmail.com}
}
