# vhcub

**Virus-Host Codon Usage Co-Adaptation Analysis**
 
<img align="right" width="250" height="250" src="https://github.com/AliYoussef96/vhcub/blob/master/logo.png">

vhcub is an R package to analyze the co-adaptation of codon usage between a virus and its host. 

## Installation

vhcub was developed using R and available on CRAN:

```R
install.packages("vhcub")
```

## The following measures are implemented in the package:

* ENC, effective number of codons ([Novembre, 2002](https://www.ncbi.nlm.nih.gov/pubmed/12140252)).

* SCUO, synonymous codon usage orderliness ([Wan et al., 2004](https://www.ncbi.nlm.nih.gov/pubmed/15222899)).

* Codon Adaptation Index, CAI ([Sharp and Li, 1987](https://www.ncbi.nlm.nih.gov/pubmed/3547335)).

* Relative  Codon Deoptimization Index, RCDI [(Puigbò et al, 2010)](https://bmcresnotes.biomedcentral.com/articles/10.1186/1756-0500-3-87).

* Relative Synonymous Codon Usage, RSCU ([Sharp and Li, 1987](https://www.ncbi.nlm.nih.gov/pubmed/3547335)).

**Using vhcub to study the CUB of a virus, its host and the co-adaptation between them is straightforward.**

For example;

coding sequences for both Escherichia virus T4 and its host Escherichia coli were downloaded in fasta format from the NCBI database.

```R
fasta <- fasta.read("EscherichiavirusT4.fasta","Escherichiacoli.fasta")
fasta.virus <- fasta[[1]]
fasta.host <- fasta[[2]]
gc.df <- GC.content(fasta.virus)
syncodon <- dinuc.syncodon(fasta.virus,permutations=100)
base <- dinuc.base(fasta.virus,permutations=100)
codon <- dinuc.codon(fasta.virus,permutations=100)
enc.df.virus <- ENc.values(fasta.virus)
enc.df.host <- ENc.values(fasta.host)
SCUO.df <- SCUO.values(fasta.virus)
cai.df <- CAI.values(fasta.virus, enc.df.host, fasta.host, genetic.code="11")
RSCU.virus <- RSCU.values(fasta.virus)
RSCU.host <- RSCU.values(fasta.host)
SiD <- SiD.value(RSCU.host,RSCU.virus)
rcdi.df <- RCDI.values(fasta.virus,fasta.host, enc.df.host)
```

vhcub package aiming to cover a complete workflow of CUB viruses analysis and its co-evolution with their host. Hence, further versions of vhcub aiming to improve, add more CUB measures and plots as well as statistical analysis to study the CUB. Hence, **contributions to the software are welcome**.



