# vhcub

## Virus-Host Codon Usage Co-Adaptation Analysis
 
<img align="right" width="250" height="250" src="https://github.com/AliYoussef96/vhcub/blob/master/logo.png">

vhcub is an R package to analyze the co-adaptation of codon usage between a virus and its host. 

### Installation

vhcub was developed using R and available on CRAN:

```R
install.packages("vhcub")
```

### The following measures are implemented in the package:

* ENC, effective number of codons ([Novembre, 2002](https://www.ncbi.nlm.nih.gov/pubmed/12140252)),

* SCUO, synonymous codon usage orderliness ([Wan et al., 2004](https://www.ncbi.nlm.nih.gov/pubmed/15222899)),

* Codon Adaptation Index, CAI ([Sharp and Li, 1987](https://www.ncbi.nlm.nih.gov/pubmed/3547335)),

* Relative  Codon Deoptimization Index, RCDI [(Puigbò et al, 2010)](https://bmcresnotes.biomedcentral.com/articles/10.1186/1756-0500-3-87),

* Relative Synonymous Codon Usage, RSCU ([Sharp and Li, 1987](https://www.ncbi.nlm.nih.gov/pubmed/3547335)),

* Also, it provides a statistical dinucleotide over- and underrepresentation with three different models.

### Usage

Using vhcub to study the CUB of a virus, its host and the co-adaptation between them is straightforward.

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

Furthermore, vhcub uses ggplot to visualize two important plots named ENc-GC3 plot and PR2plot, which help to explain what are the factors influence a virus's evolution concerning its CUB.

```R
ENc.GC3plot(enc.df.virus,gc.df)
PR2.plot(fasta.virus)
```
### Contribution Guidelines

**Contributions to the package are welcome**

For bugs and suggestions, the most effective way is by raising an issue on the github issue tracker. 
Github allows you to classify your issues so that we know if it is a bug report, feature request or feedback to the authors.

If you wish to contribute some changes to the code then you should submit a [pull request](https://github.com/AliYoussef96/BCAW-Tool/pulls)
How to create a Pull Request? [documentation on pull requests](https://help.github.com/en/articles/about-pull-requests)

## Citation

