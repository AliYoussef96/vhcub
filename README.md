# vhcub

[![](https://www.r-pkg.org/badges/version/vhcub?color=green)](https://cran.r-project.org/package=vhcub)
[![](http://cranlogs.r-pkg.org/badges/last-week/vhcub?color=green)](https://cran.r-project.org/package=vhcub)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/0de2a9a9a9e74f19ab64de2419bf1cc4)](https://www.codacy.com/manual/AliYoussef96/vhcub?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=AliYoussef96/vhcub&amp;utm_campaign=Badge_Grade)
[![Build Status](https://travis-ci.com/AliYoussef96/vhcub.svg?branch=master)](https://travis-ci.com/AliYoussef96/vhcub)
[![](https://img.shields.io/badge/doi-https%3A%2F%2Fdoi.org%2F10.1016%2Fj.jprot.2019.103613-red)](https://doi.org/10.12688/f1000research.21763.1)
[![CRAN RStudio mirror downloads](https://cranlogs.r-pkg.org/badges/grand-total/vhcub?color=blue)](https://cran.r-project.org/package=vhcub) 
[![](https://badgen.net/badge/Citations/1/:color?icon=github)](https://f1000research.com/articles/8-2137/v1)

## Virus-Host Codon Usage Co-Adaptation Analysis
 
<img align="right" width="250" height="250" src="https://github.com/AliYoussef96/vhcub/blob/master/logo2.png">

vhcub is an R package to analyze the co-adaptation of codon usage between a virus and its host. 



### Installation

1. Using devtools:

```R
if (!requireNamespace("devtools", quietly=TRUE)){
        install.packages("devtools")}
devtools::install_github('AliYoussef96/vhcub')
```


2. vhcub was developed using R and available on CRAN:

```R
install.packages("vhcub")
```

3. If 1 and 2 installation instructions failed, please try:

```R
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("coRdon")
BiocManager::install("Biostrings")


install.packages("vhcub")
```

### The following measures are implemented in the package

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
# read virus and host fasta files
fasta <- fasta.read("EscherichiavirusT4.fasta","Escherichiacoli.fasta")
fasta.virus <- fasta[[1]]
fasta.host <- fasta[[2]]

# Calculate the GC overall all content as well as GC at first, second and third codon positions for the virus
gc.df <- GC.content(fasta.virus)

# Calculate zscore using syncodon model for statistical dinucleotide over- and underrepresentation
syncodon <- dinuc.syncodon(fasta.virus,permutations=100)

# Calculate zscore using base model for statistical dinucleotide over- and underrepresentation
base <- dinuc.base(fasta.virus,permutations=100)

# Calculate zscore using codon model for statistical dinucleotide over- and underrepresentation
codon <- dinuc.codon(fasta.virus,permutations=100)

# Calculate ENc values for the virus and its host
enc.df.virus <- ENc.values(fasta.virus)
enc.df.host <- ENc.values(fasta.host)

# Calculate SCUO values for the virus
scuo.df <- SCUO.values(fasta.virus)

# Calculate CAI values for the virus using the host sequences as a reference genes set
cai.df <- CAI.values(fasta.virus, enc.df.host, fasta.host, genetic.code="11")

# Calculate RSCU values for the virus and its host
rscu.virus <- RSCU.values(fasta.virus) 
rscu.host <- RSCU.values(fasta.host)

# Calculate SiD value for the virus 
SiD <- SiD.value(rscu.host,rscu.virus)

# Calculate RCDI values for the virus
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

Anwar AM, Soudy M and Mohamed R. vhcub: Virus-host codon usage co-adaptation analysis \[version 1; peer review: 2 approved]. F1000Research 2019, 8:2137 (<https://doi.org/10.12688/f1000research.21763.1>)
