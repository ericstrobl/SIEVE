# SIEVE: Locus-Anchored Drug Prioritization for Complex Disorders

SIEVE prioritizes drugs for complex disorders by first decomposing broad clinical phenotypes into locus-anchored subphenotypes, then using genotype-coupled reference expression to infer a genetically calibrated mechanism vector for each subphenotype. It projects both mechanism vectors and perturbational signatures away from negative-control expression programs and aggregates concordance across experimental contexts to produce robust, locus-specific drug rankings.

# Installation

> library(devtools)

> install_github("ericstrobl/SIEVE")

> library(SIEVE)

# Requirements

* Outcome-item GWAS z-scores for each clinical item
* Variant genomic positions
* A reference gene-expression matrix 'G' (individuals by genes)
* An individual-level genotype matrix 'V' measured in the same individuals as the reference gene-expression data (individuals by variants)
* A drug perturbation matrix of differential expression signatures `H^diff' (conditions by genes)
* Cell-line annotations to filter perturbation conditions by cell-line
* A negative anchor set 'A_-'

# Sample Run

