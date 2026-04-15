# SIEVE: Locus-Anchored Drug Prioritization for Complex Disorders

SIEVE prioritizes drugs for complex disorders by first decomposing broad clinical phenotypes into locus-anchored subphenotypes, then using genotype-coupled reference expression to infer a genetically calibrated mechanism vector for each subphenotype. It projects both mechanism vectors and perturbational signatures away from negative-control expression programs and aggregates concordance across experimental contexts to produce robust, locus-specific drug rankings.

# Installation

> library(devtools)

> install_github("ericstrobl/SIEVE")

> library(SIEVE)

# Requirements

* Outcome-item GWAS z-scores for each clinical item
* Variant genomic positions

* A reference gene-expression matrix (individuals by genes)
* An individual-level genotype matrix measured in the same individuals as the reference gene-expression data (individuals by variants)
  
* A drug perturbation matrix of differential expression signatures (conditions by genes)
* Cell-line annotations to filter perturbation conditions by cell-line
  
* (Optional) A negative anchor set

# Sample Run

First download simulation_inputs_small.rds from the Data folder. Then load the compact semi-synthetic 1000 Genomes input bundle:

> dat <- readRDS("simulation_inputs_small.rds")

Generate one synthetic dataset. This step can take around 2 hours:

> data = simulate_data(dat$SNPs, dat$chr_pos, dat$block_id)

Run SIEVE:

> SIEVE_mod = SIEVE(data$Zstat, data$chr_pos, dat$SNPs, data$iSNPs_ref, data$jSNPs, data$G, as_FBM(data$CMAP_diff), data$iCMAP, data$cond_meta, anchor_neg_cmap_names = data$anchor_neg_pert_ids)

Display the top-ranked perturbagens for each locus. The highest-ranked hits within a locus should often include many `TRUE_Lk_*` perturbagens sharing the same `k`, indicating enrichment of a common true therapeutic class within that locus:

> print(SIEVE_mod$drug_ranks_dt[1:20,2:4])
