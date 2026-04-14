simulate_data <- function(
    SNPs,                 # n_indiv x p_snp genotype matrix or FBM-like object.
    # Rows are individuals and columns are SNPs from the source panel.
    chr_pos,              # data.frame with columns chr and pos giving genomic coordinates.
    # Rows must align with the columns of SNPs.
    block_id,             # LD-block label for each SNP column.
    # Used to choose distinct causal loci.
    col_i1 = NULL,        # Optional subset of eligible SNP columns.
    # The Mendelianization scan is restricted to SNPs sampled from these columns.
    p_use = 3000,         # Number of SNPs retained for the Mendelianization scan.
    L = 3,                # Number of causal loci.
    n_boot = 500,         # Size of the reference-expression cohort.
    n_gwas = 100000,      # Size of the GWAS cohort used to generate summary statistics.
    gwas_nsamp_range = c(50000, 100000),
    # Per-item GWAS sample-size range.
    # One sample size is drawn for each observed phenotype item.
    q_lat = 50,           # Number of latent phenotype dimensions before mixing into observed items.
    q_obs = 25,           # Number of observed phenotype items; equals ncol(Zstat).
    n_genes = 1500,       # Number of genes in the reference-expression and perturbation space.
    conf_dim = 3,         # Dimension of latent confounding shared by expression and phenotypes.
    nuis_dim = 4,         # Number of nuisance transcriptional factors in perturbational signatures.
    n_true_per_locus = 6, # Number of therapeutic perturbagens generated for each causal locus.
    n_anchor = 80,        # Number of negative-anchor perturbagens.
    n_null = 120,         # Number of null perturbagens.
    reps_per_drug = 4,    # Number of perturbation instances per perturbagen.
    cell_lines = c("neu", "npc"),
    # Retained cellular contexts for perturbation profiling.
    use_fbm_for_cmap = FALSE
    # Reserved argument; not used by the current implementation.
) {
  
  if (is.null(col_i1)) col_i1 <- seq_len(ncol(SNPs))
  stopifnot(length(block_id) == ncol(SNPs))
  
  gwas_nsamp_range <- sort(as.integer(gwas_nsamp_range))
  stopifnot(length(gwas_nsamp_range) == 2)
  stopifnot(gwas_nsamp_range[1] >= 1)
  stopifnot(gwas_nsamp_range[2] <= n_gwas)
  
  ## ---------- helpers ----------
  # Normalize a vector to unit Euclidean norm.
  unit_norm <- function(x) {
    s <- sqrt(sum(x^2))
    if (s == 0) x else x / s
  }
  
  # Construct a sparse unit-norm gene program on a specified support.
  make_sparse_vec <- function(p, support, values = NULL) {
    x <- numeric(p)
    if (is.null(values)) values <- rnorm(length(support))
    x[support] <- values
    unit_norm(x)
  }
  
  # Marginal SNP-trait z-like association score from standardized X and y.
  z_from_xy <- function(X, y) {
    # X: n x p, y: length n
    Xs <- scale(X)
    ys <- scale(y)[, 1]
    as.numeric(crossprod(ys, Xs)) / sqrt(length(ys) - 1)
  }
  
  # First principal-component score of an item matrix.
  # Used to create the total-score phenotype.
  first_pc_score <- function(Ymat) {
    Ystd <- scale(Ymat)
    pc <- prcomp(Ystd, center = FALSE, scale. = FALSE)
    score <- pc$x[, 1]
    loading <- pc$rotation[, 1]
    
    # Orient PC1 so larger values correspond to a larger mean standardized item score.
    avg_score <- rowMeans(Ystd)
    if (cor(score, avg_score) < 0) {
      score <- -score
      loading <- -loading
    }
    
    list(
      score = as.numeric(score),
      loading = as.numeric(loading)
    )
  }
  
  ## ---------- 1) pick SNP universe and causal variants ----------
  # Sample the SNP universe used in the Mendelianization scan until it spans
  # at least L distinct LD blocks.
  blocks_avail = c()
  while (length(blocks_avail) < L) {
    col_i <- sort(sample(col_i1, p_use, replace = FALSE))
    blocks_avail <- unique(block_id[col_i])
    blocks_avail <- blocks_avail[!is.na(blocks_avail)]
  }
  causal_blocks <- sample(blocks_avail, L, replace = FALSE)
  
  # Choose one causal variant from each causal locus.
  causal_variants <- sapply(causal_blocks, function(b) {
    sample(col_i[block_id[col_i] == b], 1)
  })
  
  ## Separate bootstrap cohorts:
  ## - boot_ref  for the reference genotype-expression sample
  ## - boot_gwas for the GWAS summary-statistics sample
  boot_ref  <- sample(seq_len(nrow(SNPs)), n_boot, replace = TRUE)
  boot_gwas <- sample(seq_len(nrow(SNPs)), n_gwas, replace = TRUE)
  
  V_causal_ref  <- as.matrix(SNPs[boot_ref,  causal_variants, drop = FALSE])
  V_causal_gwas <- as.matrix(SNPs[boot_gwas, causal_variants, drop = FALSE])
  
  # Standardize the causal genotype matrices.
  V_causal_ref  <- scale(V_causal_ref)
  V_causal_gwas <- scale(V_causal_gwas)
  
  ## ---------- 2) define locus-specific gene programs ----------
  # Gamma[, l] is the sparse locus-anchored gene program for locus l.
  Gamma <- matrix(0, nrow = n_genes, ncol = L)
  
  for (l in seq_len(L)) {
    supp <- sample(seq_len(n_genes), 80)
    vals <- rnorm(length(supp))
    Gamma[, l] <- make_sparse_vec(n_genes, supp, vals)
  }
  
  colnames(Gamma) <- paste0("locus", seq_len(L))
  
  # Latent confounding loadings for gene expression.
  A_g <- matrix(rnorm(conf_dim * n_genes, sd = 0.15), conf_dim, n_genes)
  
  ## ---------- 3) simulate GWAS-sample expression and phenotypes ----------
  # Latent confounders in the GWAS cohort.
  U_gwas <- matrix(rnorm(n_gwas * conf_dim), n_gwas, conf_dim)
  
  # Signed locus-to-latent-phenotype effects.
  Delta_lat <- matrix(
    sample(c(-1, 1), L * q_lat, replace = TRUE) * runif(L * q_lat, 0.01, 0.05),
    nrow = L, ncol = q_lat
  )
  
  # Choose the gene-to-latent map so that the genetic component of the latent
  # phenotypes is induced by the locus-specific gene programs.
  ridge <- 1e-8
  B_lat <- Gamma %*% solve(crossprod(Gamma) + diag(ridge, L), Delta_lat)
  
  # Gaussian latent noise with covariance induced by the latent mapping.
  Sigma_lat <- (0.7^2) * crossprod(B_lat)
  R_lat <- chol(Sigma_lat + diag(1e-10, q_lat))
  noise_lat_gwas <- matrix(rnorm(n_gwas * q_lat), n_gwas, q_lat) %*% R_lat
  
  # Latent phenotypes in the GWAS cohort:
  #   locus-driven signal
  # + latent confounding through expression
  # + Gaussian noise
  Y_lat_gwas <- V_causal_gwas %*% crossprod(Gamma, B_lat) +
    U_gwas %*% (A_g %*% B_lat) +
    noise_lat_gwas
  
  # Mixing from latent phenotypes to observed items, plus downstream confounding.
  A_y <- matrix(rnorm(q_lat * q_obs), q_lat, q_obs)
  C_y <- matrix(rnorm(conf_dim * q_obs, sd = 0.25), conf_dim, q_obs)
  
  # Observed phenotype items.
  Y_items <- Y_lat_gwas %*% A_y +
    U_gwas %*% C_y +
    matrix(rnorm(n_gwas * q_obs, sd = 1), n_gwas, q_obs)
  
  colnames(Y_items) <- paste0("item", seq_len(q_obs))
  
  ## ---------- 4) compute Mendelianization-style Z statistics ----------
  # Per-item GWAS sample sizes.
  nsamps <- sample(gwas_nsamp_range[1]:gwas_nsamp_range[2], q_obs, replace = TRUE)
  
  # SNP-by-item summary-association matrix used as input to Mendelianization.
  Zstat <- matrix(0, nrow = length(col_i), ncol = q_obs)
  colnames(Zstat) <- colnames(Y_items)
  
  blocks_use <- unique(block_id[col_i])
  blocks_use <- blocks_use[!is.na(blocks_use)]
  
  for (j in seq_len(q_obs)) {
    for (b in blocks_use) {
      idx_local <- which(block_id[col_i] == b)
      if (length(idx_local) == 0L) next
      
      # Independent GWAS subsample for this item-block pair.
      inds_jb <- sample(seq_len(n_gwas), nsamps[j], replace = FALSE)
      y_jb <- Y_items[inds_jb, j]
      
      X_jb <- as.matrix(SNPs[boot_gwas[inds_jb], col_i[idx_local], drop = FALSE])
      Zstat[idx_local, j] <- z_from_xy(X_jb, y_jb)
    }
  }
  
  ## ---------- 4b) total-score GWAS phenotype (PC1 on overlap across items) ----------
  # Total-score phenotype from the first principal component of the observed items
  # on an explicit overlap sample.
  n_overlap_gwas <- min(nsamps)
  inds_overlap_gwas <- sample(seq_len(n_gwas), n_overlap_gwas, replace = FALSE)
  
  pc_gwas <- first_pc_score(Y_items[inds_overlap_gwas, , drop = FALSE])
  Y_tot_gwas <- pc_gwas$score
  pc1_loadings_gwas <- pc_gwas$loading
  
  Zstat_tot <- matrix(0, nrow = length(col_i), ncol = 1)
  colnames(Zstat_tot) <- "Y_tot"
  
  for (b in blocks_use) {
    idx_local <- which(block_id[col_i] == b)
    if (length(idx_local) == 0L) next
    
    X_tot_b <- as.matrix(SNPs[boot_gwas[inds_overlap_gwas], col_i[idx_local], drop = FALSE])
    Zstat_tot[idx_local, 1] <- z_from_xy(X_tot_b, Y_tot_gwas)
  }
  
  chr_pos_use <- chr_pos[col_i, , drop = FALSE]
  rownames(chr_pos_use) <- NULL
  
  ## ---------- 5) simulate reference expression G and reference outcomes ----------
  # Latent confounders in the reference-expression cohort.
  U_ref <- matrix(rnorm(n_boot * conf_dim), n_boot, conf_dim)
  
  # Reference expression matrix:
  #   locus-driven signal
  # + latent confounding
  # + Gaussian noise
  G_ref_raw <- V_causal_ref %*% t(Gamma) +
    U_ref %*% A_g +
    matrix(rnorm(n_boot * n_genes, sd = 0.7), n_boot, n_genes)
  
  # Standardized reference expression matrix used by SIEVE.
  G <- scale(G_ref_raw)
  colnames(G) <- paste0("gene", seq_len(n_genes))
  
  # Latent noise term for the reference cohort.
  noise_lat_ref <- matrix(rnorm(n_boot * q_lat), n_boot, q_lat) %*% R_lat
  
  # Reference latent phenotypes and observed items.
  Y_lat_ref <- V_causal_ref %*% crossprod(Gamma, B_lat) +
    U_ref %*% (A_g %*% B_lat) +
    noise_lat_ref
  Y_items_ref <- Y_lat_ref %*% A_y +
    U_ref %*% C_y +
    matrix(rnorm(n_boot * q_obs, sd = 1), n_boot, q_obs)
  
  colnames(Y_items_ref) <- paste0("item", seq_len(q_obs))
  
  ## ---------- 5b) total-score reference phenotype (PC1) ----------
  pc_ref <- first_pc_score(Y_items_ref)
  Y_tot_ref <- pc_ref$score
  pc1_loadings_ref <- pc_ref$loading
  
  ## ---------- 6) simulate CMAP differential perturbation signatures ----------
  # Genes touched by at least one locus-specific program.
  gamma_support <- which(rowSums(abs(Gamma)) > 0)
  
  # Build nuisance factors mostly outside the causal-program support.
  nuis_pool <- setdiff(seq_len(n_genes), gamma_support)
  if (length(nuis_pool) < 50) nuis_pool <- seq_len(n_genes)
  
  Nuis <- matrix(0, nrow = n_genes, ncol = nuis_dim)
  for (r in seq_len(nuis_dim)) {
    supp <- sample(nuis_pool, min(80, length(nuis_pool)))
    Nuis[, r] <- make_sparse_vec(n_genes, supp, rnorm(length(supp)))
  }
  
  # Context-specific effects for cell line, treatment time, and dose.
  n_cell <- length(cell_lines)
  cell_eff <- matrix(rnorm(n_genes * n_cell, sd = 0.15), n_genes, n_cell)
  time_levels <- c("6h", "24h")
  dose_levels <- c("low", "high")
  time_eff <- matrix(rnorm(n_genes * length(time_levels), sd = 0.10), n_genes, length(time_levels))
  dose_eff <- matrix(rnorm(n_genes * length(dose_levels), sd = 0.10), n_genes, length(dose_levels))
  
  # Perturbagen identities for the three simulated classes.
  true_ids <- unlist(lapply(seq_len(L), function(l) {
    paste0("TRUE_L", l, "_", seq_len(n_true_per_locus))
  }))
  anchor_ids <- paste0("ANCH_", seq_len(n_anchor))
  null_ids <- paste0("NULL_", seq_len(n_null))
  
  all_pert_ids <- c(true_ids, anchor_ids, null_ids)
  
  drug_type <- c(
    rep("true", length(true_ids)),
    rep("anchor", length(anchor_ids)),
    rep("null", length(null_ids))
  )
  
  # Locus assignment for each perturbagen.
  locus_of_drug <- c(
    rep(seq_len(L), each = n_true_per_locus),
    rep(NA_integer_, length(anchor_ids)),
    sample(seq_len(L), length(null_ids), replace = TRUE)
  )
  
  # cond_meta stores the perturbation-condition annotations.
  # Each column of H corresponds to one drug-dose-time-cell-line instance.
  cond_meta <- vector("list", length(all_pert_ids) * reps_per_drug)
  H <- matrix(0, nrow = n_genes, ncol = length(all_pert_ids) * reps_per_drug)
  
  idx <- 1L
  truth <- data.frame(
    pert_id = all_pert_ids,
    drug_type = drug_type,
    locus = locus_of_drug,
    is_true_therapeutic = (drug_type == "true"),
    stringsAsFactors = FALSE
  )
  
  for (d in seq_along(all_pert_ids)) {
    pert_id <- all_pert_ids[d]
    d_type <- drug_type[d]
    l <- locus_of_drug[d]
    
    if (d_type == "true") {
      # Therapeutic perturbagens align with a locus-specific mechanism plus weak nuisance.
      a <- sample(c(-1, 1), 1) * runif(1, 0.8, 1.2)
      h_base <- a * Gamma[, l] + Nuis %*% rnorm(nuis_dim, sd = 0.15)
    } else if (d_type == "anchor") {
      # Negative anchors are dominated by nuisance transcriptional structure.
      h_base <- Nuis %*% rnorm(nuis_dim, sd = 1.0)
    } else {
      # Null perturbagens may weakly overlap a locus-specific mechanism but are
      # more nuisance-driven and less coherent.
      a <- rnorm(1, sd = 0.15)
      h_base <- a * Gamma[, l] + Nuis %*% rnorm(nuis_dim, sd = 0.35)
    }
    
    for (rep_i in seq_len(reps_per_drug)) {
      cell_i <- sample(seq_len(n_cell), 1)
      time_i <- sample(seq_along(time_levels), 1)
      dose_i <- sample(seq_along(dose_levels), 1)
      
      # Condition-level differential perturbation signature:
      #   base perturbagen effect
      # + cell-line effect
      # + treatment-time effect
      # + dose effect
      # + residual noise
      h <- h_base +
        cell_eff[, cell_i] +
        time_eff[, time_i] +
        dose_eff[, dose_i] +
        rnorm(n_genes, sd = 0.35)
      
      H[, idx] <- h
      cond_meta[[idx]] <- data.frame(
        cond_id = idx,
        pert_id = pert_id,
        cmap_name = pert_id,
        cell_iname = cell_lines[cell_i],
        pert_itime = time_levels[time_i],
        pert_idose = dose_levels[dose_i],
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
  }
  
  cond_meta <- do.call(rbind, cond_meta)
  rownames(cond_meta) <- NULL
  
  rownames(H) <- colnames(G)
  colnames(H) <- paste0("cond_", seq_len(ncol(H)))
  
  # Differential perturbation matrix used by SIEVE.
  CMAP_diff <- H
  
  ## ---------- outputs ----------
  list(
    Zstat = Zstat,                     # SNP-by-item summary-association matrix
    Zstat_tot = Zstat_tot,             # SNP-by-total-score summary-association matrix
    chr_pos = chr_pos_use,             # coordinates for the retained scan SNPs
    iSNPs = boot_gwas,                 # legacy GWAS row index output
    iSNPs_ref = boot_ref,              # reference-cohort row indices
    iSNPs_gwas = boot_gwas,            # GWAS-cohort row indices
    jSNPs = col_i,                     # SNP columns retained for the scan
    nsamps = nsamps,                   # per-item GWAS sample sizes
    n_overlap_gwas = n_overlap_gwas,   # overlap size used for total-score PC1
    inds_overlap_gwas = inds_overlap_gwas,
    Y_tot_gwas = Y_tot_gwas,           # GWAS total-score phenotype
    Y_tot_ref = Y_tot_ref,             # reference total-score phenotype
    pc1_loadings_gwas = pc1_loadings_gwas,
    pc1_loadings_ref = pc1_loadings_ref,
    G = G,                             # reference-expression matrix
    CMAP_diff = CMAP_diff,             # differential perturbation matrix
    iCMAP = seq_len(n_genes),          # identity mapping into perturbation gene space
    cond_meta = cond_meta,             # perturbation-condition metadata
    anchor_neg_pert_ids = anchor_ids,  # negative-anchor perturbagen IDs
    causal_variants = causal_variants, # causal SNPs
    causal_blocks = causal_blocks,     # causal LD blocks
    Y_items = Y_items,                 # GWAS observed phenotype items
    Gamma = Gamma,                     # locus-specific gene programs
    truth = truth,                     # perturbagen truth labels
    Y_lat_gwas = Y_lat_gwas,           # GWAS latent phenotypes
    Y_items_gwas = Y_items,            # alias of GWAS observed phenotype items
    Y_lat_ref = Y_lat_ref,             # reference latent phenotypes
    Y_items_ref = Y_items_ref          # reference observed phenotype items
  )
}
