SIEVE <- function(
    Zstat, chr_pos, SNPs, iSNPs, jSNPs,
    G, CMAP_diff, iCMAP,
    cond_meta, cell_lines=c("neu","npc"),
    # --- anchors ---
    anchor_neg_cmap_names = NULL,
    anchor_neg_pert_ids   = NULL,
    anchor_r_neg = 60,
    # --- NES / GSEA ---
    nmin = 1,
    nperm = 10000,
    gsea_param = 1
) {
  require(bigstatsr)
  require(data.table)
  require(fgsea)
  
  # ===========================================================================
  # Inputs
  # ------
  # Zstat : matrix
  #   Outcome-item GWAS summary association matrix used by Mendelianization.
  #   Rows correspond to variants; columns correspond to phenotype items; entries
  #   are typically z-statistics.
  #
  # chr_pos : data.frame
  #   Must have two columns: chromosome and base-pair position, where rows match 
  #   those of Zstat.
  #
  # SNPs : file-backed matrix (FBM, see bigstatsr)
  #   Individual-level genotype matrix. Rows are individuals and columns are
  #   variants available for the genotype-expression calibration step.
  #
  # iSNPs : integer vector
  #   Row indices selecting the individuals in SNPs used in this analysis.
  #
  # jSNPs : integer vector
  #   Column map from the Mendelianization / GWAS variant ordering into the SNPs
  #   matrix columns.
  #
  # G : numeric matrix
  #   Reference-expression matrix for the same individuals used in the genotype
  #   matrix. Rows are individuals; columns are genes.
  #
  # CMAP_diff : file-backed matrix (FBM)
  #   Drug-versus-control differential perturbation matrix H^{diff}. Rows are
  #   perturbation-space genes and columns are perturbation conditions
  #   (drug-dose-time-cell-line instances).
  #
  # iCMAP : integer vector
  #   Mapping from the genes in G to rows of CMAP_diff. NA indicates that a gene
  #   in G is absent from the perturbation gene universe.
  #
  # cond_meta : data.frame
  #   Metadata aligned to columns of CMAP_diff. Expected fields include cond_id,
  #   pert_id, cmap_name, cell_iname, pert_itime, and pert_idose.
  #
  # cell_lines : character vector
  #   Retained cell lines. The function restricts scoring and FGSEA aggregation to
  #   perturbation conditions from these cell lines.
  #
  # anchor_neg_cmap_names : character vector or NULL
  #   Optional negative-anchor drug names. Matching conditions define a nuisance
  #   expression subspace to project away.
  #
  # anchor_neg_pert_ids : character vector or NULL
  #   Optional negative-anchor perturbagen IDs, used alongside or instead of names.
  #
  # anchor_r_neg : integer
  #   Rank of the negative-anchor basis U_- built from anchor perturbations 
  #   (default = 60).
  #
  # nmin : integer
  #   Minimum number of retained perturbation instances required for a drug to
  #   enter FGSEA aggregation (default = 1).
  #
  # nperm : integer
  #   Number of FGSEA permutations (default = 1000).
  #
  # gsea_param : numeric
  #   GSEA weighting exponent passed to fgsea::fgsea / fgsea::fgseaSimple (default = 1).
  #
  # Returns
  # -------
  # A list containing:
  #   drug_ranks_dt            : wide-format locus-specific drug ranking table
  #   drug_nes_long            : long-format FGSEA results across loci
  #   deltaLs                  : condition-level cosine-similarity scores by locus
  #   lead_info                : metadata for Mendelianization lead variants
  #   cond_cols_used           : retained perturbation-condition columns
  #   anchor_neg_cols          : retained negative-anchor condition columns
  #   anchor_r_neg_used        : realized negative-anchor basis rank
  #   cmap_rows_used           : perturbation-space genes used in scoring
  #   excluded_anchor_pert_ids : anchors removed from final candidate outputs
  # ===========================================================================
  
  # Standardize and order perturbation metadata so that metadata rows align with
  # columns of the perturbation matrix.
  cm <- as.data.table(cond_meta)
  setorder(cm, cond_id)
  
  # Normalize cell-line labels and add explicit column indices.
  cm[, `:=`(
    cell_iname2 = tolower(trimws(cell_iname)),
    col_idx     = .I
  )]
  
  # Retain only the perturbation conditions from the requested cell lines.
  cm_cns <- cm[cell_iname2 %chin% cell_lines]
  
  # Column indices of the retained perturbation conditions.
  cond_idx <- sort(unique(cm_cns$col_idx))
  n_cond <- length(cond_idx)
  
  # Precompute Euclidean norms of the retained perturbation signatures. These are
  # later adjusted if the negative-anchor projection is active.
  cs <- big_colstats(CMAP_diff, ind.row = iCMAP, ind.col = cond_idx)  # data.frame with sum, var 
  p  <- length(iCMAP)
  
  sumsq <- (p - 1) * cs$var + (cs$sum^2) / p
  cond_norm_sub <- sqrt(pmax(sumsq, 0))   # guard tiny negatives
  
  
  # ---- helpers ----
  
  # First non-missing / non-empty string, used when collapsing annotations to one
  # name per perturbagen.
  first_nonempty <- function(x) {
    x <- as.character(x)
    x <- x[!is.na(x) & nzchar(x)]
    if (length(x)) x[1] else NA_character_
  }
  
  # big_cprodMat() can return either condition-by-basis or basis-by-condition;
  # reorient to condition-by-basis for downstream row-wise calculations.
  as_cond_by_basis <- function(M, n_cond) {
    M <- as.matrix(M)
    if (nrow(M) == n_cond) return(M)
    if (ncol(M) == n_cond) return(t(M))
    stop("Unexpected dimensions from big_cprodMat().")
  }
  
  # Orthogonally project a vector away from the span of U.
  proj_out  <- function(U, x) x - drop(U %*% crossprod(U, x))
  
  # Resolve retained negative-anchor conditions from metadata by drug name and/or
  # perturbagen ID.
  get_anchor_cols <- function(dt_base, cmap_names = NULL, pert_ids = NULL) {
    keep <- rep(FALSE, nrow(dt_base))
    if (!is.null(cmap_names) && length(cmap_names)) {
      keep <- keep | (tolower(trimws(dt_base$cmap_name)) %in% tolower(trimws(cmap_names)))
    }
    if (!is.null(pert_ids) && length(pert_ids)) {
      keep <- keep | (dt_base$pert_id %chin% as.character(pert_ids))
    }
    unique(dt_base[keep, cond_col])
  }
  
  # Build the negative-anchor basis U_- from anchor perturbations using a
  # truncated SVD on the retained perturbation gene universe.
  build_basis <- function(ind.col, r, gene_rows_all) {
    if (is.null(ind.col)) return(NULL)
    k_cap <- min(length(ind.col), length(gene_rows_all)) - 1L
    if (!is.finite(k_cap) || k_cap < 1L) return(NULL)
    k_use <- min(as.integer(r), as.integer(k_cap))
    sv <- bigstatsr::big_randomSVD(
      CMAP_diff,
      ind.row = gene_rows_all,
      ind.col = ind.col,
      k = k_use
    )
    if (is.null(sv$u)) stop("big_randomSVD did not return $u; check bigstatsr.")
    sv$u[, seq_len(min(k_use, ncol(sv$u))), drop = FALSE]
  }
  
  # Aggregate condition-level cosine-similarity scores to the perturbagen level. Each
  # drug is treated as the set of its retained condition instances and evaluated
  # by FGSEA over the ranked condition list.
  run_fgsea_signed <- function(dt, nperm, gsea_param) {
    # dt must have: inst_id, pert_id, score
    dt <- dt[is.finite(score)]
    if (!nrow(dt)) return(data.table())
    
    # Ranked condition-level statistics used by FGSEA.
    stats <- dt$score
    names(stats) <- dt$inst_id
    
    # Break exact ties with tiny jitter to stabilize ranking behavior.
    stats <- stats + 1e-12 * stats::rnorm(length(stats))
    stats <- sort(stats, decreasing = TRUE)
    
    # One pathway per perturbagen, defined by its retained condition instances.
    pathways <- split(dt$inst_id, dt$pert_id)
    pathways <- pathways[lengths(pathways) >= nmin]
    if (!length(pathways)) return(data.table())
    
    if ("fgseaSimple" %in% getNamespaceExports("fgsea")) {
      out <- fgsea::fgseaSimple(
        pathways = pathways,
        stats = stats,
        nperm = nperm,
        gseaParam = gsea_param
      )
    } else {
      out <- fgsea::fgsea(
        pathways = pathways,
        stats = stats,
        nperm = nperm,
        gseaParam = gsea_param
      )
    }
    
    out <- as.data.table(out)
    if (!nrow(out)) return(data.table())
    
    setnames(out, "pathway", "pert_id")
    out[]
  }
  
  # ---- Mendelianization ----
  
  # Learn composite outcomes / subphenotypes and record their lead variants.
  mod <- Mendelianization_SIEVE(Zstat, chr_pos = chr_pos)
  ord <- mod$ord
  
  chr_pos2 <- mod$chr_pos
  names(chr_pos2)[1:2] <- c("chr", "pos")
  
  lead_info <- data.table(
    locus    = seq_along(mod$leads),
    lead_idx = mod$leads,
    chr      = chr_pos2$chr[mod$leads],
    pos      = chr_pos2$pos[mod$leads]
  )
  
  # ---- perturbation gene universe ----
  
  # Retain the perturbation-space rows that are mapped from the reference
  # expression genes.
  gene_rows_all <- sort(unique(iCMAP[!is.na(iCMAP)]))
  if (!length(gene_rows_all)) stop("No non-NA rows found in iCMAP.")
  
  # ---- condition metadata for retained perturbations ----
  
  dt_base <- cm[cond_idx, .(pert_id, cmap_name, cell_iname, pert_itime, pert_idose)]
  dt_base[, `:=`(
    pert_id = as.character(pert_id),
    cmap_name = as.character(cmap_name)
  )]
  dt_base[, cond_col := cond_idx]
  
  # Assign a unique ID to each retained perturbation condition for FGSEA.
  dt_base[, inst_id := paste0("inst_", seq_len(.N))]
  
  # ---- negative-anchor precomputation ----
  
  # Identify retained perturbation conditions belonging to the negative-anchor set.
  neg_cols <- get_anchor_cols(dt_base, anchor_neg_cmap_names, anchor_neg_pert_ids)
  
  # Build the negative-anchor basis once and reuse it across all learned outcomes.
  U_neg <- build_basis(neg_cols, anchor_r_neg, gene_rows_all)
  
  # Precompute condition coordinates in the negative-anchor basis. These are used
  # to update the denominator after projection of perturbation signatures.
  H_neg_coef <- NULL
  cond_norm_use <- cond_norm_sub
  
  if (!is.null(U_neg)) {
    H_neg_coef <- as_cond_by_basis(
      bigstatsr::big_cprodMat(
        CMAP_diff,
        U_neg,
        ind.row = gene_rows_all,
        ind.col = cond_idx
      ),
      n_cond = n_cond
    )
  }
  
  # Replace raw perturbation norms with the norms of the projected signatures:
  # || (I - U_- U_-^T) H_{:,c} ||_2
  if (!is.null(U_neg)) {
    cond_norm_use <- sqrt(pmax(cond_norm_sub^2 - rowSums(H_neg_coef^2), 0))
  }
  
  # ---- exclude anchor perturbagens from final candidate outputs ----
  
  # Anchors are used only to define the nuisance subspace and are removed from the
  # downstream ranked candidate set.
  anchor_names_all <- unique(tolower(trimws(anchor_neg_cmap_names)))
  anchor_ids_all   <- unique(as.character(anchor_neg_pert_ids))
  
  exclude_pert_ids <- character(0)
  if (length(anchor_ids_all))  exclude_pert_ids <- c(exclude_pert_ids, anchor_ids_all)
  if (length(neg_cols)) exclude_pert_ids <- c(exclude_pert_ids, dt_base[cond_col %in% neg_cols, unique(pert_id)])
  if (length(anchor_names_all)) exclude_pert_ids <- c(
    exclude_pert_ids,
    dt_base[tolower(trimws(cmap_name)) %chin% anchor_names_all, unique(pert_id)]
  )
  exclude_pert_ids <- unique(exclude_pert_ids[!is.na(exclude_pert_ids) & nzchar(exclude_pert_ids)])
  
  # ---- gene-wise scaling of reference expression ----
  
  # Scale reference-expression genes before computing Cov(V, G/s). A lower-tail
  # floor prevents division by zero or near-zero standard deviations.
  sd_G <- apply(G, 2, sd)
  sd_G[sd_G == 0] <- 1
  sd_floor <- as.numeric(quantile(sd_G, 0.01, na.rm = TRUE))
  
  sJ_all <- pmax(sd_G, sd_floor)
  G_scaled_all <- sweep(G, 2, sJ_all, "/")
  
  # ---- map reference genes into perturbation space ----
  
  ok_all <- !is.na(iCMAP)
  r_ok_all <- iCMAP[ok_all]
  if (!any(ok_all)) stop("No non-NA gene mapping in iCMAP for supplied G columns.")
  
  # ---- per-locus condition-level cosine-similarity scores ----
  
  # One column per learned outcome / locus-anchored subphenotype.
  L0 <- length(mod$leads)
  deltaLs <- matrix(NA_real_, nrow = n_cond, ncol = L0)
  
  for (i in seq_along(mod$leads)) {
    # Select the lead variants from loci associated with learned outcome i.
    V <- find_leads(chr_pos2, mod$pZ[, i])
    V_cols <- jSNPs[ord[V]]
    Vmat <- as.matrix(SNPs[iSNPs, V_cols, drop = FALSE])
    
    # Genetically driven embedding:
    # C_i = Cov(V_i, G/s), then SVD. The right singular vectors define the
    # outcome-specific mechanism space.
    C  <- cov(Vmat, G_scaled_all)
    sv <- svd(C)
    k  <- min(ncol(Vmat), ncol(G_scaled_all))
    Wk <- sv$v[, 1:k, drop = FALSE]
    M  <- G_scaled_all %*% Wk
    
    # Ridge-stabilized calibration toward the Mendelianized target using the
    # local lead-variant association pattern z_L.
    Vstd <- scale(Vmat, center = TRUE, scale = TRUE)
    Mctr <- scale(M,    center = TRUE, scale = TRUE)
    B <- cov(Vstd, Mctr)
    
    z_L <- mod$Z[V, i, drop = FALSE]
    BtB <- t(B) %*% B
    theta <- solve(BtB + diag(1e-8, nrow(BtB)), t(B) %*% z_L)
    
    # Map the calibrated direction back to gene space and undo the earlier
    # gene-wise scaling.
    w_gene <- drop(Wk %*% theta) / sJ_all
    
    # Align the mechanism vector to the perturbation-gene universe. If multiple
    # reference-expression genes map to the same perturbation row, sum them.
    w_ok <- w_gene[ok_all]
    if (anyDuplicated(r_ok_all)) {
      agg <- tapply(w_ok, r_ok_all, sum)
      rows_u <- as.integer(names(agg))
      w_u    <- as.numeric(agg)
    } else {
      rows_u <- r_ok_all
      w_u    <- w_ok
    }
    
    w_full <- numeric(length(gene_rows_all))
    pos_u  <- match(rows_u, gene_rows_all)
    pos_ok <- !is.na(pos_u)
    w_full[pos_u[pos_ok]] <- w_u[pos_ok]
    
    # Anchor-guided projection: remove the component of the mechanism vector that
    # lies in the negative-anchor nuisance subspace.
    w_use <- w_full
    if (!is.null(U_neg)) w_use <- proj_out(U_neg, w_use)
    
    # Normalize the refined mechanism direction before condition-level scoring.
    nrm <- sqrt(sum(w_use^2))
    if (!is.finite(nrm) || nrm < 1e-10) next
    w_use <- w_use / (nrm + 1e-12)
    
    # Condition-level cosine-similarity score:
    # numerator = inner product between the perturbation signature and the
    # refined mechanism vector;
    # denominator = norm of the projected perturbation signature.
    sc <- drop(bigstatsr::big_cprodMat(
      CMAP_diff,
      matrix(w_use, ncol = 1),
      ind.row = gene_rows_all,
      ind.col = cond_idx
    ))
    
    # Store the per-condition statistics that will later be ranked and passed to
    # FGSEA for perturbagen-level aggregation.
    deltaLs[, i] <- sc / (cond_norm_use + 1e-8)
  }
  
  
  # ---- drug-level aggregation via FGSEA ----
  
  nes_list <- vector("list", L0)
  
  for (l in seq_len(L0)) {
    sc <- deltaLs[, l]
    if (all(is.na(sc))) next
    
    # Start from the retained condition instances for learned outcome l.
    dt <- copy(dt_base)
    dt[, score := sc]
    dt <- dt[is.finite(score)]
    
    # Remove perturbagens used as negative anchors.
    if (length(exclude_pert_ids)) {
      dt <- dt[!(pert_id %chin% exclude_pert_ids)]
    }
    if (!nrow(dt)) next
    
    # Require at least nmin retained instances before drug-level aggregation.
    counts <- dt[, .(
      n_cond  = .N,
      n_cells = uniqueN(cell_iname)
    ), by = pert_id]
    
    keep_ids <- counts[n_cond >= nmin, pert_id]
    dt <- dt[pert_id %chin% keep_ids]
    if (!nrow(dt)) next
    
    # Test whether each drug's condition instances are enriched near one end of
    # the ranked condition list.
    fg <- run_fgsea_signed(dt, nperm = nperm, gsea_param = gsea_param)
    if (!nrow(fg)) next
    
    # Attach annotations and retained-instance counts.
    ann <- dt[, .(cmap_name = first_nonempty(cmap_name)), by = pert_id]
    
    out <- merge(fg, ann, by = "pert_id", all.x = TRUE)
    out <- merge(out, counts, by = "pert_id", all.x = TRUE)
    
    # Final within-locus ranking.
    out[, abs_NES := abs(NES)]
    setorder(out, pval, padj, -abs_NES, pert_id)
    out[, `:=`(rank = .I, locus = l)]
    
    nes_list[[l]] <- out[, .(
      rank, locus, pert_id, cmap_name,
      NES, abs_NES, pval, padj, size, n_cond, n_cells
    )]
  }
  
  # Stack locus-specific FGSEA results into one long table.
  drug_nes_long <- rbindlist(nes_list, use.names = TRUE, fill = TRUE)
  
  # Also build a wide-format ranking table for easier inspection of top hits
  # across learned outcomes.
  drug_ranks_dt <- NULL
  if (nrow(drug_nes_long)) {
    drug_ranks_dt <- dcast(
      drug_nes_long, rank ~ locus,
      value.var = c("cmap_name", "pert_id", "NES", "abs_NES", "pval", "padj", "size", "n_cond", "n_cells")
    )
    setnames(
      drug_ranks_dt,
      old = setdiff(names(drug_ranks_dt), "rank"),
      new = sub("_(\\d+)$", "_rank_\\1", setdiff(names(drug_ranks_dt), "rank"))
    )
    setkey(drug_ranks_dt, rank)
  }
  
  list(
    drug_ranks_dt  = drug_ranks_dt,
    drug_nes_long  = drug_nes_long,
    deltaLs        = deltaLs,  # per-condition cosine-similarity scores used as FGSEA stats
    lead_info      = lead_info,
    cond_cols_used = cond_idx,
    anchor_neg_cols = neg_cols,
    anchor_r_neg_used = if (is.null(U_neg)) 0L else ncol(U_neg),
    cmap_rows_used = gene_rows_all,
    excluded_anchor_pert_ids = exclude_pert_ids
  )
}
