#' Add meta data to the `miRTS_score` output
#'
#' @param miR_TS.output the output from miRTS_score()
#' @param meta_data the meta data for the `Input_df` to miRTS_score();
#' should include:
#' Col_names = Sample names used in the `Input_df` (e.g., 'S001', 'S002')
#'sample_type = grouping of the sample (e.g.,'control', 'disease')
#'
#' @return Matrix: samples (rows) ×  tissues+sample_type (columns).
#' @export
#'
#' @examples
#' dim(miR_TS.output.bkup$proportions); dim(hepatitis_C.meta)
#' hepatitis_C.output <- miRTS::Add_metadata(miR_TS.output.bkup$proportions, hepatitis_C.meta)
#' dim(hepatitis_C.output)

Add_metadata <- function(miR_TS.output, meta_data) {

  req <- c("Col_names", "sample_type")
  if (!all(req %in% names(meta_data))) {
    stop("meta_data must contain columns: Col_names, sample_type")
  }

  # proportions -> data.frame with Col_names
  df <- as.data.frame(miR_TS.output$proportions, check.names = FALSE)
  df$Col_names <- rownames(df)
  rownames(df) <- NULL

  # keep only needed meta cols + inner join via merge()
  md <- meta_data[, req, drop = FALSE]
  md$Col_names <- as.character(md$Col_names)

  out <- merge(df, md, by = "Col_names", all = FALSE, sort = FALSE)

  # restore rownames + add miR_using
  rownames(out) <- out$Col_names
  out$Col_names <- NULL
  out$miR_using <- nrow(miR_TS.output$mix)

  print(table(out$sample_type))
  out
}

#' Compute MCP-style tissue scores from miRNA counts
#'
#' @param expr_tpm   Matrix or data.frame of miRNA expression (e.g., TPM),
#'                   rows = miRNAs, cols = samples.
#' @param ts_df      Data frame with columns:
#'                     - miRNA: miRNA ID matching rownames(expr_tpm)
#'                     - enriched_organ: tissue/organ name
#'                   Default to ts_miRNA, which records the enriched tissues for each miRNAs in the `miRTS_signature_v1`
#' @param log2_transform  Logical; log2(TPM + pseudo_count) if TRUE.
#' @param pseudo_count    Pseudo-count added before log2 transform.
#' @param center_genes    Logical; median-center each miRNA across samples.
#' @param min_markers     Minimum number of markers required for a tissue.
#'                        If fewer, scores are set to NA for that tissue.
#' @param summary_fun     "median" (default) or "mean" for marker aggregation.
#'
#' @return Matrix: tissues (rows) × samples (columns) of MCP-style scores.
compute_mcp_mirna_scores <- function(expr_tpm,
                                     ts_df,
                                     log2_transform = TRUE,
                                     pseudo_count   = 1,
                                     center_genes   = TRUE,
                                     min_markers    = 2,
                                     summary_fun    = c("median", "mean")) {

  summary_fun <- match.arg(summary_fun)

  # --- checks --------------------------------------------------------------
  expr_mat <- as.matrix(expr_tpm)

  if (is.null(rownames(expr_mat))) {
    stop("expr_tpm must have rownames = miRNA IDs.")
  }
  if (!all(c("miRNA", "enriched_organ") %in% colnames(ts_df))) {
    stop('ts_df must have columns "miRNA" and "enriched_organ".')
  }

  ts_df$miRNA         <- as.character(ts_df$miRNA)
  ts_df$enriched_organ <- as.character(ts_df$enriched_organ)

  # --- build marker list: one vector of miRNAs per tissue ------------------
  marker_list <- split(ts_df$miRNA, ts_df$enriched_organ)

  # --- transform expression: log2(TPM + pseudo) ----------------------------
  if (log2_transform) {
    expr_mat <- log2(expr_mat + pseudo_count)
  }

  # --- optional median centering per miRNA (gene) --------------------------
  if (center_genes) {
    # subtract row median from each row (miRNA)
    row_meds <- apply(expr_mat, 1, stats::median, na.rm = TRUE)
    expr_mat <- expr_mat - row_meds
  }

  # --- compute scores for each tissue --------------------------------------
  tissues <- names(marker_list)
  n_samp  <- ncol(expr_mat)

  score_mat <- matrix(NA_real_, nrow = length(tissues), ncol = n_samp,
                      dimnames = list(tissues, colnames(expr_mat)))

  for (i in seq_along(tissues)) {
    tissue  <- tissues[i]
    markers <- unique(marker_list[[tissue]])

    # keep only markers present in expression matrix
    markers <- intersect(markers, rownames(expr_mat))

    if (length(markers) < min_markers) {
      warning(sprintf(
        "Tissue '%s' has only %d markers in expr_tpm (< %d); scores set to NA.",
        tissue, length(markers), min_markers
      ))
      next
    }

    m_sub <- expr_mat[markers, , drop = FALSE]

    if (summary_fun == "median") {
      score_mat[tissue, ] <- apply(m_sub, 2, stats::median, na.rm = TRUE)
    } else {
      score_mat[tissue, ] <- colMeans(m_sub, na.rm = TRUE)
    }
  }

  return(score_mat)
}

# 3. Key function of miR-TS: ####

#' @title miRNA-based Tissue Signal (miR-TS) scores
#' @description Compute miRNA-based Tissue Signal (miR-TS) scores from extracellular miRNAs in circulation
#'
#' @param Input_df Input miRNAs data with miRNA names on the row
#'  and sample names on the column.
#'  Read counts (or similar on an exponential scale) should be used.
#' @param signature_matrix matrix: miRNAs x tissues,
#'  defaults to the signature matrix optimized from miRNATissueAtlas2 (PMID: 34570238)
#' @param method character: deconvolution method ("cibersort" default)
#' @param Dtct_cutoff Detect rate cutoff to include a subset of
#' miRNAs to be considered as reliably measured, default is 10%;
#'  note that a miRNA with read count being 0 still inform the miR-TS,
#'  so it is advised to keep this number low (10% or 20%) to include
#'  as many signature miRNAs as possible.
#' @param outlier_sample.cutoff Cutoff of the MAD folds of the library sizes
#' to define outlier samples. Default to NA.
#' @param useAbsolute Use Absolute mode or not. Default is TRUE.
#' @param ... additional arguments passed to the deconvolution engine
#'
#' @param Norm Method for normalizing `Input_df`, default to TMM.
#' @param ts_miR.df a parameter required when using the MCP-counter method;
#' Default to ts_miRNA, which records the enriched tissues for each miRNAs in the `miRTS_signature_v1`
#' @param Use_blood Use plasma or serum samples or not. Default is TRUE.
#'
#' @return A list, including: proportions (the estimated tissue scores),
#'   mix (the `Input_df`), signatures (the `signature_matrix`)
#' @export
#'
#' @examples
#' miR_TS.output <- miRTS_score(Input_df = example_counts)
#' all.equal(miR_TS.output, miR_TS.output.bkup) # miR_TS.output.bkup is the ultimate output saved in the package for testing.
#' dim(miR_TS.output$proportions)

miRTS_score <- function(
    Input_df,
    signature_matrix = miRTS::miRTS_signature_v1,
    method = c("cibersort", "xCell2", "MCP-counter"),
    Dtct_cutoff = 0.1,
    Norm = c("TMM", "RLE", "none"),
    ts_miR.df = ts_miRNA,
    outlier_sample.cutoff = NA,
    useAbsolute = TRUE,
    Use_blood = TRUE,
    ...
) {
  method <- match.arg(method)
  Norm <- match.arg(Norm)
  #Sanity check:
  ## Consistent miRNA naming on the row:
  rownames(Input_df) <- gsub("_","-",rownames(Input_df))
  if (!all(grepl("miR-|let-", rownames(Input_df)))) {
    stop("Error: all row names must be miRNA names (e.g., hsa-miR-1-3p)!")
  }
  ## in non-log scale:
  if (max(as.matrix(Input_df)) < 20) {
    warning("Max read count detected is less than 20! \n -The input data should be in non-log space.")
  }


  # apply sample filtering:
  if (!is.na(outlier_sample.cutoff)){
    libSizes <- colSums(as.matrix(Input_df))
    filtersamples <- function(filterParam, times= outlier_sample.cutoff ){
      samplesToRemove <- which(filterParam > median(filterParam) + times * mad(filterParam) | filterParam < median(filterParam) - times * mad(filterParam) )
      samplesToRemove
    }

    samplesToRemove <- unique(unlist(lapply(list(libSizes = libSizes), filtersamples)))
    tt <- colnames(Input_df)[samplesToRemove]
    if (length(tt)>0) {print("Samples removed:"); print(paste( tt,collapse = ", "))}
    Input_df <- Input_df[,-samplesToRemove]
  }

  # detection rate cutoff:
  keep <- which(Matrix::rowSums(Input_df > 0) >= round( Dtct_cutoff * ncol(Input_df)))
  Input_df = Input_df[keep,]

  # Normalize the input data with TMM or alternatives; this will NOT affect the miR-TS output.
  dge <- edgeR::DGEList(counts=as.matrix(Input_df))
  dge <- edgeR::calcNormFactors(dge, method = Norm)
  Y <- edgeR::cpm(dge)

  # check % of miRNAs overlapping with signature matrix:

  X = signature_matrix
  rownames(X) <- gsub("-","_",rownames(X))
  rownames(Y) <- gsub("-","_",rownames(Y))
  Num_inSig <- sum(rownames(Y) %in% rownames(X))
  Prop_inSig <- round(Num_inSig/nrow(X),digits = 2)
  print(paste("Number of miRNAs using:", Num_inSig, "(", Prop_inSig,")"))

  # Deconvolution:####
  if (method == "cibersort"){
    if (!exists("CIBERSORT")) {
      stop("CIBERSORT function not detected! Run: \n1) ?CIBERSORT_download; \n2) navigate: Menu -> CS Archive -> CS Download -> Download CIBERSORT source code; \n3) source('CIBERSORT.R')")
    }
    X <- X[order(rownames(X)),,drop=FALSE]
    Y <- Y[order(rownames(Y)),,drop=FALSE]

    # Due to the set-up of the CIBERSORT.R, I need:
    # write X and Y to temporary .txt files,
    # call CIBERSORT() with those file paths,
    # then delete the files.
    tmp_sig <- tempfile("miRTS_sig_", fileext = ".txt")
    tmp_mix <- tempfile("miRTS_mix_", fileext = ".txt")

    # ensure cleanup no matter what
    on.exit({
      if (file.exists(tmp_sig)) unlink(tmp_sig)
      if (file.exists(tmp_mix)) unlink(tmp_mix)
    }, add = TRUE)

    # write X and Y in the format CIBERSORT expects
    write.table(X, file = tmp_sig, sep = "\t", quote = FALSE, col.names = NA)
    write.table(Y, file = tmp_mix, sep = "\t", quote = FALSE, col.names = NA)

    # now call the unmodified CIBERSORT() that reads from file paths
    CIBERSORT_output <- CIBERSORT(sig_matrix = tmp_sig, mixture_file = tmp_mix, perm = 0, QN = FALSE)

    deconv.res <- list(proportions=CIBERSORT_output, mix = Y, signatures = X)

    # if you are using CIBERSORT V1.03, comment out the above and uncomment the following line:
    # deconv.res <- CIBERSORT(sig_matrix = X,mixture_file = Y, perm = 0, QN = F)

  } else if (method == "xCell2" & requireNamespace("xCell2", quietly = T)){
    min_overlap_miR = (length(intersect(rownames(Y), rownames(X))) - 1)/nrow(X)
    xcell2_results <- xCell2::xCell2Analysis(mix = Input_df,
                                             xcell2object = TA2_miR.raw_qc.257...DICE_demo.xCell2Ref,
                                             minSharedGenes = min_overlap_miR)
    deconv.res <- list(proportions=as.matrix(t(xcell2_results)), mix = Y, signatures = X)
  } else if (method == "MCP-counter" & exists("compute_mcp_mirna_scores")){
    ts_df = try( ts_miR.df[, c("miRNA", "enriched_organ")] )
    if (inherits(ts_df, 'try-error')) { stop("ts_miR.df is not correctly specified!")}
    mcp_scores <- compute_mcp_mirna_scores(
      expr_tpm    = Input_df,
      ts_df       = ts_df
    )
    deconv.res <- list(proportions=mcp_scores, mix = Y, signatures = X)
  }

  # abs mode: ####
  if (method == "cibersort" &  useAbsolute == TRUE){
    X=X[rownames(X) %in% rownames(Y), ]

    ### !!! use Absolute Mode??? (RECOMMENDED)
    print("Using absolute mode!")
    # Convert the output to the absolute mode:
    CIBERSORT.conv_ABS <- function(Input_mix,
                                   Input_Signature,
                                   Output_rel,
                                   Suppress_note = FALSE){
      if (Suppress_note==FALSE){
        print("Input mixture: rownames= miRs, colnames= sample names!")
        print("Input Signature: rownames= miRs, colnames= tissue types! -- note: these miRs should overlap with Input mixture!")
        print("Cibersort relative mode output to be converted: rownames= sample names, colnames= tissue types!")
      }
      Output_rel <- as.data.frame(Output_rel)
      # drop extra columns if present
      drop_cols <- intersect(c("P-value", "Correlation", "RMSE"), colnames(Output_rel))
      if (length(drop_cols) > 0) {
        Output_rel <- Output_rel[, setdiff(colnames(Output_rel), drop_cols), drop = FALSE]
      }
      Input_Signature=Input_Signature[, colnames(Output_rel)]
      if (!identical(colnames(Input_Signature), colnames(Output_rel)) ||
          !identical(rownames(Output_rel), colnames(Input_mix))) {
        stop("Error: Check sample names (mix vs CIBERSORT output) and tissue names (sig vs output).")
      }

      miR_overlap <- intersect(rownames(Input_Signature), rownames(Input_mix))
      if (length(miR_overlap) == 0) stop("Error: no overlapping miRs!")
      sample_mean <- apply(Input_mix[miR_overlap, ], 2, function(x) mean(x))
      sample_avg = ifelse(median(as.matrix(Input_mix))==0, # To avoid 0 in the denominator.
                          mean(as.matrix(Input_mix)),
                          median(as.matrix(Input_mix)))
      scaling_factor <- as.numeric(sample_mean/sample_avg)
      sweep(Output_rel, MARGIN = 1, scaling_factor, `*`)
    }

    res.abs <- CIBERSORT.conv_ABS(Input_mix=Y,
                                  Input_Signature= X,
                                  Output_rel=deconv.res$proportions,
                                  Suppress_note = TRUE)

    orig <- as.data.frame(deconv.res$proportions)
    keep_diag <- intersect(c("P-value", "Correlation", "RMSE"), colnames(orig))
    diag <- if (length(keep_diag) > 0) orig[, keep_diag, drop = FALSE] else NULL
    deconv.res$proportions <- cbind(res.abs, diag)
  }
  # replace the 0 with min/2:
  Transform_cibersort <- function(df) {
    df <- as.data.frame(df)

    # drop all-zero columns
    nz <- colSums(as.matrix(df), na.rm = TRUE) != 0
    df <- df[, nz, drop = FALSE]

    # drop diagnostics if present
    drop_cols <- intersect(c("P-value","Correlation","RMSE"), colnames(df))
    if (length(drop_cols) > 0) {
      df <- df[, setdiff(colnames(df), drop_cols), drop = FALSE]
    }

    # replace zeros with min/2 (second-smallest if min is 0)
    min_tissue <- vapply(df, function(x) {
      x <- x[is.finite(x)]
      if (!length(x)) return(NA_real_)
      ux <- sort(unique(x))
      if (length(ux) == 1) return(ux[1])
      if (ux[1] == 0) ux[2] else ux[1]
    }, numeric(1))

    for (i in seq_len(ncol(df))) {
      if (!is.na(min_tissue[i])) df[df[[i]] == 0, i] <- min_tissue[i] / 2
    }
    df
  }


  if (method == "cibersort"){
    deconv.res$proportions <- Transform_cibersort(as.data.frame(deconv.res$proportions))
    # output only the 17 tissues with detectable signals in plasma:
    if (Use_blood == TRUE) {
      blood_detect <- intersect(colnames(deconv.res$proportions), miRTS::Tissues.blood_detect_17)
      deconv.res$proportions <- deconv.res$proportions[, blood_detect]
    }
  }
  deconv.res

}
