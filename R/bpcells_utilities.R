#' Compute SCTransform model parameters without full residuals
#'
#' This function computes SCTransform model parameters and saves them to the
#' misc slot of a Seurat object without computing residuals over the entire
#' matrix. This is useful for large datasets where computing full residuals
#' would be memory-intensive.
#'
#' @param object A Seurat object
#' @param assay Name of assay to use (default "RNA")
#' @param layer Name of layer to use (default "counts")
#' @param ncells Number of cells to use for initial parameter estimation (default 5000)
#' @param n_genes Number of genes to use for initial parameter estimation (default 2000)
#' @param vst.flavor Flavor of VST to use: "v1" or "v2" (default "v2")
#' @param min_variance Minimum variance setting (default "umi_median")
#' @param verbose Whether to print messages (default TRUE)
#' @param seed.use Seed to set
#' @param ... Additional arguments passed to vst()
#'
#' @return A Seurat object with SCTransform model parameters saved in
#'   object@@misc$RNAFingerprinting$SCT_model
#'
#' @importFrom Seurat GetAssayData
#' @importFrom SeuratObject as.sparse
#' @importFrom Matrix colSums rowSums rowMeans
#' @importFrom stats median
#' @importFrom sctransform vst get_residual_var
#'
#' @export
ComputeSCTModel <- function(
    object,
    assay = "RNA",
    layer = "counts",
    ncells = 5000,
    n_genes = 2000,
    vst.flavor = "v2",
    min_variance = "umi_median",
    verbose = TRUE,
    seed.use = 1448145,
    ...
) {
  if (!is.null(seed.use)) {
    set.seed(seed = seed.use)
  }
  
  # Extract count matrix
  umi <- GetAssayData(object, assay = assay, layer = layer)
  
  is_iterable <- inherits(umi, "IterableMatrix")
  if (is_iterable && verbose) {
    message("Detected BPCells/IterableMatrix - using streaming operations")
  }
  
  # Compute clip range based on full matrix dimensions
  clip.range <- c(
    -sqrt(ncol(umi) / 30), 
    sqrt(ncol(umi) / 30)
  )
  
  # Prepare cell attributes for ALL cells
  if (verbose) {
    message("Computing cell attributes for all cells")
  }
  
  if (is_iterable) {
    # Use BPCells streaming operations for efficient computation
    col_stats <- BPCells::matrix_stats(umi, col_stats = "mean")
    cell_umi_counts <- col_stats$col_stats["mean", ] * nrow(umi)
    names(cell_umi_counts) <- colnames(umi)
  } else {
    # For regular matrix, compute column sums
    cell_umi_counts <- Matrix::colSums(umi)
  }
  
  # Create cell_attr for ALL cells
  all_cell_attr <- data.frame(
    log_umi = log10(cell_umi_counts),
    row.names = names(cell_umi_counts)
  )
  
  # Sample cells for model fitting
  sampled_cells <- sample.int(
    n = ncol(umi), 
    size = min(ncells, ncol(umi))
  )
  sampled_cell_names <- colnames(umi)[sampled_cells]
  
  # Extract subset for vst
  if (is_iterable) {
    umi_subset <- as.sparse(umi[, sampled_cell_names])
  } else {
    umi_subset <- umi[, sampled_cell_names]
  }
  
  cell_attr_subset <- all_cell_attr[sampled_cell_names, , drop = FALSE]
  
  # Compute gene attributes for ALL genes using full matrix
  if (verbose) {
    message("Computing gene attributes for all genes")
  }
  
  if (is_iterable) {
    # Use BPCells streaming operations
    row_stats_mean <- BPCells::matrix_stats(umi, row_stats = "mean")
    row_stats_var <- BPCells::matrix_stats(umi, row_stats = "variance")
    row_stats_nonzero <- BPCells::matrix_stats(umi, row_stats = "nonzero")
    
    gene_means <- row_stats_mean$row_stats["mean", ]
    gene_vars <- row_stats_var$row_stats["variance", ]
    gene_nonzero <- row_stats_nonzero$row_stats["nonzero", ]
    
    names(gene_means) <- rownames(umi)
    names(gene_vars) <- rownames(umi)
    names(gene_nonzero) <- rownames(umi)
    
    detection_rate <- gene_nonzero / ncol(umi)
    
    # Compute geometric mean using streaming (approximate via log-space mean)
    # For BPCells, we'll compute this on the subset used for fitting
    # to avoid materializing the full matrix
    umi_subset_mat <- as.matrix(umi_subset)
    gmean_eps <- 1
    genes_gmean_subset <- exp(rowMeans(log(umi_subset_mat + gmean_eps))) - gmean_eps
    
  } else {
    # For regular matrix
    detection_rate <- Matrix::rowSums(umi > 0) / ncol(umi)
    gene_means <- Matrix::rowMeans(umi)
    
    # Compute variance
    gene_vars <- apply(as.matrix(umi), 1, var)
    
    # Compute geometric mean
    gmean_eps <- 1
    umi_mat <- as.matrix(umi)
    genes_gmean <- exp(rowMeans(log(umi_mat + gmean_eps))) - gmean_eps
  }
  
  # Prepare vst arguments
  vst.args <- list(...)
  vst.args[['umi']] <- umi_subset
  vst.args[['cell_attr']] <- cell_attr_subset
  vst.args[['n_cells']] <- min(ncells, ncol(umi))
  vst.args[['n_genes']] <- n_genes
  vst.args[['vst.flavor']] <- vst.flavor
  vst.args[['min_variance']] <- min_variance
  vst.args[['verbosity']] <- as.numeric(verbose) * 2
  vst.args[['return_cell_attr']] <- TRUE
  vst.args[['return_gene_attr']] <- TRUE
  vst.args[['return_corrected_umi']] <- FALSE
  vst.args[['residual_type']] <- 'none'  # KEY: Don't compute residuals
  
  # Set res_clip_range for internal calculations
  vst.args[['res_clip_range']] <- c(-sqrt(ncol(umi)), sqrt(ncol(umi)))
  
  if (verbose) {
    message("Computing SCTransform model parameters (without full residuals)")
  }
  
  # Run vst to get model parameters
  vst_out <- do.call(what = vst, args = vst.args)
  
  # Update gene_attr with full matrix statistics
  genes <- rownames(vst_out$gene_attr)
  vst_out$gene_attr$detection_rate <- detection_rate[genes]
  vst_out$gene_attr$amean <- gene_means[genes]
  vst_out$gene_attr$variance <- gene_vars[genes]
  
  if (is_iterable) {
    # For IterableMatrix, use gmean from subset
    vst_out$gene_attr$gmean <- genes_gmean_subset[genes]
  } else {
    vst_out$gene_attr$gmean <- genes_gmean[genes]
  }
  
  # Compute residual mean and variance on the SUBSET used for fitting
  if (verbose) {
    message("Computing residual statistics on subset of cells used for fitting")
  }
  
  feature.variance <- get_residual_var(
    vst_out = vst_out,
    umi = umi_subset,
    residual_type = 'pearson',
    res_clip_range = vst.args[['res_clip_range']]
  )
  
  vst_out$gene_attr$residual_variance <- NA_real_
  vst_out$gene_attr[names(feature.variance), 'residual_variance'] <- feature.variance
  
  residuals_subset <- sctransform::get_residuals(
    vst_out = vst_out,
    umi = umi_subset,
    min_variance = min_variance,
    res_clip_range = vst.args[['res_clip_range']],
    verbosity = 0
  )
  
  vst_out$gene_attr$residual_mean <- NA_real_
  vst_out$gene_attr[rownames(residuals_subset), 'residual_mean'] <- rowMeans(residuals_subset)
  
  # Replace cell_attr with ALL cells
  vst_out$cell_attr <- all_cell_attr
  
  # Add sct.clip.range to arguments
  vst_out$arguments$sct.clip.range <- clip.range
  vst_out$arguments$sct.method <- "model_only"
  
  # Save to misc slot
  if (is.null(object@misc$RNAFingerprinting)) {
    object@misc$RNAFingerprinting <- list()
  }
  
  object@misc$RNAFingerprinting$SCT_model <- list(
    vst_out = vst_out
  )
  
  if (verbose) {
    message("SCTransform model parameters saved to object@misc$RNAFingerprinting$SCT_model")
  }
  
  return(object)
}

#' @keywords internal
FetchResidualsFromMisc <- function(
    object, 
    features, 
    cells
) {
  # Extract model and parameters
  sct_misc <- object@misc$RNAFingerprinting$SCT_model
  
  if (!is.null(sct_misc$vst_out)) {
    vst_out <- sct_misc$vst_out
  } else {
    model <- sct_misc$model_list[[1]]
    vst_out <- Seurat:::SCTModel_to_vst(model)
  }
  
  umi <- Seurat::GetAssayData(object, assay = "RNA", layer = "counts")
  umi_sub <- SeuratObject::as.sparse(umi[features, cells, drop = FALSE])
  clip.range <- vst_out$arguments$sct.clip.range
  clip.min <- min(clip.range)
  clip.max <- max(clip.range)
  
  cell_attr_sub <- vst_out$cell_attr[colnames(umi_sub), , drop = FALSE]
  
  # Resolve min_variance
  min_var <- vst_out$arguments$min_variance
  if (identical(min_var, "umi_median")) {
    if (is.null(object@misc$RNAFingerprinting$SCT_model$min_var_cached)) {
      chunk_size <- 2000L 
      all_cells <- colnames(Seurat::GetAssayData(object, assay = "RNA", layer = "counts"))
      cells0 <- all_cells[seq_len(min(chunk_size, length(all_cells)))]
      umi0 <- SeuratObject::as.sparse(
        Seurat::GetAssayData(object, assay = "RNA", layer = "counts")[, cells0, drop = FALSE]
      )
      nz_median <- stats::median(umi0@x)
      object@misc$RNAFingerprinting$SCT_model$min_var_cached <- (nz_median / 5)^2
    }
    min_var <- object@misc$RNAFingerprinting$SCT_model$min_var_cached
  }
  
  # Compute residuals
  vst_tmp <- vst_out
  vst_tmp$gene_attr <- vst_tmp$gene_attr[rownames(umi_sub), , drop = FALSE]
  vst_tmp$model_pars_fit <- vst_tmp$model_pars_fit[rownames(umi_sub), , drop = FALSE]
  vst_tmp$cell_attr <- cell_attr_sub
  res <- sctransform::get_residuals(
    vst_out = vst_tmp,
    umi = umi_sub,
    min_variance = min_var,
    res_clip_range = c(clip.min, clip.max),
    verbosity = 0
  )
  
  res <- as.matrix(res)
  mu <- vst_out$gene_attr[rownames(res), "residual_mean"]
  res <- sweep(res, 1, mu, FUN = "-")
  return(res)
}