#' Choose features
#'
#' Choose features for fingerprint estimation.
#'
#' @param object A Seurat object
#' @param method "de" to choose features with differential expression; "expression" to choose based on expression levels
#' @param perturbation_meta If method = "de": the name of the metadata column containing the perturbation identities
#' @param control_class If method = "de": the name of the perturbation identity corresponding to the control condition
#' @param max_per_perturbation If method = "de": the maximum number of DEGs to use per perturbation (default 500)
#' @param mean_expression If method = "expression": the minimum average expression (default 0.01)
#'
#' @return A set of features to use for fingerprint estimation
#'
#' @import Seurat
#' @importFrom Matrix rowMeans
#'
#' @export
#'
#' @examples
#' \dontrun{
#' features <- ChooseFeatures(
#' object,
#' method = 'de',
#' perturbation_meta = 'perturbation_id',
#' control_class = 'NT'
#' )
#'
#' features <- ChooseFeatures(
#' object,
#' method = 'expression'
#' )
#' }
#'
ChooseFeatures <- function(
    object,
    method,
    perturbation_meta = NULL,
    control_class = NULL,
    max_per_perturbation = 500,
    mean_expression = 0.01
) {
  # Verify parameters
  if (!method %in% c('de','expression')) {
    stop('Method must be one of de, expression')
  }

  # DE-based feature selection
  if (method == 'de') {
    if (is.null(perturbation_meta) | is.null(control_class)) {
      stop('Please specify perturbation_meta and control_class')
    }
    if (!perturbation_meta %in% colnames(object@meta.data)) {
      stop('Perturbation metadata column not found.')
    }
    if (sum(object@meta.data[,perturbation_meta]==control_class)==0) {
      stop('No cells with control class label found.')
    }
    features <- NULL
    perturbations <- unique(object@meta.data[,perturbation_meta][
      object@meta.data[,perturbation_meta]!=control_class
    ])
    for (p in perturbations) {
      markers <- FindMarkers(
        object,
        ident.1 = p,
        ident.2 = control_class,
        group.by = perturbation_meta
      )
      num <- min(sum(markers$p_val<0.05),max_per_perturbation)
      features <- c(features,rownames(markers)[1:num])
    }
    features <- unique(features)
  }

  # Expression-based feature selection
  if (method == 'expression') {
    counts <- GetAssayData(
      object,
      assay = 'RNA',
      layer = 'counts'
    )
    features <- rownames(counts)[Matrix::rowMeans(counts) > mean_expression]
  }

  return(features)
}

#' Convert genes
#'
#' Convert genes in a dictionary object
#'
#' @param dictionary Dictionary obtained from LearnDictionary()
#' @param from Key type to convert from (default 'ENSEMBL')
#' @param to Key type to convert to (default 'SYMBOL')
#'
#' @return Dictionary with re-named genes; genes that couldn't be converted are dropped
#'
#' @import org.Hs.eg.db
#' @importFrom AnnotationDbi mapIds
#'
#' @export
#'
#' @examples
#' \dontrun{
#' dictionary <- ConvertGenes(
#' dictionary
#' )
#' }
#'
ConvertGenes <- function(
    dictionary,
    from='ENSEMBL',
    to='SYMBOL'
) {
  # Convert genes
  genes <- mapIds(
    org.Hs.eg.db,
    keys = rownames(dictionary$fingerprints$Lambdas),
    keytype = from,
    column = to
  )
  to_keep <- which(!is.na(genes) & !duplicated(genes))

  # Rename
  rownames(dictionary$fingerprints$Lambdas) <- genes
  dictionary$fingerprints$Lambdas <- dictionary$fingerprints$Lambdas[to_keep,]
  for (p in 1:length(dictionary$fingerprints$dists)) {
    rownames(dictionary$fingerprints$dists[[p]]) <- genes
    dictionary$fingerprints$dists[[p]] <- dictionary$fingerprints$dists[[p]][to_keep,]
  }

  if (!is.null(dictionary$data)) {
    genes <- mapIds(
      org.Hs.eg.db,
      keys = rownames(dictionary$data),
      keytype = from,
      column = to
    )
    to_keep <- which(!is.na(genes) & !duplicated(genes))
    rownames(dictionary$data) <- genes
    dictionary$data <- dictionary$data[to_keep, ]
  }

  return(dictionary)
}

#' @keywords internal
split_and_pad <- function(x, sep, col.name) {
  split <- strsplit(x, sep)
  max_len <- max(lengths(split))
  split_padded <- do.call(rbind, lapply(split, function(x) {
    length(x) <- max_len
    x
  }))
  colnames(split_padded) <- paste0(col.name, 1:ncol(split_padded))
  split_padded[split_padded=='NA'] <- NA
  return(split_padded)
}

#' Summarize results
#'
#' Summarize results after mapping cells to fingerprints.
#'
#' @param object A Seurat object
#' @param group_level Whether fingerprinting was run at the group level (TRUE) or not (FALSE)
#' @param condition_meta The name of the metadata column containing the condition identities
#' @param control_class The name of the perturbation identity corresponding to the control condition
#' @param suffix The suffix of the output metadata columns; defaults to NULL if cell-level or the specified condition name if group-level
#' @param tidy Whether results should be in tidy format (default FALSE)
#'
#' @return A dataframe summarizing the mapping results
#'
#' @import Seurat
#'
#' @export
#'
#' @examples
#' \dontrun{
#' summary <- SummarizeResults(
#' object,
#' group_level = TRUE,
#' condition_meta = 'perturbation_id',
#' control_class = 'NT'
#' )
#'
#' summary <- SummarizeResults(
#' object,
#' group_level = FALSE,
#' condition_meta = 'perturbation_id',
#' control_class = 'NT'
#' )
#'
#' summary <- SummarizeResults(
#' object,
#' group_level = FALSE,
#' condition_meta = 'perturbation_id',
#' control_class = 'NT',
#' tidy = TRUE
#' )
#' }
#'
SummarizeResults <- function(
    object,
    group_level,
    condition_meta,
    control_class,
    suffix = NULL,
    tidy = FALSE
) {
  # Set parameters
  if (is.null(suffix) & group_level == TRUE) {
    suffix <- condition_meta
  }

  # Extract results
  metadata_cols <- c('top_credible_set', 'top_lbf', 'all_credible_sets',
                     'all_lbfs', 'all_signs')
  if (!is.null(suffix)) {
    metadata_cols <- paste0(metadata_cols, '.', suffix)
  }
  summary <- object@meta.data[,c(condition_meta, metadata_cols)]
  summary <- summary[summary[,1] != control_class, ]
  if (group_level) {
    num_conditions <- length(unique(summary[,1]))
    summary <- summary[!duplicated(summary),]
    if (nrow(summary) != num_conditions) {
      stop('Distinct results found for the same condition -- should group_level be FALSE?')
    }
    rownames(summary) <- summary[,1]
    summary <- summary[, 2:ncol(summary)]
  }
  summary <- summary[order(summary[, metadata_cols[2]],decreasing=T), ]

  # Put into tidy-format if requested
  if (tidy) {
    split_summary <- cbind(split_and_pad(summary[, metadata_cols[3]], ';', 'cs'),
                           split_and_pad(summary[, metadata_cols[4]], ';', 'lbf'),
                           split_and_pad(summary[, metadata_cols[5]], ';', 'sign'))
    rownames(split_summary) <- rownames(summary)

    tidy_summary <- do.call(rbind, lapply(1:(ncol(split_summary)/3),function(x) {
      split_summary[,paste0(c('cs','lbf','sign'),x),drop=FALSE]
    }))
    tidy_summary <- cbind(rownames(tidy_summary), tidy_summary)
    tidy_summary <- data.frame(tidy_summary[!is.na(tidy_summary[,2]),])
    colnames(tidy_summary) <- c('query','cs','lbf','sign')

    tidy_summary <- tidy_summary[
      order(factor(tidy_summary$query, levels = unique(tidy_summary$query))),]
    rownames(tidy_summary) <- NULL
    return(tidy_summary)
  } else {
    return(summary)
  }
}

#' @keywords internal
get_model_attribute <- function(
    data,
    query,
    condition_meta,
    fingerprints,
    group_level,
    attribute = c('betas','bfs','correlation','residuals'),
    match = NULL
) {
  Lambdas <- fingerprints$Lambdas
  features <- rownames(Lambdas)

  ctrl.pdev <- data[features,condition_meta=='control']
  ctrl.mus <- rowMeans(ctrl.pdev)
  PCs <- prcomp_irlba(t(ctrl.pdev),n=30)

  if (group_level) {
    w.pdev <- data[features,condition_meta==query,drop=F]
  } else {
    w.pdev <- data[features,query,drop=F]
  }
  w.reconstruct <- sweep(PCs$rotation%*%crossprod(PCs$rotation,sweep(w.pdev,1,ctrl.mus)),
                         1,ctrl.mus,'+')
  w.residual <- w.pdev - w.reconstruct

  if (attribute == 'residuals') {
    return(w.residual)
  } else {
    w.residual <- rowSums(w.residual)
  }

  if (attribute == 'correlation') {
    return(cor(w.residual,Lambdas)[1,])
  }

  if (!group_level) {
    top_genes <- order(abs(w.residual),decreasing=T)[1:5]
    top_genes <- unique(c(top_genes,which(abs(w.residual)>=10)))
    dummy <- array(0,dim=c(nrow(Lambdas),length(top_genes)))
    for (col in 1:ncol(dummy)) {
      dummy[top_genes[col],col] <- 1
    }
    dummy_vars <- (0.5*abs(w.residual[top_genes])/qnorm(.95))^2
  } else {
    dummy <- dummy_vars <- NULL
  }
  mod <- susie_plus(X=Lambdas,y=w.residual,S=dummy,varS=dummy_vars,
                    standardize=F,L=10,intercept=F,max_iter=5000,
                    estimate_prior_method='optim')

  if (attribute == 'betas') {
    betas <- sapply(1:ncol(mod$mu),function(x)
      mod$mu[which.max(mod$lbf_variable[,x]),x])
    names(betas) <- colnames(Lambdas)
    return(betas)
  }
  if (attribute == 'bfs') {
    cs <- which.max(mod$lbf_variable[,match])
    fingerprint <- Lambdas[,match,drop=F]
    bf <- get_bf(mod,fingerprint,
                 colSums(fingerprint*fingerprint),
                 w.residual,cs)
    bf_prime <- sapply(1:length(w.residual),function(x) {
      get_bf(mod,fingerprint[-x,,drop=F],
             colSums(fingerprint[-x,,drop=F]*fingerprint[-x,,drop=F]),
             w.residual[-x],cs)-bf
    })
    rank <- order(bf_prime,decreasing=F)
    names(bf_prime) <- rownames(Lambdas)
    return(bf_prime[rank])
  }
}
