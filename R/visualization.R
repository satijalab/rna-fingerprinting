#' Long Island City plot
#'
#' Visualize assignment with a Long Island City plot.
#'
#' @param object A Seurat object
#' @param query Either a cell name or a group name
#' @param condition_meta The name of the metadata column containing the condition identities
#' @param control_class The name of the condition identity corresponding to the control condition
#' @param dictionary The dictionary of fingerprints that the data was mapped to
#' @param group_level Whether fingerprinting was run at the group level (TRUE) or not (FALSE)
#' @param suffix The suffix of the output metadata columns; defaults to NULL if cell-level or the specified condition name if group-level
#' @param approx If TRUE, plot marginal correlations with fingerprints; if FALSE (default), plot conditional posterior mean coefficients
#' @param order Whether to order fingerprints by clustering (default) or alphabetically; ignored if order is stored in fingerprints
#' @param cs The maximum number of top credible sets to visualize (default 1)
#'
#' @return A Long Island City plot visualizing assignment of the query
#'
#' @import ggplot2
#' @importFrom ggrepel geom_label_repel
#'
#' @export
#'
#' @examples
#' \dontrun{
#' LongIslandCityPlot(
#' object,
#' query = 'R1-L1_AAAGGATCACCAGCTG',
#' condition_meta = 'condition_id',
#' control_class = 'untreated',
#' dictionary = dictionary,
#' group_level = FALSE,
#' cs = 1
#' )
#'
#' LongIslandCityPlot(
#' object,
#' query = 'drug1',
#' condition_meta = 'condition_id',
#' control_class = 'untreated',
#' dictionary = dictionary,
#' group_level = TRUE,
#' cs = 10
#' )
#' }
#'
LongIslandCityPlot <- function(
    object,
    query,
    condition_meta,
    control_class,
    dictionary,
    group_level,
    suffix = NULL,
    approx = FALSE,
    order = 'clustering',
    cs = 10
) {
  # Verify parameters
  if (group_level & !query %in% object@meta.data[,condition_meta]) {
    stop('Query not found in condition metadata')
  }
  if (!group_level & !query %in% colnames(object)) {
    stop('Query not found among cells')
  }
  if (!approx & is.null(getOption("fingerprinting.internal.LICplot"))) {
    message('This can be slow in large settings; you can speed up this function with approx = TRUE.')
    options(fingerprinting.internal.LICplot = TRUE)
  }
  if (!order %in% c('clustering','alphabetical')) {
    stop('Order must be one of clustering, alphabetical')
  }
  if (is.null(suffix) & group_level == TRUE) {
    suffix <- condition_meta
  }

  # Set up data and features
  sct_data <- GetAssayData(
    object,
    assay = 'SCT',
    layer = 'scale.data'
  )
  features <- rownames(dictionary$fingerprints$Lambdas)
  subset_features <- intersect(features, rownames(sct_data))
  dictionary$fingerprints$Lambdas <- dictionary$fingerprints$Lambdas[subset_features,]

  condition_metadata <- as.character(object@meta.data[,condition_meta])
  condition_metadata[condition_metadata==control_class] <- 'control'
  names(condition_metadata) <- rownames(object@meta.data)

  # Get plot components
  if (approx) {
    attribute <- 'correlation'
  } else {
    attribute <- 'betas'
  }
  betas <- get_model_attribute(
    data = sct_data,
    query = query,
    condition_meta = condition_metadata,
    fingerprints = dictionary$fingerprints,
    group_level = group_level,
    attribute = attribute
  )

  Lambdas <- dictionary$fingerprints$Lambdas
  if (!is.null(dictionary$order)) {
    hc <- dictionary$order
  } else if (order == 'clustering') {
    hc <- order.dendrogram(as.dendrogram(hclust(dist(cor(Lambdas)))))
  } else {
    hc <- order(colnames(Lambdas))
  }
  dictionary <- factor(colnames(Lambdas), levels=colnames(Lambdas)[hc])

  all_credible_sets_metadata <- 'all_credible_sets'
  if (!is.null(suffix)) {
    all_credible_sets_metadata <- paste0(all_credible_sets_metadata,'.',suffix)
  }
  all_credible_sets <- object@meta.data[,all_credible_sets_metadata]
  if (group_level) {
    assign <- all_credible_sets[condition_metadata==query][1]
  } else {
    assign <- all_credible_sets[colnames(object)==query]
  }
  prtbs <- strsplit(assign,';')[[1]]
  prtbs <- prtbs[1:min(c(length(prtbs),cs))]

  # Make plot
  plot_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
              "#FFFF33", "#A65628", "#F781BF", "#999999", "#66C2A5")
  df <- data.frame(Beta=betas,Dictionary=dictionary)
  plot <- ggplot()+
    geom_point(data=df,
               aes(x=Dictionary,y=Beta))+
    theme_classic(base_size=12)+
    theme(axis.text.x=element_blank())+
    ggtitle(query)
  for (cs_ind in 1:length(prtbs)) {
    cs_members <- strsplit(prtbs[cs_ind],',')[[1]]
    plot <- plot+
      geom_point(data=df[df$Dictionary %in% cs_members,],
                 aes(x=Dictionary,y=Beta),color=plot_colors[cs_ind])+
      geom_label_repel(data=df[df$Dictionary %in% cs_members,],
        aes(x=Dictionary,y=Beta,label=Dictionary),color=plot_colors[cs_ind])
  }
  if (approx) {
    plot <- plot + ylab('Correlation')
  }

  return(plot)
}

#' Explain match
#'
#' Explain an assignment by identifying the top genes driving a match.
#'
#' @param object A Seurat object
#' @param query Either a cell name or a group name
#' @param match The perturbation match to be explained; should be the name of one fingerprint in the dictionary
#' @param condition_meta The name of the metadata column containing the condition identities
#' @param control_class The name of the condition identity corresponding to the control condition
#' @param dictionary The dictionary of fingerprints that were used to map the data
#' @param group_level Whether fingerprinting was run at the group level (TRUE) or not (FALSE)
#' @param ref The Seurat object for the reference data from which the fingerprints were learned (only needed if plot is TRUE); if NULL, will use pseudobulked data stored in dictionary
#' @param perturbation_meta The name of the metadata column containing the perturbation identities in the reference data (only needed if plot is TRUE)
#' @param control_class_ref The name of the perturbation identity corresponding to the control condition in the reference data (only needed if plot is TRUE)
#' @param plot If TRUE (default for group-level), will return a heatmap comparison; otherwise, or for cell-level assignments, returns the top genes
#' @param num_genes Number of top genes to show (default 10)
#' @param include_de Number of DEGs from the reference and query each to show (default 0), only needed if plot is TRUE
#' @param query_key If different gene key types, the key type of the query object; NULL (default) if no conversion needed (only needed if plot is TRUE)
#' @param ref_key If different gene key types, the key type of the reference object; NULL (default) if no conversion needed (only needed if plot is TRUE)
#' @param downsample_control Whether the number of control cells should be downsampled (default TRUE), only needed if plot is TRUE
#'
#' @return If plot = TRUE, two heatmaps corresponding to the query and reference data explaining the match that was made; otherwise, a dataframe of the top genes responsible for a given match
#'
#' @import ggplot2
#' @import org.Hs.eg.db
#' @import Seurat
#' @importFrom AnnotationDbi mapIds
#'
#' @export
#'
ExplainMatch <- function(
    object,
    query,
    match,
    condition_meta,
    control_class,
    dictionary,
    group_level,
    ref = NULL,
    perturbation_meta = NULL,
    control_class_ref = NULL,
    plot = TRUE,
    num_genes = 25,
    include_de = 0,
    query_key = NULL,
    ref_key = NULL,
    downsample_control = TRUE
) {
  # Set up reference pseudobulk data if precomputed dictionary
  if (is.null(ref) && plot == TRUE) {
    ref_pseudobulk <- dictionary$data
  }

  # Determine if conversion is needed
  conversion_needed <- FALSE
  if (!is.null(query_key) & !is.null(ref_key) & inherits(ref, 'Seurat')) {
    conversion_needed <- TRUE
  }

  # Set up data and features
  sct_data <- GetAssayData(
    object,
    assay = 'SCT',
    layer = 'scale.data'
  )
  features <- rownames(dictionary$fingerprints$Lambdas)
  subset_features <- intersect(features, rownames(sct_data))
  dictionary$fingerprints$Lambdas <- dictionary$fingerprints$Lambdas[subset_features,]

  condition_metadata <- as.character(object@meta.data[,condition_meta])
  condition_metadata[condition_metadata==control_class] <- 'control'
  names(condition_metadata) <- rownames(object@meta.data)

  # Get top genes driving the match
  ranked_genes <- get_model_attribute(
    data = sct_data,
    query = query,
    condition_meta = condition_metadata,
    fingerprints = dictionary$fingerprints,
    group_level = group_level,
    attribute = 'bfs',
    match = match
  )
  top_genes <- names(ranked_genes)[1:num_genes]

  # Return if plot not requested, or if cell-level assignment
  if ((!plot) | (!group_level)) {
    top_genes_df <- data.frame(genes = top_genes,
                               delta_BF = ranked_genes[1:num_genes])
    return(top_genes_df)
  }

  # Re-order top genes based on up/down
  top_genes_signs <- sign(dictionary$fingerprints$Lambdas[top_genes, match])
  top_genes <- c(top_genes[top_genes_signs == 1],
                 top_genes[top_genes_signs == (-1)])

  # Get DE genes if requested
  if (include_de > 0) {
    query_features <- rownames(GetAssayData(object, layer = 'scale.data'))
    if (inherits(ref, 'Seurat')) {
      ref_features <- rownames(GetAssayData(ref, layer = 'scale.data'))
      if (conversion_needed) {
        ref_features <- mapIds(
          org.Hs.eg.db,
          keys = ref_features,
          keytype = ref_key,
          column = query_key
        )
      }
    } else {
      ref_features <- rownames(ref_pseudobulk)
    }
    common_features <- intersect(query_features, ref_features)
    if (length(common_features) < (num_genes + 2*include_de)) {
      stop('Not enough common features between reference and query -- do you need to specify ref_key, query_key to convert features?')
    }

    de_query <- FindMarkers(
      object,
      ident.1 = query,
      ident.2 = control_class,
      group.by = condition_meta,
      features = common_features
    )
    de_query <- de_query[which(!rownames(de_query) %in% top_genes)[1:include_de],]
    de_query <- rownames(de_query)[order(de_query$avg_log2FC)]

    if (inherits(ref, 'Seurat')) {
      de_ref <- FindMarkers(
        ref,
        ident.1 = match,
        ident.2 = control_class_ref,
        group.by = perturbation_meta,
        features = if (conversion_needed) {
          mapIds(org.Hs.eg.db, keys = common_features,
                 from = ref_key, to = query_key)
        } else {
          common_features
        }
      )
      if (conversion_needed) {
        rownames(de_ref) <- mapIds(
          org.Hs.eg.db,
          keys = rownames(de_ref),
          keytype = query_key,
          column = ref_key
        )
        de_ref <- de_ref[!is.na(rownames(de_ref)),]
      }
      de_ref <- de_ref[which(!rownames(de_ref) %in% c(top_genes, de_query))[1:include_de],]
      de_ref <- rownames(de_ref)[order(de_ref$avg_log2FC)]
    } else {
      values <- GetAssayData(ref_pseudobulk, layer = 'scale.data')[, c('control', match)]
      to_test <- common_features[which(!common_features %in% c(top_genes, de_query))]
      de_ref <- values[to_test, 2] - values[to_test, 1]
      de_ref <- de_ref[order(abs(de_ref),decreasing = T)[1:include_de]]
      de_ref <- names(de_ref)[order(de_ref)]
    }
  } else {
    de_query <- NULL
    de_ref <- NULL
  }

  # Prepare to plot
  features_to_plot <- c(top_genes, de_query, de_ref)
  query_condition_cells <- colnames(object)[object@meta.data[, condition_meta] == query]
  query_control_cells <- colnames(object)[object@meta.data[, condition_meta] == control_class]
  if (downsample_control) {
    query_control_cells <- sample(query_control_cells,
                                  min(c(length(query_control_cells),
                                        length(query_condition_cells))))
  }
  object_sub <- object[,c(query_condition_cells, query_control_cells)]
  object_sub@meta.data[, condition_meta] <- factor(object_sub@meta.data[, condition_meta],
                                                   levels = c(query, control_class))

  if (inherits(ref, 'Seurat')) {
    ref_perturbed_cells <- colnames(ref)[ref@meta.data[, perturbation_meta] == match]
    ref_control_cells <- colnames(ref)[ref@meta.data[, perturbation_meta] == control_class_ref]
    if (downsample_control) {
      ref_control_cells <- sample(ref_control_cells,
                                  min(c(length(ref_control_cells),
                                        length(ref_perturbed_cells))))
    }
    ref_sub <- ref[,c(ref_perturbed_cells, ref_control_cells)]
    ref_sub@meta.data[, perturbation_meta] <- factor(ref_sub@meta.data[, perturbation_meta],
                                                     levels = c(match, control_class_ref))
  } else {
    ref_sub <- ref_pseudobulk[,c('control', match)]
    ref_sub$gene <- factor(colnames(ref_sub), levels = c(match, 'control'))
    perturbation_meta <- 'gene'
  }

  # Make plots
  h1 <- DoHeatmap(
    object_sub,
    features = features_to_plot,
    group.by = condition_meta,
    raster = FALSE,
    label = F
  ) + ggtitle('Query')

  h2_title <- if (inherits(ref, 'Seurat')) { 'Reference' } else { 'Pseudobulked Reference' }
  h2 <- DoHeatmap(
    ref_sub,
    features = if (conversion_needed) {
      mapIds(org.Hs.eg.db, keys = features_to_plot,
             from = query_key, to = ref_key)
    } else {
      features_to_plot
    },
    group.by = perturbation_meta,
    raster = FALSE,
    label = F,
    draw.lines = if (!inherits(ref, 'Seurat')) { FALSE } else { TRUE },
  ) + ggtitle(h2_title)

  return(h1+h2)
}

#' Project fingerprints
#'
#' Project assigned cells to the fingerprints and produce a UMAP visualization.
#'
#' @param object A Seurat object
#' @param condition_meta The name of the metadata column containing the condition identities
#' @param control_class The name of the condition identity corresponding to the control condition
#' @param dictionary The dictionary of fingerprints that were used to map the data
#' @param suffix The suffix of the output metadata columns; defaults to NULL
#'
#' @return A Seurat object with a 'proj' assay containing projected values and 'projumap' reduction
#'
#' @import Seurat
#'
#' @export
#'
#' @examples
#' \dontrun{
#' object <- ProjectFingerprints(
#' object,
#' condition_meta = 'condition_id',
#' control_class = 'untreated',
#' dictionary = dictionary
#' )
#' }
#'
ProjectFingerprints <- function(
    object,
    condition_meta,
    control_class,
    dictionary,
    suffix = NULL
) {
  # Set up data and features
  sct_data <- GetAssayData(
    object,
    assay = 'SCT',
    layer = 'scale.data'
  )
  features <- rownames(dictionary$fingerprints$Lambdas)
  subset_features <- intersect(features, rownames(sct_data))
  dictionary$fingerprints$Lambdas <- dictionary$fingerprints$Lambdas[subset_features,]

  condition_metadata <- as.character(object@meta.data[,condition_meta])
  condition_metadata[condition_metadata==control_class] <- 'control'
  names(condition_metadata) <- rownames(object@meta.data)

  # Get residualized values and project
  top_credible_set_metadata <- 'top_credible_set'
  if (!is.null(suffix)) {
    top_credible_set_metadata <- paste0(top_credible_set_metadata,'.',suffix)
  }
  top_credible_set_metadata <- object@meta.data[,top_credible_set_metadata]
  residuals <- get_model_attribute(
    data = sct_data,
    query = colnames(sct_data)[which(top_credible_set_metadata != 'unassigned')],
    condition_meta = condition_metadata,
    fingerprints = dictionary$fingerprints,
    group_level = FALSE,
    attribute = 'residuals'
  )
  projected <- t(residuals) %*% dictionary$fingerprints[[1]]
  object[["proj"]] <- CreateAssayObject(data = t(projected))

  # Create UMAP reduction
  object <- RunUMAP(
    object,
    assay = 'proj',
    features = rownames(object[['proj']]),
    reduction.name = 'projUMAP',
    reduction.key = 'projUMAP_'
  )

  return(object)
}
