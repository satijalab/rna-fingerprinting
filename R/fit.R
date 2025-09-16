#' @importFrom irlba prcomp_irlba
#' @importFrom pbmcapply pbmclapply
#' @importFrom presto wilcoxauc
#' @importFrom stats cor mad median wilcox.test
#'
#' @keywords internal
fit <- function(
    data,
    perturbation_meta,
    features,
    min_cells = 10,
    override_filter = FALSE,
    mc.cores = 1
) {

  message('Processing data...')

  # Run PCA on control cells to capture baseline signal
  ctrl.pdev <- data[features,perturbation_meta=='control']
  ctrl.mus <- rowMeans(ctrl.pdev)
  PCs <- prcomp_irlba(t(ctrl.pdev),n=30)

  # Compute residuals for (a subset of) control cells after baseline subtraction
  if (ncol(ctrl.pdev)<1000) {
    samp <- 1:ncol(ctrl.pdev)
  } else {
    samp <- sample(1:ncol(ctrl.pdev),1000)
  }
  ctrl.residual <- ctrl.pdev[,samp]-sweep(PCs$rotation%*%t(PCs$x[samp,]),1,ctrl.mus,'+')

  # Identify perturbations to compute fingerprints for
  tab <- table(perturbation_meta[perturbation_meta!='control'])
  perturbations <- names(tab[tab>=min_cells])
  if (length(perturbations)==0) {
    message('No perturbations found with more than minimum number of cells')
    return(NULL)
  }
  if (length(perturbations)<length(tab)) {
    message(paste('Too few cells for perturbations: ',
                  paste(names(tab[tab<min_cells]),collapse=', ')))
  }

  message('Learning fingerprints...')

  # Learn Lambdas for each perturbation
  fits <- pbmclapply(perturbations,function(p) {

    # Perform baseline subtraction on cells that received perturbation p
    pdev <- data[features,perturbation_meta==p]
    reconstruct <- sweep(PCs$rotation%*%crossprod(PCs$rotation,sweep(pdev,1,ctrl.mus)),
                         1,ctrl.mus,'+')
    residual <- pdev-reconstruct

    # Learn initial Lambda and check for outlier cells
    PCs2 <- prcomp_irlba(t(residual),center=F,n=1)
    ell <- PCs2$x[,1]
    mads <- (ell-median(ell))/mad(ell,constant=1)
    outlier <- (max(abs(mads))>=10)

    # Iteratively discard outliers, if present, and re-compute Lambda each time
    pr <- wilcoxauc(cbind(ctrl.residual,residual),
                    y=c(rep('Control',ncol(ctrl.residual)),
                        rep('Perturbed',ncol(residual))))
    pr <- pr[pr$group=='Control',]
    track_cors <- c(cor(PCs2$rotation[,1],pr$statistic))
    store <- list(list(PCs2$rotation[,1],PCs2$x[,1],colnames(residual)))
    while (outlier) {
      pts <- which(abs(mads)>=10)
      cells_to_keep <- store[[length(store)]][[3]][-pts]
      residual_sub <- residual[,cells_to_keep]
      if (length(cells_to_keep)<min_cells) {
        break
      }
      PCs2 <- prcomp_irlba(t(residual_sub),center=F,n=1)
      ell <- PCs2$x[,1]
      mads <- (ell-median(ell))/mad(ell,constant=1)
      outlier <- (max(abs(mads))>=10)
      track_cors <- c(track_cors,cor(PCs2$rotation[,1],pr$statistic))
      store[[length(store)+1]] <- list(PCs2$rotation[,1],PCs2$x[,1],
                                       colnames(residual_sub))
    }

    # Choose the Lambda with maximum correlation to the Wilcox statistics
    pick <- which.max(abs(track_cors))
    Lam <- store[[pick]][[1]]
    ell <- store[[pick]][[2]]
    residual <- residual[,store[[pick]][[3]]]

    # Correct the sign of the Lambda, if needed
    if (track_cors[pick]>0) {
      Lam <- (-1)*Lam
      ell <- (-1)*ell
    }

    # Check if Lambda should be filtered out
    score <- wilcox.test(ell,t(Lam)%*%ctrl.residual)$p.value
    count <- sum(pr[pr$feature%in%features[order(abs(Lam),decreasing=T)[1:10]],
                    'pval']<0.05)
    keep <- ((score<0.05)&count>=2)

    # If Lambda will be kept, estimate its means and standard deviations
    if (keep | override_filter) {
      XtX <- sum(ell^2)
      XtY <- residual%*%ell
      means <- (XtY/XtX)
      sigma_squared <- colSums((t(residual)-t(means%*%ell))^2)/(length(ell)-1)
      sds <- sqrt((1/XtX)*sigma_squared)
      dists1 <- cbind(means[,1],sds)
      rownames(dists1) <- features
      return(list(Lam,dists1))
    } else {
      return(NULL)
    }
  },mc.cores=mc.cores)

  # Set perturbation names
  perturbation_names <- gsub("[,;_]","-",perturbations)
  if (!all(perturbation_names==perturbations)) {
    message('Replacing commas, underscores, and semicolons in perturbation names with dashes')
  }
  if (length(unique(perturbation_names))!=length(perturbations)) {
    perturbation_names <- make.unique(perturbation_names)
    message('Making perturbation names unique')
  }
  names(fits) <- perturbation_names

  # Put together non-filtered Lambdas and their distributions
  keep_perturbs <- names(fits)[sapply(fits,function(x) !is.null(x))]
  if (length(keep_perturbs)==0) {
    message('No perturbations with sufficiently high signal found; to learn all possible perturbations, even if low-signal, set override_filter=TRUE')
    return(NULL)
  }
  Lambdas <- do.call(cbind,lapply(keep_perturbs,function(x) fits[[x]][[1]]))
  dists <- lapply(keep_perturbs,function(x) fits[[x]][[2]])
  colnames(Lambdas) <- names(dists) <- keep_perturbs
  rownames(Lambdas) <- features

  message(paste('Learned',ncol(Lambdas),'of',length(perturbations),
                'possible perturbations; to learn all possible perturbations, even if low-signal, set override_filter=TRUE'))
  message(paste('Perturbations kept: ',
                paste(colnames(Lambdas),collapse=', ')))
  message(paste('Perturbations filtered out: ',
                paste(setdiff(perturbation_names,colnames(Lambdas)),collapse=', ')))

  return(list(Lambdas=Lambdas,dists=dists))
}

#' Learn dictionary
#'
#' Learn dictionary of fingerprints corresponding to each perturbation.
#'
#' @param object A Seurat object
#' @param perturbation_meta The name of the metadata column containing the perturbation identities
#' @param control_class The name of the perturbation identity corresponding to the control condition
#' @param features The features to be used for fingerprint estimation
#' @param min_cells The minimum number of cells that must be present to estimate a fingerprint (default 10)
#' @param override_filter If TRUE, low-signal fingerprints will not be filtered out (default FALSE)
#' @param mc.cores The number of cores to use for parallelization (default 1)
#' @param seed.use Set a random seed (default 1448145); setting NULL will not set a seed
#'
#' @return A list of fingerprints and their distributions for downstream fingerprinting analysis:
#' \describe{
#'    \item{Lambdas}{The matrix of fingerprints corresponding to each perturbation}
#'    \item{dists}{A list of means and standard deviations for each entry of each fingerprint}
#' }
#'
#' @import Seurat
#'
#' @export
#'
#' @examples
#' \dontrun{
#' dictionary <- LearnDictionary(
#' object,
#' perturbation_meta = 'perturbation_id',
#' control_class = 'NT',
#' features = features
#' )
#' }
#'
LearnDictionary <- function(
    object,
    perturbation_meta,
    control_class,
    features,
    min_cells = 10,
    override_filter = FALSE,
    mc.cores = 1,
    seed.use = 1448145
) {
  # Set up and verify parameters
  if (!is.null(x = seed.use)) {
    set.seed(seed = seed.use)
  }
  if (!"SCT" %in% Assays(object)) {
    stop("SCTransform() must be run before calling this function (recommended to set return.only.var.genes=FALSE). No 'SCT' assay found.")
  }
  sct_data <- GetAssayData(
    object,
    assay = 'SCT',
    layer = 'scale.data'
  )
  missing <- setdiff(features, rownames(sct_data))
  if (length(missing) > 0) {
    message(paste0("The following ", length(missing), " features are not present in the SCT assay:",
                   paste(head(missing, 50), collapse=', '), '...'))
    warning(paste0(length(missing), " features were not present in the SCT assay. If larger than expected, ensure SCTransform() was run with return.only.var.genes=FALSE."))
    features <- intersect(features, rownames(sct_data))
  }
  if (!perturbation_meta %in% colnames(object@meta.data)) {
    stop('Perturbation metadata column not found.')
  }
  perturbation_metadata <- as.character(object@meta.data[,perturbation_meta])
  if (sum(perturbation_metadata==control_class)==0) {
    stop('No cells with control class label found.')
  }
  perturbation_metadata[perturbation_metadata==control_class] <- 'control'

  # Learn fingerprints
  fingerprints <- fit(
    data = sct_data,
    perturbation_meta = perturbation_metadata,
    features = features,
    min_cells = min_cells,
    override_filter = override_filter,
    mc.cores = mc.cores
  )

  # Put into dictionary object
  dictionary <- structure(
    list(
      fingerprints = fingerprints,
      order = NULL,
      data = NULL
      ),
    class = "dictionary"
  )

  return(dictionary)
}
