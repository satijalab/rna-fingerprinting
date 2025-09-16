#' @keywords internal
get_bf <- function(
    mod,
    Lambdas,
    Sxx,
    y,
    l
) {
  sig2 <- mod$sigma2
  sig02 <- mod$V[l]

  Sxy <- crossprod(Lambdas,y)
  s21 <- sig2/Sxx
  bhat1 <- (1/Sxx)*Sxy
  z1 <- bhat1/sqrt(s21)
  as.numeric(log(sqrt(s21/(s21+sig02)))+(z1^2/2)*sig02/(sig02+s21))
}

#' @keywords internal
get_credible_sets <- function(
    alpha,
    store_cor,
    coverage
) {
  cs_sets <- list()
  for (l in 1:nrow(alpha)) {
    cs <- cumsum(alpha[l,][order(alpha[l,],decreasing=T)])
    cs2 <- names(cs)[1:as.numeric(which(cs>=coverage)[1])]
    cor <- store_cor[cs2,cs2,drop=F]
    if (nrow(cor)>1) {
      diag(cor) <- NA
    }
    if (median(abs(cor),na.rm=T)>=0.4) {
      cs_sets <- c(cs_sets,list(cs2))
      names(cs_sets)[[length(cs_sets)]] <- rownames(alpha)[l]
    }
  }
  return(cs_sets)
}

#' @importFrom stats rnorm
#'
#' @keywords internal
calibrate <- function(
    res0,
    mod,
    keep_cs,
    Lam,
    filter,
    newLambdas,
    Sxx,
    hold=NULL
) {
  precompute <- Lam%*%t(mod$alpha*mod$mu)
  alpha <- lapply(1:50,function(i) {
    newL <- newLambdas[[i]][,filter]
    newSxx <- Sxx[[i]][filter]
    if (!is.null(hold)) {
      for (h in 1:length(hold)) {
        dist <- hold[[h]]
        draw <- rnorm(ncol(dist),dist[1,],dist[2,])
        newL <- cbind(draw,newL)
        colnames(newL)[1] <- names(hold)[h]
      }
      newSxx <- colSums(newL*newL)
    }
    newL <- newL[,colnames(Lam)]
    newSxx <- newSxx[colnames(Lam)]
    t(sapply(keep_cs,function(l) {
      res <- res0-(mod$Xr-precompute[,l])
      bf <- get_bf(mod,newL,newSxx,res,l)
      bf <- bf-max(bf)
      exp(bf)/sum(exp(bf))
    }))
  })
  alpha <- Reduce("+", alpha) / length(alpha)
  alpha <- array(alpha,dim=c(length(keep_cs),ncol(Lam)))
  colnames(alpha) <- colnames(Lam)
  rownames(alpha) <- keep_cs
  return(alpha)
}

#' @importFrom stats setNames
#'
#' @keywords internal
post_process_cs <- function(
    cs_sets,
    mu,
    lbf
) {
  signs <- lapply(names(cs_sets),function(x) {
    mu_sign <- ifelse(sign(mu[x,cs_sets[[x]],drop=F])==1,'+','-')
    mu_sign_split <- unlist(lapply(colnames(mu_sign),function(name) {
      split_names <- unlist(strsplit(name,'_'))
      setNames(rep(mu_sign[1,name],length(split_names)),split_names)
    }))
    mu_sign_split[unique(names(mu_sign_split))]
  })

  cs_sets <- lapply(signs,function(x) names(x))
  final_cs_sets <- final_signs <- list()
  seen <- final_lbf <- NULL
  for (i in 1:length(cs_sets)) {
    cs <- cs_sets[[i]]
    if (!any(cs %in% seen)) {
      final_cs_sets[[length(final_cs_sets)+1]] <- cs
      final_signs[[length(final_signs)+1]] <- unname(signs[[i]])
      final_lbf <- c(final_lbf,lbf[i])
      seen <- c(seen,cs)
    }
  }

  final_cs_sets_str <- paste(sapply(final_cs_sets,function(x)
    paste(x,sep=',',collapse=',')),sep=';',collapse=';')
  final_signs_str <- paste(sapply(final_signs,function(x)
    paste(x,sep=',',collapse=',')),sep=';',collapse=';')
  return(list(final_cs_sets_str,unname(final_lbf),final_signs_str))
}

#' @importFrom stats rnorm qnorm
#'
#' @keywords internal
classify <- function(
    data,
    condition_meta,
    fingerprints,
    group_level = FALSE,
    top = 100,
    mc.cores = 1
) {

  message('Processing data...')

  # Extract Lambdas, distributions from fingerprints
  Lambdas <- fingerprints$Lambdas
  dists <- fingerprints$dists
  if (is.null(top)) {
    top <- ncol(Lambdas)
  }
  top <- min(c(top,ncol(Lambdas)))
  features <- rownames(Lambdas)

  # Run PCA on control cells to capture baseline signal
  ctrl.pdev <- data[features,condition_meta=='control']
  ctrl.mus <- rowMeans(ctrl.pdev)
  PCs <- prcomp_irlba(t(ctrl.pdev),n=30)

  # Perform baseline subtraction on all perturbed cells
  w.pdev <- data[features,condition_meta!='control']
  w.reconstruct <- sweep(PCs$rotation%*%crossprod(PCs$rotation,sweep(w.pdev,1,ctrl.mus)),
                         1,ctrl.mus,'+')
  w.residual <- w.pdev-w.reconstruct

  # If performing group-level analysis, aggregate cells by perturbation ID
  if (group_level) {
    labels <- condition_meta[colnames(w.residual)]
    grouped_indices <- split(1:ncol(w.residual),labels)
    w.residual <- sapply(unique(labels), function(id) {
      cols <- grouped_indices[[id]]
      rowSums(w.residual[, cols, drop = FALSE])
    })
  }

  message('Preparing...')

  # Pre-compute quantities needed for calibration and credible set identification
  newLambdas <- lapply(1:50,function(x) {
    newL <- sapply(1:ncol(Lambdas),function(k) {
      rnorm(length(features),dists[[colnames(Lambdas)[k]]][,1],
            dists[[colnames(Lambdas)[k]]][,2])
    })
    colnames(newL) <- colnames(Lambdas)
    newL
  })
  Sxx <- lapply(1:50,function(x) {
    newSxx <- colSums(newLambdas[[x]]*newLambdas[[x]])
    names(newSxx) <- colnames(Lambdas)
    newSxx
  })
  store_cor <- cor(Lambdas)

  message('Classifying...')

  # Perform classification for each cell or group
  sus <- pbmclapply(1:ncol(w.residual),function(x) {

    # Limit regression to top most correlated Lambdas to control cost
    crs <- c(abs(cor(w.residual[,x],Lambdas)))
    filter <- colnames(Lambdas)[order(crs,decreasing=T)[1:top]]

    # If single-cell classification, set up dummy variables for outliers
    if (!group_level) {
      top_genes <- order(abs(w.residual[,x]),decreasing=T)[1:5]
      top_genes <- unique(c(top_genes,which(abs(w.residual[,x])>=10)))
      dummy <- array(0,dim=c(nrow(Lambdas),length(top_genes)))
      for (col in 1:ncol(dummy)) {
        dummy[top_genes[col],col] <- 1
      }
      dummy_vars <- (0.5*abs(w.residual[top_genes,x])/qnorm(.95))^2
    } else {
      dummy <- dummy_vars <- NULL
    }

    # Run initial SuSiE regression and identify credible sets
    mod <- susie_plus(X=Lambdas[,filter],y=w.residual[,x],S=dummy,varS=dummy_vars,
                      standardize=F,L=10,intercept=F,max_iter=5000,
                      estimate_prior_method='optim')
    rownames(mod$alpha) <- rownames(mod$mu) <- names(mod$lbf) <-
      names(mod$V) <- paste0('cs',1:nrow(mod$alpha))
    keep_cs <- names(mod$lbf)[which(mod$lbf>=0.1*max(mod$lbf) & mod$V>0)]
    if (length(keep_cs)>0) {
      cs_sets <- get_credible_sets(mod$alpha[keep_cs,,drop=F],store_cor,0.95)
    } else {
      cs_sets <- list()
    }

    # If no credible sets with sufficient evidence, return 'unassigned'; else proceed
    if (length(cs_sets)==0) {
      list('unassigned')
    } else {

      # Calibrate SuSiE PIPs and re-compute credible sets
      alpha <- calibrate(w.residual[,x],mod,keep_cs,Lambdas[,filter],
                         filter,newLambdas,Sxx)
      cs_sets <- get_credible_sets(alpha,store_cor,0.95)

      # Assess if any joint Lambdas should be proposed
      prt_sets <- list()
      keep_prt <- c()
      lbf <- mod$lbf[names(cs_sets)]
      if (length(cs_sets)>1) {
        cs_sets <- cs_sets[order(lbf,decreasing=T)]
        for (cs in cs_sets) {
          if (length(prt_sets)==0) {
            prt_sets[[1]] <- cs
          } else {
            unplaced <- TRUE
            cnter <- 0
            while (unplaced & cnter < length(prt_sets)) {
              cnter <- cnter+1
              cr <- median(abs(store_cor[prt_sets[[cnter]],cs]))
              if (cr>=0.1) {
                prt_sets[[cnter]] <- c(prt_sets[[cnter]],cs)
                keep_prt <- c(keep_prt,cnter)
                unplaced <- FALSE
              }
            }
            if (unplaced) {
              prt_sets[[length(prt_sets)+1]] <- cs
            }
          }
        }
        keep_prt <- unique(keep_prt)
        keep_prt <- keep_prt[sapply(keep_prt,function(f) {
          length(unique(prt_sets[[f]]))>1
        })]
        prt_sets <- lapply(keep_prt,function(f) unique(prt_sets[[f]]))
      }

      # If no joint Lambdas to propose, return answer; else proceed
      if (length(cs_sets)==0) {
        return(list('unassigned'))
      } else if (length(prt_sets)==0) {
        return(post_process_cs(cs_sets,mod$mu[names(cs_sets),,drop=F],
                               lbf[names(cs_sets)]))
      } else {

        # Create joint Lambda
        hold <- list()
        Lam2 <- Lambdas[,filter]
        for (prt in prt_sets) {
          Lj <- Lambdas[,prt]
          lamj <- rowMeans(Lj)
          sj <- sqrt(rowSums(sapply(prt,function(x) dists[[x]][,2]^2)))/length(prt)
          lam <- rbind(lamj,sj)
          hold[[length(hold)+1]] <- lam
          Lam2 <- cbind(lam[1,],Lam2)
          colnames(Lam2)[1] <- names(hold)[[length(hold)]] <- paste(prt,collapse='_')
        }

        # Run SuSiE regression with joint Lambda and calibrate
        mod <- susie_plus(X=Lam2,y=w.residual[,x],S=dummy,varS=dummy_vars,
                          standardize=F,L=10,intercept=F,max_iter=5000,
                          estimate_prior_method='optim')
        rownames(mod$alpha) <- rownames(mod$mu) <- names(mod$lbf) <-
          names(mod$V) <- paste0('cs',1:nrow(mod$alpha))
        keep_cs <- names(mod$lbf)[which(mod$lbf>=0.1*max(mod$lbf) & mod$V>0)]
        alpha <- calibrate(w.residual[,x],mod,keep_cs,Lam2,
                           filter,newLambdas,Sxx,hold)

        # Compute credible sets
        store_cor2 <- array(0,dim=rep(length(filter)+length(prt_sets),2))
        colnames(store_cor2) <- rownames(store_cor2) <- colnames(Lam2)
        store_cor2[colnames(Lambdas[,filter]),colnames(Lambdas[,filter])] <-
          store_cor[colnames(Lambdas[,filter]),colnames(Lambdas[,filter])]
        extra_cor <- cor(Lam2[,1:length(prt_sets),drop=F],Lam2)
        store_cor2[rownames(extra_cor),colnames(extra_cor)] <- extra_cor
        store_cor2[colnames(extra_cor),rownames(extra_cor)] <- t(extra_cor)
        cs_sets <- get_credible_sets(alpha,store_cor2,0.95)
        lbf <- mod$lbf[names(cs_sets)]

        # Return answer
        if(length(cs_sets)==0) {
          return(list('unassigned'))
        } else {
          cs_sets <- cs_sets[order(lbf,decreasing=T)]
          return(post_process_cs(cs_sets,mod$mu[names(cs_sets),,drop=F],
                                 lbf[names(cs_sets)]))
        }
      }
    }
  },mc.cores=mc.cores)

  # Extract quantities to return
  assign <- sapply(sus,function(x) x[[1]])
  first <- sub(';.*$','',assign)
  maxp <- sub(',.*$','',first)
  full_lbf <- lapply(1:length(sus),function(x) {
    if (first[x]=='unassigned') { NA } else { sus[[x]][[2]] }
  })
  top_lbf <- sapply(full_lbf,function(x) round(x[1],2))
  full_lbf <- sapply(full_lbf,function(x) paste(round(x,2),sep=';',collapse=';'))
  full_signs <- sapply(1:length(sus),function(x) {
    if (first[x]=='unassigned') { NA } else { sus[[x]][[3]] }
  })
  names(maxp) <- names(first) <- names(assign) <- colnames(w.residual)
  names(full_lbf) <- names(top_lbf) <- names(full_signs) <- colnames(w.residual)
  return(list(top_perturbation=maxp,top_credible_set=first,
              full_assignment=assign,top_lbf=top_lbf,full_lbf=full_lbf,
              full_signs=full_signs))
}

#' Fingerprint cells
#'
#' Map cells to fingerprints.
#'
#' @param object A Seurat object
#' @param condition_meta The name of the metadata column containing the condition identities
#' @param control_class The name of the condition identity corresponding to the control condition
#' @param dictionary The dictionary of fingeprrints to map the data to
#' @param group_level Whether classification should be done for each cell, or for each condition in aggregate
#' @param top To speed up classification, maximum number of top fingerprints to include in regression (default 100); set NULL to use all
#' @param mc.cores The number of cores to use for parallelization (default 1)
#' @param seed.use Set a random seed (default 1448145); setting NULL will not set a seed
#' @param suffix Suffix for metadata columns; defaults to NULL under cell-level classification or the specified condition name if group-level
#'
#' @return A Seurat object with the following new metadata columns, appended by the specified suffix if any:
#' \describe{
#'    \item{top_credible_set}{The top credible set assigned to each cell (or the condition to which the cell belongs).
#'    Distinct members are separated by commas.}
#'    \item{top_lbf}{The log Bayes factor associated with the top credible set.}
#'    \item{all_credible_sets}{All credible sets assigned to each cell (or the condition to which the cell belongs).
#'    Distinct sets are separated by semicolons, and distinct members within a credible set are separated by commas.}
#'    \item{all_lbfs}{The log Bayes factors associated with each credible set, separated by semicolons.}
#'    \item{all_signs}{The signs of association for each member (separated by commas) of each credible set (separated by semicolons).}
#' }
#'
#' @import Seurat
#'
#' @export
#'
#' @examples
#' \dontrun{
#' object <- FingerprintCells(
#' object,
#' condition_meta = 'condition_id',
#' control_class = 'untreated',
#' dictionary = dictionary,
#' group_level = FALSE
#' )
#'
#' object <- FingerprintCells(
#' object,
#' condition_meta = 'condition_id',
#' control_class = 'untreated',
#' dictionary = dictionary,
#' group_level = TRUE
#' )
#' }
#'
FingerprintCells <- function(
    object,
    condition_meta,
    control_class,
    dictionary,
    group_level,
    top = 100,
    mc.cores = 1,
    seed.use = 1448145,
    suffix = NULL
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
  fingerprints <- dictionary$fingerprints
  features <- rownames(fingerprints$Lambdas)
  missing <- setdiff(features, rownames(sct_data))
  if (length(missing) > 0) {
    message(paste0("The following ", length(missing), " features are not present in the SCT assay:",
                paste(head(missing, 50), collapse=', '), '...'))
    if (length(missing) > 0.95*length(features)) {
      stop("More than 95% of features are not present. If the fingerprints were learned with a different naming system, use ConvertGenes to update them.")
    }
    if (length(missing) > 0.25*length(features)) {
      warning("More than 25% of features are not present. Ensure SCTransform() was run with return.only.var.genes=FALSE or consider re-fitting fingerprints without missing features.")
    }

    subset_features <- intersect(features, rownames(sct_data))
    fingerprints$Lambdas <- fingerprints$Lambdas[subset_features,]
    for (p in 1:length(fingerprints$dists)) {
      fingerprints$dists[[p]] <- fingerprints$dists[[p]][subset_features,]
    }
  }
  if (!condition_meta %in% colnames(object@meta.data)) {
    stop('Condition metadata column not found.')
  }
  condition_metadata <- as.character(object@meta.data[,condition_meta])
  if (sum(condition_metadata==control_class)<30) {
    stop('Not enough cells with control class label found.')
  }
  condition_metadata[condition_metadata==control_class] <- 'control'
  names(condition_metadata) <- rownames(object@meta.data)
  if (is.null(suffix) & group_level == TRUE) {
    suffix <- condition_meta
  }

  # Map cells
  out <- classify(
    data = sct_data,
    condition_meta = condition_metadata,
    fingerprints = fingerprints,
    group_level = group_level,
    top = top,
    mc.cores = mc.cores
  )

  # Assign to metadata
  if (!group_level) {
    new_meta <- data.frame(
      top_credible_set = out$top_credible_set,
      top_lbf = out$top_lbf,
      all_credible_sets = out$full_assignment,
      all_lbfs = out$full_lbf,
      all_signs = out$full_signs
    )
    if (!is.null(suffix)) {
      colnames(new_meta) <- paste0(colnames(new_meta),'.',suffix)
    }
    object <- AddMetaData(object, new_meta)
  } else {
    propagate <- function(attr) { ifelse(condition_metadata != 'control',
                                         out[[attr]][condition_metadata],
                                         NA) }
    new_meta <- data.frame(
      top_credible_set = propagate('top_credible_set'),
      top_lbf = propagate('top_lbf'),
      all_credible_sets = propagate('full_assignment'),
      all_lbfs = propagate('full_lbf'),
      all_signs = propagate('full_signs')
    )
    if (!is.null(suffix)) {
      colnames(new_meta) <- paste0(colnames(new_meta),'.',suffix)
    }
    object <- AddMetaData(object, new_meta)
  }

  return(object)
}
