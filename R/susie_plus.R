#' @importFrom matrixStats colSds
#'
#' @keywords internal
compute_colSds = function(
    X
) {
  return(matrixStats::colSds(X))
}

#' @keywords internal
compute_colstats = function(
    X,
    center = TRUE,
    scale = TRUE
) {
  n = nrow(X)
  p = ncol(X)
  if (center)
    cm = colMeans(X,na.rm = TRUE)
  else
    cm = rep(0,p)
  if (scale) {
    csd = compute_colSds(X)
    csd[csd == 0] = 1
  } else
    csd = rep(1,p)

  d = n*colMeans(X)^2 + (n-1)*compute_colSds(X)^2
  d = (d - n*cm^2)/csd^2

  return(list(cm = cm,csd = csd,d = d))
}

#' @keywords internal
init_setup = function(
    n,
    p,
    L,
    S,
    scaled_prior_variance,
    residual_variance,
    prior_weights,
    null_weight,
    varY,
    standardize
) {
  if (!is.numeric(scaled_prior_variance) || scaled_prior_variance < 0)
    stop("Scaled prior variance should be positive number")
  if (scaled_prior_variance > 1 && standardize)
    stop("Scaled prior variance should be no greater than 1 when ",
         "standardize = TRUE")
  if(is.null(residual_variance))
    residual_variance = varY
  if(is.null(prior_weights)){
    prior_weights = rep(1/p,p)
  }else{
    if(all(prior_weights == 0)){
      stop("Prior weight should greater than 0 for at least one variable.")
    }
    prior_weights = prior_weights / sum(prior_weights)
  }
  if(length(prior_weights) != p)
    stop("Prior weights must have length p")
  if (p < L)
    L = p
  s = list(alpha  = matrix(1/p,nrow = L,ncol = p),
           mu     = matrix(0,nrow = L,ncol = p),
           mu2    = matrix(0,nrow = L,ncol = p),
           Xr     = rep(0,n),
           KL     = rep(as.numeric(NA),L),
           lbf    = rep(as.numeric(NA),L),
           lbf_variable = matrix(as.numeric(NA),L,p),
           sigma2 = residual_variance,
           V      = scaled_prior_variance*varY,
           pi     = prior_weights,
           mu_S   = if (!is.null(S)) { matrix(0,nrow = ncol(S),ncol = 1) } else { NA },
           S_KL   = if (!is.null(S)) { rep(as.numeric(NA),ncol(S)) } else { 0 },
           mu2_S  = if (!is.null(S)) { matrix(0,nrow = ncol(S),ncol = 1) } else { NA })
  if (is.null(null_weight))
    s$null_index = 0
  else
    s$null_index = p
  class(s) = "susie"
  return(s)
}

#' @keywords internal
compute_Xb = function(
    X,
    b
) {
  cm = attr(X,"scaled:center")
  csd = attr(X,"scaled:scale")

  # Scale Xb.
  scaled.Xb = tcrossprod(X,t(b/csd))

  # Center Xb.
  Xb = scaled.Xb - sum(cm*b/csd)
  return(as.numeric(Xb))
}

#' @keywords internal
init_finalize = function(
    s,
    X = NULL,
    Xr = NULL
) {
  if(length(s$V) == 1)
    s$V = rep(s$V,nrow(s$alpha))

  # Check sigma2.
  if (!is.numeric(s$sigma2))
    stop("Input residual variance sigma2 must be numeric")

  # Avoid problems with dimension if input is a 1 x 1 matrix.
  s$sigma2 = as.numeric(s$sigma2)
  if (length(s$sigma2) != 1)
    stop("Input residual variance sigma2 must be a scalar")
  if (s$sigma2 <= 0)
    stop("Residual variance sigma2 must be positive (is your var(Y) zero?)")

  # check prior variance
  if (!is.numeric(s$V))
    stop("Input prior variance must be numeric")
  if (!all(s$V >= 0))
    stop("prior variance must be non-negative")
  if (!all(dim(s$mu) == dim(s$mu2)))
    stop("dimension of mu and mu2 in input object do not match")
  if (!all(dim(s$mu) == dim(s$alpha)))
    stop("dimension of mu and alpha in input object do not match")
  if (nrow(s$alpha) != length(s$V))
    stop("Input prior variance V must have length of nrow of alpha in ",
         "input object")

  # Update Xr.
  if (!missing(Xr))
    s$Xr = Xr
  if (!missing(X))
    s$Xr = compute_Xb(X,colSums(s$mu * s$alpha))

  # Reset KL and lbf.
  s$KL = rep(as.numeric(NA),nrow(s$alpha))
  s$lbf = rep(as.numeric(NA),nrow(s$alpha))
  class(s) = "susie"
  return(s)
}

#' @keywords internal
compute_Xty = function(
    X,
    y
) {
  cm = attr(X,"scaled:center")
  csd = attr(X,"scaled:scale")
  ytX = crossprod(y,X)

  scaled.Xty = t(ytX/csd)

  # Center Xty.
  centered.scaled.Xty = scaled.Xty - cm/csd * sum(y)
  return(as.numeric(centered.scaled.Xty))
}

#' @importFrom stats dnorm
#'
#' @keywords internal
loglik = function(
    V,
    betahat,
    shat2,
    prior_weights
) {

  #log(bf) for each SNP
  lbf = dnorm(betahat,0,sqrt(V+shat2),log = TRUE) -
    dnorm(betahat,0,sqrt(shat2),log = TRUE)
  lpo = lbf + log(prior_weights + sqrt(.Machine$double.eps))

  # Deal with special case of infinite shat2 (e.g., happens if X does
  # not vary).
  lbf[is.infinite(shat2)] = 0
  lpo[is.infinite(shat2)] = 0

  maxlpo = max(lpo)
  w_weighted = exp(lpo - maxlpo)
  weighted_sum_w = sum(w_weighted)
  return(log(weighted_sum_w) + maxlpo)
}

#' @keywords internal
neg.loglik.logscale = function(
    lV,
    betahat,
    shat2,
    prior_weights
) {
  -loglik(exp(lV),betahat,shat2,prior_weights)
}

#' @importFrom stats optim
#'
#' @keywords internal
optimize_prior_variance = function(
    optimize_V,
    betahat,
    shat2,
    prior_weights,
    alpha = NULL,
    post_mean2 = NULL,
    V_init = NULL,
    check_null_threshold = 0
) {
  V = V_init
  lV = optim(par = log(max(c(betahat^2-shat2,1),na.rm = TRUE)),
             fn = neg.loglik.logscale,betahat = betahat,shat2 = shat2,
             prior_weights = prior_weights,method = "Brent",lower = -30,
             upper = 15)$par
  ## if the estimated one is worse than current one, don't change it.
  if(neg.loglik.logscale(lV, betahat = betahat,shat2 = shat2,prior_weights = prior_weights) >
     neg.loglik.logscale(log(V), betahat = betahat,
                         shat2 = shat2,prior_weights = prior_weights)){
    lV = log(V)
  }
  V = exp(lV)

  if (loglik(0,betahat,shat2,prior_weights) +
      check_null_threshold >= loglik(V,betahat,shat2,prior_weights))
    V = 0
  return(V)
}

#' @importFrom stats dnorm
#'
#' @keywords internal
single_effect_regression = function(
    y,
    X,
    V,
    residual_variance = 1,
    prior_weights = NULL,
    optimize_V = c("none", "optim", "uniroot", "EM", "simple"),
    check_null_threshold = 0
) {
    optimize_V = match.arg(optimize_V)
    Xty = compute_Xty(X,y)
    betahat = (1/attr(X,"d")) * Xty
    shat2 = residual_variance/attr(X,"d")
    if (is.null(prior_weights))
      prior_weights = rep(1/ncol(X),ncol(X))
    if (optimize_V != "EM" && optimize_V != "none")
      V = optimize_prior_variance(optimize_V,betahat,shat2,prior_weights,
                                  alpha = NULL,post_mean2 = NULL,V_init = V,
                                  check_null_threshold = check_null_threshold)

    # log(po) = log(BF * prior) for each SNP
    lbf = dnorm(betahat,0,sqrt(V + shat2),log = TRUE) -
      dnorm(betahat,0,sqrt(shat2),log = TRUE)
    lpo = lbf + log(prior_weights + sqrt(.Machine$double.eps))

    # Deal with special case of infinite shat2 (e.g., happens if X does
    # not vary).
    lbf[is.infinite(shat2)] = 0
    lpo[is.infinite(shat2)] = 0
    maxlpo = max(lpo)

    # w is proportional to
    #
    #   posterior odds = BF * prior,
    #
    # but subtract max for numerical stability.
    w_weighted = exp(lpo - maxlpo)

    # Posterior prob for each SNP.
    weighted_sum_w = sum(w_weighted)
    alpha = w_weighted / weighted_sum_w
    post_var = (1/V + attr(X,"d")/residual_variance)^(-1) # Posterior variance.
    post_mean = (1/residual_variance) * post_var * Xty
    post_mean2 = post_var + post_mean^2 # Second moment.

    # BF for single effect model.
    lbf_model = maxlpo + log(weighted_sum_w)
    loglik = lbf_model + sum(dnorm(y,0,sqrt(residual_variance),log = TRUE))

    return(list(alpha = alpha,mu = post_mean,mu2 = post_mean2,lbf = lbf,
                lbf_model = lbf_model,V = V,loglik = loglik))
}

#' @keywords internal
SER_posterior_e_loglik = function(
    X,
    Y,
    s2,
    Eb,
    Eb2
) {
  n = nrow(X)
  return(-0.5*n*log(2*pi*s2) - 0.5/s2*(sum(Y*Y)
                                       - 2*sum(Y*compute_Xb(X,Eb))
                                       + sum(attr(X,"d") * Eb2)))
}

#' @keywords internal
update_each_effect = function(
    X,
    y,
    s,
    S,
    varS,
    estimate_prior_variance = FALSE,
    estimate_prior_method = "optim",
    check_null_threshold
) {
  if (!estimate_prior_variance)
    estimate_prior_method = "none"

  # Repeat for each effect to update.
  L = nrow(s$alpha)
  if (L > 0) {
    for (l in 1:L) {

      # Remove lth effect from fitted values.
      s$Xr = s$Xr - compute_Xb(X,s$alpha[l,] * s$mu[l,])

      # Compute residuals.
      R = y - s$Xr

      res = single_effect_regression(R,X,s$V[l],s$sigma2,s$pi,
                                     estimate_prior_method,
                                     check_null_threshold)
      # Update the variational estimate of the posterior mean.
      s$mu[l,]    = res$mu
      s$alpha[l,] = res$alpha
      s$mu2[l,]   = res$mu2
      s$V[l]      = res$V
      s$lbf[l]    = res$lbf_model
      s$lbf_variable[l,] = res$lbf
      s$KL[l]     = -res$loglik +
        SER_posterior_e_loglik(X,R,s$sigma2,res$alpha * res$mu,
                               res$alpha * res$mu2)

      s$Xr = s$Xr + compute_Xb(X,s$alpha[l,] * s$mu[l,])
    }
    if (!is.null(S)) {
      for (l in (L+1):(L+ncol(S))) {

        s$Xr = s$Xr - S[,l-L,drop=F]%*%s$mu_S[l-L,,drop=F]
        R = y - s$Xr
        S_pick = S[,l-L,drop=F]
        out_S_pick = compute_colstats(S_pick,center = FALSE,scale = FALSE)
        attr(S_pick,"scaled:center") = out_S_pick$cm
        attr(S_pick,"scaled:scale") = out_S_pick$csd
        attr(S_pick,"d") = out_S_pick$d

        res = single_effect_regression(R,S_pick,varS[l-L],s$sigma2,1,
                                       "none",check_null_threshold)
        s$mu_S[l-L,1] = res$mu
        s$S_KL[l-L] = -res$loglik +
          SER_posterior_e_loglik(S_pick,R,s$sigma2,res$mu,
                                 res$mu2)
        s$mu2_S[l-L,1] = res$mu2
        s$Xr = s$Xr + S[,l-L,drop=F]%*%s$mu_S[l-L,,drop=F]

      }
    }
  }
  return(s)
}

#' @keywords internal
get_objective = function(
    X,
    Y,
    s,
    S
) {
  return(Eloglik(X,Y,s,S) - sum(s$KL) - sum(s$S_KL))
}

#' @keywords internal
Eloglik = function(
    X,
    Y,
    s,
    S
) {
  n = nrow(X)
  return(-(n/2) * log(2*pi*s$sigma2) - (1/(2*s$sigma2)) * get_ER2(X,Y,s,S))
}

#' @keywords internal
compute_MXt = function(
    M,
    X
) {
  cm = attr(X,"scaled:center")
  csd = attr(X,"scaled:scale")

  return(as.matrix(t(X %*% (t(M)/csd)) - drop(M %*% (cm/csd))))
}

#' @keywords internal
get_ER2 = function(
    X,
    Y,
    s,
    S
) {
  Xr_L = compute_MXt(s$alpha * s$mu,X) # L by N matrix
  Xr_S = if (!is.null(S)) { compute_MXt(t(s$mu_S),S) } else { 0 }
  postb2 = s$alpha * s$mu2 # Posterior second moment.
  if (!is.null(S)) {
    postb2_S <- t(s$mu2_S)
    S_term <- sum(attr(S,"d") * t(postb2_S))
  } else {
    S_term <- 0
  }
  return(sum((Y - s$Xr)^2) - sum(Xr_L^2) - sum(Xr_S^2) +
           sum(attr(X,"d") * t(postb2)) + S_term)
}

#' @keywords internal
estimate_residual_variance = function(
    X,
    y,
    s,
    S
) {
  n = nrow(X)
  return((1/n)*get_ER2(X,y,s,S))
}

#' @importFrom stats var
#'
#' @keywords internal
susie_plus = function(
    X,
    y,
    L = min(10,ncol(X)),
    S = NULL,
    varS = NULL,
    scaled_prior_variance = 0.2,
    residual_variance = NULL,
    prior_weights = NULL,
    null_weight = 0,
    standardize = TRUE,
    intercept = TRUE,
    estimate_residual_variance = TRUE,
    estimate_prior_variance = TRUE,
    estimate_prior_method = c("optim", "EM", "simple"),
    check_null_threshold = 0,
    prior_tol = 1e-9,
    residual_variance_upperbound = Inf,
    s_init = NULL,
    coverage = 0.95,
    min_abs_corr = 0.5,
    compute_univariate_zscore = FALSE,
    na.rm = FALSE,
    max_iter = 100,
    tol = 1e-3,
    verbose = FALSE,
    track_fit = FALSE,
    residual_variance_lowerbound = var(drop(y))/1e4,
    refine = FALSE,
    n_purity = 100
) {

  # Process input estimate_prior_method.
  estimate_prior_method = match.arg(estimate_prior_method)

  # Check input X.
  if (!(is.double(X) & is.matrix(X)) &
      !inherits(X,"CsparseMatrix") &
      is.null(attr(X,"matrix.type")))
    stop("Input X must be a double-precision matrix, or a sparse matrix, or ",
         "a trend filtering matrix")
  if (is.numeric(null_weight) && null_weight == 0)
    null_weight = NULL
  if (!is.null(null_weight) && is.null(attr(X,"matrix.type"))) {
    if (!is.numeric(null_weight))
      stop("Null weight must be numeric")
    if (null_weight < 0 || null_weight >= 1)
      stop("Null weight must be between 0 and 1")
    if (missing(prior_weights))
      prior_weights = c(rep(1/ncol(X) * (1 - null_weight),ncol(X)),null_weight)
    else
      prior_weights = c(prior_weights * (1-null_weight),null_weight)
    X = cbind(X,0)
  }
  if (anyNA(X))
    stop("Input X must not contain missing values")
  if (anyNA(y)) {
    if (na.rm) {
      samples_kept = which(!is.na(y))
      y = y[samples_kept]
      X = X[samples_kept,]
    } else
      stop("Input y must not contain missing values")
  }
  p = ncol(X)

  # Check input y.
  n = nrow(X)
  mean_y = mean(y)

  # Center and scale input.
  if (intercept)
    y = y - mean_y

  out = compute_colstats(X,center = intercept,scale = standardize)
  attr(X,"scaled:center") = out$cm
  attr(X,"scaled:scale") = out$csd
  attr(X,"d") = out$d

  if (!is.null(S)) {
    out_S = compute_colstats(S,center = intercept,scale = standardize)
    attr(S,"scaled:center") = out_S$cm
    attr(S,"scaled:scale") = out_S$csd
    attr(S,"d") = out_S$d
  }

  # Initialize susie fit.
  s = init_setup(n,p,L,S,scaled_prior_variance,residual_variance,prior_weights,
                 null_weight,as.numeric(var(y)),standardize)
  s = init_finalize(s)

  # Initialize elbo to NA.
  elbo = rep(as.numeric(NA),max_iter + 1)
  elbo[1] = -Inf;

  for (i in 1:max_iter) {
    s = update_each_effect(X,y,s,S,varS,estimate_prior_variance,estimate_prior_method,
                           check_null_threshold)

    # Compute objective before updating residual variance because part
    # of the objective s$kl has already been computed under the
    # residual variance before the update.
    elbo[i+1] = get_objective(X,y,s,S)
    if ((elbo[i+1] - elbo[i]) < tol) {
      s$converged = TRUE
      break
    }

    if (estimate_residual_variance) {
      s$sigma2 = pmax(residual_variance_lowerbound,
                      estimate_residual_variance(X,y,s,S))
      if (s$sigma2 > residual_variance_upperbound)
        s$sigma2 = residual_variance_upperbound
    }
  }

  # Remove first (infinite) entry, and trailing NAs.
  elbo = elbo[2:(i+1)]
  s$elbo = elbo
  s$niter = i

  if (is.null(s$converged)) {
    warning(paste("IBSS algorithm did not converge in",max_iter,"iterations!"))
    s$converged = FALSE
  }

  s$intercept = 0
  s$fitted = s$Xr
  s$fitted = drop(s$fitted)
  names(s$fitted) = `if`(is.null(names(y)),rownames(X),names(y))

  if (!is.null(colnames(X))) {
    variable_names = colnames(X)
    colnames(s$alpha) = variable_names
    colnames(s$mu)    = variable_names
    colnames(s$mu2)   = variable_names
    colnames(s$lbf_variable) = variable_names
  }

  # For prediction.
  s$X_column_scale_factors = attr(X,"scaled:scale")

  return(s)
}
