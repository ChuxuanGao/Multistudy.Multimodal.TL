#' Fit ptLasso for multistudy multiview data
#'
#' Train the multiview ptLasso workflow for study-list multiview input. The
#' wrapper fits one overall multiview model, one individual model per study,
#' and one pretrained model per study.
#'
#' @param x A list of studies, each a named list of view matrices.
#' @param y A list of response vectors, one per study.
#' @param alpha Transfer-learning alpha in `[0, 1]`.
#' @param family Response family. Currently `"gaussian"` and `"binomial"` are supported.
#' @param type.measure Cross-validation metric used inside the multiview fits.
#' @param overall.lambda Lambda rule used for the stage-one overall model.
#' @param ind.lambda Lambda rule used for the individual models.
#' @param pre.lambda Lambda rule used for the pretrained models.
#' @param foldid Optional stacked fold assignment across all studies.
#' @param nfolds Number of folds used when `foldid` is not supplied.
#' @param verbose Should progress messages be printed?
#' @param fitoverall Optional pre-fit overall multiview model.
#' @param fitind Optional pre-fit list of individual multiview models.
#' @param penalty.factor Optional penalty factors across concatenated views.
#' @param group.intercepts Should study-specific stage-one baselines be used?
#' @param en.alpha Elastic-net mixing parameter passed to the multiview base learner.
#' @param parallel Logical; if `TRUE`, allow study-level parallel fits where available.
#' @param ncores Number of worker cores for study-level parallel fits.
#' @param ... Additional arguments forwarded to the multiview base fitter.
#'
#' @return A `ptLasso.multiview` object containing metadata, fitted overall,
#'   individual, and pretrained models, plus offset and support information.
#' @export
ptLasso.multiview <- function(
    x, y,
    alpha = 0.5,
    family = c("gaussian", "binomial"),
    type.measure = c("default", "mse", "auc", "deviance"),
    overall.lambda = c("lambda.1se", "lambda.min"),
    ind.lambda = c("lambda.1se", "lambda.min"),
    pre.lambda = c("lambda.1se", "lambda.min"),
    foldid = NULL,
    nfolds = 10,
    verbose = FALSE,
    fitoverall = NULL,
    fitind = NULL,
    penalty.factor = NULL,
    group.intercepts = TRUE,
    en.alpha = 1,
    parallel = FALSE,
    ncores = 1L,
    ...
) {
  this.call <- match.call()

  family <- match.arg(family)
  type.measure <- match.arg(type.measure)
  if (type.measure == "default") {
    type.measure <- if (family == "gaussian") "mse" else "deviance"
  }

  overall.lambda <- match.arg(overall.lambda, c("lambda.1se", "lambda.min"))
  ind.lambda <- match.arg(ind.lambda, c("lambda.1se", "lambda.min"))
  pre.lambda <- match.arg(pre.lambda, c("lambda.1se", "lambda.min"))

  if (!is.numeric(alpha) || length(alpha) != 1L || alpha < 0 || alpha > 1) {
    stop("alpha must be a single number between 0 and 1.")
  }
  if (!is.numeric(en.alpha) || length(en.alpha) != 1L || en.alpha < 0 || en.alpha > 1) {
    stop("en.alpha must be a single number between 0 and 1.")
  }
  if (!is.numeric(nfolds) || length(nfolds) != 1L || nfolds < 2) {
    stop("nfolds must be a single integer greater than or equal to 2.")
  }

  family_fn <- switch(family, gaussian = gaussian, binomial = binomial)
  input <- ptmv_normalize_studies(x, y)
  x <- input$x
  y <- input$y
  k <- input$k
  study_names <- input$study_names
  view_names <- input$view_names
  n_views <- input$n_views
  n_by_study <- input$n_by_study
  N_all <- input$N_all
  p_by_view <- input$p_by_view
  p <- input$p
  x_list_all <- input$x_list_all
  y_all <- input$y_all
  groups_all <- input$groups_all

  if (k == 1L) {
    message("Single-study input detected; falling back to cvar.multiview().")
    return(cvar.multiview(
      x_list = x[[1]],
      y = y[[1]],
      family = family_fn(),
      alpha = en.alpha,
      s = overall.lambda,
      nfolds = min(nfolds, length(y[[1]])),
      foldid = if (is.null(foldid)) NULL else ptmv_renumber_foldid(foldid),
      penalty.factor = penalty.factor,
      type.measure = type.measure,
      keep = TRUE,
      ...
    ))
  }

  if (n_views == 1L) {
    message("Single-view multi-study input detected; falling back to ptLasso().")
    x_single <- x_list_all[[1]]
    if (is.null(penalty.factor)) {
      penalty.factor <- rep(1, ncol(x_single))
    }
    return(ptLasso::ptLasso(
      x = x_single,
      y = y_all,
      groups = groups_all,
      alpha = alpha,
      family = family,
      type.measure = type.measure,
      use.case = "inputGroups",
      overall.lambda = overall.lambda,
      foldid = foldid,
      nfolds = nfolds,
      verbose = verbose,
      penalty.factor = penalty.factor,
      fitoverall = fitoverall,
      fitind = fitind,
      en.alpha = en.alpha,
      group.intercepts = group.intercepts,
      parallel = parallel,
      ...
    ))
  }

  if (is.null(penalty.factor)) {
    penalty.factor <- rep(1, p)
  }
  if (length(penalty.factor) != p) {
    stop(sprintf("penalty.factor must have length %d, the total number of features across views.", p))
  }

  if (!is.null(fitoverall)) {
    valid_overall <- inherits(fitoverall, "cv.multiview") ||
      inherits(fitoverall, "cvar.multiview") ||
      inherits(fitoverall, "cv.multiview.revised")
    if (!valid_overall) {
      stop("fitoverall must be a cv.multiview/cvar.multiview/cv.multiview.revised object.")
    }
  }

  if (!is.null(fitind)) {
    if (length(fitind) != k) {
      stop("fitind must contain one model per study.")
    }
    valid_ind <- vapply(
      fitind,
      function(z) inherits(z, "cv.multiview") ||
        inherits(z, "cvar.multiview") ||
        inherits(z, "cv.multiview.revised"),
      logical(1)
    )
    if (!all(valid_ind)) {
      stop("All elements of fitind must be cv.multiview/cvar.multiview/cv.multiview.revised objects.")
    }
  }

  if (is.null(foldid)) {
    foldid_all <- integer(N_all)
    start <- 1L
    for (g in seq_len(k)) {
      idx <- start:(start + n_by_study[g] - 1L)
      foldid_all[idx] <- ptmv_make_foldid(
        n = n_by_study[g],
        nfolds = nfolds,
        family = family,
        y = y[[g]]
      )
      start <- start + n_by_study[g]
    }
    foldid_all <- ptmv_renumber_foldid(foldid_all)
  } else {
    if (length(foldid) != N_all) {
      stop("Provided foldid must have length equal to total N across studies.")
    }
    foldid_all <- ptmv_renumber_foldid(foldid)
  }
  nfolds_all <- length(unique(foldid_all))

  foldid_within <- split(foldid_all, groups_all)
  foldid_within <- lapply(foldid_within, ptmv_renumber_foldid)
  foldid_within <- foldid_within[study_names]

  group_baseline <- ptmv_compute_group_baseline(y, family)
  baseline_offset_all <- NULL
  if (isTRUE(group.intercepts)) {
    baseline_offset_all <- ptmv_compute_baseline_offset(
      y_all = y_all,
      groups_all = groups_all,
      foldid_all = foldid_all,
      family_obj = family_fn()
    )
  }

  if (is.null(fitoverall)) {
    if (verbose) {
      message("Fitting overall multiview model.")
    }
    fitoverall <- cvar.multiview(
      x_list = x_list_all,
      y = y_all,
      family = family_fn(),
      alpha = en.alpha,
      s = overall.lambda,
      type.measure = type.measure,
      foldid = foldid_all,
      nfolds = nfolds_all,
      offset = baseline_offset_all,
      penalty.factor = penalty.factor,
      keep = TRUE,
      ...
    )
  }

  fitoverall_fit <- fitoverall$multiview.fit
  lamhat <- fitoverall$lambda.choice
  rho_overall <- fitoverall$rho.choice
  r_idx <- which(fitoverall$rho == rho_overall)[1]
  rho_obj <- fitoverall$cv_by_rho[[r_idx]]
  lam_idx <- which.min(abs(rho_obj$lambda - lamhat))
  preval_all <- as.numeric(rho_obj$fit.preval[, lam_idx])
  preval.offset <- ptmv_split_vector_by_study(preval_all, n_by_study, study_names)

  coef_vec <- as.numeric(coef(fitoverall_fit, s = lamhat))
  supall <- which(coef_vec[-1] != 0)

  if (is.null(fitind)) {
    if (verbose) {
      message("Fitting individual multiview models.")
    }
    fitind <- ptmv_maybe_parallel_lapply(seq_len(k), function(kk) {
      if (verbose) {
        message(sprintf("  Study %d / %d", kk, k))
      }
      cvar.multiview(
        x_list = x[[kk]],
        y = y[[kk]],
        family = family_fn(),
        alpha = en.alpha,
        foldid = foldid_within[[kk]],
        s = ind.lambda,
        type.measure = type.measure,
        keep = TRUE,
        ...
      )
    }, parallel = parallel, ncores = ncores, verbose = verbose)
    names(fitind) <- study_names
  }

  if (verbose) {
    message("Fitting pretrained multiview models.")
  }

  if (alpha == 1) {
    fitpre <- fitind
  } else {
    fitpre <- ptmv_maybe_parallel_lapply(seq_len(k), function(kk) {
      if (verbose) {
        message(sprintf("  Pretrained model %d / %d", kk, k))
      }
      alpha_eff <- max(alpha, 1e-9)
      fac <- rep(1 / alpha_eff, p)
      fac[supall] <- 1
      pf <- penalty.factor * fac
      if (alpha == 0 && length(supall) == 0L) {
        pf <- penalty.factor * rep(1e9, p)
      }
      cvar.multiview(
        x_list = x[[kk]],
        y = y[[kk]],
        family = family_fn(),
        alpha = en.alpha,
        foldid = foldid_within[[kk]],
        s = pre.lambda,
        type.measure = type.measure,
        offset = (1 - alpha) * preval.offset[[kk]],
        penalty.factor = pf,
        keep = TRUE,
        ...
      )
    }, parallel = parallel, ncores = ncores, verbose = verbose)
    names(fitpre) <- study_names
  }

  out <- list(
    call = this.call,
    k = k,
    N_all = N_all,
    n_by_study = n_by_study,
    study_names = study_names,
    group.levels = study_names,
    view_names = view_names,
    n_views = n_views,
    p_by_view = p_by_view,
    features_all = p,
    alpha = alpha,
    en.alpha = en.alpha,
    family = family,
    type.measure = type.measure,
    overall.lambda = overall.lambda,
    ind.lambda = ind.lambda,
    pre.lambda = pre.lambda,
    fitoverall.lambda = lamhat,
    fitoverall.rho = rho_overall,
    group.intercepts = group.intercepts,
    group_baseline = group_baseline,
    foldid = foldid_all,
    foldid.within = foldid_within,
    support.vars = supall,
    penalty.factor = penalty.factor,
    fitoverall = fitoverall,
    fitind = fitind,
    fitpre = fitpre,
    baseline_offset = baseline_offset_all,
    preval.offset = preval.offset,
    training_layout = input,
    parallel = parallel,
    ncores = max(1L, as.integer(ncores))
  )

  class(out) <- "ptLasso.multiview"
  out
}