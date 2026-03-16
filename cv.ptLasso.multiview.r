#' Cross-validation for multistudy multiview ptLasso
#'
#' Run `ptLasso.multiview()` across a grid of transfer-learning `alpha` values
#' and select the best value using held-out performance.
#'
#' @param x A list of studies, each a named list of view matrices.
#' @param y A list of response vectors, one per study.
#' @param alphalist Numeric vector of transfer-learning alpha values to compare.
#' @param family Response family. Currently `"gaussian"` and `"binomial"` are supported.
#' @param type.measure Cross-validation metric used to compare alphas.
#' @param nfolds Number of folds used inside each `ptLasso.multiview()` fit.
#' @param foldid Optional stacked fold assignment across all studies.
#' @param s Lambda rule used when summarizing CV performance.
#' @param alphahat.choice Whether to choose `alphahat` using overall or mean performance.
#' @param verbose Should progress messages be printed?
#' @param fitoverall Optional pre-fit overall multiview model reused across alphas.
#' @param fitind Optional pre-fit list of individual multiview models reused across alphas.
#' @param group.intercepts Should study-specific stage-one baselines be used?
#' @param parallel Logical; if `TRUE`, allow study-level parallel fits where available.
#' @param ncores Number of worker cores for study-level parallel fits.
#' @param ... Additional arguments forwarded to `ptLasso.multiview()`.
#'
#' @return A `cv.ptLasso.multiview` object containing `alphahat`,
#'   `varying.alphahat`, `alphalist`, `errpre`, `errind`, `erroverall`,
#'   and the fitted model list in `fit`.
#' @export
cv.ptLasso.multiview <- function(
    x, y,
    alphalist = seq(0, 1, length = 11),
    family = c("gaussian", "binomial"),
    type.measure = c("default", "mse", "auc", "deviance"),
    nfolds = 10,
    foldid = NULL,
    s = c("lambda.min", "lambda.1se"),
    alphahat.choice = c("overall", "mean"),
    verbose = FALSE,
    fitoverall = NULL,
    fitind = NULL,
    group.intercepts = TRUE,
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
  s <- match.arg(s)
  alphahat.choice <- match.arg(alphahat.choice)

  if (family == "binomial" && !(type.measure %in% c("auc", "deviance"))) {
    stop("For binomial family, type.measure must be 'auc' or 'deviance'.")
  }
  if (family == "gaussian" && !(type.measure %in% c("mse", "deviance"))) {
    stop("For gaussian family, type.measure must be 'mse' or 'deviance'.")
  }

  alphalist <- sort(unique(as.numeric(alphalist)))
  if (length(alphalist) < 2L || any(is.na(alphalist)) || any(alphalist < 0 | alphalist > 1)) {
    stop("alphalist must contain at least two values between 0 and 1.")
  }

  metric_rule <- ptmv_match_metric(type.measure)
  fit <- vector("list", length(alphalist))
  err_rows <- vector("list", length(alphalist))

  for (ii in seq_along(alphalist)) {
    alpha <- alphalist[ii]
    if (verbose) {
      message(sprintf("alpha = %s", format(alpha)))
    }

    fit[[ii]] <- ptLasso.multiview(
      x = x,
      y = y,
      alpha = alpha,
      family = family,
      type.measure = type.measure,
      nfolds = nfolds,
      foldid = foldid,
      verbose = verbose,
      fitoverall = fitoverall,
      fitind = fitind,
      group.intercepts = group.intercepts,
      parallel = parallel,
      ncores = ncores,
      ...
    )

    if (is.null(fitoverall)) fitoverall <- fit[[ii]]$fitoverall
    if (is.null(fitind)) fitind <- fit[[ii]]$fitind

    pred_pre <- lapply(seq_along(fit[[ii]]$fitpre), function(kk) {
      model <- fit[[ii]]$fitpre[[kk]]
      lambda <- ptmv_resolve_s(model, s)
      lam_idx <- which.min(abs(model$lambda - lambda))
      as.numeric(model$fit.preval[, lam_idx])
    })
    names(pred_pre) <- fit[[ii]]$study_names

    err_row <- c(
      overall = ptmv_metric_value(
        unlist(fit[[ii]]$training_layout$y, use.names = FALSE),
        unlist(pred_pre, use.names = FALSE),
        family = family,
        type.measure = type.measure
      ),
      mean = mean(vapply(fit[[ii]]$study_names, function(study_name) {
        ptmv_metric_value(fit[[ii]]$training_layout$y[[study_name]], pred_pre[[study_name]], family, type.measure)
      }, numeric(1))),
      setNames(
        vapply(fit[[ii]]$study_names, function(study_name) {
          ptmv_metric_value(fit[[ii]]$training_layout$y[[study_name]], pred_pre[[study_name]], family, type.measure)
        }, numeric(1)),
        paste0("group_", fit[[ii]]$study_names)
      )
    )
    err_rows[[ii]] <- err_row
  }

  errpre <- cbind(alpha = alphalist, do.call(rbind, err_rows))
  rownames(errpre) <- NULL

  base_fit <- fit[[1]]
  overall_pred <- lapply(seq_along(base_fit$study_names), function(kk) {
    study_name <- base_fit$study_names[kk]
    lambda <- ptmv_resolve_s(base_fit$fitoverall, s)
    preds <- predict(
      base_fit$fitoverall,
      newx = base_fit$training_layout$x[[study_name]],
      s = lambda,
      type = "response",
      newoffset = if (base_fit$group.intercepts) rep(base_fit$group_baseline[study_name], base_fit$n_by_study[kk]) else NULL
    )
    as.numeric(preds)
  })
  names(overall_pred) <- base_fit$study_names

  ind_pred <- lapply(base_fit$study_names, function(study_name) {
    as.numeric(predict(base_fit$fitind[[study_name]], newx = base_fit$training_layout$x[[study_name]], s = s, type = "response"))
  })
  names(ind_pred) <- base_fit$study_names

  erroverall <- ptmv_summarize_metric(overall_pred, base_fit$training_layout$y, family, type.measure, add_r2 = family == "gaussian")
  errind <- ptmv_summarize_metric(ind_pred, base_fit$training_layout$y, family, type.measure, add_r2 = family == "gaussian")

  overall_idx <- if (alphahat.choice == "mean") metric_rule$best(errpre[, "mean"]) else metric_rule$best(errpre[, "overall"])
  alphahat <- alphalist[overall_idx]
  varying.alphahat <- vapply(base_fit$study_names, function(study_name) {
    alphalist[metric_rule$best(errpre[, paste0("group_", study_name)])]
  }, numeric(1))

  out <- list(
    call = this.call,
    alphahat = alphahat,
    varying.alphahat = varying.alphahat,
    alphalist = alphalist,
    errpre = errpre,
    errind = errind,
    erroverall = erroverall,
    fitoverall = fitoverall,
    fitind = fitind,
    fit = fit,
    family = family,
    type.measure = type.measure,
    s = s,
    alphahat.choice = alphahat.choice
  )
  class(out) <- "cv.ptLasso.multiview"
  out
}
