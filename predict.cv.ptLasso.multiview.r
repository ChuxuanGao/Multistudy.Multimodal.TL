#' Predict from a cross-validated multistudy multiview ptLasso fit
#'
#' Resolve one fixed or varying `alpha` choice from a `cv.ptLasso.multiview`
#' object and generate predictions on study-list multiview test data.
#'
#' @param object A fitted `cv.ptLasso.multiview` object.
#' @param xtest A list of studies, each a named list of view matrices.
#' @param ytest Optional list of response vectors used to compute performance summaries.
#' @param alpha Optional user-specified alpha choice. May be one value or one per study.
#' @param alphatype Either `"fixed"` or `"varying"` when `alpha` is not supplied.
#' @param type Prediction scale: `"link"`, `"response"`, or `"class"` for binomial fits.
#' @param s Lambda rule used for prediction.
#' @param return.link Should link-scale predictions also be returned?
#' @param ... Reserved for future extensions.
#'
#' @return A `predict.cv.ptLasso.multiview` object containing the chosen alpha
#'   setting, study-wise predictions, support summaries, and optional error summaries.
#' @export
predict.cv.ptLasso.multiview <- function(object, xtest, ytest = NULL,
                                         alpha = NULL,
                                         alphatype = c("fixed", "varying"),
                                         type = c("link", "response", "class"),
                                         s = c("lambda.min", "lambda.1se"),
                                         return.link = FALSE, ...) {
  if (missing(xtest)) {
    stop("Please supply xtest.")
  }
  if (!inherits(object, "cv.ptLasso.multiview")) {
    stop("object must be a cv.ptLasso.multiview fit.")
  }

  this.call <- match.call()
  alphatype <- match.arg(alphatype)
  type <- match.arg(type)
  s <- match.arg(s)

  close.enough <- 1e-6
  if (is.null(alpha)) {
    alpha <- if (alphatype == "fixed") object$alphahat else object$varying.alphahat
  }

  if (length(alpha) == 1L) {
    model_idx <- which(abs(object$alphalist - alpha) < close.enough)
    if (length(model_idx) == 0L) {
      stop("Not a valid choice of alpha. Please choose alpha from object$alphalist.")
    }
    fit <- object$fit[[model_idx[1]]]
    out <- predict.ptLasso.multiview(
      fit,
      xtest = xtest,
      ytest = ytest,
      type = type,
      s = s,
      return.link = return.link,
      ...
    )
    out$call <- this.call
    out$fit <- object
    class(out) <- c("predict.cv.ptLasso.multiview", "predict.ptLasso.multiview")
    return(out)
  }

  if (is.null(names(alpha))) {
    if (length(alpha) != length(object$fit[[1]]$study_names)) {
      stop("Must have one alpha for each study.")
    }
    names(alpha) <- object$fit[[1]]$study_names
  }
  if (!setequal(names(alpha), object$fit[[1]]$study_names)) {
    stop("Alpha vector names must match the training study names exactly.")
  }
  alpha <- alpha[object$fit[[1]]$study_names]

  if (!all(vapply(alpha, function(a) any(abs(object$alphalist - a) < close.enough), logical(1)))) {
    stop("Includes at least one invalid choice of alpha. Please choose alpha from object$alphalist.")
  }

  pred_by_alpha <- lapply(object$fit, function(fit) {
    predict.ptLasso.multiview(
      fit,
      xtest = xtest,
      ytest = ytest,
      type = type,
      s = s,
      return.link = return.link,
      ...
    )
  })
  names(pred_by_alpha) <- as.character(object$alphalist)

  yhatpre <- yhatind <- yhatoverall <- vector("list", length(alpha))
  names(yhatpre) <- names(yhatind) <- names(yhatoverall) <- names(alpha)
  if (return.link) {
    linkpre <- linkind <- linkoverall <- vector("list", length(alpha))
    names(linkpre) <- names(linkind) <- names(linkoverall) <- names(alpha)
  } else {
    linkpre <- linkind <- linkoverall <- NULL
  }

  for (study_name in names(alpha)) {
    key <- as.character(alpha[[study_name]])
    pred <- pred_by_alpha[[key]]
    yhatpre[[study_name]] <- pred$yhatpre[[study_name]]
    yhatind[[study_name]] <- pred$yhatind[[study_name]]
    yhatoverall[[study_name]] <- pred$yhatoverall[[study_name]]
    if (return.link) {
      linkpre[[study_name]] <- pred$linkpre[[study_name]]
      linkind[[study_name]] <- pred$linkind[[study_name]]
      linkoverall[[study_name]] <- pred$linkoverall[[study_name]]
    }
  }

  supoverall <- pred_by_alpha[[as.character(object$alphahat)]]$supoverall
  supind <- sort(unique(unlist(lapply(names(alpha), function(study_name) {
    pred_by_alpha[[as.character(alpha[[study_name]])]]$supind
  }))))
  suppre.common <- pred_by_alpha[[as.character(object$alphahat)]]$suppre.common
  suppre.individual <- sort(unique(unlist(lapply(names(alpha), function(study_name) {
    pred_by_alpha[[as.character(alpha[[study_name]])]]$suppre.individual
  }))))

  erroverall <- errind <- errpre <- NULL
  if (!is.null(ytest)) {
    ytest_norm <- ptmv_prepare_ytest(ytest, names(alpha))
    family <- object$family
    type.measure <- object$type.measure
    erroverall <- ptmv_summarize_metric(yhatoverall, ytest_norm, family, type.measure, add_r2 = family == "gaussian")
    errind <- ptmv_summarize_metric(yhatind, ytest_norm, family, type.measure, add_r2 = family == "gaussian")
    errpre <- ptmv_summarize_metric(yhatpre, ytest_norm, family, type.measure, add_r2 = family == "gaussian")
  }

  ptmv_build_prediction_object(
    call = this.call,
    alpha = alpha,
    type.measure = object$type.measure,
    yhatoverall = yhatoverall,
    yhatind = yhatind,
    yhatpre = yhatpre,
    supoverall = supoverall,
    supind = supind,
    suppre.common = suppre.common,
    suppre.individual = setdiff(suppre.individual, suppre.common),
    linkoverall = linkoverall,
    linkind = linkind,
    linkpre = linkpre,
    erroverall = erroverall,
    errind = errind,
    errpre = errpre,
    fit = object,
    class_name = "predict.cv.ptLasso.multiview"
  )
}
