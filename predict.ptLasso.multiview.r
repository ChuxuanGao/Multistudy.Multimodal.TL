#' Predict from a multistudy multiview ptLasso fit
#'
#' Generate overall, individual, and pretrained predictions from a
#' `ptLasso.multiview` object on study-list multiview test data.
#'
#' @param object A fitted `ptLasso.multiview` object.
#' @param xtest A list of studies, each a named list of view matrices.
#' @param ytest Optional list of response vectors used to compute performance summaries.
#' @param type Prediction scale: `"link"`, `"response"`, or `"class"` for binomial fits.
#' @param s Lambda rule used for prediction.
#' @param return.link Should link-scale predictions also be returned?
#' @param ... Reserved for future extensions.
#'
#' @return A `predict.ptLasso.multiview` object containing study-wise predictions
#'   for the overall, individual, and pretrained models, support summaries,
#'   and optional error summaries.
#' @export
predict.ptLasso.multiview <- function(object, xtest, ytest = NULL,
                                      type = c("link", "response", "class"),
                                      s = c("lambda.min", "lambda.1se"),
                                      return.link = FALSE, ...) {
  if (missing(xtest)) {
    stop("Please supply xtest.")
  }
  if (!inherits(object, "ptLasso.multiview")) {
    stop("object must be a ptLasso.multiview fit.")
  }

  type <- match.arg(type)
  s <- match.arg(s)
  this.call <- match.call()

  if (object$family != "binomial" && type == "class") {
    stop("type = 'class' is only supported for binomial models.")
  }

  xtest <- ptmv_normalize_newdata(xtest, object$training_layout)
  ytest <- ptmv_prepare_ytest(ytest, object$study_names)
  test_sizes <- vapply(xtest, function(study) nrow(study[[1]]), integer(1))

  x_list_test <- ptmv_stack_by_view(xtest, object$view_names)
  overall_baseline_offset <- NULL
  if (isTRUE(object$group.intercepts)) {
    overall_baseline_offset <- ptmv_study_offsets(test_sizes, object$group_baseline[object$study_names])
  }

  overall_link_stacked <- as.numeric(predict(
    object$fitoverall,
    newx = x_list_test,
    s = s,
    type = "link",
    newoffset = overall_baseline_offset
  ))
  overall_resp_stacked <- if (type == "class") {
    as.numeric(predict(
      object$fitoverall,
      newx = x_list_test,
      s = s,
      type = "response",
      newoffset = overall_baseline_offset
    ))
  } else {
    as.numeric(predict(
      object$fitoverall,
      newx = x_list_test,
      s = s,
      type = type,
      newoffset = overall_baseline_offset
    ))
  }

  overall_link <- ptmv_split_vector_by_study(overall_link_stacked, test_sizes, object$study_names)
  overall_resp <- ptmv_split_vector_by_study(overall_resp_stacked, test_sizes, object$study_names)
  overall_pred <- if (type == "class") {
    lapply(overall_resp, function(x) ifelse(x >= 0.5, 1, 0))
  } else {
    overall_resp
  }

  stage1_link <- ptmv_split_vector_by_study(
    as.numeric(predict(
      object$fitoverall,
      newx = x_list_test,
      s = object$fitoverall.lambda,
      type = "link",
      newoffset = overall_baseline_offset
    )),
    test_sizes,
    object$study_names
  )

  pre_offsets <- lapply(stage1_link, function(x) (1 - object$alpha) * x)
  linkpre <- ptmv_predict_by_study(object$fitpre, xtest, s = s, type = "link", offsets = pre_offsets)
  yhatpre <- if (type == "class") {
    ptmv_predict_by_study(object$fitpre, xtest, s = s, type = "response", offsets = pre_offsets)
  } else {
    ptmv_predict_by_study(object$fitpre, xtest, s = s, type = type, offsets = pre_offsets)
  }
  if (type == "class") {
    yhatpre <- lapply(yhatpre, function(x) ifelse(x >= 0.5, 1, 0))
  }

  linkind <- ptmv_predict_by_study(object$fitind, xtest, s = s, type = "link")
  yhatind <- if (type == "class") {
    ptmv_predict_by_study(object$fitind, xtest, s = s, type = "response")
  } else {
    ptmv_predict_by_study(object$fitind, xtest, s = s, type = type)
  }
  if (type == "class") {
    yhatind <- lapply(yhatind, function(x) ifelse(x >= 0.5, 1, 0))
  }

  supoverall <- ptmv_get_support(object$fitoverall, s)
  supind <- ptmv_get_union_support(object$fitind, s)
  suppre.common <- ptmv_get_support(object$fitoverall, object$fitoverall.lambda)
  suppre.individual <- setdiff(ptmv_get_union_support(object$fitpre, s), suppre.common)

  erroverall <- errind <- errpre <- NULL
  if (!is.null(ytest)) {
    erroverall <- ptmv_summarize_metric(overall_resp, ytest, object$family, object$type.measure, add_r2 = object$family == "gaussian")
    errind <- ptmv_summarize_metric(yhatind, ytest, object$family, object$type.measure, add_r2 = object$family == "gaussian")
    errpre <- ptmv_summarize_metric(yhatpre, ytest, object$family, object$type.measure, add_r2 = object$family == "gaussian")
  }

  ptmv_build_prediction_object(
    call = this.call,
    alpha = object$alpha,
    type.measure = object$type.measure,
    yhatoverall = overall_pred,
    yhatind = yhatind,
    yhatpre = yhatpre,
    supoverall = supoverall,
    supind = supind,
    suppre.common = suppre.common,
    suppre.individual = suppre.individual,
    linkoverall = if (return.link) overall_link else NULL,
    linkind = if (return.link) linkind else NULL,
    linkpre = if (return.link) linkpre else NULL,
    erroverall = erroverall,
    errind = errind,
    errpre = errpre,
    class_name = "predict.ptLasso.multiview"
  )
}
