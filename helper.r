# Helper: choose whether a metric should be minimized or maximized.
# Used in `cv.ptLasso.multiview()` when selecting `alphahat`.
ptmv_match_metric <- function(type.measure) {
  if (type.measure %in% c("auc")) {
    list(best = which.max, aggregate = max)
  } else {
    list(best = which.min, aggregate = min)
  }
}

# Helper: resolve `s` into a numeric lambda value for a cv-style multiview fit.
# Used in `cv.ptLasso.multiview()`, `predict.ptLasso.multiview()`,
# and `predict.cv.ptLasso.multiview()`.
ptmv_resolve_s <- function(object, s = c("lambda.min", "lambda.1se")) {
  if (is.numeric(s)) {
    return(as.numeric(s)[1])
  }
  s <- match.arg(s)
  object[[s]]
}

# Helper: stable logit transform used for smoothed binomial baselines.
# Used in `ptmv_compute_group_baseline()` and `ptmv_compute_baseline_offset()`.
ptmv_safe_logit <- function(p) log(p / (1 - p))

# Helper: compute one baseline value per study for future prediction.
# Used in `ptLasso.multiview()` and `predict.ptLasso.multiview()`.
ptmv_compute_group_baseline <- function(y, family) {
  if (family == "gaussian") {
    return(vapply(y, mean, numeric(1)))
  }

  if (family == "binomial") {
    a <- 0.5
    b <- 0.5
    return(vapply(y, function(yi) {
      yi <- as.numeric(yi)
      n1 <- sum(yi == 1)
      n0 <- sum(yi == 0)
      ptmv_safe_logit((n1 + a) / (n1 + n0 + a + b))
    }, numeric(1)))
  }

  stop("group baseline only implemented for gaussian and binomial.")
}

# Helper: compute fold-wise study-specific offsets for stage-one training.
# Used in `ptLasso.multiview()`.
ptmv_compute_baseline_offset <- function(y_all, groups_all, foldid_all, family_obj) {
  off <- numeric(length(y_all))
  folds <- sort(unique(foldid_all))
  grps <- levels(groups_all)

  if (family_obj$family == "gaussian") {
    for (f in folds) {
      train <- foldid_all != f
      for (g in grps) {
        idx <- which(groups_all == g & !train)
        if (length(idx) == 0L) {
          next
        }
        off[idx] <- mean(y_all[train & groups_all == g])
      }
    }
    return(off)
  }

  if (family_obj$family == "binomial") {
    a <- 0.5
    b <- 0.5
    for (f in folds) {
      train <- foldid_all != f
      for (g in grps) {
        idx <- which(groups_all == g & !train)
        if (length(idx) == 0L) {
          next
        }
        y_tr_g <- as.numeric(y_all[train & groups_all == g])
        n1 <- sum(y_tr_g == 1)
        n0 <- sum(y_tr_g == 0)
        off[idx] <- ptmv_safe_logit((n1 + a) / (n1 + n0 + a + b))
      }
    }
    return(off)
  }

  stop("Baseline offset only implemented for gaussian and binomial families.")
}

# Helper: remap fold ids to consecutive integers.
# Used in `ptLasso.multiview()` and fallback calls.
ptmv_renumber_foldid <- function(fid) {
  used <- sort(unique(as.integer(fid)))
  fold_map <- setNames(seq_along(used), used)
  as.integer(fold_map[as.character(as.integer(fid))])
}

# Helper: create study-level folds with stratification for binomial outcomes.
# Used in `ptLasso.multiview()`.
ptmv_make_foldid <- function(n, nfolds, family, y = NULL) {
  nfolds_use <- min(nfolds, n)
  if (family == "binomial") {
    y01 <- as.integer(as.numeric(y) > 0)
    idx0 <- which(y01 == 0L)
    idx1 <- which(y01 == 1L)
    foldid <- integer(n)
    foldid[idx0] <- sample(rep(seq_len(nfolds_use), length.out = length(idx0)))
    foldid[idx1] <- sample(rep(seq_len(nfolds_use), length.out = length(idx1)))
    return(foldid)
  }
  sample(rep(seq_len(nfolds_use), length.out = n))
}

# Helper: switch between serial and study-level parallel lapply.
# Used in `ptLasso.multiview()`.
ptmv_maybe_parallel_lapply <- function(X, FUN, parallel = FALSE, ncores = 1L,
                                       verbose = FALSE, ...) {
  if (!isTRUE(parallel) || length(X) <= 1L) {
    return(lapply(X, FUN, ...))
  }

  ncores_use <- max(1L, as.integer(ncores))
  if (.Platform$OS.type == "windows" || ncores_use == 1L) {
    if (verbose) {
      message("parallel = TRUE requested, but using serial lapply for this platform/configuration.")
    }
    return(lapply(X, FUN, ...))
  }

  parallel::mclapply(X, FUN, ..., mc.cores = ncores_use)
}

# Helper: normalize and validate study-list training input.
# Used in `ptLasso.multiview()`.
ptmv_normalize_studies <- function(x, y) {
  if (!is.list(x) || length(x) == 0L) {
    stop("x must be a non-empty list of studies.")
  }
  if (!is.list(y) || length(y) != length(x)) {
    stop("y must be a list with one response vector per study.")
  }

  k <- length(x)
  study_names <- names(x)
  if (is.null(study_names) || any(study_names == "")) {
    study_names <- paste0("Study_", seq_len(k))
  }

  normalize_one_study <- function(study_obj, study_name) {
    if (is.matrix(study_obj)) {
      out <- list(View_1 = study_obj)
    } else if (is.list(study_obj) && length(study_obj) > 0L) {
      out <- study_obj
    } else {
      stop(sprintf("'%s' must be a matrix or a non-empty list of view matrices.", study_name))
    }

    if (is.null(names(out)) || any(names(out) == "")) {
      names(out) <- paste0("View_", seq_along(out))
    }

    for (v in names(out)) {
      if (!is.matrix(out[[v]])) {
        stop(sprintf("'%s' -> '%s' must be a matrix.", study_name, v))
      }
    }

    out
  }

  x_norm <- Map(normalize_one_study, x, study_names)
  names(x_norm) <- study_names

  y_norm <- lapply(seq_along(y), function(i) {
    yi <- y[[i]]
    if (is.null(dim(yi))) {
      yi <- drop(yi)
    }
    if (length(yi) == 0L) {
      stop(sprintf("'%s' has an empty response vector.", study_names[i]))
    }
    yi
  })
  names(y_norm) <- study_names

  view_names <- names(x_norm[[1]])
  n_views <- length(view_names)
  p_by_view <- setNames(integer(n_views), view_names)
  ref_cols <- vector("list", n_views)
  names(ref_cols) <- view_names

  for (i in seq_along(x_norm)) {
    xi <- x_norm[[i]]
    yi <- y_norm[[i]]
    if (!identical(names(xi), view_names)) {
      stop(sprintf("All studies must share identical view names and ordering. Mismatch found in '%s'.", study_names[i]))
    }

    for (v in view_names) {
      xv <- xi[[v]]
      if (nrow(xv) != length(yi)) {
        stop(sprintf(
          "Row mismatch in %s -> %s: nrow(x) = %d but length(y) = %d.",
          study_names[i], v, nrow(xv), length(yi)
        ))
      }

      if (i == 1L) {
        p_by_view[v] <- ncol(xv)
        ref_cols[[v]] <- colnames(xv)
      } else {
        if (ncol(xv) != p_by_view[v]) {
          stop(sprintf(
            "Column mismatch in %s -> %s: expected %d variables, got %d.",
            study_names[i], v, p_by_view[v], ncol(xv)
          ))
        }

        ref_has_names <- !is.null(ref_cols[[v]])
        cur_has_names <- !is.null(colnames(xv))
        if (ref_has_names != cur_has_names) {
          stop(sprintf("Colname mismatch in %s -> %s: either all studies must supply colnames or none should.", study_names[i], v))
        }
        if (ref_has_names && !identical(colnames(xv), ref_cols[[v]])) {
          stop(sprintf("Colname mismatch in %s -> %s (feature ordering/content differs from %s).", study_names[i], v, study_names[1]))
        }
      }
    }
  }

  n_by_study <- vapply(y_norm, length, integer(1))
  y_all <- unlist(y_norm, use.names = FALSE)
  x_list_all <- lapply(view_names, function(v) {
    do.call(rbind, lapply(x_norm, function(study) study[[v]]))
  })
  names(x_list_all) <- view_names
  groups_all <- factor(rep(study_names, times = n_by_study), levels = study_names)

  list(
    x = x_norm,
    y = y_norm,
    k = k,
    study_names = study_names,
    view_names = view_names,
    n_views = n_views,
    n_by_study = n_by_study,
    N_all = sum(n_by_study),
    p_by_view = p_by_view,
    p = sum(p_by_view),
    x_list_all = x_list_all,
    y_all = y_all,
    groups_all = groups_all
  )
}

# Helper: normalize and validate study-list test input against training layout.
# Used in `predict.ptLasso.multiview()` and `predict.cv.ptLasso.multiview()`.
ptmv_normalize_newdata <- function(xtest, template) {
  if (!is.list(xtest) || length(xtest) != template$k) {
    stop("xtest must be a list of studies with the same length as the training data.")
  }

  xt_names <- names(xtest)
  if (is.null(xt_names) || any(xt_names == "")) {
    names(xtest) <- template$study_names
  }

  if (!setequal(names(xtest), template$study_names)) {
    stop("xtest study names must match the training study names exactly.")
  }
  xtest <- xtest[template$study_names]

  for (study_name in template$study_names) {
    study_obj <- xtest[[study_name]]
    if (is.matrix(study_obj)) {
      study_obj <- list(View_1 = study_obj)
    }
    if (!is.list(study_obj) || length(study_obj) == 0L) {
      stop(sprintf("xtest[['%s']] must be a matrix or a non-empty list of view matrices.", study_name))
    }
    if (is.null(names(study_obj)) || any(names(study_obj) == "")) {
      names(study_obj) <- paste0("View_", seq_along(study_obj))
    }
    if (!identical(names(study_obj), template$view_names)) {
      stop(sprintf("xtest[['%s']] must use the same view names and ordering as training.", study_name))
    }

    for (view_name in template$view_names) {
      xv <- study_obj[[view_name]]
      if (!is.matrix(xv)) {
        stop(sprintf("xtest[['%s']][['%s']] must be a matrix.", study_name, view_name))
      }
      if (ncol(xv) != template$p_by_view[view_name]) {
        stop(sprintf(
          "xtest[['%s']][['%s']] has %d columns; expected %d.",
          study_name, view_name, ncol(xv), template$p_by_view[view_name]
        ))
      }

      train_cols <- colnames(template$x[[1]][[view_name]])
      test_cols <- colnames(xv)
      if (!is.null(train_cols) || !is.null(test_cols)) {
        if (!identical(train_cols, test_cols)) {
          stop(sprintf("Column names/order mismatch for xtest[['%s']][['%s']].", study_name, view_name))
        }
      }
    }

    xtest[[study_name]] <- study_obj
  }

  xtest
}

# Helper: normalize and validate study-list test outcomes.
# Used in `predict.ptLasso.multiview()` and `predict.cv.ptLasso.multiview()`.
ptmv_prepare_ytest <- function(ytest, study_names) {
  if (is.null(ytest)) {
    return(NULL)
  }
  if (!is.list(ytest) || length(ytest) != length(study_names)) {
    stop("ytest must be a list with one response vector per study.")
  }
  if (is.null(names(ytest)) || any(names(ytest) == "")) {
    names(ytest) <- study_names
  }
  if (!setequal(names(ytest), study_names)) {
    stop("ytest study names must match the training study names exactly.")
  }
  ytest <- ytest[study_names]
  out <- lapply(seq_along(ytest), function(i) {
    yi <- drop(ytest[[i]])
    if (length(yi) == 0L) {
      stop(sprintf("ytest[['%s']] is empty.", study_names[i]))
    }
    yi
  })
  names(out) <- study_names
  out
}

# Helper: expand one study-level baseline per study into one offset per sample.
# Used in `predict.ptLasso.multiview()`.
ptmv_study_offsets <- function(study_sizes, group_baseline) {
  unlist(Map(function(n, b) rep(b, n), study_sizes, as.list(group_baseline)), use.names = FALSE)
}

# Helper: stack study data into one multiview object view by view.
# Used in `predict.ptLasso.multiview()`.
ptmv_stack_by_view <- function(study_list, view_names) {
  out <- lapply(view_names, function(v) {
    do.call(rbind, lapply(study_list, function(study) study[[v]]))
  })
  names(out) <- view_names
  out
}

# Helper: split a stacked prediction vector back into studies.
# Used in `ptLasso.multiview()` and `predict.ptLasso.multiview()`.
ptmv_split_vector_by_study <- function(x, study_sizes, study_names) {
  split(x, rep(study_names, times = study_sizes))[study_names]
}

# Helper: predict one list of study-specific models and return study-wise outputs.
# Used in `predict.ptLasso.multiview()`.
ptmv_predict_by_study <- function(model_list, xtest, s, type, offsets = NULL) {
  out <- vector("list", length(model_list))
  names(out) <- names(model_list)
  for (study_name in names(model_list)) {
    model <- model_list[[study_name]]
    newoffset <- if (is.null(offsets)) NULL else offsets[[study_name]]
    pred <- predict(model, newx = xtest[[study_name]], s = s, type = type, newoffset = newoffset)
    out[[study_name]] <- as.numeric(pred)
  }
  out
}

# Helper: extract nonzero support from one cvar/cv.multiview object at `s`.
# Used in `predict.ptLasso.multiview()`.
ptmv_get_support <- function(model, s) {
  lambda <- ptmv_resolve_s(model, s)
  beta <- as.numeric(coef(model$multiview.fit, s = lambda))
  which(beta[-1] != 0)
}

# Helper: union nonzero support across multiple multiview fits.
# Used in `predict.ptLasso.multiview()`.
ptmv_get_union_support <- function(models, s) {
  sort(unique(unlist(lapply(models, ptmv_get_support, s = s))))
}

# Helper: compute one scalar performance metric for one study or stacked data.
# Used in `cv.ptLasso.multiview()`, `predict.ptLasso.multiview()`,
# and `predict.cv.ptLasso.multiview()`.
ptmv_metric_value <- function(y, pred, family, type.measure) {
  y <- as.numeric(y)
  pred <- as.numeric(pred)

  if (family == "gaussian") {
    if (type.measure %in% c("mse", "deviance")) {
      return(mean((y - pred)^2))
    }
    stop("Unsupported type.measure for gaussian.")
  }

  if (family == "binomial") {
    if (type.measure == "auc") {
      return(as.numeric(pROC::auc(
        pROC::roc(response = y, predictor = pred, levels = c(0, 1), direction = "<", quiet = TRUE)
      )))
    }
    if (type.measure == "deviance") {
      eps <- 1e-8
      pp <- pmin(pmax(pred, eps), 1 - eps)
      return(-2 * mean(y * log(pp) + (1 - y) * log(1 - pp)))
    }
    stop("Unsupported type.measure for binomial.")
  }

  stop("Unsupported family.")
}

# Helper: summarize study-wise predictions into overall and per-study metrics.
# Used in `cv.ptLasso.multiview()`, `predict.ptLasso.multiview()`,
# and `predict.cv.ptLasso.multiview()`.
ptmv_summarize_metric <- function(preds, y, family, type.measure, add_r2 = FALSE) {
  study_names <- names(y)
  study_err <- vapply(study_names, function(study_name) {
    ptmv_metric_value(y[[study_name]], preds[[study_name]], family, type.measure)
  }, numeric(1))

  all_y <- unlist(y, use.names = FALSE)
  all_pred <- unlist(preds, use.names = FALSE)
  out <- c(
    allGroups = ptmv_metric_value(all_y, all_pred, family, type.measure),
    mean = mean(study_err),
    setNames(study_err, paste0("group_", study_names))
  )

  if (add_r2 && family == "gaussian") {
    out <- c(out, "r^2" = 1 - sum((all_y - all_pred)^2) / sum((all_y - mean(all_y))^2))
  }
  out
}

# Helper: choose one prediction per study from a nested alpha-indexed list.
# Used in `predict.cv.ptLasso.multiview()` when combining varying-alpha results.
ptmv_select_predictions_by_alpha <- function(pred_list, study_names, alpha_vec) {
  out <- vector("list", length(study_names))
  names(out) <- study_names
  for (study_name in study_names) {
    out[[study_name]] <- pred_list[[study_name]][[as.character(alpha_vec[[study_name]])]]
  }
  out
}

# Helper: assemble a common prediction object for direct-fit and cv-fit methods.
# Used in `predict.ptLasso.multiview()` and `predict.cv.ptLasso.multiview()`.
ptmv_build_prediction_object <- function(call, alpha, type.measure, yhatoverall, yhatind, yhatpre,
                                         supoverall, supind, suppre.common, suppre.individual,
                                         linkoverall = NULL, linkind = NULL, linkpre = NULL,
                                         erroverall = NULL, errind = NULL, errpre = NULL,
                                         fit = NULL, class_name = "predict.ptLasso.multiview") {
  out <- list(
    call = call,
    alpha = alpha,
    yhatoverall = yhatoverall,
    yhatind = yhatind,
    yhatpre = yhatpre,
    supoverall = supoverall,
    supind = supind,
    suppre.common = suppre.common,
    suppre.individual = suppre.individual,
    type.measure = type.measure
  )
  if (!is.null(linkoverall)) out$linkoverall <- linkoverall
  if (!is.null(linkind)) out$linkind <- linkind
  if (!is.null(linkpre)) out$linkpre <- linkpre
  if (!is.null(erroverall)) out$erroverall <- erroverall
  if (!is.null(errind)) out$errind <- errind
  if (!is.null(errpre)) out$errpre <- errpre
  if (!is.null(fit)) out$fit <- fit
  if (identical(class_name, "predict.cv.ptLasso.multiview")) {
    class(out) <- c("predict.cv.ptLasso.multiview", "predict.ptLasso.multiview")
  } else {
    class(out) <- "predict.ptLasso.multiview"
  }
  out
}
