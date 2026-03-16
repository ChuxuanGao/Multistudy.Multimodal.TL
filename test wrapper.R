library(MimeSys)
library(ptLasso)
library(tidyverse)
library(multiview)

source("gen_simmba_multistudy.R")
source("cvar_multiview.R")
source("ptlasso.multiview.R")

set.seed(1234)

# Example 1: multi-study, multiview Gaussian data for the main wrapper.
test_simulated_dt <- gen_simmba_multistudy(
  nsample = 100,
  nstudy = 3,
  snr = 5,
  rho.beta = 0.5,
  sigma.alpha = 0.1,
  tau.snr = 0.1,
  outcome.type = "continuous",
  nrep = 1,
  seed = 5678
)

build_feature_lists <- function(test_simulated_dt, rep_id = "Rep_1") {
  
  studies <- names(test_simulated_dt$trainDat[[rep_id]])
  
  feature_list_by_study <- list()
  y_by_study <- list()
  
  ## store featureType levels globally (assumed consistent across studies)
  featureType_levels <- NULL
  
  for (s in studies) {
    
    feature_metadata <- test_simulated_dt$trainDat[[rep_id]][[s]]$feature_metadata
    feature_table    <- test_simulated_dt$trainDat[[rep_id]][[s]]$feature_table
    
    # ensure factor
    feature_metadata$featureType <- factor(feature_metadata$featureType)
    
    if (is.null(featureType_levels)) {
      featureType_levels <- levels(feature_metadata$featureType)
    }
    
    # split feature IDs by type
    feature_ids_by_type <- split(
      feature_metadata$featureID,
      feature_metadata$featureType
    )
    
    # slice + transpose
    feature_list <- lapply(feature_ids_by_type, function(ids) {
      ids <- intersect(ids, rownames(feature_table))
      t(feature_table[ids, , drop = FALSE])
    })
    
    # enforce consistent ordering
    feature_list <- feature_list[featureType_levels]
    
    feature_list_by_study[[s]] <- feature_list
    
    # outcome
    y_by_study[[s]] <-
      test_simulated_dt$trainDat[[rep_id]][[s]]$sample_metadata$Y
  }
  
  ## -------------------------
  ## combine all studies
  ## -------------------------
  
  feature_list_all <- lapply(featureType_levels, function(ft) {
    do.call(rbind, lapply(feature_list_by_study, function(x) x[[ft]]))
  })
  names(feature_list_all) <- featureType_levels
  
  y_all <- unlist(y_by_study, use.names = FALSE)
  
  return(list(
    feature_list_by_study = feature_list_by_study,
    y_by_study            = y_by_study,
    feature_list_all      = feature_list_all,
    y_all                 = y_all
  ))
}

build_feature_lists_split <- function(sim_obj, rep_id = "Rep_1", dataset = c("trainDat", "testDat")) {
  dataset <- match.arg(dataset)

  studies <- names(sim_obj[[dataset]][[rep_id]])
  feature_list_by_study <- list()
  y_by_study <- list()
  featureType_levels <- NULL

  for (s in studies) {
    feature_metadata <- sim_obj[[dataset]][[rep_id]][[s]]$feature_metadata
    feature_table <- sim_obj[[dataset]][[rep_id]][[s]]$feature_table

    feature_metadata$featureType <- factor(feature_metadata$featureType)
    if (is.null(featureType_levels)) {
      featureType_levels <- levels(feature_metadata$featureType)
    }

    feature_ids_by_type <- split(
      feature_metadata$featureID,
      feature_metadata$featureType
    )

    feature_list <- lapply(feature_ids_by_type, function(ids) {
      ids <- intersect(ids, rownames(feature_table))
      t(feature_table[ids, , drop = FALSE])
    })
    feature_list <- feature_list[featureType_levels]

    feature_list_by_study[[s]] <- feature_list
    y_by_study[[s]] <- sim_obj[[dataset]][[rep_id]][[s]]$sample_metadata$Y
  }

  feature_list_all <- lapply(featureType_levels, function(ft) {
    do.call(rbind, lapply(feature_list_by_study, function(x) x[[ft]]))
  })
  names(feature_list_all) <- featureType_levels

  y_all <- unlist(y_by_study, use.names = FALSE)

  list(
    feature_list_by_study = feature_list_by_study,
    y_by_study = y_by_study,
    feature_list_all = feature_list_all,
    y_all = y_all
  )
}

train_study <- build_feature_lists_split(test_simulated_dt, "Rep_1", "trainDat")
test_study <- build_feature_lists_split(test_simulated_dt, "Rep_1", "testDat")

x <- train_study$feature_list_by_study
y <- train_study$y_by_study

# Main Gaussian example used for wrapper development checks.
# This is the primary smoke test for ptLasso.multiview().

tt <- ptLasso.multiview(
  x, y,
  alpha = 0.5,
  rho = c(0, 0.5, 1),
  nlambda = 20,
  family = "gaussian",
  type.measure =  "mse",
  overall.lambda = "lambda.1se",
  ind.lambda = "lambda.1se", # NEW
  pre.lambda = "lambda.1se", # NEW
  foldid = NULL, # added
  nfolds = 3,
  verbose = TRUE,
  fitoverall = NULL,
  fitind = NULL,
  penalty.factor = NULL, # added
  group.intercepts = TRUE, # added
  en.alpha = 1# added
)

# Cross-validation wrapper for alpha selection.
cv_tt <- cv.ptLasso.multiview(
  x, y,
  alphalist = c(0, 0.5, 1),
  rho = c(0, 0.5, 1),
  nlambda = 20,
  family = "gaussian",
  type.measure = "mse",
  nfolds = 3,
  verbose = TRUE,
  group.intercepts = TRUE,
  en.alpha = 1
)

# Predict from direct fit and cv fit on held-out test data.
pred_tt <- predict(
  tt,
  xtest = test_study$feature_list_by_study,
  ytest = test_study$y_by_study,
  s = "lambda.min"
)

pred_cv_tt <- predict(
  cv_tt,
  xtest = test_study$feature_list_by_study,
  ytest = test_study$y_by_study,
  s = "lambda.min"
)

pred_cv_tt_varying <- predict(
  cv_tt,
  xtest = test_study$feature_list_by_study,
  ytest = test_study$y_by_study,
  s = "lambda.min",
  alphatype = "varying"
)

# Quick structure checks for the fitted Gaussian objects.
str(tt, max.level = 1)
str(cv_tt, max.level = 1)
str(pred_tt, max.level = 1)
str(pred_cv_tt, max.level = 1)
str(pred_cv_tt_varying, max.level = 1)
names(tt)
names(cv_tt)
names(pred_tt)
names(pred_cv_tt)

# Example 2: optional single-study multiview fallback.
# This should message and return a cvar.multiview object.
# tt_single_study <- ptLasso.multiview(
#   x = x[1],
#   y = y[1],
#   family = "gaussian",
#   type.measure = "mse",
#   rho = c(0, 0.5, 1),
#   nlambda = 20
# )

# Example 3: optional single-view multi-study fallback.
# This should message and return a ptLasso object.
# x_single_view <- lapply(x, function(study) study[1])
# tt_single_view <- ptLasso.multiview(
#   x = x_single_view,
#   y = y,
#   family = "gaussian",
#   type.measure = "mse",
#   nfolds = 3
# )

# Optional binary example kept here for future testing, but not run by default
# because it is much slower than the Gaussian smoke test.
# test_simulated_dt_bi <- gen_simmba_multistudy(
#   nsample = 100,
#   nstudy = 3,
#   snr = 5,
#   rho.beta = 0.5,
#   sigma.alpha = 0.1,
#   tau.snr = 0.1,
#   outcome.type = "binary",
#   nrep = 1,
#   seed = 91011
# )
# test_study_bi <- build_feature_lists(test_simulated_dt_bi, "Rep_1")
# x_bi <- test_study_bi$feature_list_by_study
# y_bi <- test_study_bi$y_by_study
# tt2 <- ptLasso.multiview(
#   x_bi, y_bi,
#   alpha = 0.5,
#   rho = c(0, 0.5, 1),
#   nlambda = 20,
#   family = "binomial",
#   type.measure = "auc",
#   overall.lambda = "lambda.1se",
#   ind.lambda = "lambda.min",
#   pre.lambda = "lambda.min",
#   nfolds = 5,
#   verbose = TRUE,
#   group.intercepts = TRUE,
#   en.alpha = 1
# )

# Error behavior checks.
bad_xtest_names <- test_study$feature_list_by_study
names(bad_xtest_names)[1] <- "WrongStudy"
err_bad_study <- tryCatch(
  predict(tt, xtest = bad_xtest_names, ytest = test_study$y_by_study),
  error = function(e) e$message
)

bad_xtest_view <- test_study$feature_list_by_study
names(bad_xtest_view[[1]]) <- rev(names(bad_xtest_view[[1]]))
err_bad_view <- tryCatch(
  predict(tt, xtest = bad_xtest_view, ytest = test_study$y_by_study),
  error = function(e) e$message
)

bad_xtest_dim <- test_study$feature_list_by_study
bad_xtest_dim[[1]][[1]] <- bad_xtest_dim[[1]][[1]][, -1, drop = FALSE]
err_bad_dim <- tryCatch(
  predict(tt, xtest = bad_xtest_dim, ytest = test_study$y_by_study),
  error = function(e) e$message
)

bad_ytest <- test_study$y_by_study[-1]
err_bad_ytest <- tryCatch(
  predict(tt, xtest = test_study$feature_list_by_study, ytest = bad_ytest),
  error = function(e) e$message
)

err_bad_alpha <- tryCatch(
  predict(cv_tt, xtest = test_study$feature_list_by_study, alpha = 0.33),
  error = function(e) e$message
)

err_bad_alphatype <- tryCatch(
  predict(cv_tt, xtest = test_study$feature_list_by_study, alphatype = "wrong"),
  error = function(e) e$message
)

cat("err_bad_study:", err_bad_study, "\n")
cat("err_bad_view:", err_bad_view, "\n")
cat("err_bad_dim:", err_bad_dim, "\n")
cat("err_bad_ytest:", err_bad_ytest, "\n")
cat("err_bad_alpha:", err_bad_alpha, "\n")
cat("err_bad_alphatype:", err_bad_alphatype, "\n")

