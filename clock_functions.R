#' Model Development Functions
#' 
#' Utilities for training and evaluating mortality prediction clocks.

# Required packages
library(tidyverse)
library(glmnet)
library(survival)

# Utility Functions --------------------------------------------------------

#' Process predictions and calculate age gaps
#' @param predictions Data frame of predictions
#' @param data Original data frame
#' @param loess Whether to use LOESS adjustment
#' @return Processed predictions with gaps and z-scores
process_predictions <- function(predictions, data, loess = TRUE) {
  # Join with original data
  pred_processed <- predictions |> 
    left_join(data |> select(twin_id, sex), by = "twin_id")
  
  # Apply age adjustment
  if (loess) {
    # Step 1: Fit LOWESS to get expected predicted age at each chronological age
    age_loess <- loess(.pred ~ age, data = pred_processed, span = 2/3)
    pred_processed <- pred_processed |> 
      mutate(
        loess_predicted_mean = predict(age_loess),
        # Step 2: Calculate gap as difference between predicted and expected age
        adjusted_gap = .pred - loess_predicted_mean
      )
  } else {
    # Use linear regression to model predicted age against chronological age
    age_model <- lm(.pred ~ age, data = pred_processed)
    pred_processed <- pred_processed |> 
      mutate(adjusted_gap = residuals(age_model))
  }
  
  # Calculate z-scores
  gap_mean <- mean(pred_processed$adjusted_gap)
  gap_sd <- sd(pred_processed$adjusted_gap)
  
  pred_processed |> 
    mutate(z_score = (adjusted_gap - gap_mean) / gap_sd) |> 
    select(twin_id, split, age, .pred, adjusted_gap, z_score)
}

#' Create OOF-tuned mortality clocks
#'
#' Generates out-of-fold (OOF) risk scores with pair-aware cross-fitting.
#' Models are trained separately by sex and tuned within each outer-training fold.
#'
#' @param df Data frame containing twin_id, pair_id, sex, age, predictors, time/status
#' @param predictors Named list of predictor sets (one per organ)
#' @param time_col Name of the column containing time-to-event data
#' @param status_col Name of the column containing event status
#' @param age_col Name of the column containing chronological age
#' @param loess Whether to apply LOESS smoothing for age adjustment
#' @param outer_folds Number of outer folds for cross-fitting (default: 5)
#' @param inner_folds Number of folds used by cv.glmnet within training (default: 5)
#' @param alpha_grid Grid of α values to evaluate (default: seq(0, 1, 0.1))
#' @param seed Random seed
#' @return List with `oof_pred_long`, `oof_pred_wide`, `alpha_summary`, and `alpha_by_fold`
create_oof_tuned_mortality_clock <- function(df, predictors,
                                            time_col = "time", status_col = "status",
                                            age_col = "age", loess = TRUE,
                                            outer_folds = 5, inner_folds = 5,
                                            alpha_grid = seq(0, 1, by = 0.1),
                                            seed = 350) {
  make_pair_aware_foldid <- function(pair_id, k, seed) {
    set.seed(seed)
    pairs <- unique(pair_id)
    fold_by_pair <- sample(rep(seq_len(k), length.out = length(pairs)))
    names(fold_by_pair) <- as.character(pairs)
    as.integer(fold_by_pair[as.character(pair_id)])
  }
  
  fit_cv_glmnet_cox <- function(x, y, alpha, foldid, nfolds) {
    tryCatch(
      cv.glmnet(
        x, y,
        family = "cox",
        alpha = alpha,
        nfolds = nfolds,
        foldid = foldid
      ),
      error = function(e) NULL
    )
  }
  
  oof_predict_one_organ_one_sex <- function(data, predictor_set) {
    if (nrow(data) == 0) {
      return(list(pred = tibble(), alpha_by_fold = tibble()))
    }
    
    fold_outer <- make_pair_aware_foldid(data$pair_id, k = outer_folds, seed = seed)
    
    preds <- vector("list", outer_folds)
    alpha_rows <- vector("list", outer_folds)
    
    for (k in seq_len(outer_folds)) {
      train <- data |> filter(fold_outer != k)
      test <- data |> filter(fold_outer == k)
      
      if (nrow(test) == 0 || nrow(train) < 10) {
        next
      }
      
      x_train <- as.matrix(train |> select(all_of(predictor_set)))
      scale_params <- list(
        center = apply(x_train, 2, mean),
        scale = apply(x_train, 2, sd)
      )
      scale_params$scale[scale_params$scale == 0] <- 1
      
      x_train <- scale(x_train, center = scale_params$center, scale = scale_params$scale)
      y_train <- Surv(train[[time_col]], train[[status_col]])
      
      x_test <- as.matrix(test |> select(all_of(predictor_set)))
      x_test <- scale(x_test, center = scale_params$center, scale = scale_params$scale)
      
      foldid_inner <- make_pair_aware_foldid(train$pair_id, k = inner_folds, seed = seed + k)
      
      cv_errors <- map_dbl(alpha_grid, function(a) {
        fit <- fit_cv_glmnet_cox(x_train, y_train, alpha = a, foldid = foldid_inner, nfolds = inner_folds)
        if (is.null(fit)) {
          return(Inf)
        }
        min(fit$cvm)
      })
      
      best_alpha <- alpha_grid[which.min(cv_errors)]
      final_fit <- fit_cv_glmnet_cox(x_train, y_train, alpha = best_alpha, foldid = foldid_inner, nfolds = inner_folds)
      
      if (is.null(final_fit)) {
        next
      }
      
      pred <- as.numeric(predict(final_fit, newx = x_test, s = "lambda.min", type = "link"))
      pred[is.na(pred)] <- mean(pred, na.rm = TRUE)
      
      preds[[k]] <- tibble(
        twin_id = test$twin_id,
        mortality_score = pred,
        split = "oof",
        age = test$age
      )
      
      alpha_rows[[k]] <- tibble(
        fold = k,
        best_alpha = best_alpha,
        cv_error = min(cv_errors, na.rm = TRUE)
      )
    }
    
    list(
      pred = bind_rows(preds),
      alpha_by_fold = bind_rows(alpha_rows)
    )
  }
  
  create_single_oof_clock <- function(predictor_set, organ_name) {
    male_data <- df |> filter(sex == 1) |>
      select(twin_id, pair_id, sex, age = !!sym(age_col), all_of(c(predictor_set, time_col, status_col)))
    female_data <- df |> filter(sex == 2) |>
      select(twin_id, pair_id, sex, age = !!sym(age_col), all_of(c(predictor_set, time_col, status_col)))
    
    male_oof <- oof_predict_one_organ_one_sex(male_data, predictor_set)
    female_oof <- oof_predict_one_organ_one_sex(female_data, predictor_set)
    
    pred <- bind_rows(male_oof$pred, female_oof$pred) |>
      left_join(df |> select(twin_id, sex, age = !!sym(age_col)), by = "twin_id") |>
      mutate(age = coalesce(age.x, age.y)) |>
      select(twin_id, split, age, mortality_score)
    
    pred_processed <- pred |>
      rename(.pred = mortality_score) |>
      process_predictions(df, loess)
    
    alpha_summary <- tibble(
      organ = organ_name,
      male_best_alpha_mean = mean(male_oof$alpha_by_fold$best_alpha, na.rm = TRUE),
      female_best_alpha_mean = mean(female_oof$alpha_by_fold$best_alpha, na.rm = TRUE),
      male_best_alpha_median = median(male_oof$alpha_by_fold$best_alpha, na.rm = TRUE),
      female_best_alpha_median = median(female_oof$alpha_by_fold$best_alpha, na.rm = TRUE)
    )
    
    list(
      organ = organ_name,
      pred_processed = pred_processed |> mutate(organ = organ_name),
      alpha_by_fold_male = male_oof$alpha_by_fold |> mutate(organ = organ_name, sex = "male"),
      alpha_by_fold_female = female_oof$alpha_by_fold |> mutate(organ = organ_name, sex = "female"),
      alpha_summary = alpha_summary
    )
  }
  
  clock_list <- map2(predictors, names(predictors), create_single_oof_clock)
  
  oof_pred_long <- map_dfr(clock_list, "pred_processed") |>
    select(twin_id, organ, z_score, adjusted_gap, .pred, age, split)
  
  oof_pred_wide <- oof_pred_long |>
    select(twin_id, organ, z_score) |>
    pivot_wider(names_from = organ, values_from = z_score)
  
  if ("Ovary" %in% names(oof_pred_wide)) {
    oof_pred_wide <- oof_pred_wide |>
      left_join(df |> select(twin_id, sex), by = "twin_id") |>
      mutate(Ovary = if_else(sex == 1, NA_real_, Ovary)) |>
      select(-sex)
  }
  if ("Prostate" %in% names(oof_pred_wide)) {
    oof_pred_wide <- oof_pred_wide |>
      left_join(df |> select(twin_id, sex), by = "twin_id") |>
      mutate(Prostate = if_else(sex == 2, NA_real_, Prostate)) |>
      select(-sex)
  }
  
  alpha_summary <- bind_rows(map(clock_list, "alpha_summary"))
  alpha_by_fold <- bind_rows(
    map(clock_list, "alpha_by_fold_male"),
    map(clock_list, "alpha_by_fold_female")
  )
  
  list(
    oof_pred_long = oof_pred_long,
    oof_pred_wide = oof_pred_wide,
    alpha_summary = alpha_summary,
    alpha_by_fold = alpha_by_fold
  )
}

#' Create refit tuned mortality clocks (full-data fit)
#'
#' Fits sex-specific models for each organ on all available data (within sex),
#' using provided alpha values and cv.glmnet lambda selection.
#'
#' @param df Data frame containing twin_id, pair_id, sex, age, predictors, time/status
#' @param predictors Named list of predictor sets (one per organ)
#' @param alpha_results List with per-organ `male_best_alpha` and `female_best_alpha`
#' @param time_col Name of the column containing time-to-event data
#' @param status_col Name of the column containing event status
#' @param age_col Name of the column containing chronological age
#' @param loess Whether to apply LOESS smoothing for age adjustment
#' @param nfolds Number of folds used by cv.glmnet (default: 5)
#' @param seed Random seed
#' @return List with `models`, `pred_long`, and `pred_wide`
create_refit_tuned_mortality_clock <- function(df, predictors, alpha_results,
                                              time_col = "time", status_col = "status",
                                              age_col = "age", loess = TRUE,
                                              nfolds = 5, seed = 350) {
  make_pair_aware_foldid <- function(pair_id, k, seed) {
    set.seed(seed)
    pairs <- unique(pair_id)
    fold_by_pair <- sample(rep(seq_len(k), length.out = length(pairs)))
    names(fold_by_pair) <- as.character(pairs)
    as.integer(fold_by_pair[as.character(pair_id)])
  }
  
  fit_one_sex <- function(data, predictor_set, alpha) {
    x <- as.matrix(data |> select(all_of(predictor_set)))
    scale_params <- list(
      center = apply(x, 2, mean),
      scale = apply(x, 2, sd)
    )
    scale_params$scale[scale_params$scale == 0] <- 1
    x <- scale(x, center = scale_params$center, scale = scale_params$scale)
    
    y <- Surv(data[[time_col]], data[[status_col]])
    foldid <- make_pair_aware_foldid(data$pair_id, k = nfolds, seed = seed)
    
    fit <- cv.glmnet(
      x, y,
      family = "cox",
      alpha = alpha,
      nfolds = nfolds,
      foldid = foldid
    )
    
    pred <- as.numeric(predict(fit, newx = x, s = "lambda.min", type = "link"))
    pred[is.na(pred)] <- mean(pred, na.rm = TRUE)
    
    list(
      model = fit,
      scale_params = scale_params,
      pred = tibble(twin_id = data$twin_id, mortality_score = pred, split = "refit", age = data$age)
    )
  }
  
  create_single_refit <- function(predictor_set, organ_name) {
    alphas <- alpha_results[[organ_name]]
    
    male_data <- df |> filter(sex == 1) |>
      select(twin_id, pair_id, sex, age = !!sym(age_col), all_of(c(predictor_set, time_col, status_col)))
    female_data <- df |> filter(sex == 2) |>
      select(twin_id, pair_id, sex, age = !!sym(age_col), all_of(c(predictor_set, time_col, status_col)))
    
    male_fit <- fit_one_sex(male_data, predictor_set, alpha = alphas$male_best_alpha)
    female_fit <- fit_one_sex(female_data, predictor_set, alpha = alphas$female_best_alpha)
    
    pred <- bind_rows(male_fit$pred, female_fit$pred) |>
      left_join(df |> select(twin_id, sex, age = !!sym(age_col)), by = "twin_id") |>
      mutate(age = coalesce(age.x, age.y)) |>
      select(twin_id, split, age, mortality_score)
    
    pred_processed <- pred |>
      rename(.pred = mortality_score) |>
      process_predictions(df, loess) |>
      mutate(organ = organ_name)
    
    list(
      pred_processed = pred_processed,
      male_model = male_fit$model,
      female_model = female_fit$model,
      male_scale_params = male_fit$scale_params,
      female_scale_params = female_fit$scale_params,
      male_alpha = alphas$male_best_alpha,
      female_alpha = alphas$female_best_alpha
    )
  }
  
  models <- map2(predictors, names(predictors), create_single_refit)
  
  pred_long <- map_dfr(models, "pred_processed") |>
    select(twin_id, organ, z_score, adjusted_gap, .pred, age, split)
  
  pred_wide <- pred_long |>
    select(twin_id, organ, z_score) |>
    pivot_wider(names_from = organ, values_from = z_score)
  
  if ("Ovary" %in% names(pred_wide)) {
    pred_wide <- pred_wide |>
      left_join(df |> select(twin_id, sex), by = "twin_id") |>
      mutate(Ovary = if_else(sex == 1, NA_real_, Ovary)) |>
      select(-sex)
  }
  if ("Prostate" %in% names(pred_wide)) {
    pred_wide <- pred_wide |>
      left_join(df |> select(twin_id, sex), by = "twin_id") |>
      mutate(Prostate = if_else(sex == 2, NA_real_, Prostate)) |>
      select(-sex)
  }
  
  list(
    models = models,
    pred_long = pred_long,
    pred_wide = pred_wide
  )
}

#' Extract Scaling Information from OOF Predictions (Long Format)
#'
#' @param oof_pred_long Data frame with columns organ, adjusted_gap
#' @param model_type Either "age" or "mortality"
#' @return A data frame with per-clock scaling summaries from OOF predictions
extract_scaling_info_from_oof <- function(oof_pred_long, model_type = "mortality") {
  mt <- model_type
  oof_pred_long |>
    group_by(clock = organ) |>
    summarise(
      sd_in_original_units = sd(adjusted_gap, na.rm = TRUE),
      mean_gap = mean(adjusted_gap, na.rm = TRUE),
      median_gap = median(adjusted_gap, na.rm = TRUE),
      min_gap = min(adjusted_gap, na.rm = TRUE),
      max_gap = max(adjusted_gap, na.rm = TRUE),
      q25_gap = quantile(adjusted_gap, 0.25, na.rm = TRUE),
      q75_gap = quantile(adjusted_gap, 0.75, na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(
      model_type = mt,
      interpretation = if_else(
        rep(identical(mt, "age"), length(sd_in_original_units)),
        sprintf("1 SD = %.2f years difference from chronological age", sd_in_original_units),
        sprintf("1 SD = %.2f mortality units", sd_in_original_units)
      )
    ) |>
    select(clock, model_type, sd_in_original_units, mean_gap, median_gap, min_gap, max_gap, q25_gap, q75_gap, interpretation)
}
