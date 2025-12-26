# =============================================================================
# STATISTICAL ANALYSIS AND MODELING FOR RENEWABLE ENERGY CLASSIFIER
# =============================================================================
# Author: Sahibjeet Pal Singh
# Date: December 2025
# Description: Advanced statistical analysis, hypothesis testing, and
#              machine learning model evaluation for energy classification
# =============================================================================

# =============================================================================
# PACKAGE DEPENDENCIES
# =============================================================================

#' Load analysis packages
load_analysis_packages <- function() {
  required <- c("stats", "dplyr", "tidyr", "purrr", "broom")
  
  for (pkg in required) {
    if (pkg %in% installed.packages()[, "Package"]) {
      suppressPackageStartupMessages(library(pkg, character.only = TRUE))
    }
  }
}

load_analysis_packages()

# =============================================================================
# DESCRIPTIVE STATISTICS
# =============================================================================

#' Calculate comprehensive descriptive statistics for a numeric vector
#' @param x Numeric vector
#' @param na.rm Remove NA values
#' @return Named list of statistics
descriptive_stats <- function(x, na.rm = TRUE) {
  if (na.rm) {
    x <- x[!is.na(x)]
  }
  
  n <- length(x)
  
  if (n == 0) {
    return(list(
      n = 0,
      mean = NA,
      median = NA,
      sd = NA,
      var = NA,
      min = NA,
      max = NA,
      range = NA,
      q1 = NA,
      q3 = NA,
      iqr = NA,
      skewness = NA,
      kurtosis = NA,
      se = NA,
      cv = NA
    ))
  }
  
  q <- quantile(x, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
  m <- mean(x, na.rm = TRUE)
  s <- sd(x, na.rm = TRUE)
  
  # Calculate skewness (Fisher's method)
  skewness <- if (n > 2 && s > 0) {
    sum((x - m)^3) / ((n - 1) * s^3)
  } else {
    NA
  }
  
  # Calculate excess kurtosis (Fisher's method)
  kurtosis <- if (n > 3 && s > 0) {
    (sum((x - m)^4) / ((n - 1) * s^4)) - 3
  } else {
    NA
  }
  
  list(
    n = n,
    mean = m,
    median = q[2],
    sd = s,
    var = var(x, na.rm = TRUE),
    min = min(x, na.rm = TRUE),
    max = max(x, na.rm = TRUE),
    range = max(x, na.rm = TRUE) - min(x, na.rm = TRUE),
    q1 = q[1],
    q3 = q[3],
    iqr = q[3] - q[1],
    skewness = skewness,
    kurtosis = kurtosis,
    se = s / sqrt(n),
    cv = if (m != 0) (s / abs(m)) * 100 else NA
  )
}

#' Generate summary table for multiple variables
#' @param df Dataframe
#' @param variables Vector of variable names (NULL = all numeric)
#' @return Dataframe with statistics
generate_summary_table <- function(df, variables = NULL) {
  if (is.null(variables)) {
    variables <- names(df)[sapply(df, is.numeric)]
  }
  
  stats_list <- lapply(variables, function(var) {
    stats <- descriptive_stats(df[[var]])
    data.frame(
      variable = var,
      n = stats$n,
      mean = round(stats$mean, 3),
      sd = round(stats$sd, 3),
      min = round(stats$min, 3),
      median = round(stats$median, 3),
      max = round(stats$max, 3),
      skewness = round(stats$skewness, 3),
      kurtosis = round(stats$kurtosis, 3)
    )
  })
  
  do.call(rbind, stats_list)
}

#' Calculate statistics by group
#' @param df Dataframe
#' @param group_var Grouping variable name
#' @param value_var Value variable name
#' @return Grouped statistics dataframe
grouped_stats <- function(df, group_var, value_var) {
  df %>%
    group_by(across(all_of(group_var))) %>%
    summarise(
      n = n(),
      mean = mean(.data[[value_var]], na.rm = TRUE),
      sd = sd(.data[[value_var]], na.rm = TRUE),
      se = sd / sqrt(n),
      median = median(.data[[value_var]], na.rm = TRUE),
      min = min(.data[[value_var]], na.rm = TRUE),
      max = max(.data[[value_var]], na.rm = TRUE),
      .groups = "drop"
    )
}

# =============================================================================
# HYPOTHESIS TESTING
# =============================================================================

#' Perform one-way ANOVA with post-hoc tests
#' @param df Dataframe
#' @param group_var Grouping variable
#' @param value_var Response variable
#' @return List with ANOVA results and post-hoc tests
perform_anova <- function(df, group_var, value_var) {
  # Fit ANOVA model
  formula <- as.formula(paste(value_var, "~", group_var))
  model <- aov(formula, data = df)
  
  # Get summary
  anova_summary <- summary(model)
  
  # Effect size (eta-squared)
  ss_between <- anova_summary[[1]]["Sum Sq"][1, 1]
  ss_total <- sum(anova_summary[[1]]["Sum Sq"])
  eta_squared <- ss_between / ss_total
  
  # Post-hoc test (Tukey HSD)
  tukey_result <- TukeyHSD(model)
  
  # Pairwise t-tests with Bonferroni correction
  pairwise <- pairwise.t.test(
    df[[value_var]],
    df[[group_var]],
    p.adjust.method = "bonferroni"
  )
  
  list(
    model = model,
    summary = anova_summary,
    f_statistic = anova_summary[[1]]["F value"][1, 1],
    p_value = anova_summary[[1]]["Pr(>F)"][1, 1],
    eta_squared = eta_squared,
    tukey = tukey_result,
    pairwise = pairwise
  )
}

#' Perform chi-squared test for categorical variables
#' @param df Dataframe
#' @param var1 First variable
#' @param var2 Second variable
#' @return Chi-squared test results
perform_chi_squared <- function(df, var1, var2) {
  # Create contingency table
  contingency <- table(df[[var1]], df[[var2]])
  
  # Perform test
  test_result <- chisq.test(contingency)
  
  # Calculate Cramer's V
  n <- sum(contingency)
  k <- min(nrow(contingency), ncol(contingency))
  cramers_v <- sqrt(test_result$statistic / (n * (k - 1)))
  
  list(
    contingency_table = contingency,
    chi_squared = test_result$statistic,
    df = test_result$parameter,
    p_value = test_result$p.value,
    cramers_v = as.numeric(cramers_v),
    expected = test_result$expected
  )
}

#' Perform t-test comparing two groups
#' @param df Dataframe
#' @param group_var Grouping variable
#' @param value_var Value variable
#' @param paired Whether to perform paired t-test
#' @return T-test results
perform_t_test <- function(df, group_var, value_var, paired = FALSE) {
  groups <- unique(df[[group_var]])
  
  if (length(groups) != 2) {
    stop("T-test requires exactly 2 groups")
  }
  
  group1 <- df[[value_var]][df[[group_var]] == groups[1]]
  group2 <- df[[value_var]][df[[group_var]] == groups[2]]
  
  # Perform t-test
  test_result <- t.test(group1, group2, paired = paired)
  
  # Calculate Cohen's d
  pooled_sd <- sqrt(((length(group1) - 1) * var(group1) + 
                      (length(group2) - 1) * var(group2)) / 
                     (length(group1) + length(group2) - 2))
  cohens_d <- (mean(group1) - mean(group2)) / pooled_sd
  
  list(
    groups = groups,
    t_statistic = test_result$statistic,
    df = test_result$parameter,
    p_value = test_result$p.value,
    confidence_interval = test_result$conf.int,
    mean_difference = mean(group1) - mean(group2),
    cohens_d = cohens_d
  )
}

# =============================================================================
# CORRELATION ANALYSIS
# =============================================================================

#' Calculate correlation matrix with p-values
#' @param df Dataframe
#' @param variables Variables to include (NULL = all numeric)
#' @param method Correlation method ("pearson", "spearman", "kendall")
#' @return List with correlation matrix and p-values
correlation_analysis <- function(df, variables = NULL, method = "pearson") {
  if (is.null(variables)) {
    numeric_df <- df %>% select(where(is.numeric))
  } else {
    numeric_df <- df %>% select(all_of(variables))
  }
  
  n_vars <- ncol(numeric_df)
  var_names <- names(numeric_df)
  
  # Initialize matrices
  cor_matrix <- matrix(NA, n_vars, n_vars)
  p_matrix <- matrix(NA, n_vars, n_vars)
  rownames(cor_matrix) <- colnames(cor_matrix) <- var_names
  rownames(p_matrix) <- colnames(p_matrix) <- var_names
  
  # Calculate pairwise correlations
  for (i in 1:n_vars) {
    for (j in 1:n_vars) {
      if (i == j) {
        cor_matrix[i, j] <- 1
        p_matrix[i, j] <- 0
      } else {
        test <- cor.test(
          numeric_df[[i]],
          numeric_df[[j]],
          method = method,
          use = "pairwise.complete.obs"
        )
        cor_matrix[i, j] <- test$estimate
        p_matrix[i, j] <- test$p.value
      }
    }
  }
  
  list(
    correlation = cor_matrix,
    p_values = p_matrix,
    method = method,
    n = nrow(numeric_df)
  )
}

#' Find significant correlations
#' @param cor_result Result from correlation_analysis
#' @param threshold P-value threshold
#' @param min_cor Minimum absolute correlation
#' @return Dataframe of significant correlations
find_significant_correlations <- function(cor_result, threshold = 0.05, min_cor = 0.3) {
  cor_mat <- cor_result$correlation
  p_mat <- cor_result$p_values
  
  # Extract pairs
  pairs <- expand.grid(
    var1 = rownames(cor_mat),
    var2 = colnames(cor_mat),
    stringsAsFactors = FALSE
  )
  
  pairs$correlation <- mapply(function(i, j) cor_mat[i, j], pairs$var1, pairs$var2)
  pairs$p_value <- mapply(function(i, j) p_mat[i, j], pairs$var1, pairs$var2)
  
  # Filter
  pairs %>%
    filter(
      var1 < var2,  # Remove duplicates and diagonal
      abs(correlation) >= min_cor,
      p_value < threshold
    ) %>%
    arrange(desc(abs(correlation)))
}

# =============================================================================
# REGRESSION ANALYSIS
# =============================================================================

#' Fit linear regression model
#' @param df Dataframe
#' @param response Response variable name
#' @param predictors Vector of predictor names
#' @return List with model and diagnostics
fit_linear_model <- function(df, response, predictors) {
  # Build formula
  formula <- as.formula(paste(response, "~", paste(predictors, collapse = " + ")))
  
  # Fit model
  model <- lm(formula, data = df)
  
  # Get summary
  model_summary <- summary(model)
  
  # Calculate additional diagnostics
  list(
    model = model,
    formula = formula,
    coefficients = coef(model),
    r_squared = model_summary$r.squared,
    adj_r_squared = model_summary$adj.r.squared,
    f_statistic = model_summary$fstatistic[1],
    p_value = pf(
      model_summary$fstatistic[1],
      model_summary$fstatistic[2],
      model_summary$fstatistic[3],
      lower.tail = FALSE
    ),
    residual_se = model_summary$sigma,
    aic = AIC(model),
    bic = BIC(model),
    vif = if (length(predictors) > 1) {
      tryCatch(car::vif(model), error = function(e) NULL)
    } else NULL
  )
}

#' Fit multiple regression models and compare
#' @param df Dataframe
#' @param response Response variable
#' @param predictor_sets List of predictor vectors
#' @return Comparison table
compare_models <- function(df, response, predictor_sets) {
  results <- lapply(seq_along(predictor_sets), function(i) {
    predictors <- predictor_sets[[i]]
    model_result <- fit_linear_model(df, response, predictors)
    
    data.frame(
      model = paste0("Model_", i),
      predictors = paste(predictors, collapse = " + "),
      r_squared = round(model_result$r_squared, 4),
      adj_r_squared = round(model_result$adj_r_squared, 4),
      aic = round(model_result$aic, 2),
      bic = round(model_result$bic, 2),
      p_value = format(model_result$p_value, scientific = TRUE, digits = 3)
    )
  })
  
  do.call(rbind, results)
}

# =============================================================================
# CLASSIFICATION EVALUATION
# =============================================================================

#' Calculate confusion matrix metrics
#' @param actual Actual values
#' @param predicted Predicted values
#' @param positive Positive class label (for binary)
#' @return List with metrics
confusion_matrix_metrics <- function(actual, predicted, positive = NULL) {
  # Create confusion matrix
  cm <- table(Predicted = predicted, Actual = actual)
  
  # Overall accuracy
  accuracy <- sum(diag(cm)) / sum(cm)
  
  # Per-class metrics
  classes <- unique(c(actual, predicted))
  class_metrics <- lapply(classes, function(cls) {
    tp <- cm[cls, cls]
    fp <- sum(cm[cls, ]) - tp
    fn <- sum(cm[, cls]) - tp
    tn <- sum(cm) - tp - fp - fn
    
    precision <- if ((tp + fp) > 0) tp / (tp + fp) else 0
    recall <- if ((tp + fn) > 0) tp / (tp + fn) else 0
    f1 <- if ((precision + recall) > 0) {
      2 * precision * recall / (precision + recall)
    } else 0
    specificity <- if ((tn + fp) > 0) tn / (tn + fp) else 0
    
    list(
      class = cls,
      tp = tp,
      fp = fp,
      fn = fn,
      tn = tn,
      precision = precision,
      recall = recall,
      f1 = f1,
      specificity = specificity
    )
  })
  
  # Aggregate metrics
  precision_macro <- mean(sapply(class_metrics, function(x) x$precision))
  recall_macro <- mean(sapply(class_metrics, function(x) x$recall))
  f1_macro <- mean(sapply(class_metrics, function(x) x$f1))
  
  # Weighted averages
  class_counts <- table(actual)
  weights <- class_counts / sum(class_counts)
  
  precision_weighted <- sum(sapply(class_metrics, function(x) x$precision) * 
                             weights[sapply(class_metrics, function(x) x$class)])
  recall_weighted <- sum(sapply(class_metrics, function(x) x$recall) * 
                          weights[sapply(class_metrics, function(x) x$class)])
  f1_weighted <- sum(sapply(class_metrics, function(x) x$f1) * 
                      weights[sapply(class_metrics, function(x) x$class)])
  
  # Cohen's Kappa
  expected_accuracy <- sum(rowSums(cm) * colSums(cm)) / sum(cm)^2
  kappa <- (accuracy - expected_accuracy) / (1 - expected_accuracy)
  
  list(
    confusion_matrix = cm,
    accuracy = accuracy,
    kappa = kappa,
    class_metrics = class_metrics,
    macro = list(
      precision = precision_macro,
      recall = recall_macro,
      f1 = f1_macro
    ),
    weighted = list(
      precision = precision_weighted,
      recall = recall_weighted,
      f1 = f1_weighted
    )
  )
}

#' Print formatted classification report
#' @param metrics Result from confusion_matrix_metrics
print_classification_report <- function(metrics) {
  cat("\n")
  cat("═══════════════════════════════════════════════════════════════\n")
  cat("                   CLASSIFICATION REPORT                        \n")
  cat("═══════════════════════════════════════════════════════════════\n\n")
  
  cat("CONFUSION MATRIX:\n")
  print(metrics$confusion_matrix)
  
  cat("\n\nOVERALL METRICS:\n")
  cat("───────────────────────────────────────────────────────────────\n")
  cat(sprintf("  Accuracy:      %.4f (%.2f%%)\n", 
              metrics$accuracy, metrics$accuracy * 100))
  cat(sprintf("  Cohen's Kappa: %.4f\n", metrics$kappa))
  
  cat("\n\nPER-CLASS METRICS:\n")
  cat("───────────────────────────────────────────────────────────────\n")
  cat(sprintf("%-12s %10s %10s %10s %10s\n", 
              "Class", "Precision", "Recall", "F1-Score", "Support"))
  cat("───────────────────────────────────────────────────────────────\n")
  
  for (cm in metrics$class_metrics) {
    support <- cm$tp + cm$fn
    cat(sprintf("%-12s %10.4f %10.4f %10.4f %10d\n",
                cm$class, cm$precision, cm$recall, cm$f1, support))
  }
  
  cat("───────────────────────────────────────────────────────────────\n")
  cat(sprintf("%-12s %10.4f %10.4f %10.4f\n",
              "Macro Avg", metrics$macro$precision, 
              metrics$macro$recall, metrics$macro$f1))
  cat(sprintf("%-12s %10.4f %10.4f %10.4f\n",
              "Weighted Avg", metrics$weighted$precision,
              metrics$weighted$recall, metrics$weighted$f1))
  
  cat("\n═══════════════════════════════════════════════════════════════\n\n")
}

# =============================================================================
# CROSS-VALIDATION
# =============================================================================

#' Perform k-fold cross-validation for classification
#' @param df Dataframe
#' @param classify_fn Classification function
#' @param k Number of folds
#' @param seed Random seed
#' @return Cross-validation results
cross_validate <- function(df, classify_fn, k = 5, seed = 42) {
  set.seed(seed)
  
  n <- nrow(df)
  folds <- sample(rep(1:k, length.out = n))
  
  fold_results <- lapply(1:k, function(fold) {
    # Split data
    train_idx <- folds != fold
    test_idx <- folds == fold
    
    train_data <- df[train_idx, ]
    test_data <- df[test_idx, ]
    
    # Apply classification
    predictions <- classify_fn(train_data, test_data)
    
    # Calculate metrics
    metrics <- confusion_matrix_metrics(
      test_data$best_resource,
      predictions
    )
    
    list(
      fold = fold,
      accuracy = metrics$accuracy,
      f1_macro = metrics$macro$f1,
      kappa = metrics$kappa
    )
  })
  
  # Aggregate results
  accuracies <- sapply(fold_results, function(x) x$accuracy)
  f1_scores <- sapply(fold_results, function(x) x$f1_macro)
  kappas <- sapply(fold_results, function(x) x$kappa)
  
  list(
    k = k,
    fold_results = fold_results,
    accuracy = list(
      mean = mean(accuracies),
      sd = sd(accuracies),
      values = accuracies
    ),
    f1 = list(
      mean = mean(f1_scores),
      sd = sd(f1_scores),
      values = f1_scores
    ),
    kappa = list(
      mean = mean(kappas),
      sd = sd(kappas),
      values = kappas
    )
  )
}

# =============================================================================
# FEATURE IMPORTANCE ANALYSIS
# =============================================================================

#' Analyze feature importance using permutation
#' @param df Dataframe
#' @param features Feature names
#' @param target Target variable
#' @param classify_fn Classification function
#' @param n_permutations Number of permutations
#' @return Feature importance rankings
permutation_importance <- function(df, features, target, classify_fn, 
                                    n_permutations = 10) {
  # Get baseline accuracy
  baseline_predictions <- classify_fn(df, df)
  baseline_accuracy <- mean(baseline_predictions == df[[target]])
  
  # Calculate importance for each feature
  importance_scores <- sapply(features, function(feature) {
    decreases <- sapply(1:n_permutations, function(i) {
      # Create permuted data
      permuted_df <- df
      permuted_df[[feature]] <- sample(df[[feature]])
      
      # Get predictions
      predictions <- classify_fn(df, permuted_df)
      accuracy <- mean(predictions == df[[target]])
      
      # Return decrease in accuracy
      baseline_accuracy - accuracy
    })
    
    mean(decreases)
  })
  
  # Create results dataframe
  results <- data.frame(
    feature = features,
    importance = importance_scores,
    stringsAsFactors = FALSE
  ) %>%
    arrange(desc(importance)) %>%
    mutate(
      rank = row_number(),
      normalized = importance / max(importance)
    )
  
  list(
    baseline_accuracy = baseline_accuracy,
    importance = results
  )
}

#' Calculate information gain for each feature
#' @param df Dataframe
#' @param features Feature names
#' @param target Target variable
#' @return Information gain values
information_gain <- function(df, features, target) {
  # Calculate base entropy
  target_probs <- table(df[[target]]) / nrow(df)
  base_entropy <- -sum(target_probs * log2(target_probs + 1e-10))
  
  # Calculate information gain for each feature
  gains <- sapply(features, function(feature) {
    # Bin numeric features
    if (is.numeric(df[[feature]])) {
      bins <- cut(df[[feature]], breaks = 5, labels = FALSE)
    } else {
      bins <- df[[feature]]
    }
    
    # Calculate conditional entropy
    conditional_entropy <- 0
    for (bin_val in unique(bins)) {
      subset_mask <- bins == bin_val
      if (sum(subset_mask) > 0) {
        subset_probs <- table(df[[target]][subset_mask]) / sum(subset_mask)
        subset_entropy <- -sum(subset_probs * log2(subset_probs + 1e-10))
        weight <- sum(subset_mask) / length(bins)
        conditional_entropy <- conditional_entropy + weight * subset_entropy
      }
    }
    
    base_entropy - conditional_entropy
  })
  
  data.frame(
    feature = features,
    information_gain = gains,
    stringsAsFactors = FALSE
  ) %>%
    arrange(desc(information_gain))
}

# =============================================================================
# OUTLIER DETECTION
# =============================================================================

#' Detect outliers using multiple methods
#' @param x Numeric vector
#' @param method Method: "iqr", "zscore", "mad"
#' @param threshold Threshold value
#' @return Logical vector of outliers
detect_outliers <- function(x, method = "iqr", threshold = 1.5) {
  x <- as.numeric(x)
  
  if (method == "iqr") {
    q1 <- quantile(x, 0.25, na.rm = TRUE)
    q3 <- quantile(x, 0.75, na.rm = TRUE)
    iqr <- q3 - q1
    lower <- q1 - threshold * iqr
    upper <- q3 + threshold * iqr
    outliers <- x < lower | x > upper
    
  } else if (method == "zscore") {
    z_scores <- abs(scale(x))
    outliers <- z_scores > threshold
    
  } else if (method == "mad") {
    median_x <- median(x, na.rm = TRUE)
    mad_x <- mad(x, na.rm = TRUE)
    outliers <- abs(x - median_x) / mad_x > threshold
    
  } else {
    stop("Unknown method: ", method)
  }
  
  outliers[is.na(outliers)] <- FALSE
  return(outliers)
}

#' Remove or flag outliers from dataframe
#' @param df Dataframe
#' @param columns Columns to check for outliers
#' @param method Detection method
#' @param action "flag" or "remove"
#' @return Modified dataframe
handle_outliers <- function(df, columns, method = "iqr", action = "flag") {
  if (action == "flag") {
    df$is_outlier <- FALSE
    for (col in columns) {
      outliers <- detect_outliers(df[[col]], method = method)
      df$is_outlier <- df$is_outlier | outliers
    }
  } else if (action == "remove") {
    for (col in columns) {
      outliers <- detect_outliers(df[[col]], method = method)
      df <- df[!outliers, ]
    }
  }
  
  return(df)
}

# =============================================================================
# ANALYSIS REPORT GENERATION
# =============================================================================

#' Generate comprehensive analysis report
#' @param df Classified dataframe
#' @return Analysis report as list
generate_analysis_report <- function(df) {
  report <- list()
  
  # Basic info
  report$data_summary <- list(
    n_observations = nrow(df),
    n_variables = ncol(df),
    variable_names = names(df)
  )
  
  # Descriptive statistics
  numeric_vars <- names(df)[sapply(df, is.numeric)]
  report$descriptive <- generate_summary_table(df, numeric_vars)
  
  # Classification distribution
  if ("best_resource" %in% names(df)) {
    report$classification_distribution <- df %>%
      count(best_resource) %>%
      mutate(percentage = n / sum(n) * 100)
    
    # Score statistics by resource
    report$scores_by_resource <- df %>%
      group_by(best_resource) %>%
      summarise(
        n = n(),
        solar_mean = mean(solar_score, na.rm = TRUE),
        wind_mean = mean(wind_score, na.rm = TRUE),
        hydro_mean = mean(hydro_score, na.rm = TRUE),
        .groups = "drop"
      )
  }
  
  # Correlation analysis
  score_vars <- intersect(c("solar_score", "wind_score", "hydro_score"), names(df))
  if (length(score_vars) >= 2) {
    report$correlations <- correlation_analysis(df, score_vars)
  }
  
  return(report)
}

#' Print analysis report
#' @param report Report from generate_analysis_report
print_analysis_report <- function(report) {
  cat("\n")
  cat("╔═══════════════════════════════════════════════════════════════╗\n")
  cat("║              COMPREHENSIVE ANALYSIS REPORT                    ║\n")
  cat("╚═══════════════════════════════════════════════════════════════╝\n\n")
  
  cat("📊 DATA SUMMARY\n")
  cat("───────────────────────────────────────────────────────────────\n")
  cat(sprintf("  Observations: %d\n", report$data_summary$n_observations))
  cat(sprintf("  Variables: %d\n", report$data_summary$n_variables))
  
  cat("\n📈 DESCRIPTIVE STATISTICS\n")
  cat("───────────────────────────────────────────────────────────────\n")
  print(report$descriptive, row.names = FALSE)
  
  if (!is.null(report$classification_distribution)) {
    cat("\n🎯 CLASSIFICATION DISTRIBUTION\n")
    cat("───────────────────────────────────────────────────────────────\n")
    print(report$classification_distribution, row.names = FALSE)
  }
  
  if (!is.null(report$scores_by_resource)) {
    cat("\n📊 SCORES BY RESOURCE TYPE\n")
    cat("───────────────────────────────────────────────────────────────\n")
    print(report$scores_by_resource, row.names = FALSE)
  }
  
  cat("\n═══════════════════════════════════════════════════════════════\n\n")
}

# =============================================================================
# INITIALIZATION
# =============================================================================

cat("📈 Analysis Functions v1.0.0 loaded\n")
cat("   Categories: Descriptive, Hypothesis Testing, Regression,\n")
cat("               Classification, Cross-Validation, Feature Importance\n\n")
