.HP_estimate_ps_kmeans_full <- function(
    data, y_col, covariate_cols, id_col, time_col, unit_cluster,
    ps_nstart, ps_max_iter, ps_tol, ps_seed, ps_standardize,
    ps_G_rule, ps_G_grid, ps_cv_folds, ps_cv_repeats,
    ps_select_nstart, ps_min_group_size) {
  ps_G_rule <- match.arg(ps_G_rule, c("one-se", "minimum", "n-third"))

  if (length(y_col) != 1L || !is.character(y_col) || !nzchar(y_col)) {
    stop("PS k-means requires y_col to name exactly one outcome column.")
  }
  required <- unique(c(id_col, time_col, y_col, covariate_cols))
  missing <- setdiff(required, names(data))
  if (length(missing)) {
    stop("PS k-means data are missing columns: ",
         paste(missing, collapse = ", "), ".")
  }
  if (length(covariate_cols) < 2L) {
    stop("PS k-means requires the treatment followed by at least one control in covariate_cols.")
  }
  if (anyDuplicated(data[c(id_col, time_col)])) {
    stop("PS k-means requires exactly one row for each (unit, time) pair.")
  }

  ids <- sort(unique(data[[id_col]]))
  times <- sort(unique(data[[time_col]]))
  N <- length(ids)
  TT <- length(times)
  K <- length(covariate_cols)
  if (N < 4L) stop("Full-sample PS k-means requires at least four units.")
  if (TT <= 1L) stop("Full-sample PS k-means is implemented only for T > 1.")
  if (nrow(data) != N * TT) stop("PS k-means requires a complete balanced panel.")

  if (!is.null(ps_seed) &&
      (length(ps_seed) != 1L || !is.finite(ps_seed) || ps_seed < 0 ||
       ps_seed > .Machine$integer.max || ps_seed != floor(ps_seed))) {
    stop("ps_seed must be NULL or one nonnegative integer.")
  }
  for (parameter in c("ps_nstart", "ps_max_iter", "ps_cv_folds",
                      "ps_cv_repeats", "ps_select_nstart",
                      "ps_min_group_size")) {
    value <- get(parameter)
    if (length(value) != 1L || !is.finite(value) ||
        value < 1 || value != floor(value)) {
      stop(parameter, " must be one positive integer.")
    }
  }
  if (length(ps_tol) != 1L || !is.finite(ps_tol) || ps_tol <= 0) {
    stop("ps_tol must be one positive finite number.")
  }
  if (length(ps_standardize) != 1L || is.na(ps_standardize) ||
      !is.logical(ps_standardize)) {
    stop("ps_standardize must be TRUE or FALSE.")
  }
  ps_nstart <- as.integer(ps_nstart)
  ps_max_iter <- as.integer(ps_max_iter)
  ps_cv_folds <- as.integer(ps_cv_folds)
  ps_cv_repeats <- as.integer(ps_cv_repeats)
  ps_select_nstart <- as.integer(ps_select_nstart)
  ps_min_group_size <- as.integer(ps_min_group_size)
  if (ps_cv_folds < 2L) stop("ps_cv_folds must be at least 2.")
  if (ps_cv_folds > N) stop("ps_cv_folds cannot exceed the number of units.")
  if (!is.null(ps_G_grid)) {
    if (!is.numeric(ps_G_grid) || !length(ps_G_grid) ||
        any(!is.finite(ps_G_grid)) || any(ps_G_grid < 1) ||
        any(ps_G_grid != floor(ps_G_grid))) {
      stop("ps_G_grid must be NULL or a nonempty vector of positive integers.")
    }
    ps_G_grid <- sort(unique(as.integer(ps_G_grid)))
  }

  ordered <- data[order(match(data[[id_col]], ids),
                        match(data[[time_col]], times)), , drop = FALSE]
  profile_object <- make_ps_profile(
    data = ordered, id = id_col, time = time_col,
    clustering_vars = c(y_col, covariate_cols),
    auxiliary_times = times, standardize = ps_standardize
  )
  profile <- profile_object$matrix

  seed_existed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (seed_existed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (seed_existed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  if (!is.null(ps_seed)) set.seed(as.integer(ps_seed))

  fit_seed <- function(repeat_id = 0L, fold_id = 0L, grid_id = 0L) {
    if (is.null(ps_seed)) return(sample.int(.Machine$integer.max, 1L))
    as.integer((as.double(ps_seed) + 104729 * repeat_id +
                  1009 * fold_id + 37 * grid_id) %% .Machine$integer.max)
  }
  fit_ps <- function(A, G, nstart, seed) {
    ps_kmeans(A, G = G, nstart = nstart, max_iter = ps_max_iter,
              tol = ps_tol, seed = seed)
  }

  candidate_grid <- function(max_G) {
    if (max_G < 1L) stop("No feasible PS group count remains.")
    if (!is.null(ps_G_grid)) {
      if (any(ps_G_grid > max_G)) {
        stop("ps_G_grid exceeds the feasible maximum ", max_G,
             "; reduce the grid, ps_cv_folds, or ps_min_group_size.")
      }
      return(ps_G_grid)
    }
    if (max_G <= 12L) return(seq_len(max_G))
    tail_grid <- unique(as.integer(round(exp(seq(log(9), log(max_G),
                                                  length.out = 12L)))))
    sort(unique(c(seq_len(8L), tail_grid, max_G)))
  }

  select_G <- function(A) {
    n <- nrow(A)
    designs <- list()
    design_id <- 0L
    for (repeat_id in seq_len(ps_cv_repeats)) {
      permutation <- sample.int(n)
      fold_label <- integer(n)
      fold_label[permutation] <- rep(seq_len(ps_cv_folds), length.out = n)
      for (fold in seq_len(ps_cv_folds)) {
        design_id <- design_id + 1L
        designs[[design_id]] <- list(
          repeat_id = repeat_id, fold_id = fold,
          validation = which(fold_label == fold),
          fitting = which(fold_label != fold)
        )
      }
    }
    max_G <- min(
      vapply(designs, function(x) length(x$fitting), numeric(1L)),
      vapply(designs, function(x) length(x$validation), numeric(1L)),
      vapply(designs, function(x) floor(length(x$fitting) /
                                          ps_min_group_size), numeric(1L)),
      vapply(designs, function(x) floor(length(x$validation) /
                                          ps_min_group_size), numeric(1L)),
      floor(n / ps_min_group_size)
    )
    grid <- candidate_grid(max_G)
    if (ps_G_rule == "n-third") {
      target <- max(1L, ceiling(n^(1 / 3)))
      selected <- grid[which.min(abs(grid - target))]
      return(list(
        selected_G = as.integer(selected), rule = "n-third",
        path = data.frame(G = grid, cv_loss = NA_real_, cv_se = NA_real_),
        threshold = NA_real_, minimizer_G = NA_integer_
      ))
    }

    loss_sum <- matrix(0, nrow = n, ncol = length(grid))
    loss_count <- matrix(0L, nrow = n, ncol = length(grid))
    for (design in designs) {
      for (grid_id in seq_along(grid)) {
        G <- grid[grid_id]
        fitted <- fit_ps(
          A[design$fitting, , drop = FALSE], G, ps_select_nstart,
          fit_seed(design$repeat_id, design$fold_id, grid_id)
        )
        distances <- .ps_squared_distances(
          A[design$validation, , drop = FALSE], fitted$centers
        )
        group <- .ps_balanced_assignment(
          distances, .ps_balanced_capacities(length(design$validation), G)
        )
        losses <- distances[cbind(seq_along(group), group)]
        loss_sum[design$validation, grid_id] <-
          loss_sum[design$validation, grid_id] + losses
        loss_count[design$validation, grid_id] <-
          loss_count[design$validation, grid_id] + 1L
      }
    }
    if (any(loss_count == 0L)) stop("Internal PS validation loss was not filled.")
    unit_loss <- loss_sum / loss_count
    cv_loss <- colMeans(unit_loss)
    cv_se <- apply(unit_loss, 2L, stats::sd) / sqrt(n)
    cv_se[!is.finite(cv_se)] <- 0
    minimizer_index <- which.min(cv_loss)
    threshold <- cv_loss[minimizer_index] + cv_se[minimizer_index]
    selected_index <- if (ps_G_rule == "minimum") {
      minimizer_index
    } else {
      which(cv_loss <= threshold)[1L]
    }
    list(
      selected_G = as.integer(grid[selected_index]), rule = ps_G_rule,
      path = data.frame(G = grid, cv_loss = cv_loss, cv_se = cv_se),
      threshold = threshold,
      minimizer_G = as.integer(grid[minimizer_index])
    )
  }

  max_final_G <- floor(N / ps_min_group_size)
  if (is.null(unit_cluster)) {
    selection <- select_G(profile)
    G_unit <- selection$selected_G
  } else {
    if (length(unit_cluster) != 1L || !is.finite(unit_cluster) ||
        unit_cluster < 1 || unit_cluster != floor(unit_cluster)) {
      stop("For full-sample PS k-means, unit_cluster must be one positive integer.")
    }
    G_unit <- as.integer(unit_cluster)
    selection <- list(rule = "fixed", path = NULL,
                      threshold = NA_real_, minimizer_G = NA_integer_)
  }
  if (G_unit > max_final_G) {
    stop("For full-sample PS k-means, unit_cluster must not exceed ",
         max_final_G, ", as implied by N and ps_min_group_size.")
  }
  ps_fit <- fit_ps(profile, G_unit, ps_nstart, fit_seed(repeat_id = 999L))
  unit_group <- as.integer(ps_fit$group)

  y <- matrix(NA_real_, N, TT)
  X <- array(NA_real_, c(N, TT, K))
  id_index <- match(data[[id_col]], ids)
  time_index <- match(data[[time_col]], times)
  for (row in seq_len(nrow(data))) {
    y[id_index[row], time_index[row]] <- data[[y_col]][row]
    for (k in seq_len(K)) X[id_index[row], time_index[row], k] <-
      data[[covariate_cols[k]]][row]
  }
  if (any(!is.finite(y)) || any(!is.finite(X))) {
    stop("PS k-means requires finite outcome and covariate values.")
  }

  Du <- sapply(seq_len(G_unit), function(g) as.numeric(unit_group == g))
  if (G_unit == 1L) Du <- matrix(Du, ncol = 1L)
  Mu <- diag(N) - Du %*% diag(1 / colSums(Du), G_unit) %*% t(Du)
  Z <- array(NA_real_, c(N, TT, K + 1L))
  Z[, , 1L] <- y
  for (k in seq_len(K)) Z[, , k + 1L] <- X[, , k]
  Z_projected <- array(NA_real_, dim(Z))
  for (k in seq_len(K + 1L)) Z_projected[, , k] <- Mu %*% Z[, , k]
  tY <- as.vector(Z_projected[, , 1L])
  tX <- do.call(cbind, lapply(seq_len(K), function(k) {
    as.vector(Z_projected[, , k + 1L])
  }))
  colnames(tX) <- covariate_cols
  treatment <- tX[, 1L]
  controls <- tX[, -1L, drop = FALSE]

  fit <- hdm::rlassoEffect(x = controls, y = tY, d = treatment,
                           method = "double selection")
  trans <- data.frame(y = tY, D = treatment, controls, check.names = FALSE)
  Ytilde <- hdm::rlasso(y ~ . - D - 1, data = trans)$residuals
  Dtilde <- hdm::rlasso(D ~ . - 1,
                        data = trans[, -1L, drop = FALSE])$residuals
  data_res <- data.frame(
    id = rep(ids, TT), time = rep(times, each = N),
    Ytilde = Ytilde, Dtilde = Dtilde
  )
  post_plm <- plm::plm(Ytilde ~ -1 + Dtilde, data = data_res,
                       model = "pooling", index = c("id", "time"))
  robust_se <- sqrt(diag(plm::vcovHC(post_plm, type = "HC0",
                                     method = "arellano")))
  coefs <- stats::coef(post_plm)
  df_ps <- max(1, N * TT - TT * G_unit)
  se_corrected <- robust_se * sqrt(N * TT / df_ps)
  summary_table <- data.frame(
    Estimate = coefs,
    `Std. Error corrected` = se_corrected,
    `t-value corrected` = coefs / se_corrected,
    `Pr(>|t|) corrected` = 2 * stats::pt(-abs(coefs / se_corrected), df_ps),
    check.names = FALSE
  )

  list(
    cluster_method = "ps-kmeans", cross_fitted = FALSE, cross_type = "none",
    G_unit = G_unit, unit_group = unit_group,
    fit_summary = summary(fit), post_plm_summary = summary(post_plm),
    estimate_corrected = summary_table, summary_table = summary_table,
    ps_profile_dimension = profile_object$M_d,
    ps_objective = ps_fit$objective,
    ps_group_sizes = ps_fit$group_sizes,
    ps_converged = ps_fit$converged,
    ps_iterations = ps_fit$iterations,
    ps_selection_rule = selection$rule,
    ps_selection_path = selection$path,
    ps_selection_threshold = selection$threshold,
    ps_selection_minimizer_G = selection$minimizer_G,
    ps_profile_center = profile_object$center,
    ps_profile_scale = profile_object$scale
  )
}