# Quadratic-projection K-means primitives and full-sample estimator.
#
# Unit i is represented by row i of B = Z Z' / M.  Because the shared
# balanced K-means code divides squared row distances by the number of
# columns, its loss is
#
#   N^{-1} sum_i ||B_i - center_{g_i}||_2^2 / N,
#
# which equals K-means on Phi = B / sqrt(N) and is measured in the same
# units as the quadratic aggregated-projection distance.

.pk_validate_options <- function(
    nstart, max_iter, select_nstart, select_max_iter,
    tol, seed, standardize, G_rule, G_grid,
    min_group_size, noise_multiplier, max_G = NULL,
    log_N = "log") {
  for (item in list(nstart = nstart, max_iter = max_iter,
                    select_nstart = select_nstart,
                    select_max_iter = select_max_iter,
                    min_group_size = min_group_size)) {
    if (length(item) != 1L || !is.finite(item) || item < 1L ||
        item != floor(item)) {
      stop("Projection K-means integer options must be positive integers.")
    }
  }
  if (length(tol) != 1L || !is.finite(tol) || tol <= 0) {
    stop("pk_tol must be one positive finite number.")
  }
  if (!is.null(seed) &&
      (length(seed) != 1L || !is.finite(seed) || seed < 0 ||
       seed > .Machine$integer.max || seed != floor(seed))) {
    stop("pk_seed must be NULL or one nonnegative integer.")
  }
  if (length(standardize) != 1L || is.na(standardize) ||
      !is.logical(standardize)) {
    stop("pk_standardize must be TRUE or FALSE.")
  }
  if (length(noise_multiplier) != 1L || !is.finite(noise_multiplier) ||
      noise_multiplier <= 0) {
    stop("pk_noise_multiplier must be one positive finite number.")
  }
  G_rule <- match.arg(
    G_rule[1L], c("noise-floor", "inference", "elbow", "n-third")
  )
  if (!is.null(G_grid)) {
    if (!is.numeric(G_grid) || !length(G_grid) || any(!is.finite(G_grid)) ||
        any(G_grid < 1L) || any(G_grid != floor(G_grid))) {
      stop("pk_G_grid must be NULL or a nonempty vector of positive integers.")
    }
    G_grid <- sort(unique(as.integer(G_grid)))
  }
  if (!is.null(max_G)) {
    if (!is.numeric(max_G) || length(max_G) != 1L || !is.finite(max_G) ||
        max_G < 1L || max_G != floor(max_G)) {
      stop("pk_max_G must be NULL or one positive integer.")
    }
    max_G <- as.integer(max_G)
  }
  if (is.character(log_N)) {
    if (length(log_N) != 1L || is.na(log_N) ||
        !log_N %in% c("log", "loglog")) {
      stop("pk_log_N must be 'log', 'loglog', or one positive finite number.")
    }
  } else if (!is.numeric(log_N) || length(log_N) != 1L ||
             !is.finite(log_N) || log_N <= 0) {
    stop("pk_log_N must be 'log', 'loglog', or one positive finite number.")
  }
  list(
    nstart = as.integer(nstart), max_iter = as.integer(max_iter), tol = tol,
    select_nstart = as.integer(select_nstart),
    select_max_iter = as.integer(select_max_iter),
    seed = seed, standardize = standardize, G_rule = G_rule, G_grid = G_grid,
    max_G = max_G, log_N = log_N,
    min_group_size = as.integer(min_group_size),
    noise_multiplier = noise_multiplier
  )
}


.pk_log_N_value <- function(log_N, N) {
  if (identical(log_N, "log")) {
    value <- log(N)
    rule <- "log"
  } else if (identical(log_N, "loglog")) {
    value <- log(log(N))
    rule <- "loglog"
  } else {
    value <- as.numeric(log_N)
    rule <- "user"
  }
  if (!is.finite(value) || value <= 0) {
    stop(
      "pk_log_N produces a nonpositive or non-finite value at N = ", N,
      ". Use 'log' or supply a positive finite numeric value.",
      call. = FALSE
    )
  }
  list(value = value, rule = rule)
}


.pk_candidate_grid <- function(N, min_group_size, grid = NULL, max_G = NULL) {
  feasible_max_G <- floor(N / min_group_size)
  effective_max_G <- if (is.null(max_G)) {
    feasible_max_G
  } else {
    min(feasible_max_G, as.integer(max_G))
  }
  if (effective_max_G < 1L) {
    stop("No feasible projection K-means group count remains.")
  }
  if (!is.null(grid)) {
    grid <- grid[grid <= effective_max_G]
    if (!length(grid)) {
      stop("No value in pk_G_grid is at or below the effective maximum ",
           effective_max_G, ".")
    }
    return(grid)
  }
  if (effective_max_G <= 12L) return(seq_len(effective_max_G))
  tail_grid <- unique(as.integer(round(exp(seq(
    log(9), log(effective_max_G), length.out = 12L
  )))))
  sort(unique(c(seq_len(8L), tail_grid, effective_max_G)))
}


.pk_feature_object <- function(y, X, standardize = FALSE) {
  input <- .kmedoid_profile(y, X, standardize = standardize)
  Z <- input$data
  N <- nrow(Z)
  M <- ncol(Z)
  gram <- tcrossprod(Z) / M
  Z_squared <- Z^2

  second_moment <- mean(Z_squared)
  fourth_moment <- mean(Z_squared^2)
  # Direct analogue of the original K-means noise floor.  For i != s,
  # B[i,s] is the mean of W[i,s,m] = Z[i,m] Z[s,m].  The first component
  # averages the within-m residual variation of those product means.  The
  # second component conservatively covers the N diagonal self-products,
  # whose contribution to the normalized Gram-row risk is order 1/N.
  product_residual_ss <- sum(colSums(Z_squared)^2) - M * sum(gram^2)
  diagonal_residual_ss <- sum(
    sweep(Z_squared, 1L, diag(gram), FUN = "-")^2
  )
  off_diagonal_residual_ss <- max(
    product_residual_ss - diagonal_residual_ss, 0
  )
  off_diagonal_noise <- off_diagonal_residual_ss / (N^2 * M^2)
  diagonal_noise <- sum(diag(gram)^2) / N^2
  noise_floor <- max(
    off_diagonal_noise + diagonal_noise, sqrt(.Machine$double.eps)
  )

  list(
    data = gram,
    gram = gram,
    input = input,
    N = N,
    M = M,
    second_moment = second_moment,
    fourth_moment = fourth_moment,
    off_diagonal_noise = off_diagonal_noise,
    diagonal_noise = diagonal_noise,
    noise_floor = noise_floor
  )
}


.pk_seed_for_fit <- function(seed, grid_id = 0L, fold_id = 0L) {
  if (is.null(seed)) return(sample.int(.Machine$integer.max, 1L))
  as.integer((as.double(seed) + 1009 * fold_id + 37 * grid_id) %%
               .Machine$integer.max)
}


.pk_elbow_index <- function(grid, objective) {
  if (length(grid) <= 2L || diff(range(grid)) == 0 ||
      diff(range(objective)) == 0) {
    return(which.min(objective))
  }
  x <- (grid - min(grid)) / diff(range(grid))
  y <- (objective - min(objective)) / diff(range(objective))
  line <- 1 - x
  which.max(line - y)
}


.pk_fit_features <- function(feature, unit_cluster, options, fold_id = 0L,
                             evaluation_T = feature$input$T) {
  B <- as.matrix(feature$data)
  N <- nrow(B)
  max_feasible_G <- as.integer(floor(N / options$min_group_size))

  seed_existed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (seed_existed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (seed_existed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)

  fit_one <- function(G, grid_id, nstart = options$nstart,
                      max_iter = options$max_iter) {
    ps_kmeans(
      B, G = G, nstart = nstart,
      max_iter = max_iter, tol = options$tol,
      seed = .pk_seed_for_fit(options$seed, grid_id, fold_id)
    )
  }

  if (!is.null(unit_cluster)) {
    if (length(unit_cluster) != 1L || !is.finite(unit_cluster) ||
        unit_cluster < 1L || unit_cluster != floor(unit_cluster) ||
        unit_cluster > max_feasible_G) {
      stop("unit_cluster is outside the feasible projection K-means range.")
    }
    G <- as.integer(unit_cluster)
    fit <- fit_one(G, 1L)
    selection <- list(
      selected_G = G, rule = "fixed", path = NULL,
      threshold = NA_real_, noise_floor = feature$noise_floor,
      reference_G = NA_integer_, reference_objective = NA_real_,
      inference_multiplier = NA_real_,
      complexity_threshold = NA_real_,
      log_N_rule = if (is.character(options$log_N)) options$log_N else "user",
      log_N_value = NA_real_,
      requested_max_G = options$max_G,
      effective_max_G = max_feasible_G,
      theoretical_max_G = NA_integer_,
      max_G_source = "fixed unit_cluster",
      screening_nstart = NA_integer_, screening_max_iter = NA_integer_,
      screening_objective = fit$objective, final_objective = fit$objective
    )
    return(list(fit = fit, selection = selection, feature = feature))
  }

  log_N_spec <- .pk_log_N_value(options$log_N, N)
  log_N <- log_N_spec$value
  log_K <- log(max(feature$input$K, 2L))
  theoretical_max_G <- as.integer(floor(N / (log_K * log_N)))
  if (!is.null(options$max_G)) {
    effective_max_G <- as.integer(min(max_feasible_G, options$max_G))
    max_G_source <- "user"
  } else if (options$G_rule == "inference") {
    effective_max_G <- as.integer(min(max_feasible_G, theoretical_max_G))
    max_G_source <- "theoretical"
  } else {
    effective_max_G <- max_feasible_G
    max_G_source <- "group-size feasibility"
  }
  if (effective_max_G < 1L) {
    stop(
      "The default projection K-means inference cap leaves no admissible ",
      "group count. Supply pk_max_G to override the theoretical cap.",
      call. = FALSE
    )
  }
  grid <- .pk_candidate_grid(
    N, options$min_group_size, options$G_grid, effective_max_G
  )
  noise_threshold <- options$noise_multiplier * feature$noise_floor
  # A slower-vanishing complexity bound is still o(1), but is materially less
  # restrictive at the sample sizes used in the simulations.  The relative
  # loss tolerance tends to one, so the selector stays close to the best
  # empirical fit in the admissible class instead of accepting a factor-two
  # deterioration.
  complexity_threshold <- 1 / log_N
  inference_multiplier <- 1 + 1 / N

  if (options$G_rule == "n-third") {
    target <- max(1L, ceiling((N * evaluation_T)^(1 / 3)))
    selected_G <- grid[which.min(abs(grid - target))]
    fit <- fit_one(selected_G, which(grid == selected_G)[1L])
    selection <- list(
      selected_G = as.integer(selected_G), rule = "n-third",
      path = data.frame(
        G = grid, objective = NA_real_, noise_floor = feature$noise_floor,
        threshold = NA_real_, complexity = NA_real_,
        complexity_threshold = NA_real_, theoretically_safe = NA,
        effective_max_G = effective_max_G,
        max_G_source = max_G_source, feasible = NA
      ),
      threshold = NA_real_, noise_floor = feature$noise_floor,
      reference_G = NA_integer_, reference_objective = NA_real_,
      inference_multiplier = NA_real_,
      complexity_threshold = NA_real_,
      log_N_rule = log_N_spec$rule,
      log_N_value = log_N,
      requested_max_G = options$max_G,
      effective_max_G = effective_max_G,
      theoretical_max_G = theoretical_max_G,
      max_G_source = max_G_source,
      screening_nstart = options$nstart,
      screening_max_iter = options$max_iter,
      screening_objective = fit$objective, final_objective = fit$objective
    )
    return(list(fit = fit, selection = selection, feature = feature))
  }

  threshold <- if (options$G_rule == "inference") NA_real_ else noise_threshold
  grid_id <- seq_along(grid)
  complexity <- grid * log_K / N
  theoretically_safe <- complexity <= complexity_threshold

  fits <- vector("list", length(grid))
  objective <- rep(NA_real_, length(grid))
  screening_nstart <- if (options$G_rule == "inference") {
    min(options$nstart, options$select_nstart)
  } else {
    options$nstart
  }
  screening_max_iter <- if (options$G_rule == "inference") {
    min(options$max_iter, options$select_max_iter)
  } else {
    options$max_iter
  }
  last_index <- length(grid)
  for (index in seq_along(grid)) {
    fits[[index]] <- fit_one(
      grid[index], grid_id[index], nstart = screening_nstart,
      max_iter = screening_max_iter
    )
    objective[index] <- fits[[index]]$objective
    if (options$G_rule == "noise-floor" && objective[index] <= threshold) {
      last_index <- index
      break
    }
  }
  evaluated <- seq_len(last_index)
  grid_evaluated <- grid[evaluated]
  objective_evaluated <- objective[evaluated]

  complexity_evaluated <- complexity[evaluated]
  reference_G <- NA_integer_
  reference_objective <- NA_real_
  if (options$G_rule == "inference") {
    # Compare simpler candidates with the best empirical fit under the
    # applicable maximum-G cap. By default this is the theoretical cap;
    # pk_max_G lets the user replace it for simulation and sensitivity work.
    reference_index <- which.min(objective_evaluated)
    reference_G <- as.integer(grid_evaluated[reference_index])
    reference_objective <- objective_evaluated[reference_index]
    threshold <- inference_multiplier * reference_objective
    objective_feasible <- objective_evaluated <= threshold
    # With the default cap, every evaluated candidate obeys the theoretical
    # restriction. A user-supplied pk_max_G deliberately replaces that cap;
    # theoretically_safe below reports whether each candidate nevertheless
    # remains inside the default region.
    complexity_feasible <- rep(TRUE, length(evaluated))
    feasible <- objective_feasible & complexity_feasible
    fallback_used <- rep(FALSE, length(evaluated))
    selected_index <- which(feasible)[1L]
  } else if (options$G_rule == "noise-floor") {
    objective_feasible <- objective_evaluated <= threshold
    complexity_feasible <- rep(TRUE, length(evaluated))
    feasible <- objective_feasible
    fallback_used <- rep(FALSE, length(evaluated))
    if (any(feasible)) {
      selected_index <- which(feasible)[1L]
    } else {
      diagnostic <- objective_evaluated / threshold
      selected_index <- which.min(diagnostic)
      warning(
        "No candidate satisfied the projection K-means ", options$G_rule,
        " rule; using the evaluated candidate with the smallest diagnostic.",
        call. = FALSE
      )
    }
  } else {
    selected_index <- .pk_elbow_index(grid_evaluated, objective_evaluated)
    feasible <- objective_evaluated <= threshold
    objective_feasible <- rep(NA, length(evaluated))
    complexity_feasible <- rep(NA, length(evaluated))
    fallback_used <- rep(FALSE, length(evaluated))
  }

  selection <- list(
    selected_G = as.integer(grid_evaluated[selected_index]),
    rule = options$G_rule,
    path = data.frame(
      G = grid_evaluated,
      objective = objective_evaluated,
      noise_floor = feature$noise_floor,
      threshold = threshold,
      complexity = complexity_evaluated,
      complexity_threshold = complexity_threshold,
      theoretically_safe = theoretically_safe[evaluated],
      effective_max_G = effective_max_G,
      max_G_source = max_G_source,
      reference_G = reference_G,
      reference_objective = reference_objective,
      inference_multiplier = if (options$G_rule == "inference") {
        inference_multiplier
      } else {
        NA_real_
      },
      objective_feasible = objective_feasible,
      complexity_feasible = complexity_feasible,
      feasible = feasible,
      fallback_used = fallback_used
    ),
    threshold = threshold,
    noise_floor = feature$noise_floor,
    reference_G = reference_G,
    reference_objective = reference_objective,
    inference_multiplier = if (options$G_rule == "inference") {
      inference_multiplier
    } else {
      NA_real_
    },
    complexity_threshold = complexity_threshold,
    log_N_rule = log_N_spec$rule,
    log_N_value = log_N,
    requested_max_G = options$max_G,
    effective_max_G = effective_max_G,
    theoretical_max_G = theoretical_max_G,
    max_G_source = max_G_source,
    screening_nstart = screening_nstart,
    screening_max_iter = screening_max_iter
  )

  selected_fit <- fits[[evaluated[selected_index]]]
  selection$screening_objective <- selected_fit$objective
  if (options$G_rule == "inference" &&
      (screening_nstart < options$nstart ||
       screening_max_iter < options$max_iter)) {
    refit <- fit_one(
      grid_evaluated[selected_index], grid_id[evaluated[selected_index]],
      nstart = options$nstart, max_iter = options$max_iter
    )
    # The final fit is never allowed to be worse than the screening fit.
    if (refit$objective <= selected_fit$objective) selected_fit <- refit
  }
  selection$final_objective <- selected_fit$objective
  list(
    fit = selected_fit, selection = selection,
    feature = feature
  )
}


.HP_estimate_projection_kmeans_full <- function(
    data, y_col, covariate_cols, id_col, time_col, unit_cluster,
    pk_nstart, pk_max_iter, pk_select_nstart, pk_select_max_iter,
    pk_tol, pk_seed, pk_standardize,
    pk_G_rule, pk_G_grid, pk_max_G, pk_log_N, pk_min_group_size,
    pk_noise_multiplier) {
  required <- unique(c(id_col, time_col, y_col, covariate_cols))
  missing_columns <- setdiff(required, names(data))
  if (length(missing_columns)) {
    stop("Projection K-means data are missing columns: ",
         paste(missing_columns, collapse = ", "), ".")
  }
  if (length(y_col) != 1L || length(covariate_cols) < 2L) {
    stop("Projection K-means requires one outcome and treatment followed by at least one control.")
  }
  ids <- sort(unique(data[[id_col]]))
  times <- sort(unique(data[[time_col]]))
  N <- length(ids)
  TT <- length(times)
  K <- length(covariate_cols)
  if (N < 2L || TT < 2L || nrow(data) != N * TT ||
      anyDuplicated(data[c(id_col, time_col)])) {
    stop("Projection K-means requires a complete balanced panel with N >= 2 and T >= 2.")
  }
  options <- .pk_validate_options(
    pk_nstart, pk_max_iter, pk_select_nstart, pk_select_max_iter,
    pk_tol, pk_seed, pk_standardize,
    pk_G_rule, pk_G_grid, pk_min_group_size, pk_noise_multiplier,
    max_G = pk_max_G, log_N = pk_log_N
  )

  ii <- match(data[[id_col]], ids)
  tt <- match(data[[time_col]], times)
  panel_index <- cbind(ii, tt)
  y <- matrix(NA_real_, N, TT)
  y[panel_index] <- data[[y_col]]
  X <- array(NA_real_, c(N, TT, K))
  for (k in seq_len(K)) {
    X[cbind(ii, tt, k)] <- data[[covariate_cols[k]]]
  }

  feature <- .pk_feature_object(y, X, options$standardize)
  fitted <- .pk_fit_features(feature, unit_cluster, options)
  pk_fit <- fitted$fit
  G <- pk_fit$G
  group <- pk_fit$group
  projected <- .kmedoid_project_panel(y, X, group, G)
  tY <- as.vector(projected[, , 1L])
  tX <- do.call(cbind, lapply(seq_len(K), function(k) {
    as.vector(projected[, , k + 1L])
  }))
  colnames(tX) <- covariate_cols
  treatment <- tX[, 1L]
  controls <- tX[, -1L, drop = FALSE]

  fit <- hdm::rlassoEffect(
    x = controls, y = tY, d = treatment, method = "double selection"
  )
  trans <- data.frame(y = tY, D = treatment, controls, check.names = FALSE)
  Ytilde <- hdm::rlasso(y ~ . - D - 1, data = trans)$residuals
  Dtilde <- hdm::rlasso(
    D ~ . - 1, data = trans[, -1L, drop = FALSE]
  )$residuals
  data_res <- data.frame(
    id = rep(ids, TT), time = rep(times, each = N),
    Ytilde = Ytilde, Dtilde = Dtilde
  )
  post_plm <- plm::plm(
    Ytilde ~ -1 + Dtilde, data = data_res,
    model = "pooling", index = c("id", "time")
  )
  robust_se <- sqrt(diag(plm::vcovHC(
    post_plm, type = "HC0", method = "arellano"
  )))
  coefs <- stats::coef(post_plm)
  df <- max(1, N * TT - TT * G)
  se_corrected <- robust_se * sqrt(N * TT / df)
  summary_table <- data.frame(
    Estimate = coefs,
    `Std. Error corrected` = se_corrected,
    `t-value corrected` = coefs / se_corrected,
    `Pr(>|t|) corrected` = 2 * stats::pt(-abs(coefs / se_corrected), df),
    check.names = FALSE
  )

  list(
    cluster_method = "projection-kmeans",
    cross_fitted = FALSE,
    cross_type = "none",
    G_unit = G,
    unit_group = group,
    fit_summary = summary(fit),
    post_plm_summary = summary(post_plm),
    estimate_corrected = summary_table,
    summary_table = summary_table,
    pk_input_dimension = feature$M,
    pk_feature_dimension = feature$N,
    pk_objective = pk_fit$objective,
    pk_group_sizes = pk_fit$group_sizes,
    pk_converged = pk_fit$converged,
    pk_iterations = pk_fit$iterations,
    pk_selection_rule = fitted$selection$rule,
    pk_selection_path = fitted$selection$path,
    pk_selection_threshold = fitted$selection$threshold,
    pk_reference_G = fitted$selection$reference_G,
    pk_reference_objective = fitted$selection$reference_objective,
    pk_inference_multiplier = fitted$selection$inference_multiplier,
    pk_complexity_threshold = fitted$selection$complexity_threshold,
    pk_log_N_rule = fitted$selection$log_N_rule,
    pk_log_N_value = fitted$selection$log_N_value,
    pk_requested_max_G = fitted$selection$requested_max_G,
    pk_effective_max_G = fitted$selection$effective_max_G,
    pk_theoretical_max_G = fitted$selection$theoretical_max_G,
    pk_max_G_source = fitted$selection$max_G_source,
    pk_selection_nstart = fitted$selection$screening_nstart,
    pk_selection_max_iter = fitted$selection$screening_max_iter,
    pk_screening_objective = fitted$selection$screening_objective,
    pk_final_objective = fitted$selection$final_objective,
    pk_noise_floor = fitted$selection$noise_floor,
    pk_profile_center = feature$input$center,
    pk_profile_scale = feature$input$scale
  )
}
