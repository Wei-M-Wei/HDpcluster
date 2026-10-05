.HP_estimate_hac_new_full <- function(
    data,
    y_col,
    covariate_cols,
    id_col,
    time_col,
    link,
    unit_cluster,
    hac_G_rule,
    hac_G_grid,
    hac_cv_folds,
    hac_cv_repeats,
    hac_min_group_size,
    hac_reference_fraction,
    hac_seed,
    hac_standardize) {
  if (!identical(link, "average")) {
    stop("HAC new uses average linkage; set link = 'average'.")
  }
  hac_G_rule <- match.arg(hac_G_rule, c("one-se", "minimum"))

  if (length(y_col) != 1L || !is.character(y_col) || !nzchar(y_col)) {
    stop("HAC new requires y_col to name exactly one outcome column.")
  }
  required_columns <- unique(c(id_col, time_col, y_col, covariate_cols))
  missing_columns <- setdiff(required_columns, names(data))
  if (length(missing_columns)) {
    stop(
      "HAC new data are missing columns: ",
      paste(missing_columns, collapse = ", "),
      "."
    )
  }
  if (length(covariate_cols) < 2L) {
    stop(
      "HAC new requires the treatment followed by at least one control ",
      "in covariate_cols."
    )
  }
  if (anyDuplicated(data[c(id_col, time_col)])) {
    stop("HAC new requires exactly one row for each (unit, time) pair.")
  }

  ids <- sort(unique(data[[id_col]]))
  times <- sort(unique(data[[time_col]]))
  N <- length(ids)
  T <- length(times)
  K <- length(covariate_cols)
  if (N < 4L) {
    stop("Full-sample HAC new requires at least four units.")
  }
  if (T <= 1L) {
    stop("Full-sample HAC new is implemented only for T > 1.")
  }
  if (nrow(data) != N * T) {
    stop("HAC new requires a complete balanced panel.")
  }

  if (!is.null(hac_G_grid)) {
    if (!is.numeric(hac_G_grid) || !length(hac_G_grid) ||
        any(!is.finite(hac_G_grid)) || any(hac_G_grid < 1) ||
        any(hac_G_grid != floor(hac_G_grid))) {
      stop("hac_G_grid must be NULL or a nonempty vector of positive integers.")
    }
    hac_G_grid <- sort(unique(as.integer(hac_G_grid)))
  }
  for (parameter in c(
    "hac_cv_folds", "hac_cv_repeats", "hac_min_group_size"
  )) {
    value <- get(parameter)
    if (length(value) != 1L || !is.finite(value) ||
        value < 1 || value != floor(value)) {
      stop(parameter, " must be one positive integer.")
    }
  }
  hac_cv_folds <- as.integer(hac_cv_folds)
  hac_cv_repeats <- as.integer(hac_cv_repeats)
  hac_min_group_size <- as.integer(hac_min_group_size)
  if (hac_cv_folds < 2L) {
    stop("hac_cv_folds must be at least 2.")
  }
  if (hac_cv_folds > N) {
    stop("hac_cv_folds cannot exceed the number of units.")
  }
  if (length(hac_reference_fraction) != 1L ||
      !is.finite(hac_reference_fraction) ||
      hac_reference_fraction <= 0 || hac_reference_fraction >= 1) {
    stop("hac_reference_fraction must be strictly between zero and one.")
  }
  if (!is.null(hac_seed) &&
      (length(hac_seed) != 1L || !is.finite(hac_seed) ||
       hac_seed < 0 || hac_seed > .Machine$integer.max ||
       hac_seed != floor(hac_seed))) {
    stop("hac_seed must be NULL or one nonnegative integer.")
  }
  if (length(hac_standardize) != 1L || is.na(hac_standardize) ||
      !is.logical(hac_standardize)) {
    stop("hac_standardize must be TRUE or FALSE.")
  }

  y <- matrix(NA_real_, nrow = N, ncol = T)
  X <- array(NA_real_, dim = c(N, T, K))
  id_index <- match(data[[id_col]], ids)
  time_index <- match(data[[time_col]], times)
  for (row in seq_len(nrow(data))) {
    y[id_index[row], time_index[row]] <- data[[y_col]][row]
    for (k in seq_len(K)) {
      X[id_index[row], time_index[row], k] <- data[[covariate_cols[k]]][row]
    }
  }
  if (any(!is.finite(y)) || any(!is.finite(X))) {
    stop("HAC new requires finite outcome and covariate values.")
  }

  M <- T * (K + 1L)
  profile_raw <- matrix(NA_real_, nrow = N, ncol = M)
  colnames(profile_raw) <- character(M)
  component_names <- c(y_col, covariate_cols)
  for (tt in seq_len(T)) {
    columns <- (tt - 1L) * (K + 1L) + seq_len(K + 1L)
    profile_raw[, columns] <- cbind(
      y[, tt],
      matrix(X[, tt, , drop = FALSE], nrow = N, ncol = K)
    )
    colnames(profile_raw)[columns] <- paste0(
      "t=", times[tt], "::", component_names
    )
  }

  if (hac_standardize) {
    profile_center <- apply(profile_raw, 2L, stats::median)
    profile_scale <- apply(profile_raw, 2L, stats::mad)
    bad_scale <- !is.finite(profile_scale) |
      profile_scale <= sqrt(.Machine$double.eps)
    if (any(bad_scale)) {
      profile_scale[bad_scale] <- apply(
        profile_raw[, bad_scale, drop = FALSE],
        2L,
        stats::sd
      )
    }
    positive_scale <- profile_scale[
      is.finite(profile_scale) &
        profile_scale > sqrt(.Machine$double.eps)
    ]
    reference_scale <- if (length(positive_scale)) {
      stats::median(positive_scale)
    } else {
      1
    }
    profile_scale[
      !is.finite(profile_scale) |
        profile_scale <= sqrt(.Machine$double.eps)
    ] <- reference_scale
    relative_scale <- pmin(
      pmax(profile_scale / reference_scale, 0.25),
      4
    )
    profile_scale <- reference_scale * relative_scale
    profile <- sweep(
      sweep(profile_raw, 2L, profile_center, "-"),
      2L,
      profile_scale,
      "/"
    )
  } else {
    profile_center <- rep(0, M)
    profile_scale <- rep(1, M)
    profile <- profile_raw
  }

  seed_existed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (!is.null(hac_seed)) {
    if (seed_existed) {
      old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    }
    on.exit({
      if (seed_existed) {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(as.integer(hac_seed))
  }

  split_reference <- function(n) {
    if (n < 2L) {
      stop("HAC new needs at least two fitting units before the reference split.")
    }
    permutation <- sample.int(n)
    n_reference <- max(
      1L,
      min(n - 1L, round(n * hac_reference_fraction))
    )
    list(
      reference = permutation[seq_len(n_reference)],
      auxiliary = permutation[-seq_len(n_reference)]
    )
  }

  capacities <- function(n, G) {
    result <- rep(n %/% G, G)
    remainder <- n %% G
    if (remainder > 0L) {
      result[seq_len(remainder)] <- result[seq_len(remainder)] + 1L
    }
    result
  }

  balanced_assignment <- function(distance_matrix, group_capacities) {
    n <- nrow(distance_matrix)
    G <- ncol(distance_matrix)
    if (length(group_capacities) != G ||
        sum(group_capacities) != n) {
      stop("HAC capacities do not match the assignment problem.")
    }

    if (requireNamespace("lpSolve", quietly = TRUE)) {
      solution <- lpSolve::lp.transport(
        cost.mat = distance_matrix,
        direction = "min",
        row.signs = rep("=", n),
        row.rhs = rep(1, n),
        col.signs = rep("=", G),
        col.rhs = group_capacities
      )
      if (solution$status != 0L) {
        stop("The balanced HAC transportation problem did not solve successfully.")
      }
      return(max.col(solution$solution, ties.method = "first"))
    }

    if (requireNamespace("clue", quietly = TRUE)) {
      group_slots <- rep(seq_len(G), group_capacities)
      slot_cost <- distance_matrix[, group_slots, drop = FALSE]
      assigned_slot <- as.integer(
        clue::solve_LSAP(slot_cost, maximum = FALSE)
      )
      return(group_slots[assigned_slot])
    }

    warning(
      "Neither lpSolve nor clue is installed; using a greedy balanced HAC assignment.",
      call. = FALSE
    )
    group <- integer(n)
    remaining <- group_capacities
    for (pair in order(distance_matrix)) {
      i <- ((pair - 1L) %% n) + 1L
      g <- ((pair - 1L) %/% n) + 1L
      if (group[i] == 0L && remaining[g] > 0L) {
        group[i] <- g
        remaining[g] <- remaining[g] - 1L
      }
      if (all(group > 0L)) break
    }
    if (any(group == 0L)) {
      for (i in which(group == 0L)) {
        available <- which(remaining > 0L)
        g <- available[which.min(distance_matrix[i, available])]
        group[i] <- g
        remaining[g] <- remaining[g] - 1L
      }
    }
    group
  }

  projection_scores <- function(A, reference) {
    tcrossprod(A, reference) / ncol(A)
  }

  score_distances <- function(scores, medoid_scores) {
    distance <- outer(
      rowSums(scores^2),
      rowSums(medoid_scores^2),
      "+"
    ) - 2 * tcrossprod(scores, medoid_scores)
    pmax(distance / ncol(scores), 0)
  }

  fit_hac <- function(auxiliary, reference, G) {
    auxiliary <- as.matrix(auxiliary)
    reference <- as.matrix(reference)
    n_auxiliary <- nrow(auxiliary)
    if (G < 1L || G > n_auxiliary) {
      stop(
        "The HAC group count must be between one and the auxiliary-HAC ",
        "sample size."
      )
    }
    scores <- projection_scores(auxiliary, reference)
    q_matrix <- if (n_auxiliary == 1L) {
      matrix(0, 1L, 1L)
    } else {
      as.matrix(stats::dist(scores))^2 / nrow(reference)
    }
    tree <- if (n_auxiliary == 1L) {
      NULL
    } else {
      stats::hclust(stats::as.dist(q_matrix), method = "average")
    }
    group <- if (G == 1L) {
      rep(1L, n_auxiliary)
    } else {
      stats::cutree(tree, k = G)
    }
    medoid_index <- vapply(seq_len(G), function(g) {
      members <- which(group == g)
      members[
        which.min(rowSums(q_matrix[members, members, drop = FALSE]))
      ]
    }, integer(1L))
    merge_height <- if (is.null(tree) || G == n_auxiliary) {
      0
    } else {
      tree$height[n_auxiliary - G]
    }
    list(
      G = G,
      group = as.integer(group),
      tree = tree,
      medoid_index = medoid_index,
      medoid_scores = scores[medoid_index, , drop = FALSE],
      reference = reference,
      merge_height = merge_height,
      group_sizes = tabulate(group, nbins = G)
    )
  }

  assign_hac <- function(evaluation, fit) {
    scores <- projection_scores(evaluation, fit$reference)
    distance <- score_distances(scores, fit$medoid_scores)
    group <- balanced_assignment(
      distance,
      capacities(nrow(evaluation), fit$G)
    )
    list(
      group = as.integer(group),
      distance = distance[cbind(seq_len(nrow(evaluation)), group)],
      distances = distance
    )
  }

  grid_from_max <- function(max_G) {
    if (max_G < 1L) {
      stop("No feasible HAC new group count remains.")
    }
    if (!is.null(hac_G_grid)) {
      if (any(hac_G_grid > max_G)) {
        stop(
          "hac_G_grid exceeds the feasible maximum ",
          max_G,
          "; reduce the grid, hac_cv_folds, or hac_min_group_size."
        )
      }
      return(hac_G_grid)
    }
    if (max_G <= 12L) return(seq_len(max_G))
    tail_grid <- unique(as.integer(round(exp(seq(
      log(9),
      log(max_G),
      length.out = 12L
    )))))
    sort(unique(c(seq_len(8L), tail_grid, max_G)))
  }

  select_G <- function() {
    designs <- list()
    design_id <- 0L
    for (repeat_id in seq_len(hac_cv_repeats)) {
      permutation <- sample.int(N)
      fold_label <- integer(N)
      fold_label[permutation] <- rep(
        seq_len(hac_cv_folds),
        length.out = N
      )
      for (fold in seq_len(hac_cv_folds)) {
        validation <- which(fold_label == fold)
        fitting <- which(fold_label != fold)
        split <- split_reference(length(fitting))
        design_id <- design_id + 1L
        designs[[design_id]] <- list(
          validation = validation,
          auxiliary = fitting[split$auxiliary],
          reference = fitting[split$reference]
        )
      }
    }

    max_G <- min(
      vapply(designs, function(x) length(x$auxiliary), integer(1L)),
      vapply(designs, function(x) length(x$validation), integer(1L)),
      floor(N / hac_min_group_size)
    )
    grid <- grid_from_max(max_G)
    loss_sum <- matrix(0, nrow = N, ncol = length(grid))
    loss_count <- matrix(0L, nrow = N, ncol = length(grid))
    for (design in designs) {
      for (grid_id in seq_along(grid)) {
        fit <- fit_hac(
          profile[design$auxiliary, , drop = FALSE],
          profile[design$reference, , drop = FALSE],
          grid[grid_id]
        )
        assigned <- assign_hac(
          profile[design$validation, , drop = FALSE],
          fit
        )
        loss_sum[design$validation, grid_id] <-
          loss_sum[design$validation, grid_id] + assigned$distance
        loss_count[design$validation, grid_id] <-
          loss_count[design$validation, grid_id] + 1L
      }
    }
    if (any(loss_count == 0L)) {
      stop("Internal HAC validation loss was not filled.")
    }
    unit_loss <- loss_sum / loss_count
    cv_loss <- colMeans(unit_loss)
    cv_se <- apply(unit_loss, 2L, stats::sd) / sqrt(N)
    cv_se[!is.finite(cv_se)] <- 0
    minimizer_index <- which.min(cv_loss)
    threshold <- cv_loss[minimizer_index] + cv_se[minimizer_index]
    selected_index <- if (hac_G_rule == "minimum") {
      minimizer_index
    } else {
      which(cv_loss <= threshold)[1L]
    }
    list(
      selected_G = as.integer(grid[selected_index]),
      rule = hac_G_rule,
      path = data.frame(G = grid, cv_loss = cv_loss, cv_se = cv_se),
      threshold = threshold,
      minimizer_G = as.integer(grid[minimizer_index])
    )
  }

  final_split <- split_reference(N)
  max_final_G <- min(
    length(final_split$auxiliary),
    floor(N / hac_min_group_size)
  )
  if (is.null(unit_cluster)) {
    selection <- select_G()
    G_unit <- selection$selected_G
  } else {
    if (length(unit_cluster) != 1L || !is.finite(unit_cluster) ||
        unit_cluster < 1 || unit_cluster != floor(unit_cluster)) {
      stop("For HAC new, unit_cluster must be one positive integer.")
    }
    G_unit <- as.integer(unit_cluster)
    selection <- list(
      rule = "fixed",
      path = NULL,
      threshold = NA_real_,
      minimizer_G = NA_integer_
    )
  }
  if (G_unit > max_final_G) {
    stop(
      "For full-sample HAC new, unit_cluster must not exceed ",
      max_final_G,
      ", as implied by the auxiliary split and hac_min_group_size."
    )
  }

  hac_fit <- fit_hac(
    profile[final_split$auxiliary, , drop = FALSE],
    profile[final_split$reference, , drop = FALSE],
    G_unit
  )
  assigned <- assign_hac(profile, hac_fit)
  unit_group <- assigned$group
  group_sizes <- tabulate(unit_group, nbins = G_unit)

  Du <- matrix(0, nrow = N, ncol = G_unit)
  for (g in seq_len(G_unit)) {
    Du[, g] <- as.numeric(unit_group == g)
  }
  Diagu <- diag(1 / colSums(Du), nrow = G_unit, ncol = G_unit)
  Mu <- diag(N) - Du %*% Diagu %*% t(Du)

  Z <- array(NA_real_, dim = c(N, T, K + 1L))
  Z[, , 1L] <- y
  for (k in seq_len(K)) {
    Z[, , k + 1L] <- X[, , k]
  }
  Z_projected <- array(NA_real_, dim = dim(Z))
  for (k in seq_len(K + 1L)) {
    Z_projected[, , k] <- Mu %*% Z[, , k]
  }
  tY_vector <- as.vector(Z_projected[, , 1L])
  tX_matrix <- do.call(cbind, lapply(seq_len(K), function(k) {
    as.vector(Z_projected[, , k + 1L])
  }))
  colnames(tX_matrix) <- covariate_cols
  d_vector <- tX_matrix[, 1L]
  x_matrix <- tX_matrix[, -1L, drop = FALSE]

  fit <- hdm::rlassoEffect(
    x = x_matrix,
    y = tY_vector,
    d = d_vector,
    method = "double selection"
  )
  trans <- data.frame(
    y = tY_vector,
    D = d_vector,
    x_matrix,
    check.names = FALSE
  )
  lasso_y <- hdm::rlasso(y ~ . - D - 1, data = trans)
  Ytilde <- lasso_y$residuals
  lasso_d <- hdm::rlasso(
    D ~ . - 1,
    data = trans[, -1L, drop = FALSE]
  )
  Dtilde <- lasso_d$residuals
  data_res <- data.frame(
    id = rep(ids, T),
    time = rep(times, each = N),
    Ytilde = Ytilde,
    Dtilde = Dtilde
  )
  post_plm <- plm::plm(
    Ytilde ~ -1 + Dtilde,
    data = data_res,
    model = "pooling",
    index = c("id", "time")
  )
  robust_se <- sqrt(diag(plm::vcovHC(
    post_plm,
    type = "HC0",
    method = "arellano"
  )))
  coefs <- stats::coef(post_plm)
  df_hac <- max(1, N * T - T * G_unit)
  se_corrected <- robust_se * sqrt(N * T / df_hac)
  t_values_corrected <- coefs / se_corrected
  p_values_corrected <- 2 * stats::pt(
    -abs(t_values_corrected),
    df_hac
  )
  summary_table_correct <- data.frame(
    Estimate = coefs,
    `Std. Error corrected` = se_corrected,
    `t-value corrected` = t_values_corrected,
    `Pr(>|t|) corrected` = p_values_corrected,
    check.names = FALSE
  )

  evaluation_coverage <- max(assigned$distance)
  effective_cut <- max(hac_fit$merge_height, evaluation_coverage)
  medoid_auxiliary_index <-
    final_split$auxiliary[hac_fit$medoid_index]

  list(
    cluster_method = "HAC new",
    cross_fitted = FALSE,
    cross_type = "none",
    G_unit = G_unit,
    unit_group = unit_group,
    fit_summary = summary(fit),
    post_plm_summary = summary(post_plm),
    estimate_corrected = summary_table_correct,
    summary_table = summary_table_correct,
    hac_profile_dimension = M,
    hac_reference_size = length(final_split$reference),
    hac_auxiliary_size = length(final_split$auxiliary),
    hac_reference_unit_id = ids[final_split$reference],
    hac_auxiliary_unit_id = ids[final_split$auxiliary],
    hac_merge_height = hac_fit$merge_height,
    hac_evaluation_coverage = evaluation_coverage,
    hac_effective_cut = effective_cut,
    hac_stochastic_rate = log(N) / M +
      log(N) / length(final_split$reference),
    hac_auxiliary_group_sizes = hac_fit$group_sizes,
    hac_group_sizes = group_sizes,
    hac_medoid_auxiliary_index = medoid_auxiliary_index,
    hac_medoid_unit_id = ids[medoid_auxiliary_index],
    hac_selection_rule = selection$rule,
    hac_selection_path = selection$path,
    hac_selection_threshold = selection$threshold,
    hac_selection_minimizer_G = selection$minimizer_G,
    hac_profile_center = profile_center,
    hac_profile_scale = profile_scale
  )
}
