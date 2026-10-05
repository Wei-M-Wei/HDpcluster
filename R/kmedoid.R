.kmedoid_profile <- function(y, X, standardize = TRUE,
                             center = NULL, scale = NULL) {
  y <- as.matrix(y)
  X <- as.array(X)
  N <- nrow(y)
  TT <- ncol(y)
  K <- dim(X)[3L]
  if (!identical(dim(X)[1:2], c(N, TT))) {
    stop("The k-medoids outcome and covariate arrays are not conformable.")
  }

  profile_raw <- matrix(NA_real_, nrow = N, ncol = TT * (K + 1L))
  for (tt in seq_len(TT)) {
    columns <- (tt - 1L) * (K + 1L) + seq_len(K + 1L)
    profile_raw[, columns] <- cbind(
      y[, tt],
      matrix(X[, tt, , drop = FALSE], nrow = N, ncol = K)
    )
  }
  if (any(!is.finite(profile_raw))) {
    stop("K-medoids requires finite outcome and covariate values.")
  }

  if (!standardize) {
    return(list(
      data = profile_raw,
      raw = profile_raw,
      center = rep(0, ncol(profile_raw)),
      scale = rep(1, ncol(profile_raw)),
      M = ncol(profile_raw),
      T = TT,
      K = K
    ))
  }

  if (is.null(center) || is.null(scale)) {
    center <- apply(profile_raw, 2L, stats::median)
    scale <- apply(profile_raw, 2L, stats::mad)
    bad <- !is.finite(scale) | scale <= sqrt(.Machine$double.eps)
    if (any(bad)) {
      scale[bad] <- apply(profile_raw[, bad, drop = FALSE], 2L, stats::sd)
    }
    positive <- scale[is.finite(scale) & scale > sqrt(.Machine$double.eps)]
    fallback <- if (length(positive)) stats::median(positive) else 1
    scale[!is.finite(scale) | scale <= sqrt(.Machine$double.eps)] <- fallback
  }
  if (length(center) != ncol(profile_raw) ||
      length(scale) != ncol(profile_raw) || any(scale <= 0)) {
    stop("Invalid k-medoids profile centering or scaling values.")
  }
  profile <- sweep(sweep(profile_raw, 2L, center, "-"), 2L, scale, "/")
  list(data = profile, raw = profile_raw, center = center, scale = scale,
       M = ncol(profile), T = TT, K = K)
}

.kmedoid_projection_distance <- function(profile) {
  profile <- as.matrix(profile)
  N <- nrow(profile)
  M <- ncol(profile)
  if (N < 3L) {
    stop("The leave-two-out projection distance requires at least three units.")
  }
  if (M < 1L || any(!is.finite(profile))) {
    stop("The k-medoids profile must be a finite matrix with at least one column.")
  }

  gram <- tcrossprod(profile) / M
  # For every pair (i,j), first compute the squared Euclidean distance
  # between rows i and j of the Gram matrix. This equals the sum over all
  # reference units s of <Z_i-Z_j,Z_s>_M^2. Subtract the s=i and s=j terms
  # to obtain the leave-two-out distance. Matrix multiplication moves the
  # O(N^3) arithmetic from interpreted R loops into optimized BLAS code.
  gram_norm <- rowSums(gram^2)
  all_reference_sum <- outer(gram_norm, gram_norm, "+") -
    2 * tcrossprod(gram)
  gram_diagonal_by_row <- matrix(diag(gram), nrow = N, ncol = N)
  own_i_term <- (gram_diagonal_by_row - gram)^2
  own_j_term <- (gram - t(gram_diagonal_by_row))^2
  distance <- (all_reference_sum - own_i_term - own_j_term) / (N - 2L)
  distance <- pmax(distance, 0)
  diag(distance) <- 0
  distance <- (distance + t(distance)) / 2
  dimnames(distance) <- list(rownames(profile), rownames(profile))
  distance
}

.kmedoid_capacities <- function(N, G) {
  capacities <- rep(N %/% G, G)
  remainder <- N %% G
  if (remainder > 0L) capacities[seq_len(remainder)] <-
    capacities[seq_len(remainder)] + 1L
  capacities
}

.kmedoid_capacity_variants <- function(N, G) {
  small <- N %/% G
  remainder <- N %% G
  if (remainder == 0L) return(list(rep(small, G)))
  large_groups <- utils::combn(G, remainder, simplify = FALSE)
  lapply(large_groups, function(index) {
    result <- rep(small, G)
    result[index] <- result[index] + 1L
    result
  })
}

.kmedoid_balanced_assignment <- function(cost, capacities, medoids) {
  cost <- as.matrix(cost)
  N <- nrow(cost)
  G <- ncol(cost)
  if (length(capacities) != G || sum(capacities) != N ||
      any(capacities < 1L) || length(medoids) != G) {
    stop("Invalid balanced k-medoids assignment problem.")
  }
  nonmedoid <- setdiff(seq_len(N), medoids)
  remaining <- capacities - 1L
  group <- integer(N)
  group[medoids] <- seq_len(G)
  if (!length(nonmedoid)) return(group)
  reduced_cost <- cost[nonmedoid, , drop = FALSE]

  if (requireNamespace("lpSolve", quietly = TRUE)) {
    solution <- lpSolve::lp.transport(
      cost.mat = reduced_cost,
      direction = "min",
      row.signs = rep("=", length(nonmedoid)),
      row.rhs = rep(1, length(nonmedoid)),
      col.signs = rep("=", G),
      col.rhs = remaining
    )
    if (solution$status == 0L) {
      group[nonmedoid] <- max.col(solution$solution, ties.method = "first")
      return(group)
    }
  }

  if (requireNamespace("clue", quietly = TRUE)) {
    slots <- rep(seq_len(G), remaining)
    assigned_slot <- as.integer(clue::solve_LSAP(
      reduced_cost[, slots, drop = FALSE], maximum = FALSE
    ))
    group[nonmedoid] <- slots[assigned_slot]
    return(group)
  }

  warning(
    "Neither lpSolve nor clue is available; using a greedy balanced assignment.",
    call. = FALSE
  )
  available <- remaining
  for (index in order(reduced_cost)) {
    row <- ((index - 1L) %% nrow(reduced_cost)) + 1L
    g <- ((index - 1L) %/% nrow(reduced_cost)) + 1L
    unit <- nonmedoid[row]
    if (group[unit] == 0L && available[g] > 0L) {
      group[unit] <- g
      available[g] <- available[g] - 1L
    }
    if (all(group > 0L)) break
  }
  group
}

.kmedoid_bounded_assignment <- function(cost, medoids) {
  cost <- as.matrix(cost)
  N <- nrow(cost)
  G <- ncol(cost)
  nonmedoid <- setdiff(seq_len(N), medoids)
  group <- integer(N)
  group[medoids] <- seq_len(G)
  if (!length(nonmedoid)) {
    return(list(
      group = group, capacities = rep(1L, G), objective = 0
    ))
  }

  lower_size <- N %/% G
  upper_size <- ceiling(N / G)
  lower_remaining <- rep.int(max(lower_size - 1L, 0L), G)
  upper_remaining <- rep.int(max(upper_size - 1L, 0L), G)
  reduced_cost <- cost[nonmedoid, , drop = FALSE]

  # A square assignment representation gives every group upper_size - 1
  # nonmedoid slots. When N is not divisible by G, dummy rows occupy exactly
  # G - (N %% G) optional slots, leaving each group at either floor(N/G) or
  # ceiling(N/G). One LSAP therefore replaces choose(G, N %% G) separate
  # capacity-pattern solves.
  slots_per_group <- upper_remaining[1L]
  slot_group <- rep(seq_len(G), each = slots_per_group)
  n_real <- length(nonmedoid)
  n_slots <- length(slot_group)
  n_dummy <- n_slots - n_real

  if (requireNamespace("clue", quietly = TRUE)) {
    real_slot_cost <- reduced_cost[, slot_group, drop = FALSE]
    assignment_cost <- real_slot_cost
    if (n_dummy > 0L) {
      optional_slot <- seq.int(
        from = slots_per_group, by = slots_per_group, length.out = G
      )
      big_M <- (max(real_slot_cost) + 1) * (n_real + 1)
      dummy_cost <- matrix(big_M, nrow = n_dummy, ncol = n_slots)
      dummy_cost[, optional_slot] <- 0
      assignment_cost <- rbind(real_slot_cost, dummy_cost)
    }
    selected_slot <- as.integer(
      clue::solve_LSAP(assignment_cost, maximum = FALSE)
    )
    group[nonmedoid] <- slot_group[selected_slot[seq_len(n_real)]]
  } else if (requireNamespace("lpSolve", quietly = TRUE)) {
    n_variables <- n_real * G
    unit_constraints <- matrix(0, nrow = n_real, ncol = n_variables)
    unit_constraints[cbind(rep(seq_len(n_real), G), seq_len(n_variables))] <- 1
    group_constraints <- matrix(0, nrow = G, ncol = n_variables)
    for (g in seq_len(G)) {
      columns <- (g - 1L) * n_real + seq_len(n_real)
      group_constraints[g, columns] <- 1
    }
    solution <- lpSolve::lp(
      direction = "min",
      objective.in = as.vector(reduced_cost),
      const.mat = rbind(
        unit_constraints, group_constraints, group_constraints
      ),
      const.dir = c(
        rep("=", n_real), rep(">=", G), rep("<=", G)
      ),
      const.rhs = c(
        rep(1, n_real), lower_remaining, upper_remaining
      )
    )
    if (solution$status != 0L) {
      stop("The bounded balanced-assignment problem could not be solved.")
    }
    assignment_matrix <- matrix(solution$solution, nrow = n_real, ncol = G)
    group[nonmedoid] <- max.col(assignment_matrix, ties.method = "first")
  } else {
    warning(
      "Neither clue nor lpSolve is available; using a greedy balanced ",
      "assignment. Install clue for the fast exact assignment step.",
      call. = FALSE
    )
    real_slot_cost <- reduced_cost[, slot_group, drop = FALSE]
    available_slot <- rep(TRUE, n_slots)
    if (n_dummy > 0L) {
      optional_slot <- seq.int(
        from = slots_per_group, by = slots_per_group, length.out = G
      )
      available_slot[optional_slot[seq_len(n_dummy)]] <- FALSE
    }
    assigned_real <- rep(FALSE, n_real)
    for (index in order(real_slot_cost)) {
      row <- ((index - 1L) %% n_real) + 1L
      slot <- ((index - 1L) %/% n_real) + 1L
      if (!assigned_real[row] && available_slot[slot]) {
        group[nonmedoid[row]] <- slot_group[slot]
        assigned_real[row] <- TRUE
        available_slot[slot] <- FALSE
      }
      if (all(assigned_real)) break
    }
  }

  capacities <- tabulate(group, nbins = G)
  if (any(capacities < lower_size) || any(capacities > upper_size)) {
    stop("The bounded assignment did not produce balanced group sizes.")
  }
  objective <- mean(cost[cbind(seq_len(N), group)])
  list(group = group, capacities = capacities, objective = objective)
}

.kmedoid_assign_for_medoids <- function(distance, medoids) {
  cost <- distance[, medoids, drop = FALSE]
  .kmedoid_bounded_assignment(cost, medoids)
}

.kmedoid_initial_medoids <- function(distance, G) {
  N <- nrow(distance)
  medoids <- sample.int(N, 1L)
  while (length(medoids) < G) {
    nearest <- apply(distance[, medoids, drop = FALSE], 1L, min)
    nearest[medoids] <- -Inf
    candidates <- which(nearest == max(nearest))
    # Index into candidates explicitly. In base R, sample(x, 1) treats a
    # length-one numeric x as 1:x, which can accidentally repeat a medoid.
    medoids <- c(medoids, candidates[sample.int(length(candidates), 1L)])
  }
  medoids
}

.kmedoid_fit <- function(distance, G, nstart = 30L, max_iter = 100L,
                         tol = 1e-8, seed = NULL,
                         solver = c("heuristic", "auto", "exact"),
                         exact_limit = 50000L) {
  solver <- match.arg(solver)
  distance <- as.matrix(distance)
  N <- nrow(distance)
  if (ncol(distance) != N || any(!is.finite(distance)) ||
      any(distance < 0) || max(abs(distance - t(distance))) > 1e-8) {
    stop("distance must be a finite, nonnegative, symmetric matrix.")
  }
  if (length(G) != 1L || !is.finite(G) || G < 1L ||
      G > N || G != floor(G)) {
    stop("G must be one integer between one and the number of units.")
  }
  G <- as.integer(G)
  nstart <- as.integer(nstart)
  max_iter <- as.integer(max_iter)

  seed_existed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (!is.null(seed)) {
    if (seed_existed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit({
      if (seed_existed) {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(as.integer(seed))
  }

  remainder <- N %% G
  log_problems <- lchoose(N, G) + if (remainder > 0L) lchoose(G, remainder) else 0
  exact_feasible <- is.finite(log_problems) &&
    log_problems <= log(max(1, exact_limit))
  exact_assignment_available <-
    requireNamespace("lpSolve", quietly = TRUE) ||
    requireNamespace("clue", quietly = TRUE)
  use_exact <- solver == "exact" ||
    (solver == "auto" && exact_feasible && exact_assignment_available)
  if (solver == "exact" && !exact_feasible) {
    stop("The exact balanced k-medoids search exceeds kmedoid_exact_limit.")
  }
  if (solver == "exact" && !exact_assignment_available) {
    stop("The exact solver requires either the lpSolve or clue package.")
  }

  if (use_exact) {
    medoid_sets <- utils::combn(N, G, simplify = FALSE)
    best <- NULL
    for (medoids in medoid_sets) {
      candidate <- .kmedoid_assign_for_medoids(distance, medoids)
      if (is.null(best) || candidate$objective < best$objective) {
        best <- c(candidate, list(medoids = medoids))
      }
    }
    best$iterations <- 1L
    best$converged <- TRUE
    best$solver <- "exact"
    best$global_optimum <- TRUE
  } else {
    fits <- vector("list", nstart)
    for (start in seq_len(nstart)) {
      medoids <- .kmedoid_initial_medoids(distance, G)
      previous_objective <- Inf
      converged <- FALSE
      for (iteration in seq_len(max_iter)) {
        assigned <- .kmedoid_assign_for_medoids(distance, medoids)
        updated <- vapply(seq_len(G), function(g) {
          members <- which(assigned$group == g)
          members[which.min(rowSums(distance[members, members, drop = FALSE]))]
        }, integer(1L))
        relative_change <- abs(previous_objective - assigned$objective) /
          max(abs(previous_objective), .Machine$double.eps)
        if (identical(updated, medoids) ||
            (is.finite(relative_change) && relative_change < tol)) {
          converged <- TRUE
          medoids <- updated
          break
        }
        medoids <- updated
        previous_objective <- assigned$objective
      }
      assigned <- .kmedoid_assign_for_medoids(distance, medoids)
      fits[[start]] <- c(assigned, list(
        medoids = medoids, iterations = iteration, converged = converged,
        solver = "heuristic", global_optimum = FALSE
      ))
    }
    best <- fits[[which.min(vapply(fits, `[[`, numeric(1L), "objective"))]]
  }

  best$G <- G
  best$group_sizes <- tabulate(best$group, nbins = G)
  best
}

.kmedoid_candidate_grid <- function(N, min_group_size, grid = NULL) {
  max_G <- floor(N / min_group_size)
  if (max_G < 1L) stop("No feasible k-medoids group count remains.")
  if (!is.null(grid)) {
    if (any(grid > max_G)) {
      stop("kmedoid_G_grid exceeds the feasible maximum ", max_G, ".")
    }
    return(sort(unique(as.integer(grid))))
  }
  if (max_G <= 12L) return(seq_len(max_G))
  tail <- unique(as.integer(round(exp(seq(log(9), log(max_G), length.out = 12L)))))
  sort(unique(c(seq_len(8L), tail, max_G)))
}

.kmedoid_select_G <- function(distance, grid, rule, nstart, max_iter, tol,
                              seed, solver, exact_limit, TT, K) {
  N <- nrow(distance)
  if (length(TT) != 1L || !is.finite(TT) || TT < 1L || TT != floor(TT)) {
    stop("The k-medoids clustering-fold T must be one positive integer.")
  }
  if (length(K) != 1L || !is.finite(K) || K < 1L || K != floor(K)) {
    stop("The k-medoids covariate count K must be one positive integer.")
  }
  effective_n <- N * as.integer(TT)
  inference_threshold <- 1
  theoretical_threshold <- 1 / log(effective_n)
  complexity <- grid * log(max(as.integer(K), 2L)) / N
  # An elbow path generally contains combinatorially large candidate problems.
  # When G is data-driven, "exact" therefore means exact whenever the stated
  # limit permits it, with an explicit heuristic fallback for the remaining
  # candidates. A fixed unit_cluster retains strict exact-solver semantics.
  selection_solver <- if (solver == "exact") "auto" else solver
  if (rule == "n-third") {
    target <- max(1L, ceiling(effective_n^(1 / 3)))
    selected <- grid[which.min(abs(grid - target))]
    fit_seed <- if (is.null(seed)) NULL else
      as.integer((as.double(seed) + 104729) %% .Machine$integer.max)
    fit <- .kmedoid_fit(
      distance, selected, nstart, max_iter, tol, fit_seed,
      selection_solver, exact_limit
    )
    if (solver == "exact" && !isTRUE(fit$global_optimum)) {
      warning(
        "The selected G exceeds kmedoid_exact_limit; using the multi-start ",
        "heuristic for automatic G selection. Supply a fixed unit_cluster ",
        "and raise kmedoid_exact_limit if a certified exact fit is required.",
        call. = FALSE
      )
    }
    objective_path <- ifelse(grid == selected, fit$objective, NA_real_)
    root_remainder <- sqrt(effective_n) * objective_path
    group_remainder <- sqrt(complexity + objective_path)
    inference_index <- pmax(root_remainder, group_remainder)
    return(list(
      selected_G = as.integer(selected), rule = rule,
      path = data.frame(
        G = grid,
        objective = objective_path,
        elbow_score = NA_real_,
        complexity = complexity,
        root_remainder = root_remainder,
        group_remainder = group_remainder,
        inference_index = inference_index,
        inference_threshold = inference_threshold,
        inference_feasible = inference_index <= inference_threshold,
        theoretical_threshold = theoretical_threshold,
        theoretical_feasible = inference_index <= theoretical_threshold,
        solver = ifelse(grid == selected, fit$solver, NA_character_),
        global_optimum = ifelse(grid == selected, fit$global_optimum, NA)
      ),
      fit = fit
    ))
  }
  if (rule %in% c("inference", "theoretical")) {
    selected_threshold <- if (rule == "theoretical") {
      theoretical_threshold
    } else {
      inference_threshold
    }
    evaluated_fits <- vector("list", length(grid))
    evaluated_objective <- rep(NA_real_, length(grid))
    evaluated_solver <- rep(NA_character_, length(grid))
    evaluated_global <- rep(NA, length(grid))
    selected_index <- NA_integer_

    # Both rules select the smallest feasible G. Evaluating candidates in
    # increasing order and stopping at the first success is therefore exactly
    # equivalent to fitting the entire grid and selecting afterward.
    for (index in seq_along(grid)) {
      fit_seed <- if (is.null(seed)) NULL else
        as.integer((as.double(seed) + 104729 * index) %% .Machine$integer.max)
      candidate_fit <- .kmedoid_fit(
        distance, grid[index], nstart, max_iter, tol, fit_seed,
        selection_solver, exact_limit
      )
      evaluated_fits[[index]] <- candidate_fit
      evaluated_objective[index] <- candidate_fit$objective
      evaluated_solver[index] <- candidate_fit$solver
      evaluated_global[index] <- candidate_fit$global_optimum
      candidate_root <- sqrt(effective_n) * candidate_fit$objective
      candidate_group <- sqrt(complexity[index] + candidate_fit$objective)
      candidate_index <- max(candidate_root, candidate_group)
      if (is.finite(candidate_index) && candidate_index <= selected_threshold) {
        selected_index <- index
        break
      }
    }

    last_evaluated <- if (is.na(selected_index)) length(grid) else selected_index
    evaluated <- seq_len(last_evaluated)
    objective <- evaluated_objective[evaluated]
    fit_solver <- evaluated_solver[evaluated]
    global_optimum <- evaluated_global[evaluated]
    evaluated_complexity <- complexity[evaluated]
    root_remainder <- sqrt(effective_n) * objective
    group_remainder <- sqrt(evaluated_complexity + objective)
    inference_index <- pmax(root_remainder, group_remainder)
    inference_feasible <- is.finite(inference_index) &
      inference_index <= inference_threshold
    theoretical_feasible <- is.finite(inference_index) &
      inference_index <= theoretical_threshold

    if (solver == "exact" && any(!global_optimum)) {
      warning(
        "Some evaluated G values exceed kmedoid_exact_limit; automatic G ",
        "selection uses the multi-start heuristic for those candidates. ",
        "The returned kmedoid_global_optimum flag identifies whether the ",
        "selected fit is certified exact.",
        call. = FALSE
      )
    }
    if (is.na(selected_index)) {
      stop(
        "No candidate k-medoids group count satisfies the '", rule,
        "' rule. The smallest evaluated value of max{sqrt(N*T) Q(G), ",
        "sqrt[G*log(K)/N + Q(G)]} is ",
        signif(min(inference_index), 4L), " but must not exceed ",
        signif(selected_threshold, 4L),
        if (rule == "theoretical") " = 1/log(N*T). " else ". ",
        "Enlarge kmedoid_G_grid or use 'elbow'/'n-third' for descriptive ",
        "clustering without the inference-feasibility guarantee."
      )
    }
    return(list(
      selected_G = as.integer(grid[selected_index]), rule = rule,
      path = data.frame(
        G = grid[evaluated], objective = objective,
        elbow_score = NA_real_, complexity = evaluated_complexity,
        root_remainder = root_remainder,
        group_remainder = group_remainder,
        inference_index = inference_index,
        inference_threshold = inference_threshold,
        inference_feasible = inference_feasible,
        theoretical_threshold = theoretical_threshold,
        theoretical_feasible = theoretical_feasible,
        solver = fit_solver, global_optimum = global_optimum
      ),
      fit = evaluated_fits[[selected_index]]
    ))
  }
  fits <- lapply(seq_along(grid), function(index) {
    fit_seed <- if (is.null(seed)) NULL else
      as.integer((as.double(seed) + 104729 * index) %% .Machine$integer.max)
    .kmedoid_fit(distance, grid[index], nstart, max_iter, tol, fit_seed,
                 selection_solver, exact_limit)
  })
  objective <- vapply(fits, `[[`, numeric(1L), "objective")
  fit_solver <- vapply(fits, `[[`, character(1L), "solver")
  global_optimum <- vapply(fits, `[[`, logical(1L), "global_optimum")
  if (solver == "exact" && any(!global_optimum)) {
    warning(
      "Some candidate G values exceed kmedoid_exact_limit; automatic G ",
      "selection uses the multi-start heuristic for those candidates. ",
      "The returned kmedoid_global_optimum flag identifies whether the ",
      "selected fit is certified exact.",
      call. = FALSE
    )
  }
  root_remainder <- sqrt(effective_n) * objective
  group_remainder <- sqrt(complexity + objective)
  inference_index <- pmax(root_remainder, group_remainder)
  inference_feasible <- is.finite(inference_index) &
    inference_index <= inference_threshold
  theoretical_feasible <- is.finite(inference_index) &
    inference_index <= theoretical_threshold
  if (length(grid) <= 2L || diff(range(objective)) <= .Machine$double.eps) {
    selected_index <- 1L
    elbow_score <- rep(0, length(grid))
  } else {
    x <- (grid - min(grid)) / diff(range(grid))
    y <- (objective - min(objective)) / diff(range(objective))
    line <- 1 - x
    elbow_score <- line - y
    selected_index <- which.max(elbow_score)
  }
  list(
    selected_G = as.integer(grid[selected_index]), rule = "elbow",
    path = data.frame(G = grid, objective = objective,
                      elbow_score = elbow_score,
                      complexity = complexity,
                      root_remainder = root_remainder,
                      group_remainder = group_remainder,
                      inference_index = inference_index,
                      inference_threshold = inference_threshold,
                      inference_feasible = inference_feasible,
                      theoretical_threshold = theoretical_threshold,
                      theoretical_feasible = theoretical_feasible,
                      solver = fit_solver,
                      global_optimum = global_optimum),
    fit = fits[[selected_index]]
  )
}

.kmedoid_validate_options <- function(
    nstart, max_iter, tol, seed, standardize, G_rule, G_grid,
    min_group_size, solver, exact_limit) {
  for (item in list(nstart = nstart, max_iter = max_iter,
                    min_group_size = min_group_size,
                    exact_limit = exact_limit)) {
    if (length(item) != 1L || !is.finite(item) || item < 1L ||
        item != floor(item)) stop("K-medoids integer options must be positive integers.")
  }
  if (length(tol) != 1L || !is.finite(tol) || tol <= 0) {
    stop("kmedoid_tol must be one positive finite number.")
  }
  if (!is.null(seed) && (length(seed) != 1L || !is.finite(seed) || seed < 0 ||
                         seed > .Machine$integer.max || seed != floor(seed))) {
    stop("kmedoid_seed must be NULL or one nonnegative integer.")
  }
  if (length(standardize) != 1L || is.na(standardize) || !is.logical(standardize)) {
    stop("kmedoid_standardize must be TRUE or FALSE.")
  }
  G_rule <- match.arg(
    G_rule[1L], c("inference", "theoretical", "elbow", "n-third")
  )
  solver <- match.arg(solver[1L], c("heuristic", "auto", "exact"))
  if (!is.null(G_grid)) {
    if (!is.numeric(G_grid) || !length(G_grid) || any(!is.finite(G_grid)) ||
        any(G_grid < 1L) || any(G_grid != floor(G_grid))) {
      stop("kmedoid_G_grid must be NULL or positive integers.")
    }
    G_grid <- sort(unique(as.integer(G_grid)))
  }
  list(
    nstart = as.integer(nstart), max_iter = as.integer(max_iter), tol = tol,
    seed = seed, standardize = standardize, G_rule = G_rule, G_grid = G_grid,
    min_group_size = as.integer(min_group_size), solver = solver,
    exact_limit = as.integer(exact_limit)
  )
}

.kmedoid_fit_profile <- function(profile, unit_cluster, options) {
  distance <- .kmedoid_projection_distance(profile$data)
  grid <- .kmedoid_candidate_grid(
    nrow(profile$data), options$min_group_size, options$G_grid
  )
  if (is.null(unit_cluster)) {
    selection <- .kmedoid_select_G(
      distance, grid, options$G_rule, options$nstart, options$max_iter,
      options$tol, options$seed, options$solver, options$exact_limit,
      profile$T, profile$K
    )
    G <- selection$selected_G
    fit <- selection$fit
    if (is.null(fit)) {
      fit <- .kmedoid_fit(
        distance, G, options$nstart, options$max_iter, options$tol,
        options$seed, options$solver, options$exact_limit
      )
    }
  } else {
    max_feasible <- floor(nrow(profile$data) / options$min_group_size)
    if (length(unit_cluster) != 1L || !is.finite(unit_cluster) ||
        unit_cluster < 1L || unit_cluster != floor(unit_cluster) ||
        unit_cluster > max_feasible) {
      stop("unit_cluster is outside the feasible balanced k-medoids range.")
    }
    G <- as.integer(unit_cluster)
    fit <- .kmedoid_fit(
      distance, G, options$nstart, options$max_iter, options$tol,
      options$seed, options$solver, options$exact_limit
    )
    selection <- list(selected_G = G, rule = "fixed", path = NULL, fit = fit)
  }
  list(fit = fit, selection = selection, distance = distance)
}

.kmedoid_project_panel <- function(y, X, group, G) {
  N <- nrow(y)
  TT <- ncol(y)
  K <- dim(X)[3L]
  Du <- sapply(seq_len(G), function(g) as.numeric(group == g))
  if (G == 1L) Du <- matrix(Du, ncol = 1L)
  Mu <- diag(N) - Du %*% diag(1 / colSums(Du), G) %*% t(Du)
  Z <- array(NA_real_, c(N, TT, K + 1L))
  Z[, , 1L] <- y
  for (k in seq_len(K)) Z[, , k + 1L] <- X[, , k]
  projected <- array(NA_real_, dim(Z))
  for (k in seq_len(K + 1L)) projected[, , k] <- Mu %*% Z[, , k]
  projected
}

.HP_estimate_kmedoid_full <- function(
    data, y_col, covariate_cols, id_col, time_col, unit_cluster,
    kmedoid_nstart, kmedoid_max_iter, kmedoid_tol, kmedoid_seed,
    kmedoid_standardize, kmedoid_G_rule, kmedoid_G_grid,
    kmedoid_min_group_size, kmedoid_solver, kmedoid_exact_limit) {
  required <- unique(c(id_col, time_col, y_col, covariate_cols))
  missing_columns <- setdiff(required, names(data))
  if (length(missing_columns)) {
    stop("K-medoids data are missing columns: ",
         paste(missing_columns, collapse = ", "), ".")
  }
  if (length(y_col) != 1L || length(covariate_cols) < 2L) {
    stop("K-medoids requires one outcome and treatment followed by at least one control.")
  }
  ids <- sort(unique(data[[id_col]]))
  times <- sort(unique(data[[time_col]]))
  N <- length(ids)
  TT <- length(times)
  K <- length(covariate_cols)
  if (N < 3L || TT < 2L || nrow(data) != N * TT ||
      anyDuplicated(data[c(id_col, time_col)])) {
    stop("K-medoids requires a complete balanced panel with N >= 3 and T >= 2.")
  }
  options <- .kmedoid_validate_options(
    kmedoid_nstart, kmedoid_max_iter, kmedoid_tol, kmedoid_seed,
    kmedoid_standardize, kmedoid_G_rule, kmedoid_G_grid,
    kmedoid_min_group_size, kmedoid_solver, kmedoid_exact_limit
  )

  y <- matrix(NA_real_, N, TT)
  X <- array(NA_real_, c(N, TT, K))
  ii <- match(data[[id_col]], ids)
  tt <- match(data[[time_col]], times)
  for (row in seq_len(nrow(data))) {
    y[ii[row], tt[row]] <- data[[y_col]][row]
    for (k in seq_len(K)) X[ii[row], tt[row], k] <- data[[covariate_cols[k]]][row]
  }
  profile <- .kmedoid_profile(y, X, options$standardize)
  fitted <- .kmedoid_fit_profile(profile, unit_cluster, options)
  medoid_fit <- fitted$fit
  G <- medoid_fit$G
  group <- medoid_fit$group
  projected <- .kmedoid_project_panel(y, X, group, G)
  tY <- as.vector(projected[, , 1L])
  tX <- do.call(cbind, lapply(seq_len(K), function(k) {
    as.vector(projected[, , k + 1L])
  }))
  colnames(tX) <- covariate_cols
  treatment <- tX[, 1L]
  controls <- tX[, -1L, drop = FALSE]

  fit <- hdm::rlassoEffect(x = controls, y = tY, d = treatment,
                           method = "double selection")
  trans <- data.frame(y = tY, D = treatment, controls, check.names = FALSE)
  Ytilde <- hdm::rlasso(y ~ . - D - 1, data = trans)$residuals
  Dtilde <- hdm::rlasso(D ~ . - 1, data = trans[, -1L, drop = FALSE])$residuals
  data_res <- data.frame(
    id = rep(ids, TT), time = rep(times, each = N),
    Ytilde = Ytilde, Dtilde = Dtilde
  )
  post_plm <- plm::plm(Ytilde ~ -1 + Dtilde, data = data_res,
                       model = "pooling", index = c("id", "time"))
  robust_se <- sqrt(diag(plm::vcovHC(post_plm, type = "HC0", method = "arellano")))
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
    cluster_method = "kmedoid", cross_fitted = FALSE, cross_type = "none",
    G_unit = G, unit_group = group,
    fit_summary = summary(fit), post_plm_summary = summary(post_plm),
    estimate_corrected = summary_table, summary_table = summary_table,
    kmedoid_profile_dimension = profile$M,
    kmedoid_objective = medoid_fit$objective,
    kmedoid_group_sizes = medoid_fit$group_sizes,
    kmedoid_medoid_index = medoid_fit$medoids,
    kmedoid_medoid_unit_id = ids[medoid_fit$medoids],
    kmedoid_solver = medoid_fit$solver,
    kmedoid_global_optimum = medoid_fit$global_optimum,
    kmedoid_converged = medoid_fit$converged,
    kmedoid_iterations = medoid_fit$iterations,
    kmedoid_selection_rule = fitted$selection$rule,
    kmedoid_selection_path = fitted$selection$path,
    kmedoid_distance = "quadratic-projection",
    kmedoid_profile_center = profile$center,
    kmedoid_profile_scale = profile$scale
  )
}
