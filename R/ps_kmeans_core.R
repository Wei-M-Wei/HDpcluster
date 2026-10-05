# Shared profile-sieve (PS) k-means primitives for full-sample and cross-fitted estimators.
#
# Install at least one assignment solver:
#   install.packages("lpSolve")  # compact exact solver
#   install.packages("clue")     # faster exact solver for larger N
#
# Long-data convention:
#   one row per (unit, time), one column per component of z_it.
#
# For evaluation fold d, the profile matrix uses only T_{-d}. Its ith row is
#
#   a_i^d = (z_itk : t in T_{-d}, k in C)',
#
# so its dimensions are N x M_d, where M_d = T_{-d} * p_c.

make_ps_profile <- function(data,
                            id,
                            time,
                            clustering_vars,
                            auxiliary_times,
                            evaluation_times = NULL,
                            standardize = TRUE,
                            scale_bounds = c(0.25, 4)) {
  required <- c(id, time, clustering_vars)
  missing_columns <- setdiff(required, names(data))
  if (length(missing_columns) > 0L) {
    stop("Missing columns: ", paste(missing_columns, collapse = ", "))
  }
  if (length(clustering_vars) < 1L) {
    stop("clustering_vars must contain at least one component of z_it.")
  }
  if (length(auxiliary_times) < 1L) {
    stop("auxiliary_times must be nonempty.")
  }
  if (!is.null(evaluation_times) &&
      length(intersect(auxiliary_times, evaluation_times)) > 0L) {
    stop("Auxiliary and evaluation dates overlap. This would violate time-fold cross-fitting.")
  }
  if (length(scale_bounds) != 2L ||
      any(!is.finite(scale_bounds)) ||
      scale_bounds[1L] <= 0 ||
      scale_bounds[1L] > scale_bounds[2L]) {
    stop("scale_bounds must be two positive ordered numbers.")
  }

  ids <- unique(data[[id]])
  auxiliary_times <- unique(auxiliary_times)
  keep <- data[[time]] %in% auxiliary_times
  auxiliary_data <- data[keep, required, drop = FALSE]

  unit_index <- match(auxiliary_data[[id]], ids)
  time_index <- match(auxiliary_data[[time]], auxiliary_times)
  unit_time_index <- cbind(unit_index, time_index)
  if (anyDuplicated(data.frame(unit_index, time_index))) {
    stop("There is more than one row for at least one (unit, auxiliary-time) pair.")
  }

  N <- length(ids)
  T_minus_d <- length(auxiliary_times)
  p_c <- length(clustering_vars)
  M_d <- T_minus_d * p_c
  A_raw <- matrix(NA_real_, nrow = N, ncol = M_d)

  for (row in seq_len(nrow(auxiliary_data))) {
    i <- unit_index[row]
    tt <- time_index[row]
    columns <- (tt - 1L) * p_c + seq_len(p_c)
    A_raw[i, columns] <- as.numeric(auxiliary_data[row, clustering_vars, drop = TRUE])
  }

  coordinate_names <- unlist(
    lapply(
      auxiliary_times,
      function(tt) paste0("t=", tt, "::", clustering_vars)
    ),
    use.names = FALSE
  )
  rownames(A_raw) <- as.character(ids)
  colnames(A_raw) <- coordinate_names

  incomplete <- which(!is.finite(A_raw), arr.ind = TRUE)
  if (nrow(incomplete) > 0L) {
    first_bad <- incomplete[1L, ]
    stop(
      "The auxiliary panel is incomplete or non-finite. First problem: unit ",
      ids[first_bad[1L]], ", coordinate ", coordinate_names[first_bad[2L]], "."
    )
  }

  center <- rep(0, M_d)
  scale <- rep(1, M_d)
  A <- A_raw

  if (standardize) {
    center <- apply(A_raw, 2L, stats::median)
    raw_scale <- apply(A_raw, 2L, stats::mad)

    bad_scale <- !is.finite(raw_scale) | raw_scale <= sqrt(.Machine$double.eps)
    if (any(bad_scale)) {
      raw_scale[bad_scale] <- apply(A_raw[, bad_scale, drop = FALSE], 2L, stats::sd)
    }

    positive_scale <- raw_scale[is.finite(raw_scale) &
                                  raw_scale > sqrt(.Machine$double.eps)]
    reference_scale <- if (length(positive_scale)) stats::median(positive_scale) else 1
    raw_scale[!is.finite(raw_scale) |
                raw_scale <= sqrt(.Machine$double.eps)] <- reference_scale

    relative_scale <- raw_scale / reference_scale
    relative_scale <- pmin(
      pmax(relative_scale, scale_bounds[1L]),
      scale_bounds[2L]
    )
    scale <- reference_scale * relative_scale
    A <- sweep(sweep(A_raw, 2L, center, "-"), 2L, scale, "/")
  }

  structure(
    list(
      matrix = A,
      raw_matrix = A_raw,
      unit_id = ids,
      coordinates = coordinate_names,
      auxiliary_times = auxiliary_times,
      clustering_vars = clustering_vars,
      center = center,
      scale = scale,
      N = N,
      T_minus_d = T_minus_d,
      p_c = p_c,
      M_d = M_d
    ),
    class = "ps_profile"
  )
}


.ps_squared_distances <- function(A, centers, A_norm2 = NULL) {
  if (is.null(A_norm2)) A_norm2 <- rowSums(A^2)
  distances <- outer(A_norm2, rowSums(centers^2), "+") -
    2 * tcrossprod(A, centers)
  pmax(distances / ncol(A), 0)
}


.ps_kmeans_plus_plus <- function(A, G, A_norm2 = NULL) {
  N <- nrow(A)
  if (is.null(A_norm2)) A_norm2 <- rowSums(A^2)
  chosen <- integer(G)
  chosen[1L] <- sample.int(N, 1L)
  distance_to <- function(index) {
    pmax(
      A_norm2 + A_norm2[index] -
        2 * drop(A %*% A[index, ]),
      0
    )
  }
  closest <- distance_to(chosen[1L])

  if (G >= 2L) {
    for (g in 2:G) {
      available <- setdiff(seq_len(N), chosen[seq_len(g - 1L)])
      probabilities <- closest
      probabilities[-available] <- 0
      if (sum(probabilities) <= 0) {
        chosen[g] <- sample(available, 1L)
      } else {
        chosen[g] <- sample.int(N, 1L, prob = probabilities)
      }
      new_distance <- distance_to(chosen[g])
      closest <- pmin(closest, new_distance)
    }
  }
  A[chosen, , drop = FALSE]
}


.ps_balanced_capacities <- function(N, G) {
  capacities <- rep(N %/% G, G)
  remainder <- N %% G
  if (remainder > 0L) {
    capacities[seq_len(remainder)] <- capacities[seq_len(remainder)] + 1L
  }
  capacities
}


.ps_balanced_assignment <- function(distance_matrix, capacities) {
  distance_matrix <- as.matrix(distance_matrix)
  storage.mode(distance_matrix) <- "double"
  G <- ncol(distance_matrix)
  N <- nrow(distance_matrix)
  if (length(capacities) != G || sum(capacities) != N) {
    stop("The cluster capacities do not match the distance matrix.")
  }
  if (any(capacities < 1L) || any(!is.finite(distance_matrix))) {
    stop("The balanced assignment requires positive capacities and finite distances.")
  }

  # The projection features can produce very large squared distances when the
  # original variables are not standardized.  lpSolve is sensitive to this
  # scale even though the transportation problem is feasible.  Subtracting a
  # row constant and dividing all costs by one positive constant leave the
  # minimizing assignment unchanged because every row must be assigned once.
  row_minimum <- distance_matrix[, 1L]
  if (G > 1L) {
    for (g in 2:G) row_minimum <- pmin(row_minimum, distance_matrix[, g])
  }
  cost <- sweep(distance_matrix, 1L, row_minimum, "-")
  cost_scale <- max(cost)
  if (!is.finite(cost_scale)) {
    stop("The balanced assignment costs overflowed; rescale the clustering variables.")
  }
  if (cost_scale <= sqrt(.Machine$double.eps)) {
    return(rep(seq_len(G), capacities))
  }
  cost <- cost / cost_scale

  valid_assignment <- function(group) {
    length(group) == N && all(group %in% seq_len(G)) &&
      identical(as.integer(tabulate(group, nbins = G)), as.integer(capacities))
  }

  solve_with_clue <- function() {
    group_slots <- rep(seq_len(G), capacities)
    tryCatch({
      assigned_slot <- as.integer(clue::solve_LSAP(
        cost[, group_slots, drop = FALSE], maximum = FALSE
      ))
      group_slots[assigned_slot]
    }, error = function(error) NULL)
  }

  # Replicating group slots turns the transportation problem into an exact
  # linear-sum assignment. Benchmarking over the simulation range shows that
  # its compiled solver is faster than the compact transportation solver even
  # for small panels, so use it first whenever it is installed.
  clue_available <- requireNamespace("clue", quietly = TRUE)
  if (clue_available) {
    group <- solve_with_clue()
    if (!is.null(group) && valid_assignment(group)) return(group)
  }

  # This N x G transportation problem has an integer solution because its
  # constraint matrix is totally unimodular. It is also the lower-memory exact
  # fallback when clue is unavailable or its N x N expansion fails.
  if (requireNamespace("lpSolve", quietly = TRUE)) {
    solution <- tryCatch(
      lpSolve::lp.transport(
        cost.mat = cost,
        direction = "min",
        row.signs = rep("=", N),
        row.rhs = rep(1, N),
        col.signs = rep("=", G),
        col.rhs = capacities
      ),
      error = function(error) NULL
    )
    if (!is.null(solution) && solution$status == 0L) {
      group <- max.col(solution$solution, ties.method = "first")
      if (valid_assignment(group)) return(group)
    }
  }

  # Last-resort feasible assignment.  This is used only if both exact solvers
  # are unavailable or fail.  It preserves the balance restriction and keeps
  # a single candidate/start from aborting the complete G-selection path.
  warning(
    "Exact balanced assignment failed; using a greedy balanced fallback.",
    call. = FALSE
  )
  group <- integer(N)
  available <- as.integer(capacities)
  for (index in order(cost)) {
    i <- ((index - 1L) %% N) + 1L
    g <- ((index - 1L) %/% N) + 1L
    if (group[i] == 0L && available[g] > 0L) {
      group[i] <- g
      available[g] <- available[g] - 1L
    }
    if (all(group > 0L)) break
  }
  if (!valid_assignment(group)) {
    stop(
      "The balanced assignment failed after normalized exact and greedy attempts."
    )
  }
  group
}


.ps_update_centers <- function(A, group, G) {
  group_sizes <- tabulate(group, nbins = G)
  group_sums <- rowsum(A, group, reorder = TRUE)
  unname(sweep(group_sums, 1L, group_sizes, "/"))
}


ps_kmeans <- function(profile,
                      G = ceiling(profile$N^(1 / 3)),
                      nstart = 30L,
                      max_iter = 200L,
                      tol = 1e-8,
                      seed = 12345L) {
  A <- if (inherits(profile, "ps_profile")) profile$matrix else as.matrix(profile)
  storage.mode(A) <- "double"
  N <- nrow(A)

  if (any(!is.finite(A))) stop("The profile matrix contains non-finite values.")
  if (G < 1L || G > N) stop("G must be between 1 and N.")
  if (nstart < 1L || max_iter < 1L || tol <= 0) {
    stop("nstart, max_iter, and tol must be positive.")
  }

  G <- as.integer(G)
  nstart <- as.integer(nstart)
  capacities <- .ps_balanced_capacities(N, G)
  A_norm2 <- rowSums(A^2)
  set.seed(seed)
  fits <- vector("list", nstart)

  for (start in seq_len(nstart)) {
    centers <- .ps_kmeans_plus_plus(A, G, A_norm2 = A_norm2)
    previous_group <- rep(NA_integer_, N)
    previous_objective <- Inf
    converged <- FALSE
    distances <- .ps_squared_distances(A, centers, A_norm2 = A_norm2)

    for (iteration in seq_len(max_iter)) {
      group <- .ps_balanced_assignment(distances, capacities)
      centers <- .ps_update_centers(A, group, G)
      updated_distances <- .ps_squared_distances(
        A, centers, A_norm2 = A_norm2
      )
      objective <- mean(updated_distances[cbind(seq_len(N), group)])

      same_assignment <- identical(group, previous_group)
      relative_change <- if (is.finite(previous_objective)) {
        abs(previous_objective - objective) /
          max(abs(previous_objective), .Machine$double.eps)
      } else {
        Inf
      }

      if (same_assignment || relative_change < tol) {
        converged <- TRUE
        break
      }
      previous_group <- group
      previous_objective <- objective
      distances <- updated_distances
    }

    fits[[start]] <- list(
      group = group,
      centers = centers,
      objective = objective,
      iterations = iteration,
      converged = converged
    )
  }

  objectives <- vapply(fits, `[[`, numeric(1L), "objective")
  best_start <- which.min(objectives)
  best <- fits[[best_start]]
  unit_id <- if (inherits(profile, "ps_profile")) {
    profile$unit_id
  } else if (!is.null(rownames(A))) {
    rownames(A)
  } else {
    seq_len(N)
  }

  names(best$group) <- as.character(unit_id)
  best$G <- G
  best$capacities <- capacities
  best$group_sizes <- tabulate(best$group, nbins = G)
  best$best_start <- best_start
  best$objectives_by_start <- objectives
  best$unit_id <- unit_id
  best$profile <- profile
  class(best) <- "ps_kmeans"
  best
}


