#' Estimator Function with Optional Cross-Fitting
#'@title Inference in high-dimensional data after discretizing unobserved heterogeneity
#'
#' @description R package 'HDpcluster' is dedicated to do inference for high-dimensional linear panel data model with unkown functions of fixed effects.
#'
#' @param y outcome variable
#' @param D treatment variable
#' @param X control variables
#' @param T panel size of the time period
#' @param groups_init number of the unit clusters
#' @param index index name
#' @param data data which contains the correct index of the outcome variable or the treatment variable
#' @param cluster_type different cluster means for the data, 'unit kmeans' is default, allow for 'unit pesudo'
#' @param cluster_method clustering method. Options include 'kmeans', 'hierarchical', 'pesudo', 'ps-kmeans', 'HAC new', 'kmedoid', and 'projection-kmeans'. In HP_estimate_cf, 'kmedoid' and 'projection-kmeans' use two time folds.
#' @param ps_nstart number of random k-means++ starts for each PS training fold.
#' @param ps_max_iter maximum balanced Lloyd iterations for each PS start.
#' @param ps_tol relative objective tolerance for PS convergence.
#' @param ps_seed optional random seed for reproducible PS initialization. The caller's random-number state is restored after fitting.
#' @param ps_standardize whether each profile coordinate is robustly standardized using training-fold statistics.
#' @param ps_G_rule group-count rule used when unit_cluster is NULL. "one-se" selects the smallest candidate whose cross-validated profile loss is within one standard error of the minimum; "minimum" minimizes that loss; "n-third" retains the former ceiling(N_train^(1/3)) rule.
#' @param ps_G_grid optional positive-integer candidate group counts. If NULL, a computationally sparse grid is constructed up to the largest count compatible with ps_min_group_size and the inner validation folds.
#' @param ps_cv_folds number of inner auxiliary-unit folds used to select the PS group count.
#' @param ps_cv_repeats number of repeated inner fold partitions.
#' @param ps_select_nstart number of balanced k-means++ starts used for each candidate in inner validation.
#' @param ps_min_group_size minimum training- and held-out-fold group size used to cap the candidate grid.
#' @param kmedoid_nstart number of multi-start balanced k-medoids fits when the exact search is not feasible.
#' @param kmedoid_max_iter maximum medoid-update iterations per start.
#' @param kmedoid_tol relative objective tolerance for k-medoids convergence.
#' @param kmedoid_seed optional nonnegative integer for reproducible starts.
#' @param kmedoid_standardize whether clustering coordinates are robustly standardized within the clustering block.
#' @param kmedoid_G_rule group-count rule when unit_cluster is NULL. 'inference' (the default) uses cutoff one for the combined inference diagnostic; 'theoretical' uses the vanishing cutoff 1/log(N times T); 'elbow' uses the objective elbow; and 'n-third' selects the candidate closest to ceiling((N times T)^(1/3)).
#' @param kmedoid_G_grid optional positive-integer candidate group counts.
#' @param kmedoid_min_group_size minimum balanced group size used to cap the candidate grid.
#' @param kmedoid_solver 'heuristic' (the default), 'auto', or 'exact'. Auto uses exact enumeration when it is below kmedoid_exact_limit. With data-driven G, 'exact' fits feasible candidates exactly and uses a reported heuristic fallback for candidates above the limit; with fixed unit_cluster, 'exact' remains strict.
#' @param kmedoid_exact_limit maximum number of medoid/capacity configurations considered by the exact solver.
#' @param pk_nstart number of balanced K-means starts for projection K-means.
#' @param pk_max_iter maximum balanced Lloyd iterations for each projection K-means start.
#' @param pk_select_nstart number of starts used to screen candidate group counts under the projection K-means inference rule. The selected count is refitted using pk_nstart.
#' @param pk_select_max_iter maximum Lloyd iterations used while screening candidate group counts under the projection K-means inference rule. The selected count is refitted using pk_max_iter.
#' @param pk_tol relative objective tolerance for projection K-means.
#' @param pk_seed optional nonnegative integer for reproducible starts.
#' @param pk_standardize whether the time-covariate coordinates are robustly standardized before the Gram matrix is formed.
#' @param pk_G_rule group-count rule for projection K-means when unit_cluster is NULL. 'noise-floor' (the default) is the direct analogue of the first K-means rule. 'inference' selects the smallest candidate whose objective is no more than 1 + 1/N times the smallest objective under the applicable maximum group count. 'elbow' uses the objective elbow, and 'n-third' selects the candidate closest to ceiling((N times T)^(1/3)).
#' @param pk_G_grid optional positive-integer candidate group counts.
#' @param pk_max_G optional positive-integer maximum candidate group count for projection K-means. If NULL, the inference rule uses the theoretical cap determined by pk_log_N; the other rules use the largest count allowed by pk_min_group_size. A supplied value replaces the theoretical cap and is itself limited only by group-size feasibility.
#' @param pk_log_N logarithmic factor used in the default projection K-means inference cap. Use "log" (the default) for log(N), "loglog" for log(log(N)), or supply one positive finite number to use directly. Writing this factor as L_N, the inference rule uses `theoretical_max_G = floor(N / (log(K) * L_N))` and reports the matching complexity threshold `1 / L_N`. A supplied pk_max_G overrides this cap, and a supplied unit_cluster bypasses group-count selection.
#' @param pk_min_group_size minimum balanced group size used to cap the candidate grid.
#' @param pk_noise_multiplier positive multiplier on the estimated projection-feature noise floor.
#' @param hac_G_rule group-count rule for 'HAC new' when unit_cluster is NULL: 'one-se' or 'minimum'.
#' @param hac_G_grid optional positive-integer candidate counts for 'HAC new'.
#' @param hac_cv_folds number of inner auxiliary-unit folds used to select the HAC group count.
#' @param hac_cv_repeats number of repeated inner HAC validation partitions.
#' @param hac_min_group_size minimum held-out group size used to cap the HAC candidate grid.
#' @param hac_reference_fraction fraction of nonvalidation auxiliary units used as pseudo-distance references.
#' @param hac_seed optional nonnegative integer controlling reference and validation splits.
#' @param hac_standardize whether each HAC profile coordinate is robustly standardized using outer-training statistics.
#' @param pesudo_type if cluster_type = 'unit pesudo', choice of the pesudo type
#' @param link if cluster_type = 'unit pesudo', choice of different links of pesudo type
#' @param optimal_index if cluster_type = 'unit pesudo', different ways to compute the optimal number of clusters
#' @param all_targets if TRUE, return DML estimates for every variable in covariate_cols by treating each variable in turn as the target regressor. The default FALSE preserves the original first-target behavior.
#'
#' @returns A list of fitted results is returned.
#' Within this outputted list, the following elements can be found:
#'     \item{res}{regression model.}
#'     \item{G}{number of unit clusters.}
#'     \item{estimate_correct}{corrected standard error, t value, and p value.}
#'     \item{summary_table}{summary of the model together with corrected estimates in the 'Coefficients'.}
#'
#' @import hdm
#' @import plm
#' @useDynLib HDpcluster
#' @export
HP_estimate <- function(data, y_col = NULL, covariate_cols = NULL,
                        id_col = "id", time_col = "time", cluster_type = c('two way'),
                        index = c("id", "time"), method_auto = "dynamicTreeCut", cluster_method = c('kmeans', 'hierarchical', 'pesudo', 'pseudo', 'ps-kmeans', 'HAC new', 'kmedoid', 'projection-kmeans'), link = 'average', pre_cluster = FALSE, Du_pre, Dv_pre, unit_cluster = NULL, time_cluster = NULL, gamma = 1, cc = 0, dim_moment = 1,  dist_type = "proj",
                        ps_nstart = 30L, ps_max_iter = 100L, ps_tol = 1e-8,
                        ps_seed = NULL, ps_standardize = TRUE,
                        ps_G_rule = c("one-se", "minimum", "n-third"),
                        ps_G_grid = NULL, ps_cv_folds = 2L,
                        ps_cv_repeats = 3L, ps_select_nstart = 3L,
                        ps_min_group_size = 2L,
                        kmedoid_nstart = 30L, kmedoid_max_iter = 100L,
                        kmedoid_tol = 1e-8, kmedoid_seed = NULL,
                        kmedoid_standardize = TRUE,
                        kmedoid_G_rule = c("inference", "theoretical", "elbow", "n-third"),
                        kmedoid_G_grid = NULL,
                        kmedoid_min_group_size = 2L,
                        kmedoid_solver = c("heuristic", "auto", "exact"),
                        kmedoid_exact_limit = 50000L,
                        pk_nstart = 30L, pk_max_iter = 100L,
                        pk_select_nstart = 3L, pk_select_max_iter = 25L,
                        pk_tol = 1e-8, pk_seed = NULL,
                        pk_standardize = FALSE,
                        pk_G_rule = c("noise-floor", "inference", "elbow", "n-third"),
                        pk_G_grid = NULL, pk_max_G = NULL,
                        pk_log_N = "log",
                        pk_min_group_size = 2L,
                        pk_noise_multiplier = 1,
                        hac_G_rule = c("one-se", "minimum"),
                        hac_G_grid = NULL, hac_cv_folds = 2L,
                        hac_cv_repeats = 3L, hac_min_group_size = 2L,
                        hac_reference_fraction = 0.5, hac_seed = NULL,
                        hac_standardize = TRUE,
                        all_targets = FALSE) {


  if (is.null(covariate_cols)) {
    covariate_cols <- setdiff(colnames(data), c(id_col, time_col, y_col))
  }

  if (length(all_targets) != 1L || is.na(all_targets) || !is.logical(all_targets)) {
    stop("all_targets must be TRUE or FALSE.")
  }

  if (isTRUE(all_targets)) {
    if (!length(covariate_cols)) {
      stop("all_targets = TRUE requires at least one covariate.")
    }

    # Reuse the existing, tested single-target code. For target k, move the
    # kth regressor to the first position and treat the remaining regressors
    # as high-dimensional controls. Resetting the RNG before each call keeps
    # stochastic clustering identical across target regressions whenever the
    # clustering criterion itself is invariant to the covariate ordering
    # (in particular, projection K-means).
    caller_env <- parent.frame()
    base_call <- match.call()
    base_call$all_targets <- FALSE

    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      invisible(runif(1L))
    }
    rng_start <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    rng_after_first <- NULL

    target_fits <- vector("list", length(covariate_cols))
    for (k in seq_along(covariate_cols)) {
      assign(".Random.seed", rng_start, envir = .GlobalEnv)
      call_k <- base_call
      call_k$covariate_cols <- c(covariate_cols[k], covariate_cols[-k])
      target_fits[[k]] <- eval(call_k, envir = caller_env)
      if (k == 1L) {
        rng_after_first <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
      }
    }
    if (!is.null(rng_after_first)) {
      assign(".Random.seed", rng_after_first, envir = .GlobalEnv)
    }

    summary_all <- do.call(rbind, lapply(seq_along(covariate_cols), function(k) {
      tab <- target_fits[[k]]$summary_table
      if (is.null(tab) || nrow(tab) < 1L) {
        stop("No summary-table row was returned for target ", covariate_cols[k], ".")
      }
      data.frame(
        Variable = covariate_cols[k],
        tab[1L, , drop = FALSE],
        row.names = NULL,
        check.names = FALSE
      )
    }))
    rownames(summary_all) <- covariate_cols

    result <- target_fits[[1L]]
    result$summary_table <- summary_all
    result$estimate_corrected <- summary_all
    result$all_targets <- TRUE
    result$target_names <- covariate_cols
    result$G_unit_all <- setNames(
      vapply(target_fits, function(z) as.numeric(z$G_unit)[1L], numeric(1L)),
      covariate_cols
    )
    return(result)
  }

  ids <- sort(unique(data[[id_col]]))
  times <- sort(unique(data[[time_col]]))
  N <- length(ids)
  T <- length(times)
  K <- length(covariate_cols)
  cluster_method <- match.arg(cluster_method)
  if (cluster_method == "pseudo") {
    cluster_method <- "pesudo"
  }

  if (cluster_method == "ps-kmeans") {
    if (pre_cluster) {
      stop("Full-sample PS k-means requires pre_cluster = FALSE.")
    }
    return(.HP_estimate_ps_kmeans_full(
      data = data,
      y_col = y_col,
      covariate_cols = covariate_cols,
      id_col = id_col,
      time_col = time_col,
      unit_cluster = unit_cluster,
      ps_nstart = ps_nstart,
      ps_max_iter = ps_max_iter,
      ps_tol = ps_tol,
      ps_seed = ps_seed,
      ps_standardize = ps_standardize,
      ps_G_rule = ps_G_rule,
      ps_G_grid = ps_G_grid,
      ps_cv_folds = ps_cv_folds,
      ps_cv_repeats = ps_cv_repeats,
      ps_select_nstart = ps_select_nstart,
      ps_min_group_size = ps_min_group_size
    ))
  }

  if (cluster_method == "HAC new") {
    if (pre_cluster) {
      stop("Full-sample HAC new requires pre_cluster = FALSE.")
    }
    return(.HP_estimate_hac_new_full(
      data = data,
      y_col = y_col,
      covariate_cols = covariate_cols,
      id_col = id_col,
      time_col = time_col,
      link = link,
      unit_cluster = unit_cluster,
      hac_G_rule = hac_G_rule,
      hac_G_grid = hac_G_grid,
      hac_cv_folds = hac_cv_folds,
      hac_cv_repeats = hac_cv_repeats,
      hac_min_group_size = hac_min_group_size,
      hac_reference_fraction = hac_reference_fraction,
      hac_seed = hac_seed,
      hac_standardize = hac_standardize
    ))
  }

  if (cluster_method == "kmedoid") {
    if (pre_cluster) {
      stop("Full-sample k-medoids requires pre_cluster = FALSE.")
    }
    return(.HP_estimate_kmedoid_full(
      data = data,
      y_col = y_col,
      covariate_cols = covariate_cols,
      id_col = id_col,
      time_col = time_col,
      unit_cluster = unit_cluster,
      kmedoid_nstart = kmedoid_nstart,
      kmedoid_max_iter = kmedoid_max_iter,
      kmedoid_tol = kmedoid_tol,
      kmedoid_seed = kmedoid_seed,
      kmedoid_standardize = kmedoid_standardize,
      kmedoid_G_rule = kmedoid_G_rule,
      kmedoid_G_grid = kmedoid_G_grid,
      kmedoid_min_group_size = kmedoid_min_group_size,
      kmedoid_solver = kmedoid_solver,
      kmedoid_exact_limit = kmedoid_exact_limit
    ))
  }

  if (cluster_method == "projection-kmeans") {
    if (pre_cluster) {
      stop("Full-sample projection K-means requires pre_cluster = FALSE.")
    }
    return(.HP_estimate_projection_kmeans_full(
      data = data,
      y_col = y_col,
      covariate_cols = covariate_cols,
      id_col = id_col,
      time_col = time_col,
      unit_cluster = unit_cluster,
      pk_nstart = pk_nstart,
      pk_max_iter = pk_max_iter,
      pk_select_nstart = pk_select_nstart,
      pk_select_max_iter = pk_select_max_iter,
      pk_tol = pk_tol,
      pk_seed = pk_seed,
      pk_standardize = pk_standardize,
      pk_G_rule = pk_G_rule,
      pk_G_grid = pk_G_grid,
      pk_max_G = pk_max_G,
      pk_log_N = pk_log_N,
      pk_min_group_size = pk_min_group_size,
      pk_noise_multiplier = pk_noise_multiplier
    ))
  }

  # Initialize y and X
  y <- matrix(NA_real_, nrow = N, ncol = T)
  X <- array(NA_real_, dim = c(N, T, K))

  for (i in seq_len(nrow(data))) {
    id_idx <- which(ids == data[[id_col]][i])
    time_idx <- which(times == data[[time_col]][i])
    y[id_idx, time_idx] <- data[[y_col]][i]
    for (k in seq_along(covariate_cols)) {
      X[id_idx, time_idx, k] <- data[[covariate_cols[k]]][i]
    }
  }

  # Cluster (output indicators assumed one-hot encoded)
  X_norm <- array(NA_real_, dim = dim(X))  # Same shape: (N, T, K)
  for (k in 1:dim(X)[3]) {
    X_slice <- X[,,k]
    X_norm[,,k] <- scale(X_slice)
  }
  y_norm = scale(y)

  if (cluster_method == 'pesudo'){
    if (N < 2) {
      stop("cluster_method = 'pesudo' requires at least two units.")
    }
    if (K < 1) {
      stop("cluster_method = 'pesudo' requires at least one covariate.")
    }

    # Stack y and X into Z: N x T x (K+1)
    Z <- array(NA_real_, dim = c(N, T, K + 1))
    Z[,,1] <- y
    for (k in 1:K) Z[,,k+1] <- X[,,k]

    Z_norm <- array(NA_real_, dim = c(N, T, K + 1))
    Z_norm[,,1] <- y_norm
    for (k in 1:K) Z_norm[,,k+1] <- X_norm[,,k]

    dist_mat <- pseudo_dist_unit(Z_norm)
    off_diag_dist <- dist_mat[row(dist_mat) != col(dist_mat)]
    quantile_prob <- sqrt(log(N) / K)
    quantile_prob <- min(1, max(0, quantile_prob))
    dist_cutoff <- as.numeric(quantile(off_diag_dist, probs = quantile_prob,
                                       na.rm = TRUE, names = FALSE, type = 7))

    pseudo_clusters <- vector("list", N)
    Z_proj <- array(NA_real_, dim = c(N, T, K + 1))
    for (i in seq_len(N)) {
      cluster_i <- which(dist_mat[i, ] < dist_cutoff)
      cluster_i <- sort(unique(c(i, cluster_i)))
      pseudo_clusters[[i]] <- cluster_i

      for (k in seq_len(K + 1)) {
        group_values <- matrix(Z[cluster_i, , k], nrow = length(cluster_i), ncol = T)
        group_mean <- colMeans(group_values, na.rm = TRUE)
        Z_proj[i, , k] <- Z[i, , k] - group_mean
      }
    }

    y_trans <- Z_proj[,,1]
    X_trans <- Z_proj[,,2:(K+1)]

    tY_vector <- as.vector(y_trans)
    if (T == 1) {
      X_trans_mat <- if (length(dim(X_trans)) == 3) X_trans[, 1, , drop = TRUE] else X_trans
      tX_matrix <- X_trans_mat
    } else if (length(dim(X_trans)) == 2) {
      tX_matrix <- matrix(as.vector(t(X_trans)), nrow = N * T, ncol = K)
    } else {
      tX_matrix <- matrix(aperm(X_trans, c(1, 2, 3)), nrow = N * T, ncol = K)
    }

    d_vector <- tX_matrix[, 1]
    x_matrix <- tX_matrix[, -1, drop = FALSE]

    fit <- rlassoEffect(x = x_matrix, y = tY_vector, d = d_vector, method = "double selection")

    trans <- data.frame(y = tY_vector,
                        D = d_vector,
                        x_matrix)

    lasso.Y <- rlasso(y ~ . - D - 1, data = trans)
    Ytilde <- lasso.Y$residuals

    lasso.D <- rlasso(D ~ . - 1, data = trans[, -1, drop = FALSE])
    Dtilde <- lasso.D$residuals

    data_res <- data.frame(id = rep(ids, T),
                           time = rep(times, each = N),
                           Ytilde = Ytilde,
                           Dtilde = Dtilde)

    if (T == 1) {
      Post_plm <- lm(Ytilde ~ -1 + Dtilde, data = data_res)
      robust_se <- sqrt(diag(vcovHC(Post_plm, type = "HC0")))
    } else {
      Post_plm <- plm(Ytilde ~ -1 + Dtilde, data = data_res, model = "pooling", index = c("id", "time"))
      robust_se <- sqrt(diag(vcovHC(Post_plm, type = "HC0", method = "arellano")))
    }

    coefs <- coef(Post_plm)
    G_unit <- N
    se_corrected <- robust_se
    t_values_corrected <- coefs / se_corrected

    df <- Post_plm$df.residual
    p_values_corrected <- 2 * pt(-abs(t_values_corrected), df)

    summary_table_correct <- data.frame(
      Estimate = coefs,
      `Std. Error corrected` = se_corrected,
      `t-value corrected` = t_values_corrected,
      `Pr(>|t|) corrected` = p_values_corrected
    )

    colnames(summary_table_correct) = c('Estimate', 'Std. Error corrected', 't-value corrected', 'Pr(>|t|) corrected')

    return(list(
      fit_summary = summary(fit),
      G_unit = G_unit,
      unit_group = seq_len(N),
      pseudo_clusters = pseudo_clusters,
      pseudo_cluster_size = lengths(pseudo_clusters),
      pseudo_quantile_prob = quantile_prob,
      pseudo_dist_cutoff = dist_cutoff,
      post_plm_summary = summary(Post_plm),
      estimate_corrected = summary_table_correct,
      summary_table = summary_table_correct
    ))
  }


  if (cluster_type == 'two way'){
    if (cluster_method == 'kmeans'){

      # kmeans
      X_list = lapply(seq_len(dim(X)[3]), function(k) X[ , , k])
      if (pre_cluster == FALSE){
        # cluster
        if (is.null(unit_cluster) == 1 & is.null(time_cluster) == 1){
          clusteri <- cluster_general(y, X_list, N, T, init = 50, type = "long", gamma = gamma, cc = cc, dim_moment = dim_moment  )
          G_unit <- clusteri$clusters
          klong <- clusteri$res


          clustert <- cluster_general(y, X_list, N, T, init = 50, type = "tall", gamma = gamma, cc = cc, dim_moment = dim_moment  )
          G_time <- clustert$clusters
          ktall <- clustert$res
          newcluster_T = fix_time_clusters(cluster = ktall$cluster, data = clustert$data, m = 3)
          ktall$cluster =   newcluster_T$cluster
          G_time =  newcluster_T$G_time
          #   if (G_time > floor(T/2)){
          #     clustert <- cluster_general(y, X_list, N, T, init = 50, type = "tall", groups = 1 )
          #     G_time <- clustert$clusters
          #     ktall <- clustert$res
          #   }
        }else{
          clusteri <- cluster_general(y, X_list, N, T, init = 50, type = "long" , groups = c(floor(unit_cluster)) )
          G_unit <- clusteri$clusters
          klong <- clusteri$res


          clustert <- cluster_general(y, X_list, N, T, init = 50, type = "tall" , groups = c(floor(time_cluster)))
          G_time <- clustert$clusters
          ktall <- clustert$res
        }


        Du <- matrix(0, N, G_unit)
        Dv <- matrix(0, T, G_time)
        Diagu<- matrix(0, G_unit, G_unit)
        Diagv<- matrix(0, G_time, G_time)

        for (j in seq_len(G_unit)) {
          Du[, j] <- as.numeric(klong$cluster == j)
          Diagu[j,j] <- 1/sum(Du[,j])
        }

        for (j in seq_len(G_time)) {
          Dv[, j] <- as.numeric(ktall$cluster == j)
          Diagv[j,j] <- 1/sum(Dv[,j])
        }
      }else if (pre_cluster == TRUE){
        Du = Du_pre
        Dv = Dv_pre
        Diagu<- matrix(0, G_unit, G_unit)
        Diagv<- matrix(0, G_time, G_time)

        for (j in seq_len(G_unit)) {
          Diagu[j,j] <- 1/sum(Du[,j])
        }

        for (j in seq_len(G_time)) {
          Diagv[j,j] <- 1/sum(Dv[,j])
        }
        G_unit = dim(Du)[2]
        G_time = dim(Dv)[2]
      }
    }else if (cluster_method == 'hierarchical'){

      # heriachical
      if (is.null(unit_cluster) == 1 & is.null(time_cluster) == 1){
        res_unit <- cluster_Hierarchical(y_norm, X_norm, link = link, type = "unit", method_auto = method_auto, dist_type = dist_type)
        res_time <- cluster_Hierarchical(y_norm, X_norm, link = link, type = "time", method_auto = method_auto, dist_type = dist_type)
      }else{
        res_unit <- cluster_Hierarchical(y_norm, X_norm, link = link, type = "unit", cluster = unit_cluster, dist_type = dist_type)
        res_time <- cluster_Hierarchical(y_norm, X_norm, link = link, type = "time", cluster = time_cluster, dist_type = dist_type)
      }

      Du <- res_unit$indicator    # N x G_unit
      Dv <- res_time$indicator    # T x G_time
      G_unit = res_unit$G
      G_time = res_time$G
      Diagu<- matrix(0, G_unit, G_unit)
      Diagv<- matrix(0, G_time, G_time)

      for (j in seq_len(G_unit)) {
        Diagu[j,j] <- 1/sum(Du[,j])
      }

      for (j in seq_len(G_time)) {
        Diagv[j,j] <- 1/sum(Dv[,j])
      }
    }

    unit_group <- max.col(Du)
    time_group <- max.col(Dv)

    # Projection matrices
    Mu <- diag(N) - Du %*% Diagu %*% t(Du)
    Mv <- diag(T) - Dv %*% Diagv %*% t(Dv)
    # Mu <- diag(N) - Du %*% solve(t(Du) %*% Du) %*% t(Du)
    # Mv <- diag(T) - Dv %*% solve(t(Dv) %*% Dv) %*% t(Dv)

    # Stack y and X into Z: N x T x (K+1)
    Z <- array(NA_real_, dim = c(N, T, K + 1))
    Z[,,1] <- y
    for (k in 1:K) Z[,,k+1] <- X[,,k]

    # Projection method: demean unit and time clusters slice-wise
    # When T=1, Mv becomes zero matrix, so skip time demeaning
    Z_proj <- array(NA_real_, dim = c(N, T, K + 1))
    for (k in 1:(K+1)) {
      if (T == 1) {
        Z_proj[,,k] <- Mu %*% Z[,,k]
      } else {
        Z_proj[,,k] <- Mu %*% Z[,,k] %*% Mv
      }
    }

    # Separate transformed y and X
    y_trans <- Z_proj[,,1] # Z_proj_final_array[,,1]
    X_trans <- Z_proj[,,2:(K+1)] # Z_proj_final_array[,,2:(K+1)]

    tY_vector <- as.vector(y_trans)  # (N*T)
    # Keep the same T=1 flattening convention as Niave/TWFE
    if (T == 1) {
      X_trans_mat <- if (length(dim(X_trans)) == 3) X_trans[, 1, , drop = TRUE] else X_trans
      tX_matrix <- X_trans_mat
    } else if (length(dim(X_trans)) == 2) {
      tX_matrix <- matrix(as.vector(t(X_trans)), nrow = N * T, ncol = K)
    } else {
      tX_matrix <- matrix(aperm(X_trans, c(1, 2, 3)), nrow = N * T, ncol = K)
    }

    d_vector <- tX_matrix[, 1]
    x_matrix <- tX_matrix[, -1, drop = FALSE]

    # Run rlassoEffect with correct inputs
    fit <- rlassoEffect(x = x_matrix, y = tY_vector, d = d_vector, method = "double selection")

    trans <- data.frame(y = tY_vector,
                        D = d_vector,
                        x_matrix)

    lasso.Y <- rlasso(y ~ . - D - 1, data = trans)
    Ytilde <- lasso.Y$residuals

    lasso.D <- rlasso(D ~ . - 1, data = trans[, -1, drop = FALSE])
    Dtilde <- lasso.D$residuals

    # Create proper panel structure for residuals
    # as.vector() on N x T matrix goes column-wise: unit 1 time 1, unit 2 time 1, ..., unit N time 1, unit 1 time 2, ...
    data_res <- data.frame(id = rep(ids, T),
                           time = rep(times, each = N),
                           Ytilde = Ytilde,
                           Dtilde = Dtilde)

    # For T=1, use lm() since there's no panel structure; for T>1 use plm with pooling
    if (T == 1) {
      Post_plm <- lm(Ytilde ~ -1 + Dtilde, data = data_res)
    } else {
      Post_plm <- plm(Ytilde ~ -1 + Dtilde, data = data_res, model = "pooling", index = c("id", "time"))
    }

    coefs <- coef(Post_plm)
    if (N*T - N*G_time - T*G_unit > 0){
      if (T == 1) {
        se_corrected <- sqrt(diag(vcovHC(Post_plm, type = "HC0"))) * sqrt(N * T / (N*T - N*G_time - T*G_unit))
      } else {
        se_corrected <- sqrt(diag(vcovHC(Post_plm, type = "HC0", method = "arellano"))) * sqrt(N * T / (N*T - N*G_time - T*G_unit))
      }
    }else{
      if (T == 1) {
        se_corrected <- sqrt(diag(vcovHC(Post_plm, type = "HC0"))) * sqrt(N * T / (1))
      } else {
        se_corrected <- sqrt(diag(vcovHC(Post_plm, type = "HC0", method = "arellano"))) * sqrt(N * T / (1))
      }
    }
    t_values_corrected <- coefs / se_corrected

    df <- Post_plm$df.residual
    p_values_corrected <- 2 * pt(-abs(t_values_corrected), df)

    summary_table_correct <- data.frame(
      Estimate = coefs,
      `Std. Error corrected` = se_corrected,
      `t-value corrected` = t_values_corrected,
      `Pr(>|t|) corrected` = p_values_corrected
    )

    colnames(summary_table_correct) = c('Estimate', 'Std. Error corrected', 't-value corrected', 'Pr(>|t|) corrected')

    # Return cluster counts as well
    return(list(
      fit_summary = summary(fit),
      G_unit = G_unit,
      G_time = G_time,
      unit_group = unit_group,
      time_group = time_group,
      post_plm_summary = summary(Post_plm),
      estimate_corrected = summary_table_correct,
      summary_table = summary_table_correct
    ))
  }else if (cluster_type == 'one way'){
    if (cluster_method == 'kmeans'){

      # kmeans
      X_list = lapply(seq_len(dim(X)[3]), function(k) X[ , , k])
      if (pre_cluster == FALSE){
        # cluster
        if (is.null(unit_cluster) == 1 & is.null(time_cluster) == 1){
          clusteri <- cluster_general(y, X_list, N, T, init = 50, type = "long", gamma = gamma, cc = cc, dim_moment = dim_moment  )
          G_unit <- clusteri$clusters
          klong <- clusteri$res


        }else{
          clusteri <- cluster_general(y, X_list, N, T, init = 50, type = "long" , groups = c(floor(unit_cluster)) )
          G_unit <- clusteri$clusters
          klong <- clusteri$res
        }
        Du <- matrix(0, N, G_unit)
        Diagu<- matrix(0, G_unit, G_unit)


        for (j in seq_len(G_unit)) {
          Du[, j] <- as.numeric(klong$cluster == j)
          Diagu[j,j] <- 1/sum(Du[,j])
        }


      }else if (pre_cluster == TRUE){
        Du = Du_pre
        G_unit = dim(Du)[2]
        Diagu<- matrix(0, G_unit, G_unit)

        for (j in seq_len(G_unit)) {
          Diagu[j,j] <- 1/sum(Du[,j])
        }
      }
    }else if (cluster_method == 'hierarchical'){

      # heriachical
      if (is.null(unit_cluster) == 1 & is.null(time_cluster) == 1){
        res_unit <- cluster_Hierarchical(y_norm, X_norm, link = link, type = "unit", method_auto = method_auto, dist_type = dist_type)
      }else{
        res_unit <- cluster_Hierarchical(y_norm, X_norm, link = link, type = "unit", cluster = unit_cluster, dist_type = dist_type)
      }

      Du <- res_unit$indicator    # N x G_unit
      G_unit = res_unit$G
      Diagu<- matrix(0, G_unit, G_unit)

      for (j in seq_len(G_unit)) {
        Diagu[j,j] <- 1/sum(Du[,j])
      }

    }

    unit_group <- max.col(Du)

    # Projection matrices
    Mu <- diag(N) - Du %*% Diagu %*% t(Du)

    # Stack y and X into Z: N x T x (K+1)
    Z <- array(NA_real_, dim = c(N, T, K + 1))
    Z[,,1] <- y
    for (k in 1:K) Z[,,k+1] <- X[,,k]

    # Projection method: demean unit and time clusters slice-wise
    Z_proj <- array(NA_real_, dim = c(N, T, K + 1))
    for (k in 1:(K+1)) {
      Z_proj[,,k] <- Mu %*% Z[,,k]
    }

    # Separate transformed y and X
    y_trans <- Z_proj[,,1] # Z_proj_final_array[,,1]
    X_trans <- Z_proj[,,2:(K+1)] # Z_proj_final_array[,,2:(K+1)]

    tY_vector <- as.vector(y_trans)  # (N*T)
    if (T == 1) {
      X_trans_mat <- if (length(dim(X_trans)) == 3) X_trans[, 1, , drop = TRUE] else X_trans
      tX_matrix <- X_trans_mat
    } else if (length(dim(X_trans)) == 2) {
      tX_matrix <- matrix(as.vector(t(X_trans)), nrow = N * T, ncol = K)
    } else {
      tX_matrix <- matrix(aperm(X_trans, c(1, 2, 3)), nrow = N * T, ncol = K)
    }
    # debiased lasso
    # debiased_lasso = lasso.proj(tX_matrix, tY_vector) # delas(tX_matrix, tY_vector, c(1))
    # debiased_lasso$se = debiased_lasso$se * sqrt(N * T / (N*T - T*G_unit))

    #
    d_vector <- tX_matrix[, 1]
    x_matrix <- tX_matrix[, -1, drop = FALSE]

    # Run rlassoEffect with correct inputs
    fit <- rlassoEffect(x = x_matrix, y = tY_vector, d = d_vector, method = "double selection")

    trans <- data.frame(y = tY_vector,
                        D = d_vector,
                        x_matrix)

    lasso.Y <- rlasso(y ~ . - D - 1, data = trans)
    Ytilde <- lasso.Y$residuals

    lasso.D <- rlasso(D ~ . - 1, data = trans[, -1, drop = FALSE])
    Dtilde <- lasso.D$residuals

    data_res <- data.frame(id = rep(ids, T),
                           time = rep(times, each = N),
                           Ytilde = Ytilde,
                           Dtilde = Dtilde)

    if (T == 1) {
      Post_plm <- lm(Ytilde ~ -1 + Dtilde, data = data_res)
      robust_se <- sqrt(diag(vcovHC(Post_plm, type = "HC0")))
    } else {
      Post_plm <- plm(Ytilde ~ -1 + Dtilde, data = data_res, model = "pooling", index = c("id", "time"))
      robust_se <- sqrt(diag(vcovHC(Post_plm, type = "HC0", method = "arellano")))
    }

    coefs <- coef(Post_plm)
    se_corrected <- robust_se * sqrt(N * T / max(1, N * T - T * G_unit))
    t_values_corrected <- coefs / se_corrected

    df <- Post_plm$df.residual
    p_values_corrected <- 2 * pt(-abs(t_values_corrected), df)

    summary_table_correct <- data.frame(
      Estimate = coefs,
      `Std. Error corrected` = se_corrected,
      `t-value corrected` = t_values_corrected,
      `Pr(>|t|) corrected` = p_values_corrected
    )

    colnames(summary_table_correct) = c('Estimate', 'Std. Error corrected', 't-value corrected', 'Pr(>|t|) corrected')

    # Return cluster counts as well
    return(list(
      fit_summary = summary(fit),
      G_unit = G_unit,
      unit_group = unit_group,
      post_plm_summary = summary(Post_plm),
      estimate_corrected = summary_table_correct,
      summary_table = summary_table_correct
    ))

  }else if (cluster_type == 'one way T moment'){
    if (cluster_method == 'kmeans'){

      # kmeans
      X_list = lapply(seq_len(dim(X)[3]), function(k) X[ , , k])
      if (pre_cluster == FALSE){
        # cluster
        if (is.null(unit_cluster) == 1 & is.null(time_cluster) == 1){
          clusteri <- cluster_general(y, X_list, N, T, init = 50, type = "long_T_moment", gamma = gamma, cc = cc, dim_moment = dim_moment  )
          G_unit <- clusteri$clusters
          klong <- clusteri$res


        }else{
          clusteri <- cluster_general(y, X_list, N, T, init = 50, type = "long_T_moment" , groups = c(floor(unit_cluster)) )
          G_unit <- clusteri$clusters
          klong <- clusteri$res
        }
        Du <- matrix(0, N, G_unit)
        Diagu<- matrix(0, G_unit, G_unit)


        for (j in seq_len(G_unit)) {
          Du[, j] <- as.numeric(klong$cluster == j)
          Diagu[j,j] <- 1/sum(Du[,j])
        }


      }else if (pre_cluster == TRUE){
        Du = Du_pre
        G_unit = dim(Du)[2]
        Diagu<- matrix(0, G_unit, G_unit)

        for (j in seq_len(G_unit)) {
          Diagu[j,j] <- 1/sum(Du[,j])
        }
      }
    }else if (cluster_method == 'hierarchical'){

      # heriachical
      if (is.null(unit_cluster) == 1 & is.null(time_cluster) == 1){
        res_unit <- cluster_Hierarchical(y_norm, X_norm, link = link, type = "unit", method_auto = method_auto, dist_type = dist_type)
      }else{
        res_unit <- cluster_Hierarchical(y_norm, X_norm, link = link, type = "unit", cluster = unit_cluster, dist_type = dist_type)
      }

      Du <- res_unit$indicator    # N x G_unit
      G_unit = res_unit$G
      Diagu<- matrix(0, G_unit, G_unit)

      for (j in seq_len(G_unit)) {
        Diagu[j,j] <- 1/sum(Du[,j])
      }

    }

    unit_group <- max.col(Du)

    # Projection matrices
    Mu <- diag(N) - Du %*% Diagu %*% t(Du)

    # Stack y and X into Z: N x T x (K+1)
    Z <- array(NA_real_, dim = c(N, T, K + 1))
    Z[,,1] <- y
    for (k in 1:K) Z[,,k+1] <- X[,,k]

    # Projection method: demean unit and time clusters slice-wise
    Z_proj <- array(NA_real_, dim = c(N, T, K + 1))
    for (k in 1:(K+1)) {
      Z_proj[,,k] <- Mu %*% Z[,,k]
    }

    # Separate transformed y and X
    y_trans <- Z_proj[,,1] # Z_proj_final_array[,,1]
    X_trans <- Z_proj[,,2:(K+1)] # Z_proj_final_array[,,2:(K+1)]

    tY_vector <- as.vector(y_trans)  # (N*T)
    if (T == 1) {
      X_trans_mat <- if (length(dim(X_trans)) == 3) X_trans[, 1, , drop = TRUE] else X_trans
      tX_matrix <- X_trans_mat
    } else if (length(dim(X_trans)) == 2) {
      tX_matrix <- matrix(as.vector(t(X_trans)), nrow = N * T, ncol = K)
    } else {
      tX_matrix <- matrix(aperm(X_trans, c(1, 2, 3)), nrow = N * T, ncol = K)
    }

    # debiased lasso
    # debiased_lasso = lasso.proj(tX_matrix, tY_vector) # delas(tX_matrix, tY_vector, c(1))
    # debiased_lasso$se = debiased_lasso$se * sqrt(N * T / (N*T - T*G_unit))

    #
    d_vector <- tX_matrix[, 1]
    x_matrix <- tX_matrix[, -1, drop = FALSE]

    # Run rlassoEffect with correct inputs
    fit <- rlassoEffect(x = x_matrix, y = tY_vector, d = d_vector, method = "double selection")

    trans <- data.frame(y = tY_vector,
                        D = d_vector,
                        x_matrix)

    lasso.Y <- rlasso(y ~ . - D - 1, data = trans)
    Ytilde <- lasso.Y$residuals

    lasso.D <- rlasso(D ~ . - 1, data = trans[, -1, drop = FALSE])
    Dtilde <- lasso.D$residuals

    data_res <- data.frame(id = rep(ids, T),
                           time = rep(times, each = N),
                           Ytilde = Ytilde,
                           Dtilde = Dtilde)

    if (T == 1) {
      Post_plm <- lm(Ytilde ~ -1 + Dtilde, data = data_res)
      robust_se <- sqrt(diag(vcovHC(Post_plm, type = "HC0")))
    } else {
      Post_plm <- plm(Ytilde ~ -1 + Dtilde, data = data_res, model = "pooling", index = c("id", "time"))
      robust_se <- sqrt(diag(vcovHC(Post_plm, type = "HC0", method = "arellano")))
    }

    coefs <- coef(Post_plm)
    se_corrected <- robust_se * sqrt(N * T / max(1, N * T - T * G_unit))
    t_values_corrected <- coefs / se_corrected

    df <- Post_plm$df.residual
    p_values_corrected <- 2 * pt(-abs(t_values_corrected), df)

    summary_table_correct <- data.frame(
      Estimate = coefs,
      `Std. Error corrected` = se_corrected,
      `t-value corrected` = t_values_corrected,
      `Pr(>|t|) corrected` = p_values_corrected
    )

    colnames(summary_table_correct) = c('Estimate', 'Std. Error corrected', 't-value corrected', 'Pr(>|t|) corrected')

    # Return cluster counts as well
    return(list(
      fit_summary = summary(fit),
      G_unit = G_unit,
      unit_group = unit_group,
      post_plm_summary = summary(Post_plm),
      estimate_corrected = summary_table_correct,
      summary_table = summary_table_correct
    ))

  }else if (cluster_type == 'one way more moments'){
    if (cluster_method == 'kmeans'){

      # kmeans
      X_list = lapply(seq_len(dim(X)[3]), function(k) X[ , , k])
      if (pre_cluster == FALSE){
        # cluster
        if (is.null(unit_cluster) == 1 & is.null(time_cluster) == 1){
          clusteri <- cluster_general(y, X_list, N, T, init = 50, type = "long multiple moments", gamma = gamma, cc = cc, dim_moment = dim_moment  )
          G_unit <- clusteri$clusters
          klong <- clusteri$res


        }else{
          clusteri <- cluster_general(y, X_list, N, T, init = 50, type = "long multiple moments" , groups = c(floor(unit_cluster)) )
          G_unit <- clusteri$clusters
          klong <- clusteri$res
        }
        Du <- matrix(0, N, G_unit)
        Diagu<- matrix(0, G_unit, G_unit)


        for (j in seq_len(G_unit)) {
          Du[, j] <- as.numeric(klong$cluster == j)
          Diagu[j,j] <- 1/sum(Du[,j])
        }


      }else if (pre_cluster == TRUE){
        Du = Du_pre
        G_unit = dim(Du)[2]
        Diagu<- matrix(0, G_unit, G_unit)

        for (j in seq_len(G_unit)) {
          Diagu[j,j] <- 1/sum(Du[,j])
        }
      }
    }else if (cluster_method == 'hierarchical'){

      # heriachical
      if (is.null(unit_cluster) == 1 & is.null(time_cluster) == 1){
        res_unit <- cluster_Hierarchical(y_norm, X_norm, link = link, type = "unit", method_auto = method_auto, dist_type = dist_type)
      }else{
        res_unit <- cluster_Hierarchical(y_norm, X_norm, link = link, type = "unit", cluster = unit_cluster, dist_type = dist_type)
      }

      Du <- res_unit$indicator    # N x G_unit
      G_unit = res_unit$G
      Diagu<- matrix(0, G_unit, G_unit)

      for (j in seq_len(G_unit)) {
        Diagu[j,j] <- 1/sum(Du[,j])
      }

    }

    unit_group <- max.col(Du)

    # Projection matrices
    Mu <- diag(N) - Du %*% Diagu %*% t(Du)

    # Stack y and X into Z: N x T x (K+1)
    Z <- array(NA_real_, dim = c(N, T, K + 1))
    Z[,,1] <- y
    for (k in 1:K) Z[,,k+1] <- X[,,k]

    # Projection method: demean unit and time clusters slice-wise
    Z_proj <- array(NA_real_, dim = c(N, T, K + 1))
    for (k in 1:(K+1)) {
      Z_proj[,,k] <- Mu %*% Z[,,k]
    }

    # Separate transformed y and X
    y_trans <- Z_proj[,,1] # Z_proj_final_array[,,1]
    X_trans <- Z_proj[,,2:(K+1)] # Z_proj_final_array[,,2:(K+1)]

    tY_vector <- as.vector(y_trans)  # (N*T)
    if (T == 1) {
      X_trans_mat <- if (length(dim(X_trans)) == 3) X_trans[, 1, , drop = TRUE] else X_trans
      tX_matrix <- X_trans_mat
    } else if (length(dim(X_trans)) == 2) {
      tX_matrix <- matrix(as.vector(t(X_trans)), nrow = N * T, ncol = K)
    } else {
      tX_matrix <- matrix(aperm(X_trans, c(1, 2, 3)), nrow = N * T, ncol = K)
    }

    # debiased lasso
    # debiased_lasso = lasso.proj(tX_matrix, tY_vector) # delas(tX_matrix, tY_vector, c(1))
    # debiased_lasso$se = debiased_lasso$se * sqrt(N * T / (N*T - T*G_unit))

    #
    d_vector <- tX_matrix[, 1]
    x_matrix <- tX_matrix[, -1, drop = FALSE]

    # Run rlassoEffect with correct inputs
    fit <- rlassoEffect(x = x_matrix, y = tY_vector, d = d_vector, method = "double selection")

    trans <- data.frame(y = tY_vector,
                        D = d_vector,
                        x_matrix)

    lasso.Y <- rlasso(y ~ . - D - 1, data = trans)
    Ytilde <- lasso.Y$residuals

    lasso.D <- rlasso(D ~ . - 1, data = trans[, -1, drop = FALSE])
    Dtilde <- lasso.D$residuals

    data_res <- data.frame(id = rep(ids, T),
                           time = rep(times, each = N),
                           Ytilde = Ytilde,
                           Dtilde = Dtilde)

    if (T == 1) {
      Post_plm <- lm(Ytilde ~ -1 + Dtilde, data = data_res)
      robust_se <- sqrt(diag(vcovHC(Post_plm, type = "HC0")))
    } else {
      Post_plm <- plm(Ytilde ~ -1 + Dtilde, data = data_res, model = "pooling", index = c("id", "time"))
      robust_se <- sqrt(diag(vcovHC(Post_plm, type = "HC0", method = "arellano")))
    }

    coefs <- coef(Post_plm)
    se_corrected <- robust_se * sqrt(N * T / max(1, N * T - T * G_unit))
    t_values_corrected <- coefs / se_corrected

    df <- Post_plm$df.residual
    p_values_corrected <- 2 * pt(-abs(t_values_corrected), df)

    summary_table_correct <- data.frame(
      Estimate = coefs,
      `Std. Error corrected` = se_corrected,
      `t-value corrected` = t_values_corrected,
      `Pr(>|t|) corrected` = p_values_corrected
    )

    colnames(summary_table_correct) = c('Estimate', 'Std. Error corrected', 't-value corrected', 'Pr(>|t|) corrected')

    # Return cluster counts as well
    return(list(
      fit_summary = summary(fit),
      G_unit = G_unit,
      unit_group = unit_group,
      post_plm_summary = summary(Post_plm),
      estimate_corrected = summary_table_correct,
      summary_table = summary_table_correct
    ))

  }
}

#' @export
HP_estimate_cf <- function(data, y_col = NULL, covariate_cols = NULL,
                           id_col = "id", time_col = "time", cluster_type = c('one way', 'one way T moment', 'one way more moments'),
                           index = c("id", "time"), method_auto = "dynamicTreeCut", cluster_method = c('kmeans', 'hierarchical', 'ps-kmeans', 'HAC new', 'kmedoid', 'projection-kmeans'), link = 'average', pre_cluster = FALSE, Du_pre, Dv_pre, unit_cluster = NULL, time_cluster = NULL, gamma = 1, cc = 0, dim_moment = 1,  dist_type = 'L2',
                           cross_type = c("unit", 'time'), ps_nstart = 30L,
                           ps_max_iter = 100L, ps_tol = 1e-8,
                           ps_seed = NULL, ps_standardize = TRUE,
                           ps_G_rule = c("one-se", "minimum", "n-third"),
                           ps_G_grid = NULL, ps_cv_folds = 2L,
                           ps_cv_repeats = 3L, ps_select_nstart = 3L,
                           ps_min_group_size = 2L,
                           kmedoid_nstart = 30L, kmedoid_max_iter = 100L,
                           kmedoid_tol = 1e-8, kmedoid_seed = NULL,
                           kmedoid_standardize = TRUE,
                           kmedoid_G_rule = c("inference", "theoretical", "elbow", "n-third"),
                           kmedoid_G_grid = NULL,
                           kmedoid_min_group_size = 2L,
                           kmedoid_solver = c("heuristic", "auto", "exact"),
                           kmedoid_exact_limit = 50000L,
                           pk_nstart = 30L, pk_max_iter = 100L,
                           pk_select_nstart = 3L, pk_select_max_iter = 25L,
                           pk_tol = 1e-8, pk_seed = NULL,
                           pk_standardize = FALSE,
                           pk_G_rule = c("noise-floor", "inference", "elbow", "n-third"),
                           pk_G_grid = NULL, pk_max_G = NULL,
                           pk_log_N = "log",
                           pk_min_group_size = 2L,
                           pk_noise_multiplier = 1,
                           hac_G_rule = c("one-se", "minimum"),
                           hac_G_grid = NULL, hac_cv_folds = 2L,
                           hac_cv_repeats = 3L, hac_min_group_size = 2L,
                           hac_reference_fraction = 0.5, hac_seed = NULL,
                           hac_standardize = TRUE,
                           all_targets = FALSE) {


  cross_type_was_missing <- missing(cross_type)
  if (is.null(covariate_cols)) {
    covariate_cols <- setdiff(colnames(data), c(id_col, time_col, y_col))
  }

  if (length(all_targets) != 1L || is.na(all_targets) || !is.logical(all_targets)) {
    stop("all_targets must be TRUE or FALSE.")
  }

  if (isTRUE(all_targets)) {
    if (!length(covariate_cols)) {
      stop("all_targets = TRUE requires at least one covariate.")
    }

    caller_env <- parent.frame()
    base_call <- match.call()
    base_call$all_targets <- FALSE

    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      invisible(runif(1L))
    }
    rng_start <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    rng_after_first <- NULL

    target_fits <- vector("list", length(covariate_cols))
    for (k in seq_along(covariate_cols)) {
      assign(".Random.seed", rng_start, envir = .GlobalEnv)
      call_k <- base_call
      call_k$covariate_cols <- c(covariate_cols[k], covariate_cols[-k])
      target_fits[[k]] <- eval(call_k, envir = caller_env)
      if (k == 1L) {
        rng_after_first <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
      }
    }
    if (!is.null(rng_after_first)) {
      assign(".Random.seed", rng_after_first, envir = .GlobalEnv)
    }

    summary_all <- do.call(rbind, lapply(seq_along(covariate_cols), function(k) {
      tab <- target_fits[[k]]$summary_table
      if (is.null(tab) || nrow(tab) < 1L) {
        stop("No summary-table row was returned for target ", covariate_cols[k], ".")
      }
      data.frame(
        Variable = covariate_cols[k],
        tab[1L, , drop = FALSE],
        row.names = NULL,
        check.names = FALSE
      )
    }))
    rownames(summary_all) <- covariate_cols

    result <- target_fits[[1L]]
    result$summary_table <- summary_all
    result$estimate_corrected <- summary_all
    result$all_targets <- TRUE
    result$target_names <- covariate_cols
    result$G_unit_all <- setNames(
      vapply(target_fits, function(z) as.numeric(z$G_unit)[1L], numeric(1L)),
      covariate_cols
    )
    if (!is.null(target_fits[[1L]]$G_unit_fold1)) {
      result$G_unit_fold1_all <- setNames(
        vapply(target_fits, function(z) as.numeric(z$G_unit_fold1)[1L], numeric(1L)),
        covariate_cols
      )
      result$G_unit_fold2_all <- setNames(
        vapply(target_fits, function(z) as.numeric(z$G_unit_fold2)[1L], numeric(1L)),
        covariate_cols
      )
    }
    return(result)
  }

  ids <- sort(unique(data[[id_col]]))
  times <- sort(unique(data[[time_col]]))
  N <- length(ids)
  T <- length(times)
  K <- length(covariate_cols)
  cluster_type <- match.arg(cluster_type)
  cluster_method <- match.arg(cluster_method)
  if (cluster_method %in% c("kmedoid", "projection-kmeans") &&
      cross_type_was_missing) {
    cross_type <- "time"
  }
  cross_type <- match.arg(cross_type)

  if (cluster_method %in% c("ps-kmeans", "HAC new") && cross_type != "unit") {
    stop("The cross-fitted profile methods are implemented with cross_type = 'unit'.")
  }
  if (cluster_method %in% c("kmedoid", "projection-kmeans") &&
      cross_type != "time") {
    stop("Cross-fitted k-medoids and projection K-means require cross_type = 'time'.")
  }
  if (cluster_method %in% c("ps-kmeans", "HAC new", "kmedoid",
                            "projection-kmeans") && pre_cluster) {
    stop("The cross-fitted profile methods require pre_cluster = FALSE.")
  }
  if (cluster_method %in% c("ps-kmeans", "HAC new", "kmedoid",
                            "projection-kmeans")) {
    if (length(y_col) != 1L || !is.character(y_col) || !nzchar(y_col)) {
      stop("The profile clustering methods require y_col to name exactly one outcome column.")
    }
    required_profile_columns <- unique(c(id_col, time_col, y_col, covariate_cols))
    missing_profile_columns <- setdiff(required_profile_columns, names(data))
    if (length(missing_profile_columns)) {
      stop("The profile clustering data are missing columns: ",
           paste(missing_profile_columns, collapse = ", "), ".")
    }
    if (K < 1L) {
      stop("The profile clustering methods require at least one covariate; list treatment first in covariate_cols.")
    }
    if (anyDuplicated(data[c(id_col, time_col)])) {
      stop("The profile clustering methods require exactly one row for each (unit, time) pair.")
    }
    if (nrow(data) != N * T) {
      stop("The profile clustering methods require a complete balanced panel.")
    }
  }
  if (cluster_method == "ps-kmeans") {
    if (!is.null(ps_seed) &&
        (length(ps_seed) != 1L || !is.finite(ps_seed) ||
         ps_seed < 0 || ps_seed > .Machine$integer.max ||
         ps_seed != floor(ps_seed))) {
      stop("ps_seed must be NULL or one nonnegative integer.")
    }
    if (length(ps_nstart) != 1L || !is.finite(ps_nstart) ||
        ps_nstart < 1 || ps_nstart != floor(ps_nstart)) {
      stop("ps_nstart must be one positive integer.")
    }
    if (length(ps_max_iter) != 1L || !is.finite(ps_max_iter) ||
        ps_max_iter < 1 || ps_max_iter != floor(ps_max_iter)) {
      stop("ps_max_iter must be one positive integer.")
    }
    if (length(ps_tol) != 1L || !is.finite(ps_tol) || ps_tol <= 0) {
      stop("ps_tol must be one positive finite number.")
    }
    if (length(ps_standardize) != 1L || is.na(ps_standardize) ||
        !is.logical(ps_standardize)) {
      stop("ps_standardize must be TRUE or FALSE.")
    }
    ps_G_rule <- match.arg(ps_G_rule)
    if (!is.null(ps_G_grid)) {
      if (!is.numeric(ps_G_grid) || !length(ps_G_grid) ||
          any(!is.finite(ps_G_grid)) || any(ps_G_grid < 1) ||
          any(ps_G_grid != floor(ps_G_grid))) {
        stop("ps_G_grid must be NULL or a nonempty vector of positive integers.")
      }
      ps_G_grid <- sort(unique(as.integer(ps_G_grid)))
    }
    for (parameter in c("ps_cv_folds", "ps_cv_repeats",
                        "ps_select_nstart", "ps_min_group_size")) {
      value <- get(parameter)
      if (length(value) != 1L || !is.finite(value) ||
          value < 1 || value != floor(value)) {
        stop(parameter, " must be one positive integer.")
      }
    }
    if (ps_cv_folds < 2L) {
      stop("ps_cv_folds must be at least 2.")
    }
    ps_nstart <- as.integer(ps_nstart)
    ps_max_iter <- as.integer(ps_max_iter)
    ps_cv_folds <- as.integer(ps_cv_folds)
    ps_cv_repeats <- as.integer(ps_cv_repeats)
    ps_select_nstart <- as.integer(ps_select_nstart)
    ps_min_group_size <- as.integer(ps_min_group_size)
  }
  if (cluster_method == "HAC new") {
    hac_G_rule <- match.arg(hac_G_rule)
    if (!identical(link, "average")) {
      stop("HAC new uses average linkage; set link = 'average'.")
    }
    if (!is.null(hac_G_grid)) {
      if (!is.numeric(hac_G_grid) || !length(hac_G_grid) ||
          any(!is.finite(hac_G_grid)) || any(hac_G_grid < 1) ||
          any(hac_G_grid != floor(hac_G_grid))) {
        stop("hac_G_grid must be NULL or a nonempty vector of positive integers.")
      }
      hac_G_grid <- sort(unique(as.integer(hac_G_grid)))
    }
    for (parameter in c("hac_cv_folds", "hac_cv_repeats",
                        "hac_min_group_size")) {
      value <- get(parameter)
      if (length(value) != 1L || !is.finite(value) ||
          value < 1 || value != floor(value)) {
        stop(parameter, " must be one positive integer.")
      }
    }
    if (hac_cv_folds < 2L) stop("hac_cv_folds must be at least 2.")
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
    hac_cv_folds <- as.integer(hac_cv_folds)
    hac_cv_repeats <- as.integer(hac_cv_repeats)
    hac_min_group_size <- as.integer(hac_min_group_size)
  }
  if (cluster_method == "kmedoid") {
    kmedoid_options <- .kmedoid_validate_options(
      kmedoid_nstart, kmedoid_max_iter, kmedoid_tol, kmedoid_seed,
      kmedoid_standardize, kmedoid_G_rule, kmedoid_G_grid,
      kmedoid_min_group_size, kmedoid_solver, kmedoid_exact_limit
    )
    if (N < 3L) {
      stop("Cross-fitted k-medoids requires at least three units.")
    }
  } else {
    kmedoid_options <- NULL
  }
  if (cluster_method == "projection-kmeans") {
    pk_options <- .pk_validate_options(
      pk_nstart, pk_max_iter, pk_select_nstart, pk_select_max_iter,
      pk_tol, pk_seed, pk_standardize,
      pk_G_rule, pk_G_grid, pk_min_group_size, pk_noise_multiplier,
      max_G = pk_max_G, log_N = pk_log_N
    )
    if (N < 2L) {
      stop("Cross-fitted projection K-means requires at least two units.")
    }
  } else {
    pk_options <- NULL
  }

  if (T <= 1) {
    stop("HP_estimate_cf is implemented only for T > 1.")
  }
  cf_long_type <- if (cluster_type == "one way") {
    "long"
  } else if (cluster_type == "one way T moment") {
    "long_T_moment"
  } else {
    "long multiple moments"
  }

  # Initialize y and X as balanced panel containers.
  id_index <- match(data[[id_col]], ids)
  time_index <- match(data[[time_col]], times)
  panel_index <- cbind(id_index, time_index)
  y <- matrix(NA_real_, nrow = N, ncol = T)
  y[panel_index] <- data[[y_col]]
  X <- array(NA_real_, dim = c(N, T, K))
  for (k in seq_along(covariate_cols)) {
    X[cbind(id_index, time_index, k)] <- data[[covariate_cols[k]]]
  }

  unit_feature_data <- function(y_sub, X_sub, center = NULL, scale = NULL) {
    N_sub <- nrow(y_sub)
    X_list_sub <- c(list(y_sub), lapply(seq_len(dim(X_sub)[3]), function(k) X_sub[, , k]))

    if (cf_long_type == "long") {
      mom_i <- matrix(0, nrow = N_sub, ncol = dim_moment)
      covar_groups <- if (dim_moment == 1) {
        list(seq_len(K + 1))
      } else {
        split(seq_len(K + 1), cut(seq_len(K + 1), dim_moment, labels = FALSE))
      }

      for (p in seq_len(dim_moment)) {
        X_pow <- matrix(0, nrow = N_sub, ncol = T)
        for (tt in seq_len(T)) {
          X_t <- do.call(cbind, lapply(X_list_sub, function(X_item) X_item[, tt]))
          X_pow[, tt] <- rowMeans(X_t[, covar_groups[[p]], drop = FALSE])
        }
        mom_i[, p] <- rowMeans(X_pow)
      }
    } else if (cf_long_type == "long_T_moment") {
      dm <- if (T == 1) 1 else dim_moment
      mom_i <- matrix(0, nrow = N_sub, ncol = dm)
      time_groups <- if (T == 1 || dm <= 1) {
        list(seq_len(T))
      } else {
        split(seq_len(T), cut(seq_len(T), dm, labels = FALSE))
      }
      X_mean_t <- matrix(0, nrow = N_sub, ncol = T)

      for (tt in seq_len(T)) {
        X_t <- do.call(cbind, lapply(X_list_sub, function(X_item) X_item[, tt]))
        X_mean_t[, tt] <- rowMeans(X_t)
      }
      for (p in seq_len(dm)) {
        mom_i[, p] <- rowMeans(X_mean_t[, time_groups[[p]], drop = FALSE])
      }
    } else {
      mom_i <- matrix(0, nrow = N_sub, ncol = dim_moment)
      for (p in seq_len(dim_moment)) {
        X_pow <- matrix(0, nrow = N_sub, ncol = T)
        for (tt in seq_len(T)) {
          X_t <- do.call(cbind, lapply(X_list_sub, function(X_item) X_item[, tt]))
          if (p == 1) {
            X_pow[, tt] <- rowMeans(X_t)
          } else if (p == 2) {
            X_pow[, tt] <- rowMeans(exp(X_t))
          } else if (p == 3) {
            X_pow[, tt] <- rowMeans(X_t / (1 + abs(X_t)))
          } else {
            X_pow[, tt] <- rowMeans(X_t^p)
          }
        }
        mom_i[, p] <- rowMeans(X_pow)
      }
    }

    if (is.null(center) || is.null(scale)) {
      center <- colMeans(mom_i, na.rm = TRUE)
      scale <- apply(mom_i, 2, sd, na.rm = TRUE)
      scale[!is.finite(scale)] <- 0
    }

    data_scaled <- mom_i
    for (j in seq_len(ncol(mom_i))) {
      if (scale[j] > 0) {
        data_scaled[, j] <- (mom_i[, j] - center[j]) / scale[j]
      }
    }

    list(data = data_scaled, center = center, scale = scale)
  }

  nearest_center <- function(data_test, centers) {
    centers <- as.matrix(centers)
    if (nrow(centers) == 1) {
      return(rep(1L, nrow(data_test)))
    }
    as.integer(apply(data_test, 1, function(row_i) {
      dists <- rowSums(sweep(centers, 2, row_i, "-")^2)
      which.min(dists)
    }))
  }

  profile_data <- function(y_sub, X_sub, center = NULL, scale = NULL,
                           standardize = ps_standardize) {
    N_sub <- nrow(y_sub)
    M <- T * (K + 1L)
    profile_raw <- matrix(NA_real_, nrow = N_sub, ncol = M)
    colnames(profile_raw) <- character(M)
    component_names <- c(y_col, covariate_cols)

    for (tt in seq_len(T)) {
      columns <- (tt - 1L) * (K + 1L) + seq_len(K + 1L)
      X_tt <- matrix(X_sub[, tt, , drop = FALSE], nrow = N_sub, ncol = K)
      profile_raw[, columns] <- cbind(y_sub[, tt], X_tt)
      colnames(profile_raw)[columns] <- paste0("t=", times[tt], "::", component_names)
    }

    if (any(!is.finite(profile_raw))) {
      stop("The profile clustering methods require finite y and X values.")
    }

    if (is.null(center) || is.null(scale)) {
      if (standardize) {
        center <- apply(profile_raw, 2L, stats::median)
        raw_scale <- apply(profile_raw, 2L, stats::mad)
        bad_scale <- !is.finite(raw_scale) |
          raw_scale <= sqrt(.Machine$double.eps)
        if (any(bad_scale)) {
          raw_scale[bad_scale] <- apply(
            profile_raw[, bad_scale, drop = FALSE],
            2L,
            stats::sd
          )
        }
        positive_scale <- raw_scale[
          is.finite(raw_scale) & raw_scale > sqrt(.Machine$double.eps)
        ]
        reference_scale <- if (length(positive_scale)) {
          stats::median(positive_scale)
        } else {
          1
        }
        raw_scale[!is.finite(raw_scale) |
                    raw_scale <= sqrt(.Machine$double.eps)] <- reference_scale
        relative_scale <- raw_scale / reference_scale
        relative_scale <- pmin(pmax(relative_scale, 0.25), 4)
        scale <- reference_scale * relative_scale
      } else {
        center <- rep(0, M)
        scale <- rep(1, M)
      }
    }

    if (length(center) != M || length(scale) != M ||
        any(!is.finite(center)) || any(!is.finite(scale)) ||
        any(scale <= 0)) {
      stop("Invalid training-fold center or scale for the profile.")
    }

    profile <- if (standardize) {
      sweep(sweep(profile_raw, 2L, center, "-"), 2L, scale, "/")
    } else {
      profile_raw
    }

    list(
      data = profile,
      raw_data = profile_raw,
      center = center,
      scale = scale,
      M = M
    )
  }

  ps_squared_distances <- function(A, centers) {
    distances <- outer(rowSums(A^2), rowSums(centers^2), "+") -
      2 * tcrossprod(A, centers)
    pmax(distances / ncol(A), 0)
  }

  ps_balanced_capacities <- function(n, G) {
    capacities <- rep(n %/% G, G)
    remainder <- n %% G
    if (remainder > 0L) {
      capacities[seq_len(remainder)] <- capacities[seq_len(remainder)] + 1L
    }
    capacities
  }

  ps_warned_solver <- FALSE
  ps_balanced_assignment <- function(distance_matrix, capacities) {
    n <- nrow(distance_matrix)
    G <- ncol(distance_matrix)
    if (length(capacities) != G || sum(capacities) != n) {
      stop("PS capacities do not match the assignment problem.")
    }

    if (requireNamespace("lpSolve", quietly = TRUE)) {
      solution <- lpSolve::lp.transport(
        cost.mat = distance_matrix,
        direction = "min",
        row.signs = rep("=", n),
        row.rhs = rep(1, n),
        col.signs = rep("=", G),
        col.rhs = capacities
      )
      if (solution$status != 0L) {
        stop("The balanced PS transportation problem did not solve successfully.")
      }
      return(max.col(solution$solution, ties.method = "first"))
    }

    if (requireNamespace("clue", quietly = TRUE)) {
      group_slots <- rep(seq_len(G), capacities)
      slot_cost <- distance_matrix[, group_slots, drop = FALSE]
      assigned_slot <- as.integer(clue::solve_LSAP(slot_cost, maximum = FALSE))
      return(group_slots[assigned_slot])
    }

    if (!ps_warned_solver) {
      warning(
        "Neither lpSolve nor clue is installed; using a greedy balanced PS assignment.",
        call. = FALSE
      )
      ps_warned_solver <<- TRUE
    }

    group <- integer(n)
    remaining <- capacities
    ordered_pairs <- order(distance_matrix)
    for (pair in ordered_pairs) {
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

  ps_kmeans_plus_plus <- function(A, G) {
    n <- nrow(A)
    chosen <- integer(G)
    chosen[1L] <- sample.int(n, 1L)
    closest <- rowSums((A - A[chosen[1L], ])^2)

    if (G >= 2L) {
      for (g in 2:G) {
        available <- setdiff(seq_len(n), chosen[seq_len(g - 1L)])
        probabilities <- closest
        probabilities[-available] <- 0
        if (sum(probabilities) <= 0) {
          chosen[g] <- sample(available, 1L)
        } else {
          chosen[g] <- sample.int(n, 1L, prob = probabilities)
        }
        new_distance <- rowSums((A - A[chosen[g], ])^2)
        closest <- pmin(closest, new_distance)
      }
    }
    A[chosen, , drop = FALSE]
  }

  ps_update_centers <- function(A, group, G) {
    centers <- matrix(NA_real_, nrow = G, ncol = ncol(A))
    for (g in seq_len(G)) {
      centers[g, ] <- colMeans(A[group == g, , drop = FALSE])
    }
    centers
  }

  ps_fit_balanced <- function(A, G, nstart = ps_nstart, seed = ps_seed) {
    A <- as.matrix(A)
    storage.mode(A) <- "double"
    n <- nrow(A)
    if (any(!is.finite(A))) {
      stop("The PS profile contains non-finite values.")
    }
    if (G < 1L || G > n) {
      stop("The PS group count must be between 1 and the training-fold size.")
    }
    if (length(nstart) != 1L || !is.finite(nstart) ||
        nstart < 1 || nstart != floor(nstart)) {
      stop("The PS number of starts must be one positive integer.")
    }
    nstart <- as.integer(nstart)

    if (!is.null(seed)) {
      if (length(seed) != 1L || !is.finite(seed)) {
        stop("The PS fit seed must be NULL or one finite number.")
      }
      seed_existed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
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
      set.seed(as.integer(seed))
    }

    capacities <- ps_balanced_capacities(n, G)
    fits <- vector("list", nstart)
    for (start in seq_len(nstart)) {
      centers <- ps_kmeans_plus_plus(A, G)
      previous_group <- rep(NA_integer_, n)
      previous_objective <- Inf
      converged <- FALSE

      for (iteration in seq_len(ps_max_iter)) {
        distances <- ps_squared_distances(A, centers)
        group <- ps_balanced_assignment(distances, capacities)
        centers <- ps_update_centers(A, group, G)
        updated_distances <- ps_squared_distances(A, centers)
        objective <- mean(updated_distances[cbind(seq_len(n), group)])
        relative_change <- if (is.finite(previous_objective)) {
          abs(previous_objective - objective) /
            max(abs(previous_objective), .Machine$double.eps)
        } else {
          Inf
        }

        if (identical(group, previous_group) || relative_change < ps_tol) {
          converged <- TRUE
          break
        }
        previous_group <- group
        previous_objective <- objective
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
    best <- fits[[which.min(objectives)]]
    best$G <- G
    best$capacities <- capacities
    best$group_sizes <- tabulate(best$group, nbins = G)
    best$objectives_by_start <- objectives
    best
  }

  ps_seed_for_fit <- function(repeat_id, fold_id, grid_id) {
    if (is.null(ps_seed)) return(NULL)
    modulus <- .Machine$integer.max
    as.integer((as.double(ps_seed) + 104729 * repeat_id +
                  1009 * fold_id + 37 * grid_id) %% modulus)
  }

  ps_candidate_grid <- function(n_train, n_test) {
    validation_limit <- floor(n_train / ps_cv_folds)
    max_G <- min(
      n_train,
      n_test,
      validation_limit,
      floor(n_train / ps_min_group_size),
      floor(n_test / ps_min_group_size)
    )
    if (max_G < 1L) {
      stop("The PS folds are too small for ps_min_group_size and ps_cv_folds.")
    }

    if (!is.null(ps_G_grid)) {
      invalid <- ps_G_grid > max_G
      if (any(invalid)) {
        stop(
          "ps_G_grid contains values above the feasible maximum ", max_G,
          " for this outer fold. Reduce ps_G_grid, ps_cv_folds, or ",
          "ps_min_group_size."
        )
      }
      return(ps_G_grid)
    }

    if (max_G <= 12L) return(seq_len(max_G))
    tail_grid <- unique(as.integer(round(exp(seq(
      log(9), log(max_G), length.out = 12L
    )))))
    sort(unique(c(seq_len(8L), tail_grid, max_G)))
  }

  ps_select_G <- function(A, n_test) {
    A <- as.matrix(A)
    n <- nrow(A)
    if (ps_cv_folds > n) {
      stop("ps_cv_folds cannot exceed the outer training-fold size.")
    }
    grid <- ps_candidate_grid(n, n_test)

    if (ps_G_rule == "n-third") {
      target <- max(1L, ceiling(n^(1 / 3)))
      selected <- grid[which.min(abs(grid - target))]
      return(list(
        selected_G = as.integer(selected),
        rule = "n-third",
        path = data.frame(G = grid, cv_loss = NA_real_, cv_se = NA_real_),
        threshold = NA_real_,
        minimizer_G = NA_integer_
      ))
    }

    seed_existed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (!is.null(ps_seed)) {
      if (seed_existed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
      on.exit({
        if (seed_existed) {
          assign(".Random.seed", old_seed, envir = .GlobalEnv)
        } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
          rm(".Random.seed", envir = .GlobalEnv)
        }
      }, add = TRUE)
      set.seed(as.integer((as.double(ps_seed) + 7919) %% .Machine$integer.max))
    }

    loss_sum <- matrix(0, nrow = n, ncol = length(grid))
    loss_count <- matrix(0L, nrow = n, ncol = length(grid))
    for (repeat_id in seq_len(ps_cv_repeats)) {
      permutation <- sample.int(n)
      fold_id <- integer(n)
      fold_id[permutation] <- rep(seq_len(ps_cv_folds), length.out = n)

      for (fold in seq_len(ps_cv_folds)) {
        validation_index <- which(fold_id == fold)
        fitting_index <- which(fold_id != fold)
        for (grid_id in seq_along(grid)) {
          G <- grid[grid_id]
          fit <- ps_fit_balanced(
            A[fitting_index, , drop = FALSE],
            G,
            nstart = ps_select_nstart,
            seed = ps_seed_for_fit(repeat_id, fold, grid_id)
          )
          distances <- ps_squared_distances(
            A[validation_index, , drop = FALSE], fit$centers
          )
          capacities <- ps_balanced_capacities(length(validation_index), G)
          group <- ps_balanced_assignment(distances, capacities)
          losses <- distances[cbind(seq_along(validation_index), group)]
          loss_sum[validation_index, grid_id] <-
            loss_sum[validation_index, grid_id] + losses
          loss_count[validation_index, grid_id] <-
            loss_count[validation_index, grid_id] + 1L
        }
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
      selected_G = as.integer(grid[selected_index]),
      rule = ps_G_rule,
      path = data.frame(G = grid, cv_loss = cv_loss, cv_se = cv_se),
      threshold = threshold,
      minimizer_G = as.integer(grid[minimizer_index])
    )
  }

  hac_seed_value <- function(repeat_id = 0L, fold_id = 0L, stage_id = 0L) {
    if (is.null(hac_seed)) return(NULL)
    modulus <- .Machine$integer.max
    as.integer((as.double(hac_seed) + 130363 * repeat_id +
                  2017 * fold_id + 53 * stage_id) %% modulus)
  }

  hac_split_reference <- function(n, seed = NULL) {
    if (n < 2L) stop("HAC new needs at least two auxiliary units before the reference split.")
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
    permutation <- sample.int(n)
    n_reference <- max(1L, min(n - 1L, round(n * hac_reference_fraction)))
    list(
      reference = permutation[seq_len(n_reference)],
      auxiliary = permutation[-seq_len(n_reference)]
    )
  }

  hac_projection_scores <- function(A, reference) {
    tcrossprod(A, reference) / ncol(A)
  }

  hac_score_distances <- function(scores, center_scores) {
    distances <- outer(rowSums(scores^2), rowSums(center_scores^2), "+") -
      2 * tcrossprod(scores, center_scores)
    pmax(distances / ncol(scores), 0)
  }

  hac_fit_quadratic <- function(auxiliary, reference, G) {
    auxiliary <- as.matrix(auxiliary)
    reference <- as.matrix(reference)
    n_auxiliary <- nrow(auxiliary)
    if (G < 1L || G > n_auxiliary) {
      stop("The HAC group count must be between one and the auxiliary-HAC sample size.")
    }
    scores <- hac_projection_scores(auxiliary, reference)
    if (n_auxiliary == 1L || G == 1L) {
      group <- rep(1L, n_auxiliary)
      tree <- if (n_auxiliary == 1L) NULL else
        stats::hclust(stats::as.dist(as.matrix(stats::dist(scores))^2 /
                                       nrow(reference)), method = link)
    } else {
      q_matrix <- as.matrix(stats::dist(scores))^2 / nrow(reference)
      tree <- stats::hclust(stats::as.dist(q_matrix), method = link)
      group <- stats::cutree(tree, k = G)
    }
    q_matrix <- if (n_auxiliary == 1L) matrix(0, 1L, 1L) else
      as.matrix(stats::dist(scores))^2 / nrow(reference)
    medoid_index <- vapply(seq_len(G), function(g) {
      members <- which(group == g)
      members[which.min(rowSums(q_matrix[members, members, drop = FALSE]))]
    }, integer(1L))
    merge_height <- if (is.null(tree) || G == n_auxiliary) {
      0
    } else {
      tree$height[n_auxiliary - G]
    }
    list(
      G = G, group = as.integer(group), tree = tree,
      medoid_index = medoid_index, medoid_profile = auxiliary[medoid_index, , drop = FALSE],
      medoid_scores = scores[medoid_index, , drop = FALSE],
      reference = reference, merge_height = merge_height,
      group_sizes = tabulate(group, nbins = G)
    )
  }

  hac_assign_quadratic <- function(profile, fit, balanced = TRUE) {
    scores <- hac_projection_scores(profile, fit$reference)
    distances <- hac_score_distances(scores, fit$medoid_scores)
    if (balanced) {
      capacities <- ps_balanced_capacities(nrow(profile), fit$G)
      group <- ps_balanced_assignment(distances, capacities)
    } else {
      group <- max.col(-distances, ties.method = "first")
    }
    list(
      group = as.integer(group),
      distance = distances[cbind(seq_len(nrow(profile)), group)],
      distances = distances
    )
  }

  hac_grid_from_max <- function(max_G) {
    if (max_G < 1L) stop("No feasible HAC new group count remains.")
    if (!is.null(hac_G_grid)) {
      if (any(hac_G_grid > max_G)) {
        stop("hac_G_grid exceeds the feasible maximum ", max_G,
             "; reduce the grid, hac_cv_folds, or hac_min_group_size.")
      }
      return(hac_G_grid)
    }
    if (max_G <= 12L) return(seq_len(max_G))
    tail_grid <- unique(as.integer(round(exp(seq(
      log(9), log(max_G), length.out = 12L
    )))))
    sort(unique(c(seq_len(8L), tail_grid, max_G)))
  }

  hac_select_G <- function(profile, n_test) {
    profile <- as.matrix(profile)
    n <- nrow(profile)
    if (hac_cv_folds > n) stop("hac_cv_folds cannot exceed the outer training-fold size.")

    seed_existed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (!is.null(hac_seed)) {
      if (seed_existed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
      on.exit({
        if (seed_existed) {
          assign(".Random.seed", old_seed, envir = .GlobalEnv)
        } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
          rm(".Random.seed", envir = .GlobalEnv)
        }
      }, add = TRUE)
      set.seed(hac_seed_value(stage_id = 701L))
    }

    designs <- list()
    design_id <- 0L
    for (repeat_id in seq_len(hac_cv_repeats)) {
      permutation <- sample.int(n)
      fold_label <- integer(n)
      fold_label[permutation] <- rep(seq_len(hac_cv_folds), length.out = n)
      for (fold in seq_len(hac_cv_folds)) {
        validation <- which(fold_label == fold)
        fitting <- which(fold_label != fold)
        split <- hac_split_reference(
          length(fitting), hac_seed_value(repeat_id, fold, 809L)
        )
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
      floor(n_test / hac_min_group_size)
    )
    grid <- hac_grid_from_max(max_G)
    loss_sum <- matrix(0, nrow = n, ncol = length(grid))
    loss_count <- matrix(0L, nrow = n, ncol = length(grid))
    for (design in designs) {
      for (grid_id in seq_along(grid)) {
        fit <- hac_fit_quadratic(
          profile[design$auxiliary, , drop = FALSE],
          profile[design$reference, , drop = FALSE],
          grid[grid_id]
        )
        assigned <- hac_assign_quadratic(
          profile[design$validation, , drop = FALSE], fit, balanced = TRUE
        )
        loss_sum[design$validation, grid_id] <-
          loss_sum[design$validation, grid_id] + assigned$distance
        loss_count[design$validation, grid_id] <-
          loss_count[design$validation, grid_id] + 1L
      }
    }
    if (any(loss_count == 0L)) stop("Internal HAC validation loss was not filled.")
    unit_loss <- loss_sum / loss_count
    cv_loss <- colMeans(unit_loss)
    cv_se <- apply(unit_loss, 2L, stats::sd) / sqrt(n)
    cv_se[!is.finite(cv_se)] <- 0
    minimizer_index <- which.min(cv_loss)
    threshold <- cv_loss[minimizer_index] + cv_se[minimizer_index]
    selected_index <- if (hac_G_rule == "minimum") minimizer_index else
      which(cv_loss <= threshold)[1L]
    list(
      selected_G = as.integer(grid[selected_index]),
      rule = hac_G_rule,
      path = data.frame(G = grid, cv_loss = cv_loss, cv_se = cv_se),
      threshold = threshold,
      minimizer_G = as.integer(grid[minimizer_index])
    )
  }

  estimate_time_fold <- function(train_t_idx, test_t_idx, outer_fold_id) {
    y_train <- y[, train_t_idx, drop = FALSE]
    X_train <- X[, train_t_idx, , drop = FALSE]
    T_train <- length(train_t_idx)
    T_test <- length(test_t_idx)
    kmedoid_diagnostics <- list(
      profile_dimension = NA_integer_, objective = NA_real_,
      group_sizes = NULL, medoid_index = NULL, medoid_unit_id = NULL,
      solver = NULL, global_optimum = NA, converged = NA,
      iterations = NA_integer_, selection_rule = NULL,
      selection_path = NULL, distance = NULL
    )
    pk_diagnostics <- list(
      input_dimension = NA_integer_, feature_dimension = NA_integer_,
      objective = NA_real_, group_sizes = NULL, converged = NA,
      iterations = NA_integer_, selection_rule = NULL,
      selection_path = NULL, selection_threshold = NA_real_,
      noise_floor = NA_real_
    )

    if (pre_cluster) {
      Du <- Du_pre
      if (is.null(Du) || !is.matrix(Du) || nrow(Du) != N || ncol(Du) < 1) {
        stop("When pre_cluster = TRUE, Du_pre must be an N x G_unit matrix with G_unit >= 1.")
      }
      G_unit <- ncol(Du)
    } else if (cluster_method == "kmeans") {
      X_list_train <- lapply(seq_len(dim(X_train)[3]), function(k) X_train[, , k])
      if (is.null(unit_cluster)) {
        clusteri <- cluster_general(y_train, X_list_train, N, T_train, init = 50,
                                    type = cf_long_type, gamma = gamma, cc = cc, dim_moment = dim_moment)
      } else {
        clusteri <- cluster_general(y_train, X_list_train, N, T_train, init = 50,
                                    type = cf_long_type, groups = c(floor(unit_cluster)),
                                    gamma = gamma, cc = cc, dim_moment = dim_moment)
      }
      G_unit <- clusteri$clusters
      klong <- clusteri$res
      Du <- matrix(0, N, G_unit)
      for (j in seq_len(G_unit)) {
        Du[, j] <- as.numeric(klong$cluster == j)
      }
    } else if (cluster_method == "hierarchical") {
      X_norm_train <- array(NA_real_, dim = dim(X_train))
      for (k in seq_len(dim(X_train)[3])) {
        X_norm_train[, , k] <- X_train[, , k]
      }
      y_norm_train <- y_train

      if (is.null(unit_cluster)) {
        res_unit <- cluster_Hierarchical(y_norm_train, X_norm_train, link = link,
                                         type = "unit", method_auto = method_auto,
                                         dist_type = dist_type)
      } else {
        res_unit <- cluster_Hierarchical(y_norm_train, X_norm_train, link = link,
                                         type = "unit", cluster = unit_cluster,
                                         dist_type = dist_type)
      }
      Du <- res_unit$indicator
      G_unit <- res_unit$G
    } else if (cluster_method == "kmedoid") {
      fold_options <- kmedoid_options
      if (!is.null(fold_options$seed)) {
        fold_options$seed <- as.integer(
          (as.double(fold_options$seed) + 1009 * outer_fold_id) %%
            .Machine$integer.max
        )
      }
      profile <- .kmedoid_profile(
        y_train, X_train, standardize = fold_options$standardize
      )
      fitted <- .kmedoid_fit_profile(profile, unit_cluster, fold_options)
      medoid_fit <- fitted$fit
      G_unit <- medoid_fit$G
      unit_group <- medoid_fit$group
      Du <- sapply(seq_len(G_unit), function(g) as.numeric(unit_group == g))
      if (G_unit == 1L) Du <- matrix(Du, ncol = 1L)
      kmedoid_diagnostics <- list(
        profile_dimension = profile$M,
        objective = medoid_fit$objective,
        group_sizes = medoid_fit$group_sizes,
        medoid_index = medoid_fit$medoids,
        medoid_unit_id = ids[medoid_fit$medoids],
        solver = medoid_fit$solver,
        global_optimum = medoid_fit$global_optimum,
        converged = medoid_fit$converged,
        iterations = medoid_fit$iterations,
        selection_rule = fitted$selection$rule,
        selection_path = fitted$selection$path,
        distance = "quadratic-projection"
      )
    } else if (cluster_method == "projection-kmeans") {
      fold_options <- pk_options
      if (!is.null(fold_options$seed)) {
        fold_options$seed <- as.integer(
          (as.double(fold_options$seed) + 1009 * outer_fold_id) %%
            .Machine$integer.max
        )
      }
      feature <- .pk_feature_object(
        y_train, X_train, standardize = fold_options$standardize
      )
      fitted <- .pk_fit_features(
        feature, unit_cluster, fold_options, fold_id = outer_fold_id,
        evaluation_T = T_test
      )
      projection_fit <- fitted$fit
      G_unit <- projection_fit$G
      unit_group <- projection_fit$group
      Du <- sapply(seq_len(G_unit), function(g) {
        as.numeric(unit_group == g)
      })
      if (G_unit == 1L) Du <- matrix(Du, ncol = 1L)
      pk_diagnostics <- list(
        input_dimension = feature$M,
        feature_dimension = feature$N,
        objective = projection_fit$objective,
        group_sizes = projection_fit$group_sizes,
        converged = projection_fit$converged,
        iterations = projection_fit$iterations,
        selection_rule = fitted$selection$rule,
        selection_path = fitted$selection$path,
        selection_threshold = fitted$selection$threshold,
        reference_G = fitted$selection$reference_G,
        reference_objective = fitted$selection$reference_objective,
        inference_multiplier = fitted$selection$inference_multiplier,
        complexity_threshold = fitted$selection$complexity_threshold,
        log_N_rule = fitted$selection$log_N_rule,
        log_N_value = fitted$selection$log_N_value,
        requested_max_G = fitted$selection$requested_max_G,
        effective_max_G = fitted$selection$effective_max_G,
        theoretical_max_G = fitted$selection$theoretical_max_G,
        max_G_source = fitted$selection$max_G_source,
        selection_nstart = fitted$selection$screening_nstart,
        selection_max_iter = fitted$selection$screening_max_iter,
        screening_objective = fitted$selection$screening_objective,
        final_objective = fitted$selection$final_objective,
        noise_floor = fitted$selection$noise_floor
      )
    } else {
      stop("Unsupported cluster_method for time cross-fitting. Use 'kmeans', 'hierarchical', 'kmedoid', or 'projection-kmeans'.")
    }

    Diagu <- matrix(0, G_unit, G_unit)
    for (j in seq_len(G_unit)) {
      Diagu[j, j] <- 1 / sum(Du[, j])
    }
    unit_group <- max.col(Du)

    y_test <- y[, test_t_idx, drop = FALSE]
    X_test <- X[, test_t_idx, , drop = FALSE]

    Mu <- diag(N) - Du %*% Diagu %*% t(Du)

    Z <- array(NA_real_, dim = c(N, T_test, K + 1))
    Z[, , 1] <- y_test
    for (k in seq_len(K)) {
      Z[, , k + 1] <- X_test[, , k]
    }

    Z_proj <- array(NA_real_, dim = c(N, T_test, K + 1))
    for (k in seq_len(K + 1)) {
      Z_proj[, , k] <- Mu %*% Z[, , k]
    }

    y_trans <- Z_proj[, , 1]
    X_trans <- Z_proj[, , 2:(K + 1)]

    tY_vector <- as.vector(y_trans)
    if (T_test == 1) {
      X_trans_mat <- if (length(dim(X_trans)) == 3) X_trans[, 1, , drop = TRUE] else X_trans
      tX_matrix <- X_trans_mat
    } else if (length(dim(X_trans)) == 2) {
      tX_matrix <- matrix(as.vector(t(X_trans)), nrow = N * T_test, ncol = K)
    } else {
      tX_matrix <- matrix(aperm(X_trans, c(1, 2, 3)), nrow = N * T_test, ncol = K)
    }

    d_vector <- tX_matrix[, 1]
    x_matrix <- tX_matrix[, -1, drop = FALSE]

    list(
      tY_vector = tY_vector,
      d_vector = d_vector,
      x_matrix = x_matrix,
      T_test = T_test,
      G_unit = G_unit,
      unit_group = unit_group,
      kmedoid_profile_dimension = kmedoid_diagnostics$profile_dimension,
      kmedoid_objective = kmedoid_diagnostics$objective,
      kmedoid_group_sizes = kmedoid_diagnostics$group_sizes,
      kmedoid_medoid_index = kmedoid_diagnostics$medoid_index,
      kmedoid_medoid_unit_id = kmedoid_diagnostics$medoid_unit_id,
      kmedoid_solver = kmedoid_diagnostics$solver,
      kmedoid_global_optimum = kmedoid_diagnostics$global_optimum,
      kmedoid_converged = kmedoid_diagnostics$converged,
      kmedoid_iterations = kmedoid_diagnostics$iterations,
      kmedoid_selection_rule = kmedoid_diagnostics$selection_rule,
      kmedoid_selection_path = kmedoid_diagnostics$selection_path,
      kmedoid_distance = kmedoid_diagnostics$distance,
      pk_input_dimension = pk_diagnostics$input_dimension,
      pk_feature_dimension = pk_diagnostics$feature_dimension,
      pk_objective = pk_diagnostics$objective,
      pk_group_sizes = pk_diagnostics$group_sizes,
      pk_converged = pk_diagnostics$converged,
      pk_iterations = pk_diagnostics$iterations,
      pk_selection_rule = pk_diagnostics$selection_rule,
      pk_selection_path = pk_diagnostics$selection_path,
      pk_selection_threshold = pk_diagnostics$selection_threshold,
      pk_reference_G = pk_diagnostics$reference_G,
      pk_reference_objective = pk_diagnostics$reference_objective,
      pk_inference_multiplier = pk_diagnostics$inference_multiplier,
      pk_complexity_threshold = pk_diagnostics$complexity_threshold,
      pk_log_N_rule = pk_diagnostics$log_N_rule,
      pk_log_N_value = pk_diagnostics$log_N_value,
      pk_requested_max_G = pk_diagnostics$requested_max_G,
      pk_effective_max_G = pk_diagnostics$effective_max_G,
      pk_theoretical_max_G = pk_diagnostics$theoretical_max_G,
      pk_max_G_source = pk_diagnostics$max_G_source,
      pk_selection_nstart = pk_diagnostics$selection_nstart,
      pk_selection_max_iter = pk_diagnostics$selection_max_iter,
      pk_screening_objective = pk_diagnostics$screening_objective,
      pk_final_objective = pk_diagnostics$final_objective,
      pk_noise_floor = pk_diagnostics$noise_floor
    )
  }

  estimate_unit_fold <- function(train_i_idx, test_i_idx, outer_fold_id) {
    if (pre_cluster) {
      stop("cross_type = 'unit' estimates held-out memberships from training-fold centers; use pre_cluster = FALSE.")
    }

    N_train <- length(train_i_idx)
    N_test <- length(test_i_idx)
    y_train <- y[train_i_idx, , drop = FALSE]
    X_train <- X[train_i_idx, , , drop = FALSE]
    train_features <- unit_feature_data(y_train, X_train)
    test_features <- unit_feature_data(y[test_i_idx, , drop = FALSE],
                                       X[test_i_idx, , , drop = FALSE],
                                       center = train_features$center,
                                       scale = train_features$scale)
    hac_diagnostics <- list(
      profile_dimension = NA_integer_, reference_size = NA_integer_,
      auxiliary_size = NA_integer_, merge_height = NA_real_,
      evaluation_coverage = NA_real_, effective_cut = NA_real_,
      stochastic_rate = NA_real_,
      auxiliary_group_sizes = NULL, test_group_sizes = NULL,
      medoid_auxiliary_index = NULL, medoid_unit_id = NULL,
      reference_unit_id = NULL, auxiliary_unit_id = NULL,
      selection_rule = NULL,
      selection_path = NULL, selection_threshold = NA_real_,
      selection_minimizer_G = NA_integer_
    )
    ps_diagnostics <- list(
      profile_dimension = NA_integer_,
      train_objective = NA_real_,
      train_group_sizes = NULL,
      test_group_sizes = NULL,
      converged = NA,
      iterations = NA_integer_,
      selection_rule = NULL,
      selection_path = NULL,
      selection_threshold = NA_real_,
      selection_minimizer_G = NA_integer_
    )
    test_cluster <- NULL

    if (cluster_method == "kmeans") {
      X_list_train <- lapply(seq_len(dim(X_train)[3]), function(k) X_train[, , k])

      if (is.null(unit_cluster)) {
        clusteri <- cluster_general(y_train, X_list_train, N_train, T, init = 50,
                                    type = cf_long_type, gamma = gamma, cc = cc, dim_moment = dim_moment)
      } else {
        clusteri <- cluster_general(y_train, X_list_train, N_train, T, init = 50,
                                    type = cf_long_type, groups = c(floor(unit_cluster)),
                                    gamma = gamma, cc = cc, dim_moment = dim_moment)
      }

      G_unit <- clusteri$clusters
      train_centers <- clusteri$res$centers
    } else if (cluster_method == "hierarchical") {
      if (is.null(unit_cluster)) {
        res_unit <- cluster_Hierarchical(y_train, X_train, link = link,
                                         type = "unit", method_auto = method_auto,
                                         dist_type = dist_type)
      } else {
        res_unit <- cluster_Hierarchical(y_train, X_train, link = link,
                                         type = "unit", cluster = unit_cluster,
                                         dist_type = dist_type)
      }

      G_unit <- res_unit$G
      train_centers <- do.call(rbind, lapply(seq_len(G_unit), function(j) {
        colMeans(train_features$data[res_unit$clusters == j, , drop = FALSE])
      }))
    } else if (cluster_method == "ps-kmeans") {
      if (N_train < 2L || N_test < 2L) {
        stop("Unit-cross-fitted PS k-means requires at least two units in each fold.")
      }

      train_profile <- profile_data(y_train, X_train)
      test_profile <- profile_data(
        y[test_i_idx, , drop = FALSE],
        X[test_i_idx, , , drop = FALSE],
        center = train_profile$center,
        scale = train_profile$scale
      )

      max_G <- min(
        N_train, N_test,
        floor(N_train / ps_min_group_size),
        floor(N_test / ps_min_group_size)
      )
      if (max_G < 1L) {
        stop("The outer unit folds are too small for ps_min_group_size.")
      }

      if (is.null(unit_cluster)) {
        ps_selection <- ps_select_G(train_profile$data, N_test)
        G_unit <- ps_selection$selected_G
      } else {
        G_unit <- floor(unit_cluster)
        ps_selection <- list(
          rule = "fixed", path = NULL, threshold = NA_real_,
          minimizer_G = NA_integer_
        )
      }
      if (length(G_unit) != 1L || !is.finite(G_unit) ||
          G_unit < 1L || G_unit > max_G) {
        stop(
          "For unit-cross-fitted PS k-means, unit_cluster must be an integer ",
          "between 1 and the feasible maximum ", max_G,
          " implied by the fold sizes and ps_min_group_size."
        )
      }
      G_unit <- as.integer(G_unit)

      ps_fit <- ps_fit_balanced(train_profile$data, G_unit)
      train_centers <- ps_fit$centers

      # Held-out memberships use only the transformation and centers learned
      # on the training units. The capacity constraint also keeps every
      # held-out group nonempty and approximately equal in size.
      test_distances <- ps_squared_distances(test_profile$data, train_centers)
      test_capacities <- ps_balanced_capacities(N_test, G_unit)
      test_cluster <- ps_balanced_assignment(test_distances, test_capacities)

      ps_diagnostics <- list(
        profile_dimension = train_profile$M,
        train_objective = ps_fit$objective,
        train_group_sizes = ps_fit$group_sizes,
        test_group_sizes = tabulate(test_cluster, nbins = G_unit),
        converged = ps_fit$converged,
        iterations = ps_fit$iterations,
        selection_rule = ps_selection$rule,
        selection_path = ps_selection$path,
        selection_threshold = ps_selection$threshold,
        selection_minimizer_G = ps_selection$minimizer_G
      )
    } else if (cluster_method == "HAC new") {
      if (N_train < 4L || N_test < 2L) {
        stop("HAC new requires at least four training units and two held-out units per outer fold.")
      }
      train_profile <- profile_data(
        y_train, X_train, standardize = hac_standardize
      )
      test_profile <- profile_data(
        y[test_i_idx, , drop = FALSE],
        X[test_i_idx, , , drop = FALSE],
        center = train_profile$center, scale = train_profile$scale,
        standardize = hac_standardize
      )
      outer_split <- hac_split_reference(
        N_train, hac_seed_value(fold_id = outer_fold_id, stage_id = 907L)
      )
      max_G <- min(
        length(outer_split$auxiliary),
        floor(N_test / hac_min_group_size)
      )
      if (is.null(unit_cluster)) {
        hac_selection <- hac_select_G(train_profile$data, N_test)
        G_unit <- hac_selection$selected_G
      } else {
        G_unit <- floor(unit_cluster)
        hac_selection <- list(
          rule = "fixed", path = NULL, threshold = NA_real_,
          minimizer_G = NA_integer_
        )
      }
      if (length(G_unit) != 1L || !is.finite(G_unit) ||
          G_unit < 1L || G_unit > max_G) {
        stop("For HAC new, unit_cluster must be an integer between 1 and ",
             max_G, " for this outer fold.")
      }
      G_unit <- as.integer(G_unit)
      hac_fit <- hac_fit_quadratic(
        train_profile$data[outer_split$auxiliary, , drop = FALSE],
        train_profile$data[outer_split$reference, , drop = FALSE],
        G_unit
      )
      assigned <- hac_assign_quadratic(test_profile$data, hac_fit, balanced = TRUE)
      test_cluster <- assigned$group
      evaluation_coverage <- max(assigned$distance)
      hac_diagnostics <- list(
        profile_dimension = train_profile$M,
        reference_size = length(outer_split$reference),
        auxiliary_size = length(outer_split$auxiliary),
        merge_height = hac_fit$merge_height,
        evaluation_coverage = evaluation_coverage,
        effective_cut = max(hac_fit$merge_height, evaluation_coverage),
        stochastic_rate = log(N) / train_profile$M +
          log(N) / length(outer_split$reference),
        auxiliary_group_sizes = hac_fit$group_sizes,
        test_group_sizes = tabulate(test_cluster, nbins = G_unit),
        medoid_auxiliary_index = outer_split$auxiliary[hac_fit$medoid_index],
        medoid_unit_id = ids[
          train_i_idx[outer_split$auxiliary[hac_fit$medoid_index]]
        ],
        reference_unit_id = ids[train_i_idx[outer_split$reference]],
        auxiliary_unit_id = ids[train_i_idx[outer_split$auxiliary]],
        selection_rule = hac_selection$rule,
        selection_path = hac_selection$path,
        selection_threshold = hac_selection$threshold,
        selection_minimizer_G = hac_selection$minimizer_G
      )
    } else {
      stop("Unsupported cluster_method. Use 'kmeans', 'hierarchical', 'ps-kmeans', or 'HAC new'.")
    }

    if (is.null(test_cluster)) {
      test_cluster <- nearest_center(test_features$data, train_centers)
    }
    Du <- matrix(0, N_test, G_unit)
    for (j in seq_len(G_unit)) {
      Du[, j] <- as.numeric(test_cluster == j)
    }

    used_clusters <- colSums(Du) > 0
    Du_used <- Du[, used_clusters, drop = FALSE]
    G_unit_used <- ncol(Du_used)
    Diagu <- matrix(0, G_unit_used, G_unit_used)
    for (j in seq_len(G_unit_used)) {
      Diagu[j, j] <- 1 / sum(Du_used[, j])
    }

    Mu <- diag(N_test) - Du_used %*% Diagu %*% t(Du_used)

    Z <- array(NA_real_, dim = c(N_test, T, K + 1))
    Z[, , 1] <- y[test_i_idx, , drop = FALSE]
    for (k in seq_len(K)) {
      Z[, , k + 1] <- X[test_i_idx, , k]
    }

    Z_proj <- array(NA_real_, dim = c(N_test, T, K + 1))
    for (k in seq_len(K + 1)) {
      Z_proj[, , k] <- Mu %*% Z[, , k]
    }

    y_trans <- Z_proj[, , 1]
    X_trans <- Z_proj[, , 2:(K + 1)]

    tY_vector <- as.vector(y_trans)
    if (length(dim(X_trans)) == 2) {
      tX_matrix <- matrix(as.vector(t(X_trans)), nrow = N_test * T, ncol = K)
    } else {
      tX_matrix <- matrix(aperm(X_trans, c(1, 2, 3)), nrow = N_test * T, ncol = K)
    }

    list(
      tY_vector = tY_vector,
      d_vector = tX_matrix[, 1],
      x_matrix = tX_matrix[, -1, drop = FALSE],
      N_test = N_test,
      G_unit = G_unit,
      G_unit_used = G_unit_used,
      unit_group = test_cluster,
      ps_profile_dimension = ps_diagnostics$profile_dimension,
      ps_objective = ps_diagnostics$train_objective,
      ps_train_group_sizes = ps_diagnostics$train_group_sizes,
      ps_test_group_sizes = ps_diagnostics$test_group_sizes,
      ps_converged = ps_diagnostics$converged,
      ps_iterations = ps_diagnostics$iterations,
      ps_selection_rule = ps_diagnostics$selection_rule,
      ps_selection_path = ps_diagnostics$selection_path,
      ps_selection_threshold = ps_diagnostics$selection_threshold,
      ps_selection_minimizer_G = ps_diagnostics$selection_minimizer_G,
      hac_profile_dimension = hac_diagnostics$profile_dimension,
      hac_reference_size = hac_diagnostics$reference_size,
      hac_auxiliary_size = hac_diagnostics$auxiliary_size,
      hac_merge_height = hac_diagnostics$merge_height,
      hac_evaluation_coverage = hac_diagnostics$evaluation_coverage,
      hac_effective_cut = hac_diagnostics$effective_cut,
      hac_stochastic_rate = hac_diagnostics$stochastic_rate,
      hac_auxiliary_group_sizes = hac_diagnostics$auxiliary_group_sizes,
      hac_test_group_sizes = hac_diagnostics$test_group_sizes,
      hac_medoid_auxiliary_index = hac_diagnostics$medoid_auxiliary_index,
      hac_medoid_unit_id = hac_diagnostics$medoid_unit_id,
      hac_reference_unit_id = hac_diagnostics$reference_unit_id,
      hac_auxiliary_unit_id = hac_diagnostics$auxiliary_unit_id,
      hac_selection_rule = hac_diagnostics$selection_rule,
      hac_selection_path = hac_diagnostics$selection_path,
      hac_selection_threshold = hac_diagnostics$selection_threshold,
      hac_selection_minimizer_G = hac_diagnostics$selection_minimizer_G
    )
  }

  if (cross_type == "time") {
    split_pt <- floor(T / 2)
    if (split_pt < 1 || split_pt >= T) {
      stop("Cannot split time dimension into two non-empty halves.")
    }
    first_half <- seq_len(split_pt)
    second_half <- seq.int(split_pt + 1, T)

    # Fold 1: train on first half, estimate on second half.
    fold1 <- estimate_time_fold(
      train_t_idx = first_half, test_t_idx = second_half, outer_fold_id = 1L
    )
    # Fold 2: train on second half, estimate on first half.
    fold2 <- estimate_time_fold(
      train_t_idx = second_half, test_t_idx = first_half, outer_fold_id = 2L
    )

    if (!is.numeric(fold1$T_test) || !is.numeric(fold2$T_test)) {
      stop("Fold metadata T_test must be numeric for pooled cross-fit assembly.")
    }

    # Combine fold-level projected design objects into one cross-fitted pooled design.
    y_cf_mat <- matrix(NA_real_, nrow = N, ncol = T)
    d_cf_mat <- matrix(NA_real_, nrow = N, ncol = T)
    x_cf_arr <- array(NA_real_, dim = c(N, T, ncol(fold1$x_matrix)))

    y_cf_mat[, second_half] <- matrix(fold1$tY_vector, nrow = N, ncol = fold1$T_test)
    d_cf_mat[, second_half] <- matrix(fold1$d_vector, nrow = N, ncol = fold1$T_test)
    for (k in seq_len(ncol(fold1$x_matrix))) {
      x_cf_arr[, second_half, k] <- matrix(fold1$x_matrix[, k], nrow = N, ncol = fold1$T_test)
    }

    y_cf_mat[, first_half] <- matrix(fold2$tY_vector, nrow = N, ncol = fold2$T_test)
    d_cf_mat[, first_half] <- matrix(fold2$d_vector, nrow = N, ncol = fold2$T_test)
    for (k in seq_len(ncol(fold2$x_matrix))) {
      x_cf_arr[, first_half, k] <- matrix(fold2$x_matrix[, k], nrow = N, ncol = fold2$T_test)
    }

    df_cf <- (N * fold1$T_test - fold1$T_test * fold1$G_unit) +
      (N * fold2$T_test - fold2$T_test * fold2$G_unit)
  } else {
    split_pt <- floor(N / 2)
    if (split_pt < 1 || split_pt >= N) {
      stop("Cannot split unit dimension into two non-empty halves.")
    }
    first_half <- seq_len(split_pt)
    second_half <- seq.int(split_pt + 1, N)

    # Fold 1: train on first unit fold, estimate on second unit fold.
    fold1 <- estimate_unit_fold(train_i_idx = first_half, test_i_idx = second_half, outer_fold_id = 1L)
    # Fold 2: train on second unit fold, estimate on first unit fold.
    fold2 <- estimate_unit_fold(train_i_idx = second_half, test_i_idx = first_half, outer_fold_id = 2L)

    y_cf_mat <- matrix(NA_real_, nrow = N, ncol = T)
    d_cf_mat <- matrix(NA_real_, nrow = N, ncol = T)
    x_cf_arr <- array(NA_real_, dim = c(N, T, ncol(fold1$x_matrix)))

    y_cf_mat[second_half, ] <- matrix(fold1$tY_vector, nrow = fold1$N_test, ncol = T)
    d_cf_mat[second_half, ] <- matrix(fold1$d_vector, nrow = fold1$N_test, ncol = T)
    for (k in seq_len(ncol(fold1$x_matrix))) {
      x_cf_arr[second_half, , k] <- matrix(fold1$x_matrix[, k], nrow = fold1$N_test, ncol = T)
    }

    y_cf_mat[first_half, ] <- matrix(fold2$tY_vector, nrow = fold2$N_test, ncol = T)
    d_cf_mat[first_half, ] <- matrix(fold2$d_vector, nrow = fold2$N_test, ncol = T)
    for (k in seq_len(ncol(fold2$x_matrix))) {
      x_cf_arr[first_half, , k] <- matrix(fold2$x_matrix[, k], nrow = fold2$N_test, ncol = T)
    }

    df_cf <- (fold1$N_test * T - T * fold1$G_unit_used) +
      (fold2$N_test * T - T * fold2$G_unit_used)
  }

  if (anyNA(y_cf_mat) || anyNA(d_cf_mat) || anyNA(x_cf_arr)) {
    stop("Cross-fitted projected design contains missing values; check fold assembly.")
  }

  tY_vector_cf <- as.vector(y_cf_mat)
  d_vector_cf <- as.vector(d_cf_mat)
  x_matrix_cf <- matrix(NA_real_, nrow = N * T, ncol = dim(x_cf_arr)[3])
  for (k in seq_len(dim(x_cf_arr)[3])) {
    x_matrix_cf[, k] <- as.vector(x_cf_arr[, , k])
  }

  # Run pooled DML on the combined cross-fitted projected data.
  fit <- rlassoEffect(x = x_matrix_cf, y = tY_vector_cf, d = d_vector_cf, method = "double selection")

  trans <- data.frame(y = tY_vector_cf,
                      D = d_vector_cf,
                      x_matrix_cf)

  lasso.Y <- rlasso(y ~ . - D - 1, data = trans)
  Ytilde <- lasso.Y$residuals

  lasso.D <- rlasso(D ~ . - 1, data = trans[, -1, drop = FALSE])
  Dtilde <- lasso.D$residuals

  data_res <- data.frame(id = rep(ids, T),
                         time = rep(times, each = N),
                         Ytilde = Ytilde,
                         Dtilde = Dtilde)

  Post_plm <- plm(Ytilde ~ -1 + Dtilde, data = data_res, model = "pooling", index = c("id", "time"))
  robust_se <- sqrt(diag(vcovHC(Post_plm, type = "HC0", method = "arellano")))

  coefs <- coef(Post_plm)

  df_cf <- max(1, df_cf)

  se_corrected <- robust_se * sqrt(N * T / df_cf)
  t_values_corrected <- coefs / se_corrected
  p_values_corrected <- 2 * pt(-abs(t_values_corrected), df_cf)

  summary_table_correct <- data.frame(
    Estimate = coefs,
    `Std. Error corrected` = se_corrected,
    `t-value corrected` = t_values_corrected,
    `Pr(>|t|) corrected` = p_values_corrected
  )
  colnames(summary_table_correct) <- c("Estimate", "Std. Error corrected", "t-value corrected", "Pr(>|t|) corrected")


  result <- list(
    cluster_method = cluster_method,
    cross_type = cross_type,
    G_unit_fold1 = fold1$G_unit,
    G_unit_fold2 = fold2$G_unit,
    unit_group_fold1 = fold1$unit_group,
    unit_group_fold2 = fold2$unit_group,
    fit_summary = summary(fit),
    post_plm_summary = summary(Post_plm),
    G_unit = (fold1$G_unit + fold2$G_unit) / 2,
    estimate_corrected = summary_table_correct,
    summary_table = summary_table_correct
  )

  if (cluster_method == "ps-kmeans") {
    result$ps_profile_dimension_fold1 <- fold1$ps_profile_dimension
    result$ps_profile_dimension_fold2 <- fold2$ps_profile_dimension
    result$ps_objective_fold1 <- fold1$ps_objective
    result$ps_objective_fold2 <- fold2$ps_objective
    result$ps_train_group_sizes_fold1 <- fold1$ps_train_group_sizes
    result$ps_train_group_sizes_fold2 <- fold2$ps_train_group_sizes
    result$ps_test_group_sizes_fold1 <- fold1$ps_test_group_sizes
    result$ps_test_group_sizes_fold2 <- fold2$ps_test_group_sizes
    result$ps_converged_fold1 <- fold1$ps_converged
    result$ps_converged_fold2 <- fold2$ps_converged
    result$ps_iterations_fold1 <- fold1$ps_iterations
    result$ps_iterations_fold2 <- fold2$ps_iterations
    result$ps_selection_rule_fold1 <- fold1$ps_selection_rule
    result$ps_selection_rule_fold2 <- fold2$ps_selection_rule
    result$ps_selection_path_fold1 <- fold1$ps_selection_path
    result$ps_selection_path_fold2 <- fold2$ps_selection_path
    result$ps_selection_threshold_fold1 <- fold1$ps_selection_threshold
    result$ps_selection_threshold_fold2 <- fold2$ps_selection_threshold
    result$ps_selection_minimizer_G_fold1 <- fold1$ps_selection_minimizer_G
    result$ps_selection_minimizer_G_fold2 <- fold2$ps_selection_minimizer_G
  }

  if (cluster_method == "HAC new") {
    result$hac_profile_dimension_fold1 <- fold1$hac_profile_dimension
    result$hac_profile_dimension_fold2 <- fold2$hac_profile_dimension
    result$hac_reference_size_fold1 <- fold1$hac_reference_size
    result$hac_reference_size_fold2 <- fold2$hac_reference_size
    result$hac_auxiliary_size_fold1 <- fold1$hac_auxiliary_size
    result$hac_auxiliary_size_fold2 <- fold2$hac_auxiliary_size
    result$hac_merge_height_fold1 <- fold1$hac_merge_height
    result$hac_merge_height_fold2 <- fold2$hac_merge_height
    result$hac_evaluation_coverage_fold1 <- fold1$hac_evaluation_coverage
    result$hac_evaluation_coverage_fold2 <- fold2$hac_evaluation_coverage
    result$hac_effective_cut_fold1 <- fold1$hac_effective_cut
    result$hac_effective_cut_fold2 <- fold2$hac_effective_cut
    result$hac_stochastic_rate_fold1 <- fold1$hac_stochastic_rate
    result$hac_stochastic_rate_fold2 <- fold2$hac_stochastic_rate
    result$hac_auxiliary_group_sizes_fold1 <- fold1$hac_auxiliary_group_sizes
    result$hac_auxiliary_group_sizes_fold2 <- fold2$hac_auxiliary_group_sizes
    result$hac_test_group_sizes_fold1 <- fold1$hac_test_group_sizes
    result$hac_test_group_sizes_fold2 <- fold2$hac_test_group_sizes
    result$hac_medoid_auxiliary_index_fold1 <- fold1$hac_medoid_auxiliary_index
    result$hac_medoid_auxiliary_index_fold2 <- fold2$hac_medoid_auxiliary_index
    result$hac_medoid_unit_id_fold1 <- fold1$hac_medoid_unit_id
    result$hac_medoid_unit_id_fold2 <- fold2$hac_medoid_unit_id
    result$hac_reference_unit_id_fold1 <- fold1$hac_reference_unit_id
    result$hac_reference_unit_id_fold2 <- fold2$hac_reference_unit_id
    result$hac_auxiliary_unit_id_fold1 <- fold1$hac_auxiliary_unit_id
    result$hac_auxiliary_unit_id_fold2 <- fold2$hac_auxiliary_unit_id
    result$hac_selection_rule_fold1 <- fold1$hac_selection_rule
    result$hac_selection_rule_fold2 <- fold2$hac_selection_rule
    result$hac_selection_path_fold1 <- fold1$hac_selection_path
    result$hac_selection_path_fold2 <- fold2$hac_selection_path
    result$hac_selection_threshold_fold1 <- fold1$hac_selection_threshold
    result$hac_selection_threshold_fold2 <- fold2$hac_selection_threshold
    result$hac_selection_minimizer_G_fold1 <- fold1$hac_selection_minimizer_G
    result$hac_selection_minimizer_G_fold2 <- fold2$hac_selection_minimizer_G
  }

  if (cluster_method == "kmedoid") {
    result$cross_fitted <- TRUE
    result$kmedoid_profile_dimension_fold1 <- fold1$kmedoid_profile_dimension
    result$kmedoid_profile_dimension_fold2 <- fold2$kmedoid_profile_dimension
    result$kmedoid_objective_fold1 <- fold1$kmedoid_objective
    result$kmedoid_objective_fold2 <- fold2$kmedoid_objective
    result$kmedoid_group_sizes_fold1 <- fold1$kmedoid_group_sizes
    result$kmedoid_group_sizes_fold2 <- fold2$kmedoid_group_sizes
    result$kmedoid_medoid_index_fold1 <- fold1$kmedoid_medoid_index
    result$kmedoid_medoid_index_fold2 <- fold2$kmedoid_medoid_index
    result$kmedoid_medoid_unit_id_fold1 <- fold1$kmedoid_medoid_unit_id
    result$kmedoid_medoid_unit_id_fold2 <- fold2$kmedoid_medoid_unit_id
    result$kmedoid_solver_fold1 <- fold1$kmedoid_solver
    result$kmedoid_solver_fold2 <- fold2$kmedoid_solver
    result$kmedoid_global_optimum_fold1 <- fold1$kmedoid_global_optimum
    result$kmedoid_global_optimum_fold2 <- fold2$kmedoid_global_optimum
    result$kmedoid_converged_fold1 <- fold1$kmedoid_converged
    result$kmedoid_converged_fold2 <- fold2$kmedoid_converged
    result$kmedoid_iterations_fold1 <- fold1$kmedoid_iterations
    result$kmedoid_iterations_fold2 <- fold2$kmedoid_iterations
    result$kmedoid_selection_rule_fold1 <- fold1$kmedoid_selection_rule
    result$kmedoid_selection_rule_fold2 <- fold2$kmedoid_selection_rule
    result$kmedoid_selection_path_fold1 <- fold1$kmedoid_selection_path
    result$kmedoid_selection_path_fold2 <- fold2$kmedoid_selection_path
    result$kmedoid_distance_fold1 <- fold1$kmedoid_distance
    result$kmedoid_distance_fold2 <- fold2$kmedoid_distance
  }

  if (cluster_method == "projection-kmeans") {
    result$cross_fitted <- TRUE
    result$pk_input_dimension_fold1 <- fold1$pk_input_dimension
    result$pk_input_dimension_fold2 <- fold2$pk_input_dimension
    result$pk_feature_dimension_fold1 <- fold1$pk_feature_dimension
    result$pk_feature_dimension_fold2 <- fold2$pk_feature_dimension
    result$pk_objective_fold1 <- fold1$pk_objective
    result$pk_objective_fold2 <- fold2$pk_objective
    result$pk_group_sizes_fold1 <- fold1$pk_group_sizes
    result$pk_group_sizes_fold2 <- fold2$pk_group_sizes
    result$pk_converged_fold1 <- fold1$pk_converged
    result$pk_converged_fold2 <- fold2$pk_converged
    result$pk_iterations_fold1 <- fold1$pk_iterations
    result$pk_iterations_fold2 <- fold2$pk_iterations
    result$pk_selection_rule_fold1 <- fold1$pk_selection_rule
    result$pk_selection_rule_fold2 <- fold2$pk_selection_rule
    result$pk_selection_path_fold1 <- fold1$pk_selection_path
    result$pk_selection_path_fold2 <- fold2$pk_selection_path
    result$pk_selection_threshold_fold1 <- fold1$pk_selection_threshold
    result$pk_selection_threshold_fold2 <- fold2$pk_selection_threshold
    result$pk_reference_G_fold1 <- fold1$pk_reference_G
    result$pk_reference_G_fold2 <- fold2$pk_reference_G
    result$pk_reference_objective_fold1 <- fold1$pk_reference_objective
    result$pk_reference_objective_fold2 <- fold2$pk_reference_objective
    result$pk_inference_multiplier_fold1 <- fold1$pk_inference_multiplier
    result$pk_inference_multiplier_fold2 <- fold2$pk_inference_multiplier
    result$pk_complexity_threshold_fold1 <- fold1$pk_complexity_threshold
    result$pk_complexity_threshold_fold2 <- fold2$pk_complexity_threshold
    result$pk_log_N_rule_fold1 <- fold1$pk_log_N_rule
    result$pk_log_N_rule_fold2 <- fold2$pk_log_N_rule
    result$pk_log_N_value_fold1 <- fold1$pk_log_N_value
    result$pk_log_N_value_fold2 <- fold2$pk_log_N_value
    result$pk_requested_max_G_fold1 <- fold1$pk_requested_max_G
    result$pk_requested_max_G_fold2 <- fold2$pk_requested_max_G
    result$pk_effective_max_G_fold1 <- fold1$pk_effective_max_G
    result$pk_effective_max_G_fold2 <- fold2$pk_effective_max_G
    result$pk_theoretical_max_G_fold1 <- fold1$pk_theoretical_max_G
    result$pk_theoretical_max_G_fold2 <- fold2$pk_theoretical_max_G
    result$pk_max_G_source_fold1 <- fold1$pk_max_G_source
    result$pk_max_G_source_fold2 <- fold2$pk_max_G_source
    result$pk_selection_nstart_fold1 <- fold1$pk_selection_nstart
    result$pk_selection_nstart_fold2 <- fold2$pk_selection_nstart
    result$pk_selection_max_iter_fold1 <- fold1$pk_selection_max_iter
    result$pk_selection_max_iter_fold2 <- fold2$pk_selection_max_iter
    result$pk_screening_objective_fold1 <- fold1$pk_screening_objective
    result$pk_screening_objective_fold2 <- fold2$pk_screening_objective
    result$pk_final_objective_fold1 <- fold1$pk_final_objective
    result$pk_final_objective_fold2 <- fold2$pk_final_objective
    result$pk_noise_floor_fold1 <- fold1$pk_noise_floor
    result$pk_noise_floor_fold2 <- fold2$pk_noise_floor
  }

  result
}

# HP_estimate_double <- function(data, y_col = NULL, covariate_cols = NULL,
#                                id_col = "id", time_col = "time",
#                                index = c("id", "time"), method_auto = "dynamicTreeCut", cluster_method = 'kmeans', pre_cluster = FALSE, Du_pre, Dv_pre, unit_cluster = NULL, time_cluster = NULL, gamma = 1, cc = 0, dim_moment = 1) {
#
#
#   if (is.null(covariate_cols)) {
#     covariate_cols <- setdiff(colnames(data), c(id_col, time_col, y_col))
#   }
#
#   ids <- sort(unique(data[[id_col]]))
#   times <- sort(unique(data[[time_col]]))
#   N <- length(ids)
#   T <- length(times)
#   K <- length(covariate_cols)
#
#   # Initialize y and X
#   y <- matrix(NA_real_, nrow = N, ncol = T)
#   X <- array(NA_real_, dim = c(N, T, K))
#
#   for (i in seq_len(nrow(data))) {
#     id_idx <- which(ids == data[[id_col]][i])
#     time_idx <- which(times == data[[time_col]][i])
#     y[id_idx, time_idx] <- data[[y_col]][i]
#     for (k in seq_along(covariate_cols)) {
#       X[id_idx, time_idx, k] <- data[[covariate_cols[k]]][i]
#     }
#   }
#
#   # Cluster (output indicators assumed one-hot encoded)
#   X_norm <- array(NA_real_, dim = dim(X))  # Same shape: (N, T, K)
#   for (k in 1:dim(X)[3]) {
#     X_slice <- X[,,k]
#     X_norm[,,k] <- scale(X_slice)
#   }
#   y_norm = scale(y)
#   if (cluster_method == 'kmeans'){
#
#     # kmeans
#     X_list = lapply(seq_len(dim(X)[3]), function(k) X[ , , k])
#     if (pre_cluster == FALSE){
#       # cluster
#       if (is.null(unit_cluster) == 1 & is.null(time_cluster) == 1){
#         clusteri <- cluster_general(y, X_list, N, T, init = 50, type = "long", gamma = gamma, cc = cc, dim_moment = dim_moment   )
#         G_unit <- clusteri$clusters
#         klong <- clusteri$res
#
#
#         clustert <- cluster_general(y, X_list, N, T, init = 50, type = "tall", gamma = gamma, cc = cc, dim_moment = dim_moment  )
#         G_time <- clustert$clusters
#         ktall <- clustert$res
#       }else{
#         clusteri <- cluster_general(y, X_list, N, T, init = 50, type = "long" , groups = c(floor(unit_cluster)) )
#         G_unit <- clusteri$clusters
#         klong <- clusteri$res
#
#
#         clustert <- cluster_general(y, X_list, N, T, init = 50, type = "tall" , groups = c(floor(time_cluster)))
#         G_time <- clustert$clusters
#         ktall <- clustert$res
#       }
#
#       Du <- matrix(0, N, G_unit)
#       Dv <- matrix(0, T, G_time)
#       Diagu<- matrix(0, G_unit, G_unit)
#       Diagv<- matrix(0, G_time, G_time)
#
#       for (j in seq_len(G_unit)) {
#         Du[, j] <- as.numeric(klong$cluster == j)
#         Diagu[j,j] <- 1/sum(Du[,j])
#       }
#
#       for (j in seq_len(G_time)) {
#         Dv[, j] <- as.numeric(ktall$cluster == j)
#         Diagv[j,j] <- 1/sum(Dv[,j])
#       }
#     }else if (pre_cluster == TRUE){
#       Du = Du_pre
#       Dv = Dv_pre
#       Diagu<- matrix(0, G_unit, G_unit)
#       Diagv<- matrix(0, G_time, G_time)
#
#       for (j in seq_len(G_unit)) {
#         Diagu[j,j] <- 1/sum(Du[,j])
#       }
#
#       for (j in seq_len(G_time)) {
#         Diagv[j,j] <- 1/sum(Dv[,j])
#       }
#       G_unit = dim(Du)[2]
#       G_time = dim(Dv)[2]
#     }
#   }else if (cluster_method == 'hierarchical'){
#
#     # heriachical
#     if (is.null(unit_cluster) == 1 & is.null(time_cluster) == 1){
#       res_unit <- cluster_Hierarchical(y_norm, X_norm, type = "unit", method_auto = method_auto)
#       res_time <- cluster_Hierarchical(y_norm, X_norm, type = "time", method_auto = method_auto)
#     }else{
#       res_unit <- cluster_Hierarchical(y_norm, X_norm, type = "unit", cluster = unit_cluster)
#       res_time <- cluster_Hierarchical(y_norm, X_norm, type = "time", cluster = time_cluster)
#     }
#
#     Du <- res_unit$indicator    # N x G_unit
#     Dv <- res_time$indicator    # T x G_time
#
#     G_unit = res_unit$G
#     G_time = res_time$G
#
#     Diagu<- matrix(0, G_unit, G_unit)
#     Diagv<- matrix(0, G_time, G_time)
#
#     for (j in seq_len(G_unit)) {
#       Diagu[j,j] <- 1/sum(Du[,j])
#     }
#
#     for (j in seq_len(G_time)) {
#       Diagv[j,j] <- 1/sum(Dv[,j])
#     }
#   }
#
#   unit_group <- max.col(Du)
#   time_group <- max.col(Dv)
#
#   # Projection matrices
#   Mu <- diag(N) - Du %*% Diagu %*% t(Du)
#   Mv <- diag(T) - Dv %*% Diagv %*% t(Dv)
#   # Mu <- diag(N) - Du %*% solve(t(Du) %*% Du) %*% t(Du)
#   # Mv <- diag(T) - Dv %*% solve(t(Dv) %*% Dv) %*% t(Dv)
#
#   # Stack y and X into Z: N x T x (K+1)
#   Z <- array(NA_real_, dim = c(N, T, K + 1))
#   Z[,,1] <- y
#   for (k in 1:K) Z[,,k+1] <- X[,,k]
#
#   # Projection method: demean unit and time clusters slice-wise
#   Z_proj <- array(NA_real_, dim = c(N, T, K + 1))
#   for (k in 1:(K+1)) {
#     Z_proj[,,k] <- Mu %*% Z[,,k] %*% Mv
#   }
#
#
#   y_trans <- Z_proj[,,1] # Z_proj_final_array[,,1]
#   X_trans <- Z_proj[,,2:(K+1)] # Z_proj_final_array[,,2:(K+1)]
#
#   tY_vector <- as.vector(y_trans)  # (N*T)
#   tX_matrix <- matrix(aperm(X_trans, c(1, 2, 3)), nrow = N*T, ncol = K)
#
#   d_vector <- tX_matrix[, 1]
#   x_matrix <- tX_matrix[, -1, drop = FALSE]
#   # Separate transformed y and X
#
#   trans <- data.frame(y = tY_vector,
#                       D = d_vector,
#                       x_matrix)
#
#   # Run rlassoEffect with correct inputs
#   fit <- rlassoEffect(x = x_matrix, y = tY_vector, d = d_vector, method = "double selection")
#
#   trans <- data.frame(y = tY_vector,
#                       D = d_vector,
#                       x_matrix)
#
#   lasso.Y <- rlasso(y ~ . - D - 1, data = trans)
#   Ytilde <- lasso.Y$residuals
#
#   lasso.D <- rlasso(D ~ . - 1, data = trans[, -1])
#   Dtilde <- lasso.D$residuals
#
#   data_res <- data.frame(id = data[[index[1]]],
#                          time = data[[index[2]]],
#                          Ytilde = Ytilde,
#                          Dtilde = Dtilde)
#
#   Post_plm <- plm(Ytilde ~ -1 + Dtilde, data = data_res, model = "pooling", index = c("id", "time"))
#
#   coefs <- coef(Post_plm)
#   se_corrected <- sqrt(vcovHC(Post_plm, type = "HC0", method = "arellano")) * sqrt(N * T / (N*T - N*G_time - T*G_unit))
#   t_values_corrected <- coefs / se_corrected
#
#   df <- Post_plm$df.residual
#   p_values_corrected <- 2 * pt(-abs(t_values_corrected), df)
#
#   summary_table_correct <- data.frame(
#     Estimate = coefs,
#     `Std. Error corrected` = se_corrected,
#     `t-value corrected` = t_values_corrected,
#     `Pr(>|t|) corrected` = p_values_corrected
#   )
#
#   colnames(summary_table_correct) = c('Estimate', 'Std. Error corrected', 't-value corrected', 'Pr(>|t|) corrected')
#
#   # Return cluster counts as well
#   return(list(
#     fit_summary = summary(fit),
#     G_unit = G_unit,
#     G_time = G_time,
#     unit_group = unit_group,
#     time_group = time_group,
#     post_plm_summary = summary(Post_plm),
#     estimate_corrected = summary_table_correct,
#     summary_table = summary_table_correct
#   ))
# }

#' @export
Niave_estimate <- function(data, y_col = NULL, covariate_cols = NULL,
                           id_col = "id", time_col = "time", index = c("id", "time")) {


  if (is.null(covariate_cols)) {
    covariate_cols <- setdiff(colnames(data), c(id_col, time_col, y_col))
  }

  ids <- sort(unique(data[[id_col]]))
  times <- sort(unique(data[[time_col]]))
  N <- length(ids)
  T <- length(times)
  K <- length(covariate_cols)

  # Initialize y and X
  y <- matrix(NA_real_, nrow = N, ncol = T)
  X <- array(NA_real_, dim = c(N, T, K))

  for (i in seq_len(nrow(data))) {
    id_idx <- which(ids == data[[id_col]][i])
    time_idx <- which(times == data[[time_col]][i])
    y[id_idx, time_idx] <- data[[y_col]][i]
    for (k in seq_along(covariate_cols)) {
      X[id_idx, time_idx, k] <- data[[covariate_cols[k]]][i]
    }
  }

  # Stack y and X into Z: N x T x (K+1)
  Z <- array(NA_real_, dim = c(N, T, K + 1))
  Z[,,1] <- y
  for (k in 1:K) Z[,,k+1] <- X[,,k]

  # Projection method: demean unit and time clusters slice-wise
  Z_proj <- array(NA_real_, dim = c(N, T, K + 1))
  for (k in 1:(K+1)) {
    Z_proj[,,k] <-  Z[,,k]
  }

  # Separate transformed y and X
  y_trans <- Z_proj[,,1] # Z_proj_final_array[,,1]
  X_trans <- Z_proj[,,2:(K+1)] # Z_proj_final_array[,,2:(K+1)]

  tY_vector <- as.vector(y_trans)  # (N*T)
  # Handle T=1 case where X_trans is 2D instead of 3D
  if (T == 1) {
    X_trans_mat <- if (length(dim(X_trans)) == 3) X_trans[, 1, , drop = TRUE] else X_trans
    tX_matrix <- X_trans_mat
  } else if (length(dim(X_trans)) == 2) {
    tX_matrix <- matrix(as.vector(t(X_trans)), nrow = N * T, ncol = K)
  } else {
    tX_matrix <- matrix(aperm(X_trans, c(1, 2, 3)), nrow = N * T, ncol = K)
  }

  d_vector <- tX_matrix[, 1]
  x_matrix <- tX_matrix[, -1, drop = FALSE]

  # Run rlassoEffect with correct inputs
  fit <- rlassoEffect(x = x_matrix, y = tY_vector, d = d_vector, method = "double selection")

  trans <- data.frame(y = tY_vector,
                      D = d_vector,
                      x_matrix)

  lasso.Y <- rlasso(y ~ . - D - 1, data = trans)
  Ytilde <- lasso.Y$residuals

  lasso.D <- rlasso(D ~ . - 1, data = trans[, -1, drop = FALSE])
  Dtilde <- lasso.D$residuals

  # Create proper panel structure for residuals
  # as.vector() on N x T matrix goes column-wise: unit 1 time 1, unit 2 time 1, ..., unit N time 1, unit 1 time 2, ...
  data_res <- data.frame(id = rep(ids, T),
                         time = rep(times, each = N),
                         Ytilde = Ytilde,
                         Dtilde = Dtilde)

  # For T=1, use lm() since there's no panel structure; for T>1 use plm with pooling
  if (T == 1) {
    Post_plm <- lm(Ytilde ~ -1 + Dtilde, data = data_res)
  } else {
    Post_plm <- plm(Ytilde ~ -1 + Dtilde, data = data_res, model = "pooling", index = c("id", "time"))
  }

  coefs <- coef(Post_plm)

  # Use appropriate variance estimator based on model type
  if (T == 1) {
    robust_se <- sqrt(diag(vcovHC(Post_plm, type = "HC0")))
    dof_denom <- max(1, N * T - 1)
    se_corrected <- robust_se * sqrt(N * T / dof_denom)
  } else {
    se_corrected <- sqrt(diag(vcovHC(Post_plm, type = "HC0", method = "arellano")))
  }
  t_values_corrected <- coefs / se_corrected

  df <- Post_plm$df.residual
  p_values_corrected <- 2 * pt(-abs(t_values_corrected), df)

  summary_table_correct <- data.frame(
    Estimate = coefs,
    `Std. Error corrected` = se_corrected,
    `t-value corrected` = t_values_corrected,
    `Pr(>|t|) corrected` = p_values_corrected
  )

  colnames(summary_table_correct) = c('Estimate', 'Std. Error corrected', 't-value corrected', 'Pr(>|t|) corrected')

  # Return cluster counts as well
  return(list(
    fit_summary = summary(fit),
    post_plm_summary = summary(Post_plm),
    estimate_corrected = summary_table_correct,
    summary_table = summary_table_correct
  ))
}

#' @export
# TWFE_estimate treat unit and time clusters are both 1
# TWFE_estimate <- function(data, y_col = NULL, covariate_cols = NULL,
#                           id_col = "id", time_col = "time", index = c("id", "time")) {
#   ...
# }

TWFE_estimate <- function(data, y_col = NULL, covariate_cols = NULL,
                          id_col = "id", time_col = "time", index = c("id", "time")) {

  if (is.null(covariate_cols)) {
    covariate_cols <- setdiff(colnames(data), c(id_col, time_col, y_col))
  }

  ids <- sort(unique(data[[id_col]]))
  times <- sort(unique(data[[time_col]]))
  N <- length(ids)
  T <- length(times)
  K <- length(covariate_cols)

  # Initialize y and X as balanced panel containers
  y <- matrix(NA_real_, nrow = N, ncol = T)
  X <- array(NA_real_, dim = c(N, T, K))

  for (i in seq_len(nrow(data))) {
    id_idx <- which(ids == data[[id_col]][i])
    time_idx <- which(times == data[[time_col]][i])
    y[id_idx, time_idx] <- data[[y_col]][i]
    for (k in seq_along(covariate_cols)) {
      X[id_idx, time_idx, k] <- data[[covariate_cols[k]]][i]
    }
  }

  # Residualize with respect to fixed effects before running DML.
  # When T = 1, two-way FE is not identified; fall back to no time demeaning.
  if (T == 1) {
    y_res <- y
    X_res <- X
  } else {
    Mu <- diag(N) - matrix(1 / N, nrow = N, ncol = N)
    Mv <- diag(T) - matrix(1 / T, nrow = T, ncol = T)

    y_res <- Mu %*% y %*% Mv
    X_res <- array(NA_real_, dim = c(N, T, K))
    for (k in seq_len(K)) {
      X_res[,,k] <- Mu %*% X[,,k] %*% Mv
    }
  }

  tY_vector <- as.vector(y_res)
  if (T == 1) {
    X_res_mat <- if (length(dim(X_res)) == 3) X_res[, 1, , drop = TRUE] else X_res
    tX_matrix <- X_res_mat
  } else if (length(dim(X_res)) == 2) {
    tX_matrix <- matrix(as.vector(t(X_res)), nrow = N * T, ncol = K)
  } else {
    tX_matrix <- matrix(aperm(X_res, c(1, 2, 3)), nrow = N * T, ncol = K)
  }

  d_vector <- tX_matrix[, 1]
  x_matrix <- tX_matrix[, -1, drop = FALSE]

  # DML on residualized outcome and regressors
  fit <- rlassoEffect(x = x_matrix, y = tY_vector, d = d_vector, method = "double selection")

  trans <- data.frame(
    y = tY_vector,
    D = d_vector,
    x_matrix
  )

  lasso.Y <- rlasso(y ~ . - D - 1, data = trans)
  Ytilde <- lasso.Y$residuals

  lasso.D <- rlasso(D ~ . - 1, data = trans[, -1, drop = FALSE])
  Dtilde <- lasso.D$residuals

  data_res <- data.frame(
    id = rep(ids, T),
    time = rep(times, each = N),
    Ytilde = Ytilde,
    Dtilde = Dtilde
  )

  if (T == 1) {
    Post_plm <- lm(Ytilde ~ -1 + Dtilde, data = data_res)
    robust_se <- sqrt(diag(vcovHC(Post_plm, type = "HC0")))
    dof_denom <- max(1, N * T - 1)
    G_time <- 1
  } else {
    Post_plm <- plm(Ytilde ~ -1 + Dtilde, data = data_res, model = "pooling", index = c("id", "time"))
    robust_se <- sqrt(diag(vcovHC(Post_plm, type = "HC0", method = "arellano")))
    dof_denom <- max(1, N * T - N - T + 1)
    G_time <- T
  }

  G_unit <- N
  unit_group <- seq_len(N)
  time_group <- if (T == 1) rep(1, T) else seq_len(T)

  coefs <- coef(Post_plm)
  se_corrected <- robust_se * sqrt(N * T / dof_denom)
  t_values_corrected <- coefs / se_corrected

  df <- Post_plm$df.residual
  p_values_corrected <- 2 * pt(-abs(t_values_corrected), df)

  summary_table_correct <- data.frame(
    Estimate = coefs,
    `Std. Error corrected` = se_corrected,
    `t-value corrected` = t_values_corrected,
    `Pr(>|t|) corrected` = p_values_corrected
  )

  colnames(summary_table_correct) <- c("Estimate", "Std. Error corrected", "t-value corrected", "Pr(>|t|) corrected")

  return(list(
    fit_summary = summary(fit),
    G_unit = G_unit,
    G_time = G_time,
    unit_group = unit_group,
    time_group = time_group,
    post_plm_summary = summary(Post_plm),
    estimate_corrected = summary_table_correct,
    summary_table = summary_table_correct
  ))
}

# HP_estimate_twoway <- function(data, y_col = NULL, covariate_cols = NULL,
#                         id_col = "id", time_col = "time",
#                         index = c("id", "time"), method_auto = "dynamicTreeCut", cluster = NULL) {
#
#
#   if (is.null(covariate_cols)) {
#     covariate_cols <- setdiff(colnames(data), c(id_col, time_col, y_col))
#   }
#
#   ids <- sort(unique(data[[id_col]]))
#   times <- sort(unique(data[[time_col]]))
#   N <- length(ids)
#   T <- length(times)
#   K <- length(covariate_cols)
#
#   # Initialize y and X
#   y <- matrix(NA_real_, nrow = N, ncol = T)
#   X <- array(NA_real_, dim = c(N, T, K))
#
#   for (i in seq_len(nrow(data))) {
#     id_idx <- which(ids == data[[id_col]][i])
#     time_idx <- which(times == data[[time_col]][i])
#     y[id_idx, time_idx] <- data[[y_col]][i]
#     for (k in seq_along(covariate_cols)) {
#       X[id_idx, time_idx, k] <- data[[covariate_cols[k]]][i]
#     }
#   }
#
#   # Cluster (output indicators assumed one-hot encoded)
#   X_norm <- array(NA_real_, dim = dim(X))  # Same shape: (N, T, K)
#   for (k in 1:dim(X)[3]) {
#     X_slice <- X[,,k]
#     X_norm[,,k] <- scale(X_slice)
#   }
#   y_norm = scale(y)
#   res_unit <- cluster_Hierarchical(y_norm, X_norm, type = "unit", method_auto = method_auto)
#   res_time <- cluster_Hierarchical(y_norm, X_norm, type = "time", method_auto = method_auto)
#   res_covar <- cluster_Hierarchical(y_norm, X_norm, type = "covariate", method_auto = method_auto, cluster = cluster)
#
#   Dv <- res_unit$indicator    # N x G_unit
#   Du <- res_time$indicator    # T x G_time
#   Dc <- res_covar$indicator   # (K+1) x G_covar
#
#   G_unit = res_unit$G
#   G_time = res_time$G
#   G_covar = res_covar$G
#
#   unit_group <- max.col(Dv)
#   time_group <- max.col(Du)
#   covar_group <- max.col(Dc)
#
#   # Projection matrices
#   Mv <- diag(N) - Dv %*% solve(t(Dv) %*% Dv) %*% t(Dv)
#   Mu <- diag(T) - Du %*% solve(t(Du) %*% Du) %*% t(Du)
#   Mc <- diag(K+1) - Dc %*% solve(t(Dc) %*% Dc) %*% t(Dc)
#
#   # Stack y and X into Z: N x T x (K+1)
#   Z <- array(NA_real_, dim = c(N, T, K + 1))
#   Z[,,1] <- y
#   for (k in 1:K) Z[,,k+1] <- X[,,k]
#
#   # Projection method: demean unit and time clusters slice-wise
#   Z_proj <- array(NA_real_, dim = c(N, T, K + 1))
#   for (k in 1:(K+1)) {
#     Z_proj[,,k] <- Mv %*% Z[,,k] %*% Mu
#   }
#   # Separate transformed y and X
#   y_trans <- Z_proj[,,1] # Z_proj_final_array[,,1]
#   X_trans <- Z_proj[,,2:(K+1)] # Z_proj_final_array[,,2:(K+1)]
#
#   tY_vector <- as.vector(y_trans)  # (N*T)
#   tX_matrix <- matrix(aperm(X_trans, c(1, 2, 3)), nrow = N*T, ncol = K)
#
#   d_vector <- tX_matrix[, 1]
#   x_matrix <- tX_matrix[, -1, drop = FALSE]
#
#   # Run rlassoEffect with correct inputs
#   fit <- rlassoEffect(x = x_matrix, y = tY_vector, d = d_vector, method = "partialling out")
#
#
#   trans <- data.frame(y = tY_vector,
#                       D = d_vector,
#                       x_matrix)
#
#   lasso.Y <- rlasso(y ~ . - D - 1, data = trans)
#   Ytilde <- lasso.Y$residuals
#
#   lasso.D <- rlasso(D ~ . - 1, data = trans[, -1])
#   Dtilde <- lasso.D$residuals
#
#   data_res <- data.frame(id = data[[index[1]]],
#                          time = data[[index[2]]],
#                          Ytilde = Ytilde,
#                          Dtilde = Dtilde)
#
#   Post_plm <- plm(Ytilde ~ -1 + Dtilde, data = data_res, model = "pooling", index = c("id", "time"))
#
#   coefs <- coef(Post_plm)
#   se_corrected <- sqrt(vcovHC(Post_plm, type = "HC0", method = "arellano")) * sqrt(N * T / ((N-G_unit)*(T-G_time)))
#   t_values_corrected <- coefs / se_corrected
#
#   df <- Post_plm$df.residual
#   p_values_corrected <- 2 * pt(-abs(t_values_corrected), df)
#
#   summary_table_correct <- data.frame(
#     Estimate = coefs,
#     `Std. Error corrected` = se_corrected,
#     `t-value corrected` = t_values_corrected,
#     `Pr(>|t|) corrected` = p_values_corrected
#   )
#
#   summary_table <- summary(Post_plm)
#   summary_table$coefficients <- cbind(summary_table_correct, summary_table$coefficients)
#
#   # Return cluster counts as well
#   return(list(
#     fit_summary = summary(fit),
#     G_unit = G_unit,
#     G_time = G_time,
#     unit_group = unit_group,
#     time_group = time_group,
#     post_plm_summary = summary(Post_plm),
#     estimate_corrected = summary_table_correct,
#     summary_table = summary_table
#   ))
# }
#
# HP_estimate_double_twoway <- function(data, y_col = NULL, covariate_cols = NULL,
#                                id_col = "id", time_col = "time",
#                                index = c("id", "time"), method_auto = "dynamicTreeCut", cluster = NULL) {
#
#
#   if (is.null(covariate_cols)) {
#     covariate_cols <- setdiff(colnames(data), c(id_col, time_col, y_col))
#   }
#
#   ids <- sort(unique(data[[id_col]]))
#   times <- sort(unique(data[[time_col]]))
#   N <- length(ids)
#   T <- length(times)
#   K <- length(covariate_cols)
#
#   # Initialize y and X
#   y <- matrix(NA_real_, nrow = N, ncol = T)
#   X <- array(NA_real_, dim = c(N, T, K))
#
#   for (i in seq_len(nrow(data))) {
#     id_idx <- which(ids == data[[id_col]][i])
#     time_idx <- which(times == data[[time_col]][i])
#     y[id_idx, time_idx] <- data[[y_col]][i]
#     for (k in seq_along(covariate_cols)) {
#       X[id_idx, time_idx, k] <- data[[covariate_cols[k]]][i]
#     }
#   }
#
#   # Cluster (output indicators assumed one-hot encoded)
#   X_norm <- array(NA_real_, dim = dim(X))  # Same shape: (N, T, K)
#   for (k in 1:dim(X)[3]) {
#     X_slice <- X[,,k]
#     X_norm[,,k] <- scale(X_slice)
#   }
#   y_norm = scale(y)
#   res_unit <- cluster_Hierarchical(y_norm, X_norm, type = "unit", method_auto = method_auto)
#   res_time <- cluster_Hierarchical(y_norm, X_norm, type = "time", method_auto = method_auto)
#   res_covar <- cluster_Hierarchical(y_norm, X_norm, type = "covariate", method_auto = method_auto, cluster = cluster)
#
#   Dv <- res_unit$indicator    # N x G_unit
#   Du <- res_time$indicator    # T x G_time
#   Dc <- res_covar$indicator   # (K+1) x G_covar
#
#   G_unit = res_unit$G
#   G_time = res_time$G
#   G_covar = res_covar$G
#
#   unit_group <- max.col(Dv)
#   time_group <- max.col(Du)
#   covar_group <- max.col(Dc)
#
#   # Projection matrices
#   Mv <- diag(N) - Dv %*% solve(t(Dv) %*% Dv) %*% t(Dv)
#   Mu <- diag(T) - Du %*% solve(t(Du) %*% Du) %*% t(Du)
#   Mc <- diag(K+1) - Dc %*% solve(t(Dc) %*% Dc) %*% t(Dc)
#
#   # Stack y and X into Z: N x T x (K+1)
#   Z <- array(NA_real_, dim = c(N, T, K + 1))
#   Z[,,1] <- y
#   for (k in 1:K) Z[,,k+1] <- X[,,k]
#
#   # Projection method: demean unit and time clusters slice-wise
#   Z_proj <- array(NA_real_, dim = c(N, T, K + 1))
#   for (k in 1:(K+1)) {
#     Z_proj[,,k] <- Mv %*% Z[,,k] %*% Mu
#   }
#
#   # Flatten Z_proj for covariate demeaning
#   Z_proj_mat <- matrix(NA_real_, nrow = N * T, ncol = K + 1)
#   for (k in 1:(K+1)) {
#     Z_proj_mat[,k] <- as.vector(t(Z_proj[,,k]))
#   }
#
#   # Covariate demeaning via projection
#   Z_proj_final <- Z_proj_mat %*% Mc
#
#   # Reshape back to array
#   Z_proj_final_array <- array(NA_real_, dim = c(N, T, K + 1))
#   for (k in 1:(K+1)) {
#     Z_proj_final_array[,,k] <- matrix(Z_proj_final[,k], nrow = N, ncol = T, byrow = TRUE)
#   }
#
#   first_val <- covar_group[1]
#   count_first <- sum(covar_group == first_val)
#
#   if (count_first == 1) {
#     # Find the most frequent value in covar_group
#     most_freq_val <- as.numeric(names(sort(table(covar_group), decreasing = TRUE)[1]))
#     covar_group[1] <- most_freq_val
#   }
#   #covar_group[1:length(covar_group)]=1
#   #Z_proj = group_demean_formula_cpp(Z, unit_group, time_group, covar_group)
#   y_trans <- Z_proj[,,1] # Z_proj_final_array[,,1]
#   X_trans <- Z_proj[,,2:(K+1)] # Z_proj_final_array[,,2:(K+1)]
#
#   tY_vector <- as.vector(y_trans)  # (N*T)
#   tX_matrix <- matrix(aperm(X_trans, c(1, 2, 3)), nrow = N*T, ncol = K)
#
#   d_vector <- tX_matrix[, 1]
#   x_matrix <- tX_matrix[, -1, drop = FALSE]
#   # Separate transformed y and X
#
#   trans <- data.frame(y = tY_vector,
#                       D = d_vector,
#                       x_matrix)
#
#   trans = as.data.frame(trans)
#
#   dml_data <- DoubleMLData$new(
#     data = trans,
#     y_col = "y",
#     d_cols = "D",
#     x_cols = c(colnames(trans)[-c(1,2)])
#
#   )
#
#   # Define LASSO machine learning learners for nuisance parameter estimation
#   ml_l <- lrn("regr.cv_glmnet", s = "lambda.min")  # Outcome regression model
#   ml_m <- lrn("regr.cv_glmnet", s = "lambda.min")  # Treatment model (if applicable)
#
#   # Fit the Double Machine Learning model for treatment effect estimation
#   dml_plr <- DoubleMLPLR$new(dml_data, ml_l = ml_l, ml_m = ml_m)
#
#   # Fit the model to estimate the causal effect
#   dml_plr$fit(store_predictions=TRUE)
#   g_hat <- dml_plr$predictions$ml_l
#   m_hat <- dml_plr$predictions$ml_m
#
#   # Step 2: Compute residuals
#   Ytilde <- tY_vector - g_hat
#   Dtilde <- d_vector - m_hat
#   # Double ML
#   data_res = data.frame(id = data[[index[1]]], time = data[[index[2]]], Ytilde, Dtilde)
#   Post_plm = plm(Ytilde ~ -1 + Dtilde, data = data_res, model = "pooling", index=c("id", "time"))
#
#   coefs <- coef(Post_plm)
#   se_corrected <- sqrt(vcovHC(Post_plm, type = "HC0", method = "arellano")) * sqrt(N * T / ((N-G_unit)*(T-G_time)))
#   t_values_corrected <- coefs / se_corrected
#
#   # Calculate p-values from t-distribution for each coefficient
#   df <- Post_plm$df.residual  # degrees of freedom
#   p_values_corrected <- 2 * pt(-abs(t_values_corrected), df)
#
#   summary_table_correct <- data.frame(
#     Estimate = coefs,
#     SE_Corrected = se_corrected,
#     t_value_Corrected = t_values_corrected,
#     p_value_Corrected = p_values_corrected
#   )
#   colnames(summary_table_correct) = c('Estimate', 'Std. Error corrected', 't-value corrected', 'Pr(>|t|) corrected')
#   summary_table = summary(Post_plm)
#   summary_table$coefficients = cbind(summary_table_correct, summary_table$coefficients)
#
#   # Return cluster counts as well
#   return(list(
#     fit_summary = dml_plr$summary(),
#     G_unit = G_unit,
#     G_time = G_time,
#     G_covar = G_covar,
#     unit_group = unit_group,
#     time_group = time_group,
#     covar_group = covar_group,
#     post_plm_summary = summary(Post_plm),
#     estimate_corrected = summary_table_correct,
#     summary_table = summary_table
#   ))
# }

#' @export
compute_u_hat <- function(z_array, unit_clusters, time_clusters, covar_clusters) {
  N <- dim(z_array)[1]
  T <- dim(z_array)[2]
  K <- dim(z_array)[3]

  u_hat <- array(0, dim = c(N, T, K))

  for (i in 1:N) {
    for (t in 1:T) {
      for (k in 1:K) {
        g_i <- unit_clusters[i]
        m_t <- time_clusters[t]
        l_k <- covar_clusters[k]

        # z_{itk}
        zitk <- z_array[i, t, k]

        # bar_z_{g_i t k}
        group_i_indices <- which(unit_clusters == g_i)
        bar_g_i_t_k <- mean(z_array[group_i_indices, t, k])

        # bar_z_{i m_t k}
        time_m_indices <- which(time_clusters == m_t)
        bar_i_m_t_k <- mean(z_array[i, time_m_indices, k])

        # bar_z_{i t l_k}
        covar_l_indices <- which(covar_clusters == l_k)
        bar_i_t_l_k <- mean(z_array[i, t, covar_l_indices])

        # bar_z_{g_i m_t k}
        bar_g_i_m_t_k <- mean(z_array[group_i_indices, time_m_indices, k])

        # bar_z_{g_i t l_k}
        bar_g_i_t_l_k <- mean(z_array[group_i_indices, t, covar_l_indices])

        # bar_z_{i m_t l_k}
        bar_i_m_t_l_k <- mean(z_array[i, time_m_indices, covar_l_indices])

        # Final u_hat_{itk}
        u_hat[i, t, k] <- 3 * zitk -
          2 * bar_g_i_t_k -
          2 * bar_i_m_t_k -
          2 * bar_i_t_l_k +
          bar_g_i_m_t_k +
          bar_g_i_t_l_k +
          bar_i_m_t_l_k
      }
    }
  }

  return(u_hat)
}

#' @export
cluster_Hierarchical <- function(y, X, link = "average", threshold = NULL, cluster = NULL,
                                 type = c("unit", "time", "covariate"),
                                 method_auto = c("none", "silhouette", "gap", "dynamicTreeCut"),
                                 data_for_gap = NULL, max_k = 10, deepSplit = TRUE, minClusterSize = 1, pamStage = FALSE,
                                 dist_type = c("proj", "L2")) {
  type <- match.arg(type)
  method_auto <- match.arg(method_auto)
  dist_type <- match.arg(dist_type)

  N <- nrow(y)
  T <- ncol(y)
  K <- dim(X)[3]

  # Combine y and X into N x T x (K+1)
  combined_array <- array(0, dim = c(N, T, K + 1))
  combined_array[,,1] <- y
  combined_array[,,2:(K+1)] <- X

  # Compute distance matrix based on type
  if (type == "unit") {
    if (dist_type == "proj") {
      dist_mat <- pseudo_dist_unit(combined_array)
    } else if (dist_type == "L2") {
      # d_L2(i,j) = sqrt((1/(T*(K+1))) * sum_{t,k} (z_itk - z_jtk)^2)
      unit_mat <- matrix(combined_array, nrow = N, ncol = T * (K + 1))
      dist_mat <- as.matrix(dist(unit_mat, method = "euclidean")) / sqrt(T * (K + 1))
    }
  } else if (type == "time") {
    if (dist_type == "proj") {
      dist_mat <- pseudo_dist_time(combined_array)
    } else if (dist_type == "L2") {
      # d_L2(p,q) = sqrt((1/(N*(K+1))) * sum_{i,k} (z_ipk - z_iqk)^2)
      time_mat <- matrix(aperm(combined_array, c(2, 1, 3)), nrow = T, ncol = N * (K + 1))
      dist_mat <- as.matrix(dist(time_mat, method = "euclidean")) / sqrt(N * (K + 1))
    }
  } else if (type == "covariate") {
    dist_mat <- pseudo_dist_covariate(combined_array)
  } else {
    stop("Invalid 'type' argument.")
  }

  dist_obj <- as.dist(dist_mat)
  hc <- hclust(dist_obj, method = link)

  if (!is.null(cluster)) {
    G <- cluster
    clusters <- cutree(hc, k = G)

  } else if (!is.null(threshold)) {
    clusters <- cutree(hc, h = threshold)
    G <- length(unique(clusters))

  } else if (method_auto != "none") {
    max_k <- min(max_k, ifelse(type == "unit", N, ifelse(type == "time", T, K + 1)) - 1)

    if (method_auto == "gap") {
      if (is.null(data_for_gap)) {
        stop("For method_auto = 'gap', please provide 'data_for_gap' matrix.")
      }
      gap_fun <- function(x, k) {
        dist_x <- dist(x)
        hc_x <- hclust(dist_x, method = link)
        clust <- cutree(hc_x, k = k)
        list(cluster = clust)  # must return list with $cluster
      }
      gap_stat <- cluster::clusGap(data_for_gap, FUN = gap_fun, K.max = max_k, B = 50)
      G <- maxSE(gap_stat$Tab[, "gap"], gap_stat$Tab[, "SE.sim"], method = "firstSEmax")
      clusters <- cutree(hc, k = G)

    } else if (method_auto == "silhouette") {
      sil_scores <- numeric(max_k)
      sil_scores[1] <- NA # silhouette not defined for k=1
      for (k in 2:max_k) {
        clust_try <- cutree(hc, k = k)
        ss <- silhouette(clust_try, dist_obj)
        sil_scores[k] <- mean(ss[, 3])
      }
      G <- which.max(sil_scores)
      clusters <- cutree(hc, k = G)

    } else if (method_auto == "dynamicTreeCut") {
      clusters <- cutreeDynamic(dendro = hc, distM = as.matrix(dist_obj),
                                deepSplit = deepSplit, minClusterSize = minClusterSize, pamStage = pamStage)
      G <- length(unique(clusters[clusters > 0]))  # exclude noise (0)
      noise_idx <- which(clusters == 0)
      if (length(noise_idx) > 0) {
        clusters[noise_idx] <- G + 1
        G <- G + 1
      }
    }

  } else {
    stop("Must provide either 'cluster', 'threshold', or set 'method_auto' != 'none'.")
  }

  # Format clusters and indicator matrix (like kmeans())
  clusters <- as.integer(factor(clusters))
  G <- length(unique(clusters))
  indicator <- sapply(1:G, function(g) as.integer(clusters == g))
  colnames(indicator) <- paste0("Cluster", 1:G)

  dim_cluster <- switch(type,
                        unit = N,
                        time = T,
                        covariate = K + 1)

  rownames(indicator) <- paste0(type, "_", 1:dim_cluster)

  list(G = G, clusters = clusters, indicator = indicator)
}

#' @export
cluster_general <- function(Y, X_list, N, T, init = 100, type = "long", groups = NULL, cc = 0, gamma = 1, dim_moment = 1) {
  dimtheta <- length(X_list)
  mdim <- dim_moment
  X_list <- c(list(Y), X_list)
  K = dimtheta
  if (type == "long") {
    mom_i <- c()
    ## --- Unit-side moments ---
    X_av <- matrix(0, nrow = N, ncol = dim_moment)  # N x dim_moment

    if (dim_moment == 1) {
      covar_groups <- list(seq_len(length(X_list)))
    } else {
      covar_groups <- split(
        seq_len(length(X_list)),
        cut(seq_len(length(X_list)), dim_moment, labels = FALSE)
      )
    }

    X_mat <- function(X) {
      if (is.null(dim(X))) {
        matrix(X, nrow = N, ncol = T)
      } else {
        X
      }
    }

    for (p in 1:dim_moment) {
      # For each power
      X_pow <- matrix(0, nrow = N, ncol = T)  # N x T

      for (t in 1:T) {
        # Collect all covariates at time t
        X_t <- do.call(cbind, lapply(X_list, function(X) X_mat(X)[, t]))  # N x (K+1)
        # Average the covariates in the p-th group
        cols <- covar_groups[[p]]
        X_pow[, t] <- rowMeans(X_t[, cols, drop = FALSE])
        #     if ( p == 1){
        #     X_pow[, t] <- rowMeans(X_t)             # (1/K) sum_k x_itk^p
        #     }else if (p == 2){
        #       X_pow[, t] <- rowMeans(tanh(X_t))
        #     }else if (p == 3){
        #       X_pow[, t] <- rowMeans((X_t/(1+abs(X_t))))
        #     }else{
        #       X_pow[, t] <- rowMeans(X_t^p)
        #   }
      }

      # Average over time
      X_av[, p] <- rowMeans(X_pow)               # (1/T) sum_t (1/K) sum_k x_itk^p
    }

    mom_i <- cbind(mom_i, X_av)

    # mom_micro <- c()
    mom_micro <- c()
    X_cbind <- do.call(cbind, lapply(X_list, function(X) X_mat(X)))
    mom_micro <- cbind(mom_micro, X_cbind)
    # for (p in 1:dim_moment) {
    #   if ( p == 1){
    #     mom_micro <- cbind(mom_micro, X_cbind^p)             # (1/K) sum_k x_itk^p
    #   }else if (p == 2){
    #     mom_micro <- cbind(mom_micro, tanh(X_cbind))
    #   }else if (p == 3){
    #     mom_micro <- cbind(mom_micro, (X_cbind/(1+abs(X_cbind))))
    #   }else{
    #   mom_micro <- cbind(mom_micro, X_cbind^p)
    # }
    # }
    group_lengths <- sapply(covar_groups, length)

    # new sizes
    new_lengths <- group_lengths * T

    # full sequence after extension
    full_seq <- seq_len((K+1) * T)

    # split into contiguous blocks
    col_s <- split(full_seq, rep(seq_len(dim_moment), new_lengths))
    ## --- Demean & rescale ---
    for (j in 1:mdim) {
      col_mean <- mean(mom_i[, j])
      col_sd <- sd(mom_i[, j])
      if (col_sd > 0) {
        cols <-  as.vector(col_s[[j]])
        mom_micro[, cols] <- (mom_micro[, cols] - col_mean) / col_sd
        mom_i[, j] <- (mom_i[, j] - col_mean) / col_sd
      }
    }


    # --- Variance / noise on rescaled data ---
    variance <-  sum(sapply(1:mdim, function(j) {
      norm(mom_micro[,  as.vector(col_s[[j]]) ] - mom_i[,j], type = "F")^2
    })) / (N * T^2 * (length(covar_groups[[j]]))^2)

    data <- mom_i
    dim_size <- N

  } else if (type == "tall") {
    ## --- Time-side moments ---
    mom_t <- c()
    X_av <- matrix(0, nrow = T, ncol = dim_moment)  # T x dim_moment

    for (p in 1:dim_moment) {
      # For each power
      X_pow <- matrix(0, nrow = N, ncol = T)  # N x T

      for (i in 1:N) {
        # Collect all covariates at time t
        X_i <- sapply(X_list, function(X) X[i, ])  # T x K
        # Take power p **before** averaging over K
        X_pow[i, ] <- rowMeans(X_i^p)             # (1/K) sum_k x_itk^p
      }

      # Average over unit
      X_av[, p] <- colMeans(X_pow)               # (1/N) sum_i (1/K) sum_k x_itk^p
    }

    mom_t <- cbind(mom_t, X_av)

    mom_micro2 <- c()
    X_cbind <- do.call(cbind, lapply(X_list, t))

    for (p in 1:dim_moment) {
      mom_micro2 <- cbind(mom_micro2, (X_cbind)^p)
    }

    ## --- Demean & rescale ---
    for (j in 1:mdim) {
      col_mean <- mean(mom_t[, j])
      col_sd <- sd(mom_t[, j])
      if (col_sd > 0) {
        cols <- ((j - 1) * N * (K+1) + 1):(j * N * (K+1))
        mom_micro2[, cols] <- (mom_micro2[, cols] - col_mean) / col_sd
        mom_t[, j] <- (mom_t[, j] - col_mean) / col_sd
      }
    }


    ## --- Variance / noise on rescaled data ---
    variance <- sum(sapply(1:mdim, function(i) {
      norm(mom_micro2[, ((i - 1) * N * (K+1) + 1):(i * N * (K+1))] - mom_t[,i], type = "F")^2
    })) / (T * N^2 * (K+1)^2)

    data <- mom_t
    dim_size <- T

  } else if(type == "long_T_moment") {
    mom_i <- c()
    ## --- Unit-side moments ---
    dm <- if (T == 1) 1 else dim_moment
    X_av <- matrix(0, nrow = N, ncol = dm)  # N x dim_moment

    time_groups <- if (T == 1 || dm <= 1) {
      list(seq_len(T))
    } else {
      split(seq_len(T), cut(seq_len(T), dm, labels = FALSE))
    }
    X_mean_t <- matrix(0, nrow = N, ncol = T)

    for (t in seq_len(T)) {
      # Collect all covariates at time t
      X_t <- sapply(X_list, function(X) X[, t])  # N x (K+1)
      # Average across covariates first
      X_mean_t[, t] <- rowMeans(X_t)
    }

    for (p in seq_len(dm)) {
      t_idx <- time_groups[[p]]
      X_av[, p] <- rowMeans(X_mean_t[, t_idx, drop = FALSE])
    }

    mom_i <- cbind(mom_i, X_av)

    # mom_micro <- c()
    mom_micro <- c()
    X_cbind <- do.call(cbind, X_list)
    mom_micro <- cbind(mom_micro, X_cbind)


    ## --- Demean & rescale ---
    col_s <- lapply(time_groups, function(t_idx) {
      unlist(lapply(t_idx, function(tt) {
        ((tt - 1) * (K + 1) + 1):(tt * (K + 1))
      }))
    })

    for (j in seq_len(dm)) {
      col_mean <- mean(mom_i[, j])
      col_sd <- sd(mom_i[, j])
      if (col_sd > 0) {
        cols <- as.vector(col_s[[j]])
        mom_micro[, cols] <- (mom_micro[, cols] - col_mean) / col_sd
        mom_i[, j] <- (mom_i[, j] - col_mean) / col_sd
      }
    }

    # --- Variance / noise on rescaled data ---
    variance <- sum(sapply(seq_len(dm), function(j) {
      cols <- as.vector(col_s[[j]])
      group_len <- length(time_groups[[j]])
      denom <- N * (K + 1)^2 * (group_len^2)
      norm(mom_micro[, cols] - mom_i[, j], type = "F")^2 / denom
    }))

    data <- mom_i
    dim_size <- N

  }else if (type == "long multiple moments") {
    mom_i <- c()
    ## --- Unit-side moments ---
    X_av <- matrix(0, nrow = N, ncol = dim_moment)  # N x dim_moment

    X_mat <- function(X) {
      if (is.null(dim(X))) {
        matrix(X, nrow = N, ncol = T)
      } else {
        X
      }
    }

    for (p in 1:dim_moment) {
      # For each power
      X_pow <- matrix(0, nrow = N, ncol = T)  # N x T

      for (t in 1:T) {
        # Collect all covariates at time t
        X_t <- do.call(cbind, lapply(X_list, function(X) X_mat(X)[, t]))  # N x (K+1)
        # Average the covariates in the p-th group
        if ( p == 1){
          X_pow[, t] <- rowMeans(X_t)             # (1/K) sum_k x_itk^p
        }else if (p == 2){
          X_pow[, t] <- rowMeans(exp(X_t))
        }else if (p == 3){
          X_pow[, t] <- rowMeans((X_t/(1+abs(X_t))))
        }else{
          X_pow[, t] <- rowMeans(X_t^p)
        }
      }

      # Average over time
      X_av[, p] <- rowMeans(X_pow)               # (1/T) sum_t (1/K) sum_k x_itk^p
    }

    mom_i <- cbind(mom_i, X_av)

    mom_micro <- c()
    X_cbind <- do.call(cbind, lapply(X_list, function(X) X_mat(X)))
    mom_micro <- cbind(mom_micro, X_cbind)
    for (p in 1:dim_moment) {
      if ( p == 1){
        mom_micro <- cbind(mom_micro, X_cbind^p)             # (1/K) sum_k x_itk^p
      }else if (p == 2){
        mom_micro <- cbind(mom_micro, exp(X_cbind))
      }else if (p == 3){
        mom_micro <- cbind(mom_micro, (X_cbind/(1+abs(X_cbind))))
      }else{
        mom_micro <- cbind(mom_micro, X_cbind^p)
      }
    }

    ## --- Demean & rescale ---
    for (j in 1:mdim) {
      col_mean <- mean(mom_i[, j])
      col_sd <- sd(mom_i[, j])
      if (col_sd > 0) {
        cols <- ((j - 1) * T * (K+1) + 1):(j * T * (K+1))
        mom_micro[, cols] <- (mom_micro[, cols] - col_mean) / col_sd
        mom_i[, j] <- (mom_i[, j] - col_mean) / col_sd
      }
    }


    # --- Variance / noise on rescaled data ---
    variance <-  sum(sapply(1:mdim, function(j) {
      norm(mom_micro[, ((j - 1) * T * (K+1) + 1):(j * T * (K+1))] - mom_i[,j], type = "F")^2
    })) / (N * T^2 * (K+1)^2)

    data <- mom_i
    dim_size <- N

  }else {
    stop("Invalid type. Use 'long' or 'tall'.")
  }

  ## --- Clustering ---
  if (!is.null(groups)) {
    clusters <- groups
    k_result <- kmeans(data, centers = clusters, algorithm = "Lloyd", nstart = init, iter.max = 100)
  } else {
    if (type == "long" && dim_size <= 1) {
      # With one unit there is no meaningful clustering; force one cluster.
      clusters <- 1
      k_result <- kmeans(data, centers = 1, algorithm = "Lloyd", nstart = init, iter.max = 100)
    } else {
      sd_by_col <- apply(data, 2, sd, na.rm = TRUE)
      sd_by_col[!is.finite(sd_by_col)] <- 0
      has_signal <- max(sd_by_col, na.rm = TRUE) > 0

      if (has_signal) {
        clusters <- 1
        xx = 1000
        threshold <- gamma * max(variance, .Machine$double.eps)
        while (xx >= threshold && clusters < dim_size) {
          k_result <- kmeans(data, centers = clusters, algorithm = "Lloyd", nstart = init, iter.max = 100)
          xx <- k_result$tot.withinss / dim_size
          clusters <- clusters + 1
        }
        if (type == "long") {
          clusters = max(1, min(clusters, dim_size - 1))
        } else {
          clusters = min(clusters, dim_size - 1)
        }
        k_result <- kmeans(data,centers = clusters,algorithm = "Lloyd",nstart = init,iter.max = 100)

      } else {
        clusters <- 1
        k_result <- kmeans(data, centers = 1, algorithm = "Lloyd", nstart = init, iter.max = 100)
      }
    }
  }

  ## --- Return ---
  if (type == "long") {
    list(res = k_result, clusters = clusters, data = data)
  } else {
    list(res = k_result, clusters = clusters, data = data)
  }
}

fix_time_clusters <- function(cluster, data, m = 2) {
  # cluster: raw K-means labels (length T)
  # data: T x 1 matrix (average per time period)
  # m: minimum cluster size
  #
  # Returns:
  #   cluster: cleaned labels for each time period
  #   G_time: number of clusters

  cluster <- as.integer(cluster)
  T <- nrow(data)

  # Collapse to 1 cluster if impossible
  if (T < 2 * m) {
    return(list(
      cluster = rep(1, T),
      G_time = 1
    ))
  }

  repeat {
    sizes <- table(cluster)

    small <- as.integer(names(sizes[sizes < m]))
    large <- as.integer(names(sizes[sizes >= m]))

    if (length(small) == 0) break

    # Compute centers (safe for 1-column data)
    unique_cl <- sort(unique(cluster))
    centers <- do.call(rbind, lapply(unique_cl, function(cl) {
      colMeans(data[cluster == cl, , drop = FALSE])
    }))
    rownames(centers) <- unique_cl

    # Reassign points from small clusters
    for (sc in small) {
      idx <- which(cluster == sc)
      for (i in idx) {
        target <- if (length(large) > 0) large else setdiff(unique_cl, sc)
        dists <- sapply(target, function(cl) {
          sum((data[i, ] - centers[as.character(cl), ])^2)
        })
        cluster[i] <- target[which.min(dists)]
      }
    }
  }

  # Relabel clusters 1:K
  unique_cl <- sort(unique(cluster))
  map <- setNames(seq_along(unique_cl), unique_cl)
  cluster <- map[as.character(cluster)]

  # Number of clusters
  G_time <- length(unique_cl)

  return(list(
    cluster = as.integer(cluster),
    G_time = G_time
  ))
}





#' #' Estimator Function with Optional Cross-Fitting
#' #'@title Inference in high-dimensional data after discretizing unobserved heterogeneity
#' #'
#' #' @description R package 'HDpcluster' is dedicated to do inference for high-dimensional linear panel data model with unkown functions of fixed effects.
#' #'
#' #' @param y outcome variable
#' #' @param D treatment variable
#' #' @param X control variables
#' #' @param T panel size of the time period
#' #' @param groups_init number of the unit clusters
#' #' @param index index name
#' #' @param data data which contains the correct index of the outcome variable or the treatment variable
#' #' @param cluster_type different cluster means for the data, 'unit kmeans' is default, allow for 'unit pesudo'
#' #' @param pesudo_type if cluster_type = 'unit pesudo', choice of the pesudo type
#' #' @param link if cluster_type = 'unit pesudo', choice of different links of pesudo type
#' #' @param optimal_index if cluster_type = 'unit pesudo', different ways to compute the optimal number of clusters
#' #'
#' #' @returns A list of fitted results is returned.
#' #' Within this outputted list, the following elements can be found:
#' #'     \item{res}{regression model.}
#' #'     \item{G}{number of unit clusters.}
#' #'     \item{estimate_correct}{corrected standard error, t value, and p value.}
#' #'     \item{summary_table}{summary of the model together with corrected estimates in the 'Coefficients'.}
#' #'
#' #' @import hdm
#' #' @import plm
#' #' @useDynLib HDpcluster
#' #' @export
#' HP_estimate <- function(data, y_col = NULL, covariate_cols = NULL,
#'                         id_col = "id", time_col = "time", cluster_type = c('two way'),
#'                         index = c("id", "time"), method_auto = "dynamicTreeCut", cluster_method = c('kmeans', 'hierarchical'), pre_cluster = FALSE, Du_pre, Dv_pre, unit_cluster = NULL, time_cluster = NULL, gamma = 1, cc = 0, dim_moment = 1) {
#'
#'
#'   if (is.null(covariate_cols)) {
#'     covariate_cols <- setdiff(colnames(data), c(id_col, time_col, y_col))
#'   }
#'
#'   ids <- sort(unique(data[[id_col]]))
#'   times <- sort(unique(data[[time_col]]))
#'   N <- length(ids)
#'   T <- length(times)
#'   K <- length(covariate_cols)
#'
#'   # Initialize y and X
#'   y <- matrix(NA_real_, nrow = N, ncol = T)
#'   X <- array(NA_real_, dim = c(N, T, K))
#'
#'   for (i in seq_len(nrow(data))) {
#'     id_idx <- which(ids == data[[id_col]][i])
#'     time_idx <- which(times == data[[time_col]][i])
#'     y[id_idx, time_idx] <- data[[y_col]][i]
#'     for (k in seq_along(covariate_cols)) {
#'       X[id_idx, time_idx, k] <- data[[covariate_cols[k]]][i]
#'     }
#'   }
#'
#'   # Cluster (output indicators assumed one-hot encoded)
#'   X_norm <- array(NA_real_, dim = dim(X))  # Same shape: (N, T, K)
#'   for (k in 1:dim(X)[3]) {
#'     X_slice <- X[,,k]
#'     X_norm[,,k] <- scale(X_slice)
#'   }
#'   y_norm = scale(y)
#'
#'
#'   if (cluster_type == 'two way'){
#'   if (cluster_method == 'kmeans'){
#'
#'     # kmeans
#'     X_list = lapply(seq_len(dim(X)[3]), function(k) X[ , , k])
#'     if (pre_cluster == FALSE){
#'       # cluster
#'       if (is.null(unit_cluster) == 1 & is.null(time_cluster) == 1){
#'         clusteri <- cluster_general(y, X_list, N, T, init = 50, type = "long", gamma = gamma, cc = cc, dim_moment = dim_moment  )
#'         G_unit <- clusteri$clusters
#'         klong <- clusteri$res
#'
#'
#'         clustert <- cluster_general(y, X_list, N, T, init = 50, type = "tall", gamma = gamma, cc = cc, dim_moment = dim_moment  )
#'         G_time <- clustert$clusters
#'         ktall <- clustert$res
#'         newcluster_T = fix_time_clusters(cluster = ktall$cluster, data = clustert$data, m = 3)
#'         ktall$cluster =   newcluster_T$cluster
#'         G_time =  newcluster_T$G_time
#'       #   if (G_time > floor(T/2)){
#'       #     clustert <- cluster_general(y, X_list, N, T, init = 50, type = "tall", groups = 1 )
#'       #     G_time <- clustert$clusters
#'       #     ktall <- clustert$res
#'       #   }
#'       }else{
#'         clusteri <- cluster_general(y, X_list, N, T, init = 50, type = "long" , groups = c(floor(unit_cluster)) )
#'         G_unit <- clusteri$clusters
#'         klong <- clusteri$res
#'
#'
#'         clustert <- cluster_general(y, X_list, N, T, init = 50, type = "tall" , groups = c(floor(time_cluster)))
#'         G_time <- clustert$clusters
#'         ktall <- clustert$res
#'       }
#'
#'
#'       Du <- matrix(0, N, G_unit)
#'       Dv <- matrix(0, T, G_time)
#'       Diagu<- matrix(0, G_unit, G_unit)
#'       Diagv<- matrix(0, G_time, G_time)
#'
#'       for (j in seq_len(G_unit)) {
#'         Du[, j] <- as.numeric(klong$cluster == j)
#'         Diagu[j,j] <- 1/sum(Du[,j])
#'       }
#'
#'       for (j in seq_len(G_time)) {
#'         Dv[, j] <- as.numeric(ktall$cluster == j)
#'         Diagv[j,j] <- 1/sum(Dv[,j])
#'       }
#'     }else if (pre_cluster == TRUE){
#'       Du = Du_pre
#'       Dv = Dv_pre
#'       Diagu<- matrix(0, G_unit, G_unit)
#'       Diagv<- matrix(0, G_time, G_time)
#'
#'       for (j in seq_len(G_unit)) {
#'         Diagu[j,j] <- 1/sum(Du[,j])
#'       }
#'
#'       for (j in seq_len(G_time)) {
#'         Diagv[j,j] <- 1/sum(Dv[,j])
#'       }
#'       G_unit = dim(Du)[2]
#'       G_time = dim(Dv)[2]
#'     }
#'   }else if (cluster_method == 'hierarchical'){
#'
#'     # heriachical
#'     if (is.null(unit_cluster) == 1 & is.null(time_cluster) == 1){
#'       res_unit <- cluster_Hierarchical(y_norm, X_norm, type = "unit", method_auto = method_auto)
#'       res_time <- cluster_Hierarchical(y_norm, X_norm, type = "time", method_auto = method_auto)
#'     }else{
#'       res_unit <- cluster_Hierarchical(y_norm, X_norm, type = "unit", cluster = unit_cluster)
#'       res_time <- cluster_Hierarchical(y_norm, X_norm, type = "time", cluster = time_cluster)
#'     }
#'
#'     Du <- res_unit$indicator    # N x G_unit
#'     Dv <- res_time$indicator    # T x G_time
#'     G_unit = res_unit$G
#'     G_time = res_time$G
#'     Diagu<- matrix(0, G_unit, G_unit)
#'     Diagv<- matrix(0, G_time, G_time)
#'
#'     for (j in seq_len(G_unit)) {
#'       Diagu[j,j] <- 1/sum(Du[,j])
#'     }
#'
#'     for (j in seq_len(G_time)) {
#'       Diagv[j,j] <- 1/sum(Dv[,j])
#'     }
#'   }
#'
#'   unit_group <- max.col(Du)
#'   time_group <- max.col(Dv)
#'
#'   # Projection matrices
#'   Mu <- diag(N) - Du %*% Diagu %*% t(Du)
#'   Mv <- diag(T) - Dv %*% Diagv %*% t(Dv)
#'   # Mu <- diag(N) - Du %*% solve(t(Du) %*% Du) %*% t(Du)
#'   # Mv <- diag(T) - Dv %*% solve(t(Dv) %*% Dv) %*% t(Dv)
#'
#'   # Stack y and X into Z: N x T x (K+1)
#'   Z <- array(NA_real_, dim = c(N, T, K + 1))
#'   Z[,,1] <- y
#'   for (k in 1:K) Z[,,k+1] <- X[,,k]
#'
#'   # Projection method: demean unit and time clusters slice-wise
#'   Z_proj <- array(NA_real_, dim = c(N, T, K + 1))
#'   for (k in 1:(K+1)) {
#'     Z_proj[,,k] <- Mu %*% Z[,,k] %*% Mv
#'   }
#'
#'   # Separate transformed y and X
#'   y_trans <- Z_proj[,,1] # Z_proj_final_array[,,1]
#'   X_trans <- Z_proj[,,2:(K+1)] # Z_proj_final_array[,,2:(K+1)]
#'
#'   tY_vector <- as.vector(y_trans)  # (N*T)
#'   tX_matrix <- matrix(aperm(X_trans, c(1, 2, 3)), nrow = N*T, ncol = K)
#'
#'   d_vector <- tX_matrix[, 1]
#'   x_matrix <- tX_matrix[, -1, drop = FALSE]
#'
#'   # Run rlassoEffect with correct inputs
#'   fit <- rlassoEffect(x = x_matrix, y = tY_vector, d = d_vector, method = "double selection")
#'
#'   trans <- data.frame(y = tY_vector,
#'                       D = d_vector,
#'                       x_matrix)
#'
#'   lasso.Y <- rlasso(y ~ . - D - 1, data = trans)
#'   Ytilde <- lasso.Y$residuals
#'
#'   lasso.D <- rlasso(D ~ . - 1, data = trans[, -1])
#'   Dtilde <- lasso.D$residuals
#'
#'   data_res <- data.frame(id = data[[index[1]]],
#'                          time = data[[index[2]]],
#'                          Ytilde = Ytilde,
#'                          Dtilde = Dtilde)
#'
#'   Post_plm <- plm(Ytilde ~ -1 + Dtilde, data = data_res, model = "pooling", index = c("id", "time"))
#'
#'   coefs <- coef(Post_plm)
#'   if (N*T - N*G_time - T*G_unit > 0){
#'   se_corrected <- sqrt(vcovHC(Post_plm, type = "HC0", method = "arellano")) * sqrt(N * T / (N*T - N*G_time - T*G_unit))
#'   }else{
#'     se_corrected <- sqrt(vcovHC(Post_plm, type = "HC0", method = "arellano")) * sqrt(N * T / (1))
#'   }
#'   t_values_corrected <- coefs / se_corrected
#'
#'   df <- Post_plm$df.residual
#'   p_values_corrected <- 2 * pt(-abs(t_values_corrected), df)
#'
#'   summary_table_correct <- data.frame(
#'     Estimate = coefs,
#'     `Std. Error corrected` = se_corrected,
#'     `t-value corrected` = t_values_corrected,
#'     `Pr(>|t|) corrected` = p_values_corrected
#'   )
#'
#'   colnames(summary_table_correct) = c('Estimate', 'Std. Error corrected', 't-value corrected', 'Pr(>|t|) corrected')
#'
#'   # Return cluster counts as well
#'   return(list(
#'     fit_summary = summary(fit),
#'     G_unit = G_unit,
#'     G_time = G_time,
#'     unit_group = unit_group,
#'     time_group = time_group,
#'     post_plm_summary = summary(Post_plm),
#'     estimate_corrected = summary_table_correct,
#'     summary_table = summary_table_correct
#'   ))
#'   }else if (cluster_type == 'one way'){
#'     if (cluster_method == 'kmeans'){
#'
#'       # kmeans
#'       X_list = lapply(seq_len(dim(X)[3]), function(k) X[ , , k])
#'       if (pre_cluster == FALSE){
#'         # cluster
#'         if (is.null(unit_cluster) == 1 & is.null(time_cluster) == 1){
#'           clusteri <- cluster_general(y, X_list, N, T, init = 50, type = "long", gamma = gamma, cc = cc, dim_moment = dim_moment  )
#'           G_unit <- clusteri$clusters
#'           klong <- clusteri$res
#'
#'
#'         }else{
#'           clusteri <- cluster_general(y, X_list, N, T, init = 50, type = "long" , groups = c(floor(unit_cluster)) )
#'           G_unit <- clusteri$clusters
#'           klong <- clusteri$res
#'         }
#'         Du <- matrix(0, N, G_unit)
#'         Diagu<- matrix(0, G_unit, G_unit)
#'
#'
#'         for (j in seq_len(G_unit)) {
#'           Du[, j] <- as.numeric(klong$cluster == j)
#'           Diagu[j,j] <- 1/sum(Du[,j])
#'         }
#'
#'
#'       }else if (pre_cluster == TRUE){
#'         Du = Du_pre
#'         Diagu<- matrix(0, G_unit, G_unit)
#'
#'         for (j in seq_len(G_unit)) {
#'           Diagu[j,j] <- 1/sum(Du[,j])
#'         }
#'
#'         G_unit = dim(Du)[2]
#'       }
#'     }else if (cluster_method == 'hierarchical'){
#'
#'       # heriachical
#'       if (is.null(unit_cluster) == 1 & is.null(time_cluster) == 1){
#'         res_unit <- cluster_Hierarchical(y_norm, X_norm, type = "unit", method_auto = method_auto)
#'       }else{
#'         res_unit <- cluster_Hierarchical(y_norm, X_norm, type = "unit", cluster = unit_cluster)
#'       }
#'
#'       Du <- res_unit$indicator    # N x G_unit
#'       G_unit = res_unit$G
#'       Diagu<- matrix(0, G_unit, G_unit)
#'
#'       for (j in seq_len(G_unit)) {
#'         Diagu[j,j] <- 1/sum(Du[,j])
#'       }
#'
#'     }
#'
#'     unit_group <- max.col(Du)
#'
#'     # Projection matrices
#'     Mu <- diag(N) - Du %*% Diagu %*% t(Du)
#'
#'     # Stack y and X into Z: N x T x (K+1)
#'     Z <- array(NA_real_, dim = c(N, T, K + 1))
#'     Z[,,1] <- y
#'     for (k in 1:K) Z[,,k+1] <- X[,,k]
#'
#'     # Projection method: demean unit and time clusters slice-wise
#'     Z_proj <- array(NA_real_, dim = c(N, T, K + 1))
#'     for (k in 1:(K+1)) {
#'       Z_proj[,,k] <- Mu %*% Z[,,k]
#'     }
#'
#'     # Separate transformed y and X
#'     y_trans <- Z_proj[,,1] # Z_proj_final_array[,,1]
#'     X_trans <- Z_proj[,,2:(K+1)] # Z_proj_final_array[,,2:(K+1)]
#'
#'     tY_vector <- as.vector(y_trans)  # (N*T)
#'     tX_matrix <- matrix(aperm(X_trans, c(1, 2, 3)), nrow = N*T, ncol = K)
#'
#'     # debiased lasso
#'     # debiased_lasso = lasso.proj(tX_matrix, tY_vector) # delas(tX_matrix, tY_vector, c(1))
#'     # debiased_lasso$se = debiased_lasso$se * sqrt(N * T / (N*T - T*G_unit))
#'
#'     #
#'     d_vector <- tX_matrix[, 1]
#'     x_matrix <- tX_matrix[, -1, drop = FALSE]
#'
#'     # Run rlassoEffect with correct inputs
#'     fit <- rlassoEffect(x = x_matrix, y = tY_vector, d = d_vector, method = "double selection")
#'
#'     trans <- data.frame(y = tY_vector,
#'                         D = d_vector,
#'                         x_matrix)
#'
#'     lasso.Y <- rlasso(y ~ . - D - 1, data = trans)
#'     Ytilde <- lasso.Y$residuals
#'
#'     lasso.D <- rlasso(D ~ . - 1, data = trans[, -1])
#'     Dtilde <- lasso.D$residuals
#'
#'     data_res <- data.frame(id = data[[index[1]]],
#'                            time = data[[index[2]]],
#'                            Ytilde = Ytilde,
#'                            Dtilde = Dtilde)
#'
#'     Post_plm <- plm(Ytilde ~ -1 + Dtilde, data = data_res, model = "pooling", index = c("id", "time"))
#'
#'     coefs <- coef(Post_plm)
#'     se_corrected <- sqrt(vcovHC(Post_plm, type = "HC0", method = "arellano")) * sqrt(N * T / (N*T - T*G_unit))
#'     t_values_corrected <- coefs / se_corrected
#'
#'     df <- Post_plm$df.residual
#'     p_values_corrected <- 2 * pt(-abs(t_values_corrected), df)
#'
#'     summary_table_correct <- data.frame(
#'       Estimate = coefs,
#'       `Std. Error corrected` = se_corrected,
#'       `t-value corrected` = t_values_corrected,
#'       `Pr(>|t|) corrected` = p_values_corrected
#'     )
#'
#'     colnames(summary_table_correct) = c('Estimate', 'Std. Error corrected', 't-value corrected', 'Pr(>|t|) corrected')
#'
#'     # Return cluster counts as well
#'     return(list(
#'       fit_summary = summary(fit),
#'       G_unit = G_unit,
#'       unit_group = unit_group,
#'       post_plm_summary = summary(Post_plm),
#'       estimate_corrected = summary_table_correct,
#'       summary_table = summary_table_correct
#'     ))
#'
#'   }else if (cluster_type == 'one way T moment'){
#'     if (cluster_method == 'kmeans'){
#'
#'       # kmeans
#'       X_list = lapply(seq_len(dim(X)[3]), function(k) X[ , , k])
#'       if (pre_cluster == FALSE){
#'         # cluster
#'         if (is.null(unit_cluster) == 1 & is.null(time_cluster) == 1){
#'           clusteri <- cluster_general(y, X_list, N, T, init = 50, type = "long_T_moment", gamma = gamma, cc = cc, dim_moment = dim_moment  )
#'           G_unit <- clusteri$clusters
#'           klong <- clusteri$res
#'
#'
#'         }else{
#'           clusteri <- cluster_general(y, X_list, N, T, init = 50, type = "long_T_moment" , groups = c(floor(unit_cluster)) )
#'           G_unit <- clusteri$clusters
#'           klong <- clusteri$res
#'         }
#'         Du <- matrix(0, N, G_unit)
#'         Diagu<- matrix(0, G_unit, G_unit)
#'
#'
#'         for (j in seq_len(G_unit)) {
#'           Du[, j] <- as.numeric(klong$cluster == j)
#'           Diagu[j,j] <- 1/sum(Du[,j])
#'         }
#'
#'
#'       }else if (pre_cluster == TRUE){
#'         Du = Du_pre
#'         Diagu<- matrix(0, G_unit, G_unit)
#'
#'         for (j in seq_len(G_unit)) {
#'           Diagu[j,j] <- 1/sum(Du[,j])
#'         }
#'
#'         G_unit = dim(Du)[2]
#'       }
#'     }else if (cluster_method == 'hierarchical'){
#'
#'       # heriachical
#'       if (is.null(unit_cluster) == 1 & is.null(time_cluster) == 1){
#'         res_unit <- cluster_Hierarchical(y_norm, X_norm, type = "unit", method_auto = method_auto)
#'       }else{
#'         res_unit <- cluster_Hierarchical(y_norm, X_norm, type = "unit", cluster = unit_cluster)
#'       }
#'
#'       Du <- res_unit$indicator    # N x G_unit
#'       G_unit = res_unit$G
#'       Diagu<- matrix(0, G_unit, G_unit)
#'
#'       for (j in seq_len(G_unit)) {
#'         Diagu[j,j] <- 1/sum(Du[,j])
#'       }
#'
#'     }
#'
#'     unit_group <- max.col(Du)
#'
#'     # Projection matrices
#'     Mu <- diag(N) - Du %*% Diagu %*% t(Du)
#'
#'     # Stack y and X into Z: N x T x (K+1)
#'     Z <- array(NA_real_, dim = c(N, T, K + 1))
#'     Z[,,1] <- y
#'     for (k in 1:K) Z[,,k+1] <- X[,,k]
#'
#'     # Projection method: demean unit and time clusters slice-wise
#'     Z_proj <- array(NA_real_, dim = c(N, T, K + 1))
#'     for (k in 1:(K+1)) {
#'       Z_proj[,,k] <- Mu %*% Z[,,k]
#'     }
#'
#'     # Separate transformed y and X
#'     y_trans <- Z_proj[,,1] # Z_proj_final_array[,,1]
#'     X_trans <- Z_proj[,,2:(K+1)] # Z_proj_final_array[,,2:(K+1)]
#'
#'     tY_vector <- as.vector(y_trans)  # (N*T)
#'     tX_matrix <- matrix(aperm(X_trans, c(1, 2, 3)), nrow = N*T, ncol = K)
#'
#'     # debiased lasso
#'     # debiased_lasso = lasso.proj(tX_matrix, tY_vector) # delas(tX_matrix, tY_vector, c(1))
#'     # debiased_lasso$se = debiased_lasso$se * sqrt(N * T / (N*T - T*G_unit))
#'
#'     #
#'     d_vector <- tX_matrix[, 1]
#'     x_matrix <- tX_matrix[, -1, drop = FALSE]
#'
#'     # Run rlassoEffect with correct inputs
#'     fit <- rlassoEffect(x = x_matrix, y = tY_vector, d = d_vector, method = "double selection")
#'
#'     trans <- data.frame(y = tY_vector,
#'                         D = d_vector,
#'                         x_matrix)
#'
#'     lasso.Y <- rlasso(y ~ . - D - 1, data = trans)
#'     Ytilde <- lasso.Y$residuals
#'
#'     lasso.D <- rlasso(D ~ . - 1, data = trans[, -1])
#'     Dtilde <- lasso.D$residuals
#'
#'     data_res <- data.frame(id = data[[index[1]]],
#'                            time = data[[index[2]]],
#'                            Ytilde = Ytilde,
#'                            Dtilde = Dtilde)
#'
#'     Post_plm <- plm(Ytilde ~ -1 + Dtilde, data = data_res, model = "pooling", index = c("id", "time"))
#'
#'     coefs <- coef(Post_plm)
#'     se_corrected <- sqrt(vcovHC(Post_plm, type = "HC0", method = "arellano")) * sqrt(N * T / (N*T - T*G_unit))
#'     t_values_corrected <- coefs / se_corrected
#'
#'     df <- Post_plm$df.residual
#'     p_values_corrected <- 2 * pt(-abs(t_values_corrected), df)
#'
#'     summary_table_correct <- data.frame(
#'       Estimate = coefs,
#'       `Std. Error corrected` = se_corrected,
#'       `t-value corrected` = t_values_corrected,
#'       `Pr(>|t|) corrected` = p_values_corrected
#'     )
#'
#'     colnames(summary_table_correct) = c('Estimate', 'Std. Error corrected', 't-value corrected', 'Pr(>|t|) corrected')
#'
#'     # Return cluster counts as well
#'     return(list(
#'       fit_summary = summary(fit),
#'       G_unit = G_unit,
#'       unit_group = unit_group,
#'       post_plm_summary = summary(Post_plm),
#'       estimate_corrected = summary_table_correct,
#'       summary_table = summary_table_correct
#'     ))
#'
#'   }
#' }
#'
#'
#'
#' # HP_estimate_double <- function(data, y_col = NULL, covariate_cols = NULL,
#' #                                id_col = "id", time_col = "time",
#' #                                index = c("id", "time"), method_auto = "dynamicTreeCut", cluster_method = 'kmeans', pre_cluster = FALSE, Du_pre, Dv_pre, unit_cluster = NULL, time_cluster = NULL, gamma = 1, cc = 0, dim_moment = 1) {
#' #
#' #
#' #   if (is.null(covariate_cols)) {
#' #     covariate_cols <- setdiff(colnames(data), c(id_col, time_col, y_col))
#' #   }
#' #
#' #   ids <- sort(unique(data[[id_col]]))
#' #   times <- sort(unique(data[[time_col]]))
#' #   N <- length(ids)
#' #   T <- length(times)
#' #   K <- length(covariate_cols)
#' #
#' #   # Initialize y and X
#' #   y <- matrix(NA_real_, nrow = N, ncol = T)
#' #   X <- array(NA_real_, dim = c(N, T, K))
#' #
#' #   for (i in seq_len(nrow(data))) {
#' #     id_idx <- which(ids == data[[id_col]][i])
#' #     time_idx <- which(times == data[[time_col]][i])
#' #     y[id_idx, time_idx] <- data[[y_col]][i]
#' #     for (k in seq_along(covariate_cols)) {
#' #       X[id_idx, time_idx, k] <- data[[covariate_cols[k]]][i]
#' #     }
#' #   }
#' #
#' #   # Cluster (output indicators assumed one-hot encoded)
#' #   X_norm <- array(NA_real_, dim = dim(X))  # Same shape: (N, T, K)
#' #   for (k in 1:dim(X)[3]) {
#' #     X_slice <- X[,,k]
#' #     X_norm[,,k] <- scale(X_slice)
#' #   }
#' #   y_norm = scale(y)
#' #   if (cluster_method == 'kmeans'){
#' #
#' #     # kmeans
#' #     X_list = lapply(seq_len(dim(X)[3]), function(k) X[ , , k])
#' #     if (pre_cluster == FALSE){
#' #       # cluster
#' #       if (is.null(unit_cluster) == 1 & is.null(time_cluster) == 1){
#' #         clusteri <- cluster_general(y, X_list, N, T, init = 50, type = "long", gamma = gamma, cc = cc, dim_moment = dim_moment   )
#' #         G_unit <- clusteri$clusters
#' #         klong <- clusteri$res
#' #
#' #
#' #         clustert <- cluster_general(y, X_list, N, T, init = 50, type = "tall", gamma = gamma, cc = cc, dim_moment = dim_moment  )
#' #         G_time <- clustert$clusters
#' #         ktall <- clustert$res
#' #       }else{
#' #         clusteri <- cluster_general(y, X_list, N, T, init = 50, type = "long" , groups = c(floor(unit_cluster)) )
#' #         G_unit <- clusteri$clusters
#' #         klong <- clusteri$res
#' #
#' #
#' #         clustert <- cluster_general(y, X_list, N, T, init = 50, type = "tall" , groups = c(floor(time_cluster)))
#' #         G_time <- clustert$clusters
#' #         ktall <- clustert$res
#' #       }
#' #
#' #       Du <- matrix(0, N, G_unit)
#' #       Dv <- matrix(0, T, G_time)
#' #       Diagu<- matrix(0, G_unit, G_unit)
#' #       Diagv<- matrix(0, G_time, G_time)
#' #
#' #       for (j in seq_len(G_unit)) {
#' #         Du[, j] <- as.numeric(klong$cluster == j)
#' #         Diagu[j,j] <- 1/sum(Du[,j])
#' #       }
#' #
#' #       for (j in seq_len(G_time)) {
#' #         Dv[, j] <- as.numeric(ktall$cluster == j)
#' #         Diagv[j,j] <- 1/sum(Dv[,j])
#' #       }
#' #     }else if (pre_cluster == TRUE){
#' #       Du = Du_pre
#' #       Dv = Dv_pre
#' #       Diagu<- matrix(0, G_unit, G_unit)
#' #       Diagv<- matrix(0, G_time, G_time)
#' #
#' #       for (j in seq_len(G_unit)) {
#' #         Diagu[j,j] <- 1/sum(Du[,j])
#' #       }
#' #
#' #       for (j in seq_len(G_time)) {
#' #         Diagv[j,j] <- 1/sum(Dv[,j])
#' #       }
#' #       G_unit = dim(Du)[2]
#' #       G_time = dim(Dv)[2]
#' #     }
#' #   }else if (cluster_method == 'hierarchical'){
#' #
#' #     # heriachical
#' #     if (is.null(unit_cluster) == 1 & is.null(time_cluster) == 1){
#' #       res_unit <- cluster_Hierarchical(y_norm, X_norm, type = "unit", method_auto = method_auto)
#' #       res_time <- cluster_Hierarchical(y_norm, X_norm, type = "time", method_auto = method_auto)
#' #     }else{
#' #       res_unit <- cluster_Hierarchical(y_norm, X_norm, type = "unit", cluster = unit_cluster)
#' #       res_time <- cluster_Hierarchical(y_norm, X_norm, type = "time", cluster = time_cluster)
#' #     }
#' #
#' #     Du <- res_unit$indicator    # N x G_unit
#' #     Dv <- res_time$indicator    # T x G_time
#' #
#' #     G_unit = res_unit$G
#' #     G_time = res_time$G
#' #
#' #     Diagu<- matrix(0, G_unit, G_unit)
#' #     Diagv<- matrix(0, G_time, G_time)
#' #
#' #     for (j in seq_len(G_unit)) {
#' #       Diagu[j,j] <- 1/sum(Du[,j])
#' #     }
#' #
#' #     for (j in seq_len(G_time)) {
#' #       Diagv[j,j] <- 1/sum(Dv[,j])
#' #     }
#' #   }
#' #
#' #   unit_group <- max.col(Du)
#' #   time_group <- max.col(Dv)
#' #
#' #   # Projection matrices
#' #   Mu <- diag(N) - Du %*% Diagu %*% t(Du)
#' #   Mv <- diag(T) - Dv %*% Diagv %*% t(Dv)
#' #   # Mu <- diag(N) - Du %*% solve(t(Du) %*% Du) %*% t(Du)
#' #   # Mv <- diag(T) - Dv %*% solve(t(Dv) %*% Dv) %*% t(Dv)
#' #
#' #   # Stack y and X into Z: N x T x (K+1)
#' #   Z <- array(NA_real_, dim = c(N, T, K + 1))
#' #   Z[,,1] <- y
#' #   for (k in 1:K) Z[,,k+1] <- X[,,k]
#' #
#' #   # Projection method: demean unit and time clusters slice-wise
#' #   Z_proj <- array(NA_real_, dim = c(N, T, K + 1))
#' #   for (k in 1:(K+1)) {
#' #     Z_proj[,,k] <- Mu %*% Z[,,k] %*% Mv
#' #   }
#' #
#' #
#' #   y_trans <- Z_proj[,,1] # Z_proj_final_array[,,1]
#' #   X_trans <- Z_proj[,,2:(K+1)] # Z_proj_final_array[,,2:(K+1)]
#' #
#' #   tY_vector <- as.vector(y_trans)  # (N*T)
#' #   tX_matrix <- matrix(aperm(X_trans, c(1, 2, 3)), nrow = N*T, ncol = K)
#' #
#' #   d_vector <- tX_matrix[, 1]
#' #   x_matrix <- tX_matrix[, -1, drop = FALSE]
#' #   # Separate transformed y and X
#' #
#' #   trans <- data.frame(y = tY_vector,
#' #                       D = d_vector,
#' #                       x_matrix)
#' #
#' #   # Run rlassoEffect with correct inputs
#' #   fit <- rlassoEffect(x = x_matrix, y = tY_vector, d = d_vector, method = "double selection")
#' #
#' #   trans <- data.frame(y = tY_vector,
#' #                       D = d_vector,
#' #                       x_matrix)
#' #
#' #   lasso.Y <- rlasso(y ~ . - D - 1, data = trans)
#' #   Ytilde <- lasso.Y$residuals
#' #
#' #   lasso.D <- rlasso(D ~ . - 1, data = trans[, -1])
#' #   Dtilde <- lasso.D$residuals
#' #
#' #   data_res <- data.frame(id = data[[index[1]]],
#' #                          time = data[[index[2]]],
#' #                          Ytilde = Ytilde,
#' #                          Dtilde = Dtilde)
#' #
#' #   Post_plm <- plm(Ytilde ~ -1 + Dtilde, data = data_res, model = "pooling", index = c("id", "time"))
#' #
#' #   coefs <- coef(Post_plm)
#' #   se_corrected <- sqrt(vcovHC(Post_plm, type = "HC0", method = "arellano")) * sqrt(N * T / (N*T - N*G_time - T*G_unit))
#' #   t_values_corrected <- coefs / se_corrected
#' #
#' #   df <- Post_plm$df.residual
#' #   p_values_corrected <- 2 * pt(-abs(t_values_corrected), df)
#' #
#' #   summary_table_correct <- data.frame(
#' #     Estimate = coefs,
#' #     `Std. Error corrected` = se_corrected,
#' #     `t-value corrected` = t_values_corrected,
#' #     `Pr(>|t|) corrected` = p_values_corrected
#' #   )
#' #
#' #   colnames(summary_table_correct) = c('Estimate', 'Std. Error corrected', 't-value corrected', 'Pr(>|t|) corrected')
#' #
#' #   # Return cluster counts as well
#' #   return(list(
#' #     fit_summary = summary(fit),
#' #     G_unit = G_unit,
#' #     G_time = G_time,
#' #     unit_group = unit_group,
#' #     time_group = time_group,
#' #     post_plm_summary = summary(Post_plm),
#' #     estimate_corrected = summary_table_correct,
#' #     summary_table = summary_table_correct
#' #   ))
#' # }
#'
#' #' @export
#' Niave_estimate <- function(data, y_col = NULL, covariate_cols = NULL,
#'                            id_col = "id", time_col = "time", index = c("id", "time")) {
#'
#'
#'   if (is.null(covariate_cols)) {
#'     covariate_cols <- setdiff(colnames(data), c(id_col, time_col, y_col))
#'   }
#'
#'   ids <- sort(unique(data[[id_col]]))
#'   times <- sort(unique(data[[time_col]]))
#'   N <- length(ids)
#'   T <- length(times)
#'   K <- length(covariate_cols)
#'
#'   # Initialize y and X
#'   y <- matrix(NA_real_, nrow = N, ncol = T)
#'   X <- array(NA_real_, dim = c(N, T, K))
#'
#'   for (i in seq_len(nrow(data))) {
#'     id_idx <- which(ids == data[[id_col]][i])
#'     time_idx <- which(times == data[[time_col]][i])
#'     y[id_idx, time_idx] <- data[[y_col]][i]
#'     for (k in seq_along(covariate_cols)) {
#'       X[id_idx, time_idx, k] <- data[[covariate_cols[k]]][i]
#'     }
#'   }
#'
#'   # Stack y and X into Z: N x T x (K+1)
#'   Z <- array(NA_real_, dim = c(N, T, K + 1))
#'   Z[,,1] <- y
#'   for (k in 1:K) Z[,,k+1] <- X[,,k]
#'
#'   # Projection method: demean unit and time clusters slice-wise
#'   Z_proj <- array(NA_real_, dim = c(N, T, K + 1))
#'   for (k in 1:(K+1)) {
#'     Z_proj[,,k] <-  Z[,,k]
#'   }
#'
#'   # Separate transformed y and X
#'   y_trans <- Z_proj[,,1] # Z_proj_final_array[,,1]
#'   X_trans <- Z_proj[,,2:(K+1)] # Z_proj_final_array[,,2:(K+1)]
#'
#'   tY_vector <- as.vector(y_trans)  # (N*T)
#'   tX_matrix <- matrix(aperm(X_trans, c(1, 2, 3)), nrow = N*T, ncol = K)
#'
#'   d_vector <- tX_matrix[, 1]
#'   x_matrix <- tX_matrix[, -1, drop = FALSE]
#'
#'   # Run rlassoEffect with correct inputs
#'   fit <- rlassoEffect(x = x_matrix, y = tY_vector, d = d_vector, method = "double selection")
#'
#'   trans <- data.frame(y = tY_vector,
#'                       D = d_vector,
#'                       x_matrix)
#'
#'   lasso.Y <- rlasso(y ~ . - D - 1, data = trans)
#'   Ytilde <- lasso.Y$residuals
#'
#'   lasso.D <- rlasso(D ~ . - 1, data = trans[, -1])
#'   Dtilde <- lasso.D$residuals
#'
#'   data_res <- data.frame(id = data[[index[1]]],
#'                          time = data[[index[2]]],
#'                          Ytilde = Ytilde,
#'                          Dtilde = Dtilde)
#'
#'   Post_plm <- plm(Ytilde ~ -1 + Dtilde, data = data_res, model = "pooling", index = c("id", "time"))
#'
#'   coefs <- coef(Post_plm)
#'   se_corrected <- sqrt(vcovHC(Post_plm, type = "HC0", method = "arellano"))
#'   t_values_corrected <- coefs / se_corrected
#'
#'   df <- Post_plm$df.residual
#'   p_values_corrected <- 2 * pt(-abs(t_values_corrected), df)
#'
#'   summary_table_correct <- data.frame(
#'     Estimate = coefs,
#'     `Std. Error corrected` = se_corrected,
#'     `t-value corrected` = t_values_corrected,
#'     `Pr(>|t|) corrected` = p_values_corrected
#'   )
#'
#'   colnames(summary_table_correct) = c('Estimate', 'Std. Error corrected', 't-value corrected', 'Pr(>|t|) corrected')
#'
#'   # Return cluster counts as well
#'   return(list(
#'     fit_summary = summary(fit),
#'     post_plm_summary = summary(Post_plm),
#'     estimate_corrected = summary_table_correct,
#'     summary_table = summary_table_correct
#'   ))
#' }
#'
#' #' @export
#' TWFE_estimate <- function(data, y_col = NULL, covariate_cols = NULL,
#'                           id_col = "id", time_col = "time", index = c("id", "time")) {
#'
#'
#'   if (is.null(covariate_cols)) {
#'     covariate_cols <- setdiff(colnames(data), c(id_col, time_col, y_col))
#'   }
#'
#'   ids <- sort(unique(data[[id_col]]))
#'   times <- sort(unique(data[[time_col]]))
#'   N <- length(ids)
#'   T <- length(times)
#'   K <- length(covariate_cols)
#'
#'   # Initialize y and X
#'   y <- matrix(NA_real_, nrow = N, ncol = T)
#'   X <- array(NA_real_, dim = c(N, T, K))
#'
#'   for (i in seq_len(nrow(data))) {
#'     id_idx <- which(ids == data[[id_col]][i])
#'     time_idx <- which(times == data[[time_col]][i])
#'     y[id_idx, time_idx] <- data[[y_col]][i]
#'     for (k in seq_along(covariate_cols)) {
#'       X[id_idx, time_idx, k] <- data[[covariate_cols[k]]][i]
#'     }
#'   }
#'
#'   # Cluster (output indicators assumed one-hot encoded)
#'   X_norm <- array(NA_real_, dim = dim(X))  # Same shape: (N, T, K)
#'   for (k in 1:dim(X)[3]) {
#'     X_slice <- X[,,k]
#'     X_norm[,,k] <- scale(X_slice)
#'   }
#'   y_norm = scale(y)
#'
#'   Du <- matrix(1, nrow = N, ncol = 1)    # N x G_unit
#'   Dv <- matrix(1, nrow = T, ncol = 1)    # T x G_time
#'
#'   G_unit = 1
#'   G_time = 1
#'
#'   unit_group <- max.col(Du)
#'   time_group <- max.col(Dv)
#'
#'   # Projection matrices
#'   Mu <- diag(N) - Du %*% solve(t(Du) %*% Du) %*% t(Du)
#'   Mv <- diag(T) - Dv %*% solve(t(Dv) %*% Dv) %*% t(Dv)
#'
#'   # Stack y and X into Z: N x T x (K+1)
#'   Z <- array(NA_real_, dim = c(N, T, K + 1))
#'   Z[,,1] <- y
#'   for (k in 1:K) Z[,,k+1] <- X[,,k]
#'
#'   # Projection method: demean unit and time clusters slice-wise
#'   Z_proj <- array(NA_real_, dim = c(N, T, K + 1))
#'   for (k in 1:(K+1)) {
#'     Z_proj[,,k] <- Mu %*% Z[,,k] %*% Mv
#'   }
#'
#'   # Separate transformed y and X
#'   y_trans <- Z_proj[,,1] # Z_proj_final_array[,,1]
#'   X_trans <- Z_proj[,,2:(K+1)] # Z_proj_final_array[,,2:(K+1)]
#'
#'   tY_vector <- as.vector(y_trans)  # (N*T)
#'   tX_matrix <- matrix(aperm(X_trans, c(1, 2, 3)), nrow = N*T, ncol = K)
#'
#'   d_vector <- tX_matrix[, 1]
#'   x_matrix <- tX_matrix[, -1, drop = FALSE]
#'
#'   # Run rlassoEffect with correct inputs
#'   fit <- rlassoEffect(x = x_matrix, y = tY_vector, d = d_vector, method = "double selection")
#'
#'   trans <- data.frame(y = tY_vector,
#'                       D = d_vector,
#'                       x_matrix)
#'
#'   lasso.Y <- rlasso(y ~ . - D - 1, data = trans)
#'   Ytilde <- lasso.Y$residuals
#'
#'   lasso.D <- rlasso(D ~ . - 1, data = trans[, -1])
#'   Dtilde <- lasso.D$residuals
#'
#'   data_res <- data.frame(id = data[[index[1]]],
#'                          time = data[[index[2]]],
#'                          Ytilde = Ytilde,
#'                          Dtilde = Dtilde)
#'
#'   Post_plm <- plm(Ytilde ~ -1 + Dtilde, data = data_res, model = "pooling", index = c("id", "time"))
#'
#'   coefs <- coef(Post_plm)
#'   se_corrected <- sqrt(vcovHC(Post_plm, type = "HC0", method = "arellano")) * sqrt(N * T / (N*T - N*G_time - T*G_unit))
#'   t_values_corrected <- coefs / se_corrected
#'
#'   df <- Post_plm$df.residual
#'   p_values_corrected <- 2 * pt(-abs(t_values_corrected), df)
#'
#'   summary_table_correct <- data.frame(
#'     Estimate = coefs,
#'     `Std. Error corrected` = se_corrected,
#'     `t-value corrected` = t_values_corrected,
#'     `Pr(>|t|) corrected` = p_values_corrected
#'   )
#'
#'   colnames(summary_table_correct) = c('Estimate', 'Std. Error corrected', 't-value corrected', 'Pr(>|t|) corrected')
#'
#'   # Return cluster counts as well
#'   return(list(
#'     fit_summary = summary(fit),
#'     G_unit = G_unit,
#'     G_time = G_time,
#'     unit_group = unit_group,
#'     time_group = time_group,
#'     post_plm_summary = summary(Post_plm),
#'     estimate_corrected = summary_table_correct,
#'     summary_table = summary_table_correct
#'   ))
#' }
#'
#' # HP_estimate_twoway <- function(data, y_col = NULL, covariate_cols = NULL,
#' #                         id_col = "id", time_col = "time",
#' #                         index = c("id", "time"), method_auto = "dynamicTreeCut", cluster = NULL) {
#' #
#' #
#' #   if (is.null(covariate_cols)) {
#' #     covariate_cols <- setdiff(colnames(data), c(id_col, time_col, y_col))
#' #   }
#' #
#' #   ids <- sort(unique(data[[id_col]]))
#' #   times <- sort(unique(data[[time_col]]))
#' #   N <- length(ids)
#' #   T <- length(times)
#' #   K <- length(covariate_cols)
#' #
#' #   # Initialize y and X
#' #   y <- matrix(NA_real_, nrow = N, ncol = T)
#' #   X <- array(NA_real_, dim = c(N, T, K))
#' #
#' #   for (i in seq_len(nrow(data))) {
#' #     id_idx <- which(ids == data[[id_col]][i])
#' #     time_idx <- which(times == data[[time_col]][i])
#' #     y[id_idx, time_idx] <- data[[y_col]][i]
#' #     for (k in seq_along(covariate_cols)) {
#' #       X[id_idx, time_idx, k] <- data[[covariate_cols[k]]][i]
#' #     }
#' #   }
#' #
#' #   # Cluster (output indicators assumed one-hot encoded)
#' #   X_norm <- array(NA_real_, dim = dim(X))  # Same shape: (N, T, K)
#' #   for (k in 1:dim(X)[3]) {
#' #     X_slice <- X[,,k]
#' #     X_norm[,,k] <- scale(X_slice)
#' #   }
#' #   y_norm = scale(y)
#' #   res_unit <- cluster_Hierarchical(y_norm, X_norm, type = "unit", method_auto = method_auto)
#' #   res_time <- cluster_Hierarchical(y_norm, X_norm, type = "time", method_auto = method_auto)
#' #   res_covar <- cluster_Hierarchical(y_norm, X_norm, type = "covariate", method_auto = method_auto, cluster = cluster)
#' #
#' #   Dv <- res_unit$indicator    # N x G_unit
#' #   Du <- res_time$indicator    # T x G_time
#' #   Dc <- res_covar$indicator   # (K+1) x G_covar
#' #
#' #   G_unit = res_unit$G
#' #   G_time = res_time$G
#' #   G_covar = res_covar$G
#' #
#' #   unit_group <- max.col(Dv)
#' #   time_group <- max.col(Du)
#' #   covar_group <- max.col(Dc)
#' #
#' #   # Projection matrices
#' #   Mv <- diag(N) - Dv %*% solve(t(Dv) %*% Dv) %*% t(Dv)
#' #   Mu <- diag(T) - Du %*% solve(t(Du) %*% Du) %*% t(Du)
#' #   Mc <- diag(K+1) - Dc %*% solve(t(Dc) %*% Dc) %*% t(Dc)
#' #
#' #   # Stack y and X into Z: N x T x (K+1)
#' #   Z <- array(NA_real_, dim = c(N, T, K + 1))
#' #   Z[,,1] <- y
#' #   for (k in 1:K) Z[,,k+1] <- X[,,k]
#' #
#' #   # Projection method: demean unit and time clusters slice-wise
#' #   Z_proj <- array(NA_real_, dim = c(N, T, K + 1))
#' #   for (k in 1:(K+1)) {
#' #     Z_proj[,,k] <- Mv %*% Z[,,k] %*% Mu
#' #   }
#' #   # Separate transformed y and X
#' #   y_trans <- Z_proj[,,1] # Z_proj_final_array[,,1]
#' #   X_trans <- Z_proj[,,2:(K+1)] # Z_proj_final_array[,,2:(K+1)]
#' #
#' #   tY_vector <- as.vector(y_trans)  # (N*T)
#' #   tX_matrix <- matrix(aperm(X_trans, c(1, 2, 3)), nrow = N*T, ncol = K)
#' #
#' #   d_vector <- tX_matrix[, 1]
#' #   x_matrix <- tX_matrix[, -1, drop = FALSE]
#' #
#' #   # Run rlassoEffect with correct inputs
#' #   fit <- rlassoEffect(x = x_matrix, y = tY_vector, d = d_vector, method = "partialling out")
#' #
#' #
#' #   trans <- data.frame(y = tY_vector,
#' #                       D = d_vector,
#' #                       x_matrix)
#' #
#' #   lasso.Y <- rlasso(y ~ . - D - 1, data = trans)
#' #   Ytilde <- lasso.Y$residuals
#' #
#' #   lasso.D <- rlasso(D ~ . - 1, data = trans[, -1])
#' #   Dtilde <- lasso.D$residuals
#' #
#' #   data_res <- data.frame(id = data[[index[1]]],
#' #                          time = data[[index[2]]],
#' #                          Ytilde = Ytilde,
#' #                          Dtilde = Dtilde)
#' #
#' #   Post_plm <- plm(Ytilde ~ -1 + Dtilde, data = data_res, model = "pooling", index = c("id", "time"))
#' #
#' #   coefs <- coef(Post_plm)
#' #   se_corrected <- sqrt(vcovHC(Post_plm, type = "HC0", method = "arellano")) * sqrt(N * T / ((N-G_unit)*(T-G_time)))
#' #   t_values_corrected <- coefs / se_corrected
#' #
#' #   df <- Post_plm$df.residual
#' #   p_values_corrected <- 2 * pt(-abs(t_values_corrected), df)
#' #
#' #   summary_table_correct <- data.frame(
#' #     Estimate = coefs,
#' #     `Std. Error corrected` = se_corrected,
#' #     `t-value corrected` = t_values_corrected,
#' #     `Pr(>|t|) corrected` = p_values_corrected
#' #   )
#' #
#' #   summary_table <- summary(Post_plm)
#' #   summary_table$coefficients <- cbind(summary_table_correct, summary_table$coefficients)
#' #
#' #   # Return cluster counts as well
#' #   return(list(
#' #     fit_summary = summary(fit),
#' #     G_unit = G_unit,
#' #     G_time = G_time,
#' #     unit_group = unit_group,
#' #     time_group = time_group,
#' #     post_plm_summary = summary(Post_plm),
#' #     estimate_corrected = summary_table_correct,
#' #     summary_table = summary_table
#' #   ))
#' # }
#' #
#' # HP_estimate_double_twoway <- function(data, y_col = NULL, covariate_cols = NULL,
#' #                                id_col = "id", time_col = "time",
#' #                                index = c("id", "time"), method_auto = "dynamicTreeCut", cluster = NULL) {
#' #
#' #
#' #   if (is.null(covariate_cols)) {
#' #     covariate_cols <- setdiff(colnames(data), c(id_col, time_col, y_col))
#' #   }
#' #
#' #   ids <- sort(unique(data[[id_col]]))
#' #   times <- sort(unique(data[[time_col]]))
#' #   N <- length(ids)
#' #   T <- length(times)
#' #   K <- length(covariate_cols)
#' #
#' #   # Initialize y and X
#' #   y <- matrix(NA_real_, nrow = N, ncol = T)
#' #   X <- array(NA_real_, dim = c(N, T, K))
#' #
#' #   for (i in seq_len(nrow(data))) {
#' #     id_idx <- which(ids == data[[id_col]][i])
#' #     time_idx <- which(times == data[[time_col]][i])
#' #     y[id_idx, time_idx] <- data[[y_col]][i]
#' #     for (k in seq_along(covariate_cols)) {
#' #       X[id_idx, time_idx, k] <- data[[covariate_cols[k]]][i]
#' #     }
#' #   }
#' #
#' #   # Cluster (output indicators assumed one-hot encoded)
#' #   X_norm <- array(NA_real_, dim = dim(X))  # Same shape: (N, T, K)
#' #   for (k in 1:dim(X)[3]) {
#' #     X_slice <- X[,,k]
#' #     X_norm[,,k] <- scale(X_slice)
#' #   }
#' #   y_norm = scale(y)
#' #   res_unit <- cluster_Hierarchical(y_norm, X_norm, type = "unit", method_auto = method_auto)
#' #   res_time <- cluster_Hierarchical(y_norm, X_norm, type = "time", method_auto = method_auto)
#' #   res_covar <- cluster_Hierarchical(y_norm, X_norm, type = "covariate", method_auto = method_auto, cluster = cluster)
#' #
#' #   Dv <- res_unit$indicator    # N x G_unit
#' #   Du <- res_time$indicator    # T x G_time
#' #   Dc <- res_covar$indicator   # (K+1) x G_covar
#' #
#' #   G_unit = res_unit$G
#' #   G_time = res_time$G
#' #   G_covar = res_covar$G
#' #
#' #   unit_group <- max.col(Dv)
#' #   time_group <- max.col(Du)
#' #   covar_group <- max.col(Dc)
#' #
#' #   # Projection matrices
#' #   Mv <- diag(N) - Dv %*% solve(t(Dv) %*% Dv) %*% t(Dv)
#' #   Mu <- diag(T) - Du %*% solve(t(Du) %*% Du) %*% t(Du)
#' #   Mc <- diag(K+1) - Dc %*% solve(t(Dc) %*% Dc) %*% t(Dc)
#' #
#' #   # Stack y and X into Z: N x T x (K+1)
#' #   Z <- array(NA_real_, dim = c(N, T, K + 1))
#' #   Z[,,1] <- y
#' #   for (k in 1:K) Z[,,k+1] <- X[,,k]
#' #
#' #   # Projection method: demean unit and time clusters slice-wise
#' #   Z_proj <- array(NA_real_, dim = c(N, T, K + 1))
#' #   for (k in 1:(K+1)) {
#' #     Z_proj[,,k] <- Mv %*% Z[,,k] %*% Mu
#' #   }
#' #
#' #   # Flatten Z_proj for covariate demeaning
#' #   Z_proj_mat <- matrix(NA_real_, nrow = N * T, ncol = K + 1)
#' #   for (k in 1:(K+1)) {
#' #     Z_proj_mat[,k] <- as.vector(t(Z_proj[,,k]))
#' #   }
#' #
#' #   # Covariate demeaning via projection
#' #   Z_proj_final <- Z_proj_mat %*% Mc
#' #
#' #   # Reshape back to array
#' #   Z_proj_final_array <- array(NA_real_, dim = c(N, T, K + 1))
#' #   for (k in 1:(K+1)) {
#' #     Z_proj_final_array[,,k] <- matrix(Z_proj_final[,k], nrow = N, ncol = T, byrow = TRUE)
#' #   }
#' #
#' #   first_val <- covar_group[1]
#' #   count_first <- sum(covar_group == first_val)
#' #
#' #   if (count_first == 1) {
#' #     # Find the most frequent value in covar_group
#' #     most_freq_val <- as.numeric(names(sort(table(covar_group), decreasing = TRUE)[1]))
#' #     covar_group[1] <- most_freq_val
#' #   }
#' #   #covar_group[1:length(covar_group)]=1
#' #   #Z_proj = group_demean_formula_cpp(Z, unit_group, time_group, covar_group)
#' #   y_trans <- Z_proj[,,1] # Z_proj_final_array[,,1]
#' #   X_trans <- Z_proj[,,2:(K+1)] # Z_proj_final_array[,,2:(K+1)]
#' #
#' #   tY_vector <- as.vector(y_trans)  # (N*T)
#' #   tX_matrix <- matrix(aperm(X_trans, c(1, 2, 3)), nrow = N*T, ncol = K)
#' #
#' #   d_vector <- tX_matrix[, 1]
#' #   x_matrix <- tX_matrix[, -1, drop = FALSE]
#' #   # Separate transformed y and X
#' #
#' #   trans <- data.frame(y = tY_vector,
#' #                       D = d_vector,
#' #                       x_matrix)
#' #
#' #   trans = as.data.frame(trans)
#' #
#' #   dml_data <- DoubleMLData$new(
#' #     data = trans,
#' #     y_col = "y",
#' #     d_cols = "D",
#' #     x_cols = c(colnames(trans)[-c(1,2)])
#' #
#' #   )
#' #
#' #   # Define LASSO machine learning learners for nuisance parameter estimation
#' #   ml_l <- lrn("regr.cv_glmnet", s = "lambda.min")  # Outcome regression model
#' #   ml_m <- lrn("regr.cv_glmnet", s = "lambda.min")  # Treatment model (if applicable)
#' #
#' #   # Fit the Double Machine Learning model for treatment effect estimation
#' #   dml_plr <- DoubleMLPLR$new(dml_data, ml_l = ml_l, ml_m = ml_m)
#' #
#' #   # Fit the model to estimate the causal effect
#' #   dml_plr$fit(store_predictions=TRUE)
#' #   g_hat <- dml_plr$predictions$ml_l
#' #   m_hat <- dml_plr$predictions$ml_m
#' #
#' #   # Step 2: Compute residuals
#' #   Ytilde <- tY_vector - g_hat
#' #   Dtilde <- d_vector - m_hat
#' #   # Double ML
#' #   data_res = data.frame(id = data[[index[1]]], time = data[[index[2]]], Ytilde, Dtilde)
#' #   Post_plm = plm(Ytilde ~ -1 + Dtilde, data = data_res, model = "pooling", index=c("id", "time"))
#' #
#' #   coefs <- coef(Post_plm)
#' #   se_corrected <- sqrt(vcovHC(Post_plm, type = "HC0", method = "arellano")) * sqrt(N * T / ((N-G_unit)*(T-G_time)))
#' #   t_values_corrected <- coefs / se_corrected
#' #
#' #   # Calculate p-values from t-distribution for each coefficient
#' #   df <- Post_plm$df.residual  # degrees of freedom
#' #   p_values_corrected <- 2 * pt(-abs(t_values_corrected), df)
#' #
#' #   summary_table_correct <- data.frame(
#' #     Estimate = coefs,
#' #     SE_Corrected = se_corrected,
#' #     t_value_Corrected = t_values_corrected,
#' #     p_value_Corrected = p_values_corrected
#' #   )
#' #   colnames(summary_table_correct) = c('Estimate', 'Std. Error corrected', 't-value corrected', 'Pr(>|t|) corrected')
#' #   summary_table = summary(Post_plm)
#' #   summary_table$coefficients = cbind(summary_table_correct, summary_table$coefficients)
#' #
#' #   # Return cluster counts as well
#' #   return(list(
#' #     fit_summary = dml_plr$summary(),
#' #     G_unit = G_unit,
#' #     G_time = G_time,
#' #     G_covar = G_covar,
#' #     unit_group = unit_group,
#' #     time_group = time_group,
#' #     covar_group = covar_group,
#' #     post_plm_summary = summary(Post_plm),
#' #     estimate_corrected = summary_table_correct,
#' #     summary_table = summary_table
#' #   ))
#' # }
#'
#' #' @export
#' compute_u_hat <- function(z_array, unit_clusters, time_clusters, covar_clusters) {
#'   N <- dim(z_array)[1]
#'   T <- dim(z_array)[2]
#'   K <- dim(z_array)[3]
#'
#'   u_hat <- array(0, dim = c(N, T, K))
#'
#'   for (i in 1:N) {
#'     for (t in 1:T) {
#'       for (k in 1:K) {
#'         g_i <- unit_clusters[i]
#'         m_t <- time_clusters[t]
#'         l_k <- covar_clusters[k]
#'
#'         # z_{itk}
#'         zitk <- z_array[i, t, k]
#'
#'         # bar_z_{g_i t k}
#'         group_i_indices <- which(unit_clusters == g_i)
#'         bar_g_i_t_k <- mean(z_array[group_i_indices, t, k])
#'
#'         # bar_z_{i m_t k}
#'         time_m_indices <- which(time_clusters == m_t)
#'         bar_i_m_t_k <- mean(z_array[i, time_m_indices, k])
#'
#'         # bar_z_{i t l_k}
#'         covar_l_indices <- which(covar_clusters == l_k)
#'         bar_i_t_l_k <- mean(z_array[i, t, covar_l_indices])
#'
#'         # bar_z_{g_i m_t k}
#'         bar_g_i_m_t_k <- mean(z_array[group_i_indices, time_m_indices, k])
#'
#'         # bar_z_{g_i t l_k}
#'         bar_g_i_t_l_k <- mean(z_array[group_i_indices, t, covar_l_indices])
#'
#'         # bar_z_{i m_t l_k}
#'         bar_i_m_t_l_k <- mean(z_array[i, time_m_indices, covar_l_indices])
#'
#'         # Final u_hat_{itk}
#'         u_hat[i, t, k] <- 3 * zitk -
#'           2 * bar_g_i_t_k -
#'           2 * bar_i_m_t_k -
#'           2 * bar_i_t_l_k +
#'           bar_g_i_m_t_k +
#'           bar_g_i_t_l_k +
#'           bar_i_m_t_l_k
#'       }
#'     }
#'   }
#'
#'   return(u_hat)
#' }
#'
#' #' @export
#' cluster_Hierarchical <- function(y, X, link = "average", threshold = NULL, cluster = NULL,
#'                                  type = c("unit", "time", "covariate"),
#'                                  method_auto = c("none", "silhouette", "gap", "dynamicTreeCut"),
#'                                  data_for_gap = NULL, max_k = 10, deepSplit = TRUE, minClusterSize = 1, pamStage = FALSE) {
#'   type <- match.arg(type)
#'   method_auto <- match.arg(method_auto)
#'
#'   N <- nrow(y)
#'   T <- ncol(y)
#'   K <- dim(X)[3]
#'
#'   # Combine y and X into N x T x (K+1)
#'   combined_array <- array(0, dim = c(N, T, K + 1))
#'   combined_array[,,1] <- y
#'   combined_array[,,2:(K+1)] <- X
#'
#'   # Compute distance matrix based on type
#'   if (type == "unit") {
#'     dist_mat <- pseudo_dist_unit(combined_array)
#'   } else if (type == "time") {
#'     dist_mat <- pseudo_dist_time(combined_array)
#'   } else if (type == "covariate") {
#'     dist_mat <- pseudo_dist_covariate(combined_array)
#'   } else {
#'     stop("Invalid 'type' argument.")
#'   }
#'
#'   dist_obj <- as.dist(dist_mat)
#'   hc <- hclust(dist_obj, method = link)
#'
#'   if (!is.null(cluster)) {
#'     G <- cluster
#'     clusters <- cutree(hc, k = G)
#'
#'   } else if (!is.null(threshold)) {
#'     clusters <- cutree(hc, h = threshold)
#'     G <- length(unique(clusters))
#'
#'   } else if (method_auto != "none") {
#'     max_k <- min(max_k, ifelse(type == "unit", N, ifelse(type == "time", T, K + 1)) - 1)
#'
#'     if (method_auto == "gap") {
#'       if (is.null(data_for_gap)) {
#'         stop("For method_auto = 'gap', please provide 'data_for_gap' matrix.")
#'       }
#'       gap_fun <- function(x, k) {
#'         dist_x <- dist(x)
#'         hc_x <- hclust(dist_x, method = link)
#'         clust <- cutree(hc_x, k = k)
#'         list(cluster = clust)  # must return list with $cluster
#'       }
#'       gap_stat <- cluster::clusGap(data_for_gap, FUN = gap_fun, K.max = max_k, B = 50)
#'       G <- maxSE(gap_stat$Tab[, "gap"], gap_stat$Tab[, "SE.sim"], method = "firstSEmax")
#'       clusters <- cutree(hc, k = G)
#'
#'     } else if (method_auto == "silhouette") {
#'       sil_scores <- numeric(max_k)
#'       sil_scores[1] <- NA # silhouette not defined for k=1
#'       for (k in 2:max_k) {
#'         clust_try <- cutree(hc, k = k)
#'         ss <- silhouette(clust_try, dist_obj)
#'         sil_scores[k] <- mean(ss[, 3])
#'       }
#'       G <- which.max(sil_scores)
#'       clusters <- cutree(hc, k = G)
#'
#'     } else if (method_auto == "dynamicTreeCut") {
#'       clusters <- cutreeDynamic(dendro = hc, distM = as.matrix(dist_obj),
#'                                 deepSplit = deepSplit, minClusterSize = minClusterSize, pamStage = pamStage)
#'       G <- length(unique(clusters[clusters > 0]))  # exclude noise (0)
#'       noise_idx <- which(clusters == 0)
#'       if (length(noise_idx) > 0) {
#'         clusters[noise_idx] <- G + 1
#'         G <- G + 1
#'       }
#'     }
#'
#'   } else {
#'     stop("Must provide either 'cluster', 'threshold', or set 'method_auto' != 'none'.")
#'   }
#'
#'   # Format clusters and indicator matrix (like kmeans())
#'   clusters <- as.integer(factor(clusters))
#'   G <- length(unique(clusters))
#'   indicator <- sapply(1:G, function(g) as.integer(clusters == g))
#'   colnames(indicator) <- paste0("Cluster", 1:G)
#'
#'   dim_cluster <- switch(type,
#'                         unit = N,
#'                         time = T,
#'                         covariate = K + 1)
#'
#'   rownames(indicator) <- paste0(type, "_", 1:dim_cluster)
#'
#'   list(G = G, clusters = clusters, indicator = indicator)
#' }
#'
#' #' @export
#' cluster_general <- function(Y, X_list, N, T, init = 100, type = "long", groups = NULL, cc = 0, gamma = 1, dim_moment = 1) {
#'   dimtheta <- length(X_list)
#'   mdim <- dim_moment
#'   X_list <- c(X_list, list(Y))
#'   K = dimtheta
#'   if (type == "long") {
#'     mom_i <- c()
#'     ## --- Unit-side moments ---
#'     X_av <- matrix(0, nrow = N, ncol = dim_moment)  # N x dim_moment
#'
#'     for (p in 1:dim_moment) {
#'       # For each power
#'       X_pow <- matrix(0, nrow = N, ncol = T)  # N x T
#'
#'       for (t in 1:T) {
#'         # Collect all covariates at time t
#'         X_t <- sapply(X_list, function(X) X[, t])  # N x K
#'         # Take power p **before** averaging over K
#'         if ( p == 1){
#'         X_pow[, t] <- rowMeans(X_t)             # (1/K) sum_k x_itk^p
#'         }else if (p == 2){
#'           X_pow[, t] <- rowMeans(tanh(X_t))
#'         }else if (p == 3){
#'           X_pow[, t] <- rowMeans((X_t/(1+abs(X_t))))
#'         }else{
#'           X_pow[, t] <- rowMeans(X_t^p)
#'       }
#'       }
#'
#'       # Average over time
#'       X_av[, p] <- rowMeans(X_pow)               # (1/T) sum_t (1/K) sum_k x_itk^p
#'     }
#'
#'     mom_i <- cbind(mom_i, X_av)
#'
#'     # mom_micro <- c()
#'     mom_micro <- c()
#'     X_cbind <- do.call(cbind, X_list)
#'
#'     for (p in 1:dim_moment) {
#'       if ( p == 1){
#'         mom_micro <- cbind(mom_micro, X_cbind^p)             # (1/K) sum_k x_itk^p
#'       }else if (p == 2){
#'         mom_micro <- cbind(mom_micro, tanh(X_cbind))
#'       }else if (p == 3){
#'         mom_micro <- cbind(mom_micro, (X_cbind/(1+abs(X_cbind))))
#'       }else{
#'         X_pow[, t] <- rowMeans(X_t^p)
#'       }
#'       mom_micro <- cbind(mom_micro, X_cbind^p)
#'     }
#'
#'     ## --- Demean & rescale ---
#'     for (j in 1:mdim) {
#'       col_mean <- mean(mom_i[, j])
#'       col_sd <- sd(mom_i[, j])
#'       if (col_sd > 0) {
#'         cols <- ((j - 1) * T * (K+1) + 1):(j * T * (K+1))
#'         mom_micro[, cols] <- (mom_micro[, cols] - col_mean) / col_sd
#'         mom_i[, j] <- (mom_i[, j] - col_mean) / col_sd
#'       }
#'     }
#'
#'
#'     # --- Variance / noise on rescaled data ---
#'     variance <-  sum(sapply(1:mdim, function(j) {
#'       norm(mom_micro[, ((j - 1) * T * (K+1) + 1):(j * T * (K+1))] - mom_i[,j], type = "F")^2
#'     })) / (N * T^2 * (K+1)^2)
#'
#'     data <- mom_i
#'     dim_size <- N
#'
#'   } else if (type == "tall") {
#'     ## --- Time-side moments ---
#'
#'     mom_t <- c()
#'     X_av <- matrix(0, nrow = T, ncol = dim_moment)  # T x dim_moment
#'
#'     for (p in 1:dim_moment) {
#'       # For each power
#'       X_pow <- matrix(0, nrow = N, ncol = T)  # N x T
#'
#'       for (i in 1:N) {
#'         # Collect all covariates at time t
#'         X_i <- sapply(X_list, function(X) X[i, ])  # T x K
#'         # Take power p **before** averaging over K
#'         X_pow[i, ] <- rowMeans(X_i^p)             # (1/K) sum_k x_itk^p
#'       }
#'
#'       # Average over unit
#'       X_av[, p] <- colMeans(X_pow)               # (1/N) sum_i (1/K) sum_k x_itk^p
#'     }
#'
#'     mom_t <- cbind(mom_t, X_av)
#'
#'     mom_micro2 <- c()
#'     X_cbind <- do.call(cbind, lapply(X_list, t))
#'
#'     for (p in 1:dim_moment) {
#'       mom_micro2 <- cbind(mom_micro2, (X_cbind)^p)
#'     }
#'
#'     ## --- Demean & rescale ---
#'     for (j in 1:mdim) {
#'       col_mean <- mean(mom_t[, j])
#'       col_sd <- sd(mom_t[, j])
#'       if (col_sd > 0) {
#'         cols <- ((j - 1) * N * (K+1) + 1):(j * N * (K+1))
#'         mom_micro2[, cols] <- (mom_micro2[, cols] - col_mean) / col_sd
#'         mom_t[, j] <- (mom_t[, j] - col_mean) / col_sd
#'       }
#'     }
#'
#'
#'     ## --- Variance / noise on rescaled data ---
#'     variance <- sum(sapply(1:mdim, function(i) {
#'       norm(mom_micro2[, ((i - 1) * N * (K+1) + 1):(i * N * (K+1))] - mom_t[,i], type = "F")^2
#'     })) / (T * N^2 * (K+1)^2)
#'
#'     data <- mom_t
#'     dim_size <- T
#'
#'   } else if(type == "long_T_moment") {
#'     mom_i <- c()
#'     ## --- Unit-side moments ---
#'     X_av <- matrix(0, nrow = N, ncol = T)  # N x dim_moment
#'
#'
#'     for (t in 1:T) {
#'         # Collect all covariates at time t
#'         X_t <- sapply(X_list, function(X) X[, t])  # N x K
#'         # Take power p **before** averaging over K
#'         X_av[, t] <- rowMeans(X_t)             # (1/K) sum_k x_itk^p
#'       }
#'
#'
#'     mom_i <- cbind(mom_i, X_av)
#'
#'     # mom_micro <- c()
#'     mom_micro <- c()
#'     X_cbind <- do.call(cbind, X_list)
#'     mom_micro <- cbind(mom_micro, X_cbind)
#'
#'
#'     ## --- Demean & rescale ---
#'     for (j in 1:T) {
#'       col_mean <- mean(mom_i[, j])
#'       col_sd <- sd(mom_i[, j])
#'       if (col_sd > 0) {
#'         cols <- ((j - 1)  * (K+1) + 1):(j * (K+1))
#'         mom_micro[, cols] <- (mom_micro[, cols] - col_mean) / col_sd
#'         mom_i[, j] <- (mom_i[, j] - col_mean) / col_sd
#'       }
#'     }
#'
#'
#'     # --- Variance / noise on rescaled data ---
#'     variance <-  sum(sapply(1:T, function(j) {
#'       norm(mom_micro[, ((j - 1) * (K+1) + 1):(j  * (K+1))] - mom_i[,j], type = "F")^2
#'     })) / (N * (K+1)^2)
#'
#'     data <- mom_i
#'     dim_size <- N
#'
#'   }else {
#'     stop("Invalid type. Use 'long' or 'tall'.")
#'   }
#'
#'   ## --- Clustering ---
#'   if (!is.null(groups)) {
#'     clusters <- groups
#'     k_result <- kmeans(data, centers = clusters, algorithm = "Lloyd", nstart = init, iter.max = 100)
#'   } else {
#'     if (max(apply(data, 2, sd)) > 0) {
#'       clusters <- 1
#'       xx = 1000
#'       while (xx >= gamma * variance ) {
#'         k_result <- kmeans(data, centers = clusters, algorithm = "Lloyd", nstart = init, iter.max = 100)
#'         xx <- k_result$tot.withinss / dim_size
#'         clusters <- clusters + 1
#'       }
#'       clusters = min(clusters, dim_size-1)
#'       k_result <- kmeans(data,centers = clusters,algorithm = "Lloyd",nstart = init,iter.max = 100)
#'
#'     }else{
#'       k_result <- kmeans(data, centers = 1, algorithm = "Lloyd", nstart = init, iter.max = 100)
#'     }
#'   }
#'
#'   ## --- Return ---
#'   if (type == "long") {
#'     list(res = k_result, clusters = clusters, data = data)
#'   } else {
#'     list(res = k_result, clusters = clusters, data = data)
#'   }
#' }
#'
#' fix_time_clusters <- function(cluster, data, m = 2) {
#'   # cluster: raw K-means labels (length T)
#'   # data: T x 1 matrix (average per time period)
#'   # m: minimum cluster size
#'   #
#'   # Returns:
#'   #   cluster: cleaned labels for each time period
#'   #   G_time: number of clusters
#'
#'   cluster <- as.integer(cluster)
#'   T <- nrow(data)
#'
#'   # Collapse to 1 cluster if impossible
#'   if (T < 2 * m) {
#'     return(list(
#'       cluster = rep(1, T),
#'       G_time = 1
#'     ))
#'   }
#'
#'   repeat {
#'     sizes <- table(cluster)
#'
#'     small <- as.integer(names(sizes[sizes < m]))
#'     large <- as.integer(names(sizes[sizes >= m]))
#'
#'     if (length(small) == 0) break
#'
#'     # Compute centers (safe for 1-column data)
#'     unique_cl <- sort(unique(cluster))
#'     centers <- do.call(rbind, lapply(unique_cl, function(cl) {
#'       colMeans(data[cluster == cl, , drop = FALSE])
#'     }))
#'     rownames(centers) <- unique_cl
#'
#'     # Reassign points from small clusters
#'     for (sc in small) {
#'       idx <- which(cluster == sc)
#'       for (i in idx) {
#'         target <- if (length(large) > 0) large else setdiff(unique_cl, sc)
#'         dists <- sapply(target, function(cl) {
#'           sum((data[i, ] - centers[as.character(cl), ])^2)
#'         })
#'         cluster[i] <- target[which.min(dists)]
#'       }
#'     }
#'   }
#'
#'   # Relabel clusters 1:K
#'   unique_cl <- sort(unique(cluster))
#'   map <- setNames(seq_along(unique_cl), unique_cl)
#'   cluster <- map[as.character(cluster)]
#'
#'   # Number of clusters
#'   G_time <- length(unique_cl)
#'
#'   return(list(
#'     cluster = as.integer(cluster),
#'     G_time = G_time
#'   ))
#' }


#' #' Estimator Function with Optional Cross-Fitting
#' #'@title Inference in high-dimensional data after discretizing unobserved heterogeneity
#' #'
#' #' @description R package 'HDpcluster' is dedicated to do inference for high-dimensional linear panel data model with unkown functions of fixed effects.
#' #'
#' #' @param y outcome variable
#' #' @param D treatment variable
#' #' @param X control variables
#' #' @param T panel size of the time period
#' #' @param groups_init number of the unit clusters
#' #' @param index index name
#' #' @param data data which contains the correct index of the outcome variable or the treatment variable
#' #' @param cluster_type different cluster means for the data, 'unit kmeans' is default, allow for 'unit pesudo'
#' #' @param pesudo_type if cluster_type = 'unit pesudo', choice of the pesudo type
#' #' @param link if cluster_type = 'unit pesudo', choice of different links of pesudo type
#' #' @param optimal_index if cluster_type = 'unit pesudo', different ways to compute the optimal number of clusters
#' #'
#' #' @returns A list of fitted results is returned.
#' #' Within this outputted list, the following elements can be found:
#' #'     \item{res}{regression model.}
#' #'     \item{G}{number of unit clusters.}
#' #'     \item{estimate_correct}{corrected standard error, t value, and p value.}
#' #'     \item{summary_table}{summary of the model together with corrected estimates in the 'Coefficients'.}
#' #'
#' #' @import hdm
#' #' @import plm
#' #' @useDynLib HDpcluster
#' #' @export
#' HP_estimate <- function(data, y_col = NULL, covariate_cols = NULL,
#'                         id_col = "id", time_col = "time", cluster_type = c('two way'),
#'                         index = c("id", "time"), method_auto = "dynamicTreeCut", cluster_method = c('kmeans', 'hierarchical'), pre_cluster = FALSE, Du_pre, Dv_pre, unit_cluster = NULL, time_cluster = NULL, gamma = 1, cc = 0, dim_moment = 1) {
#'
#'
#'   if (is.null(covariate_cols)) {
#'     covariate_cols <- setdiff(colnames(data), c(id_col, time_col, y_col))
#'   }
#'
#'   ids <- sort(unique(data[[id_col]]))
#'   times <- sort(unique(data[[time_col]]))
#'   N <- length(ids)
#'   T <- length(times)
#'   K <- length(covariate_cols)
#'
#'   # Initialize y and X
#'   y <- matrix(NA_real_, nrow = N, ncol = T)
#'   X <- array(NA_real_, dim = c(N, T, K))
#'
#'   for (i in seq_len(nrow(data))) {
#'     id_idx <- which(ids == data[[id_col]][i])
#'     time_idx <- which(times == data[[time_col]][i])
#'     y[id_idx, time_idx] <- data[[y_col]][i]
#'     for (k in seq_along(covariate_cols)) {
#'       X[id_idx, time_idx, k] <- data[[covariate_cols[k]]][i]
#'     }
#'   }
#'
#'   # Cluster (output indicators assumed one-hot encoded)
#'   X_norm <- array(NA_real_, dim = dim(X))  # Same shape: (N, T, K)
#'   for (k in 1:dim(X)[3]) {
#'     X_slice <- X[,,k]
#'     X_norm[,,k] <- scale(X_slice)
#'   }
#'   y_norm = scale(y)
#'
#'
#'   if (cluster_type == 'two way'){
#'     if (cluster_method == 'kmeans'){
#'
#'       # kmeans
#'       X_list = lapply(seq_len(dim(X)[3]), function(k) X[ , , k])
#'       if (pre_cluster == FALSE){
#'         # cluster
#'         if (is.null(unit_cluster) == 1 & is.null(time_cluster) == 1){
#'           clusteri <- cluster_general(y, X_list, N, T, init = 50, type = "long", gamma = gamma, cc = cc, dim_moment = dim_moment  )
#'           G_unit <- clusteri$clusters
#'           klong <- clusteri$res
#'
#'
#'           clustert <- cluster_general(y, X_list, N, T, init = 50, type = "tall", gamma = gamma, cc = cc, dim_moment = dim_moment  )
#'           G_time <- clustert$clusters
#'           ktall <- clustert$res
#'           newcluster_T = fix_time_clusters(cluster = ktall$cluster, data = clustert$data, m = 3)
#'           ktall$cluster =   newcluster_T$cluster
#'           G_time =  newcluster_T$G_time
#'           #   if (G_time > floor(T/2)){
#'           #     clustert <- cluster_general(y, X_list, N, T, init = 50, type = "tall", groups = 1 )
#'           #     G_time <- clustert$clusters
#'           #     ktall <- clustert$res
#'           #   }
#'         }else{
#'           clusteri <- cluster_general(y, X_list, N, T, init = 50, type = "long" , groups = c(floor(unit_cluster)) )
#'           G_unit <- clusteri$clusters
#'           klong <- clusteri$res
#'
#'
#'           clustert <- cluster_general(y, X_list, N, T, init = 50, type = "tall" , groups = c(floor(time_cluster)))
#'           G_time <- clustert$clusters
#'           ktall <- clustert$res
#'         }
#'
#'
#'         Du <- matrix(0, N, G_unit)
#'         Dv <- matrix(0, T, G_time)
#'         Diagu<- matrix(0, G_unit, G_unit)
#'         Diagv<- matrix(0, G_time, G_time)
#'
#'         for (j in seq_len(G_unit)) {
#'           Du[, j] <- as.numeric(klong$cluster == j)
#'           Diagu[j,j] <- 1/sum(Du[,j])
#'         }
#'
#'         for (j in seq_len(G_time)) {
#'           Dv[, j] <- as.numeric(ktall$cluster == j)
#'           Diagv[j,j] <- 1/sum(Dv[,j])
#'         }
#'       }else if (pre_cluster == TRUE){
#'         Du = Du_pre
#'         Dv = Dv_pre
#'         Diagu<- matrix(0, G_unit, G_unit)
#'         Diagv<- matrix(0, G_time, G_time)
#'
#'         for (j in seq_len(G_unit)) {
#'           Diagu[j,j] <- 1/sum(Du[,j])
#'         }
#'
#'         for (j in seq_len(G_time)) {
#'           Diagv[j,j] <- 1/sum(Dv[,j])
#'         }
#'         G_unit = dim(Du)[2]
#'         G_time = dim(Dv)[2]
#'       }
#'     }else if (cluster_method == 'hierarchical'){
#'
#'       # heriachical
#'       if (is.null(unit_cluster) == 1 & is.null(time_cluster) == 1){
#'         res_unit <- cluster_Hierarchical(y_norm, X_norm, type = "unit", method_auto = method_auto)
#'         res_time <- cluster_Hierarchical(y_norm, X_norm, type = "time", method_auto = method_auto)
#'       }else{
#'         res_unit <- cluster_Hierarchical(y_norm, X_norm, type = "unit", cluster = unit_cluster)
#'         res_time <- cluster_Hierarchical(y_norm, X_norm, type = "time", cluster = time_cluster)
#'       }
#'
#'       Du <- res_unit$indicator    # N x G_unit
#'       Dv <- res_time$indicator    # T x G_time
#'       G_unit = res_unit$G
#'       G_time = res_time$G
#'       Diagu<- matrix(0, G_unit, G_unit)
#'       Diagv<- matrix(0, G_time, G_time)
#'
#'       for (j in seq_len(G_unit)) {
#'         Diagu[j,j] <- 1/sum(Du[,j])
#'       }
#'
#'       for (j in seq_len(G_time)) {
#'         Diagv[j,j] <- 1/sum(Dv[,j])
#'       }
#'     }
#'
#'     unit_group <- max.col(Du)
#'     time_group <- max.col(Dv)
#'
#'     # Projection matrices
#'     Mu <- diag(N) - Du %*% Diagu %*% t(Du)
#'     Mv <- diag(T) - Dv %*% Diagv %*% t(Dv)
#'     # Mu <- diag(N) - Du %*% solve(t(Du) %*% Du) %*% t(Du)
#'     # Mv <- diag(T) - Dv %*% solve(t(Dv) %*% Dv) %*% t(Dv)
#'
#'     # Stack y and X into Z: N x T x (K+1)
#'     Z <- array(NA_real_, dim = c(N, T, K + 1))
#'     Z[,,1] <- y
#'     for (k in 1:K) Z[,,k+1] <- X[,,k]
#'
#'     # Projection method: demean unit and time clusters slice-wise
#'     # When T=1, Mv becomes zero matrix, so skip time demeaning
#'     Z_proj <- array(NA_real_, dim = c(N, T, K + 1))
#'     for (k in 1:(K+1)) {
#'       if (T == 1) {
#'         Z_proj[,,k] <- Mu %*% Z[,,k]
#'       } else {
#'         Z_proj[,,k] <- Mu %*% Z[,,k] %*% Mv
#'       }
#'     }
#'
#'     # Separate transformed y and X
#'     y_trans <- Z_proj[,,1] # Z_proj_final_array[,,1]
#'     X_trans <- Z_proj[,,2:(K+1)] # Z_proj_final_array[,,2:(K+1)]
#'
#'     tY_vector <- as.vector(y_trans)  # (N*T)
#'     # Handle T=1 case where X_trans is 2D instead of 3D
#'     if (length(dim(X_trans)) == 2) {
#'       tX_matrix <- matrix(as.vector(t(X_trans)), nrow = N*T, ncol = K)
#'     } else {
#'       tX_matrix <- matrix(aperm(X_trans, c(1, 2, 3)), nrow = N*T, ncol = K)
#'     }
#'
#'     d_vector <- tX_matrix[, 1]
#'     x_matrix <- tX_matrix[, -1, drop = FALSE]
#'
#'     # Run rlassoEffect with correct inputs
#'     fit <- rlassoEffect(x = x_matrix, y = tY_vector, d = d_vector, method = "double selection")
#'
#'     trans <- data.frame(y = tY_vector,
#'                         D = d_vector,
#'                         x_matrix)
#'
#'     lasso.Y <- rlasso(y ~ . - D - 1, data = trans)
#'     Ytilde <- lasso.Y$residuals
#'
#'     lasso.D <- rlasso(D ~ . - 1, data = trans[, -1])
#'     Dtilde <- lasso.D$residuals
#'
#'     # Create proper panel structure for residuals
#'     # as.vector() on N x T matrix goes column-wise: unit 1 time 1, unit 2 time 1, ..., unit N time 1, unit 1 time 2, ...
#'     data_res <- data.frame(id = rep(ids, T),
#'                            time = rep(times, each = N),
#'                            Ytilde = Ytilde,
#'                            Dtilde = Dtilde)
#'
#'     # For T=1, use lm() since there's no panel structure; for T>1 use plm with pooling
#'     if (T == 1) {
#'       Post_plm <- lm(Ytilde ~ -1 + Dtilde, data = data_res)
#'     } else {
#'       Post_plm <- plm(Ytilde ~ -1 + Dtilde, data = data_res, model = "pooling", index = c("id", "time"))
#'     }
#'
#'     coefs <- coef(Post_plm)
#'     if (N*T - N*G_time - T*G_unit > 0){
#'       if (T == 1) {
#'         se_corrected <- sqrt(diag(vcovHC(Post_plm, type = "HC0"))) * sqrt(N * T / (N*T - N*G_time - T*G_unit))
#'       } else {
#'         se_corrected <- sqrt(diag(vcovHC(Post_plm, type = "HC0", method = "arellano"))) * sqrt(N * T / (N*T - N*G_time - T*G_unit))
#'       }
#'     }else{
#'       if (T == 1) {
#'         se_corrected <- sqrt(diag(vcovHC(Post_plm, type = "HC0"))) * sqrt(N * T / (1))
#'       } else {
#'         se_corrected <- sqrt(diag(vcovHC(Post_plm, type = "HC0", method = "arellano"))) * sqrt(N * T / (1))
#'       }
#'     }
#'     t_values_corrected <- coefs / se_corrected
#'
#'     df <- Post_plm$df.residual
#'     p_values_corrected <- 2 * pt(-abs(t_values_corrected), df)
#'
#'     summary_table_correct <- data.frame(
#'       Estimate = coefs,
#'       `Std. Error corrected` = se_corrected,
#'       `t-value corrected` = t_values_corrected,
#'       `Pr(>|t|) corrected` = p_values_corrected
#'     )
#'
#'     colnames(summary_table_correct) = c('Estimate', 'Std. Error corrected', 't-value corrected', 'Pr(>|t|) corrected')
#'
#'     # Return cluster counts as well
#'     return(list(
#'       fit_summary = summary(fit),
#'       G_unit = G_unit,
#'       G_time = G_time,
#'       unit_group = unit_group,
#'       time_group = time_group,
#'       post_plm_summary = summary(Post_plm),
#'       estimate_corrected = summary_table_correct,
#'       summary_table = summary_table_correct
#'     ))
#'   }else if (cluster_type == 'one way'){
#'     if (cluster_method == 'kmeans'){
#'
#'       # kmeans
#'       X_list = lapply(seq_len(dim(X)[3]), function(k) X[ , , k])
#'       if (pre_cluster == FALSE){
#'         # cluster
#'         if (is.null(unit_cluster) == 1 & is.null(time_cluster) == 1){
#'           clusteri <- cluster_general(y, X_list, N, T, init = 50, type = "long", gamma = gamma, cc = cc, dim_moment = dim_moment  )
#'           G_unit <- clusteri$clusters
#'           klong <- clusteri$res
#'
#'
#'         }else{
#'           clusteri <- cluster_general(y, X_list, N, T, init = 50, type = "long" , groups = c(floor(unit_cluster)) )
#'           G_unit <- clusteri$clusters
#'           klong <- clusteri$res
#'         }
#'         Du <- matrix(0, N, G_unit)
#'         Diagu<- matrix(0, G_unit, G_unit)
#'
#'
#'         for (j in seq_len(G_unit)) {
#'           Du[, j] <- as.numeric(klong$cluster == j)
#'           Diagu[j,j] <- 1/sum(Du[,j])
#'         }
#'
#'
#'       }else if (pre_cluster == TRUE){
#'         Du = Du_pre
#'         Diagu<- matrix(0, G_unit, G_unit)
#'
#'         for (j in seq_len(G_unit)) {
#'           Diagu[j,j] <- 1/sum(Du[,j])
#'         }
#'
#'         G_unit = dim(Du)[2]
#'       }
#'     }else if (cluster_method == 'hierarchical'){
#'
#'       # heriachical
#'       if (is.null(unit_cluster) == 1 & is.null(time_cluster) == 1){
#'         res_unit <- cluster_Hierarchical(y_norm, X_norm, type = "unit", method_auto = method_auto)
#'       }else{
#'         res_unit <- cluster_Hierarchical(y_norm, X_norm, type = "unit", cluster = unit_cluster)
#'       }
#'
#'       Du <- res_unit$indicator    # N x G_unit
#'       G_unit = res_unit$G
#'       Diagu<- matrix(0, G_unit, G_unit)
#'
#'       for (j in seq_len(G_unit)) {
#'         Diagu[j,j] <- 1/sum(Du[,j])
#'       }
#'
#'     }
#'
#'     unit_group <- max.col(Du)
#'
#'     # Projection matrices
#'     Mu <- diag(N) - Du %*% Diagu %*% t(Du)
#'
#'     # Stack y and X into Z: N x T x (K+1)
#'     Z <- array(NA_real_, dim = c(N, T, K + 1))
#'     Z[,,1] <- y
#'     for (k in 1:K) Z[,,k+1] <- X[,,k]
#'
#'     # Projection method: demean unit and time clusters slice-wise
#'     Z_proj <- array(NA_real_, dim = c(N, T, K + 1))
#'     for (k in 1:(K+1)) {
#'       Z_proj[,,k] <- Mu %*% Z[,,k]
#'     }
#'
#'     # Separate transformed y and X
#'     y_trans <- Z_proj[,,1] # Z_proj_final_array[,,1]
#'     X_trans <- Z_proj[,,2:(K+1)] # Z_proj_final_array[,,2:(K+1)]
#'
#'     tY_vector <- as.vector(y_trans)  # (N*T)
#'     # Handle T=1 case where X_trans is 2D instead of 3D
#'     if (length(dim(X_trans)) == 2) {
#'       tX_matrix <- matrix(as.vector(t(X_trans)), nrow = N*T, ncol = K)
#'     } else {
#'       tX_matrix <- matrix(aperm(X_trans, c(1, 2, 3)), nrow = N*T, ncol = K)
#'     }
#'     # debiased lasso
#'     # debiased_lasso = lasso.proj(tX_matrix, tY_vector) # delas(tX_matrix, tY_vector, c(1))
#'     # debiased_lasso$se = debiased_lasso$se * sqrt(N * T / (N*T - T*G_unit))
#'
#'     #
#'     d_vector <- tX_matrix[, 1]
#'     x_matrix <- tX_matrix[, -1, drop = FALSE]
#'
#'     # Run rlassoEffect with correct inputs
#'     fit <- rlassoEffect(x = x_matrix, y = tY_vector, d = d_vector, method = "double selection")
#'
#'     trans <- data.frame(y = tY_vector,
#'                         D = d_vector,
#'                         x_matrix)
#'
#'     lasso.Y <- rlasso(y ~ . - D - 1, data = trans)
#'     Ytilde <- lasso.Y$residuals
#'
#'     lasso.D <- rlasso(D ~ . - 1, data = trans[, -1])
#'     Dtilde <- lasso.D$residuals
#'
#'     data_res <- data.frame(id = data[[index[1]]],
#'                            time = data[[index[2]]],
#'                            Ytilde = Ytilde,
#'                            Dtilde = Dtilde)
#'
#'     Post_plm <- plm(Ytilde ~ -1 + Dtilde, data = data_res, model = "pooling", index = c("id", "time"))
#'
#'     coefs <- coef(Post_plm)
#'     se_corrected <- sqrt(vcovHC(Post_plm, type = "HC0", method = "arellano")) * sqrt(N * T / (N*T - T*G_unit))
#'     t_values_corrected <- coefs / se_corrected
#'
#'     df <- Post_plm$df.residual
#'     p_values_corrected <- 2 * pt(-abs(t_values_corrected), df)
#'
#'     summary_table_correct <- data.frame(
#'       Estimate = coefs,
#'       `Std. Error corrected` = se_corrected,
#'       `t-value corrected` = t_values_corrected,
#'       `Pr(>|t|) corrected` = p_values_corrected
#'     )
#'
#'     colnames(summary_table_correct) = c('Estimate', 'Std. Error corrected', 't-value corrected', 'Pr(>|t|) corrected')
#'
#'     # Return cluster counts as well
#'     return(list(
#'       fit_summary = summary(fit),
#'       G_unit = G_unit,
#'       unit_group = unit_group,
#'       post_plm_summary = summary(Post_plm),
#'       estimate_corrected = summary_table_correct,
#'       summary_table = summary_table_correct
#'     ))
#'
#'   }else if (cluster_type == 'one way T moment'){
#'     if (cluster_method == 'kmeans'){
#'
#'       # kmeans
#'       X_list = lapply(seq_len(dim(X)[3]), function(k) X[ , , k])
#'       if (pre_cluster == FALSE){
#'         # cluster
#'         if (is.null(unit_cluster) == 1 & is.null(time_cluster) == 1){
#'           clusteri <- cluster_general(y, X_list, N, T, init = 50, type = "long_T_moment", gamma = gamma, cc = cc, dim_moment = dim_moment  )
#'           G_unit <- clusteri$clusters
#'           klong <- clusteri$res
#'
#'
#'         }else{
#'           clusteri <- cluster_general(y, X_list, N, T, init = 50, type = "long_T_moment" , groups = c(floor(unit_cluster)) )
#'           G_unit <- clusteri$clusters
#'           klong <- clusteri$res
#'         }
#'         Du <- matrix(0, N, G_unit)
#'         Diagu<- matrix(0, G_unit, G_unit)
#'
#'
#'         for (j in seq_len(G_unit)) {
#'           Du[, j] <- as.numeric(klong$cluster == j)
#'           Diagu[j,j] <- 1/sum(Du[,j])
#'         }
#'
#'
#'       }else if (pre_cluster == TRUE){
#'         Du = Du_pre
#'         Diagu<- matrix(0, G_unit, G_unit)
#'
#'         for (j in seq_len(G_unit)) {
#'           Diagu[j,j] <- 1/sum(Du[,j])
#'         }
#'
#'         G_unit = dim(Du)[2]
#'       }
#'     }else if (cluster_method == 'hierarchical'){
#'
#'       # heriachical
#'       if (is.null(unit_cluster) == 1 & is.null(time_cluster) == 1){
#'         res_unit <- cluster_Hierarchical(y_norm, X_norm, type = "unit", method_auto = method_auto)
#'       }else{
#'         res_unit <- cluster_Hierarchical(y_norm, X_norm, type = "unit", cluster = unit_cluster)
#'       }
#'
#'       Du <- res_unit$indicator    # N x G_unit
#'       G_unit = res_unit$G
#'       Diagu<- matrix(0, G_unit, G_unit)
#'
#'       for (j in seq_len(G_unit)) {
#'         Diagu[j,j] <- 1/sum(Du[,j])
#'       }
#'
#'     }
#'
#'     unit_group <- max.col(Du)
#'
#'     # Projection matrices
#'     Mu <- diag(N) - Du %*% Diagu %*% t(Du)
#'
#'     # Stack y and X into Z: N x T x (K+1)
#'     Z <- array(NA_real_, dim = c(N, T, K + 1))
#'     Z[,,1] <- y
#'     for (k in 1:K) Z[,,k+1] <- X[,,k]
#'
#'     # Projection method: demean unit and time clusters slice-wise
#'     Z_proj <- array(NA_real_, dim = c(N, T, K + 1))
#'     for (k in 1:(K+1)) {
#'       Z_proj[,,k] <- Mu %*% Z[,,k]
#'     }
#'
#'     # Separate transformed y and X
#'     y_trans <- Z_proj[,,1] # Z_proj_final_array[,,1]
#'     X_trans <- Z_proj[,,2:(K+1)] # Z_proj_final_array[,,2:(K+1)]
#'
#'     tY_vector <- as.vector(y_trans)  # (N*T)
#'     # Handle T=1 case where X_trans is 2D instead of 3D
#'     if (length(dim(X_trans)) == 2) {
#'       tX_matrix <- matrix(as.vector(t(X_trans)), nrow = N*T, ncol = K)
#'     } else {
#'       tX_matrix <- matrix(aperm(X_trans, c(1, 2, 3)), nrow = N*T, ncol = K)
#'     }
#'
#'     # debiased lasso
#'     # debiased_lasso = lasso.proj(tX_matrix, tY_vector) # delas(tX_matrix, tY_vector, c(1))
#'     # debiased_lasso$se = debiased_lasso$se * sqrt(N * T / (N*T - T*G_unit))
#'
#'     #
#'     d_vector <- tX_matrix[, 1]
#'     x_matrix <- tX_matrix[, -1, drop = FALSE]
#'
#'     # Run rlassoEffect with correct inputs
#'     fit <- rlassoEffect(x = x_matrix, y = tY_vector, d = d_vector, method = "double selection")
#'
#'     trans <- data.frame(y = tY_vector,
#'                         D = d_vector,
#'                         x_matrix)
#'
#'     lasso.Y <- rlasso(y ~ . - D - 1, data = trans)
#'     Ytilde <- lasso.Y$residuals
#'
#'     lasso.D <- rlasso(D ~ . - 1, data = trans[, -1])
#'     Dtilde <- lasso.D$residuals
#'
#'     data_res <- data.frame(id = data[[index[1]]],
#'                            time = data[[index[2]]],
#'                            Ytilde = Ytilde,
#'                            Dtilde = Dtilde)
#'
#'     Post_plm <- plm(Ytilde ~ -1 + Dtilde, data = data_res, model = "pooling", index = c("id", "time"))
#'
#'     coefs <- coef(Post_plm)
#'     se_corrected <- sqrt(vcovHC(Post_plm, type = "HC0", method = "arellano")) * sqrt(N * T / (N*T - T*G_unit))
#'     t_values_corrected <- coefs / se_corrected
#'
#'     df <- Post_plm$df.residual
#'     p_values_corrected <- 2 * pt(-abs(t_values_corrected), df)
#'
#'     summary_table_correct <- data.frame(
#'       Estimate = coefs,
#'       `Std. Error corrected` = se_corrected,
#'       `t-value corrected` = t_values_corrected,
#'       `Pr(>|t|) corrected` = p_values_corrected
#'     )
#'
#'     colnames(summary_table_correct) = c('Estimate', 'Std. Error corrected', 't-value corrected', 'Pr(>|t|) corrected')
#'
#'     # Return cluster counts as well
#'     return(list(
#'       fit_summary = summary(fit),
#'       G_unit = G_unit,
#'       unit_group = unit_group,
#'       post_plm_summary = summary(Post_plm),
#'       estimate_corrected = summary_table_correct,
#'       summary_table = summary_table_correct
#'     ))
#'
#'   }else if (cluster_type == 'one way more moments'){
#'     if (cluster_method == 'kmeans'){
#'
#'       # kmeans
#'       X_list = lapply(seq_len(dim(X)[3]), function(k) X[ , , k])
#'       if (pre_cluster == FALSE){
#'         # cluster
#'         if (is.null(unit_cluster) == 1 & is.null(time_cluster) == 1){
#'           clusteri <- cluster_general(y, X_list, N, T, init = 50, type = "long multiple moments", gamma = gamma, cc = cc, dim_moment = dim_moment  )
#'           G_unit <- clusteri$clusters
#'           klong <- clusteri$res
#'
#'
#'         }else{
#'           clusteri <- cluster_general(y, X_list, N, T, init = 50, type = "long multiple moments" , groups = c(floor(unit_cluster)) )
#'           G_unit <- clusteri$clusters
#'           klong <- clusteri$res
#'         }
#'         Du <- matrix(0, N, G_unit)
#'         Diagu<- matrix(0, G_unit, G_unit)
#'
#'
#'         for (j in seq_len(G_unit)) {
#'           Du[, j] <- as.numeric(klong$cluster == j)
#'           Diagu[j,j] <- 1/sum(Du[,j])
#'         }
#'
#'
#'       }else if (pre_cluster == TRUE){
#'         Du = Du_pre
#'         Diagu<- matrix(0, G_unit, G_unit)
#'
#'         for (j in seq_len(G_unit)) {
#'           Diagu[j,j] <- 1/sum(Du[,j])
#'         }
#'
#'         G_unit = dim(Du)[2]
#'       }
#'     }else if (cluster_method == 'hierarchical'){
#'
#'       # heriachical
#'       if (is.null(unit_cluster) == 1 & is.null(time_cluster) == 1){
#'         res_unit <- cluster_Hierarchical(y_norm, X_norm, type = "unit", method_auto = method_auto)
#'       }else{
#'         res_unit <- cluster_Hierarchical(y_norm, X_norm, type = "unit", cluster = unit_cluster)
#'       }
#'
#'       Du <- res_unit$indicator    # N x G_unit
#'       G_unit = res_unit$G
#'       Diagu<- matrix(0, G_unit, G_unit)
#'
#'       for (j in seq_len(G_unit)) {
#'         Diagu[j,j] <- 1/sum(Du[,j])
#'       }
#'
#'     }
#'
#'     unit_group <- max.col(Du)
#'
#'     # Projection matrices
#'     Mu <- diag(N) - Du %*% Diagu %*% t(Du)
#'
#'     # Stack y and X into Z: N x T x (K+1)
#'     Z <- array(NA_real_, dim = c(N, T, K + 1))
#'     Z[,,1] <- y
#'     for (k in 1:K) Z[,,k+1] <- X[,,k]
#'
#'     # Projection method: demean unit and time clusters slice-wise
#'     Z_proj <- array(NA_real_, dim = c(N, T, K + 1))
#'     for (k in 1:(K+1)) {
#'       Z_proj[,,k] <- Mu %*% Z[,,k]
#'     }
#'
#'     # Separate transformed y and X
#'     y_trans <- Z_proj[,,1] # Z_proj_final_array[,,1]
#'     X_trans <- Z_proj[,,2:(K+1)] # Z_proj_final_array[,,2:(K+1)]
#'
#'     tY_vector <- as.vector(y_trans)  # (N*T)
#'     # Handle T=1 case where X_trans is 2D instead of 3D
#'     if (length(dim(X_trans)) == 2) {
#'       tX_matrix <- matrix(as.vector(t(X_trans)), nrow = N*T, ncol = K)
#'     } else {
#'       tX_matrix <- matrix(aperm(X_trans, c(1, 2, 3)), nrow = N*T, ncol = K)
#'     }
#'
#'     # debiased lasso
#'     # debiased_lasso = lasso.proj(tX_matrix, tY_vector) # delas(tX_matrix, tY_vector, c(1))
#'     # debiased_lasso$se = debiased_lasso$se * sqrt(N * T / (N*T - T*G_unit))
#'
#'     #
#'     d_vector <- tX_matrix[, 1]
#'     x_matrix <- tX_matrix[, -1, drop = FALSE]
#'
#'     # Run rlassoEffect with correct inputs
#'     fit <- rlassoEffect(x = x_matrix, y = tY_vector, d = d_vector, method = "double selection")
#'
#'     trans <- data.frame(y = tY_vector,
#'                         D = d_vector,
#'                         x_matrix)
#'
#'     lasso.Y <- rlasso(y ~ . - D - 1, data = trans)
#'     Ytilde <- lasso.Y$residuals
#'
#'     lasso.D <- rlasso(D ~ . - 1, data = trans[, -1])
#'     Dtilde <- lasso.D$residuals
#'
#'     data_res <- data.frame(id = data[[index[1]]],
#'                            time = data[[index[2]]],
#'                            Ytilde = Ytilde,
#'                            Dtilde = Dtilde)
#'
#'     Post_plm <- plm(Ytilde ~ -1 + Dtilde, data = data_res, model = "pooling", index = c("id", "time"))
#'
#'     coefs <- coef(Post_plm)
#'     se_corrected <- sqrt(vcovHC(Post_plm, type = "HC0", method = "arellano")) * sqrt(N * T / (N*T - T*G_unit))
#'     t_values_corrected <- coefs / se_corrected
#'
#'     df <- Post_plm$df.residual
#'     p_values_corrected <- 2 * pt(-abs(t_values_corrected), df)
#'
#'     summary_table_correct <- data.frame(
#'       Estimate = coefs,
#'       `Std. Error corrected` = se_corrected,
#'       `t-value corrected` = t_values_corrected,
#'       `Pr(>|t|) corrected` = p_values_corrected
#'     )
#'
#'     colnames(summary_table_correct) = c('Estimate', 'Std. Error corrected', 't-value corrected', 'Pr(>|t|) corrected')
#'
#'     # Return cluster counts as well
#'     return(list(
#'       fit_summary = summary(fit),
#'       G_unit = G_unit,
#'       unit_group = unit_group,
#'       post_plm_summary = summary(Post_plm),
#'       estimate_corrected = summary_table_correct,
#'       summary_table = summary_table_correct
#'     ))
#'
#'   }
#' }
#'
#'
#'
#' # HP_estimate_double <- function(data, y_col = NULL, covariate_cols = NULL,
#' #                                id_col = "id", time_col = "time",
#' #                                index = c("id", "time"), method_auto = "dynamicTreeCut", cluster_method = 'kmeans', pre_cluster = FALSE, Du_pre, Dv_pre, unit_cluster = NULL, time_cluster = NULL, gamma = 1, cc = 0, dim_moment = 1) {
#' #
#' #
#' #   if (is.null(covariate_cols)) {
#' #     covariate_cols <- setdiff(colnames(data), c(id_col, time_col, y_col))
#' #   }
#' #
#' #   ids <- sort(unique(data[[id_col]]))
#' #   times <- sort(unique(data[[time_col]]))
#' #   N <- length(ids)
#' #   T <- length(times)
#' #   K <- length(covariate_cols)
#' #
#' #   # Initialize y and X
#' #   y <- matrix(NA_real_, nrow = N, ncol = T)
#' #   X <- array(NA_real_, dim = c(N, T, K))
#' #
#' #   for (i in seq_len(nrow(data))) {
#' #     id_idx <- which(ids == data[[id_col]][i])
#' #     time_idx <- which(times == data[[time_col]][i])
#' #     y[id_idx, time_idx] <- data[[y_col]][i]
#' #     for (k in seq_along(covariate_cols)) {
#' #       X[id_idx, time_idx, k] <- data[[covariate_cols[k]]][i]
#' #     }
#' #   }
#' #
#' #   # Cluster (output indicators assumed one-hot encoded)
#' #   X_norm <- array(NA_real_, dim = dim(X))  # Same shape: (N, T, K)
#' #   for (k in 1:dim(X)[3]) {
#' #     X_slice <- X[,,k]
#' #     X_norm[,,k] <- scale(X_slice)
#' #   }
#' #   y_norm = scale(y)
#' #   if (cluster_method == 'kmeans'){
#' #
#' #     # kmeans
#' #     X_list = lapply(seq_len(dim(X)[3]), function(k) X[ , , k])
#' #     if (pre_cluster == FALSE){
#' #       # cluster
#' #       if (is.null(unit_cluster) == 1 & is.null(time_cluster) == 1){
#' #         clusteri <- cluster_general(y, X_list, N, T, init = 50, type = "long", gamma = gamma, cc = cc, dim_moment = dim_moment   )
#' #         G_unit <- clusteri$clusters
#' #         klong <- clusteri$res
#' #
#' #
#' #         clustert <- cluster_general(y, X_list, N, T, init = 50, type = "tall", gamma = gamma, cc = cc, dim_moment = dim_moment  )
#' #         G_time <- clustert$clusters
#' #         ktall <- clustert$res
#' #       }else{
#' #         clusteri <- cluster_general(y, X_list, N, T, init = 50, type = "long" , groups = c(floor(unit_cluster)) )
#' #         G_unit <- clusteri$clusters
#' #         klong <- clusteri$res
#' #
#' #
#' #         clustert <- cluster_general(y, X_list, N, T, init = 50, type = "tall" , groups = c(floor(time_cluster)))
#' #         G_time <- clustert$clusters
#' #         ktall <- clustert$res
#' #       }
#' #
#' #       Du <- matrix(0, N, G_unit)
#' #       Dv <- matrix(0, T, G_time)
#' #       Diagu<- matrix(0, G_unit, G_unit)
#' #       Diagv<- matrix(0, G_time, G_time)
#' #
#' #       for (j in seq_len(G_unit)) {
#' #         Du[, j] <- as.numeric(klong$cluster == j)
#' #         Diagu[j,j] <- 1/sum(Du[,j])
#' #       }
#' #
#' #       for (j in seq_len(G_time)) {
#' #         Dv[, j] <- as.numeric(ktall$cluster == j)
#' #         Diagv[j,j] <- 1/sum(Dv[,j])
#' #       }
#' #     }else if (pre_cluster == TRUE){
#' #       Du = Du_pre
#' #       Dv = Dv_pre
#' #       Diagu<- matrix(0, G_unit, G_unit)
#' #       Diagv<- matrix(0, G_time, G_time)
#' #
#' #       for (j in seq_len(G_unit)) {
#' #         Diagu[j,j] <- 1/sum(Du[,j])
#' #       }
#' #
#' #       for (j in seq_len(G_time)) {
#' #         Diagv[j,j] <- 1/sum(Dv[,j])
#' #       }
#' #       G_unit = dim(Du)[2]
#' #       G_time = dim(Dv)[2]
#' #     }
#' #   }else if (cluster_method == 'hierarchical'){
#' #
#' #     # heriachical
#' #     if (is.null(unit_cluster) == 1 & is.null(time_cluster) == 1){
#' #       res_unit <- cluster_Hierarchical(y_norm, X_norm, type = "unit", method_auto = method_auto)
#' #       res_time <- cluster_Hierarchical(y_norm, X_norm, type = "time", method_auto = method_auto)
#' #     }else{
#' #       res_unit <- cluster_Hierarchical(y_norm, X_norm, type = "unit", cluster = unit_cluster)
#' #       res_time <- cluster_Hierarchical(y_norm, X_norm, type = "time", cluster = time_cluster)
#' #     }
#' #
#' #     Du <- res_unit$indicator    # N x G_unit
#' #     Dv <- res_time$indicator    # T x G_time
#' #
#' #     G_unit = res_unit$G
#' #     G_time = res_time$G
#' #
#' #     Diagu<- matrix(0, G_unit, G_unit)
#' #     Diagv<- matrix(0, G_time, G_time)
#' #
#' #     for (j in seq_len(G_unit)) {
#' #       Diagu[j,j] <- 1/sum(Du[,j])
#' #     }
#' #
#' #     for (j in seq_len(G_time)) {
#' #       Diagv[j,j] <- 1/sum(Dv[,j])
#' #     }
#' #   }
#' #
#' #   unit_group <- max.col(Du)
#' #   time_group <- max.col(Dv)
#' #
#' #   # Projection matrices
#' #   Mu <- diag(N) - Du %*% Diagu %*% t(Du)
#' #   Mv <- diag(T) - Dv %*% Diagv %*% t(Dv)
#' #   # Mu <- diag(N) - Du %*% solve(t(Du) %*% Du) %*% t(Du)
#' #   # Mv <- diag(T) - Dv %*% solve(t(Dv) %*% Dv) %*% t(Dv)
#' #
#' #   # Stack y and X into Z: N x T x (K+1)
#' #   Z <- array(NA_real_, dim = c(N, T, K + 1))
#' #   Z[,,1] <- y
#' #   for (k in 1:K) Z[,,k+1] <- X[,,k]
#' #
#' #   # Projection method: demean unit and time clusters slice-wise
#' #   Z_proj <- array(NA_real_, dim = c(N, T, K + 1))
#' #   for (k in 1:(K+1)) {
#' #     Z_proj[,,k] <- Mu %*% Z[,,k] %*% Mv
#' #   }
#' #
#' #
#' #   y_trans <- Z_proj[,,1] # Z_proj_final_array[,,1]
#' #   X_trans <- Z_proj[,,2:(K+1)] # Z_proj_final_array[,,2:(K+1)]
#' #
#' #   tY_vector <- as.vector(y_trans)  # (N*T)
#' #   tX_matrix <- matrix(aperm(X_trans, c(1, 2, 3)), nrow = N*T, ncol = K)
#' #
#' #   d_vector <- tX_matrix[, 1]
#' #   x_matrix <- tX_matrix[, -1, drop = FALSE]
#' #   # Separate transformed y and X
#' #
#' #   trans <- data.frame(y = tY_vector,
#' #                       D = d_vector,
#' #                       x_matrix)
#' #
#' #   # Run rlassoEffect with correct inputs
#' #   fit <- rlassoEffect(x = x_matrix, y = tY_vector, d = d_vector, method = "double selection")
#' #
#' #   trans <- data.frame(y = tY_vector,
#' #                       D = d_vector,
#' #                       x_matrix)
#' #
#' #   lasso.Y <- rlasso(y ~ . - D - 1, data = trans)
#' #   Ytilde <- lasso.Y$residuals
#' #
#' #   lasso.D <- rlasso(D ~ . - 1, data = trans[, -1])
#' #   Dtilde <- lasso.D$residuals
#' #
#' #   data_res <- data.frame(id = data[[index[1]]],
#' #                          time = data[[index[2]]],
#' #                          Ytilde = Ytilde,
#' #                          Dtilde = Dtilde)
#' #
#' #   Post_plm <- plm(Ytilde ~ -1 + Dtilde, data = data_res, model = "pooling", index = c("id", "time"))
#' #
#' #   coefs <- coef(Post_plm)
#' #   se_corrected <- sqrt(vcovHC(Post_plm, type = "HC0", method = "arellano")) * sqrt(N * T / (N*T - N*G_time - T*G_unit))
#' #   t_values_corrected <- coefs / se_corrected
#' #
#' #   df <- Post_plm$df.residual
#' #   p_values_corrected <- 2 * pt(-abs(t_values_corrected), df)
#' #
#' #   summary_table_correct <- data.frame(
#' #     Estimate = coefs,
#' #     `Std. Error corrected` = se_corrected,
#' #     `t-value corrected` = t_values_corrected,
#' #     `Pr(>|t|) corrected` = p_values_corrected
#' #   )
#' #
#' #   colnames(summary_table_correct) = c('Estimate', 'Std. Error corrected', 't-value corrected', 'Pr(>|t|) corrected')
#' #
#' #   # Return cluster counts as well
#' #   return(list(
#' #     fit_summary = summary(fit),
#' #     G_unit = G_unit,
#' #     G_time = G_time,
#' #     unit_group = unit_group,
#' #     time_group = time_group,
#' #     post_plm_summary = summary(Post_plm),
#' #     estimate_corrected = summary_table_correct,
#' #     summary_table = summary_table_correct
#' #   ))
#' # }
#'
#'   #' @export
#'   Niave_estimate <- function(data, y_col = NULL, covariate_cols = NULL,
#'                              id_col = "id", time_col = "time", index = c("id", "time")) {
#'
#'
#'     if (is.null(covariate_cols)) {
#'       covariate_cols <- setdiff(colnames(data), c(id_col, time_col, y_col))
#'     }
#'
#'     ids <- sort(unique(data[[id_col]]))
#'     times <- sort(unique(data[[time_col]]))
#'     N <- length(ids)
#'     T <- length(times)
#'     K <- length(covariate_cols)
#'
#'     # Initialize y and X
#'     y <- matrix(NA_real_, nrow = N, ncol = T)
#'     X <- array(NA_real_, dim = c(N, T, K))
#'
#'     for (i in seq_len(nrow(data))) {
#'       id_idx <- which(ids == data[[id_col]][i])
#'       time_idx <- which(times == data[[time_col]][i])
#'       y[id_idx, time_idx] <- data[[y_col]][i]
#'       for (k in seq_along(covariate_cols)) {
#'         X[id_idx, time_idx, k] <- data[[covariate_cols[k]]][i]
#'       }
#'     }
#'
#'     # Stack y and X into Z: N x T x (K+1)
#'     Z <- array(NA_real_, dim = c(N, T, K + 1))
#'     Z[,,1] <- y
#'     for (k in 1:K) Z[,,k+1] <- X[,,k]
#'
#'     # Projection method: demean unit and time clusters slice-wise
#'     Z_proj <- array(NA_real_, dim = c(N, T, K + 1))
#'     for (k in 1:(K+1)) {
#'       Z_proj[,,k] <-  Z[,,k]
#'     }
#'
#'     # Separate transformed y and X
#'     y_trans <- Z_proj[,,1] # Z_proj_final_array[,,1]
#'     X_trans <- Z_proj[,,2:(K+1)] # Z_proj_final_array[,,2:(K+1)]
#'
#'     tY_vector <- as.vector(y_trans)  # (N*T)
#'     # Handle T=1 case where X_trans is 2D instead of 3D
#'     if (length(dim(X_trans)) == 2) {
#'       tX_matrix <- matrix(as.vector(t(X_trans)), nrow = N*T, ncol = K)
#'     } else {
#'       tX_matrix <- matrix(aperm(X_trans, c(1, 2, 3)), nrow = N*T, ncol = K)
#'     }
#'
#'     d_vector <- tX_matrix[, 1]
#'     x_matrix <- tX_matrix[, -1, drop = FALSE]
#'
#'     # Run rlassoEffect with correct inputs
#'     fit <- rlassoEffect(x = x_matrix, y = tY_vector, d = d_vector, method = "double selection")
#'
#'     trans <- data.frame(y = tY_vector,
#'                         D = d_vector,
#'                         x_matrix)
#'
#'     lasso.Y <- rlasso(y ~ . - D - 1, data = trans)
#'     Ytilde <- lasso.Y$residuals
#'
#'     lasso.D <- rlasso(D ~ . - 1, data = trans[, -1, drop = FALSE])
#'     Dtilde <- lasso.D$residuals
#'
#'     # Create proper panel structure for residuals
#'     # as.vector() on N x T matrix goes column-wise: unit 1 time 1, unit 2 time 1, ..., unit N time 1, unit 1 time 2, ...
#'     data_res <- data.frame(id = rep(ids, T),
#'                            time = rep(times, each = N),
#'                            Ytilde = Ytilde,
#'                            Dtilde = Dtilde)
#'
#'     # For T=1, use lm() since there's no panel structure; for T>1 use plm with pooling
#'     if (T == 1) {
#'       Post_plm <- lm(Ytilde ~ -1 + Dtilde, data = data_res)
#'     } else {
#'       Post_plm <- plm(Ytilde ~ -1 + Dtilde, data = data_res, model = "pooling", index = c("id", "time"))
#'     }
#'
#'     coefs <- coef(Post_plm)
#'
#'     # Use appropriate variance estimator based on model type
#'     if (T == 1) {
#'       robust_se <- sqrt(diag(vcovHC(Post_plm, type = "HC0")))
#'       dof_denom <- max(1, N * T - 1)
#'       se_corrected <- robust_se * sqrt(N * T / dof_denom)
#'     } else {
#'       se_corrected <- sqrt(diag(vcovHC(Post_plm, type = "HC0", method = "arellano")))
#'     }
#'     t_values_corrected <- coefs / se_corrected
#'
#'     df <- Post_plm$df.residual
#'     p_values_corrected <- 2 * pt(-abs(t_values_corrected), df)
#'
#'     summary_table_correct <- data.frame(
#'       Estimate = coefs,
#'       `Std. Error corrected` = se_corrected,
#'       `t-value corrected` = t_values_corrected,
#'       `Pr(>|t|) corrected` = p_values_corrected
#'     )
#'
#'     colnames(summary_table_correct) = c('Estimate', 'Std. Error corrected', 't-value corrected', 'Pr(>|t|) corrected')
#'
#'     # Return cluster counts as well
#'     return(list(
#'       fit_summary = summary(fit),
#'       post_plm_summary = summary(Post_plm),
#'       estimate_corrected = summary_table_correct,
#'       summary_table = summary_table_correct
#'     ))
#'   }
#'
#' #' @export
#' TWFE_estimate <- function(data, y_col = NULL, covariate_cols = NULL,
#'   id_col = "id", time_col = "time", index = c("id", "time")) {
#'
#'     if (is.null(covariate_cols)) {
#'       covariate_cols <- setdiff(colnames(data), c(id_col, time_col, y_col))
#'     }
#'
#'     ids <- sort(unique(data[[id_col]]))
#'     times <- sort(unique(data[[time_col]]))
#'     N <- length(ids)
#'     T <- length(times)
#'     K <- length(covariate_cols)
#'
#'     # Initialize y and X as balanced panel containers
#'     y <- matrix(NA_real_, nrow = N, ncol = T)
#'     X <- array(NA_real_, dim = c(N, T, K))
#'
#'     for (i in seq_len(nrow(data))) {
#'       id_idx <- which(ids == data[[id_col]][i])
#'       time_idx <- which(times == data[[time_col]][i])
#'       y[id_idx, time_idx] <- data[[y_col]][i]
#'       for (k in seq_along(covariate_cols)) {
#'         X[id_idx, time_idx, k] <- data[[covariate_cols[k]]][i]
#'       }
#'     }
#'
#'     # Residualize with respect to fixed effects before running DML.
#'     # When T = 1, two-way FE is not identified; fall back to no time demeaning.
#'     if (T == 1) {
#'       y_res <- y
#'       X_res <- X
#'     } else {
#'       Mu <- diag(N) - matrix(1 / N, nrow = N, ncol = N)
#'       Mv <- diag(T) - matrix(1 / T, nrow = T, ncol = T)
#'
#'       y_res <- Mu %*% y %*% Mv
#'       X_res <- array(NA_real_, dim = c(N, T, K))
#'       for (k in seq_len(K)) {
#'         X_res[,,k] <- Mu %*% X[,,k] %*% Mv
#'       }
#'     }
#'
#'     tY_vector <- as.vector(y_res)
#'     if (T == 1) {
#'       X_res_mat <- if (length(dim(X_res)) == 3) X_res[, 1, , drop = TRUE] else X_res
#'       tX_matrix <- matrix(as.vector(t(X_res_mat)), nrow = N * T, ncol = K)
#'     } else if (length(dim(X_res)) == 2) {
#'       tX_matrix <- matrix(as.vector(t(X_res)), nrow = N * T, ncol = K)
#'     } else {
#'       tX_matrix <- matrix(aperm(X_res, c(1, 2, 3)), nrow = N * T, ncol = K)
#'     }
#'
#'     d_vector <- tX_matrix[, 1]
#'     x_matrix <- tX_matrix[, -1, drop = FALSE]
#'
#'     # DML on residualized outcome and regressors
#'     fit <- rlassoEffect(x = x_matrix, y = tY_vector, d = d_vector, method = "double selection")
#'
#'     trans <- data.frame(
#'       y = tY_vector,
#'       D = d_vector,
#'       x_matrix
#'     )
#'
#'     lasso.Y <- rlasso(y ~ . - D - 1, data = trans)
#'     Ytilde <- lasso.Y$residuals
#'
#'     lasso.D <- rlasso(D ~ . - 1, data = trans[, -1, drop = FALSE])
#'     Dtilde <- lasso.D$residuals
#'
#'     data_res <- data.frame(
#'       id = rep(ids, T),
#'       time = rep(times, each = N),
#'       Ytilde = Ytilde,
#'       Dtilde = Dtilde
#'     )
#'
#'     if (T == 1) {
#'       Post_plm <- lm(Ytilde ~ -1 + Dtilde, data = data_res)
#'       robust_se <- sqrt(diag(vcovHC(Post_plm, type = "HC0")))
#'       dof_denom <- max(1, N * T - 1)
#'       G_time <- 1
#'     } else {
#'       Post_plm <- plm(Ytilde ~ -1 + Dtilde, data = data_res, model = "pooling", index = c("id", "time"))
#'       robust_se <- sqrt(diag(vcovHC(Post_plm, type = "HC0", method = "arellano")))
#'       dof_denom <- max(1, N * T - N - T + 1)
#'       G_time <- T
#'     }
#'
#'     G_unit <- N
#'     unit_group <- seq_len(N)
#'     time_group <- if (T == 1) rep(1, T) else seq_len(T)
#'
#'     coefs <- coef(Post_plm)
#'     se_corrected <- robust_se * sqrt(N * T / dof_denom)
#'     t_values_corrected <- coefs / se_corrected
#'
#'     df <- Post_plm$df.residual
#'     p_values_corrected <- 2 * pt(-abs(t_values_corrected), df)
#'
#'     summary_table_correct <- data.frame(
#'       Estimate = coefs,
#'       `Std. Error corrected` = se_corrected,
#'       `t-value corrected` = t_values_corrected,
#'       `Pr(>|t|) corrected` = p_values_corrected
#'     )
#'
#'     colnames(summary_table_correct) <- c("Estimate", "Std. Error corrected", "t-value corrected", "Pr(>|t|) corrected")
#'
#'     return(list(
#'       fit_summary = summary(fit),
#'       G_unit = G_unit,
#'       G_time = G_time,
#'       unit_group = unit_group,
#'       time_group = time_group,
#'       post_plm_summary = summary(Post_plm),
#'       estimate_corrected = summary_table_correct,
#'       summary_table = summary_table_correct
#'     ))
#'   }
#' # TWFE_estimate <- function(data, y_col = NULL, covariate_cols = NULL,
#' #                           id_col = "id", time_col = "time", index = c("id", "time")) {
#' #
#' #
#' #   if (is.null(covariate_cols)) {
#' #     covariate_cols <- setdiff(colnames(data), c(id_col, time_col, y_col))
#' #   }
#' #
#' #   ids <- sort(unique(data[[id_col]]))
#' #   times <- sort(unique(data[[time_col]]))
#' #   N <- length(ids)
#' #   T <- length(times)
#' #   K <- length(covariate_cols)
#' #
#' #   # Initialize y and X
#' #   y <- matrix(NA_real_, nrow = N, ncol = T)
#' #   X <- array(NA_real_, dim = c(N, T, K))
#' #
#' #   for (i in seq_len(nrow(data))) {
#' #     id_idx <- which(ids == data[[id_col]][i])
#' #     time_idx <- which(times == data[[time_col]][i])
#' #     y[id_idx, time_idx] <- data[[y_col]][i]
#' #     for (k in seq_along(covariate_cols)) {
#' #       X[id_idx, time_idx, k] <- data[[covariate_cols[k]]][i]
#' #     }
#' #   }
#' #
#' #   # Cluster (output indicators assumed one-hot encoded)
#' #   X_norm <- array(NA_real_, dim = dim(X))  # Same shape: (N, T, K)
#' #   for (k in 1:dim(X)[3]) {
#' #     X_slice <- X[,,k]
#' #     X_norm[,,k] <- scale(X_slice)
#' #   }
#' #   y_norm = scale(y)
#' #
#' #   Du <- matrix(1, nrow = N, ncol = 1)    # N x G_unit
#' #   Dv <- matrix(1, nrow = T, ncol = 1)    # T x G_time
#' #
#' #   G_unit = 1
#' #   G_time = 1
#' #
#' #   unit_group <- max.col(Du)
#' #   time_group <- max.col(Dv)
#' #
#' #   # Projection matrices
#' #   Mu <- diag(N) - Du %*% solve(t(Du) %*% Du) %*% t(Du)
#' #   Mv <- diag(T) - Dv %*% solve(t(Dv) %*% Dv) %*% t(Dv)
#' #
#' #   # Stack y and X into Z: N x T x (K+1)
#' #   Z <- array(NA_real_, dim = c(N, T, K + 1))
#' #   Z[,,1] <- y
#' #   for (k in 1:K) Z[,,k+1] <- X[,,k]
#' #
#' #   # Projection method: demean unit and time clusters slice-wise
#' #   # When T=1, Mv becomes zero matrix, so skip time demeaning
#' #   Z_proj <- array(NA_real_, dim = c(N, T, K + 1))
#' #   for (k in 1:(K+1)) {
#' #     if (T == 1) {
#' #       Z_proj[,,k] <- Mu %*% Z[,,k]
#' #     } else {
#' #       Z_proj[,,k] <- Mu %*% Z[,,k] %*% Mv
#' #     }
#' #   }
#' #
#' #   # Separate transformed y and X
#' #   y_trans <- Z_proj[,,1] # Z_proj_final_array[,,1]
#' #   X_trans <- Z_proj[,,2:(K+1)] # Z_proj_final_array[,,2:(K+1)]
#' #
#' #   tY_vector <- as.vector(y_trans)  # (N*T)
#' #   # Handle T=1 case where X_trans is 2D instead of 3D
#' #   if (length(dim(X_trans)) == 2) {
#' #     tX_matrix <- matrix(as.vector(t(X_trans)), nrow = N*T, ncol = K)
#' #   } else {
#' #     tX_matrix <- matrix(aperm(X_trans, c(1, 2, 3)), nrow = N*T, ncol = K)
#' #   }
#' #
#' #   d_vector <- tX_matrix[, 1]
#' #   x_matrix <- tX_matrix[, -1, drop = FALSE]
#' #
#' #   # Run rlassoEffect with correct inputs
#' #   fit <- rlassoEffect(x = x_matrix, y = tY_vector, d = d_vector, method = "double selection")
#' #
#' #   trans <- data.frame(y = tY_vector,
#' #                       D = d_vector,
#' #                       x_matrix)
#' #
#' #   lasso.Y <- rlasso(y ~ . - D - 1, data = trans)
#' #   Ytilde <- lasso.Y$residuals
#' #
#' #   lasso.D <- rlasso(D ~ . - 1, data = trans[, -1])
#' #   Dtilde <- lasso.D$residuals
#' #
#' #   # Create proper panel structure for residuals
#' #   # as.vector() on N x T matrix goes column-wise: unit 1 time 1, unit 2 time 1, ..., unit N time 1, unit 1 time 2, ...
#' #   data_res <- data.frame(id = rep(ids, T),
#' #                          time = rep(times, each = N),
#' #                          Ytilde = Ytilde,
#' #                          Dtilde = Dtilde)
#' #
#' #   # For T=1, use lm() since there's no panel structure; for T>1 use plm with pooling
#' #   if (T == 1) {
#' #     Post_plm <- lm(Ytilde ~ -1 + Dtilde, data = data_res)
#' #   } else {
#' #     Post_plm <- plm(Ytilde ~ -1 + Dtilde, data = data_res, model = "pooling", index = c("id", "time"))
#' #   }
#' #
#' #   coefs <- coef(Post_plm)
#' #
#' #   # Use appropriate variance estimator based on model type
#' #   if (T == 1) {
#' #     se_corrected <- sqrt(diag(vcovHC(Post_plm, type = "HC0"))) * sqrt(N * T / (N*T - N*G_time - T*G_unit))
#' #   } else {
#' #     se_corrected <- sqrt(diag(vcovHC(Post_plm, type = "HC0", method = "arellano"))) * sqrt(N * T / (N*T - N*G_time - T*G_unit))
#' #   }
#' #
#' #   t_values_corrected <- coefs / se_corrected
#' #
#' #   df <- Post_plm$df.residual
#' #   p_values_corrected <- 2 * pt(-abs(t_values_corrected), df)
#' #
#' #   summary_table_correct <- data.frame(
#' #     Estimate = coefs,
#' #     `Std. Error corrected` = se_corrected,
#' #     `t-value corrected` = t_values_corrected,
#' #     `Pr(>|t|) corrected` = p_values_corrected
#' #   )
#' #
#' #   colnames(summary_table_correct) = c('Estimate', 'Std. Error corrected', 't-value corrected', 'Pr(>|t|) corrected')
#' #
#' #   # Return cluster counts as well
#' #   return(list(
#' #     fit_summary = summary(fit),
#' #     G_unit = G_unit,
#' #     G_time = G_time,
#' #     unit_group = unit_group,
#' #     time_group = time_group,
#' #     post_plm_summary = summary(Post_plm),
#' #     estimate_corrected = summary_table_correct,
#' #     summary_table = summary_table_correct
#' #   ))
#' # }
#'
#' # HP_estimate_twoway <- function(data, y_col = NULL, covariate_cols = NULL,
#' #                         id_col = "id", time_col = "time",
#' #                         index = c("id", "time"), method_auto = "dynamicTreeCut", cluster = NULL) {
#' #
#' #
#' #   if (is.null(covariate_cols)) {
#' #     covariate_cols <- setdiff(colnames(data), c(id_col, time_col, y_col))
#' #   }
#' #
#' #   ids <- sort(unique(data[[id_col]]))
#' #   times <- sort(unique(data[[time_col]]))
#' #   N <- length(ids)
#' #   T <- length(times)
#' #   K <- length(covariate_cols)
#' #
#' #   # Initialize y and X
#' #   y <- matrix(NA_real_, nrow = N, ncol = T)
#' #   X <- array(NA_real_, dim = c(N, T, K))
#' #
#' #   for (i in seq_len(nrow(data))) {
#' #     id_idx <- which(ids == data[[id_col]][i])
#' #     time_idx <- which(times == data[[time_col]][i])
#' #     y[id_idx, time_idx] <- data[[y_col]][i]
#' #     for (k in seq_along(covariate_cols)) {
#' #       X[id_idx, time_idx, k] <- data[[covariate_cols[k]]][i]
#' #     }
#' #   }
#' #
#' #   # Cluster (output indicators assumed one-hot encoded)
#' #   X_norm <- array(NA_real_, dim = dim(X))  # Same shape: (N, T, K)
#' #   for (k in 1:dim(X)[3]) {
#' #     X_slice <- X[,,k]
#' #     X_norm[,,k] <- scale(X_slice)
#' #   }
#' #   y_norm = scale(y)
#' #   res_unit <- cluster_Hierarchical(y_norm, X_norm, type = "unit", method_auto = method_auto)
#' #   res_time <- cluster_Hierarchical(y_norm, X_norm, type = "time", method_auto = method_auto)
#' #   res_covar <- cluster_Hierarchical(y_norm, X_norm, type = "covariate", method_auto = method_auto, cluster = cluster)
#' #
#' #   Dv <- res_unit$indicator    # N x G_unit
#' #   Du <- res_time$indicator    # T x G_time
#' #   Dc <- res_covar$indicator   # (K+1) x G_covar
#' #
#' #   G_unit = res_unit$G
#' #   G_time = res_time$G
#' #   G_covar = res_covar$G
#' #
#' #   unit_group <- max.col(Dv)
#' #   time_group <- max.col(Du)
#' #   covar_group <- max.col(Dc)
#' #
#' #   # Projection matrices
#' #   Mv <- diag(N) - Dv %*% solve(t(Dv) %*% Dv) %*% t(Dv)
#' #   Mu <- diag(T) - Du %*% solve(t(Du) %*% Du) %*% t(Du)
#' #   Mc <- diag(K+1) - Dc %*% solve(t(Dc) %*% Dc) %*% t(Dc)
#' #
#' #   # Stack y and X into Z: N x T x (K+1)
#' #   Z <- array(NA_real_, dim = c(N, T, K + 1))
#' #   Z[,,1] <- y
#' #   for (k in 1:K) Z[,,k+1] <- X[,,k]
#' #
#' #   # Projection method: demean unit and time clusters slice-wise
#' #   Z_proj <- array(NA_real_, dim = c(N, T, K + 1))
#' #   for (k in 1:(K+1)) {
#' #     Z_proj[,,k] <- Mv %*% Z[,,k] %*% Mu
#' #   }
#' #   # Separate transformed y and X
#' #   y_trans <- Z_proj[,,1] # Z_proj_final_array[,,1]
#' #   X_trans <- Z_proj[,,2:(K+1)] # Z_proj_final_array[,,2:(K+1)]
#' #
#' #   tY_vector <- as.vector(y_trans)  # (N*T)
#' #   tX_matrix <- matrix(aperm(X_trans, c(1, 2, 3)), nrow = N*T, ncol = K)
#' #
#' #   d_vector <- tX_matrix[, 1]
#' #   x_matrix <- tX_matrix[, -1, drop = FALSE]
#' #
#' #   # Run rlassoEffect with correct inputs
#' #   fit <- rlassoEffect(x = x_matrix, y = tY_vector, d = d_vector, method = "partialling out")
#' #
#' #
#' #   trans <- data.frame(y = tY_vector,
#' #                       D = d_vector,
#' #                       x_matrix)
#' #
#' #   lasso.Y <- rlasso(y ~ . - D - 1, data = trans)
#' #   Ytilde <- lasso.Y$residuals
#' #
#' #   lasso.D <- rlasso(D ~ . - 1, data = trans[, -1])
#' #   Dtilde <- lasso.D$residuals
#' #
#' #   data_res <- data.frame(id = data[[index[1]]],
#' #                          time = data[[index[2]]],
#' #                          Ytilde = Ytilde,
#' #                          Dtilde = Dtilde)
#' #
#' #   Post_plm <- plm(Ytilde ~ -1 + Dtilde, data = data_res, model = "pooling", index = c("id", "time"))
#' #
#' #   coefs <- coef(Post_plm)
#' #   se_corrected <- sqrt(vcovHC(Post_plm, type = "HC0", method = "arellano")) * sqrt(N * T / ((N-G_unit)*(T-G_time)))
#' #   t_values_corrected <- coefs / se_corrected
#' #
#' #   df <- Post_plm$df.residual
#' #   p_values_corrected <- 2 * pt(-abs(t_values_corrected), df)
#' #
#' #   summary_table_correct <- data.frame(
#' #     Estimate = coefs,
#' #     `Std. Error corrected` = se_corrected,
#' #     `t-value corrected` = t_values_corrected,
#' #     `Pr(>|t|) corrected` = p_values_corrected
#' #   )
#' #
#' #   summary_table <- summary(Post_plm)
#' #   summary_table$coefficients <- cbind(summary_table_correct, summary_table$coefficients)
#' #
#' #   # Return cluster counts as well
#' #   return(list(
#' #     fit_summary = summary(fit),
#' #     G_unit = G_unit,
#' #     G_time = G_time,
#' #     unit_group = unit_group,
#' #     time_group = time_group,
#' #     post_plm_summary = summary(Post_plm),
#' #     estimate_corrected = summary_table_correct,
#' #     summary_table = summary_table
#' #   ))
#' # }
#' #
#' # HP_estimate_double_twoway <- function(data, y_col = NULL, covariate_cols = NULL,
#' #                                id_col = "id", time_col = "time",
#' #                                index = c("id", "time"), method_auto = "dynamicTreeCut", cluster = NULL) {
#' #
#' #
#' #   if (is.null(covariate_cols)) {
#' #     covariate_cols <- setdiff(colnames(data), c(id_col, time_col, y_col))
#' #   }
#' #
#' #   ids <- sort(unique(data[[id_col]]))
#' #   times <- sort(unique(data[[time_col]]))
#' #   N <- length(ids)
#' #   T <- length(times)
#' #   K <- length(covariate_cols)
#' #
#' #   # Initialize y and X
#' #   y <- matrix(NA_real_, nrow = N, ncol = T)
#' #   X <- array(NA_real_, dim = c(N, T, K))
#' #
#' #   for (i in seq_len(nrow(data))) {
#' #     id_idx <- which(ids == data[[id_col]][i])
#' #     time_idx <- which(times == data[[time_col]][i])
#' #     y[id_idx, time_idx] <- data[[y_col]][i]
#' #     for (k in seq_along(covariate_cols)) {
#' #       X[id_idx, time_idx, k] <- data[[covariate_cols[k]]][i]
#' #     }
#' #   }
#' #
#' #   # Cluster (output indicators assumed one-hot encoded)
#' #   X_norm <- array(NA_real_, dim = dim(X))  # Same shape: (N, T, K)
#' #   for (k in 1:dim(X)[3]) {
#' #     X_slice <- X[,,k]
#' #     X_norm[,,k] <- scale(X_slice)
#' #   }
#' #   y_norm = scale(y)
#' #   res_unit <- cluster_Hierarchical(y_norm, X_norm, type = "unit", method_auto = method_auto)
#' #   res_time <- cluster_Hierarchical(y_norm, X_norm, type = "time", method_auto = method_auto)
#' #   res_covar <- cluster_Hierarchical(y_norm, X_norm, type = "covariate", method_auto = method_auto, cluster = cluster)
#' #
#' #   Dv <- res_unit$indicator    # N x G_unit
#' #   Du <- res_time$indicator    # T x G_time
#' #   Dc <- res_covar$indicator   # (K+1) x G_covar
#' #
#' #   G_unit = res_unit$G
#' #   G_time = res_time$G
#' #   G_covar = res_covar$G
#' #
#' #   unit_group <- max.col(Dv)
#' #   time_group <- max.col(Du)
#' #   covar_group <- max.col(Dc)
#' #
#' #   # Projection matrices
#' #   Mv <- diag(N) - Dv %*% solve(t(Dv) %*% Dv) %*% t(Dv)
#' #   Mu <- diag(T) - Du %*% solve(t(Du) %*% Du) %*% t(Du)
#' #   Mc <- diag(K+1) - Dc %*% solve(t(Dc) %*% Dc) %*% t(Dc)
#' #
#' #   # Stack y and X into Z: N x T x (K+1)
#' #   Z <- array(NA_real_, dim = c(N, T, K + 1))
#' #   Z[,,1] <- y
#' #   for (k in 1:K) Z[,,k+1] <- X[,,k]
#' #
#' #   # Projection method: demean unit and time clusters slice-wise
#' #   Z_proj <- array(NA_real_, dim = c(N, T, K + 1))
#' #   for (k in 1:(K+1)) {
#' #     Z_proj[,,k] <- Mv %*% Z[,,k] %*% Mu
#' #   }
#' #
#' #   # Flatten Z_proj for covariate demeaning
#' #   Z_proj_mat <- matrix(NA_real_, nrow = N * T, ncol = K + 1)
#' #   for (k in 1:(K+1)) {
#' #     Z_proj_mat[,k] <- as.vector(t(Z_proj[,,k]))
#' #   }
#' #
#' #   # Covariate demeaning via projection
#' #   Z_proj_final <- Z_proj_mat %*% Mc
#' #
#' #   # Reshape back to array
#' #   Z_proj_final_array <- array(NA_real_, dim = c(N, T, K + 1))
#' #   for (k in 1:(K+1)) {
#' #     Z_proj_final_array[,,k] <- matrix(Z_proj_final[,k], nrow = N, ncol = T, byrow = TRUE)
#' #   }
#' #
#' #   first_val <- covar_group[1]
#' #   count_first <- sum(covar_group == first_val)
#' #
#' #   if (count_first == 1) {
#' #     # Find the most frequent value in covar_group
#' #     most_freq_val <- as.numeric(names(sort(table(covar_group), decreasing = TRUE)[1]))
#' #     covar_group[1] <- most_freq_val
#' #   }
#' #   #covar_group[1:length(covar_group)]=1
#' #   #Z_proj = group_demean_formula_cpp(Z, unit_group, time_group, covar_group)
#' #   y_trans <- Z_proj[,,1] # Z_proj_final_array[,,1]
#' #   X_trans <- Z_proj[,,2:(K+1)] # Z_proj_final_array[,,2:(K+1)]
#' #
#' #   tY_vector <- as.vector(y_trans)  # (N*T)
#' #   tX_matrix <- matrix(aperm(X_trans, c(1, 2, 3)), nrow = N*T, ncol = K)
#' #
#' #   d_vector <- tX_matrix[, 1]
#' #   x_matrix <- tX_matrix[, -1, drop = FALSE]
#' #   # Separate transformed y and X
#' #
#' #   trans <- data.frame(y = tY_vector,
#' #                       D = d_vector,
#' #                       x_matrix)
#' #
#' #   trans = as.data.frame(trans)
#' #
#' #   dml_data <- DoubleMLData$new(
#' #     data = trans,
#' #     y_col = "y",
#' #     d_cols = "D",
#' #     x_cols = c(colnames(trans)[-c(1,2)])
#' #
#' #   )
#' #
#' #   # Define LASSO machine learning learners for nuisance parameter estimation
#' #   ml_l <- lrn("regr.cv_glmnet", s = "lambda.min")  # Outcome regression model
#' #   ml_m <- lrn("regr.cv_glmnet", s = "lambda.min")  # Treatment model (if applicable)
#' #
#' #   # Fit the Double Machine Learning model for treatment effect estimation
#' #   dml_plr <- DoubleMLPLR$new(dml_data, ml_l = ml_l, ml_m = ml_m)
#' #
#' #   # Fit the model to estimate the causal effect
#' #   dml_plr$fit(store_predictions=TRUE)
#' #   g_hat <- dml_plr$predictions$ml_l
#' #   m_hat <- dml_plr$predictions$ml_m
#' #
#' #   # Step 2: Compute residuals
#' #   Ytilde <- tY_vector - g_hat
#' #   Dtilde <- d_vector - m_hat
#' #   # Double ML
#' #   data_res = data.frame(id = data[[index[1]]], time = data[[index[2]]], Ytilde, Dtilde)
#' #   Post_plm = plm(Ytilde ~ -1 + Dtilde, data = data_res, model = "pooling", index=c("id", "time"))
#' #
#' #   coefs <- coef(Post_plm)
#' #   se_corrected <- sqrt(vcovHC(Post_plm, type = "HC0", method = "arellano")) * sqrt(N * T / ((N-G_unit)*(T-G_time)))
#' #   t_values_corrected <- coefs / se_corrected
#' #
#' #   # Calculate p-values from t-distribution for each coefficient
#' #   df <- Post_plm$df.residual  # degrees of freedom
#' #   p_values_corrected <- 2 * pt(-abs(t_values_corrected), df)
#' #
#' #   summary_table_correct <- data.frame(
#' #     Estimate = coefs,
#' #     SE_Corrected = se_corrected,
#' #     t_value_Corrected = t_values_corrected,
#' #     p_value_Corrected = p_values_corrected
#' #   )
#' #   colnames(summary_table_correct) = c('Estimate', 'Std. Error corrected', 't-value corrected', 'Pr(>|t|) corrected')
#' #   summary_table = summary(Post_plm)
#' #   summary_table$coefficients = cbind(summary_table_correct, summary_table$coefficients)
#' #
#' #   # Return cluster counts as well
#' #   return(list(
#' #     fit_summary = dml_plr$summary(),
#' #     G_unit = G_unit,
#' #     G_time = G_time,
#' #     G_covar = G_covar,
#' #     unit_group = unit_group,
#' #     time_group = time_group,
#' #     covar_group = covar_group,
#' #     post_plm_summary = summary(Post_plm),
#' #     estimate_corrected = summary_table_correct,
#' #     summary_table = summary_table
#' #   ))
#' # }
#'
#' #' @export
#' compute_u_hat <- function(z_array, unit_clusters, time_clusters, covar_clusters) {
#'   N <- dim(z_array)[1]
#'   T <- dim(z_array)[2]
#'   K <- dim(z_array)[3]
#'
#'   u_hat <- array(0, dim = c(N, T, K))
#'
#'   for (i in 1:N) {
#'     for (t in 1:T) {
#'       for (k in 1:K) {
#'         g_i <- unit_clusters[i]
#'         m_t <- time_clusters[t]
#'         l_k <- covar_clusters[k]
#'
#'         # z_{itk}
#'         zitk <- z_array[i, t, k]
#'
#'         # bar_z_{g_i t k}
#'         group_i_indices <- which(unit_clusters == g_i)
#'         bar_g_i_t_k <- mean(z_array[group_i_indices, t, k])
#'
#'         # bar_z_{i m_t k}
#'         time_m_indices <- which(time_clusters == m_t)
#'         bar_i_m_t_k <- mean(z_array[i, time_m_indices, k])
#'
#'         # bar_z_{i t l_k}
#'         covar_l_indices <- which(covar_clusters == l_k)
#'         bar_i_t_l_k <- mean(z_array[i, t, covar_l_indices])
#'
#'         # bar_z_{g_i m_t k}
#'         bar_g_i_m_t_k <- mean(z_array[group_i_indices, time_m_indices, k])
#'
#'         # bar_z_{g_i t l_k}
#'         bar_g_i_t_l_k <- mean(z_array[group_i_indices, t, covar_l_indices])
#'
#'         # bar_z_{i m_t l_k}
#'         bar_i_m_t_l_k <- mean(z_array[i, time_m_indices, covar_l_indices])
#'
#'         # Final u_hat_{itk}
#'         u_hat[i, t, k] <- 3 * zitk -
#'           2 * bar_g_i_t_k -
#'           2 * bar_i_m_t_k -
#'           2 * bar_i_t_l_k +
#'           bar_g_i_m_t_k +
#'           bar_g_i_t_l_k +
#'           bar_i_m_t_l_k
#'       }
#'     }
#'   }
#'
#'   return(u_hat)
#' }
#'
#' #' @export
#' cluster_Hierarchical <- function(y, X, link = "average", threshold = NULL, cluster = NULL,
#'                                  type = c("unit", "time", "covariate"),
#'                                  method_auto = c("none", "silhouette", "gap", "dynamicTreeCut"),
#'                                  data_for_gap = NULL, max_k = 10, deepSplit = TRUE, minClusterSize = 1, pamStage = FALSE) {
#'   type <- match.arg(type)
#'   method_auto <- match.arg(method_auto)
#'
#'   N <- nrow(y)
#'   T <- ncol(y)
#'   K <- dim(X)[3]
#'
#'   # Combine y and X into N x T x (K+1)
#'   combined_array <- array(0, dim = c(N, T, K + 1))
#'   combined_array[,,1] <- y
#'   combined_array[,,2:(K+1)] <- X
#'
#'   # Compute distance matrix based on type
#'   if (type == "unit") {
#'     dist_mat <- pseudo_dist_unit(combined_array)
#'   } else if (type == "time") {
#'     dist_mat <- pseudo_dist_time(combined_array)
#'   } else if (type == "covariate") {
#'     dist_mat <- pseudo_dist_covariate(combined_array)
#'   } else {
#'     stop("Invalid 'type' argument.")
#'   }
#'
#'   dist_obj <- as.dist(dist_mat)
#'   hc <- hclust(dist_obj, method = link)
#'
#'   if (!is.null(cluster)) {
#'     G <- cluster
#'     clusters <- cutree(hc, k = G)
#'
#'   } else if (!is.null(threshold)) {
#'     clusters <- cutree(hc, h = threshold)
#'     G <- length(unique(clusters))
#'
#'   } else if (method_auto != "none") {
#'     max_k <- min(max_k, ifelse(type == "unit", N, ifelse(type == "time", T, K + 1)) - 1)
#'
#'     if (method_auto == "gap") {
#'       if (is.null(data_for_gap)) {
#'         stop("For method_auto = 'gap', please provide 'data_for_gap' matrix.")
#'       }
#'       gap_fun <- function(x, k) {
#'         dist_x <- dist(x)
#'         hc_x <- hclust(dist_x, method = link)
#'         clust <- cutree(hc_x, k = k)
#'         list(cluster = clust)  # must return list with $cluster
#'       }
#'       gap_stat <- cluster::clusGap(data_for_gap, FUN = gap_fun, K.max = max_k, B = 50)
#'       G <- maxSE(gap_stat$Tab[, "gap"], gap_stat$Tab[, "SE.sim"], method = "firstSEmax")
#'       clusters <- cutree(hc, k = G)
#'
#'     } else if (method_auto == "silhouette") {
#'       sil_scores <- numeric(max_k)
#'       sil_scores[1] <- NA # silhouette not defined for k=1
#'       for (k in 2:max_k) {
#'         clust_try <- cutree(hc, k = k)
#'         ss <- silhouette(clust_try, dist_obj)
#'         sil_scores[k] <- mean(ss[, 3])
#'       }
#'       G <- which.max(sil_scores)
#'       clusters <- cutree(hc, k = G)
#'
#'     } else if (method_auto == "dynamicTreeCut") {
#'       clusters <- cutreeDynamic(dendro = hc, distM = as.matrix(dist_obj),
#'                                 deepSplit = deepSplit, minClusterSize = minClusterSize, pamStage = pamStage)
#'       G <- length(unique(clusters[clusters > 0]))  # exclude noise (0)
#'       noise_idx <- which(clusters == 0)
#'       if (length(noise_idx) > 0) {
#'         clusters[noise_idx] <- G + 1
#'         G <- G + 1
#'       }
#'     }
#'
#'   } else {
#'     stop("Must provide either 'cluster', 'threshold', or set 'method_auto' != 'none'.")
#'   }
#'
#'   # Format clusters and indicator matrix (like kmeans())
#'   clusters <- as.integer(factor(clusters))
#'   G <- length(unique(clusters))
#'   indicator <- sapply(1:G, function(g) as.integer(clusters == g))
#'   colnames(indicator) <- paste0("Cluster", 1:G)
#'
#'   dim_cluster <- switch(type,
#'                         unit = N,
#'                         time = T,
#'                         covariate = K + 1)
#'
#'   rownames(indicator) <- paste0(type, "_", 1:dim_cluster)
#'
#'   list(G = G, clusters = clusters, indicator = indicator)
#' }
#'
#' #' @export
#' cluster_general <- function(Y, X_list, N, T, init = 100, type = "long", groups = NULL, cc = 0, gamma = 1, dim_moment = 1) {
#'   dimtheta <- length(X_list)
#'   mdim <- dim_moment
#'   X_list <- c(list(Y), X_list)
#'   K = dimtheta
#'   if (type == "long") {
#'     mom_i <- c()
#'     ## --- Unit-side moments ---
#'     X_av <- matrix(0, nrow = N, ncol = dim_moment)  # N x dim_moment
#'
#'     if (dim_moment == 1) {
#'       covar_groups <- list(seq_len(length(X_list)))
#'     } else {
#'       covar_groups <- split(
#'         seq_len(length(X_list)),
#'         cut(seq_len(length(X_list)), dim_moment, labels = FALSE)
#'       )
#'     }
#'
#'     X_mat <- function(X) {
#'       if (is.null(dim(X))) {
#'         matrix(X, nrow = N, ncol = T)
#'       } else {
#'         X
#'       }
#'     }
#'
#'     for (p in 1:dim_moment) {
#'       # For each power
#'       X_pow <- matrix(0, nrow = N, ncol = T)  # N x T
#'
#'       for (t in 1:T) {
#'         # Collect all covariates at time t
#'         X_t <- do.call(cbind, lapply(X_list, function(X) X_mat(X)[, t]))  # N x (K+1)
#'         # Average the covariates in the p-th group
#'         cols <- covar_groups[[p]]
#'         X_pow[, t] <- rowMeans(X_t[, cols, drop = FALSE])
#'         #     if ( p == 1){
#'         #     X_pow[, t] <- rowMeans(X_t)             # (1/K) sum_k x_itk^p
#'         #     }else if (p == 2){
#'         #       X_pow[, t] <- rowMeans(tanh(X_t))
#'         #     }else if (p == 3){
#'         #       X_pow[, t] <- rowMeans((X_t/(1+abs(X_t))))
#'         #     }else{
#'         #       X_pow[, t] <- rowMeans(X_t^p)
#'         #   }
#'       }
#'
#'       # Average over time
#'       X_av[, p] <- rowMeans(X_pow)               # (1/T) sum_t (1/K) sum_k x_itk^p
#'     }
#'
#'     mom_i <- cbind(mom_i, X_av)
#'
#'     # mom_micro <- c()
#'     mom_micro <- c()
#'     X_cbind <- do.call(cbind, lapply(X_list, function(X) X_mat(X)))
#'     mom_micro <- cbind(mom_micro, X_cbind)
#'     # for (p in 1:dim_moment) {
#'     #   if ( p == 1){
#'     #     mom_micro <- cbind(mom_micro, X_cbind^p)             # (1/K) sum_k x_itk^p
#'     #   }else if (p == 2){
#'     #     mom_micro <- cbind(mom_micro, tanh(X_cbind))
#'     #   }else if (p == 3){
#'     #     mom_micro <- cbind(mom_micro, (X_cbind/(1+abs(X_cbind))))
#'     #   }else{
#'     #   mom_micro <- cbind(mom_micro, X_cbind^p)
#'     # }
#'     # }
#'     group_lengths <- sapply(covar_groups, length)
#'
#'     # new sizes
#'     new_lengths <- group_lengths * T
#'
#'     # full sequence after extension
#'     full_seq <- seq_len((K+1) * T)
#'
#'     # split into contiguous blocks
#'     col_s <- split(full_seq, rep(seq_len(dim_moment), new_lengths))
#'     ## --- Demean & rescale ---
#'     for (j in 1:mdim) {
#'       col_mean <- mean(mom_i[, j])
#'       col_sd <- sd(mom_i[, j])
#'       if (col_sd > 0) {
#'         cols <-  as.vector(col_s[[j]])
#'         mom_micro[, cols] <- (mom_micro[, cols] - col_mean) / col_sd
#'         mom_i[, j] <- (mom_i[, j] - col_mean) / col_sd
#'       }
#'     }
#'
#'
#'     # --- Variance / noise on rescaled data ---
#'     variance <-  sum(sapply(1:mdim, function(j) {
#'       norm(mom_micro[,  as.vector(col_s[[j]]) ] - mom_i[,j], type = "F")^2
#'     })) / (N * T^2 * (length(covar_groups[[1]]))^2)
#'
#'     data <- mom_i
#'     dim_size <- N
#'
#'   } else if (type == "tall") {
#'     ## --- Time-side moments ---
#'     mom_t <- c()
#'     X_av <- matrix(0, nrow = T, ncol = dim_moment)  # T x dim_moment
#'
#'     for (p in 1:dim_moment) {
#'       # For each power
#'       X_pow <- matrix(0, nrow = N, ncol = T)  # N x T
#'
#'       for (i in 1:N) {
#'         # Collect all covariates at time t
#'         X_i <- sapply(X_list, function(X) X[i, ])  # T x K
#'         # Take power p **before** averaging over K
#'         X_pow[i, ] <- rowMeans(X_i^p)             # (1/K) sum_k x_itk^p
#'       }
#'
#'       # Average over unit
#'       X_av[, p] <- colMeans(X_pow)               # (1/N) sum_i (1/K) sum_k x_itk^p
#'     }
#'
#'     mom_t <- cbind(mom_t, X_av)
#'
#'     mom_micro2 <- c()
#'     X_cbind <- do.call(cbind, lapply(X_list, t))
#'
#'     for (p in 1:dim_moment) {
#'       mom_micro2 <- cbind(mom_micro2, (X_cbind)^p)
#'     }
#'
#'     ## --- Demean & rescale ---
#'     for (j in 1:mdim) {
#'       col_mean <- mean(mom_t[, j])
#'       col_sd <- sd(mom_t[, j])
#'       if (col_sd > 0) {
#'         cols <- ((j - 1) * N * (K+1) + 1):(j * N * (K+1))
#'         mom_micro2[, cols] <- (mom_micro2[, cols] - col_mean) / col_sd
#'         mom_t[, j] <- (mom_t[, j] - col_mean) / col_sd
#'       }
#'     }
#'
#'
#'     ## --- Variance / noise on rescaled data ---
#'     variance <- sum(sapply(1:mdim, function(i) {
#'       norm(mom_micro2[, ((i - 1) * N * (K+1) + 1):(i * N * (K+1))] - mom_t[,i], type = "F")^2
#'     })) / (T * N^2 * (K+1)^2)
#'
#'     data <- mom_t
#'     dim_size <- T
#'
#'   } else if(type == "long_T_moment") {
#'     mom_i <- c()
#'     ## --- Unit-side moments ---
#'     X_av <- matrix(0, nrow = N, ncol = T)  # N x dim_moment
#'
#'
#'     for (t in 1:T) {
#'       # Collect all covariates at time t
#'       X_t <- sapply(X_list, function(X) X[, t])  # N x K
#'       # Take power p **before** averaging over K
#'       X_av[, t] <- rowMeans(X_t)             # (1/K) sum_k x_itk^p
#'     }
#'
#'
#'     mom_i <- cbind(mom_i, X_av)
#'
#'     # mom_micro <- c()
#'     mom_micro <- c()
#'     X_cbind <- do.call(cbind, X_list)
#'     mom_micro <- cbind(mom_micro, X_cbind)
#'
#'
#'     ## --- Demean & rescale ---
#'     for (j in 1:T) {
#'       col_mean <- mean(mom_i[, j])
#'       col_sd <- sd(mom_i[, j])
#'       if (col_sd > 0) {
#'         cols <- ((j - 1)  * (K+1) + 1):(j * (K+1))
#'         mom_micro[, cols] <- (mom_micro[, cols] - col_mean) / col_sd
#'         mom_i[, j] <- (mom_i[, j] - col_mean) / col_sd
#'       }
#'     }
#'
#'
#'     # --- Variance / noise on rescaled data ---
#'     variance <-  sum(sapply(1:T, function(j) {
#'       norm(mom_micro[, (((j - 1) * (K+1) + 1):(j  * (K+1)))] - mom_i[,j], type = "F")^2
#'     })) / (N * (K+1)^2)
#'
#'     data <- mom_i
#'     dim_size <- N
#'
#'   }else if (type == "long multiple moments") {
#'     mom_i <- c()
#'     ## --- Unit-side moments ---
#'     X_av <- matrix(0, nrow = N, ncol = dim_moment)  # N x dim_moment
#'
#'     X_mat <- function(X) {
#'       if (is.null(dim(X))) {
#'         matrix(X, nrow = N, ncol = T)
#'       } else {
#'         X
#'       }
#'     }
#'
#'     for (p in 1:dim_moment) {
#'       # For each power
#'       X_pow <- matrix(0, nrow = N, ncol = T)  # N x T
#'
#'       for (t in 1:T) {
#'         # Collect all covariates at time t
#'         X_t <- do.call(cbind, lapply(X_list, function(X) X_mat(X)[, t]))  # N x (K+1)
#'         # Average the covariates in the p-th group
#'         if ( p == 1){
#'           X_pow[, t] <- rowMeans(X_t)             # (1/K) sum_k x_itk^p
#'         }else if (p == 2){
#'           X_pow[, t] <- rowMeans(tanh(X_t))
#'         }else if (p == 3){
#'           X_pow[, t] <- rowMeans((X_t/(1+abs(X_t))))
#'         }else{
#'           X_pow[, t] <- rowMeans(X_t^p)
#'         }
#'       }
#'
#'       # Average over time
#'       X_av[, p] <- rowMeans(X_pow)               # (1/T) sum_t (1/K) sum_k x_itk^p
#'     }
#'
#'     mom_i <- cbind(mom_i, X_av)
#'
#'     mom_micro <- c()
#'     X_cbind <- do.call(cbind, lapply(X_list, function(X) X_mat(X)))
#'     mom_micro <- cbind(mom_micro, X_cbind)
#'     for (p in 1:dim_moment) {
#'       if ( p == 1){
#'         mom_micro <- cbind(mom_micro, X_cbind^p)             # (1/K) sum_k x_itk^p
#'       }else if (p == 2){
#'         mom_micro <- cbind(mom_micro, tanh(X_cbind))
#'       }else if (p == 3){
#'         mom_micro <- cbind(mom_micro, (X_cbind/(1+abs(X_cbind))))
#'       }else{
#'         mom_micro <- cbind(mom_micro, X_cbind^p)
#'       }
#'     }
#'
#'     ## --- Demean & rescale ---
#'     for (j in 1:mdim) {
#'       col_mean <- mean(mom_i[, j])
#'       col_sd <- sd(mom_i[, j])
#'       if (col_sd > 0) {
#'         cols <- ((j - 1) * T * (K+1) + 1):(j * T * (K+1))
#'         mom_micro[, cols] <- (mom_micro[, cols] - col_mean) / col_sd
#'         mom_i[, j] <- (mom_i[, j] - col_mean) / col_sd
#'       }
#'     }
#'
#'
#'     # --- Variance / noise on rescaled data ---
#'     variance <-  sum(sapply(1:mdim, function(j) {
#'       norm(mom_micro[, ((j - 1) * T * (K+1) + 1):(j * T * (K+1))] - mom_i[,j], type = "F")^2
#'     })) / (N * T^2 * (K+1)^2)
#'
#'     data <- mom_i
#'     dim_size <- N
#'
#'   }else {
#'     stop("Invalid type. Use 'long' or 'tall'.")
#'   }
#'
#'   ## --- Clustering ---
#'   if (!is.null(groups)) {
#'     clusters <- groups
#'     k_result <- kmeans(data, centers = clusters, algorithm = "Lloyd", nstart = init, iter.max = 100)
#'   } else {
#'     if (max(apply(data, 2, sd)) > 0) {
#'       clusters <- 1
#'       xx = 1000
#'       while (xx >= gamma * variance ) {
#'         k_result <- kmeans(data, centers = clusters, algorithm = "Lloyd", nstart = init, iter.max = 100)
#'         xx <- k_result$tot.withinss / dim_size
#'         clusters <- clusters + 1
#'       }
#'       clusters = min(clusters, dim_size-1)
#'       k_result <- kmeans(data,centers = clusters,algorithm = "Lloyd",nstart = init,iter.max = 100)
#'
#'     }else{
#'       k_result <- kmeans(data, centers = 1, algorithm = "Lloyd", nstart = init, iter.max = 100)
#'     }
#'   }
#'
#'   ## --- Return ---
#'   if (type == "long") {
#'     list(res = k_result, clusters = clusters, data = data)
#'   } else {
#'     list(res = k_result, clusters = clusters, data = data)
#'   }
#' }
#'
#' fix_time_clusters <- function(cluster, data, m = 2) {
#'   # cluster: raw K-means labels (length T)
#'   # data: T x 1 matrix (average per time period)
#'   # m: minimum cluster size
#'   #
#'   # Returns:
#'   #   cluster: cleaned labels for each time period
#'   #   G_time: number of clusters
#'
#'   cluster <- as.integer(cluster)
#'   T <- nrow(data)
#'
#'   # Collapse to 1 cluster if impossible
#'   if (T < 2 * m) {
#'     return(list(
#'       cluster = rep(1, T),
#'       G_time = 1
#'     ))
#'   }
#'
#'   repeat {
#'     sizes <- table(cluster)
#'
#'     small <- as.integer(names(sizes[sizes < m]))
#'     large <- as.integer(names(sizes[sizes >= m]))
#'
#'     if (length(small) == 0) break
#'
#'     # Compute centers (safe for 1-column data)
#'     unique_cl <- sort(unique(cluster))
#'     centers <- do.call(rbind, lapply(unique_cl, function(cl) {
#'       colMeans(data[cluster == cl, , drop = FALSE])
#'     }))
#'     rownames(centers) <- unique_cl
#'
#'     # Reassign points from small clusters
#'     for (sc in small) {
#'       idx <- which(cluster == sc)
#'       for (i in idx) {
#'         target <- if (length(large) > 0) large else setdiff(unique_cl, sc)
#'         dists <- sapply(target, function(cl) {
#'           sum((data[i, ] - centers[as.character(cl), ])^2)
#'         })
#'         cluster[i] <- target[which.min(dists)]
#'       }
#'     }
#'   }
#'
#'   # Relabel clusters 1:K
#'   unique_cl <- sort(unique(cluster))
#'   map <- setNames(seq_along(unique_cl), unique_cl)
#'   cluster <- map[as.character(cluster)]
#'
#'   # Number of clusters
#'   G_time <- length(unique_cl)
#'
#'   return(list(
#'     cluster = as.integer(cluster),
#'     G_time = G_time
#'   ))
#' }
#'
#'
#'
#'
#'
#' #' #' Estimator Function with Optional Cross-Fitting
#' #' #'@title Inference in high-dimensional data after discretizing unobserved heterogeneity
#' #' #'
#' #' #' @description R package 'HDpcluster' is dedicated to do inference for high-dimensional linear panel data model with unkown functions of fixed effects.
#' #' #'
#' #' #' @param y outcome variable
#' #' #' @param D treatment variable
#' #' #' @param X control variables
#' #' #' @param T panel size of the time period
#' #' #' @param groups_init number of the unit clusters
#' #' #' @param index index name
#' #' #' @param data data which contains the correct index of the outcome variable or the treatment variable
#' #' #' @param cluster_type different cluster means for the data, 'unit kmeans' is default, allow for 'unit pesudo'
#' #' #' @param pesudo_type if cluster_type = 'unit pesudo', choice of the pesudo type
#' #' #' @param link if cluster_type = 'unit pesudo', choice of different links of pesudo type
#' #' #' @param optimal_index if cluster_type = 'unit pesudo', different ways to compute the optimal number of clusters
#' #' #'
#' #' #' @returns A list of fitted results is returned.
#' #' #' Within this outputted list, the following elements can be found:
#' #' #'     \item{res}{regression model.}
#' #' #'     \item{G}{number of unit clusters.}
#' #' #'     \item{estimate_correct}{corrected standard error, t value, and p value.}
#' #' #'     \item{summary_table}{summary of the model together with corrected estimates in the 'Coefficients'.}
#' #' #'
#' #' #' @import hdm
#' #' #' @import plm
#' #' #' @useDynLib HDpcluster
#' #' #' @export
#' #' HP_estimate <- function(data, y_col = NULL, covariate_cols = NULL,
#' #'                         id_col = "id", time_col = "time", cluster_type = c('two way'),
#' #'                         index = c("id", "time"), method_auto = "dynamicTreeCut", cluster_method = c('kmeans', 'hierarchical'), pre_cluster = FALSE, Du_pre, Dv_pre, unit_cluster = NULL, time_cluster = NULL, gamma = 1, cc = 0, dim_moment = 1) {
#' #'
#' #'
#' #'   if (is.null(covariate_cols)) {
#' #'     covariate_cols <- setdiff(colnames(data), c(id_col, time_col, y_col))
#' #'   }
#' #'
#' #'   ids <- sort(unique(data[[id_col]]))
#' #'   times <- sort(unique(data[[time_col]]))
#' #'   N <- length(ids)
#' #'   T <- length(times)
#' #'   K <- length(covariate_cols)
#' #'
#' #'   # Initialize y and X
#' #'   y <- matrix(NA_real_, nrow = N, ncol = T)
#' #'   X <- array(NA_real_, dim = c(N, T, K))
#' #'
#' #'   for (i in seq_len(nrow(data))) {
#' #'     id_idx <- which(ids == data[[id_col]][i])
#' #'     time_idx <- which(times == data[[time_col]][i])
#' #'     y[id_idx, time_idx] <- data[[y_col]][i]
#' #'     for (k in seq_along(covariate_cols)) {
#' #'       X[id_idx, time_idx, k] <- data[[covariate_cols[k]]][i]
#' #'     }
#' #'   }
#' #'
#' #'   # Cluster (output indicators assumed one-hot encoded)
#' #'   X_norm <- array(NA_real_, dim = dim(X))  # Same shape: (N, T, K)
#' #'   for (k in 1:dim(X)[3]) {
#' #'     X_slice <- X[,,k]
#' #'     X_norm[,,k] <- scale(X_slice)
#' #'   }
#' #'   y_norm = scale(y)
#' #'
#' #'
#' #'   if (cluster_type == 'two way'){
#' #'   if (cluster_method == 'kmeans'){
#' #'
#' #'     # kmeans
#' #'     X_list = lapply(seq_len(dim(X)[3]), function(k) X[ , , k])
#' #'     if (pre_cluster == FALSE){
#' #'       # cluster
#' #'       if (is.null(unit_cluster) == 1 & is.null(time_cluster) == 1){
#' #'         clusteri <- cluster_general(y, X_list, N, T, init = 50, type = "long", gamma = gamma, cc = cc, dim_moment = dim_moment  )
#' #'         G_unit <- clusteri$clusters
#' #'         klong <- clusteri$res
#' #'
#' #'
#' #'         clustert <- cluster_general(y, X_list, N, T, init = 50, type = "tall", gamma = gamma, cc = cc, dim_moment = dim_moment  )
#' #'         G_time <- clustert$clusters
#' #'         ktall <- clustert$res
#' #'         newcluster_T = fix_time_clusters(cluster = ktall$cluster, data = clustert$data, m = 3)
#' #'         ktall$cluster =   newcluster_T$cluster
#' #'         G_time =  newcluster_T$G_time
#' #'       #   if (G_time > floor(T/2)){
#' #'       #     clustert <- cluster_general(y, X_list, N, T, init = 50, type = "tall", groups = 1 )
#' #'       #     G_time <- clustert$clusters
#' #'       #     ktall <- clustert$res
#' #'       #   }
#' #'       }else{
#' #'         clusteri <- cluster_general(y, X_list, N, T, init = 50, type = "long" , groups = c(floor(unit_cluster)) )
#' #'         G_unit <- clusteri$clusters
#' #'         klong <- clusteri$res
#' #'
#' #'
#' #'         clustert <- cluster_general(y, X_list, N, T, init = 50, type = "tall" , groups = c(floor(time_cluster)))
#' #'         G_time <- clustert$clusters
#' #'         ktall <- clustert$res
#' #'       }
#' #'
#' #'
#' #'       Du <- matrix(0, N, G_unit)
#' #'       Dv <- matrix(0, T, G_time)
#' #'       Diagu<- matrix(0, G_unit, G_unit)
#' #'       Diagv<- matrix(0, G_time, G_time)
#' #'
#' #'       for (j in seq_len(G_unit)) {
#' #'         Du[, j] <- as.numeric(klong$cluster == j)
#' #'         Diagu[j,j] <- 1/sum(Du[,j])
#' #'       }
#' #'
#' #'       for (j in seq_len(G_time)) {
#' #'         Dv[, j] <- as.numeric(ktall$cluster == j)
#' #'         Diagv[j,j] <- 1/sum(Dv[,j])
#' #'       }
#' #'     }else if (pre_cluster == TRUE){
#' #'       Du = Du_pre
#' #'       Dv = Dv_pre
#' #'       Diagu<- matrix(0, G_unit, G_unit)
#' #'       Diagv<- matrix(0, G_time, G_time)
#' #'
#' #'       for (j in seq_len(G_unit)) {
#' #'         Diagu[j,j] <- 1/sum(Du[,j])
#' #'       }
#' #'
#' #'       for (j in seq_len(G_time)) {
#' #'         Diagv[j,j] <- 1/sum(Dv[,j])
#' #'       }
#' #'       G_unit = dim(Du)[2]
#' #'       G_time = dim(Dv)[2]
#' #'     }
#' #'   }else if (cluster_method == 'hierarchical'){
#' #'
#' #'     # heriachical
#' #'     if (is.null(unit_cluster) == 1 & is.null(time_cluster) == 1){
#' #'       res_unit <- cluster_Hierarchical(y_norm, X_norm, type = "unit", method_auto = method_auto)
#' #'       res_time <- cluster_Hierarchical(y_norm, X_norm, type = "time", method_auto = method_auto)
#' #'     }else{
#' #'       res_unit <- cluster_Hierarchical(y_norm, X_norm, type = "unit", cluster = unit_cluster)
#' #'       res_time <- cluster_Hierarchical(y_norm, X_norm, type = "time", cluster = time_cluster)
#' #'     }
#' #'
#' #'     Du <- res_unit$indicator    # N x G_unit
#' #'     Dv <- res_time$indicator    # T x G_time
#' #'     G_unit = res_unit$G
#' #'     G_time = res_time$G
#' #'     Diagu<- matrix(0, G_unit, G_unit)
#' #'     Diagv<- matrix(0, G_time, G_time)
#' #'
#' #'     for (j in seq_len(G_unit)) {
#' #'       Diagu[j,j] <- 1/sum(Du[,j])
#' #'     }
#' #'
#' #'     for (j in seq_len(G_time)) {
#' #'       Diagv[j,j] <- 1/sum(Dv[,j])
#' #'     }
#' #'   }
#' #'
#' #'   unit_group <- max.col(Du)
#' #'   time_group <- max.col(Dv)
#' #'
#' #'   # Projection matrices
#' #'   Mu <- diag(N) - Du %*% Diagu %*% t(Du)
#' #'   Mv <- diag(T) - Dv %*% Diagv %*% t(Dv)
#' #'   # Mu <- diag(N) - Du %*% solve(t(Du) %*% Du) %*% t(Du)
#' #'   # Mv <- diag(T) - Dv %*% solve(t(Dv) %*% Dv) %*% t(Dv)
#' #'
#' #'   # Stack y and X into Z: N x T x (K+1)
#' #'   Z <- array(NA_real_, dim = c(N, T, K + 1))
#' #'   Z[,,1] <- y
#' #'   for (k in 1:K) Z[,,k+1] <- X[,,k]
#' #'
#' #'   # Projection method: demean unit and time clusters slice-wise
#' #'   Z_proj <- array(NA_real_, dim = c(N, T, K + 1))
#' #'   for (k in 1:(K+1)) {
#' #'     Z_proj[,,k] <- Mu %*% Z[,,k] %*% Mv
#' #'   }
#' #'
#' #'   # Separate transformed y and X
#' #'   y_trans <- Z_proj[,,1] # Z_proj_final_array[,,1]
#' #'   X_trans <- Z_proj[,,2:(K+1)] # Z_proj_final_array[,,2:(K+1)]
#' #'
#' #'   tY_vector <- as.vector(y_trans)  # (N*T)
#' #'   tX_matrix <- matrix(aperm(X_trans, c(1, 2, 3)), nrow = N*T, ncol = K)
#' #'
#' #'   d_vector <- tX_matrix[, 1]
#' #'   x_matrix <- tX_matrix[, -1, drop = FALSE]
#' #'
#' #'   # Run rlassoEffect with correct inputs
#' #'   fit <- rlassoEffect(x = x_matrix, y = tY_vector, d = d_vector, method = "double selection")
#' #'
#' #'   trans <- data.frame(y = tY_vector,
#' #'                       D = d_vector,
#' #'                       x_matrix)
#' #'
#' #'   lasso.Y <- rlasso(y ~ . - D - 1, data = trans)
#' #'   Ytilde <- lasso.Y$residuals
#' #'
#' #'   lasso.D <- rlasso(D ~ . - 1, data = trans[, -1])
#' #'   Dtilde <- lasso.D$residuals
#' #'
#' #'   data_res <- data.frame(id = data[[index[1]]],
#' #'                          time = data[[index[2]]],
#' #'                          Ytilde = Ytilde,
#' #'                          Dtilde = Dtilde)
#' #'
#' #'   Post_plm <- plm(Ytilde ~ -1 + Dtilde, data = data_res, model = "pooling", index = c("id", "time"))
#' #'
#' #'   coefs <- coef(Post_plm)
#' #'   if (N*T - N*G_time - T*G_unit > 0){
#' #'   se_corrected <- sqrt(vcovHC(Post_plm, type = "HC0", method = "arellano")) * sqrt(N * T / (N*T - N*G_time - T*G_unit))
#' #'   }else{
#' #'     se_corrected <- sqrt(vcovHC(Post_plm, type = "HC0", method = "arellano")) * sqrt(N * T / (1))
#' #'   }
#' #'   t_values_corrected <- coefs / se_corrected
#' #'
#' #'   df <- Post_plm$df.residual
#' #'   p_values_corrected <- 2 * pt(-abs(t_values_corrected), df)
#' #'
#' #'   summary_table_correct <- data.frame(
#' #'     Estimate = coefs,
#' #'     `Std. Error corrected` = se_corrected,
#' #'     `t-value corrected` = t_values_corrected,
#' #'     `Pr(>|t|) corrected` = p_values_corrected
#' #'   )
#' #'
#' #'   colnames(summary_table_correct) = c('Estimate', 'Std. Error corrected', 't-value corrected', 'Pr(>|t|) corrected')
#' #'
#' #'   # Return cluster counts as well
#' #'   return(list(
#' #'     fit_summary = summary(fit),
#' #'     G_unit = G_unit,
#' #'     G_time = G_time,
#' #'     unit_group = unit_group,
#' #'     time_group = time_group,
#' #'     post_plm_summary = summary(Post_plm),
#' #'     estimate_corrected = summary_table_correct,
#' #'     summary_table = summary_table_correct
#' #'   ))
#' #'   }else if (cluster_type == 'one way'){
#' #'     if (cluster_method == 'kmeans'){
#' #'
#' #'       # kmeans
#' #'       X_list = lapply(seq_len(dim(X)[3]), function(k) X[ , , k])
#' #'       if (pre_cluster == FALSE){
#' #'         # cluster
#' #'         if (is.null(unit_cluster) == 1 & is.null(time_cluster) == 1){
#' #'           clusteri <- cluster_general(y, X_list, N, T, init = 50, type = "long", gamma = gamma, cc = cc, dim_moment = dim_moment  )
#' #'           G_unit <- clusteri$clusters
#' #'           klong <- clusteri$res
#' #'
#' #'
#' #'         }else{
#' #'           clusteri <- cluster_general(y, X_list, N, T, init = 50, type = "long" , groups = c(floor(unit_cluster)) )
#' #'           G_unit <- clusteri$clusters
#' #'           klong <- clusteri$res
#' #'         }
#' #'         Du <- matrix(0, N, G_unit)
#' #'         Diagu<- matrix(0, G_unit, G_unit)
#' #'
#' #'
#' #'         for (j in seq_len(G_unit)) {
#' #'           Du[, j] <- as.numeric(klong$cluster == j)
#' #'           Diagu[j,j] <- 1/sum(Du[,j])
#' #'         }
#' #'
#' #'
#' #'       }else if (pre_cluster == TRUE){
#' #'         Du = Du_pre
#' #'         Diagu<- matrix(0, G_unit, G_unit)
#' #'
#' #'         for (j in seq_len(G_unit)) {
#' #'           Diagu[j,j] <- 1/sum(Du[,j])
#' #'         }
#' #'
#' #'         G_unit = dim(Du)[2]
#' #'       }
#' #'     }else if (cluster_method == 'hierarchical'){
#' #'
#' #'       # heriachical
#' #'       if (is.null(unit_cluster) == 1 & is.null(time_cluster) == 1){
#' #'         res_unit <- cluster_Hierarchical(y_norm, X_norm, type = "unit", method_auto = method_auto)
#' #'       }else{
#' #'         res_unit <- cluster_Hierarchical(y_norm, X_norm, type = "unit", cluster = unit_cluster)
#' #'       }
#' #'
#' #'       Du <- res_unit$indicator    # N x G_unit
#' #'       G_unit = res_unit$G
#' #'       Diagu<- matrix(0, G_unit, G_unit)
#' #'
#' #'       for (j in seq_len(G_unit)) {
#' #'         Diagu[j,j] <- 1/sum(Du[,j])
#' #'       }
#' #'
#' #'     }
#' #'
#' #'     unit_group <- max.col(Du)
#' #'
#' #'     # Projection matrices
#' #'     Mu <- diag(N) - Du %*% Diagu %*% t(Du)
#' #'
#' #'     # Stack y and X into Z: N x T x (K+1)
#' #'     Z <- array(NA_real_, dim = c(N, T, K + 1))
#' #'     Z[,,1] <- y
#' #'     for (k in 1:K) Z[,,k+1] <- X[,,k]
#' #'
#' #'     # Projection method: demean unit and time clusters slice-wise
#' #'     Z_proj <- array(NA_real_, dim = c(N, T, K + 1))
#' #'     for (k in 1:(K+1)) {
#' #'       Z_proj[,,k] <- Mu %*% Z[,,k]
#' #'     }
#' #'
#' #'     # Separate transformed y and X
#' #'     y_trans <- Z_proj[,,1] # Z_proj_final_array[,,1]
#' #'     X_trans <- Z_proj[,,2:(K+1)] # Z_proj_final_array[,,2:(K+1)]
#' #'
#' #'     tY_vector <- as.vector(y_trans)  # (N*T)
#' #'     tX_matrix <- matrix(aperm(X_trans, c(1, 2, 3)), nrow = N*T, ncol = K)
#' #'
#' #'     # debiased lasso
#' #'     # debiased_lasso = lasso.proj(tX_matrix, tY_vector) # delas(tX_matrix, tY_vector, c(1))
#' #'     # debiased_lasso$se = debiased_lasso$se * sqrt(N * T / (N*T - T*G_unit))
#' #'
#' #'     #
#' #'     d_vector <- tX_matrix[, 1]
#' #'     x_matrix <- tX_matrix[, -1, drop = FALSE]
#' #'
#' #'     # Run rlassoEffect with correct inputs
#' #'     fit <- rlassoEffect(x = x_matrix, y = tY_vector, d = d_vector, method = "double selection")
#' #'
#' #'     trans <- data.frame(y = tY_vector,
#' #'                         D = d_vector,
#' #'                         x_matrix)
#' #'
#' #'     lasso.Y <- rlasso(y ~ . - D - 1, data = trans)
#' #'     Ytilde <- lasso.Y$residuals
#' #'
#' #'     lasso.D <- rlasso(D ~ . - 1, data = trans[, -1])
#' #'     Dtilde <- lasso.D$residuals
#' #'
#' #'     data_res <- data.frame(id = data[[index[1]]],
#' #'                            time = data[[index[2]]],
#' #'                            Ytilde = Ytilde,
#' #'                            Dtilde = Dtilde)
#' #'
#' #'     Post_plm <- plm(Ytilde ~ -1 + Dtilde, data = data_res, model = "pooling", index = c("id", "time"))
#' #'
#' #'     coefs <- coef(Post_plm)
#' #'     se_corrected <- sqrt(vcovHC(Post_plm, type = "HC0", method = "arellano")) * sqrt(N * T / (N*T - T*G_unit))
#' #'     t_values_corrected <- coefs / se_corrected
#' #'
#' #'     df <- Post_plm$df.residual
#' #'     p_values_corrected <- 2 * pt(-abs(t_values_corrected), df)
#' #'
#' #'     summary_table_correct <- data.frame(
#' #'       Estimate = coefs,
#' #'       `Std. Error corrected` = se_corrected,
#' #'       `t-value corrected` = t_values_corrected,
#' #'       `Pr(>|t|) corrected` = p_values_corrected
#' #'     )
#' #'
#' #'     colnames(summary_table_correct) = c('Estimate', 'Std. Error corrected', 't-value corrected', 'Pr(>|t|) corrected')
#' #'
#' #'     # Return cluster counts as well
#' #'     return(list(
#' #'       fit_summary = summary(fit),
#' #'       G_unit = G_unit,
#' #'       unit_group = unit_group,
#' #'       post_plm_summary = summary(Post_plm),
#' #'       estimate_corrected = summary_table_correct,
#' #'       summary_table = summary_table_correct
#' #'     ))
#' #'
#' #'   }else if (cluster_type == 'one way T moment'){
#' #'     if (cluster_method == 'kmeans'){
#' #'
#' #'       # kmeans
#' #'       X_list = lapply(seq_len(dim(X)[3]), function(k) X[ , , k])
#' #'       if (pre_cluster == FALSE){
#' #'         # cluster
#' #'         if (is.null(unit_cluster) == 1 & is.null(time_cluster) == 1){
#' #'           clusteri <- cluster_general(y, X_list, N, T, init = 50, type = "long_T_moment", gamma = gamma, cc = cc, dim_moment = dim_moment  )
#' #'           G_unit <- clusteri$clusters
#' #'           klong <- clusteri$res
#' #'
#' #'
#' #'         }else{
#' #'           clusteri <- cluster_general(y, X_list, N, T, init = 50, type = "long_T_moment" , groups = c(floor(unit_cluster)) )
#' #'           G_unit <- clusteri$clusters
#' #'           klong <- clusteri$res
#' #'         }
#' #'         Du <- matrix(0, N, G_unit)
#' #'         Diagu<- matrix(0, G_unit, G_unit)
#' #'
#' #'
#' #'         for (j in seq_len(G_unit)) {
#' #'           Du[, j] <- as.numeric(klong$cluster == j)
#' #'           Diagu[j,j] <- 1/sum(Du[,j])
#' #'         }
#' #'
#' #'
#' #'       }else if (pre_cluster == TRUE){
#' #'         Du = Du_pre
#' #'         Diagu<- matrix(0, G_unit, G_unit)
#' #'
#' #'         for (j in seq_len(G_unit)) {
#' #'           Diagu[j,j] <- 1/sum(Du[,j])
#' #'         }
#' #'
#' #'         G_unit = dim(Du)[2]
#' #'       }
#' #'     }else if (cluster_method == 'hierarchical'){
#' #'
#' #'       # heriachical
#' #'       if (is.null(unit_cluster) == 1 & is.null(time_cluster) == 1){
#' #'         res_unit <- cluster_Hierarchical(y_norm, X_norm, type = "unit", method_auto = method_auto)
#' #'       }else{
#' #'         res_unit <- cluster_Hierarchical(y_norm, X_norm, type = "unit", cluster = unit_cluster)
#' #'       }
#' #'
#' #'       Du <- res_unit$indicator    # N x G_unit
#' #'       G_unit = res_unit$G
#' #'       Diagu<- matrix(0, G_unit, G_unit)
#' #'
#' #'       for (j in seq_len(G_unit)) {
#' #'         Diagu[j,j] <- 1/sum(Du[,j])
#' #'       }
#' #'
#' #'     }
#' #'
#' #'     unit_group <- max.col(Du)
#' #'
#' #'     # Projection matrices
#' #'     Mu <- diag(N) - Du %*% Diagu %*% t(Du)
#' #'
#' #'     # Stack y and X into Z: N x T x (K+1)
#' #'     Z <- array(NA_real_, dim = c(N, T, K + 1))
#' #'     Z[,,1] <- y
#' #'     for (k in 1:K) Z[,,k+1] <- X[,,k]
#' #'
#' #'     # Projection method: demean unit and time clusters slice-wise
#' #'     Z_proj <- array(NA_real_, dim = c(N, T, K + 1))
#' #'     for (k in 1:(K+1)) {
#' #'       Z_proj[,,k] <- Mu %*% Z[,,k]
#' #'     }
#' #'
#' #'     # Separate transformed y and X
#' #'     y_trans <- Z_proj[,,1] # Z_proj_final_array[,,1]
#' #'     X_trans <- Z_proj[,,2:(K+1)] # Z_proj_final_array[,,2:(K+1)]
#' #'
#' #'     tY_vector <- as.vector(y_trans)  # (N*T)
#' #'     tX_matrix <- matrix(aperm(X_trans, c(1, 2, 3)), nrow = N*T, ncol = K)
#' #'
#' #'     # debiased lasso
#' #'     # debiased_lasso = lasso.proj(tX_matrix, tY_vector) # delas(tX_matrix, tY_vector, c(1))
#' #'     # debiased_lasso$se = debiased_lasso$se * sqrt(N * T / (N*T - T*G_unit))
#' #'
#' #'     #
#' #'     d_vector <- tX_matrix[, 1]
#' #'     x_matrix <- tX_matrix[, -1, drop = FALSE]
#' #'
#' #'     # Run rlassoEffect with correct inputs
#' #'     fit <- rlassoEffect(x = x_matrix, y = tY_vector, d = d_vector, method = "double selection")
#' #'
#' #'     trans <- data.frame(y = tY_vector,
#' #'                         D = d_vector,
#' #'                         x_matrix)
#' #'
#' #'     lasso.Y <- rlasso(y ~ . - D - 1, data = trans)
#' #'     Ytilde <- lasso.Y$residuals
#' #'
#' #'     lasso.D <- rlasso(D ~ . - 1, data = trans[, -1])
#' #'     Dtilde <- lasso.D$residuals
#' #'
#' #'     data_res <- data.frame(id = data[[index[1]]],
#' #'                            time = data[[index[2]]],
#' #'                            Ytilde = Ytilde,
#' #'                            Dtilde = Dtilde)
#' #'
#' #'     Post_plm <- plm(Ytilde ~ -1 + Dtilde, data = data_res, model = "pooling", index = c("id", "time"))
#' #'
#' #'     coefs <- coef(Post_plm)
#' #'     se_corrected <- sqrt(vcovHC(Post_plm, type = "HC0", method = "arellano")) * sqrt(N * T / (N*T - T*G_unit))
#' #'     t_values_corrected <- coefs / se_corrected
#' #'
#' #'     df <- Post_plm$df.residual
#' #'     p_values_corrected <- 2 * pt(-abs(t_values_corrected), df)
#' #'
#' #'     summary_table_correct <- data.frame(
#' #'       Estimate = coefs,
#' #'       `Std. Error corrected` = se_corrected,
#' #'       `t-value corrected` = t_values_corrected,
#' #'       `Pr(>|t|) corrected` = p_values_corrected
#' #'     )
#' #'
#' #'     colnames(summary_table_correct) = c('Estimate', 'Std. Error corrected', 't-value corrected', 'Pr(>|t|) corrected')
#' #'
#' #'     # Return cluster counts as well
#' #'     return(list(
#' #'       fit_summary = summary(fit),
#' #'       G_unit = G_unit,
#' #'       unit_group = unit_group,
#' #'       post_plm_summary = summary(Post_plm),
#' #'       estimate_corrected = summary_table_correct,
#' #'       summary_table = summary_table_correct
#' #'     ))
#' #'
#' #'   }
#' #' }
#' #'
#' #'
#' #'
#' #' # HP_estimate_double <- function(data, y_col = NULL, covariate_cols = NULL,
#' #' #                                id_col = "id", time_col = "time",
#' #' #                                index = c("id", "time"), method_auto = "dynamicTreeCut", cluster_method = 'kmeans', pre_cluster = FALSE, Du_pre, Dv_pre, unit_cluster = NULL, time_cluster = NULL, gamma = 1, cc = 0, dim_moment = 1) {
#' #' #
#' #' #
#' #' #   if (is.null(covariate_cols)) {
#' #' #     covariate_cols <- setdiff(colnames(data), c(id_col, time_col, y_col))
#' #' #   }
#' #' #
#' #' #   ids <- sort(unique(data[[id_col]]))
#' #' #   times <- sort(unique(data[[time_col]]))
#' #' #   N <- length(ids)
#' #' #   T <- length(times)
#' #' #   K <- length(covariate_cols)
#' #' #
#' #' #   # Initialize y and X
#' #' #   y <- matrix(NA_real_, nrow = N, ncol = T)
#' #' #   X <- array(NA_real_, dim = c(N, T, K))
#' #' #
#' #' #   for (i in seq_len(nrow(data))) {
#' #' #     id_idx <- which(ids == data[[id_col]][i])
#' #' #     time_idx <- which(times == data[[time_col]][i])
#' #' #     y[id_idx, time_idx] <- data[[y_col]][i]
#' #' #     for (k in seq_along(covariate_cols)) {
#' #' #       X[id_idx, time_idx, k] <- data[[covariate_cols[k]]][i]
#' #' #     }
#' #' #   }
#' #' #
#' #' #   # Cluster (output indicators assumed one-hot encoded)
#' #' #   X_norm <- array(NA_real_, dim = dim(X))  # Same shape: (N, T, K)
#' #' #   for (k in 1:dim(X)[3]) {
#' #' #     X_slice <- X[,,k]
#' #' #     X_norm[,,k] <- scale(X_slice)
#' #' #   }
#' #' #   y_norm = scale(y)
#' #' #   if (cluster_method == 'kmeans'){
#' #' #
#' #' #     # kmeans
#' #' #     X_list = lapply(seq_len(dim(X)[3]), function(k) X[ , , k])
#' #' #     if (pre_cluster == FALSE){
#' #' #       # cluster
#' #' #       if (is.null(unit_cluster) == 1 & is.null(time_cluster) == 1){
#' #' #         clusteri <- cluster_general(y, X_list, N, T, init = 50, type = "long", gamma = gamma, cc = cc, dim_moment = dim_moment   )
#' #' #         G_unit <- clusteri$clusters
#' #' #         klong <- clusteri$res
#' #' #
#' #' #
#' #' #         clustert <- cluster_general(y, X_list, N, T, init = 50, type = "tall", gamma = gamma, cc = cc, dim_moment = dim_moment  )
#' #' #         G_time <- clustert$clusters
#' #' #         ktall <- clustert$res
#' #' #       }else{
#' #' #         clusteri <- cluster_general(y, X_list, N, T, init = 50, type = "long" , groups = c(floor(unit_cluster)) )
#' #' #         G_unit <- clusteri$clusters
#' #' #         klong <- clusteri$res
#' #' #
#' #' #
#' #' #         clustert <- cluster_general(y, X_list, N, T, init = 50, type = "tall" , groups = c(floor(time_cluster)))
#' #' #         G_time <- clustert$clusters
#' #' #         ktall <- clustert$res
#' #' #       }
#' #' #
#' #' #       Du <- matrix(0, N, G_unit)
#' #' #       Dv <- matrix(0, T, G_time)
#' #' #       Diagu<- matrix(0, G_unit, G_unit)
#' #' #       Diagv<- matrix(0, G_time, G_time)
#' #' #
#' #' #       for (j in seq_len(G_unit)) {
#' #' #         Du[, j] <- as.numeric(klong$cluster == j)
#' #' #         Diagu[j,j] <- 1/sum(Du[,j])
#' #' #       }
#' #' #
#' #' #       for (j in seq_len(G_time)) {
#' #' #         Dv[, j] <- as.numeric(ktall$cluster == j)
#' #' #         Diagv[j,j] <- 1/sum(Dv[,j])
#' #' #       }
#' #' #     }else if (pre_cluster == TRUE){
#' #' #       Du = Du_pre
#' #' #       Dv = Dv_pre
#' #' #       Diagu<- matrix(0, G_unit, G_unit)
#' #' #       Diagv<- matrix(0, G_time, G_time)
#' #' #
#' #' #       for (j in seq_len(G_unit)) {
#' #' #         Diagu[j,j] <- 1/sum(Du[,j])
#' #' #       }
#' #' #
#' #' #       for (j in seq_len(G_time)) {
#' #' #         Diagv[j,j] <- 1/sum(Dv[,j])
#' #' #       }
#' #' #       G_unit = dim(Du)[2]
#' #' #       G_time = dim(Dv)[2]
#' #' #     }
#' #' #   }else if (cluster_method == 'hierarchical'){
#' #' #
#' #' #     # heriachical
#' #' #     if (is.null(unit_cluster) == 1 & is.null(time_cluster) == 1){
#' #' #       res_unit <- cluster_Hierarchical(y_norm, X_norm, type = "unit", method_auto = method_auto)
#' #' #       res_time <- cluster_Hierarchical(y_norm, X_norm, type = "time", method_auto = method_auto)
#' #' #     }else{
#' #' #       res_unit <- cluster_Hierarchical(y_norm, X_norm, type = "unit", cluster = unit_cluster)
#' #' #       res_time <- cluster_Hierarchical(y_norm, X_norm, type = "time", cluster = time_cluster)
#' #' #     }
#' #' #
#' #' #     Du <- res_unit$indicator    # N x G_unit
#' #' #     Dv <- res_time$indicator    # T x G_time
#' #' #
#' #' #     G_unit = res_unit$G
#' #' #     G_time = res_time$G
#' #' #
#' #' #     Diagu<- matrix(0, G_unit, G_unit)
#' #' #     Diagv<- matrix(0, G_time, G_time)
#' #' #
#' #' #     for (j in seq_len(G_unit)) {
#' #' #       Diagu[j,j] <- 1/sum(Du[,j])
#' #' #     }
#' #' #
#' #' #     for (j in seq_len(G_time)) {
#' #' #       Diagv[j,j] <- 1/sum(Dv[,j])
#' #' #     }
#' #' #   }
#' #' #
#' #' #   unit_group <- max.col(Du)
#' #' #   time_group <- max.col(Dv)
#' #' #
#' #' #   # Projection matrices
#' #' #   Mu <- diag(N) - Du %*% Diagu %*% t(Du)
#' #' #   Mv <- diag(T) - Dv %*% Diagv %*% t(Dv)
#' #' #   # Mu <- diag(N) - Du %*% solve(t(Du) %*% Du) %*% t(Du)
#' #' #   # Mv <- diag(T) - Dv %*% solve(t(Dv) %*% Dv) %*% t(Dv)
#' #' #
#' #' #   # Stack y and X into Z: N x T x (K+1)
#' #' #   Z <- array(NA_real_, dim = c(N, T, K + 1))
#' #' #   Z[,,1] <- y
#' #' #   for (k in 1:K) Z[,,k+1] <- X[,,k]
#' #' #
#' #' #   # Projection method: demean unit and time clusters slice-wise
#' #' #   Z_proj <- array(NA_real_, dim = c(N, T, K + 1))
#' #' #   for (k in 1:(K+1)) {
#' #' #     Z_proj[,,k] <- Mu %*% Z[,,k] %*% Mv
#' #' #   }
#' #' #
#' #' #
#' #' #   y_trans <- Z_proj[,,1] # Z_proj_final_array[,,1]
#' #' #   X_trans <- Z_proj[,,2:(K+1)] # Z_proj_final_array[,,2:(K+1)]
#' #' #
#' #' #   tY_vector <- as.vector(y_trans)  # (N*T)
#' #' #   tX_matrix <- matrix(aperm(X_trans, c(1, 2, 3)), nrow = N*T, ncol = K)
#' #' #
#' #' #   d_vector <- tX_matrix[, 1]
#' #' #   x_matrix <- tX_matrix[, -1, drop = FALSE]
#' #' #   # Separate transformed y and X
#' #' #
#' #' #   trans <- data.frame(y = tY_vector,
#' #' #                       D = d_vector,
#' #' #                       x_matrix)
#' #' #
#' #' #   # Run rlassoEffect with correct inputs
#' #' #   fit <- rlassoEffect(x = x_matrix, y = tY_vector, d = d_vector, method = "double selection")
#' #' #
#' #' #   trans <- data.frame(y = tY_vector,
#' #' #                       D = d_vector,
#' #' #                       x_matrix)
#' #' #
#' #' #   lasso.Y <- rlasso(y ~ . - D - 1, data = trans)
#' #' #   Ytilde <- lasso.Y$residuals
#' #' #
#' #' #   lasso.D <- rlasso(D ~ . - 1, data = trans[, -1])
#' #' #   Dtilde <- lasso.D$residuals
#' #' #
#' #' #   data_res <- data.frame(id = data[[index[1]]],
#' #' #                          time = data[[index[2]]],
#' #' #                          Ytilde = Ytilde,
#' #' #                          Dtilde = Dtilde)
#' #' #
#' #' #   Post_plm <- plm(Ytilde ~ -1 + Dtilde, data = data_res, model = "pooling", index = c("id", "time"))
#' #' #
#' #' #   coefs <- coef(Post_plm)
#' #' #   se_corrected <- sqrt(vcovHC(Post_plm, type = "HC0", method = "arellano")) * sqrt(N * T / (N*T - N*G_time - T*G_unit))
#' #' #   t_values_corrected <- coefs / se_corrected
#' #' #
#' #' #   df <- Post_plm$df.residual
#' #' #   p_values_corrected <- 2 * pt(-abs(t_values_corrected), df)
#' #' #
#' #' #   summary_table_correct <- data.frame(
#' #' #     Estimate = coefs,
#' #' #     `Std. Error corrected` = se_corrected,
#' #' #     `t-value corrected` = t_values_corrected,
#' #' #     `Pr(>|t|) corrected` = p_values_corrected
#' #' #   )
#' #' #
#' #' #   colnames(summary_table_correct) = c('Estimate', 'Std. Error corrected', 't-value corrected', 'Pr(>|t|) corrected')
#' #' #
#' #' #   # Return cluster counts as well
#' #' #   return(list(
#' #' #     fit_summary = summary(fit),
#' #' #     G_unit = G_unit,
#' #' #     G_time = G_time,
#' #' #     unit_group = unit_group,
#' #' #     time_group = time_group,
#' #' #     post_plm_summary = summary(Post_plm),
#' #' #     estimate_corrected = summary_table_correct,
#' #' #     summary_table = summary_table_correct
#' #' #   ))
#' #' # }
#' #'
#' #' #' @export
#' #' Niave_estimate <- function(data, y_col = NULL, covariate_cols = NULL,
#' #'                            id_col = "id", time_col = "time", index = c("id", "time")) {
#' #'
#' #'
#' #'   if (is.null(covariate_cols)) {
#' #'     covariate_cols <- setdiff(colnames(data), c(id_col, time_col, y_col))
#' #'   }
#' #'
#' #'   ids <- sort(unique(data[[id_col]]))
#' #'   times <- sort(unique(data[[time_col]]))
#' #'   N <- length(ids)
#' #'   T <- length(times)
#' #'   K <- length(covariate_cols)
#' #'
#' #'   # Initialize y and X
#' #'   y <- matrix(NA_real_, nrow = N, ncol = T)
#' #'   X <- array(NA_real_, dim = c(N, T, K))
#' #'
#' #'   for (i in seq_len(nrow(data))) {
#' #'     id_idx <- which(ids == data[[id_col]][i])
#' #'     time_idx <- which(times == data[[time_col]][i])
#' #'     y[id_idx, time_idx] <- data[[y_col]][i]
#' #'     for (k in seq_along(covariate_cols)) {
#' #'       X[id_idx, time_idx, k] <- data[[covariate_cols[k]]][i]
#' #'     }
#' #'   }
#' #'
#' #'   # Stack y and X into Z: N x T x (K+1)
#' #'   Z <- array(NA_real_, dim = c(N, T, K + 1))
#' #'   Z[,,1] <- y
#' #'   for (k in 1:K) Z[,,k+1] <- X[,,k]
#' #'
#' #'   # Projection method: demean unit and time clusters slice-wise
#' #'   Z_proj <- array(NA_real_, dim = c(N, T, K + 1))
#' #'   for (k in 1:(K+1)) {
#' #'     Z_proj[,,k] <-  Z[,,k]
#' #'   }
#' #'
#' #'   # Separate transformed y and X
#' #'   y_trans <- Z_proj[,,1] # Z_proj_final_array[,,1]
#' #'   X_trans <- Z_proj[,,2:(K+1)] # Z_proj_final_array[,,2:(K+1)]
#' #'
#' #'   tY_vector <- as.vector(y_trans)  # (N*T)
#' #'   tX_matrix <- matrix(aperm(X_trans, c(1, 2, 3)), nrow = N*T, ncol = K)
#' #'
#' #'   d_vector <- tX_matrix[, 1]
#' #'   x_matrix <- tX_matrix[, -1, drop = FALSE]
#' #'
#' #'   # Run rlassoEffect with correct inputs
#' #'   fit <- rlassoEffect(x = x_matrix, y = tY_vector, d = d_vector, method = "double selection")
#' #'
#' #'   trans <- data.frame(y = tY_vector,
#' #'                       D = d_vector,
#' #'                       x_matrix)
#' #'
#' #'   lasso.Y <- rlasso(y ~ . - D - 1, data = trans)
#' #'   Ytilde <- lasso.Y$residuals
#' #'
#' #'   lasso.D <- rlasso(D ~ . - 1, data = trans[, -1])
#' #'   Dtilde <- lasso.D$residuals
#' #'
#' #'   data_res <- data.frame(id = data[[index[1]]],
#' #'                          time = data[[index[2]]],
#' #'                          Ytilde = Ytilde,
#' #'                          Dtilde = Dtilde)
#' #'
#' #'   Post_plm <- plm(Ytilde ~ -1 + Dtilde, data = data_res, model = "pooling", index = c("id", "time"))
#' #'
#' #'   coefs <- coef(Post_plm)
#' #'   se_corrected <- sqrt(vcovHC(Post_plm, type = "HC0", method = "arellano"))
#' #'   t_values_corrected <- coefs / se_corrected
#' #'
#' #'   df <- Post_plm$df.residual
#' #'   p_values_corrected <- 2 * pt(-abs(t_values_corrected), df)
#' #'
#' #'   summary_table_correct <- data.frame(
#' #'     Estimate = coefs,
#' #'     `Std. Error corrected` = se_corrected,
#' #'     `t-value corrected` = t_values_corrected,
#' #'     `Pr(>|t|) corrected` = p_values_corrected
#' #'   )
#' #'
#' #'   colnames(summary_table_correct) = c('Estimate', 'Std. Error corrected', 't-value corrected', 'Pr(>|t|) corrected')
#' #'
#' #'   # Return cluster counts as well
#' #'   return(list(
#' #'     fit_summary = summary(fit),
#' #'     post_plm_summary = summary(Post_plm),
#' #'     estimate_corrected = summary_table_correct,
#' #'     summary_table = summary_table_correct
#' #'   ))
#' #' }
#' #'
#' #' #' @export
#' #' TWFE_estimate <- function(data, y_col = NULL, covariate_cols = NULL,
#' #'                           id_col = "id", time_col = "time", index = c("id", "time")) {
#' #'
#' #'
#' #'   if (is.null(covariate_cols)) {
#' #'     covariate_cols <- setdiff(colnames(data), c(id_col, time_col, y_col))
#' #'   }
#' #'
#' #'   ids <- sort(unique(data[[id_col]]))
#' #'   times <- sort(unique(data[[time_col]]))
#' #'   N <- length(ids)
#' #'   T <- length(times)
#' #'   K <- length(covariate_cols)
#' #'
#' #'   # Initialize y and X
#' #'   y <- matrix(NA_real_, nrow = N, ncol = T)
#' #'   X <- array(NA_real_, dim = c(N, T, K))
#' #'
#' #'   for (i in seq_len(nrow(data))) {
#' #'     id_idx <- which(ids == data[[id_col]][i])
#' #'     time_idx <- which(times == data[[time_col]][i])
#' #'     y[id_idx, time_idx] <- data[[y_col]][i]
#' #'     for (k in seq_along(covariate_cols)) {
#' #'       X[id_idx, time_idx, k] <- data[[covariate_cols[k]]][i]
#' #'     }
#' #'   }
#' #'
#' #'   # Cluster (output indicators assumed one-hot encoded)
#' #'   X_norm <- array(NA_real_, dim = dim(X))  # Same shape: (N, T, K)
#' #'   for (k in 1:dim(X)[3]) {
#' #'     X_slice <- X[,,k]
#' #'     X_norm[,,k] <- scale(X_slice)
#' #'   }
#' #'   y_norm = scale(y)
#' #'
#' #'   Du <- matrix(1, nrow = N, ncol = 1)    # N x G_unit
#' #'   Dv <- matrix(1, nrow = T, ncol = 1)    # T x G_time
#' #'
#' #'   G_unit = 1
#' #'   G_time = 1
#' #'
#' #'   unit_group <- max.col(Du)
#' #'   time_group <- max.col(Dv)
#' #'
#' #'   # Projection matrices
#' #'   Mu <- diag(N) - Du %*% solve(t(Du) %*% Du) %*% t(Du)
#' #'   Mv <- diag(T) - Dv %*% solve(t(Dv) %*% Dv) %*% t(Dv)
#' #'
#' #'   # Stack y and X into Z: N x T x (K+1)
#' #'   Z <- array(NA_real_, dim = c(N, T, K + 1))
#' #'   Z[,,1] <- y
#' #'   for (k in 1:K) Z[,,k+1] <- X[,,k]
#' #'
#' #'   # Projection method: demean unit and time clusters slice-wise
#' #'   Z_proj <- array(NA_real_, dim = c(N, T, K + 1))
#' #'   for (k in 1:(K+1)) {
#' #'     Z_proj[,,k] <- Mu %*% Z[,,k] %*% Mv
#' #'   }
#' #'
#' #'   # Separate transformed y and X
#' #'   y_trans <- Z_proj[,,1] # Z_proj_final_array[,,1]
#' #'   X_trans <- Z_proj[,,2:(K+1)] # Z_proj_final_array[,,2:(K+1)]
#' #'
#' #'   tY_vector <- as.vector(y_trans)  # (N*T)
#' #'   tX_matrix <- matrix(aperm(X_trans, c(1, 2, 3)), nrow = N*T, ncol = K)
#' #'
#' #'   d_vector <- tX_matrix[, 1]
#' #'   x_matrix <- tX_matrix[, -1, drop = FALSE]
#' #'
#' #'   # Run rlassoEffect with correct inputs
#' #'   fit <- rlassoEffect(x = x_matrix, y = tY_vector, d = d_vector, method = "double selection")
#' #'
#' #'   trans <- data.frame(y = tY_vector,
#' #'                       D = d_vector,
#' #'                       x_matrix)
#' #'
#' #'   lasso.Y <- rlasso(y ~ . - D - 1, data = trans)
#' #'   Ytilde <- lasso.Y$residuals
#' #'
#' #'   lasso.D <- rlasso(D ~ . - 1, data = trans[, -1])
#' #'   Dtilde <- lasso.D$residuals
#' #'
#' #'   data_res <- data.frame(id = data[[index[1]]],
#' #'                          time = data[[index[2]]],
#' #'                          Ytilde = Ytilde,
#' #'                          Dtilde = Dtilde)
#' #'
#' #'   Post_plm <- plm(Ytilde ~ -1 + Dtilde, data = data_res, model = "pooling", index = c("id", "time"))
#' #'
#' #'   coefs <- coef(Post_plm)
#' #'   se_corrected <- sqrt(vcovHC(Post_plm, type = "HC0", method = "arellano")) * sqrt(N * T / (N*T - N*G_time - T*G_unit))
#' #'   t_values_corrected <- coefs / se_corrected
#' #'
#' #'   df <- Post_plm$df.residual
#' #'   p_values_corrected <- 2 * pt(-abs(t_values_corrected), df)
#' #'
#' #'   summary_table_correct <- data.frame(
#' #'     Estimate = coefs,
#' #'     `Std. Error corrected` = se_corrected,
#' #'     `t-value corrected` = t_values_corrected,
#' #'     `Pr(>|t|) corrected` = p_values_corrected
#' #'   )
#' #'
#' #'   colnames(summary_table_correct) = c('Estimate', 'Std. Error corrected', 't-value corrected', 'Pr(>|t|) corrected')
#' #'
#' #'   # Return cluster counts as well
#' #'   return(list(
#' #'     fit_summary = summary(fit),
#' #'     G_unit = G_unit,
#' #'     G_time = G_time,
#' #'     unit_group = unit_group,
#' #'     time_group = time_group,
#' #'     post_plm_summary = summary(Post_plm),
#' #'     estimate_corrected = summary_table_correct,
#' #'     summary_table = summary_table_correct
#' #'   ))
#' #' }
#' #'
#' #' # HP_estimate_twoway <- function(data, y_col = NULL, covariate_cols = NULL,
#' #' #                         id_col = "id", time_col = "time",
#' #' #                         index = c("id", "time"), method_auto = "dynamicTreeCut", cluster = NULL) {
#' #' #
#' #' #
#' #' #   if (is.null(covariate_cols)) {
#' #' #     covariate_cols <- setdiff(colnames(data), c(id_col, time_col, y_col))
#' #' #   }
#' #' #
#' #' #   ids <- sort(unique(data[[id_col]]))
#' #' #   times <- sort(unique(data[[time_col]]))
#' #' #   N <- length(ids)
#' #' #   T <- length(times)
#' #' #   K <- length(covariate_cols)
#' #' #
#' #' #   # Initialize y and X
#' #' #   y <- matrix(NA_real_, nrow = N, ncol = T)
#' #' #   X <- array(NA_real_, dim = c(N, T, K))
#' #' #
#' #' #   for (i in seq_len(nrow(data))) {
#' #' #     id_idx <- which(ids == data[[id_col]][i])
#' #' #     time_idx <- which(times == data[[time_col]][i])
#' #' #     y[id_idx, time_idx] <- data[[y_col]][i]
#' #' #     for (k in seq_along(covariate_cols)) {
#' #' #       X[id_idx, time_idx, k] <- data[[covariate_cols[k]]][i]
#' #' #     }
#' #' #   }
#' #' #
#' #' #   # Cluster (output indicators assumed one-hot encoded)
#' #' #   X_norm <- array(NA_real_, dim = dim(X))  # Same shape: (N, T, K)
#' #' #   for (k in 1:dim(X)[3]) {
#' #' #     X_slice <- X[,,k]
#' #' #     X_norm[,,k] <- scale(X_slice)
#' #' #   }
#' #' #   y_norm = scale(y)
#' #' #   res_unit <- cluster_Hierarchical(y_norm, X_norm, type = "unit", method_auto = method_auto)
#' #' #   res_time <- cluster_Hierarchical(y_norm, X_norm, type = "time", method_auto = method_auto)
#' #' #   res_covar <- cluster_Hierarchical(y_norm, X_norm, type = "covariate", method_auto = method_auto, cluster = cluster)
#' #' #
#' #' #   Dv <- res_unit$indicator    # N x G_unit
#' #' #   Du <- res_time$indicator    # T x G_time
#' #' #   Dc <- res_covar$indicator   # (K+1) x G_covar
#' #' #
#' #' #   G_unit = res_unit$G
#' #' #   G_time = res_time$G
#' #' #   G_covar = res_covar$G
#' #' #
#' #' #   unit_group <- max.col(Dv)
#' #' #   time_group <- max.col(Du)
#' #' #   covar_group <- max.col(Dc)
#' #' #
#' #' #   # Projection matrices
#' #' #   Mv <- diag(N) - Dv %*% solve(t(Dv) %*% Dv) %*% t(Dv)
#' #' #   Mu <- diag(T) - Du %*% solve(t(Du) %*% Du) %*% t(Du)
#' #' #   Mc <- diag(K+1) - Dc %*% solve(t(Dc) %*% Dc) %*% t(Dc)
#' #' #
#' #' #   # Stack y and X into Z: N x T x (K+1)
#' #' #   Z <- array(NA_real_, dim = c(N, T, K + 1))
#' #' #   Z[,,1] <- y
#' #' #   for (k in 1:K) Z[,,k+1] <- X[,,k]
#' #' #
#' #' #   # Projection method: demean unit and time clusters slice-wise
#' #' #   Z_proj <- array(NA_real_, dim = c(N, T, K + 1))
#' #' #   for (k in 1:(K+1)) {
#' #' #     Z_proj[,,k] <- Mv %*% Z[,,k] %*% Mu
#' #' #   }
#' #' #   # Separate transformed y and X
#' #' #   y_trans <- Z_proj[,,1] # Z_proj_final_array[,,1]
#' #' #   X_trans <- Z_proj[,,2:(K+1)] # Z_proj_final_array[,,2:(K+1)]
#' #' #
#' #' #   tY_vector <- as.vector(y_trans)  # (N*T)
#' #' #   tX_matrix <- matrix(aperm(X_trans, c(1, 2, 3)), nrow = N*T, ncol = K)
#' #' #
#' #' #   d_vector <- tX_matrix[, 1]
#' #' #   x_matrix <- tX_matrix[, -1, drop = FALSE]
#' #' #
#' #' #   # Run rlassoEffect with correct inputs
#' #' #   fit <- rlassoEffect(x = x_matrix, y = tY_vector, d = d_vector, method = "partialling out")
#' #' #
#' #' #
#' #' #   trans <- data.frame(y = tY_vector,
#' #' #                       D = d_vector,
#' #' #                       x_matrix)
#' #' #
#' #' #   lasso.Y <- rlasso(y ~ . - D - 1, data = trans)
#' #' #   Ytilde <- lasso.Y$residuals
#' #' #
#' #' #   lasso.D <- rlasso(D ~ . - 1, data = trans[, -1])
#' #' #   Dtilde <- lasso.D$residuals
#' #' #
#' #' #   data_res <- data.frame(id = data[[index[1]]],
#' #' #                          time = data[[index[2]]],
#' #' #                          Ytilde = Ytilde,
#' #' #                          Dtilde = Dtilde)
#' #' #
#' #' #   Post_plm <- plm(Ytilde ~ -1 + Dtilde, data = data_res, model = "pooling", index = c("id", "time"))
#' #' #
#' #' #   coefs <- coef(Post_plm)
#' #' #   se_corrected <- sqrt(vcovHC(Post_plm, type = "HC0", method = "arellano")) * sqrt(N * T / ((N-G_unit)*(T-G_time)))
#' #' #   t_values_corrected <- coefs / se_corrected
#' #' #
#' #' #   df <- Post_plm$df.residual
#' #' #   p_values_corrected <- 2 * pt(-abs(t_values_corrected), df)
#' #' #
#' #' #   summary_table_correct <- data.frame(
#' #' #     Estimate = coefs,
#' #' #     `Std. Error corrected` = se_corrected,
#' #' #     `t-value corrected` = t_values_corrected,
#' #' #     `Pr(>|t|) corrected` = p_values_corrected
#' #' #   )
#' #' #
#' #' #   summary_table <- summary(Post_plm)
#' #' #   summary_table$coefficients <- cbind(summary_table_correct, summary_table$coefficients)
#' #' #
#' #' #   # Return cluster counts as well
#' #' #   return(list(
#' #' #     fit_summary = summary(fit),
#' #' #     G_unit = G_unit,
#' #' #     G_time = G_time,
#' #' #     unit_group = unit_group,
#' #' #     time_group = time_group,
#' #' #     post_plm_summary = summary(Post_plm),
#' #' #     estimate_corrected = summary_table_correct,
#' #' #     summary_table = summary_table
#' #' #   ))
#' #' # }
#' #' #
#' #' # HP_estimate_double_twoway <- function(data, y_col = NULL, covariate_cols = NULL,
#' #' #                                id_col = "id", time_col = "time",
#' #' #                                index = c("id", "time"), method_auto = "dynamicTreeCut", cluster = NULL) {
#' #' #
#' #' #
#' #' #   if (is.null(covariate_cols)) {
#' #' #     covariate_cols <- setdiff(colnames(data), c(id_col, time_col, y_col))
#' #' #   }
#' #' #
#' #' #   ids <- sort(unique(data[[id_col]]))
#' #' #   times <- sort(unique(data[[time_col]]))
#' #' #   N <- length(ids)
#' #' #   T <- length(times)
#' #' #   K <- length(covariate_cols)
#' #' #
#' #' #   # Initialize y and X
#' #' #   y <- matrix(NA_real_, nrow = N, ncol = T)
#' #' #   X <- array(NA_real_, dim = c(N, T, K))
#' #' #
#' #' #   for (i in seq_len(nrow(data))) {
#' #' #     id_idx <- which(ids == data[[id_col]][i])
#' #' #     time_idx <- which(times == data[[time_col]][i])
#' #' #     y[id_idx, time_idx] <- data[[y_col]][i]
#' #' #     for (k in seq_along(covariate_cols)) {
#' #' #       X[id_idx, time_idx, k] <- data[[covariate_cols[k]]][i]
#' #' #     }
#' #' #   }
#' #' #
#' #' #   # Cluster (output indicators assumed one-hot encoded)
#' #' #   X_norm <- array(NA_real_, dim = dim(X))  # Same shape: (N, T, K)
#' #' #   for (k in 1:dim(X)[3]) {
#' #' #     X_slice <- X[,,k]
#' #' #     X_norm[,,k] <- scale(X_slice)
#' #' #   }
#' #' #   y_norm = scale(y)
#' #' #   res_unit <- cluster_Hierarchical(y_norm, X_norm, type = "unit", method_auto = method_auto)
#' #' #   res_time <- cluster_Hierarchical(y_norm, X_norm, type = "time", method_auto = method_auto)
#' #' #   res_covar <- cluster_Hierarchical(y_norm, X_norm, type = "covariate", method_auto = method_auto, cluster = cluster)
#' #' #
#' #' #   Dv <- res_unit$indicator    # N x G_unit
#' #' #   Du <- res_time$indicator    # T x G_time
#' #' #   Dc <- res_covar$indicator   # (K+1) x G_covar
#' #' #
#' #' #   G_unit = res_unit$G
#' #' #   G_time = res_time$G
#' #' #   G_covar = res_covar$G
#' #' #
#' #' #   unit_group <- max.col(Dv)
#' #' #   time_group <- max.col(Du)
#' #' #   covar_group <- max.col(Dc)
#' #' #
#' #' #   # Projection matrices
#' #' #   Mv <- diag(N) - Dv %*% solve(t(Dv) %*% Dv) %*% t(Dv)
#' #' #   Mu <- diag(T) - Du %*% solve(t(Du) %*% Du) %*% t(Du)
#' #' #   Mc <- diag(K+1) - Dc %*% solve(t(Dc) %*% Dc) %*% t(Dc)
#' #' #
#' #' #   # Stack y and X into Z: N x T x (K+1)
#' #' #   Z <- array(NA_real_, dim = c(N, T, K + 1))
#' #' #   Z[,,1] <- y
#' #' #   for (k in 1:K) Z[,,k+1] <- X[,,k]
#' #' #
#' #' #   # Projection method: demean unit and time clusters slice-wise
#' #' #   Z_proj <- array(NA_real_, dim = c(N, T, K + 1))
#' #' #   for (k in 1:(K+1)) {
#' #' #     Z_proj[,,k] <- Mv %*% Z[,,k] %*% Mu
#' #' #   }
#' #' #
#' #' #   # Flatten Z_proj for covariate demeaning
#' #' #   Z_proj_mat <- matrix(NA_real_, nrow = N * T, ncol = K + 1)
#' #' #   for (k in 1:(K+1)) {
#' #' #     Z_proj_mat[,k] <- as.vector(t(Z_proj[,,k]))
#' #' #   }
#' #' #
#' #' #   # Covariate demeaning via projection
#' #' #   Z_proj_final <- Z_proj_mat %*% Mc
#' #' #
#' #' #   # Reshape back to array
#' #' #   Z_proj_final_array <- array(NA_real_, dim = c(N, T, K + 1))
#' #' #   for (k in 1:(K+1)) {
#' #' #     Z_proj_final_array[,,k] <- matrix(Z_proj_final[,k], nrow = N, ncol = T, byrow = TRUE)
#' #' #   }
#' #' #
#' #' #   first_val <- covar_group[1]
#' #' #   count_first <- sum(covar_group == first_val)
#' #' #
#' #' #   if (count_first == 1) {
#' #' #     # Find the most frequent value in covar_group
#' #' #     most_freq_val <- as.numeric(names(sort(table(covar_group), decreasing = TRUE)[1]))
#' #' #     covar_group[1] <- most_freq_val
#' #' #   }
#' #' #   #covar_group[1:length(covar_group)]=1
#' #' #   #Z_proj = group_demean_formula_cpp(Z, unit_group, time_group, covar_group)
#' #' #   y_trans <- Z_proj[,,1] # Z_proj_final_array[,,1]
#' #' #   X_trans <- Z_proj[,,2:(K+1)] # Z_proj_final_array[,,2:(K+1)]
#' #' #
#' #' #   tY_vector <- as.vector(y_trans)  # (N*T)
#' #' #   tX_matrix <- matrix(aperm(X_trans, c(1, 2, 3)), nrow = N*T, ncol = K)
#' #' #
#' #' #   d_vector <- tX_matrix[, 1]
#' #' #   x_matrix <- tX_matrix[, -1, drop = FALSE]
#' #' #   # Separate transformed y and X
#' #' #
#' #' #   trans <- data.frame(y = tY_vector,
#' #' #                       D = d_vector,
#' #' #                       x_matrix)
#' #' #
#' #' #   trans = as.data.frame(trans)
#' #' #
#' #' #   dml_data <- DoubleMLData$new(
#' #' #     data = trans,
#' #' #     y_col = "y",
#' #' #     d_cols = "D",
#' #' #     x_cols = c(colnames(trans)[-c(1,2)])
#' #' #
#' #' #   )
#' #' #
#' #' #   # Define LASSO machine learning learners for nuisance parameter estimation
#' #' #   ml_l <- lrn("regr.cv_glmnet", s = "lambda.min")  # Outcome regression model
#' #' #   ml_m <- lrn("regr.cv_glmnet", s = "lambda.min")  # Treatment model (if applicable)
#' #' #
#' #' #   # Fit the Double Machine Learning model for treatment effect estimation
#' #' #   dml_plr <- DoubleMLPLR$new(dml_data, ml_l = ml_l, ml_m = ml_m)
#' #' #
#' #' #   # Fit the model to estimate the causal effect
#' #' #   dml_plr$fit(store_predictions=TRUE)
#' #' #   g_hat <- dml_plr$predictions$ml_l
#' #' #   m_hat <- dml_plr$predictions$ml_m
#' #' #
#' #' #   # Step 2: Compute residuals
#' #' #   Ytilde <- tY_vector - g_hat
#' #' #   Dtilde <- d_vector - m_hat
#' #' #   # Double ML
#' #' #   data_res = data.frame(id = data[[index[1]]], time = data[[index[2]]], Ytilde, Dtilde)
#' #' #   Post_plm = plm(Ytilde ~ -1 + Dtilde, data = data_res, model = "pooling", index=c("id", "time"))
#' #' #
#' #' #   coefs <- coef(Post_plm)
#' #' #   se_corrected <- sqrt(vcovHC(Post_plm, type = "HC0", method = "arellano")) * sqrt(N * T / ((N-G_unit)*(T-G_time)))
#' #' #   t_values_corrected <- coefs / se_corrected
#' #' #
#' #' #   # Calculate p-values from t-distribution for each coefficient
#' #' #   df <- Post_plm$df.residual  # degrees of freedom
#' #' #   p_values_corrected <- 2 * pt(-abs(t_values_corrected), df)
#' #' #
#' #' #   summary_table_correct <- data.frame(
#' #' #     Estimate = coefs,
#' #' #     SE_Corrected = se_corrected,
#' #' #     t_value_Corrected = t_values_corrected,
#' #' #     p_value_Corrected = p_values_corrected
#' #' #   )
#' #' #   colnames(summary_table_correct) = c('Estimate', 'Std. Error corrected', 't-value corrected', 'Pr(>|t|) corrected')
#' #' #   summary_table = summary(Post_plm)
#' #' #   summary_table$coefficients = cbind(summary_table_correct, summary_table$coefficients)
#' #' #
#' #' #   # Return cluster counts as well
#' #' #   return(list(
#' #' #     fit_summary = dml_plr$summary(),
#' #' #     G_unit = G_unit,
#' #' #     G_time = G_time,
#' #' #     G_covar = G_covar,
#' #' #     unit_group = unit_group,
#' #' #     time_group = time_group,
#' #' #     covar_group = covar_group,
#' #' #     post_plm_summary = summary(Post_plm),
#' #' #     estimate_corrected = summary_table_correct,
#' #' #     summary_table = summary_table
#' #' #   ))
#' #' # }
#' #'
#' #' #' @export
#' #' compute_u_hat <- function(z_array, unit_clusters, time_clusters, covar_clusters) {
#' #'   N <- dim(z_array)[1]
#' #'   T <- dim(z_array)[2]
#' #'   K <- dim(z_array)[3]
#' #'
#' #'   u_hat <- array(0, dim = c(N, T, K))
#' #'
#' #'   for (i in 1:N) {
#' #'     for (t in 1:T) {
#' #'       for (k in 1:K) {
#' #'         g_i <- unit_clusters[i]
#' #'         m_t <- time_clusters[t]
#' #'         l_k <- covar_clusters[k]
#' #'
#' #'         # z_{itk}
#' #'         zitk <- z_array[i, t, k]
#' #'
#' #'         # bar_z_{g_i t k}
#' #'         group_i_indices <- which(unit_clusters == g_i)
#' #'         bar_g_i_t_k <- mean(z_array[group_i_indices, t, k])
#' #'
#' #'         # bar_z_{i m_t k}
#' #'         time_m_indices <- which(time_clusters == m_t)
#' #'         bar_i_m_t_k <- mean(z_array[i, time_m_indices, k])
#' #'
#' #'         # bar_z_{i t l_k}
#' #'         covar_l_indices <- which(covar_clusters == l_k)
#' #'         bar_i_t_l_k <- mean(z_array[i, t, covar_l_indices])
#' #'
#' #'         # bar_z_{g_i m_t k}
#' #'         bar_g_i_m_t_k <- mean(z_array[group_i_indices, time_m_indices, k])
#' #'
#' #'         # bar_z_{g_i t l_k}
#' #'         bar_g_i_t_l_k <- mean(z_array[group_i_indices, t, covar_l_indices])
#' #'
#' #'         # bar_z_{i m_t l_k}
#' #'         bar_i_m_t_l_k <- mean(z_array[i, time_m_indices, covar_l_indices])
#' #'
#' #'         # Final u_hat_{itk}
#' #'         u_hat[i, t, k] <- 3 * zitk -
#' #'           2 * bar_g_i_t_k -
#' #'           2 * bar_i_m_t_k -
#' #'           2 * bar_i_t_l_k +
#' #'           bar_g_i_m_t_k +
#' #'           bar_g_i_t_l_k +
#' #'           bar_i_m_t_l_k
#' #'       }
#' #'     }
#' #'   }
#' #'
#' #'   return(u_hat)
#' #' }
#' #'
#' #' #' @export
#' #' cluster_Hierarchical <- function(y, X, link = "average", threshold = NULL, cluster = NULL,
#' #'                                  type = c("unit", "time", "covariate"),
#' #'                                  method_auto = c("none", "silhouette", "gap", "dynamicTreeCut"),
#' #'                                  data_for_gap = NULL, max_k = 10, deepSplit = TRUE, minClusterSize = 1, pamStage = FALSE) {
#' #'   type <- match.arg(type)
#' #'   method_auto <- match.arg(method_auto)
#' #'
#' #'   N <- nrow(y)
#' #'   T <- ncol(y)
#' #'   K <- dim(X)[3]
#' #'
#' #'   # Combine y and X into N x T x (K+1)
#' #'   combined_array <- array(0, dim = c(N, T, K + 1))
#' #'   combined_array[,,1] <- y
#' #'   combined_array[,,2:(K+1)] <- X
#' #'
#' #'   # Compute distance matrix based on type
#' #'   if (type == "unit") {
#' #'     dist_mat <- pseudo_dist_unit(combined_array)
#' #'   } else if (type == "time") {
#' #'     dist_mat <- pseudo_dist_time(combined_array)
#' #'   } else if (type == "covariate") {
#' #'     dist_mat <- pseudo_dist_covariate(combined_array)
#' #'   } else {
#' #'     stop("Invalid 'type' argument.")
#' #'   }
#' #'
#' #'   dist_obj <- as.dist(dist_mat)
#' #'   hc <- hclust(dist_obj, method = link)
#' #'
#' #'   if (!is.null(cluster)) {
#' #'     G <- cluster
#' #'     clusters <- cutree(hc, k = G)
#' #'
#' #'   } else if (!is.null(threshold)) {
#' #'     clusters <- cutree(hc, h = threshold)
#' #'     G <- length(unique(clusters))
#' #'
#' #'   } else if (method_auto != "none") {
#' #'     max_k <- min(max_k, ifelse(type == "unit", N, ifelse(type == "time", T, K + 1)) - 1)
#' #'
#' #'     if (method_auto == "gap") {
#' #'       if (is.null(data_for_gap)) {
#' #'         stop("For method_auto = 'gap', please provide 'data_for_gap' matrix.")
#' #'       }
#' #'       gap_fun <- function(x, k) {
#' #'         dist_x <- dist(x)
#' #'         hc_x <- hclust(dist_x, method = link)
#' #'         clust <- cutree(hc_x, k = k)
#' #'         list(cluster = clust)  # must return list with $cluster
#' #'       }
#' #'       gap_stat <- cluster::clusGap(data_for_gap, FUN = gap_fun, K.max = max_k, B = 50)
#' #'       G <- maxSE(gap_stat$Tab[, "gap"], gap_stat$Tab[, "SE.sim"], method = "firstSEmax")
#' #'       clusters <- cutree(hc, k = G)
#' #'
#' #'     } else if (method_auto == "silhouette") {
#' #'       sil_scores <- numeric(max_k)
#' #'       sil_scores[1] <- NA # silhouette not defined for k=1
#' #'       for (k in 2:max_k) {
#' #'         clust_try <- cutree(hc, k = k)
#' #'         ss <- silhouette(clust_try, dist_obj)
#' #'         sil_scores[k] <- mean(ss[, 3])
#' #'       }
#' #'       G <- which.max(sil_scores)
#' #'       clusters <- cutree(hc, k = G)
#' #'
#' #'     } else if (method_auto == "dynamicTreeCut") {
#' #'       clusters <- cutreeDynamic(dendro = hc, distM = as.matrix(dist_obj),
#' #'                                 deepSplit = deepSplit, minClusterSize = minClusterSize, pamStage = pamStage)
#' #'       G <- length(unique(clusters[clusters > 0]))  # exclude noise (0)
#' #'       noise_idx <- which(clusters == 0)
#' #'       if (length(noise_idx) > 0) {
#' #'         clusters[noise_idx] <- G + 1
#' #'         G <- G + 1
#' #'       }
#' #'     }
#' #'
#' #'   } else {
#' #'     stop("Must provide either 'cluster', 'threshold', or set 'method_auto' != 'none'.")
#' #'   }
#' #'
#' #'   # Format clusters and indicator matrix (like kmeans())
#' #'   clusters <- as.integer(factor(clusters))
#' #'   G <- length(unique(clusters))
#' #'   indicator <- sapply(1:G, function(g) as.integer(clusters == g))
#' #'   colnames(indicator) <- paste0("Cluster", 1:G)
#' #'
#' #'   dim_cluster <- switch(type,
#' #'                         unit = N,
#' #'                         time = T,
#' #'                         covariate = K + 1)
#' #'
#' #'   rownames(indicator) <- paste0(type, "_", 1:dim_cluster)
#' #'
#' #'   list(G = G, clusters = clusters, indicator = indicator)
#' #' }
#' #'
#' #' #' @export
#' #' cluster_general <- function(Y, X_list, N, T, init = 100, type = "long", groups = NULL, cc = 0, gamma = 1, dim_moment = 1) {
#' #'   dimtheta <- length(X_list)
#' #'   mdim <- dim_moment
#' #'   X_list <- c(X_list, list(Y))
#' #'   K = dimtheta
#' #'   if (type == "long") {
#' #'     mom_i <- c()
#' #'     ## --- Unit-side moments ---
#' #'     X_av <- matrix(0, nrow = N, ncol = dim_moment)  # N x dim_moment
#' #'
#' #'     for (p in 1:dim_moment) {
#' #'       # For each power
#' #'       X_pow <- matrix(0, nrow = N, ncol = T)  # N x T
#' #'
#' #'       for (t in 1:T) {
#' #'         # Collect all covariates at time t
#' #'         X_t <- sapply(X_list, function(X) X[, t])  # N x K
#' #'         # Take power p **before** averaging over K
#' #'         if ( p == 1){
#' #'         X_pow[, t] <- rowMeans(X_t)             # (1/K) sum_k x_itk^p
#' #'         }else if (p == 2){
#' #'           X_pow[, t] <- rowMeans(tanh(X_t))
#' #'         }else if (p == 3){
#' #'           X_pow[, t] <- rowMeans((X_t/(1+abs(X_t))))
#' #'         }else{
#' #'           X_pow[, t] <- rowMeans(X_t^p)
#' #'       }
#' #'       }
#' #'
#' #'       # Average over time
#' #'       X_av[, p] <- rowMeans(X_pow)               # (1/T) sum_t (1/K) sum_k x_itk^p
#' #'     }
#' #'
#' #'     mom_i <- cbind(mom_i, X_av)
#' #'
#' #'     # mom_micro <- c()
#' #'     mom_micro <- c()
#' #'     X_cbind <- do.call(cbind, X_list)
#' #'
#' #'     for (p in 1:dim_moment) {
#' #'       if ( p == 1){
#' #'         mom_micro <- cbind(mom_micro, X_cbind^p)             # (1/K) sum_k x_itk^p
#' #'       }else if (p == 2){
#' #'         mom_micro <- cbind(mom_micro, tanh(X_cbind))
#' #'       }else if (p == 3){
#' #'         mom_micro <- cbind(mom_micro, (X_cbind/(1+abs(X_cbind))))
#' #'       }else{
#' #'         X_pow[, t] <- rowMeans(X_t^p)
#' #'       }
#' #'       mom_micro <- cbind(mom_micro, X_cbind^p)
#' #'     }
#' #'
#' #'     ## --- Demean & rescale ---
#' #'     for (j in 1:mdim) {
#' #'       col_mean <- mean(mom_i[, j])
#' #'       col_sd <- sd(mom_i[, j])
#' #'       if (col_sd > 0) {
#' #'         cols <- ((j - 1) * T * (K+1) + 1):(j * T * (K+1))
#' #'         mom_micro[, cols] <- (mom_micro[, cols] - col_mean) / col_sd
#' #'         mom_i[, j] <- (mom_i[, j] - col_mean) / col_sd
#' #'       }
#' #'     }
#' #'
#' #'
#' #'     # --- Variance / noise on rescaled data ---
#' #'     variance <-  sum(sapply(1:mdim, function(j) {
#' #'       norm(mom_micro[, ((j - 1) * T * (K+1) + 1):(j * T * (K+1))] - mom_i[,j], type = "F")^2
#' #'     })) / (N * T^2 * (K+1)^2)
#' #'
#' #'     data <- mom_i
#' #'     dim_size <- N
#' #'
#' #'   } else if (type == "tall") {
#' #'     ## --- Time-side moments ---
#' #'
#' #'     mom_t <- c()
#' #'     X_av <- matrix(0, nrow = T, ncol = dim_moment)  # T x dim_moment
#' #'
#' #'     for (p in 1:dim_moment) {
#' #'       # For each power
#' #'       X_pow <- matrix(0, nrow = N, ncol = T)  # N x T
#' #'
#' #'       for (i in 1:N) {
#' #'         # Collect all covariates at time t
#' #'         X_i <- sapply(X_list, function(X) X[i, ])  # T x K
#' #'         # Take power p **before** averaging over K
#' #'         X_pow[i, ] <- rowMeans(X_i^p)             # (1/K) sum_k x_itk^p
#' #'       }
#' #'
#' #'       # Average over unit
#' #'       X_av[, p] <- colMeans(X_pow)               # (1/N) sum_i (1/K) sum_k x_itk^p
#' #'     }
#' #'
#' #'     mom_t <- cbind(mom_t, X_av)
#' #'
#' #'     mom_micro2 <- c()
#' #'     X_cbind <- do.call(cbind, lapply(X_list, t))
#' #'
#' #'     for (p in 1:dim_moment) {
#' #'       mom_micro2 <- cbind(mom_micro2, (X_cbind)^p)
#' #'     }
#' #'
#' #'     ## --- Demean & rescale ---
#' #'     for (j in 1:mdim) {
#' #'       col_mean <- mean(mom_t[, j])
#' #'       col_sd <- sd(mom_t[, j])
#' #'       if (col_sd > 0) {
#' #'         cols <- ((j - 1) * N * (K+1) + 1):(j * N * (K+1))
#' #'         mom_micro2[, cols] <- (mom_micro2[, cols] - col_mean) / col_sd
#' #'         mom_t[, j] <- (mom_t[, j] - col_mean) / col_sd
#' #'       }
#' #'     }
#' #'
#' #'
#' #'     ## --- Variance / noise on rescaled data ---
#' #'     variance <- sum(sapply(1:mdim, function(i) {
#' #'       norm(mom_micro2[, ((i - 1) * N * (K+1) + 1):(i * N * (K+1))] - mom_t[,i], type = "F")^2
#' #'     })) / (T * N^2 * (K+1)^2)
#' #'
#' #'     data <- mom_t
#' #'     dim_size <- T
#' #'
#' #'   } else if(type == "long_T_moment") {
#' #'     mom_i <- c()
#' #'     ## --- Unit-side moments ---
#' #'     X_av <- matrix(0, nrow = N, ncol = T)  # N x dim_moment
#' #'
#' #'
#' #'     for (t in 1:T) {
#' #'         # Collect all covariates at time t
#' #'         X_t <- sapply(X_list, function(X) X[, t])  # N x K
#' #'         # Take power p **before** averaging over K
#' #'         X_av[, t] <- rowMeans(X_t)             # (1/K) sum_k x_itk^p
#' #'       }
#' #'
#' #'
#' #'     mom_i <- cbind(mom_i, X_av)
#' #'
#' #'     # mom_micro <- c()
#' #'     mom_micro <- c()
#' #'     X_cbind <- do.call(cbind, X_list)
#' #'     mom_micro <- cbind(mom_micro, X_cbind)
#' #'
#' #'
#' #'     ## --- Demean & rescale ---
#' #'     for (j in 1:T) {
#' #'       col_mean <- mean(mom_i[, j])
#' #'       col_sd <- sd(mom_i[, j])
#' #'       if (col_sd > 0) {
#' #'         cols <- ((j - 1)  * (K+1) + 1):(j * (K+1))
#' #'         mom_micro[, cols] <- (mom_micro[, cols] - col_mean) / col_sd
#' #'         mom_i[, j] <- (mom_i[, j] - col_mean) / col_sd
#' #'       }
#' #'     }
#' #'
#' #'
#' #'     # --- Variance / noise on rescaled data ---
#' #'     variance <-  sum(sapply(1:T, function(j) {
#' #'       norm(mom_micro[, ((j - 1) * (K+1) + 1):(j  * (K+1))] - mom_i[,j], type = "F")^2
#' #'     })) / (N * (K+1)^2)
#' #'
#' #'     data <- mom_i
#' #'     dim_size <- N
#' #'
#' #'   }else {
#' #'     stop("Invalid type. Use 'long' or 'tall'.")
#' #'   }
#' #'
#' #'   ## --- Clustering ---
#' #'   if (!is.null(groups)) {
#' #'     clusters <- groups
#' #'     k_result <- kmeans(data, centers = clusters, algorithm = "Lloyd", nstart = init, iter.max = 100)
#' #'   } else {
#' #'     if (max(apply(data, 2, sd)) > 0) {
#' #'       clusters <- 1
#' #'       xx = 1000
#' #'       while (xx >= gamma * variance ) {
#' #'         k_result <- kmeans(data, centers = clusters, algorithm = "Lloyd", nstart = init, iter.max = 100)
#' #'         xx <- k_result$tot.withinss / dim_size
#' #'         clusters <- clusters + 1
#' #'       }
#' #'       clusters = min(clusters, dim_size-1)
#' #'       k_result <- kmeans(data,centers = clusters,algorithm = "Lloyd",nstart = init,iter.max = 100)
#' #'
#' #'     }else{
#' #'       k_result <- kmeans(data, centers = 1, algorithm = "Lloyd", nstart = init, iter.max = 100)
#' #'     }
#' #'   }
#' #'
#' #'   ## --- Return ---
#' #'   if (type == "long") {
#' #'     list(res = k_result, clusters = clusters, data = data)
#' #'   } else {
#' #'     list(res = k_result, clusters = clusters, data = data)
#' #'   }
#' #' }
#' #'
#' #' fix_time_clusters <- function(cluster, data, m = 2) {
#' #'   # cluster: raw K-means labels (length T)
#' #'   # data: T x 1 matrix (average per time period)
#' #'   # m: minimum cluster size
#' #'   #
#' #'   # Returns:
#' #'   #   cluster: cleaned labels for each time period
#' #'   #   G_time: number of clusters
#' #'
#' #'   cluster <- as.integer(cluster)
#' #'   T <- nrow(data)
#' #'
#' #'   # Collapse to 1 cluster if impossible
#' #'   if (T < 2 * m) {
#' #'     return(list(
#' #'       cluster = rep(1, T),
#' #'       G_time = 1
#' #'     ))
#' #'   }
#' #'
#' #'   repeat {
#' #'     sizes <- table(cluster)
#' #'
#' #'     small <- as.integer(names(sizes[sizes < m]))
#' #'     large <- as.integer(names(sizes[sizes >= m]))
#' #'
#' #'     if (length(small) == 0) break
#' #'
#' #'     # Compute centers (safe for 1-column data)
#' #'     unique_cl <- sort(unique(cluster))
#' #'     centers <- do.call(rbind, lapply(unique_cl, function(cl) {
#' #'       colMeans(data[cluster == cl, , drop = FALSE])
#' #'     }))
#' #'     rownames(centers) <- unique_cl
#' #'
#' #'     # Reassign points from small clusters
#' #'     for (sc in small) {
#' #'       idx <- which(cluster == sc)
#' #'       for (i in idx) {
#' #'         target <- if (length(large) > 0) large else setdiff(unique_cl, sc)
#' #'         dists <- sapply(target, function(cl) {
#' #'           sum((data[i, ] - centers[as.character(cl), ])^2)
#' #'         })
#' #'         cluster[i] <- target[which.min(dists)]
#' #'       }
#' #'     }
#' #'   }
#' #'
#' #'   # Relabel clusters 1:K
#' #'   unique_cl <- sort(unique(cluster))
#' #'   map <- setNames(seq_along(unique_cl), unique_cl)
#' #'   cluster <- map[as.character(cluster)]
#' #'
#' #'   # Number of clusters
#' #'   G_time <- length(unique_cl)
#' #'
#' #'   return(list(
#' #'     cluster = as.integer(cluster),
#' #'     G_time = G_time
#' #'   ))
#' #' }
