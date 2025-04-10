#' @import hdm
#' @import plm
#' @useDynLib HDpcluster
#' @export
hdpcluster_ds <- function(y, X, T, D, groups_covariate = NULL, groups_unit = NULL, index, data, type_cluster = 'one way kmeans', pesudo_type = 'seperate', link = 'average') {

  if(type_cluster == 'one way kmeans'){

    X = cbind(D, X)
    N <- dim(X)[1]/T
    K <- dim(X)[2]

    clusteri <- cluster_kmeans(y, X, T = T, type = "long", groups = groups_unit)
    G <- clusteri$clusters
    klong <- clusteri$res

    # Initialize cluster indicator matrices

    Dv <- matrix(0, N, G)

    # Populate cluster indicator matrices

    for (j in seq_len(G)) {
      Dv[, j] <- as.numeric(klong$cluster == j)
    }

    Mv <-  as.matrix(bdiag( replicate(T, diag(N) - Dv %*% solve(t(Dv) %*% Dv) %*% t(Dv), simplify = FALSE) ))
    trans = Mv %*% cbind(y, X)


    # cbind(y, X)[1:N,] - t(sweep( t(t(Dv) %*% cbind(y, X)[1:N,]) , 2, colSums(Dv), "/")[, c(klong$cluster)])

    # Apply projection to Y and X_combined
    tY <- trans[,1]
    tX_combined <- trans[,-1] # Vectorized projection of Y


    # Double ML
    fit = rlassoEffect( x = tX_combined[,-1], y = tY, d = tX_combined[,1], method="partialling out")
    # summary(fit)
    trans = data.frame(y = tY, D = tX_combined[,1], tX_combined[,-1])

    lasso.Y <- rlasso(y ~ . -D - 1 , data = trans )
    Ytilde <- lasso.Y$residuals
    lasso.D <- rlasso(D ~ .-1, data = trans[,-c(1)])
    Dtilde <- lasso.D$residuals
    data_res = data.frame(id = data[[index[1]]], time = data[[index[2]]], Ytilde, Dtilde)
    Post_plm = plm(Ytilde ~ -1 + Dtilde, data = data_res, model = "pooling", index=c("id", "time"))

    coefs <- coef(Post_plm)
    se_corrected <- sqrt(vcovHC(Post_plm, type = "HC0", method = "arellano")) * sqrt(N * T / ((N*T - G)))
    t_values_corrected <- coefs / se_corrected

    # Calculate p-values from t-distribution for each coefficient
    df <- Post_plm$df.residual  # degrees of freedom
    p_values_corrected <- 2 * pt(-abs(t_values_corrected), df)

    summary_table_correct <- data.frame(
      Estimate = coefs,
      SE_Corrected = se_corrected,
      t_value_Corrected = t_values_corrected,
      p_value_Corrected = p_values_corrected
    )
    colnames(summary_table_correct) = c('Estimate', 'Std. Error corrected', 't-value corrected', 'Pr(>|t|) corrected')
    summary_table = summary(Post_plm)
    summary_table$coefficients = cbind(summary_table_correct, summary_table$coefficients)
  }else if (type_cluster == 'one way pesudo'){

    X = cbind(D, X)
    N <- dim(X)[1]/T
    K <- dim(X)[2]

    if (T == 1){
      clusteri <- cluster_kmeans(y = y, X = X, T = T, type = "long", groups = groups_unit)
      G <- clusteri$clusters
      klong <- clusteri$res


      cluster_hie = cluster_pesudo(y = y, X = X, T = T, cluster = G, type = 'unit', link = link)
      G <- cluster_hie$G
      klong$cluster <- cluster_hie$res


      Dv <- matrix(0, N, G)

      # Populate cluster indicator matrices

      for (j in seq_len(G)) {
        Dv[, j] <- as.numeric(klong$cluster == j)
      }

      Mv <-  as.matrix(bdiag( replicate(T, diag(N) - Dv %*% solve(t(Dv) %*% Dv) %*% t(Dv), simplify = FALSE) ))
      trans = Mv %*% cbind(y, X)


      # cbind(y, X)[1:N,] - t(sweep( t(t(Dv) %*% cbind(y, X)[1:N,]) , 2, colSums(Dv), "/")[, c(klong$cluster)])

      # Apply projection to Y and X_combined
      tY <- trans[,1]
      tX_combined <- trans[,-1] # Vectorized projection of Y


      # Double ML
      fit = rlassoEffect( x = tX_combined[,-1], y = tY, d = tX_combined[,1], method="partialling out")
      # summary(fit)
      trans = data.frame(y = tY, D = tX_combined[,1], tX_combined[,-1])

      lasso.Y <- rlasso(y ~ . -D - 1 , data = trans )
      Ytilde <- lasso.Y$residuals
      lasso.D <- rlasso(D ~ .-1,data = trans[,-c(1)])
      Dtilde <- lasso.D$residuals
      data_res = data.frame(id = data[[index[1]]], time = data[[index[2]]], Ytilde, Dtilde)
      Post_plm = plm(Ytilde ~ -1 + Dtilde, data = data_res, model = "pooling", index=c("id", "time"))

      coefs <- coef(Post_plm)
      se_corrected <- sqrt(vcovHC(Post_plm, type = "HC0", method = "arellano")) * sqrt(N * T / ((N*T - G)))
      t_values_corrected <- coefs / se_corrected

      # Calculate p-values from t-distribution for each coefficient
      df <- Post_plm$df.residual  # degrees of freedom
      p_values_corrected <- 2 * pt(-abs(t_values_corrected), df)

      summary_table_correct <- data.frame(
        Estimate = coefs,
        SE_Corrected = se_corrected,
        t_value_Corrected = t_values_corrected,
        p_value_Corrected = p_values_corrected
      )
      colnames(summary_table_correct) = c('Estimate', 'Std. Error corrected', 't-value corrected', 'Pr(>|t|) corrected')
      summary_table = summary(Post_plm)
      summary_table$coefficients = cbind(summary_table_correct, summary_table$coefficients)
    }
    else{

      if (pesudo_type == 'seperate'){

        trans = c()
        for (t in 1:T){

          clusteri <- cluster_kmeans(y = y[(N*(t-1)+1):(N*t)], X = X[(N*(t-1)+1):(N*t),], T = 1, type = "long", groups = groups_unit)
          G <- clusteri$clusters
          klong <- clusteri$res


          cluster_hie = cluster_pesudo(y = y[(N*(t-1)+1):(N*t)], X = X[(N*(t-1)+1):(N*t),], T = 1, cluster = G, type = 'unit', link = link)
          G <- cluster_hie$G
          klong$cluster <- cluster_hie$res


          Dv <- matrix(0, N, G)

          # Populate cluster indicator matrices

          for (j in seq_len(G)) {
            Dv[, j] <- as.numeric(klong$cluster == j)
          }

          Mv <-  as.matrix(bdiag( replicate(1, diag(N) - Dv %*% solve(t(Dv) %*% Dv) %*% t(Dv), simplify = FALSE) ))
          trans = rbind(trans, Mv %*% cbind(y[(N*(t-1)+1):(N*t)], X[(N*(t-1)+1):(N*t),]))

        }
        # cbind(y, X)[1:N,] - t(sweep( t(t(Dv) %*% cbind(y, X)[1:N,]) , 2, colSums(Dv), "/")[, c(klong$cluster)])
      }else if (pesudo_type == 'average'){

        clusteri <- cluster_kmeans(y = y, X = X, T = T, type = "long", groups = groups_unit)
        G <- clusteri$clusters
        klong <- clusteri$res


        cluster_hie = cluster_pesudo(y = y, X = X, T = T, cluster = G, type = 'average', link = link)
        G <- cluster_hie$G
        klong$cluster <- cluster_hie$res


        Dv <- matrix(0, N, G)

        # Populate cluster indicator matrices

        for (j in seq_len(G)) {
          Dv[, j] <- as.numeric(klong$cluster == j)
        }

        Mv <-  as.matrix(bdiag( replicate(T, diag(N) - Dv %*% solve(t(Dv) %*% Dv) %*% t(Dv), simplify = FALSE) ))
        trans =  Mv %*% cbind(y, X)
      }

      # Apply projection to Y and X_combined
      tY <- trans[,1]
      tX_combined <- trans[,-1] # Vectorized projection of Y


      # Double ML
      fit = rlassoEffect( x = tX_combined[,-1], y = tY, d = tX_combined[,1], method="partialling out")
      # summary(fit)
      trans = data.frame(y = tY, D = tX_combined[,1], tX_combined[,-1])

      lasso.Y <- rlasso(y ~ . -D - 1 , data = trans )
      Ytilde <- lasso.Y$residuals
      lasso.D <- rlasso(D ~ .-1, data = trans[,-c(1)])
      Dtilde <- lasso.D$residuals
      data_res = data.frame(id = data[[index[1]]], time = data[[index[2]]], Ytilde, Dtilde)
      Post_plm = plm(Ytilde ~ -1 + Dtilde, data = data_res, model = "pooling", index=c("id", "time"))

      coefs <- coef(Post_plm)
      se_corrected <- sqrt(vcovHC(Post_plm, type = "HC0", method = "arellano")) * sqrt(N * T / ((N*T - G)))
      t_values_corrected <- coefs / se_corrected

      # Calculate p-values from t-distribution for each coefficient
      df <- Post_plm$df.residual  # degrees of freedom
      p_values_corrected <- 2 * pt(-abs(t_values_corrected), df)

      summary_table_correct <- data.frame(
        Estimate = coefs,
        SE_Corrected = se_corrected,
        t_value_Corrected = t_values_corrected,
        p_value_Corrected = p_values_corrected
      )
      colnames(summary_table_correct) = c('Estimate', 'Std. Error corrected', 't-value corrected', 'Pr(>|t|) corrected')
      summary_table = summary(Post_plm)
      summary_table$coefficients = cbind(summary_table_correct, summary_table$coefficients)

    }

  }

  return(list(res = summary(fit), G = G, res_seperate = summary(Post_plm), estimate_correct = summary_table_correct, summary_table = summary_table))
}

#' @import DoubleML
#' @import mlr3learners
#' @import mlr3
#' @export
hdpcluster_dml <- function(y, X, T, D, groups_covariate = NULL, groups_unit = NULL, index, data, type_cluster = 'one way kmeans', pesudo_type = "seperate", link = 'average') {

  if(type_cluster == 'one way kmeans'){

    X = cbind(D, X)
    N <- dim(X)[1]/T
    K <- dim(X)[2]

    clusteri <- cluster_kmeans(y, X, T = T, type = "long", groups = groups_unit)
    G <- clusteri$clusters
    klong <- clusteri$res

    # Initialize cluster indicator matrices

    Dv <- matrix(0, N, G)

    # Populate cluster indicator matrices

    for (j in seq_len(G)) {
      Dv[, j] <- as.numeric(klong$cluster == j)
    }

    Mv <-  as.matrix(bdiag( replicate(T, diag(N) - Dv %*% solve(t(Dv) %*% Dv) %*% t(Dv), simplify = FALSE) ))
    trans = Mv %*% cbind(y, X)


    # cbind(y, X)[1:N,] - t(sweep( t(t(Dv) %*% cbind(y, X)[1:N,]) , 2, colSums(Dv), "/")[, c(klong$cluster)])

    # Apply projection to Y and X_combined
    tY <- trans[,1]
    tX_combined <- trans[,-1] # Vectorized projection of Y

    trans = data.frame(y = tY, D = tX_combined[,1], tX_combined[,-1])

    dml_data <- DoubleMLData$new(
      data = trans,
      y_col = "y",
      d_cols = "D",
      x_cols = c(colnames(trans)[-c(1,2)])
    )

    # Define LASSO machine learning learners for nuisance parameter estimation
    ml_l <- lrn("regr.cv_glmnet", s = "lambda.min")  # Outcome regression model
    ml_m <- lrn("regr.cv_glmnet", s = "lambda.min")  # Treatment model (if applicable)

    # Fit the Double Machine Learning model for treatment effect estimation
    dml_plr <- DoubleMLPLR$new(dml_data, ml_l = ml_l, ml_m = ml_m)

    # Fit the model to estimate the causal effect
    dml_plr$fit(store_predictions=TRUE)
    g_hat <- dml_plr$predictions$ml_l
    m_hat <- dml_plr$predictions$ml_m

    # Step 2: Compute residuals
    Ytilde <- tY - g_hat
    Dtilde <- tX_combined[,1] - m_hat

    # Double ML
    data_res = data.frame(id = data[[index[1]]], time = data[[index[2]]], Ytilde, Dtilde)
    Post_plm = plm(Ytilde ~ -1 + Dtilde, data = data_res, model = "pooling", index=c("id", "time"))

    coefs <- coef(Post_plm)
    se_corrected <- sqrt(vcovHC(Post_plm, type = "HC0", method = "arellano")) * sqrt(N * T / ((N*T - G)))
    t_values_corrected <- coefs / se_corrected

    # Calculate p-values from t-distribution for each coefficient
    df <- Post_plm$df.residual  # degrees of freedom
    p_values_corrected <- 2 * pt(-abs(t_values_corrected), df)

    summary_table_correct <- data.frame(
      Estimate = coefs,
      SE_Corrected = se_corrected,
      t_value_Corrected = t_values_corrected,
      p_value_Corrected = p_values_corrected
    )
    colnames(summary_table_correct) = c('Estimate', 'Std. Error corrected', 't-value corrected', 'Pr(>|t|) corrected')
    summary_table = summary(Post_plm)
    summary_table$coefficients = cbind(summary_table_correct, summary_table$coefficients)

  }else if (type_cluster == 'one way pesudo'){

    X = cbind(D, X)
    N <- dim(X)[1]/T
    K <- dim(X)[2]

    if (T == 1){
      clusteri <- cluster_kmeans(y = y, X = X, T = T, type = "long", groups = groups_unit)
      G <- clusteri$clusters
      klong <- clusteri$res


      cluster_hie = cluster_pesudo(y = y, X = X, T = T, cluster = G, type = 'unit', link = link)
      G <- cluster_hie$G
      klong$cluster <- cluster_hie$res


      Dv <- matrix(0, N, G)

      # Populate cluster indicator matrices

      for (j in seq_len(G)) {
        Dv[, j] <- as.numeric(klong$cluster == j)
      }

      Mv <-  as.matrix(bdiag( replicate(T, diag(N) - Dv %*% solve(t(Dv) %*% Dv) %*% t(Dv), simplify = FALSE) ))
      trans = Mv %*% cbind(y, X)


      # cbind(y, X)[1:N,] - t(sweep( t(t(Dv) %*% cbind(y, X)[1:N,]) , 2, colSums(Dv), "/")[, c(klong$cluster)])

      # Apply projection to Y and X_combined
      tY <- trans[,1]
      tX_combined <- trans[,-1] # Vectorized projection of Y

      trans = data.frame(y = tY, D = tX_combined[,1], tX_combined[,-1])

      dml_data <- DoubleMLData$new(
        data = trans,
        y_col = "y",
        d_cols = "D",
        x_cols = c(colnames(trans)[-c(1,2)])
      )

      # Define LASSO machine learning learners for nuisance parameter estimation
      ml_l <- lrn("regr.cv_glmnet", s = "lambda.min")  # Outcome regression model
      ml_m <- lrn("regr.cv_glmnet", s = "lambda.min")  # Treatment model (if applicable)

      # Fit the Double Machine Learning model for treatment effect estimation
      dml_plr <- DoubleMLPLR$new(dml_data, ml_l = ml_l, ml_m = ml_m)

      # Fit the model to estimate the causal effect
      dml_plr$fit(store_predictions=TRUE)
      g_hat <- dml_plr$predictions$ml_l
      m_hat <- dml_plr$predictions$ml_m

      # Step 2: Compute residuals
      Ytilde <- tY - g_hat
      Dtilde <- tX_combined[,1] - m_hat

      # Double ML
      data_res = data.frame(id = data[[index[1]]], time = data[[index[2]]], Ytilde, Dtilde)
      Post_plm = plm(Ytilde ~ -1 + Dtilde, data = data_res, model = "pooling", index=c("id", "time"))

      coefs <- coef(Post_plm)
      se_corrected <- sqrt(vcovHC(Post_plm, type = "HC0", method = "arellano")) * sqrt(N * T / ((N*T - G)))
      t_values_corrected <- coefs / se_corrected

      # Calculate p-values from t-distribution for each coefficient
      df <- Post_plm$df.residual  # degrees of freedom
      p_values_corrected <- 2 * pt(-abs(t_values_corrected), df)

      summary_table_correct <- data.frame(
        Estimate = coefs,
        SE_Corrected = se_corrected,
        t_value_Corrected = t_values_corrected,
        p_value_Corrected = p_values_corrected
      )
      colnames(summary_table_correct) = c('Estimate', 'Std. Error corrected', 't-value corrected', 'Pr(>|t|) corrected')
      summary_table = summary(Post_plm)
      summary_table$coefficients = cbind(summary_table_correct, summary_table$coefficients)
    }else{
      if (pesudo_type == 'seperate'){

        trans = c()
        for (t in 1:T){

          clusteri <- cluster_kmeans(y = y[(N*(t-1)+1):(N*t)], X = X[(N*(t-1)+1):(N*t),], T = 1, type = "long", groups = groups_unit)
          G <- clusteri$clusters
          klong <- clusteri$res


          cluster_hie = cluster_pesudo(y = y[(N*(t-1)+1):(N*t)], X = X[(N*(t-1)+1):(N*t),], T = 1, cluster = G, type = 'unit', link = link)
          G <- cluster_hie$G
          klong$cluster <- cluster_hie$res


          Dv <- matrix(0, N, G)

          # Populate cluster indicator matrices

          for (j in seq_len(G)) {
            Dv[, j] <- as.numeric(klong$cluster == j)
          }

          Mv <-  as.matrix(bdiag( replicate(1, diag(N) - Dv %*% solve(t(Dv) %*% Dv) %*% t(Dv), simplify = FALSE) ))
          trans = rbind(trans, Mv %*% cbind(y[(N*(t-1)+1):(N*t)], X[(N*(t-1)+1):(N*t),]))

        }
        # cbind(y, X)[1:N,] - t(sweep( t(t(Dv) %*% cbind(y, X)[1:N,]) , 2, colSums(Dv), "/")[, c(klong$cluster)])
      }else if (pesudo_type == 'average'){

        clusteri <- cluster_kmeans(y = y, X = X, T = T, type = "long", groups = groups_unit)
        G <- clusteri$clusters
        klong <- clusteri$res


        cluster_hie = cluster_pesudo(y = y, X = X, T = T, cluster = G, type = 'average', link = link)
        G <- cluster_hie$G
        klong$cluster <- cluster_hie$res


        Dv <- matrix(0, N, G)

        # Populate cluster indicator matrices

        for (j in seq_len(G)) {
          Dv[, j] <- as.numeric(klong$cluster == j)
        }

        Mv <-  as.matrix(bdiag( replicate(T, diag(N) - Dv %*% solve(t(Dv) %*% Dv) %*% t(Dv), simplify = FALSE) ))
        trans =  Mv %*% cbind(y, X)
      }
      # cbind(y, X)[1:N,] - t(sweep( t(t(Dv) %*% cbind(y, X)[1:N,]) , 2, colSums(Dv), "/")[, c(klong$cluster)])

      # Apply projection to Y and X_combined
      tY <- trans[,1]
      tX_combined <- trans[,-1] # Vectorized projection of Y
      trans = data.frame(y = tY, D = tX_combined[,1], tX_combined[,-1])

      dml_data <- DoubleMLData$new(
        data = trans,
        y_col = "y",
        d_cols = "D",
        x_cols = c(colnames(trans)[-c(1,2)])
      )

      # Define LASSO machine learning learners for nuisance parameter estimation
      ml_l <- lrn("regr.cv_glmnet", s = "lambda.min")  # Outcome regression model
      ml_m <- lrn("regr.cv_glmnet", s = "lambda.min")  # Treatment model (if applicable)

      # Fit the Double Machine Learning model for treatment effect estimation
      dml_plr <- DoubleMLPLR$new(dml_data, ml_l = ml_l, ml_m = ml_m)

      # Fit the model to estimate the causal effect
      dml_plr$fit(store_predictions=TRUE)
      g_hat <- dml_plr$predictions$ml_l
      m_hat <- dml_plr$predictions$ml_m

      # Step 2: Compute residuals
      Ytilde <- tY - g_hat
      Dtilde <- tX_combined[,1] - m_hat

      # Double ML
      data_res = data.frame(id = data[[index[1]]], time = data[[index[2]]], Ytilde, Dtilde)
      Post_plm = plm(Ytilde ~ -1 + Dtilde, data = data_res, model = "pooling", index=c("id", "time"))

      coefs <- coef(Post_plm)
      se_corrected <- sqrt(vcovHC(Post_plm, type = "HC0", method = "arellano")) * sqrt(N * T / ((N*T - G)))
      t_values_corrected <- coefs / se_corrected

      # Calculate p-values from t-distribution for each coefficient
      df <- Post_plm$df.residual  # degrees of freedom
      p_values_corrected <- 2 * pt(-abs(t_values_corrected), df)

      summary_table_correct <- data.frame(
        Estimate = coefs,
        SE_Corrected = se_corrected,
        t_value_Corrected = t_values_corrected,
        p_value_Corrected = p_values_corrected
      )
      colnames(summary_table_correct) = c('Estimate', 'Std. Error corrected', 't-value corrected', 'Pr(>|t|) corrected')
      summary_table = summary(Post_plm)
      summary_table$coefficients = cbind(summary_table_correct, summary_table$coefficients)

    }

  }

  return(list(res = dml_plr$summary(), G = G, res_seperate = summary(Post_plm), estimate_correct = summary_table_correct, summary_table = summary_table))
}

#' @import hdm
#' @export
#'
double_selection_TWFE = function(Y, D, X, data){
  id = data$id
  time = data$time

  id_matrix <- model.matrix(~ factor(id) - 1)  # Remove intercept
  time_matrix <- model.matrix(~ factor(time) - 1)  # Remove intercept

  # Demeaning Y, D, and X variables
  Y_tilde <- demean_matrix(data$Y, id_matrix, time_matrix)
  D_tilde <- demean_matrix(data$D, id_matrix, time_matrix)
  x_names <- grep("^X", names(data), value = TRUE)  # Extract all columns starting with 'X'

  # Demean each control variable (X) using the same matrix transformation
  X_tilde <- apply(data[, x_names], 2, function(x) demean_matrix(x, id_matrix, time_matrix))

  # ---- Step 2: Apply LASSO Regression Using rlassoEffect ---- #

  # Run rlassoEffect with demeaned variables
  lasso_model <- rlassoEffect(x = X_tilde, y = Y_tilde, d = D_tilde, method = "partialling out")

  trans = data.frame(y = Y_tilde, D = D_tilde, X = X_tilde)

  lasso.Y <- rlasso(y ~ . -D - 1, data = trans )
  Ytilde <- lasso.Y$residuals
  lasso.D <- rlasso(D ~ .-1, data = trans[,-c(1)])
  Dtilde <- lasso.D$residuals
  data_res = data.frame(id = id, time = time, Ytilde, Dtilde)
  Post_plm = plm(Ytilde ~ -1 + Dtilde, data = data_res, model = "pooling", index=c("id", "time"))

  coefs <- coef(Post_plm)
  se_corrected <- sqrt(vcovHC(Post_plm, type = "HC0", method = "arellano"))
  t_values_corrected <- coefs / se_corrected

  # Calculate p-values from t-distribution for each coefficient
  df <- Post_plm$df.residual  # degrees of freedom
  p_values_corrected <- 2 * pt(-abs(t_values_corrected), df)

  summary_table_correct <- data.frame(
    Estimate = coefs,
    SE_Corrected = se_corrected,
    t_value_Corrected = t_values_corrected,
    p_value_Corrected = p_values_corrected
  )
  colnames(summary_table_correct) = c('Estimate', 'Std. Error corrected', 't-value corrected', 'Pr(>|t|) corrected')
  summary_table = summary(Post_plm)
  summary_table$coefficients = cbind(summary_table_correct, summary_table$coefficients)
  list(res = lasso_model, res_seperate = summary(Post_plm), estimate_correct = summary_table_correct, summary_table = summary_table)
}

cluster_kmeans <- function(y, X, T, type = 'long', groups = NULL) {
  # Compute row or column means for Y and each X in X_list
  if (type == "tall") {
    N = dim(X)[1]/T
    K = dim(X)[2]
    y = as.matrix(y)
    Y_mean <- colMeans(y)
    X_means <- colMeans(X)
    data <- c(Y_mean, X_means)
    variance <- sum(sapply(1:N, function(i) {
      norm(as.matrix(c(y[i], X[i, ]) - data), type = "F")^2
    })) / (N^2*(K+1))
    dim_size <- K + 1
  }else if (type == "long") {
    N = dim(X)[1]/T
    K = dim(X)[2]
    y = as.matrix(y)
    data <- rowMeans(cbind(y,X))
    data = matrix(data,N,T)
    variance <- (sum(sapply(1:K, function(j) {
      norm(as.matrix( matrix(X[, j],N,T) - data), type = "F")^2
    })) + norm(as.matrix(matrix(y,N,T) - data), type = "F")^2 ) / (N*(K+1)^2)
    dim_size <- N
  }else {
    stop("Invalid type. Use 'long' for rows or 'all' for columns.")
  }

  # Clustering logic
  if (!is.null(groups)) {
    # Use prespecified number of groups
    clusters <- groups
    k_result <- kmeans(as.matrix(data), centers = clusters, algorithm = "Hartigan-Wong", nstart = 30)
  } else {
    # Determine the number of groups adaptively
    clusters <- 1
    repeat {
      k_result <- kmeans(as.matrix(data), centers = clusters, algorithm = "Hartigan-Wong", nstart = 30)
      if (k_result$tot.withinss / dim_size <=  variance) break
      clusters <- clusters + 1
    }
  }
  list(res = k_result, clusters = clusters)
}

cluster_pesudo <- function(y, X, T, link = "average", threshold, cluster = NULL, type = 'unit'){
  N = dim(X)[1]/T
  K = dim(X)[2]
  data_dist = data.frame(y,X)
  if (type == 'unit'){
    dist_matrix = pseudo_dist(data_dist)
    hc <- hclust(dist_matrix, method = link)
    if (is.na(cluster) == 1){
      clusters <- cutree(hc, h = threshold)
      G = length(unique(clusters))
      res = clusters
    }else{
      clusters <- cutree(hc, k = cluster)
      G = length(unique(clusters))
      res = clusters
    }
  }else if (type == 'covariate'){
    dist_trans = pseudo_dist(t(data_dist))
    hc <- hclust(dist_trans, method = link)
    if (is.na(cluster) == 1){
      clusters <- cutree(hc, h = threshold)
      G = length(unique(clusters))
      res = clusters
    }else{
      clusters <- cutree(hc, k = cluster)
      G = length(unique(clusters))
      res = clusters
    }
  }else if(type == 'average'){
    mat_array <- array(t(data_dist), dim = c(K, N, T))  # transpose first
    avg <- apply(mat_array, c(2, 1), mean)  # result: N x K
    dist_matrix = pseudo_dist(avg)
    hc <- hclust(dist_matrix, method = link)
    if (is.na(cluster) == 1){
      clusters <- cutree(hc, h = threshold)
      G = length(unique(clusters))
      res = clusters
    }else{
      clusters <- cutree(hc, k = cluster)
      G = length(unique(clusters))
      res = clusters
    }
  }
  return(list(res = res, G = G))
}

pseudo_dist <- function(x) {
  x = as.matrix(x)

  # Initialize the distance matrix
  dist_matrix <- matrix(0, nrow(x), nrow(x))

  # Call the Fortran subroutine
  res <- .Fortran("calculate_pseudo_dist",
                  x = as.double(x),
                  dist_matrix = as.double(dist_matrix),
                  N = as.integer(nrow(x)),
                  T = as.integer(ncol(x)))

  # Return the distance matrix
  return(as.dist(matrix(res$dist_matrix, nrow(x), nrow(x))))
}

double_selection_fixed = function(Y, D, X, data){
  id = data$id
  time = data$time

  id_matrix <- model.matrix(~ factor(id) - 1)  # Remove intercept
  time_matrix <- model.matrix(~ factor(time) - 1)  # Remove intercept

  # Demeaning Y, D, and X variables
  Y_tilde <- demean_matrix(data$Y, id_matrix, time_matrix)
  D_tilde <- demean_matrix(data$D, id_matrix, time_matrix)
  x_names <- grep("^X", names(data), value = TRUE)  # Extract all columns starting with 'X'

  # Demean each control variable (X) using the same matrix transformation
  X_tilde <- apply(data[, x_names], 2, function(x) demean_matrix(x, id_matrix, time_matrix))

  # ---- Step 2: Apply LASSO Regression Using rlassoEffect ---- #

  # Run rlassoEffect with demeaned variables
  lasso_model <- rlassoEffect(x = X_tilde, y = Y_tilde, d = D_tilde, method = "partialling out")

  trans = data.frame(y = Y_tilde, D = D_tilde, X = X_tilde)

  lasso.Y <- rlasso(y ~ . -D - 1, data = trans )
  Ytilde <- lasso.Y$residuals
  lasso.D <- rlasso(D ~ .-1, data = trans[,-c(1)])
  Dtilde <- lasso.D$residuals
  data_res = data.frame(id = id, time = time, Ytilde, Dtilde)
  Post_plm = plm(Ytilde ~ -1 + Dtilde, data = data_res, model = "pooling", index=c("id", "time"))

  coefs <- coef(Post_plm)
  se_corrected <- sqrt(vcovHC(Post_plm, type = "HC0", method = "arellano"))
  t_values_corrected <- coefs / se_corrected

  # Calculate p-values from t-distribution for each coefficient
  df <- Post_plm$df.residual  # degrees of freedom
  p_values_corrected <- 2 * pt(-abs(t_values_corrected), df)

  summary_table_correct <- data.frame(
    Estimate = coefs,
    SE_Corrected = se_corrected,
    t_value_Corrected = t_values_corrected,
    p_value_Corrected = p_values_corrected
  )
  colnames(summary_table_correct) = c('Estimate', 'Std. Error corrected', 't-value corrected', 'Pr(>|t|) corrected')
  summary_table = summary(Post_plm)
  summary_table$coefficients = cbind(summary_table_correct, summary_table$coefficients)
  list(res = lasso_model, res_seperate = summary(Post_plm), estimate_correct = summary_table_correct, summary_table = summary_table)
}

# Demean Y, D, and X using matrix transformation
demean_matrix <- function(var, id_matrix, time_matrix) {
  # Subtract the group means for each group (id and time) from the variable
  group_mean_id <- id_matrix %*% solve(t(id_matrix) %*% id_matrix) %*% t(id_matrix) %*% var
  group_mean_time <- time_matrix %*% solve(t(time_matrix) %*% time_matrix) %*% t(time_matrix) %*% var
  var - group_mean_id - group_mean_time
}
