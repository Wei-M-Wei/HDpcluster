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
#' @param pesudo_type if cluster_type = 'unit pesudo', choice of the pesudo type
#' @param link if cluster_type = 'unit pesudo', choice of different links of pesudo type
#' @param optimal_index if cluster_type = 'unit pesudo', different ways to compute the optimal number of clusters
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
hdpcluster_ds <- function(y, D, X, T, groups_unit = NULL, index, data, cluster_type = 'unit kmeans', pesudo_type = 'average', link = 'average', optimal_index = NULL) {

  if (!all(c(index[1], index[2]) %in% colnames(data))) {
    stop("Data must contain names of 'id' and 'time'.")
  }

  # reconstruct the dataset
  data = data.frame(y = y, D = D, X, data[,index[1]], data[,index[2]])
  colnames(data)[(dim(data)[2]-1):dim(data)[2]] = index
  data = data[order(data[,index[2]], data[,index[1]]),]

  # prepare the variables
  y = data$y
  D = data$D
  X = data[, 3:(dim(data)[2]-2)]


  if(cluster_type == 'unit kmeans'){

    X = cbind(D, X)
    N <- dim(X)[1]/T
    K <- dim(X)[2]

    clusteri <- cluster_kmeans(y, X, T = T, type = "long", groups = groups_unit, index = index, data = data)
    G <- clusteri$clusters
    klong <- clusteri$res

    # Initialize cluster indicator matrices
    Du <- matrix(0, N, G)
    Diagu<- matrix(0, G, G)

    # Populate cluster indicator matrices
    for (j in seq_len(G)) {
      Du[, j] <- as.numeric(klong$cluster == j)
      Diagu[j,j] <- 1/sum(Du[,j])
    }
    Mu <-  as.matrix(bdiag( replicate(T, diag(N) - Du %*% Diagu %*% t(Du), simplify = FALSE) ))
    trans = Mu %*% as.matrix(cbind(y, X))

    # Apply projection to Y and X_combined
    tY <- trans[,1]
    tX_combined <- trans[,-1] # Vectorized projection of Y

    # Double ML
    # fit = rlassoEffect( x = tX_combined[,-1], y = tY, d = tX_combined[,1], method="partialling out")

    # Double ML
    trans = data.frame(y = tY, D = tX_combined[,1], tX_combined[,-1])
    lasso.Y <- rlasso(y ~ . -D - 1 , data = trans )
    Ytilde <- lasso.Y$residuals
    lasso.D <- rlasso(D ~ .-1, data = trans[,-c(1)])
    Dtilde <- lasso.D$residuals
    data_res = data.frame(id = data[[index[1]]], time = data[[index[2]]], Ytilde, Dtilde)
    Post_plm = plm(Ytilde ~ -1 + Dtilde, data = data_res, model = "pooling", index=c("id", "time"))

    coefs <- coef(Post_plm)
    se_corrected <- sqrt(vcovHC(Post_plm, type = "HC0", method = "arellano")) * sqrt(N * T / ((N*T - G*T)))
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
    colnames(summary_table_correct) = c('Estimate', 'Std. Error', 't-value', 'Pr(>|t|)')
    summary_table = summary(Post_plm)
    summary_table$coefficients = cbind(summary_table_correct, summary_table$coefficients)

  }else if (cluster_type == 'unit pesudo'){

    X = cbind(D, X)
    N <- dim(X)[1]/T
    K <- dim(X)[2]

    if (T == 1){

      # use the cluster from kmeans
      clusteri <- cluster_kmeans(y = y, X = X, T = T, type = "long", groups = groups_unit  , index = index, data = data)
      G <- clusteri$clusters
      klong <- clusteri$res

      cluster_hie = cluster_pesudo(y = y, X = X, T = T, cluster = G, type = 'unit', link = link, index = index, data = data, optimal_index = optimal_index)
      G <- cluster_hie$G
      klong$cluster <- cluster_hie$res

      # Initialize cluster indicator matrices
      Du <- matrix(0, N, G)
      Diagu<- matrix(0, G, G)

      # Populate cluster indicator matrices
      for (j in seq_len(G)) {
        Du[, j] <- as.numeric(klong$cluster == j)
        Diagu[j,j] <- 1/sum(Du[,j])
      }

      Mu <-  as.matrix(bdiag( replicate(T, diag(N) - Du %*% Diagu %*% t(Du), simplify = FALSE) ))
      trans = Mu %*% as.matrix(cbind(y, X))

      # Apply projection to Y and X_combined
      tY <- trans[,1]
      tX_combined <- trans[,-1] # Vectorized projection of Y


      # Double ML
      # fit = rlassoEffect( x = tX_combined[,-1], y = tY, d = tX_combined[,1], method="partialling out")

      # Double ML
      trans = data.frame(y = tY, D = tX_combined[,1], tX_combined[,-1])

      lasso.Y <- rlasso(y ~ . -D - 1 , data = trans )
      Ytilde <- lasso.Y$residuals
      lasso.D <- rlasso(D ~ .-1,data = trans[,-c(1)])
      Dtilde <- lasso.D$residuals
      data_res = data.frame(id = data[[index[1]]], time = data[[index[2]]], Ytilde, Dtilde)
      Post_plm = plm(Ytilde ~ -1 + Dtilde, data = data_res, model = "pooling", index=c("id", "time"))

      coefs <- coef(Post_plm)
      se_corrected <- sqrt(vcovHC(Post_plm, type = "HC0", method = "arellano")) * sqrt(N * T / ((N*T - G*T)))
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
      colnames(summary_table_correct) = c('Estimate', 'Std. Error', 't-value', 'Pr(>|t|)')
      summary_table = summary(Post_plm)
      summary_table$coefficients = cbind(summary_table_correct, summary_table$coefficients)
    }
    else{

      if (pesudo_type == 'seperate'){

        trans = c()
        for (t in 1:T){

          clusteri <- cluster_kmeans(y = y[(N*(t-1)+1):(N*t)], X = X[(N*(t-1)+1):(N*t),], T = 1, type = "long", groups = groups_unit, index = index, data = data[(N*(t-1)+1):(N*t),])
          G <- clusteri$clusters
          klong <- clusteri$res


          cluster_hie = cluster_pesudo(y = y[(N*(t-1)+1):(N*t)], X = X[(N*(t-1)+1):(N*t),], T = 1, cluster = G, type = 'unit', link = link, index = index, data = data[(N*(t-1)+1):(N*t),], optimal_index = optimal_index)
          G <- cluster_hie$G
          klong$cluster <- cluster_hie$res


          # Initialize cluster indicator matrices
          Du <- matrix(0, N, G)
          Diagu<- matrix(0, G, G)

          # Populate cluster indicator matrices

          for (j in seq_len(G)) {
            Du[, j] <- as.numeric(klong$cluster == j)
            Diagu[j,j] <- 1/sum(Du[,j])
          }

          Mu <-  as.matrix(bdiag( replicate(1, diag(N) - Du %*% Diagu %*% t(Du), simplify = FALSE) ))
          trans = rbind(trans, Mu %*% as.matrix(cbind(y[(N*(t-1)+1):(N*t)], X[(N*(t-1)+1):(N*t),])))

        }
        # cbind(y, X)[1:N,] - t(sweep( t(t(Du) %*% cbind(y, X)[1:N,]) , 2, colSums(Du), "/")[, c(klong$cluster)])
      }else if (pesudo_type == 'average'){

        clusteri <- cluster_kmeans(y = y, X = X, T = T, type = "long", groups = groups_unit, index = index, data = data)
        G <- clusteri$clusters
        klong <- clusteri$res


        cluster_hie = cluster_pesudo(y = y, X = X, T = T, cluster = G, type = 'average', link = link, index = index, data = data, optimal_index = optimal_index)
        G <- cluster_hie$G
        klong$cluster <- cluster_hie$res


        # Initialize cluster indicator matrices
        Du <- matrix(0, N, G)
        Diagu<- matrix(0, G, G)

        # Populate cluster indicator matrices

        for (j in seq_len(G)) {
          Du[, j] <- as.numeric(klong$cluster == j)
          Diagu[j,j] <- 1/sum(Du[,j])
        }
        Mu <-  as.matrix(bdiag( replicate(T, diag(N) - Du %*% Diagu %*% t(Du), simplify = FALSE) ))
        trans =  Mu %*% as.matrix(cbind(y, X))

      }else if (pesudo_type == 'max_out'){

        clusteri <- cluster_kmeans(y = y, X = X, T = T, type = "long", groups = groups_unit, index = index, data = data)
        G <- clusteri$clusters
        klong <- clusteri$res


        cluster_hie = cluster_pesudo(y = y, X = X, T = T, cluster = G, type = 'max_out', link = link, index = index, data = data, optimal_index = optimal_index)
        G <- cluster_hie$G
        klong$cluster <- cluster_hie$res


        # Initialize cluster indicator matrices
        Du <- matrix(0, N, G)
        Diagu<- matrix(0, G, G)

        # Populate cluster indicator matrices
        for (j in seq_len(G)) {
          Du[, j] <- as.numeric(klong$cluster == j)
          Diagu[j,j] <- 1/sum(Du[,j])
        }

        Mu <-  as.matrix(bdiag( replicate(T, diag(N) - Du %*% Diagu %*% t(Du), simplify = FALSE) ))
        trans =  Mu %*% as.matrix(cbind(y, X))
      }

      # Apply projection to Y and X_combined
      tY <- trans[,1]
      tX_combined <- trans[,-1] # Vectorized projection of Y


      # Double ML
      # fit = rlassoEffect( x = tX_combined[,-1], y = tY, d = tX_combined[,1], method="partialling out")

      # Double ML
      trans = data.frame(y = tY, D = tX_combined[,1], tX_combined[,-1])

      lasso.Y <- rlasso(y ~ . -D - 1 , data = trans )
      Ytilde <- lasso.Y$residuals
      lasso.D <- rlasso(D ~ .-1, data = trans[,-c(1)])
      Dtilde <- lasso.D$residuals
      data_res = data.frame(id = data[[index[1]]], time = data[[index[2]]], Ytilde, Dtilde)
      Post_plm = plm(Ytilde ~ -1 + Dtilde, data = data_res, model = "pooling", index=c("id", "time"))

      coefs <- coef(Post_plm)
      se_corrected <- sqrt(vcovHC(Post_plm, type = "HC0", method = "arellano")) * sqrt(N * T / ((N*T - G*T)))
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
      colnames(summary_table_correct) = c('Estimate', 'Std. Error', 't-value', 'Pr(>|t|)')
      summary_table = summary(Post_plm)
      summary_table$coefficients = cbind(summary_table_correct, summary_table$coefficients)

    }

  }

  return(list(res = summary(Post_plm), G = G, estimate_correct = summary_table_correct, summary_table = summary_table))
}

#' @import DoubleML
#' @import mlr3learners
#' @import mlr3
#' @export
hdpcluster_dml <- function(y, D, X, T, groups_unit = NULL, index, data, cluster_type = 'unit kmeans', pesudo_type = "seperate", link = 'average', optimal_index = NULL) {

  if (!all(c(index[1], index[2]) %in% colnames(data))) {
    stop("Data must contain names of 'id' and 'time'.")
  }

  data = data.frame(y = y, D = D, X, data[,index[1]], data[,index[2]])
  colnames(data)[(dim(data)[2]-1):dim(data)[2]] = index
  data = data[order(data[,index[2]], data[,index[1]]),]
  y = data$y
  D = data$D
  data[, 3:(dim(data)[2]-2)]

  if(cluster_type == 'unit kmeans'){

    X = cbind(D, X)
    N <- dim(X)[1]/T
    K <- dim(X)[2]

    clusteri <- cluster_kmeans(y, X, T = T, type = "long", groups = groups_unit, index = index, data = data)
    G <- clusteri$clusters
    klong <- clusteri$res

    # Initialize cluster indicator matrices
    Du <- matrix(0, N, G)
    Diagu<- matrix(0, G, G)

    # Populate cluster indicator matrices
    for (j in seq_len(G)) {
      Du[, j] <- as.numeric(klong$cluster == j)
      Diagu[j,j] <- 1/sum(Du[,j])
    }

    Mu <-  as.matrix(bdiag( replicate(T, diag(N) - Du %*% Diagu %*% t(Du), simplify = FALSE) ))
    trans = Mu %*% as.matrix(cbind(y, X))


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
    se_corrected <- sqrt(vcovHC(Post_plm, type = "HC0", method = "arellano")) * sqrt(N * T / ((N*T - G*T)))
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
    colnames(summary_table_correct) = c('Estimate', 'Std. Error', 't-value', 'Pr(>|t|)')
    summary_table = summary(Post_plm)
    summary_table$coefficients = cbind(summary_table_correct, summary_table$coefficients)

  }else if (cluster_type == 'unit pesudo'){

    X = cbind(D, X)
    N <- dim(X)[1]/T
    K <- dim(X)[2]

    if (T == 1){
      clusteri <- cluster_kmeans(y = y, X = X, T = T, type = "long", groups = groups_unit , index = index, data = data)
      G <- clusteri$clusters
      klong <- clusteri$res


      cluster_hie = cluster_pesudo(y = y, X = X, T = T, cluster = G, type = 'unit', link = link, index = index, data = data, optimal_index = optimal_index)
      G <- cluster_hie$G
      klong$cluster <- cluster_hie$res


      # Initialize cluster indicator matrices
      Du <- matrix(0, N, G)
      Diagu<- matrix(0, G, G)

      # Populate cluster indicator matrices
      for (j in seq_len(G)) {
        Du[, j] <- as.numeric(klong$cluster == j)
        Diagu[j,j] <- 1/sum(Du[,j])
      }

      Mu <-  as.matrix(bdiag( replicate(T, diag(N) - Du %*% Diagu %*% t(Du), simplify = FALSE) ))
      trans = Mu %*% as.matrix(cbind(y, X))


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
      se_corrected <- sqrt(vcovHC(Post_plm, type = "HC0", method = "arellano")) * sqrt(N * T / ((N*T - G*T)))
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
      colnames(summary_table_correct) = c('Estimate', 'Std. Error', 't-value', 'Pr(>|t|)')
      summary_table = summary(Post_plm)
      summary_table$coefficients = cbind(summary_table_correct, summary_table$coefficients)
    }else{
      if (pesudo_type == 'seperate'){

        trans = c()
        for (t in 1:T){

          clusteri <- cluster_kmeans(y = y[(N*(t-1)+1):(N*t)], X = X[(N*(t-1)+1):(N*t),], T = 1, type = "long", groups = groups_unit, index = index, data = data)
          G <- clusteri$clusters
          klong <- clusteri$res


          cluster_hie = cluster_pesudo(y = y[(N*(t-1)+1):(N*t)], X = X[(N*(t-1)+1):(N*t),], T = 1, cluster = G, type = 'unit', link = link, index = index, data = data[(N*(t-1)+1):(N*t),], optimal_index = optimal_index)
          G <- cluster_hie$G
          klong$cluster <- cluster_hie$res


          # Initialize cluster indicator matrices
          Du <- matrix(0, N, G)
          Diagu<- matrix(0, G, G)

          # Populate cluster indicator matrices
          for (j in seq_len(G)) {
            Du[, j] <- as.numeric(klong$cluster == j)
            Diagu[j,j] <- 1/sum(Du[,j])
          }

          Mu <-  as.matrix(bdiag( replicate(1, diag(N) - Du %*% Diagu %*% t(Du), simplify = FALSE) ))
          trans = rbind(trans, Mu %*% as.matrix(cbind(y[(N*(t-1)+1):(N*t)], X[(N*(t-1)+1):(N*t),])))

        }
        # cbind(y, X)[1:N,] - t(sweep( t(t(Du) %*% cbind(y, X)[1:N,]) , 2, colSums(Du), "/")[, c(klong$cluster)])
      }else if (pesudo_type == 'average'){

        clusteri <- cluster_kmeans(y = y, X = X, T = T, type = "long", groups = groups_unit, index = index, data = data)
        G <- clusteri$clusters
        klong <- clusteri$res


        cluster_hie = cluster_pesudo(y = y, X = X, T = T, cluster = G, type = 'average', link = link, index = index, data = data, optimal_index = optimal_index)
        G <- cluster_hie$G
        klong$cluster <- cluster_hie$res


        # Initialize cluster indicator matrices
        Du <- matrix(0, N, G)
        Diagu<- matrix(0, G, G)

        # Populate cluster indicator matrices
        for (j in seq_len(G)) {
          Du[, j] <- as.numeric(klong$cluster == j)
          Diagu[j,j] <- 1/sum(Du[,j])
        }

        Mu <-  as.matrix(bdiag( replicate(T, diag(N) - Du %*% Diagu %*% t(Du), simplify = FALSE) ))
        trans =  Mu %*% as.matrix(cbind(y, X))
      }else if (pesudo_type == 'max_out'){

        clusteri <- cluster_kmeans(y = y, X = X, T = T, type = "long", groups = groups_unit, index = index, data = data)
        G <- clusteri$clusters
        klong <- clusteri$res


        cluster_hie = cluster_pesudo(y = y, X = X, T = T, cluster = G, type = 'max_out', link = link, index = index, data = data, optimal_index = optimal_index)
        G <- cluster_hie$G
        klong$cluster <- cluster_hie$res

        # Initialize cluster indicator matrices

        Du <- matrix(0, N, G)
        Diagu<- matrix(0, G, G)
        # Populate cluster indicator matrices

        for (j in seq_len(G)) {
          Du[, j] <- as.numeric(klong$cluster == j)
          Diagu[j,j] <- 1/sum(Du[,j])
        }

        Mu <-  as.matrix(bdiag( replicate(T, diag(N) - Du %*% Diagu %*% t(Du), simplify = FALSE) ))
        trans =  Mu %*% as.matrix(cbind(y, X))
      }


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
      se_corrected <- sqrt(vcovHC(Post_plm, type = "HC0", method = "arellano")) * sqrt(N * T / ((N*T - G*T)))
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
      colnames(summary_table_correct) = c('Estimate', 'Std. Error', 't-value', 'Pr(>|t|)')
      summary_table = summary(Post_plm)
      summary_table$coefficients = cbind(summary_table_correct, summary_table$coefficients)

    }

  }

  return(list(G = G, res = summary(Post_plm), estimate_correct = summary_table_correct, summary_table = summary_table))
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
  colnames(summary_table_correct) = c('Estimate', 'Std. Error', 't-value', 'Pr(>|t|)')
  summary_table = summary(Post_plm)
  summary_table$coefficients = cbind(summary_table_correct, summary_table$coefficients)
  list(res = lasso_model, res_seperate = summary(Post_plm), estimate_correct = summary_table_correct, summary_table = summary_table)
}

cluster_kmeans <- function(y, X, T, type = 'long', groups = NULL, index, data) {

  data = data.frame(y = y, X, data[,index[1]], data[,index[2]])
  colnames(data)[(dim(data)[2]-1):dim(data)[2]] = index
  data = data[order(data[,index[2]], data[,index[1]]),]

  y = data$y
  X = data[, 2:(dim(data)[2]-2)]
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

cluster_pesudo <- function(y, X, T, link = "average", threshold, cluster = NULL, type = 'unit', index, data, optimal_index = NULL){

  data = data.frame(y = y, X, data[,index[1]], data[,index[2]])
  colnames(data)[(dim(data)[2]-1):dim(data)[2]] = index
  data = data[order(data[,index[2]], data[,index[1]]),]

  y = data$y
  X = data[, 2:(dim(data)[2] - 2)]

  N = dim(X)[1]/T
  K = dim(X)[2]
  data_dist = data.frame(X)
  if (is.null(optimal_index) == 1){
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
    }else if(type == 'max_out'){
      dist_matrix = compute_max_distance_matrix(data_dist, N, T)
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
  }else{

    if (type == 'unit'){
      dist_matrix = pseudo_dist(data_dist)
      clusters = NbClust(data = data_dist, diss = dist_matrix, distance = NULL, method = link, index = optimal_index, min.nc = 1, max.nc = N-1)
      G = length(unique(clusters$Best.partition))
      res = clusters$Best.partition

    }else if (type == 'covariate'){
      dist_trans = pseudo_dist(t(data_dist))
      clusters = NbClust(data = data_dist, diss = dist_trans, distance = NULL, method = link, index = optimal_index, min.nc = 1, max.nc = N-1)
      G = length(unique(clusters$Best.partition))
      res = clusters$Best.partition
    }else if(type == 'average'){
      mat_array <- array(t(data_dist), dim = c(K, N, T))  # transpose first
      avg <- apply(mat_array, c(2, 1), mean)  # result: N x K
      dist_matrix = pseudo_dist(avg)
      clusters = NbClust(data = avg, diss = dist_matrix, distance = NULL, method = link, index = optimal_index, min.nc = 1, max.nc = N-1)
      G = length(unique(clusters$Best.partition))
      res = clusters$Best.partition
    }else if(type == 'max_out'){
      dist_matrix = compute_max_distance_matrix(data_dist, N, T)
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
  colnames(summary_table_correct) = c('Estimate', 'Std. Error', 't-value', 'Pr(>|t|)')
  summary_table = summary(Post_plm)
  summary_table$coefficients = cbind(summary_table_correct, summary_table$coefficients)
  list(res = summary(Post_plm), estimate_correct = summary_table_correct, summary_table = summary_table)
}

# Demean Y, D, and X using matrix transformation
demean_matrix <- function(var, id_matrix, time_matrix) {
  # Subtract the group means for each group (id and time) from the variable
  group_mean_id <- id_matrix %*% solve(t(id_matrix) %*% id_matrix) %*% t(id_matrix) %*% var
  group_mean_time <- time_matrix %*% solve(t(time_matrix) %*% time_matrix) %*% t(time_matrix) %*% var
  var - group_mean_id - group_mean_time
}

dist_to_symmetric_matrix <- function(dist_obj) {
  # Convert the dist object to a full matrix
  full_matrix <- as.matrix(dist_obj)

  # Make the matrix symmetric by mirroring the upper triangle to the lower triangle
  full_matrix[lower.tri(full_matrix)] <- t(full_matrix)[lower.tri(full_matrix)]

  # Set the diagonal to zero (distance from a point to itself is 0)
  diag(full_matrix) <- 0

  return(full_matrix)
}

compute_max_distance_matrix <- function(data_mat, N, T) {
  # Initialize the max_dist_mat as a dist object
  max_dist_mat <- matrix(0, nrow = N, ncol = N)

  for (t in 1:T) {
    # Extract submatrix for time t
    start_idx <- (t - 1) * N + 1
    end_idx <- t * N
    sub_mat <- data_mat[start_idx:end_idx, , drop = FALSE]

    # Compute the distance matrix (N-1 x N)
    dist_mat <- pseudo_dist(sub_mat)

    # Convert dist object to matrix (since dist_mat is of class "dist")
    dist_matrix <- as.matrix(dist_mat)


   max_dist_mat = max_dist_mat + dist_matrix
  }

  # Convert the matrix back to a 'dist' object
  return(as.dist(max_dist_mat/T))
}

