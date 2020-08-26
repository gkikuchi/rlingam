
#' @title DirectLiNGAM class
#' @description R implementation of direct LiNGAM algorithm
#'
#' See reference for details of the algorithm
#' @references S. Shimizu, T. Inazumi, Y. Sogawa, A. Hyvärinen, Y. Kawahara, T. Washio, P. O. Hoyer and K. Bollen.
#' DirectLiNGAM: A direct method for learning a linear non-Gaussian structural equation model. Journal of Machine Learning Research, 12(Apr): 1225--1248, 2011.
#' @export
DirectLiNGAM <- R6::R6Class("DirectLiNGAM", inherit = BaseLiNGAM,
    public = list(

        #' @description fit DirectLiNGAM
        #' @param X (numeric matrix or data.frame) data matrix to fit
        fit = function(X) {
            self$estimate_causal_order(X)
            self$estimate_adjacency_matrix(X)
        },

        #' @description search causal ordering
        #' @param X (numerical matrix or data.frame) data matrix
        estimate_causal_order = function(X) {
            num_feat <- ncol(X)

            # Causal Discovery
            U <- c(1:num_feat)
            K <- c()
            X_ <- X

            for (f in 1:num_feat) {
                m <- self$search_exogenous_variable(X_, U)
                for (i in U) {
                    if (i != m) {
                        X_[, i] <- self$residual(X_[, i], X_[, m])
                    }
                }
                K <- c(K, m)
                U <- U[U != m]
            }
            self$causal_order <- K
        },

        #' @description search exogenous variable
        #' @param X (numerical matrix or data.frame) data matrix
        #' @param U (numeric vector) index of each columns
        #' @return index of estimated exogenous variable
        search_exogenous_variable = function(X, U) {
            Uc <- U
            Vj <- NULL

            M_list <- c()
            for (i in Uc) {
                M <- 0
                for (j in U) {
                    if (i != j) {
                        xi_std = (X[, i] - mean(X[, i])) / sd(X[, i])
                        xj_std = (X[, j] - mean(X[, j])) / sd(X[, j])

                        # did not work with ifelse...
                        if ((i %in% Vj) & (j %in% Uc)) {
                            ri_j <- xi_std
                        } else {
                            ri_j <- self$residual(xi_std, xj_std)
                        }

                        if ((j %in% Vj) & (i %in% Uc)) {
                            rj_i <- xj_std
                        } else {
                            rj_i <- self$residual(xj_std, xi_std)
                        }
                        M <- M + min(c(0, self$diff_mutual_info(xi_std, xj_std, ri_j, rj_i)))**2
                    }
                }
                M_list <- c(M_list, -M)
            }
            return(Uc[which.max(M_list)])
        },

        #' @description residual when xi is regressed on xj
        #' @param xi (numeric vector) target variable
        #' @param xj (numeric vector) explanatory variable
        #' @return resid (numeric vector) calculated residual
        residual = function(xi, xj) {
            resid <- xi - (cov(xi, xj) / var(xj)) * xj
            return(resid)
        },

        #' @description calculate the difference of the mutual information
        #' @param xi_std (numeric vector) standardized xi
        #' @param xj_std (numeric vector) standardized xj
        #' @param ri_j (numeric vector) resid of xi_std regressed on xj_std
        #' @param rj_i (numeric vector) resid of xj_std regressed on xi_std
        #' @return scalar value of the difference of mutual information
        diff_mutual_info = function(xi_std, xj_std, ri_j, rj_i) {
            term1 <- self$entropy(xj_std) + self$entropy(ri_j / sd(ri_j))
            term2 <- self$entropy(xi_std) + self$entropy(rj_i / sd(rj_i))
            return(term1 - term2)
        },

        #' @description calculate entropy using maximum entropy approximation
        #' @param u (numeric vector) vector to calculate entropy
        #' @return scalar value of entropy
        entropy = function(u) {
            k1 <- 79.047
            k2 <- 7.4129
            gamma <- 0.37457

            term1 <- (1 + log(2 * pi)) / 2
            term2 <- -k1 * (mean(log(cosh(u))) - gamma)^2
            term3 <- -k2 * (mean(u * exp((-u^2) / 2)))^2
            return(term1 + term2 + term3)
        }
    )
)
