
#' @title LiNGAM class
#' @description R implementation for LiNGAM models
#' This class is a base class for LiNGAM.
#'
#' This code is based on Python implementation:
#'
#' (The LiNGAM Project: https://sites.google.com/site/sshimizu06/lingam)
#' @importFrom lars lars
#' @export
BaseLiNGAM <- R6::R6Class("LiNGAM",
    public = list(
        #' @field random_state (integer) Random seed
        #' @field causal_order (numeric vector) Causal oreder of variables
        #' @field adjacency_matrix (numeric matrix) Estimated adjacency matrix
        #' @field intercept (numeric vector) Estimated intercept term
        #' @field lasso_engine (character) "lars" or "glmnet"
        random_state = NULL,
        causal_order = NULL,
        adjacency_matrix = NULL,
        intercept = NULL,
        lasso_engine = NULL,

        #' @description create lingam object
        #' @param random_state (integer) Random seed
        #' @param lasso_engine (character) "lars" or "glmnet"
        initialize = function(random_state = NULL, lasso_engine = "glmnet") {
            self$random_state <- random_state
            self$lasso_engine <- lasso_engine
        },

        #' @description subclasses should implement this method
        #' @param X (numeric matrix or data.frame) data matrix
        fit = function(X) {
        },

        #' @description estimate adjacency matrix based on causal order
        #' @param X (numeric matrix or data.frame) data matrix
        estimate_adjacency_matrix = function(X) {
            # A: intercept
            # B: adjacency matrix
            A <- rep(NA, ncol(X))
            B <- matrix(0, nrow = ncol(X), ncol = ncol(X))
            A[self$causal_order[1]] <- mean(X[, self$causal_order[1]])
            for (i in 2:length(self$causal_order)) {
                res <- self$predict_adaptive_lasso(X, self$causal_order[1:(i - 1)], self$causal_order[i])
                A[self$causal_order[i]] <- res$intercept
                B[self$causal_order[i], self$causal_order[1:(i - 1)]] <- res$coef
            }
            self$intercept <- A
            self$adjacency_matrix <- B
        },

        #' @description fit adaptice lasso
        #' @param X (numeric matrix or data.frame) data matrix
        #' @param predictors (numeric vector) index of explanatory variables
        #' @param target (integer) index of target variable
        #' @param gamma (numeric) data x will be weighted like x^(gamma) for adaptive lasso
        #' @return coef_ (numeric vector) estimated coefficients
        predict_adaptive_lasso = function(X, predictors, target, gamma = 1) {

            # 1st stage (OLS to determine weights)
            fml <- as.formula(paste0(names(X)[target], "~."))
            X_ <- X[names(X)[c(target, predictors)]]
            lr <- lm(fml, data = X_)
            weight <- abs(lr$coefficients[-1])^(gamma)

            # 2nd stage
            x <- as.matrix(X_[, -1, drop = FALSE])
            y <- as.matrix(X_[, 1], drop = FALSE)

            if (self$lasso_engine == "lars") {
                # adaptive lasso by lars()
                # use coefs which minimizes bic
                x <- t(t(x) * weight)
                reg <- lars::lars(x = x, y = y, type = "stepwise")
                bic <- log(nrow(x)) * reg$df + nrow(x) * log(reg$RSS/nrow(x))

                # lars does not explicitly return intercept term...
                idx <- which.min(bic)
                coef_ <- matrix(coef(reg), ncol = ncol(x))[idx ,]
                coef_ <- coef_ * weight
                eval(parse(text = paste0("intercept <- predict(reg, s=", idx, ", data.frame(", paste0(names(coef_), "=", 0, collapse = ","), "))$fit")))

            } else if (self$lasso_engine == "glmnet") {
                # adaptive lasso by glmnet()
                # glmnet() can not handle x with single column

                if (ncol(x) > 1) {
                    # specify lambda sequence (default did not search small lambdas)
                    lambda_seq <- exp(seq(2, -7, length.out = 80))
                    reg <- glmnet::cv.glmnet(x = x, y = y, penalty.factor = 1 / weight, type.measure = "mse", lambda = lambda_seq, relax = FALSE)

                    # use 1se rule
                    idx_best  <- which(reg$lambda == reg$lambda.1se)
                    coef_     <- reg$glmnet.fit$beta[, idx_best]
                    coef_     <- matrix(coef_, ncol = ncol(x))
                    intercept <- reg$glmnet.fit$a0[idx_best]
                } else {
                    coef_ <- matrix(lr$coefficients[-1], ncol = ncol(x))
                    intercept <- lr$coefficients[1]
                }
            }

            # return coefs and intercept
            res <- list()
            res$coef <- coef_
            res$intercept <- intercept
            return(res)
        }
    )
)

