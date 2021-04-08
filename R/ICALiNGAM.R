
#' @title ICALiNGAM class
#' @description R implementation of ICA based LiNGAM algorithm
#'              See reference for details of the algorithm
#' @references
#' S. Shimizu, P. O. Hoyer, A. Hyvärinen, and A. J. Kerminen.
#' A linear non-gaussian acyclic model for causal discovery.
#' Journal of Machine Learning Research, 7:2003-2030, 2006.
#' @importFrom fastICA fastICA
#' @importFrom clue solve_LSAP
#' @export
ICALiNGAM <- R6::R6Class("ICALiNGAM", inherit = BaseLiNGAM,
    public = list(

        #' @field max_iter (integer) maximum iterations for fastICA
        max_iter = integer(),

        #' @description create ICALiNGAM object
        #' @param random_state (integer) Random seed
        #' @param lasso_engine (character) "lars" or "glmnet"
        #' @param max_iter (integer) maximum iterations of fastICA
        initialize = function(random_state = NULL,
                              lasso_engine = "glmnet",
                              max_iter     = 1000) {
            super$initialize(random_state, lasso_engine)
            self$max_iter <- max_iter
        },

        #' @description fit DirectLiNGAM
        #' @param X (numeric matrix or data.frame) data matrix to fit
        fit = function(X) {

            # obtain unmixing matrix W_ica
            ica <- fastICA(X,
                           n.comp  = ncol(X),
                           alg.typ = "parallel",
                           maxit   = self$max_iter,
                           fun     = "logcosh"
            )
            W_ica <- (ica$K %*% ica$W)

            # obtain permuted W_ica
            # note that W_ica is t(W_ica) of scipy's fastICA
            linassign <- solve_LSAP(1 / abs(W_ica), maximum = FALSE)
            col_index <- unname(linassign)
            PW_ica    <- t(W_ica[, col_index])

            # obtain scaling vector
            D <- diag(PW_ica)

            # estimate adjacency matrix
            W_estimate <- PW_ica / D
            B_estimate <- diag(ncol(X)) - W_estimate

            causal_order <- self$estimate_causal_order(B_estimate)
            self$causal_order <- causal_order
            self$estimate_adjacency_matrix(X)
        },

        #' Estimate causal order based on estimated adjacency matrix
        #' @description if B is not DAG, set small elements to zero until DAG is obtained
        #' @param B (matrix) estimated adjacency matrix
        #' @return integer vector of length(ncol(B))
        estimate_causal_order = function(B) {
            d <- ncol(B)

            # set m(m + 2) smallest elements (in absolute value) of B
            pos_vec <- order(abs(B))
            initial_zero_num <- d * (d + 1) / 2
            for (i in pos_vec[1:initial_zero_num]) {
                B[i] <- 0
            }

            # set to zero until DAG is obtained
            for (i in pos_vec[(initial_zero_num + 1):length(pos_vec)]) {
                B[i] <- 0

                # if B is not DAG, null is returned
                causal_order <- self$search_causal_order(B)
                if (!is.null(causal_order)) {
                    break
                }
            }
            return(causal_order)
        },

        #' get causal order of estimated adjacency matrix
        #' @description if given matrix is not DAG, return NULL
        #' @param B (matrix) estimated adjacency matrix
        #' @return vector of length(ncol(B)) or NULL
        search_causal_order = function(B) {
            causal_order   <- c()
            row_num        <- nrow(B)
            original_index <- c(1:row_num)

            while (0 < ncol(B)) {
                # find a row all of which elements are zero
                row_index_list <- which(rowSums(B) == 0)
                if (length(row_index_list) == 0) break

                # append ith to the end
                target_index   <- row_index_list[1]
                causal_order   <- c(causal_order, original_index[target_index])
                original_index <- original_index[-target_index]

                # remove ith row and ith column
                B <- as.matrix(B[-target_index, -target_index])
            }
            if (!length(causal_order) == row_num) {
                causal_order = NULL
            }
            return(causal_order)
        }

    )
)
