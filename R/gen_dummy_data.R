
#' generate dummy data
#' @description generates (num_normal + num_anomaly) samples of 6 dimensional data.
#' @param num_samples (integer) number of samples
#' @param random_state (integer) random seed
#' @export
gen_dummy_data <- function(num_samples = 10000, random_state = 1) {
    set.seed(random_state)
    x3 <- runif(n = num_samples)
    x0 <- runif(n = num_samples) + 3 * x3
    x2 <- runif(n = num_samples) + 6 * x3
    x1 <- runif(n = num_samples) + 3 * x0 + 2 * x2
    x5 <- runif(n = num_samples) + 4 * x0
    x4 <- runif(n = num_samples) + 8 * x0 - 1 * x2
    X  <- data.frame(x0, x1, x2, x3, x4, x5)
    return(X)
}
