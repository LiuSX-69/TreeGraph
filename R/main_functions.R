#' Tree Graph Angle Analysis
#'
#' Main function for analyzing angles in tree-structured graphs with customizable
#' coverage levels and angle calculation methods.
#'
#' @param data A matrix or data frame with two columns representing coordinates.
#'   If NULL, the function will generate simulated data.
#' @param start Numeric vector of length 2 specifying the starting point coordinates.
#'   Default is c(0, 0).
#' @param coverage Numeric value between 0 and 1 specifying the sample coverage level.
#'   Default is 0.95 (95% coverage).
#' @param angle_method Character string specifying the angle calculation method.
#'   Options are: "mean", "min", "median", "equal", "neigh". Default is "mean".
#' @param simulate Logical indicating whether to run simulation. If TRUE, 
#'   simulation parameters should be provided.
#' @param rho1 First correlation coefficient for simulation. Default is 0.5.
#' @param rho2 Second correlation coefficient for simulation. Default is 0.7.
#' @param T Depth parameter for simulation. Default is 6.
#' @param N Branching factor for simulation. Default is 2.
#' @param Tij Time parameter for simulation. Default is 1.
#' @param n1 Outer repetitions for simulation. Default is 10.
#' @param n2 Inner repetitions for simulation. Default is 100.
#'
#' @return If analyzing data: the selected angle statistic.
#'   If running simulation: comparison results between rho1 and rho2.
#'
#' @examples
#' # Example 1: Analyze sample data
#' set.seed(123)
#' sample_data <- matrix(rnorm(200), ncol = 2)
#' result <- tree_graph_analysis(data = sample_data, 
#'                              coverage = 0.95,
#'                              angle_method = "mean")
#' 
#' # Example 2: Run simulation
#' sim_result <- tree_graph_analysis(simulate = TRUE,
#'                                  rho1 = 0.3, 
#'                                  rho2 = 0.7,
#'                                  coverage = 0.9,
#'                                  angle_method = "median")
#'
#' @export
#' @importFrom MASS mvrnorm
tree_graph_analysis <- function(data = NULL, start = c(0, 0), 
                                coverage = 0.95,
                                angle_method = c("mean", "min", "median", "equal", "neigh"),
                                simulate = FALSE,
                                rho1 = 0.5, rho2 = 0.7,
                                T = 6, N = 2, Tij = 1,
                                n1 = 10, n2 = 100) {
  
  # Input validation
  if (!is.numeric(coverage) || coverage <= 0 || coverage >= 1) {
    stop("Coverage must be a numeric value between 0 and 1")
  }
  
  angle_method <- match.arg(angle_method)
  
  alpha <- 1 - coverage
  
  if (simulate) {
    # Run simulation
    if (!requireNamespace("MASS", quietly = TRUE)) {
      stop("MASS package is required for simulation")
    }
    
    result <- rho12func(rho1, rho2, T, N, Tij, n1, n2, alpha, start)
    return(result)
    
  } else {
    # Analyze provided data
    if (is.null(data)) {
      stop("Please provide data or set simulate = TRUE")
    }
    
    if (!is.matrix(data) && !is.data.frame(data)) {
      stop("Data must be a matrix or data frame")
    }
    
    if (ncol(data) != 2) {
      stop("Data must have exactly 2 columns")
    }
    
    # Convert to matrix if data frame
    data <- as.matrix(data)
    
    # Calculate all angle statistics
    all_results <- anglealpha(data, alpha, start)
    
    # Return the selected method
    switch(angle_method,
           "mean" = all_results$anglemean,
           "min" = all_results$anglemin,
           "median" = all_results$anglemedian,
           "equal" = all_results$angleequal,
           "neigh" = all_results$angleneigh)
  }
}

#' Generate Sample Tree Data
#'
#' Generate sample tree-structured data for testing and demonstration.
#'
#' @param n_samples Number of samples to generate. Default is 100.
#' @param mu Mean vector. Default is c(0, 0).
#' @param sigma Variance vector. Default is c(1, 1).
#' @param rho Correlation coefficient. Default is 0.5.
#'
#' @return Matrix with generated data
#'
#' @examples
#' sample_data <- generate_tree_data(n_samples = 200, rho = 0.7)
#' plot(sample_data)
#'
#' @export
#' @importFrom MASS mvrnorm
generate_tree_data <- function(n_samples = 100, mu = c(0, 0), 
                               sigma = c(1, 1), rho = 0.5) {
  
  cov_matrix <- matrix(c(sigma[1], rho*sqrt(sigma[1]*sigma[2]),
                         rho*sqrt(sigma[1]*sigma[2]), sigma[2]), 2, 2)
  
  data <- MASS::mvrnorm(n_samples, mu, cov_matrix)
  return(data)
}