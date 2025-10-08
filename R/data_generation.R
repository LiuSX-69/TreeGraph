#' Generate covariance matrix
#'
#' @param sigma Variance vector
#' @param rho Correlation coefficient
#' @param i Level parameter
#'
#' @return Covariance matrix
#' @keywords internal
sigmamatrix <- function(sigma, rho, i){
  Sigma <- matrix(c(sigma[1], 
                    rho*(1-(i-1)/5)*sqrt(sigma[1]*sigma[2]),
                    rho*(1-(i-1)/5)*sqrt(sigma[1]*sigma[2]),
                    sigma[2]), 2, 2)
  return(Sigma)
}

#' Split data generation
#'
#' @param mu Mean vector
#' @param sigma Variance vector
#' @param rho Correlation coefficient
#' @param i Level parameter
#' @param N Branching factor
#' @param Tij Time parameter
#'
#' @return Generated data
#' @keywords internal
split <- function(mu, sigma, rho, i, N, Tij){
  data <- MASS::mvrnorm(Tij*N^i, mu, sigmamatrix(sigma, rho, i))
  return(data)
}

#' Generate tree-structured residual data
#'
#' @param mu Mean vector
#' @param sigma Variance vector
#' @param rho Correlation coefficient
#' @param T Depth
#' @param N Branching factor
#' @param Tij Time parameter
#'
#' @return Tree-structured data
#' @keywords internal
treedataresidual <- function(mu, sigma, rho, T, N, Tij){
  treeresidual1 <- vector('list', T)
  treeresidual2 <- vector('list', T)
  for(i in 1:T){
    residual1 <- split(mu, sigma, rho, i, N, Tij)
    treeresidual1[[i]] <- residual1[,1]
    treeresidual2[[i]] <- residual1[,2]
  }
  return(list(treeresidual1 = treeresidual1, treeresidual2 = treeresidual2))
}

#' Calculate mean and covariance
#'
#' @param data Input data
#'
#' @return Vector with mean, variances, and correlation
#' @keywords internal
musigma <- function(data){
  mu <- apply(data, 2, mean)
  sigma <- matrix(rep(0, (length(mu))^2), length(mu), length(mu))
  for(i in 1:(dim(data)[1])){
    x0 <- matrix(data[i,]-mu, length(mu), 1)
    sigma <- sigma + (x0) %*% t(x0)
  }
  sigma <- sigma/(dim(data)[1])
  
  return(c(mu, diag(sigma), sigma[1,2]/sqrt(sigma[1,1]*sigma[2,2])))
}

#' Normalize and adjust data
#'
#' @param data Input data
#'
#' @return Normalized data
#' @keywords internal
datauni <- function(data){
  T <- length(data[[1]])
  Sigma <- sqrt(musigma(cbind(data[[1]][[T]], data[[2]][[T]]))[3:4])
  for(i in 1:T){
    data[[1]][[i]] <- data[[1]][[i]]/Sigma[1]*sqrt(0.0001)
    data[[2]][[i]] <- data[[2]][[i]]/Sigma[2]*sqrt(0.0001)
    data[[1]][[i]] <- data[[1]][[i]] - mean(data[[1]][[i]]) + sqrt(-2*log(0.05)*0.0001 + sum(1/(1:i))*0.1)
    data[[2]][[i]] <- data[[2]][[i]] - mean(data[[2]][[i]]) + sqrt(-2*log(0.05)*0.0001 + sum(1/(1:i))*0.1)
  }
  data <- cbind(unlist(data[[1]]), unlist(data[[2]]))
  return(data)
}