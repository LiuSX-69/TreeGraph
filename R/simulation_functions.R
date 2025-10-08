#' Calculate delta angle for different rho values
#'
#' @param rho Correlation coefficient
#' @param mu Mean vector
#' @param sigma Variance vector
#' @param T Depth
#' @param N Branching factor
#' @param Tij Time parameter
#' @param alpha Alpha parameter (1 - coverage)
#' @param start Starting point
#'
#' @return Delta angle matrix
#' @keywords internal
deltaanglefunc <- function(rho, mu, sigma, T, N, Tij, alpha, start){
  data <- treedataresidual(mu, sigma, rho, T, N, Tij)
  data <- list(without = cbind(unlist(data[[1]]), unlist(data[[2]])), 
               with = datauni(data))
  deltaangle <- c()
  for(i in 1:2){
    out <- anglealpha(data[[i]], alpha, start)
    deltaangle <- rbind(deltaangle, unlist(lapply(out, function(out0){return(out0[1,1])})))
  }
  return(deltaangle)
}

#' Large-scale simulation function
#'
#' @param rho1 First correlation coefficient
#' @param rho2 Second correlation coefficient
#' @param mu Mean vector
#' @param sigma Variance vector
#' @param T Depth
#' @param N Branching factor
#' @param Tij Time parameter
#' @param n2 Number of repetitions
#' @param alpha Alpha parameter
#' @param start Starting point
#'
#' @return Comparison results
#' @keywords internal
largefunc_diff <- function(rho1, rho2, mu, sigma, T, N, Tij, n2, alpha, start){
  deltaangle1 <- vector('list', 2)
  deltaangle2 <- vector('list', 2)
  for(i in 1:n2){
    deltaangle10 <- deltaanglefunc(rho1, mu[1:2], sigma[1:2], T, N, Tij, alpha, start)
    deltaangle20 <- deltaanglefunc(rho2, mu[3:4], sigma[3:4], T, N, Tij, alpha, start)
    for(j in 1:2){
      deltaangle1[[j]] <- rbind(deltaangle1[[j]], deltaangle10[j,])
      deltaangle2[[j]] <- rbind(deltaangle2[[j]], deltaangle20[j,])
    }
  }
  out <- c()
  for(j in 1:2){
    out <- rbind(out, unlist(lapply(1:5, function(i){
      return(sum(deltaangle1[[j]][,i] >= deltaangle2[[j]][,i])/n2)
    })))
  }
  return(out)
}

largefunc_same <- function(rho1, rho2, mu, sigma, T, N, Tij, n2, alpha, start){
  deltaangle1 <- vector('list', 2)
  deltaangle2 <- vector('list', 2)
  for(i in 1:n2){
    deltaangle10 <- deltaanglefunc(rho1, mu, sigma, T, N, Tij, alpha, start)
    deltaangle20 <- deltaanglefunc(rho2, mu, sigma, T, N, Tij, alpha, start)
    for(j in 1:2){
      deltaangle1[[j]] <- rbind(deltaangle1[[j]], deltaangle10[j,])
      deltaangle2[[j]] <- rbind(deltaangle2[[j]], deltaangle20[j,])
    }
  }
  out <- c()
  for(j in 1:2){
    out <- rbind(out, unlist(lapply(1:5, function(i){
      return(sum(deltaangle1[[j]][,i] >= deltaangle2[[j]][,i])/n2)
    })))
  }
  return(out)
}

#' Compare two correlation coefficients
#'
#' @param rho1 First correlation coefficient
#' @param rho2 Second correlation coefficient
#' @param T Depth
#' @param N Branching factor
#' @param Tij Time parameter
#' @param n1 Outer repetitions
#' @param n2 Inner repetitions
#' @param alpha Alpha parameter
#' @param start Starting point
#'
#' @return Comparison results
#' @keywords internal
rho12func_diff <- function(rho1, rho2, T, N, Tij, n1, n2, alpha, start){
  out <- vector('list', 2)
  for(i in 1:n1){
    sigma <- round(1/rgamma(4, 50, 100), 2)/100
    mu <- sqrt((sigma)*6) + 0.01 + c(1, 1, 0, 0)
    out0 <- largefunc_diff(rho1, rho2, mu, sigma, T, N, Tij, n2, alpha, start)
    for(j in 1:2){
      out[[j]] <- rbind(out[[j]], out0[j,])
    }
  }
  return(out)
}

rho12func_same <- function(rho1, rho2, T, N, Tij, n1, n2, alpha, start){
  out <- vector('list', 2)
  for(i in 1:n1){
    sigma <- round(1/rgamma(2, 50, 100), 2)/100
    mu <- sqrt(sigma*6) + 0.01 + 7
    out0 <- largefunc_same(rho1, rho2, mu, sigma, T, N, Tij, n2, alpha, start)
    for(j in 1:2){
      out[[j]] <- rbind(out[[j]], out0[j,])
    }
  }
  return(out)
}