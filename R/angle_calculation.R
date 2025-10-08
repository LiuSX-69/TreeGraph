#' Calculate angles for tree graph analysis
#'
#' This is the core function that calculates various angle statistics
#' for points relative to a starting point.
#'
#' @param data A matrix or data frame with two columns representing coordinates
#' @param alpha The alpha parameter (1 - coverage)
#' @param start The starting point coordinates (x, y)
#'
#' @return A list containing various angle statistics
#' @keywords internal
anglealpha <- function(data, alpha, start){
  # 函数内容保持不变（你原来的完整代码）
  # 这里放置你原来的anglealpha函数的完整代码
  data <- data[((data[,1]>start[1])&(data[,2]>start[2])),]
  
  outnum <- floor(dim(data)[1]*(1-alpha))
  
  ab <- function(point, start){
    a <- (point[2]-start[2])/(point[1]-start[1])
    b <- (start[2]*point[1]-point[2]*start[1])/(point[1]-start[1])
    return(c(a,b))
  }
  
  abangle <- function(a, b, point, start){
    pointtrans <- point - start
    
    if(a >= 0){
      angle <- atan(a)/pi*180
    } else if((a < 0) & (pointtrans[1] < 0)){
      angle <- atan(a)/pi*180 + 180
    } else {
      angle <- atan(a)/pi*180
    }
    return(angle)
  }
  
  Tan <- c()
  Angle <- c()
  Coef <- c()
  
  for(i in 1:dim(data)[1]){
    point <- data[i,]
    coef <- ab(point, start)
    Tan <- c(Tan, coef[1])
    Coef <- rbind(Coef, coef)
    
    angle <- abangle(coef[1], coef[2], point, start)
    Angle <- c(Angle, angle)
  }
  
  comple <- function(outnum, Angle, orderup, orderdown){
    updown <- c()
    anglevs <- c()
    for(i in 1:(outnum+1)){
      up <- orderup[i]
      down <- orderdown[outnum+2-i]
      updown <- rbind(updown, c(up, down))
      anglevs <- c(anglevs, abs(Angle[up]-Angle[down]))
    }
    return(list(updown, anglevs))
  }
  
  orderup <- order(-Angle)
  orderdown <- order(Angle)
  
  if(outnum == 0){
    whereup <- orderup[1]
    wheredown <- orderdown[1]
    anglemin <- abs(Angle[whereup]-Angle[wheredown])
    tanup <- Tan[whereup]
    tandown <- Tan[wheredown]
    tan2 <- (tanup-tandown)/(1+tanup*tandown)
    pointup <- data[whereup,]
    coefup <- Coef[whereup,]
    pointdown <- data[wheredown,]
    coefdown <- Coef[wheredown,]
    
    angleminlist <- rbind(c(anglemin, tan2), pointup, pointdown, coefup, coefdown)
    angleequallist <- angleminlist
    anglemeanlist <- angleminlist
    angleneighlist <- angleminlist
    anglemedianlist <- angleminlist
  } else {
    updownanglevs <- comple(outnum, Angle, orderup, orderdown)
    updown <- updownanglevs[[1]]
    anglevs <- updownanglevs[[2]]
    
    where <- order(anglevs)[1]
    whereup <- updown[where,1]
    wheredown <- updown[where,2]
    anglemin <- anglevs[where]
    tanup <- Tan[whereup]
    tandown <- Tan[wheredown]
    tan2 <- (tanup-tandown)/(1+tanup*tandown)
    pointup <- data[whereup,]
    coefup <- Coef[whereup,]
    pointdown <- data[wheredown,]
    coefdown <- Coef[wheredown,]
    angleminlist <- rbind(c(anglemin, tan2), pointup, pointdown, coefup, coefdown)
    
    if(outnum == 1){
      angleequallist <- angleminlist
    } else if(outnum %% 2 == 0){
      where <- ceiling(outnum/2)
      whereup <- updown[where,1]
      wheredown <- updown[where,2]
      angleequal <- anglevs[where]
      tanup <- Tan[whereup]
      tandown <- Tan[wheredown]
      tan2 <- (tanup-tandown)/(1+tanup*tandown)
      pointup <- data[whereup,]
      coefup <- Coef[whereup,]
      pointdown <- data[wheredown,]
      coefdown <- Coef[wheredown,]
      angleequallist <- rbind(c(angleequal, tan2), pointup, pointdown, coefup, coefdown)
    } else {
      where <- floor(outnum/2)
      whereup <- updown[where,1]
      whereupup <- updown[where+1,1]
      wheredown <- updown[where,2]
      wheredowndown <- updown[where+1,2]
      pointup <- apply(data[c(whereup, whereupup),], 2, mean)
      coefup <- ab(pointup, start)
      pointdown <- apply(data[c(wheredown, wheredowndown),], 2, mean)
      coefdown <- ab(pointdown, start)
      angleequal <- abs(mean(Angle[c(whereup, whereupup)])-mean(Angle[c(wheredown, wheredowndown)]))
      tan2 <- tan(angleequal)
      angleequallist <- rbind(c(angleequal, tan2), pointup, pointdown, coefup, coefdown)
    }
    
    anglemean <- mean(anglevs)
    tan2 <- tan(anglemean)
    pointup <- apply(data[updown[,1],], 2, mean)
    coefup <- ab(pointup, start)
    pointdown <- apply(data[updown[,2],], 2, mean)
    coefdown <- ab(pointdown, start)
    anglemeanlist <- rbind(c(anglemean, tan2), pointup, pointdown, coefup, coefdown)
    
    anglemedian <- median(anglevs)
    tan2 <- tan(anglemedian)
    pointup <- apply(data[updown[,1],], 2, median)
    coefup <- ab(pointup, start)
    pointdown <- apply(data[updown[,2],], 2, median)
    coefdown <- ab(pointdown, start)
    anglemedianlist <- rbind(c(anglemedian, tan2), pointup, pointdown, coefup, coefdown)
    
    where <- order(anglevs)[1]
    whereup <- updown[where,1]
    upnum <- which(orderup == whereup)
    whereupup <- orderup[ifelse(upnum==1, upnum, upnum-1)]
    wheredown <- updown[where,2]
    downnum <- which(orderdown == wheredown)
    wheredowndown <- orderdown[ifelse(downnum==1, downnum, downnum-1)]
    pointup <- apply(data[c(whereup, whereupup),], 2, mean)
    coefup <- ab(pointup, start)
    pointdown <- apply(data[c(wheredown, wheredowndown),], 2, mean)
    coefdown <- ab(pointdown, start)
    angleneigh <- abs(mean(Angle[c(whereup, whereupup)])-mean(Angle[c(wheredown, wheredowndown)]))
    tan2 <- tan(angleneigh)
    angleneighlist <- rbind(c(angleneigh, tan2), pointup, pointdown, coefup, coefdown)
  }
  
  out <- list(
    anglemin = angleminlist,
    anglemedian = anglemedianlist,
    angleequal = angleequallist,
    anglemean = anglemeanlist,
    angleneigh = angleneighlist
  )
  
  return(out)
}