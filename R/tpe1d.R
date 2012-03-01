tpe1d <- function(d,verbose=TRUE) {
  hc <- hclust(d,"single")
  n <- length(hc$height)
  d <- as.matrix(d)
  ind <- list()
  x <- list()
  for (i in 1:n) {
    if (verbose) {
      cat("Iteration",i,"of",n,"\n")
    }
    dmin <- hc$height[i]
    pair <- hc$merge[i,]
    if (max(pair)<0) {
      ind[[i]] <- -pair
      x[[i]] <- rbind(0,dmin)
    } else if (min(pair)<0) {
      c1 <- -min(pair)
      c2 <- ind[[max(pair)]]
      ind[[i]] <- c(c1,c2)
      x[[i]] <- align1d(d[c1,c2,drop=FALSE],matrix(0,1,1),x[[max(pair)]],dmin)
    } else {
      c1 <- ind[[pair[1]]]
      c2 <- ind[[pair[2]]]
      ind[[i]] <- c(c1,c2)
      x[[i]] <- align1d(d[c1,c2,drop=FALSE],x[[pair[1]]],x[[pair[2]]],dmin)
    }
  }
  y <- x[[n]][match(1:(n+1),ind[[n]]),,drop=FALSE]
  rownames(y) <- hc$labels
  y
}

align1d <- function(d,x1,x2,dmin) {
  n1 <- dim(x1)[1]
  n2 <- dim(x2)[1]
  sigma <- function(x1,x2) {
    td <- as.matrix(dist(rbind(x1,x2)))[1:n1,(n1+1):(n1+n2),drop=FALSE]
    sum((d-td)^2)
  }
  x1r <- reflect1d(x1)
  x2r <- reflect1d(x2)
  t1 <- min(x2)-max(x1)-dmin
  t2 <- max(x2)-min(x1)+dmin
  val <- c(sigma(x1+t1,x2),
           sigma(x1+t2,x2),
           sigma(x1r+t1,x2),
           sigma(x1r+t2,x2),
           sigma(x1+t1,x2r),
           sigma(x1+t2,x2r),
           sigma(x1r+t1,x2r),
           sigma(x1r+t2,x2r))
  switch(which.min(val),
         rbind(x1+t1,x2),
         rbind(x1+t2,x2),
         rbind(x1r+t1,x2),
         rbind(x1r+t2,x2),
         rbind(x1+t1,x2r),
         rbind(x1+t2,x2r),
         rbind(x1r+t1,x2r),
         rbind(x1r+t2,x2r))
}

reflect1d <- function(x) {
  max(x)-x+min(x)
}
