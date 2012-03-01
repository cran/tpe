tpe2d <- function(d,maxit=1e4,verbose=TRUE) {
  hc <- hclust(d,"single")
  n <- dim(hc$merge)[1]
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
      x[[i]] <- rbind(c(0,0),c(dmin,0))
    } else if (min(pair)<0) {
      c1 <- -min(pair)
      c2 <- ind[[max(pair)]]
      ind[[i]] <- c(c1,c2)
      if (length(c2)==2) {
        x[[i]] <- cmdscale(d[ind[[i]],ind[[i]]])
      } else {
        x[[i]] <- align2d(d[c1,c2,drop=FALSE],matrix(0,1,2),x[[max(pair)]],dmin,maxit)
      }
    } else {
      c1 <- ind[[pair[1]]]
      c2 <- ind[[pair[2]]]
      ind[[i]] <- c(c1,c2)
      x[[i]] <- align2d(d[c1,c2,drop=FALSE],x[[pair[1]]],x[[pair[2]]],dmin,maxit)
    }
  }
  y <- x[[n]][match(1:(n+1),ind[[n]]),]
  rownames(y) <- hc$labels
  y
}

align2d <- function(d,x1,x2,dmin,maxit) {
  n1 <- dim(x1)[1]
  n2 <- dim(x2)[1]
  nn <- which(d==dmin,arr.ind=TRUE)
  if (!is.null(dim(nn))) {
    indnn <- nn[1,]
  } else {
    indnn <- c(1,nn[1])
  }
  if (n1>3) {
    c1m <- mask2d(x1,dmin)
  } else {
    c1m <- 1:n1
  }
  if (n2>3) {
    c2m <- mask2d(x2,dmin)
  } else {
    c2m <- 1:n2
  }
  indmask <- which(d[c1m,c2m]<=sort(d[c1m,c2m])[2],arr.ind=TRUE)
  if (length(c1m)==1) {
    m1 <- 1
    m2 <- c2m[indmask]
  } else if (length(c2m)==1) {
    m1 <- c1m[indmask]
    m2 <- 1
  } else {
    m1 <- c1m[indmask[,1]]
    m2 <- c2m[indmask[,2]]
  }
  cx1 <- t(t(x1)-x1[indnn[1],])
  cx2 <- t(t(x2)-x2[indnn[2],])
  d1 <- as.matrix(dist(cx1))[indnn[1],]
  d2 <- as.matrix(dist(cx2))[indnn[2],]
  L <- function(par) {
    sigma2d(par,cx1[m1,,drop=FALSE],cx2[m2,,drop=FALSE],d[m1,m2])
  }
  P <- function(par) {
    rho2d(par,cx1,cx2,dmin)
  }
  par0 <- c(0,0,dmin+2*(max(d1)+max(d2)),0)
  par <- sumt(par0,L,P,method="Nelder-Mead",control=list(maxit=maxit))$x
  tx1 <- t(t(cx1%*%rot2d(par[1]))+c(par[3],par[4]))
  tx2 <- cx2%*%rot2d(par[2])
  rbind(tx1,tx2)
}

mask2d <- function(x,dmin) {
  out <- voronoi.mosaic(x[,1],x[,2],duplicate="remove")
  tri <- triangles(out$tri)
  mat <- rep(-1,dim(x)[1])
  for (i in 1:dim(tri)[1]) {
    for (j in 1:3) {
      mat[tri[i,j]] <- max(mat[tri[i,j]],out$radius[i])
    }
  }
  mat[convex.hull(out$tri)$i] <- Inf
  which(mat>dmin)
}

sigma2d <- function(par,x1,x2,d) {
  n1 <- dim(x1)[1]
  n2 <- dim(x2)[1]
  tx1 <- t(t(x1%*%rot2d(par[1]))+c(par[3],par[4]))
  tx2 <- x2%*%rot2d(par[2])
  td <- as.matrix(dist(rbind(tx1,tx2)))[1:n1,(n1+1):(n1+n2)]
  sum((td-d)^2)
}

rho2d <- function(par,x1,x2,dmin) {
  n1 <- dim(x1)[1]
  n2 <- dim(x2)[1]
  tx1 <- t(t(x1%*%rot2d(par[1]))+c(par[3],par[4]))
  tx2 <- x2%*%rot2d(par[2])
  tdmin <- min(as.matrix(dist(rbind(tx1,tx2)))[1:n1,(n1+1):(n1+n2)])
  if (dmin > tdmin) {
    Inf
  } else {
    (tdmin-dmin)^2
  }
}

rot2d <- function(theta) {
  rbind(c(cos(theta),sin(theta)),c(-sin(theta),cos(theta)))
}
