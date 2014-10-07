# Load required libraries ------------------

library(data.table)
library(igraph)
library(NetIndices)

# Load functions ---------------------------

ran.unif <- function(motmat, pred = 10, prey = -1, random = F){
  newmat <- apply(motmat, c(1,2), function(x){
    if(x==1){runif(1, 0, pred)}else if(x==-1){runif(1, prey, 0)} else{0}
  })
  if(random){
    diag(newmat) <- runif(length(diag(newmat)), -1, 0)
  }else{diag(newmat) <- -1}
  
  return(newmat)
}

maxRE <- function(rmat){
  lam.max <- eigen(rmat)$values[which.max(Re(eigen(rmat)$values))]
  return(lam.max)
}

eig.analysis <- function(n, matrices, params){
  require(data.table)
  dims <- dim(matrices[[1]])
  cols <- length(matrices)
  rows <- n
  eigenMATRIX.re <- matrix(nrow = rows, ncol = cols)
  eigenMATRIX.im <- matrix(nrow = rows, ncol = cols)
  samps <- list()
  for(i in 1:n){
    ranmat <- lapply(matrices, ran.unif, pred = params[,1],
                     prey = params[,2], random = T)
    sampvals <- matrix(nrow = length(ranmat), ncol = dims[1]^2)
    for(j in 1:length(ranmat)){
      sampvals[j,] <- ranmat[[j]]
    }
    eigs <- sapply(ranmat, maxRE)
    eigenMATRIX.re[i,] <- Re(eigs)
    eigenMATRIX.im[i,] <- Im(eigs)
    samps[[i]] <- as.data.frame(sampvals) 
  }
  svals <- cbind(web = rep(1:length(matrices), n), n = rep(1:n, each = cols), rbindlist(samps))
  return(list(samples = svals, ematrix.re = eigenMATRIX.re, ematrix.im = eigenMATRIX.im))
}

conversion <- function(tm){
  for(i in 1:nrow(tm)){
    for(j in 1:ncol(tm)){
      if(tm[i,j] == 1){tm[j,i] <- -1}
    }
  }
  return(tm)
}

randomWEBS <- function(S = 10, numweb = 200, chain = 9, total = 14){
  require(NetIndices)
  require(igraph)
  mywebs <- list()
  for(j in 1:numweb){
    
    check <- 1
    while(!check == 0){
      myweb <- matrix(0, nrow = S, ncol = S)
      for(i in 1:chain){
        myweb[i,i+1] <- 1
      }
      tophalf <- which(myweb[upper.tri(myweb)] == 0)
      newones <- sample(tophalf, total-chain)
      myweb[upper.tri(myweb)][newones] <- 1
      mywebs[[j]] <- myweb
      
      indeg <- apply(myweb, 1, sum)
      outdeg <- apply(myweb, 2, sum)
      deg <- indeg + outdeg
      
      if(sum(deg == 0) >= 1){check <- 1}else{check <- 0}
      
    }
    
  }
  return(mywebs)
}

randomQSS <- function(mywebs, params){
  require(NetIndices)
  require(igraph)
  
  mywebs1 <- lapply(mywebs, conversion)
  myweb.tl <- lapply(mywebs, TrophInd)
  emat <- eig.analysis(1000, mywebs1, params) 
  
  qss <- apply(emat$ematrix.re, 2, function(x){sum(x<0)/1000})
  maxtl <- sapply(myweb.tl, function(x){max(x$TL)})
  meantl <- sapply(myweb.tl, function(x){mean(x$TL)})
  medtl <- sapply(myweb.tl, function(x){median(x$TL)})
  sdtl <- sapply(myweb.tl, function(x){sd(x$TL)})
  diam <- sapply(lapply(mywebs, graph.adjacency), diameter)
  
  web.dat <- data.frame(qss, diam, maxtl, meantl, medtl, sdtl)
  iter.dat <- cbind(par = rep(paste(params, collapse = "_"), nrow(emat$samples)),
                    emat$samples, reals = as.vector(emat$ematrix.re),
                    im = as.vector(emat$ematrix.im))
  
  return(list(web.dat, iter.dat))
}



## PARAMS --------------------------------
pars <- data.frame(pred = c(10, 10, 10, 5, 5, 5, 1, 1, 1), prey = c(-1, -5, -10, -1, -5, -10, -1, -5, -10))

