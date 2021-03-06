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
  
  for(i in 1:n){
    ranmat <- lapply(matrices, ran.unif, pred = params[,1],
                     prey = params[,2], random = F)
    
    eigs <- sapply(ranmat, maxRE)
    eigenMATRIX.re[i,] <- Re(eigs)
    eigenMATRIX.im[i,] <- Im(eigs)
    
  }
  
  return(list(ematrix.re = eigenMATRIX.re, ematrix.im = eigenMATRIX.im))
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
  
  
  return(web.dat)
}

niche.model<-function(S,C){
  require(igraph)
  connected = FALSE
  while(!connected){  
    new.mat<-matrix(0,nrow=S,ncol=S)
    ci<-vector()
    niche<-runif(S,0,1)
    r<-rbeta(S,1,((1/(2*C))-1))*niche
    
    for(i in 1:S){
      ci[i]<-runif(1,r[i]/2,niche[i])
    }
    
    r[which(niche==min(niche))]<-.00000001
    
    for(i in 1:S){
      
      for(j in 1:S){
        if(niche[j]>(ci[i]-(.5*r[i])) && niche[j]<(ci[i]+.5*r[i])){
          new.mat[j,i]<-1
        }
      }
    }
    
    new.mat<-new.mat[,order(apply(new.mat,2,sum))]
    
    connected <- is.connected(graph.adjacency(new.mat))
  }
  return(new.mat)
}


niche_maker <- function(n, S, C){
  niche.list <- list()
  for (i in 1:n){
    niche.list[[i]]<- niche.model(S, C)
  }
  return(niche.list)
}

webtest <- randomWEBS(S = 100, numweb = 25, chain = 10, total = 1000)

qsstest <- randomQSS(webtest, matrix(c(1, -.01, 1, -.1), nrow = 2, byrow = T))

nichetest <- niche_maker(50, 100, .05)
nichetest2 <- niche_maker(50, 100, .1)
nichetest3 <- niche_maker(50, 50, .1)

nicheqss <- randomQSS(nichetest, matrix(c(1, -1), nrow = 1))
nicheqss2 <- randomQSS(nichetest2, matrix(c(1, -1), nrow = 1))
nicheqss3 <- randomQSS(nichetest3, matrix(c(1, -1), nrow = 1))

nicheqss.1 <- randomQSS(nichetest, matrix(c(3, -1), nrow = 1))
nicheqss2.1 <- randomQSS(nichetest2, matrix(c(3, -1), nrow = 1))
nicheqss3.1 <- randomQSS(nichetest3, matrix(c(3, -1), nrow = 1))
