
## Functions to get eigenvalues of randomly sampled matrices
ran.unif <- function(motmat, pred = 10, prey = -1, random = F){
  newmat <- apply(motmat, c(1,2), function(x){
    if(x==1){runif(1, 0, pred)}else if(x==-1){runif(1, prey, 0)} else{0}
  })
  if(random){
    diag(newmat) <- runif(length(diag(newmat)), -1, 0)
  }else{diag(newmat) <- -1}
  
  return(newmat)
}

### using a lognormal
ran.lnorm <- function(motmat){
  newmat <- apply(motmat, c(1,2), function(x){
    if(x==1){rlnorm(1, -1, 1)}else if(x==-1){-rlnorm(1, -5, 1)} else{0}
  })
  diag(newmat) <- -rlnorm(1, -5, 1)
  return(newmat)
}

### Calculate largest eigenvalue (the real part)
maxRE <- function(rmat){
  lam.max <- max(Re(eigen(rmat)$values))
  return(lam.max)
}

### Wrap previous two functions together
eig.analysis <- function(n, matrices, params, mode = "unif"){
  cols <- length(matrices)
  rows <- n
  eigenMATRIX <- matrix(nrow = rows, ncol = cols)
  for(i in 1:n){
    if(mode == "unif"){
      ranmat <- lapply(matrices, ran.unif, pred = params[,1], prey = params[,2], random = T)
    }else if(mode == "lnorm"){
      ranmat <- lapply(matrices, ran.lnorm, pred = params[,1], prey = params[,2])
    }
    eigs <- sapply(ranmat, maxRE)
    eigenMATRIX[i,] <- eigs
  }
  return(eigenMATRIX)
}

conversion <- function(tm){
  for(i in 1:nrow(tm)){
    for(j in 1:ncol(tm)){
      if(tm[i,j] == 1){tm[j,i] <- -1}
    }
  }
  return(tm)
}


# -----------------------
randomQSS <- function(numweb = 200, chain = 9, total = 14, params){
  require(NetIndices)
  mywebs <- list()
  for(j in 1:numweb){
    
    check <- 1
    while(!check == 0){
      myweb <- matrix(0, nrow = 10, ncol = 10)
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
      # Maybe add an additional check to ensure that there is at least one species with indeg = 0
    }
    
  }
  
  mywebs1 <- lapply(mywebs, conversion)
  myweb.tl <- lapply(mywebs, TrophInd)
  emat <- eig.analysis(1000, mywebs1, params, mode = "unif")
  qss <- apply(emat, 2, function(x){sum(x<0)/1000})
  mtl <- sapply(myweb.tl, max)
  
  return(list(webs = mywebs, trophics = myweb.tl, eigs = emat, qss = qss, maxTL = mtl))
}

testLENGTH <- function(webiter = 200, maxchain = 9, totalINT = 14, params){
  qssLIST <- list()
  mtlLIST <- list()
  troLIST <- list()
  for(i in 1:maxchain){
    cat(i, "\n")
    test <- randomQSS(numweb = webiter, chain = i, total = totalINT, params = params)
    qssLIST[[i]] <- test$qss
    mtlLIST[[i]] <- test$maxTL
    troLIST[[i]] <- lapply(test$trophics, function(x){x$TL})
  }
  quas <- unlist(qssLIST)
  maxtl <- unlist(mtlLIST)
  meantl<- rapply(troLIST, mean)
  medtl <- rapply(troLIST, median)
  
  data <- data.frame(QSS = quas, MaxTL = maxtl, MeanTL = meantl, MedTL = medtl)
  
  return(data)
}


pars <- data.frame(pred = c(10, 10, 10, 5, 5, 5, 1, 1, 1), prey = c(-1, -5, -10, -1, -5, -10, -1, -5, -10))
ints <- c(12, 16, 20, 24, 28)


#totalDATA <- list()
for(i in 7:nrow(pars)){
  chainDATA <- data.frame(QSS = c(), MaxTL = c(), MeanTL = c(), MedTL = c(), ints = factor())
  for(j in 1:length(ints)){
    chain <- testLENGTH(webiter = 100, maxchain = 9, totalINT = ints[j], params = pars[i,])
    chain <- cbind(chain, ints = factor(rep(ints[j], 900)))
    chainDATA <- rbind(chainDATA, chain)
    cat("\n", ints[j], "is done", "\n")
  }
  totalDATA[[i]] <- chainDATA
  cat("\n", i, "th par done", "\n")
  # each item of the list corresponds to the row of parameters
}

totalDATA2 <- list()
for(i in 1:9){
  totalDATA2[[i]] <- cbind(totalDATA[[i]], scenario = paste(pars[i,], collapse = "/"))
}

head(totalDATA2[[1]])

totalDAT <- do.call(rbind, totalDATA2)
require(ggplot2)
ggplot(totalDAT, aes(x = MaxTL, y = QSS)) + geom_point(aes(col = ints)) + geom_smooth(aes(col = ints), method = "glm") + facet_wrap(~scenario)
ggplot(totalDAT, aes(x = MeanTL, y = QSS)) + geom_point(aes(col = ints)) + geom_smooth(aes(col = ints), method = "glm") + facet_wrap(~scenario)
ggplot(totalDAT, aes(x = MedTL, y = QSS)) + geom_point(aes(col = ints)) + geom_smooth(aes(col = ints), method = "glm") + facet_wrap(~scenario)

getwd()
save.image("chainINFO.Rdata")
