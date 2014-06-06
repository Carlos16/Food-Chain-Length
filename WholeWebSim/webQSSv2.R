randomWEBS <- function(numweb = 200, chain = 9, total = 14){
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
  
  return(webs = mywebs)
}

getQSS <- function(matrices, params){
  mywebs1 <- lapply(matrices, conversion)
  myweb.tl <- lapply(matrices, TrophInd)
  emat <- eig.analysis(1000, mywebs1, params, mode = "unif")
  qss <- apply(emat, 2, function(x){sum(x<0)/1000})
  maxtl <- sapply(myweb.tl, function(x){max(x$TL)})
  meantl <- sapply(myweb.tl, function(x){mean(x$TL)})
  medtl <- sapply(myweb.tl, function(x){median(x$TL)})
  
  data <- data.frame(qss, maxtl, meantl, medtl)
  
  return(data)
}




pars <- data.frame(pred = c(10, 10, 10, 5, 5, 5, 1, 1, 1), prey = c(-1, -5, -10, -1, -5, -10, -1, -5, -10))
ints <- c(12, 16, 20, 24, 28)


webLIST <- function(n, int, maxchain = 9){
  r <- list()
  for(i in 1:maxchain){
    r[[i]] <- randomWEBS(n, i, int)
  }
  return(r)
}

loopars <- function(matlist, pars){
  require(data.table)
  for(i in 1:nrow(pars)){
    d <- lapply(matlist, getQSS, params = pars[i,])
    d2[[i]] <- rbindlist(d)
    d2[[i]]$scenario <- factor(rep(paste(pars[i,], collapse = "/"), nrow(d2[[i]])))
  }
  data <- rbindlist(d2)
  return(data)
}

QSSwrapper <- function(n, numints, params){
  data <- list()
  for(i in 1:length(numints)){
    wl <- webLIST(n, numints[i], maxchain = 9)
    data[[i]] <- loopars(wl, pars)
    data[[i]]$ints <- factor(rep(numints[i]), nrow(data[[i]]))
    cat(i/length(numints)*100, "% done", "\n")
  }
  return(rbindlist(data))
}

webQSS<- QSSwrapper(n = 200, numints = ints, params = pars)

require(ggplot2)
ggplot(d2, aes(x = maxtl, y = qss)) + geom_point() + geom_smooth(method = "glm")


gmed <- ggplot(webQSS, aes(x = MedTL, y = QSS)) + geom_point(aes(col = ints)) 
gmed <- gmed + geom_smooth(aes(col = ints), method = "glm") + facet_wrap(~scenario)

gmean <- ggplot(webQSS, aes(x = MeanTL, y = QSS)) + geom_point(aes(col = ints)) 
gmean <- gmean + geom_smooth(aes(col = ints), method = "glm") + facet_wrap(~scenario)

gmax <- ggplot(webQSS, aes(x = MaxTL, y = QSS)) + geom_point(aes(col = ints)) 
gmax <- gmax + geom_smooth(aes(col = ints), method = "glm") + facet_wrap(~scenario) 

#ggsave("medtlPLOT2.png", width = 9, height = 7, dpi = 600)

setwd("C:/Users/borre_000/Desktop/GitHub/Food-Chain-Length/WholeWebSim/")

save.image("webQSSdata.Rdata")