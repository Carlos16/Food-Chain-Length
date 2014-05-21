## Functions to get eigenvalues of randomly sampled matrices

#' @title Fill a matrix with random uniform values
#' @description Will take any matrix and replace 1 with a random uniform value from 0:10 and replace -1 with
#'  a value drawn from a  uniform -1:0
#' @param motmat a matrix of 1s, 0s, and -1s to be filled
#' @param random whether the diagonal should be random (Default = TRUE) or set to -1
#' @seealso ran.lnorm
#' @details This function will take the matrix input and replace all ones with a value drawn from a 
#' uniform distribution between 0 and 10, replace all negative ones with a value drawn from a 
#' uniform distribution between  -1 and 0, and set the diagonal to either a random value drawn from a 
#' uniform distribution between -1 and 0 or just set to -1. If the parameter `random` is TRUE then the value 
#' will be drawn from the random uniform distribution, if FALSE, then the diagonal will be set to -1 for all.

ran.unif <- function(motmat, random = F){
  newmat <- apply(motmat, c(1,2), function(x){
    if(x==1){runif(1, 0, 10)}else if(x==-1){runif(1, -1, 0)} else{0}
  })
  if(random){
    diag(newmat) <- runif(length(diag(newmat)), -1, 0)
  }else{diag(newmat) <- -1}
  
  return(newmat)
}


#' @title Fill a matrix with random values from a lognormal distribution
#' @description Will take any matrix and replace 1 with a random uniform value from 0:10 and replace -1 with
#'  a value drawn from a  uniform -1:0
#' @param motmat a matrix of 1s, 0s, and -1s to be filled
#' @param random whether the diagonal should be random (Default = TRUE) or set to -1
#' @seealso ran.unif
#' @details This function will take the matrix input and replace all ones with a value drawn from a 
#' lognormal distribution (mean = -1, sd = 1), replace all negative ones with a value drawn from a 
#' lognormal distribution (mean = -5, sd = 1), and set the diagonal to a random value drawn from a 
#' lognormal distribution (mean = -5, sd = 1). 

ran.lnorm <- function(motmat){
  newmat <- apply(motmat, c(1,2), function(x){
    if(x==1){rlnorm(1, -1, 1)}else if(x==-1){-rlnorm(1, -5, 1)} else{0}
  })
  diag(newmat) <- -rlnorm(1, -5, 1)
  return(newmat)
}


#' @title Calculate the eigenvalue with the largest real part
#' @description Calculate the eigenvalues of a matrix, and return the largest real part.
#' @param rmat matrix from which eigenvalues will be calculated

maxRE <- function(rmat){
  lam.max <- max(Re(eigen(rmat)$values))
  return(lam.max)
}


#' @title Function to randomly fill a matrix and compute eigenvalue with largest real part
#' @description Function randomly fills a list of matrices and then computes the eigenvalue 
#' with the largest real part for any number of iterations
#' @param n the number of times to randomly fill each matrix
#' @param matrices a list of matrices
#' @param mode whether random values should come from a random uniform distribution or a lognormal; Defaults to "unif"
#' @details This is a wrapper function for ran.unif or ran.lnorm and maxRE. It allows eigenvalues to be 
#' computed on a list of matrices that are of interest.

eig.analysis <- function(n, matrices, mode = "unif"){
  cols <- length(matrices)
  rows <- n
  eigenMATRIX <- matrix(nrow = rows, ncol = cols)
  for(i in 1:n){
    if(mode == "unif"){
      ranmat <- lapply(matrices, ran.unif)
    }else if(mode == "lnorm"){
      ranmat <- lapply(matrices, ran.lnorm)
    }
    eigs <- sapply(ranmat, maxRE)
    eigenMATRIX[i,] <- eigs
  }
  return(eigenMATRIX)
}

#' @title Make a sign structure matrix
#' @description Converts an adjacency matrix into a sign structured predator prey matrix, 
#' where 1 indicates benefit and -1 indicates negative impadt
#' @param tm A binary adjacency matrix

conversion <- function(tm){
  for(i in 1:nrow(tm)){
    for(j in 1:ncol(tm)){
      if(tm[i,j] == 1){tm[j,i] <- -1}
    }
  }
  return(tm)
}


#' @title Compute quasi sign-stability of random food webs
#' @description Generates a set of random food webs and calculates quasi sign-stability for each. 
#' @param numweb the number of webs to generate
#' @param S how large the web should be (number of species)
#' @param chain how large should the basal food chain be
#' @param total total number of links in the web

randomQSS <- function(numweb = 200, S = 10, chain = 9, total = 14){
  # Required packages
  require(NetIndices)
  
  # Errors
  if(!chain < S){stop("Basal chain cannot be longer than number of species, "chain" must be smaller than S")}
  
  # Lopp to generate the list of webs
  mywebs <- list()
  for(j in 1:numweb){
    
    # While loop ensures that all species either eat another species or are eaten (fully connected web)
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
      # Maybe add an additional check to ensure that there is at least one species with indeg = 0
    }
    
  }
  
  #Convert the adjacency matrices into sign matrices
  mywebs1 <- lapply(mywebs, conversion)
  
  #Compute trophic level information 
  myweb.tl <- lapply(mywebs, TrophInd)
  
  #Compute eigenvalues (returns largest real part)
  emat <- eig.analysis(1000, mywebs1, mode = "unif")
  
  #Compute quasi sign-stability (proportion of randomizations when the web is stable)
  qss <- apply(emat, 2, function(x){sum(x<0)/1000})
  #Determine the largest trohic level
  mtl <- sapply(myweb.tl, max)
  
  # Returns a list of the webs generated, their trophic info, the eigenvalues and qss, and max trophic level
  return(list(webs = mywebs, trophics = myweb.tl, eigs = emat, qss = qss, maxTL = mtl))
}


#' @title Loop through different parameterizations of randomQSS
#' @description Run multiple iterations of randomQSS with different sizes for the basal chain
#' @param webiter number of random webs to generate
#' @param maxchain longest basal chain 
#' @param totalINT total number of interactions in the web
#' @param plot whether to generate a plot of QSS against maximum trophic level. Defaults to TRUE

testLENGTH <- function(webiter = 200, maxchain = 9, totalINT = 14, plot = TRUE){
  qssLIST <- list()
  mtlLIST <- list()
  troLIST <- list()
  for(i in 1:maxchain){
    cat(i, "\n")
    test <- randomQSS(numweb = webiter, chain = i, total = totalINT)
    qssLIST[[i]] <- test$qss
    mtlLIST[[i]] <- test$maxTL
    troLIST[[i]] <- lapply(test$trophics, function(x){x$TL})
  }
  
  quas <- unlist(qssLIST)
  maxtl <- unlist(mtlLIST)
  meantl<- rapply(spottest2$trolist, mean)
  medtl <- rapply(spottest2$trolist, median)
  
  data <- data.frame(QSS = quas, MaxTL = maxtl, MeanTL = meantl, MedTL = medtl)
  
  model <- glm(QSS~MaxTL, family = "quasibinomial", data = data)
  
  if(plot){
    require(ggplot2)
    g <- ggplot(cbind(data, fit = model$fitted.values), aes(x = MaxTL, y = QSS)) + geom_point()
    g <- g + geom_line(aes(x = MaxTL, y = fit), col = "blue")
    g <- g + geom_smooth(method = "loess", col = "orange")
    print(g)
  }
  return(list(data = data, model = model, trolist = troLIST))
}