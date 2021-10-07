# setwd('/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Functions')
# source('/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Functions/Tensors.R')
# require(dplyr)
###################################################################################################
#This file contains auxillary files for data subsetting and other miscellaneous tasks
#####################################################################################################

#####################################################################################################
#This function lags the data and returns both the lagged variables as well as its contemporaneous values
#####################################################################################################
lagdata <- function(idvar, X){
	condata <- X[duplicated(idvar),]
	lagdata <- X[(1:length(idvar))[duplicated(idvar)]-1,]
  idvar <- idvar[duplicated(idvar)]
	data <- cbind(idvar, condata, lagdata)
	return(data)
}
######################################################################################################
#This function returns the first time period data
######################################################################################################
t1data <- function(idvar, X){
	t1data <- X[!duplicated(idvar),] 
  idvar <- idvar[!duplicated(idvar)]
  data <- cbind(idvar, t1data)
	return(data)
}
###############################################################################################
#This function computes the piece-wise linear spline with the laplace specification
#in the tail intervals in the intercept for the vector of coefficients
###############################################################################################
lspline <- function(vectau, bvec, b1, bL, u){
    #Laplace Tail Parameters
    qb1 <- (matrix(rep(bvec[,1], each=length(u)), nrow=length(u))+1/b1*(log(u/vectau[1]))*cbind(1, array(0, c(length(u), (nrow(bvec)-1)))))*(0<u)*(u<vectau[1])
    qbL <- (matrix(rep(bvec[,length(vectau)], each=length(u)), nrow=length(u))-1/bL*(log((1-u)/(1-vectau[length(vectau)])))*cbind(1, array(0, c(length(u), (nrow(bvec)-1)))))*(vectau[length(vectau)]<u)*(u<1)
    #Initialization
    qpar <- (matrix(rep(bvec[,1], each=length(u)), nrow=length(u))+as.matrix((u-vectau[1])/(vectau[2]-vectau[1]))%*%t((bvec[,2]-bvec[,1])))*(vectau[1]<=u)*(u<vectau[2])
    for (q in 2:(length(vectau)-1)){
        qpar <- qpar+(matrix(rep(bvec[,q], each=length(u)), nrow=length(u))++as.matrix((u-vectau[q])/(vectau[q+1]-vectau[q]))%*%t((bvec[,q+1]-bvec[,q])))*(vectau[q]<=u)*(u<vectau[q+1])

    }
    lspline <- qb1+qpar+qbL
    return(lspline)
}
wlspline <- function(vectau, bvec, b1, bL, scoef, u){
    #Laplace Tail Parameters
    qb1 <- (matrix(rep(bvec[,1], each=length(u)), nrow=length(u))+scoef[1]+1/b1*(log(u/vectau[1]))*cbind(1, array(0, c(length(u), (nrow(bvec)-1)))))*(0<u)*(u<vectau[1])
    qbL <- (matrix(rep(bvec[,length(vectau)], each=length(u)), nrow=length(u))+scoef[2]-1/bL*(log((1-u)/(1-vectau[length(vectau)])))*cbind(1, array(0, c(length(u), (nrow(bvec)-1)))))*(vectau[length(vectau)]<u)*(u<1)
    #Initialization
    qpar <- (matrix(rep(bvec[,1], each=length(u)), nrow=length(u))+as.matrix((u-vectau[1])/(vectau[2]-vectau[1]))%*%t((bvec[,2]-bvec[,1])))*(vectau[1]<=u)*(u<vectau[2])
    for (q in 2:(length(vectau)-1)){
        qpar <- qpar+(matrix(rep(bvec[,q], each=length(u)), nrow=length(u))++as.matrix((u-vectau[q])/(vectau[q+1]-vectau[q]))%*%t((bvec[,q+1]-bvec[,q])))*(vectau[q]<=u)*(u<vectau[q+1])

    }
    lspline <- qb1+qpar+qbL
    return(lspline)
}
###############################################################################################
###############################################################################################
###############################################################################################

########################################################################################################
#boot resampling on IDs: bootstrapping on individuals (prodest.R)
#######################################################################################################
block.boot.resample <- function(idvar, R){
  unique.ids <- unique(idvar) # find the unique values of panels in order to reshape the data
  panel.time.indices <- apply(unique.ids, 1, function(x) {return(list(which(idvar == x)))}) # find the time indices for each panel
  seq.indices <- 1:length(unique.ids) # the panel.time.indices list is indexed with sequential numbers: we mimic it
  boot.panel.id <- replicate(R, sample(seq.indices, replace = TRUE)) # generate a matrix of new IDs - R times
  new.indices <- list() # generate the matrix of the new indices
  ind <- 1:length(unique.ids)
  for (r in 1:R){ # for each boot rep we generate a vector of indices with rownames equal to a new - and fake - ID
    new.indices[[r]] <- cbind(unlist(mapply(function(x,y) {
      names(panel.time.indices[[x]][[1]]) <- rep(y,length(panel.time.indices[[x]][[1]]))
      return(list(panel.time.indices[[x]][[1]]))
    }, boot.panel.id[,r], ind))) # return a fake ID (sequential number) as row name and the index referring to the true ID
  }
  return(new.indices)
}
#Function for creating the covariance matrix of an AR(1) process
ar1_cor <- function(n, rho){
exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - 
    (1:n - 1))
rho^exponent
}
#Better function for drawing multiple multivariate samples with subsetted covariance matrices
mvdraw <- function(N, nsize, tsize, sig){
  mv <- mvrnorm(N, mu=rep(0, max(tsize)), Sigma=sig)
  x <- list()
  for (i in 1:N){
    x[[i]] <- mv[i,1:tsize[i]]
  }
  return(unlist(x))
}









