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
	condata <- data.frame(idvar, X) %>% group_by(idvar) %>% slice(-1)
	lagdata <- data.frame(idvar, X) %>% group_by(idvar) %>% slice(-n())
	data <- cbind(data.frame(condata), data.frame(lagdata[,-1]))
	return(data)
}
######################################################################################################
#This function returns the first time period data
######################################################################################################
t0data <- function(idvar, X){
	data <- data.frame(idvar, X) %>% group_by(idvar) %>% slice(1) 
	return(data.frame(data))
}
######################################################################################################
#This function randomly samples a firm-level observation of size nbatch
#nbatch=1 draws one firm observation for stochastic gradient descent
#For n>nbatch>1 mini-batch stochastic gradient descent
#nbatch=n is for standard gradient descent
#Arguments are firm id, time period, observable data, unobservable data, nbatch, random seed, and
#a logical argument for the case where the dependent variable in the moment condition is an unobservable
#e.g. the specifciation for unobserved productivity. WT=2 is for the moment conditions where
#The dependent variable is observable, WT=1 is for when the dependent variable is the moment condition
#for the unobservable at t>1, WT=0 is for when the dependent variable is the moment condition for the unobservable
#at time t=1
#####################################################################################################
gradbatch <- function(idvar, X, U, nbatch, seed, WT){
	set.seed(seed)
	if (WT==2){
		data <- data.frame(cbind(idvar, X, U)) %>% subset(idvar%in%sample(unique(idvar), nbatch)) %>% select(-"idvar")
		names(data) <- c("Y", paste("X", 1:(ncol(X)-1), sep=""), paste("U", 1:ncol(U), sep=""))
	} else if (WT==1) {
		data <- lagdata(idvar, cbind(X, U)) %>% subset(idvar%in%sample(unique(idvar), nbatch)) %>% select(-"idvar")
		names(data) <- c(paste("X", 1:ncol(X), sep=""), paste("Y", 1:ncol(U), sep=""), paste("Xlag", 1:ncol(X), sep=""), paste("U", 1:ncol(U), sep=""))
		#Dropped Xlag since it is not used in estimation
		data <- data[,-grepl("Xlag", colnames(data))]
	} else if (WT==0) {
		data <- t0data(idvar, cbind(X, U)) %>% subset(idvar%in%sample(unique(idvar), nbatch)) %>% select(-"idvar")
		names(data) <- c(paste("X", 1:ncol(X), sep=""), paste("Y", 1:ncol(U), sep=""))
	}
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
###############################################################################################
###############################################################################################
###############################################################################################































