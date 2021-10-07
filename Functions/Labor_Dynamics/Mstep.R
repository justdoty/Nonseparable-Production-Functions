# setwd('/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Functions')
# require(dplyr)
#####################################################################################################
#Functions to estimate the exponential tail parameters on (0, \tau_{1}) and (\tau_{L}, 1)
#####################################################################################################
expb <- function(YX, pred){
	q1 <- pred[,1]
	qL <- pred[,ncol(pred)]
	b1 <- -mean(YX<=q1)/mean((YX-q1)*(YX<=q1))
	bL <- mean(YX>=qL)/mean((YX-qL)*(YX>=qL))
	return(list(b1=b1, bL=bL))
}
obj <- function(YX, XX, par, tau){
	obj <- mean((tau-(YX-XX%*%par)<0)*(YX-XX%*%par))
	return(obj)
}
dobj <- function(YX, XX, tau, par){
	dobj <- colMeans(XX*repmat((((YX-XX%*%par)<0)-tau), 1, ncol(XX)))
	return(dobj)
}
