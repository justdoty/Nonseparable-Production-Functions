setwd('/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Functions')
source('Tensors.R')
#####################################################################################################
#This file contains the functions that are used to evaluate the gradient descent algorithms
#Requires the file 'Tensors.R' to evaluate the tensor products
#####################################################################################################

#####################################################################################################
#These are the objective function and its jacobians for output, labor, and materials and productivity
#####################################################################################################
qmoment <- function(Y, X, U, MH, c, tau, WT){
	if (WT==2){
		Xc <- apply(U, 2, function(u) tensor.prod(MH, cbind(X, u))%*%c)
		obj <- rowMeans((Y-Xc)*(tau-(Y<Xc)))
	} else if (WT==1){
		Xc <- apply(U, 2, function(u) tensor.prod(MH, cbind(X, u))%*%c)
		obj <- rowMeans((Y-Xc)*(tau-(Y<Xc)))
	} else if (WT==0){
		Xc <- tensor.prod(MH, as.matrix(X))%*%c
		obj <- rowMeans((Y-Xc)*(tau-(Y<Xc)))
	}
	return(obj)
}
qjacobian <- function(Y, X, U, MH, c, tau, WT){
	if (WT==2){
		Xc <- apply(U, 2, function(u) tensor.prod(MH, cbind(X, u))%*%c)
		jac <- rowMeans(sapply(1:ncol(U), function(u) t(tau-(Y<Xc))[u,]%*%tensor.prod(MH, cbind(X, U[,u])))/length(Y))
	} else if (WT==1){
		Xc <- apply(U, 2, function(u) tensor.prod(MH, cbind(X, u))%*%c)
		jac <- rowMeans(sapply(1:ncol(U), function(u) t(tau-(Y<Xc))[u,]%*%tensor.prod(MH, cbind(X, U[,u])))/length(Y))
	} else if (WT==0){
		X <- tensor.prod(MH, as.matrix(X))
		jac <- colMeans((t(tau-(Y-X%*%c))%*%X)/nrow(Y))
	}
	return(jac)
}
###################################################################################################
#This is the objective function and its jacobian for the investment decision rule
###################################################################################################
emoment <- function(Y, X, U, MH, c){
	obj <- rowMeans(apply(U, 2, function(u) (Y-tensor.prod(MH, cbind(X, u))%*%c)))
	return(obj)
}
ejacobian <- function(Y, X, U, MH, c){
	jac <- rowMeans(apply(U, 2, function(u) t(-2*(Y-tensor.prod(MH, cbind(X, u))%*%c))%*%tensor.prod(MH, cbind(X, u)))/length(Y))
	return(jac)
}
###################################################################################################
###################################################################################################
###################################################################################################


















