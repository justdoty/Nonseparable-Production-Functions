setwd('/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Functions')
source('/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Functions/Tensors.R')
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
###############################################################################################
#Modification of the Cubic Spline Function in PRACMA to compute deriviatives of spline
###############################################################################################
cspline <- function(x, y, xi = NULL, endp2nd = FALSE, der = c(0, 0), d.ord=0) {
    n <- length(x)
    h <- x[2:n] - x[1:(n-1)]
    e <- 2 * c(h[1], h[1:(n-2)] + h[2:(n-1)], h[n-1])
    A <- Diag(e) + Diag(h, -1) + Diag(h, 1)
    d <- (y[2:n] - y[1:(n-1)]) / h
    rhs <- 3* (d[2:(n-1)] - d[1:(n-2)])
    der0 <- der[1]; dern <- der[2]
    if (endp2nd) {
        A[1, 1] <- 2 * h[1];   A[1, 2] <- h[1]
        A[n, n] <- 2 * h[n-1]; A[n-1, n-2] <- h[n-1]
        rhs <- c(3*(d[1] - der0), rhs, 3*(dern - d[n-1]))
    } else {
        A[1, ] <- 0; A[1, 1] <- 1
        A[n, ] <- 0; A[n, n] <- 1
        rhs <- c(der0, rhs, dern)
    }
    S <- zeros(n, 4)
    S[, 3] <- solve(A, rhs)
    for (m in 1:(n-1)) {
        S[m,4] = (S[m+1,3]-S[m,3]) / 3 / h[m]
        S[m,2] = d[m] - h[m]/3 * (S[m + 1,3] + 2*S[m,3])
        S[m,1] = y[m]
    }
    S <- S[1:(n-1), 4:1]
    if (d.ord==0){
        pp <- mkpp(x, S)
    } else if (d.ord==1){
    	S <- sweep(S[,-ncol(S)], MARGIN=2, c(3,2,1), `*`)
    	pp <- mkpp(x, S)
    }
    
    if (is.null(xi)) {
        return(pp)
    } else {
        yi <- ppval(pp,xi)
        return(yi)
    }
}
###############################################################################################
#This function calculates the inverse of the quantile function p=Q^{-1}(x) using Bisection Method
###############################################################################################
qpoly <- function(PY, PX, vectau, par, pp, steps, tol){
	#Boundary for Bisection Algorithm
	a <- rep(vectau[1], nrow(PY)); b <- rep(vectau[length(vectau)], nrow(PY))
	#Function to Use for Root-Finding
	bifun <- function(p) PY-rowSums(PX*sapply(1:nrow(par), function(l) ppval(pp[[l]], p)))
	#Bisection Algorithm
	for (i in 1:steps){
		#Midpoint
		c <- (a+b)/2
		#Stop if conditions are met
		if ((all(bifun(c)==0))||(all((b-a)/2<tol))){
      	return(c)
    	}
    	#Adjust end points if conditions are not met
    	signind <- sign(bifun(c))==sign(bifun(a))
    	a[signind] <- c[signind]
    	b[!signind] <- c[!signind]
	}
}
densfun <- function(PY, PX, vectau, par, parb, pp, ppd, steps, tol){
	#Function to Interpolate Derivatives of Cubic Spline
	ppint <- function(xs) sapply(1:nrow(par), function(xp) ppval(ppd[[xp]], xs))
	#Density values
	D <- ((rowSums(PX*ppint(qpoly(PY=PY, PX=PX, vectau=vectau, par=par, pp=pp, steps=steps, tol=tol))))^-1)*(PX%*%par[,1]<=PY)*(PX%*%par[,length(vectau)]>PY)+
	(vectau[1]*parb[1]*exp(parb[1]*(PY-PX%*%par[,1])))*(PY<PX%*%par[,1])+
	((1-vectau[length(vectau)])*parb[2]*exp(-parb[2]*(PY-PX%*%par[,length(vectau)])))*((PX%*%par[,length(vectau)])<=PY)
	return(D)
}

























