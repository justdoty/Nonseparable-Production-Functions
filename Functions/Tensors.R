setwd('/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Functions')
##########################################################################################################
#This file creates the tensor product of the user choice of basis polynomials, Hermite polynomials
#are used here, working on implementing other specifications such as different basis functions, splines, wavelets, etc
##########################################################################################################
#This is the reccurence construction of a Hermite polynomial, i.e:
#H_{0}(x)=1
#H_{1}(x)=2x
#H_{n+1}(x)=2xH_{n}(x)-2nH[n-1](x)
#########################################################################################################
hermite_rec <- function(n){
	if (n==0){
		return(1)
	} else if (n==1) {
		return(c(2,0))
	} else {
		h1 <- array(0, c(1,n+1))
		h1[1:n] <- 2*hermite_rec(n-1)

		h2 <- array(0, c(1,n+1))
		h2[3:length(h2)] <- 2*(n-1)*hermite_rec(n-2)

		return(h1-h2)

	}
}
######################################################################################################
#Function that creates a univariate hermite polynomial
#####################################################################################################
hermite <- function(n, x){
	#Compute the Hermite polynomials
	x <- x/sqrt(2)
	#Check n
	if (n<0){
		stop("The order of Hermite polynomial must be greater than or equal to 0.")
	}

	if (n-floor(n)!=0){
		stop("The order of Hermite polynomial must be an integer")
	}
	#Call the hermite recursive function
	h <- hermite_rec(n)
	#Evaluate the hermite polynomial function, given x
	if (length(h)==1){
		h <- rep(h, nrow(x))
	} else {
		y <- h[length(h)]*matrix(1, nrow=nrow(x), ncol=1)
		p <- 1
		for (i in (length(h)-1):1){
			y <- y+h[i]*x^p
			p <- p+1
		}
		h <- as.matrix(y, nrow=nrow(x), ncol=1)
	}
	h <- 2^(-n/2)*h
	return(h)
}
#####################################################################################################
#Function that creates a tensor product from the hermite polynomials
#####################################################################################################
tensor.prod <- function(M, vars){
	if (length(M)==1){
		prodlist <- sapply(0:M, function(z) hermite(z, as.matrix((vars-mean(vars))/sd(vars))))
		return(prodlist)
	} else {
		polymat <- lapply(1:length(M), function(x) sapply(0:M[x], function(z) hermite(z, as.matrix((vars[,x]-mean(vars[,x]))/sd(vars[,x])))))
		prodlist <- list()
		prodlist[[1]] <- polymat[[1]]
		for (p in 2:(length(M))){
			prodlist[[p]] <- t(sapply(1:nrow(vars), function(j) tcrossprod(prodlist[[p-1]][j,], polymat[[p]][j,])))
		}
		return(prodlist[[length(M)]])
	}
}
####################################################################################################
####################################################################################################
####################################################################################################










