# setwd('/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Functions')
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
	if (any(n<0)){
		stop("The order of Hermite polynomial must be greater than or equal to 0.")
	}

	if (any(n-floor(n)!=0)){
		stop("The order of Hermite polynomial must be an integer")
	}
	#Call the hermite recursive function
	herm <- lapply(1:length(n), function(i) hermite_rec(n[i]))
	hmat <- matrix(0, nrow=nrow(x), ncol=length(n))
	#Evaluate the hermite polynomial function, given x
	for (i in 1:length(n)){
		if (length(herm[[i]])==1){
			h <- rep(herm[[i]], nrow(x))
		} else {
			h <- herm[[i]]
			y <- h[length(h)]*matrix(1, nrow=nrow(x), ncol=1)
			p <- 1
			for (j in (length(h)-1):1){
				y <- array(y, length(x[,i]))
				y <- y+h[j]*x[,i]^p
				p <- p+1
			}
			h <- as.matrix(y, nrow=nrow(x), ncol=1)
		}
		hmat[,i] <- 2^(-n[i]/2)*h
	}
	return(hmat)
}
#####################################################################################################
#Function that creates a tensor product from the hermite polynomials
#####################################################################################################
tensor.prod <- function(M, vars, norm){
	#Option to normalize data (useful for gradient descent)
	vars <- as.matrix(vars)
	if (norm==TRUE){
		vars <- (vars-colMeans(vars))/apply(vars, 2, sd)
	}
	#If we just want a univariate hermite polynomial
	if (length(M)==1){
		if (length(M)<ncol(vars)){
			print("Error: Number of vars should be equal to one")
		}
		prodlist <- sapply(0:M, function(z) hermite(z, vars))
		return(prodlist)
	#This covers the case when mapply cant return a list when all M is the same
	} else if (length(unique(M))==1){
		if (length(M)!=ncol(vars)){
			print('Error: M needs to equal to the number of vars')
		}
		vecs <- mapply(seq, 0, M)
		tmp <- split(vecs, rep(1:ncol(vecs), each=nrow(vecs)))
		tmp <- as.matrix(do.call(expand.grid, tmp))
		prodlist <- apply(tmp, 1, function(z) apply(hermite(z, x=vars), 1, prod))
		return(prodlist)
	#The case where we want a linear regression/random coefficient model
	} else if (M=='linear') {
		prodlist <- cbind(1, vecs)
		return(prodlist)	
	#For the multivariate hermite polynomial
	} else {
		if (length(M)!=ncol(vars)){
			print('Error: M needs to equal to the number of vars')
		}
		vecs <- mapply(seq, 0, M)
		tmp <- as.matrix(do.call(expand.grid, vecs))
		prodlist <- apply(tmp, 1, function(z) apply(hermite(z, x=vars), 1, prod))
		return(prodlist)
	}
}
####################################################################################################
####################################################################################################
####################################################################################################
translog <- function(K, L, M, omega){
	omega <- as.matrix(omega)
	#translog productivity component
	prodf <- cbind(1, K, L, M, K*L, L*M, K*M, K^2, L^2, M^2)
	#translog productivity component (first-order for now)
	prodw <- omega
	#Combine translog production polynomial with productivity component
	#Includes the constant but not the neutral productivity component which needs to be subtracted off the dependent variable
	tprod <- cbind(prodf, sweep(prodf[,-1], 1, prodw, `*`))
	return(tprod)
}









