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
	if (norm==TRUE){
		vars <- (as.matrix(vars)-colMeans(vars))/apply(vars, 2, sd)
	} else {
		vars <- as.matrix(vars)
	}
	if (length(M)==1){
		prodlist <- sapply(0:M, function(z) hermite(z, vars))
		return(prodlist)
	} else if (all(M==1)){
		#Just a linear model: TO DO: allow all 1's to use a firt-order interaction model
		prodlist <- cbind(1, vars)
		return(prodlist)
	} else {
		vecs <- mapply(seq, 0, M)
		tmp <- as.matrix(do.call(expand.grid, vecs))
		prodlist <- apply(tmp, 1, function(z) apply(hermite(z, x=vars), 1, prod))
		return(prodlist)
	}
}
####################################################################################################
####################################################################################################
####################################################################################################










