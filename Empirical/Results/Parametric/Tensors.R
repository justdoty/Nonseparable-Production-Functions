##########################################################################################################
##########################################################################################################
#This file creates the tensor product of Hermite Polynomials and Specifications of the decision rules/accumulation rules
##########################################################################################################

####################################################################################################
#Specificatiion for the Production Function, Decision Rules, and Accumulation Processes
####################################################################################################
#The Production Function: (Hicks Neutral Cobb-Douglas and Translog, Non Hicks Neutral Cobb-Douglass and Translog
#and Non Hicks Neutral Tensor Product Hermite Polynomial)
PF <- function(A, K, L, M, omega, method){
	if (method=="cobbN"){
		prodf <- cbind(1, A, K, L, M, omega)
	} else if (method=="transN") {
		A <- (A-mean(A))/sd(A)
		K <- (K-mean(K))/sd(K)
		L <- (L-mean(L))/sd(K)
		M <- (M-mean(M))/sd(K)
		prodf <- cbind(1, A, K, L, M, K*L, L*M, K*M, K^2, L^2, M^2)
	} else if (method=="cobb"){
		A <- (A-mean(A))/sd(A)
		K <- (K-mean(K))/sd(K)
		L <- (L-mean(L))/sd(K)
		M <- (M-mean(M))/sd(K)
		omega <- (omega-mean(omega))/sd(omega)
		prod <- cbind(omega, A, K, L, M)
		prodf <- cbind(1, prod, sweep(prod[,-c(1:2)], 1, omega, `*`))
	} else if (method=="trans"){
		A <- (A-mean(A))/sd(A)
		K <- (K-mean(K))/sd(K)
		L <- (L-mean(L))/sd(K)
		M <- (M-mean(M))/sd(K)
		omega <- (omega-mean(omega))/sd(omega)
		prod <- cbind(omega, A, K, L, M, K*L, L*M, K*M, K^2, L^2, M^2)
		prodf <- cbind(1, prod, sweep(prod[,-c(1:2)], 1, omega, `*`))
	} else if (method=="hermite"){
		prodf <- cbind(1, (A-mean(A)/sd(A)), tensor.prod(c(2,2,2,2), cbind(K, L, M, omega), norm=TRUE)[,-1])
	}
	return(prodf)
}
#Labor Decision Rule
LX <- function(A, K, omega){
	LX <- cbind(1, (A-mean(A)/sd(A)), tensor.prod(c(3,3), cbind(K, omega), norm=TRUE)[,-1])
	return(LX)
}
#Material Input Decision Rule
MX <- function(A, K, omega){
	MX <- cbind(1, (A-mean(A)/sd(A)), tensor.prod(c(3,3), cbind(K, omega), norm=TRUE)[,-1])
	return(MX)
}
#Productivity Process for t>1
WX <- function(A, omega){
	WX <- cbind(1, (A-mean(A)/sd(A)), tensor.prod(4, as.matrix(omega), norm=TRUE)[,-1])
	return(WX)
}
#Capital Accumulation Process for t>1
KX <- function(A, K, I, omega){
	KX <- cbind(1, (A-mean(A)/sd(A)), tensor.prod(c(3,3), cbind(K, omega), norm=TRUE)[,-1])
	return(KX)
}
#Initial Condition for Productivity (t=1)
W1X <- function(A){
	W1 <- cbind(1, tensor.prod(3, as.matrix(A), norm=TRUE)[,-1])
	return(W1)
}
#Initial Condition for Capital (t=1)
K1X <- function(A, omega){
	K1 <- cbind(1, tensor.prod(c(3,3), cbind(A, omega), norm=TRUE)[,-1])
	return(K1)
}
##########################################################################################################
#Function for recurrence construction of a Hermite polynomial
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
################################################################################
########Function to Take Product of Rows for Tensor Product#####################
###############################################################################
rowProdBD <- function(x){
	ncol <- ncol(x)
	y <- x[,1]
	for(i in 2:ncol){
		y <- y*x[,i]
	}
	return(y)
}
#####################################################################################################
#Function that creates a tensor product from the hermite polynomials
#####################################################################################################
tensor.prod <- function(M, vars, norm){
	#Standardize Data (this should always be the case)
	vars <- as.matrix(vars)
	if (norm==TRUE){
		vars <- (vars-colMeans(vars))/apply(vars, 2, sd)
	}
	#If we just want a univariate hermite polynomial (tensor.prod(2, data, norm=TRUE))
	#For example, when the productivity process evolves exogeneously
	if (length(M)==1){
		if (length(M)<ncol(vars)){
			print("Error: Number of vars should be equal to one")
		}
		prodlist <- sapply(0:M, function(z) hermite(z, vars))
		return(prodlist)
	#This covers the case when mapply cant return a list when all M is the same
	#For example, tensor.prod(c(2,2,2,2), data, norm=TRUE) means that the dimension
	#of each univariate hermite polynomial is the same
	} else if (length(unique(M))==1){
		if (length(M)!=ncol(vars)){
			print('Error: M needs to equal to the number of vars')
		}
		vecs <- mapply(seq, 0, M)
		tmp <- split(vecs, rep(1:ncol(vecs), each=nrow(vecs)))
		tmp <- as.matrix(do.call(expand.grid, tmp))
		prodlist <- apply(tmp, 1, function(z) rowProdBD(hermite(z, vars)))
		return(prodlist)
	#This covers the case when the user wants to specify a tensor product hermite polynomial
	#with different dimensions of the univariate hermite polynomials
	#For example, tensor.prod(c(1,2,3,4), data, norm=TRUE)
	} else {
		if (length(M)!=ncol(vars)){
			print('Error: M needs to equal to the number of vars')
		}
		vecs <- mapply(seq, 0, M)
		tmp <- as.matrix(do.call(expand.grid, vecs))
		prodlist <- apply(tmp, 1, function(z) rowProdBD(hermite(z, vars)))
		return(prodlist)
	}
}
######################################################################################################
#Function that creates the derivative of a univariate hermite polynomial
#####################################################################################################
D.hermite <- function(n, x){
	if (n==0){
		Hprime <- matrix(0, nrow=nrow(x))
	} else {
		Hprime <- n*hermite(n=(n-1), x=x)
	}
	return(Hprime)
}
#####################################################################################################
#Function that creates a tensor product from the derivative of univariate hermite polynomials
#####################################################################################################
D.tensor.prod <- function(M, vars, norm){
	#Standardize Data (this should always be the case)
	vars <- as.matrix(vars)
	if (norm==TRUE){
		vars <- (vars-colMeans(vars))/apply(vars, 2, sd)
	}
	#If we just want a univariate hermite polynomial (D.tensor.prod(2, data, norm=TRUE))
	#For example, when the productivity process evolves exogeneously
	if (length(M)==1){
		if (length(M)<ncol(vars)){
			print("Error: Number of vars should be equal to one")
		}
		prodlist <- sapply(0:M, function(z) D.hermite(z, vars))
		return(prodlist)
	#This covers the case when mapply cant return a list when all M is the same
	#For example, tensor.prod(c(2,2,2,2), data, norm=TRUE) means that the dimension
	#of each univariate hermite polynomial is the same
	} else if (length(unique(M))==1){
		if (length(M)!=ncol(vars)){
			print('Error: M needs to equal to the number of vars')
		}
		vecs <- mapply(seq, 0, M)
		tmp <- split(vecs, rep(1:ncol(vecs), each=nrow(vecs)))
		tmp <- as.matrix(do.call(expand.grid, tmp))
		prodlist <- apply(tmp, 1, function(z) D.hermite(z[1], x=as.matrix(vars[,1]))*apply(hermite(z[-1], x=vars[,-1]), 1, prod))
		return(prodlist)
	#This covers the case when the user wants to specify a tensor product hermite polynomial
	#with different dimensions of the univariate hermite polynomials
	#For example, tensor.prod(c(1,2,3,4), data, norm=TRUE)
	} else {
		if (length(M)!=ncol(vars)){
			print('Error: M needs to equal to the number of vars')
		}
		vecs <- mapply(seq, 0, M)
		tmp <- as.matrix(do.call(expand.grid, vecs))
		prodlist <- apply(tmp, 1, function(z) D.hermite(z[1], x=as.matrix(vars[,1]))*apply(hermite(z[-1], x=vars[,-1]), 1, prod))
		return(prodlist)
	}
}


