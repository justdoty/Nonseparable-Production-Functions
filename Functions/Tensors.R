# setwd('/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Functions')
##########################################################################################################
##########################################################################################################
#This file creates the tensor product of Hermite Polynomials and Specifications of the decision rules/accumulation rules
##########################################################################################################

####################################################################################################
#Specificatiion for the Production Function, Decision Rules, and Accumulation Processes
####################################################################################################
#The Production Function: (Hicks Neutral Cobb-Douglas and Translog, Non Hicks Neutral Cobb-Douglass and Translog
#and Non Hicks Neutral Tensor Product Hermite Polynomial)
PF <- function(A, T, K, L, M, omega, method){
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
		# A <- (A-mean(A))/sd(A)
		# K <- (K-mean(K))/sd(K)
		# L <- (L-mean(L))/sd(K)
		# M <- (M-mean(M))/sd(K)
		# omega <- (omega-mean(omega))/sd(omega)
		prod <- cbind(A, omega, K, L, M, K*L, L*M, K*M, K^2, L^2, M^2)
		prodf <- cbind(1, prod, sweep(prod[,-c(1,2)], 1, omega, `*`))
	} else if (method=="hermite"){
		prodf <- cbind(1, tensor.prod(c(1,2,2,1), cbind(K, L, M, omega), norm=TRUE)[,-1])
	}
	return(prodf)
}
#Labor Decision Rule
LX <- function(A, T, K, omega){
	# LX <- cbind(1, (A-mean(A)/sd(A)), tensor.prod(c(2,2), cbind(K, omega), norm=TRUE)[,-1])
	# A <- (A-mean(A))/sd(A)
	# K <- (K-mean(K))/sd(K)
	# omega <- (omega-mean(omega))/sd(omega)
	LX <- cbind(1, A, K, omega, omega*K, K^2, omega^2, omega*K^2, K*omega^2, (K^2)*(omega^2), K^3, omega^3)
	return(LX)
}
#Material Input Decision Rule
MX <- function(A, T, K, omega){
	# MX <- cbind(1, (A-mean(A)/sd(A)), tensor.prod(c(2,2), cbind(K, omega), norm=TRUE)[,-1])
	# A <- (A-mean(A))/sd(A)
	# K <- (K-mean(K))/sd(K)
	# omega <- (omega-mean(omega))/sd(omega)
	MX <- cbind(1, A, K, omega, omega*K, K^2, omega^2, omega*K^2, K*omega^2, (K^2)*(omega^2), K^3, omega^3)
	return(MX)
}
#Investment Input Decision Rule
IX <- function(A, T, K, omega){
	# IX <- cbind(1, (A-mean(A)/sd(A)), tensor.prod(c(2,2), cbind(K, omega), norm=TRUE)[,-1])
	# A <- (A-mean(A))/sd(A)
	# K <- (K-mean(K))/sd(K)
	# omega <- (omega-mean(omega))/sd(omega)
	IX <- cbind(1, A, K, omega, omega*K, K^2, omega^2, omega*K^2, K*omega^2, (K^2)*(omega^2), K^3, omega^3)
	return(IX)
}
#Productivity Process for t>1
WX <- function(A, omega){
	# A <- (A-mean(A))/sd(A)
	omega <- (omega-mean(omega))
	# WX <- cbind(1, A, omega, omega^2, omega^3)
	# WX <- cbind(1, (A-mean(A)/sd(A)), tensor.prod(3, as.matrix(omega), norm=TRUE)[,-1])
	WX <- cbind(1, omega, omega^2, omega^3)
	return(WX)
}
#Capital Accumulation Process for t>1
KX <- function(A, T, K, I, omega){
	#Tensor Product Hermite Polynomial (Normalized)
	# KX <- cbind(1, (A-mean(A)/sd(A)), tensor.prod(c(2,2), cbind(K, omega), norm=TRUE)[,-1])
	# A <- (A-mean(A))/sd(A)
	# K <- (K-mean(K))/sd(K)
	# omega <- (omega-mean(omega))/sd(omega)
	KX <- cbind(1, K, omega, omega*K, K^2, omega^2, omega*K^2, K*omega^2, (K^2)*(omega^2), K^3, omega^3)
	return(KX)
}
#Initial Condition for Productivity (t=1)
W1X <- function(A){
	W1 <- cbind(1, A)
	# W1 <- cbind(1, A, A^2)
	return(W1)
}
#Initial Condition for Capital (t=1)
K1X <- function(A, omega){
	K1 <- cbind(1, A, omega)
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

