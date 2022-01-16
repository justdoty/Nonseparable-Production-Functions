####################################################################################################
#Specificatiion for the Production Function, Input Decision Rules, and Productivity
####################################################################################################
#The Production Function: (Hicks Neutral Cobb-Douglas and Translog, Non Hicks Neutral Cobb-Douglass and Translog
#and Non Hicks Neutral Tensor Product Hermite Polynomial)
PF <- function(K, L, M, omega){
	prod <- cbind(K, L, M, K*L, L*M, K*M, K^2, L^2, M^2)
	prodf <- cbind(1, omega, prod, sweep(prod, 1, omega, `*`))
	return(prodf)
}
#Labor Decision Rule
LX <- function(K, omega){
	Ldat <- cbind(K, omega)
	LX <- cbind(1, tensor(M=c(3,3), vars=Ldat)[,-1])
	return(LX)
}
#Derivative of Labor Rule
LXD <- function(K, omega, par, pos, sdpos){
	M <- c(3,3)
	Ldat <- cbind(K, omega)
	dten <- dtensor(M=M, vars=Ldat, pos=pos, sdpos=sdpos)
	LXD <- dten$prod%*%par[dten$dpos,]
	return(LXD)
}
#Material Input Decision Rule
MX <- function(K, L, omega){
	# Mdat <- cbind(K, omega)
	# M <- c(3,3)
	Mdat <- cbind(K, L, omega)
	M <- c(2,2,2)
	MX <- cbind(1, tensor(M=M, vars=Mdat)[,-1])
	return(MX)
}
#Derivative of Materials Rule
MXD <- function(K, L, omega, par, pos, sdpos){
	# Mdat <- cbind(K, omega)
	# M <- c(3,3)
	Mdat <- cbind(K, L, omega)
	M <- c(2,2,2)
	dten <- dtensor(M=M, vars=Mdat, pos=pos, sdpos=sdpos)
	MXD <- dten$prod%*%par[dten$dpos,]
	return(MXD)
}
#Investment Input Decision Rule
IX <- function(K, omega){
	Idat <- cbind(K, omega)
	IX <- cbind(1, tensor(M=c(3,3), vars=Idat)[,-1])
	return(IX)
}
#Derivative of Investment Rule
IXD <- function(K, omega, par, pos, sdpos){
	M <- c(3,3)
	Idat <- cbind(K, omega)
	dten <- dtensor(M=M, vars=Idat, pos=pos, sdpos=sdpos)
	IXD <- dten$prod%*%par[dten$dpos,]
	return(IXD)
}
#Productivity Process for t>1
WX <- function(omega){
	Wdat <- omega
	WX <- cbind(1, omega, omega^2, omega^3)
	return(WX)
}
WXD <- function(omega, par){
	M <- 3
	Wdat <- omega
	dten <- dtensor(M=M, vars=as.matrix(Wdat), pos=0)	
	WXD <- dten$prod%*%par[dten$dpos,]
	return(WXD)
}
#Initial Productivity
WX1 <- function(K, L, M){
	W1dat <- cbind(1, K, K^2, K^3)
	return(W1dat)
}
#Specification for probit regression
WBAR <- function(omega, K){
	wbardat <- cbind(1, omega, K, omega*K, omega^2, K^2)
	return(wbardat)
}
#########################################################################################################
# Function for recurrence construction of a Hermite polynomial
########################################################################################################
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
#This function quickly computes column-wise products of a matrix
colprod <- function(mat){
	tmp <- mat[,1]
	for (j in 2:ncol(mat)){
		tmp <- tmp*mat[,j]
	}
	return(tmp)
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
#Function that creates a tensor product from hermite polynomials
###############################################################################
tensor <- function(M, vars){
	#standardize the data (optional)
	vars <- as.matrix(vars)
	#If we just want a univariate hermite polynomial (tensor.prod(2, data))
	#For example, when the productivity process evolves exogeneously
	if (length(M)==1){
		if (length(M)<ncol(vars)){
			print("Error: Number of vars should be equal to one")
		}
		prodlist <- sapply(0:M, function(z) hermite(z, vars))
		return(prodlist)
	#This covers the case when mapply cant return a list when all M is the same
	#For example, tensor.prod(c(2,2,2,2), data) means that the dimension
	#of each univariate hermite polynomial is the same
	} else if (length(unique(M))==1){
		if (length(M)!=ncol(vars)){
			print('Error: M needs to equal to the number of vars')
		}
		vecs <- mapply(seq, 0, M)
		tmp <- split(vecs, rep(1:ncol(vecs), each=nrow(vecs)))
		tmp <- as.matrix(do.call(expand.grid, tmp))
		prodlist <- apply(tmp, 1, function(z) colprod(hermite(z, x=vars)))
		return(prodlist)
	#This covers the case when the user wants to specify a tensor product hermite polynomial
	#with different dimensions of the univariate hermite polynomials
	#For example, tensor.prod(c(1,2,3,4), data)
	} else {
		if (length(M)!=ncol(vars)){
			print('Error: M needs to equal to the number of vars')
		}
		vecs <- mapply(seq, 0, M)
		tmp <- as.matrix(do.call(expand.grid, vecs))
		prodlist <- apply(tmp, 1, function(z) colprod(hermite(z, x=vars)))
		return(prodlist)
	}
}
####################################################################
#Function to calculate derivative of a univariate hermite polynomial
####################################################################
dhermite <- function(n,x){
	dherm <- n*hermite(n=n-1, x=x)
	return(dherm)
}
########################################################################
#Function to compute derivative of a tensor product polynomial
#pos refers to the position of the variable of interest in vars
########################################################################
dtensor <- function(M, vars, pos, sdpos){
	#standardize the data (optional)
	vars <- as.matrix(vars)
	if (pos>ncol(vars)){
		print("Error: pos should be less than or equal number of column in vars")
	}
	#If we just want a univariate hermite polynomial (tensor.prod(2, data))
	#For example, when the productivity process evolves exogeneously
	if (length(M)==1){
		if (length(M)<ncol(vars)){
			print("Error: Number of vars should be equal to one")
		}
		prodlist <- sapply(1:M, function(z) dhermite(z, vars)/sdpos)
		dpos <- 2:(M+1)
		return(list(prod=prodlist, dpos=dpos))
	#This covers the case when mapply cant return a list when all M is the same
	#For example, tensor.prod(c(2,2,2,2), data) means that the dimension
	#of each univariate hermite polynomial is the same
	} else if (length(unique(M))==1){
		if (length(M)!=ncol(vars)){
			print('Error: M needs to equal to the number of vars')
		}
		vecs <- mapply(seq, 0, M)
		tmp <- split(vecs, rep(1:ncol(vecs), each=nrow(vecs)))
		tmp <- as.matrix(do.call(expand.grid, tmp))
		dpos <- which(tmp[,pos]!=0)
		tmp <- tmp[tmp[,pos]!=0,]
		prodlist <- apply(tmp, 1, function(z) colprod(cbind(dhermite(z[pos], as.matrix(vars[,pos]))/sdpos, hermite(z[-pos], as.matrix(vars[,-pos])))))
		return(list(prod=prodlist, dpos=dpos))
	#This covers the case when the user wants to specify a tensor product hermite polynomial
	#with different dimensions of the univariate hermite polynomials
	#For example, tensor.prod(c(1,2,3,4), data)
	} else {
		if (length(M)!=ncol(vars)){
			print('Error: M needs to equal to the number of vars')
		}
		vecs <- mapply(seq, 0, M)
		tmp <- as.matrix(do.call(expand.grid, vecs))
		dpos <- which(tmp[,pos]!=0)
		tmp <- tmp[tmp[,pos]!=0,]
		prodlist <-  apply(tmp, 1, function(z) colprod(cbind(dhermite(z[1], as.matrix(vars[,pos]))/sdpos, hermite(z[-pos], as.matrix(vars[,-pos])))))
		return(list(prod=prodlist, dpos=dpos))
	}
}






