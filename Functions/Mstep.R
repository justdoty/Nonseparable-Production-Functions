setwd('/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Functions')
# require(dplyr)
# source('NLPFQR/FUN/Moments.R')
# source('NLPFQR/FUN/Auxfuns.R')
#####################################################################################################
#Functions to estimate the exponential tail parameters on (0, \tau_{1}) and (\tau_{L}, 1)
#####################################################################################################
#For the translog production function
expby <- function(Y, K, L, M, omega, par1, parL){
	YX <- translog(K=K, L=L, M=M, omega=omega)
	b1 <- -mean((Y-omega)<=YX%*%par1)/mean(((Y-omega)-YX%*%par1)*((Y-omega)<=YX%*%par1))
	bL <- mean((Y-omega)>YX%*%parL)/mean(((Y-omega)-YX%*%parL)*((Y-omega)>YX%*%parL))
	return(list(b1=b1, bL=bL))
}
#For the input functions
expblx <- function(X, K, omega, par1, parL){
	XX <- LX(K=K, omega=omega)
	b1 <- -mean(X<=XX%*%par1)/mean((X-XX%*%par1)*(X<=XX%*%par1))
	bL <- mean(X>XX%*%parL)/mean((X-XX%*%parL)*(X>XX%*%parL))
	return(list(b1=b1, bL=bL))
}
expbmx <- function(X, K, omega, par1, parL){
	XX <- MX(K=K, omega=omega)
	b1 <- -mean(X<=XX%*%par1)/mean((X-XX%*%par1)*(X<=XX%*%par1))
	bL <- mean(X>XX%*%parL)/mean((X-XX%*%parL)*(X>XX%*%parL))
	return(list(b1=b1, bL=bL))
}
#For productivity t>1
expbwt <- function(omega, omegalag, par1, parL){
	WX <- WX(omega=omegalag)
	b1 <- -mean(omega<=WX%*%par1)/mean((omega-WX%*%par1)*(omega<=WX%*%par1))
	bL <- mean(omega>WX%*%parL)/mean((omega-WX%*%parL)*(omega>WX%*%parL))
	return(list(b1=b1, bL=bL))
}







