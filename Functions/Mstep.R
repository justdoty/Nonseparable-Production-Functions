# setwd('/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Functions')
# require(dplyr)
source('NLPFQR/FUN/Auxfuns.R')
#####################################################################################################
#Functions to estimate the exponential tail parameters on (0, \tau_{1}) and (\tau_{L}, 1)
#####################################################################################################
expb <- function(YX, XX, par1, parL){
	b1 <- -mean(YX<=XX%*%par1)/mean((YX-XX%*%par1)*(YX<=XX%*%par1))
	bL <- mean(YX>XX%*%parL)/mean((YX-XX%*%parL)*(YX>XX%*%parL))
	return(list(b1=b1, bL=bL))
}
