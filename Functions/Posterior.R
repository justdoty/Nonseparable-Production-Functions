setwd('/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Functions')
source('Tensors.R')
source('Auxfuns.R')
require(dplyr)
############################################################################################################
#Function that defines the (log) posterior density for which to evaluate the Metropolis Hastings algorithm (target distribution)
#Requires the file Tensors.R to evaluate the tensor products
#Requires the file Auxfuns.R for overhead computations
############################################################################################################
posterior <- function(data, MH, vectau, par){
	#Load the data, degrees of polynomials, tau vector, guess of the unobservables
	#and the current guess for unobservables
	MY <- MH$MY; ML <- MH$ML; MM <- MH$MM; MW <- MH$MW; MW1 <- MH$MW1; MI <- MH$MI
	parY <- par$resYinit; yb1 <- par$yb1init; ybL <- par$ybLinit; 
	parL <- par$resLinit; lb1 <- par$lb1init; lbL <- par$lbLinit;
	parM <- par$resMinit; mb1 <- par$mb1init; mbL <- par$mbLinit;
	parWT <- par$resWTinit; wtb1 <- par$wtb1init; wtbL <- par$wtbLinit;
	parW0 <- par$resW0init; w0b1 <- par$w0b1init; w0bL <- par$w0bLinit;
	parI <- par$resIinit
	#Some data preparation for posterior calculations
	data <- data.frame(data)
	names(data) <- c("idvar", "Y", "K", "L", "M", "A", "I", "U")
	WTdata <- lagdata(idvar, cbind(data$A, data$U))
	names(WTdata) <- c("idvar", "Acon", "Ucon", "Alag", "Ulag")
	W0data <- t0data(idvar, cbind(data$A, data$U))
	names(W0data) <- c("idvar", "A0", "U0")
	#Create Hermite Tensor Products
	YX <- tensor.prod(MY, cbind(data$K, data$L, data$M, data$A, data$U))
	LX <- tensor.prod(ML, cbind(data$K, data$A, data$U))
	MX <- tensor.prod(MM, cbind(data$K, data$A, data$U))
	IX <- tensor.prod(MI, cbind(data$K, data$A, data$U))
	WXT  <- tensor.prod(MW, cbind(WTdata$Acon, WTdata$Ulag))
	WX0  <- tensor.prod(MW1, as.matrix(W0data$A0))
	#Likelihood of Output Data#############################################################################
	densY <- as.matrix(rowSums(sapply(1:(length(vectau)-1), function(q) ((vectau[q+1]-vectau[q])/(YX%*%(parY[,q+1]-parY[,q])))*
							(YX%*%parY[,q]<data$Y)*(YX%*%parY[,q+1]>=data$Y))))+vectau[1]*yb1*exp(yb1*(data$Y-YX%*%parY[,1]))*(data$Y<=YX%*%parY[,1])+
					(1-vectau[length(vectau)])*ybL*exp(-ybL*(data$Y-YX%*%parY[,length(vectau)]))*(data$Y>YX%*%parY[,length(vectau)])
	densY <- as.matrix((data.frame(idvar, densY) %>% group_by(idvar) %>% summarise(densN=sum(log(densY)),.groups = 'drop'))$densN)
	#Likelihood of Labor Data#############################################################################
	densL <- as.matrix(rowSums(sapply(1:(length(vectau)-1), function(q) ((vectau[q+1]-vectau[q])/(LX%*%(parL[,q+1]-parL[,q])))*
							(LX%*%parL[,q]<data$L)*(LX%*%parL[,q+1]>=data$L))))+vectau[1]*lb1*exp(lb1*(data$L-LX%*%parL[,1]))*(data$L<=LX%*%parL[,1])+
					(1-vectau[length(vectau)])*lbL*exp(-lbL*(data$L-LX%*%parL[,length(vectau)]))*(data$L>LX%*%parL[,length(vectau)])
	densL <- as.matrix((data.frame(idvar, densL) %>% group_by(idvar) %>% summarise(densN=sum(log(densL)),.groups = 'drop'))$densN)
	#Likelihood of Materials Data########################################################################
	densM <- as.matrix(rowSums(sapply(1:(length(vectau)-1), function(q) ((vectau[q+1]-vectau[q])/(MX%*%(parM[,q+1]-parM[,q])))*
							(MX%*%parM[,q]<data$M)*(MX%*%parM[,q+1]>=data$M))))+vectau[1]*mb1*exp(mb1*(data$M-MX%*%parM[,1]))*(data$M<=MX%*%parM[,1])+
					(1-vectau[length(vectau)])*mbL*exp(-mbL*(data$M-MX%*%parM[,length(vectau)]))*(data$M>MX%*%parM[,length(vectau)])
	densM <- as.matrix((data.frame(idvar, densM) %>% group_by(idvar) %>% summarise(densN=sum(log(densM)),.groups = 'drop'))$densN)
	#Prior Omega_{t} Data############################################################################
	densWT <- as.matrix(rowSums(sapply(1:(length(vectau)-1), function(q) ((vectau[q+1]-vectau[q])/(WXT%*%(parWT[,q+1]-parWT[,q])))*
							(WXT%*%parWT[,q]<WTdata$Ucon)*(WXT%*%parWT[,q+1]>=WTdata$Ucon))))+vectau[1]*wtb1*exp(wtb1*(WTdata$Ucon-WXT%*%parWT[,1]))*(WTdata$Ucon<=WXT%*%parWT[,1])+
					(1-vectau[length(vectau)])*wtbL*exp(-wtbL*(WTdata$Ucon-WXT%*%parWT[,length(vectau)]))*(WTdata$Ucon>WXT%*%parWT[,length(vectau)])
	densWT <- as.matrix((data.frame(idvar=WTdata$idvar, densWT) %>% group_by(idvar) %>% summarise(densN=sum(log(densWT)),.groups = 'drop'))$densN)
    #Prior Omega_{0} Data############################################################################
	densW0 <- as.matrix(rowSums(sapply(1:(length(vectau)-1), function(q) ((vectau[q+1]-vectau[q])/(WX0%*%(parW0[,q+1]-parW0[,q])))*
							(WX0%*%parW0[,q]<W0data$U0)*(WX0%*%parW0[,q+1]>=W0data$U0))))+vectau[1]*w0b1*exp(w0b1*(W0data$U0-WX0%*%parW0[,1]))*(W0data$U0<=WX0%*%parW0[,1])+
					(1-vectau[length(vectau)])*w0bL*exp(-w0bL*(W0data$U0-WX0%*%parW0[,length(vectau)]))*(W0data$U0>WX0%*%parW0[,length(vectau)])
	densW0 <- log(densW0)
	#Likelihood of Investment Data (Log-Normal)#########################################################
	ResIX <- data$I-IX%*%parI[-length(parI)]
	densI <- 1/sqrt(parI[length(parI)])*dnorm(ResIX/sqrt(parI[length(parI)]))
	densI <- as.matrix((data.frame(idvar, densI) %>% group_by(idvar) %>% summarise(densN=sum(log(densI)),.groups = 'drop'))$densN)
	#Final posterior density########################################################################
	finaldens <- densY+densL+densM+densWT+densW0+densI
	return(finaldens)

}
######################################################################################################
######################################################################################################
######################################################################################################





