setwd('/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Functions')
source('Tensors.R')
source('Auxfuns.R')
# source('NLPFQR/FUN/Tensors.R')
# source('NLPFQR/FUN/Auxfuns.R')
require(dplyr)
############################################################################################################
#Function that defines the (log) posterior density for which to evaluate the Metropolis Hastings algorithm (target distribution)
#Requires the file Tensors.R to evaluate the tensor products
#Requires the file Auxfuns.R for overhead computations
############################################################################################################
posterior <- function(idvar, Y, K, L, M, I, omega, vectau, par){
	#Load the data, degrees of polynomials, tau vector, guess of the unobservables
	#and the current guess for unobservables
	parY <- par$resYinit; yb1 <- par$yb1init; ybL <- par$ybLinit; 
	parL <- par$resLinit; lb1 <- par$lb1init; lbL <- par$lbLinit;
	parM <- par$resMinit; mb1 <- par$mb1init; mbL <- par$mbLinit;
	parWT <- par$resWTinit; wtb1 <- par$wtb1init; wtbL <- par$wtbLinit;
	parW1 <- par$resW1init; parI <- par$resIinit
	#Some data preparation for posterior calculations
	WTdata <- lagdata(idvar=idvar, X=omega)
	names(WTdata) <- c("idvar", "Wcon", "Wlag")
	W1data <- t0data(idvar=idvar, X=omega)
	names(W1data) <- c("idvar", "W1")
	#Create Hermite Tensor Products
	YX <- translog(K=K, L=L, M=M, omega=omega)
	LX <- cbind(1, K, omega, K*omega, K^2, omega^2)
	MX <- cbind(1, K, omega, K*omega, K^2, omega^2)
	IX <- cbind(1, K, omega, K*omega, K^2, omega^2)
	Wcon <- WTdata$Wcon
	WX <- cbind(1, WTdata$Wlag, I(WTdata$Wlag^2), I(WTdata$Wlag^3))
	#Likelihood of Output Data#############################################################################
	densY <- as.matrix(rowSums(sapply(1:(length(vectau)-1), function(q) ((vectau[q+1]-vectau[q])/(YX%*%(parY[,q+1]-parY[,q])))*
									(YX%*%parY[,q]<=(Y-omega))*(YX%*%parY[,q+1]>(Y-omega)))))+vectau[1]*yb1*exp(yb1*((Y-omega)-YX%*%parY[,1]))*((Y-omega)<YX%*%parY[,1])+
									(1-vectau[length(vectau)])*ybL*exp(-ybL*((Y-omega)-YX%*%parY[,length(vectau)]))*((Y-omega)>=YX%*%parY[,length(vectau)])
	densY <- as.matrix((data.frame(idvar=idvar, densY) %>% group_by(idvar) %>% summarise(densN=sum(log(densY)),.groups = 'drop'))$densN)
	# print(summary(densY))
	#Likelihood of Labor Data#############################################################################
	densL <- as.matrix(rowSums(sapply(1:(length(vectau)-1), function(q) ((vectau[q+1]-vectau[q])/(LX%*%(parL[,q+1]-parL[,q])))*
							(LX%*%parL[,q]<=L)*(LX%*%parL[,q+1]>L))))+vectau[1]*lb1*exp(lb1*(L-LX%*%parL[,1]))*(L<LX%*%parL[,1])+
					(1-vectau[length(vectau)])*lbL*exp(-lbL*(L-LX%*%parL[,length(vectau)]))*(L>=LX%*%parL[,length(vectau)])
	densL <- as.matrix((data.frame(idvar=idvar, densL) %>% group_by(idvar) %>% summarise(densN=sum(log(densL)),.groups = 'drop'))$densN)
	#Likelihood of Materials Data########################################################################
	densM <- as.matrix(rowSums(sapply(1:(length(vectau)-1), function(q) ((vectau[q+1]-vectau[q])/(MX%*%(parM[,q+1]-parM[,q])))*
							(MX%*%parM[,q]<=M)*(MX%*%parM[,q+1]>M))))+vectau[1]*mb1*exp(mb1*(M-MX%*%parM[,1]))*(M<MX%*%parM[,1])+
					(1-vectau[length(vectau)])*mbL*exp(-mbL*(M-MX%*%parM[,length(vectau)]))*(M>=MX%*%parM[,length(vectau)])
	densM <- as.matrix((data.frame(idvar=idvar, densM) %>% group_by(idvar) %>% summarise(densN=sum(log(densM)),.groups = 'drop'))$densN)
	#Prior Omega_{t} Data############################################################################
	densWT <- as.matrix(rowSums(sapply(1:(length(vectau)-1), function(q) ((vectau[q+1]-vectau[q])/(WX%*%(parWT[,q+1]-parWT[,q])))*
							(WX%*%parWT[,q]<=Wcon)*(WX%*%parWT[,q+1]>Wcon))))+vectau[1]*wtb1*exp(wtb1*(Wcon-WX%*%parWT[,1]))*(Wcon<WX%*%parWT[,1])+
					(1-vectau[length(vectau)])*wtbL*exp(-wtbL*(Wcon-WX%*%parWT[,length(vectau)]))*(Wcon>=WX%*%parWT[,length(vectau)])
	test <- as.matrix(sapply(1:(length(vectau)-1), function(q) (WX%*%(parWT[,q+1]-parWT[,q]))))
	densWT <- as.matrix((data.frame(idvar=WTdata$idvar, densWT) %>% group_by(idvar) %>% summarise(densN=sum(log(densWT)),.groups = 'drop'))$densN)
    # print(summary(densWT))
    #Prior Omega_{0} Data (Normal)############################################################################
	densW1 <- (1/sqrt(2*pi*parW1[2]))*exp(-.5*((W1data$W1-parW1[1])/sqrt(parW1[2]))^2)
	densW1 <- log(densW1)
	# print(summary(densW1))
	#Likelihood of Investment Data (Log-Normal)#########################################################
	ResIX <- (I-IX%*%parI[-length(parI)])/parI[length(parI)]
	densI <- (1/sqrt(2*pi*parI[length(parI)]))*exp(-.5*ResIX^2)
	densI <- as.matrix((data.frame(idvar=idvar, densI) %>% group_by(idvar) %>% summarise(densN=sum(log(densI)),.groups = 'drop'))$densN)
	#Final posterior density########################################################################
	finaldens <- densY+densL+densM+densWT+densW1+densI
	return(finaldens)

}
######################################################################################################
######################################################################################################
######################################################################################################





