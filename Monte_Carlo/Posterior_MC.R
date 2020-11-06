# setwd('/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Functions')
source('NLPFQR/FUN/Tensors.R')
source('NLPFQR/FUN/Auxfuns.R')
require(dplyr)
############################################################################################################
#Function that defines the (log) posterior density for which to evaluate the Metropolis Hastings algorithm (target distribution)
#Requires the file Tensors.R to evaluate the tensor products
#Requires the file Auxfuns.R for overhead computations
############################################################################################################
posterior <- function(data, vectau, par){
	#Load the data, degrees of polynomials, tau vector, guess of the unobservables
	#and the current guess for unobservables
	parY <- par$resYinit; yb1 <- par$yb1init; ybL <- par$ybLinit; 
	parL <- par$resLinit; parM <- par$resMinit; 
	parWT <- par$resWTinit; wtb1 <- par$wtb1init; wtbL <- par$wtbLinit;
	parW1 <- par$resW1init; parI <- par$resIinit
	#Some data preparation for posterior calculations
	data <- data.frame(data)
	names(data) <- c("idvar", "Y", "K", "L", "M", "I", "U")
	WTdata <- lagdata(idvar=data$idvar, X=data$U)
	names(WTdata) <- c("idvar", "Ucon", "Ulag")
	W1data <- t0data(idvar=data$idvar, X=data$U)
	names(W1data) <- c("idvar", "U1")
	#Create Hermite Tensor Products
	YX <- cbind(1, data$K, data$L, data$M, data$U)
	LX <- cbind(1, data$K, data$U)
	MX <- cbind(1, data$K, data$U)
	IX <- cbind(1, data$K, data$U)
	WXT  <- cbind(1, WTdata$Ulag)
	#Likelihood of Output Data#############################################################################
	densY <- as.matrix(rowSums(sapply(1:(length(vectau)-1), function(q) ((vectau[q+1]-vectau[q])/(YX%*%(parY[,q+1]-parY[,q])))*
							(YX%*%parY[,q]<data$Y)*(YX%*%parY[,q+1]>=data$Y))))+vectau[1]*yb1*exp(yb1*(data$Y-YX%*%parY[,1]))*(data$Y<=YX%*%parY[,1])+
					(1-vectau[length(vectau)])*ybL*exp(-ybL*(data$Y-YX%*%parY[,length(vectau)]))*(data$Y>YX%*%parY[,length(vectau)])
	densY <- as.matrix((data.frame(idvar=data$idvar, densY=densY) %>% group_by(idvar) %>% summarise(densN=sum(log(densY)),.groups = 'drop'))$densN)
	#Likelihood of Labor Data#############################################################################
	ResLX <- data$L-LX%*%parL[-length(parL)]
	densL <- 1/sqrt(parL[length(parL)])*dnorm(ResLX/sqrt(parL[length(parL)]))
	densL <- as.matrix((data.frame(idvar=data$idvar, densL=densL) %>% group_by(idvar) %>% summarise(densN=sum(log(densL)),.groups = 'drop'))$densN)
	#Likelihood of Materials Data########################################################################
	ResMX <- data$M-MX%*%parM[-length(parM)]
	densM <- 1/sqrt(parM[length(parM)])*dnorm(ResMX/sqrt(parM[length(parM)]))
	densM <- as.matrix((data.frame(idvar=data$idvar, densM=densM) %>% group_by(idvar) %>% summarise(densN=sum(log(densM)),.groups = 'drop'))$densN)
	#Prior Omega_{t} Data############################################################################
	densWT <- as.matrix(rowSums(sapply(1:(length(vectau)-1), function(q) ((vectau[q+1]-vectau[q])/(WXT%*%(parWT[,q+1]-parWT[,q])))*
							(WXT%*%parWT[,q]<WTdata$Ucon)*(WXT%*%parWT[,q+1]>=WTdata$Ucon))))+vectau[1]*wtb1*exp(wtb1*(WTdata$Ucon-WXT%*%parWT[,1]))*(WTdata$Ucon<=WXT%*%parWT[,1])+
					(1-vectau[length(vectau)])*wtbL*exp(-wtbL*(WTdata$Ucon-WXT%*%parWT[,length(vectau)]))*(WTdata$Ucon>WXT%*%parWT[,length(vectau)])
	densWT <- as.matrix((data.frame(idvar=WTdata$idvar, densWT=densWT) %>% group_by(idvar) %>% summarise(densN=sum(log(densWT)),.groups = 'drop'))$densN)
    #Prior Omega_{0} Data (Normal)############################################################################
	densW1 <- as.matrix(dnorm(W1data$U1, mean=parW1[1], sd=sqrt(parW1[2])))
	densW1 <- log(densW1)
	#Likelihood of Investment Data (Log-Normal)#########################################################
	ResIX <- data$I-IX%*%parI[-length(parI)]
	densI <- 1/sqrt(parI[length(parI)])*dnorm(ResIX/sqrt(parI[length(parI)]))
	densI <- as.matrix((data.frame(idvar=data$idvar, densI=densI) %>% group_by(idvar) %>% summarise(densN=sum(log(densI)),.groups = 'drop'))$densN)
	#Final posterior density########################################################################
	finaldens <- densY+densL+densM+densWT+densW1+densI
	return(finaldens)

}
######################################################################################################
######################################################################################################
######################################################################################################





