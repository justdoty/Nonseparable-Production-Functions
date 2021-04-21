setwd('/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Functions')
source('Tensors.R')
source('Auxfuns.R')
# source('NLPFQR/FUN/Tensors.R')
# source('NLPFQR/FUN/Auxfuns.R')
# require(dplyr)
############################################################################################################
#Function that defines the (log) posterior density for which to evaluate the Metropolis Hastings algorithm (target distribution)
#Requires the file Tensors.R to evaluate the tensor products
#Requires the file Auxfuns.R for overhead computations
############################################################################################################
posterior <- function(idvar, Y, K, L, M, I, omega, vectau, par, pp, ppd){
	#Load the data, degrees of polynomials, tau vector, guess of the unobservables
	#and the current guess for unobservables
	parY <- par$resYinit; yb1 <- par$yb1init; ybL <- par$ybLinit; 
	parL <- par$resLinit; lb1 <- par$lb1init; lbL <- par$lbLinit;
	parM <- par$resMinit; mb1 <- par$mb1init; mbL <- par$mbLinit;
	parWT <- par$resWTinit; wtb1 <- par$wtb1init; wtbL <- par$wtbLinit;
	parW1 <- par$resW1init; parI <- par$resIinit
	#Load Cubic Polynomial Functions
	ppy <- pp$y; ppl <- pp$l; ppm <- pp$m; ppw <- pp$w
	ppyd <- ppd$y; ppld <- ppd$l; ppmd <- ppd$m; ppwd <- ppd$w
	#Some data preparation for posterior calculations
	WTdata <- lagdata(idvar=idvar, X=omega)
	names(WTdata) <- c("idvar", "Wcon", "Wlag")
	W1data <- t0data(idvar=idvar, X=omega)
	names(W1data) <- c("idvar", "W1")
	#Create Hermite Tensor Products
	YX <- translog(K=K, L=L, M=M, omega=omega)
	LX <- LX(K=K, omega=omega)
	MX <- MX(K=K, omega=omega)
	IX <- IX(K=K, omega=omega)
	WX <- WX(omega=WTdata$Wlag)
	Wcon <- as.matrix(WTdata$Wcon)
	W1 <- as.matrix(W1data$W1)
	#Likelihood of Output Data#############################################################################
	densY <- densfun(PY=(Y-omega), PX=YX, vectau=vectau, par=parY, parb=c(yb1, ybL), pp=ppy, ppd=ppyd, steps=50, tol=1e-3)
	densY <- as.matrix((data.frame(idvar=idvar, densY) %>% group_by(idvar) %>% summarise(densN=sum(log(densY)),.groups = 'drop'))$densN)
	#Likelihood of Labor Data#############################################################################
	densL <- densfun(PY=L, PX=LX, vectau=vectau, par=parL, parb=c(lb1, lbL), pp=ppl, ppd=ppld, steps=50, tol=1e-3)
	densL <- as.matrix((data.frame(idvar=idvar, densL) %>% group_by(idvar) %>% summarise(densN=sum(log(densL)),.groups = 'drop'))$densN)
	#Likelihood of Output Data#############################################################################
	densM <- densfun(PY=M, PX=MX, vectau=vectau, par=parM, parb=c(mb1, mbL), pp=ppm, ppd=ppmd, steps=50, tol=1e-3)
	densM <- as.matrix((data.frame(idvar=idvar, densM) %>% group_by(idvar) %>% summarise(densN=sum(log(densM)),.groups = 'drop'))$densN)
	#Likelihood of Omega_{t} Data#############################################################################
	densWT <- densfun(PY=Wcon, PX=WX, vectau=vectau, par=parWT, parb=c(wtb1, wtbL), pp=ppw, ppd=ppwd, steps=50, tol=1e-3)
	densWT <- as.matrix((data.frame(idvar=WTdata$idvar, densWT) %>% group_by(idvar) %>% summarise(densN=sum(log(densWT)),.groups = 'drop'))$densN)
    #Prior Omega_{0} Data (Normal)############################################################################
	densW1 <- log((1/sqrt(2*pi*parW1[2]))*exp(-.5*((W1-parW1[1])/sqrt(parW1[2]))^2))
	#Likelihood of Investment Data (Log-Normal)#########################################################
	densI <- (1/sqrt(2*pi*parI[length(parI)]))*exp(-.5*(I-IX%*%parI[-length(parI)])^2/parI[length(parI)])
	densI <- as.matrix((data.frame(idvar=idvar, densI) %>% group_by(idvar) %>% summarise(densN=sum(log(densI)),.groups = 'drop'))$densN)
	#Final Density
	finaldens <- densY+densL+densM+densWT+densW1+densI
	return(finaldens)

}
######################################################################################################
######################################################################################################
######################################################################################################

######################################################################################################
######################################################################################################
######################################################################################################





