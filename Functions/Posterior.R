# setwd('/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Functions')
# source('Tensors.R')
# source('Auxfuns.R')
source('NLPFQR/FUN/Tensors.R')
source('NLPFQR/FUN/Auxfuns.R')
# require(dplyr)
############################################################################################################
#Function that defines the (log) posterior density for which to evaluate the Metropolis Hastings algorithm (target distribution)
#Requires the file Tensors.R to evaluate the tensor products
#Requires the file Auxfuns.R for overhead computations
############################################################################################################
posterior <- function(data, vectau, par, method){
	#Density Function for the Conditional Quantiles
	densq <- function(idvar, vectau, YX, XX, par, parb1, parbL){
		densq <- ((vectau[2]-vectau[1])/(XX%*%(par[,2]-par[,1])))*(XX%*%par[,1]<YX)*(XX%*%par[,2]>=YX)
		for (q in 2:(length(vectau)-1)){ 
			densq <- densq+((vectau[q+1]-vectau[q])/(XX%*%(par[,q+1]-par[,q])))*(XX%*%par[,q]<YX)*(XX%*%par[,q+1]>=YX)
		}
		d1 <- array(0, length(idvar))
		dL <- array(0, length(idvar))
		d1[(YX<=XX%*%par[,1])] <- vectau[1]*parb1*exp(parb1*((YX-XX%*%par[,1])[(YX<=XX%*%par[,1])]))
		dL[(YX>XX%*%par[,length(vectau)])] <- (1-vectau[length(vectau)])*parbL*exp(-parbL*((YX-XX%*%par[,length(vectau)])[(YX>XX%*%par[,length(vectau)])]))
		densq <- as.matrix(d1)+densq+as.matrix(dL)
		densq <- as.matrix((data.frame(idvar=idvar, densq) %>% group_by(idvar) %>% summarise(densN=sum(log(densq)),.groups = 'drop'))$densN)
		return(densq)
	}
	#Some data preparation for posterior calculations
	lagdata <- lagdata(data$idvar, cbind(data$Y, data$A, data$K, data$L, data$M, data$I, data$omega))
	names(lagdata) <- c("idvar", "Ycon", "Acon", "Kcon", "Lcon", "Mcon", "Icon", "Wcon",
		"Ylag", "Alag", "Klag", "Llag", "Mlag", "Ilag", "Wlag")
	t1data <- t0data(data$idvar, cbind(data$Y, data$A, data$K, data$L, data$M, data$I, data$omega))
	names(t1data) <- c("idvar", "Y1", "A1", "K1", "L1", "M1", "I1", "W1")
    #Prior Omega_{0} Data (Normal)############################################################################
	densW1 <- function(W1, W1X, parW1){
		ResW1X <- (W1-W1X%*%parW1[-length(parW1)])^2/parW1[length(parW1)]
		densW1 <- (1/sqrt(2*pi*parW1[length(parW1)]))*exp(-.5*ResW1X)
		densW1 <- log(densW1)
		return(densW1)
	}
	#Prior K_{0} Data (Normal)############################################################################
	densK1 <- function(K1, K1X, parK1){
		ResK1X <- (K1-K1X%*%parK1[-length(parK1)])^2/parK1[length(parK1)]
		densK1 <- (1/sqrt(2*pi*parK1[length(parK1)]))*exp(-.5*ResK1X)
		densK1 <- log(densK1)
		return(densK1)
	}
	#Final posterior density########################################################################
	densY <- densq(idvar=data$idvar, vectau=vectau, YX=(data$Y-data$omega), XX=PF(A=data$A, K=data$K, L=data$L, M=data$M, omega=data$omega, method=method), par=par$resYinit, parb1=par$yb1init, parbL=par$ybLinit)
	densL <- densq(idvar=data$idvar, vectau=vectau, YX=data$L, XX=LX(A=data$A, K=data$K, omega=data$omega), par=par$resLinit, parb1=par$lb1init, parbL=par$lbLinit)
	densM <- densq(idvar=data$idvar, vectau=vectau, YX=data$M, XX=MX(A=data$A, K=data$K, omega=data$omega), par=par$resMinit, parb1=par$mb1init, parbL=par$mbLinit)
	densWT <- densq(idvar=lagdata$idvar, vectau=vectau, YX=lagdata$Wcon, XX=WX(A=lagdata$Acon, omega=lagdata$Wlag), par=par$resWTinit, parb1=par$wtb1init, parbL=par$wtbLinit)
	densKT <- densq(idvar=lagdata$idvar, vectau=vectau, YX=lagdata$Kcon, XX=KX(A=lagdata$Acon, K=lagdata$Klag, omega=lagdata$Wlag), par=par$resKTinit, parb1=par$ktb1init, parbL=par$ktbLinit)
	densW1 <- densW1(W1=t1data$W1, W1X=W1X(A=t1data$A1), parW1=par$resW1init)
	densK1 <- densK1(K1=t1data$K1, K1X=K1X(A=t1data$A1, omega=t1data$W1), parK1=par$resK1init)
	finaldens <- densY+densL+densM+densWT+densKT+densK1+densW1
	return(finaldens)

}
######################################################################################################
######################################################################################################
######################################################################################################


