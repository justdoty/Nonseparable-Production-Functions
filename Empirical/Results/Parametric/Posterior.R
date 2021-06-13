# setwd('/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Functions')
# source('Tensors.R')
# source('Auxfuns.R')
# source('NLPFQR/DATA/Parametric/Tensors.R')
# source('NLPFQR/DATA/Parametric/Posterior.R')
# require(dplyr)
############################################################################################################
#Function that defines the (log) posterior density for which to evaluate the Metropolis Hastings algorithm (target distribution)
#Requires the file Tensors.R to evaluate the tensor products
#Requires the file Auxfuns.R for overhead computations
############################################################################################################
#Density Function for the Conditional Quantiles
densq <- function(idvar, YX, vectau, pred, parb1, parbL){
	d1 <- (log(vectau[1])+log(parb1)+parb1*(YX-pred[,1]))*(YX<=pred[,1])
	densq <- ((vectau[2]-vectau[1])/(pred[,2]-pred[,1]))*(pred[,1]<YX)*(pred[,2]>=YX)
	for (q in 2:(length(vectau)-1)){ 
		densq <- densq+((vectau[q+1]-vectau[q])/(pred[,q+1]-pred[,q]))*(pred[,q]<YX)*(pred[,q+1]>=YX)
	}
	densq <- ifelse(densq>0, log(densq), 0)
	dL <- (log(1-vectau[length(vectau)])+log(parbL)-parbL*(YX-pred[,length(vectau)]))*(YX>pred[,length(vectau)])
	densq <- as.matrix(d1)+densq+as.matrix(dL)
	densq <- as.matrix((data.frame(idvar=idvar, densq) %>% group_by(idvar) %>% summarise(densN=sum(densq),.groups = 'drop'))$densN)
	return(densq)
	}
#For distributions that are specified as log-normal (T=1)
densN <- function(idvar, YX, XX, par){
	var <- par[length(par)]^2
	resx <- YX-XX%*%as.matrix(par[-length(par)])
	densq <- -0.5*log(2*pi)-0.5*log(var)-(0.5/var)*resx^2 
	densq <- as.matrix((data.frame(idvar=idvar, densq) %>% group_by(idvar) %>% summarise(densN=sum(densq),.groups = 'drop'))$densN)
	return(densq)
}
###############################################################################################################
posterior <- function(tdata, lagdata, t1data, vectau, par){
	#Final posterior density########################################################################
	densY <- densq(idvar=tdata$idvar, YX=tdata$Y, vectau=vectau, pred=predict(par$resYinit, tdata), parb1=par$yb1init, parbL=par$ybLinit)
	densWT <- densq(idvar=lagdata$idvar, YX=lagdata$Wcon, vectau=vectau, pred=predict(par$resWTinit, lagdata), parb1=par$wtb1init, parbL=par$wtbLinit)
	densL <- densN(idvar=tdata$idvar, YX=tdata$L, XX=LX(A=tdata$A, K=tdata$K, omega=tdata$W), par=par$resLinit)
	densM <- densN(idvar=tdata$idvar, YX=tdata$M, XX=MX(A=tdata$A, K=tdata$K, omega=tdata$W), par=par$resMinit)
	densKT <- densN(idvar=lagdata$idvar, YX=lagdata$Kcon, XX=KX(A=lagdata$Acon, K=lagdata$Klag, I=lagdata$Ilag, omega=lagdata$Wlag), par=par$resKTinit)
	densW1 <- densN(idvar=t1data$idvar, YX=t1data$W1, XX=W1X(A=t1data$A1), par=par$resW1init)
	densK1 <- densN(idvar=t1data$idvar, YX=t1data$K1, XX=K1X(A=t1data$A1, omega=t1data$W1), par=par$resK1init)
	finaldens <- densY+densL+densM+densWT+densKT+densK1+densW1
	return(finaldens)

}
######################################################################################################
######################################################################################################
######################################################################################################


