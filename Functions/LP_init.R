setwd('/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Functions')
#This function is used to do a basic ACF estimation of productivity to use as initial model estimates
#and initial draw for productivity in the stEM algorithm
source('Auxfuns.R')
#Function that estimates productivity (omega)
omega_est <- function(idvar, timevar, Y, K, L, M){
	regvars <- data.frame(reg1=K, reg2=L, reg3=M, reg4=K*L, reg5=K*M, 
    reg6=L*M, reg7=K^2, reg8=L^2, reg9=M^2)
	#First Stage 
	ACFfirststage <- lm(data$Y~as.matrix(regvars[, grepl('reg', colnames(regvars))]))
	#Initial Values for search
	ACF_LM <- lm(Y~K+L+M)
	ACFinit <- as.numeric(coef(ACF_LM)[-1])
	#2nd stage of LP
	#obtain consistent estimates of materials and capital
	phihat <- fitted(ACFfirststage)
    ACFphi <- phihat
	lagdata <- lagdata(idvar, cbind(Y, K, L, M, ACFphi))
	names(lagdata) <- c("idvar", "Ycon", "Kcon", "Lcon", "Mcon", "ACFphicon",
		"Ylag", "Klag", "Llag", "Mlag", "ACFphilag")
	#ACF Output
	ACFmY <-  as.matrix(lagdata$Ycon)
	#ACF Contemporary State Variables
    ACFmX <- cbind(lagdata$Kcon, lagdata$Lcon)
    #ACF Lagged State Variables
    ACFmlX <- cbind(lagdata$Klag, lagdata$Llag)
    #ACF Contemporary phi estimates
    ACFfitphi <- as.matrix(lagdata$ACFphicon)
    #ACF Lagged phi estimates
    ACFfitlagphi <- as.matrix(lagdata$ACFphilag)
    #ACF Instruments 
    ACFmZ <- cbind(1, lagdata1$Kcon, lagdata1$Llag)
    #ACF estimates for Capital and Labor
    ACFhat <- optim(par=ACFinit, fn=function(b) ACFobj(b, mY=ACFmY, mX=ACFmX, mlX=ACFmlX, 
    	mZ=ACFmZ, fitphi=ACFfitphi, fitlagphi=ACFfitlagphi), gr=NULL, method="L-BFGS-B", 
    lower=c(0,0,0), upper=c(1,1,1))$par
    omega <- phihat-cbind(K, L, M)%*%as.matrix(as.numeric(ACFhat))
    return(omega)
}
############################################################################################
#Functions for Estimating ACF Coefficients
###########################################################################################
#Function that defines the residuals
ACF_Lambda <- function(b, mY, mX, mlX, fitphi, fitlagphi){
  b <- as.matrix(as.numeric(b))
  A <- fitphi-mX%*%b[1:ncol(mX)]
  B <- fitlagphi-mlX%*%b[1:ncol(mX)]
  step1 <- lm(A~B-1)
  step1param <- as.numeric(coef(step1))
  xifit <- A-B*step1param
  return(xifit)
} 
#ACF GMM objective function
ACFobj <- function(b, mY, mX, mlX, mZ, fitphi, fitlagphi){
  xifit <- ACF_Lambda(b=b, mY=mY, mX=mX, mlX=mlX, fitphi=fitphi, fitlagphi=fitlagphi)
  mW <- solve(crossprod(mZ))/nrow(mZ)
  go <- t(crossprod(mZ, xifit))%*%mW%*%(crossprod(mZ, xifit))
  return(go)
}