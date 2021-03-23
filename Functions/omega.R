setwd('/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Functions')
#This function is used to do a basic LP estimation of productivity to use as initial model estimates
#and initial draw for productivity in the stEM algorithm
source('Auxfuns.R')
# source('NLPFQR/FUN/Auxfuns.R')
#Function that estimates productivity (omega)
omega_est <- function(idvar, timevar, Y, K, L, M){
	regvars <- data.frame(reg1=K, reg2=L, reg3=M, reg4=M*K, reg5=K^2, reg6=M^2, reg7=M*(K^2), reg8=(M^2)*K, reg9=K^3, reg10=M^3)
  #First Stage 
	LPfirststage <- lm(Y~as.matrix(regvars[, grepl('reg', colnames(regvars))]))
  LPLabor <-  as.numeric(coef(LPfirststage)[3])
	#Initial Values for search
	LP_LM <- lm(Y~K+L+M)
	LPinit <- as.numeric(coef(LP_LM)[-c(1,3)])
	#2nd stage of LP
	#obtain consistent estimates of materials and capital
	phihat <- fitted(LPfirststage)-as.matrix(L)%*%LPLabor
  LPphi <- phihat
	lagdata <- lagdata(idvar, cbind(Y, K, L, M, LPphi))
	names(lagdata) <- c("idvar", "Ycon", "Kcon", "Lcon", "Mcon", "LPphicon",
		"Ylag", "Klag", "Llag", "Mlag", "LPphilag")
	#LP Output
	LPmY <-  as.matrix(lagdata$Ycon-lagdata$Lcon*LPLabor)
	#LP Contemporary State Variables
  LPmX <- cbind(lagdata$Kcon, lagdata$Mcon)
  #LP Lagged State Variables
  LPmlX <- cbind(lagdata$Klag, lagdata$Mlag)
  #LP Contemporary phi estimates
  LPfitphi <- as.matrix(lagdata$LPphicon)
  #LP Lagged phi estimates
  LPfitlagphi <- as.matrix(lagdata$LPphilag)
  #LP Instruments 
  LPmZ <- cbind(lagdata$Kcon, lagdata$Mlag)
  #LP estimates for Capital and Labor
  LPhat <- optim(par=LPinit, fn=function(b) LPobj(b, mY=LPmY, mX=LPmX, mlX=LPmlX, 
    mZ=LPmZ, fitphi=LPfitphi, fitlagphi=LPfitlagphi), gr=NULL, method="L-BFGS-B", 
    lower=c(0,0), upper=c(1,1))$par
  omega <- phihat-cbind(K, M)%*%as.matrix(as.numeric(LPhat))
  return(omega)
}
############################################################################################
#Functions for Estimating LP Coefficients
###########################################################################################
#Function that defines the residuals
LP_Lambda <- function(b, mY, mX, mlX, fitphi, fitlagphi){
  b <- as.matrix(as.numeric(b))
  A <- fitphi-mX%*%b[1:ncol(mX)]
  B <- fitlagphi-mlX%*%b[1:ncol(mX)]
  step1 <- lm(A~B)
  step1param <- as.numeric(coef(step1))
  wfit <- cbind(1, B)%*%step1param
  resid <- mY-mX%*%b[1:ncol(mX)]-wfit
  return(resid)
} 
#LP GMM objective function
LPobj <- function(b, mY, mX, mlX, mZ, fitphi, fitlagphi){
  xifit <- LP_Lambda(b=b, mY=mY, mX=mX, mlX=mlX, fitphi=fitphi, fitlagphi=fitlagphi)
  mW <- solve(crossprod(mZ))/nrow(mZ)
  go <- t(crossprod(mZ, xifit))%*%mW%*%(crossprod(mZ, xifit))
  return(go)
}