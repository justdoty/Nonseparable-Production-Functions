setwd('/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Functions')
#This function is used to do a basic LP estimation of productivity to use as initial model estimates
#and initial draw for productivity in the stEM algorithm
# source('Auxfuns.R')
# source('NLPFQR/FUN/Auxfuns.R')
#Function that estimates productivity (omega)
omega_est <- function(idvar, timevar, Y, A, K, L, M){
	regvars <- data.frame(reg1=A, reg2=K, reg3=L, reg4=M, reg5=M*K, reg6=K^2, reg7=M^2, reg8=M*(K^2), reg9=(M^2)*K, reg10=K^3, reg11=M^3)
  #First Stage 
	LPfirststage <- lm(Y~as.matrix(regvars[, grepl('reg', colnames(regvars))]))
  LPLabor <-  as.numeric(coef(LPfirststage)[4])
	#Initial Values for search
	LP_LM <- lm(Y~A+K+L+M)
	LPinit <- as.numeric(coef(LP_LM)[-c(1,4)])
	#2nd stage of LP
	#obtain consistent estimates of materials and capital
	phihat <- fitted(LPfirststage)-as.matrix(L)%*%LPLabor
  LPphi <- phihat
	lagdata <- lagdata(idvar, cbind(Y, A, K, L, M, LPphi))
	names(lagdata) <- c("idvar", "Ycon", "Acon", "Kcon", "Lcon", "Mcon", "LPphicon",
		"Ylag", "Alag", "Klag", "Llag", "Mlag", "LPphilag")
	#LP Output
	LPmY <-  as.matrix(lagdata$Ycon-lagdata$Lcon*LPLabor)
	#LP Contemporary State Variables
  LPmX <- cbind(lagdata$Acon, lagdata$Kcon, lagdata$Mcon)
  #LP Lagged State Variables
  LPmlX <- cbind(lagdata$Alag, lagdata$Klag, lagdata$Mlag)
  #LP Contemporary phi estimates
  LPfitphi <- as.matrix(lagdata$LPphicon)
  #LP Lagged phi estimates
  LPfitlagphi <- as.matrix(lagdata$LPphilag)
  #LP Instruments 
  LPmZ <- cbind(lagdata$Acon, lagdata$Kcon, lagdata$Mlag)
  #LP estimates for Capital and Labor
  LPhat <- optim(par=LPinit, fn=function(b) LPobj(b, mY=LPmY, mX=LPmX, mlX=LPmlX, 
    mZ=LPmZ, fitphi=LPfitphi, fitlagphi=LPfitlagphi), gr=NULL, method="L-BFGS-B", 
    lower=c(0,0,0), upper=c(1,1,1))$par
  LPvec <- c(LPhat[1], LPhat[2], LPLabor, LPhat[3])
  omega <- Y-cbind(A, K, L, M)%*%as.matrix(as.numeric(LPvec))
  LPrho <- LPrho(bhat=LPhat, mY=LPmY, mX=LPmX, mlX=LPmlX, fitphi=LPfitphi, fitlagphi=LPfitlagphi)
  return(list(LPhat=LPvec, omega=omega, rho=LPrho))
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
LPrho <- function(bhat, mY, mX, mlX, fitphi, fitlagphi){
  b <- as.matrix(as.numeric(bhat))
  A <- fitphi-mX%*%b[1:ncol(mX)]
  B <- fitlagphi-mlX%*%b[1:ncol(mX)]
  step1 <- lm(A~B)
  step1param <- as.numeric(coef(step1))[2]
  return(step1param)
} 
