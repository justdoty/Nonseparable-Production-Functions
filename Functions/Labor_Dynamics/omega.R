#This function is used to do a basic LP estimation of productivity to use as initial model estimates
#and initial draw for productivity in the stEM algorithm
#Function that estimates productivity (omega)
omega_est <- function(idvar, timevar, Y, A, K, L, M){
  idcon <- duplicated(idvar)
  idlag <- duplicated(idvar, fromLast=TRUE)
	regvars <- data.frame(reg1=L, reg2=K, reg3=M, reg4=K^2, reg5=M^2)
  #First Stage 
	LPfirststage <- lm(Y~as.matrix(regvars[, grepl('reg', colnames(regvars))]))
  LPresid <- resid(LPfirststage)
  firstcoef <- as.numeric(coef(LPfirststage))
  LPLabor <-  firstcoef[2]
	#Initial Values for search
	LP_LM <- lm(Y~K+L+M)
	# LPinit <- as.numeric(coef(LP_LM)[-c(1,3)])
  LPinit <- c(firstcoef[3], firstcoef[4])+rnorm(2, 0, 0.01)
	#2nd stage of LP
	#obtain consistent estimates of materials and capital
	phihat <- fitted(LPfirststage)-as.matrix(L)%*%LPLabor
  LPphi <- phihat
	#LP Output
	LPmY <-  as.matrix(Y-L*LPLabor)[idcon,]
	#LP Contemporary State Variables
  LPmX <- cbind(K, M)[idcon,]
  #LP Lagged State Variables
  LPmlX <- cbind(K, M)[idlag,]
  #LP Contemporary phi estimates
  LPfitphi <- as.matrix(LPphi)[idcon,]
  #LP Lagged phi estimates
  LPfitlagphi <- as.matrix(LPphi)[idlag,]
  #LP Instruments 
  LPmZ <- cbind(K[idcon], M[idlag])
  #LP estimates for Capital and Labor
  LPhat <- optim(par=LPinit, fn=function(b) LP_Lambda(b, mY=LPmY, mX=LPmX, mlX=LPmlX, 
    fitphi=LPfitphi, fitlagphi=LPfitlagphi), method="BFGS")$par
  LPvec <- c(LPhat[1], LPLabor, LPhat[2])
  omega <- Y-cbind(K, L, M)%*%LPvec
  return(list(LPhat=LPvec, omega=omega, LPresid=LPresid))
}
############################################################################################
#Functions for Estimating LP Coefficients
###########################################################################################
#Function that defines the residuals
LP_Lambda <- function(b, mY, mX, mlX, fitphi, fitlagphi){
  b <- as.matrix(as.numeric(b))
  A <- fitphi-mX%*%b[1:ncol(mX)]
  B <- fitlagphi-mlX%*%b[1:ncol(mX)]
  omegapol <- poly(B, degree=3, raw=TRUE)
  omegapol <- cbind(1, omegapol)
  gb <- solve(crossprod(omegapol), tol=1e-100)%*%t(omegapol)%*%A
  wfit <- omegapol%*%gb
  # step1 <- lm(A~B)
  # step1param <- as.numeric(coef(step1))
  # wfit <- cbind(1, B)%*%step1param
  resid <- mY-mX%*%b[1:ncol(mX)]-wfit
  crit <- crossprod(resid)
  return(crit)
} 
#LP GMM objective function
LPobj <- function(b, mY, mX, mlX, mZ, fitphi, fitlagphi){
  xifit <- LP_Lambda(b=b, mY=mY, mX=mX, mlX=mlX, fitphi=fitphi, fitlagphi=fitlagphi)
  mW <- solve(crossprod(mZ))/nrow(mZ)
  go <- t(crossprod(mZ, xifit))%*%mW%*%(crossprod(mZ, xifit))
  return(go)
}
