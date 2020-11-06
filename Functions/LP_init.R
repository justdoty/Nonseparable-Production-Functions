setwd('/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Functions')
#This function is used to do a basic LP estimation of TFP to use as initial model estimates
#and initial draw for productivity in the stEM algorithm
source('Auxfuns.R')
#LP GMM Criterion Function
#x denotes the exogeneous and endogeneous regressors
#z denotes the instrumets
LP_GMM <- function(x, z, b){
  Omega1 <- x[,1]-b[1]*x[,2]-b[2]*x[,3]
  Omega2 <- x[,4]-b[1]*x[,5]-b[2]*x[,6]
  conc <- fitted(lm(Omega1~Omega2))
  resid <- x[,1]-b[1]*x[,2]-b[2]*x[,3]-conc
  Obj <- sum((z*array(data=resid, dim=dim(z)))^2)
  return(Obj)
}
LP_est <- function(idvar, timevar, Y, K, L, M){
	#First Stage of LP 
	#Obtain consistent estimates of labor and age
	L1st <- lm(L~K+M)
	Y1st <- lm(Y~K+M)
	Lfit <- L-fitted(L1st)
	Yfit <- Y-fitted(Y1st)
	stage1 <- as.numeric(coef(lm(Yfit~Lfit-1)))
	#2nd stage of LP
	#obtain consistent estimates of materials and capital
	phi <- Y-stage1[1]*L
	lagdata <- lagdata(idvar, cbind(Y, K, L, M, phi))
	names(lagdata) <- c("idvar", "Ycon", "Kcon", "Lcon", "Mcon", "phicon",
		"Ylag", "Klag", "Llag", "Mlag", "philag")
	X <- as.matrix(cbind(lagdata$phicon, lagdata$Kcon, lagdata$Mcon, lagdata$philag, lagdata$Klag, 
		lagdata$Mlag))
	Z <- as.matrix(cbind(lagdata$Kcon, lagdata$Mlag))
	obj.fn_LP <- function(b){
      momi <- LP_GMM(x=X, z=Z, b)
      return(momi)
    }
    LP1 <- stage1[1]
    LP2 <- optim(par=c(0.5, 0.5), obj.fn_LP)$par
    TFP <- Y-cbind(L, K, M)%*%c(LP1, LP2)
    return(TFP)

}