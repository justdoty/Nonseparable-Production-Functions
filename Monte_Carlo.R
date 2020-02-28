require(quantreg)
#This code extends the Monte Carlo experiment of Hu, Huang, and Sasaki (2019)
#Uses a random coefficient model instead of a linear model
nreps <- 1000
tauseq <- seq(0.25, 0.75, by=0.25)
yresmat <- array(0, dim=c(nreps, 6, length(tauseq)))
lresmat <- array(0, dim=c(nreps, 3, length(tauseq)))
mresmat <- array(0, dim=c(nreps, 3, length(tauseq)))
uresmat <- array(0, dim=c(nreps, 3, length(tauseq)))
iresmat <- array(0, dim=c(nreps, 3, length(tauseq)))
kresmat <- array(0, dim=c(nreps, 3, length(tauseq)))
#Number of Firms:
N <- 500
T <- 100
#How many time periods to use in estimation
totalT <- 30
starttime <- T-totalT
#Specification of the median of production function parameters:
medconstant <- 0; medl <- 0.3; medk <- 0.3; medm <- 0.2; medu <- 0.1; medomega <- 0.1
#Specification of the median for persistent productivity
medrho <- 1; medrho0 <- 0
#Specifcation of the median of investment function parameters
iotamed0 <- 0; mediotak <- -0.2; mediotaI <- -0.1; mediotarho <- 1
#Capital Constant Coefficients###########################################
delta <- 0.2; kappaI <- 0.5
#This plots the coefficient functionals for the output equation
#Takes the inflection point as the median and returns a function that is concave
#in lower tails and convex in higher tails. Consider changing the degree of
#concavity/convexity
beta <- function(u, med, tol){
	min <- med-tol
	minu <- min(u)
	max <- med+tol
	maxu <- max(u)
	A <- cbind(c(minu, maxu), c(1,1))
	b <- as.matrix(c(min-(1/6*minu^3-1/4*minu^2), max-(1/6*maxu^3-1/4*maxu^2)))
	x <- solve(A)%*%b
	beta <- 1/6*u^3-1/4*u^2+x[1]*u+x[2]
	return(beta)
	
}
sim.speed <- Sys.time()
for (n in 1:nreps){
	print(n)
	set.seed(123456+n)
	#Specification for ex-post production shocks:
	etadata <- matrix(runif(N*T, 0, 1), nrow=N, ncol=T)
	#Functional coefficients for production function
	# betal <- beta(c(t(epsdata)), medl, 0.05)
	betaconstant <- apply(etadata, 2, function(x) beta(x, medconstant, 0.05))
	betal <- apply(etadata, 2, function(x) beta(x, medl, 0.05))
	betak <- apply(etadata, 2, function(x) beta(x, medk, 0.05))
	betam <- apply(etadata, 2, function(x) beta(x, medm, 0.05))
	betau <- apply(etadata, 2, function(x) beta(x, medu, 0.05))
	betarho <- apply(etadata, 2, function(x) beta(x, medomega, 0.05))
	#Specification for productivity innovation shocks:
	xidata <- matrix(runif(N*T, 0, 1), nrow=N, ncol=T)
	rho <- apply(xidata, 2, function(x) beta(x, medrho, 0.2))
	rho0 <- apply(xidata, 2, function(x) beta(x, medrho0, 0.1))
	#Specification for investment shocks:
	zetadata <- matrix(runif(N*T, 0, 1), nrow=N, ncol=T)
	#Functional coefficients for investment equation
	iotak <- apply(zetadata, 2, function(x) beta(x, mediotak, 0.02))
	iotaw <- apply(zetadata, 2, function(x) beta(x, mediotarho, 0.1))
	iota0 <- apply(zetadata, 2, function(x) beta(x, iotamed0, 0.1))
	#Specification for Capital Demand Shocks
	upsilondata <- matrix(runif(N*T, 0, 1), nrow=N, ncol=T)
	#Specification for Constant Coefficient in Capital Equation
	alphaK <- 1; sigmaK <- 0.1
	kappa0 <- alphaK+sigmaK*qnorm(upsilondata)
	#Form productivity
	omgdata <- matrix(0, N, T)
	#Period 0 values of omega (ACF)
	#Note final productivity levels are very sensitive to initial values
	#Which will also affect final input choices and investment/capital data
	omgdata0 <- matrix(rnorm(N, 0, 0.3),nrow=N,ncol=1)
	omgdata[,1] <- rho0[,1]+rho[,1]*omgdata0
	for (t in 2:T){
		omgdata[,t] <- rho0[,t]+rho[,t]*omgdata[,t-1]
	}
	#Form Capital and Investment rules using accumulation law and reduced form equation
	kdata <- matrix(0, N, T)
	lnidata <- matrix(0, N, T)
	#Initial level of capital (GNR)
	kdata[,1] <- exp(kappa0[,1])
	lnidata[,1] <- iota0[,1]+log(kdata[,1])*iotak[,1]+omgdata[,1]*iotaw[,1]
	for (t in 2:T){
		kdata[,t] <- exp(kappa0[,t])+(1-delta)*kdata[,t-1]+kappaI*exp(lnidata[,t-1])
		lnidata[,t] <- iota0[,t]+log(kdata[,t])*iotak[,t]+omgdata[,t]*iotaw[,t]
	}
	#Log Capital
	lnkdata <- log(kdata)
	#Reduced form equations for input demand functions (modelled as location scale):
	#Labor Input Demand
	epsdataL <- matrix(runif(N*T, 0, 1), nrow=N, ncol=T)
	mukL <- 0.5; muwL <- 1; alphaL <- 1; sigmaL <- 0.1
	mu0L <- alphaL+sigmaL*qnorm(epsdataL)
	lnldata <- mukL*lnkdata+muwL*omgdata+mu0L
	#Material Input Demand
	epsdataM <- matrix(runif(N*T, 0, 1), nrow=N, ncol=T)
	mukM <- 0.5; muwM <- 1; alphaM <- 1; sigmaM <- 0.1
	mu0M <- alphaM+sigmaM*qnorm(epsdataM)
	lnmdata <- mukM*lnkdata+muwM*omgdata+mu0M
	#Energy Input Demand
	epsdataU <- matrix(runif(N*T, 0, 1), nrow=N, ncol=T)
	mukU <- 0.5; muwU <- 1; alphaU <- 1; sigmaU <- 0.1
	mu0U <- alphaU+sigmaU*qnorm(epsdataU)
	lnudata<- mukU*lnkdata+muwU*omgdata+mu0U
	#Output Equation
	lnydata <- betaconstant+betal*lnldata+betak*lnkdata+betam*lnmdata+betau*lnudata+betarho*omgdata
	#Stack data across firms (all the data)
	Output <- c(t(lnydata[,(starttime+1):T]))
	Labor <- c(t(lnldata[,(starttime+1):T]))
	Capital <- c(t(lnkdata[,(starttime+1):T]))
	Materials <- c(t(lnmdata[,(starttime+1):T]))
	Energy <- c(t(lnudata[,(starttime+1):T]))
	Productivity <- c(t(omgdata[,(starttime+1):T]))
	Investment <- c(t(lnidata[,(starttime+1):T]))

	#Stack data across firms (lagged data)
	Output_Lag1 <- c(t(lnydata[,(starttime+1):(T-1)]))
	Labor_Lag1 <- c(t(lnldata[,(starttime+1):(T-1)]))
	Capital_Lag1 <- c(t(lnkdata[,(starttime+1):(T-1)]))
	Materials_Lag1 <- c(t(lnmdata[,(starttime+1):(T-1)]))
	Energy_Lag1 <- c(t(lnudata[,(starttime+1):(T-1)]))
	Productivity_Lag1 <- c(t(omgdata[,(starttime+1):(T-1)]))
	Investment_Lag1 <- c(t(lnidata[,(starttime+1):(T-1)]))

	#Stack data across firms (contermporaneous data)
	Output_Con <- c(t(lnydata[,(starttime+2):T]))
	Labor_Con <- c(t(lnldata[,(starttime+2):T]))
	Capital_Con <- c(t(lnkdata[,(starttime+2):T]))
	Materials_Con <- c(t(lnmdata[,(starttime+2):T]))
	Energy_Con <- c(t(lnudata[,(starttime+2):T]))
	Productivity_Con <- c(t(omgdata[,(starttime+2):T]))
	Investment_Con <- c(t(lnidata[,(starttime+2):T]))
	#Estimation
	for (q in 1:length(tauseq)){
		ytest <- rq(Output~Labor+Capital+Materials+Energy+Productivity, tau=tauseq[q])
		ltest <- rq(Labor~Capital+Productivity, tau=tauseq[q])
		mtest <- rq(Materials~Capital+Productivity, tau=tauseq[q])
		utest <- rq(Materials~Capital+Productivity, tau=tauseq[q])
		itest <- rq(Investment~Capital+Productivity, tau=tauseq[q])
		ktest <- rq(Capital_Con~Capital_Lag1+Investment_Lag1, tau=tauseq[q])
		yresmat[,,q][n,] <- as.numeric(ytest$coefficients)
		lresmat[,,q][n,] <- as.numeric(ltest$coefficients)
		mresmat[,,q][n,] <- as.numeric(mtest$coefficients)
		uresmat[,,q][n,] <- as.numeric(utest$coefficients)
		iresmat[,,q][n,] <- as.numeric(itest$coefficients)
		kresmat[,,q][n,] <- as.numeric(ktest$coefficients)
	}
}
print(Sys.time()-sim.speed)
yBias <- array(0, dim=c(length(tauseq), 6))
yMSE <- array(0, dim=c(length(tauseq), 6))

lBias <- array(0, dim=c(length(tauseq), 3))
lMSE <- array(0, dim=c(length(tauseq), 3))

mBias <- array(0, dim=c(length(tauseq), 3))
mMSE <- array(0, dim=c(length(tauseq), 3))

uBias <- array(0, dim=c(length(tauseq), 3))
uMSE <- array(0, dim=c(length(tauseq), 3))

iBias <- array(0, dim=c(length(tauseq), 3))
iMSE <- array(0, dim=c(length(tauseq), 3))

kBias <- array(0, dim=c(length(tauseq), 3))
kMSE <- array(0, dim=c(length(tauseq), 3))
for (q in 1:length(tauseq)){
	ytrue <- c(quantile(c(betaconstant), tauseq[q]), quantile(c(betal), tauseq[q]),
		quantile(c(betak), tauseq[q]), quantile(c(betam), tauseq[q]),
		quantile(c(betau), tauseq[q]), quantile(c(betarho), tauseq[q]))
	ycoef <- colMeans(yresmat[,,q])
	yBias[q,] <- ycoef-as.numeric(ytrue)
	yMSE[q,] <- colMeans((yresmat[,,q]-ytrue)^2)

	ltrue <- c(quantile(c(mu0L), tauseq[q]), mukL, muwL)
	lcoef <- colMeans(lresmat[,,q])
	lBias[q,] <- lcoef-as.numeric(ltrue)
	lMSE[q,] <- colMeans((lresmat[,,q]-ltrue)^2)

	mtrue <- c(quantile(c(mu0M), tauseq[q]), mukM, muwM)
	mcoef <- colMeans(mresmat[,,q])
	mBias[q,] <- mcoef-as.numeric(mtrue)
	mMSE[q,] <- colMeans((mresmat[,,q]-mtrue)^2)

	utrue <- c(quantile(c(mu0U), tauseq[q]), mukU, muwU)
	ucoef <- colMeans(uresmat[,,q])
	uBias[q,] <- ucoef-as.numeric(utrue)
	uMSE[q,] <- colMeans((uresmat[,,q]-utrue)^2)

	itrue <- c(quantile(c(iota0), tauseq[q]), quantile(c(iotak), tauseq[q]), quantile(c(iotaw), tauseq[q]))
	icoef <- colMeans(iresmat[,,q])
	iBias[q,] <- icoef-as.numeric(itrue)
	iMSE[q,] <- colMeans((iresmat[,,q]-itrue)^2)

	ktrue <- c(quantile(c(kappa0), tauseq[q]), (1-delta), kappaI)
	kcoef <- colMeans(kresmat[,,q])
	kBias[q,] <- kcoef-as.numeric(ktrue)
	kMSE[q,] <- colMeans((kresmat[,,q]-ktrue)^2)
	
}
#Bias reduced, but MSE in investment estimation is still high
#Bias in capital equation: Constant and Investment Lag (good/ok in Capital Lag)







