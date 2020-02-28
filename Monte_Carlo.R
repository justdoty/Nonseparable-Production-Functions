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
#Number of Firms:
N <- 500
T <- 100
#How many time periods to use in estimation
totalT <- 30
starttime <- T-totalT
#Specification of the median of production function parameters:
medconstant <- 0.1; medl <- 0.3; medk <- 0.2; medm <- 0.2; medu <- 0.1; medomega <- 0.1
#Specification of the median for persistent productivity
medrho <- 1; medrho0 <- 0
#Specifcation of the median of investment function parameters
iotamed0 <- -0.7; mediotak <- -0.2; mediotaI <- -0.1; mediotarho <- 1
#Depreciation of Capital
delta <- 0.2
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
	Capital <- matrix(0, N, T)
	Investment <- matrix(0, N, T)
	#Initial level of capital (GNR)
	Capital[,1] <- runif(N)
	Investment[,1] <- iota0[,1]+Capital[,1]*iotak[,1]+omgdata[,1]*iotaw[,1]
	Capital[,2] <- (1-delta)*Capital[,1]+0.5*exp(Investment[,1])
	Investment[,2] <- iota0[,2]+Capital[,2]*iotak[,2]+omgdata[,2]*iotaw[,2]
	for (t in 3:T){
		Capital[,t] <- (1-delta)*Capital[,t-1]+0.5*exp(Investment[,t-1])+0.5*exp(Investment[,t-2])
		Investment[,t] <- iota0[,t]+Capital[,t]*iotak[,t]+omgdata[,t]*iotaw[,t]
	}
	#Log Capital
	capital <- log(Capital)
	#Reduced form equations for input demand functions:
	#Labor Input Demand
	epsdataL <- matrix(runif(N*T, 0, 1), nrow=N, ncol=T)
	mukL <- 0.5; muwL <- 1; alphaL <- 1; sigmaL <- 0.1
	beta0L <- alphaL+sigmaL*qnorm(epsdataL)
	labor <- mukL*capital+muwL*omgdata+beta0L
	#Material Input Demand
	epsdataM <- matrix(runif(N*T, 0, 1), nrow=N, ncol=T)
	mukM <- 0.5; muwM <- 1; alphaM <- 1; sigmaM <- 0.1
	beta0M <- alphaM+sigmaM*qnorm(epsdataM)
	materials <- mukM*capital+muwM*omgdata+beta0M
	#Energy Input Demand
	epsdataU <- matrix(runif(N*T, 0, 1), nrow=N, ncol=T)
	mukU <- 0.5; muwU <- 1; alphaU <- 1; sigmaU <- 0.1
	beta0U <- alphaU+sigmaU*qnorm(epsdataU)
	energy <- mukU*capital+muwU*omgdata+beta0U
	#Output Equation
	output <- betaconstant+betal*labor+betak*capital+betam*materials+betau*energy+betarho*omgdata
	#Stack data across firms (all the data)
	Output <- c(t(output[,(starttime+1):T]))
	Labor <- c(t(labor[,(starttime+1):T]))
	Capital <- c(t(capital[,(starttime+1):T]))
	Materials <- c(t(materials[,(starttime+1):T]))
	Energy <- c(t(energy[,(starttime+1):T]))
	Productivity <- c(t(omgdata[,(starttime+1):T]))
	Investment <- c(t(Investment[,(starttime+1):T]))
	#Estimation
	for (q in 1:length(tauseq)){
		ytest <- rq(Output~Labor+Capital+Materials+Energy+Productivity, tau=tauseq[q])
		ltest <- rq(Labor~Capital+Productivity, tau=tauseq[q])
		mtest <- rq(Materials~Capital+Productivity, tau=tauseq[q])
		utest <- rq(Materials~Capital+Productivity, tau=tauseq[q])
		itest <- rq(Investment~Capital+Productivity, tau=tauseq[q])
		yresmat[,,q][n,] <- as.numeric(ytest$coefficients)
		lresmat[,,q][n,] <- as.numeric(ltest$coefficients)
		mresmat[,,q][n,] <- as.numeric(mtest$coefficients)
		uresmat[,,q][n,] <- as.numeric(utest$coefficients)
		iresmat[,,q][n,] <- as.numeric(itest$coefficients)
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
for (q in 1:length(tauseq)){
	ytrue <- c(quantile(c(betaconstant), tauseq[q]), quantile(c(betal), tauseq[q]),
		quantile(c(betak), tauseq[q]), quantile(c(betam), tauseq[q]),
		quantile(c(betau), tauseq[q]), quantile(c(betarho), tauseq[q]))
	ycoef <- colMeans(yresmat[,,q])
	yBias[q,] <- ycoef-as.numeric(ytrue)
	yMSE[q,] <- colMeans((yresmat[,,q]-ytrue)^2)

	ltrue <- c(quantile(c(beta0L), tauseq[q]), mukL, muwL)
	lcoef <- colMeans(lresmat[,,q])
	lBias[q,] <- lcoef-as.numeric(ltrue)
	lMSE[q,] <- colMeans((lresmat[,,q]-ltrue)^2)

	mtrue <- c(quantile(c(beta0M), tauseq[q]), mukM, muwM)
	mcoef <- colMeans(mresmat[,,q])
	mBias[q,] <- mcoef-as.numeric(mtrue)
	mMSE[q,] <- colMeans((mresmat[,,q]-mtrue)^2)

	utrue <- c(quantile(c(beta0U), tauseq[q]), mukU, muwU)
	ucoef <- colMeans(uresmat[,,q])
	uBias[q,] <- ucoef-as.numeric(utrue)
	uMSE[q,] <- colMeans((uresmat[,,q]-utrue)^2)

	itrue <- c(quantile(c(iota0), tauseq[q]), quantile(c(iotak), tauseq[q]), quantile(c(iotaw), tauseq[q]))
	icoef <- colMeans(iresmat[,,q])
	iBias[q,] <- icoef-as.numeric(itrue)
	iMSE[q,] <- colMeans((iresmat[,,q]-itrue)^2)
}
#TO DO: FIX INVESTMENT ESTIMATION RESULTS







