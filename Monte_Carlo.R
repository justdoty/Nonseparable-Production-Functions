require(quantreg)
#This code extends the Monte Carlo experiment of Hu, Huang, and Sasaki (2019)
#Uses a random coefficient model instead of a linear model
nreps <- 1000
resmat <- array(0, dim=c(nreps, 6))
#Number of Firms:
N <- 500
T <- 200
#How many time periods to use in estimation
totalT <- 30
starttime <- T-totalT
#Specification of the median of production function parameters:
medconstant <- 0.1; medl <- 0.3; medk <- 0.2; medm <- 0.2; medu <- 0.1; medomega <- 0.1
#Specification of the median for persistent productivity
medrho <- 1; medrho0 <- 0.2
#Specifcation of the median of investment function parameters
iotamed0 <- 0.1; mediotak <- -0.05; mediotaI <- -0.04; mediotarho <- medrho 
#Depreciation of Capital
delta <- 0.1
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
	iotaI <- apply(zetadata, 2, function(x) beta(x, mediotak, 0.02))
	iotaw <- apply(zetadata, 2, function(x) beta(x, mediotarho, 0.1))
	iota0 <- apply(zetadata, 2, function(x) beta(x, iotamed0, 0.05))
	#Form productivity
	omgdata <- matrix(0, N, T)
	#Period 0 values of omega (GNR)
	omgdata0 <- matrix(runif(N, min=1, max=3),nrow=N,ncol=1)
	omgdata[,1] <- rho0[,1]+rho[,1]*omgdata0
	for (t in 2:T){
		omgdata[,t] <- rho0[,t]+rho[,t]*omgdata[,t-1]
	}
	#Form Capital and Investment rules using accumulation law and reduced form equation
	Capital <- matrix(0, N, T)
	Investment <- matrix(0, N, T)
	#Initial level of capital (GNR)
	Capital[,1] <- runif(N, min=11, max=400)
	Investment[,1] <- exp(iota0[,1]+iotak[,1]*Capital[,1]+iotaw[,1]*omgdata[,1])
	Capital[,2] <- (1-delta)*Capital[,1]+0.5*Investment[,1]
	Investment[,2] <- exp(iota0[,2]+iotak[,2]*Capital[,2]+iotaw[,2]*omgdata[,2])
	for (t in 3:T){
		Capital[,t] <- (1-delta)*Capital[,t-1]+0.5*Investment[,t-1]+0.5*Investment[,t-2]
		Investment[,t] <- exp(iota0[,t]+iotak[,t]*Capital[,t]+iotaw[,t]*omgdata[,t])
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
	#Estimation
	test <- rq(Output~Labor+Capital+Materials+Energy+Productivity, tau=0.5)
	resmat[n,] <- as.numeric(test$coefficients)
}
print(Sys.time()-sim.speed)
coef <- colMeans(resmat)
Bias <- coef-c(medconstant, medl, medk, medm, medu, medomega)
MSE <- apply(resmat, 2, var)+Bias^2











