setwd('/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Functions')
require(quantreg)
require(dplyr)
require(pracma)
source('Tensors.R')
source('Posterior.R')
source('Moments.R')
source('Mstep.R')
source('Auxfuns.R')
source('omega.R')
# source('NLPFQR/FUN/Tensors.R')
# source('NLPFQR/FUN/Posterior.R')
# source('NLPFQR/FUN/Moments.R')
# source('NLPFQR/FUN/Mstep.R')
# source('NLPFQR/FUN/Auxfuns.R')
# source('NLPFQR/FUN/omega.R')
Production_EM <- function(ntau, idvar, timevar, Y, K, L, M, I, maxiter, draws){
	seed <- 123456
	set.seed(seed)
	#Logical checks for output, input
	#Output
	Y <- as.matrix(Y)
	#Capital
	K <- as.matrix(K)
	#Labor
	L <- as.matrix(L)
	#Materials
	M <- as.matrix(M)
	#Investment
	I <- as.matrix(I)
	#Number of Firms
	N <- length(unique(idvar))
	T <- unique(timevar)
	#Grid of Taus
	vectau <- seq(1/(ntau+1), ntau/(ntau+1), by=1/(ntau+1))
	#Initial estimate of productivity from ACF
	omega <- as.matrix(omega_est(idvar=idvar, timevar=timevar, Y=Y, K=K, L=L, M=M))
	print(summary(omega))
	#Reformatting data for initial optimization
	WTdata <- lagdata(idvar=idvar, X=omega)
	names(WTdata) <- c("idvar", "Wcon", "Wlag")
	W1data <- t0data(idvar=idvar, X=omega)
	names(W1data) <- c("idvar", "W1")
	YX <- translog(K=K, L=L, M=M, omega=omega)
	LX <- cbind(1, K, omega, K*omega, K^2, omega^2)
	MX <- cbind(1, K, omega, K*omega, K^2, omega^2)
	WX <- cbind(1, WTdata$Wlag, I(WTdata$Wlag^2), I(WTdata$Wlag^3))
	IX <- cbind(1, K, omega, K*omega, K^2, omega^2)
	dims <- list(Y=ncol(YX), L=ncol(LX), M=ncol(MX), W=ncol(WX), I=(ncol(LX)+1))
	############################################################################
	############ Initialization for EM Algorithm ###############################
	############################################################################
	#Initialize Paramater Values for Output
	resYinit <- matrix(rq((Y-omega)~YX-1, tau=vectau)$coef, nrow=dims$Y, ncol=ntau)
	#Initialize Paramater Values for Labor
	resLinit <- matrix(rq(L~LX-1, tau=vectau)$coef, nrow=dims$L, ncol=ntau)
	#Initialize Parameter Values for Labor
	resMinit <- matrix(rq(M~MX-1, tau=vectau)$coef, nrow=dims$M, ncol=ntau)
	#Initialize Parameter Values for Omega_{t}
	resWTinit <- matrix(rq(WTdata$Wcon~WX-1, tau=vectau)$coef, nrow=dims$W, ncol=ntau)
	#Initialize Parameter Values for Omega_{0}
	#There is an issue with selecting these initial parameters as the initial guess of unobservables is far off
	#Therefore i use the parameters of the proposal distribution as initial parameters
	w1mu <- mean(W1data$W1)
	w1var <- var(W1data$W1)
	resW1init <- c(0, .05)
	#Initialize Paramater Values for Investment (Nonlinear Regression Model)
	#This doesnt have to be the case, we can set a monotone (mode) restriction in omega
	#On the investment equation that sets restrictions on the parameter space
	resIcoef <- lm(I~IX-1)
	#Vector of Investment coefficients and estimate of  variance
	resIinit <- c(as.numeric(coef(resIcoef)), mean(resid(resIcoef)^2))
	#Initial Parameter Values for Laplace Parameters
	yb1init <- ybLinit <- 1
	lb1init <- lbLinit <- 1
	mb1init <- mbLinit <- 1
	wtb1init <- wtbLinit <- 1
	############################################################################
	#E-Step
	############################################################################
	#Variance Random Walk proposal
	varRW <- .05
	#Matrices for Storing Results
	resY <- array(0, c(maxiter, dims$Y, ntau))
	resL <- array(0, c(maxiter, dims$L, ntau))
	resM <- array(0, c(maxiter, dims$M, ntau))
	resWT <- array(0, c(maxiter, dims$W, ntau))
	resW1 <- array(0, c(maxiter, 2))
	resI <- array(0, c(maxiter, dims$I))
	yb1 <- array(0, c(maxiter, 1)); ybL <- array(0, c(maxiter, 1))
	lb1 <- array(0, c(maxiter, 1)); lbL <- array(0, c(maxiter, 1))
	mb1 <- array(0, c(maxiter, 1)); mbL <- array(0, c(maxiter, 1))
	wtb1 <- array(0, c(maxiter, 1)); wtbL <- array(0, c(maxiter, 1))
	#Begin EM Algorithm
	for (iter in 1:maxiter){
		print(sprintf("EM Iteration %i", iter))
		#Initial Parameter Values 
		resinit <- list(resYinit=resYinit, yb1init=yb1init, ybLinit=ybLinit, resLinit=resLinit, lb1init=lb1init, 
			lbLinit=lbLinit, resMinit=resMinit, mb1init=mb1init, mbLinit=mbLinit, 
			resWTinit=resWTinit, wtb1init=wtb1init, wtbLinit=wtbLinit, resW1init=resW1init, resIinit=resIinit)
		#Initial Guess for Unobservables
		chain <- omega
		oldpost <- posterior(idvar=idvar, Y=Y, K=K, L=L, M=M, I=I, omega=chain, vectau=vectau, par=resinit)
		for (j in 1:(draws+1)){
			Mseed <- seed+j
			#RW Update
			newchain <- chain+runif(length(chain), -0.5, 0.5)
			newpost <- posterior(idvar=idvar, Y=Y, K=K, L=L, M=M, I=I, omega=newchain, vectau=vectau, par=resinit)
			#Log acceptance probability
			loga <- min(newpost-oldpost, 0)
			#Create a set of N indices for acceptance rule
			Nindex <- as.logical(log(runif(N))<loga)
			#Create sequence of vectors for accepting each time period of each accepted firm draw
			NTindex <- rep(Nindex, as.data.frame(table(idvar))$Freq)
			#Update unobservables that satisfy acceptance rule according to the RW process
			chain[NTindex] <- newchain[NTindex]
			oldpost[Nindex] <- newpost[Nindex]
		}
		#This is the draws+1 draw of the unobservable in the stEM algorithm
		mat <- chain
		#Note that the traditional EM algorithm requires multiple of these draws so that
		#Mat would be a matrix with the number of columns equal to the number of accepted draws
		#We do not implement this version here because it is computationally costly
		######################################################################################
		#M STEP
		#######################################################################################
		######################################################################################
		#The optimization algorithm can either be a stochastic EM algorithm (Mdraws=1)
		#or a Monte Carlo EM algorithm (M large). The advantage of the former is less
		#Computational burden however requiring large EM iterations and a non-smooth convex objective function
		#Large M is computationally expensive but allows us to use a gradient descent algorithm
		#########################################################################################
		#Reformatting data for initial optimization
		WTdata <- lagdata(idvar=idvar, X=mat)
		names(WTdata) <- c("idvar", "Wcon", "Wlag")
		W1data <- t0data(idvar=idvar, X=mat)
		names(W1data) <- c("idvar", "W1")
		YX <- translog(K=K, L=L, M=M, omega=mat)
		LX <- cbind(1, K, mat, K*mat, K^2, mat^2)
		MX <- cbind(1, K, mat, K*mat, K^2, mat^2)
		WX <- cbind(1, WTdata$Wlag, I(WTdata$Wlag^2), I(WTdata$Wlag^3))
		IX <- cbind(1, K, mat, K*mat, K^2, mat^2)
		for (q in 1:ntau){
			resY[,,q][iter,] <- rq((Y-mat)~YX-1, tau=vectau[q])$coef
			resL[,,q][iter,] <- rq(L~LX-1, tau=vectau[q])$coef
			resM[,,q][iter,] <- rq(M~MX-1, tau=vectau[q])$coef
			resWT[,,q][iter,] <- rq(WTdata$Wcon~WX-1, tau=vectau[q])$coef
		}
		##########################################################################################
		#Investment Estimation
		resILM <- lm(I~IX-1)
		resIcoef <- as.numeric(coef(resILM))
		resIvar <- mean(resid(resILM)^2)
		resI[iter,] <- c(resIcoef, resIvar)
		#########################################################################################
		#Initial Productivity
		resW1[iter,] <- c(mean(W1data$W1), var(W1data$W1))
		###########################################################################################
		#Tail Parameters
		#Output
		yb <- expby(Y=Y, K=K, L=L, M=M, omega=mat, par1=resY[,,1][iter,], parL=resY[,,ntau][iter,])
		yb1[iter,] <- yb$b1
		ybL[iter,] <- yb$bL
		#Labor
		lb <- expbx(X=L, K=K, omega=mat, par1=resL[,,1][iter,], parL=resL[,,ntau][iter,])
		lb1[iter,] <- lb$b1
		lbL[iter,] <- lb$bL
		#Materials
		mb <- expbx(X=M, K=K, omega=mat, par1=resM[,,1][iter,], parL=resM[,,ntau][iter,])
		mb1[iter,] <- mb$b1
		mbL[iter,] <- mb$bL
		#Productivity
		wtb <- expbwt(omega=WTdata$Wcon, omegalag=WTdata$Wlag, par1=resWT[,,1][iter,], parL=resWT[,,ntau][iter,])
		wtb1[iter,] <- wtb$b1
		wtbL[iter,] <- wtb$bL
		#############################################################################################
		# Use Estimates as parameters in the MH Algorithm
		resYinit <- as.matrix(resY[iter,,])
		resLinit <- as.matrix(resL[iter,,])
		resMinit <- as.matrix(resM[iter,,])
		resWTinit <- as.matrix(resWT[iter,,])
		resW1init <- as.matrix(resW1[iter,])
		resIinit <- as.matrix(resI[iter,])
		yb1init <- yb1[iter,]; ybLinit <- ybL[iter,]
		lb1init <- lb1[iter,]; lbLinit <- lbL[iter,]
		mb1init <- mb1[iter,]; mbLinit <- mbL[iter,]
		wtb1init <- wtb1[iter,]; wtbLinit <- wtbL[iter,]
		################################################################################################
		}
	resYmat <- sapply(1:ntau, function(x) colMeans(resY[,,x][(maxiter/2):maxiter,]))
	resyb1bLmat <- c(mean(yb1[(maxiter/2):maxiter]), mean(ybL[(maxiter/2):maxiter]))
	resLmat <- sapply(1:ntau, function(x) colMeans(resL[,,x][(maxiter/2):maxiter,]))
	reslb1bLmat <- c(mean(lb1[(maxiter/2):maxiter]), mean(lbL[(maxiter/2):maxiter]))
	resMmat <- sapply(1:ntau, function(x) colMeans(resM[,,x][(maxiter/2):maxiter,]))
	resmb1bLmat <- c(mean(mb1[(maxiter/2):maxiter]), mean(mbL[(maxiter/2):maxiter]))
	resWTmat <- sapply(1:ntau, function(x) colMeans(resWT[,,x][(maxiter/2):maxiter,]))
	reswtb1bLmat <- c(mean(wtb1[(maxiter/2):maxiter]), mean(wtbL[(maxiter/2):maxiter]))
	resImat <- colMeans(resI[(maxiter/2):maxiter,])
	resW1mat <- colMeans(resW1[(maxiter/2):maxiter,])
	return(list(dims=dims, vectau=vectau, resYmat=resYmat, resLmat=resLmat, resMmat=resMmat, resWTmat=resWTmat, resW1mat=resW1mat, resImat=resImat, resyb1bLmat=resyb1bLmat, reslb1bLmat=reslb1bLmat, 
		resmb1bLmat=resmb1bLmat, reswtb1bLmat=reswtb1bLmat))
}




