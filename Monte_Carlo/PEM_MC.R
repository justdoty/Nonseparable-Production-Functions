setwd('/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Monte_Carlo')
require(quantreg)
require(dplyr)
source('Tensors.R')
source('Posterior.R')
source('Moments.R')
source('Mstep.R')
source('Auxfuns.R')
source('LP_init.R')
# source('NLPFQR/FUN/Tensors.R')
# source('NLPFQR/FUN/Posterior.R')
# source('NLPFQR/FUN/Moments.R')
# source('NLPFQR/FUN/Mstep.R')
# source('NLPFQR/FUN/Auxfuns.R')
# source('NLPFQR/FUN/LP_init.R')
PEM_MC <- function(ntau, idvar, timevar, Y, K, L, M, I, maxiter, draws, Mdraws, seed){
	set.seed(seed)
	#We are assuming the researcher has access to a balanced panel. We can extend the case to unbalanced panel
	#By including possibly conditioning the productivity process on survival probability following OP (1994)
	#Or using the methods of Arellano and Bonhomme (2017)
	#Logical checks for output, input, age data
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
	#Grid of Taus
	vectau <- seq(1/(ntau+1), ntau/(ntau+1), by=1/(ntau+1))
	TFP <- as.matrix(LP_est(idvar=idvar, timevar=timevar, Y=Y, K=K, L=L, M=M))
	#Reformatting data for initial optimization
	WTdata <- lagdata(idvar=idvar, X=TFP)
	names(WTdata) <- c("idvar", "Ucon", "Ulag")
	W1data <- t0data(idvar=idvar, X=TFP)
	names(W1data) <- c("idvar", "U1")
	############################################################################
	############ Initialization for EM Algorithm ###############################
	############################################################################
	#Initialize Paramater Values for Output
	resYinit <- matrix(rq(Y~cbind(K, L, M, TFP), tau=vectau)$coef, nrow=5, ncol=ntau)
	#Initialize Parameter Values for Omega_{t}
	resWTinit <- matrix(rq(WTdata$Ucon~WTdata$Ulag, tau=vectau)$coef, nrow=2, ncol=ntau)
	#Initialize Parameter Values for Omega_{0}
	#There is an issue with selecting these initial parameters as the initial guess of unobservables is far off
	#Therefore i use the parameters of the proposal distribution as initial parameters
	w1mu <- mean(W1data$U1)
	w1var <- var(W1data$U1)
	resW1init <- c(0, .05)
	#Initialize Paramater Values for Labor
	resLcoef <- lm(L~cbind(K, TFP))
	resLinit <- c(as.numeric(coef(resLcoef)), mean(resid(resLcoef)^2))
	#Initialize Parameter Values for Materials
	resMcoef <- lm(M~cbind(K, TFP))
	resMinit <- c(as.numeric(coef(resMcoef)), mean(resid(resMcoef)^2))
	#Initialize Paramater Values for Investment 
	resIcoef <- lm(I~cbind(K, TFP))
	resIinit <- c(as.numeric(coef(resIcoef)), mean(resid(resIcoef)^2))
	#Initial Parameter Values for Laplace Parameters
	yb1init <- ybLinit <- 1
	wtb1init <- wtbLinit <- 1
	############################################################################
	############################################################################
	############################## EM Algorithm ################################
	############################################################################
	############################################################################
	#Variance Random Walk proposal
	varRW <- .05
	#Matrices for Storing Results
	resY <- array(0, c(maxiter, 5, ntau))
	resWT <- array(0, c(maxiter, 2, ntau))
	resL <- array(0, c(maxiter, 4))
	resM <- array(0, c(maxiter, 4))
	resI <- array(0, c(maxiter, 4))
	resW1 <- array(0, c(maxiter, 2))
	yb1 <- array(0, c(maxiter, 1)); ybL <- array(0, c(maxiter, 1))
	wtb1 <- array(0, c(maxiter, 1)); wtbL <- array(0, c(maxiter, 1))
	#Begin EM Algorithm
	for (iter in 1:maxiter){
		print(sprintf("EM Iteration %i", iter))
		#This sets how many Metropolis-Hastings "burn in" samples
		r <- -draws
		#Initialize Matrix for Averaging over Accepted Draws
		mat <- matrix(0, length(Y), Mdraws)
		#Initial Parameter Values 
		resinit <- list(resYinit=resYinit, yb1init=yb1init, ybLinit=ybLinit, resLinit=resLinit, resMinit=resMinit, 
			resWTinit=resWTinit, wtb1init=wtb1init, wtbLinit=wtbLinit, resW1init=resW1init, resIinit=resIinit)
		#Initial Guess for Unobservables
		un <- as.matrix(rnorm(length(Y), 0, sqrt(varRW)))
		#Calculate posterior density at this guess
		un_dens <- posterior(data=cbind(idvar, Y, K, L, M, I, un), vectau=vectau, par=resinit)
		for (j in r:Mdraws){
			Mseed <- seed+j
			set.seed(Mseed)
			#Proposal Distribution
			try_un <- as.matrix(rnorm(length(Y), 0, sqrt(varRW)))
			#Random Walk update
			y <- un+try_un
			loglik <- posterior(data=cbind(idvar, Y, K, L, M, I, y), vectau=vectau, par=resinit)
			#Log acceptance probability
			loga <- loglik-un_dens
			#Create a set of N indices for acceptance rule
			Nindex <- as.logical(log(runif(N))<loga)
			#Create sequence of vectors for accepting each time period of each accepted firm draw
			NTindex <- rep(Nindex, as.data.frame(table(idvar))$Freq)
			#Update unobservables that satisfy acceptance rule according to the RW process
			un[NTindex] <- y[NTindex]
			un_dens <- loglik
			if (j>0){
				mat[,j] <- un
			}
		}
		#The optimization algorithm can either be a stochastic EM algorithm (Mdraws=1)
		#or a Monte Carlo EM algorithm (M large). The advantage of the former is less
		#Computational burden however requiring large EM iterations and a non-smooth convex objective function
		#Large M is computationally expensive but allows us to use a gradient descent algorithm
		#Reformatting data for initial optimization
		WTdata <- lagdata(idvar=idvar, X=mat)
		names(WTdata) <- c("idvar", "Ucon", "Ulag")
		W1data <- t0data(idvar=idvar, X=mat)
		names(W1data) <- c("idvar", "U1")
		for (q in 1:ntau){
			print(sprintf("Quantile Mstep %i", q))
			if (Mdraws==1){
				resY[,,q][iter,] <- rq(Y~cbind(K, L, M, mat), tau=vectau[q])$coef
				resWT[,,q][iter,] <- rq(WTdata$Ucon~WTdata$Ulag, tau=vectau[q])$coef
			} else {
				print('Error: Still working on this part')
			}
		}
		##########################################################################################
		#Investment Estimation
		if (Mdraws==1){
			#Labor
			resLLM <- lm(L~cbind(K, mat))
			resLcoef <- as.numeric(coef(resLLM))
			resLvar <- mean(resid(resLLM)^2)
			#Materials
			resMLM <- lm(M~cbind(K, mat))
			resMcoef <- as.numeric(coef(resMLM))
			resMvar <- mean(resid(resMLM)^2)
			#Investment
			resILM <- lm(I~cbind(K, mat))
			resIcoef <- as.numeric(coef(resILM))
			resIvar <- mean(resid(resILM)^2)
			
		} else {
			print('Error: Still working on this part :)')
			}
		resL[iter,] <- c(resLcoef, resLvar)
		resM[iter,] <- c(resMcoef, resMvar)
		resI[iter,] <- c(resIcoef, resIvar)
		#########################################################################################
		#Initial Productivity
		resW1[iter,] <- c(mean(W1data$U1), var(W1data$U1))
		###########################################################################################
		#Tail Parameters
		#Output
		yb <- expb(X=cbind(Y, K, L, M), U=mat, par1=resY[,,1][iter,], parL=resY[,,ntau][iter,], WT=2)
		yb1[iter,] <- yb$b1
		ybL[iter,] <- yb$bL
		#Productivity
		wtb <- expb(idvar=idvar, U=mat, par1=resWT[,,1][iter,], parL=resWT[,,ntau][iter,], WT=1)
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
		wtb1init <- wtb1[iter,]; wtbLinit <- wtbL[iter,]
		################################################################################################
		}
	resYmat <- sapply(1:ntau, function(x) colMeans(resY[,,x][(maxiter/2):maxiter,]))
	resyb1bLmat <- c(mean(yb1[(maxiter/2):maxiter]), mean(ybL[(maxiter/2):maxiter]))
	resLmat <- colMeans(resL[(maxiter/2):maxiter,])
	resMmat <- colMeans(resM[(maxiter/2):maxiter,])
	resWTmat <- sapply(1:ntau, function(x) colMeans(resWT[,,x][(maxiter/2):maxiter,]))
	reswtb1bLmat <- c(mean(wtb1[(maxiter/2):maxiter]), mean(wtbL[(maxiter/2):maxiter]))
	resImat <- colMeans(resI[(maxiter/2):maxiter,])
	resW1mat <- colMeans(resW1[(maxiter/2):maxiter,])
	return(list(vectau=vectau, resYmat=resYmat, resLmat=resLmat, resMmat=resMmat, resWTmat=resWTmat, resW1mat=resW1mat, resImat=resImat, resyb1bLmat=resyb1bLmat, reswtb1bLmat=reswtb1bLmat))
}




