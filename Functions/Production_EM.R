setwd('/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Functions')
require(quantreg)
require(dplyr)
source('Tensors.R')
source('Posterior.R')
source('Moments.R')
source('Mstep.R')
source('Auxfuns.R')
source('LP_init.R')
Production_EM <- function(ntau, idvar, timevar, Y, K, L, M, I, maxiter, draws, Mdraws, seed, MH){
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
	WTdata <- lagdata(idvar, TFP)
	names(WTdata) <- c("idvar", "Ucon", "Ulag")
	W1data <- t0data(idvar, TFP)
	names(W1data) <- c("idvar", "U1")
	############################################################################
	############ Initialization for EM Algorithm ###############################
	############################################################################
	#Initialize Paramater Values for Output
	resYinit <- matrix(rq(Y~tensor.prod(MH$MY, cbind(K, L, M, TFP), norm=0)-1, tau=vectau)$coef, nrow=prod(MH$MY+1), ncol=ntau)
	#Initialize Paramater Values for Labor
	resLinit <- matrix(rq(L~tensor.prod(MH$ML, cbind(K, TFP), norm=0)-1, tau=vectau)$coef, nrow=prod(MH$ML+1), ncol=ntau)
	#Initialize Parameter Values for Labor
	resMinit <- matrix(rq(M~tensor.prod(MH$MM, cbind(K, TFP), norm=0)-1, tau=vectau)$coef, nrow=prod(MH$MM+1), ncol=ntau)
	#Initialize Parameter Values for Omega_{t}
	resWTinit <- matrix(rq(WTdata$Ucon~tensor.prod(MH$MW, WTdata$Ulag, norm=0)-1, tau=vectau)$coef, nrow=prod(MH$MW+1), ncol=ntau)
	#Initialize Parameter Values for Omega_{0}
	w1mu <- mean(W1data$U1)
	w1var <- var(W1data$U1)
	resW1init <- c(0, .05^2)
	#Initialize Paramater Values for Investment (Nonlinear Regression Model)
	resIcoef <- lm(I~tensor.prod(MH$MI, cbind(K, TFP), norm=0)-1)
	#Vector of Investment coefficients and estimate of  variance
	resIinit <- c(as.numeric(coef(resIcoef)), mean(resid(resIcoef)^2))
	#Initial Parameter Values for Laplace Parameters
	yb1init <- ybLinit <- 1
	lb1init <- lbLinit <- 1
	mb1init <- mbLinit <- 1
	wtb1init <- wtbLinit <- 1
	############################################################################
	############################################################################
	############################## EM Algorithm ################################
	############################################################################
	############################################################################
	#Variance Random Walk proposal
	varRW <- .05
	#Matrices for Storing Results
	resY <- array(0, c(maxiter, prod(MH$MY+1), ntau))
	resL <- array(0, c(maxiter, prod(MH$ML+1), ntau))
	resM <- array(0, c(maxiter, prod(MH$MM+1), ntau))
	resWT <- array(0, c(maxiter, prod(MH$MW+1), ntau))
	resW1 <- array(0, c(maxiter, 2))
	resI <- array(0, c(maxiter, (prod(MI+1)+1)))
	yb1 <- array(0, c(maxiter, 1)); ybL <- array(0, c(maxiter, 1))
	lb1 <- array(0, c(maxiter, 1)); lbL <- array(0, c(maxiter, 1))
	mb1 <- array(0, c(maxiter, 1)); mbL <- array(0, c(maxiter, 1))
	wtb1 <- array(0, c(maxiter, 1)); wtbL <- array(0, c(maxiter, 1))
	#Begin EM Algorithm
	for (iter in 1:maxiter){
		print(sprintf("EM Iteration %i", iter))
		#This sets how many Metropolis-Hastings "burn in" samples
		r <- -draws
		#Initialize Matrix for Averaging over Accepted Draws
		mat <- matrix(0, length(Y), Mdraws)
		#Initial Parameter Values 
		resinit <- list(resYinit=resYinit, yb1init=yb1init, ybLinit=ybLinit, resLinit=resLinit, lb1init=lb1init, 
			lbLinit=lbLinit, resMinit=resMinit, mb1init=mb1init, mbLinit=mbLinit, 
			resWTinit=resWTinit, wtb1init=wtb1init, wtbLinit=wtbLinit, resW1init=resW1init, resIinit=resIinit)
		#Initial Guess for Unobservables
		un <- as.matrix(rnorm(length(Y), 0, sqrt(varRW)))
		#Calculate posterior density at this guess
		un_dens <- posterior(data=cbind(idvar, Y, K, L, M, I, un), MH=MH, vectau=vectau, par=resinit)
		for (j in r:Mdraws){
			Mseed <- seed+j
			set.seed(Mseed)
			#Proposal Distribution
			try_un <- as.matrix(rnorm(length(Y), 0, sqrt(varRW)))
			#Random Walk update
			y <- un+try_un
			loglik <- posterior(data=cbind(idvar, Y, K, L, M, I, y), MH=MH, vectau=vectau, par=resinit)
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
		WTdata <- lagdata(idvar, mat)
		names(WTdata) <- c("idvar", "Ucon", "Ulag")
		W1data <- t0data(idvar, mat)
		names(W1data) <- c("idvar", "U1")
		for (q in 1:ntau){
			print(sprintf("Quantile Mstep %i", q))
			if (Mdraws==1){
				resY[,,q][iter,] <- rq(Y~tensor.prod(MH$MY, cbind(K, L, M, mat), norm=0)-1, tau=vectau[q])$coef
				resL[,,q][iter,] <- rq(L~tensor.prod(MH$ML, cbind(K, mat), norm=0)-1, tau=vectau[q])$coef
				resM[,,q][iter,] <- rq(M~tensor.prod(MH$MM, cbind(K, mat), norm=0)-1, tau=vectau[q])$coef
				resWT[,,q][iter,] <- rq(WTdata$Ucon~tensor.prod(MH$MW, WTdata$Ulag, norm=0)-1, tau=vectau[q])$coef
			} else {
				print('Error: Still working on this part :)')
				#TO DO: ALLOW FOR DIFFERENT M SIZE TO USE AS DIFFERENT INITIAL VALUES 
				#FOR THE EM ALGORITHM CHOOSING THE VALUES THAT CORRESPOND TO THE HIGHEST
				#LIKELIHOOD VALUE
				#Initialize Gradient Descent Algorithms
			    #Tolerance for gradient (close to zero)
			 #    tol <- 1e-2
			 #    #Stepsize for gradient
			 #    stepsize <- 100
			 #    #Number of gradient descent iterations
			 #    nsteps <- 800
			 #    nbatch <- 1
				# resY[,,q][iter,] <- sgradfunq(idvar=idvar, X=cbind(Y, K, L, M, A), U=mat, MH=MH$MY, 
				# 	init=resYinit[,q], tau=vectau[q], tol=tol, stepsize=stepsize, nsteps=nsteps, nbatch=nbatch, seed=seed, WT=2)$coef
				# resL[,,q][iter,] <- sgradfunq(idvar=idvar, X=cbind(L, K, A), U=mat, MH=MH$ML, 
				# 	init=resLinit[,q], tau=vectau[q], tol=tol, stepsize=stepsize, nsteps=nsteps, nbatch=nbatch, seed=seed, WT=2)$coef
				# resM[,,q][iter,] <- sgradfunq(idvar=idvar, X=cbind(M, K, A), U=mat, MH=MH$MM, 
				# 	init=resMinit[,q], tau=vectau[q], tol=tol, stepsize=stepsize, nsteps=nsteps, nbatch=nbatch, seed=seed, WT=2)$coef
				# resWT[,,q][iter,] <- sgradfunq(idvar=idvar, X=A, U=mat, MH=MH$MW, init=resWTinit[,q], tau=vectau[q], 
				# 	tol=tol, stepsize=stepsize, nsteps=nsteps, nbatch=nbatch, seed=seed, WT=1)$coef
			}
		}
		##########################################################################################
		#Investment Estimation
		if (Mdraws==1){
			resILM <- lm(I~tensor.prod(MI, cbind(K, mat), norm=0)-1)
			resIcoef <- as.numeric(coef(resILM))
			resIvar <- mean(resid(resILM)^2)
		} else {
			print('Error: Still working on this part :)')
			# resIcoef <- sgradfune(idvar=idvar, X=cbind(I, K, A), U=mat, MH=MI, 
			# 		init=resIinit[-length(resIinit)], tol=tol, stepsize=stepsize, nsteps=nsteps, nbatch=1, seed=seed)$coef
			# #Estimate of Variance
			# resIvar <- mean((apply(mat, 2, function(u) sum(I-tensor.prod(MI, cbind(K, A, u))%*%resIcoef)))^2)
			}
		resI[iter,] <- c(resIcoef, resIvar)
		#########################################################################################
		#Initial Productivity
		resW1[iter,] <- c(mean(W1data$U1), var(W1data$U1))
		###########################################################################################
		#Tail Parameters
		#Output
		yb <- expb(X=cbind(Y, K, L, M), U=mat, MH=MH$MY, par1=resY[,,1][iter,], parL=resY[,,ntau][iter,], WT=2)
		yb1[iter,] <- yb$b1
		ybL[iter,] <- yb$bL
		#Labor
		lb <- expb(X=cbind(L, K), U=mat, MH=MH$ML, par1=resL[,,1][iter,], parL=resL[,,ntau][iter,], WT=2)
		lb1[iter,] <- lb$b1
		lbL[iter,] <- lb$bL
		#Materials
		mb <- expb(X=cbind(M, K), U=mat, MH=MH$MM, par1=resL[,,1][iter,], parL=resL[,,ntau][iter,], WT=2)
		mb1[iter,] <- mb$b1
		mbL[iter,] <- mb$bL
		#Productivity
		wtb <- expb(idvar=idvar, U=mat, MH=MH$MW, par1=resWT[,,1][iter,], parL=resWT[,,ntau][iter,], WT=1)
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
	return(list(resYmat=resYmat, resLmat=resLmat, resMmat=resMmat, resWTmat=resWTmat, resW1mat=resW1mat, resImat=resImat, resyb1bLmat=resyb1bLmat, reslb1bLmat=reslb1bLmat, 
		resmb1bLmat=resmb1bLmat, reswtb1bLmat=reswtb1bLmat))
}




