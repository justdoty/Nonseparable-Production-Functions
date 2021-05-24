# setwd('/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Functions')
require(quantreg)
require(dplyr)
require(pracma)
# source('Tensors.R')
# source('Posterior.R')
# source('Mstep.R')
# source('Auxfuns.R')
# source('omega.R')
source('NLPFQR/FUN/Tensors.R')
source('NLPFQR/FUN/Posterior.R')
source('NLPFQR/FUN/Mstep.R')
source('NLPFQR/FUN/Auxfuns.R')
source('NLPFQR/FUN/omega.R')
Production_EM <- function(ntau, idvar, timevar, Y, A, K, L, M, I, method="trans", maxiter, draws){
	if (all(method!=c("cobbN", "transN", "cobb", "trans", "hermite"))){
		stop("Production Function must be specified as either cobbN (Hicks-Neutral Cobb Douglas),
			transN (Hicks-Neutral Translog), cobb (NonHicks-Neutral Cobb Douglass), trans (NonHicks-Neutral Translog),
			or hermite (tensor product hermite polynomial)")
	}
	seed <- 123456
	set.seed(seed)
	#Number of Firms
	N <- length(unique(idvar))
	T <- unique(timevar)
	#Grid of Taus
	vectau <- seq(1/(ntau+1), ntau/(ntau+1), by=1/(ntau+1))
	#Initial estimate of productivity from ACF plus random noise
	# omegainit <- jitter(as.matrix(omega_est(idvar=idvar, timevar=timevar, Y=Y, A=A, K=K, L=L, M=M)$omega), 50)
	# omegainit <- rnorm(length(idvar), mean=0, sd=sd(omegainit))
	omegainit <- as.matrix(rnorm(length(idvar)))
	#Reformatting data for initial optimization
	lagdata <- lagdata(idvar, cbind(Y, A, K, L, M, I, omegainit))
	names(lagdata) <- c("idvar", "Ycon", "Acon", "Kcon", "Lcon", "Mcon", "Icon", "Wcon",
		"Ylag", "Alag", "Klag", "Llag", "Mlag", "Ilag", "Wlag")
	t1data <- t0data(idvar, cbind(Y, A, K, L, M, I, omegainit))
	names(t1data) <- c("idvar", "Y1", "A1", "K1", "L1", "M1", "I1", "W1")
	############################################################################
	############ Initialization for EM Algorithm ###############################
	############################################################################
	#Initialize Paramater Values for Output
	resYinit <- matrix(rq((Y-omegainit)~PF(A=A, K=K, L=L, M=M, omega=omegainit, method=method)-1, tau=vectau)$coef, ncol=ntau)
	#Initialize Paramater Values for Labor
	resLinit <- matrix(rq(L~LX(A=A, K=K, omega=omegainit)-1, tau=vectau)$coef, ncol=ntau)
	#Initialize Parameter Values for Labor
	resMinit <- matrix(rq(M~MX(A=A, K=K, omega=omegainit)-1, tau=vectau)$coef, ncol=ntau)
	#Initialize Parameter Values for Omega_{t}
	resWTinit <- matrix(rq(lagdata$Wcon~WX(A=lagdata$Acon, omega=lagdata$Wlag)-1, tau=vectau)$coef, ncol=ntau)
	#Initialize Parameter Values for Capital
	resKTinit <- matrix(rq(lagdata$Kcon~KX(A=lagdata$Acon, K=lagdata$Klag, omega=lagdata$Wlag)-1, tau=vectau)$coef, ncol=ntau)
	#Initial Parameter Values for Initial Productivity
	#I let initial productivity be normally-distributed
	resW1coef <- lm(rnorm(N)~W1X(A=t1data$A1)-1)
	resW1init <- c(as.numeric(coef(resW1coef)), mean(resid(resW1coef)^2))
	#Initial Parameter Values for Initial Capital
	resK1coef <- lm(t1data$K1~K1X(A=t1data$A1, omega=t1data$W1)-1)
	resK1init <- c(as.numeric(coef(resK1coef)), mean(resid(resK1coef)^2))
	#Initial Parameter Values for Laplace Parameters
	#Initial Parameter Values for Laplace Parameters
	yb1init <- 20; ybLinit <- 10;
	lb1init <- 3; lbLinit <- 4;
	mb1init <- mbLinit <- 2
	wtb1init <- wtbLinit <- 5
	ktb1init <- ktbLinit <- 5
	############################################################################
	#E-Step
	############################################################################
	#Variance Random Walk proposal
	varRW <- .05
	#Matrices for Storing Results
	# dims <- list(Y=nrow(resYinit), L=nrow(resLinit), M=nrow(resMinit), W=nrow(resWTinit), I=length(resIinit))
	dims <- list(Y=nrow(resYinit), L=nrow(resLinit), M=nrow(resMinit), W=nrow(resWTinit), K=nrow(resKTinit), W1=length(resW1init), K1=length(resK1init))
	resY <- array(0, c(maxiter, dims$Y, ntau))
	resL <- array(0, c(maxiter, dims$L, ntau))
	resM <- array(0, c(maxiter, dims$M, ntau))
	resWT <- array(0, c(maxiter, dims$W, ntau))
	resKT <- array(0, c(maxiter, dims$K, ntau))
	resW1 <- array(0, c(maxiter, dims$W1))
	resK1 <- array(0, c(maxiter, dims$K1))
	yb1 <- array(0, c(maxiter, 1)); ybL <- array(0, c(maxiter, 1))
	lb1 <- array(0, c(maxiter, 1)); lbL <- array(0, c(maxiter, 1))
	mb1 <- array(0, c(maxiter, 1)); mbL <- array(0, c(maxiter, 1))
	wtb1 <- array(0, c(maxiter, 1)); wtbL <- array(0, c(maxiter, 1))
	ktb1 <- array(0, c(maxiter, 1)); ktbL <- array(0, c(maxiter, 1))
	#Diagnostic Matrices
	#Matrix for Average Acceptance Rates
	acc <- array(0, c(draws, maxiter))
	#Matrix for posterior distribution
	post <- array(0, c(N, draws, maxiter))
	#Begin EM Algorithm
	for (iter in 1:maxiter){
		print(sprintf("EM Iteration %i", iter))
		mseed <- seed+iter
		set.seed(mseed)
		resinit <- list(resYinit=resYinit, yb1init=yb1init, ybLinit=ybLinit, resLinit=resLinit, lb1init=lb1init, 
			lbLinit=lbLinit, resMinit=resMinit, mb1init=mb1init, mbLinit=mbLinit, 
			resWTinit=resWTinit, wtb1init=wtb1init, wtbLinit=wtbLinit, resW1init=resW1init, resKTinit=resKTinit, ktb1init=ktb1init, ktbLinit=ktbLinit, resK1init=resK1init)
		print(resinit$resYinit)
		print(resinit$resWTinit)
		#Initial Guess for Posterior density using initial unobservables
		omega <- omegainit
		#Posterior density
		oldpost <- posterior(data=data.frame(idvar=idvar, Y=Y, A=A, K=K, L=L, M=M, I=I, omega=omega), vectau=vectau, par=resinit, method=method)
		for (j in 1:draws){
			dseed <- mseed+j
			set.seed(dseed)
			newomega <- omega+rnorm(length(idvar), mean=0, sd=sqrt(varRW))
			#Calculate Posterior for Proposed Chain
			newpost <- posterior(data=data.frame(idvar=idvar, Y=Y, A=A, K=K, L=L, M=M, I=I, omega=newomega), vectau=vectau, par=resinit, method=method)
			#Log acceptance probability
			loga <- apply(newpost-oldpost, 1, function(x) min(x, 0))
			#Create a set of N indices for acceptance rule
			Nindex <- log(runif(N))<loga
			#Create sequence of vectors for accepting each time period of each accepted firm draw
			NTindex <- rep(Nindex, as.data.frame(table(idvar))$Freq)
			#Update unobservables that satisfy acceptance rule according to the RW process
			oldpost[Nindex] <- newpost[Nindex]
			omega[NTindex] <- newomega[NTindex]
			#Average Acceptance Rate
			acc[j,iter] <- mean(Nindex)
			post[,,iter][,j] <- oldpost
		}
		accrate <- colMeans(acc)
		print(sprintf("Acceptance Rate %.3f", accrate[iter]))
		post[,,iter] <- posterior(data=data.frame(idvar=idvar, Y=Y, A=A, K=K, L=L, M=M, I=I, omega=omega), vectau=vectau, par=resinit, method=method)
		print(sprintf("Average Posterior %.3f", mean(post[,,iter])))
		#This is the draws+1 draw of the unobservable in the stEM algorithm
		mat <- omega
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
		lagdata <- lagdata(idvar, cbind(Y, A, K, L, M, I, mat))
		names(lagdata) <- c("idvar", "Ycon", "Acon", "Kcon", "Lcon", "Mcon", "Icon", "Wcon",
		"Ylag", "Alag", "Klag", "Llag", "Mlag", "Ilag", "Wlag")
		t1data <- t0data(idvar, cbind(Y, A, K, L, M, I, mat))
		names(t1data) <- c("idvar", "Y1", "A1", "K1", "L1", "M1", "I1", "W1")
		for (q in 1:ntau){
			resY[,,q][iter,] <- matrix(rq((Y-mat)~PF(A=A, K=K, L=L, M=M, omega=mat, method=method)-1, tau=vectau[q], method="fn")$coef, nrow=dims$Y)
			resL[,,q][iter,] <- matrix(rq(L~LX(A=A, K=K, omega=mat)-1, tau=vectau[q], method="fn")$coef, nrow=dims$L)
			resM[,,q][iter,] <- matrix(rq(M~MX(A=A, K=K, omega=mat)-1, tau=vectau[q], method="fn")$coef, nrow=dims$M)
			resWT[,,q][iter,] <- matrix(rq(lagdata$Wcon~WX(A=lagdata$Acon, omega=lagdata$Wlag)-1, tau=vectau[q], method="fn")$coef, nrow=dims$W)
			resKT[,,q][iter,] <- matrix(rq(lagdata$Kcon~KX(A=lagdata$Acon, K=lagdata$Klag, omega=lagdata$Wlag)-1, tau=vectau[q], method="fn")$coef, nrow=dims$K)
		}
		##########################################################################################
		#Initial Productivity Estimation
		resW1LM <- lm(t1data$W1~W1X(A=t1data$A1)-1)
		resW1coef <- as.numeric(coef(resW1LM))
		resW1var <- mean(resid(resW1LM)^2)
		resW1[iter,] <- c(resW1coef, resW1var)
		#Initial Capital Estimation
		resK1LM <- lm(t1data$K1~K1X(A=t1data$A1, omega=t1data$W1)-1)
		resK1coef <- as.numeric(coef(resK1LM))
		resK1var <- mean(resid(resK1LM)^2)
		resK1[iter,] <- c(resK1coef, resK1var)
		###########################################################################################
		#Tail Parameters
		#Output
		yb <- expb(YX=(Y-mat), XX=PF(A=A, K=K, L=L, M=M, omega=mat, method=method), par1=resY[,,1][iter,], parL=resY[,,ntau][iter,])
		yb1[iter,] <- yb$b1
		ybL[iter,] <- yb$bL
		#Labor
		lb <- expb(YX=L, XX=LX(A=A, K=K, omega=mat), par1=resL[,,1][iter,], parL=resL[,,ntau][iter,])
		lb1[iter,] <- lb$b1
		lbL[iter,] <- lb$bL
		#Materials
		mb <- expb(YX=M, XX=MX(A=A, K=K, omega=mat), par1=resM[,,1][iter,], parL=resM[,,ntau][iter,])
		mb1[iter,] <- mb$b1
		mbL[iter,] <- mb$bL
		#Productivity
		wtb <- expb(YX=lagdata$Wcon, XX=WX(A=lagdata$Acon, omega=lagdata$Wlag), par1=resWT[,,1][iter,], parL=resWT[,,ntau][iter,])
		wtb1[iter,] <- wtb$b1
		wtbL[iter,] <- wtb$bL
		#Capital
		ktb <- expb(YX=lagdata$Kcon, XX=KX(A=lagdata$Acon, K=lagdata$Klag, omega=lagdata$Wlag), par1=resKT[,,1][iter,], parL=resKT[,,ntau][iter,])
		ktb1[iter,] <- ktb$b1
		ktbL[iter,] <- ktb$bL
		#############################################################################################
		# Use Estimates as parameters in the MH Algorithm
		resYinit <- as.matrix(resY[iter,,])
		resLinit <- as.matrix(resL[iter,,])
		resMinit <- as.matrix(resM[iter,,])
		resWTinit <- as.matrix(resWT[iter,,])
		resKTinit <- as.matrix(resKT[iter,,])
		resW1init <- as.matrix(resW1[iter,])
		resK1init <- as.matrix(resK1[iter,])
		yb1init <- yb1[iter,]; ybLinit <- ybL[iter,]
		lb1init <- lb1[iter,]; lbLinit <- lbL[iter,]
		mb1init <- mb1[iter,]; mbLinit <- mbL[iter,]
		wtb1init <- wtb1[iter,]; wtbLinit <- wtbL[iter,]
		ktb1init <- ktb1[iter,]; ktbLinit <- ktbL[iter,]
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
	resKTmat <- sapply(1:ntau, function(x) colMeans(resKT[,,x][(maxiter/2):maxiter,]))
	resktb1bLmat <- c(mean(ktb1[(maxiter/2):maxiter]), mean(ktbL[(maxiter/2):maxiter]))
	resK1mat <- colMeans(resK1[(maxiter/2):maxiter,])
	resW1mat <- colMeans(resW1[(maxiter/2):maxiter,])
	return(list(dims=dims, draws=draws, maxiter=maxiter, vectau=vectau, resYmat=resYmat, 
		resLmat=resLmat, resMmat=resMmat, resWTmat=resWTmat, resW1mat=resW1mat, 
		resKTmat=resKTmat, resK1mat=resK1mat, resyb1bLmat=resyb1bLmat, reslb1bLmat=reslb1bLmat, 
		resmb1bLmat=resmb1bLmat, reswtb1bLmat=reswtb1bLmat, accrate=accrate, 
		posterior=post, resY=resY, method=method))
}



