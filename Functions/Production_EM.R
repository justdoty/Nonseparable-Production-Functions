require(quantreg)
require(dplyr)
require(pracma)
require(MASS)
require(data.table)
source('NLPFQR/DATA/Nonparametric/Tensors.R')
source('NLPFQR/DATA/Nonparametric/Posterior.R')
source('NLPFQR/FUN/Mstep.R')
source('NLPFQR/FUN/Auxfuns.R')
source('NLPFQR/FUN/omega.R')
Production_EM <- function(ntau, idvar, timevar, Y, A, K, L, M, I, method="trans", maxiter, draws){
	ctime <- proc.time()
	if (all(method!=c("cobbN", "transN", "cobb", "trans", "hermite"))){
		stop("Production Function must be specified as either cobbN (Hicks-Neutral Cobb Douglas),
			transN (Hicks-Neutral Translog), cobb (NonHicks-Neutral Cobb Douglass), trans (NonHicks-Neutral Translog),
			or hermite (tensor product hermite polynomial)")
	}
	seed <- 123456
	set.seed(seed)
	#Number of Firms
	N <- length(unique(idvar))
	#Length of Panel
	T <- unique(timevar)
	#Number of observations per firms
	tsize <- as.data.frame(table(idvar))$Freq
	iddup <- duplicated(idvar)
	#Initial Unobservables 
	omegainit <- as.matrix(omega_est(idvar=idvar, timevar=timevar, Y=Y, A=A, K=K, L=L, M=M)$omega)
	# omegainit <- as.matrix(rnorm(length(idvar)))
	#Grid of Taus
	vectau <- seq(1/(ntau+1), ntau/(ntau+1), by=1/(ntau+1))
	#Reformatting data for initial optimization
	tdata <- data.frame(idvar=idvar, Y=Y, A=A, K=K, L=L, M=M, I=I, W=omegainit)
	lagdata <- data.frame(lagdata(idvar, cbind(Y, A, K, L, M, I, omegainit)))
	names(lagdata) <- c("idvar", "Ycon", "Acon", "Kcon", "Lcon", "Mcon", "Icon", "Wcon",
		"Ylag", "Alag", "Klag", "Llag", "Mlag", "Ilag", "Wlag")
	t1data <- data.frame(t1data(idvar, cbind(Y, A, K, L, M, I, omegainit)))
	names(t1data) <- c("idvar", "Y1", "A1", "K1", "L1", "M1", "I1", "W1")
	############################################################################
	############ Initialization for EM Algorithm ###############################
	############################################################################
	#Initialize Paramater Values for Output
	resYinit <- rq(Y~PF(A=A, K=K, L=L, M=M, omega=W, method=method)-1, tau=vectau, data=tdata)
	#Initialize Paramater Values for Labor
	resLinit <- rq(L~LX(A=A, K=K, omega=W)-1, tau=vectau, data=tdata)
	#Initialize Paramater Values for Materials
	resMinit <- rq(M~MX(A=A, K=K, omega=W)-1, tau=vectau, data=tdata)
	#Initialize Parameter Values for Omega_{t}
	resWTinit <- rq(Wcon~WX(A=Acon, omega=Wlag)-1, tau=vectau, data=lagdata)
	#Initialize Parameter Values for Capital (Using Normal)
	resKTinit <- lm(Kcon~KX(A=Acon, K=Klag, I=Ilag, omega=Wlag)-1, data=lagdata)
	resKTinit <- c(resKTinit$coef, sigma(resKTinit))
	#Initial Parameter Values for Initial Productivity
	resW1init <- lm(W1~W1X(A=A1)-1, data=t1data)
	resW1init <- c(resW1init$coef, sigma(resW1init))
	#Initial Parameter Values for Initial Capital
	resK1init <- lm(K1~K1X(A=A1, omega=W1)-1, data=t1data)
	resK1init <- c(resK1init$coef, sigma(resK1init))
	#Initial Parameter Values for Laplace Parameters
	yb1init <- 20; ybLinit <- 10;
	lb1init <- 3; lbLinit <- 4;
	mb1init <- mbLinit <- 2
	wtb1init <- wtbLinit <- 5
	#Number of Parameters
	dims <- list(Y=nrow(resYinit$coef), L=nrow(resLinit$coef), M=nrow(resMinit$coef), W=nrow(resWTinit$coef), K=length(resKTinit), W1=length(resW1init), K1=length(resK1init))
	#############################################################################################
	#Normalization of Production Function (Conditional Mean)
	#############################################################################################
	# if (method=="cobbN"){
	# 	#Centers the mean of the constant at 0
	# 	resYinit[1,] <- resYinit[1,]-mean(resYinit[1,])-((1-vectau[ntau])/yb1init-vectau[1]/ybLinit)
	# } else if (method=="transN"){
	# 	#Centers the mean of the constant at 0
	# 	resYinit[1,] <- resYinit[1,]-mean(resYinit[1,])-((1-vectau[ntau])/yb1init-vectau[1]/ybLinit)
	# } else if (method=="cobb"){
	# 	#In progress
	# } else if (method=="trans"){
	# 	#We can assume strict monotonicity in productivity for various functions. The verify, below calculates
	# 	#various functionals of the conditional output distribution and prints the fraction of samples that do not
	# 	#satisfy the monotonicity property
	# 	#Mean
	# 	print(sum((cbind(1, K, L, M, K*L, L*M, K*M, K^2, L^2, M^2)%*%as.matrix(apply(resYinit[c(2,13:21),], 1, mean)))<=0)/length(idvar))
	# 	#Median
	# 	print(sum((cbind(1, K, L, M, K*L, L*M, K*M, K^2, L^2, M^2)%*%as.matrix(resYinit[c(2,13:21),6]))<=0)/length(idvar))
	# 	#1st Quantile
	# 	print(sum((cbind(1, K, L, M, K*L, L*M, K*M, K^2, L^2, M^2)%*%as.matrix(resYinit[c(2,13:21),1]))<=0)/length(idvar))
	# 	#Last Quantile
	# 	print(sum((cbind(1, K, L, M, K*L, L*M, K*M, K^2, L^2, M^2)%*%as.matrix(resYinit[c(2,13:21),ntau]))<=0)/length(idvar))
	# }
	############################################################################
	#E-Step
	############################################################################
	#Variance Random Walk proposal
	varRW <- 0.05
	#Matrices for Storing Results
	resY <- array(0, c(maxiter, dims$Y, ntau))
	resL <- array(0, c(maxiter, dims$L, ntau))
	resM <- array(0, c(maxiter, dims$M, ntau))
	resWT <- array(0, c(maxiter, dims$W, ntau))
	resKT <- array(0, c(maxiter, dims$K))
	resW1 <- array(0, c(maxiter, dims$W1))
	resK1 <- array(0, c(maxiter, dims$K1))
	yb1 <- array(0, c(maxiter, 1)); ybL <- array(0, c(maxiter, 1))
	lb1 <- array(0, c(maxiter, 1)); lbL <- array(0, c(maxiter, 1))
	mb1 <- array(0, c(maxiter, 1)); mbL <- array(0, c(maxiter, 1))
	wtb1 <- array(0, c(maxiter, 1)); wtbL <- array(0, c(maxiter, 1))
	#Save maximum and minimum values for productivity (useful for bounds)
	maxwt <- array(0, c(maxiter, 1))
	minwt <- array(0, c(maxiter, 1))
	#Begin EM Algorithm
	# omega <- omegainit
	for (iter in 1:maxiter){
		omega <- omegainit
		resinit <- list(resYinit=resYinit, yb1init=yb1init, ybLinit=ybLinit, resLinit=resLinit, 
			lb1init=lb1init, lbLinit=lbLinit, resMinit=resMinit, mb1init=mb1init, mbLinit=mbLinit, 
			resWTinit=resWTinit, wtb1init=wtb1init, wtbLinit=wtbLinit, 
			resKTinit=resKTinit, resK1init=resK1init, resW1init=resW1init)
		#Acceptance Rates for MCMC
		acc <- array(0, c(draws, N))
		#Posterior density
		oldpost <- posterior(tdata=tdata, lagdata=lagdata, t1data=t1data, vectau=vectau, par=resinit)
		#Time the E-step
		etime <- proc.time()
		print(sprintf("Variance of Inital Draw %.3f", var(omega)))
		for (j in 1:draws){
			# newomega <- omega+rnorm(length(idvar), mean=0, sd=sqrt(varRW))
			newomega <- omega+runif(length(idvar), -2, 2)
			tdata$W <- newomega
			lagdata$Wcon <- newomega[iddup]
			lagdata$Wlag <- newomega[(1:length(idvar))[iddup]-1]
			t1data$W1 <- newomega[!iddup]
			#Calculate Posterior for Proposed Chain
			newpost <- posterior(tdata=tdata, lagdata=lagdata, t1data=t1data, vectau=vectau, par=resinit)
			#acceptance probability
			loga <- pmin(newpost-oldpost, 0)
			#Create a set of N indices for acceptance rule
			Nindex <- log(runif(N))<loga
			#Create sequence of vectors for accepting each time period of each accepted firm draw
			NTindex <- rep(Nindex, tsize)
			#Update unobservables that satisfy acceptance rule according to the RW process
			oldpost[Nindex] <- newpost[Nindex]
			omega[NTindex] <- newomega[NTindex]
			#Acceptance Rate
			acc[j,] <- Nindex
		}
		print(sprintf("Variance of Final Draw %.3f", var(omega)))
		print(summary(omega))
		print(sprintf("EM Iteration %i", iter))
		print(proc.time()-etime)
		#Average Acceptance Rate for Each Firm
		accrate <- colMeans(acc)
		#Average Acceptance Rate
		print(sprintf("Average Acceptance Rate %.3f", mean(accrate)))
		#Median Acceptance Rate
		print(sprintf("Median Acceptance Rate %.3f", median(accrate)))
		#Average Posterior
		print(sprintf("Average Posterior %.3f", mean(oldpost)))
		#This is the draws+1 draw of the unobservable in the stEM algorithm
		mat <- omega
		######################################################################################
		#M STEP
		#######################################################################################
		#Reformatting data for initial optimization
		#Reformatting data for initial optimization
		tdata$W <- mat
		lagdata$Wcon <- mat[iddup]
		lagdata$Wlag <- mat[(1:length(idvar))[iddup]-1]
		t1data$W <- mat[!iddup]
		#Time the M-Step
		mtime <- proc.time()
		resYinit <- rq(Y~PF(A=A, K=K, L=L, M=M, omega=W, method=method)-1, tau=vectau, data=tdata)
		resLinit <- rq(L~LX(A=A, K=K, omega=W)-1, tau=vectau, data=tdata)
		resMinit <- rq(M~MX(A=A, K=K, omega=W)-1, tau=vectau, data=tdata)
		resWTinit <- rq(Wcon~WX(A=Acon, omega=Wlag)-1, tau=vectau, data=lagdata)
		for (q in 1:ntau){
			resY[,,q][iter,] <- resYinit$coef[,q]
			resL[,,q][iter,] <- resLinit$coef[,q]
			resM[,,q][iter,] <- resMinit$coef[,q]
			resWT[,,q][iter,] <- resWTinit$coef[,q]
		}
		##########################################################################################
		#Capital Accumulation Estimation
		resKTinit <- lm(Kcon~KX(A=Acon, K=Klag, I=Ilag, omega=Wlag)-1, data=lagdata)
		resKTinit <- c(resKTinit$coef, sigma(resKTinit))
		resKT[iter,] <- resKTinit
		#Initial Productivity Estimation
		resW1init <- lm(W1~W1X(A=A1)-1, data=t1data)
		resW1init <- c(resW1init$coef, sigma(resW1init))
		resW1[iter,] <- resW1init
		#Initial Capital Estimation
		resK1init <- lm(K1~K1X(A=A1, omega=W1)-1, data=t1data)
		resK1init <- c(resK1init$coef, sigma(resK1init))
		resK1[iter,] <- resK1init
		############################################################################################
		#Tail Parameters
		#Output
		yb <- expb(YX=tdata$Y, resinit=predict(resYinit, tdata))
		yb1[iter,] <- yb$b1
		ybL[iter,] <- yb$bL
		#Labor
		lb <- expb(YX=tdata$L, resinit=predict(resLinit, tdata))
		lb1[iter,] <- lb$b1
		lbL[iter,] <- lb$bL
		#Labor
		mb <- expb(YX=tdata$M, resinit=predict(resMinit, tdata))
		mb1[iter,] <- mb$b1
		mbL[iter,] <- mb$bL
		#Productivity
		wtb <- expb(YX=lagdata$Wcon, resinit=predict(resWTinit, lagdata))
		wtb1[iter,] <- wtb$b1
		wtbL[iter,] <- wtb$bL
		################################################################################################
		# Use Estimates as parameters in the MH Algorithm
		resKTinit <- as.matrix(resKT[iter,])
		resW1init <- as.matrix(resW1[iter,])
		resK1init <- as.matrix(resK1[iter,])
		yb1init <- yb1[iter,]; ybLinit <- ybL[iter,]
		lb1init <- lb1[iter,]; lbLinit <- lbL[iter,]
		mb1init <- mb1[iter,]; mbLinit <- mbL[iter,]
		wtb1init <- wtb1[iter,]; wtbLinit <- wtbL[iter,]
		maxwt[iter,] <- max(mat)
		minwt[iter,] <- min(mat)
		#############################################################################################
		#Normalization of Production Function (Conditional Mean)
		#############################################################################################
		# if (method=="cobbN"){
		# 	#Centers the mean of the constant at 0
		# 	resYinit[1,] <- resYinit[1,]-mean(resYinit[1,])-((1-vectau[ntau])/yb1init-vectau[1]/ybLinit)
		# } else if (method=="transN"){
		# 	#Centers the mean of the constant at 0
		# 	resYinit[1,] <- resYinit[1,]-mean(resYinit[1,])-((1-vectau[ntau])/yb1init-vectau[1]/ybLinit)
		# } else if (method=="cobb"){
		# 	#In progress
		# } else if (method=="trans"){
		# 	#We can assume strict monotonicity in productivity for various functions. The verify, below calculates
		# 	#various functionals of the conditional output distribution and prints the fraction of samples that do not
		# 	#satisfy the monotonicity property
		# 	#Mean
		# 	print(sum((cbind(1, K, L, M, K*L, L*M, K*M, K^2, L^2, M^2)%*%as.matrix(apply(resYinit[c(2,13:21),], 1, mean)))<=0)/length(idvar))
		# 	#Median
		# 	print(sum((cbind(1, K, L, M, K*L, L*M, K*M, K^2, L^2, M^2)%*%as.matrix(resYinit[c(2,13:21),6]))<=0)/length(idvar))
		# 	#1st Quantile
		# 	print(sum((cbind(1, K, L, M, K*L, L*M, K*M, K^2, L^2, M^2)%*%as.matrix(resYinit[c(2,13:21),1]))<=0)/length(idvar))
		# 	#Last Quantile
		# 	print(sum((cbind(1, K, L, M, K*L, L*M, K*M, K^2, L^2, M^2)%*%as.matrix(resYinit[c(2,13:21),ntau]))<=0)/length(idvar))
		# }
		#Print Production Function Elasticities and Average Persistence
		print(matrix(as.numeric(resYinit$coef), ncol=ntau))
		print(resKTinit)
		print(resK1init)
		print(resW1init)
		#Computation of Average Persistence
		apar <- colMeans(sweep(WX(A=lagdata$Acon, omega=lagdata$Wlag)[,-c(2,dims$W)], 2, c(1:3), "*"))%*%matrix(as.numeric(resWTinit$coef), ncol=ntau)[-c(1,2),]
		print(as.matrix(apar))
		print(proc.time()-mtime)
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
	maxminwtmat <- c(mean(maxwt[(maxiter/2):maxiter]), mean(minwt[(maxiter/2):maxiter]))
	resKTmat <- colMeans(resKT[(maxiter/2):maxiter,])
	resK1mat <- colMeans(resK1[(maxiter/2):maxiter,])
	resW1mat <- colMeans(resW1[(maxiter/2):maxiter,])
	print(proc.time()-ctime)
	return(list(dims=dims, draws=draws, maxiter=maxiter, vectau=vectau, resYmat=resYmat, 
		resLmat=resLmat, resMmat=resMmat, resKTmat=resKTmat, resWTmat=resWTmat, resK1mat=resK1mat, resW1mat=resW1mat, resyb1bLmat=resyb1bLmat, reslb1bLmat=reslb1bLmat, 
		resmb1bLmat=resmb1bLmat, reswtb1bLmat=reswtb1bLmat, maxminwtmat=maxminwtmat, resY=resY, method=method))

}


























