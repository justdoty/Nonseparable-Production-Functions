require(quantreg)
require(dplyr)
require(pracma)
source('NLPFQR/DATA/Parametric/Tensors.R')
source('NLPFQR/DATA/Parametric/Posterior.R')
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
	omegainit <- as.matrix(rnorm(length(idvar)))
	# omegainit <- as.matrix(omega_est(idvar=idvar, timevar=timevar, Y=Y, A=A, K=K, L=L, M=M)$omega)
	#Grid of Taus
	vectau <- seq(1/(ntau+1), ntau/(ntau+1), by=1/(ntau+1))
	#Reformatting data for initial optimization
	tdata <- data.frame(idvar=idvar, Y=Y, A=A, K=K, L=L, M=M, I=I, W=omegainit)
	lagdata <- data.frame(lagdata(idvar, cbind(Y, A, K, L, M, I, omegainit)))
	names(lagdata) <- c("idvar", "Ycon", "Acon", "Kcon", "Lcon", "Mcon", "Icon", "Wcon",
		"Ylag", "Alag", "Klag", "Llag", "Mlag", "Ilag", "Wlag")
	t1data <- data.frame(t1data(idvar, cbind(Y, A, K, L, M, I, omegainit)))
	names(t1data) <- c("idvar", "Y1", "A1", "K1", "L1", "M1", "I1", "W1")
	lm.predict <- function(lm.obj, data){
		return(predict(lm.obj, data))
	}
	############################################################################
	############################################################################
	############ Initialization for EM Algorithm ###############################
	############################################################################
	#Initialize Paramater Values for Output
	resYinit <- rq(Y~PF(A=A, K=K, L=L, M=M, omega=W, method=method)-1, tau=vectau, data=tdata)
	#Initialize Parameter Values for Omega_{t}
	resWTinit <- rq(Wcon~WX(A=Acon, omega=Wlag)-1, tau=vectau, data=lagdata)
	#Initialize Paramater Values for Labor
	resLinit <- lm(L~LX(A=A, K=K, omega=omegainit)-1, data=tdata)
	resLinit <- c(resLinit$coef, sigma(resLinit))
	#Initialize Parameter Values for Materials
	resMinit <- lm(M~MX(A=A, K=K, omega=omegainit)-1, data=tdata)
	resMinit <- c(resMinit$coef, sigma(resMinit))
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
	wtb1init <- wtbLinit <- 5
	#Number of Parameters
	dims <- list(Y=nrow(resYinit$coef), L=length(resLinit), M=length(resMinit), W=nrow(resWTinit$coef), K=length(resKTinit), W1=length(resW1init), K1=length(resK1init))
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
	# 	print(sum((PF(A=A, K=K, L=L, M=M, omega=omegainit, method=method)[,c(2,13:21)]%*%as.matrix(apply(resYinit[c(2,13:21),], 1, mean)))<=0)/length(idvar))
	# 	#Median
	# 	print(sum((PF(A=A, K=K, L=L, M=M, omega=omegainit, method=method)[,c(2,13:21)]%*%as.matrix(resYinit[c(2,13:21),6]))<=0)/length(idvar))
	# 	#1st Quantile
	# 	print(sum((PF(A=A, K=K, L=L, M=M, omega=omegainit, method=method)[,c(2,13:21)]%*%as.matrix(resYinit[c(2,13:21),1]))<=0)/length(idvar))
	# 	#Last Quantile
	# 	print(sum((PF(A=A, K=K, L=L, M=M, omega=omegainit, method=method)[,c(2,13:21)]%*%as.matrix(resYinit[c(2,13:21),ntau]))<=0)/length(idvar))
	# }
	############################################################################
	#E-Step
	############################################################################
	#Variance Random Walk proposal
	varRW <- 0.05
	#Matrices for Storing Results
	resY <- array(0, c(maxiter, dims$Y, ntau))
	resWT <- array(0, c(maxiter, dims$W, ntau))
	resL <- array(0, c(maxiter, dims$L))
	resM <- array(0, c(maxiter, dims$M))
	resKT <- array(0, c(maxiter, dims$K))
	resW1 <- array(0, c(maxiter, dims$W1))
	resK1 <- array(0, c(maxiter, dims$K1))
	yb1 <- array(0, c(maxiter, 1)); ybL <- array(0, c(maxiter, 1))
	wtb1 <- array(0, c(maxiter, 1)); wtbL <- array(0, c(maxiter, 1))
	#Save maximum and minimum values for productivity (useful for bounds)
	maxwt <- array(0, c(maxiter, 1))
	minwt <- array(0, c(maxiter, 1))
	#Begin EM Algorithm
	omega <- omegainit
	for (iter in 1:maxiter){
		resinit <- list(resYinit=resYinit, yb1init=yb1init, ybLinit=ybLinit, resLinit=resLinit, 
			resMinit=resMinit, resWTinit=resWTinit, wtb1init=wtb1init, wtbLinit=wtbLinit, 
			resKTinit=resKTinit, resK1init=resK1init, resW1init=resW1init)
		#Acceptance Rates for MCMC
		acc <- array(0, c(draws, N))
		#Posterior density
		oldpost <- posterior(tdata=tdata, lagdata=lagdata, t1data=t1data, vectau=vectau, par=resinit)
		#Time the E-step
		etime <- proc.time()
		print(sprintf("Variance of Inital Draw %.3f", var(omega)))
		for (j in 1:draws){
			newomega <- omega+rnorm(length(idvar), mean=0, sd=sqrt(varRW))
			tdata$W <- newomega
			lagdata$Wcon <- newomega[iddup]
			lagdata$Wlag <- newomega[(1:length(idvar))[iddup]-1]
			t1data$W1 <- newomega[!iddup]
			#Calculate Posterior for Proposed Chain
			newpost <- posterior(tdata=tdata, lagdata=lagdata, t1data=t1data, vectau=vectau, par=resinit)
			# print(system.time(posterior(tdata=tdata, lagdata=lagdata, t1data=t1data, vectau=vectau, par=resinit)))
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
		tdata$W <- mat
		lagdata$Wcon <- mat[iddup]
		lagdata$Wlag <- mat[(1:length(idvar))[iddup]-1]
		t1data$W <- mat[!iddup]
		#Time the M-Step
		mtime <- proc.time()
		resYinit <- rq(Y~PF(A=A, K=K, L=L, M=M, omega=W, method=method)-1, tau=vectau, method="fn", data=tdata)
		resWTinit <- rq(Wcon~WX(A=Acon, omega=Wlag)-1, tau=vectau, method="fn", data=lagdata)
		for (q in 1:length(vectau)){
			resY[,,q][iter,] <- resYinit$coef[,q]
			resWT[,,q][iter,] <- resWTinit$coef[,q]
		}
		##########################################################################################
		#Labor Esimation
		resLinit <- lm(L~LX(A=A, K=K, omega=W)-1, data=tdata)
		resLinit <- c(resLinit$coef, sigma(resLinit))
		resL[iter,] <- resLinit
		#Materials Estimation
		resMinit <- lm(M~MX(A=A, K=K, omega=W)-1, data=tdata)
		resMinit <- c(resMinit$coef, sigma(resMinit))
		resM[iter,] <- resMinit
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
		#Productivity
		wtb <- expb(YX=lagdata$Wcon, resinit=predict(resWTinit, lagdata))
		wtb1[iter,] <- wtb$b1
		wtbL[iter,] <- wtb$bL
		################################################################################################
		# Use Estimates as parameters in the MH Algorithm
		yb1init <- yb1[iter,]; ybLinit <- ybL[iter,]
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
		# 	print(sum((PF(A=A, K=K, L=L, M=M, omega=mat, method=method)[,c(2,13:21)]%*%as.matrix(apply(resYinit[c(2,13:21),], 1, mean)))<=0)/length(idvar))
		# 	#Median
		# 	print(sum((PF(A=A, K=K, L=L, M=M, omega=mat, method=method)[,c(2,13:21)]%*%as.matrix(resYinit[c(2,13:21),6]))<=0)/length(idvar))
		# 	#1st Quantile
		# 	print(sum((PF(A=A, K=K, L=L, M=M, omega=mat, method=method)[,c(2,13:21)]%*%as.matrix(resYinit[c(2,13:21),1]))<=0)/length(idvar))
		# 	#Last Quantile
		# 	print(sum((PF(A=A, K=K, L=L, M=M, omega=mat, method=method)[,c(2,13:21)]%*%as.matrix(resYinit[c(2,13:21),ntau]))<=0)/length(idvar))
		# }
		#Print Production Function Elasticities and Average Persistence
		print(matrix(as.numeric(resYinit$coef), nrow=dims$Y))
		# apar <- mean(colMeans(apply(resWTinit$coef, 2,  function(x) cbind(1, (A-mean(A)/sd(A)), D.tensor.prod(4, as.matrix(omega), norm=TRUE)[,-1])%*%as.matrix(x))))
		# print(sprintf("Average Persistence %.3f", apar))
		print(proc.time()-mtime)
		################################################################################################
		}
	resYmat <- sapply(1:ntau, function(x) colMeans(resY[,,x][(maxiter/2):maxiter,]))
	resyb1bLmat <- c(mean(yb1[(maxiter/2):maxiter]), mean(ybL[(maxiter/2):maxiter]))
	resWTmat <- sapply(1:ntau, function(x) colMeans(resWT[,,x][(maxiter/2):maxiter,]))
	reswtb1bLmat <- c(mean(wtb1[(maxiter/2):maxiter]), mean(wtbL[(maxiter/2):maxiter]))
	resLmat <- colMeans(resL[(maxiter/2):maxiter,])
	resMmat <- colMeans(resM[(maxiter/2):maxiter,])
	resKTmat <- colMeans(resKT[(maxiter/2):maxiter,])
	resK1mat <- colMeans(resK1[(maxiter/2):maxiter,])
	resW1mat <- colMeans(resW1[(maxiter/2):maxiter,])
	maxminwtmat <- c(mean(maxwt[(maxiter/2):maxiter]), mean(minwt[(maxiter/2):maxiter]))
	print(proc.time()-ctime)
	return(list(dims=dims, draws=draws, maxiter=maxiter, vectau=vectau, resYmat=resYmat, 
		resLmat=resLmat, resMmat=resMmat, resKTmat=resKTmat, resWTmat=resWTmat, resK1mat=resK1mat, 
		resW1mat=resW1mat, resyb1bLmat=resyb1bLmat, reswtb1bLmat=reswtb1bLmat, maxminwtmat=maxminwtmat, resY=resY, method=method))

}























