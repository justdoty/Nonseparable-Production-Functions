require(quantreg)
require(dplyr)
require(pracma)
require(MASS)
require(data.table)
require(truncnorm)
source('NLPFQR/FUN/Tensors.R')
source('NLPFQR/FUN/Posterior.R')
source('NLPFQR/FUN/Mstep.R')
source('NLPFQR/FUN/Auxfuns.R')
source('NLPFQR/FUN/omega.R')
Production_EM <- function(ntau, idvar, timevar, Y, A, K, L, M, I, RD, method="trans", maxiter, draws){
	#Start timer
	ctime <- proc.time()
	#Random Seed
	seed <- 123456
	set.seed(seed)
	#Number of Firms
	N <- length(unique(idvar))
	#Length of Panel
	T <- length(unique(timevar))
	#Number of observations per firm
	tsize <- as.data.frame(table(idvar))$Freq
	#Number of observations per year
	nsize <- as.data.frame(table(timevar))$Freq
	#Logicals for taking contemporary and lagged data
	idcon <- duplicated(idvar)
	idlag <- duplicated(idvar, fromLast=TRUE)
	#Create an indicator for staying in operation=1, and exit=0
	opind <- !(idlag+(timevar==max(timevar)))[idcon]
	#Grid of Taus
	vectau <- seq(1/(ntau+1), ntau/(ntau+1), by=1/(ntau+1))
	###################################################################################################
	#Initializations
	###################################################################################################
	#Initial LP Estimates:
	LP <- omega_est(idvar=idvar, timevar=timevar, Y=Y, A=A, K=K, L=L, M=M)
	#Elasticities
	LPhat <- LP$LPhat
	print(LPhat)
	#TFP
	LPTFP <- LP$omega
	print(summary(LPTFP))
	#First Stage Residuals
	LPresid <- LP$LPresid
	print(summary(LPresid))
	#Just the productivity component
	LPomega <- LPTFP-LPresid
	print(summary(LPomega))
	#Estimate of Average Persistence
	print(ar.ols(LPTFP, order.max=1))
	#Different initialization attempts
	omegainit <- jitter(LPomega, amount=0.5)
	# omegainit <- LPTFP
	#Estimate of Average Persistence
	print(ar.ols(LPomega, order.max=1))
	#Fill in data
	W <- omegainit-mean(omegainit)
	# Wcon <- omegainit[idcon]
	Wcon <- omegainit[idcon]-mean(omegainit[idcon])
	# Wlag <- omegainit[idlag]
	Wlag <- omegainit[idlag]-mean(omegainit[idlag])
	W1 <- omegainit[!idcon]-mean(omegainit[!idcon])
	# W1 <- omegainit[!idcon]
	#Store in data.frames
	KY <- K-mean(K); LY <- L-mean(L); MY <- M-mean(M)
	KW <- K[idcon]-mean(K[idcon])
	K1 <- K[!idcon]-mean(K[!idcon]); L1 <- L[!idcon]-mean(L[!idcon]); M1 <- M[!idcon]-mean(!idcon);
	ydata <- data.frame(idvar=idvar, year=timevar, Y=Y-mean(Y), KY=KY, LY=LY, MY=MY, WY=W)
	inpdata <- data.frame(idvar=idvar, year=timevar, KX=KY, LX=L-mean(L), LXY=LY, MX=M-mean(M), IX=I-mean(I), WX=W)
	wtdata <- data.frame(idvar=idvar[idcon], year=timevar[idcon], Wcon=Wcon, Wlag=Wlag, Kcon=KW)
	w1data <- data.frame(idvar=idvar[!idcon], year=timevar[!idcon], W1=W1, K1=K1, L1=L1, M1=M1)
	#Initialize Paramater Values for Output
	resYinit <- rq(Y~PF(K=KY, L=LY, M=MY, omega=WY, method=method)-1, tau=vectau, data=ydata)
	#Initialize Paramater Values for Labor
	resLinit <- rq(LX~LX(K=KX, omega=WX)-1, tau=vectau, data=inpdata)
	#Initialize Paramater Values for Materials
	resMinit <- rq(MX~MX(K=KX, omega=WX)-1, tau=vectau, data=inpdata)
	#With Labor
	# resMinit <- rq(MX~MX(K=KX, L=LXY, omega=WX)-1, tau=vectau, data=inpdata)
	#Initialize Parameter Values for Investment (Optional for testing)
	resIinit <- rq(IX~IX(K=KX, omega=WX)-1, tau=vectau, data=inpdata)
	#Initialize Parameter Values for Omega_{t}
	resWTinit <- rq(Wcon~WX(omega=Wlag)-1, tau=vectau, data=wtdata)
	#Initialize Parameter Values for Omega_{1}
	# resW1init <- matrix(c(mean(W1), var(W1)), ncol=2)
	#Attempt with initial input dependence
	resW1init <- rq(W1~WX1(K=K1, L=L1 , M=M1)-1, tau=vectau, data=w1data)
	#Survival Probability
	# omegabar <- glm(opind~WBAR(omega=Wlag, K=Kcon)-1, family=binomial(link="probit"), data=wtdata)
	# barcoef <- as.numeric(omegabar$coefficients)
	# Phat <- predict(omegabar, type="response")
	# XX <- WX(omega=Wlag)
	# for (j in 1:ntau){
	# 	Ghat <- (vectau[j]-Phat)/(1-Phat)
	# 	RHS <- t(XX)%*%(1-Ghat)
	# 	restry <- rq.fit.fnb(XX, wtdata$Wcon, tau=vectau[j], rhs=RHS)$coefficients
	# 	print(restry)
	# }
	# print(resWTinit)
	#Initial Parameter Values for Laplace Parameters
	#Output
	# ybinit <- expb(YX=ydata$Y, pred=predict(resYinit, ydata))
	ybinit <- list(b1=10, bL=15)
	#Labor
	# lbinit <- expb(YX=inpdata$L, pred=predict(resLinit, inpdata))
	lbinit <- list(b1=10, bL=15)
	#Materials
	# mbinit <- expb(YX=inpdata$M, pred=predict(resMinit, inpdata))
	mbinit <- list(b1=10, bL=15)
	#Investment
	# ibinit <- expb(YX=inpdata$I, pred=predict(resIinit, inpdata))
	ibinit <- list(b1=10, bL=15)
	#Productivity t>1
	# wtbinit <- expb(YX=wtdata$Wcon, pred=predict(resWTinit, wtdata))
	wtbinit <- list(b1=10, bL=15)
	# w1binit <- expb(YX=w1data$W1, pred=predict(resW1init, w1data))
	w1binit <- list(b1=10, bL=15)
	#Number of Parameters
	# dims <- list(Y=nrow(resYinit$coef), L=nrow(resLinit$coef), M=nrow(resMinit$coef), W=nrow(resWTinit$coef), I=nrow(resIinit$coef), W1=ncol(resW1init))
	dims <- list(Y=nrow(resYinit$coef), L=nrow(resLinit$coef), M=nrow(resMinit$coef), W=nrow(resWTinit$coef), I=nrow(resIinit$coef), W1=nrow(resW1init$coef))
	#############################################################################################
	#Normalization of Production Function 
	#############################################################################################
	if (method=="cobbN"){
		#Centers the mean of the constant at 0
		resYinit$coefficients[1,] <- resYinit$coefficients[1,]-mean(resYinit$coefficients[1,])-((1-vectau[ntau])/ybinit$b1-vectau[1]/ybinit$bL)
	} 
	else if (method=="transN"){
		#Centers the mean of the constant at 0
		resYinit$coefficients[1,] <- resYinit$coefficients[1,]-mean(resYinit$coefficients[1,])-((1-vectau[ntau])/ybinit$b1-vectau[1]/ybinit$bL)
	} 
	else if (method=="cobb"){
		#Centers the mean of the constant at 0
		resYinit$coefficients[1,] <- resYinit$coefficients[1,]-mean(resYinit$coefficients[1,])-((1-vectau[ntau])/ybinit$b1-vectau[1]/ybinit$bL)
	} 
	else if (method=="trans"){
		# #Centers the mean of the constant at 0
		resYinit$coefficients[1,] <- resYinit$coefficients[1,]-mean(resYinit$coefficients[1,])-((1-vectau[ntau])/ybinit$b1-vectau[1]/ybinit$bL)
		resYinit$coefficients[2,] <- resYinit$coefficients[2,]-mean(resYinit$coefficients[2,])+1
	} else if (method=="hermite"){
		#Centers the mean of the constant at 0
		resYinit$coefficients[1,] <- resYinit$coefficients[1,]-mean(resYinit$coefficients[1,])-((1-vectau[ntau])/ybinit$b1-vectau[1]/ybinit$bL)
		#Centers the coefficient on omega at 1
		# resYinit$coefficients[2,] <- resYinit$coefficients[2,]-mean(resYinit$coefficients[2,])+1
	} 
	############################################################################
	#E-Step
	############################################################################
	#Variance for Random Walk
	#Firm-Specific
	# varRW <- rep((2.38^2/tsize), tsize)
	# varRW <- rep((2.38/tsize), tsize)
	# varRW <- rep((2/tsize), tsize)
	# varRW <- 0.1
	# varRW <- 0.05
	# varRW <- 0.01
	#For general proposals
	# eps <- 0.5
	#Good for the truncated normal initial productivity and normalized consant
	eps <- 0.15
	# eps <- rep((2/tsize), tsize)
	# eps <- rep((2.38^2/tsize), tsize)
	#Matrices for Storing Results
	resY <- array(0, c(maxiter, dims$Y, ntau))
	resL <- array(0, c(maxiter, dims$L, ntau))
	resM <- array(0, c(maxiter, dims$M, ntau))
	resI <- array(0, c(maxiter, dims$I, ntau))
	resWT <- array(0, c(maxiter, dims$W, ntau))
	# resW1 <- array(0, c(maxiter, dims$W1))
	resW1 <- array(0, c(maxiter, dims$W1, ntau))
	yb1 <- array(0, c(maxiter, 1)); ybL <- array(0, c(maxiter, 1))
	lb1 <- array(0, c(maxiter, 1)); lbL <- array(0, c(maxiter, 1))
	mb1 <- array(0, c(maxiter, 1)); mbL <- array(0, c(maxiter, 1))
	ib1 <- array(0, c(maxiter, 1)); ibL <- array(0, c(maxiter, 1))
	wtb1 <- array(0, c(maxiter, 1)); wtbL <- array(0, c(maxiter, 1))
	w1b1 <- array(0, c(maxiter, 1)); w1bL <- array(0, c(maxiter, 1))
	#Save maximum and minimum values for productivity (useful for bounds)
	maxwt <- array(0, c(maxiter, 1))
	minwt <- array(0, c(maxiter, 1))
	#Acceptance Rates for MCMC
	acc <- array(0, c(draws, N))
	#Begin EM Algorithm
	# omegainit <- rnorm(length(idvar))
	# omega <- omegainit
	for (iter in 1:maxiter){
		omega <- omegainit
		resinit <- list(resYinit=resYinit, ybinit=ybinit, resLinit=resLinit, lbinit=lbinit, 
			resMinit=resMinit, mbinit=mbinit, resIinit=resIinit, ibinit=ibinit, 
			resWTinit=resWTinit, wtbinit=wtbinit, resW1init=resW1init, w1binit=w1binit)
		#Posterior density
		#Store Values into data frames
		ydata$WY <- omega-mean(omega)
		inpdata$WX <- omega-mean(omega)
		wtdata$Wcon <- omega[idcon]-mean(omega[idcon])
		wtdata$Wlag <- omega[idlag]-mean(omega[idlag])
		w1data$W1 <- omega[!idcon]-mean(omega[!idcon])
		oldpost <- posterior(ydata=ydata, inpdata=inpdata, wtdata=wtdata, w1data=w1data, vectau=vectau, par=resinit, method=method)
		#Time the E-step
		etime <- proc.time()
		for (j in 1:draws){
			#RW Process with Gaussian Proposal (no covariance structure)
			# newomega <- rnorm(length(idvar), mean=omega, sd=sqrt(varRW))
			#RW Process with Uniform Proposal
			newomega <- omega+runif(length(idvar), min=-eps, max=eps)
			#Store Values into data frames
			ydata$WY <- newomega-mean(newomega)
			inpdata$WX <- newomega-mean(newomega)
			wtdata$Wcon <- newomega[idcon]-mean(newomega[idcon])
			wtdata$Wlag <- newomega[idlag]-mean(newomega[idlag])
			w1data$W1 <- newomega[!idcon]-mean(newomega[!idcon])
			#Calculate Posterior for Proposed Chain
			newpost <- posterior(ydata=ydata, inpdata=inpdata, wtdata=wtdata, w1data=w1data, vectau=vectau, par=resinit, method=method)
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
		#This is the last draw of the unobservable in the stEM algorithm
		mat <- omega
		print(summary(mat))
		######################################################################################
		#M STEP
		#######################################################################################
		#Reformatting data for optimization
		ydata$WY <- mat-mean(mat)
		inpdata$WX <- mat-mean(mat)
		wtdata$Wcon <- mat[idcon]-mean(mat[idcon])
		wtdata$Wlag <- mat[idlag]-mean(mat[idlag])
		w1data$W1 <- mat[!idcon]-mean(mat[!idcon])
		#Time the M-Step
		mtime <- proc.time()
		resYinit <- rq(Y~PF(K=KY, L=LY, M=MY, omega=WY, method=method)-1, tau=vectau, data=ydata)
		resLinit <- rq(LX~LX(K=KX, omega=WX)-1, tau=vectau, data=inpdata)
		resMinit <- rq(MX~MX(K=KX, omega=WX)-1, tau=vectau, data=inpdata)
		# resMinit <- rq(MX~MX(K=KX, L=LXY, omega=WX)-1, tau=vectau, data=inpdata)
		resIinit <- rq(IX~IX(K=KX, omega=WX)-1, tau=vectau, data=inpdata)
		resWTinit <- rq(Wcon~WX(omega=Wlag)-1, tau=vectau, data=wtdata)
		#Productivity t=1
		# resW1init <- matrix(c(mean(w1data$W1), var(w1data$W1)), ncol=2)
		# resW1[iter,] <- resW1init
		resW1init <- rq(W1~WX1(K=K1, L=L1 , M=M1)-1, tau=vectau, data=w1data)
		#Values for Laplace Parameters
		#Output
		ybinit <- expb(YX=ydata$Y, pred=predict(resYinit, ydata))
		yb1[iter,] <- ybinit$b1; ybL[iter,] <- ybinit$bL
		#Labor
		lbinit <- expb(YX=inpdata$LX, pred=predict(resLinit, inpdata))
		lb1[iter,] <- lbinit$b1; lbL[iter,] <- lbinit$bL
		#Materials
		mbinit <- expb(YX=inpdata$MX, pred=predict(resMinit, inpdata))
		mb1[iter,] <- mbinit$b1; mbL[iter,] <- mbinit$bL
		#Investment
		ibinit <- expb(YX=inpdata$IX, pred=predict(resIinit, inpdata))
		ib1[iter,] <- ibinit$b1; ibL[iter,] <- ibinit$bL
		#Productivity t>1
		wtbinit <- expb(YX=wtdata$Wcon, pred=predict(resWTinit, wtdata))
		wtb1[iter,] <- wtbinit$b1; wtbL[iter,] <- wtbinit$bL
		#Productivity t=1
		w1binit <- expb(YX=w1data$W1, pred=predict(resW1init, w1data))
		w1b1[iter,] <- w1binit$b1; w1bL[iter,] <- w1binit$bL
		#Bounds on the simulated productivity
		maxwt[iter,] <- max(mat)
		minwt[iter,] <- min(mat)
		#############################################################################################
		#Normalization of Production Function 
		#############################################################################################
		if (method=="cobbN"){
			#Centers the mean of the constant at 0
			resYinit$coefficients[1,] <- resYinit$coefficients[1,]-mean(resYinit$coefficients[1,])-((1-vectau[ntau])/ybinit$b1-vectau[1]/ybinit$bL)
		} 
		else if (method=="transN"){
			#Centers the mean of the constant at 0
			resYinit$coefficients[1,] <- resYinit$coefficients[1,]-mean(resYinit$coefficients[1,])-((1-vectau[ntau])/ybinit$b1-vectau[1]/ybinit$bL)
		} 
		else if (method=="cobb"){
			#Centers the mean of the constant at 0
			resYinit$coefficients[1,] <- resYinit$coefficients[1,]-mean(resYinit$coefficients[1,])-((1-vectau[ntau])/ybinit$b1-vectau[1]/ybinit$bL)
			resYinit$coefficients[2,] <- resYinit$coefficients[2,]-mean(resYinit$coefficients[2,])+1
		} 
		else if (method=="trans"){
			#Centers the mean of the constant at 0
			resYinit$coefficients[1,] <- resYinit$coefficients[1,]-mean(resYinit$coefficients[1,])-((1-vectau[ntau])/ybinit$b1-vectau[1]/ybinit$bL)
			resYinit$coefficients[2,] <- resYinit$coefficients[2,]-mean(resYinit$coefficients[2,])+1
		} else if (method=="hermite"){
			#Centers the mean of the constant at 0
			resYinit$coefficients[1,] <- resYinit$coefficients[1,]-mean(resYinit$coefficients[1,])-((1-vectau[ntau])/ybinit$b1-vectau[1]/ybinit$bL)
			#Centers the coefficient on omega at 1
			resYinit$coefficients[2,] <- resYinit$coefficients[2,]-mean(resYinit$coefficients[2,])+1
		} 
		#Store results in matrices
		for (q in 1:ntau){
			resY[,,q][iter,] <- resYinit$coef[,q]
			resL[,,q][iter,] <- resLinit$coef[,q]
			resM[,,q][iter,] <- resMinit$coef[,q]
			resI[,,q][iter,] <- resIinit$coef[,q]
			resWT[,,q][iter,] <- resWTinit$coef[,q]
			resW1[,,q][iter,] <- resW1init$coef[,q]
		}
		#Summarize A Few Marginal Effects
		kq <- (ydata$KY)
		lq <- (ydata$LY)
		mq <- (ydata$MY)
		wq <- (ydata$WY)
		#Output Elasticities
		print(matrix(as.numeric(resYinit$coef), ncol=ntau))
		# Capital
		kdat <- colMeans(cbind(1, lq, mq, 2*kq, wq, wq*lq, wq*mq, 2*kq*wq)/1)
		kpost <- c(2,5,7,8,11,14,16,17)+1
		ky <- kdat%*%matrix(as.numeric(resYinit$coef), ncol=ntau)[kpost,]
		#Labor
		ldat <- colMeans(cbind(1, kq, mq, 2*lq, wq, wq*kq, wq*mq, 2*lq*wq)/1)
		lpost <- c(3,5,6,9,12,14,15,18)+1
		ly <- ldat%*%matrix(as.numeric(resYinit$coef), ncol=ntau)[lpost,]
		#Materials
		mdat <- colMeans(cbind(1, lq, kq, 2*mq, wq, wq*lq, wq*kq, 2*mq*wq)/1)
		mpost <- c(4,6,7,10,13,15,16,19)+1
		my <- mdat%*%matrix(as.numeric(resYinit$coef), ncol=ntau)[mpost,]
		print("Output Elasticities")
		print(ky)
		print(ly)
		print(my)
		#Computation of Average Persistence
		aparWT <- matrix(as.numeric(resWTinit$coef), ncol=ntau)
		apar <- colMeans(sweep(WX(omega=wtdata$Wlag)[,-dims$W], 2, c(1:3), "*"))%*%aparWT[-1,]
		# apar <- matrix(colMeans(WXD(omega=wtdata$Wlag, par=matrix(as.numeric(resWTinit$coef), ncol=ntau))), ncol=ntau)
		print("Persistence")
		print(aparWT)
		print(apar)
		# print("Initial Productivity")
		# print(as.numeric(resW1init))
		#Average Input Response to Productivity
		lw <- colMeans(LXD(K=inpdata$KX, omega=inpdata$WX, par=matrix(as.numeric(resLinit$coef), ncol=ntau), pos=2, sdpos=1))
		mw <- colMeans(MXD(K=inpdata$KX, L=inpdata$LXY, omega=inpdata$WX, par=matrix(as.numeric(resMinit$coef), ncol=ntau), pos=2, sdpos=1))
		#With Labor
		# mw <- colMeans(MXD(K=inpdata$KX, L=inpdata$LXY, omega=inpdata$WX, par=matrix(as.numeric(resMinit$coef), ncol=ntau), pos=3, sdpos=sd(Y)))
		iw <- colMeans(IXD(K=inpdata$KX, omega=inpdata$WX, par=matrix(as.numeric(resIinit$coef), ncol=ntau), pos=2, sdpos=1))
		print("Inputs-Productivity")
		print(matrix(lw, ncol=ntau))
		print(matrix(mw, ncol=ntau))
		print(matrix(iw, ncol=ntau))
		#Update resinit list
		resinit$resYinit <- resYinit; resinit$ybinit <- ybinit
		resinit$resLinit <- resLinit; resinit$lbinit <- lbinit
		resinit$resMinit <- resMinit; resinit$mbinit <- mbinit
		resinit$resIinit <- resIinit; resinit$ibinit <- ibinit
		resinit$resWTinit <- resWTinit; resinit$wtbinit <- wtbinit
		resinit$resW1init <- resW1init; resinit$w1binit <- w1binit
		# print(resinit)
		print(proc.time()-mtime)
		################################################################################################
		}
	resYmat <- sapply(1:ntau, function(x) colMeans(resY[,,x][(maxiter/2):maxiter,]))
	resyb1bLmat <- c(mean(yb1[(maxiter/2):maxiter]), mean(ybL[(maxiter/2):maxiter]))
	resLmat <- sapply(1:ntau, function(x) colMeans(resL[,,x][(maxiter/2):maxiter,]))
	reslb1bLmat <- c(mean(lb1[(maxiter/2):maxiter]), mean(lbL[(maxiter/2):maxiter]))
	resMmat <- sapply(1:ntau, function(x) colMeans(resM[,,x][(maxiter/2):maxiter,]))
	resmb1bLmat <- c(mean(mb1[(maxiter/2):maxiter]), mean(mbL[(maxiter/2):maxiter]))
	resImat <- sapply(1:ntau, function(x) colMeans(resI[,,x][(maxiter/2):maxiter,]))
	resib1bLmat <- c(mean(ib1[(maxiter/2):maxiter]), mean(ibL[(maxiter/2):maxiter]))
	resWTmat <- sapply(1:ntau, function(x) colMeans(resWT[,,x][(maxiter/2):maxiter,]))
	reswtb1bLmat <- c(mean(wtb1[(maxiter/2):maxiter]), mean(wtbL[(maxiter/2):maxiter]))
	# resW1mat <- colMeans(resW1[(maxiter/2):maxiter,])
	resW1mat <- sapply(1:ntau, function(x) colMeans(resW1[,,x][(maxiter/2):maxiter,]))
	resw1b1bLmat <- c(mean(w1b1[(maxiter/2):maxiter]), mean(w1bL[(maxiter/2):maxiter]))
	maxminwtmat <- c(mean(maxwt[(maxiter/2):maxiter]), mean(minwt[(maxiter/2):maxiter]))
	print(proc.time()-ctime)
	return(list(dims=dims, draws=draws, maxiter=maxiter, vectau=vectau, resYmat=resYmat, 
		resLmat=resLmat, resMmat=resMmat, resImat=resImat, resWTmat=resWTmat, resW1mat=resW1mat, resw1b1bLmat=resw1b1bLmat, resyb1bLmat=resyb1bLmat, reslb1bLmat=reslb1bLmat, 
		resmb1bLmat=resmb1bLmat, reswtb1bLmat=reswtb1bLmat, resib1bLmat=resib1bLmat, maxminwtmat=maxminwtmat, resY=resY, method=method))

}











