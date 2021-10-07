require(quantreg)
require(dplyr)
require(pracma)
require(MASS)
require(data.table)
require(truncnorm)
source('NLPFQR/DATA/Selection/Tensors.R')
source('NLPFQR/DATA/Selection/Posterior.R')
source('NLPFQR/DATA/Selection/Mstep.R')
source('NLPFQR/DATA/Selection/Auxfuns.R')
source('NLPFQR/DATA/Selection/omega.R')
Production_EM <- function(ntau, idvar, timevar, Y, A, K, L, M, I, RD, method="trans", maxiter, draws){
	#Start timer
	ctime <- proc.time()
	#Set Seed
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
	#Exit=1, Stay=0
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
	#First Stage Residuals
	LPresid <- LP$LPresid
	#Just the productivity component
	LPomega <- LPTFP-LPresid
	#Initialization
	omegainit <- jitter(LPomega, amount=0.5)
	#De-mean data
	Y <- Y-mean(Y)
	K <- K-mean(K)
	L <- L-mean(L)
	M <- M-mean(M)
	I <- I-mean(I)
	#De-mean productivity
	omegainitdat <- omegainit-mean(omegainit)
	# omegainitdat <- omegainit-mean(omegainit)
	#Store Values into data framesa
	W <- omegainitdat
	Wcon <- omegainitdat[idcon]
	Wlag <- omegainitdat[idlag]
	W1 <- omegainitdat[!idcon]
	KY <- K; LY <- L; MY <- M
	KW <- K[idcon]
	K1 <- K[!idcon]; L1 <- L[!idcon]; M1 <- M[!idcon];
	ydata <- data.frame(idvar=idvar, year=timevar, Y=Y, KY=KY, LY=LY, MY=MY, WY=W)
	inpdata <- data.frame(idvar=idvar, year=timevar, KX=KY, LX=L, LXY=LY, MX=M, IX=I, WX=W)
	wtdata <- data.frame(idvar=idvar[idcon], year=timevar[idcon], Wcon=Wcon, Wlag=Wlag, Kcon=KW)
	w1data <- data.frame(idvar=idvar[!idcon], year=timevar[!idcon], W1=W1, K1=K1, L1=L1, M1=M1)
	#Initialize Paramater Values for Output
	resYinit <- rq(Y~PF(K=KY, L=LY, M=MY, omega=WY, method=method)-1, tau=vectau, data=ydata)
	#Initialize Paramater Values for Labor
	resLinit <- rq(LX~LX(K=KX, omega=WX)-1, tau=vectau, data=inpdata)
	#Initialize Paramater Values for Materials
	resMinit <- rq(MX~MX(K=KX, L=LY, omega=WX)-1, tau=vectau, data=inpdata)
	#Initialize Parameter Values for Investment 
	resIinit <- rq(IX~IX(K=KX, omega=WX)-1, tau=vectau, data=inpdata)
	#Attempt with initial input dependence
	resW1init <- rq(W1~WX1(K=K1, L=L1 , M=M1)-1, tau=vectau, data=w1data)
	########################################################################
	#Selection Correction
	#########################################################################
	#Probit Model
	omegabar <- glm(opind~WBAR(omega=Wlag, K=Kcon)-1, family=binomial(link="probit"), data=wtdata)
	barcoef <- as.numeric(omegabar$coefficients)
	#Predicted Probabilities of Exit
	Phat <- predict(omegabar, type="response")
	XX <- WX(omega=Wlag)
	resWTinit <- matrix(0, nrow=ncol(XX), ncol=ntau)
	#Selection Correction is solved using QR via Interior Point Methods with corrected tau
	for (q in 1:ntau){
		#The new "tau"
		Ghat <- (vectau[q]-Phat)/(1-Phat)
		#For firms who receive a shock above vectau[j] stay in market
	    Gind <- Ghat>0
	    Ghat <- Ghat[Gind]
		RHS <- t(XX[Gind,])%*%(1-Ghat)
		resWTinit[,q] <- as.numeric(rq.fit.fnb(XX[Gind,], wtdata$Wcon[Gind], tau=vectau[q], rhs=RHS)$coefficients)
	}
	#Selection Correction in Tails
	#For first quantile
	e1 <- wtdata$Wcon-XX%*%resWTinit[,1]
	g1 <- lm(e1~Phat+Wlag+I(Wlag*Phat)+I(Phat^2)+I(Wlag^2), data=wtdata)
	#For last quantile
	eL <- wtdata$Wcon-XX%*%resWTinit[,ntau]
	gL <- lm(eL~Phat+Wlag+I(Wlag*Phat)+I(Phat^2)+I(Wlag^2), data=wtdata)
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
	#Initial Productivity
	# w1binit <- expb(YX=w1data$W1, pred=predict(resW1init, w1data))
	w1binit <- list(b1=10, bL=15)
	#Number of Parameters
	dims <- list(Y=nrow(resYinit$coef), L=nrow(resLinit$coef), M=nrow(resMinit$coef), W=ncol(XX), I=nrow(resIinit$coef), W1=nrow(resW1init$coef), WBAR=length(barcoef))
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
	} 
	############################################################################
	#E-Step
	############################################################################
	#Set jump-size for various Metropolis-Hastings proposal distributions
	#For Normally distributed proposal distributions
	# varRW <- rep((2.38^2/tsize), tsize)
	# varRW <- 0.05
	varRW <- 0.01
	#For Uniformly distributed proposal distributions
	eps <- 0.1
	#Matrices for Storing Results
	resY <- array(0, c(maxiter, dims$Y, ntau))
	resL <- array(0, c(maxiter, dims$L, ntau))
	resM <- array(0, c(maxiter, dims$M, ntau))
	resI <- array(0, c(maxiter, dims$I, ntau))
	resWT <- array(0, c(maxiter, dims$W, ntau))
	resBAR <- array(0, c(maxiter, dims$WBAR))
	resW1 <- array(0, c(maxiter, dims$W1, ntau))
	yb1 <- array(0, c(maxiter, 1)); ybL <- array(0, c(maxiter, 1))
	lb1 <- array(0, c(maxiter, 1)); lbL <- array(0, c(maxiter, 1))
	mb1 <- array(0, c(maxiter, 1)); mbL <- array(0, c(maxiter, 1))
	ib1 <- array(0, c(maxiter, 1)); ibL <- array(0, c(maxiter, 1))
	wtb1 <- array(0, c(maxiter, 1)); wtbL <- array(0, c(maxiter, 1))
	w1b1 <- array(0, c(maxiter, 1)); w1bL <- array(0, c(maxiter, 1))
	#Save maximum and minimum values for productivity
	maxwt <- array(0, c(maxiter, 1))
	minwt <- array(0, c(maxiter, 1))
	#Acceptance Rates for MCMC
	acc <- array(0, c(draws, N))
	#Begin EM Algorithm
	# omegainit <- rnorm(length(idvar))
	omega <- omegainit
	for (iter in 1:maxiter){
		omega <- omegainit-mean(omegainit)
		resinit <- list(resYinit=resYinit, ybinit=ybinit, resLinit=resLinit, lbinit=lbinit, 
			resMinit=resMinit, mbinit=mbinit, resIinit=resIinit, ibinit=ibinit, 
			resWTinit=resWTinit, wtbinit=wtbinit, resW1init=resW1init, w1binit=w1binit, omegabar=omegabar, g1=g1, gL=gL)
		#De-mean productivity
		omegadat <- omega
		# omegadat <- omega-mean(omega)
		#Store Values into data frames
		ydata$WY <- omegadat
		inpdata$WX <- omegadat
		wtdata$Wcon <- omegadat[idcon]
		wtdata$Wlag <- omegadat[idlag]
		w1data$W1 <- omegadat[!idcon]
		#Posterior density
		oldpost <- posterior(ydata=ydata, inpdata=inpdata, wtdata=wtdata, w1data=w1data, vectau=vectau, par=resinit, method=method)
		#Time the E-step
		etime <- proc.time()
		for (j in 1:draws){
			#RW Process with Gaussian Proposal
			newomega <- rnorm(length(idvar), mean=omega, sd=sqrt(varRW))
			#RW Process with Uniform Proposal
			# newomega <- omega+runif(length(idvar), min=-eps, max=eps)
			#De-mean productivity
			newomegadat <- newomega
			# newomegadat <- newomega-mean(newomega)
			#Store Values into data frames
			ydata$WY <- newomegadat
			inpdata$WX <- newomegadat
			wtdata$Wcon <- newomegadat[idcon]
			wtdata$Wlag <- newomegadat[idlag]
			w1data$W1 <- newomegadat[!idcon]
			#Calculate Posterior for Proposed Chain
			newpost <- posterior(ydata=ydata, inpdata=inpdata, wtdata=wtdata, w1data=w1data, vectau=vectau, par=resinit, method=method)
			#Acceptance probability
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
		#De-mean productivity
		matdat <- mat
		# matdat <- mat-mean(mat)
		#Store Values into data frames
		ydata$WY <- matdat
		inpdata$WX <- matdat
		wtdata$Wcon <- matdat[idcon]
		wtdata$Wlag <- matdat[idlag]
		w1data$W1 <- matdat[!idcon]
		#Time the M-Step
		mtime <- proc.time()
		#Output
		resYinit <- rq(Y~PF(K=KY, L=LY, M=MY, omega=WY, method=method)-1, tau=vectau, data=ydata)
		#Labor
		resLinit <- rq(LX~LX(K=KX, omega=WX)-1, tau=vectau, data=inpdata)
		#Materials
		resMinit <- rq(MX~MX(K=KX, L=LY, omega=WX)-1, tau=vectau, data=inpdata)
		#Investment
		resIinit <- rq(IX~IX(K=KX, omega=WX)-1, tau=vectau, data=inpdata)
		# For Productivity at t>2
		# Probit Model
		omegabar <- glm(opind~WBAR(omega=Wlag, K=Kcon)-1, family=binomial(link="probit"), data=wtdata)
		barcoef <- as.numeric(omegabar$coefficients)
		resBAR[iter,] <- barcoef
		Phat <- predict(omegabar, type="response")
		XX <- WX(omega=wtdata$Wlag)
		for (q in 1:ntau){
			#The new "tau"
			Ghat <- (vectau[q]-Phat)/(1-Phat)
			#For firms who receive a shock above vectau[j] stay in market
		    Gind <- Ghat>0
		    Ghat <- Ghat[Gind]
			RHS <- t(XX[Gind,])%*%(1-Ghat)
			resWTinit[,q] <- as.numeric(rq.fit.fnb(XX[Gind,], wtdata$Wcon[Gind], tau=vectau[q], rhs=RHS)$coefficients)
		}
		#Initial Productivity
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
		#Selection Correction
		#For first quantile
		e1 <- wtdata$Wcon-XX%*%resWTinit[,1]
		g1 <- lm(e1~Phat+Wlag+I(Wlag*Phat)+I(Phat^2)+I(Wlag^2), data=wtdata)
		fitw1 <- XX%*%resWTinit[,1]+fitted(g1)
		#For last quantile
		eL <- wtdata$Wcon-XX%*%resWTinit[,ntau]
		gL <- lm(eL~Phat+Wlag+I(Wlag*Phat)+I(Phat^2)+I(Wlag^2), data=wtdata)
		fitwL <- XX%*%resWTinit[,ntau]+fitted(gL)
		wtbinit <- expb(YX=wtdata$Wcon, pred=cbind(fitw1, fitwL))
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
		} 
		#Store results in matrices
		for (q in 1:ntau){
			resY[,,q][iter,] <- resYinit$coef[,q]
			resL[,,q][iter,] <- resLinit$coef[,q]
			resM[,,q][iter,] <- resMinit$coef[,q]
			resI[,,q][iter,] <- resIinit$coef[,q]
			resW1[,,q][iter,] <- resW1init$coef[,q]
			resWT[,,q][iter,] <- resWTinit[,q]
		}
		#Summarize A Few Marginal Effects
		kq <- (ydata$KY)
		lq <- (ydata$LY)
		mq <- (ydata$MY)
		wq <- (ydata$WY)
		#Capital
		kdat <- colMeans(cbind(1, lq, mq, 2*kq, wq, wq*lq, wq*mq, 2*kq*wq))
		kpost <- c(2,5,7,8,11,14,16,17)+1
		ky <- kdat%*%matrix(as.numeric(resYinit$coef), ncol=ntau)[kpost,]
		#Labor
		ldat <- colMeans(cbind(1, kq, mq, 2*lq, wq, wq*kq, wq*mq, 2*lq*wq))
		lpost <- c(3,5,6,9,12,14,15,18)+1
		ly <- ldat%*%matrix(as.numeric(resYinit$coef), ncol=ntau)[lpost,]
		#Materials
		mdat <- colMeans(cbind(1, lq, kq, 2*mq, wq, wq*lq, wq*kq, 2*mq*wq))
		mpost <- c(4,6,7,10,13,15,16,19)+1
		my <- mdat%*%matrix(as.numeric(resYinit$coef), ncol=ntau)[mpost,]
		print("Output Elasticities")
		print(ky)
		print(ly)
		print(my)
		#Computation of Average Persistence
		aparWT <- resWTinit
		apar <- colMeans(sweep(WX(omega=wtdata$Wlag)[,-dims$W], 2, c(1:3), "*"))%*%aparWT[-1,]
		print("Persistence")
		print(apar)
		#Average Input Response to Productivity
		#Labor
		lw <- colMeans(LXD(K=inpdata$KX, omega=inpdata$WX, par=matrix(as.numeric(resLinit$coef), ncol=ntau), pos=2, sdpos=1))
		#Materials
		mw <- colMeans(MXD(K=inpdata$KX, L=inpdata$LXY, omega=inpdata$WX, par=matrix(as.numeric(resMinit$coef), ncol=ntau), pos=3, sdpos=1))
		#Investment
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
		resinit$omegabar <- omegabar
		resinit$g1 <- g1; resinit$gL <- gL
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
	resBARmat <- colMeans(resBAR[(maxiter/2):maxiter,])
	resW1mat <- sapply(1:ntau, function(x) colMeans(resW1[,,x][(maxiter/2):maxiter,]))
	resw1b1bLmat <- c(mean(w1b1[(maxiter/2):maxiter]), mean(w1bL[(maxiter/2):maxiter]))
	maxminwtmat <- c(mean(maxwt[(maxiter/2):maxiter]), mean(minwt[(maxiter/2):maxiter]))
	print(proc.time()-ctime)
	return(list(dims=dims, draws=draws, maxiter=maxiter, vectau=vectau, resYmat=resYmat, 
		resLmat=resLmat, resMmat=resMmat, resImat=resImat, resWTmat=resWTmat, resW1mat=resW1mat, resw1b1bLmat=resw1b1bLmat, resyb1bLmat=resyb1bLmat, reslb1bLmat=reslb1bLmat, 
		resmb1bLmat=resmb1bLmat, reswtb1bLmat=reswtb1bLmat, resBARmat=resBARmat, resib1bLmat=resib1bLmat, maxminwtmat=maxminwtmat, resY=resY, method=method))

}











