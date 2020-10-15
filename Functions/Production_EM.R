setwd('/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Functions')
require(quantreg)
require(dplyr)
source('Tensors.R')
source('Posterior.R')
source('Moments.R')
source('Mstep.R')
source('Auxfuns.R')
source('LP_init.R')
Production_EM <- function(ntau, idvar, timevar, Y, K, L, M, I, A, maxiter, draws, Mdraws, seed){
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
	#Age
	A <- as.matrix(A)
	#Number of Firms
	N <- length(unique(idvar))
	#Grid of Taus
	vectau <- seq(1/(ntau+1), ntau/(ntau+1), by=1/(ntau+1))
	#Degrees of the Univariate Hermite Polynomials
	#Output
	MY <- as.numeric(c(MYK=1, MYL=2, MYM=1, MYA=1, MYW=2))
	#Labor
	ML <- as.numeric(c(MLK=1, MLA=1, MLW=2))
	#Materials
	MM <- as.numeric(c(MMK=1, MMA=1, MMW=2))
	#Productivity t>1
	MW <- as.numeric(c(MWA=1, MWW=2))
	#Productivity t=1
	MW1 <- as.numeric(c(MW1A=1))
	#Investment
	MI <- as.numeric(c(MIK=2, MIW=2, MIA=1))
	#Combined
	MH <- list(MY=MY, ML=ML, MM=MM, MW=MW, MW1=MW1, MI=MI)
	TFP <- as.matrix(LP_est(idvar=idvar, timevar=timevar, Y=Y, K=K, L=L, M=M, A=A))
	print(var(TFP))
	#Reformatting data for initial optimization
	WTdata <- lagdata(idvar, cbind(A, TFP))
	names(WTdata) <- c("idvar", "Acon", "Ucon", "Alag", "Ulag")
	W0data <- t0data(idvar, cbind(A, TFP))
	names(W0data) <- c("idvar", "A0", "U0")
	############################################################################
	############ Initialization for EM Algorithm ###############################
	############################################################################
	#Initialize Paramater Values for Output
	resYinit <- matrix(rq(Y~tensor.prod(MY, cbind(K, L, M, A, TFP))-1, tau=vectau)$coef, nrow=prod(MY+1), ncol=ntau)
	#Initialize Paramater Values for Labor
	resLinit <- matrix(rq(L~tensor.prod(ML, cbind(K, A, TFP))-1, tau=vectau)$coef, nrow=prod(ML+1), ncol=ntau)
	#Initialize Parameter Values for Labor
	resMinit <- matrix(rq(M~tensor.prod(MM, cbind(K, A, TFP))-1, tau=vectau)$coef, nrow=prod(MM+1), ncol=ntau)
	#Initialize Parameter Values for Omega_{t}
	resWTinit <- matrix(rq(WTdata$Ucon~tensor.prod(MW, cbind(WTdata$Acon, WTdata$Ulag))-1, tau=vectau)$coef, nrow=prod(MW+1), ncol=ntau)
	#Initialize Parameter Values for Omega_{0}
	resW0init <- matrix(rq(W0data$U0~tensor.prod(MW1, W0data$A0)-1, tau=vectau)$coef, nrow=prod(MW1+1), ncol=ntau)
	#Initialize Paramater Values for Investment (Nonlinear Regression Model)
	resIcoef <- lm(I~tensor.prod(MI, cbind(K, A, TFP))-1)
	#Vector of Investment coefficients and estimate of  variance
	resIinit <- c(as.numeric(coef(resIcoef)), mean(resid(resIcoef)^2))
	#Initial Parameter Values for Laplace Parameters
	yb1init <- ybLinit <- 1
	lb1init <- lbLinit <- 1
	mb1init <- mbLinit <- 1
	wtb1init <- wtbLinit <- 1
	w0b1init <- w0bLinit <- 1
	############################################################################
	############################################################################
	############################## EM Algorithm ################################
	############################################################################
	############################################################################
	#Variance Random Walk proposal
	varRW <- .05
	#Matrices for Storing Results
	resY <- array(0, c(maxiter, prod(MY+1), ntau))
	resL <- array(0, c(maxiter, prod(ML+1), ntau))
	resM <- array(0, c(maxiter, prod(MM+1), ntau))
	resWT <- array(0, c(maxiter, prod(MW+1), ntau))
	resW0 <- array(0, c(maxiter, prod(MW1+1), ntau))
	resI <- array(0, c(maxiter, (prod(MI+1)+1)))
	yb1 <- array(0, c(maxiter, 1)); ybL <- array(0, c(maxiter, 1))
	lb1 <- array(0, c(maxiter, 1)); lbL <- array(0, c(maxiter, 1))
	mb1 <- array(0, c(maxiter, 1)); mbL <- array(0, c(maxiter, 1))
	wtb1 <- array(0, c(maxiter, 1)); wtbL <- array(0, c(maxiter, 1))
	w0b1 <- array(0, c(maxiter, 1)); w0bL <- array(0, c(maxiter, 1))
	#Begin EM Algorithm
	for (iter in 1:maxiter){
		print(iter)
		#This sets how many Metropolis-Hastings "burn in" samples
		r <- -draws
		#Initialize Matrix for Averaging over Accepted Draws
		mat <- matrix(0, length(Y), Mdraws)
		#Initial Parameter Values 
		resinit <- list(resYinit=resYinit, yb1init=yb1init, ybLinit=ybLinit, resLinit=resLinit, lb1init=lb1init, 
			lbLinit=lbLinit, resMinit=resMinit, mb1init=mb1init, mbLinit=mbLinit, 
			resWTinit=resWTinit, wtb1init=wtb1init, wtbLinit=wtbLinit, resW0init=resW0init, w0b1init=w0b1init, 
			w0bLinit=w0bLinit, resIinit=resIinit)
		#Initial Guess for Unobservables
		un <- as.matrix(rnorm(length(Y), 0, sqrt(varRW)))
		#Calculate posterior density at this guess
		un_dens <- posterior(data=cbind(idvar, Y, K, L, M, A, I, un), MH=MH, vectau=vectau, par=resinit)
		for (j in r:Mdraws){
			Mseed <- seed+j
			set.seed(Mseed)
			print(j)
			#Proposal Distribution
			try_un <- as.matrix(rnorm(length(Y), 0, sqrt(varRW)))
			#Random Walk update
			y <- un+try_un
			loglik <- posterior(data=cbind(idvar, Y, K, L, M, A, I, y), MH=MH, vectau=vectau, par=resinit)
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
		WTdata <- lagdata(idvar, cbind(A, mat))
		names(WTdata) <- c("idvar", "Acon", "Ucon", "Alag", "Ulag")
		W0data <- t0data(idvar, cbind(A, mat))
		names(W0data) <- c("idvar", "A0", "U0")
		for (q in 1:ntau){
			print(q)
			if (Mdraws==1){
				resY[,,q][iter,] <- rq(Y~tensor.prod(MY, cbind(K, L, M, A, mat))-1, tau=vectau[q])$coef
				resL[,,q][iter,] <- rq(L~tensor.prod(ML, cbind(K, A, mat))-1, tau=vectau[q])$coef
				resM[,,q][iter,] <- rq(M~tensor.prod(MM, cbind(K, A, mat))-1, tau=vectau[q])$coef
				resWT[,,q][iter,] <- rq(WTdata$Ucon~tensor.prod(MW, cbind(WTdata$Acon, WTdata$Ulag))-1, tau=vectau[q])$coef
				resW0[,,q][iter,] <- rq(W0data$U0~tensor.prod(MW1, W0data$A0)-1, tau=vectau[q])$coef
			} else {
				#Initialize Gradient Descent Algorithms
			    #Tolerance for gradient (close to zero)
			    tol <- 1e-2
			    #Stepsize for gradient
			    stepsize <- 100
			    #Number of gradient descent iterations
			    nsteps <- 800
			    nbatch <- 1
				resY[,,q][iter,] <- sgradfunq(idvar=idvar, X=cbind(Y, K, L, M, A), U=mat, MH=MY, 
					init=resYinit[,q], tau=vectau[q], tol=tol, stepsize=stepsize, nsteps=nsteps, nbatch=nbatch, seed=seed, WT=2)$coef
				resL[,,q][iter,] <- sgradfunq(idvar=idvar, X=cbind(L, K, A), U=mat, MH=ML, 
					init=resLinit[,q], tau=vectau[q], tol=tol, stepsize=stepsize, nsteps=nsteps, nbatch=nbatch, seed=seed, WT=2)$coef
				resM[,,q][iter,] <- sgradfunq(idvar=idvar, X=cbind(M, K, A), U=mat, MH=MM, 
					init=resMinit[,q], tau=vectau[q], tol=tol, stepsize=stepsize, nsteps=nsteps, nbatch=nbatch, seed=seed, WT=2)$coef
				resWT[,,q][iter,] <- sgradfunq(idvar=idvar, X=A, U=mat, MH=MW, init=resWTinit[,q], tau=vectau[q], 
					tol=tol, stepsize=stepsize, nsteps=nsteps, nbatch=nbatch, seed=seed, WT=1)$coef
				resW0[,,q][iter,] <- sgradfunq(idvar=idvar, X=A, U=mat, MH=MW1, init=resW0init[,q], tau=vectau[q], tol=tol, 
					stepsize=.1, nsteps=nsteps, nbatch=nbatch, seed=seed, WT=0)$coef
			}
		}
		##########################################################################################
		#Investment Estimation
		if (Mdraws==1){
			resILM <- lm(I~tensor.prod(MI, cbind(K, A, mat))-1)
			resIcoef <- as.numeric(coef(resILM))
			resIvar <- mean(resid(resILM)^2)
		} else {
			resIcoef <- sgradfune(idvar=idvar, X=cbind(I, K, A), U=mat, MH=MI, 
					init=resIinit[-length(resIinit)], tol=tol, stepsize=stepsize, nsteps=nsteps, nbatch=1, seed=seed)$coef
			#Estimate of Variance
			resIvar <- mean((apply(mat, 2, function(u) sum(I-tensor.prod(MI, cbind(K, A, u))%*%resIcoef)))^2)
			}
		resI[iter,] <- c(resIcoef, resIvar)
		###########################################################################################
		#Tail Parameters
		#Output
		yb <- expb(X=cbind(Y, K, L, M, A), U=mat, MH=MY, par1=resY[,,1][iter,], parL=resY[,,ntau][iter,], WT=2)
		yb1[iter,] <- yb$b1
		ybL[iter,] <- yb$bL
		#Labor
		lb <- expb(X=cbind(L, K, A), U=mat, MH=ML, par1=resL[,,1][iter,], parL=resL[,,ntau][iter,], WT=2)
		lb1[iter,] <- lb$b1
		lbL[iter,] <- lb$bL
		#Materials
		mb <- expb(X=cbind(M, K, A), U=mat, MH=MM, par1=resL[,,1][iter,], parL=resL[,,ntau][iter,], WT=2)
		mb1[iter,] <- mb$b1
		mbL[iter,] <- mb$bL
		#Productivity

		wtb <- expb(idvar=idvar, X=A, U=mat, MH=MW, par1=resWT[,,1][iter,], parL=resWT[,,ntau][iter,], WT=1)
		wtb1[iter,] <- wtb$b1
		wtbL[iter,] <- wtb$bL
		#Initial Productivity
		w0b <- expb(idvar=idvar, X=A, U=mat, MH=MW1, par1=resW0[,,1][iter,], parL=resW0[,,ntau][iter,], WT=0)
		w0b1[iter,] <- w0b$b1
		w0bL[iter,] <- w0b$bL
		#############################################################################################
		# Use Estimates as parameters in the MH Algorithm
		resYinit <- as.matrix(resY[iter,,])
		resLinit <- as.matrix(resL[iter,,])
		resMinit <- as.matrix(resM[iter,,])
		resWTinit <- as.matrix(resWT[iter,,])
		resW0init <- as.matrix(resW0[iter,,])
		resIinit <- as.matrix(resI[iter,])
		yb1init <- yb1[iter,]; ybLinit <- ybL[iter,]
		lb1init <- lb1[iter,]; lbLinit <- lbL[iter,]
		mb1init <- mb1[iter,]; mbLinit <- mbL[iter,]
		wtb1init <- wtb1[iter,]; wtbLinit <- wtbL[iter,]
		w0b1init <- w0b1[iter]; w0bLinit <- w0bL[iter]
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
	resW0mat <- sapply(1:ntau, function(x) colMeans(resW0[,,x][(maxiter/2):maxiter,]))
	resw0b1bLmat <- c(mean(w0b1[(maxiter/2):maxiter]), mean(w0bL[(maxiter/2):maxiter]))
	resImat <- colMeans(resI[(maxiter/2):maxiter,])
	return(list(resYmat, resLmat, resMmat, resWTmat, resW0mat, resImat, resyb1bLmat, reslb1bLmat, 
		resmb1bLmat, reswtb1bLmat, resw0b1bLmat))
}
#For testing purposes:
#EM Parameters
seed <- 123456
set.seed(seed)
ntau <- 11
maxiter <- 2
draws <- 50
Mdraws <- 1
#Monte Carlo Parameters
#Standard deviation of log wage process
siglnw <- 0.1
#Standard deviation of log material price process
siglnm <- 0.1
#Standard deviation of optimization error for labor
sigoptl <- 0.37
#Standard deviation of optimization error for materials
sigoptm <- 0.37
#Standard deviation of random shock to capital accumulation
sigoptk <- 0.2
#Standard deviation of random shock to investment
sigopti <- 0.37
# Number of Firms
n <- 1000
# Number of Time Periods
overallt <- 100
starttime <- 90
t <- overallt - starttime
#Production Function Parameters
beta0 <- 0
betal <- 0.3
betak <- 0.4
betam <- 0.3
betaa <- .01
betaw <- 1
#Epsilons and omega process
sigeps <- 0.1
sigomg <- 0.3 #standard deviation of omega
rho <- 0.7 #AR(1) coefficient for omega
sigxi <- sqrt((1-rho^2)*sigomg^2)
#Matrices to store data
lnkdata <- matrix(0, n, overallt) #ln(capital)
lnldata <- matrix(0, n, overallt) #ln(labor)
lnmdata <- matrix(0, n, overallt) #ln(intermediate input)
lnydata <- matrix(0, n, overallt) #ln(output)
omgdata <- matrix(0, n, overallt) #omega(t)
agedata <- matrix(0, n, overallt) #age(t)
#Location Scale Parameters for Output
eta0 <- 1
etak <- 0.3
etal <- -0.4
etam <- 0.4
#Location Scale Parameters for Labor
hetl0 <- 1
hetlk <- -0.2
hetlw <- 0.3
#Location Scale Parameters for Materials
hetm0 <- 1
hetmk <- -0.2
hetmw <- 0.3
#Specification for Production Shock Distribution
etadata <- matrix(rnorm(n*overallt, 0, sigeps), nrow=n, ncol=overallt)
#Specification for Labor Input Shock Distribution
varepsdata <- matrix(rnorm(n*overallt, 0, sigoptl), nrow=n, ncol=overallt)
#Specification for Materials Input Shock Distribution
epsdata <- matrix(rnorm(n*overallt, 0, sigoptm), nrow=n, ncol=overallt)
#Specification for Investment Input Shock Distribution
iotadata <- matrix(rnorm(n*overallt, 0, sigopti), nrow=n, ncol=overallt)
#Specification for labor cost (not serially correlated)
lnwdata <- matrix(rnorm(n*overallt, 0, siglnw), nrow=n, ncol=overallt)
#Specification for materials cost (not serially correlated)
lnpmdata <- matrix(rnorm(n*overallt, 0, siglnm), nrow=n, ncol=overallt)
#Period 0 values of omega 
omgdata0 <- matrix(rnorm(n,0,sigomg),nrow=n,ncol=1)
#Period 0 values of age
agedata0 <- floor(as.matrix(runif(n, 1, 10)))
agedata[,1] <- agedata0+1
#Period 1 values of omega 
omgdata[,1] <- betaa*agedata0[,1]+rho*omgdata0[,1]+matrix(rnorm(n,0,sigomg),nrow=n,ncol=1)
#Simulate values of omega for rest of time periods
for (s in 2:overallt){
  agedata[,s] <- agedata[,s-1]+1
  omgdata[,s] <- betaa*agedata[,s]+rho*omgdata[,s-1] + matrix(rnorm(n,0,sigomg),nrow=n,ncol=1)

}
#Intital ln(capital) level (close to 0)
lnkdata[,1] <- matrix(-100,n,1)
#Discount Rate for DP problem
disc <- 0.95
#Depreciation rate of capital
delta <- 0.2
#Capital Adjustment Cost Parameter
badj <- 1
#Simplifying components of the optimal investment rule
#Square bracket component
squarebracketterm <- exp(0.5*betal^2*sigoptl^2)*exp(0.5*betam^2*sigoptm^2)
#Constant term in front of sum including squarebracketterm
const1 <- disc*betak*betal^(betal/betak)*betam^(betam/betak)*(exp(beta0)^(1/(betak)))*
squarebracketterm/badj
vec1 <- (disc*(1-delta))^seq(100)
vec2 <- cbind(sigxi^2 * 0,cumsum(rho^(2*(seq(100)-1))))
expterm1 <- exp(0.5*(betaw/(betak))^2*rho^2*
((sigxi^2)*rho^(2*seq(100))+vec2))
#Compute Optimal Investment and Capital stock for all firms over time
investmat <- matrix(NA, n, overallt)
for (i in 1:n){
    for (s in 1:overallt){
      expterm2 <- exp((betaw/betak)*omgdata[i,s]*rho^(seq(100)))
  	  expterm3 <- exp(0.5*betal^2*((1+hetlk*lnkdata[i,s]+hetlw*omgdata[i,s])*sigoptl)^2)
  	  expterm4 <- exp(0.5*betam^2*((1+hetmk*lnkdata[i,s]+hetmw*omgdata[i,s])*sigoptm)^2)
      investmat[i,s] <- const1*sum(vec1*expterm1*expterm2*expterm3*expterm4)
      if (s >= 2){
        lnkdata[i,s] <- log((1-delta)*exp(lnkdata[i,s-1])+investmat[i,s-1])+rnorm(1, mean=0, sd=sigoptk)
      }
    }
  }

#Generate levels of labor input
lnldata <- (beta0+(1-betam)*log(betal)+betak*lnkdata+betam*log(betam)-betam*lnpmdata-(1-betam)*lnwdata+betaw*omgdata)/(1-betam-betal)
#Adding Optimization Error
hetl <- hetl0+hetlk*lnkdata+hetlw*omgdata
lnldata <- lnldata + hetl*matrix(rnorm(n*overallt,0,sigoptl),n,overallt)
#Generate levels of material input
lnmdata <- (beta0+log(betam)+betak*lnkdata+betal*(log(betal)-lnwdata-log(betam))+betaw*omgdata-(1-betal)*lnpmdata)/(1-betal-betam)
hetm <- hetm0+hetmk*lnkdata+hetmw*omgdata
lnmdata <- lnmdata + hetl*matrix(rnorm(n*overallt,0,sigoptl),n,overallt)
#Add the random shock to investmnet
investmat <- log(investmat)+iotadata
#Output
#Specifies the form of heteroskedasticity
hety <- eta0+etal*lnldata+etak*lnkdata
lnydata <- beta0 + betal*lnldata + betak*lnkdata + betam*lnmdata + omgdata + hety*etadata

#Stack data across firms (all the data)
Capital <- c(t(lnkdata[,(starttime+1):overallt]))
Labor <- c(t(lnldata[,(starttime+1):overallt]))
Materials <- c(t(lnmdata[,(starttime+1):overallt]))
Investment <- c(t(investmat[,(starttime+1):overallt]))
Wage <- c(t(lnwdata[,(starttime+1):overallt]))
Output <- c(t(lnydata[,(starttime+1):overallt]))
Productivity <- c(t(omgdata[,(starttime+1):overallt]))
Epsilon <- c(t(etadata[,(starttime+1):overallt]))

#Stack data across firms (lagged data)
Capital_Lag_1 <- c(t(lnkdata[,(starttime+1):(overallt-1)]))
Labor_Lag_1 <- c(t(lnldata[,(starttime+1):(overallt-1)]))
Labor_Lag_2 <- c(t(lnldata[,(starttime):(overallt-2)]))
Materials_Lag_1 <- c(t(lnmdata[,(starttime+1):(overallt-1)]))
Investment_Lag_1 <- c(t(investmat[,(starttime+1):overallt-1]))
Wage_Lag_1 <- c(t(lnwdata[,(starttime+1):(overallt-1)]))
Output_Lag_1 <- c(t(lnydata[,(starttime+1):(overallt-1)]))
Productivity_Lag_1 <- c(t(omgdata[,(starttime+1):(overallt-1)]))

#Stack data across firms (contemporaneous data)
Capital_Con <- c(t(lnkdata[,(starttime+2):(overallt)]))
Labor_Con <- c(t(lnldata[,(starttime+2):(overallt)]))
Materials_Con <- c(t(lnmdata[,(starttime+2):(overallt)]))
Investment_Con <- c(t(investmat[,(starttime+2):overallt]))
Wage_Con <- c(t(lnwdata[,(starttime+2):(overallt)]))
Output_Con <- c(t(lnydata[,(starttime+2):(overallt)]))
Productivity_Con <- c(t(omgdata[,(starttime+2):(overallt)]))
Age_Con <- c(t(agedata[,(starttime+2):(overallt)]))
#Create ID and Time variables
idvar <- rep(1:n, each=(overallt-(starttime+1)))
timevar <- rep(1:(overallt-(starttime+1)), n)


overall.start.time <- Sys.time()
results <- Production_EM(ntau=ntau, idvar=idvar, timevar=timevar, Y=Output_Con, K=Capital_Con, L=Labor_Con, M=Materials_Con, I=Investment_Con, A=Age_Con, maxiter=maxiter, draws=draws, Mdraws=Mdraws, seed=seed)
# save(test, file='test.Rdata')
print(Sys.time()-overall.start.time)



