setwd('/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Functions')
require(quantreg)
source('Hermite.R')
source('Posterior.R')
source('Moment.R')
Production_EM <- function(ntau, Y, K, L, M, I, A, N, T, maxiter, draws, Mdraws){
	set.seed(123456)
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
	#Grid of Taus
	vectau <- seq(1/(ntau+1), ntau/(ntau+1), by=1/(ntau+1))
	#Degrees of the Univariate Hermite Polynomials
	MK <- 2
	ML <- 2
	MM <- 2
	MA <- 1
	MW <- 1
	MI <- 1
	MH <- c(MK, ML, MM, MI, MA, MW)
	TFPreg <- lm(Y~K+L+M+A)
	TFP <- exp(resid(TFPreg))
	TFPNoise <- as.matrix(TFP+rlnorm(nrow(Y)))
	############################################################################
	############ Initialization for EM Algorithm ###############################
	############################################################################
	#Initialize Paramater Values for Output
	Yinit <- tensor.prod(c(MK, ML, MM, MA, MW), cbind(K, L, M, A, I))
	resYinit <- matrix(rq(Y~Yinit-1, tau=vectau)$coef, nrow=(MK+1)*(ML+1)*(MM+1)*(MA+1)*(MW+1), ncol=ntau)
	#Initialize Paramater Values for Labor
	Linit <- tensor.prod(c(MK, MA, MW), cbind(K, A, I))
	resLinit <- matrix(rq(L~Linit-1, tau=vectau)$coef, nrow=(MK+1)*(MA+1)*(MW+1), ncol=ntau)
	#Initialize Parameter Values for Labor
	Minit <- tensor.prod(c(MK, MA, MW), cbind(K, A, I))
	resMinit <- matrix(rq(M~Minit-1, tau=vectau)$coef, nrow=(MK+1)*(MA+1)*(MW+1), ncol=ntau)
	#Initialize Parameter Values for Omega_{t}
	Wtinit <- tensor.prod(c(MA, MW), cbind(A[-seq(1, N*T, by=T)], TFP[-seq(T, N*T, by=T)]))
	resWTinit <- matrix(rq(TFP[-seq(1, N*T, by=T)]~Wtinit-1, tau=vectau)$coef, nrow=(MA+1)*(MW+1), ncol=ntau)
	#Initialize Parameter Values for Omega_{0}
	W0init <- tensor.prod(MA, as.matrix(A[seq(1, N*T, by=T)]))
	resW0init <- matrix(rq(TFP[seq(1, N*T, by=T)]~W0init-1, tau=vectau)$coef, nrow=(MA+1), ncol=ntau)
	#Initialize Paramater Values for Investment (Nonlinear Regression Model)
	Iinit <- tensor.prod(c(MK, MA, MW), cbind(K, A, M))
	resIcoef <- lm(I~Iinit-1)
	#Vector of Investment coefficients and estimate of  variance
	resIinit <- c(as.numeric(coef(resIcoef)), mean(resid(resIcoef)^2))
	#Initial Parameter Values for Laplace Parameters
	yb1init <- ybLinit <- 10
	lb1init <- lbLinit <- 10
	mb1init <- mbLinit <- 10
	wtb1init <- wtbLinit <- 10
	w0b1init <- w0bLinit <- 10
	############################################################################
	############################################################################
	############################## EM Algorithm ################################
	############################################################################
	############################################################################
	#Variance Random Walk proposal
	varRW <- .05
	#Matrices for Storing Results
	resY <- array(0, c(maxiter, (MK+1)*(ML+1)*(MM+1)*(MA+1)*(MW+1), ntau))
	yb1 <- array(0, c(maxiter, 1))
	ybL <- array(0, c(maxiter, 1))
	resL <- array(0, c(maxiter, (MK+1)*(MA+1)*(MW+1), ntau))
	lb1 <- array(0, c(maxiter, 1))
	lbL <- array(0, c(maxiter, 1))
	resM <- array(0, c(maxiter, (MK+1)*(MA+1)*(MW+1), ntau))
	mb1 <- array(0, c(maxiter, 1))
	mbL <- array(0, c(maxiter, 1))
	resI <- array(0, c(maxiter, (MK+1)*(MA+1)*(MW+1)+1))
	resWT <- array(0, c(maxiter, (MA+1)*(MW+1), ntau))
	wtb1 <- array(0, c(maxiter, 1))
	wtbL <- array(0, c(maxiter, 1))
	resW0 <- array(0, c(maxiter, (MA+1), ntau))
	w0b1 <- array(0, c(maxiter, 1))
	w0bL <- array(0, c(maxiter, 1))
	#Begin EM Algorithm
	for (iter in 1:maxiter){
		print(iter)
		print(Sys.time()-overall.start.time)
		#This sets how many MH "burn in" samples
		r <- -draws
		#Initialize Matrix for Averaging over Accepted Draws
		mat <- matrix(0, N*T, Mdraws)
		#Initial Parameter Values 
		resinit <- list(resYinit=resYinit, yb1init=yb1init, ybLinit=ybLinit, resLinit=resLinit, lb1init=lb1init, 
			lbLinit=lbLinit, resMinit=resMinit, mb1init=mb1init, mbLinit=mbLinit, 
			resWTinit=resWTinit, wtb1init=wtb1init, wtbLinit=wtbLinit, resW0init=resW0init, w0b1init=w0b1init, 
			w0bLinit=w0bLinit, resIinit=resIinit)
		#Initial Guess for Unobservables
		un <- as.matrix(rnorm(N*T, 0, sqrt(varRW)))
		un_dens <- posterior(Y=Y, K=K, L=L, M=M, A=A, I=I, MH=MH, ntau=ntau, vectau=vectau, N=N, T=T, matdraw=un, resinit=resinit)
		for (j in r:Mdraws){
			print(j)
			#Proposal Distribution
			try_un <- as.matrix(rnorm(N*T, 0, sqrt(varRW)))
			y <- un+try_un
			#Acceptance Probability 
			loglik <- posterior(Y=Y, K=K, L=L, M=M, A=A, I=I, MH=MH, ntau=ntau, vectau=vectau, N=N, T=T, matdraw=y, resinit=resinit)
			loga <- loglik-un_dens
			#Create a set of N indices for acceptance rule
			Nindex <- which(log(runif(N))<loga)
			#Create sequence of vectors for accepting each time period of each accepted firm draw
			NTindex <- c(sapply(1:length(Nindex), function(n) seq((Nindex[n]-1)*T+1, Nindex[n]*T)))
			#Update unobservables that satisfy acceptance rule according to the RW process
			un[NTindex] <- y[NTindex]
			un_dens <- loglik
			if (j>0){
				mat[,j] <- un
			}
		}
		#Some initializations, reformatting for optimization steps
		#Create contemporary and lag values
		WTCon <- as.matrix(mat[-seq(1, N*T, by=T),])
		WTLag <- as.matrix(mat[-seq(T, N*T, by=T),])
		WT0 <- as.matrix(mat[seq(1, N*T, by=T),])
		ACon <- as.matrix(A[-seq(1, N*T, by=T)])
		A0 <- as.matrix(A[seq(1, N*T, by=T)])
		#Initialize Gradient Descent Algorithms
		#Tolerance for gradient (close to zero)
		tol <- 1e-2
		#Stepsize for gradient
		stepsize <- 100
		#Number of gradient descent iterations
		nsteps <- 500
		for (tau in 1:ntau){
				resY[,,tau][iter,] <- gradfun(Y=Y, X=cbind(K, L, M, A), U=mat, MH=c(MK, ML, MM, MA, MW), 
					init=resYinit[,tau], tau=tau, tol=tol, stepsize=stepsize, nsteps=nsteps, N=N, T=T)$coef
				resL[,,tau][iter,] <- gradfun(Y=L, X=cbind(K, A), U=mat, MH=c(MK, MA, MW), 
					init=resLinit[,tau], tau=tau, tol=tol, stepsize=stepsize, nsteps=nsteps, N=N, T=T)$coef
				resM[,,tau][iter,] <- gradfun(Y=M, X=cbind(K, A), U=mat, MH=c(MK, MA, MW), 
					init=resMinit[,tau], tau=tau, tol=tol, stepsize=stepsize, nsteps=nsteps, N=N, T=T)$coef
				resWT[,,tau][iter,] <- gradfunWT(X=A, U=mat, MH=c(MA, MW), init=resWTinit[,tau], tau=tau, 
					tol=tol, stepsize=stepsize, nsteps=nsteps, N=N, T=T)$coef
				resW0[,,tau][iter,] <- gradfunW0(X=A, U=mat, MH=MA, init=resW0init[,tau], tau=tau, tol=tol, 
					stepsize=stepsize, nsteps=nsteps, N=N, T=T)$coef
			}
			#Investment Estimation
			resIcoef <- gradfun(Y=I, X=cbind(K, A), U=mat, MH=c(MK, MA, MW), 
					init=resIinit[-length(resInit)], tol=tol, stepsize=stepsize, nsteps=nsteps, N=N, T=T)$coef
			#Estimate of Variance
			resIvar <- mean((apply(mat, 2, function(u) sum(I-tensor.prod(c(MK, MA, MW), cbind(K, A, u))%*%resIcoef)))^2)
			resI[iter,] <- c(resIcoef, resIvar)
			#Tail Parameters
			#Output
			yb1[iter,] <- -mean(sapply(1:Mdraws, function(j) sum(Y<=tensor.prod(c(MK, ML, MM, MW, MA), 
				cbind(K, L, M, mat[,j], A))%*%resY[,,1][iter,])/sum((Y-tensor.prod(c(MK, ML, MM, MW, MA), 
					cbind(K, L, M, mat[,j], A))%*%resY[,,1][iter,])*(Y<=tensor.prod(c(MK, ML, MM, MW, MA), 
					cbind(K, L, M, mat[,j], A))%*%resY[,,1][iter,]))))
			ybL[iter,] <- mean(sapply(1:Mdraws, function(j) sum(Y>tensor.prod(c(MK, ML, MM, MW, MA), 
				cbind(K, L, M, mat[,j], A))%*%resY[,,ntau][iter,])/sum((Y-tensor.prod(c(MK, ML, MM, MW, MA), 
					cbind(K, L, M, mat[,j], A))%*%resY[,,ntau][iter,])*(Y>tensor.prod(c(MK, ML, MM, MW, MA), 
					cbind(K, L, M, mat[,j], A))%*%resY[,,ntau][iter,]))))
			#Labor
			lb1[iter,] <- -mean(sapply(1:Mdraws, function(j) sum(L<=tensor.prod(c(MK, MW, MA), 
				cbind(K, mat[,j], A))%*%resL[,,1][iter,])/sum((L-tensor.prod(c(MK, MW, MA), 
					cbind(K, mat[,j], A))%*%resL[,,1][iter,])*(L<=tensor.prod(c(MK, MW, MA), 
					cbind(K, mat[,j], A))%*%resL[,,1][iter,]))))
			lbL[iter,] <- mean(sapply(1:Mdraws, function(j) sum(L>tensor.prod(c(MK, MW, MA), 
				cbind(K, mat[,j], A))%*%resL[,,ntau][iter,])/sum((L-tensor.prod(c(MK, MW, MA), 
					cbind(K, mat[,j], A))%*%resL[,,ntau][iter,])*(L>tensor.prod(c(MK, MW, MA), 
					cbind(K, mat[,j], A))%*%resL[,,ntau][iter,]))))
			#Materials
			mb1[iter,] <- -mean(sapply(1:Mdraws, function(j) sum(M<=tensor.prod(c(MK, MW, MA), 
				cbind(K, mat[,j], A))%*%resM[,,1][iter,])/sum((M-tensor.prod(c(MK, MW, MA), 
					cbind(M, mat[,j], A))%*%resM[,,1][iter,])*(M<=tensor.prod(c(MK, MW, MA), 
					cbind(M, mat[,j], A))%*%resM[,,1][iter,]))))
			mbL[iter,] <- mean(sapply(1:Mdraws, function(j) sum(M>tensor.prod(c(MK, MW, MA), 
				cbind(K, mat[,j], A))%*%resM[,,ntau][iter,])/sum((M-tensor.prod(c(MK, MW, MA), 
					cbind(K, mat[,j], A))%*%resM[,,ntau][iter,])*(M>tensor.prod(c(MK, MW, MA), 
					cbind(K, mat[,j], A))%*%resM[,,ntau][iter,]))))
			#Productivity
			wtb1[iter,] <- -mean(sapply(1:Mdraws, function(j) sum(WTCon[,j]<=tensor.prod(c(MW, MA), 
				cbind(WTLag[,j], ACon))%*%resWT[,,1][iter,])/sum((WTCon[,j]-tensor.prod(c(MW, MA), 
					cbind(WTLag[,j], ACon))%*%resWT[,,1][iter,])*(WTCon[,j]<=tensor.prod(c(MW, MA), 
					cbind(WTLag[,j], ACon))%*%resWT[,,1][iter,]))))
			wtbL[iter,] <- mean(sapply(1:Mdraws, function(j) sum(WTCon[,j]>tensor.prod(c(MW, MA), 
				cbind(WTLag[,j], ACon))%*%resWT[,,ntau][iter,])/sum((WTCon[,j]-tensor.prod(c(MW, MA), 
					cbind(WTLag[,j], ACon))%*%resWT[,,ntau][iter,])*(WTCon[,j]>tensor.prod(c(MW, MA), 
					cbind(WTLag[,j], ACon))%*%resWT[,,ntau][iter,]))))
			#Initial Productivity
			w0b1[iter,] <- -mean(sapply(1:Mdraws, function(j) sum(WT0[,j]<=tensor.prod(MA, A0)%*%resW0[,,1][iter,])/
					sum((WT0[,j]-tensor.prod(MA, A0)%*%resW0[,,1][iter,])*(WT0[,j]<=tensor.prod(MA, A0)%*%resW0[,,1][iter,]))))
			w0bL[iter,] <- mean(sapply(1:Mdraws, function(j) sum(WT0[,j]>tensor.prod(MA, A0)%*%resW0[,,ntau][iter,])))/
					mean(sapply(1:Mdraws, function(j) sum((WT0[,j]-tensor.prod(MA, A0)%*%resW0[,,ntau][iter,])*(WT0[,j]>tensor.prod(MA, A0)%*%resW0[,,ntau][iter,]))))
			#Use Estimates as parameters in the MH Algorithm
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
			w0b1init <- w0b1[iter]; w0bLinit <- w0bL[iter,]

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
ntau <- 11
maxiter <- 2
draws <- 50
Mdraws <- 5
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
#Period 0 values of omega and ln(wage)
omgdata0 <- matrix(rnorm(n,0,sigomg),nrow=n,ncol=1)
#Period 1 values of omega 
omgdata[,1] <- rho*omgdata0[,1]+matrix(rnorm(n,0,sigomg),nrow=n,ncol=1)
#Simulate values of omega for rest of time periods
for (s in 2:overallt){
  omgdata[,s] <- rho*omgdata[,s] + matrix(rnorm(n,0,sigomg),nrow=n,ncol=1)

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
Initial_Age <- floor(as.matrix(runif(n, 1, 10)))
Age_Con <- c(sapply(1:n, function(i) seq(Initial_Age[i,], Initial_Age[i,]+8, by=1)))


overall.start.time <- Sys.time()
test <- Production_EM(ntau=ntau, Y=Output_Con, K=Capital_Con, L=Labor_Con, M=Materials_Con, I=Investment_Con, A=Age_Con, N=n, T=t-1, maxiter=maxiter, draws=draws, Mdraws=Mdraws)
# save(test, file='test.Rdata')
print(Sys.time()-overall.start.time)









