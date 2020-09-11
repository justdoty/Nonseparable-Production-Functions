setwd('/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Functions')
require(quantreg)
source('Hermite.R')
source('Posterior.R')
source('Moment.R')
set.seed(123456)
Production_EM <- function(ntau, Y, K, L, M, I, A, N, T, maxiter, draws, Mdraws){
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
	TFPreg <- lm(Y~K+L+M+A)
	TFP <- exp(resid(lm))
	TFPNoise <- TFP+rlnorm(nrow(Y))
	############################################################################
	############ Initialization for EM Algorithm ###############################
	############################################################################
	#Initialize Paramater Values for Output
	H.Yinit <- tensor.prod(c(MK, ML, MM, MW, MA), cbind(K, L, M, TFPNoise, A))
	resYinit <- matrix(rq(Y~H.Yinit-1, tau=vectau)$coef, nrow=ntau, ncol=(MK+1)*(ML+1)*(MM+1)*(MW+1)*(MA+1))
	#Initialize Paramater Values for Labor
	H.Linit <- tensor.prod(c(MK, MW, MA), cbind(K, TFP, A))
	resLinit <- matrix(rq(L~H.Linit-1, tau=vectau)$coef, nrow=ntau, ncol=(MK+1)*(MW+1)*(MA+1))
	#Initialize Parameter Values for Labor
	H.Minit <- tensor.prod(c(MK, MW, MA), cbind(K, TFP, A))
	resMinit <- matrix(rq(M~H.Minit-1, tau=vectau)$coef, nrow=ntau, ncol=(MK+1)*(MW+1)*(MA+1))
	#Initialize Parameter Values for Omega_{t}
	H.Wtinit <- tensor.prod(c(MW, MA), cbind(TFP[-seq(T, N*T, by=N)], A))
	resWTinit <- matrix(rq(TFP[-seq(1, N*T, by=N)]~H.Wtinit-1, tau=vectau)$coef, nrow=ntau, ncol=(MW+1)*(MA+1))
	#Initialize Parameter Values for Omega_{0}
	H.W0init <- tensor.prod(MA, A)
	resW0init <- matrix(rq(TFP[seq(1:N*T, by=N),]~H.W0init-1, tau=vectau)$coef, nrow=ntau, ncol=(MA+1))
	#Initialize Paramater Values for Investment (Nonlinear Regression Model)
	H.Iinit <- tensor.prod(c(MK, MW, MA), cbind(K, TFP, A))
	resIcoef <- lm(L~H.Linit-1)
	#Vector of Investment coefficients and estimate of  variance
	resIinit <- as.matrix(c(as.numeric(coef(resIcoef)), mean(resid(resIcoef)^2)))
	#Initial Parameter Values for Laplace Parameters
	yb1 <- ybL <- 10
	lb1 <- lbL <- 10
	mb1 <- mbL <- 10
	wtb1 <- wtbL <- 10
	w0b1 <- w0bL <- 10
	############################################################################
	############################################################################
	############################## EM Algorithm ################################
	############################################################################
	############################################################################
	#Variance Random Walk proposal
	varRW <- .05
	#Matrices for Storing Results
	resY <- array(0, c(ntau, (MK+1)*(ML+1)*(MM+1)*(MW+1)*(MA+1), maxiter))
	yb1 <- array(0, c((MK+1)*(ML+1)*(MM+1)*(MW+1)*(MA+1), maxiter))
	ybL <- array(0, c((MK+1)*(ML+1)*(MM+1)*(MW+1)*(MA+1), maxiter))
	resL <- array(0, c(ntau, (MK+1)*(MW+1)*(MA+1), maxiter))
	lb1 <- array(0, c((MK+1)*(MW+1)*(MA+1), maxiter))
	lbL <- array(0, c((MK+1)*(MW+1)*(MA+1), maxiter))
	resM <- array(0, c(ntau, (MK+1)*(MW+1)*(MA+1), maxiter))
	mb1 <- array(0, c((MK+1)*(MW+1)*(MA+1), maxiter))
	mbL <- array(0, c((MK+1)*(MW+1)*(MA+1), maxiter))
	resI <- array(0, c((MK+1)*(MW+1)*(MA+1)+1, maxiter))
	resWT <- array(0, c(ntau, (MW+1)*(MA+1), maxiter))
	wtb1 <- array(0, c((MW+1)*(MA+1), maxiter))
	wtbL <- array(0, c((MW+1)*(MA+1), maxiter))
	resW0 <- array(0, c(ntau, (MA+1), maxiter))
	w0b1 <- array(0, c(MA, maxiter))
	w0bL <- array(0, c(MA, maxiter))
	#Begin EM Algorithm
	for (iter in 1:maxiter){
		#This sets how many MH "burn in" samples
		r <- -draws
		#Initialize Matrix for Averaging over Accepted Draws
		Mat <- matrix(0, nrow(N*T), Mdraws)
		#Initial Guess for Unobservables
		un <- rnorm(N*T, 0, sqrt(varRW))
		for (j in r:Mdraws){
			#Proposal Distribution
			try_un <- rnorm(N*T, 0, sqrt(varRW))
			#Acceptance Probability 
			try_dens <- posterior(try_un, reslist)/posterior(un, reslist)
			#Create a set of indices for acceptance rule
			aindex <- which(runif(N*T)<try_dens)
			#Update unobservables that satisfy acceptance rule according to the RW process
			un[aindex] <- un[aindex]+try_un[aindex]
			if (j>0){
				Mat[,j] <- un
			}
		}
		WTCon <- un[-seq(1, N*T, by=N),]
		WTLag <- un[-seq(T, N*T, by=N),]
		WT0 <- un[seq(1:N*T, by=N),]
		ACon <- A[-seq(1, N*T, by=N)]
		A0 <- A[seq(1:N*T, by=N)]
		for (tau in 1:ntau){
				resY[,,i][j,] <- optim(resYinit[j,], function(c) mean(sapply(1:Mdraws, function(j) 
					qmoment(Y=Y, X=tensor.prod(c(MK, ML, MM, MW, MA), cbind(K, L, M, mat[,j], A)), c, vectau[tau]))))
				resL[,,i][j,] <- optim(resLinit[j,], function(c) mean(sapply(1:Mdraws, function(j) 
					qmoment(Y=L, X=tensor.prod(c(MK, MW, MA), cbind(K, mat[,j], A)), c, vectau[tau]))))
				resM[,,i][j,] <- optim(resMinit[j,], function(c) mean(sapply(1:Mdraws, function(j) 
					qmoment(Y=M, X=tensor.prod(c(MK, MW, MA), cbind(K, mat[,j], A)), c, vectau[tau]))))
				resWT[,,i][j,] <- optim(resWTinit[j,], function(c) mean(sapply(1:Mdraws, function(j) 
					qmoment(Y=WTCon[,j], X=tensor.prod(c(MW, MA), cbind(WTLag[,j], ACon)), c, vectau[tau]))))
				resW0[,,i][j,] <- optim(resW0init[j,], function(c) mean(sapply(1:Mdraws, function(j) 
					qmoment(Y=WT0[,j], X=tensor.prod(MA, A0), c, vectau[tau]))))
			}
			#Investment Estimation
			resIcoef <- optim(resIinit[-length(resIinit)], function(c) mean(sapply(1:Mdraws, function(j) 
					qmoment(Y=I, X=tensor.prod(c(K, MW, MA), cbind(K, mat[,j], A)), c))))
			#Estimate of Variance
			resIvar <- mean(sapply(1:Mdraws, function(j) mean((I-tensor.prod(c(K, MW, MA), cbind(K, mat[,j], A))%*%resIcoef))))^2
			resImat[,i] <- c(resIcoef, resIvar)
			#Tail Parameters
			#Output
			yb1[,i] <- -mean(sapply(1:Mdraws, function(j) sum(Y<=tensor.prod(c(MK, ML, MM, MW, MA), 
				cbind(K, L, M, mat[,j], A))%*%resY[,,i][1,])/sum((Y-tensor.prod(c(MK, ML, MM, MW, MA), 
					cbind(K, L, M, mat[,j], A))%*%resY[,,i][1,])*(Y<=tensor.prod(c(MK, ML, MM, MW, MA), 
					cbind(K, L, M, mat[,j], A))%*%resY[,,i][1,]))))
			ybL[,i] <- mean(sapply(1:Mdraws, function(j) sum(Y>tensor.prod(c(MK, ML, MM, MW, MA), 
				cbind(K, L, M, mat[,j], A))%*%resY[,,i][ntau,])/sum((Y-tensor.prod(c(MK, ML, MM, MW, MA), 
					cbind(K, L, M, mat[,j], A))%*%resY[,,i][ntau,])*(Y>tensor.prod(c(MK, ML, MM, MW, MA), 
					cbind(K, L, M, mat[,j], A))%*%resY[,,i][ntau,]))))
			#Labor
			lb1[,i] <- -mean(sapply(1:Mdraws, function(j) sum(L<=tensor.prod(c(MK, MW, MA), 
				cbind(K, mat[,j], A))%*%resL[,,i][1,])/sum((L-tensor.prod(c(MK, MW, MA), 
					cbind(K, mat[,j], A))%*%resL[,,i][1,])*(L<=tensor.prod(c(MK, MW, MA), 
					cbind(K, mat[,j], A))%*%resL[,,i][1,]))))
			lbL[,i] <- mean(sapply(1:Mdraws, function(j) sum(L>tensor.prod(c(MK, MW, MA), 
				cbind(K, mat[,j], A))%*%resL[,,i][ntau,])/sum((L-tensor.prod(c(MK, MW, MA), 
					cbind(K, mat[,j], A))%*%resL[,,i][ntau,])*(L>tensor.prod(c(MK, MW, MA), 
					cbind(K, mat[,j], A))%*%resL[,,i][ntau,]))))
			#Materials
			mb1[,i] <- -mean(sapply(1:Mdraws, function(j) sum(M<=tensor.prod(c(MK, MW, MA), 
				cbind(K, mat[,j], A))%*%resM[,,i][1,])/sum((M-tensor.prod(c(MK, MW, MA), 
					cbind(M, mat[,j], A))%*%resM[,,i][1,])*(M<=tensor.prod(c(MK, MW, MA), 
					cbind(M, mat[,j], A))%*%resM[,,i][1,]))))
			mbL[,i] <- mean(sapply(1:Mdraws, function(j) sum(M>tensor.prod(c(MK, MW, MA), 
				cbind(K, mat[,j], A))%*%resM[,,i][ntau,])/sum((M-tensor.prod(c(MK, MW, MA), 
					cbind(K, mat[,j], A))%*%resM[,,i][ntau,])*(M>tensor.prod(c(MK, MW, MA), 
					cbind(K, mat[,j], A))%*%resM[,,i][ntau,]))))
			#Productivity
			wtb1[,i] <- -mean(sapply(1:Mdraws, function(j) sum(WTCon[,j]<=tensor.prod(c(MW, MA), 
				cbind(WTlag[,j], A))%*%resWT[,,i][1,])/sum((WTCon[,j]-tensor.prod(c(MW, MA), 
					cbind(WTlag[,j], A))%*%resWT[,,i][1,])*(WTCon[,j]<=tensor.prod(c(MW, MA), 
					cbind(WTLag[,j], A))%*%resWT[,,i][1,]))))
			wtbL[,i] <- mean(sapply(1:Mdraws, function(j) sum(WTCon[,j]>tensor.prod(c(MW, MA), 
				cbind(WTlag[,j], A))%*%resWT[,,i][ntau,])/sum((WTCon[,j]-tensor.prod(c(MW, MA), 
					cbind(WTlag[,j], A))%*%resWT[,,i][ntau,])*(WTCon[,j]>tensor.prod(c(MW, MA), 
					cbind(WTLag[,j], A))%*%resWT[,,i][ntau,]))))
			#Initial Productivity
			w0b1[,i] <- -mean(sapply(1:Mdraws, function(j) sum(WT0[,j]<=tensor.prod(MA, A)%*%resW0[,,i][1,])/
					sum((WT0[,j]-tensor.prod(MA, A)%*%resW0[,,i][1,])*(WT0[,j]<=tensor.prod(MA, A)%*%resW0[,,i][1,]))))
			wtbL[,i] <- mean(sapply(1:Mdraws, function(j) sum(WT0[,j]>tensor.prod(MA, A)%*%resW0[,,i][ntau,])/
					sum((WT0[,j]-tensor.prod(MA, A)%*%resW0[,,i][ntau,])*(WT0[,j]>tensor.prod(MA, A)%*%resW0[,,i][ntau,]))))
			#Use Estimates as parameters in the MH Algorithm
			resYinit <- resY[,,i]
			resLinit <- resL[,,i]
			resMinit <- resM[,,i]
			resWTinit <- resWT[,,i]
			resW0init <- resW0[,,i]
			resIinit <- resImat[,i]
			yb1 <- yb1[,i]; ybL <- ybL[,i]
			lb1 <- lb1[,i]; lbL <- lbL[,i]
			mb1 <- mb1[,i]; mbL <- mbL[,i]
			wtb1 <- wtb1[,i]; wtbL <- wtbL[,i]
			w0b1 <- w0b1[,i]; w0bL <- w0bL[,i]
			#Put them in a list to be used in Posterior.R in the next iteration
			reslist <- list(resYinit, resLinit, resMinit, resWTinit, resW0init, resIinit, yb1, ybL, lb1, lbL, 
				mb1, mbL, wtb1, wtbL, w0b1, w0bL)
		}
	resYmat <- rbind(rowMeans(yb1[,maxiter/2:maxiter]), rowMeans(resY[,,maxiter/2:maxiter], dims=2), rowMeans(ybL[,maxiter/2:maxiter]))
	resLmat <- rbind(rowMeans(lb1[,maxiter/2:maxiter]), rowMeans(resL[,,maxiter/2:maxiter], dims=2), rowMeans(lbL[,maxiter/2:maxiter]))
	resMmat <- rbind(rowMeans(mb1[,maxiter/2:maxiter]), rowMeans(resM[,,maxiter/2:maxiter], dims=2), rowMeans(mbL[,maxiter/2:maxiter]))
	resWTmat <- rbind(rowMeans(wtb1[,maxiter/2:maxiter]), rowMeans(resWT[,,maxiter/2:maxiter], dims=2), rowMeans(wtbL[,maxiter/2:maxiter]))
	resW0mat <- rbind(rowMeans(w0b1[,maxiter/2:maxiter]), rowMeans(resW0[,,maxiter/2:maxiter], dims=2), rowMeans(w0bL[,maxiter/2:maxiter]))
	resImat <- rowMeans(resI)
	return(list(resYmat, resLmat, resMmat, resWTmat, resW0mat, resImat))
}
#For testing purposes:
ntau <- 11
maxiter <- 500
draws <- 500
Mdraws <- 1

set.seed(123456)
#Standard deviation of log wage process
siglnw <- 0
#Labor chosen at time timeb
timeb <- 0
#Standard deviation of optimization error
sigoptl <- 0.37
# Number of Firms
n <- 1000
# Number of Time Periods
overallt <- 100
starttime <- 90
t <- overallt - starttime
#Production Function Parameters
alpha0 <- 0
alphal <- 0.6
alphak <- 0.4
#Epsilons and omega ln(wage) process
sigeps <- 0.1
sigomg <- 0.3 #standard deviation of omega
rho <- 0.7 #AR(1) coefficient for omega
sigxi <- sqrt((1-rho^2)*sigomg^2)
rholnw <- 0.3 # AR(1) coefficient for ln(wage)
#Matrices to store data
lnkdata <- matrix(0, n, overallt) #ln(capital)
lnldata <- matrix(0, n, overallt) #ln(labor)
lnmdata <- matrix(0, n, overallt) #ln(intermediate input)
lnwdata <- matrix(0, n, overallt) #ln(wage)
lnpdata <- matrix(0, n, overallt) #ln(output price)
lnydata <- matrix(0, n, overallt) #ln(output)
omgdata <- matrix(0, n, overallt) #omega(t)
omgdataminusb <- matrix(0, n, overallt) #omega(t-b)

#Location Scale Parameters
eta0 <- 1
etak <- 0.7
etal <- -0.6


#Specification for Error Distribution for DGPs
epsdata <- matrix(rnorm(n*overallt, 0, sigeps), nrow=n, ncol=overallt)


#subdividing the AR(1) process
rhofirst <- rho^(1-timeb)
rhosecond <- rho^(timeb)
sigxifirst <- sqrt((1-rhofirst^2)*sigomg^2)
sigxisecond <- sqrt((1-rhosecond^2)*sigomg^2) #Standard deviation of innovation in omega

sigxilnw <- sqrt((1-rholnw^2)*siglnw^2) #standard deviation of innovation in lnw(wage)

#Period 0 values of omega and ln(wage)
omgdata0 <- matrix(rnorm(n,0,sigomg),nrow=n,ncol=1)
lnwdata0 <- matrix(rnorm(n,0,siglnw),nrow=n,ncol=1)

#Period 1-b values of omega and period 1 values of omega and ln(wage)
omgdataminusb[,1] <- rhofirst*omgdata0+matrix(rnorm(n,0,sigxifirst),nrow=n,ncol=1)
omgdata[,1] <- rhosecond*omgdataminusb[,1]+matrix(rnorm(n,0,sigxisecond),nrow=n,ncol=1)
lnwdata[,1] <- rholnw*lnwdata0 + matrix(rnorm(n,0,sigxilnw),nrow=n,ncol=1)


#Simulate values of omega and ln(wage) for rest of time periods
for (s in 2:overallt){
  omgdataminusb[,s] <- rhofirst*omgdata[,s-1] + matrix(rnorm(n,0,sigxifirst),nrow=n,ncol=1)
  omgdata[,s] <- rhosecond*omgdataminusb[,s] + matrix(rnorm(n,0,sigxisecond),nrow=n,ncol=1)
  lnwdata[,s] <- rholnw*lnwdata[,s-1] + matrix(rnorm(n,0,sigxilnw),nrow=n,ncol=1)
}


#Intital ln(capital) level (close to 0)
lnkdata[,1] <- matrix(-100,n,1)
#Discount Rate for DP problem
disc <- 0.95
#Depreciation rate of capital
delta <- 0.2
#Variation in capital adjustment costs across firms
sigb <- 0.6
#See page 36 ACF 1/Phi(i) is distributed lognormally 
#across firms but constant overtime with sd 0.6
oneoverbiadj <- exp(rnorm(n,0,sigb))

#Simplifying components of the optimal investment rule
#Square bracket component
squarebracketterm <- (alphal^(alphal/(1-alphal)))*exp(0.5*alphal^2*sigoptl^2) -
(alphal^(1/(1-alphal)))*exp(0.5*sigoptl^2)
#Constant term in front of sum including squarebracketterm
const1 <- disc*(alphak/(1-alphal))*(exp(alpha0)^(1/(1-alphal)))*
squarebracketterm
vec1 <- (disc*(1-delta))^seq(100)
vec2 <- cumsum(rholnw^(2*seq(100)))
vec3 <- (sigxi^2)*cumsum(rho^(2*seq(0,99,1)))
vec3 <- cbind(sigxi^2 * 0,cumsum(rho^(2*(seq(100)-1))) )
expterm3 <- exp(0.5*((-alphal)/(1-alphal))^2*((sigxilnw^2)*vec2))
expterm4 <- exp(0.5*(1/(1-alphal))^2*rhosecond^2*
((sigxifirst^2)*rho^(2*seq(100))+vec3))
expterm5 <- exp((1/(1-alphal))*(1/2)*sigxisecond^2)
#Compute Optimal Investment and Capital stock for all firms over time
investmat <- matrix(NA, n, overallt)

for (i in 1:n){
    for (s in 1:overallt){
      expterm1 <- exp((1/(1-alphal))*omgdata[i,s]*rho^(seq(100)))
      expterm2 <- exp(((-alphal)/(1-alphal))*lnwdata[i,s]*(rholnw^seq(100)))
      investmat[i,s] <- oneoverbiadj[i]*const1*expterm5*sum(vec1*expterm1*
        expterm2*expterm3*expterm4)

      if (s >= 2){
        lnkdata[i,s] <- log((1-delta)*exp(lnkdata[i,s-1])+
          (1-0*runif(1))*investmat[i,s-1])
      }
    }
  }


#Generate levels of labor input
for (s in 1:overallt){
  lnldata[,s] <- ((sigxisecond^2)/2+log(alphal)+alpha0+
    rhosecond*omgdataminusb[,s]-lnwdata[,s]+lnpdata[,s]+
    (alphak)*lnkdata[,s])/(1-alphal)
}

#Potential Optimization Error
truelnldata <- lnldata
lnldata <- lnldata + matrix(rnorm(n*overallt,0,sigoptl),n,overallt)

#Output and Materials
#Specifies the form of heteroskedasticity
het <- etal*lnldata+etak*lnkdata
lnydata <- alpha0 + alphal*lnldata + alphak*lnkdata + omgdata + het*epsdata
lnmdata <- alpha0 + alphal*truelnldata + alphak*lnkdata + omgdata

#Stack data across firms (all the data)
Capital <- c(t(lnkdata[,(starttime+1):overallt]))
Labor <- c(t(lnldata[,(starttime+1):overallt]))
Materials <- c(t(lnmdata[,(starttime+1):overallt]))
Investment <- c(t(investmat[,(starttime+1):overallt]))
Wage <- c(t(lnwdata[,(starttime+1):overallt]))
Price <- c(t(lnpdata[,(starttime+1):overallt]))
Output <- c(t(lnydata[,(starttime+1):overallt]))
Productivity <- c(t(omgdata[,(starttime+1):overallt]))
Productivity_t_minus_b <- c(t(omgdataminusb[,(starttime+1):overallt]))
True_Labor <- c(t(truelnldata[,(starttime+1):overallt]))
Epsilon <- c(t(epsdata[,(starttime+1):overallt]))

#Stack data across firms (lagged data)
Capital_Lag_1 <- c(t(lnkdata[,(starttime+1):(overallt-1)]))
Labor_Lag_1 <- c(t(lnldata[,(starttime+1):(overallt-1)]))
Labor_Lag_2 <- c(t(lnldata[,(starttime):(overallt-2)]))
Materials_Lag_1 <- c(t(lnmdata[,(starttime+1):(overallt-1)]))
Investment_Lag_1 <- c(t(investmat[,(starttime+1):overallt-1]))
Wage_Lag_1 <- c(t(lnwdata[,(starttime+1):(overallt-1)]))
Price_Lag_1 <- c(t(lnpdata[,(starttime+1):(overallt-1)]))
Output_Lag_1 <- c(t(lnydata[,(starttime+1):(overallt-1)]))
Productivity_Lag_1 <- c(t(omgdata[,(starttime+1):(overallt-1)]))
Productivity_t_minus_b_Lag_1 <- c(t(omgdataminusb[,(starttime+1):(overallt-1)]))
True_Labor_Lag_1 <- c(t(truelnldata[,(starttime+1):(overallt-1)]))

#Stack data across firms (contemporaneous data)
Capital_Con <- c(t(lnkdata[,(starttime+2):(overallt)]))
Labor_Con <- c(t(lnldata[,(starttime+2):(overallt)]))
Materials_Con <- c(t(lnmdata[,(starttime+2):(overallt)]))
Investment_Con <- c(t(investmat[,(starttime+2):overallt]))
Wage_Con <- c(t(lnwdata[,(starttime+2):(overallt)]))
Price_Con <- c(t(lnpdata[,(starttime+2):(overallt)]))
Output_Con <- c(t(lnydata[,(starttime+2):(overallt)]))
Productivity_Con <- c(t(omgdata[,(starttime+2):(overallt)]))
Productivity_t_minus_b_Con <- c(t(omgdataminusb[,(starttime+2):(overallt)]))
True_Labor_Con <- c(t(truelnldata[,(starttime+2):(overallt)]))
Initial_Age <- floor(as.matrix(runif(n, 1, 10)))
Age_Con <- c(sapply(1:n, function(i) seq(Initial_Age[i,], Initial_Age[i,]+8, by=1)))


overall.start.time <- Sys.time()
test <- Production_EM(ntau=ntau, Y=Output_Con, K=Capital_Con, L=Labor_Con, M=Materials_Con, I=Investment_Con, A=Age_Con, N=n, T=t, maxiter=maxiter, draws=draws, Mdraws=Mdraws)
print(Sys.time()-overall.start.time)








