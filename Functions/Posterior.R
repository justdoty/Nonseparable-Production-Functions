setwd('/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Functions')
source('Hermite.R')
posterior <- function(Y, K, L, M, A, I, MH, ntau, vectau, N, T, matdraw, resinit){
	#Rename resinit list objects
	MK <- MH[1]; ML <- MH[2]; MM <- MH[3]; MA <- MH[4]; MW <- MH[5]; MI <- MH[6]
	resYinit <- resinit$resYinit; yb1init <- resinit$yb1init; ybLinit <- resinit$ybLinit; 
	resLinit <- resinit$resLinit; lb1init <- resinit$lb1init; lbLinit <- resinit$lbLinit;
	resMinit <- resinit$resMinit; mb1init <- resinit$mb1init; mbLinit <- resinit$mbLinit;
	resWTinit <- resinit$resWTinit; wtb1init <- resinit$wtb1init; wtbLinit <- resinit$wtbLinit;
	resW0init <- resinit$resW0init; w0b1init <- resinit$w0b1init; w0bLinit <- resinit$w0bLinit;
	resIinit <- resinit$resIinit
	#Create Hermite Tensor Products
	YX <- tensor.prod(c(MK, ML, MM, MW, MA), cbind(K, L, M, matdraw, A))
	LX <- tensor.prod(c(MK, MW, MA), cbind(K, matdraw, A))
	MX <- tensor.prod(c(MK, MW, MA), cbind(K, matdraw, A))
	IX <- tensor.prod(c(MK, MW, MA), cbind(K, matdraw, A))
	WXT  <- tensor.prod(c(MW, MA), cbind(matdraw[-seq(T, N*T, by=T)], A[-seq(1, N*T, by=T)]))
	WX0  <- tensor.prod(MA, as.matrix(A[seq(1, N*T, by=T)]))
	#Likelihood of Output Data
	densY <- as.matrix(rowSums(sapply(1:(ntau-1), function(tau) ((vectau[tau+1]-vectau[tau])/(YX%*%(resYinit[,tau+1]-resYinit[,tau])))*
							(YX%*%resYinit[,tau]<Y)*(YX%*%resYinit[,tau+1]>=Y))))+vectau[1]*yb1init*exp(yb1init*(Y-YX%*%resYinit[,1]))*(Y<=YX%*%resYinit[,1])+
					(1-vectau[ntau])*ybLinit*exp(-ybLinit*(Y-YX%*%resYinit[,ntau]))*(Y>YX%*%resYinit[,ntau])
	densY <- as.matrix(apply(t(matrix(log(densY), T, N)), 1, sum))
	# # #Likelihood of Labor Data
	densL <- as.matrix(rowSums(sapply(1:(ntau-1), function(tau) ((vectau[tau+1]-vectau[tau])/(LX%*%(resLinit[,tau+1]-resLinit[,tau])))*
							(LX%*%resLinit[,tau]<L)*(LX%*%resLinit[,tau+1]>=L))))+vectau[1]*lb1init*exp(lb1init*(L-LX%*%resLinit[,1]))*(L<=LX%*%resLinit[,1])+
					(1-vectau[ntau])*lbLinit*exp(-lbLinit*(L-LX%*%resLinit[,ntau]))*(L>LX%*%resLinit[,ntau])
	densL <- as.matrix(apply(t(matrix(log(densL), T, N)), 1, sum))
	# # # #Likelihood of Materials Data
	densM <- as.matrix(rowSums(sapply(1:(ntau-1), function(tau) ((vectau[tau+1]-vectau[tau])/(MX%*%(resMinit[,tau+1]-resMinit[,tau])))*
							(MX%*%resMinit[,tau]<M)*(MX%*%resMinit[,tau+1]>=M))))+vectau[1]*mb1init*exp(mb1init*(M-MX%*%resMinit[,1]))*(M<=MX%*%resMinit[,1])+
					(1-vectau[ntau])*mbLinit*exp(-mbLinit*(M-MX%*%resMinit[,ntau]))*(M>MX%*%resMinit[,ntau])
	densM <- as.matrix(apply(t(matrix(log(densM), T, N)), 1, sum))
	# # # #Likelihood of Investment Data (Log-Normal)
	Res.IX <- I-IX%*%resIinit[-length(resIinit)]
	densI <- 1/sqrt(resIinit[length(resIinit)])*dnorm(Res.IX/sqrt(resIinit[length(resIinit)]))
	densI <- as.matrix(apply(t(matrix(log(densI), T, N)), 1, sum))
	# # # # #Prior Omega_{t} Data
	matdraw.T <- matdraw[-seq(T, N*T, by=T)]
	A.T <- A[-seq(1, N*T, by=T)]
	densWT <- as.matrix(rowSums(sapply(1:(ntau-1), function(tau) ((vectau[tau+1]-vectau[tau])/(WXT%*%(resWTinit[,tau+1]-resWTinit[,tau])))*
							(WXT%*%resWTinit[,tau]<matdraw.T)*(WXT%*%resWTinit[,tau+1]>=matdraw.T))))+vectau[1]*wtb1init*exp(wtb1init*(matdraw.T-WXT%*%resWTinit[,1]))*(matdraw.T<=WXT%*%resWTinit[,1])+
					(1-vectau[ntau])*wtbLinit*exp(-wtbLinit*(matdraw.T-WXT%*%resWTinit[,ntau]))*(matdraw.T>WXT%*%resWTinit[,ntau])
	densWT <- as.matrix(apply(t(matrix(log(densWT), (T-1), N)), 1, sum))
	# # # # #Prior Omega_{0} Data
	matdraw0 <- as.matrix(matdraw[seq(1, N*T, by=T)])
	densW0 <- as.matrix(rowSums(sapply(1:(ntau-1), function(tau) ((vectau[tau+1]-vectau[tau])/(WX0%*%(resW0init[,tau+1]-resW0init[,tau])))*
							(WX0%*%resW0init[,tau]<matdraw0)*(WX0%*%resW0init[,tau+1]>=matdraw0))))+vectau[1]*w0b1init*exp(w0b1init*(matdraw0-WX0%*%resW0init[,1]))*(matdraw0<=WX0%*%resW0init[,1])+
					(1-vectau[ntau])*w0bLinit*exp(-w0bLinit*(matdraw0-WX0%*%resW0init[,ntau]))*(matdraw0>WX0%*%resW0init[,ntau])
	densW0 <- log(densW0)
	# # # # #Final posterior density
	finaldens <- densY+densL+densM+densWT+densI+densW0
	return(finaldens)

}