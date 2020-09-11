setwd('/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Functions')
source('Hermite.R')
posterior <- function(matdraw, reslist){
	YX <- tensor.prod(c(MK, ML, MM, MW, MA), cbind(K, L, M, matdraw, A))
	LX <- tensor.prod(c(MK, MW, MA), cbind(K, matdraw, A))
	MX <- tensor.prod(c(MK, MW, MA), cbind(K, matdraw, A))
	IX <- tensor.prod(c(MK, MW, MA), cbind(K, matdraw, A))
	WXT  <- tensor.prod(c(MW, MA), cbind(K, matdraw[-seq(T, N*T, by=N)], A[-seq(1, N*T, by=N)]))
	WX0  <- tensor.prod(MA, A[seq(1:N*T, by=N)])
	#Rename reslist list objects
	reslist$resYinit <- resY; reslist$yb1 <- yb1; reslist$ybL <- ybL; 
	reslist$resLinit <- resL; reslist$lb1 <- lb1; reslist$lbL <- lbL;
	reslist$resMinit <- resM; reslist$mb1 <- mb1; reslist$mbL <- mbL;
	reslist$resWTinit <- resWT; reslist$wtb1 <- wtb1; reslist$wtbL <- wtbL;
	reslist$resW0init <- resW0; reslist$w0b1 <- w0b1; reslist$w0bL <- w0bL;
	reslist$resIinit <- resI
	#Likelihood of Output Data
	densY <- sum(sapply(1:ntau-1), function(tau) ((vectau[lead(tau)]-vectau[tau])/YX%*%(resY[lead(tau)]-resY[tau]))*
					(YX%*%resY[tau]<Y)*(YX%*%resY[tau]>=Y))+vectau[1]*yb1*exp(yb1*(Y-YX%*%resY[1]))*(Y<=YX%*%resY[1])+
					(1-vectau[ntau])*ybL*exp(ybL*(Y-YX%*%resY[ntau]))*(Y>YX%*%resY[ntau])
	densY <- apply(t(matrix(densY, T, N)), 1, prod)
	#Likelihood of Labor Data
	densL <- sum(sapply(1:ntau-1), function(tau) ((vectau[lead(tau)]-vectau[tau])/LX%*%(resL[lead(tau)]-resL[tau]))*
					(LX%*%resL[tau]<L)*(LX%*%resL[tau]>=L))+vectau[1]*lb1*exp(lb1*(L-LX%*%resL[1]))*(L<=LX%*%resL[1])+
					(1-vectau[ntau])*lbL*exp(lbL*(L-LX%*%resL[ntau]))*(L>LX%*%resL[ntau])
	densL <- apply(t(matrix(densL, T, N)), 1, prod)
	#Likelihood of Materials Data
	densM <- sum(sapply(1:ntau-1), function(tau) ((vectau[lead(tau)]-vectau[tau])/MX%*%(resM[lead(tau)]-resM[tau]))*
					(MX%*%resM[tau]<M)*(MX%*%resM[tau]>=M))+vectau[1]*mb1*exp(mb1*(M-MX%*%resM[1]))*(M<=MX%*%resM[1])+
					(1-vectau[ntau])*mbL*exp(mbL*(M-MX%*%resM[ntau]))*(M>MX%*%resM[ntau])
	densM <- apply(t(matrix(densM, T, N)), 1, prod)
	#Likelihood of Investment Data (Log-Normal)
	Res.IX <- I-IX%*%resI[-length(resI)]
	densI <- 1/sqrt(resI[length(resI)])*dnorm(Res.IX/sqrt(resI[length(resI)]))
	densI <- apply(t(matrix(densI, T, N)), 1, prod)
	#Prior Omega_{t} Data
	matdraw.T <- matdraw[-seq(T, N*T, by=N)]
	A.T <- A[-seq(1, N*T, by=N)]
	densWT <- sum(sapply(1:ntau-1), function(tau) ((vectau[lead(tau)]-vectau[tau])/WXT%*%(resWT[lead(tau)]-resWT[tau]))*
					(WXT%*%resWT[tau]<matdraw.T)*(WXT%*%resWT[tau]>=matdraw.T))+vectau[1]*wtb1*exp(wtb1*(matdraw.T-WXT%*%resWT[1]))*(matdraw.T<=WXT%*%resWT[1])+
					(1-vectau[ntau])*wtbL*exp(wtbL*(matdraw.T-WXT%*%resWT[ntau]))*(matdraw.T>WXT%*%resWT[ntau])
	densWT <- apply(t(matrix(densWT, T-1, N)), 1, prod)
	#Prior Omega_{0} Data
	matdraw0 <- matdraw[seq(1:N*T, by=N)]
	densW0 <- sum(sapply(1:ntau-1), function(tau) ((vectau[lead(tau)]-vectau[tau])/WX0%*%(resW0[lead(tau)]-resW0[tau]))*
					(WX0%*%resW0[tau]<matdraw0)*(WX0%*%resW0[tau]>=matdraw0))+vectau[1]*w0b1*exp(w0b1*(matdraw0-WX0%*%resW0[1]))*(matdraw0<=WX0%*%resW0[1])+
					(1-vectau[ntau])*w0bL*exp(w0bL*(matdraw0-WX0%*%resW0[ntau]))*(matdraw0>WX0%*%resW0[ntau])
	#Final posterior density
	finaldens <- densY*densL*densM*densI*densWT*densW0
	return(finaldes)

}