############################################################################################################
#Function that defines the (log) posterior density for which to evaluate the Metropolis Hastings algorithm (target distribution)
############################################################################################################
#Faster function for ifELSE calculations
ifELSE <- function(test, yes, no) {
  out <- rep(NA, length(test))
  out[test] <- yes[test]
  out[!test] <- no[!test]
  return(out)
}
# Density Function for the Conditional Quantiles
densq <- function(idvar, YX, vectau, pred, parb1, parbL){
	densq <- suppressWarnings(ifELSE((pred[,1]<YX)&(pred[,2]>=YX), log(vectau[2]-vectau[1])-log(pred[,2]-pred[,1]), 0))
	for (q in 2:(length(vectau)-1)){ 
		densq <- suppressWarnings(ifELSE((pred[,q]<YX)&(pred[,q+1]>=YX), log(vectau[q+1]-vectau[q])-log(pred[,q+1]-pred[,q]), densq))
	}
	densq <- suppressWarnings(ifELSE(YX<=pred[,1], log(vectau[1])+log(parb1)+parb1*(YX-pred[,1]), densq))
	densq <- suppressWarnings(ifELSE(YX>pred[,length(vectau)], log(1-vectau[length(vectau)])+log(parbL)-parbL*(YX-pred[,length(vectau)]), densq))
	densq <- data.table(id=idvar, densq=densq)[,sum(densq),keyby=id]$V1
	return(densq)
	}
#Distribution for Productivity at t=1
densN1 <- function(YX, pred, var){
	resx <- YX-pred
	#Log-Normal
	densq <- -0.5*log(2*pi)-0.5*log(var)-(0.5/var)*resx^2 
	return(densq)
}
###############################################################################################################
posterior <- function(ydata, inpdata, wtdata, rdata, w1data, vectau, par, method){
	#Final posterior density########################################################################
	#For R&D Firms
	densY <- densq(idvar=ydata$idvar, YX=ydata$Y, vectau=vectau, pred=predict(par$resYinit, ydata), parb1=(par$ybinit)$b1, parbL=(par$ybinit)$bL)
	densL <- densq(idvar=inpdata$idvar, YX=inpdata$LX, vectau=vectau, pred=predict(par$resLinit, inpdata), parb1=(par$lbinit)$b1, parbL=(par$lbinit)$bL)
	densM <- densq(idvar=inpdata$idvar, YX=inpdata$MX, vectau=vectau, pred=predict(par$resMinit, inpdata), parb1=(par$mbinit)$b1, parbL=(par$mbinit)$bL)
	densI <- densq(idvar=inpdata$idvar, YX=inpdata$IX, vectau=vectau, pred=predict(par$resIinit, inpdata), parb1=(par$ibinit)$b1, parbL=(par$ibinit)$bL)
	Rind <- unique(ydata$idvar)%in%unique(ydata$idvar[ydata$R>0])
	densR <- densI
	densR[Rind] <- densq(idvar=rdata$idvar, YX=rdata$RX, vectau=vectau, pred=predict(par$resRinit, rdata), parb1=(par$rbinit)$b1, parbL=(par$rbinit)$bL)
	densR[!Rind] <- 0
	densWT <- densq(idvar=wtdata$idvar, YX=wtdata$Wcon, vectau=vectau, pred=predict(par$resWTinit, wtdata), parb1=(par$wtbinit)$b1, parbL=(par$wtbinit)$bL)
	#Initial Productivity
	densW1 <- densq(idvar=w1data$idvar, YX=w1data$W1, vectau=vectau, pred=predict(par$resW1init, w1data), parb1=(par$w1binit)$b1, parbL=(par$w1binit)$bL)
	#Final posterior
	finaldens <- densY+densL+densM+densI+densR+densWT+densW1
}
######################################################################################################
######################################################################################################
######################################################################################################

