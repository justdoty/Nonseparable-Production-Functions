# setwd('/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Functions')
# require(dplyr)
source('NLPFQR/FUN/Moments.R')
source('NLPFQR/FUN/Auxfuns.R')
#####################################################################################################
#This file contains the functions which perform all M-step computations
#Requires Moment.R to evaluate objective function and jacobian calculations
#Required Auxfuns.R to draw data for stochastic gradient descent and other overhead tasks
#####################################################################################################

#####################################################################################################
#This is the optimization routine for a gradient descent algorithm for output, labor, materials, and productivity (Quantiles)
#We reccomend a stochastic gradient descent algorithm due to high dimensional calculations of the gradient
#when N, T, M, and number of parameters in the polynomial terms are very high
###################################################################################################
sgradfunq <- function(idvar, X, U, MH, init, tau, tol, stepsize, nsteps, nbatch, WT, seed){
	Bseed <- seed
	data.init <- gradbatch(idvar=idvar, X=X, U=U, nbatch=nbatch, WT=WT, seed=Bseed)
	grad.init <- qjacobian(Y=data.init[, grepl("Y", names(data.init))], X=data.init[, grepl("X", names(data.init))], U=data.init[, grepl("U", names(data.init))], MH=MH, c=init, tau=tau, WT=WT)
	theta <- init-stepsize*grad.init
	loss <- c()
	opttime <- Sys.time()
	for (s in 1:nsteps){
	  Bseed <- Bseed+s
	  grad.data <- gradbatch(idvar=idvar, X=X, U=U, nbatch=nbatch, seed=Bseed, WT=WT)
	  loss <- sqrt(sum(qmoment(Y=grad.data[, grepl("Y", names(grad.data))], X=grad.data[, grepl("X", names(grad.data))], U=grad.data[, grepl("U", names(grad.data))], MH=MH, c=theta, tau=tau, WT=WT)^2))
	  grad <- qjacobian(Y=grad.data[, grepl("Y", names(grad.data))], X=data.init[, grepl("X", names(grad.data))], U=data.init[, grepl("U", names(grad.data))], MH=MH, c=theta, tau=tau, WT=WT)
	  theta <- theta-stepsize*grad
	  if(sqrt(sum(grad^2))<=tol){
	  	  print('Algorithm converged')
	      break
	   }
	}
  print(sprintf("Final gradient norm is %g and time was %g", round(sqrt(sum(grad^2)), digits=abs(log10(tol))), round(Sys.time()-opttime, digits=2)))
  values <- list("coef" = theta, "loss" = loss)
  return(values)
}
#####################################################################################################
#This is the optimization routine for a gradient descent algorithm for investment
#We reccomend a stochastic gradient descent algorithm due to high dimensional calculations of the gradient
#when N, T, M, and number of parameters in the polynomial terms are very high
#####################################################################################################
sgradfune <- function(idvar, X, U, MH, init, tol, stepsize, nsteps, nbatch, seed){
	Bseed <- seed
	data.init <- gradbatch(idvar=idvar, X=X, U=U, nbatch=nbatch, seed=Bseed, WT=2)
	grad.init <- ejacobian(Y=data.init[, grepl("Y", names(data.init))], X=data.init[, grepl("X", names(data.init))], U=data.init[, grepl("U", names(data.init))], MH=MH, c=init)
	theta <- init-stepsize*grad.init
	loss <- c()
	for (s in 1:nsteps){
	  Bseed <- Bseed+s
	  #Draw a new firm-level sample
	  grad.data <- gradbatch(idvar=idvar, X=X, U=U, nbatch=nbatch, seed=Bseed, WT=2)
	  loss <- sqrt(sum(emoment(Y=grad.data[, grepl("Y", names(grad.data))], X=grad.data[, grepl("X", names(grad.data))], U=grad.data[, grepl("U", names(grad.data))], MH=MH, c=theta)^2))
	  grad <- ejacobian(Y=grad.data[, grepl("Y", names(grad.data))], X=data.init[, grepl("X", names(grad.data))], U=data.init[, grepl("U", names(grad.data))], MH=MH, c=theta)
	  theta <- theta-stepsize*grad
	  if(sqrt(sum(grad^2))<=tol){
	  	  print('Algorithm converged')
	      break
	   }
	}
  print(sprintf("Final gradient norm is %g and time was %g", round(sqrt(sum(grad^2)), digits=abs(log10(tol))), round(Sys.time()-opttime, digits=2)))
  values <- list("coef" = theta, "loss" = loss)
  return(values)
}
#####################################################################################################
#Functions to estimate the exponential tail parameters on (0, \tau_{1}) and (\tau_{L}, 1)
#####################################################################################################
#For the translog production function
expby <- function(Y, K, L, M, omega, par1, parL){
	YX <- translog(K=K, L=L, M=M, omega=omega)
	b1 <- -mean((Y-omega)<=YX%*%par1)/mean(((Y-omega)-YX%*%par1)*((Y-omega)<=YX%*%par1))
	bL <- mean((Y-omega)>YX%*%parL)/mean(((Y-omega)-YX%*%parL)*((Y-omega)>YX%*%parL))
	return(list(b1=b1, bL=bL))
}
#For the input functions
expbx <- function(X, K, omega, par1, parL){
	XX <- cbind(1, K, omega, K*omega, K^2, omega^2)
	b1 <- -mean(X<=XX%*%par1)/mean((X-XX%*%par1)*(X<=XX%*%par1))
	bL <- mean(X>XX%*%parL)/mean((X-XX%*%parL)*(X>XX%*%parL))
	return(list(b1=b1, bL=bL))
}
#For productivity t>1
expbwt <- function(omega, omegalag, par1, parL){
	WX <- cbind(1, omegalag, omegalag^2, omegalag^3)
	b1 <- -mean(omega<=WX%*%par1)/mean((omega-WX%*%par1)*(omega<=WX%*%par1))
	bL <- mean(omega>WX%*%parL)/mean((omega-WX%*%parL)*(omega>WX%*%parL))
	return(list(b1=b1, bL=bL))
}







