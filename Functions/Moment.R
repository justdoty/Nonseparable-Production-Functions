setwd('/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Functions')
source('Hermite.R')
#####################################################################################################
#These are the objective function and its jacobian for output, labor, and materials
#####################################################################################################
qmoment <- function(Y, X, U, MH, c, tau){
	dimM <- prod(MH+1)
	obj <- rowSums(apply(U, 2, function(u) (Y-tensor.prod(MH, cbind(X, u))%*%c)*
			(tau-(Y<tensor.prod(MH, cbind(X, u))%*%c))))
	return(obj)
}
qjacobian <- function(Y, X, U, MH, c, tau, N, T){
	dimM <- prod(MH+1)
	NT <- N*T
	jac <- rowSums(apply(U, 2, function(u) -tensor.prod(MH, cbind(X, u))*
		array((tau-(Y<tensor.prod(MH, cbind(X, u))%*%c)), dim=c(NT, dimM))))
	jac <- unname(tapply(jac, (seq_along(jac)-1) %/% NT, sum))
	return(jac)
}
#####################################################################################################
#This is the objective function and its jacobian for productivity for t>1
#####################################################################################################
qmomentWT <- function(X, U, MH, c, tau, N, T){
	dimM <- prod(MH+1)
	XT <- as.matrix(X[-seq(1, N*T, by=T),])
	obj <- rowSums(apply(U, 2, function(u) (as.matrix(u[-seq(1, N*T, by=T)])-tensor.prod(MH, cbind(XT, as.matrix(u[-seq(T, N*T, by=T)])))%*%c)*
	(tau-(as.matrix(u[-seq(1, N*T, by=T)])<tensor.prod(MH, cbind(XT, as.matrix(u[-seq(T, N*T, by=T)])))%*%c))))
	return(obj)
}
qjacobianWT <- function(X, U, MH, c, tau, N, T){
	dimM <- prod(MH+1)
	XT <- as.matrix(X[-seq(1, N*T, by=T),])
	jac <- rowSums(apply(U, 2, function(u) tensor.prod(MH, cbind(XT, as.matrix(u[-seq(T, N*T, by=T)])))*array(tau-(as.matrix(u[-seq(1, N*T, by=T)])<tensor.prod(MH, cbind(XT, as.matrix(u[-seq(T, N*T, by=T)])))%*%c), dim=c(N*(T-1), dimM))))
	jac <- unname(tapply(jac, (seq_along(jac)-1)%/%(N*(T-1)), sum))
	return(jac)
}
#####################################################################################################
#This is the objective function and its jacobian for productivity for t=1
#####################################################################################################
qmomentW0 <- function(X, U, MH, c, tau, N, T){
	dimM <- prod(MH+1)
	X0 <- as.matrix(X[seq(1, N*T, by=T),])
	obj <- rowSums(apply(U, 2, function(u) (as.matrix(u[seq(1, N*T, by=T)])-
		tensor.prod(MH, X0)%*%c)*(tau-(as.matrix(u[seq(1, N*T, by=T)])<
			tensor.prod(MH, X0)%*%c))))
	return(obj)
}
qjacobianW0 <- function(X, U, MH, c, tau, N, T){
	dimM <- prod(MH+1)
	X0 <- as.matrix(X[seq(1, N*T, by=T),])
	jac <- rowSums(apply(U, 2, function(u) -tensor.prod(MH, X0)*
		array((tau-(as.matrix(u[seq(1, N*T, by=T)])<tensor.prod(MH, X0)%*%c)), dim=c(N, dimM))))
	jac <- unname(tapply(jac, (seq_along(jac)-1) %/% N, sum))
	return(jac)
}
###################################################################################################
#This is the objective function and its jacobian for the investment decision rule
###################################################################################################
emoment <- function(Y, X, U, MH, c){
	obj <- rowSums(apply(U, 2, function(u) (Y-tensor.prod(MH, cbind(X, u))%*%c)))
	return(obj)
}
ejacobian <- function(Y, X, U, MH, c, N, T){
	dimM <- prod(MH+1)
	NT <- N*T
	jac <- rowSums(apply(U, 2, function(u) -tensor.prod(MH, cbind(X, u))*
		array((2*(Y-tensor.prod(MH, cbind(X, u))%*%c)), dim=c(NT, dimM))))
	jac <- unname(tapply(jac, (seq_along(jac)-1) %/% NT, sum))
	return(jac)
}
#####################################################################################################
#This is the optimization routine for a gradient descent algorithm for output, labor, materials, and investment
#####################################################################################################
gradfun <- function(Y, X, U, MH, init, tau, tol, stepsize, nsteps, N, T){
	if (!is.null(tau)){
	  grad.init <- qjacobian(Y=Y, X=X, U=U, MH=MH, c=init, tau=tau, N=N, T=T)
	  } else {
	  grad.init <- ejacobian(Y=Y, X=X, U=U, MH=MH, c=init, N=N, T=T)
	  }
	theta <- init-stepsize*grad.init
	loss <- c()
	for (i in 1:nsteps){
		print(i)
	  if (!is.null(tau)){
	  	loss <- sqrt(sum(qmoment(Y=Y, X=X, U=X, MH=MH, c=theta, tau=tau)^2))
	  	grad <- qjacobian(Y=Y, X=X, U=U, MH=MH, c=theta, tau=tau, N=N, T=T)
	  } else {
	  	loss <- sqrt(sum(emoment(Y=Y, X=X, U=X, MH=MH, c=theta)^2))
	  	grad <- ejacobian(Y=Y, X=X, U=U, MH=MH, c=theta, N=N, T=T)
	  }
	  theta <- theta-stepsize*grad
	    if(sqrt(sum(grad^2))<=tol){
	      break
	    }
	}
  print("Algorithm converged")
  print(paste("Final gradient norm is",sqrt(sum(grad^2))))
  values <- list("coef" = theta, "loss" = loss)
  return(values)
}
#####################################################################################################
#This is the optimization routine for a gradient descent algorithm for productivity at t>1
#####################################################################################################
gradfunWT <- function(X, U, MH, init, tau, tol, stepsize, nsteps, N, T){
	grad.init <- qjacobianWT(X=X, U=U, MH=MH, c=init, tau=tau, N=N, T=T)
	theta <- init-stepsize*grad.init
	loss <- c()
	for (i in 1:nsteps){
		print(i)
	  	loss <- sqrt(sum(qmomentWT(X=X, U=X, MH=MH, c=theta, tau=tau, N=N, T=T)^2))
	  	grad <- qjacobianWT(X=X, U=U, MH=MH, c=theta, tau=tau, N=N, T=T)
		theta <- theta-stepsize*grad
	    if(sqrt(sum(grad^2))<=tol){
	      break
	    }
	}
  print("Algorithm converged")
  print(paste("Final gradient norm is",sqrt(sum(grad^2))))
  values <- list("coef" = theta, "loss" = loss)
  return(values)
}
#####################################################################################################
#This is the optimization routine for a gradient descent algorithm for productivity at t=1
#####################################################################################################
gradfunW0 <- function(X, U, MH, init, tau, tol, stepsize, nsteps, N, T){
	grad.init <- qjacobianW0(X=X, U=U, MH=MH, c=init, tau=tau, N=N, T=T)
	theta <- init-stepsize*grad.init
	loss <- c()
	for (i in 1:nsteps){
		print(i)
	  	loss <- sqrt(sum(qmomentW0(X=X, U=X, MH=MH, c=theta, tau=tau, N=N, T=T)^2))
	  	grad <- qjacobianW0(X=X, U=U, MH=MH, c=theta, tau=tau, N=N, T=T)
		theta <- theta-stepsize*grad
	    if(sqrt(sum(grad^2))<=tol){
	      break
	    }
	}
  print("Algorithm converged")
  print(paste("Final gradient norm is",sqrt(sum(grad^2))))
  values <- list("coef" = theta, "loss" = loss)
  return(values)
}
#####################################################################################################
#####################################################################################################
#####################################################################################################

















