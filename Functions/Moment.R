qmoment <- function(Y, X, c, tau){
	obj <- mean((Y-X%*%c[tau])*(tau-(Y<X%*%c[tau])))
	return(obj)
}
qjacobian <- function(Y, X, c, tau){
	jac <- mean((tau-(Y<X%*%c[tau])))
	return(jac)
}
emoment <- function(Y, X, c){
	obj <- colMeans(X*array(Y-X%*%c), dim=dim(X))
	return(obj)
}
ejacobian <- function(Y, X, c){
	return(colMeans(-X))
}