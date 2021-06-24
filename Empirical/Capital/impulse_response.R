require(cowplot)
require(dplyr)
require(purrr)
require(xtable)
require(stringr)
require(truncnorm)
require(RColorBrewer)
require(reshape2)
source('/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Functions/Auxfuns.R')
source('/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Functions/Tensors.R')
source('/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Functions/omega.R')
#Load US Dataset
US <- read.csv('/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Data/USdata.csv') %>% 
select(id, year, Y, K, K2, L, M, I, age) %>% transmute(id=id, year=year, Y=log(Y), K=log(K), L=log(L), M=log(M), I=log(I), A=age)
#Productivity from LP
omegainit <- omega_est(idvar=US$id, timevar=US$year, Y=US$Y, A=US$A, K=US$K, L=US$L, M=US$M)$omega
########################################################################################################
##########################################Load Results############################################
########################################################################################################
set.seed(123456)
load("/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Empirical/Capital/Capital_trans.RData")
vectau <- results$vectau
ntau <- length(vectau)
dims <- results$dims
#Load Parameter Estimates
parY <- results$resYmat
parL <- results$resLmat
parM <- results$resMmat
parKT <- results$resKTmat
parWT <- results$resWTmat
parW1 <- results$resW1mat
parK1 <- results$resK1mat
parYb <- results$resyb1bLmat
parLb <- results$reslb1bLmat
parMb <- results$resmb1bLmat
parWTb <- results$reswtb1bLmat
parW1b <- results$resw1tb1bLmat
WTminmax <- results$maxminwtmat
#Load the method used (e.g. Cobb, Translog, Hermite)
method <- results$method
resY <- results$resY
###############################################################################
#Simulate Productivity and Capital Evolution Given Model Parameters
##############################################################################
#We simuluate a balanced panel even though the original model is estimated from
#an unbalanced panel. In a later version, we consider adding a selection bias
#correction to the productivity equation and drop firms according to this rule
#in the simulated model
#############################################################################
N <- 10000
T <- length(unique(US$year))
#Age##################################################################
adata <- matrix(0, N, T)
tauinp <- c(0.1, 0.5, 0.9)
#Labor###############################################################
lnldata <- matrix(0, N, T)
labpath <- matrix(0, T, length(tauinp))
#High Labor Path
highlab <- labpath
#Median Labor Path
medlab <- labpath
#Low Labor Path
lowlab <- labpath 
#Unobservable Shock to Labor
epsdata <- matrix(runif(N*T), nrow=N, ncol=T)
#Intermediate Input#####################################################
lnmdata <- matrix(0, N, T)
matpath <- matrix(0, T, length(tauinp))
#High Materials Path
highmat <- matpath
#Median Labor Path
medmat <- matpath
#Low Labor Path
lowmat <- matpath
#Unobservable Shock to Intermediate Inputs
varepsdata <- matrix(runif(N*T), nrow=N, ncol=T)
#Capital##################################################################
lnkdata <- matrix(0, N, T)
#High Capital Path
highcap <- lnkdata
#Median Labor Path
medcap <- lnkdata
#Low Labor Path
lowcap <- lnkdata
#Productivity#############################################################
omgdata <- matrix(0, N, T)
#High Productivity Path
highomg <- omgdata
#Median Productivity Path
medomg <- omgdata
#Low Productivity Path
lowomg <- omgdata
#Unobservable Shock to Productivity
xidata <- matrix(runif(N*T), nrow=N, ncol=T)
#########################################################################################################
#Data at T=1
t1data <- US %>% group_by(id) %>% slice(1)
#Simulate Initial Age-Distribution from Log-Normal
A1 <- log(t1data$A)
A <- exp(rtruncnorm(n=N, a=min(A1), b=max(A1), mean=mean(A1), sd=sd(A1)))
adata[,1] <- A
#Initial Productivity (Give each firm the same rank initial  productivity specified by tauinit)
tauinit <- 0.5
xidata[,1] <- tauinit
omgdata[,1] <- rowSums(W1X(A=A)*lspline(vectau=vectau, bvec=parW1, b1=parW1b[1], bL=parW1b[2], u=xidata[,1]))
highomg[,1] <- omgdata[,1]
medomg[,1] <- omgdata[,1]
lowomg[,1] <- omgdata[,1]
#Initial Capital
lnkdata[,1] <- rnorm(N, mean=K1X(A=A, omega=omgdata[,1])%*%parK1[1:(length(parK1)-1)], sd=sqrt(parK1[length(parK1)]))
#Restrict Support of Initial Capital
lnkdata[,1] <- (lnkdata[,1]>max(t1data$K))*max(t1data$K)+(lnkdata[,1]<min(t1data$K))*min(t1data$K)+(lnkdata[,1]<=max(t1data$K))*(lnkdata[,1]>=min(t1data$K))*lnkdata[,1]
#Paths for Initial Capital
highcap[,1] <- lnkdata[,1]
medcap[,1] <- lnkdata[,1]
lowcap[,1] <- lnkdata[,1]
#Evolution of Productivity and Capital
for (t in 2:T){
	ttdata <- US %>% group_by(id) %>% slice(t)
	#Age
	A <- A+1
	adata[,t] <- A
	#Capital
	highcap[,t] <- rnorm(N, mean=KX(A=A, K=highcap[,t-1], omega=highomg[,t-1])%*%parKT[1:(length(parKT)-1)], sd=sqrt(parKT[length(parKT)]))
	medcap[,t] <- rnorm(N, mean=KX(A=A, K=medcap[,t-1], omega=medomg[,t-1])%*%parKT[1:(length(parKT)-1)], sd=sqrt(parKT[length(parKT)]))
	lowcap[,t] <- rnorm(N, mean=KX(A=A, K=lowcap[,t-1], omega=lowomg[,t-1])%*%parKT[1:(length(parKT)-1)], sd=sqrt(parKT[length(parKT)]))
	#Restrict the Support of Capital 
	highcap[,t] <- (highcap[,t]>max(ttdata$K))*max(ttdata$K)+(highcap[,t]<min(ttdata$K))*min(ttdata$K)+(highcap[,t]<=max(ttdata$K))*(highcap[,t]>=min(ttdata$K))*highcap[,t]
	medcap[,t] <- (medcap[,t]>max(ttdata$K))*max(ttdata$K)+(medcap[,t]<min(ttdata$K))*min(ttdata$K)+(medcap[,t]<=max(ttdata$K))*(medcap[,t]>=min(ttdata$K))*medcap[,t]
	lowcap[,t] <- (lowcap[,t]>max(ttdata$K))*max(ttdata$K)+(lowcap[,t]<min(ttdata$K))*min(ttdata$K)+(lowcap[,t]<=max(ttdata$K))*(lowcap[,t]>=min(ttdata$K))*lowcap[,t]
	#Productivity
	#At t=2, simulate productivity paths when hit by different sized innovation shocks
	if (t==2){
		#High Innovation Shock
		highxi <- 0.9
		xidata[,t] <- highxi
		highomg[,t] <- rowSums(WX(A=A, omega=highomg[,t-1])*lspline(vectau=vectau, bvec=parWT, b1=parWTb[1], bL=parWTb[2], u=xidata[,t]))
		#Restrict the Support of High Productivity
		highomg[,t] <- (highomg[,t]>WTminmax[1])*WTminmax[1]+(highomg[,t]<WTminmax[2])*WTminmax[2]+(highomg[,t]<=WTminmax[1])*(highomg[,t]>=WTminmax[2])*highomg[,t]
		#Median Innovation Shock
		medxi <- 0.5
		xidata[,t] <- medxi
		medomg[,t] <- rowSums(WX(A=A, omega=medomg[,t-1])*lspline(vectau=vectau, bvec=parWT, b1=parWTb[1], bL=parWTb[2], u=xidata[,t]))
		#Restrict the Support of Median Productivity
		medomg[,t] <- (medomg[,t]>WTminmax[1])*WTminmax[1]+(medomg[,t]<WTminmax[2])*WTminmax[2]+(medomg[,t]<=WTminmax[1])*(medomg[,t]>=WTminmax[2])*medomg[,t]
		#Low Innovation Shock
		lowxi <- 0.1
		xidata[,t] <- lowxi
		lowomg[,t] <- rowSums(WX(A=A, omega=lowomg[,t-1])*lspline(vectau=vectau, bvec=parWT, b1=parWTb[1], bL=parWTb[2], u=xidata[,t]))
		#Restrict the Support of Low Productivity
		lowomg[,t] <- (lowomg[,t]>WTminmax[1])*WTminmax[1]+(lowomg[,t]<WTminmax[2])*WTminmax[2]+(lowomg[,t]<=WTminmax[1])*(lowomg[,t]>=WTminmax[2])*lowomg[,t]
	} else {
		#High Productivity
		highomg[,t] <- rowSums(WX(A=A, omega=highomg[,t-1])*lspline(vectau=vectau, bvec=parWT, b1=parWTb[1], bL=parWTb[2], u=xidata[,t]))
		#Restrict the Support of High Productivity
		highomg[,t] <- (highomg[,t]>WTminmax[1])*WTminmax[1]+(highomg[,t]<WTminmax[2])*WTminmax[2]+(highomg[,t]<=WTminmax[1])*(highomg[,t]>=WTminmax[2])*highomg[,t]
		#Median Productivity
		medomg[,t] <- rowSums(WX(A=A, omega=medomg[,t-1])*lspline(vectau=vectau, bvec=parWT, b1=parWTb[1], bL=parWTb[2], u=xidata[,t]))
		#Restrict the Support of Median Productivity
		medomg[,t] <- (medomg[,t]>WTminmax[1])*WTminmax[1]+(medomg[,t]<WTminmax[2])*WTminmax[2]+(medomg[,t]<=WTminmax[1])*(medomg[,t]>=WTminmax[2])*medomg[,t]
		#Low Productivity
		lowomg[,t] <- rowSums(WX(A=A, omega=lowomg[,t-1])*lspline(vectau=vectau, bvec=parWT, b1=parWTb[1], bL=parWTb[2], u=xidata[,t]))
		#Restrict the Support of Low Productivity
		lowomg[,t] <- (lowomg[,t]>WTminmax[1])*WTminmax[1]+(lowomg[,t]<WTminmax[2])*WTminmax[2]+(lowomg[,t]<=WTminmax[1])*(lowomg[,t]>=WTminmax[2])*lowomg[,t]
	}
}
for (t in 1:T){
	ttdata <- US %>% group_by(id) %>% slice(t)
	for (q in 1:length(tauinp)){
		#Labor Paths
		epsdata[,t] <- tauinp[q]
		lab1 <- rowSums(LX(A=adata[,t], K=highcap[,t], omega=highomg[,t])*lspline(vectau=vectau, bvec=parL, b1=parLb[1], bL=parLb[2], u=epsdata[,t]))
		lab2 <- rowSums(LX(A=adata[,t], K=medcap[,t], omega=medomg[,t])*lspline(vectau=vectau, bvec=parL, b1=parLb[1], bL=parLb[2], u=epsdata[,t]))
		lab3 <- rowSums(LX(A=adata[,t], K=lowcap[,t], omega=lowomg[,t])*lspline(vectau=vectau, bvec=parL, b1=parLb[1], bL=parLb[2], u=epsdata[,t]))
		#Restrict the Supports
		highlab[t,q] <- mean((lab1>max(ttdata$L))*max(ttdata$L)+(lab1<min(ttdata$L))*min(ttdata$L)+(lab1<=max(ttdata$L))*(lab1>=min(ttdata$L))*lab1)
		medlab[t,q] <- mean((lab2>max(ttdata$L))*max(ttdata$L)+(lab2<min(ttdata$L))*min(ttdata$L)+(lab2<=max(ttdata$L))*(lab2>=min(ttdata$L))*lab2)
		lowlab[t,q] <- mean((lab3>max(ttdata$L))*max(ttdata$L)+(lab3<min(ttdata$L))*min(ttdata$L)+(lab3<=max(ttdata$L))*(lab3>=min(ttdata$L))*lab3)
		#Material Paths
		varepsdata[,t] <- tauinp[q]
		mat1 <- rowSums(MX(A=adata[,t], K=highcap[,t], omega=highomg[,t])*lspline(vectau=vectau, bvec=parM, b1=parMb[1], bL=parMb[2], u=varepsdata[,t]))
		mat2 <- rowSums(MX(A=adata[,t], K=medcap[,t], omega=medomg[,t])*lspline(vectau=vectau, bvec=parM, b1=parMb[1], bL=parMb[2], u=varepsdata[,t]))
		mat3 <- rowSums(MX(A=adata[,t], K=lowcap[,t], omega=lowomg[,t])*lspline(vectau=vectau, bvec=parM, b1=parMb[1], bL=parMb[2], u=varepsdata[,t]))
		#Restrict the Supports
		highmat[t,q] <- mean((mat1>max(ttdata$M))*max(ttdata$M)+(mat1<min(ttdata$M))*min(ttdata$M)+(mat1<=max(ttdata$M))*(mat1>=min(ttdata$M))*mat1)
		medmat[t,q] <- mean((mat2>max(ttdata$M))*max(ttdata$M)+(mat2<min(ttdata$M))*min(ttdata$M)+(mat2<=max(ttdata$M))*(mat2>=min(ttdata$M))*mat2)
		lowmat[t,q] <- mean((mat3>max(ttdata$M))*max(ttdata$M)+(mat3<min(ttdata$M))*min(ttdata$M)+(mat3<=max(ttdata$M))*(mat3>=min(ttdata$M))*mat3)
	}
}
#Paths for Labor
highlabpath <- highlab-medlab
lowlabpath <- lowlab-medlab
#Paths for Materials
highmatpath <- highmat-medmat
lowmatpath <- lowmat-medmat

plot(1:T, lowmatpath[,1], type="l")








