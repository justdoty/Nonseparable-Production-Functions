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
#Vector or rank of input demand functions (small, medium, large shocks)
tauinp <- c(0.1, 0.5, 0.9)
#Rank of initial productivity
tauinit <- 0.5
#Age##################################################################
adata <- matrix(0, N, T)
#Labor###############################################################
lnldata <- matrix(0, N, T)
# labpath <- matrix(0, T, length(tauinp))
#High Labor Path
highlab <- lnldata
#Median Labor Path
medlab <- lnldata
#Low Labor Path
lowlab <- lnldata
#Unobservable Shock to Labor
epsdata <- matrix(runif(N*T), nrow=N, ncol=T)
#Intermediate Input#####################################################
lnmdata <- matrix(0, N, T)
# matpath <- matrix(0, T, length(tauinp))
#High Materials Path
highmat <- lnmdata
#Median Labor Path
medmat <- lnmdata
#Low Labor Path
lowmat <- lnmdata
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
highxi <- xidata; medxi <- xidata; lowxi <- xidata
#########################################################################################################
#Data at T=1
t1data <- US %>% group_by(id) %>% slice(1)
#Initial Age
A1 <- t1data$A
#Initial Productivity (Give each firm the same rank initial  productivity specified by tauinit)
omgdata[,1] <- rowSums(W1X(A=adata[,1])*lspline(vectau=vectau, bvec=parW1, b1=parW1b[1], bL=parW1b[2], u=xidata[,1]))
omgdata[,1] <- quantile(omgdata[,1], tauinit)
#All Productivity Paths start at this inital value
highomg[,1] <- omgdata[,1]
medomg[,1] <- omgdata[,1]
lowomg[,1] <- omgdata[,1]
#Initial Capital
lnkdata[,1] <- rnorm(N, mean=K1X(A=adata[,1], omega=omgdata[,1])%*%parK1[1:(length(parK1)-1)], sd=sqrt(parK1[length(parK1)]))
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
	adata[,t] <- adata[,t-1]+1
	#Capital
	highcap[,t] <- rnorm(N, mean=KX(A=adata[,t], K=highcap[,t-1], omega=highomg[,t-1])%*%parKT[1:(length(parKT)-1)], sd=sqrt(parKT[length(parKT)]))
	medcap[,t] <- rnorm(N, mean=KX(A=adata[,t], K=medcap[,t-1], omega=medomg[,t-1])%*%parKT[1:(length(parKT)-1)], sd=sqrt(parKT[length(parKT)]))
	lowcap[,t] <- rnorm(N, mean=KX(A=adata[,t], K=lowcap[,t-1], omega=lowomg[,t-1])%*%parKT[1:(length(parKT)-1)], sd=sqrt(parKT[length(parKT)]))
	#Restrict the Support of Capital 
	highcap[,t] <- (highcap[,t]>max(ttdata$K))*max(ttdata$K)+(highcap[,t]<min(ttdata$K))*min(ttdata$K)+(highcap[,t]<=max(ttdata$K))*(highcap[,t]>=min(ttdata$K))*highcap[,t]
	medcap[,t] <- (medcap[,t]>max(ttdata$K))*max(ttdata$K)+(medcap[,t]<min(ttdata$K))*min(ttdata$K)+(medcap[,t]<=max(ttdata$K))*(medcap[,t]>=min(ttdata$K))*medcap[,t]
	lowcap[,t] <- (lowcap[,t]>max(ttdata$K))*max(ttdata$K)+(lowcap[,t]<min(ttdata$K))*min(ttdata$K)+(lowcap[,t]<=max(ttdata$K))*(lowcap[,t]>=min(ttdata$K))*lowcap[,t]
	#At t=2, simulate productivity paths when hit by different sized innovation shocks
	if (t==2){
		#High Innovation Shock
		highxi[,t] <- array(0.9, N)
		#Medium Innovation Shock
		medxi[,t] <- array(0.5, N)
		#Low Innovation Shock
		lowxi[,t] <- array(0.1, N)
	} 
	#Productivity
	highomg[,t] <- rowSums(WX(A=adata[,t], omega=highomg[,t-1])*lspline(vectau=vectau, bvec=parWT, b1=parWTb[1], bL=parWTb[2], u=highxi[,t]))
	#Restrict the Support of High Productivity
	highomg[,t] <- (highomg[,t]>WTminmax[1])*WTminmax[1]+(highomg[,t]<WTminmax[2])*WTminmax[2]+(highomg[,t]<=WTminmax[1])*(highomg[,t]>=WTminmax[2])*highomg[,t]
	#Median Productivity
	medomg[,t] <- rowSums(WX(A=adata[,t], omega=medomg[,t-1])*lspline(vectau=vectau, bvec=parWT, b1=parWTb[1], bL=parWTb[2], u=medxi[,t]))
	#Restrict the Support of Median Productivity
	medomg[,t] <- (medomg[,t]>WTminmax[1])*WTminmax[1]+(medomg[,t]<WTminmax[2])*WTminmax[2]+(medomg[,t]<=WTminmax[1])*(medomg[,t]>=WTminmax[2])*medomg[,t]
	#Low Productivity
	lowomg[,t] <- rowSums(WX(A=adata[,t], omega=lowomg[,t-1])*lspline(vectau=vectau, bvec=parWT, b1=parWTb[1], bL=parWTb[2], u=lowxi[,t]))
	#Restrict the Support of Low Productivity
	lowomg[,t] <- (lowomg[,t]>WTminmax[1])*WTminmax[1]+(lowomg[,t]<WTminmax[2])*WTminmax[2]+(lowomg[,t]<=WTminmax[1])*(lowomg[,t]>=WTminmax[2])*lowomg[,t]
}
for (t in 1:T){
	ttdata <- US %>% group_by(id) %>% slice(t)
	#Labor Paths
	lab1 <- rowSums(LX(A=adata[,t], K=highcap[,t], omega=highomg[,t])*lspline(vectau=vectau, bvec=parL, b1=parLb[1], bL=parLb[2], u=epsdata[,t]))
	lab2 <- rowSums(LX(A=adata[,t], K=medcap[,t], omega=medomg[,t])*lspline(vectau=vectau, bvec=parL, b1=parLb[1], bL=parLb[2], u=epsdata[,t]))
	lab3 <- rowSums(LX(A=adata[,t], K=lowcap[,t], omega=lowomg[,t])*lspline(vectau=vectau, bvec=parL, b1=parLb[1], bL=parLb[2], u=epsdata[,t]))
	#Restrict the Supports
	highlab[,t] <- (lab1>max(ttdata$L))*max(ttdata$L)+(lab1<min(ttdata$L))*min(ttdata$L)+(lab1<=max(ttdata$L))*(lab1>=min(ttdata$L))*lab1
	medlab[,t] <- (lab2>max(ttdata$L))*max(ttdata$L)+(lab2<min(ttdata$L))*min(ttdata$L)+(lab2<=max(ttdata$L))*(lab2>=min(ttdata$L))*lab2
	lowlab[,t] <- (lab3>max(ttdata$L))*max(ttdata$L)+(lab3<min(ttdata$L))*min(ttdata$L)+(lab3<=max(ttdata$L))*(lab3>=min(ttdata$L))*lab3
	#Material Paths
	mat1 <- rowSums(MX(A=adata[,t], K=highcap[,t], omega=highomg[,t])*lspline(vectau=vectau, bvec=parM, b1=parMb[1], bL=parMb[2], u=varepsdata[,t]))
	mat2 <- rowSums(MX(A=adata[,t], K=medcap[,t], omega=medomg[,t])*lspline(vectau=vectau, bvec=parM, b1=parMb[1], bL=parMb[2], u=varepsdata[,t]))
	mat3 <- rowSums(MX(A=adata[,t], K=lowcap[,t], omega=lowomg[,t])*lspline(vectau=vectau, bvec=parM, b1=parMb[1], bL=parMb[2], u=varepsdata[,t]))
	#Restrict the Supports
	highmat[,t] <- (mat1>max(ttdata$M))*max(ttdata$M)+(mat1<min(ttdata$M))*min(ttdata$M)+(mat1<=max(ttdata$M))*(mat1>=min(ttdata$M))*mat1
	medmat[,t] <- (mat2>max(ttdata$M))*max(ttdata$M)+(mat2<min(ttdata$M))*min(ttdata$M)+(mat2<=max(ttdata$M))*(mat2>=min(ttdata$M))*mat2
	lowmat[,t] <- (mat3>max(ttdata$M))*max(ttdata$M)+(mat3<min(ttdata$M))*min(ttdata$M)+(mat3<=max(ttdata$M))*(mat3>=min(ttdata$M))*mat3
}
highlab <- t(apply(highlab, 2, function(x) quantile(x, probs=tauinp)))
medlab <- t(apply(medlab, 2, function(x) quantile(x, probs=tauinp)))
lowlab <- t(apply(lowlab, 2, function(x) quantile(x, probs=tauinp)))
highmat <- t(apply(highmat, 2, function(x) quantile(x, probs=tauinp)))
medmat <- t(apply(medmat, 2, function(x) quantile(x, probs=tauinp)))
lowmat <- t(apply(lowmat, 2, function(x) quantile(x, probs=tauinp)))
# Paths for Labor at different ranks of labor shock
labpath <- data.frame(1:T, lowlab-medlab, highlab-medlab)
names(labpath) <- c("Time", "LowLow", "LowMed", "LowHigh", "HighLow", "HighMed", "HighHigh")
#Paths for Materials at different ranks of materials shock
matpath <- data.frame(1:T, lowmat-medmat, highmat-medmat)
names(matpath) <- c("Time", "LowLow", "LowMed", "LowHigh", "HighLow", "HighMed", "HighHigh")
#Plotting
titles <- c(0.1, 0.1, 0.1, 0.9, 0.9, 0.9)
Lplot <- list()
Mplot <- list()
for (i in 1:6){
	if (i<=3){
		ldat <- data.frame(Time=labpath$Time, Y=labpath[,i+1])
		Lplot[[i]] <- ggplot(ldat, aes(x=Time, y=Y)) + geom_line() + xlab("Time") + ylab("Labor") + coord_cartesian(ylim=c(min(ldat$Y), max(ldat$Y))*2) + ggtitle(expression(paste(tau, "-innovation=0.1")))+ geom_hline(yintercept=0, linetype='dashed', color='red')
		mdat <- data.frame(Time=matpath$Time, Y=matpath[,i+1])
		Mplot[[i]] <- ggplot(mdat, aes(x=Time, y=Y)) + geom_line() + xlab("Time") + ylab("Materials")+ coord_cartesian(ylim=c(min(mdat$Y), max(mdat$Y))*2) + ggtitle(expression(paste(tau, "-innovation=0.1")))+ geom_hline(yintercept=0, linetype='dashed', color='red')
	} else {
		ldat <- data.frame(Time=labpath$Time, Y=labpath[,i+1])
		Lplot[[i]] <- ggplot(ldat, aes(x=Time, y=Y)) + geom_line() + xlab("Time") + ylab("Labor") + coord_cartesian(ylim=c(min(ldat$Y), max(ldat$Y))*2) + ggtitle(expression(paste(tau, "-innovation=0.9")))+ geom_hline(yintercept=0, linetype='dashed', color='red')
		mdat <- data.frame(Time=matpath$Time, Y=matpath[,i+1])
		Mplot[[i]] <- ggplot(mdat, aes(x=Time, y=Y)) + geom_line() + xlab("Time") + ylab("Materials") + coord_cartesian(ylim=c(min(mdat$Y), max(mdat$Y))*2) + ggtitle(expression(paste(tau, "-innovation=0.9")))+ geom_hline(yintercept=0, linetype='dashed', color='red')
	}
}
#Labor
names(Lplot) <- c("LowLow", "LowMed", "LowHigh", "HighLow", "HighMed", "HighHigh")
Ltitle1 <- ggdraw() + draw_label(expression(paste(tau, "-labor=0.1")))
Lrow1 <- plot_grid(Ltitle1, plot_grid(Lplot$LowLow, Lplot$HighLow), ncol=1, rel_heights = c(0.1,1))
Ltitle2 <- ggdraw() + draw_label(expression(paste(tau, "-labor=0.5")))
Lrow2 <- plot_grid(Ltitle2, plot_grid(Lplot$LowMed, Lplot$HighMed), ncol=1, rel_heights = c(0.1,1))
Ltitle3 <- ggdraw() + draw_label(expression(paste(tau, "-labor=0.9")))
Lrow3 <- plot_grid(Ltitle3, plot_grid(Lplot$LowHigh, Lplot$HighHigh), ncol=1, rel_heights = c(0.1,1))
impulseL <- plot_grid(Lrow1, Lrow2, Lrow3, nrow=3)
save_plot("/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Empirical/Capital/Plots/Inputs/impulseL.png", impulseL, base_height = 13, base_width = 10)
#Materials
names(Mplot) <- c("LowLow", "LowMed", "LowHigh", "HighLow", "HighMed", "HighHigh")
Mtitle1 <- ggdraw() + draw_label(expression(paste(tau, "-materials=0.1")))
Mrow1 <- plot_grid(Mtitle1, plot_grid(Mplot$LowLow, Mplot$HighLow), ncol=1, rel_heights = c(0.1,1))
Mtitle2 <- ggdraw() + draw_label(expression(paste(tau, "-materials=0.5")))
Mrow2 <- plot_grid(Mtitle2, plot_grid(Mplot$LowMed, Mplot$HighMed), ncol=1, rel_heights = c(0.1,1))
Mtitle3 <- ggdraw() + draw_label(expression(paste(tau, "-materials=0.9")))
Mrow3 <- plot_grid(Mtitle3, plot_grid(Mplot$LowHigh, Mplot$HighHigh), ncol=1, rel_heights = c(0.1,1))
impulseM <- plot_grid(Mrow1, Mrow2, Mrow3, nrow=3)
save_plot("/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Empirical/Capital/Plots/Inputs/impulseM.png", impulseM, base_height = 13, base_width = 10)
#Capital
highcappath <- data.frame(Time=1:T, High=apply(highcap-medcap, 2, median))
lowcappath <- data.frame(Time=1:T, Low=apply(lowcap-medcap, 2, median))
highcapplot <- ggplot(highcappath, aes(x=Time, y=High)) + geom_line() + xlab("Time") + ylab("Capital") + coord_cartesian(ylim=c(min(highcappath$High), max(highcappath$High))*2)+ geom_hline(yintercept=0, linetype='dashed', color='red')
lowcapplot <- ggplot(lowcappath, aes(x=Time, y=Low)) + geom_line() + xlab("Time") + ylab("Capital") + coord_cartesian(ylim=c(min(lowcappath$Low), max(lowcappath$Low))*2)+ geom_hline(yintercept=0, linetype='dashed', color='red')
impulseK <- plot_grid(lowcapplot, highcapplot, ncol=2)
save_plot("/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Empirical/Capital/Plots/Inputs/impulseK.png", impulseK, base_height=5, base_width = 10)




