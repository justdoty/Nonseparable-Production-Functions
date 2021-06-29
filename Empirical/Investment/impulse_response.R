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
select(id, year, Y, K, K2, L, M, I, age, drate) %>% transmute(id=id, year=year, Y=log(Y), K=log(K), L=log(L), M=log(M), I=log(I), A=age, dp=drate)
#Productivity from LP
omegainit <- omega_est(idvar=US$id, timevar=US$year, Y=US$Y, A=US$A, K=US$K, L=US$L, M=US$M)$omega
#Merge with US Data
US <- US %>% mutate(omega=omegainit)
########################################################################################################
##########################################Load Results############################################
########################################################################################################
set.seed(123456)
load("/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Empirical/Investment/Investment_trans.RData")
vectau <- results$vectau
ntau <- length(vectau)
dims <- results$dims
#Load Parameter Estimates
parY <- results$resYmat
parL <- results$resLmat
parM <- results$resMmat
parI <- results$resImat
parWT <- results$resWTmat
parW1 <- results$resW1mat
parYb <- results$resyb1bLmat
parLb <- results$reslb1bLmat
parMb <- results$resmb1bLmat
parWTb <- results$reswtb1bLmat
parIb <- results$resib1bLmat
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
#Expand the sample
Nsim <- 20
N <- length(unique(US$id))*Nsim
T <- length(unique(US$year))
#Vector of ranks of input demand functions (small, medium, large shocks)
tauinp <- c(0.1, 0.5, 0.9)
#Shocks to productivity
tauxi <- c(0.1, 0.5, 0.9)
#Rank of initial productivity
tauinit <- 0.5
#Age##################################################################
adata <- matrix(0, N, T)
#Labor###############################################################
lnldata <- array(0, c(T, length(tauxi), length(tauinp)))
epsdata <- array(0, c(N, T, length(tauinp)))
#Intermediate Input#####################################################
lnmdata <- array(0, c(T, length(tauxi), length(tauinp)))
varepsdata <- array(0, c(N, T, length(tauinp)))
#Investment###############################################################
lnidata <- array(0, c(N, T, length(tauxi)))
iotadata <- matrix(0.5, nrow=N, ncol=T)
#Capital##################################################################
lnkdata <- array(0, c(N, T, length(tauxi)))
#Productivity#############################################################
omgdata <- array(0, c(N, T, length(tauxi)))
xi <- matrix(runif(N*T), nrow=N, ncol=T)
xidata <- replicate(length(tauxi), xi, simplify="array")
#Ranks for Input Shocks
for (q2 in 1:length(tauinp)){
	epsdata[,,q2] <- array(tauinp[q2], c(N,T))
	varepsdata[,,q2] <- array(tauinp[q2], c(N,T))
}
#########################################################################################################
#Data at T=1
t1data <- US %>% group_by(id) %>% slice(1)
#Initial Age
A1 <- t1data$A
#Give each firm the same age at T=1
adata[,1] <- round(median(A1))
#Initial Productivity (Give each firm the same rank initial  productivity specified by tauinit)
omg1 <- quantile(rnorm(N, mean=parW1[1], sd=sqrt(parW1[2])), tauinit)
#Initial Capital
k1 <- kronecker(array(1, c(Nsim,1)), t1data$K)
#Initial Investment (Same for Each Firm)
i1 <- rowSums(IX(A=adata[,1], K=k1, omega=omg1)*lspline(vectau=vectau, bvec=parI, b1=parIb[1], bL=parIb[2], u=iotadata[,1]))
#Restrict Support of Initial Investment
i1 <- (i1>max(t1data$I))*max(t1data$I)+(i1<min(t1data$I))*min(t1data$I)+(i1<=max(t1data$I))*(i1>=min(t1data$I))*i1
#Initial Inputs
for (q1 in 1:length(tauxi)){
	omgdata[,,q1][,1] <- omg1
	lnkdata[,,q1][,1] <- k1
	lnidata[,,q1][,1] <- i1
	xidata[,,q1][,2] <- tauxi[q1]
	for (q2 in 1:length(tauinp)){
		l1 <- rowSums(LX(A=adata[,1], K=lnkdata[,,q1][,1], omega=omgdata[,,q1][,1])*lspline(vectau=vectau, bvec=parL, b1=parLb[1], bL=parLb[2], u=epsdata[,,q2][,1]))
		m1 <- rowSums(MX(A=adata[,1], K=lnkdata[,,q1][,1], omega=omgdata[,,q1][,1])*lspline(vectau=vectau, bvec=parM, b1=parMb[1], bL=parMb[2], u=varepsdata[,,q2][,1]))
		#Restricting the Supports
		l1 <- (l1>max(t1data$L))*max(t1data$L)+(l1<min(t1data$L))*min(t1data$L)+(l1<=max(t1data$L))*(l1>=min(t1data$L))*l1
		m1 <- (m1>max(t1data$M))*max(t1data$M)+(m1<min(t1data$M))*min(t1data$M)+(m1<=max(t1data$M))*(m1>=min(t1data$M))*m1
		#Median Across Firms
		lnldata[,,q1][1,q2] <- median(l1)
		lnmdata[,,q1][1,q2] <- median(m1)
	}
}
#Evolution of Productivity, Capital, Investment, and Input Decisions for Labor and Materials
#This Loop is Pretty Slow
wmin <- -10
wmax <- 10
for (t in 2:T){
	ttdata <- US %>% group_by(id) %>% slice(t)
	#Age
	adata[,t] <- adata[,t-1]+1
	for (q1 in 1:length(tauxi)){
		omgdata[,,q1][,t] <- rowSums(WX(A=adata[,t], omega=omgdata[,,q1][,t-1])*lspline(vectau=vectau, bvec=parWT, b1=parWTb[1], bL=parWTb[2], u=xidata[,,q1][,t]))
		#Restricting the Supports
		omgdata[,,q1][,t] <- (omgdata[,,q1][,t]>wmax)*wmax+(omgdata[,,q1][,t]<wmin)*wmin+(omgdata[,,q1][,t]<=wmax)*(omgdata[,,q1][,t]>=wmin)*omgdata[,,q1][,t]
		#Generate Capital According to Accumulation Process with Industry-Average Depreciation Rates
		lnkdata[,,q1][,t] <- log(mean(ttdata$dp)*exp(lnkdata[,,q1][,t-1])+exp(lnidata[,,q1][,t-1]))
		#Generate Investment
		lnidata[,,q1][,t] <- rowSums(IX(A=adata[,t], K=lnkdata[,,q1][,t], omega=omgdata[,,q1][,t])*lspline(vectau=vectau, bvec=parI, b1=parIb[1], bL=parIb[2], u=iotadata[,t]))
		#Restricting Supports
		lnidata[,,q1][,t] <- (lnidata[,,q1][,t]>max(ttdata$I))*max(ttdata$I)+(lnidata[,,q1][,t]<min(ttdata$I))*min(ttdata$I)+(lnidata[,,q1][,t]<=max(ttdata$I))*(lnidata[,,q1][,t]>=min(ttdata$I))*lnidata[,,q1][,t]
		#Generate Optimal Input Decisions Following Innovation Shocks for Various Ranks of Demand Functions
		for (q2 in 1:length(tauinp)){
			labdata <- rowSums(LX(A=adata[,t], K=lnkdata[,,q1][,t], omega=omgdata[,,q1][,t])*lspline(vectau=vectau, bvec=parL, b1=parLb[1], bL=parLb[2], u=epsdata[,,q2][,t]))
			matdata <- rowSums(MX(A=adata[,t], K=lnkdata[,,q1][,t], omega=omgdata[,,q1][,t])*lspline(vectau=vectau, bvec=parM, b1=parMb[1], bL=parMb[2], u=varepsdata[,,q2][,t]))
			#Restricting the Supports
			labdata <- (labdata>max(ttdata$L))*max(ttdata$L)+(labdata<min(ttdata$L))*min(ttdata$L)+(labdata<=max(ttdata$L))*(labdata>=min(ttdata$L))*labdata
			matdata <- (matdata>max(ttdata$M))*max(ttdata$M)+(matdata<min(ttdata$M))*min(ttdata$M)+(matdata<=max(ttdata$M))*(matdata>=min(ttdata$M))*matdata
			#Median Across Firms
			lnldata[,,q1][t,q2] <- median(labdata)
			lnmdata[,,q1][t,q2] <- median(matdata)
		}
	}
}
#For the arrays, lnldata and lnmdata, the 1st dimension is time, 2nd is rank of input, 3rd is rank of productivity innovation
#So "Low-Low" represents low labor shock and low innovation shock
labpath <- data.frame(1:T, lnldata[,,1]-lnldata[,,2], lnldata[,,3]-lnldata[,,2])
names(labpath) <- c("Time", "LowLow", "LowMed", "LowHigh", "HighLow", "HighMed", "HighHigh")
#Paths for Materials at different ranks of materials shock
matpath <- data.frame(1:T, lnmdata[,,1]-lnmdata[,,2], lnmdata[,,3]-lnmdata[,,2])
names(matpath) <- c("Time", "LowLow", "LowMed", "LowHigh", "HighLow", "HighMed", "HighHigh")
#Plotting
titles <- c(0.1, 0.1, 0.1, 0.9, 0.9, 0.9)
Lplot <- list()
Mplot <- list()
for (i in 1:6){
	if (i<=3){
		ldat <- data.frame(Time=labpath$Time, Y=labpath[,i+1])
		Lplot[[i]] <- ggplot(ldat, aes(x=Time, y=Y)) + geom_line() + xlab("Time") + ylab("Labor") + coord_cartesian(ylim=c(min(ldat$Y), max(ldat$Y)*2))+ ggtitle(expression(paste(tau, "-innovation=0.1")))+ geom_hline(yintercept=0, linetype='dashed', color='red')
		mdat <- data.frame(Time=matpath$Time, Y=matpath[,i+1])
		Mplot[[i]] <- ggplot(mdat, aes(x=Time, y=Y)) + geom_line() + xlab("Time") + ylab("Materials")+ coord_cartesian(ylim=c(min(mdat$Y), max(mdat$Y)*2)) + ggtitle(expression(paste(tau, "-innovation=0.1")))+ geom_hline(yintercept=0, linetype='dashed', color='red')
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
save_plot("/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Empirical/Investment/Plots/Inputs/impulseL.png", impulseL, base_height = 13, base_width = 10)
#Materials
names(Mplot) <- c("LowLow", "LowMed", "LowHigh", "HighLow", "HighMed", "HighHigh")
Mtitle1 <- ggdraw() + draw_label(expression(paste(tau, "-materials=0.1")))
Mrow1 <- plot_grid(Mtitle1, plot_grid(Mplot$LowLow, Mplot$HighLow), ncol=1, rel_heights = c(0.1,1))
Mtitle2 <- ggdraw() + draw_label(expression(paste(tau, "-materials=0.5")))
Mrow2 <- plot_grid(Mtitle2, plot_grid(Mplot$LowMed, Mplot$HighMed), ncol=1, rel_heights = c(0.1,1))
Mtitle3 <- ggdraw() + draw_label(expression(paste(tau, "-materials=0.9")))
Mrow3 <- plot_grid(Mtitle3, plot_grid(Mplot$LowHigh, Mplot$HighHigh), ncol=1, rel_heights = c(0.1,1))
impulseM <- plot_grid(Mrow1, Mrow2, Mrow3, nrow=3)
save_plot("/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Empirical/Investment/Plots/Inputs/impulseM.png", impulseM, base_height = 13, base_width = 10)
#Capital
cappath <- array(0, c(T, length(tauxi)))
for (q1 in 1:length(tauxi)){
	cappath[,q1] <- apply(lnkdata[,,q1], 2, median)
}
highcappath <- data.frame(Time=1:T, High=cappath[,3]-cappath[,2])
lowcappath <- data.frame(Time=1:T, Low=cappath[,1]-cappath[,2])
highcapplot <- ggplot(highcappath, aes(x=Time, y=High)) + geom_line() + xlab("Time") + ylab("Capital") + coord_cartesian(ylim=c(min(highcappath$High), max(highcappath$High))*2) + geom_hline(yintercept=0, linetype='dashed', color='red')
lowcapplot <- ggplot(lowcappath, aes(x=Time, y=Low)) + geom_line() + xlab("Time") + ylab("Capital") + coord_cartesian(ylim=c(min(lowcappath$Low), max(lowcappath$Low))*2) + geom_hline(yintercept=0, linetype='dashed', color='red')
impulseK <- plot_grid(lowcapplot, highcapplot, ncol=2)
save_plot("/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Empirical/Investment/Plots/Inputs/impulseK.png", impulseK, base_height=5, base_width = 10)









