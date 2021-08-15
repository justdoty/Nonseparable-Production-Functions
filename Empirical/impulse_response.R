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
select(id, year, lny, lnk1, lnk2, lnl, lnm, lni, age) %>% transmute(id=id, year=year, Y=lny, K=lnk1, L=lnl, M=lnm, I=lni, A=age)
#Productivity from LP
omegainit <- omega_est(idvar=US$id, timevar=US$year, Y=US$Y, A=US$A, K=US$K, L=US$L, M=US$M)$omega
US <- US %>% mutate(omega=omegainit) %>% filter(year>=1997) %>% group_by(id) %>% filter(n()>=3)
########################################################################################################
##########################################Load Results############################################
########################################################################################################
set.seed(123456)
load("/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Empirical_trans.RData")
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
wmin <- WTminmax[2]
wmax <- WTminmax[1]
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
tauinit <- c(0.1, 0.5, 0.9)
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
iotadata <- matrix(runif(N*T), nrow=N, ncol=T)
qlnidata <- array(0, c(N, T, length(tauxi), length(tauinp)))
qmedidata <- array(0, c(T, length(tauxi), length(tauinp)))
qiotadata <- array(0, c(N, T, length(tauinp)))
#Capital##################################################################
lnkdata <- array(0, c(N, T, length(tauxi)))
qlnkdata <- array(0, c(N, T, length(tauxi), length(tauinp)))
#Productivity#############################################################
omgdata <- array(0, c(N, T, length(tauxi)))
omgqdata <- array(0, c(N, T, length(tauxi), length(tauinit)))
omgqmat <- array(0, c(T, length(tauxi), length(tauinit)))
xi <- matrix(runif(N*T), nrow=N, ncol=T)
xidata <- replicate(length(tauxi), xi, simplify="array")
#Ranks for Input Shocks
for (q2 in 1:length(tauinp)){
	epsdata[,,q2] <- array(tauinp[q2], c(N,T))
	varepsdata[,,q2] <- array(tauinp[q2], c(N,T))
	qiotadata[,,q2] <- array(tauinp[q2], c(N,T))
}
#########################################################################################################
#Data at T=1
t1data <- US %>% group_by(id) %>% slice(1)
#Initial Age
A1 <- t1data$A
#Give each firm the same age at T=1
adata[,1] <- round(median(A1))
#Initial Productivity (Give each firm the same rank initial  productivity specified by tauinit)
omgq <- quantile(rtruncnorm(N, mean=parW1[1], sd=sqrt(parW1[2])), tauinit)
#Median Shock
omg1 <- omgq[2]
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
		omgqdata[,,,q1][,,q2][,1] <- omgq[q2]
		omgqmat[,,q1][1,q2] <- median(omgqdata[,,,q1][,,q2][,1])
	}
}
#Evolution of Productivity, Capital, Investment, and Input Decisions for Labor and Materials
#This Loop is Pretty Slow
for (t in 2:T){
	ttdata <- US %>% group_by(id) %>% slice(t)
	#Age
	adata[,t] <- adata[,t-1]+1
	for (q1 in 1:length(tauxi)){
		omgdata[,,q1][,t] <- rowSums(WX(A=adata[,t], omega=omgdata[,,q1][,t-1])*lspline(vectau=vectau, bvec=parWT, b1=parWTb[1], bL=parWTb[2], u=xidata[,,q1][,t]))
		#Restricting the Supports
		omgdata[,,q1][,t] <- (omgdata[,,q1][,t]>wmax)*wmax+(omgdata[,,q1][,t]<wmin)*wmin+(omgdata[,,q1][,t]<=wmax)*(omgdata[,,q1][,t]>=wmin)*omgdata[,,q1][,t]
		omgqdata[,,,q1][,,q2][,t] <- (omgqdata[,,,q1][,,q2][,t]>wmax)*wmax+(omgqdata[,,,q1][,,q2][,t]<wmin)*wmin+(omgqdata[,,,q1][,,q2][,t]<=wmax)*(omgqdata[,,,q1][,,q2][,t]>=wmin)*omgqdata[,,,q1][,,q2][,t]
		#Generate Capital According to Accumulation Process with Industry-Average Depreciation Rates
		lnkdata[,,q1][,t] <- log(0.92*exp(lnkdata[,,q1][,t-1])+exp(lnidata[,,q1][,t-1]))
		#Restricting Supports
		lnkdata[,,q1][,t] <- (lnkdata[,,q1][,t]>max(ttdata$K))*max(ttdata$K)+(lnkdata[,,q1][,t]<min(ttdata$K))*min(ttdata$K)+(lnkdata[,,q1][,t]<=max(ttdata$K))*(lnkdata[,,q1][,t]>=min(ttdata$K))*lnkdata[,,q1][,t]
		#Generate Investment
		lnidata[,,q1][,t] <- rowSums(IX(A=adata[,t], K=lnkdata[,,q1][,t], omega=omgdata[,,q1][,t])*lspline(vectau=vectau, bvec=parI, b1=parIb[1], bL=parIb[2], u=iotadata[,t]))
		#Restricting Supports
		lnidata[,,q1][,t] <- (lnidata[,,q1][,t]>max(ttdata$I))*max(ttdata$I)+(lnidata[,,q1][,t]<min(ttdata$I))*min(ttdata$I)+(lnidata[,,q1][,t]<=max(ttdata$I))*(lnidata[,,q1][,t]>=min(ttdata$I))*lnidata[,,q1][,t]
		#Generate Optimal Input Decisions Following Innovation Shocks for Various Ranks of Demand Functions
		for (q2 in 1:length(tauinp)){
			omgqdat <- rowSums(WX(A=adata[,t], omega=omgqdata[,,,q1][,,q2][,t-1])*lspline(vectau=vectau, bvec=parWT, b1=parWTb[1], bL=parWTb[2], u=xidata[,,q1][,t]))
			labdata <- rowSums(LX(A=adata[,t], K=lnkdata[,,q1][,t], omega=omgdata[,,q1][,t])*lspline(vectau=vectau, bvec=parL, b1=parLb[1], bL=parLb[2], u=epsdata[,,q2][,t]))
			matdata <- rowSums(MX(A=adata[,t], K=lnkdata[,,q1][,t], omega=omgdata[,,q1][,t])*lspline(vectau=vectau, bvec=parM, b1=parMb[1], bL=parMb[2], u=varepsdata[,,q2][,t]))
			#Restricting the Supports
			omgqdata[,,,q1][,,q2][,t] <- (omgqdat>wmax)*wmax+(omgqdat<wmin)*wmin+(omgqdat<=wmax)*(omgqdat>=wmin)*omgqdat
			labdata <- (labdata>max(ttdata$L))*max(ttdata$L)+(labdata<min(ttdata$L))*min(ttdata$L)+(labdata<=max(ttdata$L))*(labdata>=min(ttdata$L))*labdata
			matdata <- (matdata>max(ttdata$M))*max(ttdata$M)+(matdata<min(ttdata$M))*min(ttdata$M)+(matdata<=max(ttdata$M))*(matdata>=min(ttdata$M))*matdata
			#Median Across Firms
			omgqmat[,,q1][t,q2] <- median(omgqdata[,,,q1][,,q2][,t])
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
#Paths for Productivity at different ranks of initial productivity shock
omegapath <- data.frame(1:T, omgqmat[,,1]-omgqmat[,,2], omgqmat[,,3]-omgqmat[,,2])
names(omegapath) <- c("Time", "LowLow", "LowMed", "LowHigh", "HighLow", "HighMed", "HighHigh")
#Plotting
titles <- c(0.1, 0.1, 0.1, 0.9, 0.9, 0.9)
Wplot <- list()
Lplot <- list()
Mplot <- list()
for (i in 1:6){
	if (i<=3){
		wdat <- data.frame(Time=omegapath$Time, Y=omegapath[,i+1])
		Wplot[[i]] <- ggplot(wdat, aes(x=Time, y=Y)) + geom_line() + xlab("Time") + ylab("Productivity") + coord_cartesian(ylim=c(min(wdat$Y), max(wdat$Y)*2))+ ggtitle(expression(paste(tau, "-innovation=0.1")))+ geom_hline(yintercept=0, linetype='dashed', color='red') 
		ldat <- data.frame(Time=labpath$Time, Y=labpath[,i+1])
		Lplot[[i]] <- ggplot(ldat, aes(x=Time, y=Y)) + geom_line() + xlab("Time") + ylab("Labor") + coord_cartesian(ylim=c(min(ldat$Y), max(ldat$Y)*2))+ ggtitle(expression(paste(tau, "-innovation=0.1")))+ geom_hline(yintercept=0, linetype='dashed', color='red')
		mdat <- data.frame(Time=matpath$Time, Y=matpath[,i+1])
		Mplot[[i]] <- ggplot(mdat, aes(x=Time, y=Y)) + geom_line() + xlab("Time") + ylab("Materials")+ coord_cartesian(ylim=c(min(mdat$Y), max(mdat$Y)*2)) + ggtitle(expression(paste(tau, "-innovation=0.1")))+ geom_hline(yintercept=0, linetype='dashed', color='red')
	} else {
		wdat <- data.frame(Time=omegapath$Time, Y=omegapath[,i+1])
		Wplot[[i]] <- ggplot(wdat, aes(x=Time, y=Y)) + geom_line() + xlab("Time") + ylab("Productivity") + coord_cartesian(ylim=c(min(wdat$Y), max(wdat$Y))*2) + ggtitle(expression(paste(tau, "-innovation=0.9")))+ geom_hline(yintercept=0, linetype='dashed', color='red')
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
save_plot("/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Empirical/Plots/Inputs/impulseL.png", impulseL, base_height = 13, base_width = 10)
#Materials
names(Mplot) <- c("LowLow", "LowMed", "LowHigh", "HighLow", "HighMed", "HighHigh")
Mtitle1 <- ggdraw() + draw_label(expression(paste(tau, "-materials=0.1")))
Mrow1 <- plot_grid(Mtitle1, plot_grid(Mplot$LowLow, Mplot$HighLow), ncol=1, rel_heights = c(0.1,1))
Mtitle2 <- ggdraw() + draw_label(expression(paste(tau, "-materials=0.5")))
Mrow2 <- plot_grid(Mtitle2, plot_grid(Mplot$LowMed, Mplot$HighMed), ncol=1, rel_heights = c(0.1,1))
Mtitle3 <- ggdraw() + draw_label(expression(paste(tau, "-materials=0.9")))
Mrow3 <- plot_grid(Mtitle3, plot_grid(Mplot$LowHigh, Mplot$HighHigh), ncol=1, rel_heights = c(0.1,1))
impulseM <- plot_grid(Mrow1, Mrow2, Mrow3, nrow=3)
save_plot("/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Empirical/Plots/Inputs/impulseM.png", impulseM, base_height = 13, base_width = 10)
#Productivity
names(Wplot) <- c("LowLow", "LowMed", "LowHigh", "HighLow", "HighMed", "HighHigh")
Wtitle1 <- ggdraw() + draw_label(expression(paste(tau, "-initial=0.1")))
Wrow1 <- plot_grid(Wtitle1, plot_grid(Wplot$LowLow, Wplot$HighLow), ncol=1, rel_heights = c(0.1,1))
Wtitle2 <- ggdraw() + draw_label(expression(paste(tau, "-initial=0.5")))
Wrow2 <- plot_grid(Wtitle2, plot_grid(Wplot$LowMed, Wplot$HighMed), ncol=1, rel_heights = c(0.1,1))
Wtitle3 <- ggdraw() + draw_label(expression(paste(tau, "-initial=0.9")))
Wrow3 <- plot_grid(Wtitle3, plot_grid(Wplot$LowHigh, Wplot$HighHigh), ncol=1, rel_heights = c(0.1,1))
impulseW <- plot_grid(Wrow1, Wrow2, Wrow3, nrow=3)
save_plot("/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Empirical/Plots/Inputs/impulseW.png", impulseW, base_height = 13, base_width = 10)
#Capital##################################################################################################
#I do the same procedure for capital except at different levels of shock to investment
#Initial Capital and Investment
for (q1 in 1:length(tauxi)){
	for (q2 in 1:length(tauinp)){
		qlnkdata[,,,q2][,,q1][,1] <- k1
		qi1 <- rowSums(IX(A=adata[,1], K=k1, omega=omg1)*lspline(vectau=vectau, bvec=parI, b1=parIb[1], bL=parIb[2], u=qiotadata[,,q2][,1]))
		#Restrict the Support
		qi1 <- (qi1>max(t1data$I))*max(t1data$I)+(qi1<min(t1data$I))*min(t1data$I)+(qi1<=max(t1data$I))*(qi1>=min(t1data$I))*qi1
		#Initial Investment
		qlnidata[,,,q2][,,q1][,1] <- qi1
		qmedidata[,,q1][1,q2] <- median(qi1)
	}
}
for (t in 2:T){
	#Data at time t
	ttdata <- US %>% group_by(id) %>% slice(t)
	for (q1 in 1:length(tauxi)){
		for (q2 in 1:length(tauinp)){
			qlnkdata[,,,q1][,,q2][,t] <- log(0.92*exp(qlnkdata[,,,q1][,,q2][,t-1])+exp(qlnidata[,,,q1][,,q2][,t-1]))
			qit <- rowSums(IX(A=adata[,t], K=qlnkdata[,,,q1][,,q2][,t], omega=omgdata[,,q1][,t])*lspline(vectau=vectau, bvec=parI, b1=parIb[1], bL=parIb[2], u=qiotadata[,,q2][,t]))
			#Restrict the Support of Investment
			qit <- (qit>max(ttdata$I))*max(ttdata$I)+(qit<min(ttdata$I))*min(ttdata$I)+(qit<=max(ttdata$I))*(qit>=min(ttdata$I))*qit
			qmedidata[,,q1][t,q2] <- median(qit)
		}
	}
}
cappath <- data.frame(1:T, qmedidata[,,1]-qmedidata[,,2], qmedidata[,,3]-qmedidata[,,2])
names(cappath) <- c("Time", "LowLow", "LowMed", "LowHigh", "HighLow", "HighMed", "HighHigh")
#Plotting
titles <- c(0.1, 0.1, 0.1, 0.9, 0.9, 0.9)
Kplot <- list()
for (i in 1:6){
	if (i<=3){
		kdat <- data.frame(Time=cappath$Time, Y=cappath[,i+1])
		Kplot[[i]] <- ggplot(kdat, aes(x=Time, y=Y)) + geom_line() + xlab("Time") + ylab("Capital") + coord_cartesian(ylim=c(min(kdat$Y), max(kdat$Y)*2))+ ggtitle(expression(paste(tau, "-innovation=0.1")))+ geom_hline(yintercept=0, linetype='dashed', color='red')
	} else {
		kdat <- data.frame(Time=cappath$Time, Y=cappath[,i+1])
		Kplot[[i]] <- ggplot(kdat, aes(x=Time, y=Y)) + geom_line() + xlab("Time") + ylab("Capital") + coord_cartesian(ylim=c(min(kdat$Y), max(kdat$Y))*2) + ggtitle(expression(paste(tau, "-innovation=0.9")))+ geom_hline(yintercept=0, linetype='dashed', color='red')
	}
}
names(Kplot) <- c("LowLow", "LowMed", "LowHigh", "HighLow", "HighMed", "HighHigh")
Ktitle1 <- ggdraw() + draw_label(expression(paste(tau, "-investment=0.1")))
Krow1 <- plot_grid(Ktitle1, plot_grid(Kplot$LowLow, Kplot$HighLow), ncol=1, rel_heights = c(0.1,1))
Ktitle2 <- ggdraw() + draw_label(expression(paste(tau, "-investment=0.5")))
Krow2 <- plot_grid(Ktitle2, plot_grid(Kplot$LowMed, Kplot$HighMed), ncol=1, rel_heights = c(0.1,1))
Ktitle3 <- ggdraw() + draw_label(expression(paste(tau, "-investment=0.9")))
Krow3 <- plot_grid(Ktitle3, plot_grid(Kplot$LowHigh, Kplot$HighHigh), ncol=1, rel_heights = c(0.1,1))
impulseK <- plot_grid(Krow1, Krow2, Krow3, nrow=3)
save_plot("/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Empirical/Plots/Inputs/impulseK.png", impulseK, base_height = 13, base_width = 10)









