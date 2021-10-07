require(cowplot)
require(dplyr)
require(purrr)
require(xtable)
require(stringr)
require(truncnorm)
require(RColorBrewer)
require(reshape2)
require(plotly)
library(listviewer)
source('/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Functions/Main/Auxfuns.R')
source('/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Functions/Main/Tensors.R')
source('/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Functions/Main/omega.R')
#Load US Dataset
US <- read.csv('/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Data/USdata.csv')  %>%
	select(id, year, lny, lnk1, lnk2, lnl, lnm, lni, rd) %>% transmute(id=id, time=year, Y=lny, K=lnk1, L=lnl, M=lnm, I=lni, R=rd)
#Productivity from LP
omegainit <- omega_est(idvar=US$id, timevar=US$year, Y=US$Y, A=US$A, K=US$K, L=US$L, M=US$M)$omega
#Add LP Estimate of productivity
US <- US %>% mutate(omega=omegainit) 
idcon <- duplicated(US$id)
idlag <- duplicated(US$id, fromLast=TRUE)
id1 <- !idcon
########################################################################################################
#Load Results
########################################################################################################
set.seed(123456)
load("/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Environments/main.RData")
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
parIb <- results$resib1bLmat
parWTb <- results$reswtb1bLmat
parW1b <- results$resw1b1bLmat
WTminmax <- results$maxminwtmat
#De-mean
US <- US %>% mutate(Y=Y-mean(Y), K=K-mean(K), L=L-mean(L), M=M-mean(M), I=I-mean(I))
wmin <- min(US$Y)
wmax <- max(US$Y)
#Load the method used (e.g. Cobb, Translog, Hermite)
method <- results$method
resY <- results$resY
###############################################################################
#Simulate Productivity and Capital Evolution Given Model Parameters
##############################################################################
#Expand the sample
Nsim <- 20
N <- length(unique(US$id))*Nsim
T <- length(unique(US$time))
#Vector of ranks of input demand functions (small, medium, large)
tauinp <- c(0.1, 0.5, 0.9)
#Shocks to productivity
tauxi <- c(0.1, 0.5, 0.9)
#Rank of initial productivity
tauinit <- c(0.1, 0.5, 0.9)
#Labor###############################################################
lnldata <- array(0, c(N, T, length(tauxi), length(tauinp)))
#At the median
lnlmed <- array(0, c(T, length(tauxi), length(tauinp)))
#Shocks to Labor
epsl <- array(0, c(N, T, length(tauinp)))
#Intermediate Input#####################################################
lnmdata <- array(0, c(N, T, length(tauxi), length(tauinp)))
#At the median
lnmmed <- array(0, c(T, length(tauxi), length(tauinp)))
#Shocks to materials
epsm <- array(0, c(N, T, length(tauinp)))
#Investment###############################################################
lnidata <- array(0, c(N, T, length(tauxi)))
iota <- matrix(runif(N*T), nrow=N, ncol=T)
#This is for when i vary the shocks to investment demand
qlnidata <- array(0, c(N, T, length(tauxi), length(tauinp)))
qlnimed <- array(0, c(T, length(tauxi), length(tauinp)))
qiota <- array(0, c(N, T, length(tauinp)))
#Capital##################################################################
lnkdata <- array(0, c(N, T, length(tauxi)))
#This is for when i vary the shocks to investment demand
qlnkdata <- array(0, c(N, T, length(tauxi), length(tauinp)))
qlnkmed <- array(0, c(T, length(tauxi), length(tauinp)))
#Productivity#############################################################
#This is used for labor, materials, and investment
omgdata <- array(0, c(N, T, length(tauxi)))
#This is for when i vary productivity paths for different initial productivity levels
omgqdata <- array(0, c(N, T, length(tauxi), length(tauinit)))
omgqmed <- array(0, c(T, length(tauxi), length(tauinp)))
#Innovation Shocks
xi <- matrix(runif(N*T), nrow=N, ncol=T)
xidata <- replicate(length(tauxi), xi, simplify="array")
#Ranks for Input Shocks
for (q2 in 1:length(tauinp)){
	epsl[,,q2] <- array(tauinp[q2], c(N,T))
	epsm[,,q2] <- array(tauinp[q2], c(N,T))
	qiota[,,q2] <- array(tauinp[q2], c(N,T))
}
#########################################################################################################
#Data at T=1
t1data <- US %>% group_by(id) %>% slice(1)
#Initial Capital
k1 <- kronecker(array(1, c(Nsim,1)), t1data$K)
#Initial Productivity
omg1 <- rowSums(WX1(K=k1)*lspline(vectau=vectau, bvec=parW1, b1=parW1b[1], bL=parW1b[2], u=xi[,1]))
#Restrict the Support of Initial Productivity
omg1 <- (omg1>wmax)*wmax+(omg1<wmin)*wmin+(omg1<=wmax)*(omg1>=wmin)*omg1
omg1q <- quantile(omg1, tauinit)
#Median Shock
omg1m <- omg1q[2]
#Initial Investment (at median productivity)
i1 <- rowSums(IX(K=k1, omega=omg1m)*lspline(vectau=vectau, bvec=parI, b1=parIb[1], bL=parIb[2], u=iota[,1]))
#Restrict Support of Initial Investment
i1 <- (i1>max(t1data$I))*max(t1data$I)+(i1<min(t1data$I))*min(t1data$I)+(i1<=max(t1data$I))*(i1>=min(t1data$I))*i1
#Initial Inputs
for (q1 in 1:length(tauxi)){
	omgdata[,,q1][,1] <- omg1m
	lnkdata[,,q1][,1] <- k1
	lnidata[,,q1][,1] <- i1
	#At t=2 give varying shock to productivity given by tauxi
	xidata[,,q1][,2] <- tauxi[q1]
	for (q2 in 1:length(tauinp)){
		k1 <- rep(mean(k1), length(k1))
		l1 <- rowSums(LX(K=k1, omega=omg1m)*lspline(vectau=vectau, bvec=parL, b1=parLb[1], bL=parLb[2], u=epsl[,,q2][,1]))
		m1 <- rowSums(MX(K=k1, L=l1, omega=omg1m)*lspline(vectau=vectau, bvec=parM, b1=parMb[1], bL=parMb[2], u=epsm[,,q2][,1]))
		# #Restricting the Supports
		l1 <- (l1>max(t1data$L))*max(t1data$L)+(l1<min(t1data$L))*min(t1data$L)+(l1<=max(t1data$L))*(l1>=min(t1data$L))*l1
		m1 <- (m1>max(t1data$M))*max(t1data$M)+(m1<min(t1data$M))*min(t1data$M)+(m1<=max(t1data$M))*(m1>=min(t1data$M))*m1
		# #Median Across Firms
		lnldata[,,,q2][,,q1][,1] <- l1
		lnmdata[,,,q2][,,q1][,1] <- m1
		#Different levels of initial productivity
		omgqdata[,,,q2][,,q1][,1] <- omg1q[q2]
	}
}
for (t in 2:T){
	ttdata <- US %>% group_by(id) %>% slice(t)
	for (q1 in 1:length(tauxi)){
		omg <- rowSums(WX(omega=omgdata[,,q1][,t-1])*lspline(vectau=vectau, bvec=parWT, b1=parWTb[1], bL=parWTb[2], u=xidata[,,q1][,t]))
		#Restricting the Supports
		omgdata[,,q1][,t] <- (omg>wmax)*wmax+(omg<wmin)*wmin+(omg<=wmax)*(omg>=wmin)*omg
		#Generate Capital According to Accumulation Process with Depreciation
		lnk <- log(0.98*exp(lnkdata[,,q1][,t-1])+exp(lnidata[,,q1][,t-1]))
		#Restricting Supports
		lnkdata[,,q1][,t] <- (lnk>max(ttdata$K))*max(ttdata$K)+(lnk<min(ttdata$K))*min(ttdata$K)+(lnk<=max(ttdata$K))*(lnk>=min(ttdata$K))*lnk
		#Generate Investment
		lni <- rowSums(IX(K=lnkdata[,,q1][,t], omega=omgdata[,,q1][,t])*lspline(vectau=vectau, bvec=parI, b1=parIb[1], bL=parIb[2], u=iota[,t]))
		#Restricting Supports
		lnidata[,,q1][,t] <- (lni>max(ttdata$I))*max(ttdata$I)+(lni<min(ttdata$I))*min(ttdata$I)+(lni<=max(ttdata$I))*(lni>=min(ttdata$I))*lni
		#Generate Productivity
		for (q2 in 1:length(tauinp)){
			omgq <- rowSums(WX(omega=omgqdata[,,,q2][,,q1][,t-1])*lspline(vectau=vectau, bvec=parWT, b1=parWTb[1], bL=parWTb[2], u=xidata[,,q1][,t]))
			#Restricting the Supports
			omgqdata[,,,q2][,,q1][,t] <- (omgq>wmax)*wmax+(omgq<wmin)*wmin+(omgq<=wmax)*(omgq>=wmin)*omgq
		}
	}
}
#Generate Labor
for (q1 in 1:length(tauxi)){
	for (q2 in 1:length(tauinp)){
		for (t in 2:T){
			ttdata <- US %>% group_by(id) %>% slice(t)
				lab  <- rowSums(LX(K=mean(lnkdata[,,q1]), omega=omgdata[,,q1][,t])*lspline(vectau=vectau, bvec=parL, b1=parLb[1], bL=parLb[2], u=epsl[,,q2][,t]))
				#Restricting the Supports
				lnldata[,,,q2][,,q1][,t] <- (lab>max(ttdata$L))*max(ttdata$L)+(lab<min(ttdata$L))*min(ttdata$L)+(lab<=max(ttdata$L))*(lab>=min(ttdata$L))*lab
		}
	}
}
#Generate Materials
for (q1 in 1:length(tauxi)){
	for (q2 in 1:length(tauinp)){
		for (t in 2:T){
			ttdata <- US %>% group_by(id) %>% slice(t)
			mat <- rowSums(MX(K=mean(lnkdata[,,q1]), L=mean(lnldata[,,,q2][,,q1]), omega=omgdata[,,q1][,t])*lspline(vectau=vectau, bvec=parM, b1=parMb[1], bL=parMb[2], u=epsm[,,q2][,t]))
			#Restricting the Supports
			lnmdata[,,,q2][,,q1][,t] <- (mat>max(ttdata$M))*max(ttdata$M)+(mat<min(ttdata$M))*min(ttdata$M)+(mat<=max(ttdata$M))*(mat>=min(ttdata$M))*mat
		}
	}
}
#Median across firms
for (q1 in 1:length(tauxi)){
	for (q2 in 1:length(tauinp)){
	omgqmed[,,q2][,q1] <- apply(omgqdata[,,,q2][,,q1], 2, median)
	lnlmed[,,q2][,q1] <- apply(lnldata[,,,q2][,,q1], 2, median)
	lnmmed[,,q2][,q1] <- apply(lnmdata[,,,q2][,,q1], 2, median)
	}
}
#For the arrays, lnldata and lnmdata, the 1st dimension is time, 2nd is rank of innovation shock, 3rd is rank of input shock
#So "Low-Low" represents low labor shock and low innovation shock
labpath <- data.frame(1:T, lnlmed[,1,]-lnlmed[,2,], lnlmed[,3,]-lnlmed[,2,])
names(labpath) <- c("Time", "LowLow", "LowMed", "LowHigh", "HighLow", "HighMed", "HighHigh")
#Save to Selection_Bias Folder
write.csv(labpath, "/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Functions/Selection_Bias/labpath.csv")
#Paths for Materials at different ranks of materials shock
matpath <- data.frame(1:T, lnmmed[,1,]-lnmmed[,2,], lnmmed[,3,]-lnmmed[,2,])
names(matpath) <- c("Time", "LowLow", "LowMed", "LowHigh", "HighLow", "HighMed", "HighHigh")
#Save to Selection_Bias Folder
write.csv(matpath, "/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Functions/Selection_Bias/matpath.csv")
#Paths for Productivity at different ranks of initial productivity shock
omegapath <- data.frame(1:T, omgqmed[,1,]-omgqmed[,2,], omgqmed[,3,]-omgqmed[,2,])
names(omegapath) <- c("Time", "LowLow", "LowMed", "LowHigh", "HighLow", "HighMed", "HighHigh")
#Save to Selection_Bias Folder
write.csv(omegapath, "/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Functions/Selection_Bias/omegapath.csv")
#Plotting
taup <- c("0.1", "0.5", "0.9", "0.1", "0.5", "0.9")
taualp <- c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)")
ltitles <- lapply(1:6, function(i) parse(text = paste0("tau[paste(l, \"=\", ", taup[i], ")]")))
mtitles <- lapply(1:6, function(i) parse(text = paste0("tau[paste(m, \"=\", ", taup[i], ")]")))
wtitles <- lapply(1:6, function(i) parse(text = paste0("tau[paste(omega[1], \"=\", ", taup[i], ")]")))
# wann <- lapply(1:6, function(i) list(text=expression(tau), xref="paper", yref="paper", yanchor="bottom", xanchor="center", align="center", x=0.5, y=1, showarrow=FALSE))
Wplot <- list()
Wplotly <- list()
Lplot <- list()
Lplotly <- list()
Mplot <- list()
Mplotly <- list()
for (i in 1:6){
	wdat <- data.frame(Time=omegapath$Time, Y=omegapath[,i+1])
	Wplot[[i]] <- ggplot(wdat, aes(x=Time, y=Y)) + geom_line() + xlab("Years") + ylab("Productivity") + coord_cartesian(ylim=c(min(omegapath[,-1]), max(omegapath[,-1]))) +  geom_hline(yintercept=0, linetype='dashed', color='red') + labs(title=taualp[i], subtitle=wtitles[[i]]) + theme(plot.title = element_text(size = 20), plot.subtitle=element_text(size = 20, hjust=0.5))
	Wplotly[[i]] <- plot_ly(wdat, x=~Time, y=~Y, type = 'scatter', mode = 'lines', showlegend=F, line=list(color="black")) %>% add_trace(y = 0, mode = 'lines', showlegend=F,line=list(color="red", dash="dash")) %>% layout(xaxis=list(title="Years"), yaxis=list(title="Productivity", range=list(min(omegapath[,-1]), max(omegapath[,-1]))))
	ldat <- data.frame(Time=labpath$Time, Y=labpath[,i+1])
	Lplot[[i]] <- ggplot(ldat, aes(x=Time, y=Y)) + geom_line() + xlab("Years") + ylab("Labor") + coord_cartesian(ylim=c(min(labpath[,-1]), max(labpath[,-1]))) + geom_hline(yintercept=0, linetype='dashed', color='red')+ labs(title=taualp[i], subtitle=ltitles[[i]]) + theme(plot.title = element_text(size = 20), plot.subtitle=element_text(size = 20, hjust=0.5))
	Lplotly[[i]] <- plot_ly(ldat, x=~Time, y=~Y, type = 'scatter', mode = 'lines', showlegend=F, line=list(color="black")) %>% add_trace(y = 0, mode = 'lines', showlegend=F,line=list(color="red", dash="dash")) %>% layout(xaxis=list(title="Years"), yaxis=list(title="Labor", range=list(min(labpath[,-1]), max(labpath[,-1]))))
	mdat <- data.frame(Time=matpath$Time, Y=matpath[,i+1])
	Mplot[[i]] <- ggplot(mdat, aes(x=Time, y=Y)) + geom_line() + xlab("Years") + ylab("Materials") + coord_cartesian(ylim=c(min(matpath[,-1]), max(matpath[,-1]))) + geom_hline(yintercept=0, linetype='dashed', color='red')+ labs(title=taualp[i], subtitle=mtitles[[i]]) + theme(plot.title = element_text(size = 20), plot.subtitle=element_text(size = 20, hjust=0.5))
	Mplotly[[i]] <- plot_ly(mdat, x=~Time, y=~Y, type = 'scatter', mode = 'lines', showlegend=F, line=list(color="black")) %>% add_trace(y = 0, mode = 'lines', showlegend=F,line=list(color="red", dash="dash")) %>% layout(xaxis=list(title="Years"), yaxis=list(title="Materials", range=list(min(matpath[,-1]), max(matpath[,-1]))))
}
#Labor
names(Lplot) <- c("LowLow", "LowMed", "LowHigh", "HighLow", "HighMed", "HighHigh")
Ltitle1 <- ggdraw() + draw_label(expression(paste(tau[xi], "=0.1")), fontface="bold", x=0.53, size=22)
Lrow1 <-  plot_grid(Ltitle1, plot_grid(Lplot$LowLow, Lplot$LowMed, Lplot$LowHigh, nrow=1), ncol=1, rel_heights=c(0.5,1))
Ltitle2 <- ggdraw() + draw_label(expression(paste(tau[xi], "=0.9")), fontface="bold", x=0.53, size=22)
Lrow2 <- plot_grid(Ltitle2, plot_grid(Lplot$HighLow, Lplot$HighMed, Lplot$HighHigh, nrow=1), ncol=1, rel_heights=c(0.5, 1))
impulseL <- plot_grid(Lrow1, Lrow2, nrow=2)
save_plot("/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Figures/Main/impulseL.png", impulseL, base_height = 10, base_width = 13)
#Materials
names(Mplot) <- c("LowLow", "LowMed", "LowHigh", "HighLow", "HighMed", "HighHigh")
Mtitle1 <- ggdraw() + draw_label(expression(paste(tau[xi], "=0.1")), fontface="bold", x=0.53, size=22)
Mrow1 <-  plot_grid(Mtitle1, plot_grid(Mplot$LowLow, Mplot$LowMed, Mplot$LowHigh, nrow=1), ncol=1, rel_heights=c(0.5,1))
Mtitle2 <- ggdraw() + draw_label(expression(paste(tau[xi], "=0.9")), fontface="bold", x=0.53, size=22)
Mrow2 <- plot_grid(Mtitle2, plot_grid(Mplot$HighLow, Mplot$HighMed, Mplot$HighHigh, nrow=1), ncol=1, rel_heights=c(0.5, 1))
impulseM <- plot_grid(Mrow1, Mrow2, nrow=2)
save_plot("/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Figures/Main/impulseM.png", impulseM, base_height = 10, base_width = 13)
#Productivity
names(Wplot) <- c("LowLow", "LowMed", "LowHigh", "HighLow", "HighMed", "HighHigh")
Wtitle1 <- ggdraw() + draw_label(expression(paste(tau[xi], "=0.1")), fontface="bold", x=0.53, size=22)
Wrow1 <-  plot_grid(Wtitle1, plot_grid(Wplot$LowLow, Wplot$LowMed, Wplot$LowHigh, nrow=1), ncol=1, rel_heights=c(0.5,1))
Wtitle2 <- ggdraw() + draw_label(expression(paste(tau[xi], "=0.9")), fontface="bold", x=0.53, size=22)
Wrow2 <- plot_grid(Wtitle2, plot_grid(Wplot$HighLow, Wplot$HighMed, Wplot$HighHigh, nrow=1), ncol=1, rel_heights=c(0.5, 1))
impulseW <- plot_grid(Wrow1, Wrow2, nrow=2)
save_plot("/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Figures/Main/impulseW.png", impulseW, base_height = 10, base_width = 13)
# Capital##################################################################################################
# I do the same procedure for capital except at different levels of shock to investment
#Initial Capital and Investment
for (q1 in 1:length(tauxi)){
	for (q2 in 1:length(tauinp)){
		qlnkdata[,,,q2][,,q1][,1] <- k1
		qlnkmed[,,q2][1,q1] <- median(qlnkdata[,,,q2][,,q1][,1])
		qi1 <- rowSums(IX(K=k1, omega=omg1m)*lspline(vectau=vectau, bvec=parI, b1=parIb[1], bL=parIb[2], u=qiota[,,q2][,1]))
		#Restrict the Support
		qi1 <- (qi1>max(t1data$I))*max(t1data$I)+(qi1<min(t1data$I))*min(t1data$I)+(qi1<=max(t1data$I))*(qi1>=min(t1data$I))*qi1
		#Initial Investment
		qlnidata[,,,q2][,,q1][,1] <- qi1
		qlnimed[,,q2][1,q1] <- median(qi1)
	}
}
for (t in 2:T){
	#Data at time t
	ttdata <- US %>% group_by(id) %>% slice(t)
	for (q1 in 1:length(tauxi)){
		for (q2 in 1:length(tauinp)){
			kit <- log(0.98*exp(qlnkdata[,,,q2][,,q1][,t-1])+exp(qlnidata[,,,q2][,,q1][,t-1]))
			#Restrict the Support of Capital
			qlnkdata[,,,q2][,,q1][,t] <- (kit>max(ttdata$K))*max(ttdata$K)+(kit<min(ttdata$K))*min(ttdata$K)+(kit<=max(ttdata$K))*(kit>=min(ttdata$K))*kit
			qit <- rowSums(IX(K=mean(qlnkdata[,,,q2][,,q1][,t]), omega=omgdata[,,q1][,t])*lspline(vectau=vectau, bvec=parI, b1=parIb[1], bL=parIb[2], u=qiota[,,q2][,t]))
			#Restrict the Support of Investment
			qlnidata[,,,q2][,,q1][,t] <- (qit>max(ttdata$I))*max(ttdata$I)+(qit<min(ttdata$I))*min(ttdata$I)+(qit<=max(ttdata$I))*(qit>=min(ttdata$I))*qit
		}
	}
}
#Median across firms
for (q1 in 1:length(tauxi)){
	for (q2 in 1:length(tauinp)){
	qlnimed[,,q2][,q1] <- apply(qlnidata[,,,q2][,,q1], 2, median)
	qlnkmed[,,q2][,q1] <- apply(qlnkdata[,,,q2][,,q1], 2, median)
	}
}
ipath <- data.frame(1:T, qlnimed[,1,]-qlnimed[,2,], qlnimed[,3,]-qlnimed[,2,])
names(ipath) <- c("Time", "LowLow", "LowMed", "LowHigh", "HighLow", "HighMed", "HighHigh")
write.csv(ipath, "/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Functions/Selection_Bias/ipath.csv")
ititles <- lapply(1:6, function(i) parse(text = paste0("tau[paste(zeta, \"=\", ", taup[i], ")]")))
#Plotting
Iplot <- list()
Iplotly <- list()
for (i in 1:6){
	idat <- data.frame(Time=ipath$Time, Y=ipath[,i+1])
	Iplot[[i]] <- ggplot(idat, aes(x=Time, y=Y)) + geom_line() + xlab("Years") + ylab("Investment")  + coord_cartesian(ylim=c(min(ipath[,-1]), max(ipath[,-1]))) + geom_hline(yintercept=0, linetype='dashed', color='red')+ labs(title=taualp[i], subtitle=ititles[[i]]) + theme(plot.title = element_text(size = 20), plot.subtitle=element_text(size = 20, hjust=0.5))
	Iplotly[[i]] <- plot_ly(idat, x=~Time, y=~Y, type = 'scatter', mode = 'lines', showlegend=F, line=list(color="black")) %>% add_trace(y = 0, mode = 'lines', showlegend=F,line=list(color="red", dash="dash")) %>% layout(xaxis=list(title="Years"), yaxis=list(title="Investment", range=list(min(ipath[,-1]), max(ipath[,-1]))))
}
names(Iplot) <- c("LowLow", "LowMed", "LowHigh", "HighLow", "HighMed", "HighHigh")
Ititle1 <- ggdraw() + draw_label(expression(paste(tau[xi], "=0.1")), fontface="bold", x=0.53, size=22)
Irow1 <-  plot_grid(Ititle1, plot_grid(Iplot$LowLow, Iplot$LowMed, Iplot$LowHigh, nrow=1), ncol=1, rel_heights=c(0.5,1))
Ititle2 <- ggdraw() + draw_label(expression(paste(tau[xi], "=0.9")), fontface="bold", x=0.53, size=22)
Irow2 <- plot_grid(Ititle2, plot_grid(Iplot$HighLow, Iplot$HighMed, Iplot$HighHigh, nrow=1), ncol=1, rel_heights=c(0.5, 1))
impulseI <- plot_grid(Irow1, Irow2, nrow=2)
save_plot("/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Figures/Main/impulseI.png", impulseI, base_height = 10, base_width = 13)
#Capital
kpath <- data.frame(1:T, qlnkmed[,1,]-qlnkmed[,2,], qlnkmed[,3,]-qlnkmed[,2,])
names(kpath) <- c("Time", "LowLow", "LowMed", "LowHigh", "HighLow", "HighMed", "HighHigh")
write.csv(kpath, "/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Functions/Selection_Bias/kpath.csv")
Ktitles <- lapply(1:6, function(i) parse(text = paste0("tau[paste(i, \"=\", ", taup[i], ")]")))
#Plotting
Kplot <- list()
Kplotly <- list()
for (i in 1:6){
	kdat <- data.frame(Time=kpath$Time, Y=kpath[,i+1])
	Kplot[[i]] <- ggplot(kdat, aes(x=Time, y=Y)) + geom_line() + xlab("Years") + ylab("Capital")  + coord_cartesian(ylim=c(min(kpath[,-1]), max(kpath[,-1]))) + geom_hline(yintercept=0, linetype='dashed', color='red')+ labs(title=taualp[i], subtitle=Ktitles[[i]]) + theme(plot.title = element_text(size = 20), plot.subtitle=element_text(size = 20, hjust=0.5))
	Kplotly[[i]] <- plot_ly(kdat, x=~Time, y=~Y, type = 'scatter', mode = 'lines', showlegend=F, line=list(color="black")) %>% add_trace(y = 0, mode = 'lines', showlegend=F,line=list(color="red", dash="dash")) %>% layout(xaxis=list(title="Years"), yaxis=list(title="Capital", range=list(min(kpath[,-1]), max(kpath[,-1]))))
}
names(Kplot) <- c("LowLow", "LowMed", "LowHigh", "HighLow", "HighMed", "HighHigh")
Ktitle1 <- ggdraw() + draw_label(expression(paste(tau[xi], "=0.1")), fontface="bold", x=0.53, size=22)
Krow1 <-  plot_grid(Ktitle1, plot_grid(Kplot$LowLow, Kplot$LowMed, Kplot$LowHigh, nrow=1), ncol=1, rel_heights=c(0.5,1))
Ktitle2 <- ggdraw() + draw_label(expression(paste(tau[xi], "=0.9")), fontface="bold", x=0.53, size=22)
Krow2 <- plot_grid(Ktitle2, plot_grid(Kplot$HighLow, Kplot$HighMed, Kplot$HighHigh, nrow=1), ncol=1, rel_heights=c(0.5, 1))
impulseK <- plot_grid(Krow1, Krow2, nrow=2)
save_plot("/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Figures/Main/impulseK.png", impulseK, base_height = 10, base_width = 13)
##############################################################
#For Plot.ly
###############################################################
#Productivity
annotationsW <- list(list(x=0.11, y=0.9, text=TeX("(a) \\, \\tau_{\\xi}=0.1, \\tau_{\\omega_{1}}=0.1"), font=list(size=16), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
	list(x=0.495, y=0.9, text=TeX("(b)\\,\\tau_{\\xi}=0.1, \\tau_{\\omega_{1}}=0.5"), font=list(size=16), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.88, y=0.9, text=TeX("(c)\\,\\tau_{\\xi}=0.1, \\tau_{\\omega_{1}}=0.9"), font=list(size=16), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.11, y=0.45, text=TeX("(d)\\,\\tau_{\\xi}=0.9, \\tau_{\\omega_{1}}=0.1"), font=list(size=16), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.495, y=0.45, text=TeX("(e)\\,\\tau_{\\xi}=0.9, \\tau_{\\omega_{1}}=0.5"), font=list(size=16), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.88, y=0.45, text=TeX("(f)\\,\\tau_{\\xi}=0.9, \\tau_{\\omega_{1}}=0.9"), font=list(size=16), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE))
W <- subplot(Wplotly[[1]], Wplotly[[2]], Wplotly[[3]], Wplotly[[4]], Wplotly[[5]], Wplotly[[6]], shareX=TRUE, shareY=TRUE, nrows=2)
W <- W %>% layout(annotations=annotationsW) %>% config(mathjax = 'cdn')
W
W <- plotly_json(W, FALSE)
write(W, "/Users/justindoty/Documents/Home/My_Website/static/jmp/main/impulseW.json")
#Labor 
annotationsL <- list(list(x=0.11, y=0.9, text=TeX("(a) \\, \\tau_{\\xi}=0.1, \\tau_{l}=0.1"), font=list(size=16), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
	list(x=0.495, y=0.9, text=TeX("(b)\\,\\tau_{\\xi}=0.1, \\tau_{l}=0.5"), font=list(size=16), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.88, y=0.9, text=TeX("(c)\\,\\tau_{\\xi}=0.1, \\tau_{l}=0.9"), font=list(size=16), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.11, y=0.45, text=TeX("(d)\\,\\tau_{\\xi}=0.9, \\tau_{l}=0.1"), font=list(size=16), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.495, y=0.45, text=TeX("(e)\\,\\tau_{\\xi}=0.9, \\tau_{l}=0.5"), font=list(size=16), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.88, y=0.45, text=TeX("(f)\\,\\tau_{\\xi}=0.9, \\tau_{l}=0.9"), font=list(size=16), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE))
L <- subplot(Lplotly[[1]], Lplotly[[2]], Lplotly[[3]], Lplotly[[4]], Lplotly[[5]], Lplotly[[6]], shareX=TRUE, shareY=TRUE, nrows=2)
L <- L %>% layout(annotations=annotationsL) %>% config(mathjax = 'cdn')
L
L <- plotly_json(L, FALSE)
write(L, "/Users/justindoty/Documents/Home/My_Website/static/jmp/main/impulseL.json")
#Materials
annotationsM <- list(list(x=0.11, y=0.9, text=TeX("(a) \\, \\tau_{\\xi}=0.1, \\tau_{m}=0.1"), font=list(size=16), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
	list(x=0.495, y=0.9, text=TeX("(b)\\,\\tau_{\\xi}=0.1, \\tau_{m}=0.5"), font=list(size=16), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.88, y=0.9, text=TeX("(c)\\,\\tau_{\\xi}=0.1, \\tau_{m}=0.9"), font=list(size=16), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.11, y=0.45, text=TeX("(d)\\,\\tau_{\\xi}=0.9, \\tau_{m}=0.1"), font=list(size=16), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.495, y=0.45, text=TeX("(e)\\,\\tau_{\\xi}=0.9, \\tau_{m}=0.5"), font=list(size=16), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.88, y=0.45, text=TeX("(f)\\,\\tau_{\\xi}=0.9, \\tau_{m}=0.9"), font=list(size=16), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE))
M <- subplot(Mplotly[[1]], Mplotly[[2]], Mplotly[[3]], Mplotly[[4]], Mplotly[[5]], Mplotly[[6]], shareX=TRUE, shareY=TRUE, nrows=2)
M <- M %>% layout(annotations=annotationsM) %>% config(mathjax = 'cdn')
M
M <- plotly_json(M, FALSE)
write(M, "/Users/justindoty/Documents/Home/My_Website/static/jmp/main/impulseM.json")
#Investment
annotationsI <- list(list(x=0.11, y=0.9, text=TeX("(a) \\, \\tau_{\\xi}=0.1, \\tau_{i}=0.1"), font=list(size=16), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
	list(x=0.495, y=0.9, text=TeX("(b)\\,\\tau_{\\xi}=0.1, \\tau_{i}=0.5"), font=list(size=16), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.88, y=0.9, text=TeX("(c)\\,\\tau_{\\xi}=0.1, \\tau_{i}=0.9"), font=list(size=16), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.11, y=0.45, text=TeX("(d)\\,\\tau_{\\xi}=0.9, \\tau_{i}=0.1"), font=list(size=16), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.495, y=0.45, text=TeX("(e)\\,\\tau_{\\xi}=0.9, \\tau_{i}=0.5"), font=list(size=16), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.88, y=0.45, text=TeX("(f)\\,\\tau_{\\xi}=0.9, \\tau_{i}=0.9"), font=list(size=16), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE))
I <- subplot(Iplotly[[1]], Iplotly[[2]], Iplotly[[3]], Iplotly[[4]], Iplotly[[5]], Iplotly[[6]], shareX=TRUE, shareY=TRUE, nrows=2)
I <- I %>% layout(annotations=annotationsI) %>% config(mathjax = 'cdn')
I
I <- plotly_json(I, FALSE)
write(I, "/Users/justindoty/Documents/Home/My_Website/static/jmp/main/impulseI.json")
#Capital
K <- subplot(Kplotly[[1]], Kplotly[[2]], Kplotly[[3]], Kplotly[[4]], Kplotly[[5]], Kplotly[[6]], shareX=TRUE, shareY=TRUE, nrows=2)
K <- K %>% layout(annotations=annotationsI) %>% config(mathjax = 'cdn')
K
K <- plotly_json(K, FALSE)
write(K, "/Users/justindoty/Documents/Home/My_Website/static/jmp/main/impulseK.json")






