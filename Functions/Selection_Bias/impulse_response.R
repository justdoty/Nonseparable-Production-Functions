require(cowplot)
require(dplyr)
require(purrr)
require(xtable)
require(stringr)
require(truncnorm)
require(RColorBrewer)
require(reshape2)
require(plotly)
source('/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Functions/Selection_Bias/Auxfuns.R')
source('/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Functions/Selection_Bias/Tensors.R')
source('/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Functions/Selection_Bias/omega.R')
#Load US Dataset
US <- read.csv('/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Data/USdata.csv') %>% 
select(id, year, lny, lnk1, lnk2, lnl, lnm, lni, age) %>% transmute(id=id, time=year, Y=lny, K=lnk1, L=lnl, M=lnm, I=lni, A=age)
#Productivity from LP
omegainit <- omega_est(idvar=US$id, timevar=US$year, Y=US$Y, A=US$A, K=US$K, L=US$L, M=US$M)$omega
US <- US %>% mutate(omega=omegainit)
#De-mean
stdY <- sd(US$Y)
stdK <- sd(US$K)
stdL <- sd(US$L)
stdM <- sd(US$M)
stdI <- sd(US$I)
#De-mean
US <- US %>% mutate(Y=(Y-mean(Y))/stdY, K=(K-mean(K))/stdK, L=(L-mean(L))/stdL, M=(M-mean(M))/stdM, I=(I-mean(I))/stdI)
########################################################################################################
##########################################Load Results############################################
########################################################################################################
set.seed(123456)
load("/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Environments/selection_bias.RData")
vectau <- results$vectau
ntau <- length(vectau)
dims <- results$dims
#Load Parameter Estimates
parY <- results$resYmat
parL <- results$resLmat
parM <- results$resMmat
parI <- results$resImat
parWT <- results$resWTmat
parBAR <- results$resBARmat
parW1 <- results$resW1mat
parYb <- results$resyb1bLmat
parLb <- results$reslb1bLmat
parMb <- results$resmb1bLmat
parIb <- results$resib1bLmat
parWTb <- results$reswtb1bLmat
parW1b <- results$resw1b1bLmat
WTminmax <- results$maxminwtmat
#Load the method used (e.g. Cobb, Translog, Hermite)
method <- results$method
resY <- results$resY
wmin <- -5
wmax <- 5
###############################################################################
#Simulate Productivity and Capital Evolution Given Model Parameters
##############################################################################
#Expand the sample
Nsim <- 1
N <- length(unique(US$id))*Nsim
T <- length(unique(US$time))
#Vector of ranks of input demand functions (small, medium, large)
tauinp <- c(0.1, 0.5, 0.9)
#Shocks to productivity
tauxi <- c(0.1, 0.5, 0.9)
#Rank of initial productivity
tauinit <- c(0.1, 0.5, 0.9)
#Output
lnydata <- array(0, c(N, T, length(tauxi), length(tauinit)))
qlnydata <- array(0, c(N, T, length(tauxi), length(tauinp), length(tauinit)))
qlnymed <- array(0, c(T, length(tauxi), length(tauinp), length(tauinit)))
qeta <- array(0, c(N, T, length(tauinp)))
eta <- matrix(runif(N*T), nrow=N, ncol=T)
#Labor###############################################################
lnldata <- array(0, c(N, T, length(tauxi), length(tauinit)))
qlnldata <- array(0, c(N, T, length(tauxi), length(tauinp), length(tauinit)))
#At the median
lnlmed <- array(0, c(T, length(tauxi), length(tauinit)))
qlnlmed <- array(0, c(T, length(tauxi), length(tauinp), length(tauinit)))
#Shocks to Labor
epsl <- matrix(runif(N*T), nrow=N, ncol=T)
qepsl <- array(0, c(N, T, length(tauinp)))
#Intermediate Input#####################################################
lnmdata <- array(0, c(N, T, length(tauxi), length(tauinit)))
qlnmdata <- array(0, c(N, T, length(tauxi), length(tauinp), length(tauinit)))
#At the median
lnmmed <- array(0, c(T, length(tauxi), length(tauinit)))
qlnmmed <- array(0, c(T, length(tauxi), length(tauinp), length(tauinit)))
#Shocks to materials
epsm <- matrix(runif(N*T), nrow=N, ncol=T)
qepsm <- array(0, c(N, T, length(tauinp)))
#Investment###############################################################
lnidata <- array(0, c(N, T, length(tauxi), length(tauinit)))
iota <- matrix(runif(N*T), nrow=N, ncol=T)
#This is for when i vary the shocks to investment demand
qlnidata <- array(0, c(N, T, length(tauxi), length(tauinp), length(tauinit)))
qlnimed <- array(0, c(T, length(tauxi), length(tauinp), length(tauinit)))
qiota <- array(0, c(N, T, length(tauinp)))
#Capital##################################################################
lnkdata <- array(0, c(N, T, length(tauxi), length(tauinit)))
#This is for when i vary the shocks to investment demand
qlnkdata <- array(0, c(N, T, length(tauxi), length(tauinp), length(tauinit)))
qlnkmed <- array(0, c(T, length(tauxi), length(tauinp), length(tauinit)))
#Productivity#############################################################
#This is used for labor, materials, and investment
omgdata <- array(0, c(N, T, length(tauxi), length(tauinit)))
omgmed <- array(0, c(T, length(tauxi), length(tauinit)))
#Innovation Shocks
xi <- matrix(runif(N*T), nrow=N, ncol=T)
xidata <- replicate(length(tauxi), xi, simplify="array")
#Misallocation
mpkmed <- array(0, c(T, length(tauxi), length(tauinit)))
mplmed <- array(0, c(T, length(tauxi), length(tauinit)))
mpmmed <- array(0, c(T, length(tauxi), length(tauinit)))
#Ranks for Input Shocks
for (q2 in 1:length(tauinp)){
	xidata[,,q2][,2] <- tauxi[q2]
	qeta[,,q2] <- array(tauinp[q2], c(N,T))
	qepsl[,,q2] <- array(tauinp[q2], c(N,T))
	qepsm[,,q2] <- array(tauinp[q2], c(N,T))
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
#Initialization
for (q1 in 1:length(tauinit)){
	#Initial Labor
	l1 <- rowSums(LX(K=k1, omega=omg1q[q1])*lspline(vectau=vectau, bvec=parL, b1=parLb[1], bL=parLb[2], u=epsl[,1]))
	l1 <- (l1>max(t1data$L))*max(t1data$L)+(l1<min(t1data$L))*min(t1data$L)+(l1<=max(t1data$L))*(l1>=min(t1data$L))*l1
	#Initial Materials
	m1 <- rowSums(MX(K=k1, L=l1, omega=omg1q[q1])*lspline(vectau=vectau, bvec=parM, b1=parMb[1], bL=parMb[2], u=epsm[,1]))
	m1 <- (m1>3*max(t1data$M))*3*max(t1data$M)+(m1<3*min(t1data$M))*3*min(t1data$M)+(m1<=3*max(t1data$M))*(m1>=3*min(t1data$M))*m1
	#Initial Investment (at median productivity)
	i1 <- rowSums(IX(K=k1, omega=omg1q[q1])*lspline(vectau=vectau, bvec=parI, b1=parIb[1], bL=parIb[2], u=iota[,1]))
	#Restrict Support of Initial Investment
	i1 <- (i1>3*max(t1data$I))*3*max(t1data$I)+(i1<3*min(t1data$I))*3*min(t1data$I)+(i1<=3*max(t1data$I))*(i1>=3*min(t1data$I))*i1
	#Output
	y1 <- rowSums(PF(K=k1, L=l1, M=m1, omega=omg1q[q1])*lspline(vectau=vectau, bvec=parY, b1=parYb[1], bL=parYb[2], u=eta[,1]))
	y1 <- (y1>3*max(t1data$Y))*3*max(t1data$Y)+(y1<3*min(t1data$Y))*3*min(t1data$Y)+(y1<=3*max(t1data$Y))*(y1>=3*min(t1data$Y))*y1
	#Initial Inputs
	for (q2 in 1:length(tauxi)){
		omgdata[,,,q2][,,q1][,1] <- omg1q[q1]
		lnkdata[,,,q2][,,q1][,1] <- k1
		lnidata[,,,q2][,,q1][,1] <- i1
		lnldata[,,,q2][,,q1][,1] <- l1
		lnmdata[,,,q2][,,q1][,1] <- m1
		lnydata[,,,q2][,,q1][,1] <- y1
		for (q3 in 1:length(tauinp)){
			ql1 <- rowSums(LX(K=k1, omega=omg1q[q1])*lspline(vectau=vectau, bvec=parL, b1=parLb[1], bL=parLb[2], u=qepsl[,,q3][,1]))
			qm1 <- rowSums(MX(K=k1, L=l1, omega=omg1q[q1])*lspline(vectau=vectau, bvec=parM, b1=parMb[1], bL=parMb[2], u=qepsm[,,q3][,1]))
			qi1 <- rowSums(IX(K=k1, omega=omg1q[q1])*lspline(vectau=vectau, bvec=parI, b1=parIb[1], bL=parIb[2], u=qiota[,,q3][,1]))
			# #Restricting the Supports
			qlnldata[,,,,q3][,,,q2][,,q1][,1] <- (ql1>3*max(t1data$L))*3*max(t1data$L)+(ql1<3*min(t1data$L))*3*min(t1data$L)+(ql1<=3*max(t1data$L))*(ql1>=3*min(t1data$L))*ql1
			qlnmdata[,,,,q3][,,,q2][,,q1][,1] <- (qm1>3*max(t1data$M))*3*max(t1data$M)+(qm1<3*min(t1data$M))*3*min(t1data$M)+(qm1<=3*max(t1data$M))*(qm1>=3*min(t1data$M))*qm1
			qlnidata[,,,,q3][,,,q2][,,q1][,1] <- (qi1>max(t1data$I))*max(t1data$I)+(qi1<min(t1data$I))*min(t1data$I)+(qi1<=max(t1data$I))*(qi1>=min(t1data$I))*qi1
			qlnkdata[,,,,q3][,,,q2][,,q1][,1] <- k1
		}
	}
	for (t in 2:T){
		ttdata <- US %>% group_by(id) %>% slice(t)
		for (q2 in 1:length(tauxi)){
			omg <- rowSums(WX(omega=omgdata[,,,q2][,,q1][,t-1])*lspline(vectau=vectau, bvec=parWT, b1=parWTb[1], bL=parWTb[2], u=xidata[,,q2][,t]))
			#Restricting the Supports
			omgdata[,,,q2][,,q1][,t] <- (omg>wmax)*wmax+(omg<wmin)*wmin+(omg<=wmax)*(omg>=wmin)*omg
			#Generate Capital According to Accumulation Process with Depreciation
			lnk <- log(0.98*exp(lnkdata[,,,q2][,,q1][,t-1])+exp(lnidata[,,,q2][,,q1][,t-1]))
			#Restricting Supports
			lnkdata[,,,q2][,,q1][,t] <- (lnk>max(ttdata$K))*max(ttdata$K)+(lnk<min(ttdata$K))*min(ttdata$K)+(lnk<=max(ttdata$K))*(lnk>=min(ttdata$K))*lnk
			#Generate Investment
			lni <- rowSums(IX(K=lnkdata[,,,q2][,,q1][,t], omega=omgdata[,,,q2][,,q1][,t])*lspline(vectau=vectau, bvec=parI, b1=parIb[1], bL=parIb[2], u=iota[,t]))
			#Restricting Supports
			lnidata[,,,q2][,,q1][,t] <- (lni>max(ttdata$I))*max(ttdata$I)+(lni<min(ttdata$I))*min(ttdata$I)+(lni<=max(ttdata$I))*(lni>=min(ttdata$I))*lni
		}
	}
	#Generate Labor, Materials, and Output
	for (t in 2:T){
		ttdata <- US %>% group_by(id) %>% slice(t)
		for (q2 in 1:length(tauxi)){
			lab  <- rowSums(LX(K=lnkdata[,,,q2][,,q1][,t], omega=omgdata[,,,q2][,,q1][,t])*lspline(vectau=vectau, bvec=parL, b1=parLb[1], bL=parLb[2], u=epsl[,t]))
			#Restricting the Supports
			lnldata[,,,q2][,,q1][,t] <- (lab>max(ttdata$L))*max(ttdata$L)+(lab<min(ttdata$L))*min(ttdata$L)+(lab<=max(ttdata$L))*(lab>=min(ttdata$L))*lab
			#Materials
			mat <- rowSums(MX(K=lnkdata[,,,q2][,,q1][,t], L=lnldata[,,,q2][,,q1][,t], omega=omgdata[,,,q2][,,q1][,t])*lspline(vectau=vectau, bvec=parM, b1=parMb[1], bL=parMb[2], u=epsm[,t]))
			#Restricting the Supports
			lnmdata[,,,q2][,,q1][,t] <- (mat>3*max(ttdata$M))*3*max(ttdata$M)+(mat<3*min(ttdata$M))*3*min(ttdata$M)+(mat<=3*max(ttdata$M))*(mat>=3*min(ttdata$M))*mat
			#Generate Output
			lny <- rowSums(PF(K=lnkdata[,,,q2][,,q1][,t], L=lnldata[,,,q2][,,q1][,t], M=lnmdata[,,,q2][,,q1][,t], omega=omgdata[,,,q2][,,q1][,t])*lspline(vectau=vectau, bvec=parY, b1=parYb[1], bL=parYb[2], u=eta[,t]))
			#Restricting the Supports
			lnydata[,,,q2][,,q1][,t] <- (lny>3*max(ttdata$Y))*3*max(ttdata$Y)+(lny<3*min(ttdata$Y))*3*min(ttdata$Y)+(lny<=3*max(ttdata$Y))*(lny>=3*min(ttdata$Y))*lny
			for (q3 in 1:length(tauinp)){
						qlab  <- rowSums(LX(K=lnkdata[,,,q2][,,q1][,t], omega=omgdata[,,,q2][,,q1][,t])*lspline(vectau=vectau, bvec=parL, b1=parLb[1], bL=parLb[2], u=qepsl[,,q3][,t]))
						#Restricting the Supports
						qlnldata[,,,,q3][,,,q2][,,q1][,t] <- (qlab>3*max(ttdata$L))*3*max(ttdata$L)+(qlab<3*min(ttdata$L))*3*min(ttdata$L)+(qlab<=3*max(ttdata$L))*(qlab>=3*min(ttdata$L))*qlab
						#Materials
						qmat <- rowSums(MX(K=lnkdata[,,,q2][,,q1][,t], L=lnldata[,,,q2][,,q1][,t], omega=omgdata[,,,q2][,,q1][,t])*lspline(vectau=vectau, bvec=parM, b1=parMb[1], bL=parMb[2], u=qepsm[,,q3][,t]))
						#Restricting the Supports
						qlnmdata[,,,,q3][,,,q2][,,q1][,t] <- (qmat>3*max(ttdata$M))*3*max(ttdata$M)+(qmat<3*min(ttdata$M))*3*min(ttdata$M)+(qmat<=3*max(ttdata$M))*(qmat>=3*min(ttdata$M))*qmat
						#Generate Output
						qlny <- rowSums(PF(K=lnkdata[,,,q2][,,q1][,t], L=lnldata[,,,q2][,,q1][,t], M=lnmdata[,,,q2][,,q1][,t], omega=omgdata[,,,q2][,,q1][,t])*lspline(vectau=vectau, bvec=parY, b1=parYb[1], bL=parYb[2], u=qeta[,,q3][,t]))
						#Restricting the Supports
						qlnydata[,,,,q3][,,,q2][,,q1][,t] <- (qlny>3*max(ttdata$Y))*3*max(ttdata$Y)+(qlny<3*min(ttdata$Y))*3*min(ttdata$Y)+(qlny<=3*max(ttdata$Y))*(qlny>=3*min(ttdata$Y))*qlny
						#Capital 
						qlnk <- log(0.98*exp(qlnkdata[,,,,q3][,,,q2][,,q1][,t-1])+exp(qlnidata[,,,,q3][,,,q2][,,q1][,t-1]))
						qlnkdata[,,,,q3][,,,q2][,,q1][,t] <- (qlnk>max(ttdata$K))*max(ttdata$K)+(qlnk<min(ttdata$K))*min(ttdata$K)+(qlnk<=max(ttdata$K))*(qlnk>=min(ttdata$K))*qlnk
						#Investment
						qlni <- rowSums(IX(K=qlnkdata[,,,,q3][,,,q2][,,q1][,t], omega=omgdata[,,,q2][,,q1][,t])*lspline(vectau=vectau, bvec=parI, b1=parIb[1], bL=parIb[2], u=qiota[,,q3][,t]))
						qlnidata[,,,,q3][,,,q2][,,q1][,t] <- (qlni>3*max(ttdata$I))*3*max(ttdata$I)+(qlni<3*min(ttdata$I))*3*min(ttdata$I)+(qlni<=3*max(ttdata$I))*(qlni>=3*min(ttdata$I))*qlni
			}
		}
	}
	for (q2 in 1:length(tauxi)){
		#Dispersion in Misallocation: Cross-sectional standard deviations
		mpkmed[,,q2][,q1] <- apply(lnydata[,,,q2][,,q1]-lnkdata[,,,q2][,,q1], 2, sd)
		mplmed[,,q2][,q1] <- apply(lnydata[,,,q2][,,q1]-lnldata[,,,q2][,,q1], 2, sd)
		mpmmed[,,q2][,q1] <- apply(lnydata[,,,q2][,,q1]-lnmdata[,,,q2][,,q1], 2, sd)
		omgmed[,,q2][,q1] <- apply(omgdata[,,,q2][,,q1], 2, median)
		for (q3 in 1:length(tauinp)){
		#Median across firms
		qlnlmed[,,,q3][,,q2][,q1] <- apply(qlnldata[,,,,q3][,,,q2][,,q1], 2, median)
		qlnmmed[,,,q3][,,q2][,q1] <- apply(qlnmdata[,,,,q3][,,,q2][,,q1], 2, median)
		qlnimed[,,,q3][,,q2][,q1] <- apply(qlnidata[,,,,q3][,,,q2][,,q1], 2, median)
		qlnkmed[,,,q3][,,q2][,q1] <- apply(qlnkdata[,,,,q3][,,,q2][,,q1], 2, median)
		qlnymed[,,,q3][,,q2][,q1] <- apply(qlnydata[,,,,q3][,,,q2][,,q1], 2, median)
		}
	}
}
#Labor Paths
labpath <- list(Time=1:T, low_xi_low_L=qlnlmed[,,,1][,,1]-qlnlmed[,,,1][,,2], low_xi_med_L=qlnlmed[,,,2][,,1]-qlnlmed[,,,2][,,2], low_xi_high_L=qlnlmed[,,,3][,,1]-qlnlmed[,,,3][,,2],
	high_xi_low_L=qlnlmed[,,,1][,,3]-qlnlmed[,,,1][,,2], high_xi_med_L=qlnlmed[,,,2][,,3]-qlnlmed[,,,2][,,2], high_xi_high_L=qlnlmed[,,,3][,,3]-qlnlmed[,,,3][,,2])
#Materials Path
matpath <- list(Time=1:T, low_xi_low_M=qlnmmed[,,,1][,,1]-qlnmmed[,,,1][,,2], low_xi_med_M=qlnmmed[,,,2][,,1]-qlnmmed[,,,2][,,2], low_xi_high_M=qlnmmed[,,,3][,,1]-qlnmmed[,,,3][,,2],
	high_xi_low_M=qlnmmed[,,,1][,,3]-qlnmmed[,,,1][,,2], high_xi_med_M=qlnmmed[,,,2][,,3]-qlnmmed[,,,2][,,2], high_xi_high_M=qlnmmed[,,,3][,,3]-qlnmmed[,,,3][,,2])
#Investment Paths
ipath <- list(Time=1:T, low_xi_low_I=qlnimed[,,,1][,,1]-qlnimed[,,,1][,,2], low_xi_med_I=qlnimed[,,,2][,,1]-qlnimed[,,,2][,,2], low_xi_high_I=qlnimed[,,,3][,,1]-qlnimed[,,,3][,,2],
	high_xi_low_I=qlnimed[,,,1][,,3]-qlnimed[,,,1][,,2], high_xi_med_I=qlnimed[,,,2][,,3]-qlnimed[,,,2][,,2], high_xi_high_I=qlnimed[,,,3][,,3]-qlnimed[,,,3][,,2])
#Output Paths
ypath <- list(Time=1:T, low_xi_low_Y=qlnymed[,,,1][,,1]-qlnymed[,,,1][,,2], low_xi_med_Y=qlnymed[,,,2][,,1]-qlnymed[,,,2][,,2], low_xi_high_Y=qlnymed[,,,3][,,1]-qlnymed[,,,3][,,2],
	high_xi_low_Y=qlnymed[,,,1][,,3]-qlnymed[,,,1][,,2], high_xi_med_Y=qlnymed[,,,2][,,3]-qlnymed[,,,2][,,2], high_xi_high_Y=qlnymed[,,,3][,,3]-qlnymed[,,,3][,,2])
#Misallocation Paths
mispath <- list(Time=1:T, low_xi_MPK=mpkmed[,,1]-mpkmed[,,2], low_xi_MPL=mplmed[,,1]-mplmed[,,2], low_xi_MPM=mpmmed[,,1]-mpmmed[,,2],
	high_xi_MPK=mpkmed[,,3]-mpkmed[,,2], high_xi_MPL=mplmed[,,3]-mplmed[,,2], high_xi_MPM=mpmmed[,,3]-mpmmed[,,2])
#Productivity Paths
omgpath <- list(Time=1:T, low_xi_low_W=omgmed[,,1][,1]-omgmed[,,2][,1], low_xi_med_W=omgmed[,,1][,2]-omgmed[,,2][,2], low_xi_high_W=omgmed[,,1][,3]-omgmed[,,2][,3],
	high_xi_low_W=omgmed[,,3][,1]-omgmed[,,2][,1], high_xi_med_W=omgmed[,,3][,2]-omgmed[,,2][,2], high_xi_high_W=omgmed[,,3][,3]-omgmed[,,2][,3])
#Save as R environments
mainpaths <- list(labpath, matpath, ipath, ypath, mispath, omgpath)
#Save as R environments
# save(mainpaths, file="/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Functions/Selection_Bias/mainpaths.Rdata")
#Plotting
hline <- function(y = 0, color = "red"){list(type = "line", y0 = y, y1 = y, xref = "paper", x0 = 0, x1 = 1, line = list(color = color, dash="dot"))}
Lplotly <- list()
Mplotly <- list()
Iplotly <- list()
Yplotly <- list()
Misplotly <- list()
Wplotly <- list()
#Labor
ltitles <- list(low_xi_low_L=list("Negative Shock<br>Low Labor<br>Low Productivity", "Negative Shock<br>Low Labor<br>Median Productivity", "Negative Shock<br>Low Labor<br>High Productivity"), 
	low_xi_med_L=list("Negative Shock<br>Median Labor<br>Low Productivity", "Negative Shock<br>Median Labor<br>Median Productivity", "Negative Shock<br>Median Labor<br>High Productivity"),
	low_xi_high_L=list("Negative Shock<br>High Labor<br>Low Productivity", "Negative Shock<br>High Labor<br>Median Productivity", "Negative Shock<br>High Labor<br>High Productivity"),
	high_xi_low_L=list("Positive Shock<br>Low Labor<br>Low Productivity", "Positive Shock<br>Low Labor<br>Median Productivity", "Positive Shock<br>Low Labor<br>High Productivity"), 
	high_xi_med_L=list("Positive Shock<br>Median Labor<br>Low Productivity", "Positive Shock<br>Median Labor<br>Median Productivity", "Positive Shock<br>Median Labor<br>High Productivity"),
	high_xi_high_L=list("Positive Shock<br>High Labor<br>Low Productivity", "Positive Shock<br>High Labor<br>Median Productivity", "Positive Shock<br>High Labor<br>High Productivity"))
#Materials
mtitles <- list(low_xi_low_M=list("Negative Shock<br>Low Materials<br>Low Productivity", "Negative Shock<br>Low Materials<br>Median Productivity", "Negative Shock<br>Low Materials<br>High Productivity"), 
	low_xi_med_M=list("Negative Shock<br>Median Materials<br>Low Productivity", "Negative Shock<br>Median Materials<br>Median Productivity", "Negative Shock<br>Median Materials<br>High Productivity"),
	low_xi_high_M=list("Negative Shock<br>High Materials<br>Low Productivity", "Negative Shock<br>High Materials<br>Median Productivity", "Negative Shock<br>High Materials<br>High Productivity"),
	high_xi_low_M=list("Positive Shock<br>Low Materials<br>Low Productivity", "Positive Shock<br>Low Materials<br>Median Productivity", "Positive Shock<br>Low Materials<br>High Productivity"), 
	high_xi_med_M=list("Positive Shock<br>Median Materials<br>Low Productivity", "Positive Shock<br>Median Materials<br>Median Productivity", "Positive Shock<br>Median Materials<br>High Productivity"),
	high_xi_high_M=list("Positive Shock<br>High Materials<br>Low Productivity", "Positive Shock<br>High Materials<br>Median Productivity", "Positive Shock<br>High Materials<br>High Productivity"))
#Capital
ititles <- list(low_xi_low_K=list("Negative Shock<br>Low Investment<br>Low Productivity", "Negative Shock<br>Low Investment<br>Median Productivity", "Negative Shock<br>Low Investment<br>High Productivity"), 
	low_xi_med_K=list("Negative Shock<br>Median Investment<br>Low Productivity", "Negative Shock<br>Median Investment<br>Median Productivity", "Negative Shock<br>Median Investment<br>High Productivity"),
	low_xi_high_K=list("Negative Shock<br>High Investment<br>Low Productivity", "Negative Shock<br>High Investment<br>Median Productivity", "Negative Shock<br>High Investment<br>High Productivity"),
	high_xi_low_K=list("Positive Shock<br>Low Investment<br>Low Productivity", "Positive Shock<br>Low Investment<br>Median Productivity", "Positive Shock<br>Low Investment<br>High Productivity"), 
	high_xi_med_K=list("Positive Shock<br>Median Investment<br>Low Productivity", "Positive Shock<br>Median Investment<br>Median Productivity", "Positive Shock<br>Median Investment<br>High Productivity"),
	high_xi_high_K=list("Positive Shock<br>High Investment<br>Low Productivity", "Positive Shock<br>High Investment<br>Median Productivity", "Positive Shock<br>High Investment<br>High Productivity"))
#Output
ytitles <- list(low_xi_low_Y=list("Negative Shock<br>Low Output<br>Low Productivity", "Negative Shock<br>Low Output<br>Median Productivity", "Negative Shock<br>Low Output<br>High Productivity"), 
	low_xi_med_Y=list("Negative Shock<br>Median Output<br>Low Productivity", "Negative Shock<br>Median Output<br>Median Productivity", "Negative Shock<br>Median Output<br>High Productivity"),
	low_xi_high_Y=list("Negative Shock<br>High Output<br>Low Productivity", "Negative Shock<br>High Output<br>Median Productivity", "Negative Shock<br>High Output<br>High Productivity"),
	high_xi_low_Y=list("Positive Shock<br>Low Output<br>Low Productivity", "Positive Shock<br>Low Output<br>Median Productivity", "Positive Shock<br>Low Output<br>High Productivity"), 
	high_xi_med_Y=list("Positive Shock<br>Median Output<br>Low Productivity", "Positive Shock<br>Median Output<br>Median Productivity", "Positive Shock<br>Median Output<br>High Productivity"),
	high_xi_high_Y=list("Positive Shock<br>High Output<br>Low Productivity", "Positive Shock<br>High Output<br>Median Productivity", "Positive Shock<br>High Output<br>High Productivity"))
#Misallocation
mistitles <- list(low_xi_MPK=list("Negative Shock<br>MPK<br>Low Productivity", "Negative Shock<br>MPK<br>Median Productivity", "Negative Shock<br>MPK<br>High Productivity"), 
	low_xi_MPL=list("Negative Shock<br>MPL<br>Low Productivity", "Negative Shock<br>MPL<br>Median Productivity", "Negative Shock<br>MPL<br>High Productivity"),
	low_xi_MPM=list("Negative Shock<br>MPM<br>Low Productivity", "Negative Shock<br>MPM<br>Median Productivity", "Negative Shock<br>MPM<br>High Productivity"),
	high_xi_MPK=list("Positive Shock<br>MPK<br>Low Productivity", "Positive Shock<br>MPK<br>Median Productivity", "Positive Shock<br>MPK<br>High Productivity"), 
	high_xi_MPL=list("Positive Shock<br>MPL<br>Low Productivity", "Positive Shock<br>MPL<br>Median Productivity", "Positive Shock<br>MPL<br>High Productivity"),
	high_xi_MPM=list("Positive Shock<br>MPM<br>Low Productivity", "Positive Shock<br>MPM<br>Median Productivity", "Positive Shock<br>MPM<br>High Productivity"))
#Productivity
wtitles <- list("Negative Shock<br>Low Productivity", "Negative Shock<br>Medium Productivity", "Negative Shock<br>High Productivity", "Positive Shock<br>Low Productivity", "Positive Shock<br>Medium Productivity", "Positive Shock<br>High Productivity")
wcolors <- list("red", "green", "blue", "red", "green", "blue")
for (i in 1:6){
	#Labor
	ldat <- data.frame(Time=labpath[[1]], low_W=labpath[[i+1]][,1], med_W=labpath[[i+1]][,2], high_W=labpath[[i+1]][,3])
	Lplotly[[i]] <- plot_ly(ldat, x=~Time, y=~low_W, type = 'scatter', mode = 'lines', showlegend=F, line=list(color="red"), name=ltitles[[i]][[1]], hovertemplate = paste("<i>Year<i>: %{x}", "<br>Labor Change: %{y:.3f}")) %>%  add_trace(y=~med_W, line=list(color="green"), name=ltitles[[i]][[2]]) %>% add_trace(y=~high_W, line=list(color="blue"), name=ltitles[[i]][[3]]) %>% layout(xaxis=list(title="Years", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="Labor", titlefont=list(size=18), tickfont=list(size=14), range=list(min(unlist(lapply(labpath[names(labpath)!="Time"], min))), max(unlist(lapply(labpath[names(labpath)!="Time"], max))))), shapes=list(hline(y=0)))
	#Materials
	mdat <- data.frame(Time=matpath[[1]], low_W=matpath[[i+1]][,1], med_W=matpath[[i+1]][,2], high_W=matpath[[i+1]][,3])
	Mplotly[[i]] <- plot_ly(mdat, x=~Time, y=~low_W, type = 'scatter', mode = 'lines', showlegend=F, line=list(color="red"), name=mtitles[[i]][[1]], hovertemplate = paste("<i>Year<i>: %{x}", "<br>Materials Change: %{y:.3f}")) %>%  add_trace(y=~med_W, line=list(color="green"), name=mtitles[[i]][[2]]) %>% add_trace(y=~high_W, line=list(color="blue"), name=mtitles[[i]][[3]]) %>% layout(xaxis=list(title="Years", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="Materials", titlefont=list(size=18), tickfont=list(size=14), range=list(min(unlist(lapply(matpath[names(matpath)!="Time"], min))), max(unlist(lapply(matpath[names(matpath)!="Time"], max))))), shapes=list(hline(y=0)))
	#Capital
	idat <- data.frame(Time=ipath[[1]], low_W=ipath[[i+1]][,1], med_W=ipath[[i+1]][,2], high_W=ipath[[i+1]][,3])
	Iplotly[[i]] <- plot_ly(idat, x=~Time, y=~low_W, type = 'scatter', mode = 'lines', showlegend=F, line=list(color="red"), name=ititles[[i]][[1]], hovertemplate = paste("<i>Year<i>: %{x}", "<br>Investment Change: %{y:.3f}")) %>%  add_trace(y=~med_W, line=list(color="green"), name=ititles[[i]][[2]]) %>% add_trace(y=~high_W, line=list(color="blue"), name=ititles[[i]][[3]]) %>% layout(xaxis=list(title="Years", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="Investment", titlefont=list(size=18), tickfont=list(size=14), range=list(min(unlist(lapply(ipath[names(ipath)!="Time"], min))), max(unlist(lapply(ipath[names(ipath)!="Time"], max))))), shapes=list(hline(y=0)))
	#Output
	ydat <- data.frame(Time=ypath[[1]], low_W=ypath[[i+1]][,1], med_W=ypath[[i+1]][,2], high_W=ypath[[i+1]][,3])
	Yplotly[[i]] <- plot_ly(ydat, x=~Time, y=~low_W, type = 'scatter', mode = 'lines', showlegend=F, line=list(color="red"), name=ytitles[[i]][[1]], hovertemplate = paste("<i>Year<i>: %{x}", "<br>Output Change: %{y:.3f}")) %>%  add_trace(y=~med_W, line=list(color="green"), name=ytitles[[i]][[2]]) %>% add_trace(y=~high_W, line=list(color="blue"), name=ytitles[[i]][[3]]) %>% layout(xaxis=list(title="Years", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="Output", titlefont=list(size=18), tickfont=list(size=14), range=list(min(unlist(lapply(ypath[names(ypath)!="Time"], min))), max(unlist(lapply(ypath[names(ypath)!="Time"], max))))), shapes=list(hline(y=0)))
	#Misallocation
	misdat <- data.frame(Time=mispath[[1]], low_W=mispath[[i+1]][,1], med_W=mispath[[i+1]][,2], high_W=mispath[[i+1]][,3])
	Misplotly[[i]] <- plot_ly(misdat, x=~Time, y=~low_W, type = 'scatter', mode = 'lines', showlegend=F, line=list(color="red"), name=mistitles[[i]][[1]], hovertemplate = paste("<i>Year<i>: %{x}", "<br>Misallocation Change: %{y:.3f}")) %>%  add_trace(y=~med_W, line=list(color="green"), name=mistitles[[i]][[2]]) %>% add_trace(y=~high_W, line=list(color="blue"), name=mistitles[[i]][[3]]) %>% layout(xaxis=list(title="Years", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="Misallocation", titlefont=list(size=18), tickfont=list(size=14), range=list(min(unlist(lapply(mispath[names(mispath)!="Time"], min))), max(unlist(lapply(mispath[names(mispath)!="Time"], max))))), shapes=list(hline(y=0)))
	#Productivity
	wdat <- data.frame(Time=omgpath[[1]], Y=omgpath[[i+1]])
	Wplotly[[i]] <- plot_ly(wdat, x=~Time, y=~Y, type = 'scatter', mode = 'lines', showlegend=F, line=list(color=wcolors[[i]]), name=wtitles[[i]], hovertemplate = paste("<i>Year<i>: %{x}", "<br>Productivity Change: %{y:.3f}")) %>% layout(xaxis=list(title="Years", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="Productivity", titlefont=list(size=18), tickfont=list(size=14), range=list(min(unlist(lapply(omgpath[names(omgpath)!="Time"], min))), max(unlist(lapply(omgpath[names(omgpath)!="Time"], max))))), shapes=list(hline(y=0)))
}
##############################################################
#For Plot.ly
###############################################################
#Productivity
annotationsW <- list(list(x=0.11, y=0.9, text=TeX("\\boldsymbol{(a) \\, \\tau_{\\xi}=0.1, \\tau_{\\omega_{1}}=0.1}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
	list(x=0.495, y=0.9, text=TeX("\\boldsymbol{(b)\\,\\tau_{\\xi}=0.1, \\tau_{\\omega_{1}}=0.5}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.88, y=0.9, text=TeX("\\boldsymbol{(c)\\,\\tau_{\\xi}=0.1, \\tau_{\\omega_{1}}=0.9}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.11, y=0.45, text=TeX("\\boldsymbol{(d)\\,\\tau_{\\xi}=0.9, \\tau_{\\omega_{1}}=0.1}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.495, y=0.45, text=TeX("\\boldsymbol{(e)\\,\\tau_{\\xi}=0.9, \\tau_{\\omega_{1}}=0.5}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.88, y=0.45, text=TeX("\\boldsymbol{(f)\\,\\tau_{\\xi}=0.9, \\tau_{\\omega_{1}}=0.9}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE))
W <- subplot(Wplotly[[1]], Wplotly[[2]], Wplotly[[3]], Wplotly[[4]], Wplotly[[5]], Wplotly[[6]], shareX=TRUE, shareY=TRUE, nrows=2)
Wplot <- W %>% layout(annotations=annotationsW) %>% config(mathjax = 'cdn')
Wplot
Wjson <- plotly_json(W, FALSE)
write(Wjson, "/Users/justindoty/Documents/Home/My_Website/static/jmp/selection/impulseW.json")
#Labor 
annotationsL <- list(list(x=0.11, y=0.99, text=TeX("\\boldsymbol{(a) \\, \\tau_{\\xi}=0.1, \\tau_{l}=0.1}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
	list(x=0.495, y=0.99, text=TeX("\\boldsymbol{(b)\\,\\tau_{\\xi}=0.1, \\tau_{l}=0.5}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.88, y=0.99, text=TeX("\\boldsymbol{(c)\\,\\tau_{\\xi}=0.1, \\tau_{l}=0.9}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.11, y=0.45, text=TeX("\\boldsymbol{(d)\\,\\tau_{\\xi}=0.9, \\tau_{l}=0.1}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.495, y=0.45, text=TeX("\\boldsymbol{(e)\\,\\tau_{\\xi}=0.9, \\tau_{l}=0.5}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.88, y=0.45, text=TeX("\\boldsymbol{(f)\\,\\tau_{\\xi}=0.9, \\tau_{l}=0.9}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE))
L <- subplot(Lplotly[[1]], Lplotly[[2]], Lplotly[[3]], Lplotly[[4]], Lplotly[[5]], Lplotly[[6]], shareX=TRUE, shareY=TRUE, nrows=2)
Lplot <- L %>% layout(annotations=annotationsL) %>% config(mathjax = 'cdn')
Lplot
Ljson <- plotly_json(L, FALSE)
write(Ljson, "/Users/justindoty/Documents/Home/My_Website/static/jmp/selection/impulseL.json")
#Materials
annotationsM <- list(list(x=0.11, y=0.9, text=TeX("\\boldsymbol{(a) \\, \\tau_{\\xi}=0.1, \\tau_{m}=0.1}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
	list(x=0.495, y=0.9, text=TeX("\\boldsymbol{(b)\\,\\tau_{\\xi}=0.1, \\tau_{m}=0.5}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.88, y=0.9, text=TeX("\\boldsymbol{(c)\\,\\tau_{\\xi}=0.1, \\tau_{m}=0.9}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.11, y=0.45, text=TeX("\\boldsymbol{(d)\\,\\tau_{\\xi}=0.9, \\tau_{m}=0.1}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.495, y=0.45, text=TeX("\\boldsymbol{(e)\\,\\tau_{\\xi}=0.9, \\tau_{m}=0.5}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.88, y=0.45, text=TeX("\\boldsymbol{(f)\\,\\tau_{\\xi}=0.9, \\tau_{m}=0.9}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE))
M <- subplot(Mplotly[[1]], Mplotly[[2]], Mplotly[[3]], Mplotly[[4]], Mplotly[[5]], Mplotly[[6]], shareX=TRUE, shareY=TRUE, nrows=2)
Mplot <- M %>% layout(annotations=annotationsM) %>% config(mathjax = 'cdn')
Mplot
Mjson <- plotly_json(M, FALSE)
write(Mjson, "/Users/justindoty/Documents/Home/My_Website/static/jmp/selection/impulseM.json")
#Investment
annotationsI <- list(list(x=0.13, y=0.9, text=TeX("\\boldsymbol{(a) \\, \\tau_{\\xi}=0.1, \\tau_{i}=0.1}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
	list(x=0.495, y=0.9, text=TeX("\\boldsymbol{(b)\\,\\tau_{\\xi}=0.1, \\tau_{i}=0.5}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.88, y=0.9, text=TeX("\\boldsymbol{(c)\\,\\tau_{\\xi}=0.1, \\tau_{i}=0.9}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.13, y=0.46, text=TeX("\\boldsymbol{(d)\\,\\tau_{\\xi}=0.9, \\tau_{i}=0.1}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.495, y=0.46, text=TeX("\\boldsymbol{(e)\\,\\tau_{\\xi}=0.9, \\tau_{i}=0.5}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.88, y=0.46, text=TeX("\\boldsymbol{(f)\\,\\tau_{\\xi}=0.9, \\tau_{i}=0.9}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE))
#Investment
I <- subplot(Iplotly[[1]], Iplotly[[2]], Iplotly[[3]], Iplotly[[4]], Iplotly[[5]], Iplotly[[6]], shareX=TRUE, shareY=TRUE, nrows=2)
Iplot <- I %>% layout(annotations=annotationsI) %>% config(mathjax = 'cdn')
Iplot
Ijson <- plotly_json(I, FALSE)
write(Ijson, "/Users/justindoty/Documents/Home/My_Website/static/jmp/selection/impulseI.json")
annotationsY <- list(list(x=0.13, y=0.9, text=TeX("\\boldsymbol{(a) \\, \\tau_{\\xi}=0.1, \\tau_{y}=0.1}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
	list(x=0.495, y=0.9, text=TeX("\\boldsymbol{(b)\\,\\tau_{\\xi}=0.1, \\tau_{y}=0.5}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.88, y=0.9, text=TeX("\\boldsymbol{(c)\\,\\tau_{\\xi}=0.1, \\tau_{y}=0.9}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.13, y=0.46, text=TeX("\\boldsymbol{(d)\\,\\tau_{\\xi}=0.9, \\tau_{y}=0.1}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.495, y=0.46, text=TeX("\\boldsymbol{(e)\\,\\tau_{\\xi}=0.9, \\tau_{y}=0.5}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.88, y=0.46, text=TeX("\\boldsymbol{(f)\\,\\tau_{\\xi}=0.9, \\tau_{y}=0.9}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE))
#Output
Y <- subplot(Yplotly[[1]], Yplotly[[2]], Yplotly[[3]], Yplotly[[4]], Yplotly[[5]], Yplotly[[6]], shareX=TRUE, shareY=TRUE, nrows=2)
Yplot <- Y %>% layout(annotations=annotationsY) %>% config(mathjax = 'cdn')
Yplot
Yjson <- plotly_json(Y, FALSE)
write(Yjson, "/Users/justindoty/Documents/Home/My_Website/static/jmp/selection/impulseY.json")
annotationsMis <- list(list(x=0.13, y=0.9, text=TeX("\\boldsymbol{(a) \\, MP_{k}, \\tau_{\\xi}=0.1}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
	list(x=0.495, y=0.9, text=TeX("\\boldsymbol{(b)\\, MP_{l}, \\tau_{\\xi}=0.1}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.88, y=0.9, text=TeX("\\boldsymbol{(c)\\, MP_{m}, \\tau_{\\xi}=0.1}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.13, y=0.46, text=TeX("\\boldsymbol{(d)\\, MP_{k}, \\tau_{\\xi}=0.9}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.495, y=0.46, text=TeX("\\boldsymbol{(e)\\, MP_{l}, \\tau_{\\xi}=0.9}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.88, y=0.46, text=TeX("\\boldsymbol{(f)\\, MP_{m}, \\tau_{\\xi}=0.9}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE))
#Output
MIS <- subplot(Misplotly[[1]], Misplotly[[2]], Misplotly[[3]], Misplotly[[4]], Misplotly[[5]], Misplotly[[6]], shareX=TRUE, shareY=TRUE, nrows=2)
Misplot <- MIS %>% layout(annotations=annotationsMis) %>% config(mathjax = 'cdn')
Misplot
Misjson <- plotly_json(MIS, FALSE)
write(Misjson, "/Users/justindoty/Documents/Home/My_Website/static/jmp/selection/impulseMis.json")





