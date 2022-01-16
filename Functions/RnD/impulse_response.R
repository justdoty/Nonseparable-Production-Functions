require(cowplot)
require(dplyr)
require(purrr)
require(xtable)
require(stringr)
require(truncnorm)
require(RColorBrewer)
require(reshape2)
require(plotly)
source('/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Functions/RnD/Auxfuns.R')
source('/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Functions/RnD/Tensors.R')
source('/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Functions/RnD/omega.R')
#Load US Dataset
US <- read.csv('/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Data/USdata.csv') %>% 
	mutate(rdB=ifelse(rd>0, 1, 0)) %>% group_by(id) %>% filter(all(rdB==1)|all(rdB==0)) %>% ungroup() %>%
	select(id, year, lny, lnk1, lnk2, lnl, lnm, lni, rd, rdB) %>% transmute(id=id, time=year, Y=lny, K=lnk1, L=lnl, M=lnm, I=lni, R=rd, RB=rdB)
#Productivity from LP
omegainit <- omega_est(idvar=US$id, timevar=US$year, Y=US$Y, A=US$A, K=US$K, L=US$L, M=US$M)$omega
#Add LP Estimate of productivity
US <- US %>% mutate(omega=omegainit) 
idcon <- duplicated(US$id)
idlag <- duplicated(US$id, fromLast=TRUE)
id1 <- !idcon
#Percentage of RnD Firms
rfrac <- sum((US$R[duplicated(US$id)]>0))/sum(duplicated(US$id))
########################################################################################################
#Load Results
########################################################################################################
set.seed(123456)
load("/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Environments/rnd.RData")
vectau <- results$vectau
ntau <- length(vectau)
dims <- results$dims
#Load Parameter Estimates
parY <- results$resYmat
parL <- results$resLmat
parM <- results$resMmat
parI <- results$resImat
parWT <- results$resWTmat
parR <- results$resRmat
parW1 <- results$resW1mat
parYb <- results$resyb1bLmat
parLb <- results$reslb1bLmat
parMb <- results$resmb1bLmat
parIb <- results$resib1bLmat
parWTb <- results$reswtb1bLmat
parRb <- results$resrb1bLmat
parW1b <- results$resw1b1bLmat
WTminmax <- results$maxminwtmat
#De-mean
stdY <- sd(US$Y)
stdK <- sd(US$K)
stdL <- sd(US$L)
stdM <- sd(US$M)
stdI <- sd(US$I)
#De-mean
US <- US %>% mutate(Y=(Y-mean(Y))/stdY, K=(K-mean(K))/stdK, L=(L-mean(L))/stdL, M=(M-mean(M))/stdM, I=(I-mean(I))/stdI)
lnr <- log(US$R[US$R>0])
lnr <- lnr-mean(lnr)
US$R[US$R>0] <- lnr
wmin <- -5
wmax <- 5
#Load the method used (e.g. Cobb, Translog, Hermite)
method <- results$method
resY <- results$resY
###############################################################################
#Simulate Productivity and Capital Evolution Given Model Parameters
#Expand the sample
NR <- 1
NNR <- 1
#Simulate Data at T=1
t1data <- US %>% group_by(id) %>% slice(1)
#Capital is not estimated in the main model: Use Capital Accumulation Process
#For K=1
#Capital for RnD Firms
rlnk1 <- kronecker(array(1, c(NR,1)), t1data$K[t1data$RB==1])
nr <- length(rlnk1)
#Capital for Non RnD Firms
lnk1 <- kronecker(array(1, c(NNR,1)), t1data$K[t1data$RB==0])
nnr <- length(lnk1)
#Matrices
N <- nr+nnr
T <- length(unique(US$time))
#Vector of ranks of input demand functions (small, medium, large)
tauinp <- c(0.1, 0.5, 0.9)
#Shocks to productivity
tauxi <- c(0.1, 0.5, 0.9)
#Rank of initial productivity
tauinit <- c(0.1, 0.5, 0.9)
#Position in vector for RnD Firms
pos <- c((nnr+1):N)
Nsim <- 1
k1 <- kronecker(array(1, c(Nsim,1)), t1data$K)
############################################################################
#Initialize matrices 
############################################################################
#Output##################################################################
lnydata <- array(0, c(N, T, length(tauinit), length(tauxi)))
qlnydata <- array(0, c(N, T, length(tauinit), length(tauxi), length(tauinp)))
#At the median
qlnymed <- array(0, c(T, length(tauinit), length(tauxi), length(tauinp)))
rqlnymed <- qlnymed
#Shocks to Output
qeta <- array(0, c(N, T, length(tauinp)))
eta <- matrix(runif(N*T), nrow=N, ncol=T)
#Labor###############################################################
lnldata <- array(0, c(N, T, length(tauinit), length(tauxi)))
qlnldata <- array(0, c(N, T, length(tauinit), length(tauxi), length(tauinp)))
#At the median
qlnlmed <- array(0, c(T, length(tauinit), length(tauxi), length(tauinp)))
rqlnlmed <- qlnlmed
#Shocks to Labor
epsl <- matrix(runif(N*T), nrow=N, ncol=T)
qepsl <- array(0, c(N, T, length(tauinp)))
#Materials#####################################################
lnmdata <- array(0, c(N, T, length(tauinit), length(tauxi)))
qlnmdata <- array(0, c(N, T, length(tauinit), length(tauxi), length(tauinp)))
#At the median
qlnmmed <- array(0, c(T, length(tauinit), length(tauxi), length(tauinp)))
rqlnmmed <- qlnmmed
#Shocks to materials
epsm <- matrix(runif(N*T), nrow=N, ncol=T)
qepsm <- array(0, c(N, T, length(tauinp)))
#Investment###############################################################
lnidata <- array(0, c(N, T, length(tauinit), length(tauxi)))
iota <- matrix(runif(N*T), nrow=N, ncol=T)
#This is for when i vary the shocks to investment demand
qlnidata <- array(0, c(N, T, length(tauinit), length(tauxi), length(tauinp)))
qlnimed <- array(0, c(T, length(tauinit), length(tauxi), length(tauinp)))
rqlnimed <- qlnimed
qiota <- array(0, c(N, T, length(tauinp)))
#Capital##################################################################
lnkdata <- array(0, c(N, T, length(tauinit), length(tauxi)))
#This is for when i vary the shocks to investment demand
qlnkdata <- array(0, c(N, T, length(tauinit), length(tauxi), length(tauinp)))
qlnkmed <- array(0, c(T, length(tauinit), length(tauxi), length(tauinp)))
rqlnkmed <- qlnkmed
#R&D#######################################################################
rdata <- array(0, c(N, T, length(tauinit), length(tauxi)))
rho <- matrix(runif(N*T), nrow=N, ncol=T)
#Productivity#############################################################
omgdata <- array(0, c(N, T, length(tauinit), length(tauxi)))
omgmed <- array(0, c(T, length(tauinit), length(tauxi)))
romgmed <- omgmed
#Innovation Shocks
xi <- matrix(runif(N*T), nrow=N, ncol=T)
xidata <- replicate(length(tauxi), xi, simplify="array")
#Misallocation
mpkmed <- array(0, c(T, length(tauinit), length(tauxi)))
rmpkmed <- mpkmed
mplmed <- array(0, c(T, length(tauinit), length(tauxi)))
rmplmed <- mplmed
mpmmed <- array(0, c(T, length(tauinit), length(tauxi)))
rmpmmed <- mpmmed
#Productivity Shock
for (q2 in 1:length(tauxi)){
	xidata[,,q2][,2] <- tauxi[q2]
}
#Ranks for Input Shocks
for (q3 in 1:length(tauinp)){
	qeta[,,q3] <- array(tauinp[q3], c(N,T))
	qepsl[,,q3] <- array(tauinp[q3], c(N,T))
	qepsm[,,q3] <- array(tauinp[q3], c(N,T))
	qiota[,,q3] <- array(tauinp[q3], c(N,T))
}
#I order Non RnD Firms First
k1 <- c(lnk1, rlnk1)
# Nsim <- 1
# k1 <- kronecker(array(1, c(Nsim,1)), t1data$K)
#Initial Productivity
omg1 <- rowSums(WX1(K=k1)*lspline(vectau=vectau, bvec=parW1, b1=parW1b[1], bL=parW1b[2], u=xi[,1]))
#Restrict the Support of Initial Productivity
omg1 <- (omg1>wmax)*wmax+(omg1<wmin)*wmin+(omg1<=wmax)*(omg1>=wmin)*omg1
omg1q <- quantile(omg1, tauinit)

for (q1 in 1:length(tauinit)){
	#Initial Labor
	l1 <- rowSums(LX(K=k1, omega=omg1q[q1])*lspline(vectau=vectau, bvec=parL, b1=parLb[1], bL=parLb[2], u=epsl[,1]))
	l1 <- (l1>max(t1data$L))*max(t1data$L)+(l1<min(t1data$L))*min(t1data$L)+(l1<=max(t1data$L))*(l1>=min(t1data$L))*l1
	#Initial Materials
	m1 <- rowSums(MX(K=k1, L=l1, omega=omg1q[q1])*lspline(vectau=vectau, bvec=parM, b1=parMb[1], bL=parMb[2], u=epsm[,1]))
	m1 <- (m1>max(t1data$M))*max(t1data$M)+(m1<min(t1data$M))*min(t1data$M)+(m1<=max(t1data$M))*(m1>=min(t1data$M))*m1
	#Initial Investment (at median productivity)
	i1 <- rowSums(IX(K=k1, omega=omg1q[q1])*lspline(vectau=vectau, bvec=parI, b1=parIb[1], bL=parIb[2], u=iota[,1]))
	#Restrict Support of Initial Investment
	i1 <- (i1>max(t1data$I))*max(t1data$I)+(i1<min(t1data$I))*min(t1data$I)+(i1<=max(t1data$I))*(i1>=min(t1data$I))*i1
	#Output
	y1 <- rowSums(PF(K=k1, L=l1, M=m1, omega=omg1q[q1])*lspline(vectau=vectau, bvec=parY, b1=parYb[1], bL=parYb[2], u=eta[,1]))
	y1 <- (y1>max(t1data$Y))*max(t1data$Y)+(y1<min(t1data$Y))*min(t1data$Y)+(y1<=max(t1data$Y))*(y1>=min(t1data$Y))*y1
	#Initial R&D
	r1 <- rowSums(RX(K=k1[pos], omega=omg1q[q1])*lspline(vectau=vectau, bvec=parR, b1=parRb[1], bL=parRb[2], u=rho[pos,1]))
	#Restrict Support of Initial R&D
	r1 <- (r1>max(t1data$R))*max(t1data$R)+(r1<min(t1data$R[t1data$RB==1]))*min(t1data$R[t1data$RB==1])+(r1<=max(t1data$R))*(r1>=min(t1data$R[t1data$RB==1]))*r1
	#Initial Inputs
	for (q2 in 1:length(tauxi)){
		omgdata[,,,q2][,,q1][,1] <- omg1q[q1]
		lnkdata[,,,q2][,,q1][,1] <- k1
		lnidata[,,,q2][,,q1][,1] <- i1
		lnldata[,,,q2][,,q1][,1] <- l1
		lnmdata[,,,q2][,,q1][,1] <- m1
		lnydata[,,,q2][,,q1][,1] <- y1
		rdata[,,,q2][,,q1][pos,1] <- r1
		for (q3 in 1:length(tauinp)){
			ql1 <- rowSums(LX(K=k1, omega=omg1q[q1])*lspline(vectau=vectau, bvec=parL, b1=parLb[1], bL=parLb[2], u=qepsl[,,q3][,1]))
			qm1 <- rowSums(MX(K=k1, L=l1, omega=omg1q[q1])*lspline(vectau=vectau, bvec=parM, b1=parMb[1], bL=parMb[2], u=qepsm[,,q3][,1]))
			qi1 <- rowSums(IX(K=k1, omega=omg1q[q1])*lspline(vectau=vectau, bvec=parI, b1=parIb[1], bL=parIb[2], u=qiota[,,q3][,1]))
			# #Restricting the Supports
			qlnldata[,,,,q3][,,,q2][,,q1][,1] <- (ql1>max(t1data$L))*max(t1data$L)+(ql1<min(t1data$L))*min(t1data$L)+(ql1<=max(t1data$L))*(ql1>=min(t1data$L))*ql1
			qlnmdata[,,,,q3][,,,q2][,,q1][,1] <- (qm1>max(t1data$M))*max(t1data$M)+(qm1<min(t1data$M))*min(t1data$M)+(qm1<=max(t1data$M))*(qm1>=min(t1data$M))*qm1
			qlnidata[,,,,q3][,,,q2][,,q1][,1] <- (qi1>max(t1data$I))*max(t1data$I)+(qi1<min(t1data$I))*min(t1data$I)+(qi1<=max(t1data$I))*(qi1>=min(t1data$I))*qi1
			qlnkdata[,,,,q3][,,,q2][,,q1][,1] <- k1

		}
	}
	#Generate Productivity
	for (t in 2:T){
		ttdata <- US %>% group_by(id) %>% slice(t)
		for (q2 in 1:length(tauxi)){
			rind <- rdata[,,,q2][,,q1][,t-1]!=0
			omg <- rowSums(WX(omega=omgdata[,,,q2][,,q1][,t-1], R=rdata[,,,q2][,,q1][,t-1], Rind=rind)*lspline(vectau=vectau, bvec=parWT, b1=parWTb[1], bL=parWTb[2], u=xidata[,,q2][,t]))
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
			#R&D
			lnr <- rowSums(RX(K=lnkdata[,,,q2][,,q1][pos,t], omega=omgdata[,,,q2][,,q1][pos,t])*lspline(vectau=vectau, bvec=parR, b1=parRb[1], bL=parRb[2], u=rho[pos,t]))
			#Restricting Supports
			rdata[,,,q2][,,q1][pos,t] <- (lnr>max(ttdata$R))*max(ttdata$R)+(lnr<min(ttdata$R[ttdata$RB==1]))*min(ttdata$R[ttdata$RB==1])+(lnr<=max(ttdata$R))*(lnr>=min(ttdata$R[ttdata$RB==1]))*lnr
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
			lnmdata[,,,q2][,,q1][,t] <- (mat>max(ttdata$M))*max(ttdata$M)+(mat<min(ttdata$M))*min(ttdata$M)+(mat<=max(ttdata$M))*(mat>=min(ttdata$M))*mat
			#Generate Output
			lny <- rowSums(PF(K=lnkdata[,,,q2][,,q1][,t], L=lnldata[,,,q2][,,q1][,t], M=lnmdata[,,,q2][,,q1][,t], omega=omgdata[,,,q2][,,q1][,t])*lspline(vectau=vectau, bvec=parY, b1=parYb[1], bL=parYb[2], u=eta[,t]))
			#Restricting the Supports
			lnydata[,,,q2][,,q1][,t] <- (lny>max(ttdata$Y))*max(ttdata$Y)+(lny<min(ttdata$Y))*min(ttdata$Y)+(lny<=max(ttdata$Y))*(lny>=min(ttdata$Y))*lny
			for (q3 in 1:length(tauinp)){
						qlab  <- rowSums(LX(K=lnkdata[,,,q2][,,q1][,t], omega=omgdata[,,,q2][,,q1][,t])*lspline(vectau=vectau, bvec=parL, b1=parLb[1], bL=parLb[2], u=qepsl[,,q3][,t]))
						#Restricting the Supports
						qlnldata[,,,,q3][,,,q2][,,q1][,t] <- (qlab>max(ttdata$L))*max(ttdata$L)+(qlab<min(ttdata$L))*min(ttdata$L)+(qlab<=max(ttdata$L))*(qlab>=min(ttdata$L))*qlab
						#Materials
						qmat <- rowSums(MX(K=lnkdata[,,,q2][,,q1][,t], L=lnldata[,,,q2][,,q1][,t], omega=omgdata[,,,q2][,,q1][,t])*lspline(vectau=vectau, bvec=parM, b1=parMb[1], bL=parMb[2], u=qepsm[,,q3][,t]))
						#Restricting the Supports
						qlnmdata[,,,,q3][,,,q2][,,q1][,t] <- (qmat>max(ttdata$M))*max(ttdata$M)+(qmat<min(ttdata$M))*min(ttdata$M)+(qmat<=max(ttdata$M))*(qmat>=min(ttdata$M))*qmat
						#Generate Output
						qlny <- rowSums(PF(K=lnkdata[,,,q2][,,q1][,t], L=lnldata[,,,q2][,,q1][,t], M=lnmdata[,,,q2][,,q1][,t], omega=omgdata[,,,q2][,,q1][,t])*lspline(vectau=vectau, bvec=parY, b1=parYb[1], bL=parYb[2], u=qeta[,,q3][,t]))
						#Restricting the Supports
						qlnydata[,,,,q3][,,,q2][,,q1][,t] <- (qlny>max(ttdata$Y))*max(ttdata$Y)+(qlny<min(ttdata$Y))*min(ttdata$Y)+(qlny<=max(ttdata$Y))*(qlny>=min(ttdata$Y))*qlny
						#Capital 
						qlnk <- log(0.98*exp(qlnkdata[,,,,q3][,,,q2][,,q1][,t-1])+exp(qlnidata[,,,,q3][,,,q2][,,q1][,t-1]))
						qlnkdata[,,,,q3][,,,q2][,,q1][,t] <- qlnk
						# qlnkdata[,,,,q3][,,,q2][,,q1][,t] <- (qlnk>max(ttdata$K))*max(ttdata$K)+(qlnk<min(ttdata$K))*min(ttdata$K)+(qlnk<=max(ttdata$K))*(qlnk>=min(ttdata$K))*qlnk
						#Investment
						qlni <- rowSums(IX(K=qlnkdata[,,,,q3][,,,q2][,,q1][,t], omega=omgdata[,,,q2][,,q1][,t])*lspline(vectau=vectau, bvec=parI, b1=parIb[1], bL=parIb[2], u=qiota[,,q3][,t]))
						qlnidata[,,,,q3][,,,q2][,,q1][,t] <- (qlni>max(ttdata$I))*max(ttdata$I)+(qlni<min(ttdata$I))*min(ttdata$I)+(qlni<=max(ttdata$I))*(qlni>=min(ttdata$I))*qlni
			}
		}
	}
	for (q2 in 1:length(tauxi)){
		#Dispersion in Misallocation: Cross-sectional standard deviations
		#For Non R&D Firms
		mpkmed[,,q2][,q1] <- apply(lnydata[,,,q2][,,q1][-pos,]-lnkdata[,,,q2][,,q1][-pos,], 2, sd)
		mplmed[,,q2][,q1] <- apply(lnydata[,,,q2][,,q1][-pos,]-lnldata[,,,q2][,,q1][-pos,], 2, sd)
		mpmmed[,,q2][,q1] <- apply(lnydata[,,,q2][,,q1][-pos,]-lnmdata[,,,q2][,,q1][-pos,], 2, sd)
		omgmed[,,q2][,q1] <- apply(omgdata[,,,q2][,,q1][-pos,], 2, median)
		#For R&D Firms
		rmpkmed[,,q2][,q1] <- apply(lnydata[,,,q2][,,q1][pos,]-lnkdata[,,,q2][,,q1][pos,], 2, sd)
		rmplmed[,,q2][,q1] <- apply(lnydata[,,,q2][,,q1][pos,]-lnldata[,,,q2][,,q1][pos,], 2, sd)
		rmpmmed[,,q2][,q1] <- apply(lnydata[,,,q2][,,q1][pos,]-lnmdata[,,,q2][,,q1][pos,], 2, sd)
		romgmed[,,q2][,q1] <- apply(omgdata[,,,q2][,,q1][pos,], 2, median)
		for (q3 in 1:length(tauinp)){
		#For Non R&D Firms
		qlnlmed[,,,q3][,,q2][,q1] <- apply(qlnldata[,,,,q3][,,,q2][,,q1][-pos,], 2, median)
		qlnmmed[,,,q3][,,q2][,q1] <- apply(qlnmdata[,,,,q3][,,,q2][,,q1][-pos,], 2, median)
		qlnimed[,,,q3][,,q2][,q1] <- apply(qlnidata[,,,,q3][,,,q2][,,q1][-pos,], 2, median)
		qlnkmed[,,,q3][,,q2][,q1] <- apply(qlnkdata[,,,,q3][,,,q2][,,q1][-pos,], 2, median)
		qlnymed[,,,q3][,,q2][,q1] <- apply(qlnydata[,,,,q3][,,,q2][,,q1][-pos,], 2, median)
		#For R&D Firms
		rqlnlmed[,,,q3][,,q2][,q1] <- apply(qlnldata[,,,,q3][,,,q2][,,q1][pos,], 2, median)
		rqlnmmed[,,,q3][,,q2][,q1] <- apply(qlnmdata[,,,,q3][,,,q2][,,q1][pos,], 2, median)
		rqlnimed[,,,q3][,,q2][,q1] <- apply(qlnidata[,,,,q3][,,,q2][,,q1][pos,], 2, median)
		rqlnkmed[,,,q3][,,q2][,q1] <- apply(qlnkdata[,,,,q3][,,,q2][,,q1][pos,], 2, median)
		rqlnymed[,,,q3][,,q2][,q1] <- apply(qlnydata[,,,,q3][,,,q2][,,q1][pos,], 2, median)
		}
	}
}
#Labor Paths
labpath <- list(Time=1:T, low_xi_low_L=qlnlmed[,,,1][,,1]-qlnlmed[,,,1][,,2], low_xi_med_L=qlnlmed[,,,2][,,1]-qlnlmed[,,,2][,,2], low_xi_high_L=qlnlmed[,,,3][,,1]-qlnlmed[,,,3][,,2],
	high_xi_low_L=qlnlmed[,,,1][,,3]-qlnlmed[,,,1][,,2], high_xi_med_L=qlnlmed[,,,2][,,3]-qlnlmed[,,,2][,,2], high_xi_high_L=qlnlmed[,,,3][,,3]-qlnlmed[,,,3][,,2])
rlabpath <- list(Time=1:T, low_xi_low_L=rqlnlmed[,,,1][,,1]-rqlnlmed[,,,1][,,2], low_xi_med_L=rqlnlmed[,,,2][,,1]-rqlnlmed[,,,2][,,2], low_xi_high_L=rqlnlmed[,,,3][,,1]-rqlnlmed[,,,3][,,2],
	high_xi_low_L=rqlnlmed[,,,1][,,3]-rqlnlmed[,,,1][,,2], high_xi_med_L=rqlnlmed[,,,2][,,3]-rqlnlmed[,,,2][,,2], high_xi_high_L=rqlnlmed[,,,3][,,3]-rqlnlmed[,,,3][,,2])
#Materials Path
matpath <- list(Time=1:T, low_xi_low_M=qlnmmed[,,,1][,,1]-qlnmmed[,,,1][,,2], low_xi_med_M=qlnmmed[,,,2][,,1]-qlnmmed[,,,2][,,2], low_xi_high_M=qlnmmed[,,,3][,,1]-qlnmmed[,,,3][,,2],
	high_xi_low_M=qlnmmed[,,,1][,,3]-qlnmmed[,,,1][,,2], high_xi_med_M=qlnmmed[,,,2][,,3]-qlnmmed[,,,2][,,2], high_xi_high_M=qlnmmed[,,,3][,,3]-qlnmmed[,,,3][,,2])
rmatpath <- list(Time=1:T, low_xi_low_M=rqlnmmed[,,,1][,,1]-rqlnmmed[,,,1][,,2], low_xi_med_M=rqlnmmed[,,,2][,,1]-rqlnmmed[,,,2][,,2], low_xi_high_M=rqlnmmed[,,,3][,,1]-rqlnmmed[,,,3][,,2],
	high_xi_low_M=rqlnmmed[,,,1][,,3]-rqlnmmed[,,,1][,,2], high_xi_med_M=rqlnmmed[,,,2][,,3]-rqlnmmed[,,,2][,,2], high_xi_high_M=rqlnmmed[,,,3][,,3]-rqlnmmed[,,,3][,,2])
#Investment Paths
ipath <- list(Time=1:T, low_xi_low_I=qlnimed[,,,1][,,1]-qlnimed[,,,1][,,2], low_xi_med_I=qlnimed[,,,2][,,1]-qlnimed[,,,2][,,2], low_xi_high_I=qlnimed[,,,3][,,1]-qlnimed[,,,3][,,2],
	high_xi_low_I=qlnimed[,,,1][,,3]-qlnimed[,,,1][,,2], high_xi_med_I=qlnimed[,,,2][,,3]-qlnimed[,,,2][,,2], high_xi_high_I=qlnimed[,,,3][,,3]-qlnimed[,,,3][,,2])
ripath <- list(Time=1:T, low_xi_low_I=rqlnimed[,,,1][,,1]-rqlnimed[,,,1][,,2], low_xi_med_I=rqlnimed[,,,2][,,1]-rqlnimed[,,,2][,,2], low_xi_high_I=rqlnimed[,,,3][,,1]-rqlnimed[,,,3][,,2],
	high_xi_low_I=rqlnimed[,,,1][,,3]-rqlnimed[,,,1][,,2], high_xi_med_I=rqlnimed[,,,2][,,3]-rqlnimed[,,,2][,,2], high_xi_high_I=rqlnimed[,,,3][,,3]-rqlnimed[,,,3][,,2])
#Output Paths
ypath <- list(Time=1:T, low_xi_low_Y=qlnymed[,,,1][,,1]-qlnymed[,,,1][,,2], low_xi_med_Y=qlnymed[,,,2][,,1]-qlnymed[,,,2][,,2], low_xi_high_Y=qlnymed[,,,3][,,1]-qlnymed[,,,3][,,2],
	high_xi_low_Y=qlnymed[,,,1][,,3]-qlnymed[,,,1][,,2], high_xi_med_Y=qlnymed[,,,2][,,3]-qlnymed[,,,2][,,2], high_xi_high_Y=qlnymed[,,,3][,,3]-qlnymed[,,,3][,,2])
rypath <- list(Time=1:T, low_xi_low_Y=rqlnymed[,,,1][,,1]-rqlnymed[,,,1][,,2], low_xi_med_Y=rqlnymed[,,,2][,,1]-rqlnymed[,,,2][,,2], low_xi_high_Y=rqlnymed[,,,3][,,1]-rqlnymed[,,,3][,,2],
	high_xi_low_Y=rqlnymed[,,,1][,,3]-rqlnymed[,,,1][,,2], high_xi_med_Y=rqlnymed[,,,2][,,3]-rqlnymed[,,,2][,,2], high_xi_high_Y=rqlnymed[,,,3][,,3]-rqlnymed[,,,3][,,2])
#Misallocation Paths
mispath <- list(Time=1:T, low_xi_MPK=mpkmed[,,1]-mpkmed[,,2], low_xi_MPL=mplmed[,,1]-mplmed[,,2], low_xi_MPM=mpmmed[,,1]-mpmmed[,,2],
	high_xi_MPK=mpkmed[,,3]-mpkmed[,,2], high_xi_MPL=mplmed[,,3]-mplmed[,,2], high_xi_MPM=mpmmed[,,3]-mpmmed[,,2])
rmispath <- list(Time=1:T, low_xi_MPK=rmpkmed[,,1]-rmpkmed[,,2], low_xi_MPL=rmplmed[,,1]-rmplmed[,,2], low_xi_MPM=rmpmmed[,,1]-rmpmmed[,,2],
	high_xi_MPK=rmpkmed[,,3]-rmpkmed[,,2], high_xi_MPL=rmplmed[,,3]-rmplmed[,,2], high_xi_MPM=rmpmmed[,,3]-rmpmmed[,,2])
#Productivity Paths
omgpath <- list(Time=1:T, low_xi_low_W=omgmed[,,1][,1]-omgmed[,,2][,1], low_xi_med_W=omgmed[,,1][,2]-omgmed[,,2][,2], low_xi_high_W=omgmed[,,1][,3]-omgmed[,,2][,3],
	high_xi_low_W=omgmed[,,3][,1]-omgmed[,,2][,1], high_xi_med_W=omgmed[,,3][,2]-omgmed[,,2][,2], high_xi_high_W=omgmed[,,3][,3]-omgmed[,,2][,3])
romgpath <- list(Time=1:T, low_xi_low_W=romgmed[,,1][,1]-romgmed[,,2][,1], low_xi_med_W=romgmed[,,1][,2]-romgmed[,,2][,2], low_xi_high_W=romgmed[,,1][,3]-romgmed[,,2][,3],
	high_xi_low_W=romgmed[,,3][,1]-romgmed[,,2][,1], high_xi_med_W=romgmed[,,3][,2]-romgmed[,,2][,2], high_xi_high_W=romgmed[,,3][,3]-romgmed[,,2][,3])
#Plotting
hline <- function(y = 0, color = "red"){list(type = "line", y0 = y, y1 = y, xref = "paper", x0 = 0, x1 = 1, line = list(color = color, dash="dot"))}
Lplotly <- list(); rLplotly <- list()
Mplotly <- list(); rMplotly <- list()
Iplotly <- list(); rIplotly <- list()
Yplotly <- list(); rYplotly <- list()
Misplotly <- list(); rMisplotly <- list()
Wplotly <- list(); rWplotly <- list()
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
	Lplotly[[i]] <- plot_ly(ldat, x=~Time, y=~low_W, type = 'scatter', mode = 'lines', showlegend=F, line=list(color="red"), name=ltitles[[i]][[1]], hovertemplate = paste("<i>Year<i>: %{x}", "<br>Labor Change: %{y:.3f}")) %>%  add_trace(y=~med_W, line=list(color="green"), name=ltitles[[i]][[2]]) %>% add_trace(y=~high_W, line=list(color="blue"), name=ltitles[[i]][[3]]) %>% layout(xaxis=list(title="Years", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="Labor (No R&D)", titlefont=list(size=18), tickfont=list(size=14), range=list(min(unlist(lapply(labpath[names(labpath)!="Time"], min))), max(unlist(lapply(labpath[names(labpath)!="Time"], max))))), shapes=list(hline(y=0)))
	rldat <- data.frame(Time=rlabpath[[1]], low_W=rlabpath[[i+1]][,1], med_W=rlabpath[[i+1]][,2], high_W=rlabpath[[i+1]][,3])
	rLplotly[[i]] <- plot_ly(rldat, x=~Time, y=~low_W, type = 'scatter', mode = 'lines', showlegend=F, line=list(color="red"), name=ltitles[[i]][[1]], hovertemplate = paste("<i>Year<i>: %{x}", "<br>Labor Change: %{y:.3f}")) %>%  add_trace(y=~med_W, line=list(color="green"), name=ltitles[[i]][[2]]) %>% add_trace(y=~high_W, line=list(color="blue"), name=ltitles[[i]][[3]]) %>% layout(xaxis=list(title="Years", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="Labor (R&D)", titlefont=list(size=18), tickfont=list(size=14), range=list(min(unlist(lapply(rlabpath[names(rlabpath)!="Time"], min))), max(unlist(lapply(rlabpath[names(rlabpath)!="Time"], max))))), shapes=list(hline(y=0)))
	#Materials
	mdat <- data.frame(Time=matpath[[1]], low_W=matpath[[i+1]][,1], med_W=matpath[[i+1]][,2], high_W=matpath[[i+1]][,3])
	Mplotly[[i]] <- plot_ly(mdat, x=~Time, y=~low_W, type = 'scatter', mode = 'lines', showlegend=F, line=list(color="red"), name=mtitles[[i]][[1]], hovertemplate = paste("<i>Year<i>: %{x}", "<br>Materials Change: %{y:.3f}")) %>%  add_trace(y=~med_W, line=list(color="green"), name=mtitles[[i]][[2]]) %>% add_trace(y=~high_W, line=list(color="blue"), name=mtitles[[i]][[3]]) %>% layout(xaxis=list(title="Years", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="Materials (No R&D)", titlefont=list(size=18), tickfont=list(size=14), range=list(min(unlist(lapply(matpath[names(matpath)!="Time"], min))), max(unlist(lapply(matpath[names(matpath)!="Time"], max))))), shapes=list(hline(y=0)))
	rmdat <- data.frame(Time=rmatpath[[1]], low_W=rmatpath[[i+1]][,1], med_W=rmatpath[[i+1]][,2], high_W=rmatpath[[i+1]][,3])
	rMplotly[[i]] <- plot_ly(rmdat, x=~Time, y=~low_W, type = 'scatter', mode = 'lines', showlegend=F, line=list(color="red"), name=mtitles[[i]][[1]], hovertemplate = paste("<i>Year<i>: %{x}", "<br>Materials Change: %{y:.3f}")) %>%  add_trace(y=~med_W, line=list(color="green"), name=mtitles[[i]][[2]]) %>% add_trace(y=~high_W, line=list(color="blue"), name=mtitles[[i]][[3]]) %>% layout(xaxis=list(title="Years", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="Materials (R&D)", titlefont=list(size=18), tickfont=list(size=14), range=list(min(unlist(lapply(rmatpath[names(rmatpath)!="Time"], min))), max(unlist(lapply(rmatpath[names(rmatpath)!="Time"], max))))), shapes=list(hline(y=0)))
	#Investment
	idat <- data.frame(Time=ipath[[1]], low_W=ipath[[i+1]][,1], med_W=ipath[[i+1]][,2], high_W=ipath[[i+1]][,3])
	Iplotly[[i]] <- plot_ly(idat, x=~Time, y=~low_W, type = 'scatter', mode = 'lines', showlegend=F, line=list(color="red"), name=ititles[[i]][[1]], hovertemplate = paste("<i>Year<i>: %{x}", "<br>Investment Change: %{y:.3f}")) %>%  add_trace(y=~med_W, line=list(color="green"), name=ititles[[i]][[2]]) %>% add_trace(y=~high_W, line=list(color="blue"), name=ititles[[i]][[3]]) %>% layout(xaxis=list(title="Years", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="Investment (No R&D)", titlefont=list(size=18), tickfont=list(size=14), range=list(min(unlist(lapply(ipath[names(ipath)!="Time"], min))), max(unlist(lapply(ipath[names(ipath)!="Time"], max))))), shapes=list(hline(y=0)))
	ridat <- data.frame(Time=ripath[[1]], low_W=ripath[[i+1]][,1], med_W=ripath[[i+1]][,2], high_W=ripath[[i+1]][,3])
	rIplotly[[i]] <- plot_ly(ridat, x=~Time, y=~low_W, type = 'scatter', mode = 'lines', showlegend=F, line=list(color="red"), name=ititles[[i]][[1]], hovertemplate = paste("<i>Year<i>: %{x}", "<br>Investment Change: %{y:.3f}")) %>%  add_trace(y=~med_W, line=list(color="green"), name=ititles[[i]][[2]]) %>% add_trace(y=~high_W, line=list(color="blue"), name=ititles[[i]][[3]]) %>% layout(xaxis=list(title="Years", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="Investment (R&D)", titlefont=list(size=18), tickfont=list(size=14), range=list(min(unlist(lapply(ripath[names(ripath)!="Time"], min))), max(unlist(lapply(ripath[names(ripath)!="Time"], max))))), shapes=list(hline(y=0)))
	#Output
	ydat <- data.frame(Time=ypath[[1]], low_W=ypath[[i+1]][,1], med_W=ypath[[i+1]][,2], high_W=ypath[[i+1]][,3])
	Yplotly[[i]] <- plot_ly(ydat, x=~Time, y=~low_W, type = 'scatter', mode = 'lines', showlegend=F, line=list(color="red"), name=ytitles[[i]][[1]], hovertemplate = paste("<i>Year<i>: %{x}", "<br>Output Change: %{y:.3f}")) %>%  add_trace(y=~med_W, line=list(color="green"), name=ytitles[[i]][[2]]) %>% add_trace(y=~high_W, line=list(color="blue"), name=ytitles[[i]][[3]]) %>% layout(xaxis=list(title="Years", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="Output (No R&D)", titlefont=list(size=18), tickfont=list(size=14), range=list(min(unlist(lapply(ypath[names(ypath)!="Time"], min))), max(unlist(lapply(ypath[names(ypath)!="Time"], max))))), shapes=list(hline(y=0)))
	rydat <- data.frame(Time=rypath[[1]], low_W=rypath[[i+1]][,1], med_W=rypath[[i+1]][,2], high_W=rypath[[i+1]][,3])
	rYplotly[[i]] <- plot_ly(rydat, x=~Time, y=~low_W, type = 'scatter', mode = 'lines', showlegend=F, line=list(color="red"), name=ytitles[[i]][[1]], hovertemplate = paste("<i>Year<i>: %{x}", "<br>Output Change: %{y:.3f}")) %>%  add_trace(y=~med_W, line=list(color="green"), name=ytitles[[i]][[2]]) %>% add_trace(y=~high_W, line=list(color="blue"), name=ytitles[[i]][[3]]) %>% layout(xaxis=list(title="Years", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="Output (R&D)", titlefont=list(size=18), tickfont=list(size=14), range=list(min(unlist(lapply(rypath[names(rypath)!="Time"], min))), max(unlist(lapply(rypath[names(rypath)!="Time"], max))))), shapes=list(hline(y=0)))
	#Misallocation
	misdat <- data.frame(Time=mispath[[1]], low_W=mispath[[i+1]][,1], med_W=mispath[[i+1]][,2], high_W=mispath[[i+1]][,3])
	Misplotly[[i]] <- plot_ly(misdat, x=~Time, y=~low_W, type = 'scatter', mode = 'lines', showlegend=F, line=list(color="red"), name=mistitles[[i]][[1]], hovertemplate = paste("<i>Year<i>: %{x}", "<br>Misallocation Change: %{y:.3f}")) %>%  add_trace(y=~med_W, line=list(color="green"), name=mistitles[[i]][[2]]) %>% add_trace(y=~high_W, line=list(color="blue"), name=mistitles[[i]][[3]]) %>% layout(xaxis=list(title="Years", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="Misallocation (No R&D)", titlefont=list(size=18), tickfont=list(size=14), range=list(min(unlist(lapply(mispath[names(mispath)!="Time"], min))), max(unlist(lapply(mispath[names(mispath)!="Time"], max))))), shapes=list(hline(y=0)))
	rmisdat <- data.frame(Time=rmispath[[1]], low_W=rmispath[[i+1]][,1], med_W=rmispath[[i+1]][,2], high_W=rmispath[[i+1]][,3])
	rMisplotly[[i]] <- plot_ly(rmisdat, x=~Time, y=~low_W, type = 'scatter', mode = 'lines', showlegend=F, line=list(color="red"), name=mistitles[[i]][[1]], hovertemplate = paste("<i>Year<i>: %{x}", "<br>Misallocation Change: %{y:.3f}")) %>%  add_trace(y=~med_W, line=list(color="green"), name=mistitles[[i]][[2]]) %>% add_trace(y=~high_W, line=list(color="blue"), name=mistitles[[i]][[3]]) %>% layout(xaxis=list(title="Years", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="Misallocation (R&D)", titlefont=list(size=18), tickfont=list(size=14), range=list(min(unlist(lapply(rmispath[names(rmispath)!="Time"], min))), max(unlist(lapply(rmispath[names(rmispath)!="Time"], max))))), shapes=list(hline(y=0)))
	#Productivity
	wdat <- data.frame(Time=omgpath[[1]], Y=omgpath[[i+1]])
	Wplotly[[i]] <- plot_ly(wdat, x=~Time, y=~Y, type = 'scatter', mode = 'lines', showlegend=F, line=list(color=wcolors[[i]]), name=wtitles[[i]], hovertemplate = paste("<i>Year<i>: %{x}", "<br>Productivity Change: %{y:.3f}")) %>% layout(xaxis=list(title="Years", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="Productivity (No R&D)", titlefont=list(size=18), tickfont=list(size=14), range=list(min(unlist(lapply(omgpath[names(omgpath)!="Time"], min))), max(unlist(lapply(omgpath[names(omgpath)!="Time"], max))))), shapes=list(hline(y=0)))
	rwdat <- data.frame(Time=romgpath[[1]], Y=romgpath[[i+1]])
	rWplotly[[i]] <- plot_ly(rwdat, x=~Time, y=~Y, type = 'scatter', mode = 'lines', showlegend=F, line=list(color=wcolors[[i]]), name=wtitles[[i]], hovertemplate = paste("<i>Year<i>: %{x}", "<br>Productivity Change: %{y:.3f}")) %>% layout(xaxis=list(title="Years", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="Productivity (R&D)", titlefont=list(size=18), tickfont=list(size=14), range=list(min(unlist(lapply(romgpath[names(romgpath)!="Time"], min))), max(unlist(lapply(romgpath[names(romgpath)!="Time"], max))))), shapes=list(hline(y=0)))
}
# ##############################################################
# #For Plot.ly
# ###############################################################
#Productivity
annotationsW <- list(list(x=0.07, y=0.9, text=TeX("\\boldsymbol{(a) \\, \\tau_{\\xi}=0.1, \\tau_{\\omega_{1}}=0.1}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
	list(x=0.495, y=0.9, text=TeX("\\boldsymbol{(b)\\,\\tau_{\\xi}=0.1, \\tau_{\\omega_{1}}=0.5}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.88, y=0.9, text=TeX("\\boldsymbol{(c)\\,\\tau_{\\xi}=0.1, \\tau_{\\omega_{1}}=0.9}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.07, y=0.45, text=TeX("\\boldsymbol{(d)\\,\\tau_{\\xi}=0.9, \\tau_{\\omega_{1}}=0.1}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.495, y=0.45, text=TeX("\\boldsymbol{(e)\\,\\tau_{\\xi}=0.9, \\tau_{\\omega_{1}}=0.5}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.88, y=0.45, text=TeX("\\boldsymbol{(f)\\,\\tau_{\\xi}=0.9, \\tau_{\\omega_{1}}=0.9}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE))
W <- subplot(Wplotly[[1]], Wplotly[[2]], Wplotly[[3]], Wplotly[[4]], Wplotly[[5]], Wplotly[[6]], shareX=TRUE, shareY=TRUE, nrows=2)
Wplot <- W %>% layout(annotations=annotationsW) %>% config(mathjax = 'cdn')
Wplot
Wjson <- plotly_json(W, FALSE)
rW <- subplot(rWplotly[[1]], rWplotly[[2]], rWplotly[[3]], rWplotly[[4]], rWplotly[[5]], rWplotly[[6]], shareX=TRUE, shareY=TRUE, nrows=2)
rWplot <- rW %>% layout(annotations=annotationsW) %>% config(mathjax = 'cdn')
rWplot
rWjson <- plotly_json(rW, FALSE)
write(Wjson, "/Users/justindoty/Documents/Home/My_Website/static/jmp/rnd/impulseW.json")
write(rWjson, "/Users/justindoty/Documents/Home/My_Website/static/jmp/rnd/rimpulseW.json")
# #Labor 
annotationsL <- list(list(x=0.11, y=0.92, text=TeX("\\boldsymbol{(a) \\, \\tau_{\\xi}=0.1, \\tau_{l}=0.1}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
	list(x=0.495, y=0.9, text=TeX("\\boldsymbol{(b)\\,\\tau_{\\xi}=0.1, \\tau_{l}=0.5}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.9, y=0.9, text=TeX("\\boldsymbol{(c)\\,\\tau_{\\xi}=0.1, \\tau_{l}=0.9}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.11, y=0.45, text=TeX("\\boldsymbol{(d)\\,\\tau_{\\xi}=0.9, \\tau_{l}=0.1}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.495, y=0.45, text=TeX("\\boldsymbol{(e)\\,\\tau_{\\xi}=0.9, \\tau_{l}=0.5}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.92, y=0.45, text=TeX("\\boldsymbol{(f)\\,\\tau_{\\xi}=0.9, \\tau_{l}=0.9}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE))
L <- subplot(Lplotly[[1]], Lplotly[[2]], Lplotly[[3]], Lplotly[[4]], Lplotly[[5]], Lplotly[[6]], shareX=TRUE, shareY=TRUE, nrows=2)
Lplot <- L %>% layout(annotations=annotationsL) %>% config(mathjax = 'cdn')
Lplot
Ljson <- plotly_json(L, FALSE)
rL <- subplot(rLplotly[[1]], rLplotly[[2]], rLplotly[[3]], rLplotly[[4]], rLplotly[[5]], rLplotly[[6]], shareX=TRUE, shareY=TRUE, nrows=2)
rLplot <- rL %>% layout(annotations=annotationsL) %>% config(mathjax = 'cdn')
rLplot
rLjson <- plotly_json(rL, FALSE)
write(Ljson, "/Users/justindoty/Documents/Home/My_Website/static/jmp/rnd/impulseL.json")
write(rLjson, "/Users/justindoty/Documents/Home/My_Website/static/jmp/rnd/rimpulseL.json")
#Materials
annotationsM <- list(list(x=0.07, y=0.9, text=TeX("\\boldsymbol{(a) \\, \\tau_{\\xi}=0.1, \\tau_{m}=0.1}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
	list(x=0.495, y=0.9, text=TeX("\\boldsymbol{(b)\\,\\tau_{\\xi}=0.1, \\tau_{m}=0.5}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.88, y=0.9, text=TeX("\\boldsymbol{(c)\\,\\tau_{\\xi}=0.1, \\tau_{m}=0.9}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.07, y=0.45, text=TeX("\\boldsymbol{(d)\\,\\tau_{\\xi}=0.9, \\tau_{m}=0.1}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.495, y=0.45, text=TeX("\\boldsymbol{(e)\\,\\tau_{\\xi}=0.9, \\tau_{m}=0.5}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.88, y=0.45, text=TeX("\\boldsymbol{(f)\\,\\tau_{\\xi}=0.9, \\tau_{m}=0.9}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE))
M <- subplot(Mplotly[[1]], Mplotly[[2]], Mplotly[[3]], Mplotly[[4]], Mplotly[[5]], Mplotly[[6]], shareX=TRUE, shareY=TRUE, nrows=2)
Mplot <- M %>% layout(annotations=annotationsM) %>% config(mathjax = 'cdn')
Mplot
Mjson <- plotly_json(M, FALSE)
rM <- subplot(rMplotly[[1]], rMplotly[[2]], rMplotly[[3]], rMplotly[[4]], rMplotly[[5]], rMplotly[[6]], shareX=TRUE, shareY=TRUE, nrows=2)
rMplot <- rM %>% layout(annotations=annotationsM) %>% config(mathjax = 'cdn')
rMplot
rMjson <- plotly_json(rM, FALSE)
write(Mjson, "/Users/justindoty/Documents/Home/My_Website/static/jmp/rnd/impulseM.json")
write(rMjson, "/Users/justindoty/Documents/Home/My_Website/static/jmp/rnd/rimpulseM.json")
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
rI <- subplot(rIplotly[[1]], rIplotly[[2]], rIplotly[[3]], rIplotly[[4]], rIplotly[[5]], rIplotly[[6]], shareX=TRUE, shareY=TRUE, nrows=2)
rIplot <- rI %>% layout(annotations=annotationsI) %>% config(mathjax = 'cdn')
rIplot
rIjson <- plotly_json(rI, FALSE)
write(Ijson, "/Users/justindoty/Documents/Home/My_Website/static/jmp/rnd/impulseI.json")
write(rIjson, "/Users/justindoty/Documents/Home/My_Website/static/jmp/rnd/rimpulseI.json")
annotationsY <- list(list(x=0.07, y=0.9, text=TeX("\\boldsymbol{(a) \\, \\tau_{\\xi}=0.1, \\tau_{y}=0.1}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
	list(x=0.495, y=0.9, text=TeX("\\boldsymbol{(b)\\,\\tau_{\\xi}=0.1, \\tau_{y}=0.5}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.88, y=0.9, text=TeX("\\boldsymbol{(c)\\,\\tau_{\\xi}=0.1, \\tau_{y}=0.9}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.07, y=0.46, text=TeX("\\boldsymbol{(d)\\,\\tau_{\\xi}=0.9, \\tau_{y}=0.1}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.495, y=0.46, text=TeX("\\boldsymbol{(e)\\,\\tau_{\\xi}=0.9, \\tau_{y}=0.5}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.88, y=0.46, text=TeX("\\boldsymbol{(f)\\,\\tau_{\\xi}=0.9, \\tau_{y}=0.9}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE))
# #Output
Y <- subplot(Yplotly[[1]], Yplotly[[2]], Yplotly[[3]], Yplotly[[4]], Yplotly[[5]], Yplotly[[6]], shareX=TRUE, shareY=TRUE, nrows=2)
Yplot <- Y %>% layout(annotations=annotationsY) %>% config(mathjax = 'cdn')
Yplot
Yjson <- plotly_json(Y, FALSE)
rY <- subplot(rYplotly[[1]], rYplotly[[2]], rYplotly[[3]], rYplotly[[4]], rYplotly[[5]], rYplotly[[6]], shareX=TRUE, shareY=TRUE, nrows=2)
rYplot <- rY %>% layout(annotations=annotationsY) %>% config(mathjax = 'cdn')
rYplot
rYjson <- plotly_json(rY, FALSE)
write(Yjson, "/Users/justindoty/Documents/Home/My_Website/static/jmp/rnd/impulseY.json")
write(rYjson, "/Users/justindoty/Documents/Home/My_Website/static/jmp/rnd/rimpulseY.json")
annotationsMis <- list(list(x=0.13, y=0.94, text=TeX("\\boldsymbol{(a) \\, MP_{k}, \\tau_{\\xi}=0.1}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
	list(x=0.495, y=0.94, text=TeX("\\boldsymbol{(b)\\, MP_{l}, \\tau_{\\xi}=0.1}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.88, y=0.94, text=TeX("\\boldsymbol{(c)\\, MP_{m}, \\tau_{\\xi}=0.1}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.13, y=0.46, text=TeX("\\boldsymbol{(d)\\, MP_{k}, \\tau_{\\xi}=0.9}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.495, y=0.46, text=TeX("\\boldsymbol{(e)\\, MP_{l}, \\tau_{\\xi}=0.9}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.88, y=0.46, text=TeX("\\boldsymbol{(f)\\, MP_{m}, \\tau_{\\xi}=0.9}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE))
# #Misallocation
MIS <- subplot(Misplotly[[1]], Misplotly[[2]], Misplotly[[3]], Misplotly[[4]], Misplotly[[5]], Misplotly[[6]], shareX=TRUE, shareY=TRUE, nrows=2)
Misplot <- MIS %>% layout(annotations=annotationsMis) %>% config(mathjax = 'cdn')
Misplot
Misjson <- plotly_json(MIS, FALSE)
rMIS <- subplot(rMisplotly[[1]], rMisplotly[[2]], rMisplotly[[3]], rMisplotly[[4]], rMisplotly[[5]], rMisplotly[[6]], shareX=TRUE, shareY=TRUE, nrows=2)
rMisplot <- rMIS %>% layout(annotations=annotationsMis) %>% config(mathjax = 'cdn')
rMisplot
rMisjson <- plotly_json(rMIS, FALSE)
write(Misjson, "/Users/justindoty/Documents/Home/My_Website/static/jmp/rnd/impulseMis.json")
write(rMisjson, "/Users/justindoty/Documents/Home/My_Website/static/jmp/rnd/rimpulseMis.json")














