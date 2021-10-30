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
omegainit <- as.numeric(omega_est(idvar=US$id, timevar=US$year, Y=US$Y, A=US$A, K=US$K, L=US$L, M=US$M, R=US$R)$omega)
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
US <- US %>% mutate(Y=Y-mean(Y), K=K-mean(K), L=L-mean(L), M=M-mean(M), I=I-mean(I))
lnr <- log(US$R[US$R>0])
lnr <- lnr-mean(lnr)
US$R[US$R>0] <- lnr
wmin <- min(US$Y)
wmax <- max(US$Y)
#Load the method used (e.g. Cobb, Translog, Hermite)
method <- results$method
resY <- results$resY
###############################################################################
#Simulate Productivity and Capital Evolution Given Model Parameters
#############################################################################
#Vector of ranks of input demand functions (small, medium, large)
tauinp <- c(0.1, 0.5, 0.9)
#Shocks to productivity
tauxi <- c(0.1, 0.5, 0.9)
#Rank of initial productivity
tauinit <- c(0.1, 0.5, 0.9)
#########################################################################################################
NR <- 15
NNR <- 10
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
#Position in vector for RnD Firms
pos <- c((nnr+1):N)
#Labor###############################################################
lnldata <- array(0, c(N, T, length(tauxi), length(tauinp)))
#At the median
#For Non R&D Firms
lnlmed <- array(0, c(T, length(tauxi), length(tauinp)))
#For R&D Firms
rlnlmed <- array(0, c(T, length(tauxi), length(tauinp)))
#Shocks to Labor
epsl <- array(0, c(N, T, length(tauinp)))
#Intermediate Input#####################################################
lnmdata <- array(0, c(N, T, length(tauxi), length(tauinp)))
#At the median
#For Non R&D Firms
lnmmed <- array(0, c(T, length(tauxi), length(tauinp)))
#For R&D Firms
rlnmmed <- array(0, c(T, length(tauxi), length(tauinp)))
#Shocks to Materials
epsm <- array(0, c(N, T, length(tauinp)))
#Investment###############################################################
lnidata <- array(0, c(N, T, length(tauxi)))
iota <- matrix(0.5, nrow=N, ncol=T)
#This is for when i vary the shocks to investment demand
qlnidata <- array(0, c(N, T, length(tauxi), length(tauinp)))
#For Non R&D Firms
qlnimed <- array(0, c(T, length(tauxi), length(tauinp)))
#For R&D Firms
rqlnimed <- array(0, c(T, length(tauxi), length(tauinp)))
#Shocks to Investment
qiota <- array(0, c(N, T, length(tauinp)))
#Capital##################################################################
lnkdata <- array(0, c(N, T, length(tauxi)))
#This is for when i vary the shocks to investment demand
qlnkdata <- array(0, c(N, T, length(tauxi), length(tauinp)))
#For Non R&D Firms
qlnkmed <- array(0, c(T, length(tauxi), length(tauinp)))
#For R&D Firms
rqlnkmed <- array(0, c(T, length(tauxi), length(tauinp)))
#Productivity#############################################################
#This is used for labor, materials, and investment
omgdata <- array(0, c(N, T, length(tauxi)))
#This is for when i vary productivity paths for different initial productivity levels
omgqdata <- array(0, c(N, T, length(tauxi), length(tauinit)))
#For when I vary productivity paths for different R&D sizes
omgrsize <- array(0, c(N, T, length(tauxi), length(tauinit)))
omgrsizemed <- array(0, c(T, length(tauxi), length(tauinit)))
#For Non R&D Firms
omgqmed <- array(0, c(T, length(tauxi), length(tauinp)))
#For R&D Firms
romgqmed <- array(0, c(T, length(tauxi), length(tauinp)))
#Innovation Shocks
xi <- matrix(runif(N*T), nrow=N, ncol=T)
xidata <- replicate(length(tauxi), xi, simplify="array")
#R&D
rdata <- array(0, c(N, T, length(tauxi)))
#For R&D firms of different initial productivity levels 
rqdata <- array(0, c(N, T, length(tauxi), length(tauinp)))
#For R&D firms of different investment levels
rqidata <- array(0, c(N, T, length(tauxi), length(tauinp)))
#Different sized R&D Firms
rsizedata <- array(0, c(N, T, length(tauxi), length(tauinp)))
#R&D Response to Productivity
rsizemed <- array(0, c(T, length(tauxi), length(tauinp)))
#Shocks to R&D
rhodata <- matrix(runif(N*T), nrow=N, ncol=T)
rhoq <- array(0, c(N, T, length(tauinp)))
#Ranks for Input Shocks
for (q2 in 1:length(tauinp)){
	epsl[,,q2] <- array(tauinp[q2], c(N,T))
	epsm[,,q2] <- array(tauinp[q2], c(N,T))
	qiota[,,q2] <- array(tauinp[q2], c(N,T))
	rhoq[,,q2] <- array(tauinp[q2], c(N,T))
}
#I order Non RnD Firms First
k1 <- c(lnk1, rlnk1)
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
#Initial R&D
r1 <- rowSums(RX(K=k1[pos], omega=omg1m)*lspline(vectau=vectau, bvec=parR, b1=parRb[1], bL=parRb[2], u=rhodata[pos,1]))
#Restrict Support of Initial R&D
r1 <- (r1>max(t1data$R))*max(t1data$R)+(r1<min(t1data$R[t1data$RB==1]))*min(t1data$R[t1data$RB==1])+(r1<=max(t1data$R))*(r1>=min(t1data$R[t1data$RB==1]))*r1
#Initial Inputs
for (q1 in 1:length(tauxi)){
	omgdata[,,q1][,1] <- omg1m
	lnkdata[,,q1][,1] <- k1
	lnidata[,,q1][,1] <- i1
	rdata[,,q1][pos,] <- r1
	#At t=2 give varying shock to productivity given by tauxi
	xidata[,,q1][,2] <- tauxi[q1]
	for (q2 in 1:length(tauinp)){
		k1 <- rep(mean(k1), length(k1))
		l1 <- rowSums(LX(K=k1, omega=omg1m)*lspline(vectau=vectau, bvec=parL, b1=parLb[1], bL=parLb[2], u=epsl[,,q2][,1]))
		m1 <- rowSums(MX(K=k1, L=l1, omega=omg1m)*lspline(vectau=vectau, bvec=parM, b1=parMb[1], bL=parMb[2], u=epsm[,,q2][,1]))
		rq <- rowSums(RX(K=k1[pos], omega=omg1q[q2])*lspline(vectau=vectau, bvec=parR, b1=parRb[1], bL=parRb[2], u=rhodata[pos,1]))
		rqs <- rowSums(RX(K=k1[pos], omega=omg1m)*lspline(vectau=vectau, bvec=parR, b1=parRb[1], bL=parRb[2], u=rhoq[,,q2][pos,1]))
		#Restricting the Supports
		l1 <- (l1>max(t1data$L))*max(t1data$L)+(l1<min(t1data$L))*min(t1data$L)+(l1<=max(t1data$L))*(l1>=min(t1data$L))*l1
		m1 <- (m1>max(t1data$M))*max(t1data$M)+(m1<min(t1data$M))*min(t1data$M)+(m1<=max(t1data$M))*(m1>=min(t1data$M))*m1
		rq <- (rq>max(t1data$R))*max(t1data$R)+(rq<min(t1data$R[t1data$RB==1]))*min(t1data$R[t1data$RB==1])+(rq<=max(t1data$R))*(rq>=min(t1data$R[t1data$RB==1]))*rq
		rqs <- (rqs>max(t1data$R))*max(t1data$R)+(rqs<min(t1data$R[t1data$RB==1]))*min(t1data$R[t1data$RB==1])+(rqs<=max(t1data$R))*(rqs>=min(t1data$R[t1data$RB==1]))*rqs
		#Median Across Firms
		lnldata[,,,q2][,,q1][,1] <- l1
		lnmdata[,,,q2][,,q1][,1] <- m1
		rqdata[,,,q2][,,q1][pos,1] <- rq
		rsizedata[,,,q2][,,q1][pos,1] <- rqs
		#Different levels of initial productivity
		omgqdata[,,,q2][,,q1][,1] <- omg1q[q2]
		omgrsize[,,,q2][,,q1][,1] <- omg1m

	}
}
#Generate Productivity
for (t in 2:T){
	ttdata <- US %>% group_by(id) %>% slice(t)
	for (q1 in 1:length(tauxi)){
		rind <- rdata[,,q1][,t-1]!=0
		omg <- rowSums(WX(omega=omgdata[,,q1][,t-1], R=rdata[,,q1][,t-1], Rind=rind)*lspline(vectau=vectau, bvec=parWT, b1=parWTb[1], bL=parWTb[2], u=xidata[,,q1][,t]))
		#Restricting the Supports
		omgdata[,,q1][,t] <- (omg>wmax)*wmax+(omg<wmin)*wmin+(omg<=wmax)*(omg>=wmin)*omg
		#Generate Capital According to Accumulation Process with Industry-Average Depreciation Rates
		lnk <- log(0.98*exp(lnkdata[,,q1][,t-1])+exp(lnidata[,,q1][,t-1]))
		#Restricting Supports
		lnkdata[,,q1][,t] <- (lnk>max(ttdata$K))*max(ttdata$K)+(lnk<min(ttdata$K))*min(ttdata$K)+(lnk<=max(ttdata$K))*(lnk>=min(ttdata$K))*lnk
		#Generate Investment
		lni <- rowSums(IX(K=lnkdata[,,q1][,t], omega=omgdata[,,q1][,t])*lspline(vectau=vectau, bvec=parI, b1=parIb[1], bL=parIb[2], u=iota[,t]))
		#Restricting Supports
		lnidata[,,q1][,t] <- (lni>max(ttdata$I))*max(ttdata$I)+(lni<min(ttdata$I))*min(ttdata$I)+(lni<=max(ttdata$I))*(lni>=min(ttdata$I))*lni
		#R&D
		lnr <- rowSums(RX(K=mean(lnkdata[,,q1][pos,t]), omega=omgdata[,,q1][pos,t])*lspline(vectau=vectau, bvec=parR, b1=parRb[1], bL=parRb[2], u=rhodata[pos,t]))
		#Restricting Supports
		rdata[,,q1][pos,t] <- (lnr>max(ttdata$R))*max(ttdata$R)+(lnr<min(ttdata$R[ttdata$RB==1]))*min(ttdata$R[ttdata$RB==1])+(lnr<=max(ttdata$R))*(lnr>=min(ttdata$R[ttdata$RB==1]))*lnr
		#Generate Optimal Input Decisions Following Innovation Shocks for Various Ranks of Demand Functions
		for (q2 in 1:length(tauinp)){
			rindq <- rqdata[,,,q2][,,q1][,t-1]!=0
			#Productivity
			omgq <- rowSums(WX(omega=omgqdata[,,,q2][,,q1][,t-1], R=rqdata[,,,q2][,,q1][,t-1], Rind=rindq)*lspline(vectau=vectau, bvec=parWT, b1=parWTb[1], bL=parWTb[2], u=xidata[,,q1][,t]))
			omgrqs <- rowSums(WX(omega=omgrsize[,,,q2][,,q1][,t-1], R=mean(rsizedata[,,,q2][,,q1][,t-1]), Rind=rindq)*lspline(vectau=vectau, bvec=parWT, b1=parWTb[1], bL=parWTb[2], u=xidata[,,q1][,t]))
			#R&D
			lnrq <- rowSums(RX(K=mean(lnkdata[,,q1][pos,t]), omega=omgqdata[,,,q2][,,q1][pos,t])*lspline(vectau=vectau, bvec=parR, b1=parRb[1], bL=parRb[2], u=rhodata[pos,t]))
			lnrqs <- rowSums(RX(K=mean(lnkdata[,,q1][pos,t]), omega=omgrsize[,,,q2][,,q1][pos,t])*lspline(vectau=vectau, bvec=parR, b1=parRb[1], bL=parRb[2], u=rhoq[,,q2][pos,t]))
			#Restricting the Supports
			omgqdata[,,,q2][,,q1][,t] <- (omgq>wmax)*wmax+(omgq<wmin)*wmin+(omgq<=wmax)*(omgq>=wmin)*omgq
			omgrsize[,,,q2][,,q1][,t] <- (omgrqs>wmax)*wmax+(omgrqs<wmin)*wmin+(omgrqs<=wmax)*(omgrqs>=wmin)*omgrqs
			rqdata[,,,q2][,,q1][pos,t] <- (lnrq>max(ttdata$R[ttdata$RB==1]))*max(ttdata$R[ttdata$RB==1])+(lnrq<min(ttdata$R[ttdata$RB==1]))*min(ttdata$R[ttdata$RB==1])+(lnrq<=max(ttdata$R[ttdata$RB==1]))*(lnrq>=min(ttdata$R[ttdata$RB==1]))*lnrq
			rsizedata[,,,q2][,,q1][pos,t] <- (lnrqs>max(ttdata$R[ttdata$RB==1]))*max(ttdata$R[ttdata$RB==1])+(lnrqs<min(ttdata$R[ttdata$RB==1]))*min(ttdata$R[ttdata$RB==1])+(lnrqs<=max(ttdata$R[ttdata$RB==1]))*(lnrqs>=min(ttdata$R[ttdata$RB==1]))*lnrqs
		}
	}
}
#Generate Labor and R&D at average levels of capital
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
#Generate Materials at average levels of capital and labor
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
#Take the median across firms
for (q1 in 1:length(tauxi)){
	for (q2 in 1:length(tauinp)){
		#For Non R&D Firms
		omgqmed[,,q2][,q1] <- apply(omgqdata[,,,q2][,,q1][-pos,], 2, median)
		lnlmed[,,q2][,q1] <- apply(lnldata[,,,q2][,,q1][-pos,], 2, median)
		lnmmed[,,q2][,q1] <- apply(lnmdata[,,,q2][,,q1][-pos,], 2, median)
		#For R&D Firms
		romgqmed[,,q2][,q1] <- apply(omgqdata[,,,q2][,,q1][pos,], 2, median)
		omgrsizemed[,,q2][,q1] <- apply(omgrsize[,,,q2][,,q1][pos,], 2, median)
		rlnlmed[,,q2][,q1] <- apply(lnldata[,,,q2][,,q1][pos,], 2, median)
		rlnmmed[,,q2][,q1] <- apply(lnmdata[,,,q2][,,q1][pos,], 2, median)
		rsizemed[,,q2][,q1] <- apply(rsizedata[,,,q2][,,q1][pos,], 2, median)
	}
}
#For the arrays, lnldata and lnmdata, the 1st dimension is time, 2nd is rank of innovation shock, 3rd is rank of input shock
#So "Low-Low" represents low labor shock and low innovation shock
#Labor Path for Non R&D Firms
labpath <- data.frame(1:T, lnlmed[,1,]-lnlmed[,2,], lnlmed[,3,]-lnlmed[,2,])
names(labpath) <- c("Time", "LowLow", "LowMed", "LowHigh", "HighLow", "HighMed", "HighHigh")
#Labor Path for R&D Firms
rlabpath <- data.frame(1:T, rlnlmed[,1,]-rlnlmed[,2,], rlnlmed[,3,]-rlnlmed[,2,])
names(rlabpath) <- c("Time", "LowLow", "LowMed", "LowHigh", "HighLow", "HighMed", "HighHigh")
#Materials Path for Non R&D Firms
matpath <- data.frame(1:T, lnmmed[,1,]-lnmmed[,2,], lnmmed[,3,]-lnmmed[,2,])
names(matpath) <- c("Time", "LowLow", "LowMed", "LowHigh", "HighLow", "HighMed", "HighHigh")
#Materials Path for R&D Firms
rmatpath <- data.frame(1:T, rlnmmed[,1,]-rlnmmed[,2,], rlnmmed[,3,]-rlnmmed[,2,])
names(rmatpath) <- c("Time", "LowLow", "LowMed", "LowHigh", "HighLow", "HighMed", "HighHigh")
#Productivity Path for Non R&D Firms
omegapath <- data.frame(1:T, omgqmed[,1,]-omgqmed[,2,], omgqmed[,3,]-omgqmed[,2,])
names(omegapath) <- c("Time", "LowLow", "LowMed", "LowHigh", "HighLow", "HighMed", "HighHigh")
#Productivity Path for R&D Firms
omgrpath <- data.frame(1:T, omgrsizemed[,1,]-omgrsizemed[,2,], omgrsizemed[,3,]-omgrsizemed[,2,])
names(omgrpath) <- c("Time", "LowLow", "LowMed", "LowHigh", "HighLow", "HighMed", "HighHigh")
#R&D Paths
rsizepath <- data.frame(1:T, rsizemed[,1,]-rsizemed[,2,], rsizemed[,3,]-rsizemed[,2,])
names(rsizepath) <- c("Time", "LowLow", "LowMed", "LowHigh", "HighLow", "HighMed", "HighHigh")
#Plotting
hline <- function(y = 0, color = "red") {
  list(
    type = "line",
    y0 = y,
    y1 = y,
    xref = "paper",
    x0 = 0,
    x1 = 1,
    line = list(color = color, dash="dot")
  )
}
ltitlesNOR <- list("No R&D<br>Negative Shock<br>Low Labor", "No R&D<br>Negative Shock<br>Medium Labor", "No R&D<br>Negative Shock<br>High Labor", "No R&D<br>Positive Shock<br>Low Labor", "No R&D<br>Positive Shock<br>Medium Labor", "No R&D<br>Positive Shock<br>High Labor")
mtitlesNOR <- list("No R&D<br>Negative Shock<br>Low Materials", "No R&D<br>Negative Shock<br>Medium Materials", "No R&D<br>Negative Shock<br>High Materials", "No R&D<br>Positive Shock<br>Low Materials", "No R&D<br>Positive Shock<br>Medium Materials", "No R&D<br>Positive Shock<br>High Materials")
wtitlesNOR <- list("No R&D<br>Negative Shock<br>Low Productivity", "No R&D<br>Negative Shock<br>Medium Productivity", "No R&D<br>Negative Shock<br>High Productivity", "No R&D<br>Positive Shock<br>Low Productivity", "No R&D<br>Positive Shock<br>Medium Productivity", "No R&D<br>Positive Shock<br>High Productivity")
ltitlesR <- list("R&D<br>Negative Shock<br>Low Labor", "R&D<br>Negative Shock<br>Medium Labor", "R&D<br>Negative Shock<br>High Labor", "R&D<br>Positive Shock<br>Low Labor", "R&D<br>Positive Shock<br>Medium Labor", "R&D<br>Positive Shock<br>High Labor")
mtitlesR <- list("R&D<br>Negative Shock<br>Low Materials", "R&D<br>Negative Shock<br>Medium Materials", "R&D<br>Negative Shock<br>High Materials", "R&D<br>Positive Shock<br>Low Materials", "R&D<br>Positive Shock<br>Medium Materials", "R&D<br>Positive Shock<br>High Materials")
wtitlesR <- list("R&D<br>Negative Shock<br>Low Productivity", "R&D<br>Negative Shock<br>Medium Productivity", "R&D<br>Negative Shock<br>High Productivity", "R&D<br>Positive Shock<br>Low Productivity", "R&D<br>Positive Shock<br>Medium Productivity", "R&D<br>Positive Shock<br>High Productivity")
rtitles <- list("Negative Shock<br>Low R&D", "Negative Shock<br>Medium R&D", "Negative Shock<br>High R&D", "Positive Shock<br>Low R&D", "Positive Shock<br>Medium R&D", "Positive Shock<br>High R&D")
Wplotly <- list()
Lplotly <- list()
Mplotly <- list()
Rplotly <- list()
RWplotly <- list()
for (i in 1:6){
	wdat <- data.frame(Time=omegapath$Time, Y=omegapath[,i+1], Z=omgrpath[,i+1])
	Wplotly[[i]] <- plot_ly(wdat, x=~Time, y=~Y, type = 'scatter', mode = 'lines', showlegend=F, line=list(color="black"), name=wtitlesNOR[[i]], hovertemplate = paste("<i>Year<i>: %{x}", "<br>Productivity Change: %{y:.2f}")) %>% add_trace(y = ~Z, mode = 'lines', showlegend=F, line=list(color="green", dash="dash"), name=wtitlesR[[i]], hovertemplate = paste("<i>Year<i>: %{x}", "<br>Productivity Change: %{y:.2f}")) %>% layout(xaxis=list(title="Years", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="Productivity", titlefont=list(size=18), tickfont=list(size=14), range=list(min(min(omegapath[,-1]), min(omgrpath[,-1])), max(max(omegapath[,-1]), max(omgrpath[,-1])))), shapes=list(hline(y=0)))
	ldat <- data.frame(Time=labpath$Time, Y=labpath[,i+1], Z=rlabpath[,i+1])
	Lplotly[[i]] <- plot_ly(ldat, x=~Time, y=~Y, type = 'scatter', mode = 'lines', showlegend=F, line=list(color="black"), name=ltitlesNOR[[i]], hovertemplate = paste("<i>Year<i>: %{x}", "<br>Labor Change: %{y:.2f}")) %>% add_trace(y = ~Z, mode = 'lines', showlegend=F,line=list(color="green", dash="dash"), name=ltitlesR[[i]], hovertemplate = paste("<i>Year<i>: %{x}", "<br>Labor Change: %{y:.2f}")) %>% layout(xaxis=list(title="Years", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="Labor", titlefont=list(size=18), tickfont=list(size=14), range=list(min(min(labpath[,-1]), min(rlabpath[,-1])), max(max(labpath[,-1]), max(rlabpath[,-1])))), shapes=list(hline(y=0)))
	mdat <- data.frame(Time=matpath$Time, Y=matpath[,i+1], Z=rmatpath[,i+1])
	Mplotly[[i]] <- plot_ly(mdat, x=~Time, y=~Y, type = 'scatter', mode = 'lines', showlegend=F, line=list(color="black"), name=mtitlesNOR[[i]], hovertemplate = paste("<i>Year<i>: %{x}", "<br>Materials Change: %{y:.2f}")) %>% add_trace(y = ~Z, mode = 'lines', showlegend=F,line=list(color="green", dash="dash"), name=mtitlesR[[i]], hovertemplate = paste("<i>Year<i>: %{x}", "<br>Materials Change: %{y:.2f}")) %>% layout(xaxis=list(title="Years", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="Materials", titlefont=list(size=18), tickfont=list(size=14), range=list(min(min(matpath[,-1]), min(rmatpath[,-1])), max(max(matpath[,-1]), max(rmatpath[,-1])))), shapes=list(hline(y=0)))
	rdat <- data.frame(Time=omegapath$Time, Y=omgrpath[,i+1], Z=rsizepath[,i+1])
	RWplotly[[i]] <- plot_ly(rdat, x=~Time, y=~Y, type = 'scatter', mode = 'lines', showlegend=F, line=list(color="black"), name=rtitles[[i]], hovertemplate = paste("<i>Year<i>: %{x}", "<br>Productivity Change: %{y:.2f}")) %>% layout(xaxis=list(title="Years"), yaxis=list(title="Productivity", range=list(min(omgrpath[,-1]), max(omgrpath[,-1]))), shapes=list(hline(y=0)))

	
}
# Capital##################################################################################################
# I do the same procedure for capital except at different levels of shock to investment
#Initial Capital and Investment
for (q1 in 1:length(tauxi)){
	for (q2 in 1:length(tauinp)){
		qlnkdata[,,,q2][,,q1][,1] <- k1
		qlnkmed[,,q2][1,q1] <- median(qlnkdata[,,,q2][,,q1][,1])
		qi1 <- rowSums(IX(K=k1, omega=omg1)*lspline(vectau=vectau, bvec=parI, b1=parIb[1], bL=parIb[2], u=qiota[,,q2][,1]))
		rq <- rowSums(RX(K=k1[pos], omega=omg1[pos])*lspline(vectau=vectau, bvec=parR, b1=parRb[1], bL=parRb[2], u=rhodata[pos,1])) 
		#Restrict the Support
		qi1 <- (qi1>max(t1data$I))*max(t1data$I)+(qi1<min(t1data$I))*min(t1data$I)+(qi1<=max(t1data$I))*(qi1>=min(t1data$I))*qi1
		rq <- (rq>max(t1data$R))*max(t1data$R)+(rq<min(t1data$R[t1data$RB==1]))*min(t1data$R[t1data$RB==1])+(rq<=max(t1data$R))*(rq>=min(t1data$R[t1data$RB==1]))*rq
		#Initial Investment
		qlnidata[,,,q2][,,q1][,1] <- qi1
		qlnimed[,,q2][1,q1] <- median(qi1)
		#Initial R&D
		rqidata[,,,q2][,,q1][pos,1] <- rq

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
			rq <- rowSums(RX(K=mean(qlnkdata[,,,q2][,,q1][pos,t]), omega=omgdata[,,q1][pos,t])*lspline(vectau=vectau, bvec=parR, b1=parRb[1], bL=parRb[2], u=rhodata[pos,t])) 
			#Restrict the Supports
			qlnidata[,,,q2][,,q1][,t] <- (qit>max(ttdata$I))*max(ttdata$I)+(qit<min(ttdata$I))*min(ttdata$I)+(qit<=max(ttdata$I))*(qit>=min(ttdata$I))*qit
			rqidata[,,,q2][,,q1][pos,t] <- (rq>max(ttdata$R))*max(ttdata$R)+(rq<min(ttdata$R[ttdata$RB==1]))*min(ttdata$R[ttdata$RB==1])+(rq<=max(ttdata$R))*(rq>=min(ttdata$R[ttdata$RB==1]))*rq
			#Medians across firms
			#For Non R&D Firms
			qlnimed[,,q2][t,q1] <- median(qlnidata[,,,q2][,,q1][-pos,t])
			qlnkmed[,,q2][t,q1] <- median(qlnkdata[,,,q2][,,q1][-pos,t])
			#For R&D Firms
			rqlnimed[,,q2][t,q1] <- median(qlnidata[,,,q2][,,q1][pos,t])
			rqlnkmed[,,q2][t,q1] <- median(qlnkdata[,,,q2][,,q1][pos,t])
		}
	}
}
##Capital Path for Non R&D Firms
kpath <- data.frame(1:T, qlnkmed[,1,]-qlnkmed[,2,], qlnkmed[,3,]-qlnkmed[,2,])
names(kpath) <- c("Time", "LowLow", "LowMed", "LowHigh", "HighLow", "HighMed", "HighHigh")
#Capital Path for R&D Firms
rkpath <- data.frame(1:T, rqlnkmed[,1,]-rqlnkmed[,2,], rqlnkmed[,3,]-rqlnkmed[,2,])
names(rkpath) <- c("Time", "LowLow", "LowMed", "LowHigh", "HighLow", "HighMed", "HighHigh")
ktitlesNOR <- list("No R&D<br>Negative Shock<br>Low Investment", "No R&D<br>Negative Shock<br>Medium Investment", "No R&D<br>Negative Shock<br>High Investment", "No R&D<br>Positive Shock<br>Low Investment", "No R&D<br>Positive Shock<br>Medium Investment", "No R&D<br>Positive Shock<br>High Investment")
ktitlesR <- list("R&D<br>Negative Shock<br>Low Investment", "R&D<br>Negative Shock<br>Medium Investment", "R&D<br>Negative Shock<br>High Investment", "R&D<br>Positive Shock<br>Low Investment", "R&D<br>Positive Shock<br>Medium Investment", "R&D<br>Positive Shock<br>High Investment")
#Plotting
Kplotly <- list()
for (i in 1:6){
	kdat <- data.frame(Time=kpath$Time, Y=kpath[,i+1], Z=rkpath[,i+1])
	Kplotly[[i]] <- plot_ly(kdat, x=~Time, y=~Y, type = 'scatter', mode = 'lines', showlegend=F, line=list(color="black"), name=ktitlesNOR[[i]], hovertemplate = paste("<i>Year<i>: %{x}", "<br>Capital Change: %{y:.2f}")) %>% add_trace(y = ~Z, mode = 'lines', showlegend=F,line=list(color="green", dash="dash"), name=ktitlesR[[i]], hovertemplate = paste("<i>Year<i>: %{x}", "<br>Capital Change: %{y:.2f}")) %>% layout(xaxis=list(title="Years", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="Capital", titlefont=list(size=18), tickfont=list(size=14), range=list(min(min(kpath[,-1]), min(rkpath[,-1])), max(max(kpath[,-1]), max(rkpath[,-1])))), shapes=list(hline(y=0)))
}
# ##############################################################
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
# write(Wjson, "/Users/justindoty/Documents/Home/My_Website/static/jmp/rnd/impulseW.json")
#Labor 
annotationsL <- list(list(x=0.11, y=0.9, text=TeX("\\boldsymbol{(a) \\, \\tau_{\\xi}=0.1, \\tau_{l}=0.1}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
	list(x=0.495, y=0.9, text=TeX("\\boldsymbol{(b)\\,\\tau_{\\xi}=0.1, \\tau_{l}=0.5}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.88, y=0.9, text=TeX("\\boldsymbol{(c)\\,\\tau_{\\xi}=0.1, \\tau_{l}=0.9}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.11, y=0.45, text=TeX("\\boldsymbol{(d)\\,\\tau_{\\xi}=0.9, \\tau_{l}=0.1}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.495, y=0.45, text=TeX("\\boldsymbol{(e)\\,\\tau_{\\xi}=0.9, \\tau_{l}=0.5}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.88, y=0.45, text=TeX("\\boldsymbol{(f)\\,\\tau_{\\xi}=0.9, \\tau_{l}=0.9}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE))
L <- subplot(Lplotly[[1]], Lplotly[[2]], Lplotly[[3]], Lplotly[[4]], Lplotly[[5]], Lplotly[[6]], shareX=TRUE, shareY=TRUE, nrows=2)
Lplot <- L %>% layout(annotations=annotationsL) %>% config(mathjax = 'cdn')
Lplot
Ljson <- plotly_json(L, FALSE)
# write(Ljson, "/Users/justindoty/Documents/Home/My_Website/static/jmp/rnd/impulseL.json")
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
# write(Mjson, "/Users/justindoty/Documents/Home/My_Website/static/jmp/rnd/impulseM.json")
#Capital
annotationsI <- list(list(x=0.11, y=0.9, text=TeX("\\boldsymbol{(a) \\, \\tau_{\\xi}=0.1, \\tau_{i}=0.1}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
	list(x=0.495, y=0.9, text=TeX("\\boldsymbol{(b)\\,\\tau_{\\xi}=0.1, \\tau_{i}=0.5}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.88, y=0.9, text=TeX("\\boldsymbol{(c)\\,\\tau_{\\xi}=0.1, \\tau_{i}=0.9}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.11, y=0.46, text=TeX("\\boldsymbol{(d)\\,\\tau_{\\xi}=0.9, \\tau_{i}=0.1}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.495, y=0.46, text=TeX("\\boldsymbol{(e)\\,\\tau_{\\xi}=0.9, \\tau_{i}=0.5}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.88, y=0.46, text=TeX("\\boldsymbol{(f)\\,\\tau_{\\xi}=0.9, \\tau_{i}=0.9}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE))
K <- subplot(Kplotly[[1]], Kplotly[[2]], Kplotly[[3]], Kplotly[[4]], Kplotly[[5]], Kplotly[[6]], shareX=TRUE, shareY=TRUE, nrows=2)
Kplot <- K %>% layout(annotations=annotationsI) %>% config(mathjax = 'cdn')
Kplot
Kjson <- plotly_json(K, FALSE)
# write(Kjson, "/Users/justindoty/Documents/Home/My_Website/static/jmp/rnd/impulseK.json")
#R&D Productivity
annotationsR <- list(list(x=0.11, y=0.9, text=TeX("\\boldsymbol{(a) \\, \\tau_{\\xi}=0.1, \\tau_{r}=0.1}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
	list(x=0.495, y=0.9, text=TeX("\\boldsymbol{(b)\\,\\tau_{\\xi}=0.1, \\tau_{r}=0.5}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.88, y=0.9, text=TeX("\\boldsymbol{(c)\\,\\tau_{\\xi}=0.1, \\tau_{r}=0.9}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.11, y=0.46, text=TeX("\\boldsymbol{(d)\\,\\tau_{\\xi}=0.9, \\tau_{r}=0.1}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.495, y=0.46, text=TeX("\\boldsymbol{(e)\\,\\tau_{\\xi}=0.9, \\tau_{r}=0.5}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.88, y=0.46, text=TeX("\\boldsymbol{(f)\\,\\tau_{\\xi}=0.9, \\tau_{r}=0.9}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE))
RW <- subplot(RWplotly[[1]], RWplotly[[2]], RWplotly[[3]], RWplotly[[4]], RWplotly[[5]], RWplotly[[6]], shareX=TRUE, shareY=TRUE, nrows=2)
RWplot <- RW %>% layout(annotations=annotationsR) %>% config(mathjax = 'cdn')
RWplot
RWploty <- plotly_json(RW, FALSE)
# write(RWploty, "/Users/justindoty/Documents/Home/My_Website/static/jmp/rnd/impulseRW.json")



















