require(cowplot)
require(dplyr)
require(purrr)
require(xtable)
require(stringr)
require(truncnorm)
require(RColorBrewer)
require(reshape2)
require(gridGraphics)
require(plotly)
source('/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Functions/Selection_Bias/Auxfuns.R')
source('/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Functions/Selection_Bias/Tensors.R')
source('/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Functions/Selection_Bias/omega.R')
#Load US Dataset
US <- read.csv('/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Data/USdata.csv') %>% 
select(id, year, lny, lnk1, lnk2, lnl, lnm, lni, age) %>% transmute(id=id, time=year, Y=lny, K=lnk1, L=lnl, M=lnm, I=lni, A=age)
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
#De-mean
stdY <- sd(US$Y)
stdK <- sd(US$K)
stdL <- sd(US$L)
stdM <- sd(US$M)
stdI <- sd(US$I)
#De-mean
US <- US %>% mutate(Y=(Y-mean(Y))/stdY, K=(K-mean(K))/stdK, L=(L-mean(L))/stdL, M=(M-mean(M))/stdM, I=(I-mean(I))/stdI)
wmin <- min(US$Y)
wmax <- max(US$Y)
#Load the method used (e.g. Cobb, Translog, Hermite)
method <- results$method
resY <- results$resY
###############################################################################
#Simulate Productivity and Capital Evolution Given Model Parameters
##############################################################################
Nsim <- 5
N <- length(unique(US$id))*Nsim
T <- length(unique(US$time))
#Labor
lnldata <- matrix(0, N, T) 
#Unobservable Shock to Labor
epsdata <- matrix(runif(N*T), nrow=N, ncol=T)
#Intermediate Input
lnmdata <- matrix(0, N, T)
#Unobservable Shock to Intermediate Inputs
varepsdata <- matrix(runif(N*T), nrow=N, ncol=T)
#Investment
lnidata <- matrix(0, N, T)
#Unobservable Shock to Intermediate Inputs
iotadata <- matrix(runif(N*T), nrow=N, ncol=T)
#Capital
lnkdata <- matrix(0, N, T)
#Output
lnydata <- matrix(0, N, T)
#Unobservable Shock to Output
etadata <- matrix(runif(N*T), nrow=N, ncol=T)
#Productivity
omgdata <- matrix(0, N, T)
#Unobservable Shock to Productivity
xidata <- matrix(runif(N*T), nrow=N, ncol=T)
#Exit Probabilities
pdata <- matrix(0, N, T)
#Simulate Data at T=1
t1data <- US %>% group_by(id) %>% slice(1)
#########################################################################################################
#Capital is not estimated in the Selection_Bias model: Use Capital Accumulation Process
#For K=1
lnkdata[,1] <- kronecker(array(1, c(Nsim,1)), t1data$K)
#Initial Productivity
omgdata[,1] <- rowSums(WX1(K=lnkdata[,1], L=lnldata[,1], M=lnmdata[,1])*lspline(vectau=vectau, bvec=parW1, b1=parW1b[1], bL=parW1b[2], u=xidata[,1]))
#Restrict the Support of Initial Productivity
omgdata[,1] <- (omgdata[,1]>wmax)*wmax+(omgdata[,1]<wmin)*wmin+(omgdata[,1]<=wmax)*(omgdata[,1]>=wmin)*omgdata[,1]
#Exit Probabilities
pdata[,1] <- pnorm(WBAR(omega=omgdata[,1], K=lnkdata[,1])%*%parBAR)
#Initial Investment
lnidata[,1] <- rowSums(IX(K=lnkdata[,1], omega=omgdata[,1])*lspline(vectau=vectau, bvec=parI, b1=parIb[1], bL=parIb[2], u=iotadata[,1]))
#Restrict Support of Initial Investment
lnidata[,1] <- (lnidata[,1]>max(t1data$I))*max(t1data$I)+(lnidata[,1]<min(t1data$I))*min(t1data$I)+(lnidata[,1]<=max(t1data$I))*(lnidata[,1]>=min(t1data$I))*lnidata[,1]
#Evolution of Productivity and Capital
for (t in 2:T){
	ttdata <- US %>% group_by(id) %>% slice(t)
	#Capital
	lnkdata[,t] <- log(0.98*exp(lnkdata[,t-1])+exp(lnidata[,t-1]))
	#Restrict the Support of Capital 
	lnkdata[,t] <- (lnkdata[,t]>max(ttdata$K))*max(ttdata$K)+(lnkdata[,t]<min(ttdata$K))*min(ttdata$K)+(lnkdata[,t]<=max(ttdata$K))*(lnkdata[,t]>=min(ttdata$K))*lnkdata[,t]
	#Productivity
	omgdata[,t] <- rowSums(WX(omega=omgdata[,t-1])*lspline(vectau=vectau, bvec=parWT, b1=parWTb[1], bL=parWTb[2], u=xidata[,t]))
	#Restrict the Support of Productivity
	omgdata[,t] <- (omgdata[,t]>wmax)*wmax+(omgdata[,t]<wmin)*wmin+(omgdata[,t]<=wmax)*(omgdata[,t]>=wmin)*omgdata[,t]
	#Exit Probabilities
	pdata[,t] <- pnorm(WBAR(omega=omgdata[,t-1], K=lnkdata[,t])%*%parBAR)
	#Investment
	lnidata[,t] <- rowSums(IX(K=lnkdata[,t], omega=omgdata[,t])*lspline(vectau=vectau, bvec=parI, b1=parIb[1], bL=parIb[2], u=iotadata[,t]))
	#Restrict the Support of Investment
	lnidata[,t] <- (lnidata[,t]>max(ttdata$I))*max(ttdata$I)+(lnidata[,t]<min(ttdata$I))*min(ttdata$I)+(lnidata[,t]<=max(ttdata$I))*(lnidata[,t]>=min(ttdata$I))*lnidata[,t]
}
for (t in 1:T){
	ttdata <- US %>% group_by(id) %>% slice(t)
	#Labor
	lnldata[,t] <- rowSums(LX(K=lnkdata[,t], omega=omgdata[,t])*lspline(vectau=vectau, bvec=parL, b1=parLb[1], bL=parLb[2], u=epsdata[,t]))
	#Restrict Support of Labor
	lnldata[,t] <- (lnldata[,t]>max(ttdata$L))*max(ttdata$L)+(lnldata[,t]<min(ttdata$L))*min(ttdata$L)+(lnldata[,t]<=max(ttdata$L))*(lnldata[,t]>=min(ttdata$L))*lnldata[,t]
	#Materials
	lnmdata[,t] <- rowSums(MX(K=lnkdata[,t], L=lnldata[,t], omega=omgdata[,t])*lspline(vectau=vectau, bvec=parM, b1=parMb[1], bL=parMb[2], u=varepsdata[,t]))
	#Restrict Support of Materials
	lnmdata[,t] <- (lnmdata[,t]>max(ttdata$M))*max(ttdata$M)+(lnmdata[,t]<min(ttdata$M))*min(ttdata$M)+(lnmdata[,t]<=max(ttdata$M))*(lnmdata[,t]>=min(ttdata$M))*lnmdata[,t]
	#Output
	lnydata[,t] <- rowSums(PF(K=lnkdata[,t], L=lnldata[,t], M=lnmdata[,t], omega=omgdata[,t], method=method)*lspline(vectau=vectau, bvec=parY, b1=parYb[1], bL=parYb[2], u=etadata[,t]))
	#Restrict Support of Output
	lnydata[,t] <- (lnydata[,t]>max(ttdata$Y))*max(ttdata$Y)+(lnydata[,t]<min(ttdata$Y))*min(ttdata$Y)+(lnydata[,t]<=max(ttdata$Y))*(lnydata[,t]>=min(ttdata$Y))*lnydata[,t]
}
# #Matrix of logical for exit (TRUE) or stay (False)
# pindmat <- pdata>xidata
# #Logical for which firm enters/exits
# pindi <- apply(pindmat, 1, function(x) any(x==TRUE))
# #Matrix of logical for ALL firms that exit (which year)
# pmat <- pindmat[pindi,]
# #Which year the firm exits
# pind <- apply(pmat, 1, function(x) which(x==TRUE)[1])
# for (i in 1:nrow(pmat)){
# 	pmat[i,(pind[i]:T)] <- TRUE 
# }
# pindmat[pindi,] <- pmat
#Create new data frame
data <- data.frame(id=rep(seq(1:N), each=T), time=rep(seq(min(US$time), max(US$time)), N), Y=c(t(lnydata)), 
	K=c(t(lnkdata)), L=c(t(lnldata)), M=c(t(lnmdata)), I=c(t(lnidata)), W=c(t(omgdata)), 
	xi=c(t(xidata)), epsl=c(t(epsdata)), epsm=c(t(varepsdata)), iota=c(t(iotadata)))
#Select Firms up until exit
# data <- exdata[!exdata$exitind,]
idcon <- duplicated(data$id)
idlag <- duplicated(data$id, fromLast=TRUE)
#Vectorize
#Output
out <- data$Y
outcon <- out[idcon]
outlag <- out[idlag]
#Capital
cap <- data$K
capcon <- cap[idcon]
caplag <- cap[idlag]
#Investment
inv <- data$I
invcon <- inv[idcon]
invlag <- inv[idlag]
iota <- data$iota
#Labor
lab <- data$L
labcon <- lab[idcon]
lablag <- lab[idlag]
epsL <- data$epsl
#Materials
mat <- data$M
matcon <- mat[idcon]
matlag <- mat[idlag]
epsM <- data$epsm
#Productivity
omg <- data$W
omgcon <- omg[idcon]
omglag <- omg[idlag]
xi <- data$xi
xicon <- xi[idcon]
#Capital
kpost <- c(3, 6, 8, 9, 12, 15, 17, 18)
#Labor
lpost <- c(4, 6, 7, 10, 13, 15, 16, 19)
#Materials
mpost <- c(5, 7, 8, 11, 14, 16, 17, 20)
#Hicks-Capital
hkpost <- c(12, 15, 17, 18)
#Hicks-Labor
hlpost <- c(13, 15, 16, 19)
#Hicks-Materials
hmpost <- c(14, 16, 17, 20)
#Conditional Skewness, Dispersion, and Kurtosis evaluated at percentiles of productivity
skomega <- array(0, c(ntau))
dispomega <- array(0, c(ntau))
kurpomega <- array(0, c(ntau))
for (q in 1:ntau){
	omgqlag <- as.numeric(quantile(omglag, probs=vectau[q]))
	skomega[q] <-  (WX(omega=omgqlag)%*%parWT[,ntau]+WX(omega=omgqlag)%*%parWT[,1]-2*WX(omega=omgqlag)%*%parWT[,6])/(WX(omega=omgqlag)%*%parWT[,ntau]-WX(omega=omgqlag)%*%parWT[,1])
	dispomega[q] <- WX(omega=omgqlag)%*%parWT[,ntau]-WX(omega=omgqlag)%*%parWT[,1]
	kurpomega[q] <- (WX(omega=omgqlag)%*%parWT[,(ntau-1)]-WX(omega=omgqlag)%*%parWT[,1])/(WX(omega=omgqlag)%*%parWT[,(ntau-2)]-WX(omega=omgqlag)%*%parWT[,2])
}
displot <- plot_ly(data.frame(x=vectau, y=dispomega), x=~x, y=~y, type="scatter", mode="lines", line=list(color="black"), showlegend=F, name="", hovertemplate = paste("<i>𝛕-productivity<i>: %{x}", "<br>Conditional Dispersion: %{y:.3f}")) %>% layout(xaxis=list(title="𝛕-productivity", titlefont=list(size=30), tickfont=list(size=25)), yaxis=list(title="Conditional Dispersion", titlefont=list(size=18), tickfont=list(size=25)))
skewplot <- plot_ly(data.frame(x=vectau, y=skomega), x=~x, y=~y, type="scatter", mode="lines", line=list(color="black"), showlegend=F, name="", hovertemplate = paste("<i>𝛕-productivity<i>: %{x}", "<br>Conditional Skewness: %{y:.3f}")) %>% layout(xaxis=list(title="𝛕-productivity", titlefont=list(size=30), tickfont=list(size=25)), yaxis=list(title="Conditional Skewness", titlefont=list(size=18), tickfont=list(size=25)))
kurplot <- plot_ly(data.frame(x=vectau, y=kurpomega), x=~x, y=~y, type="scatter", mode="lines", line=list(color="black"), showlegend=F, name="", hovertemplate = paste("<i>𝛕-productivity<i>: %{x}", "<br>Conditional Kurtosis: %{y:.3f}")) %>% layout(xaxis=list(title="𝛕-productivity", titlefont=list(size=30), tickfont=list(size=25)), yaxis=list(title="Conditional Kurtosis", titlefont=list(size=18), tickfont=list(size=25)))
#Annotate
latexannotationsdist <- list(list(x=0.05, y=0.95, text="(a) Conditional Dispersion", font=list(size=35, family="Times New Roman"), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
	list(x=0.5, y=0.95, text="(b) Conditional Skewness", font=list(size=35, family="Times New Roman"), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
	list(x=0.95, y=0.95, text="(c) Conditional Kurtosis", font=list(size=35, family="Times New Roman"), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE))
annotationsdist <- list(list(x=0.05, y=0.95, text="(a) Conditional Dispersion", font=list(size=18, family="Times New Roman"), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
	list(x=0.5, y=0.95, text="(b) Conditional Skewness", font=list(size=18, family="Times New Roman"), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
	list(x=0.95, y=0.95, text="(c) Conditional Kurtosis", font=list(size=18, family="Times New Roman"), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE))
omgdist <- subplot(displot, skewplot, kurplot, nrows=1, shareX=TRUE) %>% layout(scene=list(aspectratio=list(x=.6, y=.6, z=.6), xaxis=list(title="𝛕-productivity", titlefont=list(size=25), tickfont=list(size=25)), yaxis=list(title="Conditional Dispersion", titlefont=list(size=18), tickfont=list(size=20))), 
	scene2=list(aspectratio=list(x=.6, y=.6, z=.6), xaxis=list(title="𝛕-productivity", titlefont=list(size=25), tickfont=list(size=25)), yaxis=list(title="Conditional Skewness", titlefont=list(size=18), tickfont=list(size=25))),
	scene3=list(aspectratio=list(x=.6, y=.6, z=.6), xaxis=list(title="𝛕-productivity", titlefont=list(size=25), tickfont=list(size=25)), yaxis=list(title="Conditional Kurtosis", titlefont=list(size=18), tickfont=list(size=25))))
latexomgdist <- omgdist %>% layout(annotations=latexannotationsdist)
latexomgdist
omgdist <- omgdist %>% layout(annotations=annotationsdist)
omgdist   <- plotly_json(omgdist, FALSE)
write(omgdist, "/Users/justindoty/Documents/Home/My_Website/static/jmp/selection/omgdist.json")
########################################################################################
#Individual Quantile Marginal Effects####################################################
#########################################################################################
#Commands for Colors
nrz <- length(vectau)
ncz <- length(vectau)
jet.colors <-  colorRampPalette(c("midnightblue", "blue", "cyan","green", "yellow","orange","red", "darkred"))
nbcol <- 64
color <- jet.colors(nbcol)
#Plots
#For output elasticities that vary over percentiles of inputs
k3d <- array(0, c(ntau, ntau))
l3d <- k3d; m3d <- k3d
#For output elasticities that vary over productivity
kwq3d <- k3d; lwq3d <- k3d; mwq3d <- k3d
#For input demand functions
iw3d <- k3d; lw3d <- k3d; mw3d <- k3d
hk3d <- array(0, c(ntau, ntau))
hl3d <- hk3d; hm3d <- hk3d
for (q in 1:ntau){
	#For fixed quantiles of capital
	capkq <- rep(as.numeric(quantile(cap, probs=vectau[q])), length(cap))
	wmean <- rep(mean(omg), length(omg))
	labk <- rowSums(LX(K=capkq, omega=wmean)*lspline(vectau=vectau, bvec=parL, b1=parLb[1], bL=parLb[2], u=epsL))
	matk <- rowSums(MX(K=capkq, L=mean(labk), omega=wmean)*lspline(vectau=vectau, bvec=parM, b1=parMb[1], bL=parMb[2], u=epsM))
	# #For fixed quantiles of labor
	kmean <- rep(mean(cap), length(cap))
	lablq <- quantile(rowSums(LX(K=kmean, omega=wmean)*lspline(vectau=vectau, bvec=parL, b1=parLb[1], bL=parLb[2], u=epsL)), probs=vectau[q])
	matl <- rowSums(MX(K=kmean, L=lablq, omega=wmean)*lspline(vectau=vectau, bvec=parM, b1=parMb[1], bL=parMb[2], u=epsM))
	# #For fixed quantiles of materials
	labm <- rowSums(LX(K=kmean, omega=wmean)*lspline(vectau=vectau, bvec=parL, b1=parLb[1], bL=parLb[2], u=epsL))
	matmq <- quantile(rowSums(MX(K=kmean, L=mean(labm), omega=wmean)*lspline(vectau=vectau, bvec=parM, b1=parMb[1], bL=parMb[2], u=epsM)), probs=vectau[q])
	#For fixed quantiles of productivity
	wq <- rep(as.numeric(quantile(omg, probs=vectau[q])), length(omg))
	iwq <- rowSums(IX(K=kmean, omega=wq)*lspline(vectau=vectau, bvec=parI, b1=parIb[1], bL=parIb[2], u=iota))
	labw <- rowSums(LX(K=kmean, omega=wq)*lspline(vectau=vectau, bvec=parL, b1=parLb[1], bL=parLb[2], u=epsL))
	matw <- rowSums(MX(K=kmean, L=mean(labw), omega=wq)*lspline(vectau=vectau, bvec=parM, b1=parMb[1], bL=parMb[2], u=epsM))
	for (qq in 1:ntau){
		#Production Elasticities
		k3d[q,] <- colMeans(cbind(1, labk, matk, 2*capkq, wmean, wmean*labk, wmean*matk, 2*capkq*wmean))%*%parY[kpost,]
		l3d[q,] <- colMeans(cbind(1, kmean, matl, 2*lablq, wmean, wmean*kmean, wmean*matl, 2*lablq*wmean))%*%parY[lpost,]
		m3d[q,] <- colMeans(cbind(1, labm, kmean, 2*matmq, wmean, wmean*labm, wmean*kmean, 2*matmq*wmean))%*%parY[mpost,]
		#Production Elasticities (fixed percentiles of productivity)
		kwq3d[q,] <- colMeans(cbind(1, labw, matw, 2*kmean, wq, wq*labw, wq*matw, 2*kmean*wq))%*%parY[kpost,]
		lwq3d[q,] <- colMeans(cbind(1, kmean, matw, 2*labw, wq, wq*kmean, wq*matw, 2*labw*wq))%*%parY[lpost,]
		mwq3d[q,] <- colMeans(cbind(1, labw, kmean, 2*matw, wq, wq*labw, wq*kmean, 2*matw*wq))%*%parY[mpost,]
		#Production Efficiencies
		hk3d[q,] <- colMeans(cbind(1, labk, matk, 2*capkq))%*%parY[hkpost,]
		hl3d[q,] <- colMeans(cbind(1, kmean, matl, 2*lablq))%*%parY[hlpost,]
		hm3d[q,] <- colMeans(cbind(1, labm, kmean, 2*matmq))%*%parY[hmpost,]
		#Marginal Productivities
		iw3d[q,] <- colMeans(IXD(K=kmean, omega=wq, par=parI, pos=2, sdpos=1))
		lw3d[q,] <- colMeans(LXD(K=kmean, omega=wq, par=parL, pos=2, sdpos=1))
		mw3d[q,] <- colMeans(MXD(K=kmean, L=labw, omega=wq, par=parM, pos=3, sdpos=1))
	}

}
###################################################
#Capital Estimates
##################################################
#Elasticity (over capital)
kplot <- plot_ly(x=vectau, y=vectau, z=k3d, colorscale="Jet", type="surface", showscale=FALSE, scene="scene1", name=" ", hovertemplate = paste("<i>𝛕-output<i>: %{x:.2f}", "<br>𝛕-capital: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% layout(scene1=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-output"), yaxis=list(title="𝛕-capital"), zaxis=list(title="Capital Elasticity")))
# k3dplotly
# Elasticity (over productivity)
kwqplot <- plot_ly(x=vectau, y=vectau, z=kwq3d, colorscale="Jet", type="surface", showscale=FALSE, scene="scene1", name=" ", hovertemplate = paste("<i>𝛕-output<i>: %{x:.2f}", "<br>𝛕-productivity: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% layout(scene1=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-output"), yaxis=list(title="𝛕-productivity"), zaxis=list(title="Capital Elasticity")))
# kwq3dplotly
#Efficiency (over capital)
hkplot <- plot_ly(x=vectau, y=vectau, z=hk3d, colorscale="Jet", type="surface", showscale=FALSE, scene="scene1", name=" ", hovertemplate = paste("<i>𝛕-output<i>: %{x:.2f}", "<br>𝛕-capital: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% layout(scene1=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-output"), yaxis=list(title="𝛕-capital"), zaxis=list(title="Capital Efficiency"))) 
# hk3dplotly
#Investment Response to productivity (over productivity)
iwplot <- plot_ly(x=vectau, y=vectau, z=iw3d, colorscale="Jet", type="surface", showscale=FALSE, scene="scene1", name=" ", hovertemplate = paste("<i>𝛕-investment<i>: %{x:.2f}", "<br>𝛕-productivity: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% layout(scene1=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-investment"), yaxis=list(title="𝛕-productivity"), zaxis=list(title="Investment Productivity"))) 
# iw3dplotly
# ##########################################################
# #Labor Estimates
##########################################################
#Elasticity (over labor)
lplot <- plot_ly(x=vectau, y=vectau, z=l3d, colorscale="Jet", type="surface", showscale=FALSE, scene="scene2", name=" ", hovertemplate = paste("<i>𝛕-output<i>: %{x:.2f}", "<br>𝛕-labor: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% layout(scene2=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-output"), yaxis=list(title="𝛕-labor"), zaxis=list(title="Labor Elasticity")))  
# l3dplotly 
#Elasticity (over productivity)
lwqplot <- plot_ly(x=vectau, y=vectau, z=lwq3d, colorscale="Jet", type="surface", showscale=FALSE, scene="scene2", name=" ", hovertemplate = paste("<i>𝛕-output<i>: %{x:.2f}", "<br>𝛕-productivity: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% layout(scene2=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-output"), yaxis=list(title="𝛕-productivity"), zaxis=list(title="Labor Elasticity")))  
# lwq3dplotly
#Efficiency (over labor)
hlplot <- plot_ly(x=vectau, y=vectau, z=hl3d, colorscale="Jet", type="surface", showscale=FALSE, scene="scene2", name=" ", hovertemplate = paste("<i>𝛕-output<i>: %{x:.2f}", "<br>𝛕-labor: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% layout(scene2=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-output"), yaxis=list(title="𝛕-labor"), zaxis=list(title="Labor Efficiency"))) 
# hl3dplotly
#Labor response to productivity (over productivity)
lwplot <- plot_ly(x=vectau, y=vectau, z=lw3d, colorscale="Jet", type="surface", showscale=FALSE, scene="scene2", name=" ", hovertemplate = paste("<i>𝛕-labor<i>: %{x:.2f}", "<br>𝛕-productivity: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% layout(scene2=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-labor"), yaxis=list(title="𝛕-productivity"), zaxis=list(title="Labor Productivity"))) 
# lw3dplotly
#########################################################
#Materials Estimates
#########################################################
#Elasticity (over materials)
mplot <- plot_ly(x=vectau, y=vectau, z=m3d, colorscale="Jet", type="surface", showscale=FALSE, scene="scene3", name=" ", hovertemplate = paste("<i>𝛕-output<i>: %{x:.2f}", "<br>𝛕-materials: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% layout(scene3=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-output"), yaxis=list(title="𝛕-materials"), zaxis=list(title="Materials Elasticity"))) 
# m3dplotly 
#Elasticity (over productivity)
mwqplot <- plot_ly(x=vectau, y=vectau, z=mwq3d, colorscale="Jet", type="surface", showscale=FALSE, scene="scene3", name=" ", hovertemplate = paste("<i>𝛕-output<i>: %{x:.2f}", "<br>𝛕-productivity: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% layout(scene3=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-output"), yaxis=list(title="𝛕-productivity"), zaxis=list(title="Materials Elasticity"))) 
# mwq3dplotly 
#Efficiency (over materials)
hmplot <- plot_ly(x=vectau, y=vectau, z=hm3d, colorscale="Jet", type="surface", showscale=FALSE, scene="scene3", name=" ", hovertemplate = paste("<i>𝛕-output<i>: %{x:.2f}", "<br>𝛕-materials: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% layout(scene3=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-output"), yaxis=list(title="𝛕-materials"), zaxis=list(title="Materials Efficiency"))) 
# hm3dplotly 
#Material response to productivity (over productivity)
mwplot <- plot_ly(x=vectau, y=vectau, z=mw3d, colorscale="Jet", type="surface", showscale=FALSE, scene="scene3", name=" ", hovertemplate = paste("<i>𝛕-materials<i>: %{x:.2f}", "<br>𝛕-productivity: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% layout(scene3=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-materials"), yaxis=list(title="𝛕-productivity"), zaxis=list(title="Materials Productivity"))) 
# mw3dplotly
#Combined Plot.ly Elasticities
#Elasticities (over percentiles of inputs)
klmplot <- subplot(kplot, lplot, mplot, shareX=TRUE) %>% layout(scene=list(aspectratio=list(x=.6, y=.6, z=.6), xaxis=list(title="𝛕-output", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="𝛕-capital", titlefont=list(size=18), tickfont=list(size=14)), zaxis=list(title="Capital", titlefont=list(size=18), tickfont=list(size=14))), 
	scene2=list(aspectratio=list(x=.6, y=.6, z=.6), xaxis=list(title="𝛕-output", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="𝛕-labor", titlefont=list(size=18), tickfont=list(size=14)), zaxis=list(title="Labor", titlefont=list(size=18), tickfont=list(size=14))),
	scene3=list(aspectratio=list(x=.6, y=.6, z=.6), xaxis=list(title="𝛕-output", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="𝛕-materials", titlefont=list(size=18), tickfont=list(size=14)), zaxis=list(title="Materials", titlefont=list(size=18), tickfont=list(size=14))))
# Elasticities (over percentiles of productivity)
klmwqplot <- subplot(kwqplot, lwqplot, mwqplot, shareX=TRUE) %>% layout(scene=list(aspectratio=list(x=.6, y=.6, z=.6), xaxis=list(title="𝛕-output", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="𝛕-productivity", titlefont=list(size=18), tickfont=list(size=14)), zaxis=list(title="Capital", titlefont=list(size=18), tickfont=list(size=14))), 
	scene2=list(aspectratio=list(x=.6, y=.6, z=.6), xaxis=list(title="𝛕-output", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="𝛕-productivity", titlefont=list(size=18), tickfont=list(size=14)), zaxis=list(title="Labor", titlefont=list(size=18), tickfont=list(size=14))),
	scene3=list(aspectratio=list(x=.6, y=.6, z=.6), xaxis=list(title="𝛕-output", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="𝛕-productivity", titlefont=list(size=18), tickfont=list(size=14)), zaxis=list(title="Materials", titlefont=list(size=18), tickfont=list(size=14))))
#Annotations
#I use two different text sizes (one for latex and one for my website)
latexannotationsklm <- list(list(x=0.09, y=0.75, text="(a) Capital Elasticity", font=list(size=30, family="Times New Roman"), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
	list(x=0.5, y=0.75, text="(b) Labor Elasticity", font=list(size=30, family="Times New Roman"), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
	list(x=0.95, y=0.75, text="(c) Materials Elasticity", font=list(size=30, family="Times New Roman"), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE))
annotationsklm <- list(list(x=0.1, y=0.75, text="(a) Capital Elasticity", font=list(size=18), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
	list(x=0.5, y=0.75, text="(b) Labor Elasticity", font=list(size=18), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
	list(x=0.9, y=0.75, text="(c) Materials Elasticity", font=list(size=18), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE))
# Add
klmplot_latex <- klmplot %>% layout(annotations=latexannotationsklm)
klmwqplot_latex <- klmwqplot %>% layout(annotations=latexannotationsklm)
#Plot
klmplot_latex 
klmwqplot_latex 
# Save to JSON
klmplot <- klmplot %>% layout(annotations=annotationsklm)
klmwqplot <- klmwqplot %>% layout(annotations=annotationsklm)
#Json
klmplot   <- plotly_json(klmplot, FALSE)
klmwqplot  <- plotly_json(klmwqplot, FALSE)
write(klmplot, "/Users/justindoty/Documents/Home/My_Website/static/jmp/selection/klm3dplotly.json")
write(klmwqplot , "/Users/justindoty/Documents/Home/My_Website/static/jmp/selection/klmwq3dplotly.json")
#Combined Plot.ly Efficiecies
hklmplot <- subplot(hkplot, hlplot, hmplot, shareX=TRUE) %>% layout(scene=list(aspectratio=list(x=.6, y=.6, z=.6), xaxis=list(title="𝛕-output", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="𝛕-capital", titlefont=list(size=18), tickfont=list(size=14)), zaxis=list(title="Capital", titlefont=list(size=18), tickfont=list(size=14))), 
	scene2=list(aspectratio=list(x=.6, y=.6, z=.6), xaxis=list(title="𝛕-output", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="𝛕-labor", titlefont=list(size=18), tickfont=list(size=14)), zaxis=list(title="Labor", titlefont=list(size=18), tickfont=list(size=14))),
	scene3=list(aspectratio=list(x=.6, y=.6, z=.6), xaxis=list(title="𝛕-output", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="𝛕-materials", titlefont=list(size=18), tickfont=list(size=14)), zaxis=list(title="Materials", titlefont=list(size=18), tickfont=list(size=14))))
#Annotations
latexannotationshklm <- list(list(x=0.09, y=0.75, text="(a) Capital Efficiency", font=list(size=30, family="Times New Roman"), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
	list(x=0.5, y=0.75, text="(b) Labor Efficiency", font=list(size=30, family="Times New Roman"), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
	list(x=0.95, y=0.75, text="(c) Materials Efficiency", font=list(size=30, family="Times New Roman"), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE))
annotationshklm <- list(list(x=0.1, y=0.75, text="(a) Capital Efficiency", font=list(size=18), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
	list(x=0.5, y=0.75, text="(b) Labor Efficiency", font=list(size=18), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
	list(x=0.9, y=0.75, text="(c) Materials Efficiency", font=list(size=18), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE))

#Add
hklmplot_latex <- hklmplot %>% layout(annotations=latexannotationshklm)
#Plot
hklmplot_latex 
#JSON
hklmplot <- hklmplot %>% layout(annotations=annotationshklm)
#Save
hklmplot   <- plotly_json(hklmplot, FALSE)
write(hklmplot, "/Users/justindoty/Documents/Home/My_Website/static/jmp/selection/hklm3dplotly.json")
#Combined Productivities
lmiwplot <- subplot(iwplot, lwplot, mwplot) %>% layout(scene=list(aspectratio=list(x=0.6, y=0.6, z=0.6), xaxis=list(title="𝛕-investment", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="𝛕-productivity", titlefont=list(size=18), tickfont=list(size=14)), zaxis=list(title="Investment", titlefont=list(size=18), tickfont=list(size=14))), 
	scene2=list(aspectratio=list(x=0.6, y=0.6, z=0.6), xaxis=list(title="𝛕-labor", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="𝛕-productivity", titlefont=list(size=18), tickfont=list(size=14)), zaxis=list(title="Labor", titlefont=list(size=18), tickfont=list(size=14))),
	scene3=list(aspectratio=list(x=0.6, y=0.6, z=0.6), xaxis=list(title="𝛕-materials", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="𝛕-productivity", titlefont=list(size=18), tickfont=list(size=14)), zaxis=list(title="Materials", titlefont=list(size=18), tickfont=list(size=14))))
#Annotations
latexannotationsw <- list(list(x=0.08, y=0.75, text="(a) Investment Response", font=list(size=30, family="Times New Roman"), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
	list(x=0.5, y=0.75, text="(b) Labor Response", font=list(size=30, family="Times New Roman"), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
	list(x=0.92, y=0.75, text="(c) Materials Response", font=list(size=30, family="Times New Roman"), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE))
annotationsw <- list(list(x=0.07, y=0.75, text="(a) Investment Response", font=list(size=18), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
	list(x=0.5, y=0.75, text="(b) Labor Response", font=list(size=18), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
	list(x=0.92, y=0.75, text="(c) Materials Response", font=list(size=18), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE))
#Add
lmiwplot_latex <- lmiwplot %>% layout(annotations=latexannotationsw)
#Plot
lmiwplot_latex
#JSON
lmiwplot <- lmiwplot %>% layout(annotations=annotationsw)
#Save
lmiwplot   <- plotly_json(lmiwplot, FALSE)
write(lmiwplot, "/Users/justindoty/Documents/Home/My_Website/static/jmp/selection/lmiw3dplotly.json")
######################################################################################################################################################################################
#Productivity
#####################################################################################################################################################################################
omg3dq <- array(0, c(ntau, ntau))
for (q in 1:ntau){
	#Persistence of Productivity
	#Fix quantiles of previous period productivity
	omgq <- rep(as.numeric(quantile(omglag, probs=vectau[q])), length(omglag))
	#Generate productivity at fixed levels of productivity for Non R&D Firms
	omgx <- rowSums(WX(omega=omgq)*lspline(vectau=vectau, bvec=parWT, b1=parWTb[1], bL=parWTb[2], u=xicon))
	omg3dq[q,] <- colMeans(sweep(WX(omega=omgx)[,-c(dims$W)], 2, c(1:3), "*"))%*%parWT[-1,]
	
}
#Productivity Persistence
omgplotly <- plot_ly(x=vectau, y=vectau, z=omg3dq, colorscale="Jet", type="surface", showscale=FALSE, scene="scene1", name=" ", hovertemplate = paste("<i>𝛕-innovation<i>: %{x:.2f}", "<br>𝛕-productivity: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% layout(scene1=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-innovation", titlefont=list(size=20), tickfont=list(size=16)), yaxis=list(title="𝛕-productivity", titlefont=list(size=20), tickfont=list(size=16)), zaxis=list(title="Persistence", titlefont=list(size=20), tickfont=list(size=16)))) 
omgplotly
omgplotly <- plotly_json(omgplotly, FALSE)
write(omgplotly, "/Users/justindoty/Documents/Home/My_Website/static/jmp/selection/omgplotly.json")