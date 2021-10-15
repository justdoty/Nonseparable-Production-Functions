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
library(listviewer)
source('/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Functions/Main/Auxfuns.R')
source('/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Functions/Main/Tensors.R')
source('/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Functions/Main/omega.R')
#Load US Dataset
US <- read.csv('/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Data/USdata.csv') %>% 
select(id, year, lny, lnk1, lnk2, lnl, lnm, lni, age) %>% transmute(id=id, year=year, Y=lny, K=lnk1, L=lnl, M=lnm, I=lni, A=age)
#Productivity from LP
omegainit <- omega_est(idvar=US$id, timevar=US$year, Y=US$Y, A=US$A, K=US$K, L=US$L, M=US$M)$omega
#Add LP Estimate of productivity
US <- US %>% mutate(omega=omegainit) 
idcon <- duplicated(US$id)
idlag <- duplicated(US$id, fromLast=TRUE)
id1 <- !idcon
########################################################################################################
##########################################Summary Statistics############################################
########################################################################################################
#Create table for all relevant summary statistics
sumALL <- summarise_at(US, c("Y", "K", "L", "M", "I"), list(Q1=~quantile(., 0.25), med=median, Q3=~quantile(.,0.75), mean=mean, sd=sd), na.rm=TRUE)
sumstat <- matrix(round(sumALL, 2), nrow=5, ncol=5)
summary_table <- cbind(c("Output", "Capital", "Labor", "Materials", "Investment"), sumstat)
colnames(summary_table) <- c(" ", "1st Qu.", 'Median', "3rd Qu.", 'Mean', "sd")

summary_table <- xtable(summary_table, digits=c(0,0,2,2,2,2,2), type="latex")
align(summary_table) <- rep('c', 7)
addtorow <- list()
addtorow$pos <- list(-1)
addtorow$command <- '\\hline\\hline '
#For copy pasting into latex
print(summary_table, hline.after=c(0,nrow(summary_table)), add.to.row=addtorow, auto=FALSE, include.rownames=FALSE, sanitize.text.function=function(x) x, table.placement="H")
#Saves to file
# print(summary_table, hline.after=c(0,nrow(summary_table)), add.to.row=addtorow, auto=FALSE, include.rownames=FALSE, sanitize.text.function=function(x) x, table.placement="H", file="/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Empirical/US/Results/Translog/Tex/US_Summary.tex")
#Average number of firms per year
avg_firms <- round(as.numeric(colMeans(group_by(US, year) %>% summarise(firms=n()))[2]))
print(avg_firms)
unique_firms <- length(unique(US$id))
print(unique_firms)
########################################################################################################
##########################################Load Results############################################
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
#We simuluate a balanced panel even though the original model is estimated from
#an unbalanced panel. In a later version, we consider adding a selection bias
#correction to the productivity equation and drop firms according to this rule
#in the simulated model
#############################################################################
Nsim <- 5
N <- length(unique(US$id))*Nsim
T <- length(unique(US$year))
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
#Simulate Data at T=1
t1data <- US %>% group_by(id) %>% slice(1)
#########################################################################################################
#Capital is not estimated in the main model: Use Capital Accumulation Process
#For K=1
lnkdata[,1] <- kronecker(array(1, c(Nsim,1)), t1data$K)
#Initial Productivity
omgdata[,1] <- rowSums(WX1(K=lnkdata[,1], L=lnldata[,1], M=lnmdata[,1])*lspline(vectau=vectau, bvec=parW1, b1=parW1b[1], bL=parW1b[2], u=xidata[,1]))
#Restrict the Support of Initial Productivity
omgdata[,1] <- (omgdata[,1]>wmax)*wmax+(omgdata[,1]<wmin)*wmin+(omgdata[,1]<=wmax)*(omgdata[,1]>=wmin)*omgdata[,1]
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
#Vectorize
#Output
out <- c(t(lnydata))
outcon <- c(t(lnydata[,2:T]))
outlag <- c(t(lnydata[,1:(T-1)]))
#Capital
cap <- c(t(lnkdata))
capcon <- c(t(lnkdata[,2:T]))
caplag <- c(t(lnkdata[,1:(T-1)]))
#Investment
inv <- c(t(lnidata))
invcon <- c(t(lnidata[,2:T]))
invlag <- c(t(lnidata[,1:(T-1)]))
iota <- c(t(iotadata))
#Labor
lab <- c(t(lnldata))
labcon <- c(t(lnldata[,2:T]))
lablag <- c(t(lnldata[,1:(T-1)]))
epsL <- c(t(epsdata))
#Materials
mat <- c(t(lnmdata))
matcon <- c(t(lnmdata[,2:T]))
matlag <- c(t(lnmdata[,1:(T-1)]))
epsM <- c(t(varepsdata))
#Productivity
omg <- c(t(omgdata))
omgcon <- c(t(omgdata[,2:T]))
omglag <- c(t(omgdata[,1:(T-1)]))
xi <- c(t(xidata))
xicon <- c(t(xidata[,2:T]))
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
#For output elasticities that vary over percentiles of other inputs
klq3d <- k3d; kmq3d <- k3d
lkq3d <- k3d; lmq3d <- k3d
mkq3d <- k3d; mlq3d <- k3d
#For input demand functions varying over percentiles of productivity
iw3d <- k3d; lw3d <- k3d; mw3d <- k3d
#For input demand functions varying over percentiles of capital
iwkq3d <- k3d; lwkq3d <- k3d; mwkq3d <- k3d
#For non-Hicks neutral effects varying over percentiles of inputs
hk3d <- array(0, c(ntau, ntau))
hl3d <- hk3d; hm3d <- hk3d
#For non-Hicks neutral effects that vary over percentiles of other inputs
hklq3d <- k3d; hkmq3d <- k3d
hlkq3d <- k3d; hlmq3d <- k3d
hmkq3d <- k3d; hmlq3d <- k3d
for (q in 1:ntau){
	#For fixed quantiles of capital
	capkq <- rep(as.numeric(quantile(cap, probs=vectau[q])), length(cap))
	#For fixed quantiles of productivity
	wq <- rep(as.numeric(quantile(omg, probs=vectau[q])), length(omg))
	#Mean of Capital
	kmean <- rep(mean(cap), length(cap))
	#Mean of productivity
	wmean <- rep(mean(omg), length(omg))
	#Labor at fixed percentiles of capital
	labk <- rowSums(LX(K=capkq, omega=wmean)*lspline(vectau=vectau, bvec=parL, b1=parLb[1], bL=parLb[2], u=epsL))
	#Materials at fixed percentiles of capital and average labor
	matk <- rowSums(MX(K=capkq, L=mean(labk), omega=wmean)*lspline(vectau=vectau, bvec=parM, b1=parMb[1], bL=parMb[2], u=epsM))
	#Fixed quantiles of labor
	lablq <- quantile(rowSums(LX(K=kmean, omega=wmean)*lspline(vectau=vectau, bvec=parL, b1=parLb[1], bL=parLb[2], u=epsL)), probs=vectau[q])
	#Materials at Fixed Percentiles of Labor
	matl <- rowSums(MX(K=kmean, L=lablq, omega=wmean)*lspline(vectau=vectau, bvec=parM, b1=parMb[1], bL=parMb[2], u=epsM))
	#Labor
	labm <- rowSums(LX(K=kmean, omega=wmean)*lspline(vectau=vectau, bvec=parL, b1=parLb[1], bL=parLb[2], u=epsL))
	#Materials
	matmq <- quantile(rowSums(MX(K=kmean, L=mean(labm), omega=wmean)*lspline(vectau=vectau, bvec=parM, b1=parMb[1], bL=parMb[2], u=epsM)), probs=vectau[q])
	#For input productivity effects
	iwq <- rowSums(IX(K=kmean, omega=wq)*lspline(vectau=vectau, bvec=parI, b1=parIb[1], bL=parIb[2], u=iota))
	labw <- rowSums(LX(K=kmean, omega=wq)*lspline(vectau=vectau, bvec=parL, b1=parLb[1], bL=parLb[2], u=epsL))
	matw <- rowSums(MX(K=kmean, L=mean(labw), omega=wq)*lspline(vectau=vectau, bvec=parM, b1=parMb[1], bL=parMb[2], u=epsM))
	for (qq in 1:ntau){
		#Production Elasticities (fixed percentiles of inputs)
		k3d[q,] <- colMeans(cbind(1, labk, matk, 2*capkq, wmean, wmean*labk, wmean*matk, 2*capkq*wmean))%*%parY[kpost,]
		l3d[q,] <- colMeans(cbind(1, kmean, matl, 2*lablq, wmean, wmean*kmean, wmean*matl, 2*lablq*wmean))%*%parY[lpost,]
		m3d[q,] <- colMeans(cbind(1, labm, kmean, 2*matmq, wmean, wmean*labm, wmean*kmean, 2*matmq*wmean))%*%parY[mpost,]
		#Production Elasticities (fixed percentiles of productivity)
		kwq3d[q,] <- colMeans(cbind(1, labw, matw, 2*kmean, wq, wq*labw, wq*matw, 2*kmean*wq))%*%parY[kpost,]
		lwq3d[q,] <- colMeans(cbind(1, kmean, matw, 2*labw, wq, wq*kmean, wq*matw, 2*labw*wq))%*%parY[lpost,]
		mwq3d[q,] <- colMeans(cbind(1, labw, kmean, 2*matw, wq, wq*labw, wq*kmean, 2*matw*wq))%*%parY[mpost,]
		#Production Elasticities (fixed percentiles of other inputs)
		#Capital
		klq3d[q,] <- colMeans(cbind(1, lablq, matl, 2*kmean, wmean, wmean*lablq, wmean*matl, 2*kmean*wmean))%*%parY[kpost,]
		kmq3d[q,] <- colMeans(cbind(1, labm, matmq, 2*kmean, wmean, wmean*labm, wmean*matmq, 2*kmean*wmean))%*%parY[kpost,]
		#Labor
		lkq3d[q,] <- colMeans(cbind(1, capkq, matk, 2*labk, wmean, wmean*capkq, wmean*matk, 2*labk*wmean))%*%parY[lpost,]
		lmq3d[q,] <- colMeans(cbind(1, kmean, matmq, 2*labm, wmean, wmean*kmean, wmean*matmq, 2*labm*wmean))%*%parY[lpost,]
		#Materials
		mkq3d[q,] <- colMeans(cbind(1, labk, capkq, 2*matk, wmean, wmean*labk, wmean*capkq, 2*matk*wmean))%*%parY[mpost,]
		mlq3d[q,] <- colMeans(cbind(1, lablq, kmean, 2*matl, wmean, wmean*matl, wmean*kmean, 2*matl*wmean))%*%parY[mpost,]
		#Production Efficiencies (fixed percentiles of inputs)
		hk3d[q,] <- colMeans(cbind(1, labk, matk, 2*capkq))%*%parY[hkpost,]
		hl3d[q,] <- colMeans(cbind(1, kmean, matl, 2*lablq))%*%parY[hlpost,]
		hm3d[q,] <- colMeans(cbind(1, labm, kmean, 2*matmq))%*%parY[hmpost,]
		#Production Efficiencies (fixed percentiles of other inputs)
		#Capital
		hklq3d[q,] <- colMeans(cbind(1, lablq, matl, 2*kmean))%*%parY[hkpost,]
		hkmq3d[q,] <- colMeans(cbind(1, labm, matmq, 2*kmean))%*%parY[hkpost,]
		#Labor
		hlkq3d[q,] <- colMeans(cbind(1, capkq, matk, 2*labk))%*%parY[hlpost,]
		hlmq3d[q,] <- colMeans(cbind(1, kmean, matmq, 2*labm))%*%parY[hlpost,]
		#Materials
		hmkq3d[q,] <- colMeans(cbind(1, labk, capkq, 2*matk))%*%parY[hmpost,]
		hmlq3d[q,] <- colMeans(cbind(1, lablq, kmean, 2*matl))%*%parY[hmpost,]
		#Marginal Productivities (fixed percentiles of productivity)
		iw3d[q,] <- colMeans(IXD(K=kmean, omega=wq, par=parI, pos=2, sdpos=1))
		lw3d[q,] <- colMeans(LXD(K=kmean, omega=wq, par=parL, pos=2, sdpos=1))
		mw3d[q,] <- colMeans(MXD(K=kmean, L=labw, omega=wq, par=parM, pos=3, sdpos=1))
		#Marginal Productivities (fixed percentiles of capital)
		iwkq3d[q,] <- colMeans(IXD(K=capkq, omega=wmean, par=parI, pos=2, sdpos=1))
		lwkq3d[q,] <- colMeans(LXD(K=capkq, omega=wmean, par=parL, pos=2, sdpos=1))
		mwkq3d[q,] <- colMeans(MXD(K=capkq, L=labk, omega=wmean, par=parM, pos=3, sdpos=1))
	}

}
###################################################
#Capital Estimates
##################################################
#Elasticity (over capital)
k3dplotly <- plot_ly(x=vectau, y=vectau, z=k3d, colorscale="Jet", type="surface", showscale=FALSE, scene="scene1", name=" ", hovertemplate = paste("<i>𝛕-output<i>: %{x:.2f}", "<br>𝛕-capital: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% layout(scene1=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-output"), yaxis=list(title="𝛕-capital"), zaxis=list(title="Capital Elasticity")))
# k3dplotly
# Elasticity (over productivity)
kwq3dplotly <- plot_ly(x=vectau, y=vectau, z=kwq3d, colorscale="Jet", type="surface", showscale=FALSE, scene="scene1", name=" ", hovertemplate = paste("<i>𝛕-output<i>: %{x:.2f}", "<br>𝛕-productivity: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% layout(scene1=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-output"), yaxis=list(title="𝛕-productivity"), zaxis=list(title="Capital Elasticity")))
# kwq3dplotly
# Elasticity (over labor)
klq3dplotly <- plot_ly(x=vectau, y=vectau, z=klq3d, colorscale="Jet", type="surface", showscale=FALSE, scene="scene1", name=" ", hovertemplate = paste("<i>𝛕-output<i>: %{x:.2f}", "<br>𝛕-labor: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% layout(scene1=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-output"), yaxis=list(title="𝛕-labor"), zaxis=list(title="Capital Elasticity")))
# klq3dplotly
# Elasticity (over materials)
kmq3dplotly <- plot_ly(x=vectau, y=vectau, z=kmq3d, colorscale="Jet", type="surface", showscale=FALSE, scene="scene2", name=" ", hovertemplate = paste("<i>𝛕-output<i>: %{x:.2f}", "<br>𝛕-materials: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% layout(scene2=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-output"), yaxis=list(title="𝛕-materials"), zaxis=list(title="Capital Elasticity")))
# kmq3dplotly
#Efficiency (over capital)
hk3dplotly <- plot_ly(x=vectau, y=vectau, z=hk3d, colorscale="Jet", type="surface", showscale=FALSE, scene="scene1", name=" ", hovertemplate = paste("<i>𝛕-output<i>: %{x:.2f}", "<br>𝛕-capital: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% layout(scene1=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-output"), yaxis=list(title="𝛕-capital"), zaxis=list(title="Capital Efficiency"))) 
# hk3dplotly
#Efficiency (over labor)
hklq3dplotly <- plot_ly(x=vectau, y=vectau, z=hklq3d, colorscale="Jet", type="surface", showscale=FALSE, scene="scene1", name=" ", hovertemplate = paste("<i>𝛕-output<i>: %{x:.2f}", "<br>𝛕-labor: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% layout(scene1=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-output"), yaxis=list(title="𝛕-labor"), zaxis=list(title="Capital Efficiency"))) 
# hklq3dplotly
#Efficiency (over materials)
hkmq3dplotly <- plot_ly(x=vectau, y=vectau, z=hkmq3d, colorscale="Jet", type="surface", showscale=FALSE, scene="scene2", name=" ", hovertemplate = paste("<i>𝛕-output<i>: %{x:.2f}", "<br>𝛕-materials: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% layout(scene2=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-output"), yaxis=list(title="𝛕-materials"), zaxis=list(title="Capital Efficiency"))) 
# hklq3dplotly
#Investment Response to productivity (over productivity)
iw3dplotly <- plot_ly(x=vectau, y=vectau, z=iw3d, colorscale="Jet", type="surface", showscale=FALSE, scene="scene1", name=" ", hovertemplate = paste("<i>𝛕-investment<i>: %{x:.2f}", "<br>𝛕-productivity: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% layout(scene1=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-investment"), yaxis=list(title="𝛕-productivity"), zaxis=list(title="Investment Productivity"))) 
# iw3dplotly
#Investment Response to productivity (over capital)
iwkq3dplotly <- plot_ly(x=vectau, y=vectau, z=iwkq3d, colorscale="Jet", type="surface", showscale=FALSE, scene="scene1", name=" ", hovertemplate = paste("<i>𝛕-investment<i>: %{x:.2f}", "<br>𝛕-capital: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% layout(scene1=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-investment"), yaxis=list(title="𝛕-capital"), zaxis=list(title="Investment Productivity"))) 
# iwkq3dplotly
# ##########################################################
# #Labor Estimates
##########################################################
#Elasticity (over labor)
l3dplotly <- plot_ly(x=vectau, y=vectau, z=l3d, colorscale="Jet", type="surface", showscale=FALSE, scene="scene2", name=" ", hovertemplate = paste("<i>𝛕-output<i>: %{x:.2f}", "<br>𝛕-labor: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% layout(scene2=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-output"), yaxis=list(title="𝛕-labor"), zaxis=list(title="Labor Elasticity")))  
# l3dplotly 
#Elasticity (over productivity)
lwq3dplotly <- plot_ly(x=vectau, y=vectau, z=lwq3d, colorscale="Jet", type="surface", showscale=FALSE, scene="scene2", name=" ", hovertemplate = paste("<i>𝛕-output<i>: %{x:.2f}", "<br>𝛕-productivity: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% layout(scene2=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-output"), yaxis=list(title="𝛕-productivity"), zaxis=list(title="Labor Elasticity")))  
# lwq3dplotly
#Elasticity (over capital)
lkq3dplotly <- plot_ly(x=vectau, y=vectau, z=lkq3d, colorscale="Jet", type="surface", showscale=FALSE, scene="scene1", name=" ", hovertemplate = paste("<i>𝛕-output<i>: %{x:.2f}", "<br>𝛕-capital: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% layout(scene1=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-output"), yaxis=list(title="𝛕-capital"), zaxis=list(title="Labor Elasticity")))  
# lkq3dplotly
#Elasticity (over materials)
lmq3dplotly <- plot_ly(x=vectau, y=vectau, z=lmq3d, colorscale="Jet", type="surface", showscale=FALSE, scene="scene2", name=" ", hovertemplate = paste("<i>𝛕-output<i>: %{x:.2f}", "<br>𝛕-materials: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% layout(scene2=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-output"), yaxis=list(title="𝛕-materials"), zaxis=list(title="Labor Elasticity")))  
# lmq3dplotly
#Efficiency (over labor)
hl3dplotly <- plot_ly(x=vectau, y=vectau, z=hl3d, colorscale="Jet", type="surface", showscale=FALSE, scene="scene2", name=" ", hovertemplate = paste("<i>𝛕-output<i>: %{x:.2f}", "<br>𝛕-labor: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% layout(scene2=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-output"), yaxis=list(title="𝛕-labor"), zaxis=list(title="Labor Efficiency"))) 
# hl3dplotly
#Efficiency (over capital)
hlkq3dplotly <- plot_ly(x=vectau, y=vectau, z=hlkq3d, colorscale="Jet", type="surface", showscale=FALSE, scene="scene1", name=" ", hovertemplate = paste("<i>𝛕-output<i>: %{x:.2f}", "<br>𝛕-capital: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% layout(scene1=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-output"), yaxis=list(title="𝛕-capital"), zaxis=list(title="Labor Efficiency"))) 
# hlkq3dplotly
#Efficiency (over materials)
hlmq3dplotly <- plot_ly(x=vectau, y=vectau, z=hmkq3d, colorscale="Jet", type="surface", showscale=FALSE, scene="scene2", name=" ", hovertemplate = paste("<i>𝛕-output<i>: %{x:.2f}", "<br>𝛕-materials: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% layout(scene2=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-output"), yaxis=list(title="𝛕-materials"), zaxis=list(title="Labor Efficiency"))) 
# hlmq3dplotly
#Labor response to productivity (over productivity)
lw3dplotly <- plot_ly(x=vectau, y=vectau, z=lw3d, colorscale="Jet", type="surface", showscale=FALSE, scene="scene2", name=" ", hovertemplate = paste("<i>𝛕-labor<i>: %{x:.2f}", "<br>𝛕-productivity: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% layout(scene2=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-labor"), yaxis=list(title="𝛕-productivity"), zaxis=list(title="Labor Productivity"))) 
# lw3dplotly
#Labor response to productivity (over capital)
lwkq3dplotly <- plot_ly(x=vectau, y=vectau, z=lwkq3d, colorscale="Jet", type="surface", showscale=FALSE, scene="scene2", name=" ", hovertemplate = paste("<i>𝛕-labor<i>: %{x:.2f}", "<br>𝛕-capital: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% layout(scene2=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-labor"), yaxis=list(title="𝛕-capital"), zaxis=list(title="Labor Productivity"))) 
# lwkq3dplotly
#########################################################
#Materials Estimates
#########################################################
#Elasticity (over materials)
m3dplotly <- plot_ly(x=vectau, y=vectau, z=m3d, colorscale="Jet", type="surface", showscale=FALSE, scene="scene3", name=" ", hovertemplate = paste("<i>𝛕-output<i>: %{x:.2f}", "<br>𝛕-materials: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% layout(scene3=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-output"), yaxis=list(title="𝛕-materials"), zaxis=list(title="Materials Elasticity"))) 
# m3dplotly 
#Elasticity (over productivity)
mwq3dplotly <- plot_ly(x=vectau, y=vectau, z=mwq3d, colorscale="Jet", type="surface", showscale=FALSE, scene="scene3", name=" ", hovertemplate = paste("<i>𝛕-output<i>: %{x:.2f}", "<br>𝛕-productivity: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% layout(scene3=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-output"), yaxis=list(title="𝛕-productivity"), zaxis=list(title="Materials Elasticity"))) 
# mwq3dplotly 
#Elasticity (over capital)
mkq3dplotly <- plot_ly(x=vectau, y=vectau, z=mkq3d, colorscale="Jet", type="surface", showscale=FALSE, scene="scene1", name=" ", hovertemplate = paste("<i>𝛕-output<i>: %{x:.2f}", "<br>𝛕-capital: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% layout(scene1=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-output"), yaxis=list(title="𝛕-capital"), zaxis=list(title="Materials Elasticity"))) 
# mkq3dplotly 
#Elasticity (over labor)
mlq3dplotly <- plot_ly(x=vectau, y=vectau, z=mlq3d, colorscale="Jet", type="surface", showscale=FALSE, scene="scene2", name=" ", hovertemplate = paste("<i>𝛕-output<i>: %{x:.2f}", "<br>𝛕-labor: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% layout(scene2=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-output"), yaxis=list(title="𝛕-labor"), zaxis=list(title="Materials Elasticity"))) 
# mlq3dplotly 
#Efficiency (over materials)
hm3dplotly <- plot_ly(x=vectau, y=vectau, z=hm3d, colorscale="Jet", type="surface", showscale=FALSE, scene="scene3", name=" ", hovertemplate = paste("<i>𝛕-output<i>: %{x:.2f}", "<br>𝛕-materials: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% layout(scene3=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-output"), yaxis=list(title="𝛕-materials"), zaxis=list(title="Materials Efficiency"))) 
# hm3dplotly 
#Efficiency (over capital)
hmkq3dplotly <- plot_ly(x=vectau, y=vectau, z=hmkq3d, colorscale="Jet", type="surface", showscale=FALSE, scene="scene1", name=" ", hovertemplate = paste("<i>𝛕-output<i>: %{x:.2f}", "<br>𝛕-capital: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% layout(scene1=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-output"), yaxis=list(title="𝛕-capital"), zaxis=list(title="Materials Efficiency"))) 
# hmkq3dplotly 
#Efficiency (over labor)
hmlq3dplotly <- plot_ly(x=vectau, y=vectau, z=hmlq3d, colorscale="Jet", type="surface", showscale=FALSE, scene="scene2", name=" ", hovertemplate = paste("<i>𝛕-output<i>: %{x:.2f}", "<br>𝛕-labor: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% layout(scene2=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-output"), yaxis=list(title="𝛕-labor"), zaxis=list(title="Materials Efficiency"))) 
# hmlq3dplotly 
#Material response to productivity (over productivity)
mw3dplotly <- plot_ly(x=vectau, y=vectau, z=mw3d, colorscale="Jet", type="surface", showscale=FALSE, scene="scene3", name=" ", hovertemplate = paste("<i>𝛕-materials<i>: %{x:.2f}", "<br>𝛕-productivity: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% layout(scene3=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-materials"), yaxis=list(title="𝛕-productivity"), zaxis=list(title="Materials Productivity"))) 
# mw3dplotly
#Material response to productivity (over capital)
mwkq3dplotly <- plot_ly(x=vectau, y=vectau, z=mwkq3d, colorscale="Jet", type="surface", showscale=FALSE, scene="scene3", name=" ", hovertemplate = paste("<i>𝛕-materials<i>: %{x:.2f}", "<br>𝛕-capital: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% layout(scene3=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-materials"), yaxis=list(title="𝛕-capital"), zaxis=list(title="Materials Productivity"))) 
# mw3dplotly
#Combined Plot.ly Elasticities
#Elasticities (over percentiles of inputs)
klm3dplotly <- subplot(k3dplotly, l3dplotly, m3dplotly, shareX=TRUE) %>% layout(scene=list(aspectratio=list(x=.6, y=.6, z=.6), xaxis=list(title="𝛕-output"), yaxis=list(title="𝛕-capital"), zaxis=list(title="Capital")), 
	scene2=list(aspectratio=list(x=.6, y=.6, z=.6), xaxis=list(title="𝛕-output"), yaxis=list(title="𝛕-labor"), zaxis=list(title="Labor")),
	scene3=list(aspectratio=list(x=.6, y=.6, z=.6), xaxis=list(title="𝛕-output"), yaxis=list(title="𝛕-materials"), zaxis=list(title="Materials")))
#Elasticities (over percentiles of productivity)
klmwq3dplotly <- subplot(kwq3dplotly, lwq3dplotly, mwq3dplotly, shareX=TRUE) %>% layout(scene=list(aspectratio=list(x=.6, y=.6, z=.6), xaxis=list(title="𝛕-output"), yaxis=list(title="𝛕-productivity"), zaxis=list(title="Capital")), 
	scene2=list(aspectratio=list(x=.6, y=.6, z=.6), xaxis=list(title="𝛕-output"), yaxis=list(title="𝛕-productivity"), zaxis=list(title="Labor")),
	scene3=list(aspectratio=list(x=.6, y=.6, z=.6), xaxis=list(title="𝛕-output"), yaxis=list(title="𝛕-productivity"), zaxis=list(title="Materials")))
#Elasticities (over percentiles of other inputs)
#Capital
kinpq3dplotly <- subplot(klq3dplotly, kmq3dplotly, shareX=TRUE) %>% layout(scene=list(aspectratio=list(x=.9, y=.9, z=.9), xaxis=list(title="𝛕-output"), yaxis=list(title="𝛕-labor"), zaxis=list(title="Capital")), 
	scene2=list(aspectratio=list(x=.9, y=.9, z=.9), xaxis=list(title="𝛕-output"), yaxis=list(title="𝛕-materials"), zaxis=list(title="Capital")))
#Labor
linpq3dplotly <- subplot(lkq3dplotly, lmq3dplotly, shareX=TRUE) %>% layout(scene=list(aspectratio=list(x=.9, y=.9, z=.9), xaxis=list(title="𝛕-output"), yaxis=list(title="𝛕-capital"), zaxis=list(title="Labor")), 
	scene2=list(aspectratio=list(x=.9, y=.9, z=.9), xaxis=list(title="𝛕-output"), yaxis=list(title="𝛕-materials"), zaxis=list(title="Labor")))
#Materials
minpq3dplotly <- subplot(mkq3dplotly, mlq3dplotly, shareX=TRUE) %>% layout(scene=list(aspectratio=list(x=.9, y=.9, z=.9), xaxis=list(title="𝛕-output"), yaxis=list(title="𝛕-capital"), zaxis=list(title="Materials")), 
	scene2=list(aspectratio=list(x=.9, y=.9, z=.9), xaxis=list(title="𝛕-output"), yaxis=list(title="𝛕-labor"), zaxis=list(title="Materials")))
#Annotations
annotationsklm <- list(list(x=0.13, y=0.75, text="(a) Capital Elasticity", font=list(size=30), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
	list(x=0.5, y=0.75, text="(b) Labor Elasticity", font=list(size=30), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
	list(x=0.9, y=0.75, text="(c) Materials Elasticity", font=list(size=30), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE))
annotationsab <- list(list(x=0.25, y=0.8, text="(a)", font=list(size=30), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
	list(x=0.75, y=0.8, text="(b)", font=list(size=30), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE))
#Add
klm3dplotly <- klm3dplotly %>% layout(annotations=annotationsklm)
klmwq3dplotly <- klmwq3dplotly %>% layout(annotations=annotationsklm)
kinpq3dplotly <- kinpq3dplotly %>% layout(annotations=annotationsab)
linpq3dplotly <- linpq3dplotly %>% layout(annotations=annotationsab)
minpq3dplotly <- minpq3dplotly %>% layout(annotations=annotationsab)
#Plot
klm3dplotly
klmwq3dplotly 
kinpq3dplotly
linpq3dplotly
minpq3dplotly
#Json
klm3dplotly   <- plotly_json(klm3dplotly, FALSE)
klmwq3dplotly   <- plotly_json(klmwq3dplotly, FALSE)
kinpq3dplotly   <- plotly_json(kinpq3dplotly, FALSE)
linpq3dplotly   <- plotly_json(linpq3dplotly, FALSE)
minpq3dplotly   <- plotly_json(minpq3dplotly, FALSE)
write(klm3dplotly, "/Users/justindoty/Documents/Home/My_Website/static/jmp/main/klm3dplotly.json")
write(klmwq3dplotly, "/Users/justindoty/Documents/Home/My_Website/static/jmp/main/klmwq3dplotly.json")
write(kinpq3dplotly, "/Users/justindoty/Documents/Home/My_Website/static/jmp/mainext/kinpq3dplotly.json")
write(linpq3dplotly, "/Users/justindoty/Documents/Home/My_Website/static/jmp/mainext/linpq3dplotly.json")
write(minpq3dplotly, "/Users/justindoty/Documents/Home/My_Website/static/jmp/mainext/minpq3dplotly.json")
#Combined Plot.ly Efficiecies
hklm3dplotly <- subplot(hk3dplotly, hl3dplotly, hm3dplotly, shareX=TRUE) %>% layout(scene=list(aspectratio=list(x=.6, y=.6, z=.6), xaxis=list(title="𝛕-output"), yaxis=list(title="𝛕-capital"), zaxis=list(title="Capital")), 
	scene2=list(aspectratio=list(x=.6, y=.6, z=.6), xaxis=list(title="𝛕-output"), yaxis=list(title="𝛕-labor"), zaxis=list(title="Labor")),
	scene3=list(aspectratio=list(x=.6, y=.6, z=.6), xaxis=list(title="𝛕-output"), yaxis=list(title="𝛕-materials"), zaxis=list(title="Materials")))
#Over other inputs
#Capital
hkinpq3dplotly <- subplot(hklq3dplotly, hkmq3dplotly, shareX=TRUE) %>% layout(scene=list(aspectratio=list(x=.9, y=.9, z=.9), xaxis=list(title="𝛕-output"), yaxis=list(title="𝛕-labor"), zaxis=list(title="Capital")), 
	scene2=list(aspectratio=list(x=.9, y=.9, z=.9), xaxis=list(title="𝛕-output"), yaxis=list(title="𝛕-materials"), zaxis=list(title="Capital")))
#Labor
hlinpq3dplotly <- subplot(hlkq3dplotly, hlmq3dplotly, shareX=TRUE) %>% layout(scene=list(aspectratio=list(x=.9, y=.9, z=.9), xaxis=list(title="𝛕-output"), yaxis=list(title="𝛕-capital"), zaxis=list(title="Labor")), 
	scene2=list(aspectratio=list(x=.9, y=.9, z=.9), xaxis=list(title="𝛕-output"), yaxis=list(title="𝛕-materials"), zaxis=list(title="Labor")))
#Materials
hminpq3dplotly <- subplot(hmkq3dplotly, hmlq3dplotly, shareX=TRUE) %>% layout(scene=list(aspectratio=list(x=.9, y=.9, z=.9), xaxis=list(title="𝛕-output"), yaxis=list(title="𝛕-capital"), zaxis=list(title="Materials")), 
	scene2=list(aspectratio=list(x=.9, y=.9, z=.9), xaxis=list(title="𝛕-output"), yaxis=list(title="𝛕-labor"), zaxis=list(title="Materials")))
#Annotations
annotationshklm <- list(list(x=0.13, y=0.75, text="(a) Capital Efficiency", font=list(size=30), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
	list(x=0.5, y=0.75, text="(b) Labor Efficiency", font=list(size=30), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
	list(x=0.9, y=0.75, text="(c) Materials Efficiency", font=list(size=30), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE))

#Add
hklm3dplotly <- hklm3dplotly %>% layout(annotations=annotationshklm)
hkinpq3dplotly <- hkinpq3dplotly %>% layout(annotations=annotationsab)
hlinpq3dplotly <- hlinpq3dplotly %>% layout(annotations=annotationsab)
hminpq3dplotly <- hminpq3dplotly %>% layout(annotations=annotationsab)
hklm3dplotly
hkinpq3dplotly
hlinpq3dplotly
hminpq3dplotly
hklm3dplotly   <- plotly_json(hklm3dplotly, FALSE)
hkinpq3dplotly   <- plotly_json(hkinpq3dplotly, FALSE)
hlinpq3dplotly   <- plotly_json(hlinpq3dplotly, FALSE)
hminpq3dplotly   <- plotly_json(hminpq3dplotly, FALSE)
write(hklm3dplotly, "/Users/justindoty/Documents/Home/My_Website/static/jmp/main/hklm3dplotly.json")
write(hkinpq3dplotly, "/Users/justindoty/Documents/Home/My_Website/static/jmp/mainext/hkinpq3dplotly.json")
write(hlinpq3dplotly, "/Users/justindoty/Documents/Home/My_Website/static/jmp/mainext/hlinpq3dplotly.json")
write(hminpq3dplotly, "/Users/justindoty/Documents/Home/My_Website/static/jmp/mainext/hminpq3dplotly.json")
#Combined Productivities
lmiw3dplotly <- subplot(iw3dplotly, lw3dplotly, mw3dplotly) %>% layout(scene=list(aspectratio=list(x=0.6, y=0.6, z=0.6), xaxis=list(title="𝛕-investment"), yaxis=list(title="𝛕-productivity"), zaxis=list(title="Investment")), 
	scene2=list(aspectratio=list(x=0.6, y=0.6, z=0.6), xaxis=list(title="𝛕-labor"), yaxis=list(title="𝛕-productivity"), zaxis=list(title="Labor")),
	scene3=list(aspectratio=list(x=0.6, y=0.6, z=0.6), xaxis=list(title="𝛕-materials"), yaxis=list(title="𝛕-productivity"), zaxis=list(title="Materials")))
lmiwkq3dplotly <- subplot(iwkq3dplotly, lwkq3dplotly, mwkq3dplotly) %>% layout(scene=list(aspectratio=list(x=0.6, y=0.6, z=0.6), xaxis=list(title="𝛕-investment"), yaxis=list(title="𝛕-capital"), zaxis=list(title="Investment")), 
	scene2=list(aspectratio=list(x=0.6, y=0.6, z=0.6), xaxis=list(title="𝛕-labor"), yaxis=list(title="𝛕-capital"), zaxis=list(title="Labor")),
	scene3=list(aspectratio=list(x=0.6, y=0.6, z=0.6), xaxis=list(title="𝛕-materials"), yaxis=list(title="𝛕-capital"), zaxis=list(title="Materials")))
#Annotations
annotationsw <- list(list(x=0.10, y=0.75, text="(a) Investment Productivity", font=list(size=30), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
	list(x=0.5, y=0.75, text="(b) Labor Productivity", font=list(size=30), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
	list(x=0.9, y=0.75, text="(c) Materials Productivity", font=list(size=30), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE))
#Add
lmiw3dplotly <- lmiw3dplotly %>% layout(annotations=annotationsw)
lmiwkq3dplotly <- lmiwkq3dplotly %>% layout(annotations=annotationsw)
lmiw3dplotly
lmiwkq3dplotly
lmiw3dplotly   <- plotly_json(lmiw3dplotly, FALSE)
lmiwkq3dplotly   <- plotly_json(lmiwkq3dplotly, FALSE)
write(lmiw3dplotly, "/Users/justindoty/Documents/Home/My_Website/static/jmp/main/lmiw3dplotly.json")
write(lmiwkq3dplotly, "/Users/justindoty/Documents/Home/My_Website/static/jmp/mainext/lmiwkq3dplotly.json")
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
omgplotly <- plot_ly(x=vectau, y=vectau, z=omg3dq, colorscale="Jet", type="surface", showscale=FALSE, scene="scene1", name=" ", hovertemplate = paste("<i>𝛕-innovation<i>: %{x:.2f}", "<br>𝛕-productivity: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% layout(scene1=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-innovation"), yaxis=list(title="𝛕-productivity"), zaxis=list(title="Persistence"))) 
# omgplotly
omgplotly <- plotly_json(omgplotly, FALSE)
write(omgplotly, "/Users/justindoty/Documents/Home/My_Website/static/jmp/main/omgplotly.json")
