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
stdY <- sd(US$Y)
stdK <- sd(US$K)
stdL <- sd(US$L)
stdM <- sd(US$M)
stdI <- sd(US$I)
#De-mean
US <- US %>% mutate(Y=(Y-mean(Y))/stdY, K=(K-mean(K))/stdK, L=(L-mean(L))/stdL, M=(M-mean(M))/stdM, I=(I-mean(I))/stdI)
wmin <- min(US$Y)
wmax <- max(US$Y)
###############################################################################
#Simulate Productivity and Capital Evolution Given Model Parameters
##############################################################################
#############################################################################
Nsim <- 1
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
	lnydata[,t] <- rowSums(PF(K=lnkdata[,t], L=lnldata[,t], M=lnmdata[,t], omega=omgdata[,t])*lspline(vectau=vectau, bvec=parY, b1=parYb[1], bL=parYb[2], u=etadata[,t]))
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
#Scale and Ratio Effects
ak <- parY[12,]; al <- parY[13,]; am <- parY[14,]
akk <- parY[18,]; all <- parY[19,]; amm <- parY[20,]
akl <- parY[15,]; akm <- parY[17,]; alm <- parY[16,]
hspar <- cbind((ak+al+am)/3, 2*(akk+all+amm)/3)
hklpar <- cbind((ak-al)/3, 2*(akk-(3/2)*akl+all)/3)
hkmpar <- cbind((ak-am)/3, 2*(akk-(3/2)*akm+amm)/3)
hlmpar <- cbind((al-am)/3, 2*(all-(3/2)*alm+amm)/3)
########################################################################################
#Productivity Dynamics####################################################
#########################################################################################
#Linear AR(1) process
lmomega <- as.numeric(lm(omega[idcon]~omega[idlag], data=US)$coefficients)[2]
lmomega <- matrix(lmomega, nrow=ntau, ncol=ntau)
lmomega_plot <- plot_ly(x=vectau, y=vectau, z=lmomega, colorscale="Jet", type="surface", showscale=FALSE, scene="scene1", name=" ", hovertemplate = paste("<i>𝛕-innovation<i>: %{x:.2f}", "<br>𝛕-productivity: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% layout(scene1=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-innovation", titlefont=list(size=20), tickfont=list(size=16)), yaxis=list(title="𝛕-productivity", titlefont=list(size=20), tickfont=list(size=16)), zaxis=list(title="Persistence", titlefont=list(size=20), tickfont=list(size=16)))) 
#Non-linear AR(1) process
nlomega <- as.numeric(lm(omega[idcon]~omega[idlag]+I(omega[idlag]^2)+I(omega[idlag]^3), data=US)$coefficients)[-1]
nlomegamat <- matrix(0, nrow=ntau, ncol=ntau)
for (q in 1:ntau){
	omgq <- rep(as.numeric(quantile(US$omega, probs=vectau[q])), nrow(US))
	nlomegamat[q,] <- colMeans(cbind(1, 2*omgq, 3*omgq^2)%*%nlomega)
}
nlomegamat_plot <- plot_ly(x=vectau, y=vectau, z=nlomegamat , colorscale="Jet", type="surface", showscale=FALSE, scene="scene2", name=" ", hovertemplate = paste("<i>𝛕-innovation<i>: %{x:.2f}", "<br>𝛕-productivity: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% layout(scene2=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-innovation", titlefont=list(size=20), tickfont=list(size=16)), yaxis=list(title="𝛕-productivity", titlefont=list(size=20), tickfont=list(size=16)), zaxis=list(title="Persistence", titlefont=list(size=20), tickfont=list(size=16))))
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
annotationsdist <- list(list(x=0.05, y=0.9, text="(a) Conditional Dispersion", font=list(size=35, family="Times New Roman"), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
	list(x=0.5, y=0.9, text="(b) Conditional Skewness", font=list(size=35, family="Times New Roman"), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
	list(x=0.95, y=0.9, text="(c) Conditional Kurtosis", font=list(size=35, family="Times New Roman"), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE))
omgdist <- subplot(displot, skewplot, kurplot, nrows=1, shareX=TRUE) %>% layout(scene=list(aspectratio=list(x=.6, y=.6, z=.6), xaxis=list(title="𝛕-productivity", titlefont=list(size=25), tickfont=list(size=25)), yaxis=list(title="Conditional Dispersion", titlefont=list(size=18), tickfont=list(size=20))), 
	scene2=list(aspectratio=list(x=.6, y=.6, z=.6), xaxis=list(title="𝛕-productivity", titlefont=list(size=25), tickfont=list(size=25)), yaxis=list(title="Conditional Skewness", titlefont=list(size=18), tickfont=list(size=25))),
	scene3=list(aspectratio=list(x=.6, y=.6, z=.6), xaxis=list(title="𝛕-productivity", titlefont=list(size=25), tickfont=list(size=25)), yaxis=list(title="Conditional Kurtosis", titlefont=list(size=18), tickfont=list(size=25))))
omgdist <- omgdist %>% layout(annotations=annotationsdist)
omgdist
########################################################################################
#Individual Quantile Marginal Effects########################x############################
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
#For scale and ratio effects
hs <- k3d; hkl <- k3d; hkm <- k3d; hlm <- k3d
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
	#Scale and Ratio Effects
	hs[q,] <- c(1, quantile(cap+lab+mat, probs=vectau[q]))%*%t(hspar)
	hkl[q,] <- c(1, quantile(cap-lab, probs=vectau[q]))%*%t(hklpar)
	hkm[q,] <- c(1, quantile(cap-mat, probs=vectau[q]))%*%t(hkmpar)
	hlm[q,] <- c(1, quantile(lab-mat, probs=vectau[q]))%*%t(hlmpar)
	#Marginal Productivities (fixed percentiles of productivity)
	iw3d[q,] <- colMeans(IXD(K=kmean, omega=wq, par=parI, pos=2, sdpos=1))
	lw3d[q,] <- colMeans(LXD(K=kmean, omega=wq, par=parL, pos=2, sdpos=1))
	mw3d[q,] <- colMeans(MXD(K=kmean, L=labw, omega=wq, par=parM, pos=3, sdpos=1))
	#Marginal Productivities (fixed percentiles of capital)
	iwkq3d[q,] <- colMeans(IXD(K=capkq, omega=wmean, par=parI, pos=2, sdpos=1))
	lwkq3d[q,] <- colMeans(LXD(K=capkq, omega=wmean, par=parL, pos=2, sdpos=1))
	mwkq3d[q,] <- colMeans(MXD(K=capkq, L=labk, omega=wmean, par=parM, pos=3, sdpos=1))

}
###################################################
#Capital Estimates
##################################################
#Elasticity (over capital)
kplot <- plot_ly(x=vectau, y=vectau, z=k3d, colorscale="Jet", type="surface", showscale=FALSE, scene="scene1", name=" ", hovertemplate = paste("<i>𝛕-output<i>: %{x:.2f}", "<br>𝛕-capital: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% layout(scene1=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-output"), yaxis=list(title="𝛕-capital"), zaxis=list(title="Capital Elasticity")))
# kplot
# k3dplotly
# Elasticity (over productivity)
kwqplot <- plot_ly(x=vectau, y=vectau, z=kwq3d, colorscale="Jet", type="surface", showscale=FALSE, scene="scene1", name=" ", hovertemplate = paste("<i>𝛕-output<i>: %{x:.2f}", "<br>𝛕-productivity: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% layout(scene1=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-output"), yaxis=list(title="𝛕-productivity"), zaxis=list(title="Capital Elasticity")))
# kwq3dplotly
# Elasticity (over labor)
klqplot <- plot_ly(x=vectau, y=vectau, z=klq3d, colorscale="Jet", type="surface", showscale=FALSE, scene="scene1", name=" ", hovertemplate = paste("<i>𝛕-output<i>: %{x:.2f}", "<br>𝛕-labor: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% layout(scene1=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-output"), yaxis=list(title="𝛕-labor"), zaxis=list(title="Capital Elasticity")))
# klq3dplotly
# Elasticity (over materials)
kmqplot <- plot_ly(x=vectau, y=vectau, z=kmq3d, colorscale="Jet", type="surface", showscale=FALSE, scene="scene2", name=" ", hovertemplate = paste("<i>𝛕-output<i>: %{x:.2f}", "<br>𝛕-materials: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% layout(scene2=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-output"), yaxis=list(title="𝛕-materials"), zaxis=list(title="Capital Elasticity")))
# kmq3dplotly
#Efficiency (over capital)
hkplot <- plot_ly(x=vectau, y=vectau, z=hk3d, colorscale="Jet", type="surface", showscale=FALSE, scene="scene1", name=" ", hovertemplate = paste("<i>𝛕-output<i>: %{x:.2f}", "<br>𝛕-capital: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% layout(scene1=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-output"), yaxis=list(title="𝛕-capital"), zaxis=list(title="Capital Efficiency"))) 
# hk3dplotly
#Efficiency (over labor)
hklqplot <- plot_ly(x=vectau, y=vectau, z=hklq3d, colorscale="Jet", type="surface", showscale=FALSE, scene="scene1", name=" ", hovertemplate = paste("<i>𝛕-output<i>: %{x:.2f}", "<br>𝛕-labor: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% layout(scene1=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-output"), yaxis=list(title="𝛕-labor"), zaxis=list(title="Capital Efficiency"))) 
# hklq3dplotly
#Efficiency (over materials)
hkmqplot <- plot_ly(x=vectau, y=vectau, z=hkmq3d, colorscale="Jet", type="surface", showscale=FALSE, scene="scene2", name=" ", hovertemplate = paste("<i>𝛕-output<i>: %{x:.2f}", "<br>𝛕-materials: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% layout(scene2=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-output"), yaxis=list(title="𝛕-materials"), zaxis=list(title="Capital Efficiency"))) 
# hklq3dplotly
#Investment Response to productivity (over productivity)
iwplot <- plot_ly(x=vectau, y=vectau, z=iw3d, colorscale="Jet", type="surface", showscale=FALSE, scene="scene1", name=" ", hovertemplate = paste("<i>𝛕-investment<i>: %{x:.2f}", "<br>𝛕-productivity: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% layout(scene1=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-investment"), yaxis=list(title="𝛕-productivity"), zaxis=list(title="Investment Productivity"))) 
# iw3dplotly
#Investment Response to productivity (over capital)
iwkqplot <- plot_ly(x=vectau, y=vectau, z=iwkq3d, colorscale="Jet", type="surface", showscale=FALSE, scene="scene1", name=" ", hovertemplate = paste("<i>𝛕-investment<i>: %{x:.2f}", "<br>𝛕-capital: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% layout(scene1=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-investment"), yaxis=list(title="𝛕-capital"), zaxis=list(title="Investment Productivity"))) 
# iwkq3dplotly
# ##########################################################
# #Labor Estimates
##########################################################
#Elasticity (over labor)
lplot <- plot_ly(x=vectau, y=vectau, z=l3d, colorscale="Jet", type="surface", showscale=FALSE, scene="scene2", name=" ", hovertemplate = paste("<i>𝛕-output<i>: %{x:.2f}", "<br>𝛕-labor: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% layout(scene2=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-output"), yaxis=list(title="𝛕-labor"), zaxis=list(title="Labor Elasticity")))  
# l3dplotly 
#Elasticity (over productivity)
lwqplot <- plot_ly(x=vectau, y=vectau, z=lwq3d, colorscale="Jet", type="surface", showscale=FALSE, scene="scene2", name=" ", hovertemplate = paste("<i>𝛕-output<i>: %{x:.2f}", "<br>𝛕-productivity: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% layout(scene2=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-output"), yaxis=list(title="𝛕-productivity"), zaxis=list(title="Labor Elasticity")))  
# lwq3dplotly
#Elasticity (over capital)
lkqplot <- plot_ly(x=vectau, y=vectau, z=lkq3d, colorscale="Jet", type="surface", showscale=FALSE, scene="scene1", name=" ", hovertemplate = paste("<i>𝛕-output<i>: %{x:.2f}", "<br>𝛕-capital: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% layout(scene1=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-output"), yaxis=list(title="𝛕-capital"), zaxis=list(title="Labor Elasticity")))  
# lkq3dplotly
#Elasticity (over materials)
lmqplot <- plot_ly(x=vectau, y=vectau, z=lmq3d, colorscale="Jet", type="surface", showscale=FALSE, scene="scene2", name=" ", hovertemplate = paste("<i>𝛕-output<i>: %{x:.2f}", "<br>𝛕-materials: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% layout(scene2=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-output"), yaxis=list(title="𝛕-materials"), zaxis=list(title="Labor Elasticity")))  
# lmq3dplotly
#Efficiency (over labor)
hlplot <- plot_ly(x=vectau, y=vectau, z=hl3d, colorscale="Jet", type="surface", showscale=FALSE, scene="scene2", name=" ", hovertemplate = paste("<i>𝛕-output<i>: %{x:.2f}", "<br>𝛕-labor: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% layout(scene2=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-output"), yaxis=list(title="𝛕-labor"), zaxis=list(title="Labor Efficiency"))) 
# hl3dplotly
#Efficiency (over capital)
hlkqplot <- plot_ly(x=vectau, y=vectau, z=hlkq3d, colorscale="Jet", type="surface", showscale=FALSE, scene="scene1", name=" ", hovertemplate = paste("<i>𝛕-output<i>: %{x:.2f}", "<br>𝛕-capital: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% layout(scene1=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-output"), yaxis=list(title="𝛕-capital"), zaxis=list(title="Labor Efficiency"))) 
# hlkq3dplotly
#Efficiency (over materials)
hlmqplot <- plot_ly(x=vectau, y=vectau, z=hlmq3d, colorscale="Jet", type="surface", showscale=FALSE, scene="scene2", name=" ", hovertemplate = paste("<i>𝛕-output<i>: %{x:.2f}", "<br>𝛕-materials: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% layout(scene2=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-output"), yaxis=list(title="𝛕-materials"), zaxis=list(title="Labor Efficiency"))) 
# hlmq3dplotly
#Labor response to productivity (over productivity)
lwplot <- plot_ly(x=vectau, y=vectau, z=lw3d, colorscale="Jet", type="surface", showscale=FALSE, scene="scene2", name=" ", hovertemplate = paste("<i>𝛕-labor<i>: %{x:.2f}", "<br>𝛕-productivity: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% layout(scene2=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-labor"), yaxis=list(title="𝛕-productivity"), zaxis=list(title="Labor Productivity"))) 
# lw3dplotly
#Labor response to productivity (over capital)
lwkqplot <- plot_ly(x=vectau, y=vectau, z=lwkq3d, colorscale="Jet", type="surface", showscale=FALSE, scene="scene2", name=" ", hovertemplate = paste("<i>𝛕-labor<i>: %{x:.2f}", "<br>𝛕-capital: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% layout(scene2=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-labor"), yaxis=list(title="𝛕-capital"), zaxis=list(title="Labor Productivity"))) 
# lwkq3dplotly
#########################################################
#Materials Estimates
#########################################################
#Elasticity (over materials)
mplot <- plot_ly(x=vectau, y=vectau, z=m3d, colorscale="Jet", type="surface", showscale=FALSE, scene="scene3", name=" ", hovertemplate = paste("<i>𝛕-output<i>: %{x:.2f}", "<br>𝛕-materials: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% layout(scene3=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-output"), yaxis=list(title="𝛕-materials"), zaxis=list(title="Materials Elasticity"))) 
# m3dplotly 
#Elasticity (over productivity)
mwqplot <- plot_ly(x=vectau, y=vectau, z=mwq3d, colorscale="Jet", type="surface", showscale=FALSE, scene="scene3", name=" ", hovertemplate = paste("<i>𝛕-output<i>: %{x:.2f}", "<br>𝛕-productivity: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% layout(scene3=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-output"), yaxis=list(title="𝛕-productivity"), zaxis=list(title="Materials Elasticity"))) 
# mwq3dplotly 
#Elasticity (over capital)
mkqplot <- plot_ly(x=vectau, y=vectau, z=mkq3d, colorscale="Jet", type="surface", showscale=FALSE, scene="scene1", name=" ", hovertemplate = paste("<i>𝛕-output<i>: %{x:.2f}", "<br>𝛕-capital: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% layout(scene1=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-output"), yaxis=list(title="𝛕-capital"), zaxis=list(title="Materials Elasticity"))) 
# mkq3dplotly 
#Elasticity (over labor)
mlqplot <- plot_ly(x=vectau, y=vectau, z=mlq3d, colorscale="Jet", type="surface", showscale=FALSE, scene="scene2", name=" ", hovertemplate = paste("<i>𝛕-output<i>: %{x:.2f}", "<br>𝛕-labor: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% layout(scene2=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-output"), yaxis=list(title="𝛕-labor"), zaxis=list(title="Materials Elasticity"))) 
# mlq3dplotly 
#Efficiency (over materials)
hmplot <- plot_ly(x=vectau, y=vectau, z=hm3d, colorscale="Jet", type="surface", showscale=FALSE, scene="scene3", name=" ", hovertemplate = paste("<i>𝛕-output<i>: %{x:.2f}", "<br>𝛕-materials: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% layout(scene3=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-output"), yaxis=list(title="𝛕-materials"), zaxis=list(title="Materials Efficiency"))) 
# hm3dplotly 
#Efficiency (over capital)
hmkqplot <- plot_ly(x=vectau, y=vectau, z=hmkq3d, colorscale="Jet", type="surface", showscale=FALSE, scene="scene1", name=" ", hovertemplate = paste("<i>𝛕-output<i>: %{x:.2f}", "<br>𝛕-capital: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% layout(scene1=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-output"), yaxis=list(title="𝛕-capital"), zaxis=list(title="Materials Efficiency"))) 
# hmkq3dplotly 
#Efficiency (over labor)
hmlqplot <- plot_ly(x=vectau, y=vectau, z=hmlq3d, colorscale="Jet", type="surface", showscale=FALSE, scene="scene2", name=" ", hovertemplate = paste("<i>𝛕-output<i>: %{x:.2f}", "<br>𝛕-labor: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% layout(scene2=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-output"), yaxis=list(title="𝛕-labor"), zaxis=list(title="Materials Efficiency"))) 
# hmlq3dplotly 
#Material response to productivity (over productivity)
mwplot <- plot_ly(x=vectau, y=vectau, z=mw3d, colorscale="Jet", type="surface", showscale=FALSE, scene="scene3", name=" ", hovertemplate = paste("<i>𝛕-materials<i>: %{x:.2f}", "<br>𝛕-productivity: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% layout(scene3=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-materials"), yaxis=list(title="𝛕-productivity"), zaxis=list(title="Materials Productivity"))) 
# mw3dplotly
#Material response to productivity (over capital)
mwkqplot <- plot_ly(x=vectau, y=vectau, z=mwkq3d, colorscale="Jet", type="surface", showscale=FALSE, scene="scene3", name=" ", hovertemplate = paste("<i>𝛕-materials<i>: %{x:.2f}", "<br>𝛕-capital: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% layout(scene3=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-materials"), yaxis=list(title="𝛕-capital"), zaxis=list(title="Materials Productivity"))) 
# mw3dplotly
###############################################
#Scale and Ratio Effects
#Scale
hsplot <- plot_ly(x=vectau, y=vectau, z=hs, colorscale="Jet", type="surface", showscale=FALSE, scene="scene1", name=" ", hovertemplate = paste("<i>𝛕-output<i>: %{x:.2f}", "<br>𝛕-scale: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% layout(scene1=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-output", titlefont=list(size=30), tickfont=list(size=16)), yaxis=list(title="𝛕-scale", titlefont=list(size=30), tickfont=list(size=16)), zaxis=list(title="Scale Effect", titlefont=list(size=30), tickfont=list(size=16)))) 
hsplot
#Capital-Labor
hklplot <- plot_ly(x=vectau, y=vectau, z=hkl, colorscale="Jet", type="surface", showscale=FALSE, scene="scene1", name=" ", hovertemplate = paste("<i>𝛕-output<i>: %{x:.2f}", "<br>𝛕-k/l: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% layout(scene1=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-output"), yaxis=list(title="𝛕-k/l"), zaxis=list(title="Capital/Labor Effect"))) 
# hklplot
#Capital-Materials
hkmplot <- plot_ly(x=vectau, y=vectau, z=hkm, colorscale="Jet", type="surface", showscale=FALSE, scene="scene2", name=" ", hovertemplate = paste("<i>𝛕-output<i>: %{x:.2f}", "<br>𝛕-k/m: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% layout(scene2=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-output"), yaxis=list(title="𝛕-k/m"), zaxis=list(title="Capital/Materials Effect"))) 
# hkmplot
#Labor-Materials
hlmplot <- plot_ly(x=vectau, y=vectau, z=hlm, colorscale="Jet", type="surface", showscale=FALSE, scene="scene3", name=" ", hovertemplate = paste("<i>𝛕-output<i>: %{x:.2f}", "<br>𝛕-l/m: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% layout(scene3=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-output"), yaxis=list(title="𝛕-l/m"), zaxis=list(title="Labor/Materials Effect"))) 
# hlmplot
#Combined Plot.ly Elasticities
#Elasticities (over percentiles of inputs)
klmplot <- subplot(kplot, lplot, mplot, shareX=TRUE) %>% layout(scene=list(aspectratio=list(x=.6, y=.6, z=.6), xaxis=list(title="𝛕-output", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="𝛕-capital", titlefont=list(size=18), tickfont=list(size=14)), zaxis=list(title="Capital", titlefont=list(size=18), tickfont=list(size=14))), 
	scene2=list(aspectratio=list(x=.6, y=.6, z=.6), xaxis=list(title="𝛕-output", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="𝛕-labor", titlefont=list(size=18), tickfont=list(size=14)), zaxis=list(title="Labor", titlefont=list(size=18), tickfont=list(size=14))),
	scene3=list(aspectratio=list(x=.6, y=.6, z=.6), xaxis=list(title="𝛕-output", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="𝛕-materials", titlefont=list(size=18), tickfont=list(size=14)), zaxis=list(title="Materials", titlefont=list(size=18), tickfont=list(size=14))))
# Elasticities (over percentiles of productivity)
klmwqplot <- subplot(kwqplot, lwqplot, mwqplot, shareX=TRUE) %>% layout(scene=list(aspectratio=list(x=.6, y=.6, z=.6), xaxis=list(title="𝛕-output", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="𝛕-productivity", titlefont=list(size=18), tickfont=list(size=14)), zaxis=list(title="Capital", titlefont=list(size=18), tickfont=list(size=14))), 
	scene2=list(aspectratio=list(x=.6, y=.6, z=.6), xaxis=list(title="𝛕-output", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="𝛕-productivity", titlefont=list(size=18), tickfont=list(size=14)), zaxis=list(title="Labor", titlefont=list(size=18), tickfont=list(size=14))),
	scene3=list(aspectratio=list(x=.6, y=.6, z=.6), xaxis=list(title="𝛕-output", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="𝛕-productivity", titlefont=list(size=18), tickfont=list(size=14)), zaxis=list(title="Materials", titlefont=list(size=18), tickfont=list(size=14))))
#Elasticities (over percentiles of other inputs)
#Capital
kinpqplot <- subplot(klqplot, kmqplot, shareX=TRUE) %>% layout(scene=list(aspectratio=list(x=.9, y=.9, z=.9), xaxis=list(title="𝛕-output", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="𝛕-labor", titlefont=list(size=18), tickfont=list(size=14)), zaxis=list(title="Capital", titlefont=list(size=18), tickfont=list(size=14))), 
	scene2=list(aspectratio=list(x=.9, y=.9, z=.9), xaxis=list(title="𝛕-output", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="𝛕-materials", titlefont=list(size=18), tickfont=list(size=14)), zaxis=list(title="Capital", titlefont=list(size=18), tickfont=list(size=14))))
#Labor
linpqplot <- subplot(lkqplot, lmqplot, shareX=TRUE) %>% layout(scene=list(aspectratio=list(x=.9, y=.9, z=.9), xaxis=list(title="𝛕-output", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="𝛕-capital", titlefont=list(size=18), tickfont=list(size=14)), zaxis=list(title="Labor", titlefont=list(size=18), tickfont=list(size=14))), 
	scene2=list(aspectratio=list(x=.9, y=.9, z=.9), xaxis=list(title="𝛕-output", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="𝛕-materials", titlefont=list(size=18), tickfont=list(size=14)), zaxis=list(title="Labor", titlefont=list(size=18), tickfont=list(size=14))))
#Materials
minpqplot <- subplot(mkqplot, mlqplot, shareX=TRUE) %>% layout(scene=list(aspectratio=list(x=.9, y=.9, z=.9), xaxis=list(title="𝛕-output", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="𝛕-capital", titlefont=list(size=18), tickfont=list(size=14)), zaxis=list(title="Materials", titlefont=list(size=18), tickfont=list(size=14))), 
	scene2=list(aspectratio=list(x=.9, y=.9, z=.9), xaxis=list(title="𝛕-output", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="𝛕-labor", titlefont=list(size=18), tickfont=list(size=14)), zaxis=list(title="Materials", titlefont=list(size=18), tickfont=list(size=14))))
#Annotations
#I use two different text sizes (one for latex and one for my website)
latexannotationsklm <- list(list(x=0.09, y=0.75, text="(a) Capital Elasticity", font=list(size=30, family="Times New Roman"), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
	list(x=0.5, y=0.75, text="(b) Labor Elasticity", font=list(size=30, family="Times New Roman"), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
	list(x=0.95, y=0.75, text="(c) Materials Elasticity", font=list(size=30, family="Times New Roman"), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE))
latexannotationsab <- list(list(x=0.25, y=0.8, text="<b>(a)<b>", font=list(size=30, family="Times New Roman"), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
	list(x=0.75, y=0.8, text="<b>(b)<b>", font=list(size=30, family="Times New Roman"), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE))
annotationsklm <- list(list(x=0.1, y=0.75, text="(a) Capital Elasticity", font=list(size=18), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
	list(x=0.5, y=0.75, text="(b) Labor Elasticity", font=list(size=18), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
	list(x=0.9, y=0.75, text="(c) Materials Elasticity", font=list(size=18), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE))
annotationsab <- list(list(x=0.25, y=0.8, text="(a)", font=list(size=18), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
	list(x=0.75, y=0.8, text="(b)", font=list(size=18), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE))
#Add
klmplot_latex <- klmplot %>% layout(annotations=latexannotationsklm)
klmwqplot_latex <- klmwqplot %>% layout(annotations=latexannotationsklm)
kinpqplot_latex <- kinpqplot %>% layout(annotations=latexannotationsab)
linpqplot_latex <- linpqplot %>% layout(annotations=latexannotationsab)
minpqplot_latex <- minpqplot %>% layout(annotations=latexannotationsab)
#Plot
klmplot_latex 
klmwqplot_latex 
kinpqplot_latex 
linpqplot_latex  
minpqplot_latex 
# Save to JSON
klmplot <- klmplot %>% layout(annotations=annotationsklm)
klmwqplot <- klmwqplot %>% layout(annotations=annotationsklm)
kinpqplot <- kinpqplot %>% layout(annotations=annotationsab)
linpqplot <- linpqplot %>% layout(annotations=annotationsab)
minpqplot <- minpqplot %>% layout(annotations=annotationsab)
#Json
klmplot   <- plotly_json(klmplot, FALSE)
klmwqplot  <- plotly_json(klmwqplot, FALSE)
kinpqplot  <- plotly_json(kinpqplot, FALSE)
linpqplot   <- plotly_json(linpqplot, FALSE)
minpqplot   <- plotly_json(minpqplot, FALSE)
# write(klmplot, "/Users/justindoty/Documents/Home/My_Website/static/jmp/main/klm3dplotly.json")
# write(klmwqplot , "/Users/justindoty/Documents/Home/My_Website/static/jmp/main/klmwq3dplotly.json")
# write(kinpqplot, "/Users/justindoty/Documents/Home/My_Website/static/jmp/mainext/kinpq3dplotly.json")
# write(linpqplot, "/Users/justindoty/Documents/Home/My_Website/static/jmp/mainext/linpq3dplotly.json")
# write(minpqplot, "/Users/justindoty/Documents/Home/My_Website/static/jmp/mainext/minpq3dplotly.json")
#Combined Plot.ly Efficiecies
hklmplot <- subplot(hkplot, hlplot, hmplot, shareX=TRUE) %>% layout(scene=list(aspectratio=list(x=.6, y=.6, z=.6), xaxis=list(title="𝛕-output", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="𝛕-capital", titlefont=list(size=18), tickfont=list(size=14)), zaxis=list(title="Capital", titlefont=list(size=18), tickfont=list(size=14))), 
	scene2=list(aspectratio=list(x=.6, y=.6, z=.6), xaxis=list(title="𝛕-output", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="𝛕-labor", titlefont=list(size=18), tickfont=list(size=14)), zaxis=list(title="Labor", titlefont=list(size=18), tickfont=list(size=14))),
	scene3=list(aspectratio=list(x=.6, y=.6, z=.6), xaxis=list(title="𝛕-output", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="𝛕-materials", titlefont=list(size=18), tickfont=list(size=14)), zaxis=list(title="Materials", titlefont=list(size=18), tickfont=list(size=14))))
#Over other inputs
#Capital
hkinpqplot <- subplot(hklqplot, hkmqplot, shareX=TRUE) %>% layout(scene=list(aspectratio=list(x=.9, y=.9, z=.9), xaxis=list(title="𝛕-output", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="𝛕-labor", titlefont=list(size=18), tickfont=list(size=14)), zaxis=list(title="Capital", titlefont=list(size=18), tickfont=list(size=14))), 
	scene2=list(aspectratio=list(x=.9, y=.9, z=.9), xaxis=list(title="𝛕-output", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="𝛕-materials", titlefont=list(size=18), tickfont=list(size=14)), zaxis=list(title="Capital", titlefont=list(size=18), tickfont=list(size=14))))
#Labor
hlinpqplot <- subplot(hlkqplot, hlmqplot, shareX=TRUE) %>% layout(scene=list(aspectratio=list(x=.9, y=.9, z=.9), xaxis=list(title="𝛕-output", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="𝛕-capital", titlefont=list(size=18), tickfont=list(size=14)), zaxis=list(title="Labor", titlefont=list(size=18), tickfont=list(size=14))), 
	scene2=list(aspectratio=list(x=.9, y=.9, z=.9), xaxis=list(title="𝛕-output", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="𝛕-materials", titlefont=list(size=18), tickfont=list(size=14)), zaxis=list(title="Labor", titlefont=list(size=18), tickfont=list(size=14))))
#Materials
hminpqplot <- subplot(hmkqplot, hmlqplot, shareX=TRUE) %>% layout(scene=list(aspectratio=list(x=.9, y=.9, z=.9), xaxis=list(title="𝛕-output", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="𝛕-capital", titlefont=list(size=18), tickfont=list(size=14)), zaxis=list(title="Materials", titlefont=list(size=18), tickfont=list(size=14))), 
	scene2=list(aspectratio=list(x=.9, y=.9, z=.9), xaxis=list(title="𝛕-output", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="𝛕-labor", titlefont=list(size=18), tickfont=list(size=14)), zaxis=list(title="Materials", titlefont=list(size=18), tickfont=list(size=14))))
#Annotations
latexannotationshklm <- list(list(x=0.09, y=0.75, text="(a) Capital Efficiency", font=list(size=30, family="Times New Roman"), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
	list(x=0.5, y=0.75, text="(b) Labor Efficiency", font=list(size=30, family="Times New Roman"), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
	list(x=0.95, y=0.75, text="(c) Materials Efficiency", font=list(size=30, family="Times New Roman"), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE))
annotationshklm <- list(list(x=0.1, y=0.75, text="(a) Capital Efficiency", font=list(size=18), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
	list(x=0.5, y=0.75, text="(b) Labor Efficiency", font=list(size=18), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
	list(x=0.9, y=0.75, text="(c) Materials Efficiency", font=list(size=18), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE))

#Add
hklmplot_latex <- hklmplot %>% layout(annotations=latexannotationshklm)
hkinpqplot_latex <- hkinpqplot %>% layout(annotations=latexannotationsab)
hlinpqplot_latex <- hlinpqplot %>% layout(annotations=latexannotationsab)
hminpqplot_latex <- hminpqplot %>% layout(annotations=latexannotationsab)
#Plot
hklmplot_latex 
hkinpqplot_latex 
hlinpqplot_latex 
hminpqplot_latex 
#JSON
hklmplot <- hklmplot %>% layout(annotations=annotationshklm)
hkinpqplot <- hkinpqplot %>% layout(annotations=annotationsab)
hlinpqplot <- hlinpqplot %>% layout(annotations=annotationsab)
hminpqplot <- hminpqplot %>% layout(annotations=annotationsab)
#Save
hklmplot   <- plotly_json(hklmplot, FALSE)
hkinpqplot   <- plotly_json(hkinpqplot, FALSE)
hlinpqplot   <- plotly_json(hlinpqplot, FALSE)
hminpqplot   <- plotly_json(hminpqplot, FALSE)
# write(hklmplot, "/Users/justindoty/Documents/Home/My_Website/static/jmp/main/hklm3dplotly.json")
# write(hkinpqplot, "/Users/justindoty/Documents/Home/My_Website/static/jmp/mainext/hkinpq3dplotly.json")
# write(hlinpqplot, "/Users/justindoty/Documents/Home/My_Website/static/jmp/mainext/hlinpq3dplotly.json")
# write(hminpqplot, "/Users/justindoty/Documents/Home/My_Website/static/jmp/mainext/hminpq3dplotly.json")
#Scale and Ratio Effects
#Annotations
ratioannotations <- list(list(x=0.03, y=0.75, text="(a) Capital-Labor Efficiency", font=list(size=30, family="Times New Roman"), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
	list(x=0.5, y=0.75, text="(b) Capital-Materials Efficiency", font=list(size=30, family="Times New Roman"), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
	list(x=0.95, y=0.75, text="(c) Labor-Materials Efficiency", font=list(size=30, family="Times New Roman"), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE))
ratioplots <- subplot(hklplot, hkmplot, hlmplot, shareX=TRUE) %>% layout(scene=list(aspectratio=list(x=.6, y=.6, z=.6), xaxis=list(title="𝛕-output", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="𝛕-(k-l)", titlefont=list(size=18), tickfont=list(size=14)), zaxis=list(title="Capital-Labor", titlefont=list(size=18), tickfont=list(size=14))), 
	scene2=list(aspectratio=list(x=.6, y=.6, z=.6), xaxis=list(title="𝛕-output", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="𝛕-(k-m)", titlefont=list(size=18), tickfont=list(size=14)), zaxis=list(title="Capital-Materials", titlefont=list(size=18), tickfont=list(size=14))),
	scene3=list(aspectratio=list(x=.6, y=.6, z=.6), xaxis=list(title="𝛕-output", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="𝛕-(l-m)", titlefont=list(size=18), tickfont=list(size=14)), zaxis=list(title="Labor-Materials", titlefont=list(size=18), tickfont=list(size=14))))
ratioplots <- ratioplots %>% layout(annotations=ratioannotations)
ratioplots






#Combined Productivities
lmiwplot <- subplot(iwplot, lwplot, mwplot) %>% layout(scene=list(aspectratio=list(x=0.6, y=0.6, z=0.6), xaxis=list(title="𝛕-investment", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="𝛕-productivity", titlefont=list(size=18), tickfont=list(size=14)), zaxis=list(title="Investment", titlefont=list(size=18), tickfont=list(size=14))), 
	scene2=list(aspectratio=list(x=0.6, y=0.6, z=0.6), xaxis=list(title="𝛕-labor", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="𝛕-productivity", titlefont=list(size=18), tickfont=list(size=14)), zaxis=list(title="Labor", titlefont=list(size=18), tickfont=list(size=14))),
	scene3=list(aspectratio=list(x=0.6, y=0.6, z=0.6), xaxis=list(title="𝛕-materials", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="𝛕-productivity", titlefont=list(size=18), tickfont=list(size=14)), zaxis=list(title="Materials", titlefont=list(size=18), tickfont=list(size=14))))
lmiwkqplot <- subplot(iwkqplot, lwkqplot, mwkqplot) %>% layout(scene=list(aspectratio=list(x=0.6, y=0.6, z=0.6), xaxis=list(title="𝛕-investment", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="𝛕-capital", titlefont=list(size=18), tickfont=list(size=14)), zaxis=list(title="Investment", titlefont=list(size=18), tickfont=list(size=14))), 
	scene2=list(aspectratio=list(x=0.6, y=0.6, z=0.6), xaxis=list(title="𝛕-labor", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="𝛕-capital", titlefont=list(size=18), tickfont=list(size=14)), zaxis=list(title="Labor", titlefont=list(size=18), tickfont=list(size=14))),
	scene3=list(aspectratio=list(x=0.6, y=0.6, z=0.6), xaxis=list(title="𝛕-materials", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="𝛕-capital", titlefont=list(size=18), tickfont=list(size=14)), zaxis=list(title="Materials", titlefont=list(size=18), tickfont=list(size=14))))
#Annotations
latexannotationsw <- list(list(x=0.08, y=0.75, text="(a) Investment Response", font=list(size=30, family="Times New Roman"), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
	list(x=0.5, y=0.75, text="(b) Labor Response", font=list(size=30, family="Times New Roman"), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
	list(x=0.92, y=0.75, text="(c) Materials Response", font=list(size=30, family="Times New Roman"), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE))
annotationsw <- list(list(x=0.07, y=0.75, text="(a) Investment Response", font=list(size=18), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
	list(x=0.5, y=0.75, text="(b) Labor Response", font=list(size=18), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
	list(x=0.92, y=0.75, text="(c) Materials Response", font=list(size=18), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE))
#Add
lmiwplot_latex <- lmiwplot %>% layout(annotations=latexannotationsw)
lmiwkqplot_latex <- lmiwkqplot %>% layout(annotations=latexannotationsw)
#Plot
lmiwplot_latex
lmiwkqplot_latex
#JSON
lmiwplot <- lmiwplot %>% layout(annotations=annotationsw)
lmiwkqplot <- lmiwkqplot %>% layout(annotations=annotationsw)
#Save
lmiwplot   <- plotly_json(lmiwplot, FALSE)
lmiwkqplot   <- plotly_json(lmiwkqplot, FALSE)
# write(lmiwplot, "/Users/justindoty/Documents/Home/My_Website/static/jmp/main/lmiw3dplotly.json")
# write(lmiwkqplot, "/Users/justindoty/Documents/Home/My_Website/static/jmp/mainext/lmiwkq3dplotly.json")
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
# write(omgplotly, "/Users/justindoty/Documents/Home/My_Website/static/jmp/main/omgplotly.json")
#Side by Side Persistence Plots
omgplotly3 <- plot_ly(x=vectau, y=vectau, z=omg3dq, colorscale="Jet", type="surface", showscale=FALSE, scene="scene3", name=" ", hovertemplate = paste("<i>𝛕-innovation<i>: %{x:.2f}", "<br>𝛕-productivity: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% layout(scene3=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-innovation", titlefont=list(size=20), tickfont=list(size=16)), yaxis=list(title="𝛕-productivity", titlefont=list(size=20), tickfont=list(size=16)), zaxis=list(title="Persistence", titlefont=list(size=20), tickfont=list(size=16)))) 

persannotate <- list(list(x=0.09, y=0.75, text="(a) Linear Model", font=list(size=30, family="Times New Roman"), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
	list(x=0.5, y=0.75, text="(b) Nonlinear Model", font=list(size=30, family="Times New Roman"), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
	list(x=0.92, y=0.75, text="(c) Nonseparable Model", font=list(size=30, family="Times New Roman"), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE))

perspanel <- subplot(lmomega_plot, nlomegamat_plot, omgplotly3, shareX=TRUE) %>% layout(scene=list(aspectratio=list(x=.6, y=.6, z=.6), xaxis=list(title="𝛕-innovation", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="𝛕-productivity", titlefont=list(size=18), tickfont=list(size=14)), zaxis=list(title="Persistence", titlefont=list(size=18), tickfont=list(size=14))), 
	scene2=list(aspectratio=list(x=.6, y=.6, z=.6), xaxis=list(title="𝛕-innovation", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="𝛕-productivity", titlefont=list(size=18), tickfont=list(size=14)), zaxis=list(title="Persistence", titlefont=list(size=18), tickfont=list(size=14))),
	scene3=list(aspectratio=list(x=.6, y=.6, z=.6), xaxis=list(title="𝛕-innovation", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="𝛕-productivity", titlefont=list(size=18), tickfont=list(size=14)), zaxis=list(title="Persistence", titlefont=list(size=18), tickfont=list(size=14))))
perspanel %>% layout(annotations=persannotate)


















