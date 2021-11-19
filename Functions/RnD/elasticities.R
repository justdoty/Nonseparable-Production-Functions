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
wmin <- min(US$Y)
wmax <- max(US$Y)
#Load the method used (e.g. Cobb, Translog, Hermite)
method <- results$method
resY <- results$resY
###############################################################################
#Simulate Productivity and Capital Evolution Given Model Parameters
##############################################################################
NR <- 15
NNR <- 10
#Simulate Data at T=1
t1data <- US %>% group_by(id) %>% slice(1)
#########################################################################################################
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
#Investment
lnidata <- matrix(0, nrow=N, ncol=T)
#Unobservable Shock to Intermediate Inputs
iotadata <- matrix(runif(N*T), nrow=N, ncol=T)
#Capital
lnkdata <- matrix(0, nrow=N, ncol=T)
#Productivity
omgdata <- matrix(0, nrow=N, ncol=T)
#RnD
rdata <- matrix(0, nrow=N, ncol=T)
#Unobservable Shock to Productivity
xidata <- matrix(runif(N*T), nrow=N, ncol=T)
#Unobservable Shocks to RnD
rhodata <- matrix(runif(N*T), nrow=N, ncol=T)
#I order Non RnD Firms First
lnkdata[,1] <- c(lnk1, rlnk1)
#Initial Productivity
omgdata[,1] <- rowSums(WX1(K=lnkdata[,1])*lspline(vectau=vectau, bvec=parW1, b1=parW1b[1], bL=parW1b[2], u=xidata[,1]))
#Restrict the Support of Initial Productivity
omgdata[,1] <- (omgdata[,1]>wmax)*wmax+(omgdata[,1]<wmin)*wmin+(omgdata[,1]<=wmax)*(omgdata[,1]>=wmin)*omgdata[,1]
#Initial Investment
lnidata[,1] <- rowSums(IX(K=lnkdata[,1], omega=omgdata[,1])*lspline(vectau=vectau, bvec=parI, b1=parIb[1], bL=parIb[2], u=iotadata[,1]))
#Restrict Support of Initial Investment
lnidata[,1] <- (lnidata[,1]>max(t1data$I))*max(t1data$I)+(lnidata[,1]<min(t1data$I))*min(t1data$I)+(lnidata[,1]<=max(t1data$I))*(lnidata[,1]>=min(t1data$I))*lnidata[,1]
#Initial RnD
rdata[pos,1] <- rowSums(RX(K=lnkdata[pos,1], omega=omgdata[pos,1])*lspline(vectau=vectau, bvec=parR, b1=parRb[1], bL=parRb[2], u=rhodata[pos,1]))
#Restrict Support of Initial Investment
rdata[pos,1] <- (rdata[pos,1]>max(t1data$R))*max(t1data$R)+(rdata[pos,1]<min(t1data$R[t1data$RB==1]))*min(t1data$R[t1data$RB==1])+(rdata[pos,1]<=max(t1data$R))*(rdata[pos,1]>=min(t1data$R[t1data$RB==1]))*rdata[pos,1]
#Evolution of Productivity and Capital
for (t in 2:T){
	ttdata <- US %>% group_by(id) %>% slice(t)
	rind <- rdata[,t-1]!=0
	#Capital
	lnkdata[,t] <- log(0.98*exp(lnkdata[,t-1])+exp(lnidata[,t-1]))
	#Restrict the Support of Capital 
	lnkdata[,t] <- (lnkdata[,t]>max(ttdata$K))*max(ttdata$K)+(lnkdata[,t]<min(ttdata$K))*min(ttdata$K)+(lnkdata[,t]<=max(ttdata$K))*(lnkdata[,t]>=min(ttdata$K))*lnkdata[,t]
	#Productivity
	omgdata[,t] <- rowSums(WX(omega=omgdata[,t-1], R=rdata[,t-1], Rind=rind)*lspline(vectau=vectau, bvec=parWT, b1=parWTb[1], bL=parWTb[2], u=xidata[,t]))
	#Restrict the Support of Productivity
	omgdata[,t] <- (omgdata[,t]>wmax)*wmax+(omgdata[,t]<wmin)*wmin+(omgdata[,t]<=wmax)*(omgdata[,t]>=wmin)*omgdata[,t]
	#Investment
	lnidata[,t] <- rowSums(IX(K=lnkdata[,t], omega=omgdata[,t])*lspline(vectau=vectau, bvec=parI, b1=parIb[1], bL=parIb[2], u=iotadata[,t]))
	#Restrict the Support of Investment
	lnidata[,t] <- (lnidata[,t]>max(ttdata$I))*max(ttdata$I)+(lnidata[,t]<min(ttdata$I))*min(ttdata$I)+(lnidata[,t]<=max(ttdata$I))*(lnidata[,t]>=min(ttdata$I))*lnidata[,t]
	#RnD
	rdata[pos,t] <- rowSums(RX(K=lnkdata[pos,t], omega=omgdata[pos,t])*lspline(vectau=vectau, bvec=parR, b1=parRb[1], bL=parRb[2], u=rhodata[pos,t]))
	#Restrict the Support of RnD
	rdata[pos,t] <- (rdata[pos,t]>max(ttdata$R))*max(ttdata$R)+(rdata[pos,t]<min(ttdata$R[ttdata$RB==1]))*min(ttdata$R[ttdata$RB==1])+(rdata[pos,t]<=max(ttdata$R))*(rdata[pos,t]>=min(ttdata$R[ttdata$RB==1]))*rdata[pos,t]
}
#Capital
cap <- c(t(lnkdata))
capcon <- c(t(lnkdata[,2:T]))
caplag <- c(t(lnkdata[,1:(T-1)]))
#Investment
inv <- c(t(lnidata))
invcon <- c(t(lnidata[,2:T]))
invlag <- c(t(lnidata[,1:(T-1)]))
#Productivity (all)
omg <- c(t(omgdata))
omgcon <- c(t(omgdata[,2:T]))
omglag <- c(t(omgdata[,1:(T-1)]))
#Unobservable Shock to Productivity
xi <- c(t(xidata))
xicon <- c(t(xidata[,2:T]))
#For Non RnD Firms
omgnr <- c(t(omgdata[-pos,]))
#For R&D Firms
omgr <- c(t(omgdata[pos,]))
#Productivity (RnD)
omgdatar <- omgdata[pos,]
omgr <- c(t(omgdatar))
omgrcon <- c(t(omgdatar[,2:T]))
omgrlag <- c(t(omgdatar[,1:(T-1)]))
#RnD (all)
rnd <- c(t(rdata))
rcon <- c(t(rdata[,2:T]))
rlag <- c(t(rdata[,1:(T-1)]))
rind <- rlag!=0
rho <- c(t(rhodata))
rhocon <- c(t(rhodata[,2:T]))
rholag <- c(t(rhodata[,1:(T-1)])) 
#Commands for Colors
nrz <- length(vectau)
ncz <- length(vectau)
jet.colors <-  colorRampPalette(c("midnightblue", "blue", "cyan","green", "yellow","orange","red", "darkred"))
nbcol <- 64
color <- jet.colors(nbcol)
######################################################################################################################################################################################
#Productivity
#####################################################################################################################################################################################
omg3dq <- array(0, c(ntau, ntau))
romg3dq <- omg3dq
omgr3dq <- omg3dq
rdw <- omg3dq
rdr <- omg3dq
for (q in 1:ntau){
	#Persistence of Productivity for Non-R&D Firms
	#Fix quantiles of previous period productivity and previous R&D levels
	omgq <- rep(as.numeric(quantile(omglag, probs=vectau[q])), length(omglag))
	rlagq <- rep(as.numeric(quantile(rlag, probs=vectau[q])), length(rlag))
	#Generate productivity at fixed levels of productivity for Non R&D Firms
	omgx <- rowSums(WX(omega=omgq, R=rlag, Rind=rind)*lspline(vectau=vectau, bvec=parWT, b1=parWTb[1], bL=parWTb[2], u=xicon))
	omg3dq[q,] <- colMeans(WXRD(omega=omgx, R=rlag, Rind=rind, par=parWT, pos=1, sdpos=1)$d1)
	#Persistence of Productivity for R&D Firms (at percentiles of productivity)
	#Generate productivity at fixed levels of productivity for R&D Firms
	rq <- rowSums(RX(omega=omgq, K=caplag)*lspline(vectau=vectau, bvec=parR, b1=parRb[1], bL=parRb[2], u=rholag))
	rlagm <- rep(mean(rq), length(rq))
	omgrx <- rowSums(WX(omega=omgq, R=rlagm, Rind=rind)*lspline(vectau=vectau, bvec=parWT, b1=parWTb[1], bL=parWTb[2], u=xicon))
	romg3dq[q,] <- colMeans(WXRD(omega=omgrx, R=rlagm, Rind=rind, par=parWT, pos=1, sdpos=1)$d2)
	#Persistence of Productivity for R&D Firms (at percentiles of R&D)
	wmean <- rep(mean(omglag), length(omglag))
	romgx <- rowSums(WX(omega=wmean, R=rlagq, Rind=rind)*lspline(vectau=vectau, bvec=parWT, b1=parWTb[1], bL=parWTb[2], u=xicon))
	omgr3dq[q,] <- colMeans(WXRD(omega=wmean, R=rlagq, Rind=rind, par=parWT, pos=1, sdpos=1)$d2)
	#Returns to R&D for fixed percentiles of productivity
	rdw[q,] <- colMeans(WXRD(omega=omgrx, R=rlagm, Rind=rind, par=parWT, pos=2, sdpos=1)$d2)
	#Returns to R&D for fixed percentiles of R&D
	rdr[q,] <- colMeans(WXRD(omega=wmean, R=rlagq, Rind=rind, par=parWT, pos=2, sdpos=1)$d2)
	
}
#Persistence for Non-R&D Firms
omg3d <- plot_ly(x=vectau, y=vectau, z=omg3dq, colorscale="Jet", type="surface", showscale=FALSE, scene="scene1", name=" ", hovertemplate = paste("No R&D<br><i>𝛕-innovation<i>: %{x:.2f}", "<br>𝛕-productivity: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% layout(scene1=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-innovation", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="𝛕-productivity", titlefont=list(size=18), tickfont=list(size=14)), zaxis=list(title="Persistence (No R&D)", titlefont=list(size=18), tickfont=list(size=14)))) 
#Persistence for R&D Firms (fixed percentiles of productivity)
romg3d <- plot_ly(x=vectau, y=vectau, z=romg3dq, colorscale="Jet", type="surface", showscale=FALSE, scene="scene2", name=" ", hovertemplate = paste("R&D<br><i>𝛕-innovation<i>: %{x:.2f}", "<br>𝛕-productivity: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% layout(scene2=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-innovation", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="𝛕-productivity", titlefont=list(size=18), tickfont=list(size=14)), zaxis=list(title="Persistence (R&D)", titlefont=list(size=18), tickfont=list(size=14)))) 
#Persistence for R&D Firms (fixed percentiles of R&D)
omgr3d <- plot_ly(x=vectau, y=vectau, z=omgr3dq, colorscale="Jet", type="surface", showscale=FALSE, scene="scene3", name=" ", hovertemplate = paste("R&D<br><i>𝛕-innovation<i>: %{x:.2f}", "<br>𝛕-R&D: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% layout(scene3=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-innovation", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="𝛕-R&D", titlefont=list(size=18), tickfont=list(size=14)), zaxis=list(title="Persistence (R&D)", titlefont=list(size=18), tickfont=list(size=14)))) 
#Combined Plots
persplotly <- subplot(omg3d, romg3d, omgr3d) %>% layout(scene=list(aspectratio=list(x=0.6, y=0.6, z=0.6), xaxis=list(title="𝛕-innovation", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="𝛕-productivity", titlefont=list(size=18), tickfont=list(size=14)), zaxis=list(title="Persistence", titlefont=list(size=18), tickfont=list(size=14))), 
	scene2=list(aspectratio=list(x=0.6, y=0.6, z=0.6), xaxis=list(title="𝛕-innovation", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="𝛕-productivity", titlefont=list(size=18), tickfont=list(size=14)), zaxis=list(title="Persistence", titlefont=list(size=18), tickfont=list(size=14))),
	scene3=list(aspectratio=list(x=0.6, y=0.6, z=0.6), xaxis=list(title="𝛕-innovation", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="𝛕-R&D", titlefont=list(size=18), tickfont=list(size=14)), zaxis=list(title="Persistence", titlefont=list(size=18), tickfont=list(size=14))))
#Annotations
latexannotations <- list(list(x=0.10, y=0.75, text="(a) Non-performers", font=list(size=30, family="Times New Roman"), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
	list(x=0.5, y=0.75, text="(b) Performers", font=list(size=30, family="Times New Roman"), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
	list(x=0.88, y=0.75, text="(c) Performers", font=list(size=30, family="Times New Roman"), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE))
annotations <- list(list(x=0.11, y=0.75, text="(a) Non-performers", font=list(size=18), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
	list(x=0.5, y=0.75, text="(b) Performers", font=list(size=18), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
	list(x=0.88, y=0.75, text="(c) Performers", font=list(size=18), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE))
latexannotationsab <- list(list(x=0.25, y=0.8, text="<b>(a)<b>", font=list(size=30, family="Times New Roman"), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
	list(x=0.75, y=0.8, text="<b>(b)<b>", font=list(size=30, family="Times New Roman"), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE))
annotationsab <- list(list(x=0.25, y=0.8, text="(a)", font=list(size=18), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
	list(x=0.75, y=0.8, text="(b)", font=list(size=18), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE))

#Add
persplotly_latex <- persplotly %>% layout(annotations=latexannotations)
persplotly_latex
persplotly <- persplotly %>% layout(annotations=annotations)
persplotly <- plotly_json(persplotly, FALSE)
write(persplotly, "/Users/justindoty/Documents/Home/My_Website/static/jmp/rnd/persplotly.json")
#Returns to R&D (fixed percentiles of productivity)
rdwq <- plot_ly(x=vectau, y=vectau, z=rdw, colorscale="Jet", type="surface", showscale=FALSE, scene="scene1", name=" ", hovertemplate = paste("<i>𝛕-innovation<i>: %{x:.2f}", "<br>𝛕-productivity: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% layout(scene1=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-innovation", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="𝛕-productivity", titlefont=list(size=18), tickfont=list(size=14)), zaxis=list(title="Returns to R&D", titlefont=list(size=18), tickfont=list(size=14)))) 
#Returns to R&D (fixed percentiles of R&D)
rdrq <- plot_ly(x=vectau, y=vectau, z=rdr, colorscale="Jet", type="surface", showscale=FALSE, scene="scene2", name=" ", hovertemplate = paste("<i>𝛕-innovation<i>: %{x:.2f}", "<br>𝛕-R&D: %{y:.2f}<br>", "Estimate: %{z:.3f}")) %>% layout(scene2=list(camera=list(eye=list(x=-1.5, y=-1.5, z=0.5)), aspectratio=list(x=1, y=1, z=1), xaxis=list(title="𝛕-innovation", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="𝛕-R&D", titlefont=list(size=18), tickfont=list(size=14)), zaxis=list(title="Returns to R&D", titlefont=list(size=18), tickfont=list(size=14)))) 
rdplot <- subplot(rdwq, rdrq) %>% layout(scene=list(aspectratio=list(x=0.9, y=0.9, z=0.9), xaxis=list(title="𝛕-innovation", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="𝛕-productivity", titlefont=list(size=18), tickfont=list(size=14)), zaxis=list(title="Returns", titlefont=list(size=18), tickfont=list(size=14))), 
	scene2=list(aspectratio=list(x=0.9, y=0.9, z=0.9), xaxis=list(title="𝛕-innovation"), yaxis=list(title="𝛕-R&D"), zaxis=list(title="Returns")))
rdplot_latex <- rdplot %>% layout(annotations=latexannotationsab)
rdplot_latex
rdplot <- rdplot %>% layout(annotations=annotationsab)
rdplot    <- plotly_json(rdplot, FALSE)
write(rdplot, "/Users/justindoty/Documents/Home/My_Website/static/jmp/rnd/rdplot.json")



