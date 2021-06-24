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
select(id, year, Y, K, K2, L, M, I, age) %>% transmute(id=id, year=year, Y=log(Y), K=log(K2), L=log(L), M=log(M), I=log(I), A=age)
omegainit <- as.matrix(omega_est(idvar=US$id, timevar=US$year, Y=US$Y, A=US$A, K=US$K, L=US$L, M=US$M))

########################################################################################################
##########################################Summary Statistics############################################
########################################################################################################
#Create table for all relevant summary statistics
# sumNAICS <- group_by(USdata, naics2) %>% summarise_at(c("Y", "K", "L", "M", "I"), list(Q1=~quantile(., 0.25), med=median, Q3=~quantile(.,0.75), mean=mean, sd=sd), na.rm=TRUE) 
# sumALL <- cbind("All", summarise_at(USdata, c("Y", "K", "L", "M", "I"), list(Q1=~quantile(., 0.25), med=median, Q3=~quantile(.,0.75), mean=mean, sd=sd), na.rm=TRUE))
# colnames(sumALL)[1] <- "naics2"
# sizeNAICS <- group_by(USdata, naics2) %>% summarise(Firms=length(unique(id)), Total=n())
# sizeALL <- c("All", sum(sizeNAICS$Firms), sum(sizeNAICS$Total))
# size <- rbind(sizeNAICS, sizeALL)
# sumstat <- round(matrix(as.numeric(as.matrix(rbind(sumNAICS, sumALL))[,-1]), nrow=20, ncol=5), 2)
# #Some pretty formatting
# NAICS_labels <- array(NA, 5*length(NAICS)); NAICS_labels[seq(1, 5*length(NAICS), by=5)] <- paste(NAICS, paste("(Total=", size$Total, ")", sep=""))
# NAICS_labels[is.na(NAICS_labels)] <- ""
# summary_table <- cbind(NAICS_labels, rep(c("Output", "Capital", "Labor", "Materials", "Capital"), 4), sumstat)
# colnames(summary_table) <- c("Industry (NAICS code)", " ", "1st Qu.", 'Median', "3rd Qu.", 'Mean', "sd")

# summary_table <- xtable(summary_table, digits=c(2,2,0,4,4,2,2,2), type="latex")
# align(summary_table) <- rep('c', 8)
# addtorow <- list()
# addtorow$pos <- list(-1)
# addtorow$command <- '\\hline\\hline '
# #For copy pasting into latex
# print(summary_table, hline.after=c(0,nrow(summary_table)), add.to.row=addtorow, auto=FALSE, include.rownames=FALSE, sanitize.text.function=function(x) x, table.placement="H")
# #Saves to file
# print(summary_table, hline.after=c(0,nrow(summary_table)), add.to.row=addtorow, auto=FALSE, include.rownames=FALSE, sanitize.text.function=function(x) x, table.placement="H", file="/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Empirical/US/Results/Translog/Tex/US_Summary.tex")
########################################################################################################
##########################################Load Results############################################
########################################################################################################
set.seed(123456)
load("/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Empirical/Capital/Capital_trans.RData")
vectau <- results$vectau
ntau <- length(vectau)
dims <- results$dims
#Load Parameter Estimates
parY <- results$resYmat
parL <- results$resLmat
parM <- results$resMmat
parKT <- results$resKTmat
parWT <- results$resWTmat
parW1 <- results$resW1mat
parK1 <- results$resK1mat
parYb <- results$resyb1bLmat
parLb <- results$reslb1bLmat
parMb <- results$resmb1bLmat
parWTb <- results$reswtb1bLmat
parW1b <- results$resw1tb1bLmat
WTminmax <- results$maxminwtmat
#Load the method used (e.g. Cobb, Translog, Hermite)
method <- results$method
resY <- results$resY
##############################################################################
#EM Chain Diagnostics
##############################################################################
#Plot Parameters for K, L, and M over EM iterations for each tau
EMcolour <- brewer.pal(n=length(vectau), "Spectral")
#Capital
K_EM_dat <- melt(resY[,3,])
K_EM_dat$Var2 <- as.factor(K_EM_dat$Var2)
K_EM <- ggplot(K_EM_dat, aes(x=Var1, y=value, group=Var2)) + geom_line(aes(colour=Var2)) + xlab("Iteration") + ylab("Capital") + scale_colour_manual(name="", labels=round(vectau, digits=2), values=EMcolour)
save_plot("/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Empirical/Capital/Plots/Diagnostics/K_EM.png", K_EM) 
#Labor
L_EM_dat <- melt(resY[,4,])
L_EM_dat$Var2 <- as.factor(L_EM_dat$Var2)
L_EM <- ggplot(L_EM_dat, aes(x=Var1, y=value, group=Var2)) + geom_line(aes(colour=Var2)) + xlab("Iteration") + ylab("Labor") + scale_colour_manual(name="", labels=round(vectau, digits=2), values=EMcolour)
save_plot("/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Empirical/Capital/Plots/Diagnostics/L_EM.png", L_EM) 
#Materials
M_EM_dat <- melt(resY[,5,])
M_EM_dat$Var2 <- as.factor(M_EM_dat$Var2)
M_EM <- ggplot(M_EM_dat, aes(x=Var1, y=value, group=Var2)) + geom_line(aes(colour=Var2)) + xlab("Iteration") + ylab("Materials") + scale_colour_manual(name="", labels=round(vectau, digits=2), values=EMcolour)
save_plot("/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Empirical/Capital/Plots/Diagnostics/M_EM.png", M_EM)
###############################################################################
#Simulate Productivity and Capital Evolution Given Model Parameters
##############################################################################
#We simuluate a balanced panel even though the original model is estimated from
#an unbalanced panel. In a later version, we consider adding a selection bias
#correction to the productivity equation and drop firms according to this rule
#in the simulated model
#############################################################################
N <- 10000
T <- length(unique(US$year))
#Age
adata <- matrix(0, N, T)
#Capital
lnkdata <- matrix(0, N, T)
#Labor
lnldata <- matrix(0, N, T) 
#Unobservable Shock to Labor
epsdata <- matrix(runif(N*T), nrow=N, ncol=T)
#Intermediate Input
lnmdata <- matrix(0, N, T)
#Unobservable Shock to Intermediate Inputs
varepsdata <- matrix(runif(N*T), nrow=N, ncol=T)
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
A1 <- log(t1data$A)
A <- exp(rtruncnorm(n=N, a=min(A1), b=max(A1), mean=mean(A1), sd=sd(A1)))
adata[,1] <- A
#Initial Productivity
omgdata[,1] <- rowSums(W1X(A=A)*lspline(vectau=vectau, bvec=parW1, b1=parW1b[1], bL=parW1b[2], u=xidata[,1]))
#Restrict the Support of Initial Productivity
omgdata[,1] <- (omgdata[,1]>WTminmax[1])*WTminmax[1]+(omgdata[,1]<WTminmax[2])*WTminmax[2]+(omgdata[,1]<=WTminmax[1])*(omgdata[,1]>=WTminmax[2])*omgdata[,1]
#Initial Capital
lnkdata[,1] <- rnorm(N, mean=K1X(A=A, omega=omgdata[,1])%*%parK1[1:(length(parK1)-1)], sd=sqrt(parK1[length(parK1)]))
#Restrict Support of Initial Capital
lnkdata[,1] <- (lnkdata[,1]>max(t1data$K))*max(t1data$K)+(lnkdata[,1]<min(t1data$K))*min(t1data$K)+(lnkdata[,1]<=max(t1data$K))*(lnkdata[,1]>=min(t1data$K))*lnkdata[,1]
#Labor
lnldata[,1] <- rowSums(LX(A=A, K=lnkdata[,1], omega=omgdata[,1])*lspline(vectau=vectau, bvec=parL, b1=parLb[1], bL=parLb[2], u=epsdata[,1]))
#Restrict Support of Initial Labor
lnldata[,1] <- (lnldata[,1]>max(t1data$L))*max(t1data$L)+(lnldata[,1]<min(t1data$L))*min(t1data$L)+(lnldata[,1]<=max(t1data$L))*(lnldata[,1]>=min(t1data$L))*lnldata[,1]
#Materials
lnmdata[,1] <- rowSums(MX(A=A, K=lnkdata[,1], omega=omgdata[,1])*lspline(vectau=vectau, bvec=parM, b1=parMb[1], bL=parMb[2], u=varepsdata[,1]))
#Restrict Support of Initial Materials
lnmdata[,1] <- (lnmdata[,1]>max(t1data$M))*max(t1data$M)+(lnmdata[,1]<min(t1data$M))*min(t1data$M)+(lnmdata[,1]<=max(t1data$M))*(lnmdata[,1]>=min(t1data$M))*lnmdata[,1]
#Initial Output
lnydata[,1] <- rowSums(PF(A=A, K=lnkdata[,1], L=lnldata[,1], M=lnmdata[,1], omega=omgdata[,1], method=method)*lspline(vectau=vectau, bvec=parY, b1=parYb[1], bL=parYb[2], u=etadata[,1]))
#Restrict Support of Initial Output
lnydata[,1] <- (lnydata[,1]>max(t1data$Y))*max(t1data$Y)+(lnydata[,1]<min(t1data$Y))*min(t1data$Y)+(lnydata[,1]<=max(t1data$Y))*(lnydata[,1]>=min(t1data$Y))*lnydata[,1]
for (t in 2:T){
	ttdata <- US %>% group_by(id) %>% slice(t)
	#Age
	A <- A+1
	adata[,t] <- A
	#Productivity
	omgdata[,t] <- rowSums(WX(A=A, omega=omgdata[,t-1])*lspline(vectau=vectau, bvec=parWT, b1=parWTb[1], bL=parWTb[2], u=xidata[,t]))
	#Restrict the Support of Productivity
	omgdata[,t] <- (omgdata[,t]>WTminmax[1])*WTminmax[1]+(omgdata[,t]<WTminmax[2])*WTminmax[2]+(omgdata[,t]<=WTminmax[1])*(omgdata[,t]>=WTminmax[2])*omgdata[,t]
	#Capital
	lnkdata[,t] <- rnorm(N, mean=KX(A=A, K=lnkdata[,t-1], omega=omgdata[,t-1])%*%parKT[1:(length(parKT)-1)], sd=sqrt(parKT[length(parKT)]))
	#Restrict Support of Capital
	lnkdata[,t] <- (lnkdata[,t]>max(ttdata$K))*max(ttdata$K)+(lnkdata[,t]<min(ttdata$K))*min(ttdata$K)+(lnkdata[,t]<=max(ttdata$K))*(lnkdata[,t]>=min(ttdata$K))*lnkdata[,t]
	#Labor
	lnldata[,t] <- rowSums(LX(A=A, K=lnkdata[,t], omega=omgdata[,t])*lspline(vectau=vectau, bvec=parL, b1=parLb[1], bL=parLb[2], u=epsdata[,t]))
	#Restrict Support of Labor
	lnldata[,t] <- (lnldata[,t]>max(ttdata$L))*max(ttdata$L)+(lnldata[,t]<min(ttdata$L))*min(ttdata$L)+(lnldata[,t]<=max(ttdata$L))*(lnldata[,t]>=min(ttdata$L))*lnldata[,t]
	#Materials
	lnmdata[,t] <- rowSums(MX(A=A, K=lnkdata[,t], omega=omgdata[,t])*lspline(vectau=vectau, bvec=parM, b1=parMb[1], bL=parMb[2], u=varepsdata[,t]))
	#Restrict Support of Materials
	lnmdata[,t] <- (lnmdata[,t]>max(ttdata$M))*max(ttdata$M)+(lnmdata[,t]<min(ttdata$M))*min(ttdata$M)+(lnmdata[,t]<=max(ttdata$M))*(lnmdata[,t]>=min(ttdata$M))*lnmdata[,t]
	#Output
	lnydata[,t] <- rowSums(PF(A=A, K=lnkdata[,t], L=lnldata[,t], M=lnmdata[,t], omega=omgdata[,t], method=method)*lspline(vectau=vectau, bvec=parY, b1=parYb[1], bL=parYb[2], u=etadata[,t]))
	#Restrict Support of Output
	lnydata[,t] <- (lnydata[,t]>max(ttdata$Y))*max(ttdata$Y)+(lnydata[,t]<min(ttdata$Y))*min(ttdata$Y)+(lnydata[,t]<=max(ttdata$Y))*(lnydata[,t]>=min(ttdata$Y))*lnydata[,t]
}

#Vectorize
out <- c(t(lnydata))
cap <- c(t(lnkdata))
lab <- c(t(lnldata))
mat <- c(t(lnmdata))
omg <- c(t(omgdata))
age <- c(t(adata))
##################################################
#Calculate Average QMEs of the Output Elasticities
##################################################
#Capital
kpost <- c(3,6,8,9,12,15,17,18)+1
kdat <- cbind(1, lab, mat, 2*cap, omg, omg*lab, omg*mat, 2*cap*omg)
kaqme_data <- data.frame(tau=vectau, kaqme=colMeans(apply(parY[kpost,], 2, function(x) kdat%*%x)))
kaqme_plot <- ggplot(kaqme_data, aes(x=tau, y=kaqme)) + xlab(expression('percentile-'*tau)) + ylab("Capital") + geom_line(aes(y=kaqme))
save_plot("/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Empirical/Capital/Plots/Elasticities/K_AQME.png", kaqme_plot)
#Labor
lpost <- c(4,6,7,10,13,15,16,19)+1
ldat <- cbind(1, cap, mat, 2*lab, omg, omg*cap, omg*mat, 2*lab*omg)
laqme_data <- data.frame(tau=vectau, laqme=colMeans(apply(parY[lpost,], 2, function(x) ldat%*%x)))
laqme_plot <- ggplot(laqme_data, aes(x=tau, y=laqme)) + xlab(expression('percentile-'*tau)) + ylab("Labor") + geom_line(aes(y=laqme))
save_plot("/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Empirical/Capital/Plots/Elasticities/L_AQME.png", laqme_plot)
#Materials
mpost <- c(5,7,8,11,14,16,17,20)+1
mdat <- cbind(1, lab, cap, 2*mat, omg, omg*lab, omg*cap, 2*mat*omg)
maqme_data <- data.frame(tau=vectau, maqme=colMeans(apply(parY[mpost,], 2, function(x) mdat%*%x)))
maqme_plot <- ggplot(maqme_data, aes(x=tau, y=maqme)) + xlab(expression('percentile-'*tau)) + ylab("Materials") + geom_line(aes(y=maqme))
save_plot("/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Empirical/Capital/Plots/Elasticities/M_AQME.png", maqme_plot)
#Hicks-Capital
hkpost <- c(12,16,17,18)+1
hkdat <- cbind(1, lab, mat, 2*cap)
hkaqme_data <- data.frame(tau=vectau, hkaqme=colMeans(apply(parY[hkpost,], 2, function(x) hkdat%*%x)))
hkaqme_plot <- ggplot(hkaqme_data, aes(x=tau, y=hkaqme)) + xlab(expression('percentile-'*tau)) + ylab("Capital-Productivity") + geom_line(aes(y=hkaqme))
save_plot("/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Empirical/Capital/Plots/Elasticities/HICKS_K_AQME.png", hkaqme_plot)
#Hicks-Labor
hlpost <- c(13,15,16,19)+1
hldat <- cbind(1, cap, mat, 2*lab)
hlaqme_data <- data.frame(tau=vectau, hlaqme=colMeans(apply(parY[hlpost,], 2, function(x) hldat%*%x)))
hlaqme_plot <- ggplot(hlaqme_data, aes(x=tau)) + xlab(expression('percentile-'*tau)) + ylab("Labor-Productivity") + geom_line(aes(y=hlaqme))
save_plot("/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Empirical/Capital/Plots/Elasticities/HICKS_L_AQME.png", hlaqme_plot)
#Hicks-Materials
hmpost <- c(14,16,17,20)+1
hmdat <- cbind(1, lab, cap, 2*mat)
hmaqme_data <- data.frame(tau=vectau, hmaqme=colMeans(apply(parY[hmpost,], 2, function(x) hmdat%*%x)))
hmaqme_plot <- ggplot(hmaqme_data, aes(x=tau)) + xlab(expression('percentile-'*tau)) + ylab("Materials-Productivity") + geom_line(aes(y=hmaqme))
save_plot("/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Empirical/Capital/Plots/Elasticities/HICKS_M_AQME.png", hmaqme_plot)
#Hicks-Productivity
hwpost <- c(2,12:20)+1
hwdat <- cbind(1, cap, lab, mat, cap*lab, lab*mat, cap*mat, cap^2, lab^2, mat^2)
hwaqme_data <- data.frame(tau=vectau, hwaqme=colMeans(apply(parY[hwpost,], 2, function(x) hwdat%*%x)))
hwaqme_plot <- ggplot(hwaqme_data, aes(x=tau)) + xlab(expression('percentile-'*tau)) + ylab("Productivity") + geom_line(aes(y=hwaqme))
save_plot("/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Empirical/Capital/Plots/Elasticities/HICKS_W_AQME.png", hwaqme_plot)
#3D Plots#########################################################################################
#Commands for Colors
nrz <- length(vectau)
ncz <- length(vectau)
jet.colors <-  colorRampPalette(c("midnightblue", "blue", "cyan","green", "yellow","orange","red", "darkred"))
nbcol <- 64
color <- jet.colors(nbcol)
#Capital##############################################################################################
k3d <- function(x){
	return(colMeans(cbind(1, lab, mat, 2*x, omg, omg*lab, omg*mat, 2*x*omg)))
}
k3dq <- t(sapply(vectau, function(q) k3d(quantile(cap, probs=q))))%*%parY[kpost,]
kfacet <- k3dq[-1,-1]+k3dq[-1,-ncz]+k3dq[-nrz,-1]+k3dq[-nrz,-ncz]
facetcol <- cut(kfacet, nbcol)
png("/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Empirical/Capital/Plots/Elasticities/3dK.png")
persp(x=vectau, y=vectau, z=k3dq, xlab="percentile-capital", ylab="percentile-output", zlab="Capital", col=color[facetcol], ticktype="detailed", phi=20,theta=-60)
dev.off()
#Labor##############################################################################################
l3d <- function(x){
	return(colMeans(cbind(1, cap, mat, 2*x, omg, omg*cap, omg*mat, 2*x*omg)))
}
l3dq <- t(sapply(vectau, function(q) l3d(quantile(lab, probs=q))))%*%parY[lpost,]
lfacet <- l3dq[-1,-1]+l3dq[-1,-ncz]+l3dq[-nrz,-1]+l3dq[-nrz,-ncz]
facetcol <- cut(lfacet, nbcol)
png("/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Empirical/Capital/Plots/Elasticities/3dL.png")
persp(x=vectau, y=vectau, z=l3dq, xlab="percentile-labor", ylab="percentile-output", zlab="Labor", col=color[facetcol], ticktype="detailed", phi=20,theta=-60)
dev.off()
#Materials##############################################################################################
m3d <- function(x){
	return(colMeans(cbind(1, lab, cap, 2*x, omg, omg*lab, omg*cap, 2*x*omg)))
}
m3dq <- t(sapply(vectau, function(q) m3d(quantile(mat, probs=q))))%*%parY[mpost,]
mfacet <- m3dq[-1,-1]+m3dq[-1,-ncz]+m3dq[-nrz,-1]+m3dq[-nrz,-ncz]
facetcol <- cut(mfacet, nbcol)
png("/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Empirical/Capital/Plots/Elasticities/3dM.png")
persp(x=vectau, y=vectau, z=m3dq, xlab="percentile-materials", ylab="percentile-output", zlab="Materials", col=color[facetcol], ticktype="detailed", phi=20,theta=-60)
dev.off()
#Hicks-Capital###########################################################################################
hk3d <- function(x){
	return(colMeans(cbind(1, lab, mat, 2*x)))
}
hk3dq <- t(sapply(vectau, function(q) hk3d(quantile(cap, probs=q))))%*%parY[hkpost,]
hkfacet <- hk3dq[-1,-1]+hk3dq[-1,-ncz]+hk3dq[-nrz,-1]+hk3dq[-nrz,-ncz]
facetcol <- cut(hkfacet, nbcol)
png("/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Empirical/Capital/Plots/Elasticities/3dhk.png")
persp(x=vectau, y=vectau, z=hk3dq, xlab="percentile-capital", ylab="percentile-output", zlab="Hicks-Capital", col=color[facetcol], ticktype="detailed", phi=20,theta=-60)
dev.off()
#Hicks-Labor###########################################################################################
hl3d <- function(x){
	return(colMeans(cbind(1, cap, mat, 2*x)))
}
hl3dq <- t(sapply(vectau, function(q) hl3d(quantile(lab, probs=q))))%*%parY[hlpost,]
hlfacet <- hl3dq[-1,-1]+hl3dq[-1,-ncz]+hl3dq[-nrz,-1]+hl3dq[-nrz,-ncz]
facetcol <- cut(hlfacet, nbcol)
png("/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Empirical/Capital/Plots/Elasticities/3dhl.png")
persp(x=vectau, y=vectau, z=hl3dq, xlab="percentile-labor", ylab="percentile-output", zlab="Hicks-Labor", col=color[facetcol], ticktype="detailed", phi=20,theta=-60)
dev.off()
#Hicks-Materials###########################################################################################
hm3d <- function(x){
	return(colMeans(cbind(1, lab, cap, 2*x)))
}
hm3dq <- t(sapply(vectau, function(q) hm3d(quantile(mat, probs=q))))%*%parY[hmpost,]
hmfacet <- hm3dq[-1,-1]+hm3dq[-1,-ncz]+hm3dq[-nrz,-1]+hm3dq[-nrz,-ncz]
facetcol <- cut(hmfacet, nbcol)
png("/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Empirical/Capital/Plots/Elasticities/3dhm.png")
persp(x=vectau, y=vectau, z=hm3dq, xlab="percentile-materials", ylab="percentile-output", zlab="Hicks-Materials", col=color[facetcol], ticktype="detailed", phi=20,theta=-60)
dev.off()
######################################################################################################################################################################################
#Productivity
#####################################################################################################################################################################################
#Productivity Persistence
omglag <- c(t(omgdata[,1:(T-1)]))
omgcon <- c(t(omgdata[,2:T]))
acon <- c(t(adata[,2:T]))
omgx <- sweep(WX(A=acon, omega=omglag)[,-c(2,dims$W)], 2, c(1:3), "*")
omgaqme_data <- data.frame(tau=vectau, omgaqme=colMeans(apply(parWT[-c(1,2),], 2, function(x) omgx%*%x)))
omgaqme_plot <- ggplot(omgaqme_data, aes(x=tau)) + xlab(expression('percentile-'*tau)) + ylab("Average Persistence") + geom_line(aes(y=omgaqme))
save_plot("/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Empirical/Capital/Plots/TFP/OMG_AQME.png", omgaqme_plot)
#3D Productivity
omg3dq <- apply(omgx, 2, function(q) quantile(q, probs=vectau))%*%parWT[-c(1,2),]
omgfacet <- omg3dq[-1,-1]+omg3dq[-1,-ncz]+omg3dq[-nrz,-1]+omg3dq[-nrz,-ncz]
facetcol <- cut(omgfacet, nbcol)
png("/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Empirical/Capital/Plots/TFP/3dpers.png")
persp(x=vectau, y=vectau, z=omg3dq, xlab="percentile-productivity", ylab="percentile-innovation", zlab="Persistence", col=color[facetcol], ticktype="detailed", phi=20,theta=-60)
dev.off()
################
##Density Plots
###############
omgdat <- data.frame(omg)
omgdens_plot <- ggplot(omgdat, aes(x=omg)) + geom_density() + xlab("Productivity") + ylab("")
save_plot("/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Empirical/Capital/Plots/TFP/OMG_DENS.png", omgdens_plot)













