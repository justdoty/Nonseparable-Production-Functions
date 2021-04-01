require(cowplot)
require(dplyr)
require(purrr)
require(xtable)
require(stringr)
source('/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Functions/Auxfuns.R')
source('/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Functions/Tensors.R')
#Load US Dataset
USdata <- read.csv('/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Data/USdata.csv') %>% 
select(id, year, Y, K, L, M, I, dep, naics3) %>% transmute(id=id, year=year, Y=log(Y), K=log(K), L=log(L), M=log(M), 
	I=log(I), dp=dep, naics3=as.character(naics3), naics2=as.numeric(str_extract(as.character(naics3), "^.{2}"))) %>% 
	group_by(id) %>% filter(n()>=3) %>% ungroup()
#Industries for Analysis
NAICS <- c("31", "32", "33", "All")
########################################################################################################
##########################################Summary Statistics############################################
########################################################################################################
#Create table for all relevant summary statistics
sumNAICS <- group_by(USdata, naics2) %>% summarise_at(c("Y", "K", "L", "M", "I"), list(Q1=~quantile(., 0.25), med=median, Q3=~quantile(.,0.75), mean=mean, sd=sd), na.rm=TRUE) 
sumALL <- cbind("All", summarise_at(USdata, c("Y", "K", "L", "M", "I"), list(Q1=~quantile(., 0.25), med=median, Q3=~quantile(.,0.75), mean=mean, sd=sd), na.rm=TRUE))
colnames(sumALL)[1] <- "naics2"
sizeNAICS <- group_by(USdata, naics2) %>% summarise(Firms=length(unique(id)), Total=n())
sizeALL <- c("All", sum(sizeNAICS$Firms), sum(sizeNAICS$Total))
size <- rbind(sizeNAICS, sizeALL)
sumstat <- round(matrix(as.numeric(as.matrix(rbind(sumNAICS, sumALL))[,-1]), nrow=20, ncol=5), 2)
#Some pretty formatting
NAICS_labels <- array(NA, 5*length(NAICS)); NAICS_labels[seq(1, 5*length(NAICS), by=5)] <- paste(NAICS, paste("(Total=", size$Total, ")", sep=""))
NAICS_labels[is.na(NAICS_labels)] <- ""
summary_table <- cbind(NAICS_labels, rep(c("Output", "Capital", "Labor", "Materials", "Investment"), 4), sumstat)
colnames(summary_table) <- c("Industry (NAICS code)", " ", "1st Qu.", 'Median', "3rd Qu.", 'Mean', "sd")

summary_table <- xtable(summary_table, digits=c(2,2,0,4,4,2,2,2), type="latex")
align(summary_table) <- rep('c', 8)
addtorow <- list()
addtorow$pos <- list(-1)
addtorow$command <- '\\hline\\hline '
#For copy pasting into latex
print(summary_table, hline.after=c(0,nrow(summary_table)), add.to.row=addtorow, auto=FALSE, include.rownames=FALSE, sanitize.text.function=function(x) x, table.placement="H")
#Saves to file
print(summary_table, hline.after=c(0,nrow(summary_table)), add.to.row=addtorow, auto=FALSE, include.rownames=FALSE, sanitize.text.function=function(x) x, table.placement="H", file="/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Empirical/US/Results/Tex/US_Summary.tex")
########################################################################################################
##########################################Load Results############################################
########################################################################################################
NAICS <- c("31", "32", "33", "^3")
industries <- c("31", "32", "33", "All")
for (i in 1:length(NAICS)){
	load(sprintf("/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Empirical/US/Environments/NLPFQR_NAICS_%s.RData", i))
	US <- filter(USdata, str_detect(naics2, NAICS[i]))
	vectau <- results$vectau
	ntau <- length(vectau)
	dims <- results$dims
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
	##############################################################################
	#Simulate Productivity Given Model Parameters
	##############################################################################
	#We simuluate a balanced panel even though the original model is estimated from
	#an unbalanced panel. In a later version, we consider adding a selection bias
	#correction to the productivity equation and drop firms according to this rule
	#in the simulated model
	#############################################################################
	N <- 1000
	T <- length(unique(US$year))
	#Capital
	lnkdata <- matrix(0, N, T)
	#Measurement Error in Capital
	upsdata <- matrix(rnorm(N*T), nrow=N, ncol=T)
	#Investment
	lnidata <- matrix(0, N, T)
	#Labor
    lnldata <- matrix(0, N, T) 
    #Shocks to Labor
    epsdata <- matrix(runif(N*T), nrow=N, ncol=T)
    #Intermediate Input
    lnmdata <- matrix(0, N, T)
    #Shock to Intermediate Input
    varepsdata <- matrix(runif(N*T), nrow=N, ncol=T)
    #Output
    lnydata <- matrix(0, N, T)
    #Shocks to Output 
    etadata <- matrix(runif(N*T), nrow=N, ncol=T)
    #Productivity
    omgdata <- matrix(0, N, T)
    #Shocks to Productivity
    xidata <- matrix(runif(N*T), nrow=N, ncol=T)
    #Aggregate depreciations rates from data
    dp <- US %>% group_by(year) %>% summarise(dp=mean(dp))
	#Initial Productivity
	omgdata[,1] <- rnorm(N, mean=parW1[1], sd=sqrt(parW1[2]))
	#Initial Capital
	lnkdata[,1] <- log(runif(N, min=2, max=4))
	#Initial Investment
	lnidata[,1] <- cbind(1, lnkdata[,1], omgdata[,1], lnkdata[,1]*omgdata[,1], lnkdata[,1]^2, omgdata[,1]^2)%*%as.matrix(parI[-length(parI)])
	for (t in 2:T){
		omgdata[,t] <- rowSums(cbind(1, omgdata[,t-1], omgdata[,t-1]^2, omgdata[,t-1]^3)*lspline(vectau=vectau, bvec=parWT, b1=parWTb[1], bL=parWTb[2], u=xidata[,t]))
		lnidata[,t] <- cbind(1, lnkdata[,t], omgdata[,t], lnkdata[,t]*omgdata[,t], lnkdata[,t]^2, omgdata[,t]^2)%*%as.matrix(parI[-length(parI)])+rnorm(N, 0, sd=sqrt(parI[length(parI)]))
		lnkdata[,t] <- log(dp$dp[t-1]*exp(lnkdata[,t-1])+exp(lnidata[,t-1])+exp(upsdata[,t-1]))
	}
	#Labor and Materials
	for (s in 1:T){
		lnldata[,s] <- rowSums(cbind(1, lnkdata[,s], omgdata[,s], lnkdata[,s]*omgdata[,s], lnkdata[,s]^2, omgdata[,s]^2)*lspline(vectau=vectau, bvec=parL, b1=parLb[1], bL=parLb[2], u=epsdata[,s]))
		lnmdata[,s] <- rowSums(cbind(1, lnkdata[,s], omgdata[,s], lnkdata[,s]*omgdata[,s], lnkdata[,s]^2, omgdata[,s]^2)*lspline(vectau=vectau, bvec=parM, b1=parMb[1], bL=parMb[2], u=varepsdata[,s]))
		lnydata[,s] <- rowSums(translog(K=lnkdata[,s], L=lnldata[,s], M=lnmdata[,s], omega=omgdata[,s])*lspline(vectau=vectau, bvec=parY, b1=parYb[1], bL=parYb[2], u=etadata[,s]))+omgdata[,s]
	}
	#Vectorize
	out <- c(t(lnydata))
	cap <- c(t(lnkdata))
	lab <- c(t(lnldata))
	mat <- c(t(lnmdata))
	inv <- c(t(lnidata))
	omg <- c(t(omgdata))
	#######################
	#Calculate Average QMEs
	######################
	#Capital
	kpost <- c(2,5,7,8,11,14,16,17)
	kdat <- cbind(1, lab, mat, 2*cap, omg, omg*lab, omg*mat, 2*cap*omg)
	kaqme_data <- data.frame(tau=vectau, kaqme=apply(parY[kpost,], 2, function(x) mean(kdat%*%as.matrix(x))))
	kaqme_plot <- ggplot(kaqme_data, aes(x=tau)) + xlab(expression('percentile-'*tau)) + ylab("Capital") + geom_line(aes(y=kaqme))
	save_plot(paste("/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Empirical/US/Results/Plots/Elasticities/K_AQME_NAICS_", industries[i], ".png", sep=""), kaqme_plot)
	#Labor
	lpost <- c(3,5,6,9,12,14,15,18)
	ldat <- cbind(1, cap, mat, 2*lab, omg, omg*cap, omg*mat, 2*lab*omg)
	laqme_data <- data.frame(tau=vectau, laqme=apply(parY[lpost,], 2, function(x) mean(ldat%*%as.matrix(x))))
	laqme_plot <- ggplot(laqme_data, aes(x=tau)) + xlab(expression('percentile-'*tau)) + ylab("Labor") + geom_line(aes(y=laqme))
	save_plot(paste("/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Empirical/US/Results/Plots/Elasticities/L_AQME_NAICS_", industries[i], ".png", sep=""), laqme_plot)
	#Materials
	mpost <- c(4,6,7,10,13,15,16,19)
	mdat <- cbind(1, lab, cap, 2*mat, omg, omg*lab, omg*cap, 2*mat*omg)
	maqme_data <- data.frame(tau=vectau, maqme=apply(parY[mpost,], 2, function(x) mean(mdat%*%as.matrix(x))))
	maqme_plot <- ggplot(maqme_data, aes(x=tau)) + xlab(expression('percentile-'*tau)) + ylab("Materials") + geom_line(aes(y=maqme))
	save_plot(paste("/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Empirical/US/Results/Plots/Elasticities/M_AQME_NAICS_", industries[i], ".png", sep=""), maqme_plot)
	#Hicks-Capital
	hkpost <- c(11,14,16,17)
	hkdat <- cbind(1, lab, mat, 2*cap)
	hkaqme_data <- data.frame(tau=vectau, hkaqme=apply(parY[hkpost,], 2, function(x) mean(hkdat%*%as.matrix(x))))
	hkaqme_plot <- ggplot(hkaqme_data, aes(x=tau)) + xlab(expression('percentile-'*tau)) + ylab("Capital-Productivity") + geom_line(aes(y=hkaqme))
	save_plot(paste("/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Empirical/US/Results/Plots/Hicks/HICKS_K_AQME_NAICS_", industries[i], ".png", sep=""), hkaqme_plot)
	#Hicks-Labor
	hlpost <- c(12,14,15,18)
	hldat <- cbind(1, cap, mat, 2*lab)
	hlaqme_data <- data.frame(tau=vectau, hlaqme=apply(parY[hlpost,], 2, function(x) mean(hldat%*%as.matrix(x))))
	hlaqme_plot <- ggplot(hlaqme_data, aes(x=tau)) + xlab(expression('percentile-'*tau)) + ylab("Labor-Productivity") + geom_line(aes(y=hlaqme))
	save_plot(paste("/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Empirical/US/Results/Plots/Hicks/HICKS_L_AQME_NAICS_", industries[i], ".png", sep=""), hlaqme_plot)
	#Hicks-Labor
	hmpost <- c(13,15,16,19)
	hmdat <- cbind(1, lab, cap, 2*mat)
	hmaqme_data <- data.frame(tau=vectau, hmaqme=apply(parY[hmpost,], 2, function(x) mean(hmdat%*%as.matrix(x))))
	hmaqme_plot <- ggplot(hmaqme_data, aes(x=tau)) + xlab(expression('percentile-'*tau)) + ylab("Materials-Productivity") + geom_line(aes(y=hmaqme))
	save_plot(paste("/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Empirical/US/Results/Plots/Hicks/HICKS_M_AQME_NAICS_", industries[i], ".png", sep=""), hmaqme_plot)
	#Productivity Persistence
	omglag <- c(t(omgdata[,2:T]))
	omgdata <- cbind(1, 2*omglag, 3*omglag^2)
	omgaqme_data <- data.frame(tau=vectau, omgaqme=apply(parWT[-1,], 2, function(x) mean(omgdata%*%as.matrix(x))))
	omgaqme_plot <- ggplot(omgaqme_data, aes(x=tau)) + xlab(expression('percentile-'*tau)) + ylab("Average Persistence") + geom_line(aes(y=omgaqme))
	save_plot(paste("/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Empirical/US/Results/Plots/TFP/OMG_AQME_NAICS_", industries[i], ".png", sep=""), omgaqme_plot)
	###############
	#Density Plots
	##############
	omgdens_plot <- ggplot(data.frame(omg), aes(x=omg)) + geom_density() + xlab("Productivity") + ylab("")
	save_plot(paste("/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Empirical/US/Results/Plots/TFP/OMG_DENS_NAICS_", industries[i], ".png", sep=""), omgdens_plot)
}




















