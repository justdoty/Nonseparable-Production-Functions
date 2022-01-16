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
load("/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Environments/bootmain_IRF.RData")
#Load US Dataset
US <- read.csv('/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Data/USdata.csv')  %>%
	select(id, year, lny, lnk1, lnk2, lnl, lnm, lni, rd) %>% transmute(id=id, time=year, Y=lny, K=lnk1, L=lnl, M=lnm, I=lni, R=rd)
T <- length(unique(US$time))	
alp <- 0.05
ntau <- 11
vectau <- seq(1/(ntau+1), ntau/(ntau+1), by=1/(ntau+1))
#Vector of ranks of input demand functions (small, medium, large)
tauinp <- c(0.1, 0.5, 0.9)
#Shocks to productivity
tauxi <- c(0.1, 0.5, 0.9)
#Rank of initial productivity
tauinit <- c(0.1, 0.5, 0.9)
##########################
#Lower Bound
##########################
#Misallocation
mpklowLB <- array(0, c(T, length(tauinit))); mpkhighLB <- mpklowLB
mpllowLB <- mpklowLB; mplhighLB <- mpklowLB
mpmlowLB <- mpklowLB; mpmhighLB <- mpklowLB
#Productivity
omglowLB <- mpklowLB; omghighLB <- mpklowLB
#Labor
lnllowLB <- array(0, c(T, length(tauinit), length(tauinp))); lnlhighLB <- lnllowLB
#Materials
lnmlowLB <- lnllowLB; lnmhighLB <- lnllowLB
#Capital
lnklowLB <- lnllowLB; lnkhighLB <- lnllowLB
#Investment
lnilowLB <- lnllowLB; lnihighLB <- lnllowLB
#Output
lnylowLB <- lnllowLB; lnyhighLB <- lnllowLB
##########################
#Upper Bound
##########################
mpklowUB <- array(0, c(T, length(tauinit))); mpkhighUB <- mpklowUB
mpllowUB <- mpklowUB; mplhighUB <- mpklowUB
mpmlowUB <- mpklowUB; mpmhighUB <- mpklowUB
#Productivity
omglowUB <- mpklowUB; omghighUB <- mpklowUB
#Labor
lnllowUB <- array(0, c(T, length(tauinit), length(tauinp))); lnlhighUB <- lnllowUB
#Materials
lnmlowUB <- lnllowUB; lnmhighUB <- lnllowUB
#Capital
lnklowUB <- lnllowLB; lnkhighUB <- lnllowLB
#Investment
lnilowUB <- lnllowLB; lnihighUB <- lnllowLB
#Output
lnylowUB <- lnllowLB; lnyhighUB <- lnllowLB
######################
#Confidence Bands
#####################
for (t in 1:T){
	for (q1 in 1:length(tauinit)){
		#MPK Low Innovation Shock
		mpklowLB[t,q1] <- as.numeric(quantile((bootresults$mpkboot[,,1,]-bootresults$mpkboot[,,2,])[t,,][q1,], probs=alp))
		mpklowUB[t,q1] <- as.numeric(quantile((bootresults$mpkboot[,,1,]-bootresults$mpkboot[,,2,])[t,,][q1,], probs=1-alp))
		#MPK High Innovation Shock
		mpkhighLB[t,q1] <- as.numeric(quantile((bootresults$mpkboot[,,3,]-bootresults$mpkboot[,,2,])[t,,][q1,], probs=alp))
		mpkhighUB[t,q1] <- as.numeric(quantile((bootresults$mpkboot[,,3,]-bootresults$mpkboot[,,2,])[t,,][q1,], probs=1-alp))
		#MPL Low Innovation Shock
		mpllowLB[t,q1] <- as.numeric(quantile((bootresults$mplboot[,,1,]-bootresults$mplboot[,,2,])[t,,][q1,], probs=alp))
		mpllowUB[t,q1] <- as.numeric(quantile((bootresults$mplboot[,,1,]-bootresults$mplboot[,,2,])[t,,][q1,], probs=1-alp))
		#MPL High Innovation Shock
		mplhighLB[t,q1] <- as.numeric(quantile((bootresults$mplboot[,,3,]-bootresults$mplboot[,,2,])[t,,][q1,], probs=alp))
		mplhighUB[t,q1] <- as.numeric(quantile((bootresults$mplboot[,,3,]-bootresults$mplboot[,,2,])[t,,][q1,], probs=1-alp))
		#MPM Low Innovation Shock
		mpmlowLB[t,q1] <- as.numeric(quantile((bootresults$mpmboot[,,1,]-bootresults$mpmboot[,,2,])[t,,][q1,], probs=alp))
		mpmlowUB[t,q1] <- as.numeric(quantile((bootresults$mpmboot[,,1,]-bootresults$mpmboot[,,2,])[t,,][q1,], probs=1-alp))
		#MPL High Innovation Shock
		mpmhighLB[t,q1] <- as.numeric(quantile((bootresults$mpmboot[,,3,]-bootresults$mpmboot[,,2,])[t,,][q1,], probs=alp))
		mpmhighUB[t,q1] <- as.numeric(quantile((bootresults$mpmboot[,,3,]-bootresults$mpmboot[,,2,])[t,,][q1,], probs=1-alp))
		#Productivity Low Innovation Shock
		omglowLB[t,q1] <- as.numeric(quantile((bootresults$omgboot[,,1,]-bootresults$omgboot[,,2,])[t,,][q1,], probs=alp))
		omglowUB[t,q1] <- as.numeric(quantile((bootresults$omgboot[,,1,]-bootresults$omgboot[,,2,])[t,,][q1,], probs=1-alp))
		#Productivity High Innovation Shock
		omghighLB[t,q1] <- as.numeric(quantile((bootresults$omgboot[,,3,]-bootresults$omgboot[,,2,])[t,,][q1,], probs=alp))
		omghighUB[t,q1] <- as.numeric(quantile((bootresults$omgboot[,,3,]-bootresults$omgboot[,,2,])[t,,][q1,], probs=1-alp))
		for (q2 in 1:length(tauinp)){
			#Labor Low Innovation Shock
			lnllowLB[,,q2][t,q1] <- as.numeric(quantile(bootresults$qlnlboot[,,1,,][,,q2,][t,,][q1,]-bootresults$qlnlboot[,,2,,][,,q2,][t,,][q1,], probs=alp))
			lnllowUB[,,q2][t,q1] <- as.numeric(quantile(bootresults$qlnlboot[,,1,,][,,q2,][t,,][q1,]-bootresults$qlnlboot[,,2,,][,,q2,][t,,][q1,], probs=1-alp))
			#Labor High Innovation Shock
			lnlhighLB[,,q2][t,q1] <- as.numeric(quantile(bootresults$qlnlboot[,,3,,][,,q2,][t,,][q1,]-bootresults$qlnlboot[,,2,,][,,q2,][t,,][q1,], probs=alp))
			lnlhighUB[,,q2][t,q1] <- as.numeric(quantile(bootresults$qlnlboot[,,3,,][,,q2,][t,,][q1,]-bootresults$qlnlboot[,,2,,][,,q2,][t,,][q1,], probs=1-alp))
			#Materials Low Innovation Shock
			lnmlowLB[,,q2][t,q1] <- as.numeric(quantile(bootresults$qlnmboot[,,1,,][,,q2,][t,,][q1,]-bootresults$qlnmboot[,,2,,][,,q2,][t,,][q1,], probs=alp))
			lnmlowUB[,,q2][t,q1] <- as.numeric(quantile(bootresults$qlnmboot[,,1,,][,,q2,][t,,][q1,]-bootresults$qlnmboot[,,2,,][,,q2,][t,,][q1,], probs=1-alp))
			#Materials High Innovation Shock
			lnmhighLB[,,q2][t,q1] <- as.numeric(quantile(bootresults$qlnmboot[,,3,,][,,q2,][t,,][q1,]-bootresults$qlnmboot[,,2,,][,,q2,][t,,][q1,], probs=alp))
			lnmhighUB[,,q2][t,q1] <- as.numeric(quantile(bootresults$qlnmboot[,,3,,][,,q2,][t,,][q1,]-bootresults$qlnmboot[,,2,,][,,q2,][t,,][q1,], probs=1-alp))
			#Capital Low Innovation Shock
			lnklowLB[,,q2][t,q1] <- as.numeric(quantile(bootresults$qlnkboot[,,1,,][,,q2,][t,,][q1,]-bootresults$qlnkboot[,,2,,][,,q2,][t,,][q1,], probs=alp))
			lnklowUB[,,q2][t,q1] <- as.numeric(quantile(bootresults$qlnkboot[,,1,,][,,q2,][t,,][q1,]-bootresults$qlnkboot[,,2,,][,,q2,][t,,][q1,], probs=1-alp))
			#Capital High Innovation Shock
			lnkhighLB[,,q2][t,q1] <- as.numeric(quantile(bootresults$qlnkboot[,,3,,][,,q2,][t,,][q1,]-bootresults$qlnkboot[,,2,,][,,q2,][t,,][q1,], probs=alp))
			lnkhighUB[,,q2][t,q1] <- as.numeric(quantile(bootresults$qlnkboot[,,3,,][,,q2,][t,,][q1,]-bootresults$qlnkboot[,,2,,][,,q2,][t,,][q1,], probs=1-alp))
			#Investment Low Innovation Shock
			lnilowLB[,,q2][t,q1] <- as.numeric(quantile(bootresults$qlniboot[,,1,,][,,q2,][t,,][q1,]-bootresults$qlniboot[,,2,,][,,q2,][t,,][q1,], probs=alp))
			lnilowUB[,,q2][t,q1] <- as.numeric(quantile(bootresults$qlniboot[,,1,,][,,q2,][t,,][q1,]-bootresults$qlniboot[,,2,,][,,q2,][t,,][q1,], probs=1-alp))
			#Investment High Innovation Shock
			lnihighLB[,,q2][t,q1] <- as.numeric(quantile(bootresults$qlniboot[,,3,,][,,q2,][t,,][q1,]-bootresults$qlniboot[,,2,,][,,q2,][t,,][q1,], probs=alp))
			lnihighUB[,,q2][t,q1] <- as.numeric(quantile(bootresults$qlniboot[,,3,,][,,q2,][t,,][q1,]-bootresults$qlniboot[,,2,,][,,q2,][t,,][q1,], probs=1-alp))
			#Output Low Innovation Shock
			lnylowLB[,,q2][t,q1] <- as.numeric(quantile(bootresults$qlnyboot[,,1,,][,,q2,][t,,][q1,]-bootresults$qlnyboot[,,2,,][,,q2,][t,,][q1,], probs=alp))
			lnylowUB[,,q2][t,q1] <- as.numeric(quantile(bootresults$qlnyboot[,,1,,][,,q2,][t,,][q1,]-bootresults$qlnyboot[,,2,,][,,q2,][t,,][q1,], probs=1-alp))
			#Output High Innovation Shock
			lnyhighLB[,,q2][t,q1] <- as.numeric(quantile(bootresults$qlnyboot[,,3,,][,,q2,][t,,][q1,]-bootresults$qlnyboot[,,2,,][,,q2,][t,,][q1,], probs=alp))
			lnyhighUB[,,q2][t,q1] <- as.numeric(quantile(bootresults$qlnyboot[,,3,,][,,q2,][t,,][q1,]-bootresults$qlnyboot[,,2,,][,,q2,][t,,][q1,], probs=1-alp))
		}

	}
}
#Misallocation
mispath <- list(Time=1:T, low_xi_MPK=list(LB=mpklowLB, UB=mpklowUB), low_xi_MPL=list(LB=mpllowLB, UB=mpllowUB), low_xi_MPM=list(LB=mpmlowLB, UB=mpmlowUB),
	high_xi_MPK=list(LB=mpkhighLB, UB=mpkhighUB), high_xi_MPL=list(LB=mplhighLB, UB=mplhighUB), high_xi_MPM=list(LB=mpmhighLB, UB=mpmhighUB))
#Productivity
omgpath <- list(Time=1:T, low_xi_low_W=list(LB=omglowLB[,1], UB=omglowUB[,1]), low_xi_med_W=list(LB=omglowLB[,2], UB=omglowUB[,2]), low_xi_high_W=list(LB=omglowLB[,3], UB=omglowUB[,3]),
	high_xi_low_W=list(LB=omghighLB[,1], UB=omghighUB[,1]), high_xi_med_W=list(LB=omghighLB[,2], UB=omghighUB[,2]), high_xi_high_W=list(LB=omghighLB[,3], UB=omghighUB[,3]))
#Labor
labpath <- list(Time=1:T, low_xi_low_L=list(LB=lnllowLB[,,1], UB=lnllowUB[,,1]), low_xi_med_L=list(LB=lnllowLB[,,2], UB=lnllowUB[,,2]), low_xi_high_L=list(LB=lnllowLB[,,3], UB=lnllowUB[,,3]),
	high_xi_low_L=list(LB=lnlhighLB[,,1], UB=lnlhighUB[,,1]), high_xi_med_L=list(LB=lnlhighLB[,,2], UB=lnlhighUB[,,2]), high_xi_high_L=list(LB=lnlhighLB[,,3], UB=lnlhighUB[,,3]))
#Materials
matpath <- list(Time=1:T, low_xi_low_L=list(LB=lnmlowLB[,,1], UB=lnmlowUB[,,1]), low_xi_med_L=list(LB=lnmlowLB[,,2], UB=lnmlowUB[,,2]), low_xi_high_L=list(LB=lnmlowLB[,,3], UB=lnmlowUB[,,3]),
	high_xi_low_L=list(LB=lnmhighLB[,,1], UB=lnmhighUB[,,1]), high_xi_med_L=list(LB=lnmhighLB[,,2], UB=lnmhighUB[,,2]), high_xi_high_L=list(LB=lnmhighLB[,,3], UB=lnmhighUB[,,3]))
#Capital
ipath <- list(Time=1:T, low_xi_low_L=list(LB=lnilowLB[,,1], UB=lnilowUB[,,1]), low_xi_med_L=list(LB=lnilowLB[,,2], UB=lnilowUB[,,2]), low_xi_high_L=list(LB=lnilowLB[,,3], UB=lnilowUB[,,3]),
	high_xi_low_L=list(LB=lnihighLB[,,1], UB=lnihighUB[,,1]), high_xi_med_L=list(LB=lnihighLB[,,2], UB=lnihighUB[,,2]), high_xi_high_L=list(LB=lnihighLB[,,3], UB=lnihighUB[,,3]))
#Output
outpath <- list(Time=1:T, low_xi_low_L=list(LB=lnylowLB[,,1], UB=lnylowUB[,,1]), low_xi_med_L=list(LB=lnylowLB[,,2], UB=lnylowUB[,,2]), low_xi_high_L=list(LB=lnylowLB[,,3], UB=lnylowUB[,,3]),
	high_xi_low_L=list(LB=lnyhighLB[,,1], UB=lnyhighUB[,,1]), high_xi_med_L=list(LB=lnyhighLB[,,2], UB=lnyhighUB[,,2]), high_xi_high_L=list(LB=lnyhighLB[,,3], UB=lnyhighUB[,,3]))
hline <- function(y = 0, color = "red"){list(type = "line", y0 = y, y1 = y, xref = "paper", x0 = 0, x1 = 1, line = list(color = color, dash="dot"))}
MISboot <- list()
Wboot <- list()
Lboot <- list()
Mboot <- list()
Iboot <- list()
Yboot <- list()
wcolors <- list("red", "green", "blue", "red", "green", "blue")
wshades <- list("rgba(100,0,0,0.3)", "rgba(0,100,0,0.3)", "rgba(0,0,100,0.3)", "rgba(100,0,0,0.3)", "rgba(0,100,0,0.3)", "rgba(0,0,100,0.3)")
for (i in 1:6){
	#Misallocation
	misdat <- data.frame(Time=mispath[[1]], low_W_LB=mispath[[i+1]]$LB[,1], low_W_UB=mispath[[i+1]]$UB[,1], med_W_LB=mispath[[i+1]]$LB[,2], med_W_UB=mispath[[i+1]]$UB[,2], high_W_LB=mispath[[i+1]]$LB[,3], high_W_UB=mispath[[i+1]]$UB[,3])
	MISboot[[i]] <- plot_ly(misdat, x=~Time, y=~low_W_LB, type = 'scatter', mode = 'lines', showlegend=F, line=list(color="red")) %>% add_trace(y=~low_W_UB, line=list(color="red"), fill="tonexty", fillcolor="rgba(100,0,0,0.3)") %>% add_trace(y=~med_W_LB, line=list(color="green")) %>% add_trace(y=~med_W_UB, line=list(color="green"), fill="tonexty", fillcolor="rgba(0,100,0,0.3)") %>% add_trace(y=~high_W_LB, line=list(color="blue")) %>% add_trace(y=~high_W_UB, line=list(color="blue"), fill="tonexty", fillcolor="rgba(0,0,100,0.3)")  %>% layout(xaxis=list(title="Years", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="Misallocation", titlefont=list(size=18), tickfont=list(size=14), range=list(min(unlist(mispath[names(mispath)!="Time"])), max(unlist(mispath[names(mispath)!="Time"])))), shapes=list(hline(y=0)))
	#Productivity
	omgdat <- data.frame(Time=omgpath[[1]], W_LB=omgpath[[i+1]]$LB, W_UB=omgpath[[i+1]]$UB)
	Wboot[[i]] <- plot_ly(omgdat, x=~Time, y=~W_LB, type = 'scatter', mode = 'lines', showlegend=F, line=list(color=wcolors[[i]]), hovertemplate = paste("<i>Year<i>: %{x}", "<br>Productivity Change: %{y:.3f}")) %>% add_trace(y=~W_UB, line=list(color=wcolors[[i]]), fill="tonexty", fillcolor=wshades[[i]]) %>% layout(xaxis=list(title="Years", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="Productivity", titlefont=list(size=18), tickfont=list(size=14), range=list(min(unlist(omgpath[names(omgpath)!="Time"])), max(unlist(omgpath[names(omgpath)!="Time"])))), shapes=list(hline(y=0)))
	#Labor
	labdat <- data.frame(Time=labpath[[1]], low_W_LB=labpath[[i+1]]$LB[,1], low_W_UB=labpath[[i+1]]$UB[,1], med_W_LB=labpath[[i+1]]$LB[,2], med_W_UB=labpath[[i+1]]$UB[,2], high_W_LB=labpath[[i+1]]$LB[,3], high_W_UB=labpath[[i+1]]$UB[,3])
	Lboot[[i]] <- plot_ly(labdat, x=~Time, y=~low_W_LB, type = 'scatter', mode = 'lines', showlegend=F, line=list(color="red")) %>% add_trace(y=~low_W_UB, line=list(color="red"), fill="tonexty", fillcolor="rgba(100,0,0,0.3)") %>% add_trace(y=~med_W_LB, line=list(color="green")) %>% add_trace(y=~med_W_UB, line=list(color="green"), fill="tonexty", fillcolor="rgba(0,100,0,0.3)") %>% add_trace(y=~high_W_LB, line=list(color="blue")) %>% add_trace(y=~high_W_UB, line=list(color="blue"), fill="tonexty", fillcolor="rgba(0,0,100,0.3)")  %>% layout(xaxis=list(title="Years", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="Labor", titlefont=list(size=18), tickfont=list(size=14), range=list(min(unlist(labpath[names(labpath)!="Time"])), max(unlist(labpath[names(labpath)!="Time"])))), shapes=list(hline(y=0)))
	#Materials
	matdat <- data.frame(Time=matpath[[1]], low_W_LB=matpath[[i+1]]$LB[,1], low_W_UB=matpath[[i+1]]$UB[,1], med_W_LB=matpath[[i+1]]$LB[,2], med_W_UB=matpath[[i+1]]$UB[,2], high_W_LB=matpath[[i+1]]$LB[,3], high_W_UB=matpath[[i+1]]$UB[,3])
	Mboot[[i]] <- plot_ly(matdat, x=~Time, y=~low_W_LB, type = 'scatter', mode = 'lines', showlegend=F, line=list(color="red")) %>% add_trace(y=~low_W_UB, line=list(color="red"), fill="tonexty", fillcolor="rgba(100,0,0,0.3)") %>% add_trace(y=~med_W_LB, line=list(color="green")) %>% add_trace(y=~med_W_UB, line=list(color="green"), fill="tonexty", fillcolor="rgba(0,100,0,0.3)") %>% add_trace(y=~high_W_LB, line=list(color="blue")) %>% add_trace(y=~high_W_UB, line=list(color="blue"), fill="tonexty", fillcolor="rgba(0,0,100,0.3)")  %>% layout(xaxis=list(title="Years", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="Materials", titlefont=list(size=18), tickfont=list(size=14), range=list(min(unlist(matpath[names(matpath)!="Time"])), max(unlist(matpath[names(matpath)!="Time"])))), shapes=list(hline(y=0)))
	#Capital
	idat <- data.frame(Time=ipath[[1]], low_W_LB=ipath[[i+1]]$LB[,1], low_W_UB=ipath[[i+1]]$UB[,1], med_W_LB=ipath[[i+1]]$LB[,2], med_W_UB=ipath[[i+1]]$UB[,2], high_W_LB=ipath[[i+1]]$LB[,3], high_W_UB=ipath[[i+1]]$UB[,3])
	Iboot[[i]] <- plot_ly(idat, x=~Time, y=~low_W_LB, type = 'scatter', mode = 'lines', showlegend=F, line=list(color="red")) %>% add_trace(y=~low_W_UB, line=list(color="red"), fill="tonexty", fillcolor="rgba(100,0,0,0.3)") %>% add_trace(y=~med_W_LB, line=list(color="green")) %>% add_trace(y=~med_W_UB, line=list(color="green"), fill="tonexty", fillcolor="rgba(0,100,0,0.3)") %>% add_trace(y=~high_W_LB, line=list(color="blue")) %>% add_trace(y=~high_W_UB, line=list(color="blue"), fill="tonexty", fillcolor="rgba(0,0,100,0.3)")  %>% layout(xaxis=list(title="Years", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="Investment", titlefont=list(size=18), tickfont=list(size=14), range=list(min(unlist(ipath[names(ipath)!="Time"])), max(unlist(ipath[names(ipath)!="Time"])))), shapes=list(hline(y=0)))
	#Output
	outdat <- data.frame(Time=outpath[[1]], low_W_LB=outpath[[i+1]]$LB[,1], low_W_UB=outpath[[i+1]]$UB[,1], med_W_LB=outpath[[i+1]]$LB[,2], med_W_UB=outpath[[i+1]]$UB[,2], high_W_LB=outpath[[i+1]]$LB[,3], high_W_UB=outpath[[i+1]]$UB[,3])
	Yboot[[i]] <- plot_ly(outdat, x=~Time, y=~low_W_LB, type = 'scatter', mode = 'lines', showlegend=F, line=list(color="red")) %>% add_trace(y=~low_W_UB, line=list(color="red"), fill="tonexty", fillcolor="rgba(100,0,0,0.3)") %>% add_trace(y=~med_W_LB, line=list(color="green")) %>% add_trace(y=~med_W_UB, line=list(color="green"), fill="tonexty", fillcolor="rgba(0,100,0,0.3)") %>% add_trace(y=~high_W_LB, line=list(color="blue")) %>% add_trace(y=~high_W_UB, line=list(color="blue"), fill="tonexty", fillcolor="rgba(0,0,100,0.3)")  %>% layout(xaxis=list(title="Years", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="Output", titlefont=list(size=18), tickfont=list(size=14), range=list(min(unlist(outpath[names(outpath)!="Time"])), max(unlist(outpath[names(outpath)!="Time"])))), shapes=list(hline(y=0)))

}
annotationsMis <- list(list(x=0.13, y=0.9, text=TeX("\\boldsymbol{(a) \\, MP_{k}, \\tau_{\\xi}=0.1}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
	list(x=0.495, y=0.9, text=TeX("\\boldsymbol{(b)\\, MP_{l}, \\tau_{\\xi}=0.1}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.88, y=0.9, text=TeX("\\boldsymbol{(c)\\, MP_{m}, \\tau_{\\xi}=0.1}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.13, y=0.48, text=TeX("\\boldsymbol{(d)\\, MP_{k}, \\tau_{\\xi}=0.9}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.495, y=0.48, text=TeX("\\boldsymbol{(e)\\, MP_{l}, \\tau_{\\xi}=0.9}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.88, y=0.48, text=TeX("\\boldsymbol{(f)\\, MP_{m}, \\tau_{\\xi}=0.9}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE))

MIS <- subplot(MISboot[[1]], MISboot[[2]], MISboot[[3]], MISboot[[4]], MISboot[[5]], MISboot[[6]], shareX=TRUE, shareY=TRUE, nrows=2)
Misplot <- MIS %>% layout(annotations=annotationsMis) %>% config(mathjax = 'cdn')
Misplot
Misjson <- plotly_json(MIS, FALSE)
# Productivity
annotationsW <- list(list(x=0.11, y=0.9, text=TeX("\\boldsymbol{(a) \\, \\tau_{\\xi}=0.1, \\tau_{\\omega_{1}}=0.1}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
	list(x=0.495, y=0.9, text=TeX("\\boldsymbol{(b)\\,\\tau_{\\xi}=0.1, \\tau_{\\omega_{1}}=0.5}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.88, y=0.9, text=TeX("\\boldsymbol{(c)\\,\\tau_{\\xi}=0.1, \\tau_{\\omega_{1}}=0.9}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.11, y=0.45, text=TeX("\\boldsymbol{(d)\\,\\tau_{\\xi}=0.9, \\tau_{\\omega_{1}}=0.1}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.495, y=0.45, text=TeX("\\boldsymbol{(e)\\,\\tau_{\\xi}=0.9, \\tau_{\\omega_{1}}=0.5}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.88, y=0.45, text=TeX("\\boldsymbol{(f)\\,\\tau_{\\xi}=0.9, \\tau_{\\omega_{1}}=0.9}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE))
W <- subplot(Wboot[[1]], Wboot[[2]], Wboot[[3]], Wboot[[4]], Wboot[[5]], Wboot[[6]], shareX=TRUE, shareY=TRUE, nrows=2)
Wplot <- W %>% layout(annotations=annotationsW) %>% config(mathjax = 'cdn')
Wplot
Wjson <- plotly_json(W, FALSE)
annotationsL <- list(list(x=0.11, y=0.95, text=TeX("\\boldsymbol{(a) \\, \\tau_{\\xi}=0.1, \\tau_{l}=0.1}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
	list(x=0.495, y=0.95, text=TeX("\\boldsymbol{(b)\\,\\tau_{\\xi}=0.1, \\tau_{l}=0.5}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.88, y=0.95, text=TeX("\\boldsymbol{(c)\\,\\tau_{\\xi}=0.1, \\tau_{l}=0.9}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.11, y=0.48, text=TeX("\\boldsymbol{(d)\\,\\tau_{\\xi}=0.9, \\tau_{l}=0.1}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.495, y=0.48, text=TeX("\\boldsymbol{(e)\\,\\tau_{\\xi}=0.9, \\tau_{l}=0.5}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.88, y=0.48, text=TeX("\\boldsymbol{(f)\\,\\tau_{\\xi}=0.9, \\tau_{l}=0.9}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE))
L <- subplot(Lboot[[1]], Lboot[[2]], Lboot[[3]], Lboot[[4]], Lboot[[5]], Lboot[[6]], shareX=TRUE, shareY=TRUE, nrows=2)
Lplot <- L %>% layout(annotations=annotationsL) %>% config(mathjax = 'cdn')
Lplot
Ljson <- plotly_json(L, FALSE)
annotationsM <- list(list(x=0.11, y=0.9, text=TeX("\\boldsymbol{(a) \\, \\tau_{\\xi}=0.1, \\tau_{m}=0.1}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
	list(x=0.495, y=0.9, text=TeX("\\boldsymbol{(b)\\,\\tau_{\\xi}=0.1, \\tau_{m}=0.5}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.88, y=0.9, text=TeX("\\boldsymbol{(c)\\,\\tau_{\\xi}=0.1, \\tau_{m}=0.9}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.11, y=0.45, text=TeX("\\boldsymbol{(d)\\,\\tau_{\\xi}=0.9, \\tau_{m}=0.1}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.495, y=0.45, text=TeX("\\boldsymbol{(e)\\,\\tau_{\\xi}=0.9, \\tau_{m}=0.5}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.88, y=0.45, text=TeX("\\boldsymbol{(f)\\,\\tau_{\\xi}=0.9, \\tau_{m}=0.9}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE))
M <- subplot(Mboot[[1]], Mboot[[2]], Mboot[[3]], Mboot[[4]], Mboot[[5]], Mboot[[6]], shareX=TRUE, shareY=TRUE, nrows=2)
Mplot <- M %>% layout(annotations=annotationsM) %>% config(mathjax = 'cdn')
Mplot
Mjson <- plotly_json(M, FALSE)
annotationsI <- list(list(x=0.13, y=0.9, text=TeX("\\boldsymbol{(a) \\, \\tau_{\\xi}=0.1, \\tau_{i}=0.1}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
	list(x=0.495, y=0.9, text=TeX("\\boldsymbol{(b)\\,\\tau_{\\xi}=0.1, \\tau_{i}=0.5}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.88, y=0.9, text=TeX("\\boldsymbol{(c)\\,\\tau_{\\xi}=0.1, \\tau_{i}=0.9}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.13, y=0.46, text=TeX("\\boldsymbol{(d)\\,\\tau_{\\xi}=0.9, \\tau_{i}=0.1}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.495, y=0.46, text=TeX("\\boldsymbol{(e)\\,\\tau_{\\xi}=0.9, \\tau_{i}=0.5}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.88, y=0.46, text=TeX("\\boldsymbol{(f)\\,\\tau_{\\xi}=0.9, \\tau_{i}=0.9}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE))
#Capital
I <- subplot(Iboot[[1]], Iboot[[2]], Iboot[[3]], Iboot[[4]], Iboot[[5]], Iboot[[6]], shareX=TRUE, shareY=TRUE, nrows=2)
Iplot <- I %>% layout(annotations=annotationsI) %>% config(mathjax = 'cdn')
Iplot
Ijson <- plotly_json(I, FALSE)
annotationsY <- list(list(x=0.13, y=0.9, text=TeX("\\boldsymbol{(a) \\, \\tau_{\\xi}=0.1, \\tau_{y}=0.1}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
	list(x=0.495, y=0.9, text=TeX("\\boldsymbol{(b)\\,\\tau_{\\xi}=0.1, \\tau_{y}=0.5}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.88, y=0.9, text=TeX("\\boldsymbol{(c)\\,\\tau_{\\xi}=0.1, \\tau_{y}=0.9}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.13, y=0.46, text=TeX("\\boldsymbol{(d)\\,\\tau_{\\xi}=0.9, \\tau_{y}=0.1}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.495, y=0.46, text=TeX("\\boldsymbol{(e)\\,\\tau_{\\xi}=0.9, \\tau_{y}=0.5}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.88, y=0.46, text=TeX("\\boldsymbol{(f)\\,\\tau_{\\xi}=0.9, \\tau_{y}=0.9}"), font=list(size=40), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE))
#Output
Y <- subplot(Yboot[[1]], Yboot[[2]], Yboot[[3]], Yboot[[4]], Yboot[[5]], Yboot[[6]], shareX=TRUE, shareY=TRUE, nrows=2)
Yplot <- Y %>% layout(annotations=annotationsY) %>% config(mathjax = 'cdn')
Yplot
Yjson <- plotly_json(Y, FALSE)



















