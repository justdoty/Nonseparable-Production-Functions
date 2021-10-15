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
source('/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Functions/Main/Auxfuns.R')
source('/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Functions/Main/Tensors.R')
source('/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Functions/Main/omega.R')
#Load US Dataset
US <- read.csv('/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Data/USdata.csv')  %>%
	select(id, year, lny, lnk1, lnk2, lnl, lnm, lni, rd) %>% transmute(id=id, time=year, Y=lny, K=lnk1, L=lnl, M=lnm, I=lni, R=rd)
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
#Expand the sample
Nsim <- 20
N <- length(unique(US$id))*Nsim
T <- length(unique(US$time))
#Vector of ranks of input demand functions (small, medium, large)
tauinp <- c(0.1, 0.5, 0.9)
#Shocks to productivity
tauxi <- c(0.1, 0.5, 0.9)
#Rank of initial productivity
tauinit <- c(0.1, 0.5, 0.9)
#Labor###############################################################
lnldata <- array(0, c(N, T, length(tauxi), length(tauinp)))
#At the median
lnlmed <- array(0, c(T, length(tauxi), length(tauinp)))
#Shocks to Labor
epsl <- array(0, c(N, T, length(tauinp)))
#Intermediate Input#####################################################
lnmdata <- array(0, c(N, T, length(tauxi), length(tauinp)))
#At the median
lnmmed <- array(0, c(T, length(tauxi), length(tauinp)))
#Shocks to materials
epsm <- array(0, c(N, T, length(tauinp)))
#Investment###############################################################
lnidata <- array(0, c(N, T, length(tauxi)))
iota <- matrix(runif(N*T), nrow=N, ncol=T)
#This is for when i vary the shocks to investment demand
qlnidata <- array(0, c(N, T, length(tauxi), length(tauinp)))
qlnimed <- array(0, c(T, length(tauxi), length(tauinp)))
qiota <- array(0, c(N, T, length(tauinp)))
#Capital##################################################################
lnkdata <- array(0, c(N, T, length(tauxi)))
#This is for when i vary the shocks to investment demand
qlnkdata <- array(0, c(N, T, length(tauxi), length(tauinp)))
qlnkmed <- array(0, c(T, length(tauxi), length(tauinp)))
#Productivity#############################################################
#This is used for labor, materials, and investment
omgdata <- array(0, c(N, T, length(tauxi)))
#This is for when i vary productivity paths for different initial productivity levels
omgqdata <- array(0, c(N, T, length(tauxi), length(tauinit)))
omgqmed <- array(0, c(T, length(tauxi), length(tauinp)))
#Innovation Shocks
xi <- matrix(runif(N*T), nrow=N, ncol=T)
xidata <- replicate(length(tauxi), xi, simplify="array")
#Ranks for Input Shocks
for (q2 in 1:length(tauinp)){
	epsl[,,q2] <- array(tauinp[q2], c(N,T))
	epsm[,,q2] <- array(tauinp[q2], c(N,T))
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
#Median Shock
omg1m <- omg1q[2]
#Initial Investment (at median productivity)
i1 <- rowSums(IX(K=k1, omega=omg1m)*lspline(vectau=vectau, bvec=parI, b1=parIb[1], bL=parIb[2], u=iota[,1]))
#Restrict Support of Initial Investment
i1 <- (i1>max(t1data$I))*max(t1data$I)+(i1<min(t1data$I))*min(t1data$I)+(i1<=max(t1data$I))*(i1>=min(t1data$I))*i1
#Initial Inputs
for (q1 in 1:length(tauxi)){
	omgdata[,,q1][,1] <- omg1m
	lnkdata[,,q1][,1] <- k1
	lnidata[,,q1][,1] <- i1
	#At t=2 give varying shock to productivity given by tauxi
	xidata[,,q1][,2] <- tauxi[q1]
	for (q2 in 1:length(tauinp)){
		k1 <- rep(mean(k1), length(k1))
		l1 <- rowSums(LX(K=k1, omega=omg1m)*lspline(vectau=vectau, bvec=parL, b1=parLb[1], bL=parLb[2], u=epsl[,,q2][,1]))
		m1 <- rowSums(MX(K=k1, L=l1, omega=omg1m)*lspline(vectau=vectau, bvec=parM, b1=parMb[1], bL=parMb[2], u=epsm[,,q2][,1]))
		# #Restricting the Supports
		l1 <- (l1>max(t1data$L))*max(t1data$L)+(l1<min(t1data$L))*min(t1data$L)+(l1<=max(t1data$L))*(l1>=min(t1data$L))*l1
		m1 <- (m1>max(t1data$M))*max(t1data$M)+(m1<min(t1data$M))*min(t1data$M)+(m1<=max(t1data$M))*(m1>=min(t1data$M))*m1
		# #Median Across Firms
		lnldata[,,,q2][,,q1][,1] <- l1
		lnmdata[,,,q2][,,q1][,1] <- m1
		#Different levels of initial productivity
		omgqdata[,,,q2][,,q1][,1] <- omg1q[q2]
	}
}
for (t in 2:T){
	ttdata <- US %>% group_by(id) %>% slice(t)
	for (q1 in 1:length(tauxi)){
		omg <- rowSums(WX(omega=omgdata[,,q1][,t-1])*lspline(vectau=vectau, bvec=parWT, b1=parWTb[1], bL=parWTb[2], u=xidata[,,q1][,t]))
		#Restricting the Supports
		omgdata[,,q1][,t] <- (omg>wmax)*wmax+(omg<wmin)*wmin+(omg<=wmax)*(omg>=wmin)*omg
		#Generate Capital According to Accumulation Process with Depreciation
		lnk <- log(0.98*exp(lnkdata[,,q1][,t-1])+exp(lnidata[,,q1][,t-1]))
		#Restricting Supports
		lnkdata[,,q1][,t] <- (lnk>max(ttdata$K))*max(ttdata$K)+(lnk<min(ttdata$K))*min(ttdata$K)+(lnk<=max(ttdata$K))*(lnk>=min(ttdata$K))*lnk
		#Generate Investment
		lni <- rowSums(IX(K=lnkdata[,,q1][,t], omega=omgdata[,,q1][,t])*lspline(vectau=vectau, bvec=parI, b1=parIb[1], bL=parIb[2], u=iota[,t]))
		#Restricting Supports
		lnidata[,,q1][,t] <- (lni>max(ttdata$I))*max(ttdata$I)+(lni<min(ttdata$I))*min(ttdata$I)+(lni<=max(ttdata$I))*(lni>=min(ttdata$I))*lni
		#Generate Productivity
		for (q2 in 1:length(tauinp)){
			omgq <- rowSums(WX(omega=omgqdata[,,,q2][,,q1][,t-1])*lspline(vectau=vectau, bvec=parWT, b1=parWTb[1], bL=parWTb[2], u=xidata[,,q1][,t]))
			#Restricting the Supports
			omgqdata[,,,q2][,,q1][,t] <- (omgq>wmax)*wmax+(omgq<wmin)*wmin+(omgq<=wmax)*(omgq>=wmin)*omgq
		}
	}
}
#Generate Labor
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
#Generate Materials
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
#Median across firms
for (q1 in 1:length(tauxi)){
	for (q2 in 1:length(tauinp)){
	omgqmed[,,q2][,q1] <- apply(omgqdata[,,,q2][,,q1], 2, median)
	lnlmed[,,q2][,q1] <- apply(lnldata[,,,q2][,,q1], 2, median)
	lnmmed[,,q2][,q1] <- apply(lnmdata[,,,q2][,,q1], 2, median)
	}
}
#For the arrays, lnldata and lnmdata, the 1st dimension is time, 2nd is rank of innovation shock, 3rd is rank of input shock
#So "Low-Low" represents low labor shock and low innovation shock
labpath <- data.frame(1:T, lnlmed[,1,]-lnlmed[,2,], lnlmed[,3,]-lnlmed[,2,])
names(labpath) <- c("Time", "LowLow", "LowMed", "LowHigh", "HighLow", "HighMed", "HighHigh")
#Save to Selection_Bias Folder
write.csv(labpath, "/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Functions/Selection_Bias/labpath.csv")
#Paths for Materials at different ranks of materials shock
matpath <- data.frame(1:T, lnmmed[,1,]-lnmmed[,2,], lnmmed[,3,]-lnmmed[,2,])
names(matpath) <- c("Time", "LowLow", "LowMed", "LowHigh", "HighLow", "HighMed", "HighHigh")
#Save to Selection_Bias Folder
write.csv(matpath, "/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Functions/Selection_Bias/matpath.csv")
#Paths for Productivity at different ranks of initial productivity shock
omegapath <- data.frame(1:T, omgqmed[,1,]-omgqmed[,2,], omgqmed[,3,]-omgqmed[,2,])
names(omegapath) <- c("Time", "LowLow", "LowMed", "LowHigh", "HighLow", "HighMed", "HighHigh")
#Save to Selection_Bias Folder
write.csv(omegapath, "/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Functions/Selection_Bias/omegapath.csv")
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
Wplotly <- list()
Lplotly <- list()
Mplotly <- list()
wtitles <- list("Negative Shock<br>Low Productivity", "Negative Shock<br>Medium Productivity", "Negative Shock<br>High Productivity", "Positive Shock<br>Low Productivity", "Positive Shock<br>Medium Productivity", "Positive Shock<br>High Productivity")
ltitles <- list("Negative Shock<br>Low Labor", "Negative Shock<br>Medium Labor", "Negative Shock<br>High Labor", "Positive Shock<br>Low Labor", "Positive Shock<br>Medium Labor", "Positive Shock<br>High Labor")
mtitles <- list("Negative Shock<br>Low Materials", "Negative Shock<br>Medium Materials", "Negative Shock<br>High Materials", "Positive Shock<br>Low Materials", "Positive Shock<br>Medium Materials", "Positive Shock<br>High Materials")
for (i in 1:6){
	wdat <- data.frame(Time=omegapath$Time, Y=omegapath[,i+1])
	Wplotly[[i]] <- plot_ly(wdat, x=~Time, y=~Y, type = 'scatter', mode = 'lines', showlegend=F, line=list(color="black"), name=wtitles[[i]], hovertemplate = paste("<i>Year<i>: %{x}", "<br>Productivity Change: %{y:.2f}")) %>% layout(xaxis=list(title="Years"), yaxis=list(title="Productivity", range=list(min(omegapath[,-1]), max(omegapath[,-1]))), shapes=list(hline(y=0)))
	ldat <- data.frame(Time=labpath$Time, Y=labpath[,i+1])
	Lplotly[[i]] <- plot_ly(ldat, x=~Time, y=~Y, type = 'scatter', mode = 'lines', showlegend=F, line=list(color="black"), name=ltitles[[i]], hovertemplate = paste("<i>Year<i>: %{x}", "<br>Labor Change: %{y:.2f}")) %>% layout(xaxis=list(title="Years"), yaxis=list(title="Labor", range=list(min(labpath[,-1]), max(labpath[,-1]))), shapes=list(hline(y=0)))
	mdat <- data.frame(Time=matpath$Time, Y=matpath[,i+1])
	Mplotly[[i]] <- plot_ly(mdat, x=~Time, y=~Y, type = 'scatter', mode = 'lines', showlegend=F, line=list(color="black"), name=mtitles[[i]], hovertemplate = paste("<i>Year<i>: %{x}", "<br>Materials Change: %{y:.2f}")) %>% add_trace(y = 0, mode = 'lines', showlegend=F,line=list(color="red", dash="dash")) %>% layout(xaxis=list(title="Years"), yaxis=list(title="Materials", range=list(min(matpath[,-1]), max(matpath[,-1]))), shapes=list(hline(y=0)))
}
# Capital##################################################################################################
# I do the same procedure for capital except at different levels of shock to investment
# Initial Capital and Investment
for (q1 in 1:length(tauxi)){
	for (q2 in 1:length(tauinp)){
		qlnkdata[,,,q2][,,q1][,1] <- k1
		qlnkmed[,,q2][1,q1] <- median(qlnkdata[,,,q2][,,q1][,1])
		qi1 <- rowSums(IX(K=k1, omega=omg1m)*lspline(vectau=vectau, bvec=parI, b1=parIb[1], bL=parIb[2], u=qiota[,,q2][,1]))
		#Restrict the Support
		qi1 <- (qi1>max(t1data$I))*max(t1data$I)+(qi1<min(t1data$I))*min(t1data$I)+(qi1<=max(t1data$I))*(qi1>=min(t1data$I))*qi1
		#Initial Investment
		qlnidata[,,,q2][,,q1][,1] <- qi1
		qlnimed[,,q2][1,q1] <- median(qi1)
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
			#Restrict the Support of Investment
			qlnidata[,,,q2][,,q1][,t] <- (qit>max(ttdata$I))*max(ttdata$I)+(qit<min(ttdata$I))*min(ttdata$I)+(qit<=max(ttdata$I))*(qit>=min(ttdata$I))*qit
		}
	}
}
#Median across firms
for (q1 in 1:length(tauxi)){
	for (q2 in 1:length(tauinp)){
	qlnimed[,,q2][,q1] <- apply(qlnidata[,,,q2][,,q1], 2, median)
	qlnkmed[,,q2][,q1] <- apply(qlnkdata[,,,q2][,,q1], 2, median)
	}
}
# Capital
kpath <- data.frame(1:T, qlnkmed[,1,]-qlnkmed[,2,], qlnkmed[,3,]-qlnkmed[,2,])
names(kpath) <- c("Time", "LowLow", "LowMed", "LowHigh", "HighLow", "HighMed", "HighHigh")
write.csv(kpath, "/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Functions/Selection_Bias/kpath.csv")
ktitles <- list("Negative Shock<br>Low Investment", "Negative Shock<br>Medium Investment", "Negative Shock<br>High Investment", "Positive Shock<br>Low Investment", "Positive Shock<br>Medium Investment", "Positive Shock<br>High Investment")
#Plotting
Kplotly <- list()
for (i in 1:6){
	kdat <- data.frame(Time=kpath$Time, Y=kpath[,i+1])
	Kplotly[[i]] <- plot_ly(kdat, x=~Time, y=~Y, type = 'scatter', mode = 'lines', showlegend=F, line=list(color="black"), name=ktitles[[i]], hovertemplate = paste("<i>Year<i>: %{x}", "<br>Capital Change: %{y:.2f}")) %>% layout(xaxis=list(title="Years"), yaxis=list(title="Capital", range=list(min(kpath[,-1]), max(kpath[,-1]))), shapes=list(hline(y=0)))
}
##############################################################
#For Plot.ly
###############################################################
#Productivity
annotationsW <- list(list(x=0.11, y=0.9, text=TeX("(a) \\, \\tau_{\\xi}=0.1, \\tau_{\\omega_{1}}=0.1"), font=list(size=22), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
	list(x=0.495, y=0.9, text=TeX("(b)\\,\\tau_{\\xi}=0.1, \\tau_{\\omega_{1}}=0.5"), font=list(size=22), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.88, y=0.9, text=TeX("(c)\\,\\tau_{\\xi}=0.1, \\tau_{\\omega_{1}}=0.9"), font=list(size=22), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.11, y=0.45, text=TeX("(d)\\,\\tau_{\\xi}=0.9, \\tau_{\\omega_{1}}=0.1"), font=list(size=22), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.495, y=0.45, text=TeX("(e)\\,\\tau_{\\xi}=0.9, \\tau_{\\omega_{1}}=0.5"), font=list(size=22), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.88, y=0.45, text=TeX("(f)\\,\\tau_{\\xi}=0.9, \\tau_{\\omega_{1}}=0.9"), font=list(size=22), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE))
W <- subplot(Wplotly[[1]], Wplotly[[2]], Wplotly[[3]], Wplotly[[4]], Wplotly[[5]], Wplotly[[6]], shareX=TRUE, shareY=TRUE, nrows=2)
Wplot <- W %>% layout(annotations=annotationsW) %>% config(mathjax = 'cdn')
# Wplot
Wjson <- plotly_json(W, FALSE)
write(Wjson, "/Users/justindoty/Documents/Home/My_Website/static/jmp/main/impulseW.json")
#Labor 
annotationsL <- list(list(x=0.11, y=0.9, text=TeX("(a) \\, \\tau_{\\xi}=0.1, \\tau_{l}=0.1"), font=list(size=22), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
	list(x=0.495, y=0.9, text=TeX("(b)\\,\\tau_{\\xi}=0.1, \\tau_{l}=0.5"), font=list(size=22), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.88, y=0.9, text=TeX("(c)\\,\\tau_{\\xi}=0.1, \\tau_{l}=0.9"), font=list(size=22), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.11, y=0.45, text=TeX("(d)\\,\\tau_{\\xi}=0.9, \\tau_{l}=0.1"), font=list(size=22), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.495, y=0.45, text=TeX("(e)\\,\\tau_{\\xi}=0.9, \\tau_{l}=0.5"), font=list(size=22), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.88, y=0.45, text=TeX("(f)\\,\\tau_{\\xi}=0.9, \\tau_{l}=0.9"), font=list(size=22), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE))
L <- subplot(Lplotly[[1]], Lplotly[[2]], Lplotly[[3]], Lplotly[[4]], Lplotly[[5]], Lplotly[[6]], shareX=TRUE, shareY=TRUE, nrows=2)
Lplot <- L %>% layout(annotations=annotationsL) %>% config(mathjax = 'cdn')
# Lplot
Ljson <- plotly_json(L, FALSE)
write(Ljson, "/Users/justindoty/Documents/Home/My_Website/static/jmp/main/impulseL.json")
#Materials
annotationsM <- list(list(x=0.11, y=0.9, text=TeX("(a) \\, \\tau_{\\xi}=0.1, \\tau_{m}=0.1"), font=list(size=22), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
	list(x=0.495, y=0.9, text=TeX("(b)\\,\\tau_{\\xi}=0.1, \\tau_{m}=0.5"), font=list(size=22), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.88, y=0.9, text=TeX("(c)\\,\\tau_{\\xi}=0.1, \\tau_{m}=0.9"), font=list(size=22), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.11, y=0.45, text=TeX("(d)\\,\\tau_{\\xi}=0.9, \\tau_{m}=0.1"), font=list(size=22), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.495, y=0.45, text=TeX("(e)\\,\\tau_{\\xi}=0.9, \\tau_{m}=0.5"), font=list(size=22), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.88, y=0.45, text=TeX("(f)\\,\\tau_{\\xi}=0.9, \\tau_{m}=0.9"), font=list(size=22), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE))
M <- subplot(Mplotly[[1]], Mplotly[[2]], Mplotly[[3]], Mplotly[[4]], Mplotly[[5]], Mplotly[[6]], shareX=TRUE, shareY=TRUE, nrows=2)
Mplot <- M %>% layout(annotations=annotationsM) %>% config(mathjax = 'cdn')
# Mplot
Mjson <- plotly_json(M, FALSE)
write(Mjson, "/Users/justindoty/Documents/Home/My_Website/static/jmp/main/impulseM.json")
#Investment
annotationsI <- list(list(x=0.11, y=0.9, text=TeX("(a) \\, \\tau_{\\xi}=0.1, \\tau_{i}=0.1"), font=list(size=22), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
	list(x=0.495, y=0.9, text=TeX("(b)\\,\\tau_{\\xi}=0.1, \\tau_{i}=0.5"), font=list(size=22), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.88, y=0.9, text=TeX("(c)\\,\\tau_{\\xi}=0.1, \\tau_{i}=0.9"), font=list(size=22), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.11, y=0.47, text=TeX("(d)\\,\\tau_{\\xi}=0.9, \\tau_{i}=0.1"), font=list(size=22), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.495, y=0.47, text=TeX("(e)\\,\\tau_{\\xi}=0.9, \\tau_{i}=0.5"), font=list(size=22), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.88, y=0.47, text=TeX("(f)\\,\\tau_{\\xi}=0.9, \\tau_{i}=0.9"), font=list(size=22), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE))
#Capital
K <- subplot(Kplotly[[1]], Kplotly[[2]], Kplotly[[3]], Kplotly[[4]], Kplotly[[5]], Kplotly[[6]], shareX=TRUE, shareY=TRUE, nrows=2)
Kplot <- K %>% layout(annotations=annotationsI) %>% config(mathjax = 'cdn')
# Kplot
Kjson <- plotly_json(K, FALSE)
write(Kjson, "/Users/justindoty/Documents/Home/My_Website/static/jmp/main/impulseK.json")






