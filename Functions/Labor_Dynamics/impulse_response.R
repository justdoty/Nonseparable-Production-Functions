require(cowplot)
require(dplyr)
require(purrr)
require(xtable)
require(stringr)
require(truncnorm)
require(RColorBrewer)
require(reshape2)
require(plotly)
source('/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Functions/Labor_Dynamics/Auxfuns.R')
source('/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Functions//Labor_Dynamics/Tensors.R')
source('/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Functions//Labor_Dynamics/omega.R')
#Load US Dataset
US <- read.csv('/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Data/USdata.csv') %>% 
select(id, year, lny, lnk1, lnk2, lnl, lnm, lni, age) %>% transmute(id=id, year=year, Y=lny, K=lnk1, L=lnl, M=lnm, I=lni, A=age)
#Productivity from LP
omegainit <- omega_est(idvar=US$id, timevar=US$year, Y=US$Y, A=US$A, K=US$K, L=US$L, M=US$M)$omega
US <- US %>% mutate(omega=omegainit)
#De-mean
US <- US %>% mutate(Y=Y-mean(Y), K=K-mean(K), L=L-mean(L), M=M-mean(M), I=I-mean(I))
########################################################################################################
##########################################Load Results############################################
########################################################################################################
set.seed(123456)
load("/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Environments/labor_dynamics.RData")
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
#Load the method used (e.g. Cobb, Translog, Hermite)
method <- results$method
resY <- results$resY
wmin <- min(US$Y)
wmax <- max(US$Y)
###############################################################################
#Simulate Productivity and Capital Evolution Given Model Parameters
#############################################################################
#Expand the sample
Nsim <- 5
N <- length(unique(US$id))*Nsim
T <- length(unique(US$year))
#Shocks to Previous Period Labor Demand
taulag <- c(0.1, 0.25, 0.5, 0.75, 0.9)
#Size of Current Labor Demand
taul <- c(0.1, 0.5, 0.9)
#Labor###############################################################
lnldata <- array(0, c(N, T, length(taulag), length(taul)))
#At the median
lnlmed <- array(0, c(T, length(taulag), length(taul)))
#Shocks to Labor
epsl <- array(0, c(N, T, length(taulag), length(taul)))
#Investment###############################################################
lnidata <- array(0, c(N, T))
iota <- matrix(runif(N*T), nrow=N, ncol=T)
#Capital##################################################################
lnkdata <- array(0, c(N, T))
#Productivity#############################################################
omgdata <- array(0, c(N, T))
#Innovation Shocks
xidata <- matrix(runif(N*T), nrow=N, ncol=T)
#Ranks for Input Shocks
for (q2 in 1:length(taul)){
	epsl[,,,q2] <- array(taul[q2], c(N,T,length(taulag)))
}
#########################################################################################################
#Data at T=1
t1data <- US %>% group_by(id) %>% slice(1)
#Initial Capital
lnkdata[,1] <- kronecker(array(1, c(Nsim,1)), t1data$K)
#Initial Productivity
omg1 <- rowSums(WX1(K=lnkdata[,1], L=median(t1data$L))*lspline(vectau=vectau, bvec=parW1, b1=parW1b[1], bL=parW1b[2], u=xidata[,1]))
#Restrict the Support of Initial Productivity
omg1 <- (omg1>wmax)*wmax+(omg1<wmin)*wmin+(omg1<=wmax)*(omg1>=wmin)*omg1
#Median Shock
omgdata[,1] <- median(omg1)
#Initial Investment (at median productivity)
i1 <- rowSums(IX(K=lnkdata[,1], omega=omgdata[,1])*lspline(vectau=vectau, bvec=parI, b1=parIb[1], bL=parIb[2], u=iota[,1]))
#Restrict Support of Initial Investment
i1 <- (i1>max(t1data$I))*max(t1data$I)+(i1<min(t1data$I))*min(t1data$I)+(i1<=max(t1data$I))*(i1>=min(t1data$I))*i1
lnidata[,1] <- i1
#Initial Inputs
for (q1 in 1:length(taulag)){
	for (q2 in 1:length(taul)){
		#At t=2 give varying shock to labor given by taulag
		epsl[,,,q2][,,q1][,2] <- taulag[q1]
		lnldata[,,,q2][,,q1][,1] <- quantile(t1data$L, probs=taul[q2]) 
		lnlmed[,,q2][1,q1] <- median(lnldata[,,,q2][,,q1][,1])
	}
}
#This Loop is Pretty Slow
for (t in 2:T){
	ttdata <- US %>% group_by(id) %>% slice(t)
	omg <- rowSums(WX(omega=omgdata[,t-1])*lspline(vectau=vectau, bvec=parWT, b1=parWTb[1], bL=parWTb[2], u=xidata[,t]))
	#Restricting the Supports
	omgdata[,t] <- (omg>wmax)*wmax+(omg<wmin)*wmin+(omg<=wmax)*(omg>=wmin)*omg
	# Generate Capital According to Accumulation Process with Industry-Average Depreciation Rates
	lnk <- log(0.98*exp(lnkdata[,t-1])+exp(lnidata[,t-1]))
	#Restricting Supports
	lnkdata[,t] <- (lnk>max(ttdata$K))*max(ttdata$K)+(lnk<min(ttdata$K))*min(ttdata$K)+(lnk<=max(ttdata$K))*(lnk>=min(ttdata$K))*lnk
	#Generate Investment
	lni <- rowSums(IX(K=lnkdata[,t], omega=omgdata[,t])*lspline(vectau=vectau, bvec=parI, b1=parIb[1], bL=parIb[2], u=iota[,t]))
	#Restricting Supports
	lnidata[,t] <- (lni>max(ttdata$I))*max(ttdata$I)+(lni<min(ttdata$I))*min(ttdata$I)+(lni<=max(ttdata$I))*(lni>=min(ttdata$I))*lni
	for (q1 in 1:length(taulag)){
		#Generate Labor Following Shocks from previous period
		for (q2 in 1:length(taul)){
			lab <- rowSums(LX(K=lnkdata[,t], L=lnldata[,,,q2][,,q1][,t-1], omega=omgdata[,t])*lspline(vectau=vectau, bvec=parL, b1=parLb[1], bL=parLb[2], u=epsl[,,,q2][,,q1][,t]))
			#Restricting the Supports
			lnldata[,,,q2][,,q1][,t] <- (lab>max(ttdata$L))*max(ttdata$L)+(lab<min(ttdata$L))*min(ttdata$L)+(lab<=max(ttdata$L))*(lab>=min(ttdata$L))*lab
			lnlmed[,,q2][t,q1] <- median(lnldata[,,,q2][,,q1][,t])
		}
	}
}
#For the arrays, lnldata and lnmdata, the 1st dimension is time, 2nd is rank of innovation shock, 3rd is rank of input shock
#So "Low-Low" represents low labor demand and low lagged labor shock
labpath <- data.frame(1:T, lnlmed[,1,]-lnlmed[,3,], lnlmed[,2,]-lnlmed[,3,], lnlmed[,4,]-lnlmed[,3,], lnlmed[,5,]-lnlmed[,3,])
names(labpath) <- c("Time", "Q1Low", "Q1Med", "Q1High", "Q2Low", "Q2Med", "Q2High","Q3Low", "Q3Med", "Q3High","Q4Low", "Q4Med", "Q4High")
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
ltitles <- list("Low Labor Shock<br>High Labor", "Low Labor Shock<br>Medium Labor", "Low Labor Shock<br>High Labor", "Medium-Low Labor Shock<br>Low Labor", "Medium-Low Labor Shock<br>Medium Labor", "Medium-Low Labor Shock<br>High Labor", "Medium-High Labor Shock<br>Low Labor", "Medium-High Labor Shock<br>Medium Labor", "Medium-High Labor Shock<br>High Labor", "High Labor Shock<br>High Labor", "High Labor Shock<br>Medium Labor", "High Labor Shock<br>High Labor")
Lplotly <- list()
for (i in 1:12){
	ldat <- data.frame(Time=labpath$Time, Y=labpath[,i+1])
	Lplotly[[i]] <- plot_ly(ldat, x=~Time, y=~Y, type = 'scatter', mode = 'lines', showlegend=F, line=list(color="black"), name=ltitles[[i]], hovertemplate = paste("<i>Year<i>: %{x}", "<br>Labor Change: %{y:.2f}")) %>% layout(xaxis=list(title="Years", titlefont=list(size=18), tickfont=list(size=14)), yaxis=list(title="Labor", titlefont=list(size=18), tickfont=list(size=14), range=list(min(labpath[,-1]), max(labpath[,-1]))), shapes=list(hline(y=0)))


}
annotationsL <- list(list(x=0.11, y=0.97, text=TeX("\\boldsymbol{(a)\\,\\tau_{l_{t-1}}=0.1, \\tau_{l}=0.1}"), font=list(size=30), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
	list(x=0.495, y=0.97, text=TeX("\\boldsymbol{(b)\\,\\tau_{l_{t-1}}=0.1, \\tau_{l}=0.5}"), font=list(size=30), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.88, y=0.97, text=TeX("\\boldsymbol{(c)\\,\\tau_{l_{t-1}}=0.1, \\tau_{l}=0.9}"), font=list(size=30), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.11, y=0.7, text=TeX("\\boldsymbol{(d)\\,\\tau_{l_{t-1}}=0.25, \\tau_{l}=0.1}"), font=list(size=30), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.495, y=0.7, text=TeX("\\boldsymbol{(e)\\,\\tau_{l_{t-1}}=0.25, \\tau_{l}=0.5}"), font=list(size=30), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.88, y=0.7, text=TeX("\\boldsymbol{(f)\\,\\tau_{l_{t-1}}=0.25, \\tau_{l}=0.9}"), font=list(size=30), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.11, y=0.45, text=TeX("\\boldsymbol{(g)\\,\\tau_{l_{t-1}}=0.75, \\tau_{l}=0.1}"), font=list(size=30), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.495, y=0.45, text=TeX("\\boldsymbol{(h)\\,\\tau_{l_{t-1}}=0.75, \\tau_{l}=0.5}"), font=list(size=30), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.88, y=0.45, text=TeX("\\boldsymbol{(i)\\,\\tau_{l_{t-1}}=0.75, \\tau_{l}=0.9}"), font=list(size=30), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.11, y=0.2, text=TeX("\\boldsymbol{(j)\\,\\tau_{l_{t-1}}=0.9, \\tau_{l}=0.1}"), font=list(size=30), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.495, y=0.2, text=TeX("\\boldsymbol{(k)\\,\\tau_{l_{t-1}}=0.9, \\tau_{l}=0.5}"), font=list(size=30), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE),
		list(x=0.88, y=0.2, text=TeX("\\boldsymbol{(l)\\,\\tau_{l_{t-1}}=0.9, \\tau_{l}=0.9}"), font=list(size=30), xref="paper", yref="paper", xanchor="center,", yanchor="bottom", showarrow=FALSE))
L <- subplot(Lplotly[[1]], Lplotly[[2]], Lplotly[[3]], Lplotly[[4]], Lplotly[[5]], Lplotly[[6]], Lplotly[[7]], Lplotly[[8]], Lplotly[[9]], Lplotly[[10]], Lplotly[[11]], Lplotly[[12]], shareX=TRUE, shareY=TRUE, nrows=4)
Lplot <- L %>% layout(annotations=annotationsL) %>% config(mathjax = 'cdn')
Lplot
Ljson <- plotly_json(L, FALSE)
# write(Ljson, "/Users/justindoty/Documents/Home/My_Website/static/jmp/labor/impulseL.json")



