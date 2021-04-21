require(stringr) 
# source('NLPFQR/FUN/Production_EM.R')
source('/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Functions/Production_EM.R')
#This function estimates the nonlinear quantile production function batching by industry
#Function calls Production_EM to perform estimation
#EM Algorithm Parameters
#Length of tau sequence (evenly divided)
ntau <- 11
#Total number of EM steps
maxiter <- 500
#Number of Metropolis-Hastings steps (burn-in)
draws <- 300
#Load US dataset
US_panel <- read.csv('/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Data/USdata.csv')
# US_panel <- read.csv('NLPFQR/DATA/US/USdata.csv')
#Convert 3 digit NAICS code to 2 digit NAICS and take natural logs
USdata <- transmute(US_panel, id=id, year=year, lny=log(Y), lnva=log(VA), lnk=log(K), lnl=log(L), 
	lnm=log(M), lni=log(I), naics3=naics3, naics2=as.numeric(str_extract(as.character(naics3), "^.{2}"))) %>% 
	group_by(id) %>% filter(n()>=3)
#Choose which industry to select
industries <- c("31", "32", "33", "^3")
id <- as.numeric(commandArgs(TRUE)[1])
industries <- industries[id]

# For testing purposes
maxiter <- 2
draws <- 5
industries <- c("32")


US <- filter(USdata, str_detect(naics2, industries))
results <- Production_EM(ntau=ntau, idvar=US$id, timevar=US$year, Y=US$lny, K=US$lnk, L=US$lnl, M=US$lnm, 
	I=US$lni, maxiter=maxiter, draws=draws)

# filename <- paste("NLPFQR/DATA/US/NLPFQR_NAICS_", id, ".RData", sep="")
# save(results, file=filename)

#HPC Job Submissions for batches: qsub -t 1:length(industries) myjob.job
#Here length(industries)=4




#Testing for Cubic Spline
# source('/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Functions/Auxfuns.R')
# require(pracma)
# ntau <- 11
# vectau <- seq(1/(ntau+1), ntau/(ntau+1), by=1/(ntau+1))
# vecint <- seq(vectau[1] , vectau[length(vectau)], length.out=20)
# par <- parY[1,]
# ytb1 <- 19.08
# ytbL <- 10.31

# d1 <- (vectau[1]*ytb1)^-1
# dN <- (ytbL*(1-vectau[length(vectau)]))^-1
# testcoef <- csplcoef(x=vectau, y=par, xi=vectau, endp2nd=TRUE, der=c(d1,dN), deriv=0)

# dmat <- cbind(c(d1, rep(0, (nrow(parY)-1))), c(dN, rep(0, (nrow(parY)-1))))

# try <- t(sapply(1:nrow(parY), function(i) csplcoef(x=vectau, y=parY[i,], xi=vectau, endp2nd=TRUE, der=dmat[i,], deriv=1)))

