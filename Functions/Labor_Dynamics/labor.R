require(stringr) 
source('NLPFQR/DATA/Labor/Production_EM.R')
# source('/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Functions/Production_EM.R')
#Function calls Production_EM to perform estimation
#EM Algorithm Parameters
#Length of tau sequence (evenly divided)
ntau <- 11
#Total number of EM steps
maxiter <- 200
#Number of Metropolis-Hastings steps (burn-in)
draws <- 500
#Specification of the Production Function
# method <- "cobbN"
# method <- "transN"
# method <- "cobb"
method <- "trans"
# method <- "orth"
# method <- "hermite"
#Load US dataset
# US_panel <- read.csv('/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Data/USdata.csv')
US <- read.csv('NLPFQR/DATA/USdata.csv') 
US <- data.frame(US)
print(head(US))
print(summary(US))
results <- Production_EM(ntau=ntau, idvar=US$id, timevar=US$year, Y=US$lny, K=US$lnk1, A=US$age, L=US$lnl, M=US$lnm, I=US$lni, RD=US$rd, method=method, maxiter=maxiter, draws=draws)
#Save Results
ext <- paste("NLPFQR/DATA/Labor/labor_", method, ".Rdata", sep="")
save(results, file=ext)

