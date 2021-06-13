require(stringr) 
source('NLPFQR/DATA/Parametric/Production_EM.R')
# source('/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Functions/Production_EM.R')
#Function calls Production_EM to perform estimation
#EM Algorithm Parameters
#Length of tau sequence (evenly divided)
ntau <- 11
#Total number of EM steps
maxiter <- 200
#Number of Metropolis-Hastings steps (burn-in)
draws <- 200
#Specification of the Production Function
# method <- "cobbN"
# method <- "transN"
# method <- "cobb"
method <- "trans"
# method <- "hermite"
#Load US dataset
# US_panel <- read.csv('/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Data/USdata.csv')
US_panel <- read.csv('NLPFQR/DATA/USdata.csv')
#Convert 3 digit NAICS code to 2 digit NAICS and take natural logs
US <- transmute(US_panel, id=id, year=year, lny=log(Y), lnva=log(VA), lnk=log(K2), lnl=log(L), 
	lnm=log(M), lni=log(I), age=age)
results <- Production_EM(ntau=ntau, idvar=US$id, timevar=US$year, Y=US$lny, K=US$lnk, A=US$age, L=US$lnl, M=US$lnm, I=US$lni, method=method, maxiter=maxiter, draws=draws)
#Save Results
ext <- paste("NLPFQR/DATA/Parametric/US_Parametric_", method, ".Rdata", sep="")
save(results, file=ext)

