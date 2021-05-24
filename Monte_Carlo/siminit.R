require(stringr) 
source('NLPFQR/FUN/Production_EM.R')
# source('/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Functions/Production_EM.R')
#Function calls Production_EM to perform estimation
#EM Algorithm Parameters
#Length of tau sequence (evenly divided)
ntau <- 11
#Total number of EM steps
maxiter <- 500
#Number of Metropolis-Hastings steps (burn-in)
draws <- 200
#Load US dataset
# US_panel <- read.csv('/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Data/USdata.csv')
US_panel <- read.csv('NLPFQR/DATA/US/USdata.csv')
#Convert 3 digit NAICS code to 2 digit NAICS and take natural logs
US <- transmute(US_panel, id=id, year=year, lny=log(Y), lnva=log(VA), lnk=log(K), lnl=log(L), 
	lnm=log(M), lni=log(I), naics3=naics3, naics2=as.numeric(str_extract(as.character(naics3), "^.{2}")))
results <- Production_EM(ntau=ntau, idvar=US$id, timevar=US$year, Y=US$lny, K=US$lnk, L=US$lnl, M=US$lnm, 
	I=US$lni, maxiter=maxiter, draws=draws)
#Save Results
save(results, file="NLPFQR/SIM/siminit.Rdata")





