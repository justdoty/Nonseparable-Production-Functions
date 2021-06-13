require(stringr) 
source('NLPFQR/DATA/Parametric/Production_EM.R')
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
US <- read.csv('NLPFQR/DATA/USdata.csv') %>% transmute(id=id, year=year, lny=log(Y), lnva=log(VA), lnk=log(K2), lnl=log(L), lnm=log(M), lni=log(I), age=age)
#How many iterations to do in one bath
Rbatch <- 20
#How many batches to do
n <- 10
#List of results
Rlist <- list()
###############
id <- as.numeric(commandArgs(TRUE)[1])
Rseed <- id+123456
set.seed(Rseed)
for (r in 1:Rbatch){
	print(r)
	Rind <- block.boot.resample(idvar=as.matrix(US$id), R=Rbatch)
	US <- US[Rind[[r]],]
	results <- Production_EM(ntau=ntau, idvar=US$id, timevar=US$year, Y=US$lny, K=US$lnk, A=US$age, L=US$lnl, M=US$lnm, I=US$lni, method=method, maxiter=maxiter, draws=draws)
	print(results)
	Rlist[[r]] <- results
}
#Save Results
ext <- paste("NLPFQR/DATA/Parametric/Bootstrap/US_PAR_Boot_", method, id, ".Rdata", sep="")
save(Rlist, file=ext)

