require(stringr)
source('Production_EM.R')
#This function estimates the nonlinear quantile production function batching by industry
#Function calls Production_EM to perform estimation
#EM Algorithm Parameters
#Length of tau sequence (evenly divided)
ntau <- 11
#Total number of EM steps
maxiter <- 500
#Number of Metropolis-Hastings steps (burn-in)
draws <- 200
#Number of draws of Metropolis Hastings kept for computation (Mdraws=1 for stEM for now)
Mdraws <- 1
#Degrees of the Univariate Hermite Polynomials
#Output
MY <- as.numeric(c(MYK=1, MYL=2, MYM=1, MYW=2))
#Labor
ML <- as.numeric(c(MLK=1, MLW=2))
#Materials
MM <- as.numeric(c(MMK=1, MMW=2))
#Productivity t>1
MW <- as.numeric(c(MWW=2))
#Investment
MI <- as.numeric(c(MIK=1, MIW=2))
#Combined
MH <- list(MY=MY, ML=ML, MM=MM, MW=MW, MI=MI)
seed <- 123456
#Load US dataset
US_panel <- read.csv("USdata.csv")
#Convert 3 digit NAICS code to 2 digit NAICS and take natural logs
USdata <- transmute(US_panel, id=id, year=year, lny=log(Y), lnva=log(VA), lnk=log(K), lnl=log(L), lnm=log(M), lni=log(I), naics3=naics3, naics2=as.numeric(str_extract(as.character(naics3), "^.{2}")))
#Choose which industry to select
All <- "^3"
industries <- c("31", "32", "33", All)
id <- as.numeric(commandArgs(TRUE)[1])
industries <- industries[id]


US <- filter(USdata, str_detect(naics2, industries))
results <- Proudction_EM(ntau=ntau, idvar=US$id, timevar=US$year, Y=US$lny, K=US$lnk, L=US$lnl, M=US$lnm, 
	I=US$lni, maxiter=maxiter, draws=draws, Mdraws=Mdraws, seed=seed, MH=MH)

filename <- paste("NLPFQR_NAICS_", id, ".RData", sep="")
save(results, file=filename)

#HPC Job Submissions for batches: qsub -t 1:length(industries) myjob.job
#Here length(industries)=3






