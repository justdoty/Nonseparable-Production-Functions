require(stringr) 
source('NLPFQR/DATA/RnD/YesRnD/Production_EM.R')
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
US <- read.csv('NLPFQR/DATA/USdata.csv') 
#How to split the data into performing vs non-performing firms?
#Option 1
#Create a binary variable for R&D status and filter observations in which firms switch from performing and not performing R&D following Doraszelski and Jaumandreu (2018)
# US <- US %>% mutate(rdB=ifelse(rd>0, 1, 0)) %>% group_by(id) %>% filter(rdB==lag(rdB)) %>% ungroup()
#Option 2
#Only select samples of firms that either perform or non-perform the entire sample period
US <- US %>% mutate(rdB=ifelse(rd>0, 1, 0)) %>% group_by(id) %>% filter(all(rdB==1)|all(rdB==0)) %>% ungroup() 
US <- data.frame(US)
print(summary(US))
print(nrow(US))
results <- Production_EM(ntau=ntau, idvar=US$id, timevar=US$year, Y=US$lny, K=US$lnk1, A=US$age, L=US$lnl, M=US$lnm, I=US$lni, R=US$rd, method=method, maxiter=maxiter, draws=draws)
#Save Results
ext <- paste("NLPFQR/DATA/RnD/YesRnD/YesRnD_", method, ".Rdata", sep="")
save(results, file=ext)

