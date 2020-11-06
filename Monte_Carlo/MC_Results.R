load('/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Monte_Carlo/MCresults.Rdata')
source('/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Functions/Auxfuns.R')
source('/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Functions/Tensors.R')
require(dplyr)
require(purrr)
#Assign names to Monte Carlo environment
MH <- MH; nreps <- nreps; ntau <- ntau; MCdata <- MCdata; MCres <- MCres
#Grid of taus
vectau <- seq(1/(ntau+1), ntau/(ntau+1), by=1/(ntau+1))
#Number of firms
N <- length(unique(MCdata$id))
#Degrees of the Univariate Hermite Polynomials
#Output
MY <- as.numeric(c(MYK=1, MYL=2, MYM=1, MYA=1, MYW=2))
#Labor
ML <- as.numeric(c(MLK=1, MLA=1, MLW=2))
#Materials
MM <- as.numeric(c(MMK=1, MMA=1, MMW=2))
#Productivity t>1
MW <- as.numeric(c(MWA=1, MWW=2))
#Productivity t=1
MW1 <- as.numeric(c(MW1A=1))
#Investment
MI <- as.numeric(c(MIK=2, MIW=2, MIA=1))
#Combined
MH <- list(MY=MY, ML=ML, MM=MM, MW=MW, MW1=MW1, MI=MI)
###########################################################################################
#Simulate productivity using model estmates
###########################################################################################
seed <- 123456
set.seed(123456)
#############################################
#Initital Productivity
#############################################
#Initial Age
age1 <- t0data(idvar=MCdata$id, X=MCdata$A)[,2]
#Tensor product polynomial for initial age
wxt1 <- tensor.prod(MH$MW1, age1, norm=0)
#Draw of initial productivity innovation shocks
xi1 <- runif(N)
#Initial Productivity
wt1 <- rowSums(wxt1*lspline(vectau=vectau, bvec=matrix(MCres$resW0sim, prod(MW1+1), ntau), b1=MCres$resw0bsim[1], bL=MCres$resw0bsim[2], u=xi1))
#############################################
#Merge productivity data with MCdata
#############################################
set.seed(seed)
MCdata <- MCdata %>% group_by(time) %>% mutate(omega=wt1, xi=runif(length(id))) %>% ungroup() %>% group_by(id) %>%
	mutate(omega=c(omega[1], sapply(2:length(time), function(t) wtsim(wtlag=omega[t-1], X=A[t], xi=xi[t], dimM=MW,
			vectau=vectau, bvec=matrix(MCres$resWTsim, prod(MW+1), ntau), b1=MCres$reswtbsim[1], bL=MCres$reswtbsim[2]))))
test <- data.frame(MCdata %>% select(c(id, time, Y, omega, xi)))



