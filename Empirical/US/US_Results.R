require(cowplot)
load('/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Empirical/US/Environments/NLPFQR_NAICS_3.Rdata')
USdata <- read.csv('/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Data/USdata.csv', header=TRUE) %>% 
	mutate(naics2=as.numeric(str_extract(as.character(naics3), "^.{2}"))) %>% group_by(id) %>% filter(n()>=3,str_detect(naics2, "33"))
USdata <- data.frame(USdata)
source('/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Functions/Auxfuns.R')
source('/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Functions/Tensors.R')
require(dplyr)
require(purrr)
ntau <- 11
#Assign names to Monte Carlo environment
MH <- MH; nreps <- nreps; ntau <- ntau; USdata <- USdata; MCres <- MCres
#Grid of taus
vectau <- seq(1/(ntau+1), ntau/(ntau+1), by=1/(ntau+1))
#Number of firms
Nid <- unique(USdata$id)
N <- length(Nid)
#Degrees of the Univariate Hermite Polynomials
#Output
#Degrees of the Univariate Hermite Polynomials
#Output
MY <- as.numeric(c(MYK=1, MYL=1, MYM=1, MYW=1))
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
###########################################################################################
#Simulate productivity using model estmates
###########################################################################################
seed <- 123456
set.seed(123456)
#############################################
#Initital Productivity
#############################################
#Initial Productivity
wt1 <- data.frame(idvar=Nid, wt1=rnorm(N, mean=results$resW1mat[1], sd=sqrt(results$resW1mat[2])))
#############################################
#Merge productivity data with USdata
#############################################
USdata <- USdata %>% group_by(id) %>% mutate(omega=rep(wt1$wt1[which(wt1$idvar==id)], each=length(year))) %>% 
	group_by(year) %>% mutate(xi=runif(length(id))) %>% ungroup() %>% group_by(id) %>%
	mutate(omega=c(omega[1], sapply(2:length(year), function(t) wtsim(wtlag=omega[t-1], xi=xi[t], dimM=MW,
			vectau=vectau, bvec=matrix(results$resWTmat, prod(MW+1), ntau), b1=results$reswtb1bLmat[1], bL=results$reswtb1bLmat[2]))))
omega <- sort(USdata$omega)[200:(length(USdata$omega)-200)]
omegadens <- plot(density(omega))
save('/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Empirical/US/wplot.png', omegadens)

per <- results$resWTmat[2,]+2*mean(omega)*results$resWTmat[3,]
plot(vectau, per, type="l", xlab="tau", ylab="Average Persitence")