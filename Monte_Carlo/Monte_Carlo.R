source('/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Functions/Production_EM.R')
#For HPC 
# source('Production_EM.R')
#EM Algorithm Parameters
#Length of tau sequence (evenly divided)
ntau <- 11
vectau <- seq(1/(ntau+1), ntau/(ntau+1), by=1/(ntau+1))
#Total number of EM steps
maxiter <- 2
#Number of Metropolis-Hastings steps (burn-in)
draws <- 50
#Number of draws of Metropolis Hastings kept for computation (Mdraws=1 for stEM for now)
Mdraws <- 1
#Degrees of the Univariate Hermite Polynomials
#Output
MY <- "linear"
#Labor
ML <- as.numeric(c(MLK=2, MLW=2))
#Materials
MM <- as.numeric(c(MMK=2, MMW=2))
#Productivity t>1
MW <- "linear"
#Investment
MI <- as.numeric(c(MIK=2, MIW=2))
MH <- list(MY=MY, ML=ML, MM=MM, MW=MW, MI=MI)
#Monte Carlo Simulation Parameters
nreps <- 2
seed <- 123456
#Monte Carlo Parameters
#Storage Matrices for Results
resYsim <- array(0, c(5, ntau, nreps))
resybsim <- array(0, c(2, nreps))
resLsim <- array(0, c(prod(ML+1), ntau, nreps))
reslbsim <- array(0, c(2, nreps))
resMsim <- array(0, c(prod(MM+1), ntau, nreps))
resmbsim <- array(0, c(2, nreps))
resWTsim <- array(0, c(2, ntau, nreps))
reswtbsim <- array(0, c(2, nreps))
resIsim <- array(0, c((prod(MI+1)+1), nreps))
#Parameters for data evolution
##########################################
#Production:
#########################################
#Constant
beta0 <- 0
#Tail parameters (Best way to choose these?)
betaexp <- c(10,10)
#Capital
klb <- 1/11
kub <- 11/44
betak <- seq(klb, kub, length.out=ntau)
#Labor
llb <- 6/11
lub <- 20/44
betal <- seq(llb, lub, length.out=ntau)
#Materials
mlb <- 4/11
mub <- 13/44
betam <- seq(mlb, mub, length.out=ntau)
beta <- rbind(beta0, betak, betal, betam)
##########################################
#Productivity:
#########################################
#Constants
#Lower Bound and Upper Bound
#Tail Parameters
rhoexp <- c(10, 10)
#Persistency Parameters
#Lower Bound and Upper Bound
rhob <- 0.5; rhou <- 0.9
#Combined
rho <- rbind(0, seq(rhob, rhou, length.out=ntau))
##########################################
#Investment:
#########################################
iota0 <- 0
iotak <- -0.02
iotaw <- .5
#########################################
#Standard deviation of optimization error for labor
sigoptl <- 0.37
#Standard deviation of optimization error for materials
sigoptm <- 0.37
#Standard deviation of random shock to capital accumulation
sigoptk <- 0.2
#Standard deviation of random shock to investment
sigopti <- 0.37
# Number of Firms
n <- 1000
# Number of Time Periods
overallt <- 10
starttime <- 8
#Epsilons and omega process
sigomg <- 0.05 #standard deviation of omega
#Matrices to store data
lnydata <- matrix(0, n, overallt) #ln(output)
omgdata <- matrix(0, n, overallt) #omega(t)
kdata <- matrix(0, n, overallt) #ln(capital)
investmat <- matrix(0, n, overallt) #Investment
lnldata <- matrix(0, n, overallt) #ln(labor)
lnmdata <- matrix(0, n, overallt) #ln(intermediate input)
############################################################
for (niter in 1:nreps){
  print(sprintf("MC Iteration %i", niter))
  seed <- seed+niter
  set.seed(seed)
  ##################################################################
  #Specification for Production Shock Distribution
  etadata <- matrix(runif(n*overallt, 0, 1), nrow=n, ncol=overallt)
  #Specification for Productivity Innovation Shock Distribution
  xidata <- matrix(runif(n*overallt, 0, 1), nrow=n, ncol=overallt)
  #Specification for Labor Input Shock Distribution
  varepsdata <- matrix(rnorm(n*overallt, 0, sigoptl), nrow=n, ncol=overallt)
  #Specification for Materials Input Shock Distribution
  epsdata <- matrix(rnorm(n*overallt, 0, sigoptm), nrow=n, ncol=overallt)
  #Specification for Investment Input Shock Distribution
  iotadata <- matrix(rnorm(n*overallt, 0, sigopti), nrow=n, ncol=overallt)
  #Specification for labor cost (not serially correlated)
  pl <- 0.3
  #Specification for materials cost (not serially correlated)
  pm <- 0.2
  ##########################################################################

  #####################################################################
  #Period 1 values of omega 
  omgdata[,1] <- matrix(rnorm(n,0,sigomg),nrow=n,ncol=1)
  #Simulate values of omega for rest of time periods
  for (s in 2:overallt){
    omgdata[,s] <- rowSums(omgdata[,s-1]*lspline(vectau=vectau, bvec=rho, b1=rhoexp[1], bL=rhoexp[2], u=xidata[,s]))

  }
  ########################################################################
  #########################################
  #Intital ln(capital) level (close to 0)
  kdata[,1] <- matrix(1,n,1)
  #Depreciation rate of capital
  delta <- 0.2
  ########################################

  ######################################################################
  #Compute Optimal Investment and Capital stock for all firms over time
  for (i in 1:n){
      for (s in 1:overallt){
        investmat[i,s] <- iota0+iotak*log(kdata[i,s])+omgdata[i,s]+iotadata[i,s]
        if (s >= 2){
          kdata[i,s] <- (1-delta)*kdata[i,s-1]+exp(investmat[i,s-1])+rlnorm(1, meanlog=0, sdlog=sigoptk)
        }
      }
    }
  lnkdata <- log(kdata)
  ############################################################################

  ###########################################################################################
  #Generate the interpolating spline for output coefficients
  betaspl <- lspline(vectau=vectau, bvec=beta, b1=betaexp[1], bL=betaexp[2], u=c(t(etadata)))
  #Generate the mean of the quantile output coefficients
  betabar <- colMeans(betaspl)
  ###########,f################################################################################
  #################################
  #Generate levels of labor input
  #################################
  betal0 <- (1/(1-betabar["betal"]-betabar["betam"]))*(log(betabar["betal"])-log(pl)+betabar["betam"]*log((pl*betabar["betam"])/pm*betabar["betal"]))
  betalk <- betabar["betak"]/(1-betabar["betal"]-betabar["betam"])
  betalw <- 1/(1-betabar["betal"]-betabar["betam"])
  lnldata <- betal0+betalk*lnkdata+betalw*omgdata+varepsdata
  #################################
  #Generate levels of material input
  #################################
  betam0 <- (1/(1-betabar["betal"]-betabar["betam"]))*(log(betabar["betam"])-log(pm)+betabar["betal"]*log((pm*betabar["betal"])/pl*betabar["betam"]))
  betamk <- betabar["betak"]/(1-betabar["betal"]-betabar["betam"])
  betamw <- 1/(1-betabar["betal"]-betabar["betam"])
  lnmdata <- betam0+betamk*lnkdata+betamw*omgdata+epsdata
  ##############################################################
  #Output
  ###############################################################
  lnydata <- rowSums(betaspl*cbind(1, c(t(lnkdata)), c(t(lnldata)), c(t(lnmdata))))

  #Stack data across firms 
  Output <- matrix(lnydata, nrow=n, ncol=overallt, byrow=TRUE)
  Capital <- matrix(lnkdata, nrow=n, ncol=overallt, byrow=TRUE)
  Labor <- matrix(lnldata, nrow=n, ncol=overallt, byrow=TRUE)
  Materials <- matrix(lnmdata, nrow=n, ncol=overallt, byrow=TRUE)
  Productivity <- matrix(omgdata, nrow=n, ncol=overallt, byrow=TRUE)
  Investment <- matrix(investmat, nrow=n, ncol=overallt, byrow=TRUE)
  #Take the last 3 time periods for estimation and stack across firms
  Output_Con <- c(t(Output[,starttime:overallt]))
  Capital_Con <- c(t(Capital[,starttime:overallt]))
  Labor_Con <- c(t(Labor[,starttime:overallt]))
  Materials_Con <- c(t(Materials[,starttime:overallt]))
  Productivity_Con <- c(t(Productivity[,starttime:overallt]))
  Investment_Con <- c(t(Investment[,starttime:overallt]))
  #Create ID and Time variables
  idvar <- rep(1:n, each=(overallt-starttime+1))
  timevar <- rep(1:(overallt-starttime+1), n)
  if (niter==nreps){
    MCdata <- data.frame(id=idvar, time=timevar, Y=Output_Con, K=Capital_Con, L=Labor_Con, M=Materials_Con, I=Investment_Con)
  }
  #Start Timer
  overall.start.time <- Sys.time()
  #Results
  # results <- Production_EM(ntau=ntau, idvar=idvar, timevar=timevar, Y=Output_Con, K=Capital_Con, L=Labor_Con, M=Materials_Con, I=Investment_Con, maxiter=maxiter, draws=draws, Mdraws=Mdraws, seed=seed, MH)
  # resYsim[,,niter] <- results$resYmat; resybsim[,niter] <- results$resyb1bLmat
  # resLsim[,,niter], <- results$resLmat; reslbsim[,niter] <- results$reslb1bLmat
  # resMsim[,,niter] <- results$resMmat; resmbsim[,niter] <- results$resmb1bLmat
  # resWTsim[,,niter] <- results$resWTmat; reswtbsim[,niter] <- results$reswtb1bLmat
  # resW1sim[,,niter] <- results$resW1mat; resIsim[,niter] <- results$resImat
  #Loop time
  print(Sys.time()-overall.start.time)
}
# MCres <- list(resYsim=resYsim, resybsim=resybsim, resLsim=resLsim, reslbsim=reslbsim, resMsim=resMsim,
#   resmbsim=resmbsim, resWTsim=resWTsim, reswtbsim=reswtbsim, resW1sim=resW1sim, resIsim=resIsim)
#Save Results
# save(MH, nreps, ntau, MCres, MCdata, file='/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Monte_Carlo/MCresults.Rdata')
#Total time of MC experiment
print(Sys.time()-overall.start.time)






