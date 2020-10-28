source('/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Functions/Production_EM.R')
#For HPC 
# source('Production_EM.R')
#EM Algorithm Parameters
#Length of tau sequence (evenly divided)
ntau <- 11
#Total number of EM steps
maxiter <- 2
#Number of Metropolis-Hastings steps (burn-in)
draws <- 50
#Number of draws of Metropolis Hastings kept for computation (Mdraws=1 for stEM for now)
Mdraws <- 1
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
#Monte Carlo Simulation Parameters
nreps <- 2
seed <- 123456
#Monte Carlo Parameters
#Storage Matrices for Results
resYsim <- array(0, c(prod(MY+1), ntau, nreps))
resybsim <- array(0, c(2, nreps))
resLsim <- array(0, c(prod(ML+1), ntau, nreps))
reslbsim <- array(0, c(2, nreps))
resMsim <- array(0, c(prod(MM+1), ntau, nreps))
resmbsim <- array(0, c(2, nreps))
resWTsim <- array(0, c(prod(MW+1), ntau, nreps))
reswtbsim <- array(0, c(2, nreps))
resW0sim <- array(0, c(prod(MW1+1), ntau, nreps))
resw0bsim <- array(0, c(2, nreps))
resIsim <- array(0, c((prod(MI+1)+1), nreps))
#Parameters for data evolution
#Standard deviation of log wage process
siglnw <- 0.1
#Standard deviation of log material price process
siglnm <- 0.1
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
overallt <- 100
starttime <- 90
t <- overallt - starttime
#Production Function Parameters
beta0 <- 0
betal <- 0.3
betak <- 0.4
betam <- 0.3
betaa <- .01
betaw <- 1
#Epsilons and omega process
sigeps <- 0.1
sigomg <- 0.3 #standard deviation of omega
rho <- 0.7 #AR(1) coefficient for omega
sigxi <- sqrt((1-rho^2)*sigomg^2)
#Matrices to store data
lnkdata <- matrix(0, n, overallt) #ln(capital)
lnldata <- matrix(0, n, overallt) #ln(labor)
lnmdata <- matrix(0, n, overallt) #ln(intermediate input)
lnydata <- matrix(0, n, overallt) #ln(output)
omgdata <- matrix(0, n, overallt) #omega(t)
agedata <- matrix(0, n, overallt) #age(t)
#Location Scale Parameters for Output
eta0 <- 1
etak <- 0.3
etal <- -0.4
etam <- 0.4
#Location Scale Parameters for Labor
hetl0 <- 1
hetlk <- -0.2
hetlw <- 0.3
#Location Scale Parameters for Materials
hetm0 <- 1
hetmk <- -0.2
hetmw <- 0.3
for (niter in 1:nreps){
  print(sprintf("MC Iteration %i", niter))
  seed <- seed+niter
  set.seed(seed)
  #Specification for Production Shock Distribution
  etadata <- matrix(rnorm(n*overallt, 0, sigeps), nrow=n, ncol=overallt)
  #Specification for Labor Input Shock Distribution
  varepsdata <- matrix(rnorm(n*overallt, 0, sigoptl), nrow=n, ncol=overallt)
  #Specification for Materials Input Shock Distribution
  epsdata <- matrix(rnorm(n*overallt, 0, sigoptm), nrow=n, ncol=overallt)
  #Specification for Investment Input Shock Distribution
  iotadata <- matrix(rnorm(n*overallt, 0, sigopti), nrow=n, ncol=overallt)
  #Specification for labor cost (not serially correlated)
  lnwdata <- matrix(rnorm(n*overallt, 0, siglnw), nrow=n, ncol=overallt)
  #Specification for materials cost (not serially correlated)
  lnpmdata <- matrix(rnorm(n*overallt, 0, siglnm), nrow=n, ncol=overallt)
  #Period 0 values of omega 
  omgdata0 <- matrix(rnorm(n,0,sigomg),nrow=n,ncol=1)
  #Period 0 values of age
  agedata0 <- floor(as.matrix(runif(n, 1, 10)))
  agedata[,1] <- agedata0+1
  #Period 1 values of omega 
  omgdata[,1] <- betaa*agedata0[,1]+rho*omgdata0[,1]+matrix(rnorm(n,0,sigomg),nrow=n,ncol=1)
  #Simulate values of omega for rest of time periods
  for (s in 2:overallt){
    agedata[,s] <- agedata[,s-1]+1
    omgdata[,s] <- betaa*agedata[,s]+rho*omgdata[,s-1] + matrix(rnorm(n,0,sigomg),nrow=n,ncol=1)

  }
  #Intital ln(capital) level (close to 0)
  lnkdata[,1] <- matrix(-100,n,1)
  #Discount Rate for DP problem
  disc <- 0.95
  #Depreciation rate of capital
  delta <- 0.2
  #Capital Adjustment Cost Parameter
  badj <- 1
  #Simplifying components of the optimal investment rule
  #Square bracket component
  squarebracketterm <- exp(0.5*betal^2*sigoptl^2)*exp(0.5*betam^2*sigoptm^2)
  #Constant term in front of sum including squarebracketterm
  const1 <- disc*betak*betal^(betal/betak)*betam^(betam/betak)*(exp(beta0)^(1/(betak)))*
  squarebracketterm/badj
  vec1 <- (disc*(1-delta))^seq(100)
  vec2 <- cbind(sigxi^2 * 0,cumsum(rho^(2*(seq(100)-1))))
  expterm1 <- exp(0.5*(betaw/(betak))^2*rho^2*
  ((sigxi^2)*rho^(2*seq(100))+vec2))
  #Compute Optimal Investment and Capital stock for all firms over time
  investmat <- matrix(NA, n, overallt)
  for (i in 1:n){
      for (s in 1:overallt){
        expterm2 <- exp((betaw/betak)*omgdata[i,s]*rho^(seq(100)))
        expterm3 <- exp(0.5*betal^2*((1+hetlk*lnkdata[i,s]+hetlw*omgdata[i,s])*sigoptl)^2)
        expterm4 <- exp(0.5*betam^2*((1+hetmk*lnkdata[i,s]+hetmw*omgdata[i,s])*sigoptm)^2)
        investmat[i,s] <- const1*sum(vec1*expterm1*expterm2*expterm3*expterm4)
        if (s >= 2){
          lnkdata[i,s] <- log((1-delta)*exp(lnkdata[i,s-1])+investmat[i,s-1])+rnorm(1, mean=0, sd=sigoptk)
        }
      }
    }

  #Generate levels of labor input
  lnldata <- (beta0+(1-betam)*log(betal)+betak*lnkdata+betam*log(betam)-betam*lnpmdata-(1-betam)*lnwdata+betaw*omgdata)/(1-betam-betal)
  #Adding Optimization Error
  hetl <- hetl0+hetlk*lnkdata+hetlw*omgdata
  lnldata <- lnldata + hetl*matrix(rnorm(n*overallt,0,sigoptl),n,overallt)
  #Generate levels of material input
  lnmdata <- (beta0+log(betam)+betak*lnkdata+betal*(log(betal)-lnwdata-log(betam))+betaw*omgdata-(1-betal)*lnpmdata)/(1-betal-betam)
  hetm <- hetm0+hetmk*lnkdata+hetmw*omgdata
  lnmdata <- lnmdata + hetl*matrix(rnorm(n*overallt,0,sigoptl),n,overallt)
  #Add the random shock to investmnet
  investmat <- log(investmat)+iotadata
  #Output
  #Specifies the form of heteroskedasticity
  hety <- eta0+etal*lnldata+etak*lnkdata
  lnydata <- beta0 + betal*lnldata + betak*lnkdata + betam*lnmdata + omgdata + hety*etadata

  #Stack data across firms (all the data)
  Capital <- c(t(lnkdata[,(starttime+1):overallt]))
  Labor <- c(t(lnldata[,(starttime+1):overallt]))
  Materials <- c(t(lnmdata[,(starttime+1):overallt]))
  Investment <- c(t(investmat[,(starttime+1):overallt]))
  Wage <- c(t(lnwdata[,(starttime+1):overallt]))
  Output <- c(t(lnydata[,(starttime+1):overallt]))
  Productivity <- c(t(omgdata[,(starttime+1):overallt]))
  Epsilon <- c(t(etadata[,(starttime+1):overallt]))

  #Stack data across firms (lagged data)
  Capital_Lag_1 <- c(t(lnkdata[,(starttime+1):(overallt-1)]))
  Labor_Lag_1 <- c(t(lnldata[,(starttime+1):(overallt-1)]))
  Labor_Lag_2 <- c(t(lnldata[,(starttime):(overallt-2)]))
  Materials_Lag_1 <- c(t(lnmdata[,(starttime+1):(overallt-1)]))
  Investment_Lag_1 <- c(t(investmat[,(starttime+1):overallt-1]))
  Wage_Lag_1 <- c(t(lnwdata[,(starttime+1):(overallt-1)]))
  Output_Lag_1 <- c(t(lnydata[,(starttime+1):(overallt-1)]))
  Productivity_Lag_1 <- c(t(omgdata[,(starttime+1):(overallt-1)]))

  #Stack data across firms (contemporaneous data)
  Capital_Con <- c(t(lnkdata[,(starttime+2):(overallt)]))
  Labor_Con <- c(t(lnldata[,(starttime+2):(overallt)]))
  Materials_Con <- c(t(lnmdata[,(starttime+2):(overallt)]))
  Investment_Con <- c(t(investmat[,(starttime+2):overallt]))
  Wage_Con <- c(t(lnwdata[,(starttime+2):(overallt)]))
  Output_Con <- c(t(lnydata[,(starttime+2):(overallt)]))
  Productivity_Con <- c(t(omgdata[,(starttime+2):(overallt)]))
  Age_Con <- c(t(agedata[,(starttime+2):(overallt)]))
  #Create ID and Time variables
  idvar <- rep(1:n, each=(overallt-(starttime+1)))
  timevar <- rep(1:(overallt-(starttime+1)), n)
  if (niter==nreps){
    MCdata <- data.frame(id=idvar, time=timevar, Y=Output_Con, K=Capital_Con, L=Labor_Con, M=Materials_Con, I=Investment_Con, A=Age_Con)
  }
  #Start Timer
  overall.start.time <- Sys.time()
  #Results
  results <- Production_EM(ntau=ntau, idvar=idvar, timevar=timevar, Y=Output_Con, K=Capital_Con, L=Labor_Con, M=Materials_Con, I=Investment_Con, A=Age_Con, maxiter=maxiter, draws=draws, Mdraws=Mdraws, seed=seed, MH)
  resYsim[,,niter] <- results$resYmat; resybsim[,niter] <- results$resyb1bLmat
  resLsim[,,niter] <- results$resLmat; reslbsim[,niter] <- results$reslb1bLmat
  resMsim[,,niter] <- results$resMmat; resmbsim[,niter] <- results$resmb1bLmat
  resWTsim[,,niter] <- results$resWTmat; reswtbsim[,niter] <- results$reswtb1bLmat
  resW0sim[,,niter] <- results$resW0mat; resw0bsim[,niter] <- results$resw0b1bLmat
  resIsim[,niter] <- results$resImat
  #Loop time
  print(Sys.time()-overall.start.time)
}
MCres <- list(resYsim=resYsim, resybsim=resybsim, resLsim=resLsim, reslbsim=reslbsim, resMsim=resMsim,
  resmbsim=resmbsim, resWTsim=resWTsim, reswtbsim=reswtbsim, resW0sim=resW0sim, resw0bsim=resw0bsim, resIsim=resIsim)
#Save Results
save(MH, nreps, ntau, MCres, MCdata, file='/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Code/Monte_Carlo/MCresults.Rdata')
#Total time of MC experiment
print(Sys.time()-overall.start.time)


