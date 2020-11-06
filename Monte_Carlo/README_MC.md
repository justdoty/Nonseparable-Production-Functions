#This folder contains the code necessary to replicate the Monte Carlo Experiment

#File Monte_Carlo.R produces the data generating processes. Uses functions PEM_SIM.R and Posterior_MC.R
#Saves output as an R environment.

#File PEM_SIM.R is the file that estimates the model. Note that this is different from the function 
#Production_EM.R in the Functions folder which generalizes to tensor product polynomials.

#Posterior_MC.R is the file that computes the posterior distribution. This is also different from Posterior.R 
#in the Functions folder which generalizes to tensor prodcut polynomial. Here, Labor, Materials, Investment,
#and initial productivity are parametrically specified.

#The file MC_Results.R takes the environment from Monte_Carlo.R and loads it. Computes all final estimates and graphs used in the paper