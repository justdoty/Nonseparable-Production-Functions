# A Dynamic Framework for Identification and Estimation of Nonseparable Production Functions
by Justin Doty (2021)
## Abstract
This paper studies identification and estimation of a nonseparable model for production functions with unobserved heterogeneity. Non-parametric identification results are established for the production function and productivity process under stationarity conditions. This framework allows for heterogeneous effects of output elasticities and factor efficiencies in addition to non-linear productivity persistence. It also allows for additional unobservables in the input demand functions, which would violate the scalar unobservability requirement in proxy variables under previous approaches. This extension is used to show firm's heterogeneous responses to productivity shocks corresponding to the size of their input demand. This paper illustrates these results in an application to US manufacturing firms where the proposed model is estimated using non-linear quantile regression. 

## Software Implementation
All code is written in R language and evaluated on an HPC system.

## Replication Files
The following folders contain replication files for the tables and figures in the paper.
1. [Environments](/Environments): A folder containing all R environments used in the analysis
2. [Figures](/Figures): A folder containing all the figures produced by the replication files
3. [Functions](/Functions): A folder containing all of the replication files which contains the following
	- [Labor Dynamics](/Functions/Labor_Dynamics): Contains replication files for producing impulse response functions for labor in the Appendix
	- [Main](/Functions/Main): Contains replication files for producing the results in the main paper
	- [RnD](/Functions/RnD): Contains replication files for producing R&D comparison estimates in the Appendix 
	- [Selection_Bias](/Functions/Selection_Bias): Reproduces the estimates in the main paper correcting for selection bias in the Appendix
	- [Data Cleaning](/Functions/Compustat_Cleaning.R): Procedure for cleaning the data from Compustat
	
In addition, each of the folders contain variations of the same files used in estimation
- Production_EM.R: Main estimation file. Uses functions below as dependencies.
- Posterior.R: Function used to compute to posterior density
- Mstep.R: Function to perform Maximization step of algorithm
- Tensors.R: Function used to specify functional forms for the production function, inputs, and productivity
- omega.R: Function used to estimate a production function using [Levinsohn and Petrin (2003)](https://doi.org/10.1111/1467-937X.00246)
- Auxfuns.R: Contains auxiliary files used in estimation procedures
- elasticities.R: Function used to construct production function estimates
- impulse_response.R: Function used to construct productivity and input responses to shocks


