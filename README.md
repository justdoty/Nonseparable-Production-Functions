# A Dynamic Framework for Identification and Estimation of Nonlinear Production Functions
by Justin Doty (2021)
## Abstract
This paper studies identification and estimation of a non-linear model for production functions with unobserved heterogeneity. Non-parametric identification results are established for the production function and productivity process under stationarity assumptions. This paper then estimates the conditional quantiles of firm production using non-linear quantile regression. This paper finds that objects of interest, such as output elasticities with respect to inputs vary considerably with respect to the rank of unobserved technology shocks. In the application to US firm-level manufacturing data, this paper considers a Translog production function with non-Hicksian neutral productivity shocks and a productivity process that features non-linear persistence.

## Software Implementation
All code is written in R language and evaluated on an HPC system.

## Replication Files
The following folders contain replication files for the tables and figures in the paper.
1. [Empirical](/Empirical): A folder containing all R environments used in the analysis and files used to produce figures and tables
2. [Functions](/Functions): A folder containing all the replication files for estimation
	- [Production_EM](/Functions/Production_EM.R): Main estimation file. Uses functions below as dependencies.
	- [Posterior](/Functions/Posterior.R): Function used to compute to posterior density
	- [Mstep](/Functions/Mstep.R): Function to perform Maximization step of algorithm
	- [Tensors](/Functions/Tensors.R): Function used to specify functional forms for the production function, inputs, and productivity
	- [omega](/Functions/omega.R): Function used to estimate a production function using [Levinsohn and Petrin (2003)](https://doi.org/10.1111/1467-937X.00246)
	- [Auxiliary Files](/Functions/Auxfuns.R): Contains auxiliary files used in estimation procedures


