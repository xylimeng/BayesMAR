# BayesMAR
Bayesian Median Autoregressive model for time series forecasting 


## Example 

Refer to the R script ```BayesMAR.R```. 

Computational time may be different due to sample size and settings insides BayesMAR.R, e.g. numbers of iterations, acceptance region of auto-tuning (binary search) and etc

```r
source('BayesMAR.R')

#------------------------------------------------------------------
#  Example 1 Simulation
#------------------------------------------------------------------

# simulated data
yt = matrix(0,202,1)
for(i in 3:202){
  yt[i] = 0.3 + 0.75*yt[i-1] - 0.35*yt[i-2] + rnorm(1,0,1)
}
yt = yt[3:202]

# apply BayesMAR, with order 2 [AR(2) with intercept]
results = BMAR(yt, 2)

# list[[1]] coefficients
results[1]
# list[[2]] acceptance rate
results[2] 
# list[[3]] burn-in chain
results[3]
# list[[4]] chain after burn-in
results[4]

# 4-step ahead prediction
BMAR_pred(yt,results[[1]][3,],4)

#------------------------------------------------------------------
#  Example 2 BIC
#------------------------------------------------------------------
# Assume we want to choose model from order 1 to order 20 and use data generated above

# load existing codes provided by Koenker 2005, solving the MAE problem 
library(quantreg)
library(zoo)

# generate matrix to storing values
# use data generated above
BIC = matrix(0,20,1)

for( i in 1:20){
  BIC[i] = BMAR_BIC(yt,i,20)
}
# find the optimum order
which.min(BIC)

```
