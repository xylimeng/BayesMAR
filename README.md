# BayesMAR
Bayesian Median Autoregressive model for time series forecasting 


## Example 

Refer to the R script ```BayesMAR.R```. 
```r
# Example 1 simulation
# computational time may be different depending on
# sample size and settings insides BayesMAR.R
# e.g. numbers of iterations, acceptance region of auto-tuning

source('BayesMAR.R')
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
```
