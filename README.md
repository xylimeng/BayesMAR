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
# Assume we want to choose model from order 1 to order 20 and use 
# data generated above, which means, sample size = 200 and max_p = 20

# load existing codes provided by Koenker 2005, solving the MAE problem 
library(quantreg)
library(zoo)

# generate matrix to storing values
# use data generated above
BIC = matrix(0,20,1)
y = yt

# loop from order 1 to order 20
# 1. transform data, satisfying the form used by <quantreg>
# 2. calculate BIC for each order
# Note: when calculate BIC, we have 180 = 200 - 20 = Sample Size - Max_p

for( i in 1:20){
# 1. 
  data = matrix( y[(i+1):200], length(y[(i+1):(200)]),(i+1))
  name = c('y')  
  for( j in 1:i){
    data[,(j+1)] = y[(i+1-j):(200-j)]
    name = cbind(name, paste('y',j,sep=''))
  }
  data = split(data, rep( 1:ncol(data), each = nrow(data)))
  names(data) = name
  qar = dynrq(reformulate(name[2:(i+1)], "y"), tau = 0.5, data = data)
  
# 2.
  b = qar$coefficients
  p = length(b) + 1
  n = length(qar$y) # adjusted obs
  
  # S for \sum | y_t - \hat{y_t} |/2
  S = sum( abs( qar$residuals) )/2
  # tau 
  MAP.tau = S/(n+2) 
  
  # likelihood
  MAP.L = (4*MAP.tau)^(-(180)) * exp( -MAP.tau^(-1)*sum( abs(qar$residuals[ (1+(20-i)):(180+(20-i)) ]))/2)
  BIC[i] = log(n)*(p) - 2*log(MAP.L)
}
# find the optimum order
which.min(BIC)
```
