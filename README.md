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
#  Example 2 Replication of TBR data
#  Including BIC, Bayes_MAP and Bayes_BMA
#------------------------------------------------------------------
# We save results of step 1. for easy reference.
# the following codes return BayesMAR prediction array
# fp[ time = 1:35, BayesMAR_order = 1:20, step_ahead = 1:4 ]
load("~/Desktop/tem/TBR_BayesMAR_forecast.RData")

# 1. BayesMAR prediction
# could be skipped by loading the .RData file above
# read data, generate matrix to store
diffr <- as.matrix(read.csv("diffr.csv"))[,2]
f <- matrix(0,35,4)
# recursive forecast, via each oder
for( BMAR_order in 1:20){
  for (k in 1:35){
    # recursively get data
    y <- diffr[1:(160+k)]
    # forecast
    r <- BMAR( y, BMAR_order)
    f[k,] <- BMAR_pred( y, r[[1]][3,] , 4)
  }
  # save data
  write.csv( f, file = paste('diffr_order_',BMAR_order,'.csv',sep=''))
}


# 2. BIC
BIC = matrix(0,35,20)
# recursive data, via order
for( k in 1:35){
  y = diffr[1:(160+k)]
  for( BMAR_order in 1:20){
    BIC[k,BMAR_order] = BMAR_BIC(y,BMAR_order,20)
  }
}

# 3. BayesMAR - MAP
# pick the optimum order by the BIC at the beginning period
p = which.min(BIC[1,])
BayesMAR_MAP = fp[,p,]

# 4. BayesMAR - BMA
# decided by the BIC at the beginning period 
# calculate weight
bic = BIC[1,]
exp_bic = exp(-bic/2)
weights = exp_bic/sum(exp_bic)

# preidciton at time i, j-step head
# = sum( preidciton_order_p * weights_order_p )
BayesMAR_BMA = matrix(0, 35, 4)
for( i in 1:35){
  for( j in 1:4){
    BayesMAR_BMA[i,j] = sum( fp[i,,j]* weights)
  }
}

# 5. RMSE and MAE 
# read raw data tr / predict target r
tr <- read.csv(paste('tr.csv',sep=''))
r <- matrix(0,35,4)
for(i in 1:4){
  r[,i] <- tr[(162+i):(196+i),2]
}

# prediction = previous raw data + predictive changes
# prediction for BayesMAR_MAP
tem <- matrix(0,35,4)
tem[,1] <- BayesMAR_MAP[,1] + tr[162:196,2]
for( j in 1:3){
  tem[,(1+j)] = BayesMAR_MAP[,(1+j)]+tem[,j]
}
# RMSE and MAE for BayesMAR_MAP
sqrt(apply((tem - r)^2,2,mean))
apply( abs(tem - r), 2, mean)

# prediction for BayesMAR_BMA
tem[,1] <- BayesMAR_BMA[,1] + tr[162:196,2]
for( j in 1:3){
  tem[,(1+j)] = BayesMAR_BMA[,(1+j)]+tem[,j]
}
# RMSE and MAE for BayesMAR_BMA
sqrt(apply((tem - r)^2,2,mean))
apply( abs(tem - r), 2, mean)
```
