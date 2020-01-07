# BayesMAR
Bayesian Median Autoregressive model for time series forecasting 


## Example 

Note: Computational time depends on the sample size and settings for MCMC, e.g. the number of iterations, acceptance region of auto-tuning (binary search), etc.


### Example 1. Simulation
```r
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

# 4-step ahead prediction
BMAR_pred(yt,results[[1]][3,],4)
# return 1:h=4 ahead predicitons
```

### Example 2. Real data application using 3-Month Treasury Bill: secondary market Rate (TBR) data


```r 
#  We demonstrate how to use BIC to implement Bayes_MAP and Bayes_BMA in this example. 

source('BayesMAR.R')

# 1. BayesMAR prediction
# We have saved results of step 1 (i.e., the variable "fp") into "TBR_BayesMAR_forecast.RData" for easy reference.
# "fp" stores the array of prediction: fp[ time = 1:35, BayesMAR_order = 1:20, step_ahead = 1:4 ]
# Skip the rest of Step 1 by loading the .RData file 
load("TBR_BayesMAR_forecast.RData")

# Otherwise, run the following code: 
# read data & generate array 'fp' to store predictions  
diffr <- as.matrix(read.csv("diffr.csv"))[,2]
fp <- array(NA, c(35, 20, 4))
# recursive forecast, via each oder
for( p in 1:20){ # order in BayesMAR
  for (k in 1:35){ # timepoints for recursive forecasting 
      # recursively get data
      y <- diffr[1:(160+k)]
      # forecast
      r <- BMAR(y, p)
      fp[k,p,] <- BMAR_pred( y, r[[1]][3,], 4)
  }
}

# 2. BIC
BIC = matrix(0,35,20)
# recursive data, via order
for( k in 1:35){ # timepoints for recursive forecasting
  y = diffr[1:(160+k)]
  for( p in 1:20){ # order in BayesMAR
    BIC[k,p] = BMAR_BIC(y,p,20)
  }
}

# 3. BayesMAR - MAP
# pick the optimum order by the BIC at the beginning period
p = which.min(BIC[1,])
BayesMAR_MAP = fp[,p,]

# 4. BayesMAR - BMA
# model weight determined by the BIC value at the beginning period 
bic = BIC[1,]
exp_bic = exp(-bic/2)
weights = exp_bic/sum(exp_bic)

BayesMAR_BMA = matrix(0, 35, 4)
for( k in 1:35){ # timepoints for recursive forecasting
  for( h in 1:4){ # h-step ahead prediction
    BayesMAR_BMA[k,h] = sum( fp[k,,h]* weights)
  }
}

# 5. RMSE and MAE 
# read raw data tr / predict target r
tr <- read.csv(paste('tr.csv',sep=''))
r <- matrix(0,35,4)
for(h in 1:4){ # h-step ahead true data
  r[,h] <- tr[(162+h):(196+h),2]
}

# prediction = previous raw data + predictive changes
# prediction for BayesMAR_MAP
tem <- matrix(0,35,4)
tem[,1] <- BayesMAR_MAP[,1] + tr[162:196,2]
for( h in 1:3){ # (h+1)-step ahead prediction
  tem[,(1+h)] = BayesMAR_MAP[,(1+h)]+tem[,h]
}
# RMSE and MAE for BayesMAR_MAP
sqrt(apply((tem - r)^2,2,mean))
apply( abs(tem - r), 2, mean)

# prediction for BayesMAR_BMA
tem[,1] <- BayesMAR_BMA[,1] + tr[162:196,2]
for( h in 1:3){ # (h+1)-step ahead prediction
  tem[,(1+h)] = BayesMAR_BMA[,(1+h)]+tem[,h]
}
# RMSE and MAE for BayesMAR_BMA
sqrt(apply((tem - r)^2,2,mean))
apply( abs(tem - r), 2, mean)
```

## Reference

Zijian Zeng and Meng Li (2020). Bayesian Median Autoregression for Robust Time Series Forecasting. <https://arxiv.org/pdf/2001.01116.pdf>
