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


# 6. CRPS
# CRPS for BayesMAR-MAP
diffp <- as.matrix(read.csv("diffr.csv"))[,2]
rt = cbind(diffp[162:196],diffp[163:197],
           diffp[164:198],diffp[165:199])
BIC_max = read.csv('BIC_max.csv')[,2:21]
order = apply(BIC_max, 1 , which.min)[1]
data_name <- 'diffp'
f1 = matrix(0,35,15000)
f2 = matrix(0,35,15000)
f3 = matrix(0,35,15000)
f4 = matrix(0,35,15000)

for( i in 1:35){
  beta = as.matrix(read.csv(paste0('rw_beta_samples_',i,'_order_',order,'.csv'))[,(2:(2+order))])
  y = get(data_name)[1:(160+i)]
  # sample tau
  a = length(y)-order
  Y = matrix(1,(1+order),(length(y)-order))
  for( j in 1:order){
    Y[(1+j),] = y[(160+i-j):(1+order-j)]
  }
  true_Y = matrix(rep(y[(160+i):(1+order)],15000), nrow = 15000, byrow = T)
  b = apply( abs( true_Y - beta%*%Y) ,1,sum)/2
  tau = matrix(0,15000,1)
  tau = 1/rgamma(15000,a,b)
  
  # f1 
  Y = matrix(1, 15000, (1+order))
  fill = matrix( rep( y[(160+i):(160+i-(order)+1)], 15000), nrow = 15000, byrow = T)
  Y[,2:(1+order)] = fill
  f1[i,] = apply(beta*Y,1,sum) + rlaplace(15000,0, (2*tau))
  # f2 
  fill = cbind(f1[i,], fill)
  Y[,2:(1+order)] = fill[,1:order]
  f2[i,] = apply(beta*Y,1,sum) + rlaplace(15000,0, (2*tau))
  # f3
  fill = cbind(f2[i,], fill)
  Y[,2:(1+order)] = fill[,1:order]
  f3[i,] = apply(beta*Y,1,sum) +  rlaplace(15000,0, (2*tau))
  # f4 
  fill = cbind(f3[i,], fill)
  Y[,2:(1+order)] = fill[,1:order]
  f4[i,] = apply(beta*Y,1,sum) +  rlaplace(15000,0, (2*tau))
}

bmarru = matrix(0,35,4)
bmarrm = matrix(0,35,4)
bmarrd = matrix(0,35,4)

for( j in 1:4){
  bmarru[,j] = apply( get(paste0('f',j)), 1, function(x){quantile(x,0.975)})
  bmarrm[,j] = apply( get(paste0('f',j)), 1, mean)
  bmarrd[,j] = apply( get(paste0('f',j)), 1, function(x){quantile(x,0.025)})
}

bmarr_crps = matrix(0, 35, 4)
for( i in 1:35){
  for( j in 1:4){
    dat = get(paste0('f',j))[i,]
    bmarr_crps[i,j] = crps_sample( y = rt[i,j], dat = dat)
  }
}
bmarr_crps = apply(bmarr_crps,2,mean)
```

## Reference

Zeng, Z. and Li, M. (2020). Bayesian Median Autoregression for Robust Time Series Forecasting. *International Journal of Forecasting*, accepted. <https://arxiv.org/pdf/2001.01116.pdf>
