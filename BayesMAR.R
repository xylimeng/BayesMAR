{
  # need package 'zoo' and 'quantreg'
  library(zoo)
  library(quantreg) 
  
  #-----------------------------------------------
  # main function -- BMAR
  #-----------------------------------------------
  # input:
  # Data -- time-series Data
  # Order -- order of ar process;
  #             intercept term will be included
  #             eg: order 1 is
  #               y_t = b1 + b2*y_{t-1}  
  
  # output:
  # list:
  # [[1]] beta
  # [[2]] accept rate
  # [[3]] burnin chain
  # [[4]] chain
  
  
  # init_f = 0
  # init_beta = runif(0,1)
  # iterations = 40,000 with 25,000 burnin
  # a = binar search, around 35% \pm 10%
  
  # A binary search algorithm is used to achieve tuning-free
  # A long-burnin period is used to ensure convergence
  # They are convenient for auto-recursive forecast, but will
  # hugely increase non-necessarily computational burnden.
  # Please manually change them if necessary 
  
  BMAR <- function( Data, order = 1){ 
    
    #-----------------------------------------------
    # check the input
    #-----------------------------------------------
    # Data
    SampleSize = length(Data)
    y = matrix(Data, SampleSize, 1)
    # order
    order = matrix( (order+1) , 1, 1)      
    
    #-----------------------------------------------
    # Burnin MCMC
    #-----------------------------------------------
    
    # init_beta 
    init_beta = runif(order,0,1)
    
    # find optim a
    a = 1
    r.b_a = 1
    l.b_a = 0
    
    tem_target = quick_accept_rate(y, order, a, 25000, init_beta)
    
    # ensure right bound of a is large enough
    while(tem_target > 0.45){
      r.b_a = r.b_a + 1
      a = r.b_a
      tem_target = quick_accept_rate(y, order, a, 25000, init_beta)
    }
    
    # find optim a by binary search
    while( abs(tem_target - 0.35) > 0.10 ){
      if( tem_target > 0.45 ){
        l.b_a = a
        a = (a + r.b_a)
        tem_target = quick_accept_rate(y, order, a, 25000, init_beta)
      }else{
        r.b_a = a
        a = (a + l.b_a)/2
        tem_target = quick_accept_rate(y, order, a, 25000, init_beta)
      }
    }
    
    # MCMC
    Burnin_chain <- MCMC(y, order, a, 25000, init_beta)
    
    #-----------------------------------------------
    # RW-MCMC
    #-----------------------------------------------
    
    # set the initial
    RW_initial <- Burnin_chain[25000,]
    
    # find optim a
    a = 1
    r.b_a = 1
    l.b_a = 0
    
    tem_target = quick_accept_rate(y, order, a, 10000, RW_initial)
    
    # ensure right bound of a is large enough
    while(tem_target > 0.45){
      r.b_a = r.b_a + 1
      a = r.b_a
      tem_target = quick_accept_rate(y, order, a, 10000, RW_initial)
    }
    
    # find optim a by binary search
    while( abs(tem_target - 0.35) > 0.10 ){
      if( tem_target > 0.45 ){
        l.b_a = a
        a = (a + r.b_a)
        tem_target = quick_accept_rate(y, order, a, 10000, RW_initial)
      }else{
        r.b_a = a
        a = (a + l.b_a)/2
        tem_target = quick_accept_rate(y, order, a, 10000, RW_initial)
      }
    }
    
    # MCMC
    RW_chain <- MCMC(y, order, a, 15000, RW_initial)
    
    #-----------------------------------------------
    # report results
    #-----------------------------------------------
    
    Accept = matrix(0, 1, 2)
    Accept[1] = 1 - mean( duplicated(Burnin_chain[,1]))
    Accept[2] = 1 - mean( duplicated(RW_chain[,1]))
    
    colnames(Accept) <- c('burn_in','random_walk')
    
    # report beta and variance
    mean_BI = apply( Burnin_chain, 2, mean)
    sd_BI = apply( Burnin_chain, 2, sd)
    mean_RW = apply( RW_chain, 2, mean)
    sd_RW = apply( RW_chain, 2, sd)
    beta <- matrix(c( mean_BI, sd_BI, mean_RW, sd_RW), 4, order, byrow = T)
    rownames(beta) <- c('BI_Coeff.','BI_S.E','RW_Coeff.','RW_S.E')
    
    output <- list(beta,Accept,Burnin_chain,RW_chain)
    
    return(output)
  }
  
  #-----------------------------------------------
  #       proposal
  #-----------------------------------------------
  #   using t-proposal
  #   beta_{proposal} = beta_{last period} + t
  #   where t ~ t(5)
  #  
  #   input row vector \beta
  #   parameter a
  #   return row vector \beta
  proposal_beta <- function(param, a){
    
    # check input
    order = length(param)
    param <- matrix(param, 1, order)
    
    # update proposal
    param = param + t(diag(a,order)%*%runif(order,-0.1,0.1))
    
    return(param)
  }
  
  #-----------------------------------------------
  #       MCMC Sampler
  #-----------------------------------------------
  MCMC <- function( y, order, a, Iterations, init_beta){
    #-----------------------------------------------
    # conclude infromations needs
    #-----------------------------------------------  
    
    # creat the matrix to store chain
    chain <- array( dim = c(Iterations, order))
    # input initial value of the chain
    chain[1,] <- matrix(init_beta , 1, order)
    
    #-----------------------------------------------
    # MCMC process
    #-----------------------------------------------
    for ( i in 2:Iterations){
      # propose
      proposal = proposal_beta( chain[i-1,], a)
      
      # accept/reject
      prob = exp( posterior(proposal, y) - posterior( chain[i-1,], y))
      
      if ( runif(1) < prob ){
        chain[i,] <- proposal
      }else{
        chain[i,] <- chain[i-1,]
      }
      
    }
    return(chain)
  } 
  
  #-----------------------------------------------
  #       quick-accept-rate
  #-----------------------------------------------
  
  quick_accept_rate = function( y, order, a, Iterations, init_beta){
    tem_chain = MCMC(y, order, a, Iterations, init_beta)
    tem_accept = 1 - mean(duplicated(tem_chain[,1]))
    return(tem_accept)
  }
  
  #-----------------------------------------------
  #       posterior function
  #-----------------------------------------------
  #    priors of beta_i : \propto 1
  #    Likelihood : Skewed-Laplace distribution ~ SL(0,\tau,\alpha)
  #    posterior is [\sum^{n}_{t=2} \frac{1}{2} \abs{ y_t - f_{t}(\beta) }]^{-n}
  posterior <- function(param, y){
    #-----------------------------------------------
    # check input
    #-----------------------------------------------
    order = length(param)
    beta = matrix(param, order, 1)
    SampleSize = length(y)
    y <- matrix(y, SampleSize, 1)
    f <- matrix(0,SampleSize,1)
    
    # set scalar to store sum
    t <- 0 
    
    #-----------------------------------------------
    # calculate posterior for given Data
    #-----------------------------------------------
    Y = matrix(1, order, SampleSize-order+1)
    for(i in 2:order){
      Y[i,] = y[(SampleSize-i+1):(order-i+1)]
    }
    
    f[SampleSize:order] = t(beta) %*% Y 
    f[SampleSize:order] = (abs(y[SampleSize:order] - f[SampleSize:order]))/2
    
    t <- (-(SampleSize-order+1))*log(sum(f[SampleSize:order]))
    
    return(t)
  }
  
  #---------------------------------------
  #   BMAR_predict
  #---------------------------------------
  
  # input： 
  #      data, beta, predict-step
  # output prediction
  
  BMAR_pred = function( Data, beta, step = 1){
    p = length(beta)
    y = as.matrix(Data[ length(Data) : (length(Data) - p)])
    tem_y = as.matrix(c(1,y))
    Forecast = matrix(0,1,step)
    
    for(i in 1:step){
      Forecast[1,i] = t(beta) %*% tem_y[1:p]
      y = rbind(Forecast[1,i],y)
      tem_y = as.matrix( c(1,y))
    }
    return(Forecast)
  }
  
  BMAR_BIC = function( Data, order, max_order = 20){
    y = Data
    size = length(Data)
    
    data = matrix( y[(order+1):size], length(y[(order+1):size]),(order+1))
    name = c('y')  
    for( j in 1:order){
      data[,(j+1)] = y[(order+1-j):(size-j)]
      name = cbind(name, paste('y',j,sep=''))
    }
    data = split(data, rep( 1:ncol(data), each = nrow(data)))
    names(data) = name
    qar = dynrq(reformulate(name[2:(order+1)], "y"), tau = 0.5, data = data)
    
    b = qar$coefficients
    p = length(b) + 1
    n = length(qar$y) # adjusted obs
    
    # S for \sum | y_t - \hat{y_t} |/2
    S = sum( abs( qar$residuals) )/2
    # tau 
    MAP.tau = S/(n+2) 
    
    # likelihood
    MAP.L = (4*MAP.tau)^(-(size - max_order)) * exp( -MAP.tau^(-1)*sum( abs(qar$residuals[ (1+(max_order-order)):(size - max_order+(max_order-order)) ]))/2)
    return( log(n)*(p) - 2*log(MAP.L) )
  }
}

