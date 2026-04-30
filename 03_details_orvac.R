library(prodlim) ## load libraries
library(dfoptim)
library(numDeriv)
library(survival)
library(survcomp)
library(simsurv)
library(coxed)
library(R.utils)

w <- 1 ## learning parameter
fsims <- 500 # Number of sims for fut (final predictive probability)
ssims <- fsims # Number of sims for suc (interim predictive probability)
success_thres <- seq(0.95, 0.995, 0.005) # Overall rule for success (several options for tuning decision thresholds)
c <- seq(1, 4, by = 0.2) ## relationship between first interim sample size and other cumulative-stage sample sizes
numT <- length(c) ## number of analyses
Smax <- 250 # Maximum first interim analysis sample size
recruit <- 50/3 ## recruitment rate
Tmax <- 40 # Median survival time (without censoring) is half of this
xknots <- c(0.00, 0.25, 0.50, 0.75, 1.00) ## splines for baseline hazard (data generation)
yknots <- c(0.1,0.25,0.3,0.035,0.05)
max_follow <- 30 ## maximum follow up time

prior_mu_super <- 0 ## prior distribution on beta in generalized posterior
prior_sd_super <- 10

## precompute summary statistics for logpl
## inputs are survival times (surv.time) and
## indicator of whether censored (surv.event)
precompute_logpl <- function(surv.time, surv.event) {
  
  ## reorder the data so that the earliest events come first
  ord <- order(surv.time, decreasing = TRUE)
  time  <- surv.time[ord]
  event <- surv.event[ord]
  
  r <- rle(time) ## computes runs of identical times
  ends <- cumsum(r$lengths) ## note when runs start and end
  starts <- c(1, head(ends, -1) + 1)
  
  list(ord = ord, starts = starts, ends = ends,
       event = event)
}

## faster function to compute log partial likelihood of
## a Cox model (O(nlogn) instead of O(n^2))
## pred: linear predictor in Cox model
## pre: preprocessed summary statistics
logpl_fast_precomp <- function(pred, pre) {
  eta   <- pred[pre$ord] ## extract predictors in decreasing order
  event <- pre$event
  
  max_eta <- max(eta) ## log-sum-exp stabilization (if all predictors are large/small)
  exp_eta <- exp(eta - max_eta)
  cum_risk <- cumsum(exp_eta)
  
  loglik <- 0
  
  starts <- pre$starts
  ends   <- pre$ends
  
  for (k in seq_along(starts)) {
    idx <- starts[k]:ends[k] ## get all obs at that time
    d_t <- sum(event[idx]) ## how many events at this time
    
    if (d_t > 0) {
      ## log denominator with correction
      log_denom <- log(cum_risk[starts[k]]) + max_eta
      
      loglik <- loglik + sum(eta[idx]*event[idx]) - d_t*log_denom
    }
  }
  
  loglik
}

## function to compute loss function based on faster 
## helper functions; inputs are as follows:
## beta: coefficient in partial likelihood
## d: treatment assignment
## pre: pre-processed data summaries
laplace_loss_fast <- function(beta, d, pre) {
  
  eta <- beta*d ## linear predictor
  
  log_prior <- sum(dnorm(beta, prior_mu_super, prior_sd_super, log = TRUE)) ## log-prior
  
  ## fast log partial likelihood
  logpl_val <- logpl_fast_precomp(pred = eta, pre = pre)[1]
  
  loss_fun <- -logpl_val ## self-information loss
  
  log_loss <- -w*loss_fun
  
  -(log_prior + log_loss)
}

# Function to get spline for baseline hazard (data generation)
ihaz <- smooth.spline(x=xknots, y=yknots, df=3)
ihaz_fun <- function(t){
  predict(ihaz,x = t/Tmax)$y
}

## simulate survival data for all stages; the inputs are
## all_x: treatment assignments
## ihaz_fun: function for the baseline hazard
sim_surv <- function(all_x, ihaz_fun){
  
  ## simulate data for both treatments separately
  all_data <- rep(-1, length(all_x))
  all_data[which(all_x==0)] <- sim.survdata_fast(sum(all_x == 0), Tmax, 
                                                 ihaz_fun, x = 0, beta = prior_mu)
  all_data[which(all_x==1)] <- sim.survdata_fast(sum(all_x == 1), Tmax, 
                                                 ihaz_fun, x = 1, beta = prior_mu)
  
  as.numeric(all_data)
  
}

# Run futility rule 
fut_super <- function(fsims, mu0_loss, sigma0_loss,
                      prior_mu, Smax, x,
                      all_x, temp_haz,
                      current_y, out){
  fut <- rep(0,fsims)
  for(j in 1:fsims){
    fbeta <- rnorm(1,mu0_loss,sqrt(sigma0_loss))
    fn <- Smax - length(current_y)
    fxi <- all_x[(1+length(current_y)):Smax]
    
    ##################################
    simdatax0 <- sim.survdata_fast(N=1000, T=Tmax, hazard.fun = temp_haz, x = 0, beta = fbeta)
    simdatax1 <- sim.survdata_fast(N=1000, T=Tmax, hazard.fun = temp_haz, x = 1, beta = fbeta)
    y01 <- rbind(simdatax0,simdatax1)
    ##################################
    
    # Simulate left censored data
    fy <- current_y
    ind <- which(out == 0)
    y_temp <- rep(0,length(ind))
    for(k in 1:length(ind)){
      pre_sim_data <- y01[x[ind[k]]+1,]
      jk <- which(pre_sim_data > current_y[ind[k]])
      if(length(jk) > 0){y_temp[k] <- sample(x=pre_sim_data[jk],size=1)}else{y_temp[k] <- current_y[ind[k]]}
    }
    fy[ind] <- y_temp
    
    ## simulate untruncated data for future observations
    fyy <- rep(0,fn)
    fyy[which(fxi==0)] <- simdatax0[sample(1:1000, sum(fxi==0), replace = TRUE)]
    fyy[which(fxi==1)] <- simdatax1[sample(1:1000, sum(fxi==1), replace = TRUE)]
    
    ## update trt assignment and censoring indicators
    fx <- c(x,fxi)
    fy <- c(fy,fyy)
    fcen <- rep(1,Smax)
    
    ind <- which(fy > max_follow)
    fcen[ind] <- 0
    
    fy <- pmin(fy, max_follow)
    
    fpre <- precompute_logpl(surv.time = fy + runif(length(fy),0,1e-7), surv.event = fcen)
    
    out_loss <- optim(par=prior_mu,fn = laplace_loss_fast,d=fx,pre = fpre,
                      hessian=TRUE,method="L-BFGS-B",lower=-10,upper=10)
    fmu0_loss <- out_loss$par
    fhess_loss <- out_loss$hessian
    fsigma0_loss <- 1/fhess_loss
    ## get complementary posterior probability based on predictive data
    fut[j] <- pnorm(q = 0, mean = fmu0_loss, sd = sqrt(fsigma0_loss), lower.tail = FALSE)
  }
  fut
}

# Run superiority rule 
suc_super <- function(ssims, mu0_loss, sigma0_loss,
                      prior_mu, x,
                      temp_haz,
                      current_y, out){
  suc <- rep(0,ssims)
  for(j in 1:ssims){
    sbeta <- rnorm(1,mu0_loss,sqrt(sigma0_loss))
    
    ##################################
    simdatax0 <- sim.survdata_fast(N=1000, T=Tmax, hazard.fun = temp_haz, x = 0, beta = sbeta)
    simdatax1 <- sim.survdata_fast(N=1000, T=Tmax, hazard.fun = temp_haz, x = 1, beta = sbeta)
    y01 <- rbind(simdatax0,simdatax1)
    ##################################
    
    # Simulate left censored data
    sy <- current_y
    ind <- which(out == 0)
    y_temp <- rep(0,length(ind))
    for(k in 1:length(ind)){
      pre_sim_data <- y01[x[ind[k]]+1,]
      jk <- which(pre_sim_data > current_y[ind[k]])
      if(length(jk) > 0){y_temp[k] <- sample(x=pre_sim_data[jk],size=1)}else{y_temp[k] <- current_y[ind[k]]}
    }
    sy[ind] <- y_temp
    
    ## no untruncated data for future stages in this case; otherwise, same process as futility
    sx <- x
    scen <- rep(1,length(sy))
    
    ind <- which(sy > max_follow)
    scen[ind] <- 0
    
    sy <- pmin(sy, max_follow)
    
    spre <- precompute_logpl(surv.time = sy + runif(length(sy),0,1e-7), surv.event = scen)
    
    out_loss <- optim(par=prior_mu,fn = laplace_loss_fast,d=sx,pre = spre,
                      hessian=TRUE,method="L-BFGS-B",lower=-10,upper=10)
    
    smu0_loss <- out_loss$par
    shess_loss <- out_loss$hessian
    ssigma0_loss <- 1/shess_loss
    suc[j] <- pnorm(q = 0, mean = smu0_loss, sd = sqrt(ssigma0_loss), lower.tail = FALSE)
  }
  suc
}

# Log-likelihood for super model (uniform/constant prior)
## inputs are as follows:
## theta: y-values for knots in spline
## yt: survival times (rounded)
## x: treatment assignments
## cen: censoring indicators
## beta: coefficent from generalized posterior
sur_loglike <- function(theta,yt,x,cen,beta){
  
  haz_like <- smooth.spline(x=xknots, y=theta, df=3)
  haz <- function(t){
    predict(haz_like, x = t/Tmax)$y
  }
  
  # From hazard to pdf/sur
  #######################
  time <- 0:Tmax
  haz0 <- haz(time)
  surv0 <- exp(-cumsum(haz0))
  failCDF0 <- 1 - surv0
  failPDF0 <- c(0, diff(failCDF0))
  
  haz1 <- haz0*exp(beta)
  surv1 <- surv0^exp(beta)
  failCDF1 <- 1 - surv1
  failPDF1 <- c(0, diff(failCDF1))
  
  #######################
  ## construct the log-likelihood
  
  ind <- which(cen==1)
  ytemp <- yt[ind]
  xtemp <- x[ind]
  pdf <- ifelse(xtemp==1, failPDF1[ytemp + 1], failPDF0[ytemp +1]) 
  
  ind <- which(cen==0)
  ytemp <- yt[ind]
  xtemp <- x[ind]
  sur <- ifelse(xtemp==1, surv1[ytemp + 1], surv0[ytemp +1])
  
  like <- c(pdf,sur)
  like[like < 1e-20] <- 1e-20
  out <- -sum(log(like))
  if(!is.finite(out)){out <- 100000}
  out
}

## this is a faster alternative to sim.surv_data()
## the inputs are as follows:
## N: number of observations
## Tmax: maximum follow up time
## hazard.fun: function for baseline hazard
## x: value for trt assignment (scalar)
## beta: regression coefficient value in partial likelihood
sim.survdata_fast <- function(N, Tmax, hazard.fun, x, beta) {
  ## Compute linear predictor for all individuals
  eta <- if (x == 1) beta else 0
  
  ## Baseline hazards
  h0 <- hazard.fun(1:Tmax)
  
  ## Adjust for proportional hazards
  H <- 1 - (1 - h0)^exp(eta)
  
  ## Simulate uniform random numbers for each individual
  U <- matrix(runif(N*Tmax), nrow = N, ncol = Tmax)
  
  ## Event occurs at first time interval where U < H
  events <- U < matrix(H, nrow = N, ncol = Tmax, byrow = TRUE)
  event_time <- max.col(events, ties.method = "first")  # vectorized
  
  ## If no event occurred (all FALSE), set to Tmax
  event_time[!rowSums(events)] <- Tmax
  
  event_time
}

## function to get the logit of the posterior probability
## the inputs are as follows:
## probs: complementary posterior probabilities returned by futility or superiority procedure
## thres: decision threshold for success at analysis t (gamma not xi)
getLogits <- function(probs, thres){
  
  qs <- log(probs) - log(1 - probs)
  
  ## get complementary posterior predictive probability using kernel density estimate
  kd.qs <- density(na.omit(qs))
  np.prob <- mean(pnorm(log(1 - thres) - log(thres), ifelse(is.finite(qs), qs,
                                                            ifelse(qs > 0,  max(na.omit(qs)) + 1,
                                                                   min(na.omit(qs)) - 1)),
                        kd.qs$bw, lower.tail = FALSE))
  ## if greater than 0.5, calculate on complementary scale to avoid rounding to 1
  if (np.prob > 0.5){
    np.prob <- mean(pnorm(log(1 - thres) - log(thres), ifelse(is.finite(qs), qs,
                                                              ifelse(qs > 0,  max(na.omit(qs)) + 1,
                                                                     min(na.omit(qs)) - 1)),
                          kd.qs$bw, lower.tail = TRUE))
    
    ## if probability is very small use normal approximation (more stable)
    if (np.prob < 0.00004){
      temp <- pnorm(log(1 - thres) - log(thres), mean(na.omit(qs)), sd(na.omit(qs)), lower.tail = TRUE)
      return(qlogis(temp))
    } else {
      return(qlogis(np.prob))
    }
    
  }
  
  ## repeat for case where original probability was small
  if (np.prob < 0.00004){
    temp <- pnorm(log(1 - thres) - log(thres), mean(na.omit(qs)), sd(na.omit(qs)), lower.tail = FALSE)
    return(-1*qlogis(temp))
  } else {
    return(-1*qlogis(np.prob))
  }
  
  ## the logit of the posterior predictive probability is returned (not on complementary scale)
}