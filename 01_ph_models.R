## BEGIN SETUP ##

## load necessary packages
require(foreach)
require(doParallel)
require(doSNOW)
require(MASS)
require(grid)
require(ggplot2)
require(cowplot)
require(ggpubr)

expit <- function(x){1/(1 + exp(-x))}
logit <- function(x){log(x) - log(1-x)}

## set up parallelization
cores=detectCores()
cl <- makeSOCKcluster(71)

## number of draws M for posterior predictive probabilities
draws <- 1000

## simulation repetitions
m <- 10000
registerDoSNOW(cl)
pb <- txtProgressBar(max = m, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

## set up prior distributions for exponential baseline
prior_mu <- rep(0, 2); prior_sd <- rep(sqrt(10), 2)

## function to compute the log-posterior for the model 
## with the exponential baseline hazard; inputs are as follows:
## par: values of logeta and theta
## precomp: precomputed summary statistics (see next function)
## prior_mu: prior means for logeta and theta
## prior_sd: prior standard deviations
laplace_exp <- function(par, precomp, prior_mu, prior_sd) {
  
  ## extract parameter values and summary statistics
  logeta <- par[1]
  theta <- par[2]
  
  ddot <- precomp[1]
  dAdot <- precomp[2]
  y0dot <- precomp[3]
  y1dot <- precomp[4]
  
  log_prior <- sum(dnorm(par, prior_mu, prior_sd, log = TRUE))
  
  ## fast log partial likelihood
  log_lik <- ddot*logeta + theta*dAdot - exp(logeta)*y0dot - exp(logeta)*y1dot*exp(theta)
  
  -(log_prior + log_lik)
}

## function to simulate data from the model with exponential baseline
## the inputs are as follows:
## n: sample size
## eta: value for baseline parameter
## theta: hazard ratio
sim_exp <- function(n, eta, theta){
  u <- runif(n)
  A <- round(runif(n))
  return(cbind(A, -log(1-u)/(eta*exp(theta*A))))
}

## function to precompute summary statistics from the data;
## this is done to avoid recomputing these within optim to get
## the Laplace approximation
## the inputs are as follows:
## dat: data frame of treatments and observed times
## delta: vector of censoring status (1 = observed, 0 = censored)
precomp_exp <- function(dat, delta){
  A <- dat[,1]
  t_obs <- dat[,2]
  
  ddot  <- sum(delta)
  dAdot <- sum(delta*A)
  y0dot <- sum(t_obs*(1 - A))
  y1dot <- sum(t_obs*A)
  
  return(c(ddot, dAdot, y0dot, y1dot))
} 

## function to simulate truncated data from the model with exponential baseline
## the inputs are as follows:
## y: the current data (left bound for truncation)
## A: current treatment assignments
## eta: value for baseline parameter
## theta: hazard ratio
sim_trunc_exp <- function(y, A, eta, theta){
  rate.temp <- eta*exp(theta*A)
  return(y + rexp(length(y), rate.temp))
}

tic <- Sys.time()
## set up exponential parameters
exp.par <- c(0.4,-0.8)
## sample sizes considered
n_vals <- seq(40,100,5)
## consider various (complementary) threshold values of posterior predictive probabilities
thres1 <- 1 - 0.975
thres2 <- 1 - 0.95
thres3 <- 1 - 0.9
thres4 <- 1 - 0.99
for (k in 1:length(n_vals)){
  n <- n_vals[k]
  sim.res <- foreach(i=1:m, .combine=rbind, .packages = c("MASS"),
                     .options.snow=opts) %dopar% {
                       
                       set.seed(k*m + i + 100000)
                       ## simulate initial data
                       dat.temp <- sim_exp(n, exp.par[1], exp.par[2])
                       
                       ## precompute summary statistics and get the Laplace approximation
                       precomp.temp <- precomp_exp(cbind(dat.temp[,1], pmin(dat.temp[, 2], 5)), 
                                                         delta = as.integer(dat.temp[,2] <= 5))
                       
                       out_loss <- optim(par=c(0,0),fn = laplace_exp, precomp= precomp.temp,
                                         prior_mu = prior_mu, prior_sd = prior_sd,
                                         hessian=TRUE,method="L-BFGS-B",lower=rep(-10, 2),upper=rep(10, 2))
                       mu0_loss <- out_loss$par
                       hess_loss <- out_loss$hessian
                       sigma0_loss <- solve(hess_loss)
                       
                       ## get posterior probability
                       p1 <- pnorm(q = 0, mean = mu0_loss[2], sd = sqrt(sigma0_loss[2,2]), lower.tail = FALSE)
                       
                       ## get the posterior mode
                       m1 <- mu0_loss[2]
                       
                       ## now consider an interim predictive probability 
                       ## get the starting time
                       t0 <- cumsum(rexp(n, 10))
                       x0 <- dat.temp[,2]
                       artif <- ifelse(t0 < max(t0) - 5, 0,
                                       ifelse(x0 + t0 > max(t0), 1, 0))
                       
                       ## compute posterior based on artificially censored data
                       delta.temp <- as.numeric(x0 < 5)
                       delta.temp <- ifelse(artif, 0, delta.temp)
                       x <- ifelse(delta.temp, x0, pmin(5, max(t0) - t0))
                       
                       precomp.temp <- precomp_exp(cbind(dat.temp[,1], x), delta = delta.temp)
                       
                       ## get the Laplace approximation
                       out_loss <- optim(par=c(0,0),fn = laplace_exp, precomp= precomp.temp,
                                         prior_mu = prior_mu, prior_sd = prior_sd,
                                         hessian=TRUE,method="L-BFGS-B",lower=rep(-10, 2),upper=rep(10, 2))
                       mu0_loss <- out_loss$par
                       hess_loss <- out_loss$hessian
                       sigma0_loss <- solve(hess_loss)
                       
                       ## get the posterior mode
                       m2 <- mu0_loss[2]
                       
                       ## generate new eta and theta values from current posterior
                       new.pars <- mvrnorm(draws, mu0_loss, sigma0_loss)
                       new.pars[,1] <- exp(new.pars[,1])
                       
                       prob.temp <- NULL
                       for (j in 1:draws){
                         
                         ## simulate truncated data and replace for artificially censored observations
                         dat.pred <- sim_trunc_exp(dat.temp[,2], dat.temp[,1], new.pars[j,1], new.pars[j,2])
                         dat.pred <- ifelse(artif, dat.pred, dat.temp[,2])
                         
                         ## get Laplace approximation for posterior based on predictive data
                         precomp.temp <- precomp_exp(cbind(dat.temp[,1], dat.pred), delta = as.integer(dat.pred <= 5))
                         
                         out_loss <- optim(par=c(0,0),fn = laplace_exp, precomp= precomp.temp,
                                           prior_mu = prior_mu, prior_sd = prior_sd,
                                           hessian=TRUE,method="L-BFGS-B",lower=rep(-10, 2),upper=rep(10, 2))
                         mu0_loss <- out_loss$par
                         hess_loss <- out_loss$hessian
                         sigma0_loss <- solve(hess_loss)
                         
                         ## get (complementary) posterior probability (better for round-off)
                         evidence <- pnorm(q = 0, mean = mu0_loss[2], sd = sqrt(sigma0_loss[2,2]), lower.tail = FALSE)
                         prob.temp <- c(prob.temp, evidence)
                       }
                       ## take the logits of posterior probabilities and get kernel density estimate of their distribution
                       lp <- logit(prob.temp)
                       lp <- ifelse(is.finite(lp), lp,
                                ifelse(lp > 0, max(na.omit(lp))+1, min(na.omit(lp))-1))
                       kd.lp <- density(na.omit(lp))
                       
                       ## compute complementary posterior predictive probabilities based on each threshold value
                       pp1 <- mean(pnorm(logit(thres1), lp, kd.lp$bw, lower.tail = FALSE))
                       pp2 <- mean(pnorm(logit(thres2), lp, kd.lp$bw, lower.tail = FALSE))
                       pp3 <- mean(pnorm(logit(thres3), lp, kd.lp$bw, lower.tail = FALSE))
                       pp4 <- mean(pnorm(logit(thres4), lp, kd.lp$bw, lower.tail = FALSE))
                       
                       ## if probabilities are very small, use normal approximation to distribution (more stable)
                       if (pp1 < 0.00004){
                         pp1 <- pnorm(logit(thres1), mean(lp), sd(lp), lower.tail = FALSE)
                       } 
                       
                       if (pp2 < 0.00004){
                         pp2 <- pnorm(logit(thres2), mean(lp), sd(lp), lower.tail = FALSE)
                       } 
                       
                       if (pp3 < 0.00004){
                         pp3 <- pnorm(logit(thres3), mean(lp), sd(lp), lower.tail = FALSE)
                       } 
                       
                       if (pp4 < 0.00004){
                         pp4 <- pnorm(logit(thres4), mean(lp), sd(lp), lower.tail = FALSE)
                       } 
                       
                       ## return posterior summaries and modes
                       c(p1, pp1, pp2, pp3, m1, m2, pp4)
                     }
  
  ## write results to a .csv file
  write.csv(sim.res, paste0("sec5p1_d12_n_", n, "_H1.csv"), row.names = FALSE)   
}

toc <- Sys.time()
toc -tic

## set up Weibull parameters for designs 3 and 4
wei.par <- c(0.08, 2.25, -1.145)

## function to simulate data with Weibull baseline hazard; the inputs are
## n: sample size
## eta: scale parameter for Weibull baseline
## rho: shape parameter for Weibull baseline
## theta: hazard ratio
sim_weibull <- function(n, eta, rho, theta){
  u <- runif(n)
  A <- round(runif(n))
  return(cbind(A, (-log(1-u)/(eta*exp(theta*A)))^(1/rho)))
}

## repeat the process for designs 3 and 4
tic <- Sys.time()
n_vals <- seq(40,100,5)
thres1 <- 1 - 0.975
thres2 <- 1 - 0.95
thres3 <- 1 - 0.9
thres4 <- 1 - 0.99
for (k in 1:length(n_vals)){
  n <- n_vals[k]
  sim.res <- foreach(i=1:m, .combine=rbind, .packages = c("MASS"),
                     .options.snow=opts) %dopar% {
                       
                       set.seed(k*m + i + 200000)
                       
                       ## simulate initial Weibull baseline data
                       dat.temp <- sim_weibull(n, wei.par[1], wei.par[2], wei.par[3])
                       
                       ## get Laplace approximation for design 3
                       precomp.temp <- precomp_exp(cbind(dat.temp[,1], pmin(dat.temp[, 2], 5)), 
                                                   delta = as.integer(dat.temp[,2] <= 5))
                       
                       out_loss <- optim(par=c(0,0),fn = laplace_exp, precomp= precomp.temp,
                                         prior_mu = prior_mu, prior_sd = prior_sd,
                                         hessian=TRUE,method="L-BFGS-B",lower=rep(-10, 2),upper=rep(10, 2))
                       mu0_loss <- out_loss$par
                       hess_loss <- out_loss$hessian
                       sigma0_loss <- solve(hess_loss)
                       
                       ## get posterior probability
                       p1 <- pnorm(q = 0, mean = mu0_loss[2], sd = sqrt(sigma0_loss[2,2]), lower.tail = FALSE)
                       
                       ## get the posterior mode
                       m1 <- mu0_loss[2]
                       
                       ## now consider an interim predictive probability 
                       ## get the starting time
                       t0 <- cumsum(rexp(n, 10))
                       x0 <- dat.temp[,2]
                       artif <- ifelse(t0 < max(t0) - 5, 0,
                                       ifelse(x0 + t0 > max(t0), 1, 0))
                       
                       ## compute posterior based on censored data
                       delta.temp <- as.numeric(x0 < 5)
                       delta.temp <- ifelse(artif, 0, delta.temp)
                       x <- ifelse(delta.temp, x0, pmin(5, max(t0) - t0))
                       
                       ## get the Laplace approximation for design 4
                       precomp.temp <- precomp_exp(cbind(dat.temp[,1], x), delta = delta.temp)
                       
                       out_loss <- optim(par=c(0,0),fn = laplace_exp, precomp= precomp.temp,
                                         prior_mu = prior_mu, prior_sd = prior_sd,
                                         hessian=TRUE,method="L-BFGS-B",lower=rep(-10, 2),upper=rep(10, 2))
                       mu0_loss <- out_loss$par
                       hess_loss <- out_loss$hessian
                       sigma0_loss <- solve(hess_loss)
                       
                       ## get the posterior mode
                       m2 <- mu0_loss[2]
                       
                       new.pars <- mvrnorm(draws, mu0_loss, sigma0_loss)
                       new.pars[,1] <- exp(new.pars[,1])
                       
                       prob.temp <- NULL
                       for (j in 1:draws){
                         
                         ## simulate truncated data from model with *exponential* baseline
                         dat.pred <- sim_trunc_exp(dat.temp[,2], dat.temp[,1], new.pars[j,1], new.pars[j,2])
                         dat.pred <- ifelse(artif, dat.pred, dat.temp[,2])
                         
                         precomp.temp <- precomp_exp(cbind(dat.temp[,1], dat.pred), delta = as.integer(dat.pred <= 5))
                         
                         out_loss <- optim(par=c(0,0),fn = laplace_exp, precomp= precomp.temp,
                                           prior_mu = prior_mu, prior_sd = prior_sd,
                                           hessian=TRUE,method="L-BFGS-B",lower=rep(-10, 2),upper=rep(10, 2))
                         mu0_loss <- out_loss$par
                         hess_loss <- out_loss$hessian
                         sigma0_loss <- solve(hess_loss)
                         
                         ## get complementary posterior probability
                         evidence <- pnorm(q = 0, mean = mu0_loss[2], sd = sqrt(sigma0_loss[2,2]), lower.tail = FALSE)
                         prob.temp <- c(prob.temp, evidence)
                       }
                       lp <- logit(prob.temp)
                       lp <- ifelse(is.finite(lp), lp,
                                    ifelse(lp > 0, max(na.omit(lp))+1, min(na.omit(lp))-1))
                       kd.lp <- density(na.omit(lp))
                       
                       ## compute various complementary interim predictive probabilities
                       pp1 <- mean(pnorm(logit(thres1), lp, kd.lp$bw, lower.tail = FALSE))
                       pp2 <- mean(pnorm(logit(thres2), lp, kd.lp$bw, lower.tail = FALSE))
                       pp3 <- mean(pnorm(logit(thres3), lp, kd.lp$bw, lower.tail = FALSE))
                       pp4 <- mean(pnorm(logit(thres4), lp, kd.lp$bw, lower.tail = FALSE))
                       
                       if (pp1 < 0.00004){
                         pp1 <- pnorm(logit(thres1), mean(lp), sd(lp), lower.tail = FALSE)
                       } 
                       
                       if (pp2 < 0.00004){
                         pp2 <- pnorm(logit(thres2), mean(lp), sd(lp), lower.tail = FALSE)
                       } 
                       
                       if (pp3 < 0.00004){
                         pp3 <- pnorm(logit(thres3), mean(lp), sd(lp), lower.tail = FALSE)
                       } 
                       
                       if (pp4 < 0.00004){
                         pp4 <- pnorm(logit(thres4), mean(lp), sd(lp), lower.tail = FALSE)
                       } 
                       
                       ## output various posterior summaries and modes
                       c(p1, pp1, pp2, pp3, m1, m2, pp4)
                     }
  
  ## return results in .csv file
  write.csv(sim.res, paste0("sec5p1_d34_n_", n, "_H1.csv"), row.names = FALSE)   
}

toc <- Sys.time()
toc -tic

## set parameters for designs 5 and 6
wei.par <- c(0.062, 2.25, -0.8)

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
## helper functions
laplace_loss_fast <- function(beta, d, pre) {
  
  eta <- beta*d ## linear predictor
  
  log_prior <- sum(dnorm(beta, prior_mu_super, prior_sd_super, log = TRUE)) ## log-prior
  
  ## fast log partial likelihood
  logpl_val <- logpl_fast_precomp(pred = eta, pre = pre)[1]
  
  loss_fun <- -logpl_val ## self-information loss
  
  log_loss <- -w*loss_fun
  
  -(log_prior + log_loss)
}

## function to approximate the cumulative baseline hazard using
## Breslow's estimator; the inputs are as follows:
## dat: data.frame with columns (trt, time, delta)
## theta: hazard ratio
breslow_baseline <- function(dat, theta) {
  
  ## extract treatment, time, and censoring status
  x <- dat$trt
  time <- dat$time
  delta <- dat$delta
  
  ## calculate risk scores based on treatment
  risk <- exp(theta*x)
  
  ## resort the vectors based on observed time
  ord <- order(time)
  time <- time[ord]
  delta <- delta[ord]
  risk <- risk[ord]
  
  ## get the unique set of event times
  event_times <- sort(unique(time[delta == 1]))
  
  ## precompute the total risk set sum
  total_risk <- sum(risk)
  
  ## initialize outputs
  H0 <- NULL; dH0 <- NULL
  
  ## running index
  j <- 1; current_risk <- total_risk
  
  for (k in 1:length(event_times)) {
    t <- event_times[k]
    
    d_j <- 0 ## count events at this time
    
    ## remove individuals with time < t from the risk set
    while (j <= length(time) && time[j] < t) {
      current_risk <- current_risk - risk[j]
      j <- j + 1
    }
    
    ## count events at time t
    idx <- j
    while (idx <= length(time) && time[idx] == t) {
      if (delta[idx] == 1) {
        d_j <- d_j + 1
      }
      idx <- idx + 1
    }
    
    ## Breslow increment
    dH0[k] <- d_j/current_risk
    
    ## output cumulative hazard
    H0[k] <- if (k == 1) dH0[k] else H0[k-1] + dH0[k]
  }
  
  ## return list of times and cumulative baseline hazard
  ## points are added for the endpoints of the time interval [0,5]
  return(list(time = c(0, event_times, 5),
              H0 = c(0, H0, pmax(2, max(H0) + 0.1))))
}

## this function will make a sampler that can be used to 
## simulate truncated data from the survival model
## the inputs are as follows:
## H0_grid: grid of cumulative baseline hazard estimates
## t_grid: grid of times
## t_max: maximum follow-up time
make_sampler <- function(H0_grid, t_grid, t_max = 5) {
  
  ## ensure that H0 is always positive and increasing
  H0_grid <- cummax(H0_grid)
  eps <- 1e-10
  H0_grid <- H0_grid + eps*seq_along(H0_grid)
  
  ## get maximum baseline hazard over follow-up
  H_max <- as.numeric(tail(H0_grid, 1))
  
  ## pre-compute inverse cumulative baseline hazard
  inv_H0 <- function(y) {
    approx(H0_grid, t_grid,
           xout = y, method = "linear",
           ties = "ordered", rule = 2)$y}
  
  ## pre-compute cumulative baseline hazard (for truncation)
  H0_fun <- function(t) {
    approx(t_grid, H0_grid,
           xout = t,
           method = "linear", rule = 2)$y}
  
  ## return a sampler with the previous functions pre-computed
  ## inputs are sample size (n), trt assignments (A), hazard ratio (theta),
  ## and L (truncation times)
  function(n, A, theta, L = NULL) {
    ## draw uniform realizations
    u <- runif(n)
    
    ## compute untruncated input into inverse cumulative baseline hazard
    y <- -log(u)/exp(A * theta)
    
    ## apply left truncation
    H_L <- H0_fun(L)
    y <- H_L + y
    
    ## determine whether events occur
    event <- y <= H_max
    
    ## cap at boundary (for inversion later)
    y_cap <- pmin(y, H_max)
    
    ## get the times by inversion
    t_sim <- inv_H0(y_cap)
    
    ## apply censoring at end of follow-up
    t_obs <- pmin(t_sim, t_max)
    
    ## return times and event indicators
    return(list(time = t_obs,
      event = as.integer(event)))
  }
}

## set prior distributions and learning rate for designs 5 and 6
prior_mu_super <- 0
prior_sd_super <- sqrt(10)
w <- 1

## repeat the process for generalized
tic <- Sys.time()
n_vals <- seq(40,100,5)
thres1 <- 1 - 0.975
thres2 <- 1 - 0.95
thres3 <- 1 - 0.9
thres4 <- 1 - 0.99
for (k in 1:length(n_vals)){
  n <- n_vals[k]
  sim.res <- foreach(i=1:m, .combine=rbind, .packages = c("MASS", "mgcv"),
                     .options.snow=opts) %dopar% {
                       
                       set.seed(k*m + i + 300000)
                       
                       ## simulate data with Weibull baseline
                       dat.temp <- sim_weibull(n, wei.par[1], wei.par[2], wei.par[3])
                       
                       current_y <- pmin(dat.temp[,2], 5)
                       cen <- ifelse(dat.temp[,2] <= 5, 1, 0)
                       
                       ## get posterior for design 5
                       pre <- precompute_logpl(surv.time = current_y, surv.event = cen)
                       
                       out_loss <- optim(par=wei.par[3],fn = laplace_loss_fast, d=dat.temp[,1],pre=pre,
                                         hessian=TRUE,method="L-BFGS-B",lower=-10,upper=10)
                       mu0_loss <- out_loss$par
                       hess_loss <- out_loss$hessian
                       sigma0_loss <- 1/hess_loss
                       
                       ## get posterior probability
                       p1 <- pnorm(q = 0, mean = mu0_loss, sd = sqrt(sigma0_loss), lower.tail = FALSE)
                       
                       ## get the posterior mode
                       m1 <- mu0_loss
                       
                       ## now consider an interim predictive probability 
                       ## get the starting time
                       t0 <- cumsum(rexp(n, 10))
                       x0 <- dat.temp[,2]
                       artif <- ifelse(t0 < max(t0) - 5, 0,
                                       ifelse(x0 + t0 > max(t0), 1, 0))
                       
                       ## compute posterior based on censored data for design 6
                       delta.temp <- as.numeric(x0 < 5)
                       delta.temp <- ifelse(artif, 0, delta.temp)
                       current_y <- ifelse(delta.temp, x0, pmin(5, max(t0) - t0))
                       
                       pre <- precompute_logpl(surv.time = current_y, surv.event = delta.temp)
                       
                       out_loss <- optim(par=wei.par[3],fn = laplace_loss_fast, d=dat.temp[,1],pre=pre,
                                         hessian=TRUE,method="L-BFGS-B",lower=-10,upper=10)
                       mu0_loss <- out_loss$par
                       hess_loss <- out_loss$hessian
                       sigma0_loss <- 1/hess_loss
                       
                       ## get the posterior mode
                       m2 <- mu0_loss
                       
                       ## construct simulation function for posterior predictive data
                       dat.temp2 <- data.frame(trt = dat.temp[,1], time = current_y, 
                                               delta = delta.temp)
                       
                       ## get approximation to cumulative baseline hazard and make sampler function
                       bres.temp <- breslow_baseline(dat.temp2, mu0_loss)
                       make_sampler.i <- make_sampler(bres.temp$H0, bres.temp$time)
                       
                       ## draw new hazard ratios
                       new.par <- rnorm(draws, mu0_loss, sigma0_loss)
                       
                       ## get data for artficially censored observations
                       n.art <- sum(artif == 1)
                       y.art <- current_y[artif == 1]
                       A.art <- dat.temp[artif==1, 1]
                       
                       prob.temp <- NULL
                       for (j in 1:draws){
                         
                         dat.pred.y <- make_sampler.i(n.art, A.art, theta = new.par[j], L = y.art)$time
                         dat.pred <- dat.temp2$time
                         dat.pred[artif==1] <- dat.pred.y
                         
                         cen.art <- ifelse(dat.pred <= 5, 1, 0)
                         
                         ## get posterior for the predictive data
                         pre <- precompute_logpl(surv.time = dat.pred, surv.event = cen.art)
                         
                         out_loss <- optim(par=new.par[j],fn = laplace_loss_fast, d=dat.temp[,1],pre=pre,
                                           hessian=TRUE,method="L-BFGS-B",lower=-10,upper=10)
                         mu0_loss <- out_loss$par
                         hess_loss <- out_loss$hessian
                         sigma0_loss <- 1/hess_loss
                         
                         ## get complementary posterior probability
                         evidence <- pnorm(q = 0, mean = mu0_loss, sd = sqrt(sigma0_loss), lower.tail = FALSE)
                         prob.temp <- c(prob.temp, evidence)
                       }
                       lp <- logit(prob.temp)
                       lp <- ifelse(is.finite(lp), lp,
                                    ifelse(lp > 0, max(na.omit(lp))+1, min(na.omit(lp))-1))
                       kd.lp <- density(na.omit(lp))
                       
                       ## compute various interim predictive probabilities
                       pp1 <- mean(pnorm(logit(thres1), lp, kd.lp$bw, lower.tail = FALSE))
                       pp2 <- mean(pnorm(logit(thres2), lp, kd.lp$bw, lower.tail = FALSE))
                       pp3 <- mean(pnorm(logit(thres3), lp, kd.lp$bw, lower.tail = FALSE))
                       pp4 <- mean(pnorm(logit(thres4), lp, kd.lp$bw, lower.tail = FALSE))
                       
                       if (pp1 < 0.00004){
                         pp1 <- pnorm(logit(thres1), mean(lp), sd(lp), lower.tail = FALSE)
                       } 
                       
                       if (pp2 < 0.00004){
                         pp2 <- pnorm(logit(thres2), mean(lp), sd(lp), lower.tail = FALSE)
                       } 
                       
                       if (pp3 < 0.00004){
                         pp3 <- pnorm(logit(thres3), mean(lp), sd(lp), lower.tail = FALSE)
                       } 
                       
                       if (pp4 < 0.00004){
                         pp4 <- pnorm(logit(thres4), mean(lp), sd(lp), lower.tail = FALSE)
                       } 
                       
                       ## return posterior summaries and modes
                       c(p1, pp1, pp2, pp3, m1, m2, pp4)
                     }
  
  ## output results in a .csv file
  write.csv(sim.res, paste0("sec5p1_d56_n_", n, "_H1.csv"), row.names = FALSE)   
}

toc <- Sys.time()
toc -tic

## this function approximates the power curve for designs 1, 3, and 5
## (i.e., using linear approximations for n); the inputs are as follows:
## m1: matrix of complementary probabilities at first sample size (from .csv files)
## m2: matrix of complementary probabilities at second sample size (from .csv files)
## n0: first sample size
## n1: second sample size
## lb: lower bound for which to compute OCs
## ub: upper bound for which to compute OCs
## by: increments between lower and upper bound
## gam: decision threshold for success
pwr_one_lin <- function(m1, m2, n0, n1, lb, ub, by, gam){
  
  ## get logits for the p-values in one tail of 
  ls <- l1(m1)
  li <- l1(m2)
  
  ## adjust any infinite logits
  ls <- ifelse(ls == -Inf, min(subset(ls, is.finite(ls))) - 1, ls)
  ls <- ifelse(ls == Inf, max(subset(ls, is.finite(ls))) + 1, ls)
  
  ## adjust any infinite logits
  li <- ifelse(li == -Inf, min(subset(li, is.finite(li))) - 1, li)
  li <- ifelse(li == Inf, max(subset(li, is.finite(li))) + 1, li)
  
  slopes <- NULL
  ints <- NULL
  
  ## construct the slopes using the order statistics of the sampling
  ## distribution estimates
  ls_s <- ls[order(ls)]
  li_s <- li[order(li)]
  
  slopes <- (li_s - ls_s)/(n1-n0)
  ints <- ls_s - slopes*n0
  
  ## create matrix to calculate power
  samps <- seq(lb,ub,by)
  res.vec <- NULL
  for (i in 1:length(samps)){
    ## check which probabilities are less than the threshold 
    stop.temp <- ints + slopes*samps[i] >= logit(gam)
    res.vec[i] <- mean(stop.temp)
  }
  
  ## return matrix with samples sizes and estimated power
  return(cbind(samps, res.vec))
  
}

## this function approximates the power curve for designs 2, 4, and 6
## (i.e., using linear approximations for n^2); the inputs are as follows:
## m1: matrix of complementary probabilities at first sample size (from .csv files)
## m2: matrix of complementary probabilities at second sample size (from .csv files)
## n0: first sample size
## n1: second sample size
## lb: lower bound for which to compute OCs
## ub: upper bound for which to compute OCs
## by: increments between lower and upper bound
## xi: decision threshold for success
pwr_one_quad <- function(m1, m2, n0, n1, lb, ub, by, xi){
  
  ## get logits for the p-values in one tail of 
  ls <- l1(m1)
  li <- l1(m2)
  
  ## adjust any infinite logits
  ls <- ifelse(ls == -Inf, min(subset(ls, is.finite(ls))) - 1, ls)
  ls <- ifelse(ls == Inf, max(subset(ls, is.finite(ls))) + 1, ls)
  
  ## adjust any infinite logits
  li <- ifelse(li == -Inf, min(subset(li, is.finite(li))) - 1, li)
  li <- ifelse(li == Inf, max(subset(li, is.finite(li))) + 1, li)
  
  slopes <- NULL
  ints <- NULL
  
  ## construct the slopes using the order statistics of the sampling
  ## distribution estimates
  ls_s <- ls[order(ls)]
  li_s <- li[order(li)]
  
  slopes <- (li_s - ls_s)/(n1^2-n0^2)
  ints <- ls_s - slopes*n0^2
  
  ## create matrix to calculate power
  samps <- seq(lb,ub,by)
  res.vec <- NULL
  for (i in 1:length(samps)){
    ## check which probabilities are less than the threshold 
    stop.temp <- ints + slopes*samps[i]^2 >= logit(xi)
    res.vec[i] <- mean(stop.temp)
  }
  
  ## return matrix with samples sizes and estimated power
  return(cbind(samps, res.vec))
  
}

## create the first figure for this section
pwr.lin <- pwr_one_lin(read.csv("sec5p1_d12_n_45_H1.csv")[,1],
                       read.csv("sec5p1_d12_n_95_H1.csv")[,1],
                       45, 95, 40, 100, 1, 0.99)

pwr.sim <- NULL
for (i in 1:length(seq(40, 100, 5))){
  pwr.temp <- mean(1 - read.csv(paste0("sec5p1_d12_n_",seq(40, 100, 5)[i],"_H1.csv"))[,1] >= 0.99)
  pwr.sim <- c(pwr.sim, pwr.temp)
}

## combine different power curve estimates into one data frame
df1 <- data.frame(n = c(n_vals, seq(40, 100, 1)),
                  power = c(pwr.sim, pwr.lin[,2]),
                  curve = c(rep("C_Simulation", length(n_vals)),
                            rep("A_Algorithm 2", length(seq(40, 100, 1)))))

## create subplot
plot1 <- ggplot(df1, aes(x=n)) + theme_bw() +
  geom_line(aes(y = power, color=as.factor(curve), linetype = as.factor(curve)), 
            alpha = 0.9, size = 1) +
  labs(title=bquote('Design 1')) +
  labs(x= '', y= '') +
  theme(plot.title = element_text(size=18,face="bold", hjust =  0.5,
                                  margin = margin(t = 0, 0, 5, 0))) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18)) +
  theme(legend.position="none") +
  scale_color_manual(name = " ", 
                     labels = c("Estimated  ", "Simulated"),
                     values = c("steelblue1", "firebrick")) +
  scale_linetype_manual(name = " ", 
                        labels = c("Estimated  ", "Simulated"),
                        values = c(1, 5)) +
  theme(legend.text=element_text(size=18)) +
  theme(legend.key.size = unit(1, "cm")) +
  ylim(0,1) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

## get the subfigure for design 2
pwr.lin <- pwr_one_quad(read.csv("sec5p1_d12_n_45_H1.csv")[,7],
                       read.csv("sec5p1_d12_n_95_H1.csv")[,7],
                       45, 95, 40, 100, 1, 0.95)

pwr.sim <- NULL
for (i in 1:length(seq(40, 100, 5))){
  pwr.temp <- mean(1 - read.csv(paste0("sec5p1_d12_n_",seq(40, 100, 5)[i],"_H1.csv"))[,7] >= 0.95)
  pwr.sim <- c(pwr.sim, pwr.temp)
}

## combine different power curve estimates into one data frame
df2 <- data.frame(n = c(n_vals, seq(40, 100, 1)),
                  power = c(pwr.sim, pwr.lin[,2]),
                  curve = c(rep("C_Simulation", length(n_vals)),
                            rep("A_Algorithm 2", length(seq(40, 100, 1)))))

## create subplot
plot2 <- ggplot(df2, aes(x=n)) + theme_bw() +
  geom_line(aes(y = power, color=as.factor(curve), linetype = as.factor(curve)), 
            alpha = 0.9, size = 1) +
  labs(title=bquote('Design 2')) +
  labs(x= bquote(italic(n)), y= '') +
  theme(plot.title = element_text(size=18,face="bold", hjust =  0.5,
                                  margin = margin(t = 0, 0, 5, 0))) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18)) +
  theme(legend.position="none") +
  scale_color_manual(name = " ", 
                     labels = c("Estimated  ", "Simulated"),
                     values = c("steelblue1", "firebrick")) +
  scale_linetype_manual(name = " ", 
                        labels = c("Estimated  ", "Simulated"),
                        values = c(1, 5)) +
  theme(legend.text=element_text(size=18)) +
  theme(legend.key.size = unit(1, "cm")) +
  ylim(0,1) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) 

## get plot for design 3
pwr.lin <- pwr_one_lin(read.csv("sec5p1_d34_n_45_H1.csv")[,1],
                       read.csv("sec5p1_d34_n_95_H1.csv")[,1],
                       45, 95, 40, 100, 1, 0.99)

pwr.sim <- NULL
for (i in 1:length(seq(40, 100, 5))){
  pwr.temp <- mean(1 - read.csv(paste0("sec5p1_d34_n_",seq(40, 100, 5)[i],"_H1.csv"))[,1] >= 0.99)
  pwr.sim <- c(pwr.sim, pwr.temp)
}

## combine different power curve estimates into one data frame
df3 <- data.frame(n = c(n_vals, seq(40, 100, 1)),
                  power = c(pwr.sim, pwr.lin[,2]),
                  curve = c(rep("C_Simulation", length(n_vals)),
                            rep("A_Algorithm 2", length(seq(40, 100, 1)))))

## create subplot
plot3 <- ggplot(df3, aes(x=n)) + theme_bw() +
  geom_line(aes(y = power, color=as.factor(curve), linetype = as.factor(curve)), 
            alpha = 0.9, size = 1) +
  labs(title=bquote('Design 3')) +
  labs(x= '', y= '') +
  theme(plot.title = element_text(size=18,face="bold", hjust =  0.5,
                                  margin = margin(t = 0, 0, 5, 0))) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18)) +
  theme(legend.position="none") +
  scale_color_manual(name = " ", 
                     labels = c("Estimated  ", "Simulated"),
                     values = c("steelblue1", "firebrick")) +
  scale_linetype_manual(name = " ", 
                        labels = c("Estimated  ", "Simulated"),
                        values = c(1, 5)) +
  theme(legend.text=element_text(size=18)) +
  theme(legend.key.size = unit(1, "cm")) +
  ylim(0,1) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

## get the subfigure for design 4
pwr.lin <- pwr_one_quad(read.csv("sec5p1_d34_n_45_H1.csv")[,7],
                        read.csv("sec5p1_d34_n_95_H1.csv")[,7],
                        45, 95, 40, 100, 1, 0.95)

pwr.sim <- NULL
for (i in 1:length(seq(40, 100, 5))){
  pwr.temp <- mean(1 - read.csv(paste0("sec5p1_d34_n_",seq(40, 100, 5)[i],"_H1.csv"))[,7] >= 0.95)
  pwr.sim <- c(pwr.sim, pwr.temp)
}

## combine different power curve estimates into one data frame
df4 <- data.frame(n = c(n_vals, seq(40, 100, 1)),
                  power = c(pwr.sim, pwr.lin[,2]),
                  curve = c(rep("C_Simulation", length(n_vals)),
                            rep("A_Algorithm 2", length(seq(40, 100, 1)))))

## create subplot
plot4 <- ggplot(df4, aes(x=n)) + theme_bw() +
  geom_line(aes(y = power, color=as.factor(curve), linetype = as.factor(curve)), 
            alpha = 0.9, size = 1) +
  labs(title=bquote('Design 4')) +
  labs(x= bquote(italic(n)), y= '') +
  theme(plot.title = element_text(size=18,face="bold", hjust =  0.5,
                                  margin = margin(t = 0, 0, 5, 0))) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18)) +
  theme(legend.position="none") +
  scale_color_manual(name = " ", 
                     labels = c("Estimated  ", "Simulated"),
                     values = c("steelblue1", "firebrick")) +
  scale_linetype_manual(name = " ", 
                        labels = c("Estimated  ", "Simulated"),
                        values = c(1, 5)) +
  theme(legend.text=element_text(size=18)) +
  theme(legend.key.size = unit(1, "cm")) +
  ylim(0,1) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

## get plot for design 5
pwr.lin <- pwr_one_lin(read.csv("sec5p1_d56_n_45_H1.csv")[,1],
                       read.csv("sec5p1_d56_n_95_H1.csv")[,1],
                       45, 95, 40, 100, 1, 0.99)

pwr.sim <- NULL
for (i in 1:length(seq(40, 100, 5))){
  pwr.temp <- mean(1 - read.csv(paste0("sec5p1_d56_n_",seq(40, 100, 5)[i],"_H1.csv"))[,1] >= 0.99)
  pwr.sim <- c(pwr.sim, pwr.temp)
}

## combine different power curve estimates into one data frame
df5 <- data.frame(n = c(n_vals, seq(40, 100, 1)),
                  power = c(pwr.sim, pwr.lin[,2]),
                  curve = c(rep("C_Simulation", length(n_vals)),
                            rep("A_Algorithm 2", length(seq(40, 100, 1)))))

## create subplot
plot5 <- ggplot(df5, aes(x=n)) + theme_bw() +
  geom_line(aes(y = power, color=as.factor(curve), linetype = as.factor(curve)), 
            alpha = 0.9, size = 1) +
  labs(title=bquote('Design 5')) +
  labs(x= '', y= '') +
  theme(plot.title = element_text(size=18,face="bold", hjust =  0.5,
                                  margin = margin(t = 0, 0, 5, 0))) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18)) +
  theme(legend.position="none") +
  scale_color_manual(name = " ", 
                     labels = c("Estimated  ", "Simulated"),
                     values = c("steelblue1", "firebrick")) +
  scale_linetype_manual(name = " ", 
                        labels = c("Estimated  ", "Simulated"),
                        values = c(1, 5)) +
  theme(legend.text=element_text(size=18)) +
  theme(legend.key.size = unit(1, "cm")) +
  ylim(0,1) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

## get the subfigure for design 6
pwr.lin <- pwr_one_quad(read.csv("sec5p1_d56_n_45_H1.csv")[,3],
                        read.csv("sec5p1_d56_n_95_H1.csv")[,3],
                        45, 95, 40, 100, 1, 0.8)

pwr.sim <- NULL
for (i in 1:length(seq(40, 100, 5))){
  pwr.temp <- mean(1 - read.csv(paste0("sec5p1_d56_n_",seq(40, 100, 5)[i],"_H1.csv"))[,3] >= 0.8)
  pwr.sim <- c(pwr.sim, pwr.temp)
}

## combine different power curve estimates into one data frame
df6 <- data.frame(n = c(n_vals, seq(40, 100, 1)),
                  power = c(pwr.sim, pwr.lin[,2]),
                  curve = c(rep("C_Simulation", length(n_vals)),
                            rep("A_Algorithm 2", length(seq(40, 100, 1)))))

## create subplot
plot6 <- ggplot(df6, aes(x=n)) + theme_bw() +
  geom_line(aes(y = power, color=as.factor(curve), linetype = as.factor(curve)), 
            alpha = 0.9, size = 1) +
  labs(title=bquote('Design 6')) +
  labs(x= bquote(italic(n)), y= '') +
  theme(plot.title = element_text(size=18,face="bold", hjust =  0.5,
                                  margin = margin(t = 0, 0, 5, 0))) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18)) +
  theme(legend.position="none") +
  scale_color_manual(name = " ", 
                     labels = c("Estimated  ", "Simulated"),
                     values = c("steelblue1", "firebrick")) +
  scale_linetype_manual(name = " ", 
                        labels = c("Estimated  ", "Simulated"),
                        values = c(1, 5)) +
  theme(legend.text=element_text(size=18)) +
  theme(legend.key.size = unit(1, "cm")) +
  ylim(0,1) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

## get plot with legend
plot6.legend <- ggplot(df6, aes(x=n)) + theme_bw() +
  geom_line(aes(y = power, color=as.factor(curve), linetype = as.factor(curve)), 
            alpha = 0.9, size = 1) +
  labs(title=bquote('Design 6')) +
  labs(x= bquote(italic(n)), y= '') +
  theme(plot.title = element_text(size=18,face="bold", hjust =  0.5,
                                  margin = margin(t = 0, 0, 5, 0))) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18)) +
  theme(legend.position="bottom") +
  scale_color_manual(name = " ", 
                     labels = c("Estimated  ", "Simulated"),
                     values = c("steelblue1", "firebrick")) +
  scale_linetype_manual(name = " ", 
                        labels = c("Estimated  ", "Simulated"),
                        values = c(1, 5)) +
  theme(legend.text=element_text(size=18)) +
  theme(legend.key.size = unit(1, "cm")) +
  ylim(0,1) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(legend.key.size = unit(1.5, "cm"))

## compile the first figure

figp.row1 <- plot_grid(plot1 + theme(plot.margin=unit(c(0.2,0.3,0,-0.55),"cm")), 
                       plot3 + theme(plot.margin=unit(c(0.2,0.3,0,-0.45),"cm")),
                       plot5 + theme(plot.margin=unit(c(0.2,0.3,0,-0.45),"cm")),
                       rel_widths = c(1.11, 1, 1), nrow = 1)
figp.row2 <- plot_grid(plot2 + theme(plot.margin=unit(c(0.2,0.3,0,-0.55),"cm")), 
                       plot4 + theme(plot.margin=unit(c(0.2,0.3,0,-0.45),"cm")),
                       plot6 + theme(plot.margin=unit(c(0.2,0.3,0,-0.45),"cm")),
                       rel_widths = c(1.11, 1, 1), nrow = 1)

figp <- plot_grid(figp.row1, figp.row2, nrow = 2, rel_heights = c(1, 1.085))

## add y-axis label
y_label <- textGrob("Success Probability", rot = 90, gp = gpar(fontsize = 18))

fig.final.top <- plot_grid(y_label, figp, ncol = 2, rel_widths = c(0.5, 12))
fig.final <- plot_grid(fig.final.top, ggpubr::get_legend(plot6.legend), nrow = 2, rel_heights = c(12, 1.2))

# output as .pdf file for the article
pdf(file = "Fig_Sec5p1.pdf",   # The directory you want to save the file in
    width = 11, # The width of the plot in inches
    height = 8*11/12) # The height of the plot in inches

fig.final

dev.off()

## get the figure for the supplement

## extract the medians of the posterior mode
mode.sim <- NULL
for (i in 1:length(seq(40, 100, 5))){
  mode.temp <- cbind(read.csv(paste0("sec5p1_d12_n_",seq(40, 100, 5)[i],"_H1.csv"))[,5:6],
                     read.csv(paste0("sec5p1_d34_n_",seq(40, 100, 5)[i],"_H1.csv"))[,5:6],
                     read.csv(paste0("sec5p1_d5_n_",seq(40, 100, 5)[i],"_H1.csv"))[,5:6])
  mode.temp <- apply(mode.temp, 2, median)
  mode.sim <- rbind(mode.sim, mode.temp)
}

dfmode <- data.frame(n = rep(seq(40, 100, 5), 6), 
                     mode = as.numeric(mode.sim), 
                     curve = rep(1:6, each = 13))

plotmode <- ggplot(dfmode, aes(x=n)) + theme_bw() +
  geom_line(aes(y = mode, color=as.factor(curve), linetype = as.factor(curve)), 
            alpha = 0.9, size = 1) +
  labs(x= bquote(italic(n)), y= 'Median of the Posterior Mode') +
  theme(plot.title = element_text(size=18,face="bold", hjust =  0.5,
                                  margin = margin(t = 0, 0, 5, 0))) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18)) +
  theme(legend.position="bottom") +
  theme(legend.text=element_text(size=18), legend.title = element_text(size = 18)) +
  theme(legend.key.size = unit(1.5, "cm")) +
  ylim(-1,-0.75) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  geom_hline(yintercept=-0.8, lty = 2) +
  guides(color=guide_legend(title="Design"), linetype = guide_legend(title="Design"))


# output as .pdf file for the article
pdf(file = "FigMode.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 7) # The height of the plot in inches

plotmode

dev.off()