require(foreach)
require(doParallel)
require(doSNOW)
require(grid)
require(ggplot2)
require(cowplot)
require(ggpubr)

## set up parallelization with 10000 simulation repetitions
cores=detectCores()
cl <- makeSOCKcluster(71)

R <- 1000
registerDoSNOW(cl)
pb <- txtProgressBar(max = R, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

## set up parallelization
clusterEvalQ(cl, {
  source("03_details_orvac.R")
})
toc <- Sys.time()

## set effect sizes
source("03_details_orvac.R")
prior_mu <- -0.2 ## prior for general Bayes
prior_sd <- 10
ns <- seq(100, 250, 15) ## first interim sample sizes considered
logit <- function(x){log(x) - log(1-x)}

for (j in 1:length(ns)){
  n <- ns[j]
  Smax <- round(max(c)*n)
    
  ## simulate data for all trials
  mat_x <- matrix(round(runif(R*Smax)), nrow = R, ncol = Smax)
  mat_y <- t(apply(mat_x, 1, sim_surv, ihaz_fun = ihaz_fun))
  mat_t <- matrix(rexp(Smax, recruit), nrow = R, ncol = Smax)
  mat_t <- t(apply(mat_t, 1, function(x){round(cumsum(x))}))
  
  for (t in 1:(numT-1)){ ## interim simulation results
    res_int <- foreach(k=1:R, .packages=c('prodlim', 'dfoptim', 'numDeriv', 'survival',
                                          'survcomp', 'simsurv', 'coxed', 'R.utils'), .combine=rbind,
                       .errorhandling = "remove", .options.snow=opts) %dopar% {
                         
                         ## extract relevant data for analysis t
                         all_x <- mat_x[k,]
                         all_y <- mat_y[k,]
                         all_enroll <- mat_t[k,]
                         
                         ## get the time for analysis t
                         Smax_temp <- round(c[t]*n)
                         time_t <- all_enroll[Smax_temp]
                         
                         x <- all_x[1:Smax_temp]
                         y <- all_y[1:Smax_temp]
                         enroll <- all_enroll[1:Smax_temp]
                         
                         avail_time <- time_t - enroll
                         
                         ## apply articifical censoring and censoring due to max follow-up (out = 1)
                         out <- rep(0, Smax_temp)
                         cen <- rep(2, Smax_temp)
                         
                         current_y <- pmin(y, max_follow, avail_time)
                         cen <- ifelse(y <= pmin(max_follow, avail_time), 1, 0)
                         out <- ifelse(enroll + pmin(y, max_follow) <= time_t, 1, 0)
                         
                         ## get posterior for the interim analysis (precompute dat summaries)
                         pre <- precompute_logpl(surv.time = current_y+runif(length(current_y),0,1e-7), surv.event = cen)
                         
                         out_loss <- optim(par=prior_mu,fn = laplace_loss_fast, d=x,pre=pre,
                                           hessian=TRUE,method="L-BFGS-B",lower=-10,upper=10)
                         mu0_loss <- out_loss$par
                         hess_loss <- out_loss$hessian
                         sigma0_loss <- 1/hess_loss
                         
                         # Fit super model (here just take posterior modes)
                         skeep <- matrix(0,nrow=10,ncol=6)
                         for(ijk in 1:10){
                           opt <- optim(runif(5),fn = sur_loglike,yt=round(current_y),x=x,cen=cen,
                                        beta = mu0_loss, method="L-BFGS-B",
                                        lower=c(rep(0,5)),upper=c(rep(0.5,5)))
                           skeep[ijk,] <- c(opt$par,opt$value)
                         }
                         ind <- which.min(skeep[,6])[1]
                         skeep_out <- skeep[ind,1:5]
                         
                         ## bound the baseline hazard to reasonable values (if needed)
                         thaz <- smooth.spline(x=xknots, y=skeep_out, df=3)
                         temp_haz <- function(t){
                           pmin(pmax(0.02, predict(thaz,x = t/Tmax)$y), 0.5)
                         }
                         
                         ## run futility analysis
                         fut <- fut_super(fsims = fsims, mu0_loss = mu0_loss, sigma0_loss = sigma0_loss,
                                          prior_mu = prior_mu, Smax = Smax, x = x,
                                          all_x = all_x, temp_haz = temp_haz,
                                          current_y = current_y, out = out)
                         
                         ## get the logit of posterior predictive probabilities
                         ## for various candidate thresholds
                         fut_fp <- NULL
                         for (kk in 1:length(success_thres)){
                           fut_fp[kk] <- getLogits(fut, thres = success_thres[kk])
                         }
                         
                         ## run superiority analysis
                         suc <- suc_super(ssims = ssims, mu0_loss = mu0_loss, sigma0_loss = sigma0_loss,
                                          prior_mu = prior_mu, x = x, temp_haz = temp_haz,
                                          current_y = current_y, out = out)
                         
                         ## get the logit of posterior predictive probabilities
                         ## for various candidate thresholds
                         suc_ip <- NULL
                         for (kk in 1:length(success_thres)){
                           suc_ip[kk] <- getLogits(suc, thres = success_thres[kk])
                         }
                        
                        c(k, fut_fp, suc_ip) 
                       }
    
    ## record results in .csv file
    write.csv(res_int[,seq(1, length(success_thres) + 1, 1)], 
              paste0("fut_n_", n, "_t_", t, ".csv"), row.names = FALSE)
    write.csv(res_int[,c(1, seq(length(success_thres) + 2, 2*length(success_thres) + 1, 1))], 
              paste0("suc_n_", n, "_t_", t, ".csv"), row.names = FALSE)
  }  
  
  ## final analysis results 
  t <- numT
  
  res_fin <- foreach(k=1:R, .packages=c('prodlim', 'dfoptim', 'numDeriv', 'survival',
                                        'survcomp', 'simsurv', 'coxed', 'R.utils'), .combine=rbind,
                     .errorhandling = "remove", .options.snow=opts) %dopar% {
                       
                       all_x <- mat_x[k,]
                       all_y <- mat_y[k,]
                       
                       x <- all_x
                       y <- all_y
                       
                       current_y <- pmin(y, max_follow)
                       cen <- ifelse(y <= max_follow, 1, 0)
                       
                       ## get posterior for the final analysis
                       pre <- precompute_logpl(surv.time = current_y+runif(length(current_y),0,1e-7), surv.event = cen)
                       
                       out_loss <- optim(par=prior_mu,fn = laplace_loss_fast, d=x,pre=pre,
                                         hessian=TRUE,method="L-BFGS-B",lower=-10,upper=10)
                       mu0_loss <- out_loss$par
                       hess_loss <- out_loss$hessian
                       sigma0_loss <- 1/hess_loss
                       
                       ## get logit of posterior probabilty
                       if (mu0_loss < 0){
                         evidence <- pnorm(q = 0, mean = mu0_loss, sd = sqrt(sigma0_loss), lower.tail = FALSE)
                         evidence <- -1*qlogis(evidence)
                       } else {
                         evidence <- pnorm(q = 0, mean = mu0_loss, sd = sqrt(sigma0_loss), lower.tail = TRUE)
                         evidence <- qlogis(evidence)
                       }
                       
                       c(k, evidence)
                     }
  
  ## record the results
  write.csv(res_fin, paste0("suc_n_", n, "_t_", t, ".csv"), row.names = FALSE)
  
}  

## now get the sampling distribution estimates under H0

prior_mu <- 0 ## type I error
prior_sd <- 10
ns <- c(115, 235)

for (j in 1:length(ns)){
  n <- ns[j]
  Smax <- round(max(c)*n)
  
  ## simulate data for all trials
  mat_x <- matrix(round(runif(R*Smax)), nrow = R, ncol = Smax)
  mat_y <- t(apply(mat_x, 1, sim_surv, ihaz_fun = ihaz_fun))
  mat_t <- matrix(rexp(Smax, recruit), nrow = R, ncol = Smax)
  mat_t <- t(apply(mat_t, 1, function(x){round(cumsum(x))}))
  
  for (t in 1:(numT-1)){ ## interim simulation results
    res_int <- foreach(k=1:R, .packages=c('prodlim', 'dfoptim', 'numDeriv', 'survival',
                                          'survcomp', 'simsurv', 'coxed', 'R.utils'), .combine=rbind,
                       .errorhandling = "remove", .options.snow=opts) %dopar% {
                         
                         ## extract relevant data for analysis t
                         all_x <- mat_x[k,]
                         all_y <- mat_y[k,]
                         all_enroll <- mat_t[k,]
                         
                         ## get the time for analysis t
                         Smax_temp <- round(c[t]*n)
                         time_t <- all_enroll[Smax_temp]
                         
                         x <- all_x[1:Smax_temp]
                         y <- all_y[1:Smax_temp]
                         enroll <- all_enroll[1:Smax_temp]
                         
                         avail_time <- time_t - enroll
                         
                         ## apply articifical censoring and censoring due to max follow-up (out = 1)
                         out <- rep(0, Smax_temp)
                         cen <- rep(2, Smax_temp)
                         
                         current_y <- pmin(y, max_follow, avail_time)
                         cen <- ifelse(y <= pmin(max_follow, avail_time), 1, 0)
                         out <- ifelse(enroll + pmin(y, max_follow) <= time_t, 1, 0)
                         
                         ## get posterior for the interim analysis (precompute dat summaries)
                         pre <- precompute_logpl(surv.time = current_y+runif(length(current_y),0,1e-7), surv.event = cen)
                         
                         out_loss <- optim(par=prior_mu,fn = laplace_loss_fast, d=x,pre=pre,
                                           hessian=TRUE,method="L-BFGS-B",lower=-10,upper=10)
                         mu0_loss <- out_loss$par
                         hess_loss <- out_loss$hessian
                         sigma0_loss <- 1/hess_loss
                         
                         # Fit super model (here just take posterior modes)
                         skeep <- matrix(0,nrow=10,ncol=6)
                         for(ijk in 1:10){
                           opt <- optim(runif(5),fn = sur_loglike,yt=round(current_y),x=x,cen=cen,
                                        beta = mu0_loss, method="L-BFGS-B",
                                        lower=c(rep(0,5)),upper=c(rep(0.5,5)))
                           skeep[ijk,] <- c(opt$par,opt$value)
                         }
                         ind <- which.min(skeep[,6])[1]
                         skeep_out <- skeep[ind,1:5]
                         
                         ## bound the baseline hazard to reasonable values (if needed)
                         thaz <- smooth.spline(x=xknots, y=skeep_out, df=3)
                         temp_haz <- function(t){
                           pmin(pmax(0.02, predict(thaz,x = t/Tmax)$y), 0.5)
                         }
                         
                         ## run futility analysis
                         fut <- fut_super(fsims = fsims, mu0_loss = mu0_loss, sigma0_loss = sigma0_loss,
                                          prior_mu = prior_mu, Smax = Smax, x = x,
                                          all_x = all_x, temp_haz = temp_haz,
                                          current_y = current_y, out = out)
                         
                         ## get the logit of posterior predictive probabilities
                         ## for various candidate thresholds
                         fut_fp <- NULL
                         for (kk in 1:length(success_thres)){
                           fut_fp[kk] <- getLogits(fut, thres = success_thres[kk])
                         }
                         
                         ## run superiority analysis
                         suc <- suc_super(ssims = ssims, mu0_loss = mu0_loss, sigma0_loss = sigma0_loss,
                                          prior_mu = prior_mu, x = x, temp_haz = temp_haz,
                                          current_y = current_y, out = out)
                         
                         ## get the logit of posterior predictive probabilities
                         ## for various candidate thresholds
                         suc_ip <- NULL
                         for (kk in 1:length(success_thres)){
                           suc_ip[kk] <- getLogits(suc, thres = success_thres[kk])
                         }
                         
                         c(k, fut_fp, suc_ip) 
                       }
    
    ## record results in .csv file
    write.csv(res_int[,seq(1, length(success_thres) + 1, 1)], 
              paste0("H0_fut_n_", n, "_t_", t, ".csv"), row.names = FALSE)
    write.csv(res_int[,c(1, seq(length(success_thres) + 2, 2*length(success_thres) + 1, 1))], 
              paste0("H0_suc_n_", n, "_t_", t, ".csv"), row.names = FALSE)
  }  
  
  ## final analysis results 
  t <- numT
  
  res_fin <- foreach(k=1:R, .packages=c('prodlim', 'dfoptim', 'numDeriv', 'survival',
                                        'survcomp', 'simsurv', 'coxed', 'R.utils'), .combine=rbind,
                     .errorhandling = "remove", .options.snow=opts) %dopar% {
                       
                       all_x <- mat_x[k,]
                       all_y <- mat_y[k,]
                       
                       x <- all_x
                       y <- all_y
                       
                       current_y <- pmin(y, max_follow)
                       cen <- ifelse(y <= max_follow, 1, 0)
                       
                       ## get posterior for the final analysis
                       pre <- precompute_logpl(surv.time = current_y+runif(length(current_y),0,1e-7), surv.event = cen)
                       
                       out_loss <- optim(par=prior_mu,fn = laplace_loss_fast, d=x,pre=pre,
                                         hessian=TRUE,method="L-BFGS-B",lower=-10,upper=10)
                       mu0_loss <- out_loss$par
                       hess_loss <- out_loss$hessian
                       sigma0_loss <- 1/hess_loss
                       
                       ## get logit of posterior probability
                       if (mu0_loss < 0){
                         evidence <- pnorm(q = 0, mean = mu0_loss, sd = sqrt(sigma0_loss), lower.tail = FALSE)
                         evidence <- -1*qlogis(evidence)
                       } else {
                         evidence <- pnorm(q = 0, mean = mu0_loss, sd = sqrt(sigma0_loss), lower.tail = TRUE)
                         evidence <- qlogis(evidence)
                       }
                       
                       c(k, evidence)
                     }
  
  ## record the results
  write.csv(res_fin, paste0("H0_suc_n_", n, "_t_", t, ".csv"), row.names = FALSE)
  
} 

## this function computes stopping probabilities for a group sequential design using 
## the sampling distribution estimates contained in .csv files generated above
stopMatData <- function(mat, gam, xi){
  
  ## the inputs are described as follows:
  ## mat: a matrix of posterior summaries (obtained from .csv files)
  ## gam: the vector of sucess thresholds for all analyses
  ## xi: the vector of failure thresholds for the first T-1 analyses
  
  res.mat <- NULL
  ## check which logits are larger than first success threshold based on
  ## the linear approximations
  j <- 1
  stop.p <- mat[,j] >= logit(gam[j])
  res.mat[j] <- mean(stop.p)
  
  ## repeat process for failure thresholds at first analysis
  stop.f <- mat[,j+1] < logit(xi[j])
  ## make sure we didn't just stop for success
  stop.f <- ifelse(stop.p, 0, stop.f)
  res.mat[j+1] <- mean(stop.f)
  
  # ## repeat this process for the other interim analyses
  for (j in 2:(0.5*ncol(mat) - 0.5)){
    ## check success thresholds
    stop.p.temp <- mat[,2*j-1] >= logit(gam[j])
    stop.p <- ifelse(stop.f, 0, ifelse(stop.p, 1, stop.p.temp))
    res.mat[2*j-1] <- mean(stop.p)
    
    ## repeat process for failure thresholds
    stop.f.temp <- mat[,2*j] < logit(xi[j])
    stop.f <- ifelse(stop.p, 0, ifelse(stop.f, 1, stop.f.temp))
    res.mat[2*j] <- mean(stop.f)
  }
  ## implement process for final analysis (success thresholds only)
  j <- 0.5*ncol(mat) + 0.5
  stop.p.temp <- mat[,2*j-1] >= logit(gam[j])
  stop.p <- ifelse(stop.f, 0, ifelse(stop.p, 1, stop.p.temp))
  res.mat[2*j-1] <- mean(stop.p)
  
  return(res.mat)
}

## this function is used to create linear approximations to posterior and posterior predictive 
## probabilities under the conditional approach
getLines <- function(m1, m2, n1s, n2s){
  
  ## the inputs are described as follows:
  # m1 is matrix of sequential probs (first joint sampling distribution estimate)
  # m2 is the matrix of sequential probs (second joint sampling distribution estimate)
  # n1s are the sample size for the first joint sampling distribution estimate
  # n2s are the sample size for the second joint sampling distribution estimate
  
  # get logits for the sequential probs in both sampling distribution estimates
  ls <- m1
  li <- m2
  
  # adjust an infinite logits
  for (j in 1:ncol(ls)){
    ls[,j] <- ifelse(ls[,j] == -Inf, min(subset(ls, is.finite(ls[,j]))[,j]) - 1, ls[,j])
    ls[,j] <- ifelse(ls[,j] == Inf, max(subset(ls, is.finite(ls[,j]))[,j]) + 1, ls[,j])
  }
  
  for (j in 1:ncol(li)){
    li[,j] <- ifelse(li[,j] == -Inf, min(subset(li, is.finite(li[,j]))[,j]) - 1, li[,j])
    li[,j] <- ifelse(li[,j] == Inf, max(subset(li, is.finite(li[,j]))[,j]) + 1, li[,j])
  }
  
  # get indexes to combine individual logits later
  ls <- cbind(ls, seq(1, nrow(ls), 1))
  li <- cbind(li, seq(1, nrow(li), 1))
  
  slopes <- NULL
  ints <- NULL
  
  ## construct the slopes separately for each analysis
  for (j in 1:ncol(m1)){
    ls_s <- ls[order(ls[,j]),j]
    li_s <- li[order(li[,j]),j]
    
    l_slope <- (li_s - ls_s)/(n2s[j]-n1s[j])
    l_int <- ls_s - l_slope*n1s[j]
    
    # reorder according to smaller sample size
    l_slope[ls[order(ls[,j]),ncol(ls)]] <- l_slope 
    l_int[ls[order(ls[,j]),ncol(ls)]] <- l_int 
    
    slopes <- cbind(slopes, l_slope)
    ints <- cbind(ints, l_int)
  }
  
  return(cbind(ints, slopes))
}

## this function uses linear approximations to create a matrix of stopping probabilities across 
## a range of sample sizes that can be plotted or used to find optimal sample sizes
stopMat <- function(lines, c_vec, lb, ub, by, gam, xi){
  
  ## the inputs are described as follows:
  ## lines: a matrix of intercepts and slopes returned by getLines()
  ## c_vec: the vector dictating the spacing of the analyses
  ## lb: the smallest sample size (n1) for which the stopping probs are estimated
  ## ub: the largest sample size (n1) for which the stopping probs are estimated
  ## by: the increments between lb and ub at which the stopping probs are estimated
  ## gam: the vector of sucess thresholds for all analyses
  ## xi: the vector of failure thresholds for the first T-1 analyses
  
  ## get vector of sample sizes and ratios
  samps <- seq(lb,ub,by)
  ratios <- c_vec/c_vec[1]
  
  ## extract the number of decisions
  tps <- 0.5*ncol(lines)
  
  ## extract intercepts and slopes from the lines matrix
  res.mat <- matrix(0, ncol = tps, nrow = length(samps))
  ints <- lines[, seq(1, tps, 1)]
  slopes <- lines[, seq(tps + 1, 2*tps, 1)]
  
  ## repeat process for each sample size
  for (i in 1:length(samps)){
    
    ## check which logits are larger than first success threshold based on
    ## the linear approximations
    j <- 1
    stop.p <- ints[,j] + (samps[i]*ratios[j])^2*slopes[,j] >= logit(gam[j])
    res.mat[i, j] <- mean(stop.p)
    
    ## repeat process for failure thresholds at first analysis
    stop.f <- ints[,j + 1] + (samps[i]*ratios[j + 1])^1*slopes[,j + 1] < logit(xi[j])
    ## make sure we didn't just stop for success
    stop.f <- ifelse(stop.p, 0, stop.f)
    res.mat[i, j+1] <- mean(stop.f)
    
    # ## repeat this process for the other interim analyses
    for (j in 2:(0.5*length(c_vec) - 0.5)){
      # print(j)
      ## check success thresholds
      stop.p.temp <- ints[,2*j-1] + (samps[i]*ratios[2*j-1])^2*slopes[,2*j-1] >= logit(gam[j])
      stop.p <- ifelse(stop.f, 0, ifelse(stop.p, 1, stop.p.temp))
      res.mat[i, 2*j-1] <- mean(stop.p)
      
      ## repeat process for failure thresholds
      stop.f.temp <- ints[,2*j] + samps[i]*ratios[2*j]*slopes[,2*j] < logit(xi[j])
      stop.f <- ifelse(stop.p, 0, ifelse(stop.f, 1, stop.f.temp))
      res.mat[i, 2*j] <- mean(stop.f)
    }
    ## implement process for final analysis (success thresholds only)
    j <- 0.5*length(c_vec) + 0.5
    stop.p.temp <- ints[,2*j-1] + samps[i]*ratios[2*j-1]*slopes[,2*j-1] >= logit(gam[j])
    stop.p <- ifelse(stop.f, 0, ifelse(stop.p, 1, stop.p.temp))
    res.mat[i, 2*j-1] <- mean(stop.p)
    
  }
  
  return(res.mat)
}

## tune the decision thresholds using sampling distribution estimates under H0
## consider type I error
ns <- c(115, 235)
gamgam <- c(seq(0.97, 0.9, length.out = 15), 0.97)
xixi <- rep(0.05, 15)
for (j in 1:length(ns)){
  temp <- matrix(0, nrow = 1000, ncol = 31)
  for (k in 1:15){ ## get posterior predictive probabilities corresponding to gamma = 0.97
    temp[,2*k-1] <- read.csv(paste0("H0_suc_n_",ns[j],"_t_", k, ".csv"))[,6]
    temp[,2*k] <- read.csv(paste0("H0_fut_n_",ns[j],"_t_", k, ".csv"))[,6]
  }
  temp[,31] <- read.csv(paste0("H0_suc_n_",ns[j],"_t_16.csv"))[,2]
  assign(paste0("jointH0", ns[j]), temp)
}

## use these sampling distribution estimates to verify the type I error rate is controlled
simH0_alln <- NULL
for (j in 1:length(ns)){
  temp <- stopMatData(get(paste0("jointH0", ns[j])), gamgam, xixi)
  simH0_alln <- rbind(simH0_alln, temp)
}

## verify that the type I error rate is reasonably controlled
simH0_alln[,31]

## construct the sampling distribution estimates under H1
ns <- seq(100, 250, 15)
for (j in 1:length(ns)){
  temp <- matrix(0, nrow = 1000, ncol = 31)
  for (k in 1:15){ ## get posterior predictive probabilities corresponding to gamma = 0.97
    temp[,2*k-1] <- read.csv(paste0("suc_n_",ns[j],"_t_", k, ".csv"))[,6]
    temp[,2*k] <- read.csv(paste0("fut_n_",ns[j],"_t_", k, ".csv"))[,6]
  }
  temp[,31] <- read.csv(paste0("suc_n_",ns[j],"_t_16.csv"))[,2]
  assign(paste0("joint", ns[j]), temp)
}

## get the sample sizes for each column of the sampling distribution estimate matrix
nlow <- c(c(rbind((115*seq(1, 3.8, 0.2))^2, 115*seq(1, 3.8, 0.2))), 115*4)
nhigh <- c(c(rbind((235*seq(1, 3.8, 0.2))^2, 235*seq(1, 3.8, 0.2))), 235*4)

## get linear approximation and resulting stopping probabilities
lines.cond <- getLines(joint115, joint235, nlow, nhigh)
stop.cond <- stopMat(lines.cond, c(rep(seq(1, 3.8, 0.2), each = 2), 4), 
                     100, 250, 1, gamgam, xixi)

cbind(100:250, stop.cond[,31]) ## recommended sample size is 195 at first interim

## use simulated sampling distribution estimates to get stopping probabilities directly
sim_alln <- NULL
for (j in 1:length(ns)){
  temp <- stopMatData(get(paste0("joint", ns[j])), gamgam, xixi)
  sim_alln <- rbind(sim_alln, temp)
}

## construct efficacy plots
css <- seq(1, 4, 0.2)
suc_lb <- rep(c(0, 0.05, 0.1, 0.2), each = 4)
suc_ub <- rep(c(0.35, 0.6, 0.75, 1), each = 4)
## for loop for all efficacy analyses
for (k in 1:16){
  pwr.sim <- sim_alln[, 2*k-1]
  ## combine different power curve estimates into one data frame
  assign(paste0("dfsuc", k), data.frame(n = css[k]*c(n_vals, seq(head(n_vals, 1), tail(n_vals, 1), 1),
                                                     n_vals, n_vals),
                    power = c(pwr.sim, stop.cond[,2*k-1],
                              pwr.sim - 1.96*sqrt(pwr.sim*(1-pwr.sim)/1000),
                              pwr.sim + 1.96*sqrt(pwr.sim*(1-pwr.sim)/1000)),
                    curve = c(rep("C_Simulation", length(n_vals)),
                              rep("A_Algorithm 2", length(seq(head(n_vals, 1), tail(n_vals, 1), 1))),
                              rep("D_CI", length(n_vals)), rep("E_CI", length(n_vals)))))
  
  ## create subplot
  assign(paste0("plotsuc", k), ggplot(get(paste0("dfsuc", k)), aes(x=n)) + theme_bw() +
    geom_line(aes(y = power, color=as.factor(curve), linetype = as.factor(curve)), 
              alpha = 0.9, size = 1) +
    labs(title='') +
    labs(x= bquote(italic(n)[.(k)]), y= bquote('')) +
    theme(plot.title = element_text(size=18,face="bold", hjust =  0.5,
                                    margin = margin(t = 0, 0, 5, 0))) +
    theme(axis.text=element_text(size=16),
          axis.title=element_text(size=18)) +
    theme(legend.position="none") +
    scale_color_manual(name = " ", 
                       labels = c("Estimated  ", "Simulated  ", "95% CI  ",  "95% CI  "),
                       values = c("steelblue1", "firebrick", "firebrick", "firebrick")) +
    scale_linetype_manual(name = " ", 
                          labels = c("Estimated  ", "Simulated  ", "95% CI  ",  "95% CI  "),
                          values = c(1, 5, 3, 3)) +
    theme(legend.text=element_text(size=18)) +
    theme(legend.key.size = unit(1, "cm")) +
    ylim(suc_lb[k],suc_ub[k]) +
    theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
    theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))))
  
  if ((k%%4) != 1){
    assign(paste0("plotsuc", k), get(paste0("plotsuc", k)) +
             theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
             labs(x= bquote(italic(n)[.(k)]), y= ''))
  }
  
  if (k %in% c(13, 14, 15)){
    assign(paste0("plotsuc", k), get(paste0("plotsuc", k)) +
             scale_x_continuous(breaks = seq(400, 850, 150)))
  }
  
  if (k == 16){
    assign(paste0("plotsuc", k), get(paste0("plotsuc", k)) +
             geom_hline(yintercept=0.8, lty = 2))
  }
}

figp.row1 <- plot_grid(plotsuc1 + theme(plot.margin=unit(c(-0.25,0.3,0,-0.55),"cm")), 
                       plotsuc2 + theme(plot.margin=unit(c(-0.25,0.3,0,-0.45),"cm")),
                       plotsuc3 + theme(plot.margin=unit(c(-0.25,0.3,0,-0.45),"cm")),
                       plotsuc4 + theme(plot.margin=unit(c(-0.25,0.3,0,-0.45),"cm")),
                       rel_widths = c(1.11, 1, 1, 1), nrow = 1)
figp.row2 <- plot_grid(plotsuc5 + theme(plot.margin=unit(c(-0.25,0.3,0,-0.55),"cm")), 
                       plotsuc6 + theme(plot.margin=unit(c(-0.25,0.3,0,-0.45),"cm")),
                       plotsuc7 + theme(plot.margin=unit(c(-0.25,0.3,0,-0.45),"cm")),
                       plotsuc8 + theme(plot.margin=unit(c(-0.25,0.3,0,-0.45),"cm")),
                       rel_widths = c(1.11, 1, 1, 1), nrow = 1)
figp.row3 <- plot_grid(plotsuc9 + theme(plot.margin=unit(c(-0.25,0.3,0,-0.55),"cm")), 
                       plotsuc10 + theme(plot.margin=unit(c(-0.25,0.3,0,-0.45),"cm")),
                       plotsuc11 + theme(plot.margin=unit(c(-0.25,0.3,0,-0.45),"cm")),
                       plotsuc12 + theme(plot.margin=unit(c(-0.25,0.3,0,-0.45),"cm")),
                       rel_widths = c(1.11, 1, 1, 1), nrow = 1)
figp.row4 <- plot_grid(plotsuc13 + theme(plot.margin=unit(c(-0.25,0.3,0,-0.55),"cm")), 
                       plotsuc14 + theme(plot.margin=unit(c(-0.25,0.3,0,-0.45),"cm")),
                       plotsuc15 + theme(plot.margin=unit(c(-0.25,0.3,0,-0.45),"cm")),
                       plotsuc16 + theme(plot.margin=unit(c(-0.25,0.3,0,-0.45),"cm")),
                       rel_widths = c(1.11, 1, 1, 1), nrow = 1)

figp <- plot_grid(figp.row1, figp.row2, figp.row3, figp.row4, nrow = 4)

## create additional plot for legend
K <- 16
pwr.sim <- sim_alln[, 2*k-1]
## combine different power curve estimates into one data frame
assign(paste0("dfsuc", k), data.frame(n = css[k]*c(n_vals, seq(head(n_vals, 1), tail(n_vals, 1), 1),
                                                   n_vals),
                                      power = c(pwr.sim, stop.cond[,2*k-1],
                                                pwr.sim - 1.96*sqrt(pwr.sim*(1-pwr.sim)/1000)),
                                      curve = c(rep("C_Simulation", length(n_vals)),
                                                rep("A_Algorithm 2", length(seq(head(n_vals, 1), tail(n_vals, 1), 1))),
                                                rep("D_CI", length(n_vals)))))

## create subplot
assign(paste0("plotsuc.legend"), ggplot(get(paste0("dfsuc", k)), aes(x=n)) + theme_bw() +
         geom_line(aes(y = power, color=as.factor(curve), linetype = as.factor(curve)), 
                   alpha = 0.9, size = 1) +
         labs(title='') +
         labs(x= bquote(italic(n)[.(k)]), y= bquote('Cumulative Success Probability')) +
         theme(plot.title = element_text(size=18,face="bold", hjust =  0.5,
                                         margin = margin(t = 0, 0, 5, 0))) +
         theme(axis.text=element_text(size=16),
               axis.title=element_text(size=18)) +
         theme(legend.position="bottom") +
         scale_color_manual(name = " ", 
                            labels = c("Estimated  ", "Simulated  ", "95% CI  "),
                            values = c("steelblue1", "firebrick", "firebrick")) +
         scale_linetype_manual(name = " ", 
                               labels = c("Estimated  ", "Simulated  ",  "95% CI  "),
                               values = c(1, 5, 3)) +
         theme(legend.text=element_text(size=18)) +
         theme(legend.key.size = unit(1.5, "cm")) +
         ylim(suc_lb[k],suc_ub[k]) +
         theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
         theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))))

## add y-axis label
y_label <- textGrob("Cumulative Success Probability", rot = 90, gp = gpar(fontsize = 18))

fig.final.top <- plot_grid(y_label, figp, ncol = 2, rel_widths = c(0.5, 12))
fig.final <- plot_grid(fig.final.top, ggpubr::get_legend(plotsuc.legend), nrow = 2, rel_heights = c(12, 0.7))

# output as .pdf file for the article
pdf(file = "Fig_ORVAC_suc.pdf",   # The directory you want to save the file in
    width = 12, # The width of the plot in inches
    height = 12) # The height of the plot in inches

fig.final

dev.off()

## construct futility plots
css <- seq(1, 4, 0.2)
fut_lb <- rep(c(0, 0,0,0), each = 4)
fut_ub <- rep(c(0.15, 0.25, 0.3, 0.35), each = 4)
## for loop for all futility analyses
for (k in 1:15){
  pwr.sim <- sim_alln[, 2*k]
  ## combine different power curve estimates into one data frame
  assign(paste0("dffail", k), data.frame(n = css[k]*c(n_vals, seq(head(n_vals, 1), tail(n_vals, 1), 1),
                                                     n_vals, n_vals),
                                        power = c(pwr.sim, stop.cond[,2*k],
                                                  pwr.sim - 1.96*sqrt(pwr.sim*(1-pwr.sim)/1000),
                                                  pwr.sim + 1.96*sqrt(pwr.sim*(1-pwr.sim)/1000)),
                                        curve = c(rep("C_Simulation", length(n_vals)),
                                                  rep("A_Algorithm 2", length(seq(head(n_vals, 1), tail(n_vals, 1), 1))),
                                                  rep("D_CI", length(n_vals)), rep("E_CI", length(n_vals)))))
  
  ## create subplot
  assign(paste0("plotfail", k), ggplot(get(paste0("dffail", k)), aes(x=n)) + theme_bw() +
           geom_line(aes(y = power, color=as.factor(curve), linetype = as.factor(curve)), 
                     alpha = 0.9, size = 1) +
           labs(title='') +
           labs(x= bquote(italic(n)[.(k)]), y= bquote('')) +
           theme(plot.title = element_text(size=18,face="bold", hjust =  0.5,
                                           margin = margin(t = 0, 0, 5, 0))) +
           theme(axis.text=element_text(size=16),
                 axis.title=element_text(size=18)) +
           theme(legend.position="none") +
           scale_color_manual(name = " ", 
                              labels = c("Estimated  ", "Simulated  ", "95% CI  ",  "95% CI  "),
                              values = c("steelblue1", "firebrick", "firebrick", "firebrick")) +
           scale_linetype_manual(name = " ", 
                                 labels = c("Estimated  ", "Simulated  ", "95% CI  ",  "95% CI  "),
                                 values = c(1, 5, 3, 3)) +
           theme(legend.text=element_text(size=18)) +
           theme(legend.key.size = unit(1, "cm")) +
           ylim(fut_lb[k],fut_ub[k]) +
           theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
           theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))))
  
  if ((k%%4) != 1){
    assign(paste0("plotfail", k), get(paste0("plotfail", k)) +
             theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
             labs(x= bquote(italic(n)[.(k)]), y= ''))
  }
  
  if (k %in% c(13, 14, 15)){
    assign(paste0("plotfail", k), get(paste0("plotfail", k)) +
             scale_x_continuous(breaks = seq(400, 850, 150)))
  }
  
}

figp.row1 <- plot_grid(plotfail1 + theme(plot.margin=unit(c(-0.25,0.3,0,-0.55),"cm")), 
                       plotfail2 + theme(plot.margin=unit(c(-0.25,0.3,0,-0.45),"cm")),
                       plotfail3 + theme(plot.margin=unit(c(-0.25,0.3,0,-0.45),"cm")),
                       plotfail4 + theme(plot.margin=unit(c(-0.25,0.3,0,-0.45),"cm")),
                       rel_widths = c(1.11, 1, 1, 1), nrow = 1)
figp.row2 <- plot_grid(plotfail5 + theme(plot.margin=unit(c(-0.25,0.3,0,-0.55),"cm")), 
                       plotfail6 + theme(plot.margin=unit(c(-0.25,0.3,0,-0.45),"cm")),
                       plotfail7 + theme(plot.margin=unit(c(-0.25,0.3,0,-0.45),"cm")),
                       plotfail8 + theme(plot.margin=unit(c(-0.25,0.3,0,-0.45),"cm")),
                       rel_widths = c(1.11, 1, 1, 1), nrow = 1)
figp.row3 <- plot_grid(plotfail9 + theme(plot.margin=unit(c(-0.25,0.3,0,-0.55),"cm")), 
                       plotfail10 + theme(plot.margin=unit(c(-0.25,0.3,0,-0.45),"cm")),
                       plotfail11 + theme(plot.margin=unit(c(-0.25,0.3,0,-0.45),"cm")),
                       plotfail12 + theme(plot.margin=unit(c(-0.25,0.3,0,-0.45),"cm")),
                       rel_widths = c(1.11, 1, 1, 1), nrow = 1)
figp.row4 <- plot_grid(plotfail13 + theme(plot.margin=unit(c(-0.25,0.3,0,-0.55),"cm")), 
                       plotfail14 + theme(plot.margin=unit(c(-0.25,0.3,0,-0.45),"cm")),
                       plotfail15 + theme(plot.margin=unit(c(-0.25,0.3,0,-0.45),"cm")),
                       NULL,
                       rel_widths = c(1.11, 1, 1, 1), nrow = 1)

figp <- plot_grid(figp.row1, figp.row2, figp.row3, figp.row4, nrow = 4)

## add y-axis label
y_label <- textGrob("Cumulative Failure Probability", rot = 90, gp = gpar(fontsize = 18))

fig.final.top <- plot_grid(y_label, figp, ncol = 2, rel_widths = c(0.5, 12))
fig.final <- plot_grid(fig.final.top, ggpubr::get_legend(plotsuc.legend), nrow = 2, rel_heights = c(12, 0.7))

# output as .pdf file for the article
pdf(file = "Fig_ORVAC_fail.pdf",   # The directory you want to save the file in
    width = 12, # The width of the plot in inches
    height = 12) # The height of the plot in inches

fig.final

dev.off()