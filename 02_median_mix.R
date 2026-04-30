## this code file is used to reproduce the results from Section 6.2

## load necessary packages
require(foreach)
require(doParallel)
require(doSNOW)
require(numDeriv)
require(ggplot2)
require(ggpubr)
require(cowplot)
require(rjags)
require(coda)

## define helper functions
expit <- function(x){1/(1 + exp(-x))}
logit <- function(x){log(x) - log(1-x)}
l1 <- function(x){-1*logit(x)}

## set up parallelization
cores=detectCores()
cl <- makeSOCKcluster(71)

m <- 10000
registerDoSNOW(cl)
pb <- txtProgressBar(max = m, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

## first conduct simulations under conditional approach for H0
probs0 <- c(1/2, 1/4, 1/4)
alphas0 <- c(11.33153, 16.832127, 22.33243)
betas0 <- c(10, 15, 20)
thres <- 1.1

n.chains = 2 ## number of chains for MCMC
n.burnin = 1000 ## number of draws to discard at beginning of chain
n.draws = 10000 ## number of draws to retain per chain

## set array of sample sizes
n_vals <- seq(60, 200, 10)
for (j in 1:length(n_vals)){
  n <- n_vals[j]
  sim.res <- foreach(i=1:m, .packages=c('rjags', 'coda'), .combine=rbind,
                     .options.snow=opts) %dopar% {
                       
                       set.seed(j*m + i)
                       
                       ## get rjags model from text file
                       model.path <- paste(getwd(), '/jags_gbi.txt', sep='')
                       
                       ## simulate data from mixture gamma
                       t11 <- rgamma(n, alphas0[1], betas0[1])
                       t12 <- rgamma(n, alphas0[2], betas0[2])
                       t13 <- rgamma(n, alphas0[3], betas0[3])
                       
                       u <- runif(n)
                       t1 <- ifelse(u < probs0[1], t11, ifelse(u < sum(probs0[1:2]), t12, t13))
                       
                       ## process the data for rjags
                       n_jags <- length(t1)
                       y_jags <- t1
                       
                       ## estimate omega
                       dens.dat <- density(t1)
                       med.dat <- median(t1)
                      
                       omega.dat <- 2*mean(dnorm(med.dat, t1, dens.dat$bw))
                       
                       ## initialize model and burn-in
                       jags.m <- jags.model(model.path,
                                            data = list(n = n_jags, y = y_jags, zero = 0, omega = omega.dat), 
                                            n.chains=n.chains, quiet = TRUE)
                       
                       update(jags.m, n.burnin, progress.bar = "none")
                       
                       ## MCMC sampling
                       sim_gbi <- coda.samples(jags.m, c("theta"), n.iter=n.draws, progress.bar = "none")
                       
                       sim_gbi.all <- as.numeric(do.call(rbind, sim_gbi))
                       
                       ## get logit of posterior probability using kernel density estimate
                       lp <- log(sim_gbi.all)
                       lp <- ifelse(is.finite(lp), lp,
                                    ifelse(lp > 0, max(na.omit(lp))+1, min(na.omit(lp))-1))
                       kd.lp <- density(na.omit(lp))
                       pp <- mean(pnorm(log(thres), lp, kd.lp$bw, lower.tail = TRUE))
                       if (pp < 0.00004){
                         pp <- pnorm(log(thres), mean(lp), sd(lp), lower.tail = TRUE)
                       }
                       
                       pp
                       
                     }
  
  write.csv(sim.res, paste0("lmed_cond_n_", n, "_H0.csv"), row.names = FALSE)   
}

## repeat with predictive approach and the informative prior
sig_pct <- 0.075
pL <- pnorm(1.15, 1.2, sig_pct)
pU <- pnorm(1.25, 1.2, sig_pct)

set.seed(1)
n_vals <- seq(60, 200, 10)
for (j in 1:length(n_vals)){
  ## simulate true parameters for data generation according to the predictive approach
  n <- n_vals[j]
  pct_dec <- sort(qnorm(runif(m, pL, pU), 1.2, sig_pct))
  probs1 <- c(1/2, 1/4, 1/4)
  betas1 <- c(10, 15, 20)
  
  alphas11 <- 11.83161 + (pct_dec - 1.15)*(12.83175 - 11.83161)/0.1
  alphas12 <- 17.58218 + (pct_dec - 1.15)*(19.08228 - 17.58218)/0.1
  alphas13 <- 23.33247 + (pct_dec - 1.15)*(25.33255 - 23.33247)/0.1
  
  ## repeat process to estimate the sampling distribution under H1
  sim.res <- foreach(i=1:m, .packages=c('rjags', 'coda'), .combine=rbind,
                     .options.snow=opts) %dopar% {
                       
                       set.seed(j*m + i)
                       
                       model.path <- paste(getwd(), '/jags_gbi.txt', sep='')
                       
                       ## simulate data
                       t11 <- rgamma(n, alphas11[i], betas1[1])
                       t12 <- rgamma(n, alphas12[i], betas1[2])
                       t13 <- rgamma(n, alphas13[i], betas1[3])
                       
                       u <- runif(n)
                       t1 <- ifelse(u < probs1[1], t11, ifelse(u < sum(probs1[1:2]), t12, t13))
                       
                       ## process the data for rjags
                       n_jags <- length(t1)
                       y_jags <- t1
                       
                       dens.dat <- density(t1)
                       med.dat <- median(t1)
                       
                       omega.dat <- 2*mean(dnorm(med.dat, t1, dens.dat$bw))
                       
                       jags.m <- jags.model(model.path,
                                            data = list(n = n_jags, y = y_jags, zero = 0, omega = omega.dat), 
                                            n.chains=n.chains, quiet = TRUE)
                       
                       update(jags.m, n.burnin, progress.bar = "none")
                       
                       sim_gbi <- coda.samples(jags.m, c("theta"), n.iter=n.draws, progress.bar = "none")
                       
                       sim_gbi.all <- as.numeric(do.call(rbind, sim_gbi))
                       
                       lp <- log(sim_gbi.all)
                       lp <- ifelse(is.finite(lp), lp,
                                    ifelse(lp > 0, max(na.omit(lp))+1, min(na.omit(lp))-1))
                       kd.lp <- density(na.omit(lp))
                       pp <- mean(pnorm(log(thres), lp, kd.lp$bw, lower.tail = TRUE))
                       if (pp < 0.00004){
                         pp <- pnorm(log(thres), mean(lp), sd(lp), lower.tail = TRUE)
                       }
                       
                       pp
                       
                     }
  
  write.csv(sim.res, paste0("lmed_pred_n_", n, "_H1.csv"), row.names = FALSE)   
  
}


## use this function to estimate the type I error rate (or power)
## under the conditional approach; the inputs are as follows:
## m1: matrix of complementary probabilities at first sample size (from .csv files)
## m2: matrix of complementary probabilities at second sample size (from .csv files)
## n0: first sample size
## n1: second sample size
## lb: lower bound for which to compute OCs
## ub: upper bound for which to compute OCs
## by: increments between lower and upper bound
## gam: decision threshold for success
pwr_one_lin <- function(m1, m2, n0, n1, lb, ub, by, gam){
  
  ## get logits for the posterior probabilities
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
  
  ## create matrix to calculate power or type I error rate
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

## use this to estimate power under the predictive approach;
## the inputs are as follows:
## m1: matrix of logits at first sample size (from .csv files)
## m2: matrix of logits at second sample size (from .csv files)
## n0: first sample size
## n1: second sample size
## lb: lower bound for which to compute OCs
## ub: upper bound for which to compute OCs
## by: increments between lower and upper bound
## gam: decision threshold for success
## M: number of bins for logits at each sample size
pwr_one_pred <- function(m1, m2, n0, n1, lb, ub, by, gam, M = 10){
  
  m <- length(m1)/M
  
  ## get logits for the posterior probs
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
  
  ## split out by bins to construct linear approximations
  for (i in 1:M){
    ls_i <- ls[seq((i-1)*m + 1, i*m, 1)]
    li_i <- li[seq((i-1)*m + 1, i*m, 1)]
    
    ## construct the slopes using the order statistics of the sampling
    ## distribution estimates
    ls_s <- ls_i[order(ls_i)]
    li_s <- li_i[order(li_i)]
    
    slopes_i <- (li_s - ls_s)/(n1-n0)
    ints_i <- ls_s - slopes_i*n0
    
    slopes <- c(slopes, slopes_i)
    ints <- c(ints, ints_i)
  }
  
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

## success threshold
gam <- 0.0575

## construct the figure for the paper
pwr.lin <- pwr_one_lin(read.csv("lmed_cond_n_180_H0.csv")[,1],
                       read.csv("lmed_cond_n_80_H0.csv")[,1],
                       180, 80, 60, 200, 1, 1 - gam)

pwr.sim <- NULL
for (i in 1:length(seq(60, 200, 10))){
  pwr.temp <- mean(1 - read.csv(paste0("lmed_cond_n_",seq(60, 200, 10)[i],"_H0.csv"))[,1] >= 1 - gam)
  pwr.sim <- c(pwr.sim, pwr.temp)
}

## combine different power curve estimates into one data frame
df1 <- data.frame(n = c(n_vals, seq(60, 200, 1)),
                  power = c(pwr.sim, pwr.lin[,2]),
                  curve = c(rep("C_Simulation", length(n_vals)),
                            rep("A_Algorithm 2", length(seq(60, 200, 1)))))

## create subplot
plot1 <- ggplot(df1, aes(x=n)) + theme_bw() +
  geom_line(aes(y = power, color=as.factor(curve), linetype = as.factor(curve)), 
            alpha = 0.9, size = 1) +
  labs(title=bquote('Type I Error Rate')) +
  labs(x= bquote(italic(n)), y= bquote('Success Probability')) +
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
  ylim(0,0.1) + 
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  geom_hline(yintercept=0.05, lty = 2) +
  scale_x_continuous(breaks = seq(80, 200, 40))

pwr.lin <- pwr_one_pred(read.csv("lmed_pred_n_180_H1.csv")[,1],
                        read.csv("lmed_pred_n_80_H1.csv")[,1],
                        180, 80, 60, 200, 1, 1 - gam)

pwr.sim <- NULL
for (i in 1:length(seq(60, 200, 10))){
  pwr.temp <- mean(1 - read.csv(paste0("lmed_pred_n_",seq(60, 200, 10)[i],"_H1.csv"))[,1] >= 1 - gam)
  pwr.sim <- c(pwr.sim, pwr.temp)
}

## combine different power curve estimates into one data frame
df2 <- data.frame(n = c(n_vals, seq(60, 200, 1)),
                  power = c(pwr.sim, pwr.lin[,2]),
                  curve = c(rep("C_Simulation", length(n_vals)),
                            rep("A_Algorithm 2", length(seq(60, 200, 1)))))

## create subplot
plot2 <- ggplot(df2, aes(x=n)) + theme_bw() +
  geom_line(aes(y = power, color=as.factor(curve), linetype = as.factor(curve)), 
            alpha = 0.9, size = 1) +
  labs(title=bquote('Power')) +
  labs(x= bquote(italic(n)), y= bquote('Success Probability')) +
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
  ylim(0.6,1) + 
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  geom_hline(yintercept=0.9, lty = 2)  +
  scale_x_continuous(breaks = seq(80, 200, 40))

## arrange for the final plot
figp.row1 <- plot_grid(plot1 + theme(plot.margin=unit(c(0.25,0.75,0.25,0.75),"cm")), 
                       plot2 + theme(plot.margin=unit(c(0.25,0.75,0.25,0.75),"cm")),
                       rel_widths = c(1, 1))

## get a common legend for the larger figure
plot1.legend <- ggplot(df2, aes(x=n)) + theme_bw() +
  geom_line(aes(y = power, color=as.factor(curve), linetype = as.factor(curve)), 
            alpha = 0.9, size = 1) +
  labs(title=bquote('Power')) +
  labs(x= bquote(italic(n)), y= bquote('Success Probability')) +
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
  ylim(0.6,1) + 
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  geom_hline(yintercept=0.9, lty = 2) +
  theme(legend.key.size = unit(1.5, "cm"))

## add legend to bottom of the plot
fig_final <- plot_grid(figp.row1, ggpubr::get_legend(plot1.legend), ncol = 1, rel_heights = c(1, .1))

# output as .pdf file for the article
pdf(file = "Figure_Median.pdf",   # The directory you want to save the file in
    width = 9.5, # The width of the plot in inches
    height = 4.5) # The height of the plot in inches

fig_final

dev.off()