library(DPpackage)
library(arm)
library(splines)
library(rjags)

X <- 1
Y <-  2 + .25 * X + .25 * rt(50000, df = 15)

n <- 50

sim.bias <- function(n, freq = TRUE) {
  
  X <- rnorm(n, sd = 3)
  Y <-  2 + .25 * X + .25 * rt(length(X), df = 15)
  enns <- sample(ceiling(runif(n, 100, 500)))
  
  samp.indi.data <- function(i) {
    
    Z <- rbinom(enns[i], 1, .5)
    mu.S <- 2.5 + X[i] * Z
    lambda.T <- -1 + Y[i] * Z
    
    S <- rnorm(enns[i], mean = mu.S)
    T <- rpois(enns[i], exp(lambda.T))
    data.frame(Z, S, T)
    
    
  }
  
  analyze.indi.data <- function(test) {
    
    s1 <- bayesglm(S ~ Z, data = test)
    t1 <- bayesglm(T ~ Z, data = test, family = "poisson")
    spost <- sim(s1, n.sims = 1000)@coef[, 2]
    tpost <- sim(t1, n.sims = 1000)@coef[, 2]
    
    c(mean(tpost),mean(spost), sd(tpost),  sd(spost))
    
  }
  
  analyze.indi.data.freq <- function(test) {
    
    s1 <- glm(S ~ Z, data = test)
    t1 <- glm(T ~ Z, data = test, family = "poisson")
    spost <- s1$coefficients[2]
    tpost <- t1$coefficients[2]
    
    c(tpost,spost, sqrt(diag(vcov(t1)))[2],  sqrt(diag(vcov(s1)))[2])
    
  }
  
  if(freq) myanalyze <- analyze.indi.data.freq else myanalyze <- analyze.indi.data
  
  indi.results <- do.call(rbind, lapply(lapply(1:n, samp.indi.data), myanalyze))
  
  y <- indi.results[, 1]
  x1 <- indi.results[, 2]
  
  chidat <- data.frame(x1 = x1, y = y)
  #ggplot(chidat, aes(x = x1, y = y)) + geom_point() + stat_smooth()
  
  w <- cbind(y, x1)#, x2, x3)
  wbar <- apply(w, 2, mean) 
  wcov <- var(w) 
  prior <- list(a0 = 10, 
                b0 = 1, 
                nu1 = 4, 
                nu2 = 4, 
                s2 = 0.5 * wcov, 
                m2 = wbar, 
                psiinv2 = 2 * solve(wcov),
                tau1 = 6.01, tau2 = 3.01)
  
  mcmc <- list(nburn=5000,
               nsave=1000,
               nskip=2,
               ndisplay=2000)
  
  ################################################
  # covariate values where the density
  # and mean function is evaluated
  ################################################
  xpred <- c(-1, 0, 1)
  ################################################
  # fiiting the model
  ################################################
  
  truy <- 2 + .25 * xpred
  
  
  fitLDDP <- DPcdensity(y = y, x = x1,mus=indi.results[,1:2],sds=indi.results[,3:4], xpred = xpred,
                        ngrid=500, 
                        compute.band=FALSE,
                        type.band="HPD",
                        prior=prior, 
                        mcmc=mcmc, 
                        state=NULL, 
                        status=TRUE)
  
  
  msedpp <- mean((fitLDDP$meanfp.m - truy)^2)
  
  get.variance <- function(fitLDDP) {
    mu <- fitLDDP$meanfp.m
    
    i <- 1
    mean(((fitLDDP$grid - mu[i])^2) * fitLDDP$densp.m[i, ])
    
  }
  
}

bayesglm <- replicate(n = 10, sim.bias(20, freq = FALSE))

library(ggplot2)
ggplot(data.frame(mse = c(bayesglm[1,], bayesglm[2,], bayesglm[3,]), type = c(rep("dpp", 10), rep("lin", 10), rep("jags", 10))), 
       aes(x = type, y = mse))  + geom_boxplot()+ geom_jitter()


freqglm <- replicate(n = 10, sim.bias(20, freq = TRUE))
ggplot(data.frame(mse = c(bayesglm[1,], bayesglm[2,], bayesglm[3,], freqglm[1, ]), 
                  type = c(rep("dpp", 10), rep("lin", 10), rep("jags", 10), rep("dpp.freq", 10))), 
       aes(x = type, y = mse))  + geom_boxplot()+ geom_jitter()

