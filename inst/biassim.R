library(DPpackage)
library(arm)
library(splines)
library(rjags)


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
  xpred <- seq(-7, 7, .1)	
  ################################################
  # fiiting the model
  ################################################
  
  truy <- 2 + .25 * xpred
  
  fitLDDP <- DPcdensity(y = y, x = x1,mus=indi.results[,1:2],sds=indi.results[,3:4], xpred = xpred,
                        ngrid=1, 
                        compute.band=TRUE,
                        type.band="HPD",
                        prior=prior, 
                        mcmc=mcmc, 
                        state=NULL, 
                        status=TRUE)
  
  msedpp <- mean((fitLDDP$meanfp.m - truy)^2)
  
  
  mselm <- mean((predict(lm(y ~ bs(x1, df = 3)), newdata = data.frame(x1 = xpred)) - truy)^2)
  
  
  
  
  ## full measurement error model
  model_string_old <- "model{

# Likelihood
for(i in 1:n){
X2[i] ~ dnorm(X[i], 1.0)
Y[i]   ~ dnorm(mu[i],inv.var)
mu[i] <- beta[1] + beta[2]*X[i]
X[i] ~ dnorm(mux, varx)

}

# Prior for beta
for(j in 1:2){
beta[j] ~ dnorm(0,0.0001)
}

# Prior for the inverse variance
inv.var   ~ dgamma(0.01, 0.01)
mux ~ dnorm(0.0, 0.001)
varx ~ dgamma(0.01, 0.01)
#inv.det   ~ dgamma(0.01, 0.01)
sigma     <- 1/sqrt(inv.var)

}"


model <- jags.model(textConnection(model_string_old), 
                    data = list(Y = y, n = length(y), X2 = x1), quiet = TRUE)


update(model, 1000, progress.bar = "none")

samp0 <- coda.samples(
  model,
  variable.names = c("beta", "X"),
  n.iter = 1000, thin = 2,
  progress.bar = "none"
)

res.full <- colMeans(samp0[[1]])

msejags <- mean((res.full[-c(1:n)][1] + xpred * res.full[-c(1:n)][2] - truy)^2)

c(msedpp, mselm, msejags)

}

bayesglm <- replicate(n = 10, sim.bias(20, freq = FALSE))

library(ggplot2)
ggplot(data.frame(mse = c(bayesglm[1,], bayesglm[2,], bayesglm[3,]), type = c(rep("dpp", 10), rep("lin", 10), rep("jags", 10))), 
       aes(x = type, y = mse))  + geom_boxplot()+ geom_jitter()


freqglm <- replicate(n = 10, sim.bias(20, freq = TRUE))
ggplot(data.frame(mse = c(bayesglm[1,], bayesglm[2,], bayesglm[3,], freqglm[1, ]), 
                  type = c(rep("dpp", 10), rep("lin", 10), rep("jags", 10), rep("dpp.freq", 10))), 
       aes(x = type, y = mse))  + geom_boxplot()+ geom_jitter()

