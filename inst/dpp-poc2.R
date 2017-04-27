library(DPpackage)
library(arm)
library(splines)


n <- 50

X <- rnorm(n, sd = 3)
Y <-  2 + cos(2 * X * pi / 5) + .25 * X + .25 * rt(length(X), df = 15)
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


indi.results <- do.call(rbind, lapply(lapply(1:n, samp.indi.data), analyze.indi.data))

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
             ndisplay=500)

################################################
# covariate values where the density
# and mean function is evaluated
################################################
xpred <- seq(min(X), max(X), .1)	
################################################
# fiiting the model
################################################


plot(Y ~ X)
lines(y ~ x1, type = "p", pch = 20)

fitLDDP <- DPcdensity(y = y, x = x1,mus=indi.results[,1:2],sds=indi.results[,3:4], xpred = xpred,
                      ngrid=1, 
                      compute.band=TRUE,
                      type.band="HPD",
                      prior=prior, 
                      mcmc=mcmc, 
                      state=NULL, 
                      status=TRUE)

lines(xpred, fitLDDP$meanfp.m, type = "l", lwd = 2, col = "red")
lines(xpred, fitLDDP$meanfp.l, type = "l", lwd = 1, lty = 2, col = "red")
lines(xpred, fitLDDP$meanfp.h, type = "l", lwd = 1, lty = 2, col = "red")

curve(2 + cos(2 * x * pi / 5) + .25 * x, add = TRUE, col = "blue")
hsy <- lm(Y ~ X)
abline(hsy$coefficients[1], hsy$coefficients[2], col = "green")

s1 <- lm(y ~ bs(x1, df = 3))
lines(predict(s1, newdata = data.frame(x1 = xpred)) ~ xpred, col = "purple")

library(rjags)



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
  n.iter = 50000, thin = 10,
  progress.bar = "none"
)

res.full <- colMeans(samp0[[1]])

abline(res.full[-c(1:n)][1], res.full[-c(1:n)][2], col = "purple")


