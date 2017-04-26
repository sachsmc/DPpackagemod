library(ggplot2)
library(DPpackage)

x1 <- c(rnorm(15, sd = 2), 2 + rt(15, df = 5))
y <- 2 + cos(2 * x1 * pi / 5) + .25 * x1 + .25 * rt(length(x1), df = 15)

n <- 200
X <- rnorm(n, mean = 0)
x1 <- rnorm(n, mean = X, sd = 0.5) ## observed

y <- 1 + 2.5 * X + rnorm(n)

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
xpred <- seq(-2, 2, .1)	
################################################
# fiiting the model
################################################

fitLDDP <- DPcdensity(y = y, x = x1,muy=y,mux=x1*.8, xpred = xpred,
                       ngrid=1, 
                       compute.band=TRUE,
                       type.band="HPD",
                       prior=prior, 
                       mcmc=mcmc, 
                       state=NULL, 
                       status=TRUE)

plot(y ~ X)
lines(y ~ x1, type = "p", pch = 20)
lines(xpred, fitLDDP$meanfp.m, type = "l", lwd = 2, col = "red")
lines(xpred, fitLDDP$meanfp.l, type = "l", lwd = 1, lty = 2, col = "red")
lines(xpred, fitLDDP$meanfp.h, type = "l", lwd = 1, lty = 2, col = "red")



library(rjags)


## full measurement error model
model_string <- "model{

# Likelihood
for(i in 1:n){
Y[i]   ~ dnorm(mu[i],inv.var)
mu[i] <- beta[1] + beta[2]*X[i]
X[i] ~ dnorm(jqe * X2[i], xvar)
#X2[i] ~ dnorm(mux, varx)
}

jqe <- .5
xvar <- pow(.5, -1)

# Prior for beta
for(j in 1:2){
beta[j] ~ dnorm(0,0.0001)
}

# Prior for the inverse variance
inv.var   ~ dgamma(0.01, 0.01)
mux ~ dnorm(0.0, 0.001)
varx ~ dgamma(0.01, 0.01)
# inv.det   ~ dgamma(0.01, 0.01)
sigma     <- 1/sqrt(inv.var)

}"


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

abline(res.full[-c(1:n)][1], res.full[-c(1:n)][2], col = "blue")

hsy <- lm(y ~ x1)
abline(hsy$coefficients[1], hsy$coefficients[2], col = "green")

abline(1, 2.5, lwd = 2, col = "purple")

plot(X ~ res.full[1:n])
abline(0, 1)
