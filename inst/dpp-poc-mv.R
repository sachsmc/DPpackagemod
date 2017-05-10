library(DPpackagemod)
library(arm)
library(splines)


n <- 350

X <- rnorm(n, sd = 3)
X2 <- rt(n, df = 2) + .25 * X
Y <-  2 + cos(2 * X * pi / 5) + .25 * X2 + .25 * rt(length(X), df = 15)
enns <- sample(ceiling(runif(n, 100, 500)))

samp.indi.data <- function(i) {
  
  Z <- rbinom(enns[i], 1, .5)
  mu.S <- 2.5 + X[i] * Z
  mu.S2 <- .5 + X2[i] * Z
  lambda.T <- -1 + Y[i] * Z
  
  S <- rnorm(enns[i], mean = mu.S)
  S2 <- rnorm(enns[i], mean = mu.S2)
  T <- rpois(enns[i], exp(lambda.T))
  data.frame(Z, S, S2, T)
  
  
}

analyze.indi.data <- function(test) {
  
  s1 <- bayesglm(S ~ Z, data = test)
  s2 <- bayesglm(S2 ~ Z, data = test)
  t1 <- bayesglm(T ~ Z, data = test, family = "poisson")
  spost <- sim(s1, n.sims = 1000)@coef[, 2]
  s2post <- sim(s2, n.sims = 1000)@coef[, 2]
  tpost <- sim(t1, n.sims = 1000)@coef[, 2]
  
  c(mean(tpost),mean(spost), mean(s2post), sd(tpost),  sd(spost), sd(s2post))
  
}


indi.results <- do.call(rbind, lapply(lapply(1:n, samp.indi.data), analyze.indi.data))

y <- indi.results[, 1]
x1 <- indi.results[, 2]
x2 <- indi.results[, 3]

w <- cbind(y, x1, x2)#, x2, x3)
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
xpred <- cbind(X, X2)
################################################
# fiiting the model
################################################



fitLDDP <- DPcdensity(y = y, x = cbind(x1, x2),mus=indi.results[,1:3],sds=indi.results[,4:6], xpred = xpred,
                      ngrid=1, 
                      compute.band=TRUE,
                      type.band="HPD",
                      prior=prior, 
                      mcmc=mcmc, 
                      state=NULL, 
                      status=TRUE)

plot(Y, fitLDDP$meanfp.m, ylim = c(-1, 6), xlim = c(-1, 6))

abline(0, 1)

s1 <- lm(y ~ bs(x1, df = 3) + bs(x2, df = 3))
lines(Y ~ predict(s1), type = "p", col = "purple", pch = 20)

res <- data.frame(fitted = c(fitLDDP$meanfp.m, predict(s1)), true = c(Y, Y), type = c(rep("DPP", length(Y)), rep("spline", length(Y))))
library(ggplot2)

ggplot(res, aes(x = fitted, y = true, color = type)) + geom_point() + stat_smooth(method = "lm", se = FALSE) + geom_abline(intercept = 0, slope = 1)
