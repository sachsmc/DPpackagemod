n <- 20
freq <- TRUE

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
             ndisplay=1e6)

xpred <- indi.results[, 2]
xpred.samp <- unlist(lapply(1:nrow(indi.results), function(i) rnorm(100, mean = indi.results[i, 2], sd = indi.results[i, 4])))
dex.list <- lapply(1:nrow(indi.results), function(i) ((i - 1) * 100) + 1:100)

fitLDDP <- DPcdensity(y = y, x = x1,mus=indi.results[,1:2],sds=indi.results[,3:4], xpred = xpred,
                      ngrid=50, 
                      compute.band=FALSE,
                      type.band="HPD",
                      prior=prior, 
                      mcmc=mcmc, 
                      state=NULL, 
                      status=TRUE)


fitLDDP.2 <- DPcdensity(y = y, x = x1,mus=indi.results[,1:2],sds=indi.results[,3:4], xpred = xpred.samp,
                      ngrid=50, 
                      compute.band=FALSE,
                      type.band="HPD",
                      prior=prior, 
                      mcmc=mcmc, 
                      state=NULL, 
                      status=TRUE)

plot(sapply(dex.list, function(ii) mean(xpred.samp[ii])) ~ xpred)

ave.density <- lapply(dex.list, function(ii) colMeans(fitLDDP.2$densp.m[ii, ]))
hat.density <- lapply(1:20, function(i) fitLDDP$densp.m[i, ])
# 1.) just using the sample from 3 
# as the distrbution for D 
# or the true value in D 
# 2.) getting the denisty at Tsj hat 
# yeah?
# 3.) True Tcjs
samp.from.density <- function(n, x, y) {
  
  sample(x, size = n, prob = y, replace = TRUE)
  
}

Tc1 <- rnorm(50, mean = indi.results[, 1], sd = indi.results[, 3])
Tc2 <- fitLDDP$densp.m
Tc3 <- Y

Tc0 <- ## doesn't exist yet