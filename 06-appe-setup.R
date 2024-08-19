## load required libraries
require(ggplot2)
require(cowplot)
require(rjags)
require(coda)
require(MASS)

## get gamma parameters for vldv using summary statistics
## pre-study  
a1 <- (1/.458)^2; b1 <- a1/24.5
a2 <- (1/.465)^2; b2 <- a2/24.9

## after 68 weeks
c1 <- (1/.458)^2; d1 <- c1/(0.78*24.5)
c2 <- (1/.465)^2; d2 <- c2/(0.93*24.9)

## generate from Gaussian copula for group A
set.seed(1)
samp1 <- mvrnorm(1306, mu = c(0,0), Sigma = rbind(c(1, 0.8), c(0.8, 1)))
u <- pnorm(samp1)
## obtain pre and post measurements with CDF inversion
pre1 <- qgamma(u[,1], a1, b1)
post1 <- qgamma(u[,2], c1, d1)
## obtain ratios
plot(density(post1/pre1))

## repeat process for group B
set.seed(2)
samp2 <- mvrnorm(655, mu = c(0,0), Sigma = rbind(c(1, 0.85), c(0.85, 1)))
v <- pnorm(samp2)
pre2 <- qgamma(v[,1], a2, b2)
post2 <- qgamma(v[,2], c2, d2)
plot(density(post2/pre2))

## save data for group A (w1) and group B (w2)
w1 <- post1/pre1
w2 <- post2/pre2

## use JAGS to obtain posterior draws via MCMC
set.seed(6)
burnin <- 1000
nchains <- 1
nthin <- 2
ndraws <- nthin*100000
n1 <- length(w1)
mu0 <- 2; tau0 <- 1; kappa0 <- 2; lambda0 <- 1
model1.fit <- jags.model(file="JAGS_weibull.txt",
                         data=list(n=n1, y = w1, 
                                   tau0 = tau0, mu0 = mu0,
                                   kappa0 = kappa0, lambda0 = lambda0), 
                         n.chains = nchains)

update(model1.fit, burnin)
model1.samples <- coda.samples(model1.fit, c("nu", "l"), n.iter=ndraws, thin=nthin)

## save posterior draws for group A
nu.1 <- unlist(model1.samples[,2])
lambda.1 <- unlist(model1.samples[,1])

## repeat process for group B
set.seed(7)
n2 <- length(w2)
model2.fit <- jags.model(file="JAGS_weibull.txt",
                         data=list(n=n2, y = w2, 
                                   tau0 = tau0, mu0 = mu0,
                                   kappa0 = kappa0, lambda0 = lambda0), 
                         n.chains = nchains)

update(model2.fit, burnin)
model2.samples <- coda.samples(model2.fit, c("nu", "l"), n.iter=ndraws, thin=nthin)

## save posterior draws for group B
nu.2 <- unlist(model2.samples[,2])
lambda.2 <- unlist(model2.samples[,1])

## find posterior means for design values (nu is shape a, lambda is scale b)
mean(nu.1) ## should be 3.09
mean(lambda.1) ## should be 0.90
mean(nu.2) ## should be 3.76
mean(lambda.2) ## should be 1.07

## save posterior draws to .csv files for reference
write.csv(nu.1, "nu1s.csv", row.names = FALSE)
write.csv(lambda.1, "lambda1s.csv", row.names = FALSE)
write.csv(nu.2, "nu2s.csv", row.names = FALSE)
write.csv(lambda.2, "lambda2s.csv", row.names = FALSE)

## read posterior draws from .csv files if necessary
alpha.1 <- read.csv("nu1s.csv")$x
beta.1 <- read.csv("lambda1s.csv")$x
alpha.2 <- read.csv("nu2s.csv")$x
beta.2 <- read.csv("lambda2s.csv")$x

## fit informative priors
## define factor by which to inflate the variance (25 here)
var_inf <- 25

## consider alpha1 = a1
## find the posterior mode and variance
a1dens <- density(alpha.1, n = 1024, bw = "SJ")
a1_mode <- a1dens$x[which.max(a1dens$y)]
a1_var <- var(alpha.1)

## solve system of equations to get desired gamma distribution
b1_star <- (-a1_mode - sqrt(a1_mode^2 + 4*var_inf*a1_var))/(-2*var_inf*a1_var); b1_star
a1_star <- 1 + a1_mode*b1_star; a1_star
informs <- c(a1_star, b1_star)

## consider beta1 = b1
## find the posterior mode and variance
a1dens <- density(beta.1, n = 1024, bw = "SJ")
a1_mode <- a1dens$x[which.max(a1dens$y)]
a1_var <- var(beta.1)

## solve system of equations to get desired gamma distribution
b1_star <- (-a1_mode - sqrt(a1_mode^2 + 4*var_inf*a1_var))/(-2*var_inf*a1_var); b1_star
a1_star <- 1 + a1_mode*b1_star; a1_star
informs <- rbind(informs, c(a1_star, b1_star))

## consider alpha2 = a2
## find the posterior mode and variance
a1dens <- density(alpha.2, n = 1024, bw = "SJ")
a1_mode <- a1dens$x[which.max(a1dens$y)]
a1_var <- var(alpha.2)

## solve system of equations to get desired gamma distribution
b1_star <- (-a1_mode - sqrt(a1_mode^2 + 4*var_inf*a1_var))/(-2*var_inf*a1_var); b1_star
a1_star <- 1 + a1_mode*b1_star; a1_star
informs <- rbind(informs, c(a1_star, b1_star))

## consider beta2 = b2
## find the posterior mode and variance
a1dens <- density(beta.2, n = 1024, bw = "SJ")
a1_mode <- a1dens$x[which.max(a1dens$y)]
a1_var <- var(beta.2)

## solve system of equations to get desired gamma distribution
b1_star <- (-a1_mode - sqrt(a1_mode^2 + 4*var_inf*a1_var))/(-2*var_inf*a1_var); b1_star
a1_star <- 1 + a1_mode*b1_star; a1_star
informs <- rbind(informs, c(a1_star, b1_star))

## output informative gamma distribution hyperparameters to .csv file
write.csv(informs, "informs_weibull.csv", row.names = FALSE)

## function to do MCMC for the posterior of theta1/theta2
WeibullPost <- function(y1, y2, mu1, tau1, kappa1, lambda1,
                        mu2, tau2, kappa2, lambda2, burnin = 1000, nchains = 1,
                        nthin = 1, ndraws = 50000){
  
  if (length(y1) <= 100){
    ndraws <- 75000
  }
  else if (length(y1) >= 200){
    ndraws <- 25000
  }
  
  ## y1 and y2 are the observations in each group (w in paper)
  ## nu_j has a Gamma(mu_j, tau_j) prior, where tau_j is a rate (a in paper)
  ## lambda_j has a Gamma(kappa_j, lambda_j) prior, where lambda_j is a rate (b in paper)
  ## tau is the threshold for the tail probability
  ## burnin is the number of MCMC iterations to discard at the start of each chain
  ## nchains is the number of chains to generate
  ## nthin is the thinning parameter for the MCMC process
  ## ndraws is the number of draws to generate (excluding burnin but including thinned draws)
  
  n1 <- length(y1)
  model1.fit <- jags.model(file="JAGS_weibull.txt",
                           data=list(n=n1, y = y1, 
                                     tau0 = tau1, mu0 = mu1,
                                     kappa0 = kappa1, lambda0 = lambda1), 
                           n.chains = nchains, quiet = TRUE)
  
  update(model1.fit, burnin, progress.bar = "none")
  model1.samples <- coda.samples(model1.fit, c("nu", "l"), n.iter=ndraws, thin=nthin, progress.bar = "none")
  
  nu.1 <- unlist(model1.samples[,2])
  lambda.1 <- unlist(model1.samples[,1])
  
  n2 <- length(y2)
  model2.fit <- jags.model(file="JAGS_weibull.txt",
                           data=list(n=n2, y = y2, 
                                     tau0 = tau2, mu0 = mu2,
                                     kappa0 = kappa2, lambda0 = lambda2), 
                           n.chains = nchains, quiet = TRUE)
  
  update(model2.fit, burnin, progress.bar = "none")
  model2.samples <- coda.samples(model2.fit, c("nu", "l"), n.iter=ndraws, thin=nthin, progress.bar = "none")
  
  nu.2 <- unlist(model2.samples[,2])
  lambda.2 <- unlist(model2.samples[,1])
  
  theta1 <- lambda.1*gamma(1 + 1/nu.1)
  theta2 <- lambda.2*gamma(1 + 1/nu.2)
  theta <- theta1/theta2
  theta <- ifelse(is.na(theta), Inf, theta)
  return(theta)
}