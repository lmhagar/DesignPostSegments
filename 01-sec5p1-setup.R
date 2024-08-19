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
set.seed(3)
burnin <- 1000
nchains <- 1
nthin <- 2
ndraws <- nthin*100000
sum_y1 <- sum(w1)
sum_logy1 <- sum(log(w1))
n1 <- length(w1)
mu0 <- 2; tau0 <- 0.25; kappa0 <- 2; lambda0 <- 0.25
model1.fit <- jags.model(file="JAGS_gamma.txt",
                         data=list(n=n1, sum_y = sum_y1, sum_logy = sum_logy1, 
                                   tau0 = tau0, mu0 = mu0, zero = 0,
                                   kappa0 = kappa0, lambda0 = lambda0), 
                         n.chains = nchains)

update(model1.fit, burnin)
model1.samples <- coda.samples(model1.fit, c("alpha", "beta"), n.iter=ndraws, thin=nthin)

## save posterior draws for group A
alpha.1 <- unlist(model1.samples[,1])
beta.1 <- unlist(model1.samples[,2])

## repeat process for group B
set.seed(6)
sum_y2 <- sum(w2)
sum_logy2 <- sum(log(w2))
n2 <- length(w2)
model2.fit <- jags.model(file="JAGS_gamma.txt",
                         data=list(n=n2, sum_y = sum_y2, sum_logy = sum_logy2, 
                                   tau0 = tau0, mu0 = mu0, zero = 0,
                                   kappa0 = kappa0, lambda0 = lambda0), 
                         n.chains = nchains)

update(model2.fit, burnin)
model2.samples <- coda.samples(model2.fit, c("alpha", "beta"), n.iter=ndraws, thin=nthin)

## save posterior draws for group B
alpha.2 <- unlist(model2.samples[,1])
beta.2 <- unlist(model2.samples[,2])

## save posterior draws to .csv files
write.csv(alpha.1, "alpha1.csv", row.names = FALSE)
write.csv(beta.1, "beta1.csv", row.names = FALSE)
write.csv(alpha.2, "alpha2.csv", row.names = FALSE)
write.csv(beta.2, "beta2.csv", row.names = FALSE)

## obtain eta values under H1 using these posterior means
eta <- round(c(mean(alpha.1), mean(beta.1), mean(alpha.2), mean(beta.2)),2)
## these values are eta <- c(10.04, 12.37, 14.46, 14.90)

## load in data from .csv files if necessary
alpha.1 <- read.csv("alpha1.csv")$x
beta.1 <- read.csv("beta1.csv")$x
alpha.2 <- read.csv("alpha2.csv")$x
beta.2 <- read.csv("beta2.csv")$x

## fit informative priors that are approximately gamma distributed
## define factor by which to inflate the variance
## we do not inflate the variance for the gamma example,
## so this factor is 1
var_inf <- 1

## consider alpha1 (a_A in the paper)
## find the posterior mode and variance
a1dens <- density(alpha.1, n = 1024, bw = "SJ")
a1_mode <- a1dens$x[which.max(a1dens$y)]
a1_var <- var(alpha.1)

## solve system of equations to get desired gamma distribution
b1_star <- (-a1_mode - sqrt(a1_mode^2 + 4*var_inf*a1_var))/(-2*var_inf*a1_var); b1_star
a1_star <- 1 + a1_mode*b1_star; a1_star
informs <- c(a1_star, b1_star)

## consider beta1 (b_A in the paper)
## find the posterior mode and variance
a1dens <- density(beta.1, n = 1024, bw = "SJ")
a1_mode <- a1dens$x[which.max(a1dens$y)]
a1_var <- var(beta.1)

## solve system of equations to get desired gamma distribution
b1_star <- (-a1_mode - sqrt(a1_mode^2 + 4*var_inf*a1_var))/(-2*var_inf*a1_var); b1_star
a1_star <- 1 + a1_mode*b1_star; a1_star
informs <- rbind(informs, c(a1_star, b1_star))

## consider alpha2 (a_B in the paper)
## find the posterior mode and variance
a1dens <- density(alpha.2, n = 1024, bw = "SJ")
a1_mode <- a1dens$x[which.max(a1dens$y)]
a1_var <- var(alpha.2)

## solve system of equations to get desired gamma distribution
b1_star <- (-a1_mode - sqrt(a1_mode^2 + 4*var_inf*a1_var))/(-2*var_inf*a1_var); b1_star
a1_star <- 1 + a1_mode*b1_star; a1_star
informs <- rbind(informs, c(a1_star, b1_star))

## consider beta2 (b_B in the paper)
## find the posterior mode and variance
a1dens <- density(beta.2, n = 1024, bw = "SJ")
a1_mode <- a1dens$x[which.max(a1dens$y)]
a1_var <- var(beta.2)

## solve system of equations to get desired gamma distribution
b1_star <- (-a1_mode - sqrt(a1_mode^2 + 4*var_inf*a1_var))/(-2*var_inf*a1_var); b1_star
a1_star <- 1 + a1_mode*b1_star; a1_star
informs <- rbind(informs, c(a1_star, b1_star))

## output informative gamma distribution hyperparameters to .csv file
write.csv(informs, "informs_gamma.csv", row.names = FALSE)

## function to do MCMC for the posterior of theta1/theta2
GammaPostMean <- function(y1, y2, mu1, tau1, kappa1, lambda1,
                      mu2, tau2, kappa2, lambda2, burnin = 1000, nchains = 1,
                      nthin = 2, ndraws = 200000){
  
  ## y1 and y2 are the observations in each group (w in paper)
  ## alpha_j has a Gamma(mu_j, tau_j) prior, where tau_j is a rate (a in paper)
  ## beta_j has a Gamma(kappa_j, lambda_j) prior, where lambda_j is a rate (b in paper)
  ## burnin is the number of MCMC iterations to discard at the start of each chain
  ## nchains is the number of chains to generate
  ## nthin is the thinning parameter for the MCMC process
  ## ndraws is the number of draws to generate (excluding burnin but including thinned draws)
  
  ## precompute summary statistics to make JAGS faster
  ## implement MCMC for group A
  sum_y1 <- sum(y1)
  sum_logy1 <- sum(log(y1))
  n1 <- length(y1)
  model1.fit <- jags.model(file="JAGS_gamma.txt",
                           data=list(n=n1, sum_y = sum_y1, sum_logy = sum_logy1, 
                                     tau0 = tau1, mu0 = mu1, zero = 0,
                                     kappa0 = kappa1, lambda0 = lambda1), 
                           n.chains = nchains, quiet = TRUE)
  
  update(model1.fit, burnin, progress.bar = "none")
  model1.samples <- coda.samples(model1.fit, c("alpha", "beta"), n.iter=ndraws, thin=nthin, progress.bar = "none")
  
  ## get posterior draws for alpha_j and beta_j
  alpha.1 <- unlist(model1.samples[,1])
  beta.1 <- unlist(model1.samples[,2])
  
  ## implement MCMC for group B
  sum_y2 <- sum(y2)
  sum_logy2 <- sum(log(y2))
  n2 <- length(y2)
  model2.fit <- jags.model(file="JAGS_gamma.txt",
                           data=list(n=n2, sum_y = sum_y2, sum_logy = sum_logy2, 
                                     tau0 = tau2, mu0 = mu2, zero = 0,
                                     kappa0 = kappa2, lambda0 = lambda2), 
                           n.chains = nchains, quiet = TRUE)
  
  update(model2.fit, burnin, progress.bar = "none")
  model2.samples <- coda.samples(model2.fit, c("alpha", "beta"), n.iter=ndraws, thin=nthin, progress.bar = "none")
  
  alpha.2 <- unlist(model2.samples[,1])
  beta.2 <- unlist(model2.samples[,2])
  
  ## obtain posterior draws for theta (the ratio of the two gamma means)
  theta1 <- alpha.1/beta.1
  theta2 <- alpha.2/beta.2
  theta <- theta1/theta2
  theta <- ifelse(is.na(theta), Inf, theta)
  return(theta)
}