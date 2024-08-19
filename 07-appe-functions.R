## functions for Weibull numerical study

## BEGIN SETUP ##

## load necessary packages
require(qrng)
require(nleqslv)

## function to solve for the posterior mode of each group
## x is the point (logalpha, logbeta) at which this function is evaluated at
## y = c(length(data), sum(log(data)), sum(data))
## alpha = a and beta = b in the paper
fn_grad <- function(x, y, mu, tau, kappa, lambda) {
  
  res1 <- -y[3]*(x[1] - y[1])*exp(y[2])^2 + y[3]*(x[2] - y[2])*exp(y[2])*(1 + digamma(1)) - exp(x[1])*tau + mu
  
  res2 <- y[3]*(x[1] - y[1])*exp(y[2])*(1 + digamma(1)) - y[3]*(x[2] - y[2])*y[4] - exp(x[2])*lambda + kappa
  
  return(c(res1, res2))
}

## function to calculate the covariance matrix for the
## approximate posterior using Laplace's method. u is log(alpha_j), and
## v is log(beta_j). yy_star is c(length(data), sum(log(data)), sum(data)),
## and hyper = c(tau_j, lambda_j) are the relevant hyperparameters.
## alpha = a and beta = b in the paper
calc_covar <- function(u, v, yy_star, hyper){
  l <- exp(u); k <- exp(v); n <- yy_star[3]
  tau <- hyper[1]; lambda <- hyper[2]
  d11 <- n*k^2 + l*tau
  d12 <- -n*k*(1 + digamma(1))
  d22 <- n*yy_star[4] + k*lambda
  mat <- rbind(c(d11, d12), c(d12, d22))
  return(solve(mat))
}

## function to calculate the logit of the posterior probability
## u is point from Sobol' sequence, params are the draw from Psi_0
## or Psi_1, deltas is the interval c(delta_L, delta_U), n_val is the 
## sample size presently explored, hyper is a matrix of the hyperparameters, 
## and q is the constant for imbalanced sample size determination
logitP <- function(u, params, deltas, n_val, hyper, q){
  
  ## return negative power if sample size is not positive
  if (n_val <= 0){return(-1.5)}
  
  wei_lambda.1 <- params[1]
  wei_k.1 <- params[2]
  wei_lambda.2 <- params[3]
  wei_k.2 <- params[4]
  
  emc <- pi^2/6 + (digamma(1))^2 + 2*digamma(1)
  ## generate approximately normal MLEs for group A using delta method
  rho1 <- (1 + digamma(1))/sqrt(1 + emc)
  mat1 <- matrix((sqrt(6)/(pi*wei_k.1))^2*c(1 + emc, wei_k.1*(1 + digamma(1)), wei_k.1*(1 + digamma(1)), wei_k.1^2), nrow = 2)
  l1 <- qnorm(u[1], log(wei_lambda.1), sqrt(mat1[1,1]/(q*n_val)))
  k1 <- qnorm(u[2], log(wei_k.1) + rho1*(l1 - log(wei_lambda.1))*(sqrt(mat1[2,2])/sqrt(mat1[1,1])), sqrt(1- rho1^2)*sqrt(mat1[2,2]/(q*n_val)))
  
  ## generate approximately normal MLEs for group B using delta method
  rho2 <- (1 + digamma(1))/sqrt(1 + emc)
  mat2 <- matrix((sqrt(6)/(pi*wei_k.2))^2*c(1 + emc, wei_k.2*(1 + digamma(1)), wei_k.2*(1 + digamma(1)), wei_k.2^2), nrow = 2)
  l2 <- qnorm(u[3], log(wei_lambda.2), sqrt(mat2[1,1]/(n_val)))
  k2 <- qnorm(u[4], log(wei_k.2) + rho2*(l2 - log(wei_lambda.2))*(sqrt(mat2[2,2])/sqrt(mat2[1,1])), sqrt(1- rho2^2)*sqrt(mat2[2,2]/(n_val)))
  
  ## summarize information from first group of data (faster computation)
  yy_star1 <- c(l1, k1, q*n_val, 1 + emc)
  ## find posterior modes for the first group (logalpha and logbeta)
  modes <- nleqslv(c(l1, k1), fn_grad, y = yy_star1, mu = hyper[1,1], tau = hyper[1,2],
                   kappa = hyper[2,1], lambda = hyper[2,2] )$x
  
  mat1_new <- calc_covar(modes[1], modes[2], yy_star1, c(hyper[1,2], hyper[2,2]))
  ## exponentiate modes to return to standard scale
  modes1 <- exp(modes)
  
  ## repeat all steps for group B
  yy_star2 <- c(l2, k2, n_val, 1 + emc)
  modes <- nleqslv(c(l2, k2), fn_grad, y = yy_star2, mu = hyper[3,1], tau = hyper[3,2],
                   kappa = hyper[4,1], lambda = hyper[4,2] )$x
  
  mat2_new <- calc_covar(modes[1], modes[2], yy_star2, c(hyper[3,2], hyper[4,2]))
  modes2 <- exp(modes)
  l1 <- modes1[1]; k1 <- modes1[2]; l2 <- modes2[1]; k2 <- modes2[2]
  
  ## ensure no modes are 0 due to underflow errors
  if(max(l1 <= 0, k1 <= 0, l2 <= 0, k2<= 0)){return(-1.5)}
  
  ## define theta metrics in terms of design values
  theta1 <- l1*gamma(1 + 1/k1)
  theta2 <- l2*gamma(1 + 1/k2)
  
  if(max(theta1 <= 0, !is.finite(theta1), theta2 <= 0, !is.finite(theta2))){return(-1.5)}
  
  ## compute partial derivatives of theta with respect to logalpha and logbeta for each group
  ## this is different from the previous function that computes the derivatives with respect
  ## to alpha and beta
  d_lambda1 <- l1*gamma(1 + 1/k1)
  
  d_k1 <- -(l1/k1)*gamma(1 + 1/k1)*digamma(1 + 1/k1)
  
  d_lambda2 <- l2*gamma(1 + 1/k2)
  
  d_k2 <- -(l2/k2)*gamma(1 + 1/k2)*digamma(1 + 1/k2)
  
  ## apply the delta method to get the limiting variance for each group log mean
  avar1 <- t(c(d_lambda1, d_k1))%*%mat1_new%*%c(d_lambda1, d_k1)
  
  avar2 <- t(c(d_lambda2, d_k2))%*%mat2_new%*%c(d_lambda2, d_k2)
  
  ## apply the delta method to get the limiting variance for logtheta
  Fish_ratio_mu <- avar1/theta1^2 + avar2/theta2^2 
  
  ## return negative power if division causes NA/Inf values
  if (is.na(Fish_ratio_mu)){return(-1.5)}
  if (Fish_ratio_mu < -10000){return(-1.5)}
  
  ## return power based on normal approximation induced by Bernstein-von Mises
  realP <- pnorm(deltas[2], log(theta1/theta2), sqrt(Fish_ratio_mu)) - pnorm(deltas[1], log(theta1/theta2), sqrt(Fish_ratio_mu))
  
  ## slight perturbation of the posterior probability if it is too close to 0 or 1
  ## this ensures that the logit of the posterior probability is finite.
  if (realP > 1 - 10^(-7)){
    realP <- 1 - 10^(-7)
  }
  else if (realP < .Machine$double.eps){
    realP <- .Machine$double.eps
  }
  
  ## return the logit of the posterior probability
  return(log(realP) - log(1 - realP))
}

## Code to implement Algorithm 3 for the Weibull example. We explore sample sizes
## and return the optimal design (i.e., the (n, gamma) combination)
findGamma <- function(eta_plus1, eta_plus0, pwr, typeI, deltas, q,
                      alphas1, betas1, alphas2, betas2, m = 8192, m0 = 512,
                      seed_H1 = 1, seed_H0 = 2, upper_init = 1000, prng = FALSE, 
                      contour = FALSE){
  
  ## eta_plus1: parameter draws from Psi_1
  ## eta_plus0: parameter draws from Psi_0
  ## pwr: target power (1 - beta)
  ## typeI: desired type I error rate (alpha)
  ## q: constant for imbalanced sample size determination
  ## deltas: interval to define H1
  ## a_A ~ Gamma(alphas1[1], alphas1[2])
  ## b_A ~ Gamma(betas1[1], betas1[2])
  ## a_B ~ Gamma(alphas2[1], alphas2[2])
  ## b_B ~ Gamma(betas2[1], betas2[2])
  ## m: number of simulation repetitions
  ## m_0: number of repetitions to explore sampling distribution segments
  ## seed_H1: seed to generate Sobol' sequence for Psi_1
  ## seed_H0: seed to generate Sobol' sequence for Psi_0
  ## upper_init: maximum sample size to explore (not an issue if this is too small)
  ## prng: TRUE is pseudorandom sequence is used (for comparison)
  ## contour: when TRUE, returns logits of posterior probabilities that construct contour plot
  
  ## generate Sobol' sequence for each sampling distribution
  ## pseudorandom option available for comparison
  if (prng == FALSE){
    sob_H1 <- sobol(m, d = 4, randomize = "digital.shift", seed = seed_H1)
    sob_H0 <- sobol(m, d = 4, randomize = "digital.shift", seed = seed_H0)
  }
  else{
    set.seed(seed_H1); sob_H1 <- matrix(runif(4*m), ncol = 4)
    
    set.seed(seed_H0); sob_H0 <- matrix(runif(4*m), ncol = 4)
  }
  
  if (m < m0){stop("m0 should be less than m")}
  lower <- 1
  n <- upper_init
  upper <- upper_init
  H1_vec <- rep(0, m0)
  H0_vec <- rep(1, m0)
  
  mid1 <- ceiling((1 - typeI)*m0)
  mid2 <- floor((1-pwr)*m0)
  
  ## find a larger upper bound if the initial upper bound provided
  ## is not large enough
  while (sort(H0_vec)[mid1] > sort(H1_vec)[mid2]){
    H1_vec <- NULL
    H0_vec <- NULL
    ## implement Algorithm 2 for each of the first m0 points
    for (i in 1:m0){
      H1_vec[i] <- logitP(n_val = n, q = q,
                                  params = as.numeric(eta_plus1),
                                  delta = deltas, u = sob_H1[i,],
                                  hyper = rbind(alphas1, betas1, alphas2, betas2))
      H0_vec[i] <- logitP(n_val = n, q = q,
                                params = as.numeric(eta_plus0),
                                delta = deltas, u = sob_H0[i,],
                                hyper = rbind(alphas1, betas1, alphas2, betas2))
    }
    ## increase upper bound if the criterion on the order statistics is not satisfied
    if (sort(H0_vec)[mid1] <= sort(H1_vec)[mid2]){
      H1_max <- H1_vec; H0_max <- H0_vec
      upper <- n
      n_max <- n
    }
    else{
      lower <- n
      H1_min <- H1_vec; H0_min <- H0_vec
      n_min <- n
      n <- 2*n
    }
  }
  
  ## get a lower point to construct linear approximations to the logits
  ## of the posterior probabilities (this is only required if the initial
  ## upper bound for n was sufficiently large)
  if(lower == 1){
    n <- ceiling(0.5*(upper + lower))
    H1_vec <- NULL
    H0_vec <- NULL
    for (i in 1:m0){
      H1_vec[i] <- logitP(n_val = n, q = q,
                                  params = as.numeric(eta_plus1),
                                  delta = deltas, u = sob_H1[i,],
                                  hyper = rbind(alphas1, betas1, alphas2, betas2))
      H0_vec[i] <- logitP(n_val = n, q = q,
                                params = as.numeric(eta_plus0),
                                delta = deltas, u = sob_H0[i,],
                                hyper = rbind(alphas1, betas1, alphas2, betas2))
    }
    ## if the criterion is still satisfied then this smaller sample size is still
    ## an upper bound
    if (sort(H0_vec)[mid1] <= sort(H1_vec)[mid2]){
      H1_min <- H1_vec; H0_min <- H0_vec
      upper <- n
      n_min <- n
    }
    else{
      lower <- n
      H1_min <- H1_vec; H0_min <- H0_vec
      n_min <- n
    }
  }
  
  ## construct the linear approximations under H1 and H0
  H1_slope <- (H1_max - H1_min)/(n_max-n_min)
  H0_slope <- (H0_max - H0_min)/(n_max-n_min)
  H1_int <- H1_min - H1_slope*n_min
  H0_int <- H0_min - H0_slope*n_min
  
  ## get an initial sample size using the linear approximations to the posterior 
  ## probabilities (no Algorithm 3)
  upper_temp <- upper
  lower_temp <- lower
  while ((upper_temp - lower_temp) > 1){
    n <- ceiling(0.5*(upper_temp + lower_temp))
    H1_vec <- NULL
    H0_vec <- NULL
    for (i in 1:m0){
      H1_vec <- H1_int + H1_slope*n
      
      H0_vec <- H0_int + H0_slope*n
    }
    if (sort(H0_vec)[mid1] <= sort(H1_vec)[mid2]){
      upper_temp <- n
    }
    else{
      lower_temp <- n
    }
  }
  
  ii <- 0
  n <- max(10,upper_temp)
  n_last <- n
  while (ii < 1 | (abs(n_last - n)/n) >= 0.1 | abs(n_last - n) > 5){
    ii <- ii + 1
    n_last <- n
    H1_vec <- NULL
    H0_vec <- NULL
    ## we next run Algorithm 3 at that initialize sample size for the first m0 points
    for (i in 1:m0){
      H1_vec[i] <- logitP(n_val = n, q = q,
                                  params = as.numeric(eta_plus1),
                                  delta = deltas, u = sob_H1[i,],
                                  hyper = rbind(alphas1, betas1, alphas2, betas2))
      H0_vec[i] <- logitP(n_val = n, q = q,
                                params = as.numeric(eta_plus0),
                                delta = deltas, u = sob_H0[i,],
                                hyper = rbind(alphas1, betas1, alphas2, betas2))
    }
    ## these posterior probabilities will be used later to create better linear approximations
    ## on the logit scale (we need to determine whether n is a lower or upper bound here with
    ## respect to the first m0 points)
    if (sort(H0_vec)[mid1] <= sort(H1_vec)[mid2]){
      upper <- n
      if (n < n_min){
        n_max <- n_min
        H1_max <- H1_min; H0_max <- H0_min
        n_min <- n
        H1_min <- H1_vec; H0_min <- H0_vec
      }
      else{
        n_max <- n
        H1_max <- H1_vec; H0_max <- H0_vec
      }
    }
    else{
      lower <- n
      if (n > n_min){
        n_min <- n
        H1_min <- H1_vec; H0_min <- H0_vec
      }
      else{
        n_max <- n_min
        H1_max <- H1_min; H0_max <- H0_min
        n_min <- n
        H1_min <- H1_vec; H0_min <- H0_vec
      }
    }
    
    ## now use the probabilities at n and the other bound for binary search to get
    ## better linear approximations on the logit scale. We then explore sample sizes
    ## again with the first m0 points using these linear approximations (no Algorithm 3)
    H1_slope <- (H1_max - H1_min)/(n_max-n_min)
    H0_slope <- (H0_max - H0_min)/(n_max-n_min)
    H1_int <- H1_min - H1_slope*n_min
    H0_int <- H0_min - H0_slope*n_min
    upper_temp <- upper
    lower_temp <- lower
    while ((upper_temp - lower_temp) > 1){
      n <- ceiling(0.5*(upper_temp + lower_temp))
      H1_vec <- NULL
      H0_vec <- NULL
      for (i in 1:m0){
        H1_vec <- H1_int + H1_slope*n
        
        H0_vec <- H0_int + H0_slope*n
      }
      if (sort(H0_vec)[mid1] <= sort(H1_vec)[mid2]){
        upper_temp <- n
      }
      else{
        lower_temp <- n
      }
    }
    n <- upper_temp
    ## we check how much the sample size recommendation changed with respect to
    ## what we started this iteration of the floor loop with; if there is not
    ## much change, then we just take this new recommendation as n(0)
  }
  
  ## compute the posterior probabilities using Algorithm 3 for all m points for
  ## each Sobol' sequence at the initial sample size n(0)
  n_init <- n
  H1_vec <- NULL
  H0_vec <- NULL
  for (i in 1:m){
    H1_vec[i] <- logitP(n_val = n, q = q,
                                params = as.numeric(eta_plus1),
                                delta = deltas, u = sob_H1[i,],
                                hyper = rbind(alphas1, betas1, alphas2, betas2))
    H0_vec[i] <- logitP(n_val = n, q = q,
                              params = as.numeric(eta_plus0),
                              delta = deltas, u = sob_H0[i,],
                              hyper = rbind(alphas1, betas1, alphas2, betas2))
  }
  
  ## save the logits of the probabilities
  H1_vec_n0 <- H1_vec
  H0_vec_n0 <- H0_vec
  n0 <- n
  
  ## these flags denote whether alpha and beta are too small to explore
  ## symmetrically around the relevant order statistic
  flag1 <- ifelse(2*typeI*m < m0, 1, 0)
  flag2 <- ifelse(2*(1 - pwr)*m < m0, 1, 0)
  
  ## if beta is too small, we take the smallest m0 order statistics
  if (flag2){
    low2 <- 1; high2 <- m0
  }
  ## otherwise we explore symmetrically around the order statistic
  else{
    low2 <- floor((1 - pwr)*m) - 0.5*m0 + 1
    high2 <- low2 + m0 - 1
  }
  
  ## if alpha is too small, we take the largest m0 order statistics
  if (flag1){
    low1 <- m - m0 + 1; high1 <- m
  }
  else{
    low1 <- floor((1 - typeI)*m) - 0.5*m0 + 1
    high1 <- low1 + m0 - 1
  }
  
  ## we update the indices of the order statistics to reflect the larger
  ## set of points (m vs. m0)
  mid1 <- ceiling((1 - typeI)*m)
  mid2 <- floor((1-pwr)*m)
  
  ## if we do not satisfy the criterion on the order statistics, we create the linear
  ## approximation using a larger second sample size 
  if (sort(H0_vec)[mid1] > sort(H1_vec)[mid2]){
    H1_low <- H1_vec; H0_low <- H0_vec
    lower <- n
    n_low <- n
    ## choose a larger sample size n(1)
    n <- max(ceiling(1.1*n), n + 5)
    H1_vec <- NULL
    H0_vec <- NULL
    ## use Algorithm 3 to approximate the posterior probabilities at all m points
    for (i in 1:m){
      H1_vec[i] <- logitP(n_val = n, q = q,
                                  params = as.numeric(eta_plus1),
                                  delta = deltas, u = sob_H1[i,],
                                  hyper = rbind(alphas1, betas1, alphas2, betas2))
      H0_vec[i] <- logitP(n_val = n, q = q,
                                params = as.numeric(eta_plus0),
                                delta = deltas, u = sob_H0[i,],
                                hyper = rbind(alphas1, betas1, alphas2, betas2))
    }
    
    ## save the logits of the posterior probabilities
    H1_vec_n1 <- H1_vec
    H0_vec_n1 <- H0_vec
    n1 <- n
    
    n_high <- n
    H1_high <- H1_vec; H0_high <- H0_vec
    
    ## obtain the linear approximations for all m points
    H1_slope <- (H1_high - H1_low)/(n_high-n_low)
    H0_slope <- (H0_high - H0_low)/(n_high-n_low)
    H1_int <- H1_low - H1_slope*n_low
    H0_int <- H0_low - H0_slope*n_low
    
    ## if the criterion on the order statistics is still not satisfied for this
    ## larger sample size, we need to find a larger upper bound for the binary search
    ## we do not do this using Algorithm 3 -- just the linear approximations to the
    ## posterior probabilities on the logit scale.
    if (sort(H0_vec)[mid1] <= sort(H1_vec)[mid2]){
      upper <- n
    }
    else {
      n <- ceiling(1.5*n_low)
      H1_vec <- H1_int + H1_slope*n
      H0_vec <- H0_int + H0_slope*n
      upper <- n
      
      while(sort(H0_vec)[mid1] > sort(H1_vec)[mid2]){
        n <- ceiling(1.5*n)
        H1_vec <- H1_int + H1_slope*n
        H0_vec <- H0_int + H0_slope*n
        if (sort(H0_vec)[mid1] <= sort(H1_vec)[mid2]){
          upper <- n
        }
      }
    }
  }
  ## if we satisfy the criterion on the order statistics, we create the linear
  ## approximation using a smaller second sample size 
  else{
    H1_high <- H1_vec; H0_high <- H0_vec
    upper <- n
    n_high <- n
    ## choose a smaller sample size n(1)
    n <- min(floor(0.9*n), n - 5)
    H1_vec <- NULL
    H0_vec <- NULL
    ## use Algorithm 3 to approximate the posterior probabilities at all m points
    for (i in 1:m){
      H1_vec[i] <- logitP(n_val = n, q = q,
                                  params = as.numeric(eta_plus1),
                                  delta = deltas, u = sob_H1[i,],
                                  hyper = rbind(alphas1, betas1, alphas2, betas2))
      H0_vec[i] <- logitP(n_val = n, q = q,
                                params = as.numeric(eta_plus0),
                                delta = deltas, u = sob_H0[i,],
                                hyper = rbind(alphas1, betas1, alphas2, betas2))
    }
    
    ## save the logits of the posterior probabilities
    H1_vec_n1 <- H1_vec
    H0_vec_n1 <- H0_vec
    n1 <- n
    
    n_low <- n
    H1_low <- H1_vec; H0_low <- H0_vec
    
    ## obtain the linear approximations for all m points
    H1_slope <- (H1_high - H1_low)/(n_high-n_low)
    H0_slope <- (H0_high - H0_low)/(n_high-n_low)
    H1_int <- H1_low - H1_slope*n_low
    H0_int <- H0_low - H0_slope*n_low
    
    ## if the criterion on the order statistics is still satisfied for this
    ## larger sample size, we need to find a smaller lower bound for the binary search
    if (sort(H0_vec)[mid1] > sort(H1_vec)[mid2]){
      lower <- n
    }
    else{
      n <- floor(0.5*n_high)
      H1_vec <- H1_int + H1_slope*n
      H0_vec <- H0_int + H0_slope*n
      lower <- n
      
      while(sort(H0_vec)[mid1] <= sort(H1_vec)[mid2]){
        n <- floor(0.5*n)
        H1_vec <- H1_int + H1_slope*n
        H0_vec <- H0_int + H0_slope*n
        if (sort(H0_vec)[mid1] > sort(H1_vec)[mid2]){
          lower <- n
        }
      }
    }
  }
  
  ## implement binary search given these bounds for the sample size
  while ((upper - lower) > 1){
    n <- ceiling(0.5*(upper + lower))
    
    ## use the linear approximations to select the points from the Sobol' 
    ## sequence in a targeted way; these points will be used to explore the 
    ## sampling distributions of posterior probabilities in a targeted way
    H1_vec <- H1_int + H1_slope*n
    H0_vec <- H0_int + H0_slope*n
    ## get the indices for these points
    sub_H1 <- which(rank(H1_vec) >= low2 & rank(H1_vec) <= high2)
   sub_H0 <- which(rank(H0_vec) >= low1 & rank(H0_vec) <= high1)
    
    H1_vec <- NULL
    H0_vec <- NULL
    for (i in 1:length(sub_H1)){
      H1_vec[i] <- logitP(n_val = n, q = q,
                                  params = as.numeric(eta_plus1),
                                  delta = deltas, u = sob_H1[sub_H1[i],],
                                  hyper = rbind(alphas1, betas1, alphas2, betas2))
      H0_vec[i] <- logitP(n_val = n, q = q,
                                params = as.numeric(eta_plus0),
                                delta = deltas, u = sob_H0[sub_H0[i],],
                                hyper = rbind(alphas1, betas1, alphas2, betas2))
    }
    
    ## for the unselected points, we just assume that their probabilities are less (greater)
    ## than the relevant order statistic if their anticipated probability based on the 
    ## linear approximations on the logit scale is less (greater) than the relevant order statistic
    H0_aug <- c(rep(c(min(H0_vec)-1), low1 - 1), H0_vec, rep(c(max(H0_vec)+1), m - high1))
    H1_aug <- c(rep(c(min(H1_vec)-1), low2 - 1), H1_vec, rep(c(max(H1_vec)+1), m - high2))
    if (sort(H0_aug)[mid1] <= sort(H1_aug)[mid2]){
      H1_max <- H1_vec; H0_max <- H0_vec
      sub_H1_max <- sub_H1; sub_H0_max <-sub_H0
      upper <- n
    }
    else{
      lower <- n
    }
  }
  
  n <- upper
  ## we don't need to reapproximate the same posterior probabilities if we have already explored
  ## this sample size
  if (n == n_high){
    
    ## do not need to do this; it is just so we have three sample sizes to construct each contour plot
    n <- ifelse(n0 == n_high, max(ceiling(1.1*n0), n0 + 5), max(ceiling(1.1*n1), n1 + 5))
    H1_vec <- NULL
    H0_vec <- NULL
    ## use Algorithm 3 to approximate the posterior probabilities at all m points
    for (i in 1:m){
      H1_vec[i] <- logitP(n_val = n, q = q,
                                  params = as.numeric(eta_plus1),
                                  delta = deltas, u = sob_H1[i,],
                                  hyper = rbind(alphas1, betas1, alphas2, betas2))
      H0_vec[i] <- logitP(n_val = n, q = q,
                                params = as.numeric(eta_plus0),
                                delta = deltas, u = sob_H0[i,],
                                hyper = rbind(alphas1, betas1, alphas2, betas2))
    }
    
    ## save the logits of the posterior probabilities
    H1_vec_n2 <- H1_vec
    H0_vec_n2 <- H0_vec
    n2 <- n
    
    n_final <- ifelse(n1 > n0, 1, 0)
    
    ## first component of list is (n, gamma, relevant order statistic of H1 probabilities, confirmatory type I error
    ## estimate, confirmatory power estimate, n0, n1, n2)
    ## components 2 to 6 of the list are only needed for the contour plot (they are the posterior probabilities under
    ## H0 and H1 for the sample sizes n0, n1, n2)
    ## we sort the sample sizes and the posterior probabilities to help with making the plot
    results <- list(c(get(paste0("n", n_final)), 1/(1+exp(-as.numeric(sort(get(paste0("H0_vec_n", n_final)))[mid1]))), 
                      1/(1+exp(-as.numeric(sort(get(paste0("H1_vec_n", n_final)))[mid2]))),
                      as.numeric(mean(get(paste0("H0_vec_n", n_final)) > sort(get(paste0("H0_vec_n", n_final)))[mid1])), 
                      as.numeric(mean(get(paste0("H1_vec_n", n_final)) > sort(get(paste0("H0_vec_n", n_final)))[mid1])),
                      sort(c(n0, n1, n2))),
                    get(paste0("H0_vec_n", which.min(c(n0, n1, n2))-1)), get(paste0("H1_vec_n", which.min(c(n0, n1, n2))-1)), 
                    get(paste0("H0_vec_n", which(sort(c(n0, n1, n2))[2] == c(n0, n1, n2))-1)), get(paste0("H1_vec_n", which(sort(c(n0, n1, n2))[2] == c(n0, n1, n2))-1)),
                    get(paste0("H0_vec_n", which.max(c(n0, n1, n2))-1)), get(paste0("H1_vec_n", which.max(c(n0, n1, n2))-1)))
    
    if (contour == TRUE){
      return(results)
    }
    else{
      return(results[[1]][1:5])
    }
  }
  
  if (n == n_low){
    
    ## do not need to do this; it is just so we have three sample sizes to construct each contour plot
    n <- ifelse(n0 == n_low, min(floor(0.9*n0), n0 - 5), min(floor(0.9*n1), n1 - 5))
    H1_vec <- NULL
    H0_vec <- NULL
    ## use Algorithm 2 to approximate the posterior probabilities at all m points
    for (i in 1:m){
      H1_vec[i] <- logitP(n_val = n, q = q,
                                  params = as.numeric(eta_plus1),
                                  delta = deltas, u = sob_H1[i,],
                                  hyper = rbind(alphas1, betas1, alphas2, betas2))
      H0_vec[i] <- logitP(n_val = n, q = q,
                                params = as.numeric(eta_plus0),
                                delta = deltas, u = sob_H0[i,],
                                hyper = rbind(alphas1, betas1, alphas2, betas2))
    }
    
    ## save the logits of the posterior probabilities
    H1_vec_n2 <- H1_vec
    H0_vec_n2 <- H0_vec
    n2 <- n
    
    n_final <- ifelse(n1 < n0, 1, 0)
    
    ## first component of list is (n, gamma, relevant order statistic of H1 probabilities, confirmatory type I error
    ## estimate, confirmatory power estimate, n0, n1, n2)
    ## components 2 to 6 of the list are only needed for the contour plot (they are the posterior probabilities under
    ## H1 and H0 for the sample sizes n0, n1, n2)
    ## we sort the sample sizes and the posterior probabilities to help with making the plot
    results <- list(c(get(paste0("n", n_final)), 1/(1+exp(-as.numeric(sort(get(paste0("H0_vec_n", n_final)))[mid1]))), 
                      1/(1+exp(-as.numeric(sort(get(paste0("H1_vec_n", n_final)))[mid2]))),
                      as.numeric(mean(get(paste0("H0_vec_n", n_final)) > sort(get(paste0("H0_vec_n", n_final)))[mid1])), 
                      as.numeric(mean(get(paste0("H1_vec_n", n_final)) > sort(get(paste0("H0_vec_n", n_final)))[mid1])),
                      sort(c(n0, n1, n2))),
                    get(paste0("H0_vec_n", which.min(c(n0, n1, n2))-1)), get(paste0("H1_vec_n", which.min(c(n0, n1, n2))-1)), 
                    get(paste0("H0_vec_n", which(sort(c(n0, n1, n2))[2] == c(n0, n1, n2))-1)), get(paste0("H1_vec_n", which(sort(c(n0, n1, n2))[2] == c(n0, n1, n2))-1)),
                    get(paste0("H0_vec_n", which.max(c(n0, n1, n2))-1)), get(paste0("H1_vec_n", which.max(c(n0, n1, n2))-1)))
    
    if (contour == TRUE){
      return(results)
    }
    else{
      return(results[[1]][1:5])
    }
  }
  
  ## otherwise, we approximate the remaining posterior probabilities that were not
  ## selected at the final sample size recommendation
  H1_vec <- NULL
  H0_vec <- NULL
  remain_H1 <- subset(1:m, ! 1:m %in% sub_H1_max)
  remain_H0 <- subset(1:m, ! 1:m %in% sub_H0_max)
  for (i in remain_H1){
    H1_vec[i] <- logitP(n_val = n, q = q,
                                params = as.numeric(eta_plus1),
                                delta = deltas, u = sob_H1[i,],
                                hyper = rbind(alphas1, betas1, alphas2, betas2))
  }
  for (i in remain_H0){
    H0_vec[i] <- logitP(n_val = n, q = q,
                              params = as.numeric(eta_plus0),
                              delta = deltas, u = sob_H0[i,],
                              hyper = rbind(alphas1, betas1, alphas2, betas2))
  }
  H1_vec[sub_H1_max] <- H1_max; H0_vec[sub_H0_max] <- H0_max
  
  ## save the logits of the posterior probabilities
  H1_vec_n2 <- H1_vec
  H0_vec_n2 <- H0_vec
  n2 <- n
  
  ## first component of list is (n, gamma, relevant order statistic of H1 probabilities, confirmatory type I error
  ## estimate, confirmatory power estimate, n0, n1, n2)
  ## components 2 to 6 of the list are only needed for the contour plot (they are the posterior probabilities under
  ## H1 and H0 for the sample sizes n0, n1, n2)
  ## we sort the sample sizes and the posterior probabilities to help with making the plot
  
  results <- list(c(n, 1/(1 + exp(-as.numeric(sort(H0_vec)[mid1]))),
                    1/(1 + exp(-as.numeric(sort(H1_vec)[mid2]))),
                    as.numeric(mean(H0_vec > sort(H0_vec)[mid1])), as.numeric(mean(H1_vec > sort(H0_vec)[mid1])),
                    sort(c(n0, n1, n2))),
                  get(paste0("H0_vec_n", which.min(c(n0, n1, n2))-1)), get(paste0("H1_vec_n", which.min(c(n0, n1, n2))-1)), 
                  get(paste0("H0_vec_n", which(sort(c(n0, n1, n2))[2] == c(n0, n1, n2))-1)), get(paste0("H1_vec_n", which(sort(c(n0, n1, n2))[2] == c(n0, n1, n2))-1)),
                  get(paste0("H0_vec_n", which.max(c(n0, n1, n2))-1)), get(paste0("H1_vec_n", which.max(c(n0, n1, n2))-1)))
  
  if (contour == TRUE){
    return(results)
  }
  else{
    return(results[[1]][1:5])
  }
}

## alternative function to explore sample sizes using all points from the 
## hypercube for each sample size that we consider (for confirmation purposes)
findGammaFull <- function(eta_plus1, eta_plus0, pwr, typeI, deltas, q,
                          alphas1, betas1, alphas2, betas2, m = 8192, m0 = 512,
                          seed_H1 = 1, seed_H0 = 2, upper_init = 1000){
  
  ## inputs are the same as findGamma(), but m0 is not used
  
  ## generate Sobol' sequence for each sampling distribution
  sob_H1 <- sobol(m, d = 4, randomize = "digital.shift", seed = seed_H1)
  sob_H0 <- sobol(m, d = 4, randomize = "digital.shift", seed = seed_H0)
  
  if (m < m0){stop("m0 should be less than m")}
  lower <- 1
  n <- upper_init
  upper <- upper_init
  H1_vec <- rep(0, m)
  H0_vec <- rep(1, m)
  
  mid1 <- ceiling((1 - typeI)*m)
  mid2 <- floor((1-pwr)*m)
  
  ## find a larger upper bound if the initial upper bound provided
  ## is not large enough
  while (sort(H0_vec)[mid1] > sort(H1_vec)[mid2]){
    H1_vec <- NULL
    H0_vec <- NULL
    ## implement Algorithm 2 for each of the first m0 points
    for (i in 1:m){
      H1_vec[i] <- logitP(n_val = n, q = q,
                                  params = as.numeric(eta_plus1),
                                  delta = deltas, u = sob_H1[i,],
                                  hyper = rbind(alphas1, betas1, alphas2, betas2))
      H0_vec[i] <- logitP(n_val = n, q = q,
                                params = as.numeric(eta_plus0),
                                delta = deltas, u = sob_H0[i,],
                                hyper = rbind(alphas1, betas1, alphas2, betas2))
    }
    ## increase upper bound if the criterion on the order statistics is not satisfied
    if (sort(H0_vec)[mid1] <= sort(H1_vec)[mid2]){
      H1_max <- H1_vec; H0_max <- H0_vec
      upper <- n
      n_max <- n
    }
    else{
      lower <- n
      H1_min <- H1_vec; H0_min <- H0_vec
      n_min <- n
      n <- 2*n
    }
  }
  
  ## implement binary search given these bounds for the sample size
  while ((upper - lower) > 1){
    n <- ceiling(0.5*(upper + lower))
    
    H1_vec <- NULL
    H0_vec <- NULL
    for (i in 1:m){
      H1_vec[i] <- logitP(n_val = n, q = q,
                                  params = as.numeric(eta_plus1),
                                  delta = deltas, u = sob_H1[i,],
                                  hyper = rbind(alphas1, betas1, alphas2, betas2))
      H0_vec[i] <- logitP(n_val = n, q = q,
                                params = as.numeric(eta_plus0),
                                delta = deltas, u = sob_H0[i,],
                                hyper = rbind(alphas1, betas1, alphas2, betas2))
    }
    
    if (sort(H0_vec)[mid1] <= sort(H1_vec)[mid2]){
      H1_max <- H1_vec; H0_max <- H0_vec
      upper <- n
    }
    else{
      lower <- n
    }
  }
  
  n <- upper
  H1_vec <- H1_max
  H0_vec <- H0_max
  
  ## list is (n, gamma, relevant order statistic of H1 probabilities, confirmatory type I error
  ## estimate, confirmatory power estimate)
  return(list(c(n, 1/(1 + exp(-as.numeric(sort(H0_vec)[mid1]))),
                1/(1 + exp(-as.numeric(sort(H1_vec)[mid2]))),
                as.numeric(mean(H0_vec > sort(H0_vec)[mid1])), as.numeric(mean(H1_vec > sort(H0_vec)[mid1])))))
  
}

## Version that finds n given fixed gamma
findn <- function(eta_plus1, eta_plus0, pwr, typeI, deltas, q,
                      alphas1, betas1, alphas2, betas2, m = 8192, m0 = 512,
                      seed_H1 = 1, seed_H0 = 2, upper_init = 1000, prng = FALSE, 
                      contour = FALSE){
  
  ## the inputs are the same as findGamma()
  
  ## however, the sampling distribution of posterior probabilities under H1 is always compared to
  ## the logit of typeI (for all sample sizes considered). Thus, the criterion for the type I error
  ## rate may not always be satisfied
  lA <- log(1-typeI) - log(typeI)
  
  ## generate Sobol' sequence for each sampling distribution
  ## pseudorandom option available for comparison
  if (prng == FALSE){
    sob_H1 <- sobol(m, d = 4, randomize = "digital.shift", seed = seed_H1)
    sob_H0 <- sobol(m, d = 4, randomize = "digital.shift", seed = seed_H0)
  }
  else{
    set.seed(seed_H1); sob_H1 <- matrix(runif(4*m), ncol = 4)
    
    set.seed(seed_H0); sob_H0 <- matrix(runif(4*m), ncol = 4)
  }
  
  if (m < m0){stop("m0 should be less than m")}
  lower <- 1
  n <- upper_init
  upper <- upper_init
  H1_vec <- rep(0, m0)
  H0_vec <- rep(1, m0)
  
  mid1 <- ceiling((1 - typeI)*m0)
  mid2 <- floor((1-pwr)*m0)
  
  ## find a larger upper bound if the initial upper bound provided
  ## is not large enough
  while (lA > sort(H1_vec)[mid2]){
    H1_vec <- NULL
    H0_vec <- NULL
    ## implement Algorithm 2 for each of the first m0 points
    for (i in 1:m0){
      H1_vec[i] <- logitP(n_val = n, q = q,
                                  params = as.numeric(eta_plus1),
                                  delta = deltas, u = sob_H1[i,],
                                  hyper = rbind(alphas1, betas1, alphas2, betas2))
      H0_vec[i] <- logitP(n_val = n, q = q,
                                params = as.numeric(eta_plus0),
                                delta = deltas, u = sob_H0[i,],
                                hyper = rbind(alphas1, betas1, alphas2, betas2))
    }
    ## increase upper bound if the criterion on the order statistics is not satisfied
    if (lA <= sort(H1_vec)[mid2]){
      H1_max <- H1_vec; H0_max <- H0_vec
      upper <- n
      n_max <- n
    }
    else{
      lower <- n
      H1_min <- H1_vec; H0_min <- H0_vec
      n_min <- n
      n <- 2*n
    }
  }
  
  ## get a lower point to construct linear approximations to the logits
  ## of the posterior probabilities (this is only required if the initial
  ## upper bound for n was sufficiently large)
  if(lower == 1){
    n <- ceiling(0.5*(upper + lower))
    H1_vec <- NULL
    H0_vec <- NULL
    for (i in 1:m0){
      H1_vec[i] <- logitP(n_val = n, q = q,
                                  params = as.numeric(eta_plus1),
                                  delta = deltas, u = sob_H1[i,],
                                  hyper = rbind(alphas1, betas1, alphas2, betas2))
      H0_vec[i] <- logitP(n_val = n, q = q,
                                params = as.numeric(eta_plus0),
                                delta = deltas, u = sob_H0[i,],
                                hyper = rbind(alphas1, betas1, alphas2, betas2))
    }
    ## if the criterion is still satisfied then this smaller sample size is still
    ## an upper bound
    if (lA <= sort(H1_vec)[mid2]){
      H1_min <- H1_vec; H0_min <- H0_vec
      upper <- n
      n_min <- n
    }
    else{
      lower <- n
      H1_min <- H1_vec; H0_min <- H0_vec
      n_min <- n
    }
  }
  
  ## construct the linear approximations under H1 and H0
  H1_slope <- (H1_max - H1_min)/(n_max-n_min)
  H0_slope <- (H0_max - H0_min)/(n_max-n_min)
  H1_int <- H1_min - H1_slope*n_min
  H0_int <- H0_min - H0_slope*n_min
  
  ## get an initial sample size using the linear approximations to the posterior 
  ## probabilities (no Algorithm 2)
  upper_temp <- upper
  lower_temp <- lower
  while ((upper_temp - lower_temp) > 1){
    n <- ceiling(0.5*(upper_temp + lower_temp))
    H1_vec <- NULL
    H0_vec <- NULL
    for (i in 1:m0){
      H1_vec <- H1_int + H1_slope*n
      
      H0_vec <- H0_int + H0_slope*n
    }
    if (lA <= sort(H1_vec)[mid2]){
      upper_temp <- n
    }
    else{
      lower_temp <- n
    }
  }
  
  ii <- 0
  n <- max(10,upper_temp)
  n_last <- n
  while (ii < 1 | (abs(n_last - n)/n) >= 0.1 | abs(n_last - n) > 5){
    ii <- ii + 1
    n_last <- n
    H1_vec <- NULL
    H0_vec <- NULL
    ## we next run Algorithm 2 at that initialize sample size for the first m0 points
    for (i in 1:m0){
      H1_vec[i] <- logitP(n_val = n, q = q,
                                  params = as.numeric(eta_plus1),
                                  delta = deltas, u = sob_H1[i,],
                                  hyper = rbind(alphas1, betas1, alphas2, betas2))
      H0_vec[i] <- logitP(n_val = n, q = q,
                                params = as.numeric(eta_plus0),
                                delta = deltas, u = sob_H0[i,],
                                hyper = rbind(alphas1, betas1, alphas2, betas2))
    }
    ## these posterior probabilities will be used later to create better linear approximations
    ## on the logit scale (we need to determine whether n is a lower or upper bound here with
    ## respect to the first m0 points)
    if (lA <= sort(H1_vec)[mid2]){
      upper <- n
      if (n < n_min){
        n_max <- n_min
        H1_max <- H1_min; H0_max <- H0_min
        n_min <- n
        H1_min <- H1_vec; H0_min <- H0_vec
      }
      else{
        n_max <- n
        H1_max <- H1_vec; H0_max <- H0_vec
      }
    }
    else{
      lower <- n
      if (n > n_min){
        n_min <- n
        H1_min <- H1_vec; H0_min <- H0_vec
      }
      else{
        n_max <- n_min
        H1_max <- H1_min; H0_max <- H0_min
        n_min <- n
        H1_min <- H1_vec; H0_min <- H0_vec
      }
    }
    
    ## now use the probabilities at n and the other bound for binary search to get
    ## better linear approximations on the logit scale. We then explore sample sizes
    ## again with the first m0 points using these linear approximations (no Algorithm 2)
    H1_slope <- (H1_max - H1_min)/(n_max-n_min)
    H0_slope <- (H0_max - H0_min)/(n_max-n_min)
    H1_int <- H1_min - H1_slope*n_min
    H0_int <- H0_min - H0_slope*n_min
    upper_temp <- upper
    lower_temp <- lower
    while ((upper_temp - lower_temp) > 1){
      n <- ceiling(0.5*(upper_temp + lower_temp))
      H1_vec <- NULL
      H0_vec <- NULL
      for (i in 1:m0){
        H1_vec <- H1_int + H1_slope*n
        
        H0_vec <- H0_int + H0_slope*n
      }
      if (lA <= sort(H1_vec)[mid2]){
        upper_temp <- n
      }
      else{
        lower_temp <- n
      }
    }
    n <- upper_temp
    ## we check how much the sample size recommendation changed with respect to
    ## what we started this iteration of the floor loop with; if there is not
    ## much change, then we just take this new recommendation as n(0)
  }
  
  ## compute the posterior probabilities using Algorithm 2 for all m points for
  ## each Sobol' sequence at the initial sample size n(0)
  n_init <- n
  H1_vec <- NULL
  H0_vec <- NULL
  for (i in 1:m){
    H1_vec[i] <- logitP(n_val = n, q = q,
                                params = as.numeric(eta_plus1),
                                delta = deltas, u = sob_H1[i,],
                                hyper = rbind(alphas1, betas1, alphas2, betas2))
    H0_vec[i] <- logitP(n_val = n, q = q,
                              params = as.numeric(eta_plus0),
                              delta = deltas, u = sob_H0[i,],
                              hyper = rbind(alphas1, betas1, alphas2, betas2))
  }
  
  ## save the logits of the probabilities
  H1_vec_n0 <- H1_vec
  H0_vec_n0 <- H0_vec
  n0 <- n
  
  ## these flags denote whether alpha and beta are too small to explore
  ## symmetrically around the relevant order statistic
  flag1 <- ifelse(2*typeI*m < m0, 1, 0)
  flag2 <- ifelse(2*(1 - pwr)*m < m0, 1, 0)
  
  ## if beta is too small, we take the smallest m0 order statistics
  if (flag2){
    low2 <- 1; high2 <- m0
  }
  ## otherwise we explore symmetrically around the order statistic
  else{
    low2 <- floor((1 - pwr)*m) - 0.5*m0 + 1
    high2 <- low2 + m0 - 1
  }
  
  ## if alpha is too small, we take the largest m0 order statistics
  if (flag1){
    low1 <- m - m0 + 1; high1 <- m
  }
  else{
    low1 <- floor((1 - typeI)*m) - 0.5*m0 + 1
    high1 <- low1 + m0 - 1
  }
  
  ## we update the indices of the order statistics to reflect the larger
  ## set of points (m vs. m0)
  mid1 <- ceiling((1 - typeI)*m)
  mid2 <- floor((1-pwr)*m)
  
  ## if we do not satisfy the criterion on the order statistics, we create the linear
  ## approximation using a larger second sample size 
  if (lA > sort(H1_vec)[mid2]){
    H1_low <- H1_vec; H0_low <- H0_vec
    lower <- n
    n_low <- n
    ## choose a larger sample size n(1)
    n <- max(ceiling(1.1*n), n + 5)
    H1_vec <- NULL
    H0_vec <- NULL
    ## use Algorithm 2 to approximate the posterior probabilities at all m points
    for (i in 1:m){
      H1_vec[i] <- logitP(n_val = n, q = q,
                                  params = as.numeric(eta_plus1),
                                  delta = deltas, u = sob_H1[i,],
                                  hyper = rbind(alphas1, betas1, alphas2, betas2))
      H0_vec[i] <- logitP(n_val = n, q = q,
                                params = as.numeric(eta_plus0),
                                delta = deltas, u = sob_H0[i,],
                                hyper = rbind(alphas1, betas1, alphas2, betas2))
    }
    
    ## save the logits of the posterior probabilities
    H1_vec_n1 <- H1_vec
    H0_vec_n1 <- H0_vec
    n1 <- n
    
    n_high <- n
    H1_high <- H1_vec; H0_high <- H0_vec
    
    ## obtain the linear approximations for all m points
    H1_slope <- (H1_high - H1_low)/(n_high-n_low)
    H0_slope <- (H0_high - H0_low)/(n_high-n_low)
    H1_int <- H1_low - H1_slope*n_low
    H0_int <- H0_low - H0_slope*n_low
    
    ## if the criterion on the order statistics is still not satisfied for this
    ## larger sample size, we need to find a larger upper bound for the binary search
    ## we do not do this using Algorithm 2 -- just the linear approximations to the
    ## posterior probabilities on the logit scale.
    if (lA <= sort(H1_vec)[mid2]){
      upper <- n
    }
    else {
      n <- ceiling(1.5*n_low)
      H1_vec <- H1_int + H1_slope*n
      H0_vec <- H0_int + H0_slope*n
      upper <- n
      
      while(lA > sort(H1_vec)[mid2]){
        n <- ceiling(1.5*n)
        H1_vec <- H1_int + H1_slope*n
        H0_vec <- H0_int + H0_slope*n
        if (lA <= sort(H1_vec)[mid2]){
          upper <- n
        }
      }
    }
  }
  ## if we satisfy the criterion on the order statistics, we create the linear
  ## approximation using a smaller second sample size 
  else{
    H1_high <- H1_vec; H0_high <- H0_vec
    upper <- n
    n_high <- n
    ## choose a smaller sample size n(1)
    n <- min(floor(0.9*n), n - 5)
    H1_vec <- NULL
    H0_vec <- NULL
    ## use Algorithm 2 to approximate the posterior probabilities at all m points
    for (i in 1:m){
      H1_vec[i] <- logitP(n_val = n, q = q,
                                  params = as.numeric(eta_plus1),
                                  delta = deltas, u = sob_H1[i,],
                                  hyper = rbind(alphas1, betas1, alphas2, betas2))
      H0_vec[i] <- logitP(n_val = n, q = q,
                                params = as.numeric(eta_plus0),
                                delta = deltas, u = sob_H0[i,],
                                hyper = rbind(alphas1, betas1, alphas2, betas2))
    }
    
    ## save the logits of the posterior probabilities
    H1_vec_n1 <- H1_vec
    H0_vec_n1 <- H0_vec
    n1 <- n
    
    n_low <- n
    H1_low <- H1_vec; H0_low <- H0_vec
    
    ## obtain the linear approximations for all m points
    H1_slope <- (H1_high - H1_low)/(n_high-n_low)
    H0_slope <- (H0_high - H0_low)/(n_high-n_low)
    H1_int <- H1_low - H1_slope*n_low
    H0_int <- H0_low - H0_slope*n_low
    
    ## if the criterion on the order statistics is still satisfied for this
    ## larger sample size, we need to find a smaller lower bound for the binary search
    if (lA > sort(H1_vec)[mid2]){
      lower <- n
    }
    else{
      n <- floor(0.5*n_high)
      H1_vec <- H1_int + H1_slope*n
      H0_vec <- H0_int + H0_slope*n
      lower <- n
      
      while(lA <= sort(H1_vec)[mid2]){
        n <- floor(0.5*n)
        H1_vec <- H1_int + H1_slope*n
        H0_vec <- H0_int + H0_slope*n
        if (lA > sort(H1_vec)[mid2]){
          lower <- n
        }
      }
    }
  }
  
  ## implement binary search given these bounds for the sample size
  while ((upper - lower) > 1){
    n <- ceiling(0.5*(upper + lower))
    
    ## use the linear approximations to select the points from the Sobol' 
    ## sequence in a targeted way; these points will be used to explore the 
    ## sampling distributions of posterior probabilities in a targeted way
    H1_vec <- H1_int + H1_slope*n
    H0_vec <- H0_int + H0_slope*n
    ## get the indices for these points
    sub_H1 <- which(rank(H1_vec) >= low2 & rank(H1_vec) <= high2)
   sub_H0 <- which(rank(H0_vec) >= low1 & rank(H0_vec) <= high1)
    
    H1_vec <- NULL
    H0_vec <- NULL
    for (i in 1:length(sub_H1)){
      H1_vec[i] <- logitP(n_val = n, q = q,
                                  params = as.numeric(eta_plus1),
                                  delta = deltas, u = sob_H1[sub_H1[i],],
                                  hyper = rbind(alphas1, betas1, alphas2, betas2))
      H0_vec[i] <- logitP(n_val = n, q = q,
                                params = as.numeric(eta_plus0),
                                delta = deltas, u = sob_H0[sub_H0[i],],
                                hyper = rbind(alphas1, betas1, alphas2, betas2))
    }
    
    ## for the unselected points, we just assume that their probabilities are less (greater)
    ## than the relevant order statistic if their anticipated probability based on the 
    ## linear approximations on the logit scale is less (greater) than the relevant order statistic
    H0_aug <- c(rep(c(min(H0_vec)-1), low1 - 1), H0_vec, rep(c(max(H0_vec)+1), m - high1))
    H1_aug <- c(rep(c(min(H1_vec)-1), low2 - 1), H1_vec, rep(c(max(H1_vec)+1), m - high2))
    if (lA <= sort(H1_aug)[mid2]){
      H1_max <- H1_vec; H0_max <- H0_vec
      sub_H1_max <- sub_H1; sub_H0_max <-sub_H0
      upper <- n
    }
    else{
      lower <- n
    }
  }
  
  n <- upper
  ## we don't need to reapproximate the same posterior probabilities if we have already explored
  ## this sample size
  if (n == n_high){
    
    ## do not need to do this; it is just so we have three sample sizes to construct each contour plot
    n <- ifelse(n0 == n_high, max(ceiling(1.1*n0), n0 + 5), max(ceiling(1.1*n1), n1 + 5))
    H1_vec <- NULL
    H0_vec <- NULL
    ## use Algorithm 2 to approximate the posterior probabilities at all m points
    for (i in 1:m){
      H1_vec[i] <- logitP(n_val = n, q = q,
                                  params = as.numeric(eta_plus1),
                                  delta = deltas, u = sob_H1[i,],
                                  hyper = rbind(alphas1, betas1, alphas2, betas2))
      H0_vec[i] <- logitP(n_val = n, q = q,
                                params = as.numeric(eta_plus0),
                                delta = deltas, u = sob_H0[i,],
                                hyper = rbind(alphas1, betas1, alphas2, betas2))
    }
    
    ## save the logits of the posterior probabilities
    H1_vec_n2 <- H1_vec
    H0_vec_n2 <- H0_vec
    n2 <- n
    
    n_final <- ifelse(n1 > n0, 1, 0)
    
    ## first component of list is (n, gamma, relevant order statistic of H1 probabilities, confirmatory type I error
    ## estimate, confirmatory power estimate, n0, n1, n2)
    ## components 2 to 6 of the list are only needed for the contour plot (they are the posterior probabilities under
    ## H1 and H0 for the sample sizes n0, n1, n2)
    ## we sort the sample sizes and the posterior probabilities to help with making the plot
    results <- list(c(get(paste0("n", n_final)), 1/(1+exp(-lA)), 
                      1/(1+exp(-as.numeric(sort(get(paste0("H1_vec_n", n_final)))[mid2]))),
                      as.numeric(mean(get(paste0("H0_vec_n", n_final)) >lA)), 
                      as.numeric(mean(get(paste0("H1_vec_n", n_final)) > lA)),
                      sort(c(n0, n1, n2))),
                    get(paste0("H0_vec_n", which.min(c(n0, n1, n2))-1)), get(paste0("H1_vec_n", which.min(c(n0, n1, n2))-1)), 
                    get(paste0("H0_vec_n", which(sort(c(n0, n1, n2))[2] == c(n0, n1, n2))-1)), get(paste0("H1_vec_n", which(sort(c(n0, n1, n2))[2] == c(n0, n1, n2))-1)),
                    get(paste0("H0_vec_n", which.max(c(n0, n1, n2))-1)), get(paste0("H1_vec_n", which.max(c(n0, n1, n2))-1)))
    
    if (contour == TRUE){
      return(results)
    }
    else{
      return(results[[1]][1:5])
    }
  }
  
  if (n == n_low){
    
    ## do not need to do this; it is just so we have three sample sizes to construct each contour plot
    n <- ifelse(n0 == n_low, min(floor(0.9*n0), n0 - 5), min(floor(0.9*n1), n1 - 5))
    H1_vec <- NULL
    H0_vec <- NULL
    ## use Algorithm 2 to approximate the posterior probabilities at all m points
    for (i in 1:m){
      H1_vec[i] <- logitP(n_val = n, q = q,
                                  params = as.numeric(eta_plus1),
                                  delta = deltas, u = sob_H1[i,],
                                  hyper = rbind(alphas1, betas1, alphas2, betas2))
      H0_vec[i] <- logitP(n_val = n, q = q,
                                params = as.numeric(eta_plus0),
                                delta = deltas, u = sob_H0[i,],
                                hyper = rbind(alphas1, betas1, alphas2, betas2))
    }
    
    ## save the logits of the posterior probabilities
    H1_vec_n2 <- H1_vec
    H0_vec_n2 <- H0_vec
    n2 <- n
    
    n_final <- ifelse(n1 < n0, 1, 0)
    
    ## first component of list is (n, gamma, relevant order statistic of H1 probabilities, confirmatory type I error
    ## estimate, confirmatory power estimate, n0, n1, n2)
    ## components 2 to 6 of the list are only needed for the contour plot (they are the posterior probabilities under
    ## H1 and H0 for the sample sizes n0, n1, n2)
    ## we sort the sample sizes and the posterior probabilities to help with making the plot
    results <- list(c(get(paste0("n", n_final)), 1/(1+exp(-lA)), 
                      1/(1+exp(-as.numeric(sort(get(paste0("H1_vec_n", n_final)))[mid2]))),
                      as.numeric(mean(get(paste0("H0_vec_n", n_final)) > lA)), 
                      as.numeric(mean(get(paste0("H1_vec_n", n_final)) > lA)),
                      sort(c(n0, n1, n2))),
                    get(paste0("H0_vec_n", which.min(c(n0, n1, n2))-1)), get(paste0("H1_vec_n", which.min(c(n0, n1, n2))-1)), 
                    get(paste0("H0_vec_n", which(sort(c(n0, n1, n2))[2] == c(n0, n1, n2))-1)), get(paste0("H1_vec_n", which(sort(c(n0, n1, n2))[2] == c(n0, n1, n2))-1)),
                    get(paste0("H0_vec_n", which.max(c(n0, n1, n2))-1)), get(paste0("H1_vec_n", which.max(c(n0, n1, n2))-1)))
    
    if (contour == TRUE){
      return(results)
    }
    else{
      return(results[[1]][1:5])
    }
  }
  
  ## otherwise, we approximate the remaining posterior probabilities that were not
  ## selected at the final sample size recommendation
  H1_vec <- NULL
  H0_vec <- NULL
  remain_H1 <- subset(1:m, ! 1:m %in% sub_H1_max)
  remain_H0 <- subset(1:m, ! 1:m %in% sub_H0_max)
  for (i in remain_H1){
    H1_vec[i] <- logitP(n_val = n, q = q,
                                params = as.numeric(eta_plus1),
                                delta = deltas, u = sob_H1[i,],
                                hyper = rbind(alphas1, betas1, alphas2, betas2))
  }
  for (i in remain_H0){
    H0_vec[i] <- logitP(n_val = n, q = q,
                              params = as.numeric(eta_plus0),
                              delta = deltas, u = sob_H0[i,],
                              hyper = rbind(alphas1, betas1, alphas2, betas2))
  }
  H1_vec[sub_H1_max] <- H1_max; H0_vec[sub_H0_max] <- H0_max
  
  ## save the logits of the posterior probabilities
  H1_vec_n2 <- H1_vec
  H0_vec_n2 <- H0_vec
  n2 <- n
  
  ## first component of list is (n, gamma, relevant order statistic of H1 probabilities, confirmatory type I error
  ## estimate, confirmatory power estimate, n0, n1, n2)
  ## components 2 to 6 of the list are only needed for the contour plot (they are the posterior probabilities under
  ## H1 and H0 for the sample sizes n0, n1, n2)
  ## we sort the sample sizes and the posterior probabilities to help with making the plot
  
  results <- list(c(n, 1/(1 + exp(-lA)),
                    1/(1 + exp(-as.numeric(sort(H1_vec)[mid2]))),
                    as.numeric(mean(H0_vec > lA)), as.numeric(mean(H1_vec > lA)),
                    sort(c(n0, n1, n2))),
                  get(paste0("H0_vec_n", which.min(c(n0, n1, n2))-1)), get(paste0("H1_vec_n", which.min(c(n0, n1, n2))-1)), 
                  get(paste0("H0_vec_n", which(sort(c(n0, n1, n2))[2] == c(n0, n1, n2))-1)), get(paste0("H1_vec_n", which(sort(c(n0, n1, n2))[2] == c(n0, n1, n2))-1)),
                  get(paste0("H0_vec_n", which.max(c(n0, n1, n2))-1)), get(paste0("H1_vec_n", which.max(c(n0, n1, n2))-1)))
  
  if (contour == TRUE){
    return(results)
  }
  else{
    return(results[[1]][1:5])
  }
}