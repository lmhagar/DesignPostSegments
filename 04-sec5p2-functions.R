## BEGIN SETUP ##

## load necessary packages
require(qrng)
require(foreach)
require(doParallel)
require(doSNOW)
require(qrng)

## function to calculate the logit of the posterior probability
## u is point from Sobol' sequence, params are the draw from Psi_0
## or Psi_1, deltas is the interval c(delta_L, delta_U), n_val is the 
## sample size presently explored, hyper is a list of the hyperparameters, 
## and q is the constant for imbalanced sample size determination
logitP <- function(u, params, delta, n_val, q, 
                       hyper = list(rep(0,3), 0.01*diag(3), 1, 1)){
  
  ## get coefficient values for data generation
  beta0 <- params[1]
  beta1 <- params[2]
  beta2 <- params[3]
  
  ## get values to generate covariates and error terms
  mu2 <- params[4]
  sigma2 <- params[5]
  sigmae <- params[6]
  
  ## save hyperparameters
  mu0 <- hyper[[1]]
  lambda0 <- hyper[[2]]
  a0 <- hyper[[3]]
  b0 <- hyper[[4]]
  
  ## obtain sample sizes for groups A (n1) and B (n2)
  n1 <- q*n_val
  n2 <- n_val

  ## generate group-specific sample means for x2
  n <- n1 + n2
  sum_x1 <- n1
  x21_bar <- qnorm(u[1], mu2, sigma2/sqrt(n1))
  x22_bar <- qnorm(u[2], mu2, sigma2/sqrt(n2))
  
  ## generate sample variance for x2 (both groups)
  ss_x2 <- sigma2^2*qchisq(u[3], n1 + n2 - 1)
  
  ## use algebra to get sufficient statistics
  sum_x2 <- n1*x21_bar + n2*x22_bar
  sum_x1x2 <- n1*x21_bar
  sum_x22 <- ss_x2 + sum_x2^2/n
  
  ## generate group-specific sample means for error terms
  eps21_bar <- qnorm(u[4], 0, sigmae/sqrt(n1))
  eps22_bar <- qnorm(u[5], 0, sigmae/sqrt(n2))
  
  ## use algebra to get sufficient statistics
  sum_eps <- n1*eps21_bar + n2*eps22_bar
  sum_x1eps <- n1*eps21_bar
  
  ## use Barlett decomposition to get remaining components of the 
  ## covariance matrix for x2 and varepsilon
  ss_epsx2 <- sigma2*sigmae*sqrt(qchisq(u[3], n1 + n2 - 1))*qnorm(u[6])
  ss_eps <- sigmae^2*(qchisq(u[7], n - 2) + qnorm(u[6])^2)
  sum_epsx2 <- ss_epsx2 + sum_x2*sum_eps/(n1 + n2)
  sum_eps2 <- ss_eps + sum_eps^2/n
  
  ## use algebra to get sufficient statistics
  sum_y <- beta0*(n1 + n2) + beta1*n1 + beta2*sum_x2 + sum_eps
  sum_yx1 <- beta0*n1 + beta1*n1 + beta2*n1*x21_bar + n1*eps21_bar
  sum_yx2 <- beta0*sum_x2 + beta1*sum_x1x2 + beta2*sum_x22 + sum_epsx2
  sum_y2 <- n*beta0^2 + beta1^2*n1 + beta2^2*sum_x22 + sum_eps2 +
    2*beta0*beta1*n1 + 2*beta0*beta2*sum_x2 + 2*beta0*sum_eps +
    2*beta1*beta2*sum_x1x2 + 2*beta1*sum_x1eps + 2*beta2*sum_epsx2
  
  ## compute the following statistics to get the location and
  ## scale parameter of the t-distribution posterior
  XtX <- rbind(c(n1 + n2, n1, sum_x2),
               c(n1, n1, sum_x1x2),
               c(sum_x2, sum_x1x2, sum_x22))
  
  Xty <- c(sum_y, sum_yx1, sum_yx2)
  lambdaN <- XtX + lambda0
  muN <- solve(lambdaN)%*%t(t(mu0)%*%lambda0 + Xty)
  
  aN <- a0 + 0.5*(n1 + n2)
  bN <- b0 + 0.5*(sum_y2 + t(mu0)%*%lambda0%*%mu0 - t(muN)%*%lambdaN%*%muN)
  
  ## compute posterior probability
  realP <- pt((delta[2]-muN[2])/sqrt(solve(lambdaN)[2,2]*bN/aN), 2*aN) - pt((delta[1]-muN[2])/sqrt(solve(lambdaN)[2,2]*bN/aN), 2*aN)
  
  ## slight perturbation of the posterior probability if it is too close to 0 or 1
  ## this ensures that the logit of the posterior probability is finite.
  if (realP > 1 - 10^(-7)){
    realP <- pt((delta[2]-muN[2])/sqrt(solve(lambdaN)[2,2]*bN/aN), 2*aN, lower.tail = FALSE) + pt((delta[1]-muN[2])/sqrt(solve(lambdaN)[2,2]*bN/aN), 2*aN)
    logitP <- log(realP) - log(1 - realP)
    logitP <- -1*logitP
  }
  else if (realP < .Machine$double.eps){
    realP <- .Machine$double.eps
    logitP <- log(realP) - log(1 - realP)
  }
  else{
    logitP <- log(realP) - log(1 - realP)
  }
  
  ## return the logit of the posterior probability
  return(logitP)
  
}

## function where we return what we need to get the contour plots
## Code to implement Algorithm 3 for the illustrative example. We explore sample sizes
## and return the optimal design (i.e., the (n, gamma) combination)
findGamma <- function(eta_plus1, eta_plus0, pwr, typeI, deltas, q, 
                      hypers, m = 8192, m0 = 512,
                      seed_H1 = 1, seed_H0 = 2, upper_init = 1000, prng = FALSE, 
                      contour = FALSE, b1.low = 9, b1.high = 12){
  
  ## eta_plus1: parameter draws from Psi_1
  ## eta_plus0: parameter draws from Psi_0
  ## pwr: target power (1 - beta)
  ## typeI: desired type I error rate (alpha)
  ## deltas: interval to define H1
  ## q: constant for imbalanced sample size determination
  ## hypers: hyperparamters for conjugate prior
  ## m: number of simulation repetitions
  ## m_0: number of repetitions to explore sampling distribution segments
  ## seed_H1: seed to generate Sobol' sequence for Psi_1
  ## seed_H0: seed to generate Sobol' sequence for Psi_0
  ## upper_init: maximum sample size to explore (not an issue if this is too small)
  ## prng: TRUE is pseudorandom sequence is used (for comparison)
  ## contour: when TRUE, returns logits of posterior probabilities that construct contour plot
  ## b1.low: lower endpoint for uniform distribution for beta1 (specific to this Psi1)
  ## b1.high: upper endpoint for uniform distribution for beta1 (specific to this Psi1)
  
  ## generate Sobol' sequence for each sampling distribution
  ## pseudorandom option available for comparison
  if (prng == FALSE){
    sob_H1 <- sobol(m, d = 8, randomize = "digital.shift", seed = seed_H1)
    sob_H0 <- sobol(m, d = 7, randomize = "digital.shift", seed = seed_H0)
  }
  else{
    set.seed(seed_H1); sob_H1 <- matrix(runif(8*m), ncol = 8)
    
    set.seed(seed_H0); sob_H0 <- matrix(runif(7*m), ncol = 7)
  }
  
  eta_plus1[,2] <- b1.low + (b1.high - b1.low)*sob_H1[,8]
  
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
    ## implement Algorithm 3 for each of the first m0 points
    for (i in 1:m0){
      H1_vec[i] <- logitP(n_val = n, q = q,
                                  params = as.numeric(eta_plus1[i,]),
                                  delta = deltas, u = sob_H1[i,1:7],
                                  hyper = hypers)
      H0_vec[i] <- logitP(n_val = n, q = q,
                                params = as.numeric(eta_plus0[i,]),
                                delta = deltas, u = sob_H0[i,],
                                hyper = hypers)
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
                                  params = as.numeric(eta_plus1[i,]),
                                  delta = deltas, u = sob_H1[i,1:7],
                                  hyper = hypers)
      H0_vec[i] <- logitP(n_val = n, q = q,
                                params = as.numeric(eta_plus0[i,]),
                                delta = deltas, u = sob_H0[i,],
                                hyper = hypers)
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
  
  ## construct the linear approximations under H0 and H1
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
    if (n == n_min){
      n <- n - 1
    }
    if (n == n_max){
      n <- n + 1
    }
    H1_vec <- NULL
    H0_vec <- NULL
    ## we next run Algorithm 3 at that initialize sample size for the first m0 points
    for (i in 1:m0){
      H1_vec[i] <- logitP(n_val = n, q = q,
                                  params = as.numeric(eta_plus1[i,]),
                                  delta = deltas, u = sob_H1[i,1:7],
                                  hyper = hypers)
      H0_vec[i] <- logitP(n_val = n, q = q,
                                params = as.numeric(eta_plus0[i,]),
                                delta = deltas, u = sob_H0[i,],
                                hyper = hypers)
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
                                params = as.numeric(eta_plus1[i,]),
                                delta = deltas, u = sob_H1[i,1:7],
                                hyper = hypers)
    H0_vec[i] <- logitP(n_val = n, q = q,
                              params = as.numeric(eta_plus0[i,]),
                              delta = deltas, u = sob_H0[i,],
                              hyper = hypers)
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
                                  params = as.numeric(eta_plus1[i,]),
                                  delta = deltas, u = sob_H1[i,1:7],
                                  hyper = hypers)
      H0_vec[i] <- logitP(n_val = n, q = q,
                                params = as.numeric(eta_plus0[i,]),
                                delta = deltas, u = sob_H0[i,],
                                hyper = hypers)
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
                                  params = as.numeric(eta_plus1[i,]),
                                  delta = deltas, u = sob_H1[i,1:7],
                                  hyper = hypers)
      H0_vec[i] <- logitP(n_val = n, q = q,
                                params = as.numeric(eta_plus0[i,]),
                                delta = deltas, u = sob_H0[i,],
                                hyper = hypers)
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
                                  params = as.numeric(eta_plus1[sub_H1[i],]),
                                  delta = deltas, u = sob_H1[sub_H1[i],],
                                  hyper = hypers)
      H0_vec[i] <- logitP(n_val = n, q = q,
                                params = as.numeric(eta_plus0[sub_H0[i],]),
                                delta = deltas, u = sob_H0[sub_H0[i],],
                                hyper = hypers)
    }
    
    ## for the unselected points, we just assume that their probabilities are less (greater)
    ## than the relevant order statistic if their anticipated probability based on the 
    ## linear approximations on the logit scale is less (greater) than the relevant order statistic
    H0_aug <- c(rep(c(min(H0_vec)-1), low1 - 1), H0_vec, rep(c(max(H0_vec)+1), m - high1))
    H1_aug <- c(rep(c(min(H1_vec)-1), low2 - 1), H1_vec, rep(c(max(H1_vec)+1), m - high2))
    sub_H1_max <- NULL; sub_H0_max <- NULL
    if (sort(H0_aug)[mid1] <= sort(H1_aug)[mid2]){
      H1_max <- H1_vec; H0_max <- H0_vec
      sub_H1_max <- sub_H1; sub_H0_max <- sub_H0
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
                                  params = as.numeric(eta_plus1[i,]),
                                  delta = deltas, u = sob_H1[i,1:7],
                                  hyper = hypers)
      H0_vec[i] <- logitP(n_val = n, q = q,
                                params = as.numeric(eta_plus0[i,]),
                                delta = deltas, u = sob_H0[i,],
                                hyper = hypers)
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
                      sort(c(n0, n1, n2)), c(n0, n1, n2)),
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
    ## use Algorithm 3 to approximate the posterior probabilities at all m points
    for (i in 1:m){
      H1_vec[i] <- logitP(n_val = n, q = q,
                                  params = as.numeric(eta_plus1[i,]),
                                  delta = deltas, u = sob_H1[i,1:7],
                                  hyper = hypers)
      H0_vec[i] <- logitP(n_val = n, q = q,
                                params = as.numeric(eta_plus0[i,]),
                                delta = deltas, u = sob_H0[i,],
                                hyper = hypers)
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
                      sort(c(n0, n1, n2)), c(n0, n1, n2)),
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
                                params = as.numeric(eta_plus1[i,]),
                                delta = deltas, u = sob_H1[i,1:7],
                                hyper = hypers)
  }
  for (i in remain_H0){
    H0_vec[i] <- logitP(n_val = n, q = q,
                              params = as.numeric(eta_plus0[i,]),
                              delta = deltas, u = sob_H0[i,],
                              hyper = hypers)
  }
  H1_vec[sub_H1_max] <- H1_max; H0_vec[sub_H0_max] <- H0_max
  
  ## save the logits of the posterior probabilities
  H1_vec_n2 <- H1_vec
  H0_vec_n2 <- H0_vec
  n2 <- n
  
  ## first component of list is (n, gamma, relevant order statistic of H1 probabilities, confirmatory type I error
  ## estimate, confirmatory power estimate, n0, n1, n2)
  ## components 2 to 6 of the list are only needed for the contour plot (they are the posterior probabilities under
  ## H0 and H1 for the sample sizes n0, n1, n2)
  ## we sort the sample sizes and the posterior probabilities to help with making the plot
  
  results <- list(c(n, 1/(1 + exp(-as.numeric(sort(H0_vec)[mid1]))),
                    1/(1 + exp(-as.numeric(sort(H1_vec)[mid2]))),
                    as.numeric(mean(H0_vec > sort(H0_vec)[mid1])), as.numeric(mean(H1_vec > sort(H0_vec)[mid1])),
                    sort(c(n0, n1, n2)), c(n0, n1, n2)),
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
                          hypers, m = 8192, m0 = 512,
                          seed_H1 = 1, seed_H0 = 2, upper_init = 1000,
                          b1.low = 9, b1.high = 12){
  
  ## inputs are the same as findGamma(), but m0 is not used
  
  ## generate Sobol' sequence for each sampling distribution
  sob_H1 <- sobol(m, d = 8, randomize = "digital.shift", seed = seed_H1)
  sob_H0 <- sobol(m, d = 7, randomize = "digital.shift", seed = seed_H0)
  
  eta_plus1[,2] <- b1.low + (b1.high - b1.low)*sob_H1[,8]
  
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
    ## implement Algorithm 3 for each of the first m0 points
    for (i in 1:m){
      H1_vec[i] <- logitP(n_val = n, q = q,
                                  params = as.numeric(eta_plus1[i,]),
                                  delta = deltas, u = sob_H1[i,1:7],
                                  hyper = hypers)
      H0_vec[i] <- logitP(n_val = n, q = q,
                                params = as.numeric(eta_plus0[i,]),
                                delta = deltas, u = sob_H0[i,],
                                hyper = hypers)
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
                                  params = as.numeric(eta_plus1[i,]),
                                  delta = deltas, u = sob_H1[i,1:7],
                                  hyper = hypers)
      H0_vec[i] <- logitP(n_val = n, q = q,
                                params = as.numeric(eta_plus0[i,]),
                                delta = deltas, u = sob_H0[i,],
                                hyper = hypers)
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
                      hypers, m = 8192, m0 = 512,
                      seed_H1 = 1, seed_H0 = 2, upper_init = 1000, prng = FALSE, 
                      contour = FALSE, b1.low = 9, b1.high = 12){
  
  ## the inputs are the same as findGamma()
  
  ## however, the sampling distribution of posterior probabilities under H1 is always compared to
  ## the logit of typeI (for all sample sizes considered). Thus, the criterion for the type I error
  ## rate may not always be satisfied
  lA <- log(1 - typeI) - log(typeI)
  
  ## generate Sobol' sequence for each sampling distribution
  ## pseudorandom option available for comparison
  if (prng == FALSE){
    sob_H1 <- sobol(m, d = 8, randomize = "digital.shift", seed = seed_H1)
    sob_H0 <- sobol(m, d = 7, randomize = "digital.shift", seed = seed_H0)
  }
  else{
    set.seed(seed_H1); sob_H1 <- matrix(runif(8*m), ncol = 8)
    
    set.seed(seed_H0); sob_H0 <- matrix(runif(7*m), ncol = 7)
  }
  
  eta_plus1[,2] <- b1.low + (b1.high - b1.low)*sob_H1[,8]
  
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
    ## implement Algorithm 3 for each of the first m0 points
    for (i in 1:m0){
      H1_vec[i] <- logitP(n_val = n, q = q,
                                  params = as.numeric(eta_plus1[i,]),
                                  delta = deltas, u = sob_H1[i,1:7],
                                  hyper = hypers)
      H0_vec[i] <- logitP(n_val = n, q = q,
                                params = as.numeric(eta_plus0[i,]),
                                delta = deltas, u = sob_H0[i,],
                                hyper = hypers)
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
                                  params = as.numeric(eta_plus1[i,]),
                                  delta = deltas, u = sob_H1[i,1:7],
                                  hyper = hypers)
      H0_vec[i] <- logitP(n_val = n, q = q,
                                params = as.numeric(eta_plus0[i,]),
                                delta = deltas, u = sob_H0[i,],
                                hyper = hypers)
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
    if (n == n_min){
      n <- n - 1
    }
    if (n == n_max){
      n <- n + 1
    }
    H1_vec <- NULL
    H0_vec <- NULL
    ## we next run Algorithm 3 at that initialize sample size for the first m0 points
    for (i in 1:m0){
      H1_vec[i] <- logitP(n_val = n, q = q,
                                  params = as.numeric(eta_plus1[i,]),
                                  delta = deltas, u = sob_H1[i,1:7],
                                  hyper = hypers)
      H0_vec[i] <- logitP(n_val = n, q = q,
                                params = as.numeric(eta_plus0[i,]),
                                delta = deltas, u = sob_H0[i,],
                                hyper = hypers)
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
  
  ## compute the posterior probabilities using Algorithm 3 for all m points for
  ## each Sobol' sequence at the initial sample size n(0)
  n_init <- n
  H1_vec <- NULL
  H0_vec <- NULL
  for (i in 1:m){
    H1_vec[i] <- logitP(n_val = n, q = q,
                                params = as.numeric(eta_plus1[i,]),
                                delta = deltas, u = sob_H1[i,1:7],
                                hyper = hypers)
    H0_vec[i] <- logitP(n_val = n, q = q,
                              params = as.numeric(eta_plus0[i,]),
                              delta = deltas, u = sob_H0[i,],
                              hyper = hypers)
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
    ## use Algorithm 3 to approximate the posterior probabilities at all m points
    for (i in 1:m){
      H1_vec[i] <- logitP(n_val = n, q = q,
                                  params = as.numeric(eta_plus1[i,]),
                                  delta = deltas, u = sob_H1[i,1:7],
                                  hyper = hypers)
      H0_vec[i] <- logitP(n_val = n, q = q,
                                params = as.numeric(eta_plus0[i,]),
                                delta = deltas, u = sob_H0[i,],
                                hyper = hypers)
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
    ## use Algorithm 3 to approximate the posterior probabilities at all m points
    for (i in 1:m){
      H1_vec[i] <- logitP(n_val = n, q = q,
                                  params = as.numeric(eta_plus1[i,]),
                                  delta = deltas, u = sob_H1[i,1:7],
                                  hyper = hypers)
      H0_vec[i] <- logitP(n_val = n, q = q,
                                params = as.numeric(eta_plus0[i,]),
                                delta = deltas, u = sob_H0[i,],
                                hyper = hypers)
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
                                  params = as.numeric(eta_plus1[sub_H1[i],]),
                                  delta = deltas, u = sob_H1[sub_H1[i],],
                                  hyper = hypers)
      H0_vec[i] <- logitP(n_val = n, q = q,
                                params = as.numeric(eta_plus0[sub_H0[i],]),
                                delta = deltas, u = sob_H0[sub_H0[i],],
                                hyper = hypers)
    }
    
    ## for the unselected points, we just assume that their probabilities are less (greater)
    ## than the relevant order statistic if their anticipated probability based on the 
    ## linear approximations on the logit scale is less (greater) than the relevant order statistic
    H0_aug <- c(rep(c(min(H0_vec)-1), low1 - 1), H0_vec, rep(c(max(H0_vec)+1), m - high1))
    H1_aug <- c(rep(c(min(H1_vec)-1), low2 - 1), H1_vec, rep(c(max(H1_vec)+1), m - high2))
    sub_H1_max <- NULL; sub_H0_max <- NULL
    if (sort(H0_aug)[mid1] <= sort(H1_aug)[mid2]){
      H1_max <- H1_vec; H0_max <- H0_vec
      sub_H1_max <- sub_H1; sub_H0_max <- sub_H0
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
                                  params = as.numeric(eta_plus1[i,]),
                                  delta = deltas, u = sob_H1[i,1:7],
                                  hyper = hypers)
      H0_vec[i] <- logitP(n_val = n, q = q,
                                params = as.numeric(eta_plus0[i,]),
                                delta = deltas, u = sob_H0[i,],
                                hyper = hypers)
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
    results <- list(c(get(paste0("n", n_final)), 1/(1 + exp(-as.numeric(lA))), 
                      1/(1+exp(-as.numeric(sort(get(paste0("H1_vec_n", n_final)))[mid2]))),
                      as.numeric(mean(get(paste0("H0_vec_n", n_final)) > lA)), 
                      as.numeric(mean(get(paste0("H1_vec_n", n_final)) > lA)),
                      sort(c(n0, n1, n2)), c(n0, n1, n2)),
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
    ## use Algorithm 3 to approximate the posterior probabilities at all m points
    for (i in 1:m){
      H1_vec[i] <- logitP(n_val = n, q = q,
                                  params = as.numeric(eta_plus1[i,]),
                                  delta = deltas, u = sob_H1[i,1:7],
                                  hyper = hypers)
      H0_vec[i] <- logitP(n_val = n, q = q,
                                params = as.numeric(eta_plus0[i,]),
                                delta = deltas, u = sob_H0[i,],
                                hyper = hypers)
    }
    
    ## save the logits of the posterior probabilities
    H1_vec_n2 <- H1_vec
    H0_vec_n2 <- H0_vec
    n2 <- n
    
    n_final <- ifelse(n1 < n0, 1, 0)
    
    ## first component of list is (n, gamma, relevant order statistic of H1 probabilities, confirmatory type I error
    ## estimate, confirmatory power estimate, n^((0)), n^((1)), n^((2)))
    ## components 2 to 6 of the list are only needed for the contour plot (they are the posterior probabilities under
    ## H1 and H0 regions for the sample sizes n0, n1, n2)
    ## we sort the sample sizes and the posterior probabilities to help with making the plot
    results <- list(c(get(paste0("n", n_final)), 1/(1 + exp(-as.numeric(lA))), 
                      1/(1+exp(-as.numeric(sort(get(paste0("H1_vec_n", n_final)))[mid2]))),
                      as.numeric(mean(get(paste0("H0_vec_n", n_final)) > lA)), 
                      as.numeric(mean(get(paste0("H1_vec_n", n_final)) > lA)),
                      sort(c(n0, n1, n2)), c(n0, n1, n2)),
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
                                params = as.numeric(eta_plus1[i,]),
                                delta = deltas, u = sob_H1[i,1:7],
                                hyper = hypers)
  }
  for (i in remain_H0){
    H0_vec[i] <- logitP(n_val = n, q = q,
                              params = as.numeric(eta_plus0[i,]),
                              delta = deltas, u = sob_H0[i,],
                              hyper = hypers)
  }
  H1_vec[sub_H1_max] <- H1_max; H0_vec[sub_H0_max] <- H0_max
  
  ## save the logits of the posterior probabilities
  H1_vec_n2 <- H1_vec
  H0_vec_n2 <- H0_vec
  n2 <- n
  
  ## first component of list is (n, gamma, relevant order statistic of H1 probabilities, confirmatory type I error
  ## estimate, confirmatory power estimate, n0, n1, n2)
  ## components 2 to 6 of the list are only needed for the contour plot (they are the posterior probabilities under
  ## H0 and H1 for the sample sizes n0, n1, n2)
  ## we sort the sample sizes and the posterior probabilities to help with making the plot
  
  results <- list(c(n, 1/(1 + exp(-as.numeric(lA))),
                    1/(1 + exp(-as.numeric(sort(H1_vec)[mid2]))),
                    as.numeric(mean(H0_vec > lA)), as.numeric(mean(H1_vec > lA)),
                    sort(c(n0, n1, n2)), c(n0, n1, n2)),
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

## this is the same function as logitP; however, we do not account for
## the x2 covariate (modifications are commented later in this algorithm)
## also this function takes different hyperparameters (i.e., 2 dimensions
## instead of 3)
logitPNoCov <- function(u, params, delta, n_val, q, 
                        hyper = list(rep(0,2), 0.01*diag(2), 1, 1)){
  beta0 <- params[1]
  beta1 <- params[2]
  beta2 <- params[3]
  
  mu2 <- params[4]
  sigma2 <- params[5]
  sigmae <- params[6]
  
  mu0 <- hyper[[1]]
  lambda0 <- hyper[[2]]
  a0 <- hyper[[3]]
  b0 <- hyper[[4]]
  
  n1 <- n_val
  n2 <- q*n_val
  
  n <- n1 + n2
  sum_x1 <- n1
  x21_bar <- qnorm(u[1], mu2, sigma2/sqrt(n1))
  x22_bar <- qnorm(u[2], mu2, sigma2/sqrt(n2))
  
  ss_x2 <- sigma2^2*qchisq(u[3], n1 + n2 - 1)
  
  sum_x2 <- n1*x21_bar + n2*x22_bar
  sum_x1x2 <- n1*x21_bar
  sum_x22 <- ss_x2 + sum_x2^2/n
  
  eps21_bar <- qnorm(u[4], 0, sigmae/sqrt(n1))
  eps22_bar <- qnorm(u[5], 0, sigmae/sqrt(n2))
  
  sum_eps <- n1*eps21_bar + n2*eps22_bar
  sum_x1eps <- n1*eps21_bar
  
  ss_epsx2 <- sigma2*sigmae*sqrt(qchisq(u[3], n1 + n2 - 1))*qnorm(u[6])
  ss_eps <- sigmae^2*(qchisq(u[7], n - 2) + qnorm(u[6])^2)
  sum_epsx2 <- ss_epsx2 + sum_x2*sum_eps/(n1 + n2)
  sum_eps2 <- ss_eps + sum_eps^2/n
  
  sum_y <- beta0*(n1 + n2) + beta1*n1 + beta2*sum_x2 + sum_eps
  sum_yx1 <- beta0*n1 + beta1*n1 + beta2*n1*x21_bar + n1*eps21_bar
  sum_yx2 <- beta0*sum_x2 + beta1*sum_x1x2 + beta2*sum_x22 + sum_epsx2
  sum_y2 <- n*beta0^2 + beta1^2*n1 + beta2^2*sum_x22 + sum_eps2 +
    2*beta0*beta1*n1 + 2*beta0*beta2*sum_x2 + 2*beta0*sum_eps +
    2*beta1*beta2*sum_x1x2 + 2*beta1*sum_x1eps + 2*beta2*sum_epsx2
  
  ## XtX and Xty are updated from logitP()
  XtX <- rbind(c(n1 + n2, n1),
               c(n1, n1))
  
  Xty <- c(sum_y, sum_yx1)
  
  lambdaN <- XtX + lambda0
  muN <- solve(lambdaN)%*%t(t(mu0)%*%lambda0 + Xty)
  
  aN <- a0 + 0.5*(n1 + n2)
  bN <- b0 + 0.5*(sum_y2 + t(mu0)%*%lambda0%*%mu0 - t(muN)%*%lambdaN%*%muN)
  
  ## compute posterior probability given this normal approximation
  realP <- pt((delta[2]-muN[2])/sqrt(solve(lambdaN)[2,2]*bN/aN), 2*aN) - pt((delta[1]-muN[2])/sqrt(solve(lambdaN)[2,2]*bN/aN), 2*aN)
  
  ## slight perturbation of the posterior probability if it is too close to 0 or 1
  ## this ensures that the logit of the posterior probability is finite.
  if (realP > 1 - 10^(-7)){
    realP <- pt((delta[2]-muN[2])/sqrt(solve(lambdaN)[2,2]*bN/aN), 2*aN, lower.tail = FALSE) + pt((delta[1]-muN[2])/sqrt(solve(lambdaN)[2,2]*bN/aN), 2*aN)
    logitPNoCov <- log(realP) - log(1 - realP)
    logitPNoCov <- -1*logitPNoCov
  }
  else if (realP < .Machine$double.eps){
    realP <- .Machine$double.eps
    logitPNoCov <- log(realP) - log(1 - realP)
  }
  else{
    logitPNoCov <- log(realP) - log(1 - realP)
  }
  
  ## return the logit of the posterior probability
  return(logitPNoCov)
  
}

findGammaNoCov <- function(eta_plus1, eta_plus0, pwr, typeI, deltas, q, 
                           hypers, m = 8192, m0 = 512,
                           seed_H1 = 1, seed_H0 = 2, upper_init = 1000, prng = FALSE, 
                           contour = FALSE, b1.low = 9, b1.high = 12){
  
  ## generate Sobol' sequence for each sampling distribution
  ## pseudorandom option available for comparison
  if (prng == FALSE){
    sob_H1 <- sobol(m, d = 8, randomize = "digital.shift", seed = seed_H1)
    sob_H0 <- sobol(m, d = 7, randomize = "digital.shift", seed = seed_H0)
  }
  else{
    set.seed(seed_H1); sob_H1 <- matrix(runif(8*m), ncol = 8)
    
    set.seed(seed_H0); sob_H0 <- matrix(runif(7*m), ncol = 7)
  }
  
  eta_plus1[,2] <- b1.low + (b1.high - b1.low)*sob_H1[,8]
  
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
    ## implement Algorithm 3 for each of the first m0 points
    for (i in 1:m0){
      H1_vec[i] <- logitPNoCov(n_val = n, q = q,
                               params = as.numeric(eta_plus1[i,]),
                               delta = deltas, u = sob_H1[i,1:7],
                               hyper = hypers)
      H0_vec[i] <- logitPNoCov(n_val = n, q = q,
                               params = as.numeric(eta_plus0[i,]),
                               delta = deltas, u = sob_H0[i,],
                               hyper = hypers)
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
      H1_vec[i] <- logitPNoCov(n_val = n, q = q,
                               params = as.numeric(eta_plus1[i,]),
                               delta = deltas, u = sob_H1[i,1:7],
                               hyper = hypers)
      H0_vec[i] <- logitPNoCov(n_val = n, q = q,
                               params = as.numeric(eta_plus0[i,]),
                               delta = deltas, u = sob_H0[i,],
                               hyper = hypers)
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
    if (n == n_min){
      n <- n - 1
    }
    if (n == n_max){
      n <- n + 1
    }
    H1_vec <- NULL
    H0_vec <- NULL
    ## we next run Algorithm 3 at that initialize sample size for the first m0 points
    for (i in 1:m0){
      H1_vec[i] <- logitPNoCov(n_val = n, q = q,
                               params = as.numeric(eta_plus1[i,]),
                               delta = deltas, u = sob_H1[i,1:7],
                               hyper = hypers)
      H0_vec[i] <- logitPNoCov(n_val = n, q = q,
                               params = as.numeric(eta_plus0[i,]),
                               delta = deltas, u = sob_H0[i,],
                               hyper = hypers)
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
    H1_vec[i] <- logitPNoCov(n_val = n, q = q,
                             params = as.numeric(eta_plus1[i,]),
                             delta = deltas, u = sob_H1[i,1:7],
                             hyper = hypers)
    H0_vec[i] <- logitPNoCov(n_val = n, q = q,
                             params = as.numeric(eta_plus0[i,]),
                             delta = deltas, u = sob_H0[i,],
                             hyper = hypers)
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
      H1_vec[i] <- logitPNoCov(n_val = n, q = q,
                               params = as.numeric(eta_plus1[i,]),
                               delta = deltas, u = sob_H1[i,1:7],
                               hyper = hypers)
      H0_vec[i] <- logitPNoCov(n_val = n, q = q,
                               params = as.numeric(eta_plus0[i,]),
                               delta = deltas, u = sob_H0[i,],
                               hyper = hypers)
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
      H1_vec[i] <- logitPNoCov(n_val = n, q = q,
                               params = as.numeric(eta_plus1[i,]),
                               delta = deltas, u = sob_H1[i,1:7],
                               hyper = hypers)
      H0_vec[i] <- logitPNoCov(n_val = n, q = q,
                               params = as.numeric(eta_plus0[i,]),
                               delta = deltas, u = sob_H0[i,],
                               hyper = hypers)
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
      H1_vec[i] <- logitPNoCov(n_val = n, q = q,
                               params = as.numeric(eta_plus1[sub_H1[i],]),
                               delta = deltas, u = sob_H1[sub_H1[i],],
                               hyper = hypers)
      H0_vec[i] <- logitPNoCov(n_val = n, q = q,
                               params = as.numeric(eta_plus0[sub_H0[i],]),
                               delta = deltas, u = sob_H0[sub_H0[i],],
                               hyper = hypers)
    }
    
    ## for the unselected points, we just assume that their probabilities are less (greater)
    ## than the relevant order statistic if their anticipated probability based on the 
    ## linear approximations on the logit scale is less (greater) than the relevant order statistic
    H0_aug <- c(rep(c(min(H0_vec)-1), low1 - 1), H0_vec, rep(c(max(H0_vec)+1), m - high1))
    H1_aug <- c(rep(c(min(H1_vec)-1), low2 - 1), H1_vec, rep(c(max(H1_vec)+1), m - high2))
    sub_H1_max <- NULL; sub_H0_max <- NULL
    if (sort(H0_aug)[mid1] <= sort(H1_aug)[mid2]){
      H1_max <- H1_vec; H0_max <- H0_vec
      sub_H1_max <- sub_H1; sub_H0_max <- sub_H0
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
      H1_vec[i] <- logitPNoCov(n_val = n, q = q,
                               params = as.numeric(eta_plus1[i,]),
                               delta = deltas, u = sob_H1[i,1:7],
                               hyper = hypers)
      H0_vec[i] <- logitPNoCov(n_val = n, q = q,
                               params = as.numeric(eta_plus0[i,]),
                               delta = deltas, u = sob_H0[i,],
                               hyper = hypers)
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
                      sort(c(n0, n1, n2)), c(n0, n1, n2)),
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
    ## use Algorithm 3 to approximate the posterior probabilities at all m points
    for (i in 1:m){
      H1_vec[i] <- logitPNoCov(n_val = n, q = q,
                               params = as.numeric(eta_plus1[i,]),
                               delta = deltas, u = sob_H1[i,1:7],
                               hyper = hypers)
      H0_vec[i] <- logitPNoCov(n_val = n, q = q,
                               params = as.numeric(eta_plus0[i,]),
                               delta = deltas, u = sob_H0[i,],
                               hyper = hypers)
    }
    
    ## save the logits of the posterior probabilities
    H1_vec_n2 <- H1_vec
    H0_vec_n2 <- H0_vec
    n2 <- n
    
    n_final <- ifelse(n1 < n0, 1, 0)
    
    ## first component of list is (n, gamma, relevant order statistic of H1 probabilities, confirmatory type I error
    ## estimate, confirmatory power estimate, n0, n1, n2)
    ## components 2 to 6 of the list are only needed for the contour plot (they are the posterior probabilities under
    ## H0 and H1 for the sample sizes n0, n1, n2)
    ## we sort the sample sizes and the posterior probabilities to help with making the plot
    results <- list(c(get(paste0("n", n_final)), 1/(1+exp(-as.numeric(sort(get(paste0("H0_vec_n", n_final)))[mid1]))), 
                      1/(1+exp(-as.numeric(sort(get(paste0("H1_vec_n", n_final)))[mid2]))),
                      as.numeric(mean(get(paste0("H0_vec_n", n_final)) > sort(get(paste0("H0_vec_n", n_final)))[mid1])), 
                      as.numeric(mean(get(paste0("H1_vec_n", n_final)) > sort(get(paste0("H0_vec_n", n_final)))[mid1])),
                      sort(c(n0, n1, n2)), c(n0, n1, n2)),
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
    H1_vec[i] <- logitPNoCov(n_val = n, q = q,
                             params = as.numeric(eta_plus1[i,]),
                             delta = deltas, u = sob_H1[i,1:7],
                             hyper = hypers)
  }
  for (i in remain_H0){
    H0_vec[i] <- logitPNoCov(n_val = n, q = q,
                             params = as.numeric(eta_plus0[i,]),
                             delta = deltas, u = sob_H0[i,],
                             hyper = hypers)
  }
  H1_vec[sub_H1_max] <- H1_max; H0_vec[sub_H0_max] <- H0_max
  
  ## save the logits of the posterior probabilities
  H1_vec_n2 <- H1_vec
  H0_vec_n2 <- H0_vec
  n2 <- n
  
  ## first component of list is (n, gamma, relevant order statistic of H1 probabilities, confirmatory type I error
  ## estimate, confirmatory power estimate, n0, n1, n2)
  ## components 2 to 6 of the list are only needed for the contour plot (they are the posterior probabilities under
  ## H0 and H1 for the sample sizes n0, n1, n2)
  ## we sort the sample sizes and the posterior probabilities to help with making the plot
  
  results <- list(c(n, 1/(1 + exp(-as.numeric(sort(H0_vec)[mid1]))),
                    1/(1 + exp(-as.numeric(sort(H1_vec)[mid2]))),
                    as.numeric(mean(H0_vec > sort(H0_vec)[mid1])), as.numeric(mean(H1_vec > sort(H0_vec)[mid1])),
                    sort(c(n0, n1, n2)), c(n0, n1, n2)),
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


