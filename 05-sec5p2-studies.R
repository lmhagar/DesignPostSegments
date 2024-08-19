## BEGIN SETUP ##

## load necessary packages
require(qrng)
require(foreach)
require(doParallel)
require(doSNOW)
require(qrng)

## set up hyperparameters and draws from Psi1 and Psi0 for numerical studies
## we overwrite the beta1 value in FindGamma()
eta_plus1s <- c(3 - 115*0.25, 10, 0.25, 115, 14.5, sqrt(10.7^2 - 0.25^2*(14.5^2)))
eta_plus1s <- matrix(rep(eta_plus1s, 4096), nrow = 4096, byrow = TRUE)

eta_plus0s <- c(3 - 115*0.25, 5, 0.25, 115, 14.5, sqrt(10.7^2 - 0.25^2*(14.5^2)))
eta_plus0s <- matrix(rep(eta_plus0s, 4096), nrow = 4096, byrow = TRUE)
hyper_ex <- list(rep(0,3), 0.01*diag(3), 1, 1)


## set up parallelization with 1000 repetitions and 72 cores
cores=detectCores()
cl <- makeSOCKcluster(cores[1]-1)

registerDoSNOW(cl)
pb <- txtProgressBar(max = 1000, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

## repeat the sample size calculation for the illustrative example 1000 times
## with the different Sobol' sequences
temp <- foreach(i=1:1000, .combine='rbind', .packages = c("qrng"),
                .options.snow=opts, .errorhandling = "remove") %dopar% {
                  temp_res <- findGamma(eta_plus1 = eta_plus1s, eta_plus0 = eta_plus0s, 
                                        pwr = 0.8, typeI = 0.05, deltas = c(5, Inf), q = 2, 
                                        hypers = hyper_ex, m = 4096, m0 = 64,
                                        seed_H1 = i, seed_H0 = i + 1000, contour = FALSE)
                  
                  c(i, as.numeric(unlist(temp_res)))
                }

## output the optimal design and summary for each simulation repetition
write.csv(temp, "sema_unif_linear_4096_summary.csv", row.names = FALSE)

## get 95% bootstrap CI for sample size n_B
quantile(temp[,2], c(0.025, 0.975))
## repeat for gamma
round(quantile(temp[,3], c(0.025, 0.975)),4)

## repeat the 1000 sample size calculations with the SAME Sobol' sequences but
## explore the entire hypercube and each sample size considered
tempFull <- foreach(i=1:1000, .combine='rbind', .packages = c("qrng"),
                    .options.snow=opts, .errorhandling = "remove") %dopar% {
                      temp_res <- findGammaFull(eta_plus1 = eta_plus1s, eta_plus0 = eta_plus0s, 
                                                pwr = 0.8, typeI = 0.05, deltas = c(5, Inf), q = 2, 
                                                hypers = hyper_ex, m = 4096, m0 = 64,
                                                seed_H1 = i, seed_H0 = i + 1000)
                      
                      c(i, as.numeric(unlist(temp_res)))
                    }

write.csv(tempFull, "sema10_unif_4096_full_summary.csv", row.names = FALSE)

## check if (n, gamma) combinations are the same
sum(temp[,2] == tempFull[,2])
sum(temp[,3] == tempFull[,3])

## sample size calculation for fixed critical value gamma
tempFixed <- foreach(i=1:1000, .combine='rbind', .packages = c("qrng"),
                     .options.snow=opts, .errorhandling = "remove") %dopar% {
                       temp_res <- findn(eta_plus1 = eta_plus1s, eta_plus0 = eta_plus0s, 
                                         pwr = 0.8, typeI = 0.05, deltas = c(5, Inf), q = 2, 
                                         hypers = hyper_ex, m = 4096, m0 = 64,
                                         seed_H1 = i, seed_H0 = i + 1000, contour = FALSE)
                       
                       c(i, as.numeric(unlist(temp_res)))
                     }

write.csv(tempFixed, "lin_4096_fixed_summary.csv", row.names = FALSE)

mean(tempFixed[,2]) ## recommended sample sizes is 35

## comparison precision of Sobol' sequences with precision of 
## longer pseudorandom ones
mm <-35000
eta_plus0s <- c(3 - 115*0.25, 5, 0.25, 115, 14.5, sqrt(10.7^2 - 0.25^2*(14.5^2)))
eta_plus0s <- matrix(rep(eta_plus0s, mm), nrow = mm, byrow = TRUE)
hyper_ex <- list(rep(0,3), 0.01*diag(3), 1, 1)

eta_plus1s <- c(3 - 115*0.25, 9, 0.25, 115, 14.5, sqrt(10.7^2 - 0.25^2*(14.5^2)),
                3 - 115*0.25, 12, 0.25, 115, 14.5, sqrt(10.7^2 - 0.25^2*(14.5^2)))
eta_plus1s <- matrix(rep(eta_plus1s, mm/2), nrow =mm, byrow = TRUE)

## repeat the sample size calculation for the gamma example 1000 times
## with the different pseudorandom sequences
tempPRNG <- foreach(i=1:1000, .combine='rbind', .packages = c("qrng"),
                    .options.snow=opts, .errorhandling = "remove") %dopar% {
                      temp_res <- findGamma(eta_plus1 = eta_plus1s, eta_plus0 = eta_plus0s, prng = TRUE,
                                            pwr = 0.8, typeI = 0.05, deltas = c(5, Inf), q = 2, 
                                            hypers = hyper_ex, m = 35000, m0 = 200,
                                            seed_H1 = i, seed_H0 = i + 1000, contour = FALSE)
                      
                      c(i, as.numeric(unlist(temp_res)))
                    }

## output the optimal design and summary for each simulation repetition
write.csv(tempPRNG, "sema_unif_linear_35K_summary_prng.csv", row.names = FALSE)

## get 95% bootstrap CI for sample size n_B
quantile(temp[,2], c(0.025, 0.975))
## repeat for gamma
round(quantile(temp[,3], c(0.025, 0.975)),4)

## now we assess the type I error rate and power for various simpler designs
## we go back to the initial m = 4096 value
eta_plus1s <- c(3 - 115*0.25, 10, 0.25, 115, 14.5, sqrt(10.7^2 - 0.25^2*(14.5^2)))
eta_plus1s <- matrix(rep(eta_plus1s, 4096), nrow = 4096, byrow = TRUE)

eta_plus0s <- c(3 - 115*0.25, 5, 0.25, 115, 14.5, sqrt(10.7^2 - 0.25^2*(14.5^2)))
eta_plus0s <- matrix(rep(eta_plus0s, 4096), nrow = 4096, byrow = TRUE)
hyper_ex <- list(rep(0,3), 0.01*diag(3), 1, 1)

## define parameters for this setting
n <- 32
lA <- log(0.95) - log(0.05)
m <- 4096
deltas <- c(5, Inf)
temp320p95 <- foreach(i=1:1000, .combine='rbind', .packages = c("qrng"),
                      .options.snow=opts, .errorhandling = "remove") %dopar% {
                        
                        sob_H1 <- sobol(m, d = 8, randomize = "digital.shift", seed = i)
                        sob_H0 <- sobol(m, d = 7, randomize = "digital.shift", seed = i + 1000)
                        
                        eta_plus1s[,2] <- 9 + (3)*sob_H1[,8]
                        
                        H1_vec <- NULL
                        H0_vec <- NULL
                        ## implement Algorithm 2 for all points
                        for (i in 1:nrow(sob_H1)){
                          H1_vec[i] <- logitP(n_val = n, q = 2,
                                              params = as.numeric(eta_plus1s[i,]),
                                              delta = deltas, u = sob_H1[i,1:7],
                                              hyper = hyper_ex)
                          H0_vec[i] <- logitP(n_val = n, q = 2,
                                              params = as.numeric(eta_plus0s[i,]),
                                              delta = deltas, u = sob_H0[i,],
                                              hyper = hyper_ex)
                        }
                        c(i, mean(H0_vec >= lA), mean(H1_vec >= lA))
                      }

write.csv(temp32p95, "lin_4096_32p95_summary.csv", row.names = FALSE)

## calculate type I error rate for this design
mean(temp32p95[,2])
## calculate power for this design
mean(temp32p95[,3])

## now consider n = 35 and gamma = 0.95
n <- 35
lA <- log(0.95) - log(0.05)
m <- 4096
deltas <- c(5, Inf)
temp350p95 <- foreach(i=1:1000, .combine='rbind', .packages = c("qrng"),
                      .options.snow=opts, .errorhandling = "remove") %dopar% {
                        
                        sob_H1 <- sobol(m, d = 8, randomize = "digital.shift", seed = i)
                        sob_H0 <- sobol(m, d = 7, randomize = "digital.shift", seed = i + 1000)
                        
                        eta_plus1s[,2] <- 9 + (3)*sob_H1[,8]
                        
                        H1_vec <- NULL
                        H0_vec <- NULL
                        ## implement Algorithm 2 for all points
                        for (i in 1:nrow(sob_H1)){
                          H1_vec[i] <- logitP(n_val = n, q = 2,
                                              params = as.numeric(eta_plus1s[i,]),
                                              delta = deltas, u = sob_H1[i,1:7],
                                              hyper = hyper_ex)
                          H0_vec[i] <- logitP(n_val = n, q = 2,
                                              params = as.numeric(eta_plus0s[i,]),
                                              delta = deltas, u = sob_H0[i,],
                                              hyper = hyper_ex)
                        }
                        c(i, mean(H0_vec >= lA), mean(H1_vec >= lA))
                      }

write.csv(temp35p95, "lin_4096_35p95_summary.csv", row.names = FALSE)

## calculate type I error rate for this design
mean(temp35p95[,2])
## calculate power for this design
mean(temp35p95[,3])


## repeat the sample size calculation for the gamma example 1000 times
## with the different Sobol' sequences (with contour = TRUE)
temp <- foreach(i=1:1000, .combine='rbind', .packages = c("qrng"),
                .options.snow=opts, .errorhandling = "remove") %dopar% {
                  temp_res <- findGamma(eta_plus1 = eta_plus1s, eta_plus0 = eta_plus0s, 
                                        pwr = 0.8, typeI = 0.05, deltas = c(5, Inf), q = 2, 
                                        hypers = hyper_ex, m = 4096, m0 = 128, upper_init = 100,
                                        seed_H1 = i, seed_H0 = i + 1000, contour = TRUE)
                  
                  c(i, as.numeric(unlist(temp_res)))
                }

## output the optimal design and summary for each simulation repetition
write.csv(temp[,1:12], "lin_sobol_4096_summary.csv", row.names = FALSE)

## output the red and green posterior probabilities for each of the three sample
## sizes n^((0)) = small, n^((1)) = mid, and n^((2)) = large
write.csv(temp[,13:4108], "lin_sobol_4096_H0_small.csv", row.names = FALSE)
write.csv(temp[,4109:8204], "lin_sobol_4096_H1_small.csv", row.names = FALSE)
write.csv(temp[,8205:12300], "lin_sobol_4096_H0_mid.csv", row.names = FALSE)
write.csv(temp[,12301:16396], "lin_sobol_4096_H1_mid.csv", row.names = FALSE)
write.csv(temp[,16397:20492], "lin_sobol_4096_H0_large.csv", row.names = FALSE)
write.csv(temp[,20493:24588], "lin_sobol_4096_H1_large.csv", row.names = FALSE)

## contour plot confirmatory simulations with 81920 repetitions
resamps <- 81920
pb <- txtProgressBar(max = resamps, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

## get parameter values from Psi1 and Psi0 to generate data from the linear model
samps <- seq(25, 45, 1)

delta_L <- 5
delta_U <- Inf
deltas <- c(delta_L, delta_U)

eta_plus1s <- c(3 - 115*0.25, 10, 0.25, 115, 14.5, sqrt(10.7^2 - 0.25^2*(14.5^2)))
eta_plus1s <- matrix(rep(eta_plus1s, resamps), nrow = resamps, byrow = TRUE)
## overwrite the beta1 value
eta_plus1s[,2] <- runif(resamps)*3 + 9

eta_plus0s <- c(3 - 115*0.25, 5, 0.25, 115, 14.5, sqrt(10.7^2 - 0.25^2*(14.5^2)))
eta_plus0s <- matrix(rep(eta_plus0s, resamps), nrow = resamps, byrow = TRUE)
hyper_ex <- list(rep(0,3), 0.01*diag(3), 1, 1)

for (k in 1:length(samps)){
  probs_in <- foreach(j=1:resamps, .combine='rbind', .packages = c("rjags", "coda"),
                      .options.snow=opts, .errorhandling = "remove") %dopar% {
                        
                        set.seed(k*resamps + j)
                        
                        temp <- logitP(n_val = samps[k], q = q,
                                       params = eta_plus1s[j,],
                                       delta = deltas, u = runif(7),
                                       hyper = hyper_ex)
                        
                        1/(1 + exp(-temp))
                        
                      }
  
  ## output results to .csv file
  write.csv(probs_in, paste0("probs_H1_lin_", samps[k], ".csv"), row.names = FALSE)
}

## get parameter values from design priors to generate
## the Weibull data for the red region
params <- rbind(eta_plus0s, eta_plus0s, eta_plus0s, eta_plus0s, eta_plus0s)

for (k in 1:length(samps)){
  
  probs_out <- foreach(j=1:resamps, .combine='rbind', .packages = c("rjags", "coda"),
                       .options.snow=opts, .errorhandling = "remove") %dopar% {
                         
                         set.seed(k*resamps + j)
                         
                         temp <- logitP(n_val = samps[k], q = q,
                                        params = eta_plus0s[j,],
                                        delta = deltas, u = runif(7),
                                        hyper = hyper_ex)
                         
                         1/(1 + exp(-temp))
                       }
  
  ## output results to .csv file
  write.csv(probs_out, paste0("probs_H0_lin_", samps[k], ".csv"), row.names = FALSE)
}

## create the two contour plots for Figure D.1
n_lows <- as.numeric(read.csv("lin_sobol_4096_summary.csv")[,7])
n_mids <- as.numeric(read.csv("lin_sobol_4096_summary.csv")[,8])
n_highs <- as.numeric(read.csv("lin_sobol_4096_summary.csv")[,9])

H1_mids <- read.csv("lin_sobol_4096_H1_mid.csv")
H0_mids <- read.csv("lin_sobol_4096_H0_mid.csv")
H1_lows <- read.csv("lin_sobol_4096_H1_small.csv")
H0_lows <- read.csv("lin_sobol_4096_H0_small.csv")
H1_highs <- read.csv("lin_sobol_4096_H1_large.csv")
H0_highs <- read.csv("lin_sobol_4096_H0_large.csv")

z_full <- matrix(0, nrow = 21, ncol = 50)
w_full <- matrix(0, nrow = 21, ncol = 50)
for (k in 1:1000){
  ## the process below repeats the process detailed for sample size calculation
  ## 1 above for all repetitions k = {1, 2, ..., 1000}
  H1_low <- unlist(H1_lows[k,]); H1_mid <- unlist(H1_mids[k,]); 
  H1_high <- unlist(H1_highs[k,])
  
  H0_low <- unlist(H0_lows[k,]); H0_mid <- unlist(H0_mids[k,]); 
  H0_high <- unlist(H0_highs[k,])
  
  n_low <- n_lows[k]; n_mid <- n_mids[k]; n_high <- n_highs[k]
  
  H1_slope <- (H1_mid - H1_low)/(n_mid-n_low)
  H0_slope <- (H0_mid - H0_low)/(n_mid-n_low)
  H1_int <- H1_low - H1_slope*n_low
  H0_int <- H0_low - H0_slope*n_low
  
  for (i in seq(25, as.numeric(n_mid))){
    assign(paste0("H1_vec_", i), H1_int + H1_slope*i)
    assign(paste0("H0_vec_", i), H0_int + H0_slope*i)
  }
  
  H1_slope <- (H1_high - H1_mid)/(n_high-n_mid)
  H0_slope <- (H0_high - H0_mid)/(n_high-n_mid)
  H1_int <- H1_mid - H1_slope*n_mid
  H0_int <- H0_mid - H0_slope*n_mid
  
  for (i in seq(as.numeric(n_mid) + 1, 45)){
    assign(paste0("H1_vec_", i), H1_int + H1_slope*i)
    assign(paste0("H0_vec_", i), H0_int + H0_slope*i)
  }
  
  gammas <- seq(log(0.93) - log(0.07), log(0.97) - log(0.03), length.out = 50)
  
  x <- seq(25, 45, 1)
  y <- 1/(1 + exp(-gammas))
  
  z_mat <- NULL
  for (i in 1:length(x)){
    z_mat <- rbind(z_mat, get(paste0("H1_vec_", x[i])))
  }
  z <- NULL
  for (j in 1:length(y)){
    z <- cbind(z, rowMeans(z_mat > gammas[j]))
  }
  
  ## now we multiply the matrix for each sample size calculation by 1/1000
  ## to get the matrices corresponding to the average
  z_full <- 0.001*z + z_full
  
  w_mat <- NULL
  for (i in 1:length(x)){
    w_mat <- rbind(w_mat, get(paste0("H0_vec_", x[i])))
  }
  w <- NULL
  for (j in 1:length(y)){
    w <- cbind(w, rowMeans(w_mat > gammas[j]))
  }
  
  w_full <- 0.001*w + w_full
  
}
write.csv(z_full, "z_full_lin.csv", row.names = FALSE)
write.csv(w_full, "w_full_lin.csv", row.names = FALSE)

## create matrices for right contour plot in Figure D.1 (based on
## linear regression data simulated above)
z_full2 <- matrix(0, nrow = 21, ncol = 50)
w_full2 <- matrix(0, nrow = 21, ncol = 50)
## convert the posterior probabilities for each approximated sampling
## distribution to the logit scale (with error checking to ensure
## to logits are finite)
for (i in seq(25,45,1)){
  assign(paste0("H1_vec_", i), 
         pmin(pmax(as.numeric(unlist(read.csv(paste0("probs_H1_lin_", i,".csv")))),
                   .Machine$double.eps), 1 - 10^(-7)))
  assign(paste0("H0_vec_", i),
         pmin(pmax(as.numeric(unlist(read.csv(paste0("probs_H0_lin_", i,".csv")))),
                   .Machine$double.eps), 1 - 10^(-7)))
}

## this process mirrors what was done to create the z and w matrices in 
## the previous two plots but with the estimates obtained by simulating data
gammas <- seq(log(0.93) - log(0.07), log(0.97) - log(0.03), length.out = 50)

x <- seq(25, 45, 1)
y <- 1/(1 + exp(-gammas))

z_mat <- NULL
for (i in 1:length(x)){
  z_mat <- rbind(z_mat, get(paste0("H1_vec_", x[i])))
}
z_mat <- log(z_mat) - log(1-z_mat)
z <- NULL
for (j in 1:length(y)){
  z <- cbind(z, rowMeans(z_mat > gammas[j]))
}

z_full2 <- z

w_mat <- NULL
for (i in 1:length(x)){
  w_mat <- rbind(w_mat, get(paste0("H0_vec_", x[i])))
}
w_mat <- log(w_mat) - log(1-w_mat)
w <- NULL
for (j in 1:length(y)){
  w <- cbind(w, rowMeans(w_mat > gammas[j]))
}

w_full2 <- w

## write output to a .csv file  
write.csv(z_full2, "z_full_lin2.csv", row.names = FALSE)
write.csv(w_full2, "w_full_lin2.csv", row.names = FALSE)

## create the three contour plots and output as .pdf file for the article
pdf(file = "FigureLinContour.pdf",   # The directory you want to save the file in
    width = 6.5, 
    height = 6.5) 

par(mfrow=c(2,2), mar = c(3.75, 3.75, 2, 0.35) + 0.1, mgp=c(2.25,1,0))

x <- seq(25,45, 1)
gammas <- seq(log(0.93) - log(0.07), log(0.97) - log(0.03), length.out = 50)
y <- 1/(1 + exp(-gammas))

z_full <- matrix(unlist(read.csv("z_full_lin.csv")), nrow = 21, ncol = 50)
w_full <- matrix(unlist(read.csv("w_full_lin.csv")), nrow = 21, ncol = 50)

contour(x, y, w_full, levels = c(seq(0.02, 0.035, 0.005), 0.045, seq(0.065, 0.08, 0.005)), 
        xlab = expression(italic("n")), xlim = c(25,45),
        ylab = expression(gamma),  main = "Type I Error Rate", labcex = 0.8, method = "edge", axes = FALSE, cex.lab = 1.25)
contour(x, y, w_full, levels = c(0.05), col = "firebrick", add = TRUE, labcex = 0.8, method = "edge")
contour(x, y, z_full, levels = c(0.8), col = "seagreen", add = TRUE, labcex = 0.8, labels = "", method = "edge")
contour(x, y, w_full, levels = c(0.055), col = "black", add = TRUE, labcex = 0.8, method = "edge")
contour(x, y, w_full, levels = c(0.06), col = "black", add = TRUE, labcex = 0.8, method = "edge")
contour(x, y, w_full, levels = c(0.04), col = "black", add = TRUE, labcex = 0.8, method = "edge")
axis(side = 1, at = seq(25, 45, 5), cex.axis = 1.15)
axis(side = 2, at = seq(0.93, 0.97, 0.01), cex.axis = 1.15)
box()

z_full2 <- matrix(unlist(read.csv("z_full_lin2.csv")), nrow = 21, ncol = 50)
w_full2 <- matrix(unlist(read.csv("w_full_lin2.csv")), nrow = 21, ncol = 50)

contour(x, y, w_full2, levels = c(seq(0.02, 0.035, 0.005), 0.045, seq(0.065, 0.08, 0.005)), 
        xlab = expression(italic("n")), xlim = c(25,45),
        ylab = expression(gamma),  main = "Type I Error Rate", labcex = 0.8, method = "edge", axes = FALSE, cex.lab = 1.25)
contour(x, y, w_full2, levels = c(0.05), col = "firebrick", add = TRUE, labcex = 0.8, method = "edge")
contour(x, y, z_full2, levels = c(0.8), col = "seagreen", add = TRUE, labcex = 0.8, labels = "", method = "edge")
contour(x, y, w_full2, levels = c(0.055), col = "black", add = TRUE, labcex = 0.8, method = "edge")
contour(x, y, w_full2, levels = c(0.06), col = "black", add = TRUE, labcex = 0.8, method = "edge")
contour(x, y, w_full2, levels = c(0.04), col = "black", add = TRUE, labcex = 0.8)
axis(side = 1, at = seq(25, 45, 5), cex.axis = 1.15)
axis(side = 2, at = seq(0.93, 0.97, 0.01), cex.axis = 1.15)
box()

contour(x, y, z_full, levels = c(seq(0.68, 0.77, 0.015), 0.80, seq(0.83, 0.89, 0.015)), 
        xlab = expression(italic("n")), xlim = c(25,45),axes = FALSE, cex.lab = 1.25,
        ylab = expression(gamma),  main = "Power", labcex = 0.8, method = "edge")
contour(x, y, w_full, levels = c(0.05), col = "firebrick", add = TRUE, labcex = 0.8,labels = "", method = "edge") 
contour(x, y, z_full, levels = c(0.8), col = "seagreen", add = TRUE, labcex = 0.8, method = "edge")
contour(x, y, z_full, levels = c(0.815), col = "black", add = TRUE, method = "edge", labcex = 0.8)
contour(x, y, z_full, levels = c(0.785), col = "black", add = TRUE, method = "edge", labcex = 0.8)
axis(side = 1, at = seq(25,45, 5), cex.axis = 1.15)
axis(side = 2, at = seq(0.93, 0.97, 0.01), cex.axis = 1.15)
box()

contour(x, y, z_full2, levels = c(seq(0.68, 0.77, 0.015), 0.80, seq(0.83, 0.89, 0.015)), 
        xlab = expression(italic("n")), xlim = c(25,45),axes = FALSE, cex.lab = 1.25,
        ylab = expression(gamma),  main = "Power", labcex = 0.8, method = "edge")
contour(x, y, w_full2, levels = c(0.05), col = "firebrick", add = TRUE, labcex = 0.8,labels = "", method = "edge") 
contour(x, y, z_full2, levels = c(0.8), col = "seagreen", add = TRUE, labcex = 0.8, method = "edge")
contour(x, y, z_full2, levels = c(0.815), col = "black", add = TRUE, method = "edge", labcex = 0.8)
contour(x, y, z_full2, levels = c(0.785), col = "black", add = TRUE, method = "edge", labcex = 0.8)
axis(side = 1, at = seq(25,45, 5), cex.axis = 1.15)
axis(side = 2, at = seq(0.93, 0.97, 0.01), cex.axis = 1.15)
box()

par(mfrow=c(1,1))
dev.off()

## implement the sample size calculation 1000 times without covariates using different
## Sobol' sequences

## set hyperparameters with proper dimension
hyper_ex <- list(rep(0,2), 0.01*diag(2), 1, 1)

tempNoCov <- foreach(i=1:1000, .combine='rbind', .packages = c("qrng"),
                     .options.snow=opts, .errorhandling = "remove") %dopar% {
                       temp_res <- findGammaNoCov(eta_plus1 = eta_plus1s, eta_plus0 = eta_plus0s, 
                                                  pwr = 0.8, typeI = 0.05, deltas = c(5, Inf), q = 2, 
                                                  hypers = hyper_ex, m = 4096, m0 = 64,
                                                  seed_H1 = i, seed_H0 = i + 1000, contour = FALSE)
                       
                       c(i, as.numeric(unlist(temp_res)))
                     }

write.csv(tempNoCov, "lin_4096_nocov_summary.csv", row.names = FALSE)

## get 95% bootstrap CI for sample size n_B
quantile(tempNoCov[,2], c(0.025, 0.975))
## repeat for gamma
round(quantile(tempNoCov[,3], c(0.025, 0.975)),4)