## BEGIN SETUP ##

## code for numerical studies with gamma example
## load necessary packages
require(qrng)
require(foreach)
require(doParallel)
require(doSNOW)
require(ggplot2)
require(ggpubr)
require(cowplot)
require(colorspace)
require(rjags)
require(coda)

## set up parallelization with 1000 repetitions on 72 cores
cores=detectCores()
cl <- makeSOCKcluster(cores[1]-1)

registerDoSNOW(cl)
pb <- txtProgressBar(max = 1000, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

## we use the same parameter draws from Psi_1 and Psi_0
## for all m simulation repetitons for this example,
## but we could change this
m <- 8192
eta_plus1s <- c(10.04, 12.37, 14.46, 14.90)
eta_plus0s <- c(10.04, 11.495, 14.46, 14.90)
eta_plus1s <- matrix(eta_plus1s, ncol = 4, nrow =m, byrow = TRUE)
eta_plus0s <- matrix(eta_plus0s, ncol = 4, nrow =m, byrow = TRUE)
informs <- read.csv("informs.csv")

## define hyperparameters for analysis priors
alphas1 <- c(informs[1,1], informs[1,2])
betas1 <- c(informs[2,1], informs[2,2])
alphas2 <- c(informs[3,1], informs[3,2])
betas2 <- c(informs[4,1], informs[4,2])

## repeat the sample size calculation for the gamma example 1000 times
## with the different Sobol' sequences (and informative priors)
temp <- foreach(i=1:1000, .combine='rbind', .packages = c("qrng", "nleqslv"),
                .options.snow=opts, .errorhandling = "remove") %dopar% {
                  temp_res <-  findGamma(eta_plus1 = eta_plus1s, eta_plus0 = eta_plus0s, 
                                         pwr = 0.9, typeI = 0.01, 
                                         deltas = c(-Inf, log(0.9)), q = 2, 
                                         alphas1, betas1, alphas2, betas2, m = 8192,
                                         seed_H1 = i, seed_H0 = i + 1000, contour = FALSE)
                  
                  c(i, as.numeric(unlist(temp_res)))
                }

## output the optimal design and summary for each simulation repetition
write.csv(temp, "gamma_sobol_8192_inf.csv", row.names = FALSE)

## get 95% bootstrap CI for sample size n_B
quantile(temp[,2], c(0.025, 0.975))
## repeat for gamma
round(quantile(temp[,3], c(0.025, 0.975)),4)

## repeat the sample size calculation for the gamma example 1000 times
## with the different Sobol' sequences (and informative priors)
tempFull <- foreach(i=1:1000, .combine='rbind', .packages = c("qrng", "nleqslv"),
                .options.snow=opts, .errorhandling = "remove") %dopar% {
                  temp_res <-  findGammaFull(eta_plus1 = eta_plus1s, eta_plus0 = eta_plus0s, 
                                             pwr = 0.9, typeI = 0.01, 
                                             deltas = c(-Inf, log(0.9)), q = 2, 
                                             alphas1, betas1, alphas2, betas2, m = 8192,
                                             seed_H1 = i, seed_H0 = i + 1000)
                  
                  c(i, as.numeric(unlist(temp_res)))
                }

## output the optimal design and summary for each simulation repetition
write.csv(tempFull, "gamma_sobol_8192_inf_full.csv", row.names = FALSE)

## check if (n, gamma) combinations are the same
sum(temp[,2] == tempFull[,2])
sum(temp[,3] == tempFull[,3])

## comparison precision of Sobol' sequences with precision of 
## longer pseudorandom ones
m <- 75000
eta_plus1s <- c(10.04, 12.37, 14.46, 14.90)
eta_plus0s <- c(10.04, 11.495, 14.46, 14.90)
eta_plus1s <- matrix(eta_plus1s, ncol = 4, nrow =m, byrow = TRUE)
eta_plus0s <- matrix(eta_plus0s, ncol = 4, nrow =m, byrow = TRUE)

## repeat the sample size calculation for the gamma example 1000 times
## with the different pseudorandom sequences
temp <- foreach(i=1:1000, .combine='rbind', .packages = c("qrng", "nleqslv"),
                .options.snow=opts, .errorhandling = "remove") %dopar% {
                  temp_res <-  findGamma(eta_plus1 = eta_plus1s, eta_plus0 = eta_plus0s, 
                                         pwr = 0.9, typeI = 0.01, prng = TRUE,
                                         deltas = c(-Inf, log(0.9)), q = 2, 
                                         alphas1, betas1, alphas2, betas2, m = 75000,
                                         seed_H1 = i + 8000, seed_H0 = i + 9000, contour = FALSE)
                  
                  c(i, as.numeric(unlist(temp_res)))
                }

## output the optimal design and summary for each simulation repetition
write.csv(temp, "gamma_prng_75k.csv", row.names = FALSE)

## get 95% bootstrap CI for sample size n_B
quantile(temp[,2], c(0.025, 0.975))
## repeat for gamma
round(quantile(temp[,3], c(0.025, 0.975)),4)

## confirmatory estimates of operating characteristics
## analytical sample size calculation is 288, so consider (288, 0.99) combination
q <- 2
n <- 288
lA <- log(0.99) - log(0.01)
m <- 8192
deltas <- c(-Inf, log(0.9))
eta_plus1 <- eta_plus1s
eta_plus0 <- eta_plus0s
temp288p99 <- foreach(i=1:1000, .combine='rbind', .packages = c("qrng"),
                      .options.snow=opts, .errorhandling = "remove") %dopar% {
                        
                        sob_H1 <- sobol(m, d = 4, randomize = "digital.shift", seed = i)
                        sob_H0 <- sobol(m, d = 4, randomize = "digital.shift", seed = i + 1000)
                        
                        H1_vec <- NULL
                        H0_vec <- NULL
                        ## implement Algorithm 2 for each point
                        for (i in 1:m){
                          H1_vec[i] <- targetPower(n_val = n, q = q,
                                                      params = as.numeric(eta_plus1[i,]),
                                                      delta = deltas, u = sob_H1[i,1:4],
                                                      hyper = rbind(alphas1, betas1, 
                                                                    alphas2, betas2))
                          H0_vec[i] <- targetPower(n_val = n, q = q,
                                                    params = as.numeric(eta_plus0[i,]),
                                                    delta = deltas, u = sob_H0[i,1:4],
                                                    hyper = rbind(alphas1, betas1, 
                                                                  alphas2, betas2))
                        }
                        c(i, mean(H0_vec >= lA), mean(H1_vec >= lA))
                      }

write.csv(temp288p99, "gamma_inf_8192_288p99_summary.csv", row.names = FALSE)

## calculate type I error rate for this design
mean(temp288p99[,2])
## calculate power for this design
mean(temp288p99[,3])

## define hyperparameters for uninformative analysis priors (for comparison)
alphas1 <- c(2, 0.25)
betas1 <- c(2, 0.25)
alphas2 <- c(2, 0.25)
betas2 <- c(2, 0.25)

temp288p99un <- foreach(i=1:1000, .combine='rbind', .packages = c("qrng"),
                        .options.snow=opts, .errorhandling = "remove") %dopar% {
                          
                          sob_H1 <- sobol(m, d = 4, randomize = "digital.shift", seed = i)
                          sob_H0 <- sobol(m, d = 4, randomize = "digital.shift", seed = i + 1000)
                          
                          H1_vec <- NULL
                          H0_vec <- NULL
                          ## implement Algorithm 3 for each of the first m0 points
                          for (i in 1:m){
                            H1_vec[i] <- targetPower(n_val = n, q = q,
                                                        params = as.numeric(eta_plus1[i,]),
                                                        delta = deltas, u = sob_H1[i,1:4],
                                                        hyper = rbind(alphas1, betas1, 
                                                                      alphas2, betas2))
                            H0_vec[i] <- targetPower(n_val = n, q = q,
                                                      params = as.numeric(eta_plus0[i,]),
                                                      delta = deltas, u = sob_H0[i,1:4],
                                                      hyper = rbind(alphas1, betas1, 
                                                                    alphas2, betas2))
                          }
                          c(i, mean(H0_vec >= lA), mean(H1_vec >= lA))
                        }

write.csv(temp288p99un, "gamma_uninf_8192_288p99_summary.csv", row.names = FALSE)

## calculate type I error rate for this design
mean(temp288p99un[,2])
## calculate power for this design
mean(temp288p99un[,3])

## switch back to informative priors to construct contour plots in Section 6
informs <- read.csv("informs.csv")
alphas1 <- c(informs[1,1], informs[1,2])
betas1 <- c(informs[2,1], informs[2,2])
alphas2 <- c(informs[3,1], informs[3,2])
betas2 <- c(informs[4,1], informs[4,2])

## repeat the sample size calculation for the gamma example 1000 times
## with the different Sobol' sequences (with contour = TRUE)
temp <- foreach(i=1:1000, .combine='rbind', .packages = c("qrng", "nleqslv"),
                .options.snow=opts, .errorhandling = "remove") %dopar% {
                  temp_res <-  findGamma(eta_plus1 = eta_plus1s, eta_plus0 = eta_plus0s, 
                                             pwr = 0.9, typeI = 0.01, 
                                             deltas = c(-Inf, log(0.9)), q = 2, 
                                             alphas1, betas1, alphas2, betas2, m = 8192,
                                             seed_H1 = i, seed_H0 = i + 1000, contour = TRUE)
                  
                  c(i, as.numeric(unlist(temp_res)))
                }

## output the optimal design and summary for each simulation repetition
write.csv(temp[,1:9], "gam_sobol_8192_summary.csv", row.names = FALSE)

## output the red and green posterior probabilities for each of the three sample
## sizes n^((0)) = small, n^((1)) = mid, and n^((2)) = large
write.csv(temp[,10:8201], "gam_sobol_8192_H0_small.csv", row.names = FALSE)
write.csv(temp[,8202:16393], "gam_sobol_8192_H1_small.csv", row.names = FALSE)
write.csv(temp[,16394:24585], "gam_sobol_8192_H0_mid.csv", row.names = FALSE)
write.csv(temp[,24586:32777], "gam_sobol_8192_H1_mid.csv", row.names = FALSE)
write.csv(temp[,32778:40969], "gam_sobol_8192_H0_large.csv", row.names = FALSE)
write.csv(temp[,40970:49161], "gam_sobol_8192_H1_large.csv", row.names = FALSE)

## ensure GammaPostMean() was loaded in "01-sec5p1-setup.R" since it is used here

## now use 81920 simulation repetitions
resamps <- 81920
pb <- txtProgressBar(max = resamps, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

## get parameter values from design priors to generate
## the gamma data under H1
samps <- seq(275, 325, 2)

delta_L <- 0
delta_U <- 0.9

## estimate sampling distributions of posterior probabilities under H1
for (k in 1:length(samps)){
  probs_in <- foreach(j=1:resamps, .combine='rbind', .packages = c("rjags", "coda"),
                      .options.snow=opts, .errorhandling = "remove") %dopar% {
                        
                        set.seed(k*resamps + j)
                        
                        y_star1 <- rgamma(2*samps[k], eta_plus1s[1], eta_plus1s[2])
                        y_star2 <- rgamma(samps[k], eta_plus1s[3], eta_plus1s[4])
                        
                        theta.diff <- GammaPostMean(y_star1, y_star2, informs[1,1], informs[1,2],
                                                    informs[2,1], informs[2,2], informs[3,1], 
                                                    informs[3,2], informs[4,1], informs[4,2],
                                                    threshold)
                        
                        mean(ifelse(theta.diff > delta_L, theta.diff < delta_U,0))
                      }
  
  ## output results to .csv file
  write.csv(probs_in, paste0("probs_H1_gamma_", samps[k], ".csv"), row.names = FALSE)
}

## get parameter values from design priors to generate
## the gamma data under H0
params <- rbind(eta_plus0s, eta_plus0s, eta_plus0s, eta_plus0s, eta_plus0s)

## estimate sampling distributions of posterior probabilities under H0
for (k in 1:length(samps)){
  
  probs_out <- foreach(j=1:resamps, .combine='rbind', .packages = c("rjags", "coda"),
                       .options.snow=opts, .errorhandling = "remove") %dopar% {
                         
                         set.seed(k*resamps + j)
                         
                         y_star1 <- rgamma(2*samps[k], eta_plus0s[1], eta_plus0s[2])
                         y_star2 <- rgamma(samps[k], eta_plus0s[3], eta_plus0s[4])
                         
                         theta.diff <- GammaPostMean(y_star1, y_star2, informs[1,1], informs[1,2],
                                                     informs[2,1], informs[2,2], informs[3,1], 
                                                     informs[3,2], informs[4,1], informs[4,2],
                                                     threshold)
                         
                         mean(ifelse(theta.diff > delta_L, theta.diff < delta_U,0))
                       }
  
  ## output results to .csv file
  write.csv(probs_out, paste0("probs_H0_gamma_", samps[k], ".csv"), row.names = FALSE)
}

## code to create contour plots
## create matrices for left contour plot in Figure 1 (based on the first
## sample size calculation)
first_rep <- as.numeric(read.csv("gam_sobol_8192_summary.csv")[1,])
n_low <- first_rep[7]; n_mid <- first_rep[8];
n_high <- first_rep[9]

## read in the posterior probabilities corresponding to n0 and n1
H1_mid <- unlist(read.csv("gam_sobol_8192_H1_mid.csv")[1,])
H0_mid <- unlist(read.csv("gam_sobol_8192_H0_mid.csv")[1,])
H1_low <- unlist(read.csv("gam_sobol_8192_H1_small.csv")[1,])
H0_low <- unlist(read.csv("gam_sobol_8192_H0_small.csv")[1,])

## get the slopes and intercepts for the linear approximations
H1_slope <- (H1_mid - H1_low)/(n_mid-n_low)
H0_slope <- (H0_mid - H0_low)/(n_mid-n_low)
H1_int <- H1_low - H1_slope*n_low
H0_int <- H0_low - H0_slope*n_low

## approximate the sampling distributions of posterior probabilities
## on the logit scale using these approximations
for (i in seq(275, as.numeric(n_mid))){
  assign(paste0("H1_vec_", i), H1_int + H1_slope*i)
  assign(paste0("H0_vec_", i), H0_int + H0_slope*i)
}

## read in the posterior probabilities corresponding to n^((2))
H1_high <- unlist(read.csv("gam_sobol_8192_H1_large.csv")[1,])
H0_high <- unlist(read.csv("gam_sobol_8192_H0_large.csv")[1,])

## get the slopes and intercepts for the linear approximations
## use n2 instead of n0 this time
H1_slope <- (H1_high - H1_mid)/(n_high-n_mid)
H0_slope <- (H0_high - H0_mid)/(n_high-n_mid)
H1_int <- H1_mid - H1_slope*n_mid
H0_int <- H0_mid - H0_slope*n_mid

## approximate the sampling distributions of posterior probabilities
## on the logit scale using these approximations
for (i in seq(as.numeric(n_mid) + 1, 325)){
  assign(paste0("H1_vec_", i), H1_int + H1_slope*i)
  assign(paste0("H0_vec_", i), H0_int + H0_slope*i)
}

## create a vector of gamma values on the logit scale to compute power
## and type I error rate estimates
opt_gamma <- as.numeric(first_rep[3])
opt_gamma <- log(opt_gamma) - log(1 - opt_gamma)

gammas <- seq(log(0.988) - log(1 - 0.988), log(0.997) - log(1 - 0.997), length.out = 50)
gammas <- sort(c(opt_gamma, gammas))

x <- seq(275, 325, 1)
y <- 1/(1 + exp(-gammas))

## z matrix is for power
z_mat <- NULL
for (i in 1:length(x)){
  z_mat <- rbind(z_mat, get(paste0("H1_vec_", x[i])))
}
z <- NULL
for (j in 1:length(y)){
  z <- cbind(z, rowMeans(z_mat > gammas[j]))
}

## w matrix is for type I error
w_mat <- NULL
for (i in 1:length(x)){
  w_mat <- rbind(w_mat, get(paste0("H0_vec_", x[i])))
}
w <- NULL
for (j in 1:length(y)){
  w <- cbind(w, rowMeans(w_mat > gammas[j]))
}

write.csv(w, "w_mat1.csv", row.names = FALSE)
write.csv(z, "z_mat1.csv", row.names = FALSE)

## create matrices for center contour plot in Figure 1 (based on the average
## of 1000 sample size calculations)

## read in the sample sizes and posterior probabilities (logit scale) for 
## all 1000 repetitions
n_lows <- as.numeric(read.csv("gam_sobol_8192_summary.csv")[,7])
n_mids <- as.numeric(read.csv("gam_sobol_8192_summary.csv")[,8])
n_highs <- as.numeric(read.csv("gam_sobol_8192_summary.csv")[,9])

H1_mids <- read.csv("gam_sobol_8192_H1_mid.csv")
H0_mids <- read.csv("gam_sobol_8192_H0_mid.csv")
H1_lows <- read.csv("gam_sobol_8192_H1_small.csv")
H0_lows <- read.csv("gam_sobol_8192_H0_small.csv")
H1_highs <- read.csv("gam_sobol_8192_H1_large.csv")
H0_highs <- read.csv("gam_sobol_8192_H0_large.csv")

z_full <- matrix(0, nrow = 51, ncol = 50)
w_full <- matrix(0, nrow = 51, ncol = 50)
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
  
  for (i in seq(275, as.numeric(n_mid))){
    assign(paste0("H1_vec_", i), H1_int + H1_slope*i)
    assign(paste0("H0_vec_", i), H0_int + H0_slope*i)
  }
  
  H1_slope <- (H1_high - H1_mid)/(n_high-n_mid)
  H0_slope <- (H0_high - H0_mid)/(n_high-n_mid)
  H1_int <- H1_mid - H1_slope*n_mid
  H0_int <- H0_mid - H0_slope*n_mid
  
  for (i in seq(as.numeric(n_mid) + 1, 325)){
    assign(paste0("H1_vec_", i), H1_int + H1_slope*i)
    assign(paste0("H0_vec_", i), H0_int + H0_slope*i)
  }
  
  gammas <- seq(log(0.988) - log(1 - 0.988), log(0.997) - log(1 - 0.997), length.out = 50)
  
  x <- seq(275, 325, 1)
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
write.csv(z_full, "z_full_mat.csv", row.names = FALSE)
write.csv(w_full, "w_full_mat.csv", row.names = FALSE)

## create matrices for right contour plot in Figure 1 (based on
## simulating gamma data)
z_full2 <- matrix(0, nrow = 26, ncol = 50)
w_full2 <- matrix(0, nrow = 26, ncol = 50)
## convert the posterior probabilities for each approximated sampling
## distribution to the logit scale (with error checking to ensure
## to logits are finite)
for (i in seq(275, 325,2)){
  assign(paste0("H1_vec_", i), 
         pmin(pmax(as.numeric(unlist(read.csv(paste0("probs_H1_gamma_", i,".csv")))),
                   .Machine$double.eps), 1 - 10^(-7)))
  assign(paste0("H0_vec_", i), 
         pmin(pmax(as.numeric(unlist(read.csv(paste0("probs_H0_gamma_", i,".csv")))),
                   .Machine$double.eps), 1 - 10^(-7)))
}

## this process mirrors what was done to create the z and w matrices in 
## the previous two plots but with the estimates obtained by simulating data
gammas <- seq(log(0.988) - log(1 - 0.988), log(0.997) - log(1 - 0.997), length.out = 50)

x <- seq(275, 325, 2)
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
write.csv(z_full2, "z_full_mat2.csv", row.names = FALSE)
write.csv(w_full2, "w_full_mat2.csv", row.names = FALSE)

## create the three contour plots and output as .pdf file for the article
pdf(file = "Figure1.pdf",   # The directory you want to save the file in
    width = 7.5, 
    height = 5) 

par(mfrow=c(2,3), mar = c(3.75, 3.75, 2, 0.35) + 0.1, mgp=c(2.25,1,0))

## read in matrices for left plot
z <- matrix(unlist(read.csv("z_mat1.csv")), nrow = 51, ncol = 51)
w <- matrix(unlist(read.csv("w_mat1.csv")), nrow =51, ncol = 51)
gammas <- seq(log(0.988) - log(1 - 0.988), log(0.997) - log(1 - 0.997), length.out = 50)
gammas <- sort(c(opt_gamma, gammas))
y <- 1/(1 + exp(-gammas))
x <- seq(275, 325, 1)

contour(x, y, w, levels = c(seq(0.0025, 0.0100, 0.0025), seq(0.0125, 0.02, 0.0025)), 
        xlab = expression(italic("n")[italic("B")]), xlim = c(275,325),
        ylab = expression(gamma),  main = "Type I Error Rate", labcex = 0.8, method = "edge", axes = FALSE, cex.lab = 1.25)
points(x = first_rep[2], y = 1/(1 + exp(-opt_gamma)), pch = 19, col = adjustcolor("grey50", 0.75))
contour(x, y, w, levels = c(0.01), col = "firebrick", add = TRUE, labcex = 0.8, method = "edge")
contour(x, y, z, levels = c(0.9), col = "seagreen", add = TRUE, labcex = 0.8, labels = "", method = "edge")
axis(side = 1, at = seq(275, 325, 10), cex.axis = 1.15)
axis(side = 2, at = seq(0.988, 0.997, 0.003), cex.axis = 1.15)
box()

gammas <- seq(log(0.988) - log(1 - 0.988), log(0.997) - log(1 - 0.997), length.out = 50)
y <- 1/(1 + exp(-gammas))

z_full <- matrix(unlist(read.csv("z_full_mat.csv")), nrow = 51, ncol = 50)
w_full <- matrix(unlist(read.csv("w_full_mat.csv")), nrow = 51, ncol = 50)

contour(x, y, w_full, levels = c(seq(0.0025, 0.0100, 0.0025), seq(0.0125, 0.02, 0.0025)), 
        xlab = expression(italic("n")[italic("B")]), 
        xlim = c(275,325), axes = FALSE, cex.lab = 1.25,
        ylab = expression(gamma),  main = "Type I Error Rate", labcex = 0.8, method = "edge")
contour(x, y, w_full, levels = c(0.01), col = "firebrick", add = TRUE, labcex = 0.8, method = "edge")
contour(x, y, z_full, levels = c(0.9), col = "seagreen", add = TRUE, labcex = 0.8, labels = "", method = "edge")
axis(side = 1, at = seq(275, 325, 10), cex.axis = 1.15)
axis(side = 2, at = seq(0.988, 0.997, 0.003), cex.axis = 1.15)
box()

z_full2 <- matrix(unlist(read.csv("z_full_mat2_long.csv")), nrow = 26, ncol = 50)
w_full2 <- matrix(unlist(read.csv("w_full_mat2_long.csv")), nrow = 26, ncol = 50)
x <- seq(275, 325, 2)

contour(x, y, w_full2, levels = c(seq(0.0025, 0.0050, 0.0025), seq(0.015, 0.02, 0.0025)), 
        xlab = expression(italic("n")[italic("B")]), xlim = c(275,325), axes = FALSE, cex.lab = 1.25,
        ylab = expression(gamma),  main = "Type I Error Rate", labcex = 0.8, method = "edge")
contour(x, y, w_full2, levels = c(0.01), col = "firebrick", add = TRUE, labcex = 0.8, method = "edge")
contour(x, y, z_full2, levels = c(0.9), col = "seagreen", add = TRUE, labcex = 0.8, labels = "", method = "edge")
contour(x, y, w_full2, levels = c(0.0075), col = "black", add = TRUE, method = "flattest", labcex = 0.8)
contour(x, y, w_full2, levels = c(0.0125), col = "black", add = TRUE, method = "flattest", labcex = 0.8)
axis(side = 1, at = seq(275, 325, 10), cex.axis = 1.15)
axis(side = 2, at = seq(0.988, 0.997, 0.003), cex.axis = 1.15)
box()

gammas <- seq(log(0.988) - log(1 - 0.988), log(0.997) - log(1 - 0.997), length.out = 50)
gammas <- sort(c(opt_gamma, gammas))
y <- 1/(1 + exp(-gammas))
x <- seq(275, 325, 1)

contour(x, y, z, levels = c(seq(0.795, 0.885, 0.015), seq(0.915, 0.96, 0.015)), 
        xlab = expression(italic("n")[italic("B")]), xlim = c(275,325),axes = FALSE, cex.lab = 1.25,
        ylab = expression(gamma),  main = "Power", labcex = 0.8, method = "edge")
points(x = first_rep[2], y = 1/(1 + exp(-opt_gamma)), pch = 19, col = adjustcolor("grey50", 0.75))
contour(x, y, w, levels = c(0.01), col = "firebrick", add = TRUE, labcex = 0.8,labels = "", method = "edge") 
contour(x, y, z, levels = c(0.9), col = "seagreen", add = TRUE, labcex = 0.8, method = "edge")
axis(side = 1, at = seq(275, 325, 10), cex.axis = 1.15)
axis(side = 2, at = seq(0.988, 0.997, 0.003), cex.axis = 1.15)
box()

gammas <- seq(log(0.988) - log(1 - 0.988), log(0.997) - log(1 - 0.997), length.out = 50)
y <- 1/(1 + exp(-gammas))

contour(x, y, z_full, levels = c(seq(0.795, 0.885, 0.015), seq(0.915, 0.96, 0.015)), 
        xlab = expression(italic("n")[italic("B")]), xlim = c(275,325), axes = FALSE,
        ylab = expression(gamma),  main = "Power", labcex = 0.8, method = "edge", cex.lab = 1.25)
axis(side = 1, at = seq(100, 120, 5), cex.axis = 1.15)
axis(side = 2, at = seq(0.91, 0.955, 0.015), cex.axis = 1.15)
box()
contour(x, y, w_full, levels = c(0.01), col = "firebrick", add = TRUE, labcex = 0.8,labels = "", method = "edge") 
contour(x, y, z_full, levels = c(0.9), col = "seagreen", add = TRUE, labcex = 0.8, method = "edge")
axis(side = 1, at = seq(275, 325, 10), cex.axis = 1.15)
axis(side = 2, at = seq(0.988, 0.997, 0.003), cex.axis = 1.15)
box()
x <- seq(275, 325, 2)

contour(x, y, z_full2, levels = c(seq(0.795, 0.885, 0.015), seq(0.915, 0.96, 0.015)), 
        xlab = expression(italic("n")[italic("B")]), xlim = c(275,325), axes = FALSE, cex.lab = 1.25,
        ylab = expression(gamma),  main = "Power", labcex = 0.8, method = "edge")
contour(x, y, w_full2, levels = c(0.01), col = "firebrick", add = TRUE, labcex = 0.8,labels = "", method = "edge") 
contour(x, y, z_full2, levels = c(0.9), col = "seagreen", add = TRUE, labcex = 0.8, method = "edge")
axis(side = 1, at = seq(275, 325, 10), cex.axis = 1.15)
axis(side = 2, at = seq(0.988, 0.997, 0.003), cex.axis = 1.15)
box()

par(mfrow=c(1,1))
dev.off()