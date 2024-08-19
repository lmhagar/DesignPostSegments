## BEGIN SETUP ##

## code for numerical studies with Weibull example
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

## save these parameter values drawn from degenerate Psi0 and Psi0
eta_plus1s <- c(0.90, 3.09, 1.07, 3.76)
eta_plus0s <- c(0.973, 3.09, 1.07, 3.76)

## define hyperparameters for analysis priors
informs <- read.csv("informs_weibull.csv") ## generated in file 06
alphas1 <- informs[2,]
betas1 <- informs[1,]
alphas2 <- informs[4,]
betas2 <- informs[3,]

## set up parallelization with 1000 repetitions on 72 cores
cores=detectCores()
cl <- makeSOCKcluster(cores[1]-1)

registerDoSNOW(cl)
pb <- txtProgressBar(max = 1000, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

## repeat the sample size calculation for the Weibull example 1000 times
## with the different Sobol' sequences
temp <- foreach(i=1:1000, .combine='rbind', .packages = c("qrng", "nleqslv"),
                .options.snow=opts, .errorhandling = "remove") %dopar% {
                  temp_res <-  findGamma(eta_plus1 = eta_plus1s, eta_plus0 = eta_plus0s, 
                                         pwr = 0.9, typeI = 0.01, deltas = c(-Inf, log(0.9)), q = 2, 
                                         alphas1, betas1, alphas2, betas2, m = 8192,
                                         seed_H1 = i + 8000, seed_H0 = i + 9000, contour = TRUE)
                  
                  c(i, as.numeric(unlist(temp_res)))
                }

## output the optimal design and summary for each simulation repetition
write.csv(temp[,1:9], "wei_sobol_8192_summary.csv", row.names = FALSE)

## get 95% bootstrap CI for sample size n_B
quantile(temp[,2], c(0.025, 0.975))
## repeat for gamma
round(quantile(temp[,3], c(0.025, 0.975)),4)

## output the red and green posterior probabilities for each of the three sample
## sizes n^((0)) = small, n^((1)) = mid, and n^((2)) = large
write.csv(temp[,10:8201], "wei_sobol_8192_H0_small.csv", row.names = FALSE)
write.csv(temp[,8202:16393], "wei_sobol_8192_H1_small.csv", row.names = FALSE)
write.csv(temp[,16394:24585], "wei_sobol_8192_H0_mid.csv", row.names = FALSE)
write.csv(temp[,24586:32777], "wei_sobol_8192_H1_mid.csv", row.names = FALSE)
write.csv(temp[,32778:40969], "wei_sobol_8192_H0_large.csv", row.names = FALSE)
write.csv(temp[,40970:49161], "wei_sobol_8192_H1_large.csv", row.names = FALSE)

## repeat the 1000 sample size calculations with the SAME Sobol' sequences but
## explore the entire hypercube and each sample size considered
tempFull <- foreach(i=1:1000, .combine='rbind', .packages = c("qrng", "nleqslv"),
                    .options.snow=opts, .errorhandling = "remove") %dopar% {
                      temp_res <- findGammaFull(eta_plus1 = eta_plus1s, eta_plus0 = eta_plus0s, 
                                                pwr = 0.9, typeI = 0.01, deltas = c(-Inf, log(0.9)), q = 2, 
                                                alphas1, betas1, alphas2, betas2, m = 8192,
                                                seed_H1 = i + 8000, seed_H0 = i + 9000)
                      
                      c(i, as.numeric(unlist(temp_res)))
                    }

write.csv(tempFull, "wei_sobol_8192_full_summary.csv", row.names = FALSE)

## check if (n, gamma) combinations are the same
sum(temp[,2] == tempFull[,2])
sum(temp[,3] == tempFull[,3])

## comparison precision of Sobol' sequences with precision of 
## longer pseudorandom ones
temp <- foreach(i=1:1000, .combine='rbind', .packages = c("qrng", "nleqslv"),
                .options.snow=opts, .errorhandling = "remove") %dopar% {
                  temp_res <-  findGamma(eta_plus1 = eta_plus1s, eta_plus0 = eta_plus0s, 
                                         pwr = 0.9, typeI = 0.01, deltas = c(-Inf, log(0.9)), q = 2, 
                                         alphas1, betas1, alphas2, betas2, m = 75000, prng = TRUE,
                                         seed_H1 = i + 10000, seed_H0 = i + 11000, contour = FALSE)
                  
                  c(i, as.numeric(unlist(temp_res)))
                }

write.csv(temp, "wei_prng_75k_summary.csv", row.names = FALSE)

## get 95% bootstrap CI for sample size n_B
quantile(temp[,2], c(0.025, 0.975))
## repeat for gamma
round(quantile(temp[,3], c(0.025, 0.975)),4)

## now we get the sample size calculations for fixed gamma
temp <- foreach(i=1:1000, .combine='rbind', .packages = c("qrng", "nleqslv"),
                .options.snow=opts, .errorhandling = "remove") %dopar% {
                  temp_res <-  findn(eta_plus1 = eta_plus1s, eta_plus0 = eta_plus0s, 
                                         pwr = 0.9, typeI = 0.01, deltas = c(-Inf, log(0.9)), q = 2, 
                                         alphas1, betas1, alphas2, betas2, m = 8192,
                                         seed_H1 = i + 8000, seed_H0 = i + 9000, contour = FALSE)
                  
                  c(i, as.numeric(unlist(temp_res)))
                }

## output the optimal design and summary for each simulation repetition
write.csv(temp, "wei_sobol_8192_summary_fixed0p99.csv", row.names = FALSE) 

## calculate sample size recommendation with fixed gamma
mean(temp[,2]) ## should be 307

## create matrices for left contour plot in Figure E.1 (based on the average
## of 1000 sample size calculations); we do not have a contour plot for the
## single sample size calculation in this figure

## read in the sample sizes and posterior probabilities (logit scale) for 
## all 1000 repetitions
n_lows <- as.numeric(read.csv("wei_sobol_8192_summary.csv")[,7])
n_mids <- as.numeric(read.csv("wei_sobol_8192_summary.csv")[,8])
n_highs <- as.numeric(read.csv("wei_sobol_8192_summary.csv")[,9])

H1_mids <- read.csv("wei_sobol_8192_H1_mid.csv")
H0_mids <- read.csv("wei_sobol_8192_H0_mid.csv")
H1_lows <- read.csv("wei_sobol_8192_H1_small.csv")
H0_lows <- read.csv("wei_sobol_8192_H0_small.csv")
H1_highs <- read.csv("wei_sobol_8192_H1_large.csv")
H0_highs <- read.csv("wei_sobol_8192_H0_large.csv")

z_full <- matrix(0, nrow = 51, ncol = 50)
w_full <- matrix(0, nrow = 51, ncol = 50)
for (k in 1:1000){
  print(k)
  ## repeat the process from 04-multinomial-study-alg3.R with Weibull model
  H1_low <- unlist(H1_lows[k,]); H1_mid <- unlist(H1_mids[k,]); 
  H1_high <- unlist(H1_highs[k,])
  
  H0_low <- unlist(H0_lows[k,]); H0_mid <- unlist(H0_mids[k,]); 
  H0_high <- unlist(H0_highs[k,])
  
  n_low <- n_lows[k]; n_mid <- n_mids[k]; n_high <- n_highs[k]
  
  H1_slope <- (H1_mid - H1_low)/(n_mid-n_low)
  H0_slope <- (H0_mid - H0_low)/(n_mid-n_low)
  H1_int <- H1_low - H1_slope*n_low
  H0_int <- H0_low - H0_slope*n_low
  
  for (i in seq(300, as.numeric(n_mid))){
    assign(paste0("H1_vec_", i), H1_int + H1_slope*i)
    assign(paste0("H0_vec_", i), H0_int + H0_slope*i)
  }
  
  H1_slope <- (H1_high - H1_mid)/(n_high-n_mid)
  H0_slope <- (H0_high - H0_mid)/(n_high-n_mid)
  H1_int <- H1_mid - H1_slope*n_mid
  H0_int <- H0_mid - H0_slope*n_mid
  
  for (i in seq(as.numeric(n_mid) + 1, 350)){
    assign(paste0("H1_vec_", i), H1_int + H1_slope*i)
    assign(paste0("H0_vec_", i), H0_int + H0_slope*i)
  }
  
  gammas <- seq(log(0.988) - log(1 - 0.988), log(0.997) - log(1 - 0.997), length.out = 50)
  
  x <- seq(300, 350, 1)
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
write.csv(z_full, "wei_z_full_mat.csv", row.names = FALSE)
write.csv(w_full, "wei_w_full_mat.csv", row.names = FALSE)

## estimate operating characteristics with 40960 simulation repetitions
## using data simulated according to Psi0 or Psi1
resamps <- 40960
pb <- txtProgressBar(max = resamps, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

## get parameter values from Psi1 (degenerate) to generate
## the Weibull data under H1
## sample size of 307 is needed to get estimates for 
## operating characteristics at that sample size
samps <- c(seq(300, 350, 5), 307)

delta_L <- 0
delta_U <- 0.9

for (k in 1:length(samps)){
  probs_in <- foreach(j=1:resamps, .combine='rbind', .packages = c("rjags", "coda"),
                      .options.snow=opts, .errorhandling = "remove") %dopar% {
                        
                        y_star1 <- rweibull(q*samps[k], eta_plus1s[2], eta_plus1s[1])
                        y_star2 <- rweibull(samps[k], eta_plus1s[4], eta_plus1s[3])
                        
                        theta.diff <- WeibullPost(y_star1, y_star2, informs[1,1], informs[1,2], 
                                                  informs[2,1], informs[2,2],
                                                  informs[3,1], informs[3,2], 
                                                  informs[4,1], informs[4,2], ndraws = 25000)
                        
                        mean(ifelse(theta.diff > delta_L, theta.diff < delta_U,0))
                      }
  
  ## output results to .csv file
  write.csv(probs_in, paste0("probs_H1_appe_", samps[k], ".csv"), row.names = FALSE)
}

## get parameter values from Psi0 (degenerate) to generate
## the Weibull data under H0
for (k in 1:length(samps)){
  
  probs_out <- foreach(j=1:resamps, .combine='rbind', .packages = c("rjags", "coda"),
                       .options.snow=opts, .errorhandling = "remove") %dopar% {
                         
                         y_star1 <- rweibull(q*samps[k], eta_plus0s[2], eta_plus0s[1])
                         y_star2 <- rweibull(samps[k], eta_plus0s[4], eta_plus0s[3])
                         
                         theta.diff <- WeibullPost(y_star1, y_star2, informs[1,1], informs[1,2], 
                                                   informs[2,1], informs[2,2],
                                                   informs[3,1], informs[3,2], 
                                                   informs[4,1], informs[4,2], ndraws = 25000)
                         
                         mean(ifelse(theta.diff > delta_L, theta.diff < delta_U,0))
                       }
  
  ## output results to .csv file
  write.csv(probs_out, paste0("probs_H0_appe_", samps[k], ".csv"), row.names = FALSE)
}

## create matrices for right contour plot in Figure E.1 (based on
## Weibull data simulated above)
z_full2 <- matrix(0, nrow = 11, ncol = 60)
w_full2 <- matrix(0, nrow = 11, ncol = 60)
## convert the posterior probabilities for each approximated sampling
## distribution to the logit scale (with error checking to ensure
## to logits are finite)
for (i in seq(300,350,5)){
  assign(paste0("H1_vec_", i), 
         pmin(pmax(as.numeric(unlist(read.csv(paste0("probs_H1_appe_", i,".csv")))),
                   .Machine$double.eps), 1 - 10^(-7)))
  assign(paste0("H0_vec_", i),
         pmin(pmax(as.numeric(unlist(read.csv(paste0("probs_H0_appe_", i,".csv")))),
                   .Machine$double.eps), 1 - 10^(-7)))
}

## this process mirrors what was done to create the z and w matrices in 
## the previous two plots but with the estimates obtained by simulating data
gammas <- seq(log(0.988) - log(1 - 0.988), log(0.997) - log(1 - 0.997), length.out = 50)

x <- seq(300, 350, 5)
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
write.csv(z_full2, "wei_z_full_mat2.csv", row.names = FALSE)
write.csv(w_full2, "wei_w_full_mat2.csv", row.names = FALSE)

## create the three contour plots and output as .pdf file for the article
pdf(file = "FigureE1.pdf",   # The directory you want to save the file in
    width = 6.5, 
    height = 6.5) 

par(mfrow=c(2,2), mar = c(3.75, 3.75, 2, 0.35) + 0.1, mgp=c(2.25,1,0))

x <- seq(300, 350, 1)
gammas <- seq(log(0.988) - log(1 - 0.988), log(0.997) - log(1 - 0.997), length.out = 50)
y <- 1/(1 + exp(-gammas))

z_full <- matrix(unlist(read.csv("wei_z_full_mat.csv")), nrow = 51, ncol = 50)
w_full <- matrix(unlist(read.csv("wei_w_full_mat.csv")), nrow = 51, ncol = 50)

contour(x, y, w_full, levels = c(seq(0.0025, 0.0100, 0.0025), seq(0.015, 0.02, 0.0025)), 
        xlab = expression(italic("n")[italic("B")]), 
        xlim = c(300,350), axes = FALSE, cex.lab = 1.25,
        ylab = expression(gamma),  main = "Type I Error Rate", labcex = 0.8, method = "edge")
contour(x, y, w_full, levels = c(0.01), col = "firebrick", add = TRUE, labcex = 0.8, method = "edge")
contour(x, y, z_full, levels = c(0.9), col = "seagreen", add = TRUE, labcex = 0.8, labels = "", method = "edge")
axis(side = 1, at = seq(300, 350, 10), cex.axis = 1.15)
axis(side = 2, at = seq(0.988, 0.997, 0.003), cex.axis = 1.15)
contour(x, y, w_full, levels = c(0.0125), col = "black", add = TRUE, method = "edge", labcex = 0.8, label = "")
box()

z_full2 <- matrix(unlist(read.csv("wei_z_full_mat2.csv")), nrow = 11, ncol = 50)
w_full2 <- matrix(unlist(read.csv("wei_w_full_mat2.csv")), nrow = 11, ncol = 50)

x <- seq(300, 350, 5)

contour(x, y, w_full2, levels = c(seq(0.0025, 0.0075, 0.0025), seq(0.0125, 0.02, 0.0025)), 
        xlab = expression(italic("n")[italic("B")]), xlim = c(300, 350), axes = FALSE, cex.lab = 1.25,
        ylab = expression(gamma),  main = "Type I Error Rate", labcex = 0.8, method = "edge")
contour(x, y, w_full2, levels = c(0.01), col = "firebrick", add = TRUE, labcex = 0.8, method = "edge")
contour(x, y, z_full2, levels = c(0.9), col = "seagreen", add = TRUE, labcex = 0.8, labels = "", method = "edge")
axis(side = 1, at = seq(300, 350, 10), cex.axis = 1.15)
axis(side = 2, at = seq(0.988, 0.997, 0.003), cex.axis = 1.15)
box()

x <- seq(300, 350, 1)
contour(x, y, z_full, levels = c(seq(0.795, 0.87, 0.015), seq(0.915, 0.96, 0.015)), 
        xlab = expression(italic("n")[italic("B")]), xlim = c(300, 350), axes = FALSE,
        ylab = expression(gamma),  main = "Power", labcex = 0.8, method = "edge", cex.lab = 1.25)
axis(side = 1, at = seq(100, 120, 5), cex.axis = 1.15)
axis(side = 2, at = seq(0.91, 0.955, 0.015), cex.axis = 1.15)
contour(x, y, w_full, levels = c(0.01), col = "firebrick", add = TRUE, labcex = 0.8,labels = "", method = "edge") 
contour(x, y, z_full, levels = c(0.9), col = "seagreen", add = TRUE, labcex = 0.8, method = "edge")
contour(x, y, z_full, levels = c(0.885), col = "black", add = TRUE, method = "edge", labcex = 0.8, label = "")
axis(side = 1, at = seq(300, 350, 10), cex.axis = 1.15)
axis(side = 2, at = seq(0.988, 0.997, 0.003), cex.axis = 1.15)
box()

x <- seq(300, 350, 5)
contour(x, y, z_full2, levels = c(seq(0.795, 0.87, 0.015), seq(0.915, 0.96, 0.015)), 
        xlab = expression(italic("n")[italic("B")]), xlim = c(300,350), axes = FALSE, cex.lab = 1.25,
        ylab = expression(gamma),  main = "Power", labcex = 0.8, method = "edge")
contour(x, y, w_full2, levels = c(0.01), col = "firebrick", add = TRUE, labcex = 0.8,labels = "", method = "edge") 
contour(x, y, z_full2, levels = c(0.9), col = "seagreen", add = TRUE, labcex = 0.8, method = "edge")
contour(x, y, z_full2, levels = c(0.885), col = "black", add = TRUE, method = "edge", labcex = 0.8, label = "")
axis(side = 1, at = seq(300, 350, 10), cex.axis = 1.15)
axis(side = 2, at = seq(0.988, 0.997, 0.003), cex.axis = 1.15)
box()

par(mfrow=c(1,1))
dev.off()