# Linear Risk Sharing in Community-Based Insurance: Ruin Reduction in the Compound Poisson Model 
# R Code to calculate the ruin probability when each participant's claim sizes are Exponentially distributed - Generate Figure 1 in Section 3.1
# Authors: Denuit, M., Flores-Contró, J. M. and Robert, C. Y.

rm(list = ls())

######################################################################################################################################
######################################################################################################################################

# We load the required packages.

library(parallel)
library(tikzDevice)
options(tikzMetricPackages = c("\\usepackage[utf8]{inputenc}", "\\usepackage[T1]{fontenc}", "\\usetikzlibrary{calc}", "\\usepackage{amssymb}"))

######################################################################################################################################
######################################################################################################################################

# First, we create a function that generates a random number which is drawn from a mixture of LogNormally distributed random varaibles.

rmixtureoflognormals <- function(n, probs, mu, sigma) {
  
  if (length(probs) != length(mu) || length(mu) != length(sigma))
    stop("probs, mu, and sigma must have the same length")
  
  if (abs(sum(probs) - 1) > 1e-8)
    stop("probs must sum to 1")
  
  components <- sample(seq_along(probs), size = n, replace = TRUE, prob = probs)
  
  x <- rlnorm(n, meanlog = mu[components], sdlog = sigma[components])
  
  list(x = x, components = components)
}

# rmixtureoflognormals(3, probs = c(0.8, 0.05, 0.05, 0.05, 0.05), mu = c(0, 0.5, -0.3, 0.7, 1), sigma = c(0.25, 0.25, 0.25, 0.25, 0.25))

# Second, we create a function that generates a random number which is drawn from a mixture of scaled LogNormally distributed random varaibles. 
# Indeed, recall that we need this scaled version for our risk processes once pooling takes place. The scale factors are given by the entries 
# in the allocation matrix (i.e., the transfer ratios).

rmixture_scaled_lognormals <- function(n, probs, mu, sigma, M_i) {
  
  if (length(probs) != length(mu) ||
      length(mu) != length(sigma) ||
      length(mu) != length(M_i))
    stop("probs, mu, sigma, and M_i must have the same length")
  
  if (abs(sum(probs) - 1) > 1e-8)
    stop("probs must sum to 1")
  
  components <- sample(seq_along(probs), size = n, replace = TRUE, prob = probs)
  
  x <- M_i[components] * rlnorm(n, meanlog = mu[components], sdlog = sigma[components])
  
  list(x = x, components = components)
}

# rmixture_scaled_lognormals(3, probs = c(0.8, 0.05, 0.05, 0.05, 0.05), mu = c(0, 0.5, -0.3, 0.7, 1), sigma = c(0.25, 0.25, 0.25, 0.25, 0.25), M_i = c(0.1, 0.3, 0.4, 0.1, 0.1))

# Finally, we create a function that employs the Euler Maruyama method, a standard numerical procedure used to simulate the trajectory 
# (or paths) of stochastic processes. In this case, we are particularly interested in generating the path of three Cramér-Lundberg 
# processes. That is, we will generate six different plots; two plots per participant: (i) one before pooling takes place and (ii) 
# another plot that shows the Cramér-Lundberg path after pooling takes place. This will allow us to illustrate how pooling affects 
# the claim experience (i.e., claim arrival times and claim sizes).

EulerMaruyamaMethod <- function(T = 20, eta = 2/5, kappa_1 = 15, lambda_1 = 0.25, mu_1 = 0, sigma_1 = 0.25, kappa_2 = 10, lambda_2 = 0.5, mu_2 = 0.5, sigma_2 = 0.25, kappa_3 = 7, lambda_3 = 0.75, mu_3 = -0.3, sigma_3 = 0.25, file = '/Users/jose/Library/CloudStorage/OneDrive-UCL/Documents/Postdoc/Linear Risk Sharing Project/R/Graphs/Latex Codes to Generate Graphs'){
  
  # First, using the expected value principle, we estimate the premium rate for each participant.
  c_1 <- (1 + eta) * lambda_1 * exp(mu_1 + sigma_1^2/2)
  c_2 <- (1 + eta) * lambda_2 * exp(mu_2 + sigma_2^2/2)
  c_3 <- (1 + eta) * lambda_3 * exp(mu_3 + sigma_3^2/2)
  
  # Now, we compute the pool intensity and the mixing probabilities.  
  lambda <- lambda_1 + lambda_2 + lambda_3
  probs <- c(lambda_1, lambda_2, lambda_3)/lambda
  
  # Then, we calculate the transfer ratios; that is, the entries of our allocation matrix. 
  # It is worth noting that, here, we are considering the mean-proportional risk-sharing rule. 
  # Thus, transfer ratios are computed accordingly.
  M_1s <- c((lambda_1 * exp(mu_1 + sigma_1^2/2))/(lambda_1 * exp(mu_1 + sigma_1^2/2) + lambda_2 * exp(mu_2 + sigma_2^2/2) + lambda_3 * exp(mu_3 + sigma_3^2/2)), (lambda_1 * exp(mu_1 + sigma_1^2/2))/(lambda_1 * exp(mu_1 + sigma_1^2/2) + lambda_2 * exp(mu_2 + sigma_2^2/2) + lambda_3 * exp(mu_3 + sigma_3^2/2)), (lambda_1 * exp(mu_1 + sigma_1^2/2))/(lambda_1 * exp(mu_1 + sigma_1^2/2) + lambda_2 * exp(mu_2 + sigma_2^2/2) + lambda_3 * exp(mu_3 + sigma_3^2/2)))
  M_2s = c((lambda_2 * exp(mu_2 + sigma_2^2/2))/(lambda_1 * exp(mu_1 + sigma_1^2/2) + lambda_2 * exp(mu_2 + sigma_2^2/2) + lambda_3 * exp(mu_3 + sigma_3^2/2)), (lambda_2 * exp(mu_2 + sigma_2^2/2))/(lambda_1 * exp(mu_1 + sigma_1^2/2) + lambda_2 * exp(mu_2 + sigma_2^2/2) + lambda_3 * exp(mu_3 + sigma_3^2/2)), (lambda_2 * exp(mu_2 + sigma_2^2/2))/(lambda_1 * exp(mu_1 + sigma_1^2/2) + lambda_2 * exp(mu_2 + sigma_2^2/2) + lambda_3 * exp(mu_3 + sigma_3^2/2)))
  M_3s = c((lambda_3 * exp(mu_3 + sigma_3^2/2))/(lambda_1 * exp(mu_1 + sigma_1^2/2) + lambda_2 * exp(mu_2 + sigma_2^2/2) + lambda_3 * exp(mu_3 + sigma_3^2/2)), (lambda_3 * exp(mu_3 + sigma_3^2/2))/(lambda_1 * exp(mu_1 + sigma_1^2/2) + lambda_2 * exp(mu_2 + sigma_2^2/2) + lambda_3 * exp(mu_3 + sigma_3^2/2)), (lambda_3 * exp(mu_3 + sigma_3^2/2))/(lambda_1 * exp(mu_1 + sigma_1^2/2) + lambda_2 * exp(mu_2 + sigma_2^2/2) + lambda_3 * exp(mu_3 + sigma_3^2/2)))
  
  # We assign line styles and color for claim origins. This will allow us to differentiate
  # in our plots which participant originated the claim.
  # Participant 1: orange, dotted (lty 3)
  # Participant 2: teal, dashed (lty 2)
  # Participant 3: pink, dot-dash (lty 4)
  jump_col <- c("#E69F00", "#009E73", "#CC79A7")
  jump_lty <- c(3, 2, 4)
  
  # Set working directory and choose a seed so that we generate reproducible plots.
  setwd(file)
  set.seed(15)
  
  # EULER-MARUYAMA METHOD: BEFORE POOLING TAKES PLACE
  
  # Participant 1 - Initial computations
  # Simulate claims
  N <- rpois(1, lambda_1 * T) # Number of claims in interval (0, T)
  claim_times_1 <- sort(runif(N, 0, T)) # The probability of a claim occurring at an specific time is uniformly distributed over the interval (0, T)
  claim_sizes_1 <- rlnorm(N, meanlog = mu_1, sdlog = sigma_1) # Generate claim severities (or sizes)
  
  # Participant 1 - Generate Plot
  tikz('PlotCramerLundbergIndividual_1.tex', standAlone = TRUE, width = 4, height = 4, packages = c("\\usepackage{tikz}", "\\usepackage[active,tightpage,psfixbb]{preview}", "\\PreviewEnvironment{pgfpicture}", "\\setlength\\PreviewBorder{0pt}", "\\usepackage{amssymb}", "\\usepackage{amsmath}"))
  par(mgp = c(2.5, 1, 0), mar = c(3.5, 3.5, 1, 1) + 0.1)
  plot(NA, xlim = c(0, T), ylim = c(0, 20), xaxs = "i", yaxs = "i", cex.lab = 1, cex.axis = 1, xlab = "$t$", ylab = "$V_{1,t}$", yaxt = "n")
  axis(side = 2, at = c(0, 5, 10, 15, 20), labels = c("0", "5", "10", "15", "20"))
  reserve <- kappa_1 # Initial reserve 
  last_time <- 0
  
  for(i in seq_along(claim_times_1)){
    reserve_before <- reserve + c_1 * (claim_times_1[i] - last_time) # Reserve immediately before the claim
    segments(last_time, reserve, claim_times_1[i], reserve_before, lwd = 1, col = "blue") # Solid blue line, premium accumulation with slope c_1
    segments(claim_times_1[i], reserve_before, claim_times_1[i], reserve_before - claim_sizes_1[i], lwd = 1, lty = jump_lty[1], col = jump_col[1]) # Jumps in the path; i.e., claims occurrences, styled by originating participant (here, participant 1) 
    reserve <- reserve_before - claim_sizes_1[i]  # Update surplus process
    last_time <- claim_times_1[i] # Update time of last occurrence
  }
  segments(last_time, reserve, T, reserve + c_1*(T-last_time), lwd = 1, col = "blue") # Final blue line segment, premium accumulation with slope c_1
  dev.off()
  
  # Participant 2 - Initial computations
  # Simulate claims
  N <- rpois(1, lambda_2 * T) # Number of claims in interval (0, T)
  claim_times_2 <- sort(runif(N, 0, T)) # The probability of a claim occurring at an specific time is uniformly distributed over the interval (0, T)
  claim_sizes_2 <- rlnorm(N, meanlog = mu_2, sdlog = sigma_2) # Generate claim severities (or sizes)
  
  # Participant 2 - Generate Plot
  tikz('PlotCramerLundbergIndividual_2.tex', standAlone = TRUE, width = 4, height = 4, packages = c("\\usepackage{tikz}", "\\usepackage[active,tightpage,psfixbb]{preview}", "\\PreviewEnvironment{pgfpicture}", "\\setlength\\PreviewBorder{0pt}", "\\usepackage{amssymb}", "\\usepackage{amsmath}"))
  par(mgp = c(2.5, 1, 0), mar = c(3.5, 3.5, 1, 1) + 0.1)
  plot(NA, xlim = c(0, T), ylim = c(0, 20), xaxs = "i", yaxs = "i", cex.lab = 1, cex.axis = 1, xlab = "$t$", ylab = "$V_{2,t}$", yaxt = "n")
  axis(side = 2, at = c(0, 5, 10, 15, 20), labels = c("0", "5", "10", "15", "20"))
  reserve <- kappa_2 # Initial reserve
  last_time <- 0
  
  for(i in seq_along(claim_times_2)){
    reserve_before <- reserve + c_2 * (claim_times_2[i] - last_time) # Reserve immediately before the claim
    segments(last_time, reserve, claim_times_2[i], reserve_before, lwd = 1, col = "blue") # Solid blue line, premium accumulation with slope c_2
    segments(claim_times_2[i], reserve_before, claim_times_2[i], reserve_before - claim_sizes_2[i], lwd = 1, lty = jump_lty[2], col = jump_col[2]) # Jumps in the path; i.e., claims occurrences, styled by originating participant (here, participant 2) 
    reserve <- reserve_before - claim_sizes_2[i] # Update surplus process
    last_time <- claim_times_2[i] # Update time of last occurrence
  }
  segments(last_time, reserve, T, reserve + c_2*(T-last_time), lwd = 1, col = "blue") # Final blue line segment, premium accumulation with slope c_2
  dev.off()
  
  # Participant 3 - Initial computations
  # Simulate claims
  N <- rpois(1, lambda_3 * T) # Number of claims in interval (0, T)
  claim_times_3 <- sort(runif(N, 0, T)) # The probability of a claim occurring at an specific time is uniformly distributed over the interval (0, T)
  claim_sizes_3 <- rlnorm(N, meanlog = mu_3, sdlog = sigma_3) # Generate claim severities (or sizes)
  
  # Participant 3 - Generate Plot
  tikz('PlotCramerLundbergIndividual_3.tex', standAlone = TRUE, width = 4, height = 4, packages = c("\\usepackage{tikz}", "\\usepackage[active,tightpage,psfixbb]{preview}", "\\PreviewEnvironment{pgfpicture}", "\\setlength\\PreviewBorder{0pt}", "\\usepackage{amssymb}", "\\usepackage{amsmath}"))
  par(mgp = c(2.5, 1, 0), mar = c(3.5, 3.5, 1, 1) + 0.1)
  plot(NA, xlim = c(0, T), ylim = c(0, 20), xaxs = "i", yaxs = "i", cex.lab = 1, cex.axis = 1, xlab = "$t$", ylab = "$V_{3,t}$", yaxt = "n")
  axis(side = 2, at = c(0, 5, 10, 15, 20), labels = c("0", "5", "10", "15", "20"))
  reserve <- kappa_3 # Initial reserve
  last_time <- 0
  
  for(i in seq_along(claim_times_3)){
    reserve_before <- reserve + c_3 * (claim_times_3[i] - last_time) # Reserve immediately before the claim
    segments(last_time, reserve, claim_times_3[i], reserve_before, lwd = 1, col = "blue") # Solid blue line, premium accumulation with slope c_3
    segments(claim_times_3[i], reserve_before, claim_times_3[i], reserve_before - claim_sizes_3[i], lwd = 1, lty = jump_lty[3], col = jump_col[3]) # Jumps in the path; i.e., claims occurrences, styled by originating participant (here, participant 3) 
    reserve <- reserve_before - claim_sizes_3[i] # Update surplus process
    last_time <- claim_times_3[i] # Update time of last occurrence
  }
  segments(last_time, reserve, T, reserve + c_3*(T-last_time), lwd = 1, col = "blue") # Final blue line segment, premium accumulation with slope c_3
  dev.off()
  
  # EULER-MARUYAMA METHOD: AFTER POOLING TAKES PLACE
  
  claim_times <- c(claim_times_1, claim_times_2, claim_times_3) # Number of claims in the pool on the interval (0, T)
  
  # Now, we need to keep track of which participant originated each claim, in the same order as claim_times,
  # so that after sorting we know which color/line type to use for each jump (claim occurrence) in the pooled 
  # surplus processes.
  claim_origin <- c(rep(1, length(claim_times_1)), rep(2, length(claim_times_2)), rep(3, length(claim_times_3)))
  
  # We sort the occurrence times and origins
  idx <- order(claim_times)
  claim_times_sorted <- claim_times[idx]
  claim_origin_sorted <- claim_origin[idx]
  
  # Pooled participant 1 - Initial computations
  # Simulate claims
  claim_sizes_1 <- c(M_1s[1] * claim_sizes_1, M_1s[2] * claim_sizes_2, M_1s[3] * claim_sizes_3) # Generate claim severities (or sizes). Recall that these are now scaled by the transfer ratios
  claim_sizes_sorted_1 <- claim_sizes_1[idx] # We sort the claim severities 
  
  # Pooled participant 1 - Generate Plot
  tikz('PlotCramerLundbergPooled_1.tex', standAlone = TRUE, width = 4, height = 4, packages = c("\\usepackage{tikz}", "\\usepackage[active,tightpage,psfixbb]{preview}", "\\PreviewEnvironment{pgfpicture}", "\\setlength\\PreviewBorder{0pt}", "\\usepackage{amssymb}", "\\usepackage{amsmath}"))
  par(mgp = c(2.5, 1, 0), mar = c(3.5, 3.5, 1, 1) + 0.1)
  plot(NA, xlim = c(0, T), ylim = c(0, 20), xaxs = "i", yaxs = "i", cex.lab = 1, cex.axis = 1, xlab = "$t$", ylab = "$V^{\\text{pool}}_{1,t}$", yaxt = "n")
  axis(side = 2, at = c(0, 5, 10, 15, 20), labels = c("0", "5", "10", "15", "20"))
  reserve <- kappa_1 # Initial reserve
  last_time <- 0
  
  for(i in seq_along(claim_times_sorted)){
    reserve_before <- reserve + c_1 * (claim_times_sorted[i] - last_time) # Reserve immediately before the claim
    segments(last_time, reserve, claim_times_sorted[i], reserve_before, lwd = 1, col = "blue") # Solid blue line, premium accumulation with slope c_1
    segments(claim_times_sorted[i], reserve_before, claim_times_sorted[i], reserve_before - claim_sizes_sorted_1[i], lwd = 1, lty = jump_lty[claim_origin_sorted[i]], col = jump_col[claim_origin_sorted[i]]) # Jumps in the path; color/line type depend on which participant originated the claim
    reserve <- reserve_before - claim_sizes_sorted_1[i] # Update surplus process
    last_time <- claim_times_sorted[i] # Update time of last occurrence
  }
  segments(last_time, reserve, T, reserve + c_1*(T-last_time), lwd = 1, col = "blue") # Final blue line segment, premium accumulation with slope c_1
  dev.off()
  
  # Pooled participant 2 - Initial computations
  # Simulate claims
  claim_sizes_2 <- c(M_2s[1] * claim_sizes_1, M_2s[2] * claim_sizes_2, M_2s[3] * claim_sizes_3) # Generate claim severities (or sizes). Recall that these are now scaled by the transfer ratios
  claim_sizes_sorted_2 <- claim_sizes_2[idx] # We sort the claim severities 
  
  # Pooled participant 2 - Generate Plot
  tikz('PlotCramerLundbergPooled_2.tex', standAlone = TRUE, width = 4, height = 4, packages = c("\\usepackage{tikz}", "\\usepackage[active,tightpage,psfixbb]{preview}", "\\PreviewEnvironment{pgfpicture}", "\\setlength\\PreviewBorder{0pt}", "\\usepackage{amssymb}", "\\usepackage{amsmath}"))
  par(mgp = c(2.5, 1, 0), mar = c(3.5, 3.5, 1, 1) + 0.1)
  plot(NA, xlim = c(0, T), ylim = c(0, 20), xaxs = "i", yaxs = "i", cex.lab = 1, cex.axis = 1, xlab = "$t$", ylab = "$V^{\\text{pool}}_{2,t}$", yaxt = "n")
  axis(side = 2, at = c(0, 5, 10, 15, 20), labels = c("0", "5", "10", "15", "20"))
  reserve <- kappa_2 # Initial reserve
  last_time <- 0
  
  for(i in seq_along(claim_times_sorted)){
    reserve_before <- reserve + c_2 * (claim_times_sorted[i] - last_time) # Reserve immediately before the claim
    segments(last_time, reserve, claim_times_sorted[i], reserve_before, lwd = 1, col = "blue") # Solid blue line, premium accumulation with slope c_2
    segments(claim_times_sorted[i], reserve_before, claim_times_sorted[i], reserve_before - claim_sizes_sorted_2[i], lwd = 1, lty = jump_lty[claim_origin_sorted[i]], col = jump_col[claim_origin_sorted[i]]) # Jumps in the path; color/line type depend on which participant originated the claim
    reserve <- reserve_before - claim_sizes_sorted_2[i] # Update surplus process
    last_time <- claim_times_sorted[i] # Update time of last occurrence
  }
  segments(last_time, reserve, T, reserve + c_2*(T-last_time), lwd = 1, col = "blue") # Final blue line segment, premium accumulation with slope c_2
  dev.off()
  
  # Pooled participant 3 - Initial computations
  # Simulate claims
  claim_sizes_3 <- c(M_3s[1] * claim_sizes_1, M_3s[2] * claim_sizes_2, M_3s[3] * claim_sizes_3) # Generate claim severities (or sizes). Recall that these are now scaled by the transfer ratios
  claim_sizes_sorted_3 <- claim_sizes_3[idx] # We sort the claim severities
  
  # Pooled participant 3 - Generate Plot
  tikz('PlotCramerLundbergPooled_3.tex', standAlone = TRUE, width = 4, height = 4, packages = c("\\usepackage{tikz}", "\\usepackage[active,tightpage,psfixbb]{preview}", "\\PreviewEnvironment{pgfpicture}", "\\setlength\\PreviewBorder{0pt}", "\\usepackage{amssymb}", "\\usepackage{amsmath}"))
  par(mgp = c(2.5, 1, 0), mar = c(3.5, 3.5, 1, 1) + 0.1)
  plot(NA, xlim = c(0, T), ylim = c(0, 20), xaxs = "i", yaxs = "i", cex.lab = 1, cex.axis = 1, xlab = "$t$", ylab = "$V^{\\text{pool}}_{3,t}$", yaxt = "n")
  axis(side = 2, at = c(0, 5, 10, 15, 20), labels = c("0", "5", "10", "15", "20"))
  reserve <- kappa_3 # Initial reserve
  last_time <- 0
  
  for(i in seq_along(claim_times_sorted)){
    reserve_before <- reserve + c_3 * (claim_times_sorted[i] - last_time) # Reserve immediately before the claim
    segments(last_time, reserve, claim_times_sorted[i], reserve_before, lwd = 1, col = "blue") # Solid blue line, premium accumulation with slope c_3
    segments(claim_times_sorted[i], reserve_before, claim_times_sorted[i], reserve_before - claim_sizes_sorted_3[i], lwd = 1, lty = jump_lty[claim_origin_sorted[i]], col = jump_col[claim_origin_sorted[i]]) # Jumps in the path; color/line type depend on which participant originated the claim
    reserve <- reserve_before - claim_sizes_sorted_3[i] # Update surplus process
    last_time <- claim_times_sorted[i] # Update time of last occurrence
  }
  segments(last_time, reserve, T, reserve + c_3*(T-last_time), lwd = 1, col = "blue") # Final blue line segment, premium accumulation with slope c_3
  dev.off()
  
}

EulerMaruyamaMethodResult <- EulerMaruyamaMethod(T = 20, eta = 2/5, kappa_1 = 5, lambda_1 = 0.3, mu_1 = 0.5, sigma_1 = 1, kappa_2 = 1, lambda_2 = 0.25, mu_2 = 0.45, sigma_2 = 1, kappa_3 = 4, lambda_3 = 0.2, mu_3 = 0.35, sigma_3 = 1, file = '/Users/jose/Library/CloudStorage/OneDrive-UCL/Documents/Postdoc/Linear Risk Sharing Project/R/Graphs/Latex Codes to Generate Graphs')


