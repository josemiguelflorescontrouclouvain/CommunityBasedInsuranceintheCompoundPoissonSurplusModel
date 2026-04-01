# Community-Based Insurance in the Compound Poisson Surplus Model
# R Code to calculate the ruin probability when each participant's claim sizes are LogNormally distributed
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

# The function below is in charge of simulating random losses that follow a mixture of LogNormal distributions.
# Recall that, under the assumption in which each participant's losses follow an Exponential distribution,
# then, the losses after pooling for each participant "i" follow an Hyperexponential distribution with 
# with parameter \beta_{ij} := \alpha_{i}/M_{ij}. The inputs for this function are:
# n: a numeric integer that provides the total number of random numbers with an Hyperexponential 
# distribution to be generated.
# probs: a vector containing the weights allocated to each of the Exponential distributions.
# mu: a vector containing the mu parameters for each of the LogNormal distributions.
# sigma: a vector containing the sigma parameters for each of the LogNormal distributions.

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

# The function below is in charge of simulating random losses that follow a scaled mixture of LogNormal distributions.
# Recall that, under the assumption in which each participant's losses follow an Exponential distribution,
# then, the losses after pooling for each participant "i" follow an Hyperexponential distribution with 
# with parameter \beta_{ij} := \alpha_{i}/M_{ij}. The inputs for this function are:
# n: a numeric integer that provides the total number of random numbers with an Hyperexponential 
# distribution to be generated.
# probs: a vector containing the weights allocated to each of the Exponential distributions.
# mu: a vector containing the mu parameters for each of the LogNormal distributions.
# sigma: a vector containing the sigma parameters for each of the LogNormal distributions.
# M_i: a vector containing the transfer rations M_i for the "i"th participant.

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

# The function below is in charge of simulating the ruin probability for the Cramér-Lundberg model 
# for participant "i" before pooling has taken place. The inputs for this function are:
# InitialCapital: a number indicating the initial capital for participant "i" (\kappa_{i})
# c: the premium rate for participant "i" (c_{i} = (1 + \eta) * E[S_{i,1}] = (1 + \eta) * \lambda_{i} * E[Y_{i,1}])
# lambda: the intensity of the compound Poission process. Recall that this parameter somehow describes us the frequency
# of the losses experienced by the participant (\lambda_{i})
# mu: the mu parameter of the LogNormal distribution.
# sigma: the sigma parameter of the LogNormal distribution.
# T_max: is the finite time horizon over which we simulate the Cramér–Lundberg risk process
# Simulations: represents the number of independent surplus paths we simulate in order to estimate the ruin probability.

SimulateRuinProbabilityIndividual <- function(InitialCapital = 0, c = 5.6, lambda = 2, mu = 0, sigma = 0.25, T_max = 1000, Simulations = 10000) {
  
  ruin_flags <- logical(Simulations) # create a vector with logical entries (TRUE or FALSE) that has lenght equal to the number of simulations
  
  # we start our simulations
  for (k in 1:Simulations) {
    n_claims <- rpois(1, lambda * T_max) # we simulate the total number of claims (frequency) over the time interval [0, T]
    if (n_claims == 0) { 
      ruin_flags[k] <- FALSE # if there is not claims, then ruin is impossible and thus we assign FALSE to that entry of the ruin_flags vector
      next
    }
    claim_times <- sort(runif(n_claims, 0, T_max))  # we generate the exact times when the n_claims occur; these times are uniformly distributed over the time interval [0, T]. That is, the have the same probability of occurring at any time in this interval.
    claim_sizes <- rlnorm(n_claims, meanlog = mu, sdlog = sigma) # we generate the claim sizes which will follow a Lognormal distribution.
    capital <- InitialCapital + c * claim_times - cumsum(claim_sizes) # we calculate the capital values immediately after each claim
    ruin_flags[k] <- any(capital <= 0) # since ruin can only occur at claim times, then we just need to now check if capital <= 0 for any of its entries
  }
  mean(ruin_flags) # lastly, the estimate of the ruin probability is given by the sum of TRUEs divided by the total number of simulations; in other words, we can also calculate this by estimating the mean of the ruin_flags vector
}

# We define general parameters
T_max <- 1000
Simulations <- 100000
eta <- 2/5
lambda_1 <- 2
mu_1 <- 0
sigma_1 <- 0.25 
c_1 <- (1 + eta) * lambda_1 * exp(mu_1 + sigma_1^2/2)
lambda_2 <- 1
mu_2 <- 0.5
sigma_2 <- 0.25 
c_2 <- (1 + eta) * lambda_2 * exp(mu_2 + sigma_2^2/2)
lambda_3 <- 3
mu_3 <- -0.3
sigma_3 <- 0.25 
c_3 <- (1 + eta) * lambda_3 * exp(mu_3 + sigma_3^2/2)
# Example 1. Estimation of ruin probability for participant 1
SimulateRuinProbabilityIndividual(InitialCapital = 1.4, c = c_1, lambda = lambda_1, mu = mu_1, sigma = sigma_1, T_max = T_max, Simulations = Simulations)
# Example 2. Estimation of ruin probability for participant 2
SimulateRuinProbabilityIndividual(InitialCapital = 0, c = c_2, lambda = lambda_2, mu = mu_2, sigma = sigma_2, T_max = T_max, Simulations = Simulations)
# Example 3. Estimation of ruin probability for participant 3
SimulateRuinProbabilityIndividual(InitialCapital = 0, c = c_3, lambda = lambda_3, mu = mu_3, sigma = sigma_3, T_max = T_max, Simulations = Simulations)

# The function below is in charge of simulating the ruin probability for the Cramér-Lundberg model 
# for participant "i" after pooling has taken place. The inputs for this function are:
# InitialCapital: a number indicating the initial capital for participant "i" (\kappa_{i})
# c: the premium rate for participant "i" (c_{i} = (1 + \eta) * E[S_{i,1}] = (1 + \eta) * \lambda_{i} * E[Y_{i,1}])
# lambda: the intensity of the compound Poission process. Recall that this parameter somehow describes us the frequency
# of the losses experienced in the pool (\lambda = \lambda_{1} + \lambda_{2} + ... + \lambda_{n})
# probs: a vector containing the weights allocated to each of the Lognormal distributions (this vector contains the values \lamda_{i}/\lambda).
# mus: a vector containing the mu parameters for each of the LogNormal distributions.
# sigmas: a vector containing the sigma parameters for each of the LogNormal distributions.
# T_max: is the finite time horizon over which we simulate the Cramér–Lundberg risk process
# Simulations: represents the number of independent surplus paths we simulate in order to estimate the ruin probability.

SimulateRuinProbabilityPool <- function(InitialCapital = 0, c = 5.6, lambda = lambda_1 + lambda_2 + lambda_3, probs = c(2/3, 1/3), mus = c(0, 0.5, -0.3), sigmas = c(0.25, 0.25, 0.25), M_is = c(8/15, 8/15, 8/15), T_max = 1000, Simulations = 10000) {
  
  ruin_flags <- logical(Simulations) # create a vector with logical entries (TRUE or FALSE) that has lenght equal to the number of simulations
  
  # we start our simulations
  for (k in 1:Simulations) {
    n_claims <- rpois(1, lambda * T_max) # we simulate the total number of claims (frequency) over the time interval [0, T]
    if (n_claims == 0) { 
      ruin_flags[k] <- FALSE # if there is not claims, then ruin is impossible and thus we assign FALSE to that entry of the ruin_flags vector
      next
    }
    claim_times <- sort(runif(n_claims, 0, T_max))  # we generate the exact times when the n_claims occur; these times are uniformly distributed over the time interval [0, T]. That is, the have the same probability of occurring at any time in this interval.
    claim_sizes <- rmixture_scaled_lognormals(n_claims, probs, mus, sigmas, M_is)$x # we generate the claim sizes which will follow a Hyperexponential distribution.
    capital <- InitialCapital + c * claim_times - cumsum(claim_sizes) # we calculate the capital values immediately after each claim
    ruin_flags[k] <- any(capital <= 0) # since ruin can only occur at claim times, then we just need to now check if capital <= 0 for any of its entries
  }
  mean(ruin_flags) # lastly, the estimate of the ruin probability is given by the sum of TRUEs divided by the total number of simulations; in other words, we can also calculate this by estimating the mean of the ruin_flags vector
}

# We define general parameters
T_max <- 1000
Simulations <- 10000
eta <- 2/5
lambda_1 <- 2
lambda_2 <- 1
lambda_3 <- 3  
lambda <- lambda_1 + lambda_2 + lambda_3
probs <- c(lambda_1, lambda_2, lambda_3)/lambda
mus <- c(0, 0.5, -0.3)
sigmas <- c(0.25, 0.25, 0.25)
c_1 <- (1 + eta) * lambda_1 * exp(mu_1 + sigma_1^2/2)
c_2 <- (1 + eta) * lambda_2 * exp(mu_2 + sigma_2^2/2)
c_3 <- (1 + eta) * lambda_3 * exp(mu_3 + sigma_3^2/2)
M_1s <- c((lambda_1 * exp(mu_1 + sigma_1^2/2))/(lambda_1 * exp(mu_1 + sigma_1^2/2) + lambda_2 * exp(mu_2 + sigma_2^2/2) + lambda_3 * exp(mu_3 + sigma_3^2/2)), (lambda_1 * exp(mu_1 + sigma_1^2/2))/(lambda_1 * exp(mu_1 + sigma_1^2/2) + lambda_2 * exp(mu_2 + sigma_2^2/2) + lambda_3 * exp(mu_3 + sigma_3^2/2)), (lambda_1 * exp(mu_1 + sigma_1^2/2))/(lambda_1 * exp(mu_1 + sigma_1^2/2) + lambda_2 * exp(mu_2 + sigma_2^2/2) + lambda_3 * exp(mu_3 + sigma_3^2/2)))
M_2s = c((lambda_2 * exp(mu_2 + sigma_2^2/2))/(lambda_1 * exp(mu_1 + sigma_1^2/2) + lambda_2 * exp(mu_2 + sigma_2^2/2) + lambda_3 * exp(mu_3 + sigma_3^2/2)), (lambda_2 * exp(mu_2 + sigma_2^2/2))/(lambda_1 * exp(mu_1 + sigma_1^2/2) + lambda_2 * exp(mu_2 + sigma_2^2/2) + lambda_3 * exp(mu_3 + sigma_3^2/2)), (lambda_2 * exp(mu_2 + sigma_2^2/2))/(lambda_1 * exp(mu_1 + sigma_1^2/2) + lambda_2 * exp(mu_2 + sigma_2^2/2) + lambda_3 * exp(mu_3 + sigma_3^2/2)))
M_3s = c((lambda_3 * exp(mu_3 + sigma_3^2/2))/(lambda_1 * exp(mu_1 + sigma_1^2/2) + lambda_2 * exp(mu_2 + sigma_2^2/2) + lambda_3 * exp(mu_3 + sigma_3^2/2)), (lambda_3 * exp(mu_3 + sigma_3^2/2))/(lambda_1 * exp(mu_1 + sigma_1^2/2) + lambda_2 * exp(mu_2 + sigma_2^2/2) + lambda_3 * exp(mu_3 + sigma_3^2/2)), (lambda_3 * exp(mu_3 + sigma_3^2/2))/(lambda_1 * exp(mu_1 + sigma_1^2/2) + lambda_2 * exp(mu_2 + sigma_2^2/2) + lambda_3 * exp(mu_3 + sigma_3^2/2)))
# Example 1. Estimation of ruin probability for participant 1
SimulateRuinProbabilityPool(InitialCapital = 0.58, c = c_1, lambda = lambda, probs = probs, mus = mus, sigmas = sigmas, M_is = M_1s, T_max = T_max, Simulations = Simulations)
# Example 2. Estimation of ruin probability for participant 2
SimulateRuinProbabilityPool(InitialCapital = 1, c = c_2, lambda = lambda, probs = probs, mus = mus, sigmas = sigmas, M_is = M_2s, T_max = T_max, Simulations = Simulations)
# Example 3. Estimation of ruin probability for participant 3
SimulateRuinProbabilityPool(InitialCapital = 1, c = c_3, lambda = lambda, probs = probs, mus = mus, sigmas = sigmas, M_is = M_3s, T_max = T_max, Simulations = Simulations)

# Example 1. Estimation of ruin probability for participant 1, 2 and 3 (with a Pool of three individuals)
# Vector of initial capitals for participant 1
initial_capitals <- seq(0, 5, by = 0.1)
# Number of simulations
Simulations <- 100000
# Finite time horizon
T_max <- 1000
# Safety loading
eta <- 2/5
# Remaining parameters for  participant 1
lambda_1 <- 2
mu_1 <- 0
sigma_1 <- sqrt(1)
c_1 <- (1 + eta) * lambda_1 * exp(mu_1 + sigma_1^2/2)
# Remaining parameters for  participant 2
lambda_2 <- 1
mu_2 <- 0.5
sigma_2 <- sqrt(1)
c_2 <- (1 + eta) * lambda_2 * exp(mu_2 + sigma_2^2/2)
# Remaining parameters for  participant 3
lambda_3 <- 3
mu_3 <- -0.3
sigma_3 <- sqrt(1) 
c_3 <- (1 + eta) * lambda_3 * exp(mu_3 + sigma_3^2/2)
# Remaining parameters for  the pool
lambdas <- c(lambda_1, lambda_2, lambda_3)
lambda <- lambda_1 + lambda_2 + lambda_3
probs <- c(lambda_1, lambda_2, lambda_3)/lambda
mus <- c(mu_1, mu_2, mu_3)
sigmas <- c(sigma_1, sigma_2, sigma_3)

MP_1s <- c((lambda_1 * exp(mu_1 + sigma_1^2/2))/(lambda_1 * exp(mu_1 + sigma_1^2/2) + lambda_2 * exp(mu_2 + sigma_2^2/2) + lambda_3 * exp(mu_3 + sigma_3^2/2)), (lambda_1 * exp(mu_1 + sigma_1^2/2))/(lambda_1 * exp(mu_1 + sigma_1^2/2) + lambda_2 * exp(mu_2 + sigma_2^2/2) + lambda_3 * exp(mu_3 + sigma_3^2/2)), (lambda_1 * exp(mu_1 + sigma_1^2/2))/(lambda_1 * exp(mu_1 + sigma_1^2/2) + lambda_2 * exp(mu_2 + sigma_2^2/2) + lambda_3 * exp(mu_3 + sigma_3^2/2)))
MP_2s = c((lambda_2 * exp(mu_2 + sigma_2^2/2))/(lambda_1 * exp(mu_1 + sigma_1^2/2) + lambda_2 * exp(mu_2 + sigma_2^2/2) + lambda_3 * exp(mu_3 + sigma_3^2/2)), (lambda_2 * exp(mu_2 + sigma_2^2/2))/(lambda_1 * exp(mu_1 + sigma_1^2/2) + lambda_2 * exp(mu_2 + sigma_2^2/2) + lambda_3 * exp(mu_3 + sigma_3^2/2)), (lambda_2 * exp(mu_2 + sigma_2^2/2))/(lambda_1 * exp(mu_1 + sigma_1^2/2) + lambda_2 * exp(mu_2 + sigma_2^2/2) + lambda_3 * exp(mu_3 + sigma_3^2/2)))
MP_3s = c((lambda_3 * exp(mu_3 + sigma_3^2/2))/(lambda_1 * exp(mu_1 + sigma_1^2/2) + lambda_2 * exp(mu_2 + sigma_2^2/2) + lambda_3 * exp(mu_3 + sigma_3^2/2)), (lambda_3 * exp(mu_3 + sigma_3^2/2))/(lambda_1 * exp(mu_1 + sigma_1^2/2) + lambda_2 * exp(mu_2 + sigma_2^2/2) + lambda_3 * exp(mu_3 + sigma_3^2/2)), (lambda_3 * exp(mu_3 + sigma_3^2/2))/(lambda_1 * exp(mu_1 + sigma_1^2/2) + lambda_2 * exp(mu_2 + sigma_2^2/2) + lambda_3 * exp(mu_3 + sigma_3^2/2)))

ALT_1s <- c(0.400, 0.100, 0.465759)
ALT_2s <- c(0.539003, 0.300, 0.034241)
ALT_3s <- c(0.060997, 0.600, 0.500)

# Use parallel computing
cores <- detectCores() - 1
cl <- makeCluster(cores)
clusterExport(cl, varlist = c("SimulateRuinProbabilityPool", "SimulateRuinProbabilityIndividual", "rmixtureoflognormals", "rmixture_scaled_lognormals", "Simulations", "T_max", "eta", "lambda_1", "mu_1", "sigma_1", "c_1", "MP_1s", "ALT_1s", "lambda_2", "mu_2", "sigma_2", "c_2", "MP_2s", "ALT_2s", "lambda_3", "mu_3", "sigma_3", "c_3", "MP_3s", "ALT_3s", "lambda", "probs", "mus", "sigmas"))

ruin_probsMP_pool_1 <- parSapply(cl, initial_capitals, function(IC) {SimulateRuinProbabilityPool(InitialCapital = IC, c = c_1, lambda = lambda, probs = probs, mus = mus, sigmas = sigmas, M_is = MP_1s, T_max = T_max, Simulations = Simulations)})
ruin_probsMP_pool_2 <- parSapply(cl, initial_capitals, function(IC) {SimulateRuinProbabilityPool(InitialCapital = IC, c = c_2, lambda = lambda, probs = probs, mus = mus, sigmas = sigmas, M_is = MP_2s, T_max = T_max, Simulations = Simulations)})
ruin_probsMP_pool_3 <- parSapply(cl, initial_capitals, function(IC) {SimulateRuinProbabilityPool(InitialCapital = IC, c = c_3, lambda = lambda, probs = probs, mus = mus, sigmas = sigmas, M_is = MP_3s, T_max = T_max, Simulations = Simulations)})

ruin_probsALT_pool_1 <- parSapply(cl, initial_capitals, function(IC) {SimulateRuinProbabilityPool(InitialCapital = IC, c = c_1, lambda = lambda, probs = probs, mus = mus, sigmas = sigmas, M_is = ALT_1s, T_max = T_max, Simulations = Simulations)})
ruin_probsALT_pool_2 <- parSapply(cl, initial_capitals, function(IC) {SimulateRuinProbabilityPool(InitialCapital = IC, c = c_2, lambda = lambda, probs = probs, mus = mus, sigmas = sigmas, M_is = ALT_2s, T_max = T_max, Simulations = Simulations)})
ruin_probsALT_pool_3 <- parSapply(cl, initial_capitals, function(IC) {SimulateRuinProbabilityPool(InitialCapital = IC, c = c_3, lambda = lambda, probs = probs, mus = mus, sigmas = sigmas, M_is = ALT_3s, T_max = T_max, Simulations = Simulations)})

ruin_probs_individual_1 <- parSapply(cl, initial_capitals, function(IC) {SimulateRuinProbabilityIndividual(InitialCapital = IC, c = c_1, lambda = lambda_1, mu = mu_1, sigma = sigma_1, T_max = T_max, Simulations = Simulations)})
ruin_probs_individual_2 <- parSapply(cl, initial_capitals, function(IC) {SimulateRuinProbabilityIndividual(InitialCapital = IC, c = c_2, lambda = lambda_2, mu = mu_2, sigma = sigma_2, T_max = T_max, Simulations = Simulations)})
ruin_probs_individual_3 <- parSapply(cl, initial_capitals, function(IC) {SimulateRuinProbabilityIndividual(InitialCapital = IC, c = c_3, lambda = lambda_3, mu = mu_3, sigma = sigma_3, T_max = T_max, Simulations = Simulations)})

stopCluster(cl)

file <- '/Users/josemiguelflorescontro/Documents/UCLouvain/Projects/Linear Risk Sharing/R/Graphs/Latex Codes to Generate Graphs'
setwd(file)

tikz('PlotInfiniteTimeRuinProbabilityExample53 (Monte Carlo) - All Participants MP.tex', standAlone = TRUE, width = 4, height = 4, packages = c("\\usepackage{tikz}", "\\usepackage[active,tightpage,psfixbb]{preview}", "\\PreviewEnvironment{pgfpicture}", "\\setlength\\PreviewBorder{0pt}", "\\usepackage{amssymb}", "\\usepackage{scalerel}"))
my.expressions <- c("$\\psi_{{\\scaleto{1}{3pt}}}^{{\\scaleto{ {\\fontfamily{qcr}\\selectfont POOL}}{3pt}}}\\left(\\kappa\\right)$", "$\\psi_{{\\scaleto{2}{3pt}}}^{{\\scaleto{ {\\fontfamily{qcr}\\selectfont POOL}}{3pt}}}\\left(\\kappa\\right)$", "$\\psi_{{\\scaleto{3}{3pt}}}^{{\\scaleto{ {\\fontfamily{qcr}\\selectfont POOL}}{3pt}}}\\left(\\kappa\\right)$", "$\\psi_{{\\scaleto{1}{3pt}}}\\left(\\kappa\\right)$", "$\\psi_{{\\scaleto{2}{3pt}}}\\left(\\kappa\\right)$", "$\\psi_{{\\scaleto{3}{3pt}}}\\left(\\kappa\\right)$")
par(mgp = c(2.5, 1, 0), mar = c(3.5, 3.5, 1, 1) + 0.1)
plot(initial_capitals, ruin_probsMP_pool_1, type = "l", lwd = 2, lty = 1, col = "blue", xaxs = "i", yaxs = "i", xlim = c(min(initial_capitals), max(initial_capitals)), ylim = c(0, 1), xlab = "$\\kappa$", ylab = "Ruin Probability")
grid(col = "gray70", lty = "dotted")
lines(initial_capitals, ruin_probsMP_pool_2, lwd = 2, lty = 1, col = "green")
lines(initial_capitals, ruin_probsMP_pool_3, lwd = 2, lty = 1, col = "red")
lines(initial_capitals, ruin_probs_individual_1, lwd = 2, lty = 2, col = "orange")
lines(initial_capitals, ruin_probs_individual_2, lwd = 2, lty = 3, col = "gray")
lines(initial_capitals, ruin_probs_individual_3, lwd = 2, lty = 4, col = "brown")
legend("topright", inset = 0.02, legend = my.expressions, col = c("blue", "green", "red", "orange", "gray", "brown"), lwd = c(2, 2, 2, 2, 2, 2), lty = c(1, 1, 1, 2, 3, 4), cex = 0.8)
dev.off()

tikz('PlotInfiniteTimeRuinProbabilityExample53 (Monte Carlo) - All Participants ALT.tex', standAlone = TRUE, width = 4, height = 4, packages = c("\\usepackage{tikz}", "\\usepackage[active,tightpage,psfixbb]{preview}", "\\PreviewEnvironment{pgfpicture}", "\\setlength\\PreviewBorder{0pt}", "\\usepackage{amssymb}", "\\usepackage{scalerel}"))
my.expressions <- c("$\\psi_{{\\scaleto{1}{3pt}}}^{{\\scaleto{ {\\fontfamily{qcr}\\selectfont POOL}}{3pt}}}\\left(\\kappa\\right)$", "$\\psi_{{\\scaleto{2}{3pt}}}^{{\\scaleto{ {\\fontfamily{qcr}\\selectfont POOL}}{3pt}}}\\left(\\kappa\\right)$", "$\\psi_{{\\scaleto{3}{3pt}}}^{{\\scaleto{ {\\fontfamily{qcr}\\selectfont POOL}}{3pt}}}\\left(\\kappa\\right)$", "$\\psi_{{\\scaleto{1}{3pt}}}\\left(\\kappa\\right)$", "$\\psi_{{\\scaleto{2}{3pt}}}\\left(\\kappa\\right)$", "$\\psi_{{\\scaleto{3}{3pt}}}\\left(\\kappa\\right)$")
par(mgp = c(2.5, 1, 0), mar = c(3.5, 3.5, 1, 1) + 0.1)
plot(initial_capitals, ruin_probsALT_pool_1, type = "l", lwd = 2, lty = 1, col = "blue", xaxs = "i", yaxs = "i", xlim = c(min(initial_capitals), max(initial_capitals)), ylim = c(0, 1), xlab = "$\\kappa$", ylab = "Ruin Probability")
grid(col = "gray70", lty = "dotted")
lines(initial_capitals, ruin_probsALT_pool_2, lwd = 2, lty = 1, col = "green")
lines(initial_capitals, ruin_probsALT_pool_3, lwd = 2, lty = 1, col = "red")
lines(initial_capitals, ruin_probs_individual_1, lwd = 2, lty = 2, col = "orange")
lines(initial_capitals, ruin_probs_individual_2, lwd = 2, lty = 3, col = "gray")
lines(initial_capitals, ruin_probs_individual_3, lwd = 2, lty = 4, col = "brown")
legend("topright", inset = 0.02, legend = my.expressions, col = c("blue", "green", "red", "orange", "gray", "brown"), lwd = c(2, 2, 2, 2, 2, 2), lty = c(1, 1, 1, 2, 3, 4), cex = 0.8)
dev.off()


