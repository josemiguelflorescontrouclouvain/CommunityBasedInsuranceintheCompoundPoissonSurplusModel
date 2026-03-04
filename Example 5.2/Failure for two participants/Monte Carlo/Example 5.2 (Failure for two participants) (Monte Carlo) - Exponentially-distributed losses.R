# Community-Based Insurance in the Compound Poisson Surplus Model
# R Code to calculate the ruin probability when each participant's claim sizes are Exponentially distributed
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

# The function below is in charge of simulating random losses that follow an Hyperexponential distribution.
# Recall that, under the assumption in which each participant's losses follow an Exponential distribution,
# then, the losses after pooling for each participant "i" follow an Hyperexponential distribution with 
# with parameter \beta_{ij} := \alpha_{i}/M_{ij}. The inputs for this function are:
# n: a numeric integer that provides the total number of random numbers with an Hyperexponential 
# distribution to be generated.
# probs: a vector containing the weights allocated to each of the Exponential distributions.
# alphas: a vector containing the parameters for each of the Exponential distributions.

rmixtureofexponentials <- function(n, probs, alphas) {
  if(length(probs) != length(alphas)) stop("probs and alphas must have the same length") # probs and alphas need to have the same length
  if(abs(sum(probs) - 1) > 1e-8) stop("probs must sum to 1") # the sum of the elements in the vector probs needs to be equal to one
  
  components <- sample(seq_along(probs), size = n, replace = TRUE, prob = probs) # we randomly (obviously, considering the probability/weight of each exponential) select one of the Exponential distributions
  x <- rexp(n, rate = alphas[components])
  
  list(x = x, components = components) # we return the generated random numbers and also a vector containing the position of the Exponential distribution that was selected to generate that random number
}

# rmixtureofexponentials(3, probs = c(0.8, 0.05, 0.05, 0.05, 0.05), alphas = c(0.5, 2, 1, 1.5, 3))

# The function below is in charge of simulating the ruin probability for the Cramér-Lundberg model 
# for participant "i" before pooling has taken place. The inputs for this function are:
# InitialCapital: a number indicating the initial capital for participant "i" (\kappa_{i})
# c: the premium rate for participant "i" (c_{i} = (1 + \eta) * E[S_{i,1}] = (1 + \eta) * \lambda_{i} * E[Y_{i,1}])
# lambda: the intensity of the compound Poission process. Recall that this parameter somehow describes us the frequency
# of the losses experienced by the participant (\lambda_{i})
# alpha: the parameter of the Exponential distribution
# T_max: is the finite time horizon over which we simulate the Cramér–Lundberg risk process
# Simulations: represents the number of independent surplus paths we simulate in order to estimate the ruin probability.

SimulateRuinProbabilityIndividual <- function(InitialCapital = 0, c = 5.6, lambda = 2, alpha = 0.5, T_max = 1000, Simulations = 10000) {
  
  ruin_flags <- logical(Simulations) # create a vector with logical entries (TRUE or FALSE) that has lenght equal to the number of simulations
  
  # we start our simulations
  for (k in 1:Simulations) {
    n_claims <- rpois(1, lambda * T_max) # we simulate the total number of claims (frequency) over the time interval [0, T]
    if (n_claims == 0) { 
      ruin_flags[k] <- FALSE # if there is not claims, then ruin is impossible and thus we assign FALSE to that entry of the ruin_flags vector
      next
    }
    claim_times <- sort(runif(n_claims, 0, T_max))  # we generate the exact times when the n_claims occur; these times are uniformly distributed over the time interval [0, T]. That is, the have the same probability of occurring at any time in this interval.
    claim_sizes <- rexp(n_claims, alpha) # we generate the claim sizes which will follow a Hyperexponential distribution.
    capital <- InitialCapital + c * claim_times - cumsum(claim_sizes) # we calculate the capital values immediately after each claim
    ruin_flags[k] <- any(capital <= 0) # since ruin can only occur at claim times, then we just need to now check if capital <= 0 for any of its entries
  }
  mean(ruin_flags) # lastly, the estimate of the ruin probability is given by the sum of TRUEs divided by the total number of simulations; in other words, we can also calculate this by estimating the mean of the ruin_flags vector
}

# We define general parameters
T_max <- 10000
Simulations <- 10000
eta <- 2/5
lambda_1 <- 2
alpha_1 <- 1/2
c_1 <- (1 + eta) * lambda_1 * 1/alpha_1
lambda_2 <- 1
alpha_2 <- 2
c_2 <- (1 + eta) * lambda_2 * 1/alpha_2
lambda_3 <- 3
alpha_3 <- 1
c_3 <- (1 + eta) * lambda_3 * 1/alpha_3
# Example 1. Estimation of ruin probability for participant 1
SimulateRuinProbabilityIndividual(InitialCapital = 0, c = c_1, lambda = lambda_1, alpha = alpha_1, T_max = T_max, Simulations = Simulations)
# Example 2. Estimation of ruin probability for participant 2
SimulateRuinProbabilityIndividual(InitialCapital = 0, c = c_2, lambda = lambda_2, alpha = alpha_2, T_max = T_max, Simulations = Simulations)
# Example 3. Estimation of ruin probability for participant 3
SimulateRuinProbabilityIndividual(InitialCapital = 0, c = c_3, lambda = lambda_3, alpha = alpha_3, T_max = T_max, Simulations = Simulations)

# The function below is in charge of simulating the ruin probability for the Cramér-Lundberg model 
# for participant "i" after pooling has taken place. The inputs for this function are:
# InitialCapital: a number indicating the initial capital for participant "i" (\kappa_{i})
# c: the premium rate for participant "i" (c_{i} = (1 + \eta) * E[S_{i,1}] = (1 + \eta) * \lambda_{i} * E[Y_{i,1}])
# lambda: the intensity of the compound Poission process. Recall that this parameter somehow describes us the frequency
# of the losses experienced in the pool (\lambda = \lambda_{1} + \lambda_{2} + ... + \lambda_{n})
# probs: a vector containing the weights allocated to each of the Exponential distributions (this vector contains the values \lamda_{i}/\lambda).
# alphas: a vector containing the parameters for each of the Exponential distributions.
# T_max: is the finite time horizon over which we simulate the Cramér–Lundberg risk process
# Simulations: represents the number of independent surplus paths we simulate in order to estimate the ruin probability.

SimulateRuinProbabilityPool <- function(InitialCapital = 0, c = 5.6, lambda = 3, probs = c(2/3, 1/3), alphas = c(1, 1/2), T_max = 1000, Simulations = 10000) {
  
  ruin_flags <- logical(Simulations) # create a vector with logical entries (TRUE or FALSE) that has length equal to the number of simulations
  
  # we start our simulations
  for (k in 1:Simulations) {
    n_claims <- rpois(1, lambda * T_max) # we simulate the total number of claims (frequency) over the time interval [0, T]
    if (n_claims == 0) { 
      ruin_flags[k] <- FALSE # if there is not claims, then ruin is impossible and thus we assign FALSE to that entry of the ruin_flags vector
      next
    }
    claim_times <- sort(runif(n_claims, 0, T_max))  # we generate the exact times when the n_claims occur; these times are uniformly distributed over the time interval [0, T]. That is, the have the same probability of occurring at any time in this interval.
    claim_sizes <- rmixtureofexponentials(n_claims, probs, alphas)$x # we generate the claim sizes which will follow a Hyperexponential distribution.
    capital <- InitialCapital + c * claim_times - cumsum(claim_sizes) # we calculate the capital values immediately after each claim
    ruin_flags[k] <- any(capital <= 0) # since ruin can only occur at claim times, then we just need to now check if capital <= 0 for any of its entries
  }
  mean(ruin_flags) # lastly, the estimate of the ruin probability is given by the sum of TRUEs divided by the total number of simulations; in other words, we can also calculate this by estimating the mean of the ruin_flags vector
}

# We define general parameters
T_max <- 1000
Simulations <- 10000
# Safety loading
eta <- 2/5
# Remaining parameters for  participant 1
lambda_1 <- 2
alpha_1 <- 1/2
c_1 <- (1 + eta) * lambda_1 * 1/alpha_1
# Remaining parameters for  participant 2
lambda_2 <- 1
alpha_2 <- 2
c_2 <- (1 + eta) * lambda_2 * 1/alpha_2
# Remaining parameters for  participant 3
lambda_3 <- 3
alpha_3 <- 1
c_3 <- (1 + eta) * lambda_3 * 1/alpha_3
# Remaining parameters for  the pool 
lambda <- lambda_1 + lambda_2 + lambda_3
probs <- c(lambda_1, lambda_2, lambda_3)/lambda
MP_1s <- c((lambda_1 * 1/alpha_1)/(lambda_1 * 1/alpha_1 + lambda_2 * 1/alpha_2 + lambda_3 * 1/alpha_3), (lambda_1 * 1/alpha_1)/(lambda_1 * 1/alpha_1 + lambda_2 * 1/alpha_2 + lambda_3 * 1/alpha_3), (lambda_1 * 1/alpha_1)/(lambda_1 * 1/alpha_1 + lambda_2 * 1/alpha_2 + lambda_3 * 1/alpha_3))
MP_2s <- c((lambda_2 * 1/alpha_2)/(lambda_1 * 1/alpha_1 + lambda_2 * 1/alpha_2 + lambda_3 * 1/alpha_3), (lambda_2 * 1/alpha_2)/(lambda_1 * 1/alpha_1 + lambda_2 * 1/alpha_2 + lambda_3 * 1/alpha_3), (lambda_2 * 1/alpha_2)/(lambda_1 * 1/alpha_1 + lambda_2 * 1/alpha_2 + lambda_3 * 1/alpha_3))
MP_3s <- c((lambda_3 * 1/alpha_3)/(lambda_1 * 1/alpha_1 + lambda_2 * 1/alpha_2 + lambda_3 * 1/alpha_3), (lambda_3 * 1/alpha_3)/(lambda_1 * 1/alpha_1 + lambda_2 * 1/alpha_2 + lambda_3 * 1/alpha_3), (lambda_3 * 1/alpha_3)/(lambda_1 * 1/alpha_1 + lambda_2 * 1/alpha_2 + lambda_3 * 1/alpha_3))
betasMP_1 <- c(alpha_1, alpha_2, alpha_3)/MP_1s
betasMP_2 <- c(alpha_1, alpha_2, alpha_3)/MP_2s
betasMP_3 <- c(alpha_1, alpha_2, alpha_3)/MP_3s
# Example 1. Estimation of ruin probability for participant 1
SimulateRuinProbabilityPool(InitialCapital = 0, c = c_1, lambda = lambda, probs = probs, alphas = betasMP_1, T_max = 1000, Simulations = 10000)
# Example 2. Estimation of ruin probability for participant 2
SimulateRuinProbabilityPool(InitialCapital = 0, c = c_2, lambda = lambda, probs = probs, alphas = betasMP_2, T_max = 1000, Simulations = 10000)
# Example 3. Estimation of ruin probability for participant 3
SimulateRuinProbabilityPool(InitialCapital = 0, c = c_3, lambda = lambda, probs = probs, alphas = betasMP_3, T_max = 1000, Simulations = 10000)

# Example 1. Estimation of ruin probability for participant 1, 2 and 3 (with a pool of three individuals)
# Vector of initial capitals for participant 1
initial_capitals <- seq(0, 5, by = 0.1)
# We define general parameters
T_max <- 1000
Simulations <- 100000
# Safety loading
eta <- 2/5
# Remaining parameters for  participant 1
lambda_1 <- 100
alpha_1 <- 1
c_1 <- (1 + eta) * lambda_1 * 1/alpha_1
# Remaining parameters for  participant 2
lambda_2 <- 2
alpha_2 <- 1/50
c_2 <- (1 + eta) * lambda_2 * 1/alpha_2
# Remaining parameters for  participant 3
lambda_3 <- 100
alpha_3 <- 1
c_3 <- (1 + eta) * lambda_3 * 1/alpha_3
# Remaining parameters for  the pool 
lambda <- lambda_1 + lambda_2 + lambda_3
probs <- c(lambda_1, lambda_2, lambda_3)/lambda

MP_1s <- c((lambda_1 * 1/alpha_1)/(lambda_1 * 1/alpha_1 + lambda_2 * 1/alpha_2 + lambda_3 * 1/alpha_3), (lambda_1 * 1/alpha_1)/(lambda_1 * 1/alpha_1 + lambda_2 * 1/alpha_2 + lambda_3 * 1/alpha_3), (lambda_1 * 1/alpha_1)/(lambda_1 * 1/alpha_1 + lambda_2 * 1/alpha_2 + lambda_3 * 1/alpha_3))
MP_2s <- c((lambda_2 * 1/alpha_2)/(lambda_1 * 1/alpha_1 + lambda_2 * 1/alpha_2 + lambda_3 * 1/alpha_3), (lambda_2 * 1/alpha_2)/(lambda_1 * 1/alpha_1 + lambda_2 * 1/alpha_2 + lambda_3 * 1/alpha_3), (lambda_2 * 1/alpha_2)/(lambda_1 * 1/alpha_1 + lambda_2 * 1/alpha_2 + lambda_3 * 1/alpha_3))
MP_3s <- c((lambda_3 * 1/alpha_3)/(lambda_1 * 1/alpha_1 + lambda_2 * 1/alpha_2 + lambda_3 * 1/alpha_3), (lambda_3 * 1/alpha_3)/(lambda_1 * 1/alpha_1 + lambda_2 * 1/alpha_2 + lambda_3 * 1/alpha_3), (lambda_3 * 1/alpha_3)/(lambda_1 * 1/alpha_1 + lambda_2 * 1/alpha_2 + lambda_3 * 1/alpha_3))

ALT_1s <- c(0.5, 0.4, 0.1)
ALT_2s <- c(0.3, 0.2, 0.5)
ALT_3s <- c(0.2, 0.4, 0.4)

betasMP_1 <- c(alpha_1, alpha_2, alpha_3)/MP_1s
betasMP_2 <- c(alpha_1, alpha_2, alpha_3)/MP_2s
betasMP_3 <- c(alpha_1, alpha_2, alpha_3)/MP_3s

betasALT_1 <- c(alpha_1, alpha_2, alpha_3)/ALT_1s
betasALT_2 <- c(alpha_1, alpha_2, alpha_3)/ALT_2s
betasALT_3 <- c(alpha_1, alpha_2, alpha_3)/ALT_3s

# Use parallel computing
cores <- detectCores() - 1
cl <- makeCluster(cores)
clusterExport(cl, varlist = c("rmixtureofexponentials", "SimulateRuinProbabilityIndividual", "SimulateRuinProbabilityPool", "T_max", "Simulations", "lambda_1", "alpha_1", "c_1", "lambda_2", "alpha_2", "c_2", "lambda_3", "alpha_3", "c_3", "lambda", "probs", "MP_1s", "MP_2s", "MP_3s", "ALT_1s", "ALT_2s", "ALT_3s", "betasMP_1", "betasMP_2", "betasMP_3", "betasALT_1", "betasALT_2", "betasALT_3"))

ruin_probsMP_pool_1 <- parSapply(cl, initial_capitals, function(IC) {SimulateRuinProbabilityPool(InitialCapital = IC, c = c_1, lambda = lambda, probs = probs, alphas = betasMP_1, T_max = T_max, Simulations = Simulations)})
ruin_probsMP_pool_2 <- parSapply(cl, initial_capitals, function(IC) {SimulateRuinProbabilityPool(InitialCapital = IC, c = c_2, lambda = lambda, probs = probs, alphas = betasMP_2, T_max = T_max, Simulations = Simulations)})
ruin_probsMP_pool_3 <- parSapply(cl, initial_capitals, function(IC) {SimulateRuinProbabilityPool(InitialCapital = IC, c = c_3, lambda = lambda, probs = probs, alphas = betasMP_3, T_max = T_max, Simulations = Simulations)})

ruin_probsALT_pool_1 <- parSapply(cl, initial_capitals, function(IC) {SimulateRuinProbabilityPool(InitialCapital = IC, c = c_1, lambda = lambda, probs = probs, alphas = betasALT_1, T_max = T_max, Simulations = Simulations)})
ruin_probsALT_pool_2 <- parSapply(cl, initial_capitals, function(IC) {SimulateRuinProbabilityPool(InitialCapital = IC, c = c_2, lambda = lambda, probs = probs, alphas = betasALT_2, T_max = T_max, Simulations = Simulations)})
ruin_probsALT_pool_3 <- parSapply(cl, initial_capitals, function(IC) {SimulateRuinProbabilityPool(InitialCapital = IC, c = c_3, lambda = lambda, probs = probs, alphas = betasALT_3, T_max = T_max, Simulations = Simulations)})

ruin_probs_individual_1 <- parSapply(cl, initial_capitals, function(IC) {SimulateRuinProbabilityIndividual(InitialCapital = IC, c = c_1, lambda = lambda_1, alpha = alpha_1, T_max = T_max, Simulations = Simulations)})
ruin_probs_individual_2 <- parSapply(cl, initial_capitals, function(IC) {SimulateRuinProbabilityIndividual(InitialCapital = IC, c = c_2, lambda = lambda_2, alpha = alpha_2, T_max = T_max, Simulations = Simulations)})
ruin_probs_individual_3 <- parSapply(cl, initial_capitals, function(IC) {SimulateRuinProbabilityIndividual(InitialCapital = IC, c = c_3, lambda = lambda_3, alpha = alpha_3, T_max = T_max, Simulations = Simulations)})

stopCluster(cl)

file <- '/Users/josemiguelflorescontro/Documents/UCLouvain/Projects/Linear Risk Sharing/R/Graphs/Latex Codes to Generate Graphs'
setwd(file)

tikz('PlotInfiniteTimeRuinProbabilityExample52 (Version 2) (Monte Carlo) - All Participants MP.tex', standAlone = TRUE, width = 4, height = 4, packages = c("\\usepackage{tikz}", "\\usepackage[active,tightpage,psfixbb]{preview}", "\\PreviewEnvironment{pgfpicture}", "\\setlength\\PreviewBorder{0pt}", "\\usepackage{amssymb}", "\\usepackage{scalerel}"))
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

tikz('PlotInfiniteTimeRuinProbabilityExample52 (Version 2) (Monte Carlo) - All Participants ALT.tex', standAlone = TRUE, width = 4, height = 4, packages = c("\\usepackage{tikz}", "\\usepackage[active,tightpage,psfixbb]{preview}", "\\PreviewEnvironment{pgfpicture}", "\\setlength\\PreviewBorder{0pt}", "\\usepackage{amssymb}", "\\usepackage{scalerel}"))
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