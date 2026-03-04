# Community-Based Insurance in the Compound Poisson Surplus Model
# R Code to calculate the ruin probability when each participant's claim sizes are Exponentially distributed
# Authors: Denuit, M., Flores-Contró, J. M. and Robert, C. Y.

rm(list = ls())

######################################################################################################################################
######################################################################################################################################

# We load the required packages.

library(ggsci)
library(tikzDevice)
options(tikzMetricPackages = c("\\usepackage[utf8]{inputenc}", "\\usepackage[T1]{fontenc}", "\\usetikzlibrary{calc}", "\\usepackage{amssymb}"))

# We define the input parameters.

lambda_i <- c(100, 2, 100) # Frequencies (intensities) for each participant
lambda <- sum(lambda_i) # Sum of all Frequencies (intensities)
alpha_i  <- c(1, 1/50, 1) # Rate parameter of the exponential distribution for each participant

mu_i <- 1/alpha_i # Mean claim sizes

eta  <- 2/5 # Safety loading factor

MP_1s <- rep((lambda_i[1] * mu_i[1])/sum(lambda_i * mu_i), 3) # Mean-proportional (MP) risk-sharing scheme for participant 1
MP_2s <- rep((lambda_i[2] * mu_i[2])/sum(lambda_i * mu_i), 3) # Mean-proportional (MP) risk-sharing scheme for participant 2
MP_3s <- rep((lambda_i[3] * mu_i[3])/sum(lambda_i * mu_i), 3) # Mean-proportional (MP) risk-sharing scheme for participant 3

ALT_1s <- c(0.5, 0.4, 0.1) # Alternative (ALT) risk-sharing scheme for participant 1
ALT_2s <- c(0.3, 0.2, 0.5) # Alternative (ALT) risk-sharing scheme for participant 2
ALT_3s <- c(0.2, 0.4, 0.4) # Alternative (ALT) risk-sharing scheme for participant 3

betaMP_1s <- alpha_i/MP_1s # Corresponding rate parameters of the mixture of exponential distributions for participant 1
betaMP_2s <- alpha_i/MP_2s # Corresponding rate parameters of the mixture of exponential distributions for participant 2
betaMP_3s <- alpha_i/MP_3s # Corresponding rate parameters of the mixture of exponential distributions for participant 3

betaALT_1s <- alpha_i/ALT_1s # Corresponding rate parameters of the mixture of exponential distributions for participant 1
betaALT_2s <- alpha_i/ALT_2s # Corresponding rate parameters of the mixture of exponential distributions for participant 2
betaALT_3s <- alpha_i/ALT_3s # Corresponding rate parameters of the mixture of exponential distributions for participant 3

# We calculate the mean of Z_{i,j} for each participant

muZMP_1s <- sum(lambda_i/lambda * 1/betaMP_1s)
muZMP_2s <- sum(lambda_i/lambda * 1/betaMP_2s)
muZMP_3s <- sum(lambda_i/lambda * 1/betaMP_3s)

muZALT_1s <- sum(lambda_i/lambda * 1/betaALT_1s)
muZALT_2s <- sum(lambda_i/lambda * 1/betaALT_2s)
muZALT_3s <- sum(lambda_i/lambda * 1/betaALT_3s)

# We define G_Zi(r), the moment generating function of the random variable Z_{i,k} (which follows a mixture of Exponential distributions with rate parameters beta_{i,j}).

G_Zi <- function(r, lambda, lambda_vec, beta_vec) {
  sum_terms <- sum((beta_vec * lambda_vec) / (beta_vec - r))
  return((1 / lambda) * sum_terms)
}

# We define the denominator function from Equation (13.6.9) from Bowers et. al (1997).

denominator_func <- function(r, eta, mu_Zi, lambda, lambda_vec, beta_vec) {
  gz_r <- G_Zi(r, lambda, lambda_vec, beta_vec)
  val <- 1 + (1 + eta) * mu_Zi * r - gz_r
  return(val)
}

# We find the roots: r_{1}, r_{2} and r_{3}.

# First, we sort beta_vec to identify the intervals between singularities.

sorted_betasMP_1s <- sort(betaMP_1s)
sorted_betasMP_2s <- sort(betaMP_2s)
sorted_betasMP_3s <- sort(betaMP_3s)

sorted_betasALT_1s <- sort(betaALT_1s)
sorted_betasALT_2s <- sort(betaALT_2s)
sorted_betasALT_3s <- sort(betaALT_3s)

# Then, we define the three intervals. Note that we use a small epsilon (1e-10) to avoid the exact singularity points.

intervalsMP_1s <- list(
  c(1e-10, sorted_betasMP_1s[1] - 1e-10), # Interval for r_1
  c(sorted_betasMP_1s[1] + 1e-10, sorted_betasMP_1s[2] - 1e-10), # Interval for r_2
  c(sorted_betasMP_1s[2] + 1e-10, sorted_betasMP_1s[3] - 1e-10)  # Interval for r_3
)

intervalsMP_2s <- list(
  c(1e-10, sorted_betasMP_2s[1] - 1e-10), # Interval for r_1
  c(sorted_betasMP_2s[1] + 1e-10, sorted_betasMP_2s[2] - 1e-10), # Interval for r_2
  c(sorted_betasMP_2s[2] + 1e-10, sorted_betasMP_2s[3] - 1e-10)  # Interval for r_3
)

intervalsMP_3s <- list(
  c(1e-10, sorted_betasMP_3s[1] - 1e-10), # Interval for r_1
  c(sorted_betasMP_3s[1] + 1e-10, sorted_betasMP_3s[2] - 1e-10), # Interval for r_2
  c(sorted_betasMP_3s[2] + 1e-10, sorted_betasMP_3s[3] - 1e-10)  # Interval for r_3
)

intervalsALT_1s <- list(
  c(1e-10, sorted_betasALT_1s[1] - 1e-10), # Interval for r_1
  c(sorted_betasALT_1s[1] + 1e-10, sorted_betasALT_1s[2] - 1e-10), # Interval for r_2
  c(sorted_betasALT_1s[2] + 1e-10, sorted_betasALT_1s[3] - 1e-10)  # Interval for r_3
)

intervalsALT_2s <- list(
  c(1e-10, sorted_betasALT_2s[1] - 1e-10), # Interval for r_1
  c(sorted_betasALT_2s[1] + 1e-10, sorted_betasALT_2s[2] - 1e-10), # Interval for r_2
  c(sorted_betasALT_2s[2] + 1e-10, sorted_betasALT_2s[3] - 1e-10)  # Interval for r_3
)

intervalsALT_3s <- list(
  c(1e-10, sorted_betasALT_3s[1] - 1e-10), # Interval for r_1
  c(sorted_betasALT_3s[1] + 1e-10, sorted_betasALT_3s[2] - 1e-10), # Interval for r_2
  c(sorted_betasALT_3s[2] + 1e-10, sorted_betasALT_3s[3] - 1e-10)  # Interval for r_3
)

# We also define the vector in which the roots will be stored.

rootsMP_1s <- numeric(3)
rootsMP_2s <- numeric(3)
rootsMP_3s <- numeric(3)

rootsALT_1s <- numeric(3)
rootsALT_2s <- numeric(3)
rootsALT_3s <- numeric(3)

# Lastly, we compute all roots for each of the three participants

for(i in 1:3) {
  resMP_1s <- uniroot(
    f = denominator_func,
    interval = intervalsMP_1s[[i]],
    tol = .Machine$double.eps^0.75, # High precision tolerance
    eta = eta,
    mu_Zi = muZMP_1s,
    lambda = lambda,
    lambda_vec = lambda_i,
    beta_vec = betaMP_1s
  )
  rootsMP_1s[i] <- resMP_1s$root
}

for(i in 1:3) {
  resMP_2s <- uniroot(
    f = denominator_func,
    interval = intervalsMP_2s[[i]],
    tol = .Machine$double.eps^0.75, # High precision tolerance
    eta = eta,
    mu_Zi = muZMP_2s,
    lambda = lambda,
    lambda_vec = lambda_i,
    beta_vec = betaMP_2s
  )
  rootsMP_2s[i] <- resMP_2s$root
}

for(i in 1:3) {
  resMP_3s <- uniroot(
    f = denominator_func,
    interval = intervalsMP_3s[[i]],
    tol = .Machine$double.eps^0.75, # High precision tolerance
    eta = eta,
    mu_Zi = muZMP_3s,
    lambda = lambda,
    lambda_vec = lambda_i,
    beta_vec = betaMP_3s
  )
  rootsMP_3s[i] <- resMP_3s$root
}

for(i in 1:3) {
  resALT_1s <- uniroot(
    f = denominator_func,
    interval = intervalsALT_1s[[i]],
    tol = .Machine$double.eps^0.75, # High precision tolerance
    eta = eta,
    mu_Zi = muZALT_1s,
    lambda = lambda,
    lambda_vec = lambda_i,
    beta_vec = betaALT_1s
  )
  rootsALT_1s[i] <- resALT_1s$root
}

for(i in 1:3) {
  resALT_2s <- uniroot(
    f = denominator_func,
    interval = intervalsALT_2s[[i]],
    tol = .Machine$double.eps^0.75, # High precision tolerance
    eta = eta,
    mu_Zi = muZALT_2s,
    lambda = lambda,
    lambda_vec = lambda_i,
    beta_vec = betaALT_2s
  )
  rootsALT_2s[i] <- resALT_2s$root
}

for(i in 1:3) {
  resALT_3s <- uniroot(
    f = denominator_func,
    interval = intervalsALT_3s[[i]],
    tol = .Machine$double.eps^0.75, # High precision tolerance
    eta = eta,
    mu_Zi = muZALT_3s,
    lambda = lambda,
    lambda_vec = lambda_i,
    beta_vec = betaALT_3s
  )
  rootsALT_3s[i] <- resALT_3s$root
}

# Output all roots for each of the three participants.

print(paste("Root 1 under the Mean-proportional (MP) risk-sharing scheme for participant 1:", rootsMP_1s[1]))
print(paste("Root 2 under the Mean-proportional (MP) risk-sharing scheme for participant 1:", rootsMP_1s[2]))
print(paste("Root 3 under the Mean-proportional (MP) risk-sharing scheme for participant 1:", rootsMP_1s[3]))

print(paste("Root 1 under the Mean-proportional (MP) risk-sharing scheme for participant 2:", rootsMP_2s[1]))
print(paste("Root 2 under the Mean-proportional (MP) risk-sharing scheme for participant 2:", rootsMP_2s[2]))
print(paste("Root 3 under the Mean-proportional (MP) risk-sharing scheme for participant 2:", rootsMP_2s[3]))

print(paste("Root 1 under the Mean-proportional (MP) risk-sharing scheme for participant 3:", rootsMP_3s[1]))
print(paste("Root 2 under the Mean-proportional (MP) risk-sharing scheme for participant 3:", rootsMP_3s[2]))
print(paste("Root 3 under the Mean-proportional (MP) risk-sharing scheme for participant 3:", rootsMP_3s[3]))

print(paste("Root 1 under the Alternative (ALT) risk-sharing scheme for participant 1:", rootsALT_1s[1]))
print(paste("Root 2 under the Alternative (ALT) risk-sharing scheme for participant 1:", rootsALT_1s[2]))
print(paste("Root 3 under the Alternative (ALT) risk-sharing scheme for participant 1:", rootsALT_1s[3]))

print(paste("Root 1 under the Alternative (ALT) risk-sharing scheme for participant 2:", rootsALT_2s[1]))
print(paste("Root 2 under the Alternative (ALT) risk-sharing scheme for participant 2:", rootsALT_2s[2]))
print(paste("Root 3 under the Alternative (ALT) risk-sharing scheme for participant 2:", rootsALT_2s[3]))

print(paste("Root 1 under the Alternative (ALT) risk-sharing scheme for participant 3:", rootsALT_3s[1]))
print(paste("Root 2 under the Alternative (ALT) risk-sharing scheme for participant 3:", rootsALT_3s[2]))
print(paste("Root 3 under the Alternative (ALT) risk-sharing scheme for participant 3:", rootsALT_3s[3]))

# We find the coefficients: C_{1}, C_{2} and C_{3}.

find_coefficients <- function(roots, eta, mu_Zi, lambda, lambda_vec, beta_vec) {
  
  dG_dr <- function(r) {
    sum((beta_vec * lambda_vec) / (beta_vec - r)^2) / lambda
  }
  
  Ci_values <- numeric(length(roots))
  
  for (i in seq_along(roots)) {
    ri <- roots[i]

    gz_ri <- G_Zi(ri, lambda, lambda_vec, beta_vec)
    num <- (eta * (gz_ri - 1)) / (1 + eta)
    
    den_prime <- ((1 + eta) * mu_Zi) - dG_dr(ri)
  
    Ci_values[i] <- (num / -den_prime) / ri 
  }
  return(Ci_values)
}
 
CiMP_1s <- find_coefficients(rootsMP_1s, eta, muZMP_1s, lambda, lambda_i, betaMP_1s) # Run for participant 1
CiMP_2s <- find_coefficients(rootsMP_2s, eta, muZMP_2s, lambda, lambda_i, betaMP_2s) # Run for participant 2
CiMP_3s <- find_coefficients(rootsMP_3s, eta, muZMP_3s, lambda, lambda_i, betaMP_3s) # Run for participant 3

CiALT_1s <- find_coefficients(rootsALT_1s, eta, muZALT_1s, lambda, lambda_i, betaALT_1s) # Run for participant 1
CiALT_2s <- find_coefficients(rootsALT_2s, eta, muZALT_2s, lambda, lambda_i, betaALT_2s) # Run for participant 2
CiALT_3s <- find_coefficients(rootsALT_3s, eta, muZALT_3s, lambda, lambda_i, betaALT_3s) # Run for participant 3

# Output all coefficients for each of the three participants.

print(paste("Coefficient 1 under the Mean-proportional (MP) risk-sharing scheme for participant 1:", CiMP_1s[1]))
print(paste("Coefficient 2 under the Mean-proportional (MP) risk-sharing scheme for participant 1:", CiMP_1s[2]))
print(paste("Coefficient 3 under the Mean-proportional (MP) risk-sharing scheme for participant 1:", CiMP_1s[3]))

print(paste("Coefficient 1 under the Mean-proportional (MP) risk-sharing scheme for participant 2:", CiMP_2s[1]))
print(paste("Coefficient 2 under the Mean-proportional (MP) risk-sharing scheme for participant 2:", CiMP_2s[2]))
print(paste("Coefficient 3 under the Mean-proportional (MP) risk-sharing scheme for participant 2:", CiMP_2s[3]))

print(paste("Coefficient 1 under the Mean-proportional (MP) risk-sharing scheme for participant 3:", CiMP_3s[1]))
print(paste("Coefficient 2 under the Mean-proportional (MP) risk-sharing scheme for participant 3:", CiMP_3s[2]))
print(paste("Coefficient 3 under the Mean-proportional (MP) risk-sharing scheme for participant 3:", CiMP_3s[3]))

print(paste("Coefficient 1 under the Alternative (ALT) risk-sharing scheme for participant 1:", CiALT_1s[1]))
print(paste("Coefficient 2 under the Alternative (ALT) risk-sharing scheme for participant 1:", CiALT_1s[2]))
print(paste("Coefficient 3 under the Alternative (ALT) risk-sharing scheme for participant 1:", CiALT_1s[3]))

print(paste("Coefficient 1 under the Alternative (ALT) risk-sharing scheme for participant 2:", CiALT_2s[1]))
print(paste("Coefficient 2 under the Alternative (ALT) risk-sharing scheme for participant 2:", CiALT_2s[2]))
print(paste("Coefficient 3 under the Alternative (ALT) risk-sharing scheme for participant 2:", CiALT_2s[3]))

print(paste("Coefficient 1 under the Alternative (ALT) risk-sharing scheme for participant 3:", CiALT_3s[1]))
print(paste("Coefficient 2 under the Alternative (ALT) risk-sharing scheme for participant 3:", CiALT_3s[2]))
print(paste("Coefficient 3 under the Alternative (ALT) risk-sharing scheme for participant 3:", CiALT_3s[3]))

# We define the functions to calculate the ruin probabilities (with and without pooling).

ruin_probability_exponential <- function(kappa, lambda, alpha, eta){
  c <- (1 + eta) * lambda * 1/alpha
  ruin_probability <- lambda/(alpha * c) * exp( -(alpha - lambda/c) * kappa)
  return(ruin_probability)
}

ruin_probability_mixture_exponentials <- function(kappa, roots, coefficients){
  ruin_probability <- coefficients[1] * exp(-roots[1] * kappa) + coefficients[2] * exp(-roots[2] * kappa) + coefficients[3] * exp(-roots[3] * kappa)
  return(ruin_probability)
}

# Then, we calculate the ruin probabilities for each of the participants.

initial_capitals <- seq(0, 60, by = 0.01)

ruin_probs_1s <- sapply(initial_capitals, function(kappa) {ruin_probability_exponential(kappa, lambda_i[1], alpha_i[1], eta)})
ruin_probs_2s <- sapply(initial_capitals, function(kappa) {ruin_probability_exponential(kappa, lambda_i[2], alpha_i[2], eta)})
ruin_probs_3s <- sapply(initial_capitals, function(kappa) {ruin_probability_exponential(kappa, lambda_i[3], alpha_i[3], eta)})

ruin_probsMP_1s <- sapply(initial_capitals, function(kappa) {ruin_probability_mixture_exponentials(kappa, rootsMP_1s, CiMP_1s)})
ruin_probsMP_2s <- sapply(initial_capitals, function(kappa) {ruin_probability_mixture_exponentials(kappa, rootsMP_2s, CiMP_2s)})
ruin_probsMP_3s <- sapply(initial_capitals, function(kappa) {ruin_probability_mixture_exponentials(kappa, rootsMP_3s, CiMP_3s)})

ruin_probsALT_1s <- sapply(initial_capitals, function(kappa) {ruin_probability_mixture_exponentials(kappa, rootsALT_1s, CiALT_1s)})
ruin_probsALT_2s <- sapply(initial_capitals, function(kappa) {ruin_probability_mixture_exponentials(kappa, rootsALT_2s, CiALT_2s)})
ruin_probsALT_3s <- sapply(initial_capitals, function(kappa) {ruin_probability_mixture_exponentials(kappa, rootsALT_3s, CiALT_3s)})

# Lastly, we plot our results.

file <- '/Users/josemiguelflorescontro/Documents/UCLouvain/Projects/Linear Risk Sharing/R/Graphs/Latex Codes to Generate Graphs'
setwd(file)

tikz('PlotInfiniteTimeRuinProbabilityExample52 (Version 2) - Participant 1.tex', standAlone = TRUE, width = 4, height = 4, packages = c("\\usepackage{tikz}", "\\usepackage[active,tightpage,psfixbb]{preview}", "\\PreviewEnvironment{pgfpicture}", "\\setlength\\PreviewBorder{0pt}", "\\usepackage{amssymb}", "\\usepackage{scalerel}"))
par(mgp = c(2.5, 1, 0), mar = c(3.5, 3.5, 1, 1) + 0.1)
MyColors <- pal_jco()(3)
MyLines <-seq(from = 1, to = 3, by = 1)
plot(initial_capitals, ruin_probs_1s, type = "l", lwd = 2, lty = MyLines[1], col = MyColors[1], xaxs = "i", yaxs = "i", xlim = c(min(initial_capitals), max(initial_capitals)), ylim = c(0, 1), xlab = "$\\kappa_{1}$", ylab = "Ruin Probability")
lines(initial_capitals, ruin_probsMP_1s, lwd = 2, lty = MyLines[2], col = MyColors[2])
lines(initial_capitals, ruin_probsALT_1s, lwd = 2, lty = MyLines[3], col = MyColors[3])
my.expressions <- c("$\\psi_{1}\\left(\\kappa_{1}\\right)$", "$\\psi^{\\scaleto{{\\fontfamily{qcr}\\selectfont pool - MP}}{3pt}}_{1}\\left(\\kappa_{1}\\right)$", "$\\psi^{\\scaleto{{\\fontfamily{qcr}\\selectfont pool - ALT}}{3pt}}_{1}\\left(\\kappa_{1}\\right)$")
legend("topright", inset = 0.02, legend = my.expressions, lty = MyLines, lwd = 2, col = MyColors, cex = 0.8)
dev.off()

tikz('PlotInfiniteTimeRuinProbabilityExample52 (Version 2) - Participant 2.tex', standAlone = TRUE, width = 4, height = 4, packages = c("\\usepackage{tikz}", "\\usepackage[active,tightpage,psfixbb]{preview}", "\\PreviewEnvironment{pgfpicture}", "\\setlength\\PreviewBorder{0pt}", "\\usepackage{amssymb}", "\\usepackage{scalerel}"))
par(mgp = c(2.5, 1, 0), mar = c(3.5, 3.5, 1, 1) + 0.1)
MyColors <- pal_jco()(3)
MyLines <-seq(from = 1, to = 3, by = 1)
plot(initial_capitals, ruin_probs_2s, type = "l", lwd = 2, lty = MyLines[1], col = MyColors[1], xaxs = "i", yaxs = "i", xlim = c(min(initial_capitals), max(initial_capitals)), ylim = c(0, 1), xlab = "$\\kappa_{2}$", ylab = "Ruin Probability")
lines(initial_capitals, ruin_probsMP_2s, lwd = 2, lty = MyLines[2], col = MyColors[2])
lines(initial_capitals, ruin_probsALT_2s, lwd = 2, lty = MyLines[3], col = MyColors[3])
my.expressions <- c("$\\psi_{2}\\left(\\kappa_{2}\\right)$", "$\\psi^{\\scaleto{{\\fontfamily{qcr}\\selectfont pool - MP}}{3pt}}_{2}\\left(\\kappa_{2}\\right)$", "$\\psi^{\\scaleto{{\\fontfamily{qcr}\\selectfont pool - ALT}}{3pt}}_{2}\\left(\\kappa_{2}\\right)$")
legend("topright", inset = 0.02, legend = my.expressions, lty = MyLines, lwd = 2, col = MyColors, cex = 0.8)
dev.off()

tikz('PlotInfiniteTimeRuinProbabilityExample52 (Version 2) - Participant 3.tex', standAlone = TRUE, width = 4, height = 4, packages = c("\\usepackage{tikz}", "\\usepackage[active,tightpage,psfixbb]{preview}", "\\PreviewEnvironment{pgfpicture}", "\\setlength\\PreviewBorder{0pt}", "\\usepackage{amssymb}", "\\usepackage{scalerel}"))
par(mgp = c(2.5, 1, 0), mar = c(3.5, 3.5, 1, 1) + 0.1)
MyColors <- pal_jco()(3)
MyLines <-seq(from = 1, to = 3, by = 1)
plot(initial_capitals, ruin_probs_3s, type = "l", lwd = 2, lty = MyLines[1], col = MyColors[1], xaxs = "i", yaxs = "i", xlim = c(min(initial_capitals), max(initial_capitals)), ylim = c(0, 1), xlab = "$\\kappa_{3}$", ylab = "Ruin Probability")
lines(initial_capitals, ruin_probsMP_3s, lwd = 2, lty = MyLines[2], col = MyColors[2])
lines(initial_capitals, ruin_probsALT_3s, lwd = 2, lty = MyLines[3], col = MyColors[3])
my.expressions <- c("$\\psi_{3}\\left(\\kappa_{3}\\right)$", "$\\psi^{\\scaleto{{\\fontfamily{qcr}\\selectfont pool - MP}}{3pt}}_{3}\\left(\\kappa_{3}\\right)$", "$\\psi^{\\scaleto{{\\fontfamily{qcr}\\selectfont pool - ALT}}{3pt}}_{3}\\left(\\kappa_{3}\\right)$")
legend("topright", inset = 0.02, legend = my.expressions, lty = MyLines, lwd = 2, col = MyColors, cex = 0.8)
dev.off()

# tikz('PlotInfiniteTimeRuinProbabilityExample52 (Version 2) - All Participants MP.tex', standAlone = TRUE, width = 4, height = 4, packages = c("\\usepackage{tikz}", "\\usepackage[active,tightpage,psfixbb]{preview}", "\\PreviewEnvironment{pgfpicture}", "\\setlength\\PreviewBorder{0pt}", "\\usepackage{amssymb}", "\\usepackage{scalerel}"))
# my.expressions <- c("$\\psi_{{\\scaleto{1}{3pt}}}^{{\\scaleto{ {\\fontfamily{qcr}\\selectfont POOL}}{3pt}}}\\left(\\kappa\\right)$", "$\\psi_{{\\scaleto{2}{3pt}}}^{{\\scaleto{ {\\fontfamily{qcr}\\selectfont POOL}}{3pt}}}\\left(\\kappa\\right)$", "$\\psi_{{\\scaleto{3}{3pt}}}^{{\\scaleto{ {\\fontfamily{qcr}\\selectfont POOL}}{3pt}}}\\left(\\kappa\\right)$", "$\\psi_{{\\scaleto{1}{3pt}}}\\left(\\kappa\\right)$", "$\\psi_{{\\scaleto{2}{3pt}}}\\left(\\kappa\\right)$", "$\\psi_{{\\scaleto{3}{3pt}}}\\left(\\kappa\\right)$")
# par(mgp = c(2.5, 1, 0), mar = c(3.5, 3.5, 1, 1) + 0.1)
# plot(initial_capitals, ruin_probsMP_1s, type = "l", lwd = 2, lty = 1, col = "blue", xaxs = "i", yaxs = "i", xlim = c(min(initial_capitals), 5), ylim = c(0, 1), xlab = "$\\kappa$", ylab = "Ruin Probability")
# grid(col = "gray70", lty = "dotted")
# lines(initial_capitals, ruin_probsMP_2s, lwd = 2, lty = 1, col = "green")
# lines(initial_capitals, ruin_probsMP_3s, lwd = 2, lty = 1, col = "red")
# lines(initial_capitals, ruin_probs_1s, lwd = 2, lty = 2, col = "orange")
# lines(initial_capitals, ruin_probs_2s, lwd = 2, lty = 2, col = "gray")
# lines(initial_capitals, ruin_probs_3s, lwd = 2, lty = 2, col = "brown")
# legend("topright", inset = 0.02, legend = my.expressions, col = c("blue", "green", "red", "orange", "gray", "brown"), lwd = c(2, 2, 2, 2, 2, 2), lty = c(1, 1, 1, 2, 3, 4), cex = 0.8)
# dev.off()
# 
# tikz('PlotInfiniteTimeRuinProbabilityExample52 (Version 2) - All Participants ALT.tex', standAlone = TRUE, width = 4, height = 4, packages = c("\\usepackage{tikz}", "\\usepackage[active,tightpage,psfixbb]{preview}", "\\PreviewEnvironment{pgfpicture}", "\\setlength\\PreviewBorder{0pt}", "\\usepackage{amssymb}", "\\usepackage{scalerel}"))
# my.expressions <- c("$\\psi_{{\\scaleto{1}{3pt}}}^{{\\scaleto{ {\\fontfamily{qcr}\\selectfont POOL}}{3pt}}}\\left(\\kappa\\right)$", "$\\psi_{{\\scaleto{2}{3pt}}}^{{\\scaleto{ {\\fontfamily{qcr}\\selectfont POOL}}{3pt}}}\\left(\\kappa\\right)$", "$\\psi_{{\\scaleto{3}{3pt}}}^{{\\scaleto{ {\\fontfamily{qcr}\\selectfont POOL}}{3pt}}}\\left(\\kappa\\right)$", "$\\psi_{{\\scaleto{1}{3pt}}}\\left(\\kappa\\right)$", "$\\psi_{{\\scaleto{2}{3pt}}}\\left(\\kappa\\right)$", "$\\psi_{{\\scaleto{3}{3pt}}}\\left(\\kappa\\right)$")
# par(mgp = c(2.5, 1, 0), mar = c(3.5, 3.5, 1, 1) + 0.1)
# plot(initial_capitals, ruin_probsALT_1s, type = "l", lwd = 2, lty = 1, col = "blue", xaxs = "i", yaxs = "i", xlim = c(min(initial_capitals), 5), ylim = c(0, 1), xlab = "$\\kappa$", ylab = "Ruin Probability")
# grid(col = "gray70", lty = "dotted")
# lines(initial_capitals, ruin_probsALT_2s, lwd = 2, lty = 1, col = "green")
# lines(initial_capitals, ruin_probsALT_3s, lwd = 2, lty = 1, col = "red")
# lines(initial_capitals, ruin_probs_1s, lwd = 2, lty = 2, col = "orange")
# lines(initial_capitals, ruin_probs_2s, lwd = 2, lty = 2, col = "gray")
# lines(initial_capitals, ruin_probs_3s, lwd = 2, lty = 2, col = "brown")
# legend("topright", inset = 0.02, legend = my.expressions, col = c("blue", "green", "red", "orange", "gray", "brown"), lwd = c(2, 2, 2, 2, 2, 2), lty = c(1, 1, 1, 2, 3, 4), cex = 0.8)
# dev.off()
