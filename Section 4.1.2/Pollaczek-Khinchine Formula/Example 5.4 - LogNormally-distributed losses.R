# Community-Based Insurance in the Compound Poisson Surplus Model
# R Code to calculate the ruin probability when each participant's claim sizes are LogNormally distributed
# Authors: Denuit, M., Flores-Contró, J. M. and Robert, C. Y.

rm(list = ls())

######################################################################################################################################
######################################################################################################################################

# We load the required packages.

library(actuar)
library(tikzDevice)
options(tikzMetricPackages = c("\\usepackage[utf8]{inputenc}", "\\usepackage[T1]{fontenc}", "\\usetikzlibrary{calc}", "\\usepackage{amssymb}"))

# We define the input parameters.

initial_capitals <- seq(0, 20, by = 0.1)
eta <- 2/5
lambda_j <- c(1, 1, 1)
lambda <- sum(lambda_j)
mu <- c(-3.2238, -0.1711, 0.2876)
sigma <- sqrt(c(4.615193, 0.342225, 0.2357102))
EYj <- exp(mu + sigma^2 / 2)

# Parameters for participant 1.

M_val <- (lambda_j[1] * EYj[1]) / sum(lambda_j * EYj)
MP_1s <- rep(M_val, 3)
ALT_1s <- c(0.08, 0.1, 0.1786667)
meanMP_pool_S1 <- sum(lambda_j * MP_1s * EYj)
meanALT_pool_S1 <- sum(lambda_j * ALT_1s * EYj)
mean_S1 <- lambda_j[1] * EYj[1]
c_1 <- (1 + eta) * lambda_j[1] * EYj[1]
rhoMP_pool_1 <- meanMP_pool_S1 / c_1
rhoALT_pool_1 <- meanALT_pool_S1 / c_1
rho_1 <- mean_S1 / c_1

# Parameters for participant 2.

M_val <- (lambda_j[2] * EYj[2]) / sum(lambda_j * EYj)
MP_2s <- rep(M_val, 3)
ALT_2s <- c(0.02, 0.85, 0.0946667)
meanMP_pool_S2 <- sum(lambda_j * MP_2s * EYj)
meanALT_pool_S2 <- sum(lambda_j * ALT_2s * EYj)
mean_S2 <- lambda_j[2] * EYj[2]
c_2 <- (1 + eta) * lambda_j[2] * EYj[2]
rhoMP_pool_2 <- meanMP_pool_S2 / c_2
rhoALT_pool_2 <- meanALT_pool_S2 / c_2
rho_2 <- mean_S2 / c_2

# Parameters for participant 3.

M_val <- (lambda_j[3] * EYj[3]) / sum(lambda_j * EYj)
MP_3s <- rep(M_val, 3)
ALT_3s <- c(0.9, 0.05, 0.72666666)
meanMP_pool_S3 <- sum(lambda_j * MP_3s * EYj)
meanALT_pool_S3 <- sum(lambda_j * ALT_3s * EYj)
mean_S3 <- lambda_j[3] * EYj[3]
c_3 <- (1 + eta) * lambda_j[3] * EYj[3]
rhoMP_pool_3 <- meanMP_pool_S3 / c_3
rhoALT_pool_3 <- meanALT_pool_S3 / c_3
rho_3 <- mean_S3 / c_3

# We double check that actuarial fairness holds.

cat("Fairness check participant 1:\n")
cat("  Pool mean:", sum(lambda_j * MP_1s * EYj), " | Individual mean:", lambda_j[1] * EYj[1], "\n")
cat("  Pool mean:", sum(lambda_j * ALT_1s * EYj), " | Individual mean:", lambda_j[1] * EYj[1], "\n")

cat("Fairness check participant 2:\n")
cat("  Pool mean:", sum(lambda_j * MP_2s * EYj), " | Individual mean:", lambda_j[2] * EYj[2], "\n")
cat("  Pool mean:", sum(lambda_j * ALT_2s * EYj), " | Individual mean:", lambda_j[2] * EYj[2], "\n")

cat("Fairness check participant 3:\n")
cat("  Pool mean:", sum(lambda_j * MP_3s * EYj), " | Individual mean:", lambda_j[3] * EYj[3], "\n")
cat("  Pool mean:", sum(lambda_j * ALT_3s * EYj), " | Individual mean:", lambda_j[3] * EYj[3], "\n")

# Now, we discretise the CDF for each of the individual insurance accounts.

step  <- 0.1
upper <- 1000

Rdist_1 <- discretize(plnorm(x, mu[1], sigma[1]), from = 0, to = upper, step = step, method = "unbiased", lev = levlnorm(x, mu[1], sigma[1]))
Rdist_2 <- discretize(plnorm(x, mu[2], sigma[2]), from = 0, to = upper, step = step, method = "unbiased", lev = levlnorm(x, mu[2], sigma[2]))
Rdist_3 <- discretize(plnorm(x, mu[3], sigma[3]), from = 0, to = upper, step = step, method = "unbiased", lev = levlnorm(x, mu[3], sigma[3]))

# We double-check that they sum up 1.

cat("\nDiscretisation mass check (should be close to 1):\n")
cat("  Rdist_1:", sum(Rdist_1), "\n")
cat("  Rdist_2:", sum(Rdist_2), "\n")
cat("  Rdist_3:", sum(Rdist_3), "\n")

# We double-check that the mean of the discretiside distribution actually matches the theoretical mean.

grid <- seq(0, upper, by = step)
cat("\nDiscretisation mean check (exp(mu[i] + sigma[i]^2/2)):\n") 
cat("  Rdist_1:", sum(grid * Rdist_1), "\n")
cat("  Rdist_2:", sum(grid * Rdist_2), "\n")
cat("  Rdist_3:", sum(grid * Rdist_3), "\n")
cat("  Theoretical means:", EYj, "\n")

# We create a function that estimates a discretised version of the integrated tail (equilibrium / size-biased) distribution.

integrated_tail <- function(Rdist, EY, step) {
  Fdist <- cumsum(Rdist)
  Sdist <- 1 - Fdist
  fe <- c(0, Sdist[-length(Sdist)]) * step / EY
  fe <- fe / sum(fe)
  return(fe)
}

# We estimate the integrated tail distribution for each participant.

integrated_tail_1 <- integrated_tail(Rdist_1, EYj[1], step)
integrated_tail_2 <- integrated_tail(Rdist_2, EYj[2], step)
integrated_tail_3 <- integrated_tail(Rdist_3, EYj[3], step)

# Using Panjer, we estimate the aggregate claim amount distribution.

aggregate_claim_amount_distribution_1 <- aggregateDist("recursive", model.freq = "geometric", model.sev = integrated_tail_1, x.scale = step, prob = 1 - rho_1)
aggregate_claim_amount_distribution_2 <- aggregateDist("recursive", model.freq = "geometric", model.sev = integrated_tail_2, x.scale = step, prob = 1 - rho_2)
aggregate_claim_amount_distribution_3 <- aggregateDist("recursive", model.freq = "geometric", model.sev = integrated_tail_3, x.scale = step, prob = 1 - rho_3)

# Then, we calculate the ruin probability for each participant.

ruin_probability_1 <- 1 - aggregate_claim_amount_distribution_1(initial_capitals)
ruin_probability_2 <- 1 - aggregate_claim_amount_distribution_2(initial_capitals)
ruin_probability_3 <- 1 - aggregate_claim_amount_distribution_3(initial_capitals)

# Now, we focus on the pooled ruin probabilities.

# First, we define a vector that contains the weights for our mixture (must sum to 1).

w <- lambda_j / lambda

# Next, we create a function that estimates the CDF for our mixture of lognormals.

pmixlnorm <- function(x, w, M_is, mu, sigma) {
  rowSums(sapply(seq_along(w), function(j) w[j] * plnorm(x/M_is[j], mu[j], sigma[j])))
}

# Then, we also create a function that calculates the limited expected value (E[min(X, x)]) for our mixture of lognormals.

levmixlnorm <- function(x, w, M_is, mu, sigma) {
  rowSums(sapply(seq_along(w), function(j) w[j] * M_is[j] * levlnorm(x/M_is[j], mu[j], sigma[j])))
}

# We discretise the CDF for each of the pooled insurance accounts.

RdistMP_pool_1 <- discretize(pmixlnorm(x, w, MP_1s, mu, sigma), from = 0, to = upper, step = step, method = "unbiased", lev = levmixlnorm(x, w, MP_1s, mu, sigma))
RdistMP_pool_2 <- discretize(pmixlnorm(x, w, MP_2s, mu, sigma), from = 0, to = upper, step = step, method = "unbiased", lev = levmixlnorm(x, w, MP_2s, mu, sigma))
RdistMP_pool_3 <- discretize(pmixlnorm(x, w, MP_3s, mu, sigma), from = 0, to = upper, step = step, method = "unbiased", lev = levmixlnorm(x, w, MP_3s, mu, sigma))

RdistALT_pool_1 <- discretize(pmixlnorm(x, w, ALT_1s, mu, sigma), from = 0, to = upper, step = step, method = "unbiased", lev = levmixlnorm(x, w, ALT_1s, mu, sigma))
RdistALT_pool_2 <- discretize(pmixlnorm(x, w, ALT_2s, mu, sigma), from = 0, to = upper, step = step, method = "unbiased", lev = levmixlnorm(x, w, ALT_2s, mu, sigma))
RdistALT_pool_3 <- discretize(pmixlnorm(x, w, ALT_3s, mu, sigma), from = 0, to = upper, step = step, method = "unbiased", lev = levmixlnorm(x, w, ALT_3s, mu, sigma))

# We double-check that the mean of the discretiside distribution actually matches the theoretical mean.

cat("\nDiscretisation mean check:\n") 
cat("  RdistMP_pool_1:", sum(grid * RdistMP_pool_1), "\n")
cat("  RdistMP_pool_2:", sum(grid * RdistMP_pool_2), "\n")
cat("  RdistMP_pool_3:", sum(grid * RdistMP_pool_3), "\n")
cat("  Theoretical mean for participant 1:", sum(lambda_j/lambda * exp(mu + sigma^2 / 2) * MP_1s), "\n")
cat("  Theoretical mean for participant 2:", sum(lambda_j/lambda * exp(mu + sigma^2 / 2) * MP_2s), "\n")
cat("  Theoretical mean for participant 3:", sum(lambda_j/lambda * exp(mu + sigma^2 / 2) * MP_3s), "\n")

cat("  RdistALT_pool_1:", sum(grid * RdistALT_pool_1), "\n")
cat("  RdistALT_pool_2:", sum(grid * RdistALT_pool_2), "\n")
cat("  RdistALT_pool_3:", sum(grid * RdistALT_pool_3), "\n")
cat("  Theoretical mean for participant 1:", sum(lambda_j/lambda * exp(mu + sigma^2 / 2) * ALT_1s), "\n")
cat("  Theoretical mean for participant 2:", sum(lambda_j/lambda * exp(mu + sigma^2 / 2) * ALT_2s), "\n")
cat("  Theoretical mean for participant 3:", sum(lambda_j/lambda * exp(mu + sigma^2 / 2) * ALT_3s), "\n")

# We estimate the integrated tail distribution for each participant.

integrated_tailMP_pool_1 <- integrated_tail(RdistMP_pool_1, sum(lambda_j/lambda * exp(mu + sigma^2 / 2) * MP_1s), step)
integrated_tailMP_pool_2 <- integrated_tail(RdistMP_pool_2, sum(lambda_j/lambda * exp(mu + sigma^2 / 2) * MP_2s), step)
integrated_tailMP_pool_3 <- integrated_tail(RdistMP_pool_3, sum(lambda_j/lambda * exp(mu + sigma^2 / 2) * MP_3s), step)

integrated_tailALT_pool_1 <- integrated_tail(RdistALT_pool_1, sum(lambda_j/lambda * exp(mu + sigma^2 / 2) * ALT_1s), step)
integrated_tailALT_pool_2 <- integrated_tail(RdistALT_pool_2, sum(lambda_j/lambda * exp(mu + sigma^2 / 2) * ALT_2s), step)
integrated_tailALT_pool_3 <- integrated_tail(RdistALT_pool_3, sum(lambda_j/lambda * exp(mu + sigma^2 / 2) * ALT_3s), step)

# Pool rho for each participant (all should be < 1 for ruin probability to make sense).

cat("rhoMP_pool_1:", rhoMP_pool_1, "\n")
cat("rhoMP_pool_2:", rhoMP_pool_2, "\n")
cat("rhoMP_pool_3:", rhoMP_pool_3, "\n")

cat("rhoALT_pool_1:", rhoALT_pool_1, "\n")
cat("rhoALT_pool_2:", rhoALT_pool_2, "\n")
cat("rhoALT_pool_3:", rhoALT_pool_3, "\n")

# Using Panjer, we estimate the aggregate claim amount distribution.

aggregate_claim_amount_distributionMP_pool_1 <- aggregateDist("recursive", model.freq = "geometric", model.sev = integrated_tailMP_pool_1, x.scale = step, prob = 1 - rhoMP_pool_1)
aggregate_claim_amount_distributionMP_pool_2 <- aggregateDist("recursive", model.freq = "geometric", model.sev = integrated_tailMP_pool_2, x.scale = step, prob = 1 - rhoMP_pool_2)
aggregate_claim_amount_distributionMP_pool_3 <- aggregateDist("recursive", model.freq = "geometric", model.sev = integrated_tailMP_pool_3, x.scale = step, prob = 1 - rhoMP_pool_3)

aggregate_claim_amount_distributionALT_pool_1 <- aggregateDist("recursive", model.freq = "geometric", model.sev = integrated_tailALT_pool_1, x.scale = step, prob = 1 - rhoALT_pool_1)
aggregate_claim_amount_distributionALT_pool_2 <- aggregateDist("recursive", model.freq = "geometric", model.sev = integrated_tailALT_pool_2, x.scale = step, prob = 1 - rhoALT_pool_2)
aggregate_claim_amount_distributionALT_pool_3 <- aggregateDist("recursive", model.freq = "geometric", model.sev = integrated_tailALT_pool_3, x.scale = step, prob = 1 - rhoALT_pool_3)

# Then, we calculate the ruin probability for each participant.

ruin_probabilityMP_pool_1 <- 1 - aggregate_claim_amount_distributionMP_pool_1(initial_capitals)
ruin_probabilityMP_pool_2 <- 1 - aggregate_claim_amount_distributionMP_pool_2(initial_capitals)
ruin_probabilityMP_pool_3 <- 1 - aggregate_claim_amount_distributionMP_pool_3(initial_capitals)

ruin_probabilityALT_pool_1 <- 1 - aggregate_claim_amount_distributionALT_pool_1(initial_capitals)
ruin_probabilityALT_pool_2 <- 1 - aggregate_claim_amount_distributionALT_pool_2(initial_capitals)
ruin_probabilityALT_pool_3 <- 1 - aggregate_claim_amount_distributionALT_pool_3(initial_capitals)

# Sanity check: psi(0) should equal rho_pool_i.

cat("psi_pool_1(0):", ruin_probabilityMP_pool_1[1], "| rho_pool_1:", rhoMP_pool_1, "\n")
cat("psi_pool_2(0):", ruin_probabilityMP_pool_2[1], "| rho_pool_2:", rhoMP_pool_2, "\n")
cat("psi_pool_3(0):", ruin_probabilityMP_pool_3[1], "| rho_pool_3:", rhoMP_pool_3, "\n")

cat("psi_pool_1(0):", ruin_probabilityALT_pool_1[1], "| rho_pool_1:", rhoALT_pool_1, "\n")
cat("psi_pool_2(0):", ruin_probabilityALT_pool_2[1], "| rho_pool_2:", rhoALT_pool_2, "\n")
cat("psi_pool_3(0):", ruin_probabilityALT_pool_3[1], "| rho_pool_3:", rhoALT_pool_3, "\n")

# Lastly, we plot our results.

file <- '/Users/josemiguelflorescontro/Documents/UCLouvain/Projects/Linear Risk Sharing/R/Graphs/Latex Codes to Generate Graphs'
setwd(file)

tikz('PlotInfiniteTimeRuinProbabilityExample54 - Participant 1.tex', standAlone = TRUE, width = 4, height = 4, packages = c("\\usepackage{tikz}", "\\usepackage[active,tightpage,psfixbb]{preview}", "\\PreviewEnvironment{pgfpicture}", "\\setlength\\PreviewBorder{0pt}", "\\usepackage{amssymb}", "\\usepackage{scalerel}"))
par(mgp = c(2.5, 1, 0), mar = c(3.5, 3.5, 1, 1) + 0.1)
MyColors <- pal_jco()(3)
MyLines <-seq(from = 1, to = 3, by = 1)
plot(initial_capitals, ruin_probability_1, type = "l", lwd = 2, lty = MyLines[1], col = MyColors[1], xaxs = "i", yaxs = "i", xlim = c(min(initial_capitals), max(initial_capitals)), ylim = c(0, 1), xlab = "$\\kappa_{1}$", ylab = "Ruin Probability")
lines(initial_capitals, ruin_probabilityMP_pool_1, lwd = 2, lty = MyLines[2], col = MyColors[2])
lines(initial_capitals, ruin_probabilityALT_pool_1, lwd = 2, lty = MyLines[3], col = MyColors[3])
my.expressions <- c("$\\psi_{1}\\left(\\kappa_{1}\\right)$", "$\\psi^{\\scaleto{{\\fontfamily{qcr}\\selectfont pool - MP}}{3pt}}_{1}\\left(\\kappa_{1}\\right)$", "$\\psi^{\\scaleto{{\\fontfamily{qcr}\\selectfont pool - ALT}}{3pt}}_{1}\\left(\\kappa_{1}\\right)$")
legend("topright", inset = 0.02, legend = my.expressions, lty = MyLines, lwd = 2, col = MyColors, cex = 0.8)
dev.off()

tikz('PlotInfiniteTimeRuinProbabilityExample54 - Participant 2.tex', standAlone = TRUE, width = 4, height = 4, packages = c("\\usepackage{tikz}", "\\usepackage[active,tightpage,psfixbb]{preview}", "\\PreviewEnvironment{pgfpicture}", "\\setlength\\PreviewBorder{0pt}", "\\usepackage{amssymb}", "\\usepackage{scalerel}"))
par(mgp = c(2.5, 1, 0), mar = c(3.5, 3.5, 1, 1) + 0.1)
MyColors <- pal_jco()(3)
MyLines <-seq(from = 1, to = 3, by = 1)
plot(initial_capitals, ruin_probability_2, type = "l", lwd = 2, lty = MyLines[1], col = MyColors[1], xaxs = "i", yaxs = "i", xlim = c(min(initial_capitals), max(initial_capitals)), ylim = c(0, 1), xlab = "$\\kappa_{2}$", ylab = "Ruin Probability")
lines(initial_capitals, ruin_probabilityMP_pool_2, lwd = 2, lty = MyLines[2], col = MyColors[2])
lines(initial_capitals, ruin_probabilityALT_pool_2, lwd = 2, lty = MyLines[3], col = MyColors[3])
my.expressions <- c("$\\psi_{2}\\left(\\kappa_{2}\\right)$", "$\\psi^{\\scaleto{{\\fontfamily{qcr}\\selectfont pool - MP}}{3pt}}_{2}\\left(\\kappa_{2}\\right)$", "$\\psi^{\\scaleto{{\\fontfamily{qcr}\\selectfont pool - ALT}}{3pt}}_{2}\\left(\\kappa_{2}\\right)$")
legend("topright", inset = 0.02, legend = my.expressions, lty = MyLines, lwd = 2, col = MyColors, cex = 0.8)
dev.off()

tikz('PlotInfiniteTimeRuinProbabilityExample54 - Participant 3.tex', standAlone = TRUE, width = 4, height = 4, packages = c("\\usepackage{tikz}", "\\usepackage[active,tightpage,psfixbb]{preview}", "\\PreviewEnvironment{pgfpicture}", "\\setlength\\PreviewBorder{0pt}", "\\usepackage{amssymb}", "\\usepackage{scalerel}"))
par(mgp = c(2.5, 1, 0), mar = c(3.5, 3.5, 1, 1) + 0.1)
MyColors <- pal_jco()(3)
MyLines <-seq(from = 1, to = 3, by = 1)
plot(initial_capitals, ruin_probability_3, type = "l", lwd = 2, lty = MyLines[1], col = MyColors[1], xaxs = "i", yaxs = "i", xlim = c(min(initial_capitals), max(initial_capitals)), ylim = c(0, 1), xlab = "$\\kappa_{3}$", ylab = "Ruin Probability")
lines(initial_capitals, ruin_probabilityMP_pool_3, lwd = 2, lty = MyLines[2], col = MyColors[2])
lines(initial_capitals, ruin_probabilityALT_pool_3, lwd = 2, lty = MyLines[3], col = MyColors[3])
my.expressions <- c("$\\psi_{3}\\left(\\kappa_{3}\\right)$", "$\\psi^{\\scaleto{{\\fontfamily{qcr}\\selectfont pool - MP}}{3pt}}_{3}\\left(\\kappa_{3}\\right)$", "$\\psi^{\\scaleto{{\\fontfamily{qcr}\\selectfont pool - ALT}}{3pt}}_{3}\\left(\\kappa_{3}\\right)$")
legend("topright", inset = 0.02, legend = my.expressions, lty = MyLines, lwd = 2, col = MyColors, cex = 0.8)
dev.off()

# tikz('PlotInfiniteTimeRuinProbabilityExample54 - All Participants MP.tex', standAlone = TRUE, width = 4, height = 4, packages = c("\\usepackage{tikz}", "\\usepackage[active,tightpage,psfixbb]{preview}", "\\PreviewEnvironment{pgfpicture}", "\\setlength\\PreviewBorder{0pt}", "\\usepackage{amssymb}", "\\usepackage{scalerel}"))
# my.expressions <- c("$\\psi_{{\\scaleto{1}{3pt}}}^{{\\scaleto{ {\\fontfamily{qcr}\\selectfont POOL}}{3pt}}}\\left(\\kappa\\right)$", "$\\psi_{{\\scaleto{2}{3pt}}}^{{\\scaleto{ {\\fontfamily{qcr}\\selectfont POOL}}{3pt}}}\\left(\\kappa\\right)$", "$\\psi_{{\\scaleto{3}{3pt}}}^{{\\scaleto{ {\\fontfamily{qcr}\\selectfont POOL}}{3pt}}}\\left(\\kappa\\right)$", "$\\psi_{{\\scaleto{1}{3pt}}}\\left(\\kappa\\right)$", "$\\psi_{{\\scaleto{2}{3pt}}}\\left(\\kappa\\right)$", "$\\psi_{{\\scaleto{3}{3pt}}}\\left(\\kappa\\right)$")
# par(mgp = c(2.5, 1, 0), mar = c(3.5, 3.5, 1, 1) + 0.1)
# plot(initial_capitals, ruin_probabilityMP_pool_1, type = "l", lwd = 2, lty = 1, col = "blue", xaxs = "i", yaxs = "i", xlim = c(min(initial_capitals), 5), ylim = c(0, 1), xlab = "$\\kappa$", ylab = "Ruin Probability")
# grid(col = "gray70", lty = "dotted")
# lines(initial_capitals, ruin_probabilityMP_pool_2, lwd = 2, lty = 1, col = "green")
# lines(initial_capitals, ruin_probabilityMP_pool_3, lwd = 2, lty = 1, col = "red")
# lines(initial_capitals, ruin_probability_1, lwd = 2, lty = 2, col = "orange")
# lines(initial_capitals, ruin_probability_2, lwd = 2, lty = 2, col = "gray")
# lines(initial_capitals, ruin_probability_3, lwd = 2, lty = 2, col = "brown")
# legend("topright", inset = 0.02, legend = my.expressions, col = c("blue", "green", "red", "orange", "gray", "brown"), lwd = c(2, 2, 2, 2, 2, 2), lty = c(1, 1, 1, 2, 3, 4), cex = 0.8)
# dev.off()
# 
# tikz('PlotInfiniteTimeRuinProbabilityExample54 - All Participants ALT.tex', standAlone = TRUE, width = 4, height = 4, packages = c("\\usepackage{tikz}", "\\usepackage[active,tightpage,psfixbb]{preview}", "\\PreviewEnvironment{pgfpicture}", "\\setlength\\PreviewBorder{0pt}", "\\usepackage{amssymb}", "\\usepackage{scalerel}"))
# my.expressions <- c("$\\psi_{{\\scaleto{1}{3pt}}}^{{\\scaleto{ {\\fontfamily{qcr}\\selectfont POOL}}{3pt}}}\\left(\\kappa\\right)$", "$\\psi_{{\\scaleto{2}{3pt}}}^{{\\scaleto{ {\\fontfamily{qcr}\\selectfont POOL}}{3pt}}}\\left(\\kappa\\right)$", "$\\psi_{{\\scaleto{3}{3pt}}}^{{\\scaleto{ {\\fontfamily{qcr}\\selectfont POOL}}{3pt}}}\\left(\\kappa\\right)$", "$\\psi_{{\\scaleto{1}{3pt}}}\\left(\\kappa\\right)$", "$\\psi_{{\\scaleto{2}{3pt}}}\\left(\\kappa\\right)$", "$\\psi_{{\\scaleto{3}{3pt}}}\\left(\\kappa\\right)$")
# par(mgp = c(2.5, 1, 0), mar = c(3.5, 3.5, 1, 1) + 0.1)
# plot(initial_capitals, ruin_probabilityALT_pool_1, type = "l", lwd = 2, lty = 1, col = "blue", xaxs = "i", yaxs = "i", xlim = c(min(initial_capitals), 5), ylim = c(0, 1), xlab = "$\\kappa$", ylab = "Ruin Probability")
# grid(col = "gray70", lty = "dotted")
# lines(initial_capitals, ruin_probabilityALT_pool_2, lwd = 2, lty = 1, col = "green")
# lines(initial_capitals, ruin_probabilityALT_pool_3, lwd = 2, lty = 1, col = "red")
# lines(initial_capitals, ruin_probability_1, lwd = 2, lty = 2, col = "orange")
# lines(initial_capitals, ruin_probability_2, lwd = 2, lty = 2, col = "gray")
# lines(initial_capitals, ruin_probability_3, lwd = 2, lty = 2, col = "brown")
# legend("topright", inset = 0.02, legend = my.expressions, col = c("blue", "green", "red", "orange", "gray", "brown"), lwd = c(2, 2, 2, 2, 2, 2), lty = c(1, 1, 1, 2, 3, 4), cex = 0.8)
# dev.off()