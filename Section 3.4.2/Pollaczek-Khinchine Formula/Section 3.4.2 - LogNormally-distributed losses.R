# Linear Risk Sharing in Community-Based Insurance: Ruin Reduction in the Compound Poisson Model 
# R Code to calculate the ruin probability when each participant's claim sizes are LogNormally distributed - Section 3.4.2
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
lambda_j <- c(2, 1, 3)
lambda <- sum(lambda_j)
mu <- c(0, 0.5, -0.3)
sigma <- c(1, 1, 1)
EYj <- exp(mu + sigma^2 / 2)

# Parameters for participant 1.

MP_val_1s <- (lambda_j[1] * EYj[1]) / sum(lambda_j * EYj)
MP_1s <- rep(MP_val_1s, 3)
ALT_1s <- c(0.400, 0.300, 0.3173894)
meanMP_pool_S1 <- sum(lambda_j * MP_1s * EYj)
meanALT_pool_S1 <- sum(lambda_j * ALT_1s * EYj)
mean_S1 <- lambda_j[1] * EYj[1]
c_1 <- (1 + eta) * lambda_j[1] * EYj[1]
rhoMP_pool_1 <- meanMP_pool_S1 / c_1
rhoALT_pool_1 <- meanALT_pool_S1 / c_1
rho_1 <- mean_S1 / c_1

# Parameters for participant 2.

MP_val_2s <- (lambda_j[2] * EYj[2]) / sum(lambda_j * EYj)
MP_2s <- rep(MP_val_2s, 3)
ALT_2s <- c(0.3741306, 0.300, 0.1826106)
meanMP_pool_S2 <- sum(lambda_j * MP_2s * EYj)
meanALT_pool_S2 <- sum(lambda_j * ALT_2s * EYj)
mean_S2 <- lambda_j[2] * EYj[2]
c_2 <- (1 + eta) * lambda_j[2] * EYj[2]
rhoMP_pool_2 <- meanMP_pool_S2 / c_2
rhoALT_pool_2 <- meanALT_pool_S2 / c_2
rho_2 <- mean_S2 / c_2

# Parameters for participant 3.

MP_val_3s <- (lambda_j[3] * EYj[3]) / sum(lambda_j * EYj)
MP_3s <- rep(MP_val_3s, 3)
ALT_3s <- c(0.2258694, 0.400, 0.500)
meanMP_pool_S3 <- sum(lambda_j * MP_3s * EYj)
meanALT_pool_S3 <- sum(lambda_j * ALT_3s * EYj)
mean_S3 <- lambda_j[3] * EYj[3]
c_3 <- (1 + eta) * lambda_j[3] * EYj[3]
rhoMP_pool_3 <- meanMP_pool_S3 / c_3
rhoALT_pool_3 <- meanALT_pool_S3 / c_3
rho_3 <- mean_S3 / c_3

# Parameters for pooled fund.

mean_S <- mean_S1 + mean_S2 + mean_S3
c <- c_1 + c_2 + c_3
rho_pooled_fund <- mean_S / c

# We double check that actuarial fairness holds.

# Participant 1
cat("Participant 1:\n")
cat("  Fairness check:\n")
cat("    Mean-proportional risk-sharing rule:\n")
cat("      Pool mean:", sum(lambda_j * MP_1s * EYj), " | Individual mean:", lambda_j[1] * EYj[1], "\n")
cat("    Alternative risk-sharing rule:\n")
cat("      Pool mean:", sum(lambda_j * ALT_1s * EYj), " | Individual mean:", lambda_j[1] * EYj[1], "\n")

# Participant 2
cat("Participant 2:\n")
cat("  Fairness check:\n")
cat("    Mean-proportional risk-sharing rule:\n")
cat("      Pool mean:", sum(lambda_j * MP_2s * EYj), " | Individual mean:", lambda_j[2] * EYj[2], "\n")
cat("    Alternative risk-sharing rule:\n")
cat("      Pool mean:", sum(lambda_j * ALT_2s * EYj), " | Individual mean:", lambda_j[2] * EYj[2], "\n")

# Participant 3
cat("Participant 3:\n")
cat("  Fairness check:\n")
cat("    Mean-proportional risk-sharing rule:\n")
cat("      Pool mean:", sum(lambda_j * MP_3s * EYj), " | Individual mean:", lambda_j[3] * EYj[3], "\n")
cat("    Alternative risk-sharing rule:\n")
cat("      Pool mean:", sum(lambda_j * ALT_3s * EYj), " | Individual mean:", lambda_j[3] * EYj[3], "\n")

# We double check that the full allocation property holds.

# Column 1
cat("Column 1:\n")
cat("  Full allocation check:\n")
cat("    Mean-proportional risk-sharing rule:\n")
cat("      Sum (needs to be equal to one):", MP_1s[1] + MP_2s[1] + MP_3s[1], "\n")
cat("    Alternative risk-sharing rule:\n")
cat("      Sum (needs to be equal to one):", ALT_1s[1] + ALT_2s[1] + ALT_3s[1], "\n")

# Column 2
cat("Column 2:\n")
cat("  Full allocation check:\n")
cat("    Mean-proportional risk-sharing rule:\n")
cat("      Sum (needs to be equal to one):", MP_1s[2] + MP_2s[2] + MP_3s[2], "\n")
cat("    Alternative risk-sharing rule:\n")
cat("      Sum (needs to be equal to one):", ALT_1s[2] + ALT_2s[2] + ALT_3s[2], "\n")

# Column 3
cat("Column 3:\n")
cat("  Full allocation check:\n")
cat("    Mean-proportional risk-sharing rule:\n")
cat("      Sum (needs to be equal to one):", MP_1s[3] + MP_2s[3] + MP_3s[3], "\n")
cat("    Alternative risk-sharing rule:\n")
cat("      Sum (needs to be equal to one):", ALT_1s[3] + ALT_2s[3] + ALT_3s[3], "\n")

# We double check that the capacity constraint holds.

# Participant 1
cat("Participant 1:\n")
cat("  Mean-proportional risk-sharing rule:\n")
cat("    Capacity Constraint check for a_11:\n")
cat("      The required inequality is: ", MP_1s[1] * exp(mu[1] + sigma[1]^2/2) <= exp(mu[1] + sigma[1]^2/2), "\n")
cat("    Capacity Constraint check for a_12:\n")
cat("      The required inequality is: ", MP_1s[2] * exp(mu[2] + sigma[2]^2/2) <= exp(mu[1] + sigma[1]^2/2), "\n")
cat("    Capacity Constraint check for a_13:\n")
cat("      The required inequality is: ", MP_1s[3] * exp(mu[3] + sigma[3]^2/2) <= exp(mu[1] + sigma[1]^2/2), "\n")
cat("  Alternative risk-sharing rule:\n")
cat("    Capacity Constraint check for a_11:\n")
cat("      The required inequality is: ", ALT_1s[1] * exp(mu[1] + sigma[1]^2/2) <= exp(mu[1] + sigma[1]^2/2), "\n")
cat("    Capacity Constraint check for a_12:\n")
cat("      The required inequality is: ", ALT_1s[2] * exp(mu[2] + sigma[2]^2/2) <= exp(mu[1] + sigma[1]^2/2), "\n")
cat("    Capacity Constraint check for a_13:\n")
cat("      The required inequality is: ", ALT_1s[3] * exp(mu[3] + sigma[3]^2/2) <= exp(mu[1] + sigma[1]^2/2), "\n")

# Participant 2
cat("Participant 2:\n")
cat("  Mean-proportional risk-sharing rule:\n")
cat("    Capacity Constraint check for a_21:\n")
cat("      The required inequality is: ", MP_2s[1] * exp(mu[1] + sigma[1]^2/2) <= exp(mu[2] + sigma[2]^2/2), "\n")
cat("    Capacity Constraint check for a_22:\n")
cat("      The required inequality is: ", MP_2s[2] * exp(mu[2] + sigma[2]^2/2) <= exp(mu[2] + sigma[2]^2/2), "\n")
cat("    Capacity Constraint check for a_23:\n")
cat("      The required inequality is: ", MP_2s[3] * exp(mu[3] + sigma[3]^2/2) <= exp(mu[2] + sigma[2]^2/2), "\n")
cat("  Alternative risk-sharing rule:\n")
cat("    Capacity Constraint check for a_21:\n")
cat("      The required inequality is: ", ALT_2s[1] * exp(mu[1] + sigma[1]^2/2) <= exp(mu[2] + sigma[2]^2/2), "\n")
cat("    Capacity Constraint check for a_22:\n")
cat("      The required inequality is: ", ALT_2s[2] * exp(mu[2] + sigma[2]^2/2) <= exp(mu[2] + sigma[2]^2/2), "\n")
cat("    Capacity Constraint check for a_23:\n")
cat("      The required inequality is: ", ALT_2s[3] * exp(mu[3] + sigma[3]^2/2) <= exp(mu[2] + sigma[2]^2/2), "\n")

# Participant 3
cat("Participant 3:\n")
cat("  Mean-proportional risk-sharing rule:\n")
cat("    Capacity Constraint check for a_31:\n")
cat("      The required inequality is: ", MP_3s[1] * exp(mu[1] + sigma[1]^2/2) <= exp(mu[3] + sigma[3]^2/2), "\n")
cat("    Capacity Constraint check for a_32:\n")
cat("      The required inequality is: ", MP_3s[2] * exp(mu[2] + sigma[2]^2/2) <= exp(mu[3] + sigma[3]^2/2), "\n")
cat("    Capacity Constraint check for a_33:\n")
cat("      The required inequality is: ", MP_3s[3] * exp(mu[3] + sigma[3]^2/2) <= exp(mu[3] + sigma[3]^2/2), "\n")
cat("  Alternative risk-sharing rule:\n")
cat("    Capacity Constraint check for a_31:\n")
cat("      The required inequality is: ", ALT_3s[1] * exp(mu[1] + sigma[1]^2/2) <= exp(mu[3] + sigma[3]^2/2), "\n")
cat("    Capacity Constraint check for a_32:\n")
cat("      The required inequality is: ", ALT_3s[2] * exp(mu[2] + sigma[2]^2/2) <= exp(mu[3] + sigma[3]^2/2), "\n")
cat("    Capacity Constraint check for a_33:\n")
cat("      The required inequality is: ", ALT_3s[3] * exp(mu[3] + sigma[3]^2/2) <= exp(mu[3] + sigma[3]^2/2), "\n")

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

# We double-check that they sum up 1.

cat("\nDiscretisation mass check (should be close to 1):\n")
cat("  RdistMP_pool_1:", sum(RdistMP_pool_1), "\n")
cat("  RdistMP_pool_2:", sum(RdistMP_pool_2), "\n")
cat("  RdistMP_pool_3:", sum(RdistMP_pool_3), "\n")

cat("  RdistALT_pool_1:", sum(RdistALT_pool_1), "\n")
cat("  RdistALT_pool_2:", sum(RdistALT_pool_2), "\n")
cat("  RdistALT_pool_3:", sum(RdistALT_pool_3), "\n")

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

# Now, we focus on the pooled fund ruin probabilities.

# First, we create a function that estimates the CDF for our non-scaled mixture of lognormals.

pmixlnorm_nonscaled <- function(x, w, mu, sigma) {
  rowSums(sapply(seq_along(w), function(j) w[j] * plnorm(x, mu[j], sigma[j])))
}

# Then, we also create a function that calculates the limited expected value (E[min(X, x)]) for our non-scaled mixture of lognormals.

levmixlnorm <- function(x, w, mu, sigma) {
  rowSums(sapply(seq_along(w), function(j) w[j] * levlnorm(x, mu[j], sigma[j])))
}

# We discretise the CDF.

Rdist_pooled_fund <- discretize(pmixlnorm_nonscaled(x, w, mu, sigma), from = 0, to = upper, step = step, method = "unbiased", lev = levmixlnorm(x, w, mu, sigma))

# We double-check that it sums up 1.

cat("\nDiscretisation mass check (should be close to 1):\n")
cat("  Rdist_pooled_fund:", sum(Rdist_pooled_fund), "\n")

# We double-check that the mean of the discretiside distribution actually matches the theoretical mean.

cat("\nDiscretisation mean check:\n") 
cat("  Rdist_pooled_fund:", sum(grid * Rdist_pooled_fund), "\n")
cat("  Theoretical mean for pooled fund:", sum(lambda_j/lambda * exp(mu + sigma^2 / 2)), "\n")

# We estimate the integrated tail distribution for each participant.

integrated_tail_pooled_fund <- integrated_tail(Rdist_pooled_fund, sum(lambda_j/lambda * exp(mu + sigma^2 / 2)), step)

# Pool rho for each participant (all should be < 1 for ruin probability to make sense).

cat("rho_pooled_fund:", rho_pooled_fund, "\n")

# Using Panjer, we estimate the aggregate claim amount distribution.

aggregate_claim_amount_distribution_pooled_fund <- aggregateDist("recursive", model.freq = "geometric", model.sev = integrated_tail_pooled_fund, x.scale = step, prob = 1 - rho_pooled_fund)

# Then, we calculate the ruin probability for each participant.

ruin_probability_pooled_fund <- 1 - aggregate_claim_amount_distribution_pooled_fund(initial_capitals)

# Sanity check: psi(0) should equal rho_pool_i.

cat("psi_pooled_fund(0):", ruin_probability_pooled_fund[1], "| rho_pooled_fund:", rho_pooled_fund, "\n")

# Lastly, we plot our results.

file <- '/Users/jose/Library/CloudStorage/OneDrive-UCL/Documents/Postdoc/Linear Risk Sharing Project/R/Graphs/Latex Codes to Generate Graphs'
setwd(file)

tikz('PlotInfiniteTimeRuinProbabilitySection342 - Participant 1.tex', standAlone = TRUE, width = 4, height = 4, packages = c("\\usepackage{tikz}", "\\usepackage[active,tightpage,psfixbb]{preview}", "\\PreviewEnvironment{pgfpicture}", "\\setlength\\PreviewBorder{0pt}", "\\usepackage{amssymb}", "\\usepackage{scalerel}", "\\usepackage{amsmath}"))
par(mgp = c(2.5, 1, 0), mar = c(3.5, 3.5, 1, 1) + 0.1)
MyColors <- pal_jco()(4)
MyLines <-seq(from = 1, to = 4, by = 1)
plot(initial_capitals, ruin_probability_1, type = "l", lwd = 2, lty = MyLines[1], col = MyColors[1], xaxs = "i", yaxs = "i", xlim = c(min(initial_capitals), max(initial_capitals)), ylim = c(0, 1), xlab = "$\\kappa, \\kappa_{1}$", ylab = "Ruin Probability")
lines(initial_capitals, ruin_probabilityMP_pool_1, lwd = 2, lty = MyLines[2], col = MyColors[2])
lines(initial_capitals, ruin_probabilityALT_pool_1, lwd = 2, lty = MyLines[3], col = MyColors[3])
lines(initial_capitals, ruin_probability_pooled_fund, lwd = 2, lty = MyLines[4], col = MyColors[4])
my.expressions <- c("$\\psi_{1}\\left(\\kappa_{1}\\right)$", "$\\psi^{\\scaleto{{\\fontfamily{qcr}\\selectfont \\text{pool-MP}}}{3pt}}_{1}\\left(\\kappa_{1}\\right)$", "$\\psi^{\\scaleto{{\\fontfamily{qcr}\\selectfont \\text{pool-ALT}}}{3pt}}_{1}\\left(\\kappa_{1}\\right)$", "$\\psi\\left(\\kappa\\right)$")
legend("topright", inset = 0.02, legend = my.expressions, lty = MyLines, lwd = 2, col = MyColors, cex = 0.8)
dev.off()

tikz('PlotInfiniteTimeRuinProbabilitySection342 - Participant 2.tex', standAlone = TRUE, width = 4, height = 4, packages = c("\\usepackage{tikz}", "\\usepackage[active,tightpage,psfixbb]{preview}", "\\PreviewEnvironment{pgfpicture}", "\\setlength\\PreviewBorder{0pt}", "\\usepackage{amssymb}", "\\usepackage{scalerel}", "\\usepackage{amsmath}"))
par(mgp = c(2.5, 1, 0), mar = c(3.5, 3.5, 1, 1) + 0.1)
MyColors <- pal_jco()(4)
MyLines <-seq(from = 1, to = 4, by = 1)
plot(initial_capitals, ruin_probability_2, type = "l", lwd = 2, lty = MyLines[1], col = MyColors[1], xaxs = "i", yaxs = "i", xlim = c(min(initial_capitals), max(initial_capitals)), ylim = c(0, 1), xlab = "$\\kappa, \\kappa_{2}$", ylab = "Ruin Probability")
lines(initial_capitals, ruin_probabilityMP_pool_2, lwd = 2, lty = MyLines[2], col = MyColors[2])
lines(initial_capitals, ruin_probabilityALT_pool_2, lwd = 2, lty = MyLines[3], col = MyColors[3])
lines(initial_capitals, ruin_probability_pooled_fund, lwd = 2, lty = MyLines[4], col = MyColors[4])
my.expressions <- c("$\\psi_{2}\\left(\\kappa_{2}\\right)$", "$\\psi^{\\scaleto{{\\fontfamily{qcr}\\selectfont \\text{pool-MP}}}{3pt}}_{2}\\left(\\kappa_{2}\\right)$", "$\\psi^{\\scaleto{{\\fontfamily{qcr}\\selectfont \\text{pool-ALT}}}{3pt}}_{2}\\left(\\kappa_{2}\\right)$", "$\\psi\\left(\\kappa\\right)$")
legend("topright", inset = 0.02, legend = my.expressions, lty = MyLines, lwd = 2, col = MyColors, cex = 0.8)
dev.off()

tikz('PlotInfiniteTimeRuinProbabilitySection342 - Participant 3.tex', standAlone = TRUE, width = 4, height = 4, packages = c("\\usepackage{tikz}", "\\usepackage[active,tightpage,psfixbb]{preview}", "\\PreviewEnvironment{pgfpicture}", "\\setlength\\PreviewBorder{0pt}", "\\usepackage{amssymb}", "\\usepackage{scalerel}", "\\usepackage{amsmath}"))
par(mgp = c(2.5, 1, 0), mar = c(3.5, 3.5, 1, 1) + 0.1)
MyColors <- pal_jco()(4)
MyLines <-seq(from = 1, to = 4, by = 1)
plot(initial_capitals, ruin_probability_3, type = "l", lwd = 2, lty = MyLines[1], col = MyColors[1], xaxs = "i", yaxs = "i", xlim = c(min(initial_capitals), max(initial_capitals)), ylim = c(0, 1), xlab = "$\\kappa, \\kappa_{3}$", ylab = "Ruin Probability")
lines(initial_capitals, ruin_probabilityMP_pool_3, lwd = 2, lty = MyLines[2], col = MyColors[2])
lines(initial_capitals, ruin_probabilityALT_pool_3, lwd = 2, lty = MyLines[3], col = MyColors[3])
lines(initial_capitals, ruin_probability_pooled_fund, lwd = 2, lty = MyLines[4], col = MyColors[4])
my.expressions <- c("$\\psi_{3}\\left(\\kappa_{3}\\right)$", "$\\psi^{\\scaleto{{\\fontfamily{qcr}\\selectfont \\text{pool-MP}}}{3pt}}_{3}\\left(\\kappa_{3}\\right)$", "$\\psi^{\\scaleto{{\\fontfamily{qcr}\\selectfont \\text{pool-ALT}}}{3pt}}_{3}\\left(\\kappa_{3}\\right)$", "$\\psi\\left(\\kappa\\right)$")
legend("topright", inset = 0.02, legend = my.expressions, lty = MyLines, lwd = 2, col = MyColors, cex = 0.8)
dev.off()

tikz('PlotInfiniteTimeRuinProbabilitySection342 - All Participants MP.tex', standAlone = TRUE, width = 4, height = 4, packages = c("\\usepackage{tikz}", "\\usepackage[active,tightpage,psfixbb]{preview}", "\\PreviewEnvironment{pgfpicture}", "\\setlength\\PreviewBorder{0pt}", "\\usepackage{amssymb}", "\\usepackage{scalerel}", "\\usepackage{amsmath}"))
my.expressions <- c("$\\psi_{{\\scaleto{1}{3pt}}}^{{\\scaleto{ {\\fontfamily{qcr}\\selectfont \\text{POOL}}}{3pt}}}\\left(\\kappa\\right)$", "$\\psi_{{\\scaleto{2}{3pt}}}^{{\\scaleto{ {\\fontfamily{qcr}\\selectfont \\text{POOL}}}{3pt}}}\\left(\\kappa\\right)$", "$\\psi_{{\\scaleto{3}{3pt}}}^{{\\scaleto{ {\\fontfamily{qcr}\\selectfont \\text{POOL}}}{3pt}}}\\left(\\kappa\\right)$", "$\\psi_{{\\scaleto{1}{3pt}}}\\left(\\kappa\\right)$", "$\\psi_{{\\scaleto{2}{3pt}}}\\left(\\kappa\\right)$", "$\\psi_{{\\scaleto{3}{3pt}}}\\left(\\kappa\\right)$", "$\\psi\\left(\\kappa\\right)$")
par(mgp = c(2.5, 1, 0), mar = c(3.5, 3.5, 1, 1) + 0.1)
plot(initial_capitals, ruin_probabilityMP_pool_1, type = "l", lwd = 2, lty = 1, col = "blue", xaxs = "i", yaxs = "i", xlim = c(min(initial_capitals), 5), ylim = c(0, 1), xlab = "$\\kappa$", ylab = "Ruin Probability")
grid(col = "gray70", lty = "dotted")
lines(initial_capitals, ruin_probabilityMP_pool_2, lwd = 2, lty = 1, col = "green")
lines(initial_capitals, ruin_probabilityMP_pool_3, lwd = 2, lty = 1, col = "red")
lines(initial_capitals, ruin_probability_1, lwd = 2, lty = 2, col = "orange")
lines(initial_capitals, ruin_probability_2, lwd = 2, lty = 3, col = "gray")
lines(initial_capitals, ruin_probability_3, lwd = 2, lty = 4, col = "brown")
lines(initial_capitals, ruin_probability_pooled_fund, lwd = 2, lty = 5, col = "purple")
legend("topright", inset = 0.02, legend = my.expressions, col = c("blue", "green", "red", "orange", "gray", "brown", "purple"), lwd = c(2, 2, 2, 2, 2, 2, 2), lty = c(1, 1, 1, 2, 3, 4, 5), cex = 0.8)
dev.off()

tikz('PlotInfiniteTimeRuinProbabilitySection342 - All Participants ALT.tex', standAlone = TRUE, width = 4, height = 4, packages = c("\\usepackage{tikz}", "\\usepackage[active,tightpage,psfixbb]{preview}", "\\PreviewEnvironment{pgfpicture}", "\\setlength\\PreviewBorder{0pt}", "\\usepackage{amssymb}", "\\usepackage{scalerel}", "\\usepackage{amsmath}"))
my.expressions <- c("$\\psi_{{\\scaleto{1}{3pt}}}^{{\\scaleto{ {\\fontfamily{qcr}\\selectfont \\text{POOL}}}{3pt}}}\\left(\\kappa\\right)$", "$\\psi_{{\\scaleto{2}{3pt}}}^{{\\scaleto{ {\\fontfamily{qcr}\\selectfont \\text{POOL}}}{3pt}}}\\left(\\kappa\\right)$", "$\\psi_{{\\scaleto{3}{3pt}}}^{{\\scaleto{ {\\fontfamily{qcr}\\selectfont \\text{POOL}}}{3pt}}}\\left(\\kappa\\right)$", "$\\psi_{{\\scaleto{1}{3pt}}}\\left(\\kappa\\right)$", "$\\psi_{{\\scaleto{2}{3pt}}}\\left(\\kappa\\right)$", "$\\psi_{{\\scaleto{3}{3pt}}}\\left(\\kappa\\right)$", "$\\psi\\left(\\kappa\\right)$")
par(mgp = c(2.5, 1, 0), mar = c(3.5, 3.5, 1, 1) + 0.1)
plot(initial_capitals, ruin_probabilityALT_pool_1, type = "l", lwd = 2, lty = 1, col = "blue", xaxs = "i", yaxs = "i", xlim = c(min(initial_capitals), 5), ylim = c(0, 1), xlab = "$\\kappa$", ylab = "Ruin Probability")
grid(col = "gray70", lty = "dotted")
lines(initial_capitals, ruin_probabilityALT_pool_2, lwd = 2, lty = 1, col = "green")
lines(initial_capitals, ruin_probabilityALT_pool_3, lwd = 2, lty = 1, col = "red")
lines(initial_capitals, ruin_probability_1, lwd = 2, lty = 2, col = "orange")
lines(initial_capitals, ruin_probability_2, lwd = 2, lty = 3, col = "gray")
lines(initial_capitals, ruin_probability_3, lwd = 2, lty = 4, col = "brown")
lines(initial_capitals, ruin_probability_pooled_fund, lwd = 2, lty = 5, col = "purple")
legend("topright", inset = 0.02, legend = my.expressions, col = c("blue", "green", "red", "orange", "gray", "brown", "purple"), lwd = c(2, 2, 2, 2, 2, 2, 2), lty = c(1, 1, 1, 2, 3, 4, 5), cex = 0.8)
dev.off()