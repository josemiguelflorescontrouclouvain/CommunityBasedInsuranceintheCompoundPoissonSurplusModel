# Linear Risk Sharing in Community-Based Insurance: Ruin Reduction in the Compound Poisson Model 
# R Code to calculate Pareto-efficient transfer ratios for a pool with two participants when claim sizes are Exponentially distributed - Section 5.2
# Authors: Denuit, M., Flores-Contró, J. M. and Robert, C. Y.

rm(list = ls())

######################################################################################################################################
######################################################################################################################################

# We load the required packages.

library(ggsci)
library(tikzDevice)
options(tikzMetricPackages = c("\\usepackage[utf8]{inputenc}", "\\usepackage[T1]{fontenc}", "\\usetikzlibrary{calc}", "\\usepackage{amssymb}"))

######################################################################################################################################
######################################################################################################################################

# First, we define the general parameters. Note that, we store these parameters in a list called "parameters"

parameters <- list(
  eta      = 2/5, # Loading factor
  alpha_1  = 0.5, # Claim-severity rate for participant 1
  alpha_2  = 0.7, # Claim-severity rate for participant 2
  lambda_1 = 0.5, # Claim-frequency intensity rate for participant 1
  lambda_2 = 0.6, # Claim-frequency intensity rate for participant 2
  w_1      = 0.5, # Weight assigned to ruin probability of participant 1
  w_2      = 0.5 # Weight assigned to ruin probability of participant 2
)

# Second, we incorporate to the "parameters" list, a few other (calculated) parameters

parameters$b_1     <- 1/parameters$alpha_1 # Expected-severity loss for participant 1
parameters$b_2     <- 1/parameters$alpha_2 # Expected-severity loss for participant 2
parameters$lambda  <- parameters$lambda_1 + parameters$lambda_2 # Claim-frequency intensity rate for the pool
parameters$p       <- parameters$lambda_1/parameters$lambda   # In this case, we have that: 1 - p = lambda_2/lambda

######################################################################################################################################
######################################################################################################################################

# The function "compute_k_and_EZ" is in charge of calculating k_11, k_12, k_21, k_22, E[Z_1] and E[Z_2]
# The input for this function is:
# parameters: a list containing all the required parameters
# Moreover, note that this will be a function of the parameter "a" that we are looking for
# The output of this function is a list containing the calculated parameters: k_11, k_12, k_21, k_22, E[Z_1] and E[Z_2]

compute_k_and_EZ <- function(a, parameters) {
  # First, all the parameters are read
  alpha_1  <- parameters$alpha_1
  alpha_2  <- parameters$alpha_2
  lambda_1 <- parameters$lambda_1
  lambda_2 <- parameters$lambda_2
  lambda <- parameters$lambda
  b_1 <- parameters$b_1
  b_2 <- parameters$b_2
  # Then, one can proceed to calculate: k_11, k_12, k_21, k_22, E[Z_1] and E[Z_2]
  # Participant 1 - Derivation of parameters 
  k_11  <- a/alpha_1
  k_12  <- (lambda_1 * (1 - a)) / (lambda_2 * alpha_1)
  E_Z_1 <- (a * lambda_1)/(alpha_1 * lambda) + (lambda_1 * (1 - a))/(alpha_1 * lambda)
  # Participant 2 - Derivation of parameters 
  k_21  <- (1 - a) / alpha_1
  k_22  <- (lambda_2 * b_2 - lambda_1 * b_1 * (1 - a))/(lambda_2 * b_2 * alpha_2)
  E_Z_2 <- ((1 - a) * lambda_1) / (alpha_1 * lambda) + (lambda_2 * (lambda_2 * b_2 - lambda_1 * b_1 * (1 - a)))/(lambda * lambda_2 * b_2 * alpha_2)
  # Finally, we store our results for k_11, k_12, k_21, k_22, E[Z_1] and E[Z_2] in a list
  list(k_11 = k_11, k_12 = k_12, E_Z_1 = E_Z_1, k_21 = k_21, k_22 = k_22, E_Z_2 = E_Z_2)
}

######################################################################################################################################
######################################################################################################################################

# The function "compute_coefficients_a" is in charge of calculating the coefficients a_2, a_1 and a_0
# The input for this function is:
# eta: loading factor
# E_Z: expected value either of Z_1 or Z_2
# k_1: parameter either k_11 or k_21
# k_2: parameter either k_12 or k_22
# The output of this function is a list containing the calculated parameters: a_2, a_1 and a_0

compute_coefficients_a <- function(eta, E_Z, k_1, k_2) {
  # First, we compute a_2
  a_2 <- (1 + eta) * E_Z * k_1 * k_2
  # Second, a_1 is computed
  a_1 <- -((1 + eta) * E_Z * (k_1 + k_2) - k_1 * k_2)
  # Lastly, the function computes a_0 
  a_0 <- eta * E_Z
  # Finally, we store our results for a_2, a_1 and a_0 in a list
  list(a_2 = a_2, a_1 = a_1, a_0 = a_0)
}

######################################################################################################################################
######################################################################################################################################

# The function "compute_r_C" is in charge of calculating r_11, r_12, C_11 and C_12 (or equivalently, r_21, r_22, C_21 and C_22)
# The input for this function is:
# eta: loading factor
# E_Z: expected value either of Z_1 or Z_2
# k_1: parameter either k_11 or k_21
# k_2: parameter either k_12 or k_22
# The output of this function is a list containing the calculated values for: r_11, r_12, C_11 and C_12 (or equivalently, if working with participant 2, r_21, r_22, C_21 and C_22)

compute_r_C <- function(eta, E_Z, k_1, k_2) {
  # First, we call our function "compute_coefficients_a", which will allow us to compute a_2, a_1 and a_0
  compute_coefficients_a  <- compute_coefficients_a(eta, E_Z, k_1, k_2)
  # Then, we save the variables a_2, a_1 and a_0
  a_2 <- compute_coefficients_a$a_2
  a_1 <- compute_coefficients_a$a_1
  a_0 <- compute_coefficients_a$a_0
  # Sanity check: verify that the values are real numbers (i.e., we check that the discriminant is not negative)
  disc <- a_1^2 - 4 * a_2 * a_0
  if (is.na(disc) || disc < 0 || a_2 == 0) {
    return(list(r_1 = NA_real_, r_2 = NA_real_, C_1 = NA_real_, C_2 = NA_real_))
  }
  sq <- sqrt(disc)
  # One can then calculate the roots r_11 and r_12 (or equivalently, if working with participant 2, r_21 and r_22)
  r_1 <- (-a_1 + sq)/(2 * a_2)
  r_2 <- (-a_1 - sq)/(2 * a_2)
  # Similarly, one proceeds to calculate the coefficients C_11 and C_12 (or equivalently, if working with participant 2, C_21 and C_22)
  C_1 <- -(eta * (E_Z - k_1 * k_2 * r_1))/(r_1 * (1 + eta) * (2 * a_2 * r_1 + a_1))
  C_2 <- -(eta * (E_Z - k_1 * k_2 * r_2))/(r_2 * (1 + eta) * (2 * a_2 * r_2 + a_1))
  # As a last step, one can store the results for r_11, r_12, C_11 and C_12 in a list (or equivalently, if working with participant 2, the results for r_21, r_22, C_21 and C_22)
  list(r_1 = r_1, r_2 = r_2, C_1 = C_1, C_2 = C_2)
}

######################################################################################################################################
######################################################################################################################################

# The function "compute_ruin_probability" is in charge of computing the ruin probability (clearly, this function can be used to
# compute the ruin probability of participant 1 or participant 2)
# The input for this function is:
# eta: loading factor
# E_Z: expected value (Z_1 or Z_2)
# k_1: parameter (k_11 or k_21)
# k_2: parameter (k_12 or k_22)
# kappa : initial reserve

compute_ruin_probability <- function(eta, E_Z, k_1, k_2, kappa) {
  # First, we call our function "compute_r_C" to compute the required values: r_11, r_12, C_11 and C_12 (or equivalently, if working with participant 2, r_21, r_22, C_21 and C_22)
  rc <- compute_r_C(eta, E_Z, k_1, k_2)
  # Then, we are now in position to compute the ruin probability
  rc$C_1 * exp(-rc$r_1 * kappa) + rc$C_2 * exp(-rc$r_2 * kappa)
}

######################################################################################################################################
######################################################################################################################################

# The function "weighted_ruin_probability" is in charge of estimating the weighted ruin probability for both participants
# The input for this function is:
# kappa_1: initial reserve for participant 1 
# kappa_2: initial reserve for participant 2
# parameters: a list containing all the required parameters
# It is worth noting that this will be a function of the transfer ratio "a". The optimal transfer ratio, "a★", is the value
# of "a" that that minimizes the weighted sum of ruin probabilities

weighted_ruin_probability <- function(a, kappa_1, kappa_2, parameters) {
  # First, we call our function "compute_k_and_EZ" to compute the required parameters: k_11, k_12, k_21, k_22, E[Z_1] and E[Z_2]
  kz <- compute_k_and_EZ(a, parameters)
  # Then, we are now in position to compute the ruin probability for each of the participants
  psi_1 <- compute_ruin_probability(parameters$eta, kz$E_Z_1, kz$k_11, kz$k_12, kappa_1)
  psi_2 <- compute_ruin_probability(parameters$eta, kz$E_Z_2, kz$k_21, kz$k_22, kappa_2)
  # Sanity check: we verify that our ruin probabilities are in the range (0,1)
  if (is.na(psi_1) || is.na(psi_2) || psi_1 < 0 || psi_1 > 1 || psi_2 < 0 || psi_2 > 1) {
    return(NA_real_)
  }
  # Lastly, we compute the weighted sum of the ruin probabilities for participant 1 and participant 2
  parameters$w_1 * psi_1 + parameters$w_2 * psi_2
}

# The function "stand_alone_ruin_probability" is in charge of estimating the ruin probability for each participant
# The input for this function is:
# kappa_1: initial reserve for participant 1 
# kappa_2: initial reserve for participant 2
# parameters: a list containing all the required parameters
# It is worth noting that this will be a function of the transfer ratio "a". The optimal transfer ratio, "a★", is the value
# of "a" that that minimizes the weighted sum of ruin probabilities

pooled_ruin_probability <- function(a, kappa_1, kappa_2, parameters) {
  # First, we call our function "compute_k_and_EZ" to compute the required parameters: k_11, k_12, k_21, k_22, E[Z_1] and E[Z_2]
  kz <- compute_k_and_EZ(a, parameters)
  # Then, we are now in position to compute the ruin probability for each of the participants
  psi_1 <- compute_ruin_probability(parameters$eta, kz$E_Z_1, kz$k_11, kz$k_12, kappa_1)
  psi_2 <- compute_ruin_probability(parameters$eta, kz$E_Z_2, kz$k_21, kz$k_22, kappa_2)
  # Sanity check: we verify that our ruin probabilities are in the range (0,1)
  if (is.na(psi_1) || is.na(psi_2) || psi_1 < 0 || psi_1 > 1 || psi_2 < 0 || psi_2 > 1) {
    return(NA_real_)
  }
  # As a last step, one can store the ruin probabilities in a list
  list(psi_1 = psi_1, psi_2 = psi_2)
}

######################################################################################################################################
######################################################################################################################################

# PARETO OPTIMALITY FOR RUIN PROBABILITIES: OPTIMIZATION PROBLEM

# First, we define some other general parameters

kappa_1 <- 3 # Initial reserve for participant 1   
kappa_2 <- 3 # Initial reserve for participant 2

# Moreover, we define our objective functions for each scenario

# SCENARIO 1: FULL WEIGHT TO RUIN PROBABILITY OF PARTICIPANT 1
objective_function_ws_1 <- function(a) {
  parameters$w_1 <- 1
  parameters$w_2 <- 0
  weighted_ruin_probability(a, kappa_1, kappa_2, parameters)
  }

# SCENARIO 2: FULL WEIGHT TO RUIN PROBABILITY OF PARTICIPANT 2
objective_function_ws_2 <- function(a) {
  parameters$w_1 <- 0
  parameters$w_2 <- 1
  weighted_ruin_probability(a, kappa_1, kappa_2, parameters)
}

# SCENARIO 3: WEIGHT BOTH PARTICIPANTS EQUALLY AT 1/2
objective_function_ws_3 <- function(a) {
  parameters$w_1 <- 0.5
  parameters$w_2 <- 0.5
  weighted_ruin_probability(a, kappa_1, kappa_2, parameters)
}

# SCENARIO 4: WEIGHT FOR PARTICIPANTS 1 AT 0.7 AND WEIGHT FOR PARTICIPANTS 2 AT 0.3
objective_function_ws_4 <- function(a) {
  parameters$w_1 <- 0.7
  parameters$w_2 <- 0.3
  weighted_ruin_probability(a, kappa_1, kappa_2, parameters)
}

# SCENARIO 5: WEIGHT FOR PARTICIPANTS 1 AT 0.3 AND WEIGHT FOR PARTICIPANTS 2 AT 0.7
objective_function_ws_5 <- function(a) {
  parameters$w_1 <- 0.3
  parameters$w_2 <- 0.7
  weighted_ruin_probability(a, kappa_1, kappa_2, parameters)
}

# Now, we derive the feasible range for the transfer ratio "a". Recall that, the range is obtained
# from considering the definition of the transfer ratio (i.e., knowing that a ∈ [0, 1]), actuarial fairness,
# full allocation and capacity constraints. All these together lead to the following feasible range for 
# the transfer ratio "a".
a_min <- max(1 - (parameters$lambda_2 * parameters$b_2)/(parameters$lambda_1 * parameters$b_1), 1 - parameters$b_2/parameters$b_1, 1 - parameters$lambda_2/parameters$lambda_1)
a_min <- max(0, a_min) + 0.0000001
a_max <- 1 - 0.0000001
cat(sprintf("\nFeasible range for a given current parameters: [%.4f, %.4f)\n", a_min, a_max)) # Print the feasible range for the transfer ratio "a"

# Now, we can find the optimal transfer ratio "a★" for each of our scenarios

# SCENARIO 1: FULL WEIGHT TO RUIN PROBABILITY OF PARTICIPANT 1
optimal_a_ws_1 <- optimize(objective_function_ws_1, interval = c(a_min, a_max))
cat("Optimal for Weights-Scenario 1 a* =", optimal_a_ws_1$minimum, " -> weighted ruin probability =", optimal_a_ws_1$objective, "\n")

# SCENARIO 2: FULL WEIGHT TO RUIN PROBABILITY OF PARTICIPANT 2
optimal_a_ws_2 <- optimize(objective_function_ws_2, interval = c(a_min, a_max))
cat("Optimal for Weights-Scenario 2 a* =", optimal_a_ws_2$minimum, " -> weighted ruin probability =", optimal_a_ws_2$objective, "\n")

# SCENARIO 3: WEIGHT BOTH PARTICIPANTS EQUALLY AT 1/2
optimal_a_ws_3 <- optimize(objective_function_ws_3, interval = c(a_min, a_max))
cat("Optimal for Weights-Scenario 3 a* =", optimal_a_ws_3$minimum, " -> weighted ruin probability =", optimal_a_ws_3$objective, "\n")

# SCENARIO 4: WEIGHT FOR PARTICIPANTS 1 AT 0.7 AND WEIGHT FOR PARTICIPANTS 2 AT 0.3
optimal_a_ws_4 <- optimize(objective_function_ws_4, interval = c(a_min, a_max))
cat("Optimal for Weights-Scenario 4 a* =", optimal_a_ws_4$minimum, " -> weighted ruin probability =", optimal_a_ws_4$objective, "\n")

# SCENARIO 5: WEIGHT FOR PARTICIPANTS 1 AT 0.3 AND WEIGHT FOR PARTICIPANTS 2 AT 0.7
optimal_a_ws_5 <- optimize(objective_function_ws_5, interval = c(a_min, a_max))
cat("Optimal for Weights-Scenario 5 a* =", optimal_a_ws_5$minimum, " -> weighted ruin probability =", optimal_a_ws_5$objective, "\n")

# Set working directory
file <- '/Users/jose/Library/CloudStorage/OneDrive-UCL/Documents/Postdoc/Linear Risk Sharing Project/R/Graphs/Latex Codes to Generate Graphs'
setwd(file)

# Set feasible range for transfer ratio "a" in plots
a_grid <- seq(a_min, 0.999, length.out = 100)

# We compute the weighted sum of the ruin probabilities throughout the feasible range for the transfer ratio "a".
# This is done for each of our scenarios.

# SCENARIO 1: FULL WEIGHT TO RUIN PROBABILITY OF PARTICIPANT 1
psi_weighted_grid_ws_1 <- sapply(a_grid, objective_function_ws_1)

# SCENARIO 2: FULL WEIGHT TO RUIN PROBABILITY OF PARTICIPANT 2
psi_weighted_grid_ws_2 <- sapply(a_grid, objective_function_ws_2)

# SCENARIO 3: WEIGHT BOTH PARTICIPANTS EQUALLY AT 1/2
psi_weighted_grid_ws_3 <- sapply(a_grid, objective_function_ws_3)

# SCENARIO 4: WEIGHT FOR PARTICIPANTS 1 AT 0.7 AND WEIGHT FOR PARTICIPANTS 2 AT 0.3
psi_weighted_grid_ws_4 <- sapply(a_grid, objective_function_ws_4)

# SCENARIO 5: WEIGHT FOR PARTICIPANTS 1 AT 0.3 AND WEIGHT FOR PARTICIPANTS 2 AT 0.7
psi_weighted_grid_ws_5 <- sapply(a_grid, objective_function_ws_5)

# For Figure 1 (a), we need to compute the pooled ruin probability for each participant 
# throughout the feasible range of the transfer ratio "a". We proceed to compute these
# ruin probabilities

# Participant 1
psi_1 <- function(a) pooled_ruin_probability(a, kappa_1, kappa_2, parameters)$psi_1
psi_1_non_optimal <- sapply(a_grid, psi_1)

# Participant 2
psi_2 <- function(a) pooled_ruin_probability(a, kappa_1, kappa_2, parameters)$psi_2
psi_2_non_optimal <- sapply(a_grid, psi_2)

# In Figure 1 (a), we also want to show the ruin probability values (for both participants)
# when considering the optimal transfer ratio "a★". Clearly, we want to do this for each
# scenario.

# SCENARIO 1: FULL WEIGHT TO RUIN PROBABILITY OF PARTICIPANT 1
# Participant 1
psi_1_optimal_ws_1 <- sapply(optimal_a_ws_1$minimum, psi_1)
# Participant 2
psi_2_optimal_ws_1 <- sapply(optimal_a_ws_1$minimum, psi_2)

# SCENARIO 2: FULL WEIGHT TO RUIN PROBABILITY OF PARTICIPANT 2
# Participant 1
psi_1_optimal_ws_2 <- sapply(optimal_a_ws_2$minimum, psi_1)
# Participant 2
psi_2_optimal_ws_2 <- sapply(optimal_a_ws_2$minimum, psi_2)

# SCENARIO 3: WEIGHT BOTH PARTICIPANTS EQUALLY AT 1/2
# Participant 1
psi_1_optimal_ws_3 <- sapply(optimal_a_ws_3$minimum, psi_1)
# Participant 2
psi_2_optimal_ws_3 <- sapply(optimal_a_ws_3$minimum, psi_2)

# SCENARIO 4: WEIGHT FOR PARTICIPANTS 1 AT 0.7 AND WEIGHT FOR PARTICIPANTS 2 AT 0.3
# Participant 1
psi_1_optimal_ws_4 <- sapply(optimal_a_ws_4$minimum, psi_1)
# Participant 2
psi_2_optimal_ws_4 <- sapply(optimal_a_ws_4$minimum, psi_2)

# SCENARIO 5: WEIGHT FOR PARTICIPANTS 1 AT 0.3 AND WEIGHT FOR PARTICIPANTS 2 AT 0.7
# Participant 1
psi_1_optimal_ws_5 <- sapply(optimal_a_ws_5$minimum, psi_1)
# Participant 2
psi_2_optimal_ws_5 <- sapply(optimal_a_ws_5$minimum, psi_2)

# Now, we compute the allocation matrix under each of our scenarios.
# Furthermore, we verify that the transfer ratios in such allocation matrix comply with
# actuarial fairness, full allocation and capacity constraints. We do these for each of
# our scenarios.

# SCENARIO 1: FULL WEIGHT TO RUIN PROBABILITY OF PARTICIPANT 1
A_ws_1 <- matrix(c(optimal_a_ws_1$minimum, (parameters$lambda_1 * parameters$b_1)/(parameters$lambda_2 * parameters$b_2) * (1 - optimal_a_ws_1$minimum),
              1 - optimal_a_ws_1$minimum, 1 - (parameters$lambda_1 * parameters$b_1)/(parameters$lambda_2 * parameters$b_2) * (1 - optimal_a_ws_1$minimum)),
            nrow = 2, byrow = TRUE)

cat("The Matrix A for Weights-Scenario 1 is:\n")
print(A_ws_1)
cat("The Sum of the Columns of Matrix A for Weights-Scenario 1 is =", colSums(A_ws_1), "\n")
cat("Check of Full Allocattion for Column 1 of Matrix A in Weights-Scenario 1 is: ", A_ws_1[1, 1] + A_ws_1[2,1]  == 1, "\n")
cat("Check of Full Allocattion for Column 2 of Matrix A in Weights-Scenario 1 is: ", A_ws_1[1, 2] + A_ws_1[2,2]  == 1, "\n")
cat("Check of Actuarial Fairness - Participant 1 - for Matrix A in Weights-Scenario 1 is: ", round(A_ws_1[1, 1] * (parameters$lambda_1 * parameters$b_1) + A_ws_1[1, 2] * (parameters$lambda_2 * parameters$b_2), 3) == round(parameters$lambda_1 * parameters$b_1, 3), "\n")
cat("Check of Actuarial Fairness - Participant 2 - for Matrix A in Weights-Scenario 1 is: ", round(A_ws_1[2, 1] * (parameters$lambda_1 * parameters$b_1) + A_ws_1[2, 2] * (parameters$lambda_2 * parameters$b_2), 3) == round(parameters$lambda_2 * parameters$b_2, 3), "\n")
cat("Check of Capacity Constraint - Participant 1 - for Matrix A in Weights-Scenario 1 is: ", A_ws_1[1, 1] * parameters$b_1 <= parameters$b_1, "\n")
cat("Check of Capacity Constraint - Participant 1 - for Matrix A in Weights-Scenario 1 is: ", A_ws_1[1, 2] * parameters$b_2 <= parameters$b_1, "\n")
cat("Check of Capacity Constraint - Participant 2 - for Matrix A in Weights-Scenario 1 is: ", A_ws_1[2, 1] * parameters$b_1 <= parameters$b_2, "\n")
cat("Check of Capacity Constraint - Participant 2 - for Matrix A in Weights-Scenario 1 is: ", A_ws_1[2, 2] * parameters$b_2 <= parameters$b_2, "\n")

# SCENARIO 2: FULL WEIGHT TO RUIN PROBABILITY OF PARTICIPANT 2
A_ws_2 <- matrix(c(optimal_a_ws_2$minimum, (parameters$lambda_1 * parameters$b_1)/(parameters$lambda_2 * parameters$b_2) * (1 - optimal_a_ws_2$minimum),
                   1 - optimal_a_ws_2$minimum, 1 - (parameters$lambda_1 * parameters$b_1)/(parameters$lambda_2 * parameters$b_2) * (1 - optimal_a_ws_2$minimum)),
                 nrow = 2, byrow = TRUE)

cat("The Matrix A for Weights-Scenario 2 is:\n")
print(A_ws_2)
cat("The Sum of the Columns of Matrix A for Weights-Scenario 2 is =", colSums(A_ws_2), "\n")
cat("Check of Full Allocattion for Column 1 of Matrix A in Weights-Scenario 2 is: ", A_ws_2[1, 1] + A_ws_2[2,1]  == 1, "\n")
cat("Check of Full Allocattion for Column 2 of Matrix A in Weights-Scenario 2 is: ", A_ws_2[1, 2] + A_ws_2[2,2]  == 1, "\n")
cat("Check of Actuarial Fairness - Participant 1 - for Matrix A in Weights-Scenario 2 is: ", round(A_ws_2[1, 1] * (parameters$lambda_1 * parameters$b_1) + A_ws_2[1, 2] * (parameters$lambda_2 * parameters$b_2), 3) == round(parameters$lambda_1 * parameters$b_1, 3), "\n")
cat("Check of Actuarial Fairness - Participant 2 - for Matrix A in Weights-Scenario 2 is: ", round(A_ws_2[2, 1] * (parameters$lambda_1 * parameters$b_1) + A_ws_2[2, 2] * (parameters$lambda_2 * parameters$b_2), 3) == round(parameters$lambda_2 * parameters$b_2, 3), "\n")
cat("Check of Capacity Constraint - Participant 1 - for Matrix A in Weights-Scenario 1 is: ", A_ws_2[1, 1] * parameters$b_1 <= parameters$b_1, "\n")
cat("Check of Capacity Constraint - Participant 1 - for Matrix A in Weights-Scenario 1 is: ", A_ws_2[1, 2] * parameters$b_2 <= parameters$b_1, "\n")
cat("Check of Capacity Constraint - Participant 2 - for Matrix A in Weights-Scenario 1 is: ", A_ws_2[2, 1] * parameters$b_1 <= parameters$b_2, "\n")
cat("Check of Capacity Constraint - Participant 2 - for Matrix A in Weights-Scenario 1 is: ", A_ws_2[2, 2] * parameters$b_2 <= parameters$b_2, "\n")

# SCENARIO 3: WEIGHT BOTH PARTICIPANTS EQUALLY AT 1/2
A_ws_3 <- matrix(c(optimal_a_ws_3$minimum, (parameters$lambda_1 * parameters$b_1)/(parameters$lambda_2 * parameters$b_2) * (1 - optimal_a_ws_3$minimum),
                   1 - optimal_a_ws_3$minimum, 1 - (parameters$lambda_1 * parameters$b_1)/(parameters$lambda_2 * parameters$b_2) * (1 - optimal_a_ws_3$minimum)),
                 nrow = 2, byrow = TRUE)

cat("The Matrix A for Weights-Scenario 3 is:\n")
print(A_ws_3)
cat("The Sum of the Columns of Matrix A for Weights-Scenario 3 is =", colSums(A_ws_3), "\n")
cat("Check of Full Allocattion for Column 1 of Matrix A in Weights-Scenario 3 is: ", A_ws_3[1, 1] + A_ws_3[2,1]  == 1, "\n")
cat("Check of Full Allocattion for Column 2 of Matrix A in Weights-Scenario 3 is: ", A_ws_3[1, 2] + A_ws_3[2,2]  == 1, "\n")
cat("Check of Actuarial Fairness - Participant 1 - for Matrix A in Weights-Scenario 3 is: ", round(A_ws_3[1, 1] * (parameters$lambda_1 * parameters$b_1) + A_ws_3[1, 2] * (parameters$lambda_2 * parameters$b_2), 3) == round(parameters$lambda_1 * parameters$b_1, 3), "\n")
cat("Check of Actuarial Fairness - Participant 2 - for Matrix A in Weights-Scenario 3 is: ", round(A_ws_3[2, 1] * (parameters$lambda_1 * parameters$b_1) + A_ws_3[2, 2] * (parameters$lambda_2 * parameters$b_2), 3) == round(parameters$lambda_2 * parameters$b_2, 3), "\n")
cat("Check of Capacity Constraint - Participant 1 - for Matrix A in Weights-Scenario 1 is: ", A_ws_3[1, 1] * parameters$b_1 <= parameters$b_1, "\n")
cat("Check of Capacity Constraint - Participant 1 - for Matrix A in Weights-Scenario 1 is: ", A_ws_3[1, 2] * parameters$b_2 <= parameters$b_1, "\n")
cat("Check of Capacity Constraint - Participant 2 - for Matrix A in Weights-Scenario 1 is: ", A_ws_3[2, 1] * parameters$b_1 <= parameters$b_2, "\n")
cat("Check of Capacity Constraint - Participant 2 - for Matrix A in Weights-Scenario 1 is: ", A_ws_3[2, 2] * parameters$b_2 <= parameters$b_2, "\n")

# SCENARIO 4: WEIGHT FOR PARTICIPANTS 1 AT 0.7 AND WEIGHT FOR PARTICIPANTS 2 AT 0.3
A_ws_4 <- matrix(c(optimal_a_ws_4$minimum, (parameters$lambda_1 * parameters$b_1)/(parameters$lambda_2 * parameters$b_2) * (1 - optimal_a_ws_4$minimum),
                   1 - optimal_a_ws_4$minimum, 1 - (parameters$lambda_1 * parameters$b_1)/(parameters$lambda_2 * parameters$b_2) * (1 - optimal_a_ws_4$minimum)),
                 nrow = 2, byrow = TRUE)

cat("The Matrix A for Weights-Scenario 4 is:\n")
print(A_ws_4)
cat("The Sum of the Columns of Matrix A for Weights-Scenario 4 is =", colSums(A_ws_4), "\n")
cat("Check of Full Allocattion for Column 1 of Matrix A in Weights-Scenario 4 is: ", A_ws_4[1, 1] + A_ws_4[2,1]  == 1, "\n")
cat("Check of Full Allocattion for Column 2 of Matrix A in Weights-Scenario 4 is: ", A_ws_4[1, 2] + A_ws_4[2,2]  == 1, "\n")
cat("Check of Actuarial Fairness - Participant 1 - for Matrix A in Weights-Scenario 4 is: ", round(A_ws_4[1, 1] * (parameters$lambda_1 * parameters$b_1) + A_ws_4[1, 2] * (parameters$lambda_2 * parameters$b_2), 3) == round(parameters$lambda_1 * parameters$b_1, 3), "\n")
cat("Check of Actuarial Fairness - Participant 2 - for Matrix A in Weights-Scenario 4 is: ", round(A_ws_4[2, 1] * (parameters$lambda_1 * parameters$b_1) + A_ws_4[2, 2] * (parameters$lambda_2 * parameters$b_2), 3) == round(parameters$lambda_2 * parameters$b_2, 3), "\n")
cat("Check of Capacity Constraint - Participant 1 - for Matrix A in Weights-Scenario 1 is: ", A_ws_4[1, 1] * parameters$b_1 <= parameters$b_1, "\n")
cat("Check of Capacity Constraint - Participant 1 - for Matrix A in Weights-Scenario 1 is: ", A_ws_4[1, 2] * parameters$b_2 <= parameters$b_1, "\n")
cat("Check of Capacity Constraint - Participant 2 - for Matrix A in Weights-Scenario 1 is: ", A_ws_4[2, 1] * parameters$b_1 <= parameters$b_2, "\n")
cat("Check of Capacity Constraint - Participant 2 - for Matrix A in Weights-Scenario 1 is: ", A_ws_4[2, 2] * parameters$b_2 <= parameters$b_2, "\n")

# SCENARIO 5: WEIGHT FOR PARTICIPANTS 1 AT 0.3 AND WEIGHT FOR PARTICIPANTS 2 AT 0.7
A_ws_5 <- matrix(c(optimal_a_ws_5$minimum, (parameters$lambda_1 * parameters$b_1)/(parameters$lambda_2 * parameters$b_2) * (1 - optimal_a_ws_5$minimum),
                   1 - optimal_a_ws_5$minimum, 1 - (parameters$lambda_1 * parameters$b_1)/(parameters$lambda_2 * parameters$b_2) * (1 - optimal_a_ws_5$minimum)),
                 nrow = 2, byrow = TRUE)

cat("The Matrix A for Weights-Scenario 5 is:\n")
print(A_ws_5)
cat("The Sum of the Columns of Matrix A for Weights-Scenario 5 is =", colSums(A_ws_5), "\n")
cat("Check of Full Allocattion for Column 1 of Matrix A in Weights-Scenario 5 is: ", A_ws_5[1, 1] + A_ws_5[2,1]  == 1, "\n")
cat("Check of Full Allocattion for Column 2 of Matrix A in Weights-Scenario 5 is: ", A_ws_5[1, 2] + A_ws_5[2,2]  == 1, "\n")
cat("Check of Actuarial Fairness - Participant 1 - for Matrix A in Weights-Scenario 5 is: ", round(A_ws_5[1, 1] * (parameters$lambda_1 * parameters$b_1) + A_ws_5[1, 2] * (parameters$lambda_2 * parameters$b_2), 3) == round(parameters$lambda_1 * parameters$b_1, 3), "\n")
cat("Check of Actuarialfairness - Participant 2 - for Matrix A in Weights-Scenario 5 is: ", round(A_ws_5[2, 1] * (parameters$lambda_1 * parameters$b_1) + A_ws_5[2, 2] * (parameters$lambda_2 * parameters$b_2), 3) == round(parameters$lambda_2 * parameters$b_2, 3), "\n")
cat("Check of Capacity Constraint - Participant 1 - for Matrix A in Weights-Scenario 1 is: ", A_ws_5[1, 1] * parameters$b_1 <= parameters$b_1, "\n")
cat("Check of Capacity Constraint - Participant 1 - for Matrix A in Weights-Scenario 1 is: ", A_ws_5[1, 2] * parameters$b_2 <= parameters$b_1, "\n")
cat("Check of Capacity Constraint - Participant 2 - for Matrix A in Weights-Scenario 1 is: ", A_ws_5[2, 1] * parameters$b_1 <= parameters$b_2, "\n")
cat("Check of Capacity Constraint - Participant 2 - for Matrix A in Weights-Scenario 1 is: ", A_ws_5[2, 2] * parameters$b_2 <= parameters$b_2, "\n")

# Lastly, we plot our results.

tikz('PlotParetoEfficientCandidatesSection52 - Figure 1.tex', standAlone = TRUE, width = 4, height = 4, packages = c("\\usepackage{tikz}", "\\usepackage[active,tightpage,psfixbb]{preview}", "\\PreviewEnvironment{pgfpicture}", "\\setlength\\PreviewBorder{0pt}", "\\usepackage{amssymb}", "\\usepackage{scalerel}", "\\usepackage{amsmath}"))
par(mgp = c(2.5, 1, 0), mar = c(3.5, 4.5, 1, 4.5) + 0.1)
my.expressions <- c(sprintf("$w_{1} = %.3f, w_{2} = %.3f, a^{\\star} = %.3f$", 1, 0, optimal_a_ws_1$minimum), sprintf("$w_{1} = %.3f, w_{2} = %.3f, a^{\\star} = %.3f$", 0, 1, optimal_a_ws_2$minimum), sprintf("$w_{1} = %.3f, w_{2} = %.3f, a^{\\star} = %.3f$", 0.5, 0.5, optimal_a_ws_3$minimum), sprintf("$w_{1} = %.3f, w_{2} = %.3f, a^{\\star} = %.3f$", 0.7, 0.3, optimal_a_ws_4$minimum), sprintf("$w_{1} = %.3f, w_{2} = %.3f, a^{\\star} = %.3f$", 0.3, 0.7, optimal_a_ws_5$minimum))
MyColors <- pal_jco()(5)
cols <- hcl.colors(100, "Oslo")
MyLines <-seq(from = 1, to = 4, by = 1)
plot(psi_1_non_optimal, psi_2_non_optimal, type = "n", lwd = 2, lty = MyLines[1], col = MyColors[1], xaxs = "i", yaxs = "i", xlim = c(0, 0.6), ylim = c(0, 0.4), xlab = "$\\psi_{{\\scaleto{1}{3pt}}}^{{\\scaleto{ {\\fontfamily{qcr}\\selectfont \\text{POOL}}}{3pt}}}\\left(\\kappa_{1};a_{1,1},a_{1,2}\\right)$", ylab = "$\\psi_{{\\scaleto{2}{3pt}}}^{{\\scaleto{ {\\fontfamily{qcr}\\selectfont \\text{POOL}}}{3pt}}}\\left(\\kappa_{2};a_{2,1},a_{2,2}\\right)$", xaxt ="n", yaxt ="n")
x_ticks <- seq(0, 0.6, by = 0.2)
y_ticks <- seq(0, 0.4, by = 0.1)
axis(1, at = x_ticks, labels = sprintf("%.1f", x_ticks))
axis(2, at = y_ticks, labels = sprintf("%.1f", y_ticks))
for (i in 1:(length(a_grid) - 1)) {
  lines(psi_1_non_optimal[i:(i+1)],
        psi_2_non_optimal[i:(i+1)],
        col = hcl.colors(length(a_grid) - 1, "Oslo")[i],
        lwd = 2)
}
points(psi_1_optimal_ws_1, psi_2_optimal_ws_1, pch = 21, bg = MyColors[1], col = MyColors[1], cex = 0.6)
points(psi_1_optimal_ws_2, psi_2_optimal_ws_2, pch = 22, bg = MyColors[2], col = MyColors[2], cex = 0.6)
points(psi_1_optimal_ws_3, psi_2_optimal_ws_3, pch = 23, bg = MyColors[3], col = MyColors[3], cex = 0.6)
points(psi_1_optimal_ws_4, psi_2_optimal_ws_4, pch = 24, bg = MyColors[4], col = MyColors[4], cex = 0.6)
points(psi_1_optimal_ws_5, psi_2_optimal_ws_5, pch = 25, bg = MyColors[5], col = MyColors[5], cex = 0.6)
legend("bottomright", inset = 0.02, legend = my.expressions, col = c(MyColors[1], MyColors[2], MyColors[3], MyColors[4], MyColors[5]), pch = c(21, 22, 23, 24, 25), pt.bg = c(MyColors[1], MyColors[2], MyColors[3], MyColors[4], MyColors[5]), cex = 0.55)
par(fig = c(0.83, 0.86, 0.21, 0.86), new = TRUE, mar = c(0, 0, 0, 0))
image(x = 1, y = a_grid, z = matrix(a_grid, nrow = 1), col = cols, axes = FALSE, xlab = "", ylab = "")
ticks <- seq(a_min, 1, length.out = 5)
axis(4, at = ticks, labels = sprintf("$%.3f$",  round(ticks, 3)))
mtext("$a$", side = 4, line = 2)
dev.off()

tikz('PlotParetoEfficientCandidatesSection52 - Figure 2.tex', standAlone = TRUE, width = 4, height = 4, packages = c("\\usepackage{tikz}", "\\usepackage[active,tightpage,psfixbb]{preview}", "\\PreviewEnvironment{pgfpicture}", "\\setlength\\PreviewBorder{0pt}", "\\usepackage{amssymb}", "\\usepackage{scalerel}", "\\usepackage{amsmath}"))
my.expressions <- c(sprintf("$w_{1} = %.3f, w_{2} = %.3f, a^{\\star} = %.3f$", 1, 0, optimal_a_ws_1$minimum), sprintf("$w_{1} = %.3f, w_{2} = %.3f, a^{\\star} = %.3f$", 0, 1, optimal_a_ws_2$minimum), sprintf("$w_{1} = %.3f, w_{2} = %.3f, a^{\\star} = %.3f$", 0.5, 0.5, optimal_a_ws_3$minimum), sprintf("$w_{1} = %.3f, w_{2} = %.3f, a^{\\star} = %.3f$", 0.7, 0.3, optimal_a_ws_4$minimum), sprintf("$w_{1} = %.3f, w_{2} = %.3f, a^{\\star} = %.3f$", 0.3, 0.7, optimal_a_ws_5$minimum))
par(mgp = c(2.5, 1, 0), mar = c(3.5, 3.5, 1, 1) + 0.1)
MyColors <- pal_jco()(5)
MyLines <-seq(from = 1, to = 5, by = 1)
plot(a_grid, psi_weighted_grid_ws_1, type = "l", lwd = 2, lty = MyLines[1], col = MyColors[1], xaxs = "i", yaxs = "i", xlim = c(max(0, round(a_min - 0.1, 1)), 1), ylim = c(0, 0.6), xlab = "$a$", ylab = "$w_{1}\\psi_{{\\scaleto{1}{3pt}}}^{{\\scaleto{ {\\fontfamily{qcr}\\selectfont \\text{POOL}}}{3pt}}}\\left(\\kappa_{1};a_{1,1},a_{1,2}\\right) + w_{2}\\psi_{{\\scaleto{2}{3pt}}}^{{\\scaleto{ {\\fontfamily{qcr}\\selectfont \\text{POOL}}}{3pt}}}\\left(\\kappa_{2};a_{2,1},a_{2,2}\\right)$")
lines(a_grid, psi_weighted_grid_ws_2, lwd = 2, lty = MyLines[2], col = MyColors[2])
lines(a_grid, psi_weighted_grid_ws_3, lwd = 2, lty = MyLines[3], col = MyColors[3])
lines(a_grid, psi_weighted_grid_ws_4, lwd = 2, lty = MyLines[4], col = MyColors[4])
lines(a_grid, psi_weighted_grid_ws_5, lwd = 2, lty = MyLines[5], col = MyColors[5])
abline(v = a_min, col = "gray70", lty = "dotted")
points(optimal_a_ws_1$minimum, optimal_a_ws_1$objective, bg = MyColors[1], col = MyColors[1], pch = 21, cex = 0.6)
points(optimal_a_ws_2$minimum, optimal_a_ws_2$objective, bg = MyColors[2], col = MyColors[2], pch = 22, cex = 0.6)
points(optimal_a_ws_3$minimum, optimal_a_ws_3$objective, bg = MyColors[3], col = MyColors[3], pch = 23, cex = 0.6)
points(optimal_a_ws_4$minimum, optimal_a_ws_4$objective, bg = MyColors[4], col = MyColors[4], pch = 24, cex = 0.6)
points(optimal_a_ws_5$minimum, optimal_a_ws_5$objective, bg = MyColors[5], col = MyColors[5], pch = 25, cex = 0.6)
legend("bottomright", inset = 0.02, legend = my.expressions, lty = c(MyLines[1], MyLines[2], MyLines[3], MyLines[4], MyLines[5]), col = c(MyColors[1], MyColors[2], MyColors[3], MyColors[4], MyColors[5]), pch = c(21, 22, 23, 24, 25), pt.bg = c(MyColors[1], MyColors[2], MyColors[3], MyColors[4], MyColors[5]), cex = 0.55)
dev.off()