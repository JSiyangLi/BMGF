############
# example 1: marginal likelihood calculation
###########
# model structure: hyperparameters (alpha, Beta) -> parameter lambda -> observation 0
# "observed" data from Poisson likelihood
y = 0

# hyperparameters
alpha = 4
Beta = 5

# moment generating function for gamma prior
mgf <- function(t, alpha, Beta) {
  Beta^alpha / (Beta - t)^alpha
}

# mgf marginalisation
marginal_likelihood <- function(y, t, alpha, Beta) {
  mgf(t, alpha, Beta) *
    prod(1/factorial(y))
}
(m <- marginal_likelihood(y, t = -1, alpha = alpha, Beta = Beta))


# verification via conjugacy + Chib's method
# posterior kernel
Lambda = 1 # choose a point in the sample space of lambda
post_kernel <- function(lambda, alpha, Beta, y) {
  dgamma(lambda, shape = alpha, rate = Beta) * prod(dpois(y, lambda = lambda))
}
(p <- post_kernel(Lambda, alpha, Beta, y) / 
  dgamma(Lambda, shape = alpha + sum(y), rate = Beta + length(y)))
round(p, 15) == round(m, 15) # TRUE

# verification via NegBin
round(dnbinom(y, size=alpha, prob=Beta/(Beta+1)), 15) == round(m, 15)

###############################################
############
# example not included: marginal likelihood calculation
###########
# model structure: hyperparameters (alpha, Beta) -> parameter lambda -> observation 3
# "observed" data from Poisson likelihood
y = 3

# hyperparameters
alpha = 4
Beta = 5

# moment generating function for gamma prior
mgf <- function(t, alpha, Beta) {
  Beta^alpha / (Beta - t)^alpha
}
mgf_first_derivative <- function(t, alpha, Beta) {
  alpha * mgf(t, alpha, Beta) * (Beta - t)^(-1)
}
mgf_second_derivative <- function(t, alpha, Beta) {
  (alpha + 1) * mgf_first_derivative(t, alpha, Beta) * (Beta - t)^(-1)
}
mgf_third_derivative <- function(t, alpha, Beta) {
  (alpha + 2) * mgf_second_derivative(t, alpha, Beta) * (Beta - t)^(-1)
}

# mgf marginalisation
marginal_likelihood <- function(y, t, alpha, Beta) {
  mgf_third_derivative(t, alpha, Beta) *
    prod(1/factorial(y))
}
(m <- marginal_likelihood(y, t = -1, alpha = alpha, Beta = Beta))


# verification via conjugacy + Chib's method
# posterior kernel
Lambda = 1 # choose a point in the sample space of lambda
post_kernel <- function(lambda, alpha, Beta, y) {
  dgamma(lambda, shape = alpha, rate = Beta) * prod(dpois(y, lambda = lambda))
}
(p <- post_kernel(Lambda, alpha, Beta, y) / 
  dgamma(Lambda, shape = alpha + sum(y), rate = Beta + length(y)))
round(p, 15) == round(m, 15) # TRUE

# verification via NegBin
round(dnbinom(y, size=alpha, prob=Beta/(Beta+1)), 15) == round(m, 15)

###########################################
############
# example 2: hierarchical model marginalisation
###########
# model structure: hyperparameters (alpha, Beta) -> parameter lambda[1] -> observation 0
#                                                              lambda[2]-> observation 1
#                                                              lambda[3]-> observation 2
#                                                              lambda[4]-> observation 3
# "observed" data from Poisson likelihood
y = 0:3

# hyperparameters
alpha = 6
Beta = 5

# moment generating function for gamma prior
mgf <- function(t, alpha, Beta) {
  Beta^alpha / (Beta - t)^alpha
}
mgf_first_derivative <- function(t, alpha, Beta) {
  alpha * mgf(t, alpha, Beta) * (Beta - t)^(-1)
}
mgf_second_derivative <- function(t, alpha, Beta) {
  (alpha + 1) * mgf_first_derivative(t, alpha, Beta) * (Beta - t)^(-1)
}
mgf_third_derivative <- function(t, alpha, Beta) {
  (alpha + 2) * mgf_second_derivative(t, alpha, Beta) * (Beta - t)^(-1)
}

# mgf-marginalisation
marginalisation <- function(y, t, alpha, Beta) {
  mgf(t, alpha, Beta) * 
    mgf_first_derivative(t, alpha, Beta) *
    mgf_second_derivative(t, alpha, Beta) * 
    mgf_third_derivative(t, alpha, Beta) *
    prod(1 / factorial(y))
}
(m <- marginalisation(y, t = -1, alpha = alpha, Beta = Beta))

# verification via conjugacy + Chib's method
# posterior kernel
Lambda = 1
post_kernel <- function(lambda, alpha, Beta, y) {
  dgamma(lambda, shape = alpha, rate = Beta) * prod(dpois(y, lambda = lambda))
}
Chib_marginalisation <- function(i) {
  post_kernel(Lambda, alpha, Beta, y[i]) / 
    dgamma(Lambda, shape = alpha + y[i], rate = Beta + 1)
} 
# applying Chib's method 4 times (once for each sub model) for the 4 marginal likelihoods
(p <- prod(sapply(1:4, Chib_marginalisation)))
round(p, 14) == round(m, 14) # TRUE

# verification via NegBin
round(prod(dnbinom(y, size=alpha, prob=Beta/(Beta+1))), 15) == round(m, 15)

###############################################
############
# example 3: marginal likelihood calculation
###########
# model structure: hyperparameters (alpha, Beta) -> parameter lambda -> observation 0
#                                                                                   0
#                                                                                   1
#                                                                                   2
# "observed" data from Poisson likelihood
y = c(0, 0:2)

# hyperparameters
alpha = 4
Beta = 6

# moment generating function for gamma prior
mgf <- function(t, alpha, Beta) {
  Beta^alpha / (Beta - t)^alpha
}
mgf_first_derivative <- function(t, alpha, Beta) {
  alpha * mgf(t, alpha, Beta) * (Beta - t)^(-1)
}
mgf_second_derivative <- function(t, alpha, Beta) {
  (alpha + 1) * mgf_first_derivative(t, alpha, Beta) * (Beta - t)^(-1)
}
mgf_third_derivative <- function(t, alpha, Beta) {
  (alpha + 2) * mgf_second_derivative(t, alpha, Beta) * (Beta - t)^(-1)
}

# mgf marginalisation
marginal_likelihood <- function(y, t, alpha, Beta) {
  mgf_third_derivative(t, alpha, Beta) *
    prod(1/factorial(y))
}
(m <- marginal_likelihood(y, t = -4, alpha = alpha, Beta = Beta))


# verification via conjugacy + Chib's method
# posterior kernel
Lambda = 1 # choose a point in the sample space of lambda
post_kernel <- function(lambda, alpha, Beta, y) {
  dgamma(lambda, shape = alpha, rate = Beta) * prod(dpois(y, lambda = lambda))
}
(p <- post_kernel(Lambda, alpha, Beta, y) / 
    dgamma(Lambda, shape = alpha + sum(y), rate = Beta + length(y)))
round(p, 16) == round(m, 16) # TRUE

################################################
###########
# example 4: linearly transformed hierarchical model marginalisation
###########
# "observed" data from Poisson likelihood
y = c(0, 1, 0, 2, 3)

# hyperparameters
o <- optim(par = rep(1, 2), 
           fn = function(ab) -sum(dgamma(x = c(1, 20/9, 10/3), shape = ab[1], rate = ab[2], log = TRUE)),
           hessian = TRUE)
(alpha = round(o$par[1], digits = 1))
(Beta = round(o$par[2]))

# mgf-marginalisation
A <- rbind(c(0.1, 0, 0),
           c(0.9, 0.1, 0),
           c(0, 0.1, 0),
           c(0, 0.8, 0.1),
           c(0, 0, 0.9))
marginalisation <- function(y, tvec, alpha, Beta, A) {
  a <- Beta - tvec %*% A[, 1]
  b <- Beta - tvec %*% A[, 2]
  c <- Beta - tvec %*% A[, 3]
  
  (1 / prod(factorial(y))) * Beta^(3 * alpha) * (-alpha) * A[5, 3]^3 *
    ((alpha^2 * (alpha + 1)^2 * (-alpha-2) * A[2, 1] * A[4, 2]^2) / (a^(alpha+1) * b^(alpha+2) * c^(alpha+3)) + 
       (2 * alpha^2 * (-alpha-1) * (-alpha-2) * (-alpha-3) * A[2, 1] * A[4, 3] * A[4, 2]) / (a^(alpha+1) * b^(alpha+1) * c^(alpha+4)) + 
       ((-alpha) * (alpha+1)^2 * (alpha+2)^2 * A[2, 2] * A[4, 2]^2) / (a^alpha * b^(alpha+3) * c^(alpha+3)) +
       ((-alpha) * (-alpha-1) * (-alpha-2) * (-alpha-3) * (-alpha-4) * A[2, 1] * A[4, 3]^2) / (a^(alpha+1) * b^alpha * c^(alpha+5)) +
       (2 * (-alpha) * (alpha+1)^2 * (-alpha-2) * (-alpha-3) * A[2, 2] * A[4, 3] * A[4, 2]) / (a^alpha * b^(alpha+2) * c^(alpha+4)) + 
       ((-alpha) * (-alpha-1) * (-alpha-2) * (-alpha-3) * (-alpha-4) * A[2, 2] * A[4, 3]^2) / (a^alpha * b^(alpha+1) * c^(alpha+5)))
}
(m <- marginalisation(y, tvec = rep(-1, 5), alpha = alpha, Beta = Beta, A = A))

# verification via simulation
set.seed(42)
iterations <- 1e6
region_rates <- matrix(rgamma(3 * iterations, shape = alpha, rate = Beta), nrow = 3)
region_obs <- matrix(rpois(5 * iterations, lambda = t(A %*% region_rates)), ncol = 5)
colMeans(region_obs) # should be around (0.225, 2.25, 0.225, 2.025, 2.025)
region_y <- subset(region_obs, region_obs[, 1] == y[1] & region_obs[, 2] == y[2] & region_obs[, 3] == y[3] & region_obs[, 4] == y[4] & region_obs[, 5] == y[5])
nrow(region_y) # number of times the observed Y happen
nrow(region_y) / iterations # frequency the observed Y happens
(freqq <- qbinom(c(0.025, 0.975), size = iterations, prob = m)) # theoratical exact confidence interval
nrow(region_y) >= freqq[1] & nrow(region_y) <= freqq[2] # the frequency Y happens is within the CI
#####################################
############
# example 5: pump failure hierarchical model marginalisation
############
# observations
ti <- c(94.32, 15.72, 62.88, 125.76, 5.24, 31.44, 1.048, 1.048, 2.096, 10.48)
y <- c(5, 1, 5, 14, 3, 19, 1, 1, 4, 22)

# mgf marginalisation
log_gamma_mgf_derivative <- function(order, alpha, beta, s) {
  lgamma(alpha + order) - lgamma(alpha) + alpha * log(beta) - (alpha + order) * log(beta - s)
}
log_marginalisation <- function(alpha, beta, ti, y) {
  sum(y * log(ti) - lfactorial(y) + sapply(1:length(ti), function(i) log_gamma_mgf_derivative(y[i], alpha, beta, -ti[i])))
}
(m <- exp(log_marginalisation(alpha = 1.27, beta = 0.82, ti, y)))

# NegBin verification
(p <- exp(sum(dnbinom(y, size = 1.27, prob = 0.82 / (0.82 + ti), log = TRUE))))
round(p, 29) == round(m, 29) # TRUE

##################################
###########
# example 6: Pareto pump failure marginal likelihood (analytical)
###########
sum(y)
sum(ti)
exp(sum(y * log(ti) - lfactorial(y)))
