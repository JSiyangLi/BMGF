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
# example 4: marginal likelihood calculation
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

