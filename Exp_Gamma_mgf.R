dexpgamma <- function(x, shape, rate, log = FALSE) {
  log_dens <- log(rate) + log(shape) + (shape - 1) * log(x) - (shape + 1) * log(rate + x)
  if (is.na(log_dens)) {
    stop("invalid shape or rate value")
  } else {
    ifelse(log, log_dens, exp(log_dens))
  }
}
dmultiexpgamma <- function(x, shape, rate, log = FALSE) {
  n <- length(y)
  log_dens <- log(rate) + lgamma(n * shape + 1) + (shape - 1) * sum(log(y)) - (n * lgamma(shape) + (n * shape + 1) * log(rate + sum(y)))
  if (is.na(log_dens)) {
    stop("invalid shape or rate value")
  } else {
    ifelse(log, log_dens, exp(log_dens))
  }
}

############
# example 1: y=3.4, shape=1, rate=1
###########
# model structure: hyperparameter lambda -> parameter Beta -> observation 3.4
# "observed" data from gamma likelihood
y = 3.4

# hyperparameters
alpha = 1
lambda = 1

# moment generating function for exponential prior
mgf <- function(t, lambda) {
  lambda / (lambda - t)
}
mgf_first_derivative <- function(t, lambda) {
  mgf(t, lambda) / (lambda - t)
}

# mgf marginalisation
marginal_likelihood <- function(y, t, alpha, lambda) {
  mgf_first_derivative(t, lambda) *
    exp(sum((alpha - 1) * log(y) - lgamma(alpha)))
}
(m <- marginal_likelihood(y, t = -y, alpha = alpha, lambda = lambda))

# verification using exponential-gamma distribution
dexpgamma(3.4, 1, 1) == m

#########################
############
# example 2: y=(0.4, 2.2), shape=(1.5, 2), rate=0.9
###########
# model structure: hyperparameter lambda -> parameter Beta[1] -> observation 0.4
#                                                     Beta[2] -> observation 2.2
# "observed" data from gamma likelihood
y = c(0.4, 2.2)

# hyperparameters
alpha = c(1.5, 2)
lambda = 0.9

# moment generating function for exponential prior
mgf <- function(t, lambda) {
  lambda / (lambda - t)
}
mgf_first_derivative <- function(t, lambda) {
  mgf(t, lambda) / (lambda - t)
}
mgf_second_derivative <- function(t, lambda) {
  mgf_first_derivative(t, lambda) * (2 / (lambda - t))
}

# using Riemann-Liouville fractional derivative with lower limit of -Inf
mgf_1.5_derivative <- function(t, lambda) {
  0.5 * 1.5 * lambda * pi * (lambda - t)^-2.5 / gamma(0.5)
}

# mgf marginalisation
marginal_likelihood <- function(y, t, alpha, lambda) {
  mgf_1.5_derivative(t[1], lambda) *
    mgf_second_derivative(t[2], lambda) *
    exp(sum((alpha - 1) * log(y) - lgamma(alpha)))
}
(m <- marginal_likelihood(y, t = -y, alpha = alpha, lambda = lambda))

# verification using exponential-gamma distribution
round(dexpgamma(y[1], alpha[1], lambda) * dexpgamma(y[2], alpha[2], lambda), 15) == round(m, 15)

#########################
############
# example 3: y=(0.7, 2.3, 3.6), shape=2.5, rate=0.9
###########
# model structure: hyperparameter lambda -> parameter Beta -> observation 2.7
#                                                                         3.3
#                                                                         3.6
# "observed" data from gamma likelihood
y = c(2.7, 3.3, 3.6)

# hyperparameters
alpha = 0.5
lambda = 1.1

# moment generating function for exponential prior
mgf_0.5_derivative <- function(t, lambda) {
  0.5 * lambda * pi * (lambda - t)^-1.5 / gamma(0.5)
}
mgf_1.5_derivative <- function(t, lambda) {
  1.5 * 0.5 * lambda * pi * (lambda - t)^-2.5 / gamma(0.5)
}

# mgf marginalisation
marginal_likelihood <- function(y, t, alpha, lambda) {
  n <- length(y)
  mgf_1.5_derivative(t, lambda) *
    exp(sum((alpha - 1) * log(y)) - n * lgamma(alpha))
}
(m <- marginal_likelihood(y, t = -sum(y), alpha = alpha, lambda = lambda))

# verification using compound gamma
round(dmultiexpgamma(y, alpha, lambda), 18) == round(m, 18)
