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

#############################
##############
# example 4: cake baking
##############
# design matrices
Xbold_mat <- cbind(matrix(rep(175 + 10*(0:5), each = 45), ncol = 6))
Xbold <- c(Xbold_mat)
rbold <- cbind(c(rep(1, 15), rep(0, 30)),
               c(rep(0, 15), rep(1, 15), rep(0, 15)),
               c(rep(0, 30), rep(1, 15)))
# responses
ybold_mat <- 
  rbind(c(42, 46, 47, 39, 53, 42), # recipe 1
        c(47, 29, 35, 47, 57, 45),
        c(32, 32, 37, 43, 45, 45),
        c(26, 32, 35, 24, 39, 26),
        c(28, 30, 31, 37, 41, 47),
        c(24, 22, 22, 29, 35, 26),
        c(26, 23, 25, 27, 33, 35),
        c(24, 33, 23, 32, 31, 34),
        c(24, 27, 28, 33, 34, 23),
        c(24, 33, 27, 31, 30, 33),
        c(33, 39, 33, 28, 33, 30),
        c(28, 31, 27, 39, 35, 43),
        c(29, 28, 31, 29, 37, 33),
        c(24, 40, 29, 40, 40, 31),
        c(26, 28, 32, 25, 37, 33),
        c(39, 46, 51, 49, 55, 42), # recipe 2
        c(35, 46, 47, 39, 52, 61),
        c(34, 30, 42, 35, 42, 35),
        c(25, 26, 28, 46, 37, 37),
        c(31, 30, 29, 35, 40, 36),
        c(24, 29, 29, 29, 24, 35),
        c(22, 25, 26, 26, 29, 36),
        c(26, 23, 24, 31, 27, 37),
        c(27, 26, 32, 28, 32, 33),
        c(21, 24, 24, 27, 37, 30),
        c(20, 27, 33, 31, 28, 33),
        c(23, 28, 31, 34, 31, 29),
        c(32, 35, 30, 27, 35, 30),
        c(23, 25, 22, 19, 21, 35),
        c(21, 21, 28, 26, 27, 20),
        c(46, 44, 45, 46, 48, 63), # recipe 3
        c(43, 43, 43, 46, 47, 58),
        c(33, 24, 40, 37, 41, 38),
        c(38, 41, 38, 30, 36, 35),
        c(21, 25, 31, 35, 33, 23),
        c(24, 33, 30, 30, 37, 35),
        c(20, 21, 31, 24, 30, 33),
        c(24, 23, 21, 24, 21, 35),
        c(24, 18, 21, 26, 28, 28),
        c(26, 28, 27, 27, 35, 35),
        c(28, 25, 26, 25, 38, 28), 
        c(24, 30, 28, 35, 33, 28),
        c(28, 29, 43, 28, 33, 37),
        c(19, 22, 27, 25, 25, 35),
        c(21, 28, 25, 25, 31, 25))
ybold <- c(ybold_mat)

# Lee (2019) construction
group_separation <- 1.25 * (-1:1)
df0 <- data.frame(cbind(ybold, Xbold + rep(group_separation, each = 15), rep(1:3, each = 15), factor(1:15), rep(1:45, 6))) # the data frame for plotting only
df1 <- data.frame(cbind(ybold, Xbold, rep(1:3, each = 15), factor(1:15), rep(1:45, 6))) # the data frame for analysis
colnames(df0) <- colnames(df1) <- c("angle", "temperature", "recipe", "replication", "rr_inter")

# plotting the observations
library(ggplot2)
ggplot(data = df1, mapping = aes(x = temperature, y = angle, colour = as.character(recipe), shape = as.character(recipe))) +
  geom_point(alpha = 0.5) + # alpha blending
  scale_colour_discrete(name = "recipe", labels = c("I", "II", "III")) +
  scale_shape_discrete(name = "recipe", labels = c("I", "II", "III")) +
  labs(title = "Breaking angles of cakes vs. recipes and temperatures",
       x = "baking temperatures (degree celcius)",
       y = "cake breaking angle (degree)") +
  geom_smooth(alpha = 0.1, linewidth = 0.4, method = "glm", method.args = list(family = Gamma(link = "log")))
ggplot(data = df1, mapping = aes(x = replication, y = angle, colour = as.character(recipe), shape = as.character(recipe))) +
  geom_point(alpha = 0.5) + # alpha blending
  scale_colour_discrete(name = "recipe", labels = c("I", "II", "III")) +
  scale_shape_discrete(name = "recipe", labels = c("I", "II", "III")) +
  labs(title = "Breaking angles of cakes vs. recipes and replication order",
       x = "replication (time order)",
       y = "cake breaking angle (degree)") +
  geom_smooth(alpha = 0.1, linewidth = 0.4, method = "glm", method.args = list(family = Gamma(link = "log")))

gplot <- 
  ggplot(data = df0, mapping = aes(x = temperature, y = angle, colour = as.character(recipe), group = factor(recipe))) +
  geom_point(alpha = 0.5) + # alpha blending
  scale_colour_discrete(name = "recipe", labels = c("I", "II", "III")) +
  scale_shape_discrete(name = "recipe", labels = c("I", "II", "III")) +
  labs(title = "Breaking angles of cakes vs. recipes and temperatures",
       x = "baking temperatures (degree celcius)",
       y = "cake breaking angle (degree)")

#################################
# no-random-interaction model
#################################
library(hglm)
noint_fit <- 
  hglm2(meanmodel = angle ~ factor(recipe) * factor(temperature) + (1|replication),
        data = df1,
        family = Gamma(link = "log"),
        rand.family = inverse.gamma(link = "log"),
        method = "EQL1",
        calc.like = TRUE,
        verbose = TRUE)
(noint_s <- summary(noint_fit))
logLik(noint_fit)
plot(noint_fit)

# marginal likelihoods
log_marginal_likelihood <- function(abold, xi, alpha, ybold, rbold, Xbold) {
  m <- length(ybold) # the number of responses
  n <- ncol(rbold) # the number of random covariates
  log_bzeta <- -Xbold %*% abold
  log_front <- m * alpha * log(alpha) - m * lgamma(alpha) + (alpha - 1) * sum(log(ybold)) + alpha * sum(log_bzeta) +
    n * lgamma((m / n) * alpha + xi + 1) - n * lgamma(xi + 1)
  log_mgf_fraction <- (xi + 1) * log(xi) - ((m/n)*alpha + xi + 1) * log(xi + alpha * t(ybold * exp(log_bzeta)) %*% rbold)
  log_front + sum(log_mgf_fraction)
}

# maximum marginal likelihood estimates for fixed effects
nointXbold <- model.matrix(angle ~ factor(recipe) * factor(temperature), data = df1)
nointRbold <- model.matrix(angle ~ -1 + factor(replication), data = df1) # or lmer(angle ~ factor(recipe) * factor(temperature) + (1 | replication), data = df1) |> getME("Z") or matrix(rep(t(diag(1, nrow = 15)), 18), ncol = 15, byrow = TRUE)
noint_mmle <-
  optim(par = noint_s$FixCoefMat[, 1], 
        fn = function(abold) -log_marginal_likelihood(abold, xi = 1/noint_s$varRanef, alpha = round(1/noint_s$varFix), ybold = ybold, rbold = nointRbold, Xbold = nointXbold))
log_marginal_likelihood(noint_s$FixCoefMat[, 1], xi = 1/noint_s$varRanef, alpha = 52, ybold = ybold, rbold = nointRbold, Xbold = nointXbold)

# testing the consistency of estimates from HGLM and MMLE
noint_confint <- noint_s$FixCoefMat[, 1] + cbind(qnorm(0.025) * noint_s$FixCoefMat[, 2], qnorm(0.975) * noint_s$FixCoefMat[, 2])
# are all the estimates from maximum hglm and maximum mgf-marginal likelihood consistent?
sapply(1:length(noint_s$FixCoefMat[, 1]), function(j) noint_confint[j, 1] <= noint_mmle$par[j] & noint_mmle$par[j] <= noint_confint[j, 2]) |> all()

# graph with fitted values
noint_recipe_1_fitted <- c(noint_s$FixCoefMat[1, 1], noint_s$FixCoefMat[1, 1] + noint_s$FixCoefMat[4:8, 1])
noint_recipe_2_fitted <- noint_recipe_1_fitted + noint_s$FixCoefMat[c(2, 9, 11, 13, 15, 17), 1]
noint_recipe_3_fitted <- noint_recipe_1_fitted + noint_s$FixCoefMat[c(2, 9, 11, 13, 15, 17) + 1, 1]

noint_mmle_recipe_1_fitted <- c(noint_mmle$par[1], noint_mmle$par[1] + noint_mmle$par[4:8])
noint_mmle_recipe_2_fitted <- noint_mmle_recipe_1_fitted + noint_mmle$par[c(2, 9, 11, 13, 15, 17)]
noint_mmle_recipe_3_fitted <- noint_mmle_recipe_1_fitted + noint_mmle$par[c(2, 9, 11, 13, 15, 17) + 1]

noint_fitted0 <- cbind(exp(c(noint_recipe_1_fitted, noint_recipe_2_fitted, noint_recipe_3_fitted)),
                            exp(c(noint_mmle_recipe_1_fitted, noint_mmle_recipe_2_fitted, noint_mmle_recipe_3_fitted)),
                            rep(175 + 10*(0:5), 3) + rep(group_separation, each = 6), rep(1:3, each = 6)) |> as.data.frame() # the data frame for plotting only
colnames(noint_fitted0) <- c("maphle_fitted", "mmle_fitted", "temperature", "recipe")

maphle_cakes <- gplot + geom_point(data = noint_fitted0, mapping = aes(x = temperature, y = maphle_fitted, fill = factor(recipe)), alpha = 0.7, color="darkred", shape = 23, size = 3, show.legend = FALSE) +
  geom_line(data = noint_fitted0, mapping = aes(x = temperature, y = maphle_fitted, colour = factor(recipe))) +
  labs(subtitle = "with the maximum h-likelihood estimates")

mmle_cakes <- gplot + geom_point(data = noint_fitted0, mapping = aes(x = temperature, y = mmle_fitted, fill = factor(recipe)), alpha = 0.7, color="darkred", shape = 23, size = 3, show.legend = FALSE) +
  geom_line(data = noint_fitted0, mapping = aes(x = temperature, y = mmle_fitted, colour = factor(recipe))) +
  labs(subtitle = "with the maximum marginal likelihood estimates")

pdf(file = "cakes.pdf", width = 5, height = 7)
maphle_cakes
mmle_cakes
dev.off()

#####################
# models in Lee et al. (2017) and Cochran & Cox (1957)
#####################
# HGLM in Lee et al. (2017)
lee_fit <- 
  hglm2(meanmodel = angle ~ factor(recipe) * factor(temperature) + (1|replication) + (1|rr_inter),
        data = df1,
        family = Gamma(link = "log"),
        rand.family = inverse.gamma(link = "log"),
        method = "EQL1",
        calc.like = TRUE,
        verbose = TRUE)
(s <- summary(lee_fit))
logLik(lee_fit)
plot(lee_fit)

designXbold <- model.matrix(angle ~ factor(recipe) * factor(temperature), data = df1)
library(lme4)
designRbold <- lmer(angle ~ factor(recipe) * factor(temperature) + (1 | replication) + (1 | recipe:replication), data = df1) |> getME("Z")

# test the marginal likelihood function
log_marginal_likelihood(abold = s$FixCoefMat[, 1], xi = 1/s$varRanef[2], alpha = 52, ybold = ybold, rbold = designRbold, Xbold = designXbold)


# maximum marginal likelihood estimates for fixed effects
lee_mmle <-
  optim(par = s$FixCoefMat[, 1], 
        fn = function(abold) log_marginal_likelihood(abold, xi = 219.1, alpha = 52, ybold = ybold, rbold = designRbold, Xbold = designXbold),
        method = "L-BFGS-B", control = list(fnscale = -1, trace = TRUE))
mmle <-
  optim(par = s$FixCoefMat[, 1], 
        fn = function(abold) log_marginal_likelihood(abold, xi = 1/s$varRanef[2], alpha = 52, ybold = ybold, rbold = designRbold, Xbold = designXbold),
        method = "L-BFGS-B", control = list(fnscale = -1, trace = TRUE))

# testing the consistency of estimates
lee_confint <- s$FixCoefMat[, 1] + cbind(qnorm(0.025) * s$FixCoefMat[, 2], qnorm(0.975) * s$FixCoefMat[, 2])
# are all the estimates from maximum hglm and maximum mgf-marginal likelihood consistent?
sapply(1:length(s$FixCoefMat[, 1]), function(j) lee_confint[j, 1] <= mmle$par[j] & mmle$par[j] <= lee_confint[j, 2])
# likelihood-ratio test for the two models
lrt(lee_fit, noint_fit)

# graph with fitted values
lee_recipe_1_fitted <- c(s$FixCoefMat[1, 1], s$FixCoefMat[1, 1] + s$FixCoefMat[4:8, 1])
lee_recipe_2_fitted <- lee_recipe_1_fitted + s$FixCoefMat[c(2, 9, 11, 13, 15, 17), 1]
lee_recipe_3_fitted <- lee_recipe_1_fitted + s$FixCoefMat[c(2, 9, 11, 13, 15, 17) + 1, 1]

mmle_recipe_1_fitted <- c(mmle$par[1], mmle$par[1] + mmle$par[4:8])
mmle_recipe_2_fitted <- mmle_recipe_1_fitted + mmle$par[c(2, 9, 11, 13, 15, 17)]
mmle_recipe_3_fitted <- mmle_recipe_1_fitted + mmle$par[c(2, 9, 11, 13, 15, 17) + 1]

lee_fitted <- cbind(exp(c(lee_recipe_1_fitted, lee_recipe_2_fitted, lee_recipe_3_fitted)),
                    exp(c(mmle_recipe_1_fitted, mmle_recipe_2_fitted, mmle_recipe_3_fitted)),
                    rep(175 + 10*(0:5), 3), rep(1:3, each = 6)) # the data frame for analysis
lee_fitted0 <- cbind(exp(c(lee_recipe_1_fitted, lee_recipe_2_fitted, lee_recipe_3_fitted)),
                     exp(c(mmle_recipe_1_fitted, mmle_recipe_2_fitted, mmle_recipe_3_fitted)),
                     rep(175 + 10*(0:5), 3) + rep(2*(-1:1), each = 6), rep(1:3, each = 6)) |> as.data.frame() # the data frame for plotting only
colnames(lee_fitted0) <- colnames(lee_fitted) <- c("maphle_fitted", "mmle_fitted", "temperature", "recipe")

gplot + geom_point(data = lee_fitted0, mapping = aes(x = temperature, y = maphle_fitted, fill = factor(recipe)), alpha = 0.7, color="darkred", shape = 23, size = 3, show.legend = FALSE) +
  geom_line(data = lee_fitted0, mapping = aes(x = temperature, y = maphle_fitted, colour = factor(recipe))) +
  labs(subtitle = "with the maximum h-likelihood estimates")

gplot + geom_point(data = lee_fitted0, mapping = aes(x = temperature, y = mmle_fitted, fill = factor(recipe)), alpha = 0.7, color="darkred", shape = 23, size = 3, show.legend = FALSE) +
  geom_line(data = lee_fitted0, mapping = aes(x = temperature, y = mmle_fitted, colour = factor(recipe))) +
  labs(subtitle = "with the maximum marginal likelihood estimates")

#########################
# the HGLM corresponding to Cochran and Cox (1957)
cox_fit <- 
  hglm2(meanmodel = angle ~ factor(recipe) * factor(temperature) + factor(replication) + (1|rr_inter),
        data = df1,
        family = Gamma(link = "log"),
        rand.family = inverse.gamma(link = "log"),
        method = "EQL1",
        calc.like = TRUE,
        verbose = TRUE)
(cox_s <- summary(cox_fit))
logLik(cox_fit) 
plot(cox_fit)
# likelihood-ratio test for the two models
lrt(cox_fit, noint_fit)

# maximum marginal likelihood estimates for fixed effects
coxXbold <- model.matrix(angle ~ factor(recipe) * factor(temperature) + factor(replication), data = df1)
coxRbold <- model.matrix(angle ~ -1 + factor(rr_inter), data = df1) # alternatively, lmer(angle ~ factor(recipe) * factor(temperature) + factor(replication) + (1 | recipe:replication), data = df1) |> getME("Z")
cox_mmle <-
  optim(par = cox_s$FixCoefMat[, 1], 
        fn = function(abold) -log_marginal_likelihood(abold, xi = 1/cox_s$varRanef, alpha = round(1/cox_s$varFix), ybold = ybold, rbold = coxRbold, Xbold = coxXbold))
log_marginal_likelihood(noint_s$FixCoefMat[, 1], xi = 1/noint_s$varRanef, alpha = round(1/cox_s$varFix), ybold = ybold, rbold = nointRbold, Xbold = nointXbold)

# testing the consistency of estimates
cox_confint <- cox_s$FixCoefMat[, 1] + cbind(qnorm(0.025) * cox_s$FixCoefMat[, 2], qnorm(0.975) * cox_s$FixCoefMat[, 2])
# are all the estimates from maximum hglm and maximum mgf-marginal likelihood consistent?
sapply(1:length(cox_s$FixCoefMat[, 1]), function(j) cox_confint[j, 1] <= cox_mmle$par[j] & cox_mmle$par[j] <= cox_confint[j, 2]) |> all()
