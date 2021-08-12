simulate_weights <- function(n_weights, coef_var) {
  # mu <- 1 / n_weights
  # sigma2 <- (coef_var / n_weights)^2
  # alpha <- ((1 - mu) / sigma2 - 1 / mu) * mu^2
  # beta <- alpha * (1 / mu - 1)

  alpha <- 1 / coef_var^2 - 1 / (n_weights * coef_var^2) - 1 / n_weights
  beta <- alpha * (n_weights - 1)

  x <- rbeta(n_weights, alpha, beta)
  if(sum(is.na(x))) return(NA)
  x / sum(x)
}
n_weights <- 10
coef_var <- 3



set.seed(200)
tmp <-
  crossing(n_weights = 10^(1:5),
         coef_var = c(0.1, 0.2, 0.3, 0.5, 1, 2, 3, 5, 10)) %>%
  slice(rep(1:n(), times = 1)) %>%
  mutate(x = map2(n_weights, coef_var, ~simulate_weights(.x, .y)))


tmp %>%
  mutate()
# k is the number of weights
k <- 8000

# coef_var is standard deviation of the weights / the mean of the weights
coef_var <- 0.2

# since we know the mean of the weights is 1 / k, the coef_var determines the variance of the weights
mu <- 1 / k
sigma2 <- (coef_var / k)^2

# Then we can sample from beta distribution with the above mean and variance
alpha <- ((1 - mu) / sigma2 - 1 / mu) * mu^2
beta <- alpha * (1 / mu - 1)

x <- rbeta(k, alpha, beta)

# since mean(x) should be close to 1/k, sum(x) should be close to 1
# but not exactly, so we can normalize
x <- x / sum(x)

# Exactly mu, (duh):
mu
mean(x)
# Very close to coef_var:
coef_var
sd(x) / mean(x)

# This only works because k is huge
k <- 10
coef_var <- 0.2
mu <- 1 / k
sigma2 <- (coef_var / k)^2
alpha <- ((1 - mu) / sigma2 - 1 / mu) * mu^2
beta <- alpha * (1 / mu - 1)

x <- rbeta(k, alpha, beta)
x <- x / sum(x)

# Exactly mu, (duh):
mu
mean(x)
# Not close to coef_var:
coef_var
sd(x) / mean(x)
