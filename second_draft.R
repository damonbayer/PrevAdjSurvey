library(extraDistr)
library(tidyverse)
n_replications <- 1000
survey_sample_size <- 100
n_weights <- 30
ICC <- 0.1
prevalence <- 0.005
alpha <- (1 - ICC) / ICC * prevalence
beta <- (1 - ICC) / ICC * (1 - prevalence)



tmp <- tibble(group = rep(1:n_replications, each = n_weights),
       weights = as.vector(t(rdirichlet(n_replications, alpha = rep(1, n_weights)))),
       true_positivity = rbeta(n_weights * n_replications, alpha, beta),
       sample_size = survey_sample_size) %>%
  mutate(sample_positives = rbinom(n = n_weights * n_replications, size = sample_size, prob = true_positivity)) %>%
  mutate(sample_positivity = sample_positives / sample_size)


tmp2 <- tmp %>%
  group_by(group) %>%
  summarize(result_wspoissonTest_midp = list(wspoissonTest(x = cur_data()[["sample_positives"]], w = cur_data()[["weights"]] / cur_data()[["sample_size"]], midp = F)),
            result_wspoissonTest_midp = list(wspoissonTest(x = cur_data()[["sample_positives"]], w = cur_data()[["weights"]] / cur_data()[["sample_size"]], midp = T)))


tmp %>%
  group_by(group) %>%
  summarize(p_hat = sum(cur_data()[["sample_positivity"]] * cur_data()[["weights"]]))

tmp2 %>%
  mutate(conf.int = wspoissonTest %>% map("conf.int") %>% map(~set_names(., c("l", "u")))) %>%
  unnest_wider(conf.int) %>%
  mutate(lower_error = l > 0.05,
         upper_error = u > 0.05,
         covered = !(lower_error | upper_error)) %>%
  summarize(mean(lower_error), mean(upper_error), mean(covered))


# y = sum s * x / n
# y = sum weigts * sample_positivity


# -------------------------------------------------------------------------



wspoissonTest(x = tmp[["sample_positives"]], w = tmp[["weights"]] / tmp[["sample_size"]], midp = T)
wspoissonTest(x = tmp[["sample_positives"]], w = tmp[["weights"]] / tmp[["sample_size"]], midp = T)

wspoissonTest(x = tmp$sample_positives, w = tmp$weights / tmp$sample_size)$conf.int
wspoissonTest(x = dat[["cases"]], w = dat[["group_weight"]] / dat[["size"]])


tmp %>%
  group_by(group) %>%
  summarize(group_positivity = sum(sample_positivity * weights)) %>%
  summarize(mean(group_positivity))


  group_by(group) %>%
  summarize(sum(weights))


map_dbl(1:n_replications, ~{
  tibble(weights = as.vector(rdirichlet(1, alpha = rep(20, n_weights))),
         true_positivity = rbeta(n_weights, alpha, beta),
         sample_size = survey_sample_size) %>%
    mutate(positives = rbinom(n = n_weights, size = sample_size, prob = true_positivity)) %>%
    summarize(a = mean(positives / sample_size)) %>%
    pull(a)
})

# Generate weights ike this:
rdirichlet(1, alpha = rep(20, n_weights))
# smaller alpha => larger variance


# Generate positivity rates like this
# ICC, Prevalence




rbeta(n_weights, alpha, beta)



tibble(psu_size = round(rgamma(n_psu, shape = 2, scale = 100)),
       psu_positives = round(psu_size * rbeta(n_psu, alpha, beta)),
       psu_prevalence = psu_positives / psu_size)}))
# Generate sample sizes like this?


tibble()

as.vector(rdirichlet(1, alpha = rep(20, n_weights)))

rdirichlet(n_replications, alpha = rep(20, n_weights)) %>%
  t() %>%
  as_tibble

rdirichlet(n_replications, alpha = rep(100, n_weights)) %>%
  apply(1, var) %>%
  mean()
1.095742e-05


rdirichlet(n_replications, alpha = rep(1, n_weights)) %>%
  apply(1, var) %>%
  mean()
0.001066417

rdirichlet(n_replications, alpha = rep(.1, n_weights)) %>%
  apply(1, var) %>%
  mean()

rdirmnom(n = 100, size = 1000, alpha = rep(0.5, 30)) %>%
  `/`(1000) %>%
  apply(1, var) %>%
  mean()
