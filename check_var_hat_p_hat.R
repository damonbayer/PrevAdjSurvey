n_replications <- 1

simulation_design_200_50[1,]
prevalence <- 0.005
ICC <- 0.0001
dirichlet_alpha <- 0.01
n_weights <- 200
samples_per_weight <- 50

simulation_design_200_50 %>%
  slice(rep(1:n(), each = n_replications)) %>%
  group_by(design) %>%
  mutate(replication = row_number()) %>%
  ungroup() %>%
  mutate(interval_results =
           future_pmap(.options = furrr_options(seed = 200),
                       list(prevalence, ICC, dirichlet_alpha, n_weights, samples_per_weight),
                       function(prevalence, ICC, dirichlet_alpha, n_weights, samples_per_weight) {
                         alpha <- (1 - ICC) / ICC * prevalence
                         beta <- (1 - ICC) / ICC * (1 - prevalence)

                         dat <- tibble(
                           weight = as.vector(rdirichlet(1, alpha = rep(dirichlet_alpha, n_weights))),
                           true_positivity = rbeta(n_weights, alpha, beta),
                           n_tested = samples_per_weight) %>%
                           mutate(sample_positives = rbinom(n = n(), size = n_tested, prob = true_positivity)) %>%
                           mutate(sample_positivity = sample_positives / n_tested)

                         mean_weight <- mean(dat$weight)
                         var_weight <- var(dat$weight)
                         coef_var <- sqrt(var_weight) / mean_weight

                         p_hat <- sum(dat$sample_positivity * dat$weight)
                         var_hat_p_hat <- sum((dat$weight / dat$n_tested)^2 * dat$sample_positives)

                         result_wspoissonTest <-
                           wspoissonTest(x = dat$sample_positives,
                                         w = dat$weight / dat$n_tested,
                                         midp = F)

                         result_wspoissonTest_midp <-
                           wspoissonTest(x = dat$sample_positives,
                                         w = dat$weight / dat$n_tested,
                                         midp = T)

                         result_AC <- AC_method(p_hat = p_hat,
                                                var_hat_p_hat = var_hat_p_hat,
                                                adjusted = F)

                         result_AC_adjusted <- AC_method(p_hat = p_hat,
                                                         var_hat_p_hat = var_hat_p_hat,
                                                         n_psu = n_weights,
                                                         adjusted = T)

                         result_CP <- CP_method(p_hat = p_hat,
                                                var_hat_p_hat = var_hat_p_hat,
                                                adjusted = F)

                         result_CP_adjusted <- CP_method(p_hat = p_hat,
                                                         var_hat_p_hat = var_hat_p_hat,
                                                         n_psu = n_weights,
                                                         adjusted = T)

                         list(mean_weight = mean_weight,
                              var_weight = var_weight,
                              coef_var = coef_var,
                              result_wspoissonTest = result_wspoissonTest,
                              result_wspoissonTest_midp = result_wspoissonTest_midp,
                              result_AC = result_AC,
                              result_AC_adjusted = result_AC_adjusted,
                              result_CP = result_CP,
                              result_CP_adjusted = result_CP_adjusted)
                       }))
