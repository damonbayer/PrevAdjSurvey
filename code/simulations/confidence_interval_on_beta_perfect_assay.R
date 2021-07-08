n_psu <- 30
n_replication <- 100
tmp<-
crossing(replication = 1:n_replication,
         ICC = 0.01,
         # ICC = c(0.005, 0.010, 0.050, 0.100, 0.150, 0.250, 0.500),
         prevalence = seq(0.005, 0.02, by = 0.005)) %>%
  mutate(psu = map2(ICC, prevalence,
                    ~{
                      alpha <- (1 - .x) / .x * .y
                      beta <- (1 - .x) / .x * (1 - .y)
                      tibble(psu_size = round(rgamma(n_psu, shape = 2, scale = 100)),
                             psu_positives = round(psu_size * rbeta(n_psu, alpha, beta)),
                             psu_prevalence = psu_positives / psu_size)}))
  # mutate(p_hat_summary = map(psu, ~{
  #   p_hat <- sum(.[["psu_size"]] * .[["psu_prevalence"]]) / sum(.[["psu_size"]])
  #   var_hat_p_hat <- sum((.[["psu_prevalence"]] - p_hat)^2) / (n_psu * (n_psu - 1))
  #   c("p_hat" = p_hat, "var_hat_p_hat" = var_hat_p_hat)
  #   })) %>%
  # unnest_wider(p_hat_summary)

tmp2<-
tmp %>%
  group_by(replication, ICC, prevalence) %>%
  unnest(psu) %>%
  summarize(wspoissonTest_result = list(wspoissonTest(x = psu_positives, w = psu_size / sum(psu_size), midp = T)),
    total_psu_size = sum(psu_size),
    .groups = "drop") %>%
  # This is wrong
  mutate(conf.int = map2(wspoissonTest_result, total_psu_size,
                         ~set_names(.x$conf.int / .y, c("low", "high")))) %>%
  unnest_wider(conf.int, names_sep = "_")


tmp2 %>%
  count(prevalence < conf.int_high)
