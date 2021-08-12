source("generate_state_sample.R")
library(asht)
library(microbenchmark)
library(kableExtra)

AC_method <- function(y, n_survey_design, m_survey_design, adj = F, conf.level = 0.95, n_strata = 51) {
  alpha <- 1 - conf.level
  z_val <- qnorm(1 - alpha / 2)

  # df is calculated as the number of primary sampling units minus the number of
  # strata.
  t_val <- qt(p = 1 - alpha / 2, df = m_survey_design - n_strata)
  c <- z_val * z_val / 2

  p_hat <- mean(y / n_survey_design)
  var_hat_p_hat <- sum((y / n_survey_design - p_hat)^2) / (m_survey_design * (m_survey_design - 1))

  n_eff <- p_hat * (1 - p_hat) / var_hat_p_hat * ifelse(adj, (z_val / t_val)^2, 1)

  x_tilde <- p_hat * n_eff + c
  n_tilde <- n_eff + 2 * c
  p_tilde <- x_tilde / n_tilde

  ci <- p_tilde + c(-1,1) * z_val * sqrt(p_tilde * (1 - p_tilde) / n_tilde)
  attr(ci, "conf.level") <- conf.level

  rval <- list(conf.int = ci, estimate = p_tilde, p_hat = p_hat, var_hat_p_hat = var_hat_p_hat)
  class(rval) <- "htest"
  return(rval)
}

n_simulations <- 10000
# n_simulations <- 5000
m_survey_design <- 12000
n_survey_design <- 1
set.seed(200)

state_samples <-
  tibble(simulation = 1:n_simulations,
         # data_few_state = map(simulation, ~generate_state_sample(scenario = few_state_scenario, m = m_survey_design, n_tilde = n_survey_design)),
         # data_many_state = map(simulation, ~generate_state_sample(scenario = many_state_scenario, m = m_survey_design, n_tilde = n_survey_design)),
         # data_two_state = map(simulation, ~generate_state_sample(scenario = two_state_scenario, m = m_survey_design, n_tilde = n_survey_design)),
         data_one_state = map(simulation, ~generate_state_sample(scenario = one_state_scenario, m = m_survey_design, n_tilde = n_survey_design))) %>%
  pivot_longer(starts_with("data_"),
               names_prefix = "data_",
               names_to = "scenario",
               values_to = "data") %>%
  mutate(prevalence = case_when(
    scenario == "few_state" ~ sum(few_state_scenario$population_cases) / sum(few_state_scenario$population),
    scenario == "many_state" ~ sum(many_state_scenario$population_cases) / sum(many_state_scenario$population),
    scenario == "two_state" ~ sum(two_state_scenario$population_cases) / sum(two_state_scenario$population),
    scenario == "one_state" ~ sum(one_state_scenario$population_cases) / sum(one_state_scenario$population))) %>%
  mutate(result_wspoissonTest = map(data, function(data) {
    dat <-
      data %>%
      group_by(state, population) %>%
      summarize(cases = sum(sample_cases),
                size = sum(sample_size),
                .groups = "drop") %>%
      mutate(group_weight = population / sum(population))

    wspoissonTest(x = dat[["cases"]], w = dat[["group_weight"]] / dat[["size"]])
  }),
  result_wspoissonTest_midp = map(data, function(data) {
    dat <-
      data %>%
      group_by(state, population) %>%
      summarize(cases = sum(sample_cases),
                size = sum(sample_size),
                .groups = "drop") %>%
      mutate(group_weight = population / sum(population))

    wspoissonTest(x = dat[["cases"]], w = dat[["group_weight"]] / dat[["size"]], midp = T)
  }),
  result_AC = map(data, ~AC_method(y = .[["sample_cases"]],
                                   n_survey_design = n_survey_design,
                                   m_survey_design = m_survey_design,
                                   adj = F,
                                   n_strata = 51)),
  result_AC_adj = map(data, ~AC_method(y = .[["sample_cases"]],
                                       n_survey_design = n_survey_design,
                                       m_survey_design = m_survey_design,
                                       adj = T,
                                       n_strata = 51))) %>%
  pivot_longer(cols = starts_with("result"),
               names_prefix = "result_",
               names_to = "method",
               values_to = "htest")


# Summarize Results -------------------------------------------------------
state_samples %>%
  mutate(conf.int = htest %>%
           map("conf.int") %>%
           map(~set_names(., c("l", "u")))) %>%
  unnest_wider(conf.int, names_sep = "_") %>%
  group_by(scenario, method) %>%
  summarize(coverage = mean(conf.int_l <= prevalence & prevalence <= conf.int_u),
            lower_error_freq = mean(conf.int_l > prevalence),
            upper_error_freq = mean(prevalence > conf.int_u),
            .groups = "drop") %>%
  arrange(scenario, method)


structure(list(scenario = c(
  "few_state", "few_state", "few_state",
  "few_state", "many_state", "many_state", "many_state", "many_state",
  "one_state", "one_state", "one_state", "one_state", "two_state",
  "two_state", "two_state", "two_state"
), method = c(
  "AC", "AC_adj",
  "wspoissonTest", "wspoissonTest_midp", "AC", "AC_adj", "wspoissonTest",
  "wspoissonTest_midp", "AC", "AC_adj", "wspoissonTest", "wspoissonTest_midp",
  "AC", "AC_adj", "wspoissonTest", "wspoissonTest_midp"
), coverage = c(
  0.957,
  0.957, 0.9629, 0.9559, 0.9515, 0.9515, 0.9597, 0.9517, 0.966,
  0.966, 0.9928, 0.9852, 0.9484, 0.9484, 0.9632, 0.9538
), lower_error_freq = c(
  0.0259,
  0.0259, 0.0216, 0.025, 0.0287, 0.0287, 0.0246, 0.0281, 0.0232,
  0.0232, 0.0039, 0.0083, 0.0256, 0.0256, 0.0193, 0.0235
), upper_error_freq = c(
  0.0171,
  0.0171, 0.0155, 0.0191, 0.0198, 0.0198, 0.0157, 0.0202, 0.0108,
  0.0108, 0.0033, 0.0065, 0.026, 0.026, 0.0175, 0.0227
)), row.names = c(
  NA,
  -16L
), class = c("tbl_df", "tbl", "data.frame")) %>%
  knitr::kable(format = "latex", booktabs = T) %>%
  kable_styling(position = "center")
