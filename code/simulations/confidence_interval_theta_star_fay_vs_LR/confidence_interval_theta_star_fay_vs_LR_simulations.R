library(tidyverse)
library(scales)
library(asht)
library(cowplot)
library(fs)
library(furrr)
source("WprevSeSp_SRS.R")

plan(multisession, workers = 40)
furrr_options(seed = 200)

n_replications <- 10000

simulation_design <-
  crossing(
    n_tested_for_prevalence = 100,
    n_tested_for_sensitivity = 60,
    n_tested_for_specificity = 300,
    prevalence = seq(0.005, 0.02, by = 0.005),
    sensitivity = seq(0.75, 1, by = 0.05),
    specificity = seq(0.75, 1, by = 0.05)) %>%
  slice(rep(1:n(), each = n_replications))

n_experiments <- nrow(simulation_design)

set.seed(200)
experiments <-
  simulation_design %>%
  # Simulate AP
  mutate(n_known_positive = map2_int(n_tested_for_prevalence, prevalence, ~rbinom(1, .x, .y))) %>%
  mutate(n_known_negative = n_tested_for_prevalence - n_known_positive) %>%
  mutate(n_true_positive = map2_int(n_known_positive, sensitivity, ~rbinom(1, .x, .y))) %>%
  mutate(n_true_negative = map2_int(n_known_negative, specificity, ~rbinom(1, .x, .y))) %>%
  mutate(n_false_positive = n_known_negative - n_true_negative,
         est_prevalence = (n_true_positive + n_false_positive) / n_tested_for_prevalence) %>%
  # Simulate Se and SP
  mutate(est_specificity = map2_dbl(n_tested_for_specificity, specificity, ~rbinom(1, .x, .y) / .x),
         est_sensitivity = map2_dbl(n_tested_for_sensitivity, sensitivity, ~rbinom(1, .x, .y) / .x)) %>%
  select(names(simulation_design), starts_with("est"))

results <-
  experiments %>%
  mutate(LR = future_pmap(.l = list(est_prevalence = est_prevalence,
                                    n_tested_for_prevalence = n_tested_for_prevalence,
                                    est_sensitivity = est_sensitivity,
                                    n_tested_for_sensitivity = n_tested_for_sensitivity,
                                    est_specificity = est_specificity,
                                    n_tested_for_specificity = n_tested_for_specificity),
                          .f = function(est_prevalence,
                                        n_tested_for_prevalence,
                                        est_sensitivity,
                                        n_tested_for_sensitivity,
                                        est_specificity,
                                        n_tested_for_specificity)
                            prevSeSp(AP = est_prevalence,
                                     nP = n_tested_for_prevalence,
                                     Se = est_sensitivity,
                                     nSe = n_tested_for_sensitivity,
                                     Sp = est_specificity,
                                     nSp = n_tested_for_specificity,
                                     conf.level = 0.95)),
         Fay = future_pmap(.l = list(est_prevalence = est_prevalence,
                                     n_tested_for_prevalence = n_tested_for_prevalence,
                                     est_sensitivity = est_sensitivity,
                                     n_tested_for_sensitivity = n_tested_for_sensitivity,
                                     est_specificity = est_specificity,
                                     n_tested_for_specificity = n_tested_for_specificity),
                           .f = function(est_prevalence,
                                         n_tested_for_prevalence,
                                         est_sensitivity,
                                         n_tested_for_sensitivity,
                                         est_specificity,
                                         n_tested_for_specificity)
                             WprevSeSp_SRS(AP = est_prevalence,
                                           nP = n_tested_for_prevalence,
                                           Se = est_sensitivity,
                                           nSe = n_tested_for_sensitivity,
                                           Sp = est_specificity,
                                           nSp = n_tested_for_specificity,
                                           conf.level = 0.95, nmc = 5e4))) %>%
  mutate(LR_estimate = LR %>% map_dbl("estimate"),
         LR_conf.int.lower = LR %>% map("conf.int") %>% map_dbl(1),
         LR_conf.int.upper = LR %>% map("conf.int") %>% map_dbl(2),
         Fay_estimate = Fay %>% map_dbl("estimate"),
         Fay_conf.int.lower = Fay %>% map("conf.int") %>% map_dbl(1),
         Fay_conf.int.upper = Fay %>% map("conf.int") %>% map_dbl(2)) %>%
  mutate(LR_covered = LR_conf.int.lower <= prevalence &
           prevalence <= LR_conf.int.upper,
         Fay_covered = Fay_conf.int.lower <= prevalence &
           prevalence <= Fay_conf.int.upper)

write_rds(results, "/dfs6/pub/bayerd/sim_results_10000.rds")
