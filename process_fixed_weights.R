library(tidyverse)
library(scales)
library(fs)
library(cowplot)
results <- read_rds("fixed_weights_results_10000.rds")

# NOTE some designs are problematic
results %>%
  filter(is.na(covered)) %>%
  count(design, sort = T)

# Provisiional becasue some intervals are degerate

results %>%
  filter(design == 246) %>%
  select_if(.predicate = ~!is.logical(.)) %>%
  select(-method, -replication) %>%
  distinct()

design_key <-
  results %>%
  select_if(.predicate = ~!is.logical(.)) %>%
  select(-replication, -method) %>%
  distinct()

results_summary <-
  results %>%
  group_by(design, method) %>%
  summarize(lower_error_freq = mean(lower_error, na.rm = T),
            upper_error_freq = mean(upper_error, na.rm = T),
            coverage = mean(covered, na.rm = T),
            .groups = "drop") %>%
  left_join(design_key)

rm(results)

coef_var_levels <- c(0.01, 0.1, 0.2, 0.3, 0.5, 1, 2, 3, 5, 10)

coef_var_key <-
  tibble(coef_var = results_summary$coef_var %>% unique() %>% sort(),
         coef_var_rounded = sapply(coef_var, function(x) coef_var_levels[which.min(abs(x - coef_var_levels) / coef_var_levels)]))
# missing 100 and 10000 ugh

results_summary <-
  results_summary %>%
  left_join(coef_var_key) %>%
  mutate(coef_var_rounded = coef_var_rounded %>% percent() %>% as_factor()) %>%
  select(-coef_var) %>%
  pivot_longer(cols = c(lower_error_freq, upper_error_freq, coverage)) %>%
  mutate(name = fct_recode(name,
                           "Lower Error Frequency" = "lower_error_freq",
                           "Upper Error Frequency" = "upper_error_freq",
                           "Coverage" = "coverage"),
         method = fct_recode(method,
                             "Agresti-Coull (Unadjusted)" = "AC",
                             "Agresti-Coull (Adjusted)" = "AC_adjusted",
                             "Clopper-Pearson (Unadjusted)" = "CP",
                             "Clopper-Pearson (Adjusted)" = "CP_adjusted",
                             "wsPoisson" = "wspoissonTest",
                             "wsPoisson with mid-p" = "wspoissonTest_midp"),
         group_distribution = group_distribution %>%
           fct_relevel("high", "uniform", "low") %>%
           fct_relabel(str_to_title))



cbp2 <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

generate_plot <- function(name_to_plot, n_groups_to_plot){
  group_size <- if_else(n_groups_to_plot == 8000, 1, 200)

  generated_plot <-
    results_summary %>%
    filter(n_groups == n_groups_to_plot,
           name == name_to_plot) %>%
    ggplot(aes(coef_var_rounded, value, color = method, group = method)) +
    facet_grid(group_distribution ~ prop_groups_with_prev,
               scales = "free_y",
               labeller = labeller(
                 prop_groups_with_prev = function(x) x %>% as.numeric() %>% percent(accuracy = 1) %>% str_c("Groups with Prev.\n", .),
                 group_distribution = function(x) str_c("Groups with Prev.\n", x, " Weights"))) +
    geom_hline(yintercept = case_when(name_to_plot == "Coverage" ~ 0.95,
                                      name_to_plot == "Upper Error Frequency" ~ 0.025,
                                      name_to_plot == "Lower Error Frequency" ~ 0.025),
               linetype = "dashed") +
    geom_point(size = 4, alpha = 0.5) +
    geom_line(size = 2, alpha = 0.25) +
    cowplot::theme_minimal_grid() +
    scale_x_discrete(name = "Weight Coefficient of Variation") +
    scale_y_continuous(name = name_to_plot, label = ~percent(., accuracy = 1)) +
    scale_color_discrete(name = "Method") +
    theme(legend.position = "bottom") +
    ggtitle(label = glue("{name_to_plot} Properties for Simulations with {comma(n_groups_to_plot)} Groups of {comma(group_size)}"),
            subtitle = "Each Point = 10,000 Replications") +
    scale_color_brewer(palette = "Dark2")
    # scale_color_manual(values = cbp2)
  generated_plot
}



results_summary %>%
  select(n_groups, name) %>%
  distinct() %>%
  mutate(plot = map2(name, n_groups, ~generate_plot(.x, .y))) %>%
  mutate(file_name = map2_chr(name, n_groups, ~path(str_c(str_replace_all(str_to_lower(.x), "\\s", "_"), .y, sep = "_"), ext = "pdf"))) %>%
  mutate(walk2(plot, file_name, ~save_plot(filename = .y, plot = .x, ncol = 5, nrow = 3, base_asp = (16 /9) / (5 / 3), base_height = 5.5)))

library(glue)



tmp +
