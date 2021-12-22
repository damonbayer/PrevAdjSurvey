library(tidyverse)
library(scales)
library(glue)
library(cowplot)
library(fs)


results_summary <-
  read_csv("code/simulations/confidence_interval_perfect_assay/combined_perfect_assay_raw.csv") %>%
  left_join(read_rds("code/simulations/confidence_interval_perfect_assay/experimental_design.rds") %>%
              select_if(~!is.list(.))) %>%
  mutate(method = str_remove(method, "_result")) %>%
  pivot_longer(cols = c(conf_int_1, conf_int_2, width, lower_error_freq, upper_error_freq, coverage)) %>%
  mutate(name = fct_recode(name,
                           "Lower Error Frequency" = "lower_error_freq",
                           "Upper Error Frequency" = "upper_error_freq",
                           "Coverage" = "coverage",
                           "Lower Confidence Limit" = "conf_int_1",
                           "Upper Confidence Limit" = "conf_int_2",
                           "Confidence Interval Width" = "width")) %>%
  # mutate(method = fct_recode(method,
  #                            "Agresti-Coull (Unadjusted)" = "AC",
  #                            "Agresti-Coull (Adjusted)" = "AC_adjusted",
  #                            "Clopper-Pearson (Unadjusted)" = "CP",
  #                            "Clopper-Pearson (Adjusted)" = "CP_adjusted",
  #                            "wsPoisson" = "wspoissonTest",
  #                            "wsPoisson with mid-p" = "wspoissonTest_midp")) %>%
  filter(method != "wspoissonTest_midp") %>%
  mutate(method = fct_recode(method,
                             "Agresti-Coull" = "AC",
                             "Korn-Graubard" = "CP",
                             "wsPoisson" = "wspoissonTest")) %>%
  mutate(group_distribution = group_distribution %>%
           fct_relevel("high", "uniform", "low") %>%
           fct_relabel(str_to_title))


generate_plot <- function(name_to_plot, n_groups_to_plot, prev_to_plot){
  group_size <- if_else(n_groups_to_plot == 8000, 1, 200)

  generated_plot <-
    results_summary %>%
    filter(n_groups == n_groups_to_plot,
           name == name_to_plot,
           prev == prev_to_plot,
           prop_groups_with_prev %in% c(0.05, 0.25, 0.75)) %>%
    mutate(
      tests = 10000,
      successes = as.integer(value * tests),
      failures = tests - successes) %>%
    ggplot(aes(x = coef_var,
               y = successes / tests,
               color = method,
               group = method,
               # linetype = method,
               successes = successes,
               failures = failures)) +
    facet_grid(group_distribution ~ prop_groups_with_prev,
               scales = "free_y",
               labeller = labeller(
                 prop_groups_with_prev = function(x) x %>% as.numeric() %>% percent(accuracy = 1) %>% str_c("Weights with Non-0 Prev.\n", .),
                 group_distribution = function(x) str_c("Weights with Prev.\n", x, " Weights"))) +
    geom_hline(yintercept = case_when(name_to_plot == "Coverage" ~ 0.95,
                                      name_to_plot == "Upper Error Frequency" ~ 0.025,
                                      name_to_plot == "Lower Error Frequency" ~ 0.025),
               linetype = "dashed") +
    geom_point(size = 1, alpha = 0.1) +
    geom_smooth(
      method="glm",
      method.args=list(family="binomial"),
      formula = cbind(successes, failures) ~ splines::bs(x),
      se = F,
      alpha = 0.5, linetype = "dashed"
    ) +
    cowplot::theme_minimal_grid() +
    # scale_x_continuous(name = "Weight Coefficient of Variation") +
    scale_x_continuous(name = "Weight Coefficient of Variation", labels = scales::percent) +
    # scale_y_continuous(name = name_to_plot, label = ~percent(., accuracy = 1)) +
    scale_y_continuous(name = name_to_plot, label = percent) +
    scale_color_discrete(name = "Method") +
    theme(legend.position = "bottom",
          panel.border = element_rect(color = "black", fill = NA, size = 1)) +
    ggtitle(label = glue("{name_to_plot} Properties for Simulations with {percent(prev_to_plot, accuracy = 0.1)} Prevalence Among {comma(n_groups_to_plot)} Groups of {comma(group_size)}"),
            subtitle = "Each Point = 10,000 Replications")
  generated_plot
}


results_summary %>%
  select(n_groups, name, prev) %>%
  distinct() %>%
  mutate(plot = pmap(
    list(name, n_groups, prev),
    ~generate_plot(name_to_plot = ..1,
                   n_groups_to_plot = ..2,
                   prev_to_plot = ..3))) %>%
  mutate(file_name = path("figures", str_c("perfect", str_replace_all(str_to_lower(name), "\\s", "_"), n_groups, "groups", str_replace_all(prev, "\\.", "_"), "prev", sep = "_"), ext = "pdf")) %>%
  # head(1) %>%
  with(walk2(file_name, plot, ~save_plot(filename = .x, plot = .y, ncol = 3, nrow = 3, base_asp = 1)))

