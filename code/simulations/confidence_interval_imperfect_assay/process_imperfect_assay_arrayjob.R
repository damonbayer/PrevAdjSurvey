library(tidyverse)
library(scales)
library(glue)
library(cowplot)
library(fs)



results_summary <-
  read_csv("code/simulations/confidence_interval_imperfect_assay/combined_imperfect_assay_summary.csv") %>%
  left_join(read_csv("code/simulations/confidence_interval_imperfect_assay/experimental_design.csv")) %>%
  mutate(method = str_remove(method, "_result")) %>%
  pivot_longer(cols = c(lower_error_freq, upper_error_freq, coverage)) %>%
  mutate(name = fct_recode(name,
                           "Lower Error Frequency" = "lower_error_freq",
                           "Upper Error Frequency" = "upper_error_freq",
                           "Coverage" = "coverage")) %>%
  mutate(method = fct_recode(method,
                             "WprevSeSp Binomial" = "WprevSeSp",
                             "WprevSeSp Poisson" = "WprevSeSp_gamma")) %>%
  mutate(group_distribution = group_distribution %>%
           fct_relevel("high", "uniform", "low") %>%
           fct_relabel(str_to_title))







generate_plot_2 <- function(name_to_plot, n_groups_to_plot, prev_to_plot){
  group_size <- if_else(n_groups_to_plot == 8000, 1, 200)

  generated_plot <-
    results_summary %>%
    filter(n_groups == n_groups_to_plot,
           name == name_to_plot,
           prev == prev_to_plot,
           specificity %in% c(0.8, 0.9, 1)) %>%
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
    facet_grid(group_distribution ~ specificity,
               scales = "free_y",
               labeller = labeller(
                 specificity = function(x) str_c("Specificity\n", scales::percent(as.numeric(x), accuracy = 1)),
                 group_distribution = function(x) str_c("Groups with Prev.\n", x, " Weights"))) +
    geom_hline(yintercept = case_when(name_to_plot == "Coverage" ~ 0.95,
                                      name_to_plot == "Upper Error Frequency" ~ 0.025,
                                      name_to_plot == "Lower Error Frequency" ~ 0.025),
               linetype = "dashed") +
    geom_point(size = 1, alpha = 0.25) +
    geom_smooth(
      method="glm",
      method.args=list(family="binomial"),
      formula = cbind(successes, failures) ~ splines::bs(x),
      se = F,
      alpha = 0.1, linetype = "dashed"
    ) +
    cowplot::theme_minimal_grid() +
    scale_x_continuous(name = "Weight Coefficient of Variation", labels = scales::percent) +
    scale_y_continuous(name = name_to_plot, label = ~percent(., accuracy = 1)) +
    scale_color_discrete(name = "Method") +
    theme(legend.position = "bottom") +
    ggtitle(label = glue("{name_to_plot} Properties for Simulations with {percent(prev_to_plot, accuracy = 0.1)} Prevalence Among {comma(n_groups_to_plot)} Groups of {comma(group_size)}"),
            subtitle = "Each Point = 10,000 Replications")
  generated_plot
}

results_summary %>%
  select(n_groups, name, prev) %>%
  distinct() %>%
  mutate(plot = pmap(
    list(name, n_groups, prev),
    ~generate_plot_2(name_to_plot = ..1,
                     n_groups_to_plot = ..2,
                     prev_to_plot = ..3))) %>%
  mutate(file_name = path("figures", str_c("imperfect", str_replace_all(str_to_lower(name), "\\s", "_"), n_groups, str_replace_all(prev, "\\.", "_"), "reduced", sep = "_"), ext = "pdf")) %>%
  with(walk2(file_name, plot, ~save_plot(filename = .x, plot = .y, ncol = 3, nrow = 2, base_width = 3.5)))
