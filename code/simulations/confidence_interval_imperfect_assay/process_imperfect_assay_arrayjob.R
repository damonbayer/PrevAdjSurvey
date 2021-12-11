library(tidyverse)
library(scales)
library(glue)
library(cowplot)
library(fs)



results_summary <-
  read_csv("code/simulations/confidence_interval_imperfect_assay/combined_imperfect_assay_raw.csv") %>%
  left_join(read_csv("code/simulations/confidence_interval_imperfect_assay/experimental_design.csv")) %>%
  mutate(method = str_remove(method, "_result")) %>%
  pivot_longer(cols = c(conf_int_1, conf_int_2, width, lower_error_freq, upper_error_freq, coverage)) %>%
  mutate(name = fct_recode(name,
                           "Lower Error Frequency" = "lower_error_freq",
                           "Upper Error Frequency" = "upper_error_freq",
                           "Coverage" = "coverage")) %>%
  mutate(method = fct_recode(method,
                             "WprevSeSp Binomial" = "WprevSeSp",
                             "WprevSeSp Poisson" = "WprevSeSp_gamma",
                             "wsPoisson" = "wspoissonTest")) %>%
  mutate(group_distribution = group_distribution %>%
           fct_relevel("high", "uniform", "low") %>%
           fct_relabel(str_to_title)) %>%
  filter(method != "wsPoisson")




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
    theme(legend.position = "bottom",
          panel.border = element_rect(color = "black", fill = NA, size = 1)) +
    ggtitle(label = glue("{name_to_plot} Properties for Simulations with {percent(prev_to_plot, accuracy = 0.1)} Prevalence Among {comma(n_groups_to_plot)} Groups of {comma(group_size)}"),
            subtitle = "Each Facet = 95% Sensitivity, Each Point = 10,000 Replications")
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
  mutate(file_name = path("figures", str_c("imperfect", str_replace_all(str_to_lower(name), "\\s", "_"), n_groups, "groups", str_replace_all(prev, "\\.", "_"), "prev", sep = "_"), ext = "pdf")) %>%
  # head(1) %>%
  with(walk2(file_name, plot, ~save_plot(filename = .x, plot = .y, ncol = 3, nrow = 3, base_asp = 1)))


# dir_ls("figures/") %>%
#   enframe(name = NULL, value = "fig_path") %>%
#   mutate(fig_name = fig_path %>% path_file() %>% path_ext_remove()) %>%
#   mutate(latex_code = map2(fig_path, fig_name,
#                            ~c("\\begin{figure}",
#                               "\\centering",
#                               str_c("\\includegraphics[width=\\textwidth]{", .x, "}"),
#                               "\\caption{Caption}",
#                               str_c("\\label{fig:", .y, "}"),
#                               "\\end{figure}", ""))) %>%
#   pull(latex_code) %>%
#   unlist() %>%
#   write_lines("~/Desktop/untitled.txt")
