library(tidyverse)
library(scales)
library(asht)
library(cowplot)
library(fs)
library(furrr)

results <- read_csv("code/simulations/confidence_interval_theta_star_fay_vs_LR/confidence_interval_theta_star_fay_vs_LR_results.csv")

# Process Results ---------------------------------------------------------
coverage_results <-
  results %>%
  group_by(n_tested_for_prevalence, prevalence, sensitivity, specificity) %>%
  summarize(LR_coverage = mean(LR_covered),
            Fay_coverage = mean(Fay_covered),
            .groups = "drop") %>%
  pivot_longer(ends_with("coverage"), names_to = "method", values_to = "coverage") %>%
  mutate(method = str_sub(method, end = -10))

coverage_distance_difference_data <-
  coverage_results %>%
  mutate(coverage_distance_from_truth = abs(0.95 - coverage)) %>%
  select(-coverage) %>%
  pivot_wider(names_from = method, values_from = coverage_distance_from_truth) %>%
  mutate(coverage_distance_difference = LR - Fay)

ci_results <-
  results %>%
  select(n_tested_for_prevalence, n_tested_for_sensitivity,
         n_tested_for_specificity, prevalence, sensitivity, specificity,
         est_prevalence, est_specificity, est_sensitivity,
         contains("conf.int")) %>%
  pivot_longer(contains("conf.int"), names_sep = "_", names_to = c("method", "ci_bound")) %>%
  mutate(ci_bound = str_sub(ci_bound, start = 10)) %>%
  pivot_wider(names_from = ci_bound, values_from = value, values_fn = list) %>%
  unnest(c(lower, upper)) %>%
  group_by(n_tested_for_prevalence, n_tested_for_sensitivity, n_tested_for_specificity, prevalence, sensitivity, specificity, method) %>%
  summarize(lower_error_freq = mean(lower > prevalence),
            upper_error_freq = mean(upper < prevalence),
            .groups = "drop")


# Plot Results ------------------------------------------------------------

coverage_comparison_plot <-
  coverage_results %>%
  ggplot(aes(specificity, sensitivity, fill = coverage)) +
  facet_grid(method ~ prevalence,
             labeller = labeller(
               prevalence = function(.) str_c("Prevalence:\n", scales::percent(as.numeric(.), accuracy = .1)),
               n_tested_for_prevalence = function(.) str_c("Tested for Prevalence Size:\n", scales::comma(as.numeric(.))))
  ) +
  # geom_raster(interpolate = F) +
  geom_tile() +
  cowplot::theme_cowplot() +
  scale_x_continuous(name = "Specificity", labels = function(.) percent(., accuracy = 1)) +
  scale_y_continuous(name = "Sensitivity", labels = function(.) percent(., accuracy = 1)) +
  scale_fill_distiller(name = "Coverage",type = "div", palette = "PiYG", direction = 1, limits = c(0.9, 1), labels = function(.) percent(., accuracy = 1)) +
  guides(fill = guide_colourbar(barwidth = 20)) +
  ggtitle("Coverage Comparison", subtitle = "10,000 replications each, nP = 100, nSp = 300, nSe = 60") +
  theme(legend.position = "bottom")

coverage_dist_limit <- max(abs(coverage_distance_difference_data$coverage_distance_difference)) + 0.01

coverage_error_comparison_plot <-
  coverage_distance_difference_data %>%
  ggplot(aes(specificity, sensitivity, fill = coverage_distance_difference)) +
  facet_grid(n_tested_for_prevalence ~ prevalence,
             labeller = labeller(
               prevalence = function(.) str_c("Prevalence:\n", scales::percent(as.numeric(.), accuracy = .1)),
               n_tested_for_prevalence = function(.) str_c("Tested for Prevalence Size:\n", scales::comma(as.numeric(.))))
  ) +
  # geom_raster(interpolate = F) +
  geom_tile() +
  cowplot::theme_cowplot() +
  scale_x_continuous(name = "Specificity", labels = function(.) percent(., accuracy = 1)) +
  scale_y_continuous(name = "Sensitivity", labels = function(.) percent(., accuracy = 1)) +
  # scale_fill_gradient2() +
  scale_fill_distiller(name = "Coverage Error Difference",type = "div", palette = "PiYG", direction = 1,
                       limits = c(-coverage_dist_limit, coverage_dist_limit), breaks = c(-coverage_dist_limit, -0.01, 0, 0.01, coverage_dist_limit),
                       labels = function(.) case_when(. == -coverage_dist_limit ~ "LR Better",
                                                      . == coverage_dist_limit ~ "Fay Better",
                                                      TRUE ~percent(., accuracy = 1))
  ) +
  guides(fill = guide_colourbar(barwidth = 20, title.position = "top")) +
  ggtitle("Coverage Error Comparison", subtitle = "10,000 replications each, nP = 100, nSp = 300, nSe = 60") +
  theme(legend.position = "bottom")

lower_error_frequency_comparison_plot <-
  ci_results %>%
  ggplot(aes(specificity, sensitivity, fill = lower_error_freq)) +
  facet_grid(method ~ prevalence,
             labeller = labeller(
               prevalence = function(.) str_c("Prevalence:\n", scales::percent(as.numeric(.), accuracy = .1)),
               n_tested_for_prevalence = function(.) str_c("Tested for Prevalence Size:\n", scales::comma(as.numeric(.))))
  ) +
  geom_tile() +
  cowplot::theme_cowplot() +
  scale_x_continuous(name = "Specificity", labels = function(.) percent(., accuracy = 1)) +
  scale_y_continuous(name = "Sensitivity", labels = function(.) percent(., accuracy = 1)) +
  scale_fill_distiller(name = "Lower Error Frequency",type = "div", limits = c(0, 0.05), palette = "PiYG", direction = -1, labels = function(.) percent(., accuracy = 1)) +
  guides(fill = guide_colourbar(barwidth = 20)) +
  ggtitle("Lower Error Frequency Comparison", subtitle = "10,000 replications each, nP = 100, nSp = 300, nSe = 60") +
  theme(legend.position = "bottom")


upper_error_frequency_comparison_plot <-
  ci_results %>%
  ggplot(aes(specificity, sensitivity, fill = upper_error_freq)) +
  facet_grid(method ~ prevalence,
             labeller = labeller(
               prevalence = function(.) str_c("Prevalence:\n", scales::percent(as.numeric(.), accuracy = .1)),
               n_tested_for_prevalence = function(.) str_c("Tested for Prevalence Size:\n", scales::comma(as.numeric(.))))
  ) +
  geom_tile() +
  cowplot::theme_cowplot() +
  scale_x_continuous(name = "Specificity", labels = function(.) percent(., accuracy = 1)) +
  scale_y_continuous(name = "Sensitivity", labels = function(.) percent(., accuracy = 1)) +
  scale_fill_distiller(name = "Upper Error Frequency",type = "div", limits = c(0, 0.05), palette = "PiYG", direction = -1, labels = function(.) percent(., accuracy = 1)) +
  guides(fill = guide_colourbar(barwidth = 20)) +
  ggtitle("Upper Error Frequency Comparison", subtitle = "10,000 replications each, nP = 100, nSp = 300, nSe = 60") +
  theme(legend.position = "bottom")

plot_save_dir <- "figures"
dir_create(plot_save_dir)

ls()[str_ends(ls(), "_plot")]

combined_plot <-
  plot_grid(coverage_comparison_plot,
            upper_error_frequency_comparison_plot,
            lower_error_frequency_comparison_plot,
            nrow = 3,
            ncol = 1,
            align = "hv")


save_plot(filename = path(plot_save_dir, "Combined_comparison_plots.pdf"), plot = combined_plot, ncol = 4, nrow = 6, base_height = 2.5)



sapply(ls()[str_ends(ls(), "_plot")],
       function(plot_name) save_plot(filename = path(plot_save_dir, plot_name, ext = "pdf"),
                                     plot = get(plot_name),
                                     base_height = 8,
                                     base_asp = 16/9))


coverage_results %>%
  group_by(method) %>%
  summarize(mean(coverage > 0.95),
            mean(coverage),
            min(coverage),
            max(coverage))

