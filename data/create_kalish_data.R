library(tidyverse)
library(fs)
source("code/simulations/confidence_interval_imperfect_assay/WprevSeSp_SRS.R")
load("C:/Users/bayerdm/Downloads/dataForDamon.RData")

kalish_data <-
  data.for.Damon %>%
  as_tibble() %>%
  mutate(y = as.integer(y),
         hispanic = hispanic == "Yes")

dir_create("data")
write_csv(kalish_data, "data/kalish_data.csv")

n_groups <- nrow(kalish_data)
tests_per_group <- 1

n_tested_for_prevalence <- n_groups * tests_per_group

n_tested_for_specificity <- 300
n_tested_for_sensitivity <- 56
# https://github.com/niaid/SARSCoV2-US-serosurvey-2020/blob/09ef4d285e2519e0501c9f80d4f9251595e97305/Figure3_Figure4/R/ForestPlots.R
est_specificity <- 1
est_sensitivity <- 1


# Full Data ---------------------------------------------------------------
apparent_positive_counts <- kalish_data$y
weights <- kalish_data$w
apparent_prevalence <- sum(apparent_positive_counts / tests_per_group * weights)
var_hat_p_hat <- sum((weights / tests_per_group)^2 * apparent_positive_counts)
stdErrPrev <- sqrt(var_hat_p_hat)
coef_var <- sd(weights) / mean(weights)
dput(coef_var)
2.52363306559155
adjusted_prevalence <- prevAdj(AP = apparent_prevalence,
                               sen = est_sensitivity,
                               spec = est_specificity)


WprevSeSp_result <-
  WprevSeSp_original(AP = apparent_prevalence,
                     stdErrPrev = stdErrPrev,
                     nP = n_tested_for_prevalence,
                     Se = est_sensitivity,
                     nSe = n_tested_for_sensitivity,
                     Sp = est_specificity,
                     nSp = n_tested_for_specificity)

dput(WprevSeSp_result)
structure(list(estimate = c(`adjusted prevalence` = 0.0455598302324494),
               statistic = c(`Sensitivity (using nSe=56)` = 1), parameter = c(`Specificity (using nSp=300)` = 1),
               conf.int = structure(c(`2.5%` = 0.0252636260762333, `97.5%` = 0.0668264610875286
               ), conf.level = 0.95), data.name = "Unadjusted prevalence=0.04556(se=0.008562)",
               method = "Prevalence Adjusted for Sensitivity and Specificity (CI by Korn-Graubard with melding with negatives set to zero)",
               nPeff = 8058), class = "htest")

w <- weights / tests_per_group
WprevSeSp_gamma_result <- WprevSeSp_gamma(apparent_positive_counts = apparent_positive_counts,
                                          w = w,
                                          nP = n_tested_for_prevalence,
                                          Se = est_sensitivity,
                                          nSe = n_tested_for_sensitivity,
                                          Sp = est_specificity,
                                          nSp = n_tested_for_specificity)
dput(WprevSeSp_gamma_result)
structure(list(estimate = c(`adjusted prevalence` = 0.0455598302324494),
               statistic = c(`Sensitivity (using nSe=56)` = 1), parameter = c(`Specificity (using nSp=300)` = 1),
               conf.int = structure(c(`2.5%` = 0.025569736751475, `97.5%` = 0.0753507157282565
               ), conf.level = 0.95), data.name = "Unadjusted prevalence=0.04556",
               method = "Prevalence Adjusted for Sensitivity and Specificity (CI by Korn-Graubard with melding with negatives set to zero)",
               nPeff = 8058), class = "htest")
WprevSeSp_result$conf.int
WprevSeSp_gamma_result$conf.int


# Hispanic Data -----------------------------------------------------------
kalish_data <-
  kalish_data %>%
  filter(hispanic)
apparent_positive_counts <- kalish_data$y
weights <- kalish_data$w / sum(kalish_data$w)
apparent_prevalence <- sum(apparent_positive_counts / tests_per_group * weights)
var_hat_p_hat <- sum((weights / tests_per_group)^2 * apparent_positive_counts)
stdErrPrev <- sqrt(var_hat_p_hat)
coef_var <- sd(weights) / mean(weights)
dput(coef_var)
3.05941523379367

adjusted_prevalence <- prevAdj(AP = apparent_prevalence,
                               sen = est_sensitivity,
                               spec = est_specificity)


WprevSeSp_result <-
  WprevSeSp_original(AP = apparent_prevalence,
                     stdErrPrev = stdErrPrev,
                     nP = n_tested_for_prevalence,
                     Se = est_sensitivity,
                     nSe = n_tested_for_sensitivity,
                     Sp = est_specificity,
                     nSp = n_tested_for_specificity)

dput(WprevSeSp_result)
structure(list(estimate = c(`adjusted prevalence` = 0.0612791931338004),
               statistic = c(`Sensitivity (using nSe=56)` = 1), parameter = c(`Specificity (using nSp=300)` = 1),
               conf.int = structure(c(`2.5%` = 0.0234604022334971, `97.5%` = 0.11753836349672
               ), conf.level = 0.95), data.name = "Unadjusted prevalence=0.06128(se=0.02045)",
               method = "Prevalence Adjusted for Sensitivity and Specificity (CI by Korn-Graubard with melding with negatives set to zero)",
               nPeff = 8058), class = "htest")



w <- weights / tests_per_group

WprevSeSp_gamma_result <- WprevSeSp_gamma(apparent_positive_counts = apparent_positive_counts,
                                          w = w,
                                          nP = n_tested_for_prevalence,
                                          Se = est_sensitivity,
                                          nSe = n_tested_for_sensitivity,
                                          Sp = est_specificity,
                                          nSp = n_tested_for_specificity)
dput(WprevSeSp_gamma_result)
structure(list(estimate = c(`adjusted prevalence` = 0.0612791931338004),
               statistic = c(`Sensitivity (using nSe=56)` = 1), parameter = c(`Specificity (using nSp=300)` = 1),
               conf.int = structure(c(`2.5%` = 0.023961917726667, `97.5%` = 0.200151985872185
               ), conf.level = 0.95), data.name = "Unadjusted prevalence=0.06128",
               method = "Prevalence Adjusted for Sensitivity and Specificity (CI by Korn-Graubard with melding with negatives set to zero)",
               nPeff = 8058), class = "htest")
WprevSeSp_result$conf.int
WprevSeSp_gamma_result$conf.int
