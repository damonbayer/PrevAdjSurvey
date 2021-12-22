library(tidyverse)

g <- function(theta_1, phi_n, phi_p) {
  if ((phi_n < phi_p) & (phi_p < theta_1)) {
    return(1)
  } else if ((phi_p >= theta_1) & (theta_1 >= phi_n)) {
    return((theta_1 - phi_n) / (phi_p - phi_n))
  } else {
    return(0)
  }
}

monotonic <- function(x) all(x == cummax(x)) | all(x == cummin(x))

tmp <-
  crossing(theta_1 = seq(0, 1, 0.005),
           phi_n = seq(0, 1, 0.01),
           phi_p = seq(0, 1, 0.01)) %>%
  filter(phi_p > phi_n) %>%
  filter(phi_p >= theta_1 & theta_1 >= phi_n) %>%
  mutate(g = pmap_dbl(list(theta_1 = theta_1, phi_n = phi_n, phi_p = phi_p),
                      function(theta_1, phi_n, phi_p) g(theta_1, phi_n, phi_p)))

tmp %>%
  group_by(theta_1, phi_n) %>%
  arrange(phi_p) %>%
  summarize(phi_p = list(phi_p), .groups = "drop") %>%
  mutate(monotonic = map_lgl(phi_p, monotonic)) %>%
  summarize(mean(monotonic))

tmp %>%
  group_by(phi_p, phi_n) %>%
  arrange(theta_1) %>%
  summarize(theta_1 = list(theta_1), .groups = "drop") %>%
  mutate(monotonic = map_lgl(theta_1, monotonic)) %>%
  summarize(mean(monotonic))

tmp %>%
  group_by(theta_1, phi_p) %>%
  arrange(phi_n) %>%
  summarize(phi_n = list(phi_n), .groups = "drop") %>%
  mutate(monotonic = map_lgl(phi_n, monotonic)) %>%
  summarize(mean(monotonic))




# -------------------------------------------------------------------------


tmp <-
  crossing(theta_1 = seq(0, 1, 0.005),
           phi_n = seq(0, 1, 0.01),
           phi_p = seq(0, 1, 0.01)) %>%
  mutate(g = pmap_dbl(list(theta_1 = theta_1, phi_n = phi_n, phi_p = phi_p),
                      function(theta_1, phi_n, phi_p) g(theta_1, phi_n, phi_p))) %>%
  replace_na(list(g = 0))

tmp %>%
  filter(phi_n == theta_1) %>%
  group_by(phi_n, theta_1) %>%
  arrange(phi_p) %>%
  summarize(monotonic(g), .groups = "drop") %>%
  summarize(mean(`monotonic(g)`))

tmp %>%
  filter(phi_n == phi_p) %>%
  group_by(phi_n, phi_p) %>%
  arrange(theta_1) %>%
  summarize(monotonic(g), .groups = "drop") %>%
  summarize(mean(`monotonic(g)`))


tmp %>%
  filter(theta_1 == phi_p) %>%
  group_by(theta_1, phi_p) %>%
  arrange(phi_n) %>%
  summarize(monotonic(g), .groups = "drop") %>%
  summarize(mean(`monotonic(g)`))

