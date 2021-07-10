source("generate_state_sample.R")

n_simulations <- 10


set.seed(200)
tmp <- tibble(simulation = 1:n_simulations,
              prevalence = sum(few_state_scenario$population_cases) / sum(few_state_scenario$population),
              data = map(simulation, ~generate_state_sample(scenario = few_state_scenario, total_people_sampled = 12000))) %>%
  mutate(wspoissonTest_result = map(data, ~wspoissonTest(x = .[["cases"]], w = .[["group_weight"]] / .[["size"]])))


dat <- tmp[[1, "data"]][[1]]

library(asht)
wspoissonTest(x = dat$cases, w = few_state_scenario$population / sum(few_state_scenario$population) / dat$size)$conf.int

library(binom)
binom.confint(x = dat$cases, n = dat$size, methods = "agresti-coull")
wspoissonTest(x = dat$cases, w = few_state_scenario$population / sum(few_state_scenario$population) / dat$size)
binom.confint

# end of section 4.1 of Dean and Pagan tells how to get var(p) estimate

AC_method <- function(x, n, adj = F, conf.level = 0.95) {
  p <- x/n
  alpha <- 1 - conf.level
  z <- qnorm(1 - alpha / 2)
  c <- z * z / 2

  .x <- x + 0.5 * z2
  .n <- n + z2
  .p <- .x/.n


  lcl <- .p - z * sqrt(.p * (1 - .p)/.n)
  ucl <- .p + z * sqrt(.p * (1 - .p)/.n)
  res.ac <- data.frame(method = rep("agresti-coull", NROW(x)),
                       xn, mean = p, lower = lcl, upper = ucl)
  res <- res.ac
  res
}

