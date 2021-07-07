# https://github.com/niaid/SARSCoV2-US-serosurvey-2020/blob/main/Figure3_Figure4/R/WprevSeSp_MF.R
### Author: Mike Fay
### apply this function to perform sensitivity and specificity adjustment

prevAdj <- function(AP, sen, spec) (AP + spec - 1) / (sen + spec - 1)

WprevSeSp_SRS <- function(AP, nP, Se, nSe, Sp, nSp, conf.level = 0.95,
                          neg.to.zero = TRUE, nmc = 1e5, seed = 49201) {
  if (!is.null(seed)) set.seed(seed)

  # This doesn't actually get used in any calculations, just the final output
  # Maybe this estimated standard error is supposed to account for Se & Sp,
  # but probably not.
  # Thus, assume nP * AP ~ Binomial(nP, true prevalence)
  # => Var(nP * AP) =  nP * true prevalence * (1 - true prevalence)
  # => Var(AP) =  (true prevalence * (1 - true prevalence)) / nP
  # => std.dev(AP) =  sqrt((true prevalence * (1 - true prevalence)) / nP)
  # => Plug-in estimate:
  stdErrPrev <- sqrt((AP * (1 - AP)) / nP)


  P <- prevAdj(AP, Se, Sp)

  xAP <- AP * nP

  # use Beta lower and upper confidence distributions
  # for apparent prevalence
  # defined using Korn-Graubard methods
  B_AP_L<- rbeta(nmc, xAP, nP - xAP + 1)
  B_AP_U<- rbeta(nmc, xAP + 1, nP - xAP)

  # use Beta lower and upper confidence distributions
  # for sensitivity (from binomial assumptions)
  B_sen_L <- rbeta(nmc, nSe * Se, nSe - nSe * Se + 1)
  B_sen_U <- rbeta(nmc, nSe * Se + 1, nSe - nSe * Se)

  # use Beta lower and upper confidence distributions
  # for specificity (from binomial assumptions)
  B_spec_L <- rbeta(nmc, nSp * Sp, nSp - nSp * Sp + 1)
  B_spec_U <- rbeta(nmc, nSp * Sp + 1, nSp - nSp * Sp)

  # Monte Carlo melding method
  lo <- prevAdj(B_AP_L, B_sen_U, B_spec_L)
  up <- prevAdj(B_AP_U, B_sen_L, B_spec_U)
  alpha <- 1 - conf.level
  LCL <- quantile(lo, prob = alpha / 2)
  UCL <- quantile(up, prob = 1 - alpha / 2)

  ci <- c(LCL, UCL)
  attr(ci, "conf.level") <- conf.level
  # when neg.to.zero=TRUE, set negative values to zero
  if (neg.to.zero) {
    if (P < 0) P <- 0
    if (ci[1] < 0) ci[1] <- 0
  }
  estimate <- P
  names(estimate) <- "adjusted prevalence"
  data <- paste0("Unadjusted prevalence=", signif(AP, 4), "(se=", signif(stdErrPrev, 4), ")")
  statistic <- Se
  names(statistic) <- paste0("Sensitivity (using nSe=", nSe, ")")
  parameter <- Sp
  names(parameter) <- paste0("Specificity (using nSp=", nSp, ")")
  method <- "Prevalence Adjusted for Sensitivity and Specificity (CI by Korn-Graubard with melding)"
  if (neg.to.zero) method <- "Prevalence Adjusted for Sensitivity and Specificity (CI by Korn-Graubard with melding with negatives set to zero)"
  output <- list(estimate = estimate, statistic = statistic, parameter = parameter, conf.int = ci, data.name = data, method = method, nPeff = nP)
  class(output) <- "htest"
  output
}

WprevSeSp_original <- function(AP, stdErrPrev, nP, Se, nSe, Sp, nSp, conf.level = 0.95,
                      neg.to.zero = TRUE, nmc = 1e5, seed = 49201) {
  # by default set the seed so that we automatically have reproducible research
  # if doing simulations, then set seed to NULL
  if (!is.null(seed)) set.seed(seed)

  P <- prevAdj(AP, Se, Sp)
  # use Korn and Graubard, 1998, Survey Methodology 24(2): 193-201
  # (note they used F distributions, but this beta distributional representation is equivalent)
  # Note for the exact 100(1- alpha)% Clopper-Pearons CI for binomial (as in binom.test)
  # for x ~ binomial(n, theta)
  # is
  # lower = qbeta(alpha/2, x, n-x+1)
  # upper = qbeta(1-alpha/2, x+1, n-x)
  # So Korn-Graubard say
  #  when x=0, use the binomial (ignoring weights)
  #  otherwise use
  #  nstar = AP*(1-AP)/stderr^2
  #  where AP=apparent prevalence
  # then use xstar = AP*nstar
  # then replace nstar for n and xstar for x in the exact binomial
  # below we give the lower=L and upper=U parameter lists for the beta distributions
  if (AP == 0) {
    L <- list(a = 0, b = nP + 1)
    U <- list(a = 1, b = nP)
  } else {
    nstar <- AP * (1 - AP) / (stdErrPrev^2)
    L <- list(a = AP * nstar, b = nstar - AP * nstar + 1)
    U <- list(a = AP * nstar + 1, b = nstar - AP * nstar)
  }

  # use Beta lower and upper confidence distributions
  # for apparent prevalence
  # defined using Korn-Graubard methods
  B_AP_L <- rbeta(nmc, L$a, L$b)
  B_AP_U <- rbeta(nmc, U$a, U$b)

  # use Beta lower and upper confidence distributions
  # for sensitivity (from binomial assumptions)
  B_sen_L <- rbeta(nmc, nSe * Se, nSe - nSe * Se + 1)
  B_sen_U <- rbeta(nmc, nSe * Se + 1, nSe - nSe * Se)

  # use Beta lower and upper confidence distributions
  # for specificity (from binomial assumptions)
  B_spec_L <- rbeta(nmc, nSp * Sp, nSp - nSp * Sp + 1)
  B_spec_U <- rbeta(nmc, nSp * Sp + 1, nSp - nSp * Sp)

  # Monte Carlo melding method
  lo <- prevAdj(B_AP_L, B_sen_U, B_spec_L)
  up <- prevAdj(B_AP_U, B_sen_L, B_spec_U)
  alpha <- 1 - conf.level
  LCL <- quantile(lo, prob = alpha / 2)
  UCL <- quantile(up, prob = 1 - alpha / 2)

  ci <- c(LCL, UCL)
  attr(ci, "conf.level") <- conf.level
  # when neg.to.zero=TRUE, set negative values to zero
  if (neg.to.zero) {
    if (P < 0) P <- 0
    if (ci[1] < 0) ci[1] <- 0
  }
  estimate <- P
  names(estimate) <- "adjusted prevalence"
  data <- paste0("Unadjusted prevalence=", signif(AP, 4), "(se=", signif(stdErrPrev, 4), ")")
  statistic <- Se
  names(statistic) <- paste0("Sensitivity (using nSe=", nSe, ")")
  parameter <- Sp
  names(parameter) <- paste0("Specificity (using nSp=", nSp, ")")
  method <- "Prevalence Adjusted for Sensitivity and Specificity (CI by Korn-Graubard with melding)"
  if (neg.to.zero) method <- "Prevalence Adjusted for Sensitivity and Specificity (CI by Korn-Graubard with melding with negatives set to zero)"
  output <- list(estimate = estimate, statistic = statistic, parameter = parameter, conf.int = ci, data.name = data, method = method, nPeff = nP)
  class(output) <- "htest"
  output
}
