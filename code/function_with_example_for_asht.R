#' @param method either "binomial" or "poisson"
#' @param x integer vector of apparent positive tests for each group
#' @param n integer vector of number of tests for each group
#' @param w numeric vector of weights for each group (must sum to 1)
#' @param cn number of positive tests for negative controls
#' @param mn number of positive tests results on negative controls
#' @param cp number of positive tests for positive controls
#' @param mp number of positive tests results on positive controls
#' @param conf.level confidence level of the interval
#' @param nmc number of Monte Carlo replications
#' @param seed seed for random number generation
#'
#' @return A list with class "htest" containing the following components:
#' `estimate` the adjusted prevalence estimate, adjusted for sensitivity and specificity
#' `statistic` the estimated sensitivity given by cp / mp
#' `parameter` the estimated specificity given by 1 - cp / mp
#' `conf.int` a confidence interval for the prevalence
#' `data.name` a character string giving the unadjusted prevalence value
#' `method` the character string describing the output.
#' @export
#'
#' @examples
WprevSeSp <- function(method = c("binomial", "poisson"),
                      x, n, w, cn, mn, cp, mp,
                      conf.level = 0.95,
                      nmc = 1e5, seed = 49201) {
  method <- match.arg(method)

  if (length(unique(sapply(list(x, n, w), length))) > 1) {
    stop("x, n, and w must be of equal length")
  }

  if (!isTRUE(all.equal(sum(w), 1))) {
    stop("Weights must sum to 1.")
  }

  if (!is.null(seed)) set.seed(seed)

  g <- function(theta_1, phi_n, phi_p) {
    result <- numeric(length(theta_1))

    ind_for_computed <- which((phi_p >= theta_1) & (theta_1 >= phi_n))
    computed <- (theta_1 - phi_n) / (phi_p - phi_n)

    result[(phi_n < phi_p) & (phi_p < theta_1)] <- 1
    result[ind_for_computed] <- computed[ind_for_computed]
    result
  }

  alpha <- 1 - conf.level
  theta_hat <- x / n
  s_w_theta_hat <- sum(w * theta_hat)

  samples_B_phi_n_L <- rbeta(n = nmc, shape1 = cn, shape2 = mn - cn + 1)
  samples_B_phi_n_U <- rbeta(n = nmc, shape1 = cn + 1, shape2 = mn - cn)

  samples_B_phi_p_L <- rbeta(n = nmc, shape1 = cp, shape2 = mp - cp + 1)
  samples_B_phi_p_U <- rbeta(n = nmc, shape1 = cp + 1, shape2 = mp - cp)

  if (method == "binomial") {
    n_eff <-
      if (s_w_theta_hat == 0) {
      sum(n)
      } else {
        (s_w_theta_hat * (1 - s_w_theta_hat)) / (sum(w^2 / n * theta_hat))
    }
    x_eff <- n_eff * s_w_theta_hat

    samples_B_KG_L <- rbeta(n = nmc, shape1 = x_eff, shape2 = n_eff - x_eff + 1)
    samples_B_KG_U <- rbeta(n = nmc, shape1 = x_eff + 1, shape2 = n_eff - x_eff)

    samples_g_L <- g(samples_B_KG_L, samples_B_phi_n_U, samples_B_phi_p_U)
    samples_g_U <- g(samples_B_KG_U, samples_B_phi_n_L, samples_B_phi_p_L)
  } else if (method == "poisson") {
    y <- sum(w / n * x)
    v <- sum((w / n)^2 * x)
    wm <- max(w / n)
    y_star <- y + wm
    v_star <- v + wm^2

    samples_G_beta_star_L <-
      if (y == 0) {
        numeric(nmc)
      } else {
        rgamma(n = nmc, y^2 / v, scale = v / y)
      }

    samples_G_beta_star_U <- rgamma(n = nmc, y_star^2 / v_star, scale = v_star / y_star)

    samples_g_L <- g(samples_G_beta_star_L, samples_B_phi_n_U, samples_B_phi_p_U)
    samples_g_U <- g(samples_G_beta_star_U, samples_B_phi_n_L, samples_B_phi_p_L)
  }

  confidence_limit_L <- quantile(samples_g_L, prob = alpha / 2)
  confidence_limit_U <- quantile(samples_g_U, prob = 1 - alpha / 2)

  ci <- c(confidence_limit_L, confidence_limit_U)
  AP <- g(s_w_theta_hat, cn / mn, cp / mp)
  Sp <- 1 - cp / mp
  Se <- cp / mp

  estimate <- c("adjusted prevalence" = g(s_w_theta_hat, cn / mn, cp / mp))
  data <- paste0("Unadjusted prevalence=", s_w_theta_hat)
  statistic <- Se
  names(statistic) <- paste0("Sensitivity (using nSe=", mp, ")")
  parameter <- Sp
  names(parameter) <- paste0("Specificity (using nSp=", mp, ")")
  method_text <- paste0(
    "Prevalence Adjusted for Sensitivity and Specificity (CI by ",
    ifelse(method == "binomial", "Korn-Graubard", "Fay and Feuer"),
    " with melding)"
  )
  output <- list(
    estimate = estimate,
    statistic = statistic,
    parameter = parameter,
    conf.int = ci,
    data.name = data,
    method = method_text
  )
  class(output) <- "htest"
  output
}


# Example -----------------------------------------------------------------
example_data_WprevSeSp <- list(
  x = c(
    53, 47, 63, 50, 54, 54, 57, 51, 66, 51, 52, 48, 37, 44, 59,
    55, 50, 58, 52, 54, 41, 45, 49, 54, 37, 53, 57, 58, 55, 55, 56,
    42, 58, 47, 49, 63, 54, 54, 54, 41, 43, 56, 44, 49, 47, 45, 62,
    53, 54, 47
  ),
  n = c(
    200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200,
    200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200,
    200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200,
    200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200
  ),
  w = c(
    0.0205962892844504, 0.0204062236737538, 0.0203843096373626,
    0.0202785701233134, 0.0202617051778543, 0.0202138087214499, 0.0201972974884707,
    0.0201818190015587, 0.0201631543739836, 0.0201560795402158, 0.0201555234250465,
    0.0201461978246263, 0.0201342022821394, 0.0201264004067009, 0.0201167314250592,
    0.0201015081093692, 0.0201003484427457, 0.0201002680000886, 0.0200817537259523,
    0.0200573433887284, 0.0200443907258367, 0.0200358187073312, 0.0200349749335002,
    0.0200264994605187, 0.0200112846914561, 0.020006219121804, 0.0199975642569458,
    0.0199649774153205, 0.0199614929059539, 0.0199426355876479, 0.0199334287088002,
    0.0199298633246975, 0.0199150015155486, 0.0199063452368827, 0.0198920051366782,
    0.0198877425787182, 0.0198679831412633, 0.0198500844815989, 0.0198381388412286,
    0.0198348595904904, 0.0198348180141822, 0.0198174510243331, 0.0197922036364436,
    0.0197821574067888, 0.0197204417557631, 0.0197004976818864, 0.019682896458092,
    0.019649677766428, 0.0196158425485035, 0.019563169292488
  ),
  cn = 77,
  cp = 58,
  mn = 300,
  mp = 60
)


WprevSeSp(
  method = "binomial",
  x = example_data_WprevSeSp$x,
  n = example_data_WprevSeSp$n,
  w = example_data_WprevSeSp$w,
  cn = example_data_WprevSeSp$cn,
  mn = example_data_WprevSeSp$mn,
  cp = example_data_WprevSeSp$cp,
  mp = example_data_WprevSeSp$mp
)

WprevSeSp(
  method = "poisson",
  x = example_data_WprevSeSp$x,
  n = example_data_WprevSeSp$n,
  w = example_data_WprevSeSp$w,
  cn = example_data_WprevSeSp$cn,
  mn = example_data_WprevSeSp$mn,
  cp = example_data_WprevSeSp$cp,
  mp = example_data_WprevSeSp$mp
)

# Check error handling:
WprevSeSp(
  method = "missing method",
  x = example_data_WprevSeSp$x,
  n = example_data_WprevSeSp$n,
  w = example_data_WprevSeSp$w,
  cn = example_data_WprevSeSp$cn,
  mn = example_data_WprevSeSp$mn,
  cp = example_data_WprevSeSp$cp,
  mp = example_data_WprevSeSp$mp
)


WprevSeSp(
  method = "poisson",
  x = example_data_WprevSeSp$x,
  n = 200,
  w = example_data_WprevSeSp$w,
  cn = example_data_WprevSeSp$cn,
  mn = example_data_WprevSeSp$mn,
  cp = example_data_WprevSeSp$cp,
  mp = example_data_WprevSeSp$mp
)

WprevSeSp(
  method = "binomial",
  x = rep(0, length(example_data_WprevSeSp$x)),
  n = example_data_WprevSeSp$n,
  w = example_data_WprevSeSp$w,
  cn = example_data_WprevSeSp$cn,
  mn = example_data_WprevSeSp$mn,
  cp = example_data_WprevSeSp$cp,
  mp = example_data_WprevSeSp$mp
)

WprevSeSp(
  method = "poisson",
  x = rep(0, length(example_data_WprevSeSp$x)),
  n = example_data_WprevSeSp$n,
  w = example_data_WprevSeSp$w,
  cn = example_data_WprevSeSp$cn,
  mn = example_data_WprevSeSp$mn,
  cp = example_data_WprevSeSp$cp,
  mp = example_data_WprevSeSp$mp
)
