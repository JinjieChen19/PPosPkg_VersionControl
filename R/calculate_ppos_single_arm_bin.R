#' Function to calculate the posterior probability of success for a single-arm trial with binary endpoint
#' based on Bayesian go no-go dual-criterion design
#'
#' @param x interim number of responses
#' @param n interim sample size
#' @param alpha_prior shape parameter 1 of Beta prior
#' @param beta_prior shape parameter 2 of Beta prior
#' @param N final sample size
#' @param cv1 benchmark for  criterion 1
#' @param cv2 benchmark for  criterion 2
#' @param sig_level significance level for criterion 1
#' @param n_sims iterations for Monte Carlo simulation
#' @importFrom VGAM dbetabinom.ab
#'
#' @return Posterior probability of success, defined as the probability that the dual-criterion is met at final given the interim data.
#' @export calculate_ppos_single_arm_bin
#'
#' @examples
#' x <- 13
#' n <- 16
#' alpha_prior <- 1
#' beta_prior <- 1
#' N <- 26
#' cv1 <- 0.4
#' cv2 <- 0.6
#' sig_level <- 0.95
#' n_sims <- 10000
#' seed <- 42
#' ppos <- calculate_ppos_single_arm_bin(x, n, alpha_prior, beta_prior, N, cv1, cv2, sig_level, n_sims)
#' print(ppos)
calculate_ppos_single_arm_bin <- function(x = 13, n = 16, alpha_prior = 1, beta_prior = 1,
                                          N = 26, cv1 = 0.4, cv2 = 0.6,
                                          sig_level = 0.95, n_sims = 10000) {

  # Validate inputs
  if (n > N) stop("Interim sample size (n) must be less than or equal to final sample size (N).")
  if (!all(c(cv1, cv2, sig_level) >= 0 & c(cv1, cv2, sig_level) <= 1)) {
    stop("Probabilities (cv1, cv2, sig_level) must be between 0 and 1.")
  }

  # Update posterior parameters based on interim data
  alpha_post <- alpha_prior + x
  beta_post <- beta_prior + n - x

  # Enumerate all possible future outcomes and their probabilities
  x_future <- 0:(N - n)
  p_future_x <- VGAM::dbetabinom.ab(x_future, N - n, alpha_post, beta_post)

  # Simulate posterior probabilities for each possible future outcome
  final_x <- x + x_future
  alpha_post_final <- alpha_prior + final_x
  beta_post_final <- beta_prior + N - final_x

  # Vectorized Monte Carlo sampling
  sims <- mapply(function(alpha, beta) rbeta(n_sims, alpha, beta),
                 alpha = alpha_post_final, beta = beta_post_final, SIMPLIFY = FALSE)

  # Calculate probabilities of meeting dual criteria for each outcome
  prob_meets_cv1 <- sapply(sims, function(s) mean(s > cv1))
  prob_meets_cv2 <- sapply(sims, function(s) mean(s > cv2))
  meets_criteria <- (prob_meets_cv1 >= sig_level) & (prob_meets_cv2 >= 0.5)

  # Calculate PPoS
  ppos <- sum(p_future_x[meets_criteria])

  return(ppos)
}
