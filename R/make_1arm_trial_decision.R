


#' Decision-making function for single arm binary endpoint trial with single interim analysis
#' based on Bayesian go no-go dual-criterion design per posteriors probability of success
#' given interim data
#'
#' @param x interim number of responses
#' @param n interim sample size
#' @param alpha_prior shape parameter 1 of Beta prior
#' @param beta_prior shape parameter 2 of Beta prior
#' @param N final sample size
#' @param cv1 benchmark for  criterion 1
#' @param cv2 benchmark for  criterion 2
#' @param sig_level significance level for criterion 1
#' @param ppos_go cutoff for PPoS to early go at interim
#' @param ppos_no_go cutoff for PPoS to early no-go at interim
#' @param orr assumed overall response rate
#' @param n_final_iterations number of samples to draw from the posterior distribution
#'
#' @return A data frame with the following columns:
#'   \describe{
#'     \item{ppos_interim}{Posterior Predictive Probability of Success at interim.}
#'     \item{interim_decision}{Decision at the interim analysis: "early go", "early no-go", or "continue".}
#'     \item{final_decision}{Decision at the final stage: "early go", "early no-go", or "indeterminate".}
#'     \item{total_sample_size}{The total sample size used in the trial.}
#'   }
#' @export make_1arm_trial_decision
#' @importFrom stats rbeta rbinom
#'
#' @examples   x = 13; n = 16; alpha_prior = 1; beta_prior = 1;
#'             N = 26; cv1 = 0.4; cv2 = 0.6; sig_level = 0.95
#'             ppos_go = 0.90; ppos_no_go = 0.10; orr = 0.4;
#'             n_final_iterations = 10000
#'             make_1arm_trial_decision(x, n, alpha_prior, beta_prior,
#'                        N, cv1, cv2, sig_level,
#'                         ppos_go, ppos_no_go, orr,
#'                         n_final_iterations)
make_1arm_trial_decision <- function(x, n, alpha_prior, beta_prior,
                                     N, cv1, cv2, sig_level,
                                     ppos_go, ppos_no_go, orr,
                                     n_final_iterations = 10000) {

  # Input validation
  if (x > n || n > N) stop("x must be <= n and n must be <= N.")
  if (!all(c(ppos_go, ppos_no_go, cv1, cv2, sig_level) >= 0 & c(ppos_go, ppos_no_go, cv1, cv2, sig_level) <= 1)) {
    stop("All probability parameters (ppos_go, ppos_no_go, cv1, cv2, sig_level) must be between 0 and 1.")
  }

  # Step 1: Calculate PPoS at interim
  ppos_interim <- calculate_ppos_single_arm_bin(x, n, alpha_prior, beta_prior, N, cv1, cv2, sig_level)

  # Initialize decisions and total sample size
  interim_decision <- NULL
  final_decision <- NULL
  total_sample_size <- n

  # Step 2: Interim decision based on PPoS thresholds
  if (ppos_interim > ppos_go) {
    interim_decision <- "early go"
    final_decision <- "go"
  } else if (ppos_interim < ppos_no_go) {
    interim_decision <- "early no-go"
    final_decision <- "no-go"
  } else {
    # Continue to final stage
    interim_decision <- "continue"

    # Step 3: Simulate remaining responses
    x_remaining <- rbinom(1, N - n, orr)

    # Update posterior parameters
    alpha_posterior <- alpha_prior + (x + x_remaining)
    beta_posterior <- beta_prior + N - (x + x_remaining)

    # Step 4: Generate posterior samples
    samples <- rbeta(n_final_iterations, alpha_posterior, beta_posterior)
    prob_meets_cv1 <- mean(samples > cv1)
    prob_meets_cv2 <- mean(samples > cv2)

    # Step 5: Final decision based on posterior probabilities
    if (prob_meets_cv1 >= sig_level && prob_meets_cv2 >= 0.5) {
      final_decision <- "go"
    } else if (prob_meets_cv1 < sig_level && prob_meets_cv2 < 0.5) {
      final_decision <- "no-go"
    } else {
      final_decision <- "indeterminate"
    }

    # Update total sample size
    total_sample_size <- N
  }

  # Return decisions and total sample size
  return(data.frame(ppos_interim = ppos_interim,
                    interim_decision = interim_decision,
                    final_decision = final_decision,
                    total_sample_size = total_sample_size))
}


