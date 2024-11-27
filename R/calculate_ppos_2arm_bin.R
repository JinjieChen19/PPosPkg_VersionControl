#' Function to calculate the posterior probability of success for a two-arm trial with binary endpoints based on Bayesian go no-go dual-criterion design
#'
#' @param n1 Number of patients in treatment arm
#' @param n2 Number of patients in control arm
#' @param x1 Number of responses in treatment arm
#' @param x2 Number of responses in control arm
#' @param alpha_trt_prior shape parameter 1 of Beta prior for treatment arm
#' @param beta_trt_prior shape parameter 2 of Beta prior for treatment arm
#' @param alpha_ctrl_prior shape parameter 1 of Beta prior for control arm
#' @param beta_ctrl_prior shape parameter 2 of Beta prior for control arm
#' @param N1 Total number of patients in treatment arm
#' @param N2 Total number of patients in control arm
#' @param cv1 Benchmark for the first criterion
#' @param cv2 Benchmark for the second criterion
#' @param sig_level Significance level for the first of dual-criterion
#'
#' @return Posterior probability of success, defined as the probability that the dual-criterion is met at final given the interim data.
#' @importFrom VGAM dbetabinom.ab
#' @importFrom stats rbeta
#' @export calculate_ppos_2arm_bin
#'
#' @examples x1 = 16; x2 = 9; n1 = 25;
#'           n2 = 25; N1 = 32; N2 = 32;
#'           alpha_trt_prior = 1; beta_trt_prior = 1;
#'           alpha_ctrl_prior = 1; beta_ctrl_prior = 1;
#'           cv1 = 0.0; cv2 = 0.4; sig_level = 0.9
#'           calculate_ppos_2arm_bin(
#'             n1, n2, x1, x2, alpha_trt_prior, beta_trt_prior,
#'             alpha_ctrl_prior, beta_ctrl_prior, N1, N2, cv1, cv2,
#'             sig_level
#'           )
calculate_ppos_2arm_bin <- function(n1, n2, x1, x2, alpha_trt_prior, beta_trt_prior, alpha_ctrl_prior, beta_ctrl_prior, N1, N2, cv1, cv2, sig_level) {
  #Calculate posterior probability of success in each arm
  alpha1 <- alpha_trt_prior + x1
  beta1 <- beta_trt_prior + n1 - x1
  alpha2 <- alpha_ctrl_prior + x2
  beta2 <- beta_ctrl_prior + n2 - x2

  # Initial PPoS
  ppos <- 0

  #Enumerate all future possible outcomes
  for (x1_future in 0: (N1-n1)) {
    for (x2_future in 0: (N2-n2)){
      # update overall response counts
      final_x1 <- x1 + x1_future
      final_x2 <- x2 + x2_future
      # update posterior parameters based on combined interim and future data
      alpha_trt_post <- alpha_trt_prior + final_x1
      beta_trt_post <- beta_trt_prior + N1 - final_x1
      alpha_ctrl_post <- alpha_ctrl_prior + final_x2
      beta_ctrl_post <- beta_ctrl_prior + N2 - final_x2

      #Monte Carlo sampling from updated beta distributions
      n_sims <- 10000
      trt_sims <- rbeta(n_sims, alpha_trt_post, beta_trt_post)
      ctrl_sims <- rbeta(n_sims, alpha_ctrl_post, beta_ctrl_post)
      delta_sims <- trt_sims - ctrl_sims

      #calculate probability of success in each arm
      ppos_trt <- mean(delta_sims > cv1)
      ppos_ctrl <- mean(delta_sims > cv2)

      #check if dual criterion is met
      if (ppos_trt >= sig_level & ppos_ctrl >= 0.5) {
        #compute probability of this outcome under the beta-binomial model
        p_future_x1 <- VGAM::dbetabinom.ab(x1_future, N1-n1, alpha1, beta1)
        p_future_x2 <- VGAM::dbetabinom.ab(x2_future, N2-n2, alpha2, beta2)
        ppos <- ppos + p_future_x1 * p_future_x2
      }

    }

  }

  #Calculate posterior probability of success in each arm
  return(ppos)
}


