



#' Title: Make 2-arm trial decision for Binary Endpoints
#'
#' @param n1 interim sample size for treatment arm
#' @param n2 interim sample size for control arm
#' @param x1 observed number of events in treatment arm
#' @param x2 observed number of events in control arm
#' @param alpha_trt_prior shape parameter 1 for beta prior for treatment arm
#' @param beta_trt_prior  shape parameter 2 for beta prior for treatment arm
#' @param alpha_ctrl_prior shape parameter 1 for beta prior for control arm
#' @param beta_ctrl_prior shape parameter 2 for beta prior for control arm
#' @param N1 final sample size for treatment arm
#' @param N2 final sample size for control arm
#' @param cv1 delta_orr benchmark for criterion 1
#' @param cv2 delta_orr benchmark for criterion 2
#' @param sig_level significance level for dual-criterion
#' @param ppos_go early go cutoff for PPoS
#' @param ppos_no_go early no go cutoff for PPoS
#' @param orr_trt assumed ORR for treatment arm
#' @param orr_ctl assumed ORR for control arm
#' @param n_final_iterations number of samples to draw from the posterior distribution
#'
#' @return a data frame with interim decision, final decision, and total sample size
#' @export make_2arm_trial_decision
#'
#' @examples n1 = 25; n2 = 25; x1 = 10; x2 = 5;
#'           alpha_trt_prior = 1; beta_trt_prior = 1;
#'           alpha_ctrl_prior = 1; beta_ctrl_prior = 1;
#'           N1 = 32; N2 = 32; cv1 = 0.0; cv2 = 0.23;
#'           sig_level = 0.90; ppos_go = 0.90; ppos_no_go = 0.10;
#'           orr_trt = 0.4; orr_ctl = 0.2;
#'           n_final_iterations = 10000
#'           make_2arm_trial_decision(n1, n2, x1, x2,
#'                         alpha_trt_prior, beta_trt_prior,
#'                         alpha_ctrl_prior, beta_ctrl_prior,
#'                         N1, N2, cv1, cv2, sig_level,
#'                         ppos_go, ppos_no_go, orr_trt, orr_ctl,
#'                         n_final_iterations)
make_2arm_trial_decision <- function(n1, n2, x1, x2,
                                     alpha_trt_prior, beta_trt_prior,
                                     alpha_ctrl_prior, beta_ctrl_prior,
                                     N1, N2, cv1, cv2, sig_level=0.90,
                                     ppos_go = 0.90, ppos_no_go = 0.10,
                                     orr_trt = 0.4, orr_ctl = 0.2,
                                     n_final_iterations = 1000
){
  # step 1: calculate ppos at interim using calculate_ppos_2arm_bin() function
  ppos_interim <- calculate_ppos_2arm_bin(n1, n2, x1, x2,
                                          alpha_trt_prior, beta_trt_prior,
                                          alpha_ctrl_prior, beta_ctrl_prior,
                                          N1, N2, cv1, cv2, sig_level)
  #Initialize decisions and final sample size
  interim_decision <- NULL
  final_decision <- NULL
  total_sample_size <- n1 + n2

  #step 2: Interim decision based on ppos_go and ppos_no_go
  if(ppos_interim > ppos_go){
    interim_decision <- "early go"
    final_decision <- "early terminate"
  } else if(ppos_interim < ppos_no_go){
    interim_decision <- "early no-go"
    final_decision <- "early terminate"
  } else {
    #if continue, proceed to final stage
    interim_decision <- "continue"

    #simulate remaining responses based on input ORR and sample size
    x1_remaining <- rbinom(1, N1-n1, orr_trt)
    x2_remaining <- rbinom(1, N2-n2, orr_ctl)

    #update posteriors based on combined interim+final data for each arm
    alpha_trt_posterior <- alpha_trt_prior + x1 + x1_remaining
    beta_trt_posterior <- beta_trt_prior + n1 + N1 - x1 + N1 - x1_remaining
    alpha_ctrl_posterior <- alpha_ctrl_prior + x2 + x2_remaining
    beta_ctrl_posterior <- beta_ctrl_prior + n2 + N2 - x2 + N2 - x2_remaining

    #step 4: generate samples from the posterior distribution
    trt_samples <- rbeta(n_final_iterations, alpha_trt_posterior, beta_trt_posterior)
    ctrl_samples <- rbeta(n_final_iterations, alpha_ctrl_posterior, beta_ctrl_posterior)
    delta_samples <- trt_samples - ctrl_samples

    #step 5: calculate the probability of the treatment being superior
    prob_final_cv1 <- mean(delta_samples > cv1)
    prob_final_cv2 <- mean(delta_samples < cv2)

    #step 6: make final decision based on the probability of meeting dual-criterion
    if(prob_final_cv1 >= sig_level & prob_final_cv2 >= 0.5){
      final_decision <- "final go"
    } else if(prob_final_cv2 < sig_level & prob_final_cv1 < 0.5){
      final_decision <- "final no-go"
    } else{
      final_decision <- "inderterminate"
    }

    total_sample_size <- N1 + N2
  }

  return(data.frame(interim_decision = interim_decision,
                    final_decision = final_decision,
                    total_sample_size = total_sample_size))
}
