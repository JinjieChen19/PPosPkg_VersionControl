
#' Title: Simulate 2-arm trial with binary endpoints per PPoS and GNG dual-criterion design
#'
#' @param n_simulation number of simulations
#' @param n1 interim sample size for treatment arm
#' @param n2 interim sample size for control arm
#' @param N1 final sample size for treatment arm
#' @param N2 final sample size for control arm
#' @param orr_trt assumed ORR for treatment arm
#' @param orr_ctl assumed ORR for control arm
#' @param ppos_go cutoff for PPoS to early go at interim
#' @param ppos_no_go cutoff for PPoS to early no-go at interim
#' @param sig_level significance level for  criterion 1
#' @param cv1 benchmark for response difference criterion 1
#' @param cv2 benchmark for response difference criterion 2
#' @param alpha_trt_prior shape parameter 1 for beta prior for treatment arm
#' @param beta_trt_prior shape parameter 2 for beta prior for treatment arm
#' @param alpha_ctrl_prior shape parameter 1 for beta prior for control arm
#' @param beta_ctrl_prior shape parameter 2 for beta prior for control arm
#' @param seed seed for reproducibility
#'
#' @return a data frame with interim decision, final decision, and expected sample size
#' @export simulate_OC_2arm_binary
#' @importFrom purrr map_dfr
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @importFrom dplyr n
#' @examples n_simulation = 10; n1 = 25; n2 = 25;
#'           N1 = 32; N2 = 32; orr_trt = 0.4; orr_ctl = 0.2;
#'           ppos_go = 0.90; ppos_no_go = 0.10; sig_level = 0.90;
#'           cv1 = 0.0; cv2 = 0.23; alpha_trt_prior = 1; beta_trt_prior = 1;
#'           alpha_ctrl_prior = 1; beta_ctrl_prior = 1;
#'           seed = 112024
#'           simulate_OC_2arm_binary(n_simulation, n1, n2, N1, N2, orr_trt, orr_ctl,
#'                        ppos_go, ppos_no_go, sig_level, cv1, cv2,
#'                        alpha_trt_prior, beta_trt_prior, alpha_ctrl_prior, beta_ctrl_prior,
#'                        seed)
simulate_OC_2arm_binary <- function(n_simulation = 1000, n1 = 25, n2 = 25,
                                    N1 = 32, N2 = 32, orr_trt = 0.4, orr_ctl = 0.2,
                                    ppos_go = 0.90, ppos_no_go = 0.10, sig_level = 0.90,
                                    cv1 = 0.0, cv2 = 0.23, alpha_trt_prior = 1, beta_trt_prior = 1,
                                    alpha_ctrl_prior = 1, beta_ctrl_prior = 1,
                                    seed = 112024) {
        #Initialize vectors to store results
  interim_decision <- rep(NA, n_simulation)
  final_decision <- rep(NA, n_simulation)
  sample_size <- rep(NA, n_simulation)
  overall_decisions <- rep(NA, n_simulation)
  # set seed for reproducibility
  set.seed(seed)
  #wrapper function to run one simulation
  one_simulation <- function(){
    x1 <- rbinom(1, n1, orr_trt)
    x2 <- rbinom(1, n2, orr_ctl)
    result <- make_2arm_trial_decision(n1, n2, x1, x2,
                                       alpha_trt_prior, beta_trt_prior,
                                       alpha_ctrl_prior, beta_ctrl_prior,
                                       N1, N2, cv1, cv2, sig_level,
                                       ppos_go, ppos_no_go, orr_trt, orr_ctl)
    return(result)
  }

  #run multiple simulations using map_dfr
  results <- purrr::map_dfr(1:n_simulation, ~one_simulation())
  #extract results
  results_df <- data.frame(interim_decision = results$interim_decision,
                           final_decision = results$final_decision,
                           sample_size = results$total_sample_size)
  # add a column to indicate overall decision
  results_df$overall_decision <- ifelse(results_df$interim_decision == "early go" | results_df$final_decision == "final go", "overall go",
                                        ifelse(results_df$final_decision == "early no-go" | results_df$final_decision == "final no-go", "overll no-go", "indeterminate"))
  #calculate the summary table, including interim and final deciions count and percentage, as well as expected sample size
  summary_table <- results_df %>%
    group_by(interim_decision, final_decision,overall_decision) %>%
    summarise(n = n(), percentage = n()/n_simulation)
  summary_table$expected_sample_size <- mean(results_df$sample_size)

  return(summary_table)
}

