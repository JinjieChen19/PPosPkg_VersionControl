#' Simulate Multiple Bayesian Adaptive Trials allowing for multiple interims, and Summarize Results
#'
#' This function simulates a Bayesian adaptive trial for a specified number of iterations
#' (e.g., 1000) and summarizes the results, including counts and percentages of each
#' decision outcome (e.g., "early go", "early no-go", "final go", "final no-go") and
#' the expected sample size across all trials.
#'
#' @param n_sim Number of simulations to perform (default: 1000).
#' @param alpha_prior Shape parameter 1 of the Beta prior for the response rate.
#' @param beta_prior Shape parameter 2 of the Beta prior for the response rate.
#' @param N Total sample size for the trial.
#' @param n_interim A vector specifying the sample sizes for each interim stage.
#'   The sum of this vector must not exceed `N`.
#' @param go_cutoffs A vector of go decision thresholds for each stage (excluding the final stage).
#'   Must be non-increasing.
#' @param no_go_cutoffs A vector of no-go decision thresholds for each stage (excluding the final stage).
#'   Must be non-decreasing.
#' @param cv1 Benchmark for the first criterion (e.g., minimum acceptable response rate).
#' @param cv2 Benchmark for the second criterion (e.g., desirable response rate for success).
#' @param sig_level Significance level for the first criterion (e.g., probability of success >= `cv1`).
#' @param orr True overall response rate assumed for the simulation.
#' @param n_final_iterations Number of Monte Carlo samples to draw from the posterior distribution
#'   for the final stage. Default is 10000.
#' @param seed An optional random seed for reproducibility. Default is 12345.
#'
#' @return A list containing:
#'   \itemize{
#'     \item \code{summary}: A data frame summarizing the results, with columns:
#'       \itemize{
#'         \item \code{Decision}: The decision type ("early go", "early no-go", "final go", "final no-go").
#'         \item \code{Count}: The number of times each decision occurred.
#'         \item \code{Percentage}: The percentage of simulations resulting in each decision.
#'         \item \code{Expected_Sample_Size}: The mean sample size used across all trials.
#'       }
#'     \item \code{results}: A list of data frames, each representing the detailed results of a single trial.
#'   }
#' @importFrom purrr  map map_dfr
#' @importFrom dplyr group_by summarise n mutate
#' @export sim_multi_IAs_1arm_bin_oc
#'
#' @examples
#' sim <- sim_multi_IAs_1arm_bin_oc(
#'   n_sim = 10,
#'   alpha_prior = 1,
#'   beta_prior = 1,
#'   N = 40,
#'   n_interim = c(10, 10, 10),
#'   go_cutoffs = c(0.95, 0.90, 0.85),
#'   no_go_cutoffs = c(0.04, 0.07, 0.10),
#'   cv1 = 0.4,
#'   cv2 = 0.6,
#'   sig_level = 0.9,
#'   orr = 0.4,
#'   n_final_iterations = 10000,
#'   seed = 12345
#' )
#' print(sim$summary)
sim_multi_IAs_1arm_bin_oc <- function(n_sim = 100, alpha_prior, beta_prior,
                                     N, n_interim, go_cutoffs, no_go_cutoffs,
                                     cv1, cv2, sig_level, orr,
                                     n_final_iterations = 10000, seed = 12345) {

  # Set the random seed for reproducibility
  set.seed(seed)

  # Simulate trials and collect results
  all_results <- map(1:n_sim, ~{
    sim_multi_IAs_1arm_bin_1run(
      alpha_prior = alpha_prior,
      beta_prior = beta_prior,
      N = N,
      n_interim = n_interim,
      go_cutoffs = go_cutoffs,
      no_go_cutoffs = no_go_cutoffs,
      cv1 = cv1,
      cv2 = cv2,
      sig_level = sig_level,
      orr = orr,
      n_final_iterations = n_final_iterations,
      seed = NULL  # Use a random seed for each simulation
    )
  })

  # Extract final decision and sample size from each trial
  final_results <- map_dfr(all_results, ~{
    data.frame(
      Decision = tail(.x$Decision, 1),
      Sample_Size = tail(.x$Cumulative_Sample_Size, 1)
    )
  })

  # Summarize results
  summary <- final_results %>%
    group_by(Decision) %>%
    summarise(
      Count = n(),
      Percentage = (n() / n_sim) * 100,
      .groups = "drop"
    ) %>%
    mutate(Expected_Sample_Size = mean(final_results$Sample_Size))

  return(list(summary = summary, results = all_results))
}


# sim <- sim_multi_IAs_1arm_bin_oc(
#   n_sim = 100,
#   alpha_prior = 1,
#   beta_prior = 1,
#   N = 40,
#   n_interim = c(10, 10, 10),
#   go_cutoffs = c(0.9, 0.9, 0.9),
#   no_go_cutoffs = c(0.1, 0.1, 0.1),
#   cv1 = 0.4,
#   cv2 = 0.6,
#   sig_level = 0.9,
#   orr = 0.4,
#   n_final_iterations = 10000,
#   seed = 42
# )
#
# # Print summary
# print(sim$summary)
#
# # Access detailed results for the first trial
# print(sim$results[[1]])

