#' Simulate an Adaptive Bayesian Trial with Multiple Interim Analyses
#'
#' This function simulates a Bayesian adaptive trial with multiple interim analyses
#' for a single-arm study with a binary endpoint. At each interim stage, the function
#' evaluates the posterior predictive probability of success (PPoS) based on cumulative
#' data and predefined decision thresholds (go/no-go cutoffs). If no early stopping
#' decision is made, the trial proceeds to the next stage or the final analysis.
#'
#' @param alpha_prior Shape parameter 1 of the Beta prior for the response rate.
#' @param beta_prior Shape parameter 2 of the Beta prior for the response rate.
#' @param N Total sample size for the trial.
#' @param n_interim A vector specifying the sample sizes for each interim stage.
#'   The sum of this vector must not exceed `N`.
#' @param go_cutoffs A vector of go decision thresholds for each stage (excluding the final stage).
#'   Must be non-increasing.
#' @param no_go_cutoffs A vector of no-go decision thresholds for each stage (excluding the final stage).
#'   Must be non-decreasing.
#' @param cv1 Benchmark for the first criterion (e.g., lower reference response rate).
#' @param cv2 Benchmark for the second criterion (e.g., minimal clinical meaningful response rate for success).
#' @param sig_level Significance level for the first criterion (e.g., probability of success >= `cv1`).
#' @param orr True overall response rate assumed for the simulation.
#' @param n_final_iterations Number of Monte Carlo samples to draw from the posterior distribution
#'   for the final stage. Default is 10000.
#' @param seed An optional random seed for reproducibility. Default is 12345.
#'
#' @return A data frame with the following columns:
#'   \itemize{
#'     \item \code{Stage}: The stage number (interim stages and final stage).
#'     \item \code{Stage_Sample_Size}: The sample size for the current stage.
#'     \item \code{Cumulative_Sample_Size}: The cumulative sample size up to the current stage.
#'     \item \code{Cumulative_Responses}: The cumulative number of responses up to the current stage.
#'     \item \code{PPoS}: The posterior predictive probability of success at interim stages (NA for the final stage).
#'     \item \code{Prob_Meets_CV1}: Probability of meeting the first criterion (only for the final stage).
#'     \item \code{Prob_Meets_CV2}: Probability of meeting the second criterion (only for the final stage).
#'     \item \code{Decision}: The decision made at the current stage ("early go", "early no-go", "continue",
#'       or "final go"/"final no-go").
#'   }
#'
#' @export sim_multi_IAs_1arm_bin_1run
#'
#' @examples
#' result <- sim_multi_IAs_1arm_bin_1run(
#'   alpha_prior = 1,
#'   beta_prior = 1,
#'   N = 40,
#'   n_interim = c(10, 10, 10),
#'   go_cutoffs = c(0.9, 0.8, 0.7),
#'   no_go_cutoffs = c(0.1, 0.2, 0.3),
#'   cv1 = 0.4,
#'   cv2 = 0.6,
#'   sig_level = 0.9,
#'   orr = 0.6,
#'   n_final_iterations = 10000,
#'   seed = sample(1e6,1)
#' )
#' print(result)
sim_multi_IAs_1arm_bin_1run <- function(alpha_prior, beta_prior,
                                    N, n_interim,
                                    go_cutoffs, no_go_cutoffs,
                                    cv1, cv2, sig_level,
                                    orr, n_final_iterations = 10000, seed = 12345) {
  # Set random seed if provided
  if (!is.null(seed)) set.seed(seed)

  # Ensure the sum of n_interim does not exceed the total sample size N
  if (sum(n_interim) > N) stop("The total interim sample size cannot exceed the total sample size N.")

  # Ensure go_cutoffs are non-increasing and no_go_cutoffs are non-decreasing
  if (!all(diff(go_cutoffs) <= 0)) stop("go_cutoffs must be non-increasing.")
  if (!all(diff(no_go_cutoffs) >= 0)) stop("no_go_cutoffs must be non-decreasing.")

  # Add final stage to n_interim (remaining sample size)
  n_interim <- c(n_interim, N - sum(n_interim))
  total_stages <- length(n_interim)
  go_cutoffs <- c(go_cutoffs, NA)  # No go_cutoff for final stage
  no_go_cutoffs <- c(no_go_cutoffs, NA)  # No no_go_cutoff for final stage

  # Initialize variables
  cumulative_n <- 0
  cumulative_x <- 0
  decisions <- list()  # Save the decision results for each stage

  # Loop through each stage
  for (stage in 1:total_stages) {
    n_stage <- n_interim[stage]
    cumulative_n <- cumulative_n + n_stage

    # Generate response data for the current stage
    x_stage <- rbinom(1, n_stage, orr)
    cumulative_x <- cumulative_x + x_stage

    if (stage < total_stages) {
      # Interim stage: Compute PPoS
      ppos <- calculate_ppos_single_arm_bin(
        x = cumulative_x,
        n = cumulative_n,
        alpha_prior = alpha_prior, # Always use the initial prior as we use cumulative data
        beta_prior = beta_prior,   # Always use the initial prior as we use cumulative data
        N = N,
        cv1 = cv1,
        cv2 = cv2,
        sig_level = sig_level
      )

      # Get stage-specific cutoffs
      go_cutoff <- go_cutoffs[stage]
      no_go_cutoff <- no_go_cutoffs[stage]

      # Record the stage
      decisions[[stage]] <- list(
        stage = stage,
        n_stage = n_stage,
        cumulative_n = cumulative_n,
        cumulative_x = cumulative_x,
        ppos = ppos
      )

      # Make decision using stage-specific cutoffs
      if (ppos > go_cutoff) {
        decisions[[stage]]$decision <- "early go"
        break
      } else if (ppos < no_go_cutoff) {
        decisions[[stage]]$decision <- "early no-go"
        break
      } else {
        decisions[[stage]]$decision <- "continue"
      }
    } else {
      # Final stage: Simulate remaining responses and update posterior
      x_remaining <- rbinom(1, N - cumulative_n, orr)
      # update posterior parameters based on all data
      alpha_post <- alpha_prior + cumulative_x + x_remaining
      beta_post <- beta_prior + N - cumulative_n - x_remaining

      # Monte Carlo sampling from the updated posterior distribution at final
      posterior_samples <- rbeta(n_final_iterations, alpha_post, beta_post)

      # Calculate probabilities of meeting the criteria
      prob_meets_cv1 <- mean(posterior_samples > cv1)
      prob_meets_cv2 <- mean(posterior_samples > cv2)

      # Record the final stage
      decisions[[stage]] <- list(
        stage = stage,
        n_stage = n_stage,
        cumulative_n = cumulative_n,
        cumulative_x = cumulative_x,
        ppos = NA,  # PPoS is not applicable at the final stage
        prob_meets_cv1 = prob_meets_cv1,
        prob_meets_cv2 = prob_meets_cv2
      )

      # Make final decision
      if (prob_meets_cv1 >= sig_level && prob_meets_cv2 >= 0.5) {
        decisions[[stage]]$decision <- "final go"
      } else if (prob_meets_cv2 < sig_level && prob_meets_cv1 < 0.5) {
        decisions[[stage]]$decision <- "final no-go"
      } else {
        decisions[[stage]]$decision <- "indeterminate"
      }
    }
  }

  # Transform the list to a data frame
  result <- do.call(rbind, lapply(decisions, function(x) data.frame(
    Stage = x$stage,
    Stage_Sample_Size = x$n_stage,
    Cumulative_Sample_Size = x$cumulative_n,
    Cumulative_Responses = x$cumulative_x,
    PPoS = x$ppos,
    Prob_Meets_CV1 = ifelse(is.null(x$prob_meets_cv1), NA, x$prob_meets_cv1),
    Prob_Meets_CV2 = ifelse(is.null(x$prob_meets_cv2), NA, x$prob_meets_cv2),
    Decision = x$decision
  )))

  return(result)
}
