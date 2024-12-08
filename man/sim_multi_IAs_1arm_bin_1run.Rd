% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate_adaptive_trial_single_arm_bin.R
\name{sim_multi_IAs_1arm_bin_1run}
\alias{sim_multi_IAs_1arm_bin_1run}
\title{Simulate an Adaptive Bayesian Trial with Multiple Interim Analyses}
\usage{
sim_multi_IAs_1arm_bin_1run(
  alpha_prior,
  beta_prior,
  N,
  n_interim,
  go_cutoffs,
  no_go_cutoffs,
  cv1,
  cv2,
  sig_level,
  orr,
  n_final_iterations = 10000,
  seed = 12345
)
}
\arguments{
\item{alpha_prior}{Shape parameter 1 of the Beta prior for the response rate.}

\item{beta_prior}{Shape parameter 2 of the Beta prior for the response rate.}

\item{N}{Total sample size for the trial.}

\item{n_interim}{A vector specifying the sample sizes for each interim stage.
The sum of this vector must not exceed \code{N}.}

\item{go_cutoffs}{A vector of go decision thresholds for each stage (excluding the final stage).
Must be non-increasing.}

\item{no_go_cutoffs}{A vector of no-go decision thresholds for each stage (excluding the final stage).
Must be non-decreasing.}

\item{cv1}{Benchmark for the first criterion (e.g., lower reference response rate).}

\item{cv2}{Benchmark for the second criterion (e.g., minimal clinical meaningful response rate for success).}

\item{sig_level}{Significance level for the first criterion (e.g., probability of success >= \code{cv1}).}

\item{orr}{True overall response rate assumed for the simulation.}

\item{n_final_iterations}{Number of Monte Carlo samples to draw from the posterior distribution
for the final stage. Default is 10000.}

\item{seed}{An optional random seed for reproducibility. Default is 12345.}
}
\value{
A data frame with the following columns:
\itemize{
\item \code{Stage}: The stage number (interim stages and final stage).
\item \code{Stage_Sample_Size}: The sample size for the current stage.
\item \code{Cumulative_Sample_Size}: The cumulative sample size up to the current stage.
\item \code{Cumulative_Responses}: The cumulative number of responses up to the current stage.
\item \code{PPoS}: The posterior predictive probability of success at interim stages (NA for the final stage).
\item \code{Prob_Meets_CV1}: Probability of meeting the first criterion (only for the final stage).
\item \code{Prob_Meets_CV2}: Probability of meeting the second criterion (only for the final stage).
\item \code{Decision}: The decision made at the current stage ("early go", "early no-go", "continue",
or "final go"/"final no-go").
}
}
\description{
This function simulates a Bayesian adaptive trial with multiple interim analyses
for a single-arm study with a binary endpoint. At each interim stage, the function
evaluates the posterior predictive probability of success (PPoS) based on cumulative
data and predefined decision thresholds (go/no-go cutoffs). If no early stopping
decision is made, the trial proceeds to the next stage or the final analysis.
}
\examples{
result <- sim_multi_IAs_1arm_bin_1run(
  alpha_prior = 1,
  beta_prior = 1,
  N = 40,
  n_interim = c(10, 10, 10),
  go_cutoffs = c(0.9, 0.8, 0.7),
  no_go_cutoffs = c(0.1, 0.2, 0.3),
  cv1 = 0.4,
  cv2 = 0.6,
  sig_level = 0.9,
  orr = 0.6,
  n_final_iterations = 10000,
  seed = sample(1e6,1)
)
print(result)
}
