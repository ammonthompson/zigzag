% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/zigzagBurnin.R
\name{burnin}
\alias{burnin}
\title{burnin}
\usage{
burnin(sample_frequency = 100, ngen = 1000, burnin_target_acceptance_rate=0.44,
progress_plot = FALSE, write_to_files = TRUE, burninprefix = "burnin", append = FALSE)
}
\arguments{
\item{sample_frequency}{Number of generations between samples from the chain}

\item{ngen}{Number of generations to run the chain}

\item{burnin_target_acceptance_rate}{proportion of proposals that are accepted. For proposal size tuning.}

\item{progress_plot}{Show plots of model statistics as MCMC progresses}

\item{write_to_files}{Write burnin samples to output files in burninprefix_burnin_output directory}

\item{burninprefix}{The prefix for the burnin output directory as well as the burnin output files within}

\item{append}{If files already exist from previous interrupted burnin, append to those files if TRUE.}
}
\description{
Run burnin and tune proposal size parameters for hieararchical bayesian mixture model
}
