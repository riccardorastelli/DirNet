\name{dirnet_posterior}
\alias{dirnet_posterior}
\title{dirnet_posterior}
\description{
Evaluates the log-posterior function for a given set of parameter values. 
}
\usage{
dirnet_posterior(edgelist, model_type, pars, hypers, verbose = F)
}
\arguments{
See \code{\link{dirnet_GS}}.
}
\value{
  \item{computing_time}{Number of seconds required for the evaluation.}
  \item{log_prior_value}{Natural logarithm of the value of the prior distribution.}
  \item{log_likelihood_value}{Natural logarithm of the value of the likelihood function.}
  \item{log_posterior_value}{Natural logarithm of the value of the posterior distribution.}
}
\details{
See \code{\link{dirnet_GS}} for more information and examples.
}
