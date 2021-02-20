\name{dirnet_GS}
\alias{dirnet_GS}
\title{dirnet_GS}
\description{
Runs the Metropolis-within-Gibbs sampler to fit the statistical model. 
}
\usage{
dirnet_GS(edgelist, model_type, pars, hypers, mcmc, proposals, verbose = F)
}
\arguments{
  \item{edgelist}{A matrix with 4 columns. In the generic row \code{(t, i, j, v)}: \code{t} is the time-frame number; \code{i} is the sender node; \code{j} is the receiver node; \code{v} is the natural logarithm of a proportion. This proportion corresponds to the relative weight carried by the edge \code{(i,j)} at time \code{t}. The edgelist must not contain self-edges, and the out-going edge proportions are assumed to sum to one. The number of nodes and number of time-frames in the network are deduced from this edgelist.}
  \item{model_type}{A three-dimensional binary vector identifying the model that should be fitted to the data. The first number determines whether the trend component should be included, the second number determines whether a sender random effect should be included, and the third number determines whether a receiver random effect should be included.}
  \item{pars}{A list indicating the initial values for the parameters. It must have the following elements: \code{mu} for the trend effect; \code{theta} for the sender effect; \code{gamma} for the receiver effect; \code{tau_eta} for \code{mu}'s precision; \code{tau_theta} for \code{theta}'s precision; \code{tau_gamma} for \code{gamma}'s precision. Values should be provided regardless of whether the effects are included in the model or not. If some of the parameters are not included, they should be initialized as zeros. The dimensionalities of each of the objects should match the dimensions derived from the edgelist.}
  \item{hypers}{A list indicating the hyper-parameters for the prior distributions. It must have the following elements: \code{sigma_mu} for the variance of \code{mu} in the first time-frame; \code{a_eta} and \code{b_eta} for the gamma prior on \code{tau_eta}; \code{a_theta}, \code{b_theta}, \code{a_gamma} and \code{b_gamma} analogously for the other gamma priors.}
  \item{mcmc}{A list setting some options regarding the sampling procedure. It must have the following elements: \code{burnin} the number of iterations to be discarded as burnin; \code{thin} to set how frequently new sampled observations should be retained (i.e. save the values every \code{thin} samples); \code{n_iterations} as the number of observations that will be available in the final sample, for each parameter. The thinning starts after the burn-in period, meaning that the total number of iterations of the algorithm is \code{burnin + thin * n_iterations}.}
  \item{proposals}{A list characterizing the proposal distributions for the parameters. It must have the following elements: \code{mu} a vector indicating the standard deviations of the Gaussian proposals individually for each \code{mu} parameter; \code{theta} and \code{gamma} are analogous for the other parameters. }
  \item{verbose}{\code{TRUE} or \code{FALSE} indicating whether a lengthy output should be printed out.}
}
\value{
  \item{computing_time}{Number of seconds required for the sampling.}
  \item{mu_sample}{Posterior samples for \code{mu} parameters.}
  \item{theta_sample}{Posterior samples for \code{theta} parameters.}
  \item{gamma_sample}{Posterior samples for \code{gamma} parameters.}
  \item{tau_eta_sample}{Posterior samples for \code{tau_eta} parameters.}
  \item{tau_theta_sample}{Posterior samples for \code{tau_theta} parameters.}
  \item{tau_gamma_sample}{Posterior samples for \code{tau_gamma} parameters.}
  \item{acceptance_mu}{Acceptance probability for each of the \code{mu} parameters.}
  \item{acceptance_theta}{Acceptance probability for each of the \code{theta} parameters.}
  \item{acceptance_gamma}{Acceptance probability for each of the \code{gamma} parameters.}
  \item{log_likelihood_values}{Stored values of the log-likelihood function at each iteration.}
  \item{log_posterior_values}{Stored values of the log-posterior function at each iteration.}
}
\examples{
data("simulated")

pars <- list(mu = rep(0,t_frames), 
             theta = rep(0,n_nodes), 
             gamma = rep(0,n_nodes),
             tau_eta = 1, 
             tau_theta = 1, 
             tau_gamma = 1)

hypers <- list(sigma_mu = 100, 
               a_eta = 0.001,
               b_eta = 0.001,
               a_theta = 0.001,
               b_theta = 0.001,
               a_gamma = 0.001,
               b_gamma = 0.001)

mcmc <- list(burnin = 0, thin = 1, n_iterations = 1000)

proposals <- list(mu = rep(1,t_frames), theta = rep(1,n_nodes), gamma = rep(1,n_nodes))

verbose <- TRUE

set.seed(12345)
res <- dirnet_GS(edgelist, model_type = c(1, 1, 1), pars, hypers, mcmc, proposals, verbose)
ts.plot(res$log_posterior_values)
}
