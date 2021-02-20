#include "dirnet.hpp"

dirnet_static::dirnet_static(arma::mat edgelist_, 
                             bool include_mu_, bool include_theta_, bool include_gamma_, 
                             arma::vec mu_, arma::vec theta_, arma::vec gamma_, double tau_eta_, double tau_theta_, double tau_gamma_, 
                             double sigma_mu_, double a_eta_, double b_eta_, double a_theta_, double b_theta_, double a_gamma_, double b_gamma_, 
                             unsigned int burnin_, unsigned int thin_, unsigned int net_n_iter_, 
                             arma::vec proposal_mu_, arma::vec proposal_theta_, arma::vec proposal_gamma_, 
                             bool verbose_)
{
  T = mu_.n_elem;
  N = theta_.n_elem;
  edgelist = edgelist_;
  L = edgelist.n_rows;
  include_mu = include_mu_;
  include_theta = include_theta_;
  include_gamma = include_gamma_;
  mu = mu_;
  theta = theta_;
  gamma = gamma_;
  gamma.at(0) = - accu(gamma) + gamma.at(0);
  tau_eta = tau_eta_;
  tau_theta = tau_theta_;
  tau_gamma = tau_gamma_;
  sigma_mu = sigma_mu_;
  a_eta = a_eta_;
  b_eta = b_eta_;
  a_theta = a_theta_;
  b_theta = b_theta_;
  a_gamma = a_gamma_;
  b_gamma = b_gamma_;
  verbose = verbose_;
  EvaluateAllValues();
  
  net_n_iter = net_n_iter_;
  burnin = burnin_;
  thin = thin_;
  total_n_iter = burnin + net_n_iter*thin;
  mu_sample.zeros(T,net_n_iter);
  theta_sample.zeros(N,net_n_iter);
  gamma_sample.zeros(N,net_n_iter);
  tau_eta_sample.zeros(net_n_iter);
  tau_theta_sample.zeros(net_n_iter);
  tau_gamma_sample.zeros(net_n_iter);
  likelihood_values.zeros(net_n_iter);
  posterior_values.zeros(net_n_iter);
  proposal_mu = proposal_mu_;
  proposal_theta = proposal_theta_;
  proposal_gamma = proposal_gamma_;
  accepted_counts_mu.zeros(T);
  accepted_counts_theta.zeros(N);
  accepted_counts_gamma.zeros(N);
  iteration_counts_mu.ones(T);// these three are overwritten just below
  iteration_counts_theta.ones(N);
  iteration_counts_gamma.ones(N);
  if (total_n_iter > 0)
  {
    iteration_counts_mu.zeros(T);
    iteration_counts_theta.zeros(N);
    iteration_counts_gamma.zeros(N);
  }
}

void dirnet_static::EvaluateAllValues()
{
  EvaluateDegrees();
  EvaluateAdjacency();
  EvaluateStatistics();
  EvaluatePrior();
  EvaluateLikelihood();
  EvaluatePosterior();
}

void dirnet_static::EvaluateDegrees()
{
  unsigned int l, t, i, j;
  n_in_edges.zeros(N,T);
  n_out_edges.zeros(N,T);
  for (l=0; l<L; ++l)
  {
    t = edgelist.at(l,0);
    i = edgelist.at(l,1);
    j = edgelist.at(l,2);
    n_out_edges.at(i,t) += 1;
    n_in_edges.at(j,t) += 1;
  }
}

void dirnet_static::EvaluateAdjacency()
{
  unsigned int l, t, i, j;
  double value;
  log_adj.set_size(N,N,T);
  log_adj.fill(-100);// Extremely low value to avoid zeros in dirichlet likelihood
  for (l=0; l<L; ++l)
  {
    t = edgelist.at(l,0);
    i = edgelist.at(l,1);
    j = edgelist.at(l,2);
    value = edgelist.at(l,3);
    log_adj.at(i,j,t) = value;
  }
}

void dirnet_static::EvaluateStatistics()
{
  unsigned int t, i, j;
  gamma.at(0) = - accu(gamma) + gamma.at(0);// this is repeated here for debug purposes
  sum_of_gamma_exp = 0;
  sum_of_squared_mu_increments = 0;
  sum_of_squared_theta = 0;
  sum_of_squared_gamma = 0;
  for (j=0; j<N; ++j) sum_of_gamma_exp += exp(gamma.at(j));
  for (t=1; t<T; ++t) sum_of_squared_mu_increments += (mu.at(t)-mu.at(t-1)) * (mu.at(t)-mu.at(t-1));
  for (i=0; i<N; ++i) sum_of_squared_theta += theta.at(i) * theta.at(i);
  for (j=1; j<N; ++j) sum_of_squared_gamma += gamma.at(j) * gamma.at(j);
}

void dirnet_static::EvaluatePrior()
{
  unsigned int t, i, j;
  prior_value = 0;
  if (include_mu) prior_value += R::dgamma(tau_eta,a_eta,1/b_eta,true);
  if (include_theta) prior_value += R::dgamma(tau_theta,a_theta,1/b_theta,true);
  if (include_gamma) prior_value += R::dgamma(tau_gamma,a_gamma,1/b_gamma,true);
  prior_value += R::dnorm(mu.at(0),0,sigma_mu,true);
  if (include_mu) for (t=1; t<T; ++t) prior_value += R::dnorm(mu.at(t)-mu.at(t-1),0,1/sqrt(tau_eta),true);
  if (include_theta) for (i=0; i<N; ++i) prior_value += R::dnorm(theta.at(i),0,1/sqrt(tau_theta),true);
  if (include_gamma) for (j=1; j<N; ++j) prior_value += R::dnorm(gamma.at(j),0,1/sqrt(tau_gamma),true);
}

void dirnet_static::EvaluateLikelihood()
{
  unsigned int t, i, j;
  likelihood_value = 0;
  for (t=0; t<T; ++t) for (i=0; i<N; ++i) if (n_out_edges.at(i,t) > 0)
  {
    likelihood_value += lgamma( exp(mu.at(t)+theta.at(i)) * (sum_of_gamma_exp-exp(gamma.at(i))) );
    for (j=0; j<N; ++j) if (i != j) 
    {
      likelihood_value += - lgamma(exp(mu.at(t)+theta.at(i)+gamma.at(j)));
      likelihood_value += (exp(mu.at(t)+theta.at(i)+gamma.at(j)) - 1) * log_adj.at(i,j,t);
    }
  }
}

void dirnet_static::EvaluatePosterior()
{
  posterior_value = prior_value + likelihood_value;
}

void dirnet_static::CheckValues()// this is a debug function for the mcmc updates
{
  double prior_check = prior_value;
  double likelihood_check = likelihood_value;
  EvaluateAllValues();
  std::cout << "DEBUG: \tError on prior = \t\t\t\t\t\t\t\t\t" << std::abs(prior_value - prior_check)  << std::endl;
  std::cout << "DEBUG: \tError on likelihood = \t\t\t\t\t\t\t\t\t" << std::abs(likelihood_value - likelihood_check) << std::endl;
}
