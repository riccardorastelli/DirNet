#include "dirnet.hpp"

double dirnet_static::UpdateMu(unsigned int t, double sd)
{
  // if (verbose) std::cout << "\nUpdate of mu (" << t << ") is starting" << std::endl;//                    DEBUG ONLY
  double res = 0;
  unsigned int i, j;
  double prior_delta, likelihood_delta;
  double increment, mu_new;
  increment = R::rnorm(0,sd);
  mu_new = mu.at(t) + increment;
  
  prior_delta = 0;
  if (t == 0) prior_delta += - (mu_new*mu_new - mu.at(t)*mu.at(t)) / (2*sigma_mu*sigma_mu);
  if (include_mu) if (t > 0) prior_delta += - 0.5 * tau_eta * (   (mu_new-mu.at(t-1))*(mu_new-mu.at(t-1)) - (mu.at(t)-mu.at(t-1))*(mu.at(t)-mu.at(t-1))   );
  if (include_mu) if (t < T-1) prior_delta += - 0.5 * tau_eta * (   (mu_new-mu.at(t+1))*(mu_new-mu.at(t+1)) - (mu.at(t)-mu.at(t+1))*(mu.at(t)-mu.at(t+1))   );
  
  if (include_mu)
  {
    likelihood_delta = 0;
    for (i=0; i<N; ++i) if (n_out_edges.at(i,t) > 0)
    {
      for (j=0; j<N; ++j) if (i != j) 
      {
        likelihood_delta += log_adj.at(i,j,t) * (  exp(mu_new+theta.at(i)+gamma.at(j))-exp(mu.at(t)+theta.at(i)+gamma.at(j))  );
        likelihood_delta += - lgamma(exp(mu_new+theta.at(i)+gamma.at(j))) + lgamma(exp(mu.at(t)+theta.at(i)+gamma.at(j)));
      }
      likelihood_delta += lgamma(exp(mu_new+theta.at(i))*(sum_of_gamma_exp-exp(gamma.at(i)))) - lgamma(exp(mu.at(t)+theta.at(i))*(sum_of_gamma_exp-exp(gamma.at(i))));
    }
  }
  
  if (!include_mu)
  {
    likelihood_delta = 0;
    for (unsigned int tt = 0; tt < T; ++tt) for (i=0; i<N; ++i) if (n_out_edges.at(i,tt) > 0)
    {
      for (j=0; j<N; ++j) if (i != j) 
      {
        likelihood_delta += log_adj.at(i,j,tt) * (  exp(mu_new+theta.at(i)+gamma.at(j))-exp(mu.at(tt)+theta.at(i)+gamma.at(j))  );
        likelihood_delta += - lgamma(exp(mu_new+theta.at(i)+gamma.at(j))) + lgamma(exp(mu.at(tt)+theta.at(i)+gamma.at(j)));
      }
      likelihood_delta += lgamma(exp(mu_new+theta.at(i))*(sum_of_gamma_exp-exp(gamma.at(i)))) - lgamma(exp(mu.at(tt)+theta.at(i))*(sum_of_gamma_exp-exp(gamma.at(i))));
    }
  }
  
  if (log(R::runif(0,1)) < prior_delta + likelihood_delta)
  {
    res += 1;
    if (include_mu) if (t > 0)
    {
      sum_of_squared_mu_increments -= (mu.at(t)-mu.at(t-1)) * (mu.at(t)-mu.at(t-1));
      sum_of_squared_mu_increments += (mu_new-mu.at(t-1)) * (mu_new-mu.at(t-1));
    }
    if (include_mu) if (t < T-1)
    {
      sum_of_squared_mu_increments -= (mu.at(t)-mu.at(t+1)) * (mu.at(t)-mu.at(t+1));
      sum_of_squared_mu_increments += (mu_new-mu.at(t+1)) * (mu_new-mu.at(t+1));
    }
    if (!include_mu) for (unsigned tt = 0; tt < T; ++tt) mu.at(tt) = mu_new;
    if (include_mu) mu.at(t) = mu_new;
    prior_value += prior_delta;
    likelihood_value += likelihood_delta;
    EvaluatePosterior();
    // CheckValues();//                    DEBUG ONLY
  }
  // if (verbose) std::cout << "Update of mu (" << t << ") has ended with res = " << res << std::endl << std::endl;//                    DEBUG ONLY
  return (res);
}

double dirnet_static::UpdateTheta(unsigned int i, double sd)
{
  // if (verbose) std::cout << "\nUpdate of theta (" << i << ") is starting" << std::endl;//                    DEBUG ONLY
  double res = 0;
  unsigned int t, j;
  double prior_delta, likelihood_delta;
  double increment, theta_new;
  increment = R::rnorm(0,sd);
  theta_new = theta.at(i) + increment;
  
  prior_delta = - 0.5 * tau_theta * (theta_new*theta_new - theta.at(i)*theta.at(i));
  
  likelihood_delta = 0;
  for (t=0; t<T; ++t) if (n_out_edges.at(i,t) > 0)
  {
    likelihood_delta += lgamma(exp(mu.at(t)+theta_new)*(sum_of_gamma_exp-exp(gamma.at(i)))) - lgamma(exp(mu.at(t)+theta.at(i))*(sum_of_gamma_exp-exp(gamma.at(i))));
    for (j=0; j<N; ++j) if (i != j) 
    {
      likelihood_delta += log_adj.at(i,j,t) * exp(mu.at(t)+gamma.at(j)) * (  exp(theta_new)-exp(theta.at(i))  );
      likelihood_delta += - lgamma(exp(mu.at(t)+theta_new+gamma.at(j))) + lgamma(exp(mu.at(t)+theta.at(i)+gamma.at(j)));
    }
  }
  
  if (log(R::runif(0,1)) < prior_delta + likelihood_delta)
  {
    res += 1;
    sum_of_squared_theta -= theta.at(i) * theta.at(i);
    sum_of_squared_theta += theta_new * theta_new;
    theta.at(i) = theta_new;
    prior_value += prior_delta;
    likelihood_value += likelihood_delta;
    EvaluatePosterior();
    // CheckValues();//                    DEBUG ONLY
  }
  // if (verbose) std::cout << "Update of theta (" << i << ") has ended with res = " << res << std::endl << std::endl;//                    DEBUG ONLY
  return (res);
}

double dirnet_static::UpdateGamma(unsigned int j, double sd)
{
  // if (verbose) std::cout << "\nUpdate of gamma (" << j << ") is starting" << std::endl;//                    DEBUG ONLY
  double res = 0;
  unsigned int t, i;
  double prior_delta, likelihood_delta;
  arma::vec gamma_new;
  double increment, sum_of_gamma_exp_new;
  increment = R::rnorm(0,sd);
  gamma_new = gamma;
  gamma_new.at(j) += increment;
  gamma_new.at(0) -= increment;
  sum_of_gamma_exp_new = sum_of_gamma_exp - exp(gamma.at(j)) + exp(gamma_new.at(j)) - exp(gamma.at(0)) + exp(gamma_new.at(0));
  
  prior_delta = - 0.5 * tau_gamma * (gamma_new.at(j)*gamma_new.at(j) - gamma.at(j)*gamma.at(j));
  
  likelihood_delta = 0;
  for (t=0; t<T; ++t) for (i=0; i<N; ++i) if (n_out_edges.at(i,t) > 0)
  {
    likelihood_delta += lgamma(exp(mu.at(t)+theta.at(i))*(sum_of_gamma_exp_new-exp(gamma_new.at(i)))) - lgamma(exp(mu.at(t)+theta.at(i))*(sum_of_gamma_exp-exp(gamma.at(i))));
    if (i != j)
    {
      likelihood_delta += - lgamma(exp(mu.at(t)+theta.at(i)+gamma_new.at(j))) + lgamma(exp(mu.at(t)+theta.at(i)+gamma.at(j)));
      likelihood_delta += log_adj.at(i,j,t) * exp(mu.at(t)+theta.at(i)) * (  exp(gamma_new.at(j))-exp(gamma.at(j))  );
    }
    if (i != 0)
    {
      likelihood_delta += - lgamma(exp(mu.at(t)+theta.at(i)+gamma_new.at(0))) + lgamma(exp(mu.at(t)+theta.at(i)+gamma.at(0)));
      likelihood_delta += log_adj.at(i,0,t) * exp(mu.at(t)+theta.at(i)) * (  exp(gamma_new.at(0))-exp(gamma.at(0))  );
    }
  }
  
  if (log(R::runif(0,1)) < prior_delta + likelihood_delta)
  {
    res += 1;
    sum_of_squared_gamma += -gamma.at(j)*gamma.at(j) + gamma_new.at(j)*gamma_new.at(j);
    gamma.at(j) = gamma_new.at(j);
    gamma.at(0) = gamma_new.at(0);
    sum_of_gamma_exp = sum_of_gamma_exp_new;
    prior_value += prior_delta;
    likelihood_value += likelihood_delta;
    EvaluatePosterior();
    // CheckValues();//                    DEBUG ONLY
  }
  // if (verbose) std::cout << "Update of gamma (" << j << ") has ended with res = " << res << std::endl << std::endl;//                    DEBUG ONLY
  return (res);
}

void dirnet_static::UpdateTauEta()
{
  // if (verbose) std::cout << "\nUpdate of tau_eta is starting" << std::endl;//                    DEBUG ONLY
  double par1 = a_eta + 0.5*((double)T-1);
  // double par1 = a_eta + (double)T - 1;
  double par2 = b_eta + sum_of_squared_mu_increments/2;
  double tau_eta_new = R::rgamma(par1, 1/par2);
  prior_value += ((double)par1-1) * (log(tau_eta_new)-log(tau_eta)) - par2 * (tau_eta_new-tau_eta);
  EvaluatePosterior();
  tau_eta = tau_eta_new;
  // CheckValues();//                    DEBUG ONLY
  // if (verbose) std::cout << "Update of tau_eta has ended" << std::endl << std::endl;//                    DEBUG ONLY
}

void dirnet_static::UpdateTauTheta()
{
  // if (verbose) std::cout << "\nUpdate of tau_theta is starting" << std::endl;//                    DEBUG ONLY
  double par1 = a_theta + 0.5*((double)N);
  double par2 = b_theta + sum_of_squared_theta/2;
  double tau_theta_new = R::rgamma(par1, 1/par2);
  prior_value += ((double)par1-1) * (log(tau_theta_new)-log(tau_theta)) - par2 * (tau_theta_new-tau_theta);
  EvaluatePosterior();
  tau_theta = tau_theta_new;
  // CheckValues();//                    DEBUG ONLY
  // if (verbose) std::cout << "Update of tau_theta has ended" << std::endl << std::endl;//                    DEBUG ONLY
}

void dirnet_static::UpdateTauGamma()
{
  // if (verbose) std::cout << "\nUpdate of tau_gamma is starting" << std::endl;//                    DEBUG ONLY
  double par1 = a_gamma + 0.5*((double)N-1);
  double par2 = b_gamma + sum_of_squared_gamma/2;
  double tau_gamma_new = R::rgamma(par1, 1/par2);
  prior_value += ((double)par1-1) * (log(tau_gamma_new)-log(tau_gamma)) - par2 * (tau_gamma_new-tau_gamma);
  EvaluatePosterior();
  tau_gamma = tau_gamma_new;
  // CheckValues();//                    DEBUG ONLY
  // if (verbose) std::cout << "Update of tau_gamma has ended" << std::endl << std::endl;//                    DEBUG ONLY
}

void dirnet_static::GibbsSampler()
{
  if (verbose) std::cout << "\nGibbs sampling has started ..." << std::endl;
  unsigned int t, i, j, index_inner, index_outer;
  arma::wall_clock timer;
  timer.tic();
  index_outer = 0;
  index_inner = 0;
  while (index_inner < net_n_iter)
  {
    if (!include_mu) 
    {
      accepted_counts_mu.at(0) += UpdateMu(0,proposal_mu.at(0));
      iteration_counts_mu.at(0) ++;
    }
    if (include_mu) for (t=0; t<T; ++t)
    {
      accepted_counts_mu.at(t) += UpdateMu(t,proposal_mu.at(t));
      iteration_counts_mu.at(t) ++;
    }
    if (include_theta) for (i=0; i<N; ++i)
    {
      accepted_counts_theta.at(i,t) += UpdateTheta(i,proposal_theta.at(i,t));
      iteration_counts_theta.at(i,t) ++;
    }
    if (include_gamma) for (j=0; j<N; ++j) if (j != 0)
    {
      accepted_counts_gamma.at(j) += UpdateGamma(j,proposal_gamma.at(j));
      iteration_counts_gamma.at(j) ++;
    }
    if (include_mu) UpdateTauEta();
    if (include_theta) UpdateTauTheta();
    if (include_gamma) UpdateTauGamma();
    if (index_outer > burnin) if (index_outer % thin == 0)
    {
      mu_sample.col(index_inner) = mu;
      theta_sample.col(index_inner) = theta;
      gamma_sample.col(index_inner) = gamma;
      tau_eta_sample.at(index_inner) = tau_eta;
      tau_theta_sample.at(index_inner) = tau_theta;
      tau_gamma_sample.at(index_inner) = tau_gamma;
      likelihood_values.at(index_inner) = likelihood_value;
      posterior_values.at(index_inner) = posterior_value;
      ++index_inner;
    }
    if (index_outer % 100 == 0) if (verbose) std::cout << "Elapsed Time " << floor(10*timer.toc())/10 << "\t\tEnd of iteration " << index_outer << " out of " << total_n_iter << "\t\tCurrent posterior value: " << posterior_value << std::endl;
    ++index_outer;
  }
  if (verbose) std::cout << "... Gibbs sampling has terminated after " << floor(10*timer.toc())/10 << " seconds\n" << std::endl;
  // CheckValues();//                DEBUG ONLY
}
