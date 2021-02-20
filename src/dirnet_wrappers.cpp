#include "dirnet.hpp"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List dirnet_posterior_cpp(arma::mat edgelist_, 
                                bool include_mu_, bool include_theta_, bool include_gamma_, 
                                arma::vec mu_, arma::vec theta_, arma::vec gamma_, double tau_eta_, double tau_theta_, double tau_gamma_, 
                                double sigma_mu_, double a_eta_, double b_eta_, double a_theta_, double b_theta_, double a_gamma_, double b_gamma_, 
                                unsigned int burnin_, unsigned int thin_, unsigned int net_n_iter_, 
                                arma::vec proposal_mu_, arma::vec proposal_theta_, arma::vec proposal_gamma_, 
                                bool verbose_) {
  double computing_time;
  arma::wall_clock timer;
  timer.tic();
  dirnet_static network(edgelist_, 
                        include_mu_, include_theta_, include_gamma_, 
                        mu_, theta_, gamma_, tau_eta_, tau_theta_, tau_gamma_, 
                        sigma_mu_, a_eta_, b_eta_, a_theta_, b_theta_, a_gamma_, b_gamma_, 
                        burnin_, thin_, net_n_iter_, 
                        proposal_mu_, proposal_theta_, proposal_gamma_, 
                        verbose_);
  network.EvaluateAllValues();
  computing_time = timer.toc();
  return (Rcpp::List::create(Rcpp::Named("computing_time") = computing_time,
                             Rcpp::Named("log_prior_value") = network.prior_value,
                             Rcpp::Named("log_likelihood_value") = network.likelihood_value,
                             Rcpp::Named("log_posterior_value") = network.posterior_value));
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List dirnet_GS_cpp(arma::mat edgelist_, 
                                bool include_mu_, bool include_theta_, bool include_gamma_, 
                                arma::vec mu_, arma::vec theta_, arma::vec gamma_, double tau_eta_, double tau_theta_, double tau_gamma_, 
                                double sigma_mu_, double a_eta_, double b_eta_, double a_theta_, double b_theta_, double a_gamma_, double b_gamma_, 
                                unsigned int burnin_, unsigned int thin_, unsigned int net_n_iter_, 
                                arma::vec proposal_mu_, arma::vec proposal_theta_, arma::vec proposal_gamma_, 
                                bool verbose_) {
  double computing_time;
  arma::wall_clock timer;
  timer.tic();
  dirnet_static network(edgelist_, 
                        include_mu_, include_theta_, include_gamma_, 
                        mu_, theta_, gamma_, tau_eta_, tau_theta_, tau_gamma_, 
                        sigma_mu_, a_eta_, b_eta_, a_theta_, b_theta_, a_gamma_, b_gamma_, 
                        burnin_, thin_, net_n_iter_, 
                        proposal_mu_, proposal_theta_, proposal_gamma_, 
                        verbose_);
  network.GibbsSampler();
  computing_time = timer.toc();
  return (Rcpp::List::create(Rcpp::Named("computing_time") = computing_time,
                             Rcpp::Named("mu_sample") = network.mu_sample,
                             Rcpp::Named("theta_sample") = network.theta_sample,
                             Rcpp::Named("gamma_sample") = network.gamma_sample,
                             Rcpp::Named("tau_eta_sample") = network.tau_eta_sample,
                             Rcpp::Named("tau_theta_sample") = network.tau_theta_sample,
                             Rcpp::Named("tau_gamma_sample") = network.tau_gamma_sample,
                             Rcpp::Named("acceptance_mu") = network.accepted_counts_mu / network.iteration_counts_mu,
                             Rcpp::Named("acceptance_theta") = network.accepted_counts_theta / network.iteration_counts_theta,
                             Rcpp::Named("acceptance_gamma") = network.accepted_counts_gamma / network.iteration_counts_gamma,
                             Rcpp::Named("log_likelihood_values") = network.likelihood_values,
                             Rcpp::Named("log_posterior_values") = network.posterior_values));
}
