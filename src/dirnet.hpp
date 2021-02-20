#ifndef DIRNET_STATIC_H
#define DIRNET_STATIC_H

#include <RcppArmadillo.h>

class dirnet_static
{
public:
  dirnet_static(arma::mat, 
                bool, bool, bool, 
                arma::vec, arma::vec, arma::vec, double, double, double, 
                double, double, double, double, double, double, double, 
                unsigned int, unsigned int, unsigned int, 
                arma::vec, arma::vec, arma::vec, 
                bool);
  bool verbose;
  
  // global
  unsigned int T;
  unsigned int N;
  unsigned int L;
  
  // data
  arma::mat edgelist;
  arma::cube log_adj;
  arma::mat n_in_edges;
  arma::mat n_out_edges;
  
  // model identification
  bool include_mu;
  bool include_theta;
  bool include_gamma;
  
  // model parameters
  arma::vec mu;
  arma::vec theta;
  arma::vec gamma;
  double tau_eta;
  double tau_theta;
  double tau_gamma;
  
  // statistics
  double sum_of_gamma_exp;
  double sum_of_squared_mu_increments;
  double sum_of_squared_theta;
  double sum_of_squared_gamma;
  
  // hyperparameters
  double sigma_mu;
  double a_eta;
  double b_eta;
  double a_theta;
  double b_theta;
  double a_gamma;
  double b_gamma;
  
  // inference
  double prior_value;
  double likelihood_value;
  double posterior_value;
  
  // global functions
  void EvaluateAllValues();
  void EvaluateDegrees();
  void EvaluateAdjacency();
  void EvaluateStatistics();
  
  // inference functions
  void EvaluatePrior();
  void EvaluateLikelihood();
  void EvaluatePosterior();
  void CheckValues();
  
  // mcmc functions
  double UpdateMu(unsigned int, double);
  double UpdateTheta(unsigned int, double);
  double UpdateGamma(unsigned int, double);
  void UpdateTauEta();
  void UpdateTauTheta();
  void UpdateTauGamma();
  
  // MCMC
  unsigned int burnin;
  unsigned int thin;
  unsigned int total_n_iter;
  unsigned int net_n_iter;
  arma::mat mu_sample;
  arma::mat theta_sample;
  arma::mat gamma_sample;
  arma::vec tau_eta_sample;
  arma::vec tau_theta_sample;
  arma::vec tau_gamma_sample;
  arma::vec likelihood_values;
  arma::vec posterior_values;
  arma::vec proposal_mu;
  arma::vec proposal_theta;
  arma::vec proposal_gamma;
  arma::vec accepted_counts_mu;
  arma::vec iteration_counts_mu;
  arma::vec accepted_counts_theta;
  arma::vec iteration_counts_theta;
  arma::vec accepted_counts_gamma;
  arma::vec iteration_counts_gamma;
  
  void GibbsSampler();
  
protected:
private:
};

#endif // DIRNET_STATIC_H
