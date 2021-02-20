dirnet_posterior <- function(edgelist, model_type, pars, hypers, verbose = F)
{
  edgelist[,1] = edgelist[,1] - 1
  edgelist[,2] = edgelist[,2] - 1
  edgelist[,3] = edgelist[,3] - 1
  dirnet_posterior_cpp(edgelist, 
                              model_type[1], model_type[2], model_type[3], 
                              pars$mu, pars$theta, pars$gamma, pars$tau_eta, pars$tau_theta, pars$tau_gamma, 
                              hypers$sigma_mu, hypers$a_eta, hypers$b_eta, hypers$a_theta, hypers$b_theta, hypers$a_gamma, hypers$b_gamma, 
                              0, 5, 0,
                              1:5, 1:5, 1:5, 
                              verbose)
}

dirnet_GS <- function(edgelist, model_type, pars, hypers, mcmc, proposals, verbose = F)
{
  edgelist[,1] = edgelist[,1] - 1
  edgelist[,2] = edgelist[,2] - 1
  edgelist[,3] = edgelist[,3] - 1
  dirnet_GS_cpp(edgelist, 
                       model_type[1], model_type[2], model_type[3], 
                       pars$mu, pars$theta, pars$gamma, pars$tau_eta, pars$tau_theta, pars$tau_gamma, 
                       hypers$sigma_mu, hypers$a_eta, hypers$b_eta, hypers$a_theta, hypers$b_theta, hypers$a_gamma, hypers$b_gamma, 
                       mcmc$burnin, mcmc$thin, mcmc$n_iterations,
                       proposals$mu, proposals$theta, proposals$gamma, 
                       verbose)
}

