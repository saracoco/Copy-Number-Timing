plotting <- function(res, input_data, K, simulation_params) {
  draws <- res$draws(format = "matrix")
  names <- paste("tau[", 1:K, "]", sep = "")
  
  areas_tau <- mcmc_areas(
    draws,
    pars = names,
    prob = 0.8,  # 80% intervals
    prob_outer = 0.99,  # 99%
    point_est = "mean"
  ) + labs(
    title = "Approximate Posterior distributions",
    subtitle = "with mean and 80% intervals"
  ) + xlim(0, 1)
  
  color_scheme_set("blue")
  intervals <- mcmc_intervals(draws, regex_pars = c("w"))
  
  # Posterior predictive check
  stanfit <- rstan::read_stan_csv(res$output_files())
  y_rep <- as.matrix(stanfit, pars = "NV_pred")
  y <- input_data$NV
  
  ppc <- ppc_dens_overlay(y = y, yrep = y_rep)
  ecdf_compare <- ppc_ecdf_overlay(y, y_rep)
  intervals_compare <- ppc_intervals(y, y_rep)
  mean_compare <- ppc_stat(y, y_rep, stat = 'mean')
  max_compare <- ppc_stat(y, y_rep, stat = 'max')
  min_compare <- ppc_stat(y, y_rep, stat = 'min')
  median_compare <- ppc_stat(y, y_rep, stat = 'median')
  
  # Simulation parameters plot
  param_plot <- ggplot(simulation_params, aes(x = factor(vector_karyo), y = vector_tau)) +
    geom_bar(stat = "identity", fill = "blue") +
    labs(title = "Simulation Parameters", x = "Karyotype", y = "Tau")
  
  
  
  
  final_plot <- (areas_tau | intervals) / ppc / (ecdf_compare | intervals_compare) /
    (mean_compare | max_compare | min_compare | median_compare) /
    param_plot + 
    plot_layout(widths = c(6, 8, 8, 8), heights = c(15, 12, 8, 8, 10)) + 
    plot_annotation(
      title = 'Title',
      subtitle = "Subtitle",
      caption = "caption"
    ) & theme(
      text = element_text(size = 8), 
      plot.title = element_text(size = 10), 
      plot.subtitle = element_text(size = 8), 
      axis.text = element_text(size = 8), 
      plot.caption = element_text(size = 5)
    )
  
  return(final_plot)
}
