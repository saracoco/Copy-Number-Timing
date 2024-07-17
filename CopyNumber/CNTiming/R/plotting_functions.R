
#' plotting Function
#'
#' This function obtains the peaks to be used for CN timing in the input data.
#' @param res result
#' @param input_data data
#' @param K number of mixture components
#' @keywords input
#' @export
#' @examples
#' plotting()


plotting <- function(res, input_data, K){

  draws <- res$draws(format = "matrix")

  #color_scheme_set("red")
  #manage to get the K from the model fit directly rather than as input
  names <- paste("tau[", 1:K, "]", sep = "")

  areas_tau <- mcmc_areas(
    draws,
    pars = names,
    prob = 0.8, # 80% intervals
    prob_outer = 0.99, # 99%
    point_est = "mean"
  )+
    labs(
      title = "Approximate Posterior distributions",
      subtitle = "with mean and 80% intervals"
    )+
    xlim(0, 1)

  color_scheme_set("blue")
  intervals <- mcmc_intervals(draws, regex_pars = c("w"))


  #posterior predictive check
  stanfit <- rstan::read_stan_csv(res$output_files())
  y_rep <- as.matrix(stanfit, pars = "NV_pred")
  y = input_data$NV
  #distribution of replicated data vs real data
  ppc <- ppc_dens_overlay(
    y = input_data$NV,
    yrep = y_rep)
  #empirical cumulative distribution
  ecdf_compare <-ppc_ecdf_overlay(y,y_rep)

  #predictive intervals vs observed values
  intervals_compare <- ppc_intervals(y,y_rep)

  #Compare statistics
  #compute the estimated Bayesian p-value
  mean_compare <- ppc_stat(y,y_rep, stat = 'mean')
  max_compare <- ppc_stat(y,y_rep, stat = 'max')
  min_compare <- ppc_stat(y,y_rep, stat = 'min')
  median_compare <- ppc_stat(y,y_rep, stat = 'median')
  #}


    final_plot <- (areas_tau | intervals) / ppc /(ecdf_compare | intervals_compare) / (mean_compare|max_compare|min_compare|median_compare) +
      plot_layout(widths = c(6, 8, 8, 8), heights = c(15, 12, 8, 8)) +
      plot_annotation(
        title = 'Title  ',
        subtitle = "Subtitle",
        caption = "caption"
      ) & theme(text = element_text(size = 8), plot.title = element_text(size = 10), plot.subtitle = element_text(size = 8), axis.text = element_text(size = 8), plot.caption = element_text(size = 5))

    # ggsave("./plots/plots_inference.pdf", plot = final_plot, width = 12, height = 24)
    # ggsave("./plots/plots_inference.png", width = 12, height = 24, device = "png")
    #



  return(final_plot)
}
