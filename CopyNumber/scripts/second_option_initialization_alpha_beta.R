

# SECOND OPTION INITIALIZATION (to be completed): extract the alpha and beta as the #mutations under the 2 theoretical peaks
#perform clustering on the events by setting the data as points with 2 variables : alpha - beta
#data_kmeans <- data.frame()

for (segment_idx in 1:n_segments) {
  alpha = .05
  
  print(segment_idx)
  
  # Segments
  segment_mutations <- all_sim %>% filter(segment_id==segment_idx)
  
  peaks <- get_clonal_peaks(data$karyo[segment_idx], purity=1)
  
  accepted_mutations <- data.frame()
  
  # Check if mutation is inside CI
  probs <- c(alpha/2, 1 - alpha/2)
  
  DP <- segment_mutations$DP
  NV <- segment_mutations$NV
  
  alpha_idx <- lapply(1:length(DP), function(i) {
    p = 1
    quantiles <- qbinom(probs, DP[i], p)
    
    if ((NV[i] >= quantiles[1]) && (NV[i] <= quantiles[2])) {
      return(i)
    }
  }) %>% unlist()
  
  beta_idx <- lapply(1:length(DP), function(i) {
    p = 2
    quantiles <- qbinom(probs, DP[i], p)
    
    if ((NV[i] >= quantiles[1]) && (NV[i] <= quantiles[2])) {
      return(i)
    }
  }) %>% unlist()
  
  # Get only good mutations
  accepted_mutations <- data.frame(DP = DP[accepted_idx], NV = NV[accepted_idx])
  id = c(id,segment_idx)
  alpha = c(alpha,)
  
}