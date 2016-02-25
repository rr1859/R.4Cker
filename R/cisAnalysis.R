#' function to get all interactions on the bait chromosome
#' @param obj 4C-ker object
#' @param k number of observed fragments in each window to build adaptive windows
cisAnalysis <- function(obj, k){
  #lower and upper quantiles to separate the three HMM states
  lower <- 0.6
  upper <- 0.9
  ############
  num_samples <- length(obj@samples)
  #build adaptive windows
  window_counts <- buildAdaptiveWindowsCis(obj@data_cis, obj@bait_coord,
                                           obj@bait_chr, obj@bait_name, obj@output_dir,k, "cis")
  num_windows <- nrow(window_counts)
  #get count for each sample for each window
  counts_results <- getWindowCounts(obj@data_cis, window_counts, num_windows, obj@samples,obj@output_dir, "cis")
  #build synthetic samples by shuffling reads
  if(length(obj@samples) > 1)
    synth_counts_results <- generateSyntheticSamples(window_counts, num_windows, obj@data_cis, "cis")
  else
    synth_counts_results <- counts_results
  ####parameter estimation####
  pars_valid <- FALSE
  ##while loop to change the quantile separation if the parameter estimation fails
  while(!pars_valid & lower > 0.1){
 
    starting_values <- startingValuesCis(window_counts, num_windows, obj@bait_coord,
                                          num_samples,synth_counts_results$norm_counts_log,
                                          synth_counts_results$dist_log,lower, upper)
    int_glm <- starting_values$int_glm
    lint_glm <- starting_values$int_glm
    nint_glm <- starting_values$int_glm

    par_est_results <- parameterEstimationCis(synth_counts_results$hmm_input,num_samples,
                                              starting_values$trstart,starting_values$respstart,
                                              starting_values$instart, ineqfunCis)
    is_valid <- validateParametersCis(starting_values$int_glm, par_est_results$pars,
                                      synth_counts_results$hmm_input)
    par_est_results$pars <- is_valid$pars
    int_glm$coefficients <- par_est_results$pars[19:20]
    lint_glm$coefficients <- par_est_results$pars[16:17]
    nint_glm$coefficients <- par_est_results$pars[13:14]
    if(all(par_est_results$pars[13:21] == starting_values$respstart) | !is_valid$valid){
      lower <- lower-0.05
      upper <- upper-0.05
    } else
      pars_valid <- TRUE
  }
  ####end par est####

  ####viterbi####
  j <- 1
  l <- 1
  for(i in 1:length(obj@conditions)){
    reps = obj@replicates[i]
    samples = obj@samples[l:(l+reps-1)]
    hmm_input <- counts_results$hmm_input[j:(j+(num_windows*reps)-1),]
    viterbi_results <- viterbi3State(hmm_input, par_est_results$mod_fit,samples,
                                        window_counts, num_windows,
                                        paste(obj@bait_name, obj@conditions[i], sep = "_"),
                                        obj@bait_chr,reps, obj@output_dir, "cis")
    j <- j+(num_windows*reps)
    l <- l+reps
  }
  print(paste("BED file of highest interacting domains for the cis chromosome are saved in ",
              obj@output_dir, sep = ""))
  ####end viterbi####
  norm_counts_avg=NULL
  j=1
  for(i in 1:length(obj@replicates)){
    reps=obj@replicates[i]
    norm_counts_avg = rbind(norm_counts_avg, cbind(rowMeans(counts_results$norm_counts[,j:(j+reps-1)]), rep(obj@conditions[i],nrow(window_counts))))
    j=j+i
  }
  norm_counts_avg = data.frame(cbind(rowMeans(window_counts[,2:3]),norm_counts_avg))
  norm_counts_avg[,1] =as.numeric(as.character(norm_counts_avg[,1]))
  norm_counts_avg[,2] =as.numeric(as.character(norm_counts_avg[,2]))
  colnames(norm_counts_avg) = c("Coord", "Count", "Condition")
  return(list(norm_counts_avg=norm_counts_avg, window_counts=counts_results$window_counts))
}




