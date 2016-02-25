#' function to get all interactions on all trans chromosomes
#' @param obj 4C-ker object
#' @param k number of observed fragments in each window to build adaptive windows
ineqfunTrans = function(x, env = globalenv()){
   vec = c((-x[15]+x[17]),(-x[13]+x[15]), x[17])
   return(vec)
}
transAnalysis <- function(obj, k){
  #lower and upper quantiles to separate the three HMM states
  lower <- 0.6
  upper <- 0.8
  region = "trans"
  ############
  num_samples <- length(obj@samples)
  window_counts <- buildAdaptiveWindowsTrans(obj@data_trans, obj@samples, obj@chrs_trans,k)
  num_windows <- nrow(window_counts)
  counts_results <- getWindowCounts(obj@data_trans, window_counts, num_windows, obj@samples,obj@output_dir, region)
  if(length(obj@samples) > 1)
    synth_counts_results <- generateSyntheticSamples(window_counts, num_windows, obj@data_trans, region)
  else
    synth_counts_results <- counts_results
  pars_valid <- FALSE
  ##while loop to change the quantile separation if the parameter estimation fails
  while(!pars_valid & lower > 0.1){
    starting_values <- startingValuesTrans(synth_counts_results$synth_hmm_input, lower, upper)
    par_est_results = parameterEstimationTrans(synth_counts_results$synth_hmm_input,obj@replicates,starting_values$trstart,starting_values$respstart, starting_values$instart, ineqfunTrans)
    if(all(par_est_results$pars[4:12] == starting_values$trstart)){
      lower = lower-0.05
      upper = upper-0.05
    }
    else{
      pars_valid = TRUE
    }
  }
  j <- 1
  l <- 1
  for(i in 1:length(obj@conditions)){
    reps = obj@replicates[i]
    samples = obj@samples[l:(l+reps-1)]
    hmm_input <- data.frame(counts_results$hmm_input[j:(j+(num_windows*reps)-1),])
    colnames(hmm_input) = "counts"
    viterbi_results <- viterbi3State(hmm_input, par_est_results$mod_fit,samples,
                                     window_counts, num_windows,
                                     paste(obj@bait_name, obj@conditions[i], sep = "_"),
                                     obj@bait_chr,reps, obj@output_dir, region)
    j <- j+(num_windows*reps)
    l <- l+reps
  }
  print(paste("BED file of highest interacting domains for trans chromosomes are saved in ",
              obj@output_dir, sep = ""))
}

