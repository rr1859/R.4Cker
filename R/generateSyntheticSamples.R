generateSyntheticSamples <- function(window_counts, num_windows, data, region){
  set.seed(2)
  n_art_samples <- length(data)
  reps <- length(data)
  chrs <- unique(window_counts[,1])
  non_overlap_windows_all <- NULL
  print("Generating synthetic samples....")
  for(chr in chrs){
    window_counts_chr <- window_counts[window_counts[,1] == chr,]
    num_windows_chr <- nrow(window_counts_chr)
    non_overlap_windows <- data.frame(t(sapply(1:(num_windows_chr-1), function(i) c(window_counts_chr[i,2], window_counts_chr[(i+1),2]))))
    non_overlap_windows <- rbind(non_overlap_windows, c(window_counts_chr[(num_windows_chr-1),3], window_counts_chr[num_windows_chr,3]))
    non_overlap_windows <- cbind(rep(window_counts_chr[1,1], nrow(non_overlap_windows)),non_overlap_windows )
    non_overlap_windows_all <- rbind(non_overlap_windows_all, non_overlap_windows)
  }
  sample_n <- sapply(1:n_art_samples, function(i) sample(1:reps, num_windows, replace = TRUE))
  all_art_samples <- vector("list", n_art_samples)
  window_counts_art <- window_counts[,1:3]
  for(i in 1:n_art_samples){
    art_sample <- NULL
    for(j in 1:reps){
      data_j <- data[[j]]
      rows_sample <- which(sample_n[,i] == j)
      art_sample <- rbind(art_sample, do.call(rbind, lapply(rows_sample, function(k) data_j[data_j[,1] == as.character(non_overlap_windows_all[k,1]) & data_j[,2] > non_overlap_windows_all[k,2] & data_j[,2] < non_overlap_windows_all[k,3],])))
    }
    art_sample <- art_sample[order(art_sample[,2]),]
    window_counts_art <- cbind(window_counts_art, unlist(lapply(1:num_windows, function(i)
      removePCR(data <- art_sample[art_sample[,1] == as.character(window_counts[i,1]) & art_sample[,2] > window_counts[i,2] & art_sample[,2] < window_counts[i,3],4]))))
    all_art_samples[[i]] <- art_sample
  }
  colnames(window_counts_art) <- c("chr", "start","end", unlist(lapply(1:n_art_samples, function(i) paste("sample",i, sep = "_"))))
  norm_counts_art <- normalizeCounts(window_counts_art[,-c(1:3)])
  norm_counts_art_log <- log(norm_counts_art+1,10)
  if(region == "cis" | region == "nearbait"){
    dist_log <- log(window_counts[,4]+1,10)
    hmm_input_art <- data.frame(counts = c(norm_counts_art_log), distance = rep(dist_log, n_art_samples), row.names = NULL)
    return(list(hmm_input = hmm_input_art, norm_counts_log=norm_counts_art_log, dist_log = dist_log))
  }
  if(region == "trans"){
    hmm_input_art <- data.frame(counts = c(norm_counts_art_log), row.names = NULL)
    return(list(hmm_input = hmm_input_art))
  }
}

