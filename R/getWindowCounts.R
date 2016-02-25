getWindowCounts <- function(data,window_counts,num_windows,samples,output_dir, region){
  for(i in 1:length(samples)){
    data_sample <- data[[i]]
    window_counts <- cbind(window_counts,unlist(lapply(1:num_windows, function(i)
                    removePCR(data = data_sample[data_sample[,1] == as.character(window_counts[i,1]) & data_sample[,2] > window_counts[i,2] & data_sample[,2] < window_counts[i,3],4]))))
  }
  if(region == "trans"){
    colnames(window_counts) <- c("chr", "start", "end", sapply(4:ncol(window_counts), function(i) paste("sample",(i-3), sep = "_")))
    print("Normalizing counts...")
    norm_counts <- normalizeCounts(window_counts[,-c(1:3)])
    norm_counts_log <- log(norm_counts+1,10)
    write.table(cbind(window_counts, norm_counts), paste(output_dir, "_", region, "_norm_counts.txt", sep = ""), quote = FALSE, row.names = FALSE)
    for(i in 1:length(samples)){
      write.table(cbind(window_counts,norm_counts[,i]),
                  paste(output_dir,samples[i], "_",region,"_norm_counts.bedGraph", sep = ""),
                  quote =FALSE, col.names = FALSE, row.names = FALSE)
    }
    hmm_input <- data.frame(counts = c(norm_counts_log),row.names = NULL)
    return(list(hmm_input = hmm_input, norm_counts_log=norm_counts_log))
  }
  else{
    colnames(window_counts) <- c("chr", "start", "end",  "distance", sapply(5:ncol(window_counts), function(i) paste("sample",(i-4), sep = "_")))
    if(length(samples) > 1){
      print("Normalizing counts...")
      norm_counts <- normalizeCounts(window_counts[,-c(1:4)])
    }
    else{
      norm_counts <- data.frame(counts=window_counts[,5])
    }
    dist_log <- log(window_counts$distance+1,10)
    norm_counts_log <- log(norm_counts+1,10)
    for(i in 1:length(samples)){
      write.table(cbind(window_counts[,1:3],norm_counts[,i]),
                  paste(output_dir,samples[i], "_",region,"_norm_counts.bedGraph", sep = ""),
                  quote =FALSE, col.names = FALSE, row.names = FALSE)
    }
    hmm_input <- data.frame(counts = c(norm_counts_log), distance = rep(dist_log, length(samples)), row.names = NULL)
    return(list(window_counts=window_counts,hmm_input = hmm_input, norm_counts = norm_counts,norm_counts_log=norm_counts_log, dist_log = dist_log))
  }
}
