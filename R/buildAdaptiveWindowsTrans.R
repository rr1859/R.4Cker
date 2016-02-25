bwMax <- function(coord,k){
  end <- length(coord)
  output <- rep(0,end)
  #i takes the value of all coord left of the bait and j right of the bait
  for(i in 1:end){
    if(i+k <=end)
      output[i] <- coord[(i+k)]-coord[i]
    else
      output[i] <- coord[end]-coord[(end-k)]
  }
  return(output)
}
overlappingWindows <- function(bw_all, s_fit){
  start <- min(bw_all[,1])
  end <- max(bw_all[,1])
  window_counts <- NULL
  window_sizes <- NULL
  pos <- start
  max_ws <- 1.5e6
  while(pos <= end){
    if(predict(s_fit,pos)$y <5){
      print("Try more fragments per window")
      break
    }
    else{
      if(is.null(window_sizes))
        size <- predict(s_fit,pos)$y
      else{
        if(predict(s_fit,pos)$y >= 1.5*max(window_sizes))
          size <- max(window_sizes)
        if(predict(s_fit,pos)$y <= median(window_sizes))
          size <- median(window_sizes)
        if(predict(s_fit,pos)$y > max_ws)
          size <- max_ws
        else
          size <- predict(s_fit,pos)$y
      }
      #print(paste("right_size",size))
      window <- c(pos,pos+size)
      window_counts <- rbind(window_counts, window)
      window_sizes <- append(window_sizes,size)
      pos <- mean(window)
    }
  }
  #remove - rows
  neg_rows <- which(window_counts[,1] <=0 | window_counts[,2] <=0)
  if(length(neg_rows) != 0)
    window_counts <- window_counts[-neg_rows,]
  window_counts <- window_counts[order(window_counts[,2]),]
  num_windows <- nrow(window_counts)
  #print (paste("Number of windows on cis chr is ", num_windows, sep = ""))
  window_counts[,1] <- round(window_counts[,1])
  window_counts[,2] <- round(window_counts[,2])
  return(window_counts)
}

buildAdaptiveWindowsTrans <- function(data, samples, chrs_trans,k){
  print("Building adaptive windows...")
  reps <- length(samples)
  bw_all <- vector("list", length(chrs_trans))
  windows_all <- NULL
  for(j in 1:length(chrs_trans)){
    bw_reps <- NULL
    for(i in 1:reps){
      trans_chr <- data[[i]]
      coord <- trans_chr[trans_chr[,1] == chrs_trans[j],2]
      #print(paste( chrs_trans[j], length(coord)))
      if(length(coord) > k){
        coord <- coord[order(coord)]
        bw <- bwMax(coord,k)
        bw_reps <- rbind(bw_reps, cbind(coord,rep(paste("sample", i, sep = ""),length(bw)),bw))
        bw_reps <- data.frame(bw_reps)
        bw_reps[,1] <- as.numeric(as.character(bw_reps[,1]))
        bw_reps[,3] <- as.numeric(as.character(bw_reps[,3]))
        bw_reps <- bw_reps[order(bw_reps[,1]),]
        s_fit <- smooth.spline(bw_reps[,1], bw_reps[,3], spar= 0.75)
        spline <- cbind(s_fit$x,rep("spline",length(s_fit$x)),s_fit$y)
        colnames(spline) <- c("coord", "V2", "bw")
        bw_all[[j]] <- bw_reps
        windows <- overlappingWindows(bw_reps, s_fit)
        windows <- cbind(rep(chrs_trans[j], nrow(windows)), windows)
        windows <- data.frame(windows, row.names = NULL)
        windows[,2] <- as.numeric(as.character(windows[,2]))
        windows[,3] <- as.numeric(as.character(windows[,3]))
        windows_all <- rbind(windows_all, windows)
      }
    }
  }
  colnames(windows_all) = c("chr", "start", "end")
  windows_all <- windows_all[order(windows_all[,1], windows_all[,2]),]
  return(windows_all)
}



