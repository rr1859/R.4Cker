buildAdaptiveWindowsCis <- function(data_all,bait_coord, bait_chr,bait_name, output_dir,k, region){
  print("Building adaptive windows...")
  bw_all <- NULL
  for(i in 1:length(data_all)){
    bait_row <- tail(which(data_all[[i]][,2] <= bait_coord),1)
    coord <- data_all[[i]][,2]
    start <- 1
    end <- length(coord)
    bw <- rep(0, end)
    #i takes the value of all coord left of the bait and j right of the bait
    #for trans bait_row is set to zero
    i <- bait_row
    j <- bait_row
    while(i >= start | j <= end){
      if(i >= start){
        if(i-k >= 1)
          bw[i] <- coord[i]-coord[(i-k)]
        else
          bw[i] <- coord[k]-coord[1]
        i <- i-1
      }
      if(j <= end){
        if (j+k <= end)
          bw[j] <- coord[(j+k)]-coord[j]
        else
          bw[j] <- coord[end]-coord[(end-k)]
        j <- j+1
      }
    }
    bw_all <- rbind(bw_all, cbind(coord,rep(paste("sample", i, sep = ""),length(bw)),bw))
  }
  bw_all <- data.frame(bw_all)
  bw_all[,1] <- as.numeric(as.character(bw_all[,1]))
  bw_all[,3] <- as.numeric(as.character(bw_all[,3]))
  bw_all <- bw_all[order(bw_all[,1]),]
  #cubic spline fit to coord, window size
  s_fit <- smooth.spline(bw_all[,1], bw_all[,3], spar = 0.75)
  spline <- cbind(s_fit$x,rep("spline",length(s_fit$x)),s_fit$y)
  colnames(spline) <- c("coord", "V2", "bw")
  #bw_all_2 for plot
  bw_all_2 <- rbind(bw_all,spline)
  bw_all_2[,1] <- as.numeric(as.character(bw_all_2[,1]))
  bw_all_2[,3] <- as.numeric(as.character(bw_all_2[,3]))
  start <- min(bw_all[,1])
  end <- max(bw_all[,1])
  left <- bait_coord
  right <- bait_coord
  window_counts <- NULL
  left_window_sizes <- NULL
  right_window_sizes <- NULL
  #max window size is 500kb
  max_ws <- 500000
  while(left >= start | right <= end){
    if(left >= start){
      #print(paste("left:",left,sep=""))
      if(predict(s_fit,left)$y <5){
        print("Try more fragments per window")
        break
      }
      else{
        if(is.null(left_window_sizes)){
          size <- predict(s_fit,left)$y
        }
        else{
          if(predict(s_fit,left)$y >= 1.5*max(left_window_sizes)){
            #print("left_1")
            size <- max(left_window_sizes)}
          if(predict(s_fit,left)$y <= median(left_window_sizes))
            size <- median(left_window_sizes)
          if(predict(s_fit,left)$y > max_ws)
            size <- max_ws
          else
            size <- predict(s_fit,left)$y
        }
        #print(paste("left_size",size))
        window_left <- c(left-size, left)
        window_counts <- rbind(window_counts, window_left)
        left_window_sizes <- append(left_window_sizes, size)
        left <- mean(window_left)
      }
    }
    if(right <= end){
      #print(paste("right:",right,sep=""))
      if(predict(s_fit,right)$y < 5){
        #print("Try more fragments per window")
        break
      }
      else{
        if(is.null(right_window_sizes))
          size <- predict(s_fit,right)$y
        else{
          if(predict(s_fit,right)$y >= 1.5*max(right_window_sizes))
            size <- max(right_window_sizes)
          if(predict(s_fit,right)$y <= median(right_window_sizes))
            size <- median(right_window_sizes)
          if(predict(s_fit,right)$y > max_ws)
            size <- max_ws
          else
            size <- predict(s_fit,right)$y
        }
        #print(paste("right_size",size))
        window_right <- c(right,right+size)
        window_counts <- rbind(window_counts, window_right)
        right_window_sizes <- append(right_window_sizes,size)
        right <- mean(window_right)
      }
    }
  }
  #remove - rows
  neg_rows <- which(window_counts[,1] <=0 | window_counts[,2] <=0)
  if(length(neg_rows) != 0)
    window_counts <- window_counts[-neg_rows,]
  window_counts <- window_counts[order(window_counts[,2]),]
  window_counts[,1] <- round(window_counts[,1])
  window_counts[,2] <- round(window_counts[,2])
  num_windows <- nrow(window_counts)
  write.table(cbind(rep(bait_chr, num_windows), window_counts), paste(output_dir,bait_name, "_",region,"_adaptive_windows.bed",sep =""), quote = FALSE, row.names = FALSE, col.names = FALSE)
  window_counts <- cbind(rep(bait_chr, nrow(window_counts)), window_counts,abs(rowMeans(window_counts[,1:2])-bait_coord))
  window_counts <- data.frame(window_counts, row.names = NULL)
  window_counts[,2] <- as.numeric(as.character(window_counts[,2]))
  window_counts[,3] <- as.numeric(as.character(window_counts[,3]))
  window_counts[,4] <- as.numeric(as.character(window_counts[,4]))
  colnames(window_counts) = c("chr", "start", "end", "distance")
  return(window_counts)
}


