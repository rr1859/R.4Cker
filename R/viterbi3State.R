trimOtherStates <- function(hi, li, ni, window_counts){
  other_states <- c(li, ni)
  left <- hi-1
  right <- hi+1
  left_conflict <- left[which(left %in% other_states)]
  right_conflict <- right[which(right %in% other_states)]
  output <- NULL
  for(i in hi){
    if((i-1) %in% left_conflict & !(i+1) %in% right_conflict)
      output <- rbind(output, data.frame(chr = window_counts[(i+1),1],start = window_counts[(i+1),2], end= window_counts[i,3]))
    if(!(i-1) %in% left_conflict & !(i+1) %in% right_conflict)
      output <- rbind(output, data.frame(chr = window_counts[i,1], start = window_counts[i,2], end = window_counts[i,3]))
  }
  if(!is.null(output))
    output[,2:3] = round(output[,2:3])
  return(output)
}
viterbi3State <- function(hmm_input,mod, samples, window_counts, num_windows, bait_name, bait_chr,reps, output_dir, region){
  pars <- c(unlist(getpars(mod)))

  ntimes <- rep(nrow(hmm_input)/reps,reps)
  if(region == "cis" | region == "nearbait"){
    glmformula <- formula(counts ~ distance)
    mod_new <- depmix(response = glmformula,
                      data = hmm_input,
                      nstates = 3,
                      trstart = pars[4:12],
                      family = gaussian(link = "identity"),
                      respstart = pars[13:21],
                      instart = pars[1:3],
                      ntimes = ntimes
    )
  }
  if(region == "trans"){
    glmformula <- formula(counts ~ 1)
    mod_new = depmix(response = glmformula,
                     data = hmm_input,
                     nstates = 3,
                     trstart = pars[4:12],
                     family = gaussian(link = "identity"),
                     respstart = pars[13:18],
                     instart = pars[1:3],
                     ntimes = ntimes
    )
  }
  v <- viterbi(mod_new)
  row <- 1
  sample_hiint_all <- NULL
  sample_liint_all <- NULL
  sample_niint_all <- NULL
  chrs = unique(window_counts[,1])
  for(rep in 1:reps){
    sample_v <- v[row:(row+num_windows-1),1]
    sample_hiint <- which(sample_v == 3)
    sample_liint <- which(sample_v == 2)
    sample_niint <- which(sample_v == 1)
    sample_hiint_all <- append(sample_hiint_all, sample_hiint)
    sample_liint_all <- append(sample_liint_all, sample_liint)
    sample_niint_all <- append(sample_niint_all, sample_niint)
    trim_sample = NULL
    for(chr in chrs){
      rows = which(window_counts[,1] == as.character(chr))
      if(length(rows) > 0){
        sample_hiint_chr = sample_hiint[which(sample_hiint >= min(rows) & sample_hiint <= max(rows))]
        sample_liint_chr = sample_liint[which(sample_liint >= min(rows) & sample_liint <= max(rows))]
        sample_niint_chr = sample_niint[which(sample_niint >= min(rows) & sample_niint <= max(rows))]
        trim = trimOtherStates(sample_hiint_chr, sample_liint_chr, sample_niint_chr,window_counts)
        trim_sample = rbind(trim_sample, trim)
      }
    }
    write.table(merge_windows(trim_sample),paste(output_dir,samples[rep], "_", region, "_highinter.bed", sep =""),
                quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
    if(length(sample_liint) > 0){
      write.table(merge_windows(window_counts[sample_liint,1:3]),paste(output_dir,samples[rep], "_", region, "_lowinter.bed", sep =""),
                  quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
    }
    if(length(sample_niint) >0){
      write.table(merge_windows(window_counts[sample_niint,1:3]),paste(output_dir,samples[rep], "_", region, "_noninter.bed", sep =""),
                  quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
    }

    row <- row+num_windows
  }


  union_hiint <- unique(sample_hiint_all)
  union_liint <- unique(sample_liint_all)
  union_niint <- unique(sample_niint_all)
  #hi-int overlap between replicates
  table_intersect <- data.frame(table(sample_hiint_all))
  intersect <- table_intersect[table_intersect[,2] == reps,1]
  intersect <- as.numeric(as.character(intersect))
  intersect <- intersect[order(intersect)]
  final_intersect_trim_hi = NULL
  for(chr in chrs){
    rows = which(window_counts[,1] == as.character(chr))
    if(length(rows) > 0){
      intersect_chr = intersect[which(intersect >= min(rows) & intersect <= max(rows))]
      union_liint_chr = union_liint[which(union_liint >= min(rows) & union_liint <= max(rows))]
      union_niint_chr = union_niint[which(union_niint >= min(rows) & union_niint <= max(rows))]
      intersect_trim = trimOtherStates(intersect_chr, union_liint_chr, union_niint_chr,window_counts)
      final_intersect_trim_hi = rbind(final_intersect_trim_hi, intersect_trim)
    }
  }
  write.table(merge_windows(final_intersect_trim_hi),paste(output_dir,bait_name, "_", region, "_highinter.bed", sep =""),
              quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  #low-int overlap between replicates
  table_intersect <- data.frame(table(sample_liint_all))
  intersect <- table_intersect[table_intersect[,2] == reps,1]
  intersect <- as.numeric(as.character(intersect))
  intersect <- intersect[order(intersect)]
  final_intersect_trim_li = NULL
  for(chr in chrs){
    rows = which(window_counts[,1] == as.character(chr))
    if(length(rows) > 0){
      intersect_chr = intersect[which(intersect >= min(rows) & intersect <= max(rows))]
      union_hiint_chr = union_hiint[which(union_hiint >= min(rows) & union_hiint <= max(rows))]
      union_niint_chr = union_niint[which(union_niint >= min(rows) & union_niint <= max(rows))]
      intersect_trim = trimOtherStates(intersect_chr, union_hiint_chr, union_niint_chr,window_counts)
      final_intersect_trim_li = rbind(final_intersect_trim_li, intersect_trim)
    }
  }
  write.table(merge_windows(final_intersect_trim_li),paste(output_dir,bait_name, "_", region, "_lowinter.bed", sep =""),
              quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  #non-int overlap between replicates
  if(length(sample_niint_all) > 0){
    table_intersect <- data.frame(table(sample_niint_all))
    intersect <- table_intersect[table_intersect[,2] == reps,1]
    if(length(intersect) >0){
      intersect <- as.numeric(as.character(intersect))
      intersect <- intersect[order(intersect)]
      final_intersect_trim_ni = NULL
      for(chr in chrs){
        rows = which(window_counts[,1] == as.character(chr))
        if(length(rows) > 0){
          intersect_chr = intersect[which(intersect >= min(rows) & intersect <= max(rows))]
          union_hiint_chr = union_hiint[which(union_hiint >= min(rows) & union_hiint <= max(rows))]
          union_liint_chr = union_liint[which(union_liint >= min(rows) & union_liint <= max(rows))]
          intersect_trim = trimOtherStates(intersect_chr, union_hiint_chr, union_liint_chr,window_counts)
          final_intersect_trim_ni = rbind(final_intersect_trim_ni, intersect_trim)
        }
      }
      if(nrow(final_intersect_trim_ni) > 0){
        write.table(merge_windows(final_intersect_trim_ni),paste(output_dir,bait_name, "_", region, "_noninter.bed", sep =""),
                    quote = FALSE, row.names = FALSE, col.names = FALSE)
      }

    }
  }
}




