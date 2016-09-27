differentialAnalysis <- function(obj,norm_counts_avg, windows,conditions, region,coordinates, pval){
  pval_options=c(0.01, 0.05,0.1)
  if(length(conditions) != 2)
    stop("Only 2 conditions can be analyzed")
  all_highinter_domains = NULL
  for(condition in conditions){
    condition_highinter=read.table(paste(obj@output_dir,obj@bait_name, "_", condition, "_", region,"_highinter.bed", sep = ""), stringsAsFactors = FALSE)
    all_highinter_domains = rbind(all_highinter_domains, condition_highinter)
  }
  windows_domains_all = NULL
  for(i in 1:nrow(all_highinter_domains)){
    windows_domain = windows[which(windows[,2] > all_highinter_domains[i,2] & windows[,2] < all_highinter_domains[i,3]),]
    windows_domains_all = rbind(windows_domains_all, windows_domain)
  }
  windows_domains_all = unique(windows_domains_all)
  if(!is.null(coordinates)){
    windows_domains_all=windows_domains_all[which(windows_domains_all[,2] >= coordinates[1] & windows_domains_all[,2] <= coordinates[2]),]
    norm_counts_avg=norm_counts_avg[which(windows_domains_all[,2] >= coordinates[1] & windows_domains_all[,2] <= coordinates[2]),]
    windows=windows[which(windows_domains_all[,2] >= coordinates[1] & windows_domains_all[,2] <= coordinates[2]),]
  }
  window_counts = data.frame(windows_domains_all[,-c(1:4)])
  
  colData_df = data.frame(unlist(lapply(1:length(obj@conditions), function(j) c(rep(obj@conditions[j], obj@replicates[j])))))
  
  rownames(colData_df) = colnames(window_counts)
  colnames(colData_df) = c("condition")
  
  condition = factor(unlist(lapply(1:length(obj@conditions), function(j) c(rep(obj@conditions[j], obj@replicates[j])))))
  #DESeq steps
  dds = DESeqDataSetFromMatrix(countData=window_counts,
                               colData = colData_df,
                               design = ~ condition)
  
  dds <- estimateSizeFactors(dds)
  dds <- estimateDispersions(dds, fitType = "local")
  dds <- nbinomWaldTest(dds)
  res <- results(dds, c("condition", conditions[1], conditions[2]))
  if(length(which(res$padj < pval)) == 0){
    print("No significant changes..trying a higher p-value...")
    pval_row=which(pval_options == pval)
    while(!length(which(res$padj < pval)) == 0 & pval_row <=3){
      pval_row=pval_row+1
    }
    if(length(which(res$padj < pval_options[pval_row])) == 0)
      stop("No significant interactions")
    else{
      pval=pval_options[pval_row]
      print("Using ",pval , "...")
    }
  }
  norm_counts = counts(dds, normalized = TRUE)
  norm_counts_log = log(norm_counts+1,10)
  condition1_row = which(obj@conditions == conditions[1])
  condition2_row = which(obj@conditions == conditions[2])
  cols_conditions=NULL
  j=1
  for(i in obj@replicates){
    cols_conditions = rbind(cols_conditions, c(j, (j+i-1)))
    j=j+i
  }
  if(region == "cis"){
    sig_rows = rep("not_sig", nrow(windows_domains_all))
    sig_rows[which(res$padj < pval)] = "sig"
    plot_df = data.frame(coord=c(rowMeans(windows_domains_all[,2:3]),rowMeans(windows_domains_all[,2:3])),
                         counts=c(rowMeans(norm_counts_log[,cols_conditions[condition1_row,1]:cols_conditions[condition1_row,2]]),
                                  rowMeans(norm_counts_log[,cols_conditions[condition2_row,1]:cols_conditions[condition2_row,2]])),
                         conditions=c(rep(conditions[1], nrow(windows_domains_all)), rep(conditions[2], nrow(windows_domains_all))),
                         sig=c(sig_rows, sig_rows))
    quartz()
    ggplot(plot_df, aes(x=coord, y=counts, colour=conditions))+
      theme_bw()+
      xlab(paste("Chromosome coordinates (", obj@bait_chr, ")", sep =""))+
      ylab("Normalized counts")+geom_point(data=subset(plot_df,sig=="not_sig"), shape=1, size=0.5)+
      geom_point(data=subset(plot_df,sig=="sig"))
  }
  if(region == "nearbait"){
    sig_rows = rep("not_sig", nrow(windows))
    sig_windows=paste(windows_domains_all[which(res$padj < pval), 1], 
                      windows_domains_all[which(res$padj < pval), 2], 
                      windows_domains_all[which(res$padj < pval), 3], sep = "_")
    sig_windows_rows=rep("not_sig", nrow(windows))
    sig_windows_rows[which(paste(windows[,1], windows[,2], windows[,3], sep = "_") %in% sig_windows)] = "sig"
    plot_df = cbind(norm_counts_avg,sig=c(sig_windows_rows, sig_windows_rows))
    plot_df_log= plot_df
    plot_df_log[,2] = log(plot_df_log[,2]+1,10)
    #left of here
    
    quartz()
    print(ggplot(plot_df, aes(x=Coord, y=Count, colour=Condition))+
      geom_line()+theme_bw()+
      xlab(paste("Chromosome coordinates (", obj@bait_chr, ")", sep =""))+
      ylab("Normalized counts")+geom_point(data=subset(plot_df,sig=="not_sig"), shape=1, size=0.5)+
      geom_point(data=subset(plot_df,sig=="sig")))
    quartz()
    print(ggplot(plot_df_log, aes(x=Coord, y=Count, colour=Condition))+
      theme_bw()+
      xlab(paste("Chromosome coordinates (", obj@bait_chr, ")", sep =""))+
      ylab("Normalized counts(log)")+geom_point(data=subset(plot_df_log,sig=="not_sig"), shape=1, size=0.5)+
      geom_point(data=subset(plot_df_log,sig=="sig")))
  }
  sig_merge_windows=merge_windows(windows_domains_all[which(res$padj < pval), 1:3])
  print(paste("BED file of significant domains saved in ", obj@output_dir, sep = ""))
  write.table(sig_merge_windows,
              paste(obj@output_dir, obj@bait_name, "_", conditions[1], "_", conditions[2], "_", region, "_pval", pval,"_diff.bed", sep = ""),
              quote=FALSE, col.names=FALSE, row.names=FALSE, sep = "\t")
}

