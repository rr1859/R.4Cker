differentialAnalysis <- function(obj,norm_counts_avg, windows,conditions, region,coordinates=NULL, pval){
  if(.Platform$OS.type=="windows"){
    quartz<-function() windows()
  }
  pval_options=c(0.01, 0.05,0.1)
  if(length(conditions) != 2)
    stop("Only 2 conditions can be analyzed")
  
  colData_df = data.frame(unlist(lapply(1:length(obj@conditions), function(j) c(rep(obj@conditions[j], obj@replicates[j])))))
  
  rownames(colData_df) = colnames(windows)[-c(1:4)]
  colnames(colData_df) = c("condition")
  
  condition = factor(unlist(lapply(1:length(obj@conditions), function(j) c(rep(obj@conditions[j], obj@replicates[j])))))
  #DESeq steps
  dds = DESeqDataSetFromMatrix(countData=windows[, -c(1:4)],
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
  sig_rows = rep("not_sig", nrow(windows))
  sig_rows[which(res$padj < pval)] = "sig"
  if(region == "nearbait"){
    plot_df = data.frame(coord=c(rowMeans(windows[,2:3]),rowMeans(windows[,2:3])),
                         counts=c(rowMeans(norm_counts[,cols_conditions[condition1_row,1]:cols_conditions[condition1_row,2]]),
                                  rowMeans(norm_counts[,cols_conditions[condition2_row,1]:cols_conditions[condition2_row,2]])),
                         conditions=c(rep(conditions[1], nrow(windows)), rep(conditions[2], nrow(windows))),
                         sig=c(sig_rows, sig_rows))
    if(!is.null(coordinates)){
      plot_df=plot_df[which(plot_df[,1] >= coordinates[1] & plot_df[,1] <= coordinates[2]),]
    }
    
    quartz()
    print(ggplot(plot_df, aes(x=coord, y=counts, colour=conditions))+geom_line()+
      theme_bw()+
      xlab(paste("Chromosome coordinates (", obj@bait_chr, ")", sep =""))+
      ylab("Normalized counts")+geom_point(data=subset(plot_df,sig=="not_sig"), shape=1, size=0.5)+
      geom_point(data=subset(plot_df,sig=="sig")))
  }
  if(region == "cis"){
    plot_df = data.frame(coord=c(rowMeans(windows[,2:3]),rowMeans(windows[,2:3])),
                         counts=c(rowMeans(norm_counts_log[,cols_conditions[condition1_row,1]:cols_conditions[condition1_row,2]]),
                                  rowMeans(norm_counts_log[,cols_conditions[condition2_row,1]:cols_conditions[condition2_row,2]])),
                         conditions=c(rep(conditions[1], nrow(windows)), rep(conditions[2], nrow(windows))),
                         sig=c(sig_rows, sig_rows))
    quartz()
    print(ggplot(plot_df, aes(x=coord, y=counts, colour=conditions))+
      theme_bw()+
      xlab(paste("Chromosome coordinates (", obj@bait_chr, ")", sep =""))+
      ylab("Normalized counts")+geom_point(data=subset(plot_df,sig=="not_sig"), shape=1, size=0.5)+
      geom_point(data=subset(plot_df,sig=="sig")))
  }
  sig_merge_windows=merge_windows(windows[which(res$padj < pval), 1:3])
  print(paste("BED file of significant domains saved in ", obj@output_dir, sep = ""))
  write.table(sig_merge_windows,
    paste(obj@output_dir, obj@bait_name, "_", conditions[1], "_", conditions[2], "_", region, "_pval", pval,"_diff.bed", sep = ""),
    quote=FALSE, col.names=FALSE, row.names=FALSE, sep = "\t")
  op.df <- cbind(windows[,1:3],data.frame(res))
  op.l <- list(op.df, plot_df)
  op.l
}

