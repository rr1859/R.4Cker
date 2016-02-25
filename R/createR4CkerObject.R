#' 4C-ker class
#' @param data_cis list containing counts per fragment on the bait chromosome for each replicate
#' @param data_trans list containing counts per fragment on trans chromosome for each replicate
#' @param data_nearbait containing counts per fragment in the region surrounding the bait for each replicate
#' @param chrs_trans names of tran chromososomes
#' @param bait_name name of the bait as provided by the user
#' @param bait_chr chromosome name that the bait is located on
#' @param bait_coord chromosomal position of bait primer
#' @param samples names of all samples for analysis (including replicates)
#' @param condition names of all conditions for analysis
#' @param replicates number of replicates for each condition
#' @param species species code
#' @param output_dir directory where all files generated will be saved

R4CkerData <- setClass("R4CkerData",
  representation(data_cis = "list",
    data_trans = "list",
    data_nearbait = "list",
    chrs_trans ="vector",
    bait_name = "character",
    bait_chr = "character",
    bait_coord = "numeric",
    primary_enz = "character",
    samples = "vector",
    conditions = "vector",
    replicates = "vector",
    species = "character",
    output_dir = "character"
  )
)
#' Create 4C-ker object
#' @param files path to files used for analysis
#' @param bait_chr chromosome name where bait is located
#' @param bait_coord chromosomal position of bait primer
#' @param bait_name name of the bait
#' @param samples list of samples used in analysis (including replicates)
#' @param condition names of all conditions for analysis
#' @param replicates number of replicates for each condition
#' @param species species code
#' @param output_dir directory where all files generated will be saved

createR4CkerObjectFromFiles <- function(files,bait_chr,bait_coord,bait_name,primary_enz,samples,conditions,
                             replicates, species, output_dir){
  if(sum(replicates) != length(files))
    stop("Number of samples does not match")
  if(substr(output_dir, nchar(output_dir), nchar(output_dir)) != "/")
    output_dir = paste(output_dir,"/", sep = "")
  dir.create(output_dir)
  data_cis <- vector("list", sum(replicates))
  data_nearbait <- vector("list",sum(replicates))
  data_trans <- vector("list",sum(replicates))
  chrs_trans <- NULL
  if(nchar(primary_enz) == 4)
    nearbait_size = 1e6
  else
    nearbait_size = 5e6
  stats <- NULL
  for(i in 1:length(files)){
    data <- read.table(files[i], stringsAsFactors = FALSE)
    data <- data[order(data[,1], data[,2]),]
    data_sample_cis <- data[data[,1] == bait_chr,]
    data_sample_trans <- data[data[,1] != bait_chr,]
    data_cis[[i]] <- data_sample_cis
    data_trans[[i]] <- data_sample_trans
    data_nearbait[[i]] <- data_sample_cis[which(data_sample_cis[,2] > (bait_coord-nearbait_size) &
                                             data_sample_cis[,2] < (bait_coord+nearbait_size)),]
    chrs_trans <- append(chrs_trans,as.character(data[data[,1] != bait_chr,1]))
    total_num_reads <- sum(data_sample_trans[,4], data_sample_cis[,4])
    total_num_sites <- length(c(data_sample_trans[,4], data_sample_cis[,4]))
    num_reads_cis <- sum(data_sample_cis[,4])
    perc_reads_cis <- (num_reads_cis/total_num_reads)*100
    num_sites_cis <- length(data_sample_cis[,4])
    perc_sites_cis <- (num_sites_cis/total_num_sites)*100
    num_reads_trans <- total_num_reads-num_reads_cis
    perc_reads_trans <- (num_reads_trans/total_num_reads)*100
    num_sites_trans <- total_num_sites-num_sites_cis
    perc_sites_trans <- (num_sites_trans/total_num_sites)*100
    stats <- cbind(stats, c(total_num_reads,total_num_sites,
                            num_reads_cis,perc_reads_cis,
                            num_sites_cis,perc_sites_cis,
                            num_reads_trans,perc_reads_trans,
                            num_sites_trans, perc_sites_trans))
  }
  chrs_trans = unique(chrs_trans)
  rownames(stats) <- c("Total_number_of_reads","Total_number_of_observed_fragments",
                       "Reads_in_cis","Percentage_of_reads_in_cis",
                       "Observed_fragments_in_cis","Percentage_of_observed_fragments_in_cis",
                       "Reads_in_trans","Percentage_of_reads_in_trans",
                       "Observed_fragments_in_trans","Percentage_of_observed_fragments_in_trans")
  colnames(stats) <- sub(".bedGraph", "",files)
  write.table(stats, paste(output_dir,bait_name,"_stats.txt", sep = ""), quote = FALSE, sep = "\t")

  obj_4Cker = R4CkerData(data_cis = data_cis,
                      data_nearbait = data_nearbait,
                      data_trans = data_trans,
                      chrs_trans = chrs_trans,
                      bait_name = bait_name,
                      bait_chr = bait_chr,
                      bait_coord = bait_coord,
                      primary_enz = primary_enz,
                      samples = samples,
                      conditions = conditions,
                      replicates = replicates,
                      species = species,
                      output_dir = output_dir)
  return(obj_4Cker)
}

createR4CkerObjectFromDFs <- function(dfs,bait_chr,bait_coord,bait_name,primary_enz,samples,conditions,
                                        replicates, species, output_dir){
  if(sum(replicates) != length(dfs))
    stop("Number of samples does not match")
  if(substr(output_dir, nchar(output_dir), nchar(output_dir)) != "/")
    output_dir = paste(output_dir,"/", sep = "")
  dir.create(output_dir)
  data_cis <- vector("list", sum(replicates))
  data_nearbait <- vector("list",sum(replicates))
  data_trans <- vector("list",sum(replicates))
  chrs_trans <- NULL
  if(nchar(primary_enz) == 4)
    nearbait_size = 1e6
  else
    nearbait_size = 5e6
  stats <- NULL
  for(i in 1:length(dfs)){
    data <- get(dfs[i])
    data <- data[order(data[,1], data[,2]),]
    data_sample_cis <- data[data[,1] == bait_chr,]
    data_sample_trans <- data[data[,1] != bait_chr,]
    data_cis[[i]] <- data_sample_cis
    data_trans[[i]] <- data_sample_trans
    data_nearbait[[i]] <- data_sample_cis[which(data_sample_cis[,2] > (bait_coord-nearbait_size) &
                                                  data_sample_cis[,2] < (bait_coord+nearbait_size)),]
    chrs_trans <- append(chrs_trans,as.character(data[data[,1] != bait_chr,1]))
    total_num_reads <- sum(data_sample_trans[,4], data_sample_cis[,4])
    total_num_sites <- length(c(data_sample_trans[,4], data_sample_cis[,4]))
    num_reads_cis <- sum(data_sample_cis[,4])
    perc_reads_cis <- (num_reads_cis/total_num_reads)*100
    num_sites_cis <- length(data_sample_cis[,4])
    perc_sites_cis <- (num_sites_cis/total_num_sites)*100
    num_reads_trans <- total_num_reads-num_reads_cis
    perc_reads_trans <- (num_reads_trans/total_num_reads)*100
    num_sites_trans <- total_num_sites-num_sites_cis
    perc_sites_trans <- (num_sites_trans/total_num_sites)*100
    stats <- cbind(stats, c(total_num_reads,total_num_sites,
                            num_reads_cis,perc_reads_cis,
                            num_sites_cis,perc_sites_cis,
                            num_reads_trans,perc_reads_trans,
                            num_sites_trans, perc_sites_trans))
  }
  chrs_trans = unique(chrs_trans)
  rownames(stats) <- c("Total_number_of_reads","Total_number_of_observed_fragments",
                       "Reads_in_cis","Percentage_of_reads_in_cis",
                       "Observed_fragments_in_cis","Percentage_of_observed_fragments_in_cis",
                       "Reads_in_trans","Percentage_of_reads_in_trans",
                       "Observed_fragments_in_trans","Percentage_of_observed_fragments_in_trans")
  colnames(stats) <- sub(".bedGraph", "",dfs)
  write.table(stats, paste(output_dir,bait_name,"_stats.txt", sep = ""), quote = FALSE, sep = "\t")
  
  obj_4Cker = R4CkerData(data_cis = data_cis,
                         data_nearbait = data_nearbait,
                         data_trans = data_trans,
                         chrs_trans = chrs_trans,
                         bait_name = bait_name,
                         bait_chr = bait_chr,
                         bait_coord = bait_coord,
                         primary_enz = primary_enz,
                         samples = samples,
                         conditions = conditions,
                         replicates = replicates,
                         species = species,
                         output_dir = output_dir)
  return(obj_4Cker)
}

