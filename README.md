# R.4Cker

4Cker is a method developed to analyze 4C-Seq (circularized chromosome conformation capture) data. The manuscript has been submitted and a preprint is available on bioarxiv (<a href="http://biorxiv.org/content/early/2015/11/03/030569">4C-ker: A method to reproducibly identify genome-wide interactions captured by 4C-Seq
experiments</a>).

For reviews on 4C-seq refer to the following articles:
<ul>
	<li>van de Werken HJ, Methods in Enzymology, 2012</li>
	<li>Splinter E, Methods, 2012</li>
	<li>Raviram R, Epigenomics, 2014</li>
</ul>
4Cker is designed to identify domains of interaction
<ol>
	<li>near the bait (1MB on either side for 4bp enzymes and 5Mb on either side for 6bp enzymes)</li> 
	<li><i>cis</i> chromosome</li>
	<li><i>trans</i> chromosomes</li>
</ol>

The required R packages to run R.4Cker are: MASS, DESeq2, psych, depmixS4, miscTools, devtools, RColorBrewer, ggplot2. To install R.4Cker:
```
library(devtools)
install_github("rr1859/R.4Cker")
library(R.4Cker)
```
The input to 4C-ker are 4 column tab-delimited count files with <i>chr</i> <i>start</i> <i>end</i> and <i>counts per observed fragment</i>. The extension to files can be named ‘.bedGraph’ to visualize the interaction profile on IGV genome browser. If you already have the count files skip to step 3, if not 4C-Seq FASTQ files can be mapped to a reduced genome (refer to steps 1-2). 

**1. Creating a reduced genome and mapping 4C-Seq reads**

4C-Seq reads are usually mapped to a reduced genome consisting of unique sequence fragments adjacent to the primary restriction enzyme sites in the genome. A script (reduced_genome.sh) has been provided to create a custom reduced genome (this requires oligoMatch (ucsc command line tools - kent) and fastaFromBed (bedtools)). Modify line 13,15,17 in reduced_genome.sh to reflect the correct size of fragment, primary enzyme and genome. In addition a FASTA(.fa) file for the primary enzyme and genome are required (files should be named according to what is provided in the .sh script - example: mm10.fa and hindiii.fa). The reduced genome can be used to map 4C-Seq single-end reads with <i>bowtie2</i>. First a bowtie2 index is created using <i>bowtie2-build</i>. With bowtie2, the -5 option can be used to trim the first ‘x’ bps that contain the barcode (if present) and the bait sequence including the RE. Below is an example to map 51bp long reads where the first 26 bps of the read contain a 6bp barcode + 20 bp of the sequence containing the bait sequence:

```
bowtie2-build mm10_hindiii_flanking_sequences_25_unique.fa mm10_hindiii_flanking_sequences_25_unique
bowtie2 -p 12 -N 0 -5 26 \
	--un sample_unaligned.sam \
	-x mm10_hindiii_flanking_sequences_25_unique \
	-U sample.fastq \
	-S sample_aligned.sam
```

**2. Creating a counts file from mapped data (SAM file output from bowtie2)**

The SAM output from bowtie2 reads can be used to create a counts file(.bedGraph) with number of reads per RE fragment using the sam_bedGraph.sh file with the argument for the directory where .SAM files can be found. This will make a .bedGraph file for ALL *aligned.sam files in the folder. Once files are generated, remove the self-ligated and undigested fragments before analyzing 4C interactions.
```
bash sam_bedGraph.sh input_dir_forSAMfiles/
```

**3. Removing the self-ligated and undigested fragments**

To remove the self-ligated and undigested fragments, first find the coordinates based on the primer sequence (Watson strand regardless of primer strand and not including the RE sequence)
```
grep -i -B 1 "CTTCTATCTGCAGAGA" mm10_hindiii_flanking_sequences_25_unique_2.fa

>chr1:100025183-100025208
CTTCTATCTGCAGAGAAATATAGCC
```
Then find the sequence before and after the bait fragment
```
grep -E -A 1 -B 1 'chr1.100025183.100025208' mm10_hindiii_flanking_sites_25_unique_2.bed

chr1	100025152	100025177
chr1	100025183	100025208
chr1	100038369	100038394
```
Remove the fragments matching the coordinates from line 1 and 3 printed above from all .bedGraph file and generate new files with suffix '_rm_self_und.bedGraph'. Remember to do this for each bait.
```
for file in bait*bedGraph; do grep -v -E 'chr1.100025152.100025177|chr1.100038369.100038394' $file > $(echo $file | sed 's/.bedGraph/_rm_self_und.bedGraph/g'); done
```

**4. Analyzing 4C interactions**

As an example, datasets with baits on CD83 (HindIII) and the Tcrb Eb enhancer(NlaIII) will be used to demonstrate the use of 4C-ker. Count files for these datasets can be found in the data folder.

Creating a 4C-ker object: A 4C-ker object can be created using exiting data frames loaded in R (<i>create4CkerObjectFromDFs</i>) or by specifying the path to the .bedGraph files (<i>createR4CkerObjectFromFiles</i>). The function <i>create4CkerObjectFromDFs</i> can be used with the example datasets provided. The arguments for the function are as follows:<br>

<i>files/dfs</i>: path to files or character string of the names of the data.frame variables<br>
<i>bait_chr</i>: chromosome name of bait location<br>
<i>bait_coord</i>: chromosomal coordinate of bait primer<br>
<i>primary_enz</i>: sequence of primary restriction enzyme<br>
<i>samples</i>: names of samples for all files<br>
<i>conditions</i>: names of conditions being analyzed (min of 1)<br>
<i>replicates</i>: number of replicates for each condition<br>
<i>species</i>: code for species (Mouse=mm, Human=hs)<br>
<i>output_dir</i>: path to directory where all output files will be saved

For example to analyze interactions for 1 condition with 2 replicates using the example data (CD83_HindIII):
```
data(CD83_HindIII)
my_obj = createR4CkerObjectFromDFs(dfs = c("CD83_HindIII_1", "CD83_HindIII_2"),
                       bait_chr="chr13",
                       bait_coord= 43773612,
                       bait_name = "CD83",
                       primary_enz = "AAGCTT",
                       samples = c("CD83_H_1", "CD83_H_2"),
                       conditions = "CD83",
                       replicates = 2,
                       species = "mm",
                       output_dir = "~/Documents/CD83_results_R4Cker/")
```

To analyze multiple conditions with the same bait:
```
data(Tcrb_Eb_DN_ImmB_WT)
my_obj = createR4CkerObjectFromDFs(dfs = c("Tcrb_Eb_DN_WT_1","Tcrb_Eb_DN_WT_2","Tcrb_Eb_ImmB_WT_1","Tcrb_Eb_ImmB_WT_2"),
                       bait_chr="chr6",
                       bait_coord= 41553025,
                       bait_name = "Eb",
                       primary_enz = "CATG",
                       samples = c("DN_Eb_1", "DN_Eb_2", "ImmB_Eb_1", "ImmB_Eb_2"),
                       conditions = c("DN", "ImmB"),
                       replicates = c(2,2),
                       species = "mm",
                       output_dir = "~/Documents/Eb_DN_ImmB_results_R4Cker/")
```

To use other datasets from the paper you can download the ZIP of R.4Cker from github and use the .bedGraph files in the data folder.
```
my_obj = createR4CkerObjectFromFiles(files = c("~/Downloads/R.4Cker-master/data/IghCg1_HindIII_1.bedGraph", 
					"~/Downloads/R.4Cker-master/data/IghCg1_HindIII_2.bedGraph"),
                       bait_chr="chr12",
                       bait_coord= 113321402,
                       bait_name = "Igh",
                       primary_enz = "AAGCTT",
                       samples = c("Igh_H_1", "Igh_H_2"),
                       conditions = "Igh",
                       replicates = 2,
                       species = "mm",
                       output_dir = "~/Documents/Igh_results_R4Cker/")
```

Once my_obj has been created, the following functions can be called to define domains of interaction with the bait: <i>nearBaitAnalysis</i>, <i>cisAnalysis</i>, <i>transAnalysis</i>. Each of these functions takes as arguments, my_obj and the value of 'k' to determine the size of the adaptive windows. The BED files with the resulting domains of interaction files (*_highinter.bed) will be saved in the output directory. There will be one file for each replicate and one that contains domains called by all replicates. In addition, a normalized counts file will be created for the adaptive windows and can be viewed on IGV genome browser.

Near the bait analysis: For 4bp cutter experiments we recommend a k=3-5 and for 6bp cutter experiments a k=10. The function returns a data frame of windows used for the analysis and the mean window counts for each conditions.
```
nb_results=nearBaitAnalysis(my_obj,k=5)
```
To plot the 4C profile for all conditions near the bait:
```
library(ggplot2)
ggplot(nb_results$norm_counts_avg, aes(x=Coord, y=Count, colour=Condition))+
	theme_bw()+
	geom_line()+xlab(paste("Chromosome coordinates (", my_obj@bait_chr, ")", sep =""))+
	ylab("Normalized counts")+
	ggtitle(paste("Near bait analysis (", my_obj@bait_name, " bait)", sep = ""))
```
The X and Y axis limits can be changed to zoom into a particular region. Add to the end of ggplot function: +ylim(c(start,end))+xlim(c(start,end))

<i>Cis</i> analysis: We recommend a k=10 for all experiments.
```
cis_results=cisAnalysis(my_obj,k=10)
```
<i>Trans</i> analysis. We recommend a k=20 for 6bp experiments and a k=50+ for 4bp cutter experiments.
```
transAnalysis(my_obj,k=20)
```

**4. Differential interactions using DESeq2**

Differential analysis is done between two conditions on the windows that are interacting with the bait in at least 2 condition. We recommend this analysis only for nearbait and cis interactions. A plot of significant interactions will be generated and a BED file of significantly different domains will be saved in the output directory. The arguments for the function are as follows:<br>
<i>obj</i>: 4CkerObject created for analysis<br>
<i>norm_counts_avg</i>: Average counts of replicates for each condition per window<br>
<i>windows</i>: Raw counts per window for each sample<br>
<i>conditions</i>: Conditions to compare - only two conditions for each comparison<br>
<i>region</i>: Region for analysis.options are: "nearbait" or "cis"<br>
<i>coordinates</i>: Coordinates for a more focused analysis. example: coordinates=c(41e6,42e6)<br>
<i>pval</i>: p-value cut off for DESeq2 analysis. Recommended values are 0.01,0.05 or 0.1. If given p-value fails, the function will try higher values upto 0.1<br>

Example for nearbait analysis with Tcrb dataset:

```
differentialAnalysis(obj=my_obj,
			norm_counts_avg=nb_results$norm_counts_avg,
			windows=nb_results$window_counts,
			conditions=c("DN", "ImmB"),
			region="nearbait",
			coordinates=NULL,
			pval=0.05)
```
