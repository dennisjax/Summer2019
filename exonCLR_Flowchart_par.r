###########################################################################
####################   Load exonCLR library  ##############################
###########################################################################
setwd("./")
library("rtracklayer")
library("stringr")
library("GenomicRanges")
library("AnnotationHub")
library("openxlsx")
library("mgcv")
library("progress")
library("dplyr")
library("reshape2")
library("ggplot2")
library("qpcR")
library("foreach")
library("doParallel")
source("exonCLR.r")

registerDoParallel(detectCores() - 2)

###########################################################################
##################### Extract Individual GTFs #############################
###########################################################################
### Input file:
### 1) GTF file from LyRICS
### It cointains all genome annotation for all the target genes
### 2) Target gene list
### The data has three columns: "Gene_stable_ID", "Transcript_stable_ID”, "HGNC_symbol"
### Output:
### Individual GTF files + a pooled Targets_raw.gtf in raw_gtfs/ folder

Extract_individual_gtfs("./annotated.dontScreenMatchAnnot.dontScreenProteinCoding.gtf","genes.txt")
# Modified to test all genes
# Extract_individual_gtfs("./annotated.dontScreenMatchAnnot.dontScreenProteinCoding.gtf")

target0_gtfs<-paste("raw_gtfs", dir("raw_gtfs"), sep = "/")
target0_gtfs<-target0_gtfs[target0_gtfs!="raw_gtfs/Targets_raw.gtf"]
target1_gtfs<-paste("raw_gtfs", dir(pattern="_remain2.0.gtf"), sep = "/")
target1_gtfs<-gsub("_remain2.0","_raw",target1_gtfs)
target_gtfs<-setdiff(target0_gtfs,target1_gtfs)

############################################################################
####################### intron_filter_by_GRanges ###########################
############################################################################
### Filter the transcripts with partial exons in comparison to the Gencode Reference Genome
### Input individual GTFs
### Output files:
### "$Gene"_filtered2.0.gtf 
### "$Gene"_potential_new_exon.gtf
### "$Gene"_remain2.0.gtf
### gencode_"$Gene".gtf
###
### temp_result contains:
### $all_exon
### $prefiltered_new_exon
### $filtered
### $prefiltered_new_exon_count
### $filtered_binary
### $remain
### $filtered_count
### $remain_binary
### $isoform_count
### $remain_count
### $prefiltered_binary

for (in_gtf in target_gtfs){
p.inv <- try(assign(paste("temp_result",strsplit(x=in_gtf,split="/|_")[[1]][3],sep="_"),intron_filter_by_GRanges(import.gff(in_gtf))), silent=TRUE)
if ('try-error' %in% class(p.inv)) next
}
save.image(file="exonCLR_Projects.RData")

#############################################################################
#####################  intron_filter_by_gff  ################################
#############################################################################
### Systemetically go through the process of intron_filter_by_GRanges of all selected genes in a gtf file.
### Beside the four types of gtf files for each gene in the same directory
### It creates a 
### 1) subdirectory called /binary_results;
### All the binary forms are stored as csv files in this subdirectory ("$Gene"_binary.csv).
### 2) filter_stats.csv
###
### The variable temp_gff_result contains:
### $analyzed_gene
### $percent_remain
### $stat_remain
### $percent_filtered
### $stat_filtered        
### $percent_prefiltered
### $stat_prefiltered

for (i in 1:length(target_gtfs)){
x.inv <- try(intron_filter_by_gff(target_gtfs[i]), silent=TRUE)
    if ('try-error' %in% class(x.inv)) next
}
save.image(file="exonCLR_Projects.RData")

############################################################################################# Cluster the isoforms based on their binary forms. ################
################################################################################
### "index" is a vector of integers that represent the group that an isoform belongs to.
### "exon_binary" represents the simplified binary forms after the clustering.
### "group" represents a list containing isoform names in each subgroup.

Files_temp<-ls(pattern="temp_")
for (result_tmp in Files_temp){
y.inv <- try(assign(paste("group_result_",strsplit(x=result_tmp,split="_")[[1]][3],sep=""),isoform_group(get(result_tmp)$remain_binary, get(result_tmp)$remain)), silent=TRUE)
if ('try-error' %in% class(y.inv)) next
}
save.image(file="exonCLR_Projects.RData")

###############################################################################
#################### generate_exon_only_binary  ###############################
###############################################################################
### It generate 
### 1) "$Gene"_exon_only_binary.csv contains exon binary matrix for the gene
### 2) Skipped_Genes.txt

generate_exon_only_binary("./binary_results", "grouped_remain_binary_results")
save.image(file="exonCLR_Projects.RData")

###############################################################################
################## Generate heatmaps of exon coverage #########################
###############################################################################
### Generate a heatmap based on the "_exon_only_binary.csv" files,
### that contain simplified binary forms of transcripts that remained.
### The heatmap is saved as png files in the /output_heatmaps/ folder

library(pheatmap)
original_dir <- getwd()
setwd("./binary_results/grouped_remain_binary_results")
file_list <- list.files(pattern = "_exon_only_binary.csv")
for (i in 1:length(file_list)){
  temp_final <- read.csv(file_list[i])
  gene = strsplit(file_list[i], split = "_exon_only_binary.csv")[[1]]
  iso_names = names(temp_final)[-1]
  iso_index <- as.numeric(temp_final[1, ][-1])
  file_path <- file.path(original_dir, paste(gene, "_remain2.0.gtf", sep = ""))
  temp = import.gff(file_path) 
  print(paste("Analyzing Gene #", i, ":", gene))
  if(nrow(temp_final) > 2){
    heatmap_coverage(temp, getwd(), iso_index)
  }
}
setwd(original_dir)
save.image(file="exonCLR_Projects.RData")

##############################################################################
##################### Simplify a gtf/gff file ################################
##############################################################################
### Simplify a gtf/gff file based on the process from previous functions which
### 1) filter transcripts with partial exons and intron retentions, and 
### 2) cluster the rest into simplified binary.
### The simplified gtf file for each gene is stored in the subdirectory /processed_gtf_files/
### The final output will be a gtf file with an additional prefix of "simplified_" in the working directory

load("exonCLR_Projects.RData")
system(paste("cat *_remain2.0.gtf >", "Targets_remain.gtf"))
gtf_remain<-dir(pattern="_remain2.0.gtf")
for (individual_remain in gtf_remain){
z.inv <- try(simplify_gff(individual_remain), silent=TRUE)
if ('try-error' %in% class(z.inv)) next
}

##############################################################################
######################## Group gtf files #####################################
##############################################################################
### Next step of grouping after simplify_gff
### Change to the directory when contains all the "_simplified.gtf" files.
### Output a csv file in the format of xlsx.
### Row names are the gene names. The column names are the group index.
### The column contains either NA or the transcript IDs in that group for a gene
### The output file name is "result.csv".

output_folder<-"simplified_gtf"
system(paste0("cat simplified_gtf/*_simplified.gtf >", "Targets_simplified.gtf"))
group_gff("Targets_simplified.gtf",custom_gene=intersect(gsub("temp_result_","",ls(pattern="temp_result_")),gsub("_simplified.gtf","",dir("simplified_gtf/"))))
save.image(file="exonCLR_Projects.RData")

################################################################################
########################### Quantify isoforms #################################
###############################################################################
### Quantify the number of isoforms based on the xlsx/csv file generated in the previous function.
### The output is a csv file named "quantify_NA.csv"
### which 0 appearance is replaced by NA to facilitate the bubbleplot in the next step.

quantify("./result.csv", "./all_samples.chained_count.txt")

###############################################################################
######################### bubblePlots of Isoforms #############################
### Create a bubble plot based on quantification for a targeted gene.
### Create a subdirectory called /bubbleplot/ and output the bubble plot there.
### Be careful about the sample selection in the script

foreach(gene_bubble = gsub("_simplified.gtf","",dir("simplified_gtf/"))) %dopar% {
m.inv <- try(bubblePlots("quantity_NA.csv",plot_gene_name=gene_bubble), silent=TRUE)
if ('try-error' %in% class(m.inv)) next
}

stopImplicitCluster()