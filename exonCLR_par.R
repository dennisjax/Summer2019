###########################################################################################
######################  Extract GTF files for individual genes  ###########################
###########################################################################################
### Input file:
### 1) GTF file from LyRICS
### It cointains all genome annotation for all the target genes
### 2) Target gene list
### The data has three columns: "Gene_stable_ID", "Transcript_stable_ID”, "HGNC_symbol"
### Output:
### Individual GTF files + a pooled Targets_raw.gtf in raw_gtfs/ folder
require(foreach)
Extract_individual_gtfs <- function(raw_gtf, target_gene_list = NULL){
dir.create("./raw_gtfs",showWarnings = TRUE)
library(rtracklayer)
## Input the gtf from LyRICS pipeline
gtf_raw<-import.gff(raw_gtf)
# Input target genes
if(is.null(target_gene_list)){
  target_genes <- as.data.frame(unique(gtf_raw$gffcompare_gene_name[!is.na(gtf_raw$gffcompare_gene_name)]))
}
else{
  target_genes <- read.table(target_gene_list)
}
print(target_genes)
colnames(target_genes)=c("Genes")
Transcripts<-unique(gtf_raw$transcript_id[gtf_raw$gffcompare_gene_name %in% target_genes$Genes])
target_gtf<-gtf_raw[gtf_raw$transcript_id %in% Transcripts]
setwd("./raw_gtfs")
export.gff(target_gtf,"Targets_raw.gtf",format="gtf")
foreach(gene_symbol = target_genes$Genes) %dopar% {
ensembl_transcript<-unique(gtf_raw$transcript_id[gtf_raw$gffcompare_gene_name %in% gene_symbol])
individual_gtf<-gtf_raw[gtf_raw$transcript_id %in% ensembl_transcript]
if(length(individual_gtf)>0){
export.gff(individual_gtf,paste(gene_symbol,"raw.gtf",sep="_"),format="gtf")
}
else{
next
}
}
setwd("..")
}

bubblePlots <- function(csv, plot_gene_name){
  require(dplyr)
  require(reshape2)
  require(ggplot2)
  dir.create("bubbleplot", showWarnings = F)
  res <- read.csv(csv, header = T, stringsAsFactors = F)

  gene <- res %>% dplyr::filter(gene_name %in% c(plot_gene_name))
  
  gene.melt <- melt(gene, id=c("gene_name", "group"))
  
  
  gene.melt$label <- paste0(gene.melt$gene_name, "_", sprintf("%2d", gene.melt$group))
########################################################################  
  targetSamples <- grep("MEL", unique(gene.melt$variable), value = T)
  gene.melt.filt <- gene.melt %>% dplyr::filter(variable %in% targetSamples)
  
  sorted_labels <- unique(gene.melt.filt$label[order(gene.melt.filt$label, decreasing = T)])
  
  gene.melt.filt$label <- factor(gene.melt.filt$label,
                                 levels = sorted_labels)
  
  ggplot(gene.melt.filt, 
         aes(x=variable, y=label, color = value, 
                             size = value)) +
    geom_point(shape = 16) +
    geom_count() +
    scale_fill_distiller(palette = "RdYlBu") +
    scale_color_distiller(palette = "RdYlBu") +
    theme_bw() +
    theme(axis.text.x = element_text(size=10, angle = 30, hjust = 1),
          axis.text.y = element_text(size=12),
          axis.title.x = element_text(face="plain", colour="black", size=12),
          axis.title.y = element_text(face="plain", colour="black", size=12),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          plot.background = element_blank(),
          legend.position="right") +
    labs(title="", x = "Samples", y="Isoform groups" )
    graph_height <- 0.7*nrow(gene)
    if(graph_height < 7){graph_height <- 6.5}
    ggsave(file.path(getwd(), "bubbleplot", paste(plot_gene_name, "_bubbleplot.png", sep = "")), height = graph_height, dpi = 400)
}

## Data Cleaning function
find_retention <- function(input_binary_mat, head_tail_test = F){
  test_output <- logical(ncol(input_binary_mat)) == F
  num_zero <- 1
  jud <- TRUE

  while (jud & num_zero < (nrow(input_binary_mat)-1)){
    pat <- c(1, rep(0, num_zero), 1)
    pat_1 <- rep(1, length(pat))
    foreach(i = 1:ncol(input_binary_mat)) %:% {
      temp <- input_binary_mat[, i]
      index_list <- search_pat(pat, temp)
      foreach(index = index_list) %dopar% {
        foreach(j = 1:ncol(input_binary_mat)) %do% {
          if(j != i){
            if(identical(input_binary_mat[, j][index:(index+length(pat)-1)],pat_1)){
              test_output[j] <- FALSE
              if(head_tail_test){
                foreach(k=(1:ncol(input_binary_mat))) %do% {
                  if(any(is.na(input_binary_mat[, k][index:(index+length(pat)-1)]))){test_output[k] = F}
                }
              }
            }
          }
        }
      }
    }
    num_zero <- num_zero+1
  }
  test_output
}
generate_exon_only_binary <- function(binary_folder_path, result_folder_name, custom_gene = NULL, file_pattern = "_binary.csv"){
  print(paste(Sys.time(), ": processing files.")); flush.console()
  init_files = list.files(binary_folder_path, pattern = file_pattern)
  dir.create(file.path(binary_folder_path, result_folder_name), showWarnings = FALSE)
  files <- init_files
  skipped_gene = NULL
  if(length(custom_gene) > 0){
    files = NULL
    foreach(i=1:length(files)) %dopar% {
      temp_gene_name <- strsplit(init_files[i], split = file_pattern)[[1]]
      if(temp_gene_name %in% custom_gene){processed_file <- c(processed_file, init_files[i])}
    }
  }

  foreach(i = 1:length(files)) %dopar% {
  #for(i in 53:53){
    gene = strsplit(files[i], split = file_pattern)[[1]]
    print(paste("Processing Files for Gene #", i,":", gene))
    file <- files[i]
    df <- read.csv(file.path(binary_folder_path, file))
    isoform <- names(df)[which(df[1, ] == "R")]
    remain_binary <- df[, which(df[1, ] == "R")]

    if(length(dim(remain_binary)) != 0){
      if(dim(remain_binary)[2] != 0){
        remain_binary <- remain_binary[-1, ]
        remain_binary_mat <- as.matrix(remain_binary)
        remain_binary_mat <- apply(remain_binary_mat, 2, as.numeric)
        result <- isoform_group(remain_binary_mat, isoform)
        file_name <- paste(gene, "_exon_only_binary.csv", sep = "")
        d <- data.frame(result$exon_binary)
        group_index <- result$index
        names(d) <- isoform
        d <- rbind(group_index, d)
        write.csv(d, file.path(binary_folder_path, result_folder_name, file_name))
      }
      else{
        skipped_gene <- c(skipped_gene, gene)
      }
    }

    else{
      skipped_gene <- c(skipped_gene, gene)
    }
  }

  write(skipped_gene, file.path(binary_folder_path, result_folder_name, "Skipped_Genes.txt"))
  print(paste(Sys.time(), ": done.")); flush.console()
  print(paste("Skipped Genes:", skipped_gene))

}

group_gff <- function(gff, custom_gene = NULL){
#  source("https://bioconductor.org/biocLite.R")
#  require("stringr"); require("rtracklayer"); require("GenomicRanges"); require("qpcR"); require("AnnotationHub"); require("openxlsx")
  annot_gr = import.gff(gff)
  gene_names = unique(annot_gr[!is.na(annot_gr$gffcompare_gene_name)]$gffcompare_gene_name)
  count_skipped = 0
  output_matrix <- list()
  analyzed_gene <- NULL
  if (length(custom_gene) == 0){
      count = 0
      foreach(i = 1:(length(gene_names))) %dopar% {
      gene = gene_names[i]
      count <- count + 1
      print(paste("Analyzing Gene # ", count,": ", gene, sep = "" ))
      gr_temp_inter = annot_gr[annot_gr$gffcompare_gene_name %in% gene]
      temp_gene_id = unique(gr_temp_inter$gffcompare_gene_id)[1]
      gr_temp = annot_gr[annot_gr$gffcompare_gene_id == temp_gene_id]
      temp_gr_exon = gr_temp[gr_temp$type == "exon"]
      transcript_count = length(unique(gr_temp$transcript_id))
      repetition_test = TRUE
      if(file.exists(file.path(getwd(),output_folder,paste(gene, "_simplified.gtf", sep = "")))){
        repetition_test = FALSE
        print(paste("Match Requirements: Skipped Completed Computation for Gene - ", gene, sep = "" ))
      }
      if (repetition_test & transcript_count > 3 & length(temp_gr_exon) > 2*transcript_count & nchar(gene) < 13){
        print(paste("Match Requirements: Computing for Gene - ", gene, sep = "" ))
        analyzed_gene <- c(analyzed_gene, gene)
        temp_data <- intron_filter_by_GRanges(gr_temp)
        if(length(temp_data$remain) <= 1){
          print("Skipped Isoform Group as the number of remaining isoforms is <= 1.")
          temp_reduce <- reduce(temp_data$all_exon[temp_data$all_exon$transcript_id == temp_data$remain])
          mcols(temp_reduce) <- cbind(paste(transcript_id, "group", 1, sep = "."), "exon", transcript_id, gene)
          names(elementMetadata(temp_reduce)) <- c("transcript_id", "type", "gffcompare_gene_id", "gffcompare_gene_name")
          temp_simplified_gr <- temp_reduce
          output_file_name <- paste(gene, "_simplified.gtf", sep = "")
          export(temp_simplified_gr, file.path(getwd(), output_folder, output_file_name))
          output_gr_list <- c(output_gr_list, temp_simplified_gr)
          } else {
          isoform_group_result <- isoform_group(temp_data$remain_binary, temp_data$remain)
          temp_simplified_gr <- simplify_GRanges(temp_data$all_exon, isoform_group_result$group, gene)
        }
        }
      }
    } else {
    count = 0
    foreach(i=1:(length(custom_gene))) %dopar% {
      gene = custom_gene[i]
      print(paste("Analyzing Gene # ", i,": ", gene, sep = "" ))
      gr_temp_inter = annot_gr[annot_gr$gffcompare_gene_name %in% gene]
      temp_gene_id = unique(gr_temp_inter$gffcompare_gene_id)[1]
      if(!is.na(temp_gene_id)){
        gr_temp = annot_gr[annot_gr$gffcompare_gene_id == temp_gene_id]
        transcript_count = length(unique(gr_temp$transcript_id))
        if (transcript_count > 3){
          analyzed_gene <- c(analyzed_gene, gene)
          print(paste("Match Requirements: Computing for Gene - ", gene, sep = "" ))
          
	#temp_data <- intron_filter_by_GRanges(gr_temp)
	temp_data <- get(paste0("temp_result_",gene))
          if(length(temp_data$remain) > 1){
            isoform_group_result <- isoform_group(temp_data$remain_binary, temp_data$remain)
            output_matrix <- qpcR:::rbind.na(output_matrix, c(gene, t(isoform_group_result$group)))
          }else{
            output_matrix <- qpcR:::rbind.na(output_matrix, c(gene, temp_data$remain))
          }
        }
      }
    }
  }
  output_matrix <- output_matrix[-1,]
  rownames(output_matrix) <- output_matrix[,1]
  output_matrix <- output_matrix[,-1]
  colnames(output_matrix) <- paste("Group_", 1:ncol(output_matrix), sep="")
  write.xlsx(data.frame(output_matrix), "result.csv", row.names = TRUE)
  list("output" = output_matrix, "computed_genes" = analyzed_gene)
}

heatmap_coverage <- function(gr_input, output_file_path = getwd(), isoform_group_index){
  ##Prepare gtf file from hg38
  if(!file.exists("Gencode_proteinCoding.gtf")){
    library(AnnotationHub)
    ah = AnnotationHub()
    ah = subset(ah,species=="Homo sapiens")
    qhs <- query(ah, c("Ensembl", "gene", "annotation", "grch38"))
    gtf <- qhs[["AH51014"]] #Homo_sapiens.GRCh38.85.gtf
    export(gtf[gtf$transcript_biotype %in% "protein_coding", ], "Gencode_proteinCoding.gtf")
  }

  gencode_gtf <- import.gff("Gencode_proteinCoding.gtf")
  gencode_gtf_gene <- gencode_gtf[gencode_gtf$gene_name %in% gr_input$gffcompare_gene_name]
  if(unique(gr_input$gffcompare_gene_name[!is.na(gr_input$gffcompare_gene_name)]) > 1){
    gencode_gtf_gene <- gencode_gtf[gencode_gtf$gene_name %in% unique(gr_input$gffcompare_gene_name[!is.na(gr_input$gffcompare_gene_name)])[1]]
  }
  gencode_exon = gencode_gtf_gene[gencode_gtf_gene$type == "exon"]
  gencode_list_input = split(gencode_exon, gencode_exon$transcript_id)
  gencode_list_ranges = reduce(gencode_list_input)

  gene_name <- unique(gr_input$gffcompare_gene_name)[1]
  dir.create(file.path(output_file_path, "output_heatmaps"), showWarnings = F)
  gr_exon <- gr_input[gr_input$type == "exon"]
  gr_list_input = split(gr_exon, gr_exon$transcript_id)
  isoform_names <- names(gr_list_input)
  isoform_count = length(unique(gr_exon$transcript_id))

  chr = unique(seqnames(gr_list_input[[1]]))
  chr_num = as.numeric(strsplit(toString(chr), split = "chr")[[1]][2])
  if(is.na(chr_num)){chr_num = 23}
  cov = base::as.vector(GenomicRanges::coverage(gr_exon)[[1]])
  extract = which(cov>0)
  extract_high = which(cov>isoform_count*0.1)

  gr_extract = reduce(GRanges(chr, range=IRanges(start=extract, end=extract)))
  gr_extract_high = reduce(GRanges(seqnames=chr, range=IRanges(extract_high, extract_high)))

  gr_base = reduce(gr_exon)
  stick = NULL
  foreach(i = 1:length(gr_base)) %dopar% {stick=c(stick, seq(start(gr_base[i]), end(gr_base[i])))}
  stick = GRanges(seqnames=chr, ranges =IRanges(start=stick, end=stick))
  mat = matrix(unlist(lapply(gr_list_input, populateMatrix, stick)), ncol=length(stick), byrow=T)
  mat = t(mat)
  coverage_out = matrix(NA, nrow = length(gr_base), ncol = length(gr_list_input))

  new_isoform <- NULL
  foreach(i = 1:length(unique(isoform_group_index))) %dopar% {
    temp_reduce_group <- reduce(unlist(gr_list_input[isoform_group_index == i]))
    seqlevels(temp_reduce_group) <- substring(seqlevels(temp_reduce_group), 4)
    temp_test <- TRUE
    j = 1
    while(j < (length(gencode_list_ranges)+1) & temp_test){
      temp_gencode <- unlist(gencode_list_ranges[j])
      temp_combined <- reduce(c(temp_reduce_group, temp_gencode))
      if(length(temp_combined) == length(temp_gencode)){
        #if(sum(start(temp_combined) == start(temp_gencode)) == length(temp_gencode) &
        #   sum(end(temp_combined) == end(temp_gencode)) == length(temp_gencode)){temp_test = FALSE}
        temp_test = FALSE
      }
      j <- j + 1
    }
    new_isoform <- c(new_isoform, temp_test)

  }

  foreach(i = 1:length(gr_base)) %:% {
    index_end <- sum(width(gr_base)[1:i])
    index_start <- index_end - width(gr_base)[i] + 1
    wid <- width(gr_base)[i]
    foreach(j = 1:length(gr_list_input)) %dopar% {
      val <- sum(mat[, j][index_start:index_end])/wid
      coverage_out[i, j] <- val
    }
  }
  coverage_out <- t(coverage_out)
  rownames(coverage_out) <- isoform_names
  colnames(coverage_out) <- width(gr_base)

  quantity <- NULL
  final_coverage_out <- matrix(NA, nrow = length(unique(iso_index)), ncol = length(gr_base))
  final_coverage_out_2 <- matrix(NA, nrow = length(unique(iso_index)), ncol = length(gr_base))
  foreach(i = 1:length(unique(iso_index))) %dopar% {
    if(length(which(iso_index == i)) == 1){
      final_coverage_out[i, ] <- coverage_out[which(iso_index == i), ]
      final_coverage_out_2[i, ] <- c(1)
      final_coverage_out_2[i, which(final_coverage_out[i, ] == 0)] <- NaN
      }
    else{
      foreach(j = 1:length(gr_base)) %dopar% {
        match_index <- which(iso_index == i)
        val2 <- sum(coverage_out[match_index, ][, j])/(max(coverage_out[match_index, ][, j])*length(coverage_out[match_index, ][, j]))
        val <- sum(coverage_out[match_index, ][, j])/length(match_index)
        final_coverage_out[i, j] <- val
        final_coverage_out_2[i, j] <- val2
        #if(is.na(val2)){final_coverage_out_2[i, j] <- 1}
      }
    }
    quantity <- c(quantity, length(which(iso_index == i)))
  }
  
  annot <- paste(" ", gene_name, "_", 1:nrow(final_coverage_out_2), sep = "")
  annot[new_isoform] <- paste(annot[new_isoform], "*")
  annot[!new_isoform] <- paste(annot[!new_isoform], " ")
  annot <- paste(annot, "  (m = ", quantity, ")", sep = "")
  final_coverage_out_2[which(is.na(final_coverage_out_2))] <- 0
  #rownames(final_coverage_out_2) <- paste("M = ", quantity, " ", gene_name, "_", 1:nrow(final_coverage_out_2), sep = "")
  rownames(final_coverage_out_2) <- annot
  colnames(final_coverage_out_2) <- width(gr_base)

  out_png <- file.path(output_file_path, "output_heatmaps", paste(gene, "_coverage_heatmap.png", sep = ""))
  myColor <- colorRampPalette(c("white", "blue", "black"))(100)
  mat <- matrix(NA, nrow = nrow(final_coverage_out), ncol = ncol(final_coverage_out))
  foreach(i = 1:nrow(mat)) %:% {
    foreach(j = 1:ncol(mat)) %dopar% {
      if(final_coverage_out[i, j] != 0){mat[i, j] <- round(final_coverage_out[i, j], 2)}
      else{mat[i, j] <- NA}
    }
  }
  cluster_r = FALSE
  if(length(unique(isoform_group_index)) > 1){cluster_r = TRUE}

  #pheatmap(final_coverage_out, cluster_rows = cluster_r, cluster_cols = FALSE, display_numbers = TRUE, color = myColor, main = paste("Isoform Groups of", gene), filename = out_png,
  #         cellwidth = 48.5, cellheight = 30, number_color = "white")
  pheatmap(final_coverage_out_2, cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = TRUE, color = myColor, main = paste("Isoform Groups of", gene), filename = out_png,
           cellwidth = 48.5, cellheight = 30, number_color = "white")
  
}
heatmap_exon_binary <- function(input_file){
  c <- read.csv(input_file)
  isoform <- names(c)
  c <- c[-1, ]
  c_mat <- apply(c, 2, as.numeric)
  c_mat <- t(c_mat[, -1])
  colnames(c_mat) <- paste("E", 1:ncol(c_mat), sep = "")
  pheatmap(c_mat, cluster_cols = FALSE)
}
#' Filter the transcripts with partial exons.
#'
#' @param gff gtf/gff file path
#' @param custom_gene a vector of targeted gene names
#' @return a list of accessible varaiables.
#' @examples
#' intron_filter_by_gff(your.gff.file)
#' intron_filter_by_gff("HumanGenome.gtf", custom_gene = c("CD5", "IL27RA"))
#' @export 
#' a csv file named "result.csv"
#' @export 
#' gtf file(s) for each gene after the process
#' 
#' @details Systemetically go through the process of intron_filter_by_GRanges of all selected genes in a gtf file. 
#' It produces three types of gtf files for each gene in the same directory. 
#' (_remain2.0.gtf, _potential_new_exon.gtf, _filtered2.0.gtf) 
#' It also creates a subdirectory called /binary_results. 
#' All the binary forms are stored as csv files in this subdirectory.
#' \describe{
#' \item{\strong{Output List}}{}
#'     \item{stat_prefiltered}{A vector of number of transcript(s) with partial exon(s) for each gene}
#'     \item{stat_remain}{A vector of number of remaining transcript(s)}
#'     \item{stat_filtered}{A vector of number of transcript(s) with intron retention(s) for each gene}
#'     \item{percent_prefiltered}{A vector of percentage of transcript(s) with partial exon(s) for each gene}
#'     \item{percent_remain}{A vector of percentage of remaining transcript(s)}
#'     \item{percent_filtered}{A vector of percentage of transcript(s) with intron retention(s) for each gene}
#'     \item{anlyzed_gene}{A vector of gene name(s) that have been analyzed}
#'     
#' }
intron_filter_by_gff <- function(gff, custom_gene = NULL){
  require("stringr"); require("rtracklayer"); require("GenomicRanges"); require("progress"); require("AnnotationHub")
  annot_gr = import.gff(gff)
  gene_names = unique(annot_gr[!is.na(annot_gr$gffcompare_gene_name)]$gffcompare_gene_name)
  stat_prefiltered = NULL
  stat_remain = NULL
  stat_filtered = NULL
  stat_isoform = NULL
  perc_prefiltered = NULL
  perc_remain = NULL
  perc_filtered = NULL
  analyzed_gene = NULL
  count_skipped = 0
  #exception <- c("PSMB9")
  #custom_gene <- custom_gene[!custom_gene %in% exception]
  dir.create(file.path(getwd(), "binary_results"), showWarnings = FALSE)

  if (length(custom_gene) == 0){
    foreach(i = 1:(length(gene_names))) %dopar% {
      #for(i in 1:10){
      gene = gene_names[i]
      gr_temp_inter = annot_gr[annot_gr$gffcompare_gene_name %in% gene]
      temp_gene_id = unique(gr_temp_inter$gene_id)[1]
      gr_temp = annot_gr[annot_gr$gene_id == temp_gene_id]
      temp_gr_exon = gr_temp[gr_temp$type == "exon"]
      transcript_count = length(unique(gr_temp$transcript_id))
      repetition_test = TRUE
      if(file.exists(paste("gencode_", gene, ".gtf", sep = ""))
         & file.exists(paste(getwd(), "/binary_results/", gene, "_binary.csv", sep = ""))
         & (file.exists(paste(gene,"_remain2.0.gtf", sep=""))
            | file.exists(paste(gene,"_potential_new_exon.gtf", sep=""))
            | file.exists(paste(gene,"_filtered2.0.gtf", sep="")))){
        repetition_test = FALSE
      }

      if (repetition_test & transcript_count > 5 & length(temp_gr_exon) > 2*transcript_count ){
        analyzed_gene <- c(analyzed_gene, gene)
        temp_data <- intron_filter_by_GRanges(gr_temp)
        stat_prefiltered <- c(stat_prefiltered, temp_data$prefiltered_new_exon_count)
        stat_remain <- c(stat_remain, temp_data$remain_count)
        stat_filtered <- c(stat_filtered, temp_data$filtered_count)
        isoform_num <- temp_data$isoform_count
        temp_data$remain
        stat_isoform <- c(stat_isoform, isoform_num)
        perc_prefiltered <- c(perc_prefiltered, temp_data$prefiltered_new_exon_count/isoform_num)
        perc_remain <- c(perc_remain, temp_data$remain_count/isoform_num)
        perc_filtered <- c(perc_filtered, temp_data$filtered_count/isoform_num)
        d <- data.frame(stat_prefiltered, perc_prefiltered, stat_filtered, perc_filtered, stat_remain, perc_remain, stat_isoform, analyzed_gene)
        names(d) <- c("Count_Prefiltered", "Percent_Prefiltered", "Count_Filtered", "Percent_Filtered", "Count_Remain", "Percent_Remain", "Count_Isoforms", "Gene_Name")
        write.csv(d, "filter_stats.csv")
        category = c(rep("P", temp_data$prefiltered_new_exon_count), rep("F", temp_data$filtered_count), rep("R", temp_data$remain_count))
        binary_data = data.frame(temp_data$prefiltered_binary, temp_data$filtered_binary, temp_data$remain_binary)
        #print(c(temp_data$prefiltered_new_exon, temp_data$filtered, temp_data$remain))
        names(binary_data) <- c(temp_data$prefiltered_new_exon, temp_data$filtered, temp_data$remain)
        binary_data <- rbind(category, binary_data)
        bin_file_name <- paste(gene, "_binary.csv", sep = "")
        write.csv(binary_data, paste(getwd(), "/binary_results/", bin_file_name, sep = ""))
      }


    }
  }

  else{
    count = 0
    foreach(i = 1:(length(gene_names))) %dopar% {
      #for(i in 1:12){
      gene = gene_names[i]
      if(gene %in% custom_gene){
        count <- count + 1
        print(paste("Analyzing Gene # ", count,": ", gene, sep = "" ))
        gr_temp_inter = annot_gr[annot_gr$gffcompare_gene_name %in% gene]
        temp_gene_id = unique(gr_temp_inter$gene_id)[1]
        gr_temp = annot_gr[annot_gr$gene_id == temp_gene_id]
        temp_gr_exon = gr_temp[gr_temp$type == "exon"]
        transcript_count = length(unique(gr_temp$transcript_id))
        repetition_test = TRUE
        if(file.exists(paste("gencode_", gene, ".gtf", sep = ""))
           & file.exists(paste(getwd(), "/binary_results/", gene, "_binary.csv", sep = ""))
           & (file.exists(paste(gene,"_remain2.0.gtf", sep=""))
              | file.exists(paste(gene,"_potential_new_exon.gtf", sep=""))
              | file.exists(paste(gene,"_filtered2.0.gtf", sep="")))){
          repetition_test = FALSE
        }

        if (repetition_test & transcript_count > 1 & length(temp_gr_exon) > 2*transcript_count ){
          analyzed_gene <- c(analyzed_gene, gene)
          temp_data <- intron_filter_by_GRanges(gr_temp)
          stat_prefiltered <- c(stat_prefiltered, temp_data$prefiltered_new_exon_count)
          stat_remain <- c(stat_remain, temp_data$remain_count)
          stat_filtered <- c(stat_filtered, temp_data$filtered_count)
          isoform_num <- temp_data$isoform_count
          stat_isoform <- c(stat_isoform, isoform_num)
          perc_prefiltered <- c(perc_prefiltered, temp_data$prefiltered_new_exon_count/isoform_num)
          perc_remain <- c(perc_remain, temp_data$remain_count/isoform_num)
          perc_filtered <- c(perc_filtered, temp_data$filtered_count/isoform_num)
          d <- data.frame(stat_prefiltered, perc_prefiltered, stat_filtered, perc_filtered, stat_remain, perc_remain, stat_isoform, analyzed_gene)
          names(d) <- c("Count_Prefiltered", "Percent_Prefiltered", "Count_Filtered", "Percent_Filtered", "Count_Remain", "Percent_Remain", "Count_Isoforms", "Gene_Name")
          #print(d)
          write.csv(d, "filter_stats.csv")
          category = c(rep("P", temp_data$prefiltered_new_exon_count), rep("F", temp_data$filtered_count), rep("R", temp_data$remain_count))
          binary_data = data.frame(temp_data$prefiltered_binary, temp_data$filtered_binary, temp_data$remain_binary)
          #print(c(temp_data$prefiltered_new_exon, temp_data$filtered, temp_data$remain))
          names(binary_data) <- c(temp_data$prefiltered_new_exon, temp_data$filtered, temp_data$remain)
          binary_data <- rbind(category, binary_data)
          bin_file_name <- paste(gene, "_binary.csv", sep = "")
          write.csv(binary_data, paste(getwd(), "/binary_results/", bin_file_name, sep = ""))
        }

      }
    }
  }

  #print(paste("Total Skipped Genes: ", count_skipped))

  list("stat_prefiltered" = stat_prefiltered, "stat_remain" = stat_remain, "stat_filtered" = stat_filtered,
       "analyzed_gene" = analyzed_gene, "percent_filtered" = perc_filtered,
       "percent_prefiltered" = perc_prefiltered, "percent_remain" = perc_remain)
}
#' Filter the transcripts with partial exons in comparison to the Gencode Reference Genome. Then filter the ones with intron retention through a self-learning process within the remaining group. Output the binary forms of all three groups, their names, count, and all exons as in GRange.
#'
#' @param gr_input GenomicRanges::GRanges Object
#' @return a list. 
#' @examples
#' intron_filter_by_GRanges(import.gff(file_name))
#' @export
#' @details This function returns a list containing the following elements.    
#' * element 1
#' * element 2 is thatl...
#' binary forms, count and names of prefiltered, filtered and remain isoforms; 
#' {"prefiltered_new_exon", "remain", "filtered", "prefiltered_binary", "remain_binary", "filtered_binary", "prefiltered_new_exon_count", "remain_count", "filtered_count","isoform_count", "all_exon"}
#' @import AnnotationHub
#' @importFrom rtracklayer export

intron_filter_by_GRanges <- function(gr_input){
  ##Prepare gtf file from hg38
  if(!file.exists("Gencode_proteinCoding.gtf")){
    ah = AnnotationHub::AnnotationHub()
    ah = AnnotationHub::subset(ah,species=="Homo sapiens")
    qhs <- AnnotationHub::query(ah, c("Ensembl", "gene", "annotation", "grch38"))
    gtf <- qhs[["AH51014"]] #Homo_sapiens.GRCh38.85.gtf
    rtracklayer::export(gtf[gtf$transcript_biotype %in% "protein_coding", ], "Gencode_proteinCoding.gtf")
  }


  gene_name = unique(gr_input[!is.na(gr_input$gffcompare_gene_name)]$gffcompare_gene_name)
  export_gtf_name <- paste("gencode_", gene_name, ".gtf", sep = "")
  if(!file.exists(export_gtf_name)){
    gencode_gtf <- import.gff("Gencode_proteinCoding.gtf")
    gr = gr_input
    gencode_gtf_gene <- gencode_gtf[gencode_gtf$gene_name %in% gr$gffcompare_gene_name]
    if(unique(gr$gffcompare_gene_name[!is.na(gr$gffcompare_gene_name)]) > 1){
      gencode_gtf_gene <- gencode_gtf[gencode_gtf$gene_name %in% unique(gr$gffcompare_gene_name[!is.na(gr$gffcompare_gene_name)])[1]]
    }
    gencode_exon = gencode_gtf_gene[gencode_gtf_gene$type == "exon"]
    export(gencode_exon, export_gtf_name)
  }

  gencode_exon = import.gff(export_gtf_name)

  gencode_list_input = split(gencode_exon, gencode_exon$transcript_id)
  gencode_base = reduce(gencode_exon)
  gr = gr_input
  gr_exon = gr[gr$type=="exon"]
  gr_list_input = split(gr_exon, gr_exon$transcript_id)
  isoform_names <- names(gr_list_input)
  isoform_count = length(unique(gr_exon$transcript_id))

  chr = unique(seqnames(gr_list_input[[1]]))
  chr_num = as.numeric(strsplit(toString(chr), split = "chr")[[1]][2])
  if(is.na(chr_num)){chr_num = 23}
  if(length(GenomicRanges::coverage(gr_exon)) > 1){
    cov = base::as.vector(GenomicRanges::coverage(gr_exon)[[chr_num]])
  }
  if(length(GenomicRanges::coverage(gr_exon)) == 1){
    cov = base::as.vector(GenomicRanges::coverage(gr_exon)[[1]])
  }
  
  extract = which(cov>0)
  extract_high = which(cov>isoform_count*0.1)

  gr_extract = reduce(GRanges(chr, range=IRanges(start=extract, end=extract)))
  gr_extract_high = reduce(GRanges(seqnames=chr, range=IRanges(extract_high, extract_high)))

  gr_base = reduce(gr_exon)
  gr_subject = gr_extract_high
  gr_intron = GenomicRanges::setdiff(GRanges(seqname=chr, ranges = IRanges(start(gr_base), end(gr_base))), gr_subject)
  gr_exons = gr_subject; gr_exons$type = "exon"
  gr_introns = gr_intron;
  if (length(gr_introns) > 0) { gr_introns$type="intron" }
  gr_tract = c(gr_exons, gr_introns)
  gr_tract = gr_tract[order(start(gr_tract))]


  new1 <- GRanges(seqnames = chr, range = IRanges(start(gencode_exon), end(gencode_exon)))
  new2 <- GRanges(seqnames = chr, range = IRanges(start(gr_exon), end(gr_exon)))
  recombined_gencode_gr <- c(new1, new2)
  recombined_gencode_gr <- recombined_gencode_gr[order(start(recombined_gencode_gr))]
  recombined_gr_base <- reduce(recombined_gencode_gr)

  reduce_intron_start = NULL
  reduce_intron_end = NULL
  if (length(recombined_gr_base) > 1){
    foreach (i = 1:(length(recombined_gr_base)-1)) %dopar% {
      reduce_intron_start <- c(reduce_intron_start, end(recombined_gr_base)[i]+1)
      reduce_intron_end <- c(reduce_intron_end, start(recombined_gr_base)[i+1]-1)
    }
    reduce_intron_start <- reduce_intron_start[!is.na(reduce_intron_start)]
    reduce_intron_end <- reduce_intron_end[!is.na(reduce_intron_end)]
    reduce_intron = GRanges(seqname=chr, ranges = IRanges(reduce_intron_start, reduce_intron_end))
    gr_combined = c(recombined_gr_base, reduce_intron)
    gr_combined = gr_combined[order(start(gr_combined))]
    piece = c(disjoin(recombined_gencode_gr), reduce_intron)
    }

  if(length(recombined_gr_base) == 1){
    gr_combined = recombined_gr_base
    piece = disjoin(recombined_gencode_gr)
    }

  
  piece = piece[order(start(piece))]
  piece_width = width(piece)
  piece_start = start(piece)
  piece_end = end(piece)



  stick = NULL
  foreach (i = 1:length(gr_combined)) %dopar% {stick=c(stick, seq(start(gr_combined[i]), end(gr_combined[i])))} #changed
  stick = GRanges(seqnames=chr, ranges =IRanges(start=stick, end=stick))
  mat = matrix(unlist(lapply(gr_list_input, populateMatrix, stick)), ncol=length(stick), byrow=T)
  mat = t(mat)

  mat_gencode = matrix(FALSE, nrow = nrow(mat), ncol = length(gencode_list_input))
  foreach (j = (1:length(gencode_list_input))) %:% {
    temp <- unlist(gencode_list_input[j])
    foreach(i = (1:length(temp))) %dopar% {
      gencode_start <- start(temp)[i] - start(gr_base)[1] + 1
      gencode_end <- gencode_start + width(temp)[i] - 1
      if(gencode_start < 0 & gencode_end > 0){gencode_start = 1}
      if (gencode_end > 0){mat_gencode[, j][gencode_start:gencode_end] = T}
    }
  }

  #for (i in (1:length(gencode_base))){
  #  gencode_start <- start(gencode_base)[i] - start(gr_base)[1]
  #  gencode_end <- gencode_start + width(gencode_base)[i]
  #  print(c(gencode_start, gencode_end))
  #  mat_gencode[gencode_start:gencode_end] = T
  #}
  out_disjoin <- matrix(NA, ncol = ncol(mat), nrow = length(piece))
  gencode_out_disjoin <- matrix(NA, ncol= ncol(mat_gencode), length(piece))


  # Based on Disjoin
  foreach(i = 1:(ncol(out_disjoin))) %:% {
    a <- mat[, i]
    foreach(j = 1:(nrow(out_disjoin))) %dopar% {
      tract_e <- sum(width(piece)[1:j])
      tract_s <- tract_e - width(piece)[j] + 1
      compare_region <- a[tract_s:tract_e]
      val <- sum(compare_region, na.rm = TRUE)/(width(piece)[j])
      out_disjoin[j, i] <- val
    }
  }
  # Separation Line


  # Calculating for gencode standard set
  foreach(i = (1:ncol(gencode_out_disjoin))) %:% {
    foreach(j = 1:(nrow(gencode_out_disjoin))) %dopar% {
      tract_e <- sum(width(piece)[1:j])
      tract_s <- tract_e - width(piece)[j] + 1
      compare_region <- mat_gencode[, i][tract_s:tract_e]
      val <- sum(compare_region, na.rm = TRUE)/(width(piece)[j])
      if(val != 0 & val != 1){val <- round(val)}
      gencode_out_disjoin[j, i] <- val
    }
  }
  # Recognize head and tail
  foreach(i = (1:ncol(gencode_out_disjoin))) %dopar% {
    temp <- gencode_out_disjoin[, i]
    first_one <- which(temp==1)[1]
    last_one <- which(temp==1)[sum(temp)]
    if(!is.na(first_one) & first_one > 1){gencode_out_disjoin[, i][1:(first_one-1)] = NA}
    if((length(last_one) != 0)){
      if(last_one < nrow(gencode_out_disjoin)){gencode_out_disjoin[, i][(last_one+1):nrow(gencode_out_disjoin)] = NA}
      }
    }
  # Separation Line


  # Getting rid of isoforms with potential new exons
  pre_test <- logical(ncol(out_disjoin)) == F
  foreach (i = (1:ncol(out_disjoin))) %dopar% {
    gencode_test <- find_retention(cbind(gencode_out_disjoin, out_disjoin[,i]), head_tail_test = T)
    if(identical(gencode_test[1:(length(gencode_test)-1)],rep(FALSE, length(gencode_test)-1))){pre_test[i] = F} # & gencode_test[length(gencode_test)] == T
  }
  # Separation Line

  isoform_names_filtered = isoform_names[pre_test] # Potential isoforms with new exons are removed
  out_disjoin_filtered = out_disjoin[, pre_test]


  if(!is.matrix(out_disjoin_filtered)){
    out_disjoin_filtered <- as.matrix(out_disjoin_filtered)
  }

  if(ncol(out_disjoin_filtered) > 0){
    test <- find_retention(out_disjoin_filtered)
    foreach(i = 1:length(isoform_names_filtered)) %dopar% {
      if (gr[gr$transcript_id == isoform_names_filtered[i]]$matchAnnot_gene[1] != gr[gr$transcript_id == isoform_names_filtered[i]]$gffcompare_gene_name[1]){
        test[i] <- FALSE
      }

    }
  }

  if(ncol(out_disjoin_filtered) == 0){
    test = NULL
  }




  if(sum(!pre_test)>0){export(gr[gr$transcript_id %in% isoform_names[!pre_test]], paste(gene_name,"_potential_new_exon.gtf", sep=""))}
  if(sum(test)>0){export(gr[gr$transcript_id %in% isoform_names_filtered[test]], paste(gene_name,"_remain2.0.gtf", sep=""))}
  if((length(test) - sum(test))>0){export(gr[gr$transcript_id %in% isoform_names_filtered[test==F]], paste(gene_name,"_filtered2.0.gtf", sep=""))}
  
  filtered_binary = NULL
  if((length(test) - sum(test))>0){filtered_binary = out_disjoin_filtered[, !test]}
  list("prefiltered_new_exon" = isoform_names[!pre_test], "remain" = isoform_names_filtered[test], "filtered" = isoform_names_filtered[test==F],
       "prefiltered_binary" = out_disjoin[, !pre_test], "remain_binary" = out_disjoin_filtered[, test], "filtered_binary" = filtered_binary,
       "prefiltered_new_exon_count" = sum(!pre_test), "remain_count" = sum(test), "filtered_count" = (length(test) - sum(test)),
       "isoform_count" = isoform_count, "all_exon" = gr_exon)
}
isoform_group <- function(input_binary_mat, isoform_names){
  require(mgcv)
  if (length(isoform_names) != ncol(input_binary_mat)){print("Warning: Isoform Names Input might not match the binary matrix.")}
  exon_start_index = NULL
  exon_end_index = NULL
  index = 1
  standard_exon_binary_list = rep(0, nrow(input_binary_mat))
  while (index < nrow(input_binary_mat)){
    if (sum(input_binary_mat[index, ]) > 0){
      temp = 0
      while((index+temp) < nrow(input_binary_mat) & sum(input_binary_mat[(index+temp), ]) > 0){
        temp <- temp + 1
      }
      exon_start_index <- c(exon_start_index, index)
      exon_end_index <- c(exon_end_index, index+temp-1)
      standard_exon_binary_list[index:(index+temp-1)] = 1
      index <- index + temp
    }

    else{index <- index + 1}

  }

  if(sum(input_binary_mat[nrow(input_binary_mat), ]) > 0 & sum(input_binary_mat[nrow(input_binary_mat)-1, ]) == 0){
    exon_start_index <- c(exon_start_index, nrow(input_binary_mat))
    exon_end_index <- c(exon_end_index, nrow(input_binary_mat))
  }

  exon_binary_mat <- matrix(0, nrow = length(exon_start_index), ncol = ncol(input_binary_mat))
  foreach(j = (1:ncol(input_binary_mat))) %:% {
    foreach(i = (1:nrow(exon_binary_mat))) %dopar% {
      if (sum(input_binary_mat[, j][exon_start_index[i]:exon_end_index[i]]) > 0){
        exon_binary_mat[i, j] <- 1
      }
    }
  }

  category <- uniquecombs(t(exon_binary_mat))
  group_index <- attr(category, "index")
  cluster = group_index
  indx <- order(cluster)
  pre_group <- list()

  foreach(i = 1:length(unique(group_index))) %dopar% {
    pos = which(group_index==i)
    pre_group[[i]] <- isoform_names[pos]
  }

  list("index"= group_index, "exon_binary" = exon_binary_mat, "group" = pre_group)

}
populateMatrix = function(query, subject){
  info = findOverlaps(query, subject)
  ind = 1:length(subject)
  return((ind %in% subjectHits(info)))
}
quantify <- function(group_csv, quant_file){
  require("openxlsx")
  quant_data <- read.table(quant_file)
  group_info <- read.xlsx(group_csv)
  output_data <- list()
  foreach(i = 1:length(group_info[, 1])) %:% {
    gene = group_info[i, 1]
    print(paste("Analyzing Gene #", i, "-", gene))
    temp_no_na <- group_info[i, ][!is.na(group_info[i,])][-1]
    foreach(j = 1:length(temp_no_na)) %dopar% {
      temp_col_name <- paste(gene, j, sep = "_")
      temp_isoform_names <- unlist(strsplit(temp_no_na[j], split = ", "))
      temp_sum <- rep(0, ncol(quant_data)-1)
      foreach(k = 1:length(temp_isoform_names)) %do% {
        temp_transcript_id = temp_isoform_names[k]
        temp_index <- which(quant_data[,1] == temp_transcript_id)
        temp_row_data<- quant_data[temp_index,][-1]
        temp_processed_row_data <- as.numeric(as.vector(as.matrix(temp_row_data)))
        temp_sum <- rowSums(cbind(temp_sum, temp_processed_row_data), na.rm = TRUE)
      }
      output_data <- rbind(output_data, c(gene, j, temp_sum))

    }
  }
  colnames(output_data) <- c("gene_name", "group", as.vector(as.matrix(quant_data[1,]))[-1])
  output_data[which(output_data == "0")] <- NA #could be removed if 0 wants to be kept
  write.csv(output_data, "quantity_NA.csv", row.names = FALSE)
}
search_pat <- function(pattern, x){
  result <- NULL
  foreach(i = 1:(length(x)-length(pattern))) %dopar% {
    compare <- x[i:(i+length(pattern)-1)]
    if (identical(pattern,compare)){
      result <- c(result,i)
    }
  }

  result
}
simplify_gff <- function(gff, custom_gene = NULL, output_folder = "simplified_gtf"){
  #source("https://bioconductor.org/biocLite.R")
  #require("stringr"); require("rtracklayer"); require("GenomicRanges"); require("progress"); require("AnnotationHub")
  annot_gr = import.gff(gff)
  gene_names = unique(annot_gr[!is.na(annot_gr$gffcompare_gene_name)]$gffcompare_gene_name)
  count_skipped = 0
  dir.create(file.path(getwd(), "processed_gtf_files"), showWarnings = FALSE)
  output_gr_list <- NULL
  analyzed_gene <- NULL

  if (length(custom_gene) == 0){
      count = 0
      foreach(i = 1:(length(gene_names))) %dopar% {
      gene = gene_names[i]
      count <- count + 1
      print(paste("Analyzing Gene # ", count,": ", gene, sep = "" ))
      gr_temp_inter = annot_gr[annot_gr$gffcompare_gene_name %in% gene]
      temp_gene_id = unique(gr_temp_inter$gene_id)[1]
      gr_temp = annot_gr[annot_gr$gene_id == temp_gene_id]
      temp_gr_exon = gr_temp[gr_temp$type == "exon"]
      transcript_count = length(unique(gr_temp$transcript_id))
      repetition_test = TRUE
      if(file.exists(file.path(getwd(), output_folder, paste(gene, "_simplified.gtf", sep = "")))){
        repetition_test = FALSE
        print(paste("Match Requirements: Skipped Completed Computation for Gene - ", gene, sep = "" ))
        output_gr_list <- c(output_gr_list, import.gff(file.path(getwd(), output_folder, paste(gene, "_simplified.gtf", sep = ""))))
      }

      if (repetition_test & transcript_count > 3){
        print(paste("Match Requirements: Computing for Gene - ", gene, sep = "" ))
        analyzed_gene <- c(analyzed_gene, gene)
        temp_data <- intron_filter_by_GRanges(gr_temp)
        if(length(temp_data$remain) <= 1){
          print("Skipped Isoform Group as the number of remaining isoforms is <= 1.")
          temp_reduce <- reduce(temp_data$all_exon[temp_data$all_exon$transcript_id == temp_data$remain])
          mcols(temp_reduce) <- cbind(paste(transcript_id, "group", 1, sep = "."), "exon", transcript_id, gene)
          names(elementMetadata(temp_reduce)) <- c("transcript_id", "type", "gffcompare_gene_id", "gffcompare_gene_name")
          temp_simplified_gr <- temp_reduce
          output_file_name <- paste(gene, "_simplified.gtf", sep = "")
          export(temp_simplified_gr, file.path(getwd(), output_folder, output_file_name))
          output_gr_list <- c(output_gr_list, temp_simplified_gr)
          } else {
          isoform_group_result <- isoform_group(temp_data$remain_binary, temp_data$remain)
          temp_simplified_gr <- simplify_GRanges(temp_data$all_exon, isoform_group_result$group, gene)
          output_gr_list <- c(output_gr_list, temp_simplified_gr)
        }
        }
      }
    } else {
    count = 0
    foreach(i = 1:(length(gene_names))) %dopar% {
      gene = gene_names[i]
      analyzed_gene <- c(analyzed_gene, gene)
      if(gene %in% custom_gene){
        count <- count + 1
        print(paste("Analyzing Gene # ", i=count,": ", gene, sep = "" ))
        gr_temp_inter = annot_gr[annot_gr$gffcompare_gene_name %in% gene]
        temp_gene_id = unique(gr_temp_inter$gene_id)[1]
        gr_temp = annot_gr[annot_gr$gene_id == temp_gene_id]
        temp_gr_exon = gr_temp[gr_temp$type == "exon"]
        transcript_count = length(unique(gr_temp$transcript_id))
        repetition_test = TRUE
        if(file.exists(file.path(getwd(), output_folder, paste(gene, "_simplified.gtf", sep = "")))){
          repetition_test = FALSE
          print(paste("Match Requirements: Skipped Completed Computation for Gene - ", gene, sep = "" ))
          output_gr_list <- c(output_gr_list, import.gff(file.path(getwd(), output_folder, paste(gene, "_simplified.gtf", sep = ""))))
        }

        if (repetition_test & transcript_count > 3){
          print(paste("Match Requirements: Computing for Gene - ", gene, sep = "" ))
          temp_data <- intron_filter_by_GRanges(gr_temp)
          if(length(temp_data$remain) <= 1){
            print("Skipped Isoform Group as the number of remaining isoforms is <= 1.")
            temp_reduce <- reduce(temp_data$all_exon[temp_data$all_exon$transcript_id == temp_data$remain])
            mcols(temp_reduce) <- cbind(paste(transcript_id, "group", 1, sep = "."), "exon", transcript_id, gene)
            names(elementMetadata(temp_reduce)) <- c("transcript_id", "type", "gffcompare_gene_id", "gffcompare_gene_name")
            temp_simplified_gr <- temp_reduce
            output_file_name <- paste(gene, "_simplified.gtf", sep = "")
            export(temp_simplified_gr, file.path(getwd(), output_folder, output_file_name))
            output_gr_list <- c(output_gr_list, import.gff(output_file_name))
            } else {
            isoform_group_result <- isoform_group(temp_data$remain_binary, temp_data$remain)
            temp_simplified_gr <- simplify_GRanges(temp_data$all_exon, isoform_group_result$group, gene)
            output_gr_list <- c(output_gr_list, temp_simplified_gr)
            }
          }
        }
      }
    }
  output_gr_list <- unlist(GRangesList(output_gr_list))
  if(length(output_gr_list) > 0){export(output_gr_list, paste("simplified_", gff, sep = ""))}
  list("output" = output_gr_list, "computed_genes" = analyzed_gene)



}
simplify_GRanges <- function(input_gr_exon, group_list, input_gene_name, output_folder = "simplified_gtf"){
  addition_gtf <- GRanges(NULL)
  transcript_id <- paste(unlist(strsplit(input_gr_exon$transcript_id[1], split = "[.]"))[1], unlist(strsplit(input_gr_exon$transcript_id[1], split = "[.]"))[2], sep = ".")
  foreach(i = 1:length(group_list)) %dopar% {
    subgroup_isoforms <- group_list[[i]]
    subgroup_exon <- input_gr_exon[input_gr_exon$transcript_id %in% subgroup_isoforms]
    temp_reduce <- reduce(subgroup_exon)
    mcols(temp_reduce) <- cbind(paste(transcript_id, "group", i, sep = "."), "exon", transcript_id, input_gene_name)
    names(elementMetadata(temp_reduce)) <- c("transcript_id", "type", "gffcompare_gene_id", "gffcompare_gene_name")
    addition_gtf <- c(addition_gtf, temp_reduce)
  }

  dir.create(file.path(getwd(), output_folder), showWarnings = FALSE)
  output_file_name <- paste(input_gene_name, "_simplified.gtf", sep = "")
  export(addition_gtf, file.path(getwd(), output_folder, output_file_name))
  addition_gtf
}
