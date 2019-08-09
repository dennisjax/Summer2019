filtered_gene_names <- unique(gtf$gffcompare_gene_name[which(gtf$transcript_id %in% filtered_genes)], na.rm = T)


for (gene_bubble in gsub("_simplified.gtf","",dir("simplified_gtf/"))){
  bubblePlots("quantity_NA.csv",plot_gene_name=gene_bubble)
}







library(foreach)
library(doParallel)
library(dplyr)
file <- read.csv("all_samples.chained_count.csv", row.names = 1)
PV_files <- file %>% dplyr::select(contains("PV"))
num_PV <- ncol(PV_files)
non_PV_files <- file %>% dplyr::select(-contains("PV"))
num_non_PV <- ncol(non_PV_files)
filtered_genes <- c()
stats <- c()

# Select the PBIDs whose sums of transcripts are non-zero
candidate_PBID <- row.names(PV_files)[which(rowSums(PV_files, na.rm = T) != 0)]
all_PBID <- row.names(PV_files)

#for (ID in candidate_PBID){
#  non_PV_score <- sum(non_PV_files[ID, ], na.rm = T)/num_non_PV
#  PV_score <- sum(PV_files[ID, ], na.rm = T)/num_PV
#  if ((PV_score/non_PV_score) >= 2){
#    filtered_genes <- c(filtered_genes, ID)
#  }
#}

registerDoParallel(detectCores() - 2)
filtered_genes <- foreach(ID = candidate_PBID, .combine = c) %dopar% {
  non_PV_score <- sum(non_PV_files[ID, ], na.rm = T)/num_non_PV
  PV_score <- sum(PV_files[ID, ], na.rm = T)/num_PV
  if ((PV_score/non_PV_score) >= 2){
    ID
  }
}
stopImplicitCluster()

respondent <- c(import.gff("PV018.gtf"), import.gff("PV032.gtf"), import.gff("PV034.gtf"), import.gff("PV042.gtf"))
nonrespondent <- c(import.gff("PV006.gtf"), import.gff("PV013.gtf"), import.gff("PV011.gtf"), import.gff("PV024.gtf"))

confirm_result <- as.data.frame(filtered_genes)

p_val <- c()
#start below
for(candidate in filtered_genes[7717:9167]){
  res_group <- as.numeric(subset(respondent, transcript_id == candidate)$TPM)
  res_group[which(is.na(res_group))] = 0
  non_group <- as.numeric(subset(nonrespondent, transcript_id == candidate)$TPM)
  non_group[which(is.na(non_group))] = 0
  print(candidate)
  if(length(non_group) == 0 | length(res_group) == 0){
    p_val <- c(p_val, NA)
  }
  else{
    p_val <- c(p_val, t.test(res_group, non_group, mu = 0, 
              alternative = "greater")$p.value)
  }
}

confirm_result$pvalue <- p_val



PV032_files <- file %>% dplyr::select(contains("PV032"))
PV034_files <- file %>% dplyr::select(contains("PV034"))
non_PV_files <- file %>% dplyr::select(-contains("PV"))
df <- data.frame(rowSums(PV032_files, na.rm = T), rowSums(PV034_files, na.rm = T), rowSums(non_PV_files, na.rm = T))
colnames(df) <- c("PV032", "PV034", "Other Healthy")
write.table(df, file = "all_samples.chained_count_rearranged.txt")

quantify("./result.csv", "./all_samples.chained_count.txt")


PV032_files <- res %>% dplyr::select(contains("PV032"))
PV034_files <- res %>% dplyr::select(contains("PV034"))
non_PV_files <- res %>% dplyr::select(-contains("PV"))
non_PV_files <- non_PV_files[, -1]
non_PV_files <- non_PV_files[, -1]
df <- data.frame(rowSums(PV032_files, na.rm = T)/ncol(PV032_files), 
                 rowSums(PV034_files, na.rm = T)/ncol(PV034_files), 
                 rowSums(non_PV_files, na.rm = T)/ncol(non_PV_files))
colnames(df) <- c("PV032", "PV034", "Other Healthy")

res <- cbind(res[, 1:2], df)
      
      
bubblePlots <- function(csv, plot_gene_name){
  require(dplyr)
  require(reshape2)
  require(ggplot2)
  dir.create("bubbleplot", showWarnings = F)
  res <- read.csv(csv, header = T, stringsAsFactors = F)
  
  ########################################### NEW
  PV032_files <- res %>% dplyr::select(contains("PV032"))
  PV034_files <- res %>% dplyr::select(contains("PV034"))
  non_PV_files <- res %>% dplyr::select(-contains("PV"))
  non_PV_files <- non_PV_files[, -1]
  non_PV_files <- non_PV_files[, -1]
  df <- data.frame(rowSums(PV032_files, na.rm = T)/ncol(PV032_files), 
                   rowSums(PV034_files, na.rm = T)/ncol(PV034_files), 
                   rowSums(non_PV_files, na.rm = T)/ncol(non_PV_files))
  colnames(df) <- c("PV032", "PV034", "Other Healthy")
  
  res <- cbind(res[, 1:2], df*100)
  res[res == 0] <- NA
  ##################################ABOVE NEW############################
  
  gene <- res %>% dplyr::filter(gene_name %in% c(plot_gene_name))
  
  gene.melt <- melt(gene, id=c("gene_name", "group"))
  
  
  gene.melt$label <- paste0(gene.melt$gene_name, "_", sprintf("%2d", gene.melt$group))
  ########################################################################  
  #targetSamples <- grep("MEL", unique(gene.melt$variable), value = T)
  targetSamples <- unique(gene.melt$variable)
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
  graph_width <- 0.3*ncol(gene)
  if(graph_height < 7){graph_height <- 6.5}
  if(graph_width < 7){graph_width <- 6.5}
  ggsave(file.path(getwd(), "bubbleplot", paste(plot_gene_name, "_bubbleplot.png", sep = "")), height = graph_height, width = graph_width, dpi = 400)
}



t <- read.csv("LRfiltered_pathway.csv")
pathway_name <- as.vector(t$Term)
for(i in 1:length(pathway_name)){
  temp <- unlist(strsplit(pathway_name[i], split=':', fixed=TRUE))[2]
  print(temp)
  pathway_name[i] <- temp
}


barplot(t$PValue, names.arg = pathway_name, cex.names = 0.5, horiz = T, las = 1)




sample_names <- scan("samples.txt", character())
sample_gtfs <- paste(sample_names, ".gtf", sep = "")
TPM_threshold = 0.1


gene_list <- c() 
for (gtf in sample_gtfs){
  gtf_data <- import.gff(gtf)
  print(gtf)
  
    #gtf_data <- import.gff(gtf)
  filtered_ID <- gtf_data$transcript_id[which(gtf_data$TPM >= TPM_threshold)]
  print(filtered_ID)
  gene_list <- c(gene_list, filtered_ID)
}


RNA_seq_filtered_gene_names <- unique(gtf$gffcompare_gene_name[which(gtf$transcript_id %in% gene_list)], na.rm = T)
write(RNA_seq_filtered_gene_names, "RNA_seq_gene_names.txt")


