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

sample_names <- scan("samples.txt", character())
sample_gtfs <- paste(sample_names, ".gtf", sep = "")
TPM_cutoff <- c(0.1, 1, 10)

registerDoParallel(detectCores() - 2)
output <- foreach(gtf = sample_gtfs, .combine = rbind) %do% {
  gtf_data <- import.gff(gtf)
  print(gtf)
  temp <- foreach(TPM_threshold = TPM_cutoff, .combine = c) %dopar%{
    #gtf_data <- import.gff(gtf)
    filtered_ID <- gtf_data$transcript_id[which(gtf_data$TPM >= TPM_threshold)]
    c(length(filtered_ID),
      #length(intersect(filtered_ID,filtered_genes)),
      length(intersect(filtered_ID,all_PBID)),
      #length(intersect(filtered_ID,filtered_genes))/length(filtered_genes))
      length(intersect(filtered_ID,all_PBID))/length(all_PBID))
  }
  print(temp)
  temp
}

stopImplicitCluster()

row.names(output) <- sample_names
colnames(output) <- as.vector(rbind(paste(TPM_cutoff, "-TPM", sep = ""), "Overlap", "%Overlap"))
output <- cbind(output, length(all_PBID))
colnames(output)[ncol(output)] <- "LRs"



#b$transcript_id[which(b$TPM > 0.1)]