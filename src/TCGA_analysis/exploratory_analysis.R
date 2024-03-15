# Load libraries for differential expression analysis
library(DESeq2)
library(dplyr)

# Load sample data
sample_data <- read.csv("TCGA_data/gdc_sample_sheet.2024-02-27.tsv", sep = "\t")

# Preping RNAseq data

# Filter out the RNAseq data
rnaseq_samples <- sample_data %>%
  filter(Data.Category == "Transcriptome Profiling" & Data.Type == "Gene Expression Quantification") %>%
  # Remove Metastaic samples
    filter(!grepl("Metastatic", Sample.Type))

# Load the RNAseq data using File.ID dir and File.Name as file
rnaseq_files <- file.path("TCGA_data", rnaseq_samples$File.ID, rnaseq_samples$File.Name)

# Add sample name as S1, S2, S3, etc to rnaseq_samples
rnaseq_samples$sample <- paste0("S", 1:nrow(rnaseq_samples))

# Loop through all the files and read the unstranded column into rnaseq_counts data
for (i in 1:length(rnaseq_files)) {
# for (i in 1:10) {
    rnaseq_data <- read.table(rnaseq_files[i], header = TRUE, sep = "\t")
    if (i == 1) {
        rnaseq_counts <- rnaseq_data[, "unstranded"]
    } else {
        rnaseq_counts <- cbind(rnaseq_counts, rnaseq_data[, "unstranded"])
    }

    # print every 100 iterations
    if (i %% 100 == 0) {
        print(i)
    }
    
    # Check dimensions of rnaseq_data
    # print(dim(rnaseq_data))
}


# Add sample names to rnaseq_counts
colnames(rnaseq_counts) <- rnaseq_samples$Sample.ID

# Add rownames to rnaseq_counts from gene_id as rnaseq_data
rownames(rnaseq_counts) <- rnaseq_data$gene_id

# Remove first four rows
# rnaseq_counts <- rnaseq_counts[-c(1:4), ]

# Filter low count genes
keep <- rowSums(rnaseq_counts) >= 100
rnaseq_counts_filtered <- rnaseq_counts[keep,]



# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = rnaseq_counts,
                              colData = rnaseq_samples,
                              design = ~ Sample.Type)

# Filter low count genes
keep <- rowSums(counts(dds)) >= 100
dds <- dds[keep,]

# Split the rownames with 
duplicated(str_split(rownames(dds), pattern =  "\\.", simplify = T)[,1]) %>% table

rownames(dds) <- str_split(rownames(dds), pattern =  "\\.", simplify = T)[,1]


# Run DESeq2
dds <- DESeq(dds, betaPrior = TRUE)

# Save the  dds object
saveRDS(dds, "TCGA_data/mRNA_dds.rds")


# Differential expression analysis
res <- results(dds, contrast = c("Sample.Type", "Primary Tumor", "Solid Tissue Normal"))

# Save the results
write.csv(as.data.frame(res), "TCGA_data/mRNA_deseq2_results.csv")

# read in microRNA data
# Filter out the RNAseq data
rnaseq_samples <- sample_data %>%
  filter(Data.Category == "Transcriptome Profiling" & Data.Type == "miRNA Expression Quantification") %>%
  # Remove Metastaic samples
    filter(!grepl("Metastatic", Sample.Type))

# Load the RNAseq data using File.ID dir and File.Name as file
rnaseq_files <- file.path("TCGA_data", rnaseq_samples$File.ID, rnaseq_samples$File.Name)

# Add sample name as S1, S2, S3, etc to rnaseq_samples
rnaseq_samples$sample <- paste0("S", 1:nrow(rnaseq_samples))

# Loop through all the files and read the unstranded column into rnaseq_counts data
for (i in 1:length(rnaseq_files)) {
# for (i in 1:10) {
    rnaseq_data <- read.table(rnaseq_files[i], header = TRUE, sep = "\t")
    if (i == 1) {
        rnaseq_counts <- rnaseq_data[, "read_count"]
    } else {
        rnaseq_counts <- cbind(rnaseq_counts, rnaseq_data[, "read_count"])
    }

    # print every 100 iterations
    if (i %% 100 == 0) {
        print(i)
    }
    
    # Check dimensions of rnaseq_data
    # print(dim(rnaseq_data))
}


# Add sample names to rnaseq_counts
colnames(rnaseq_counts) <- rnaseq_samples$Sample.ID

# Add rownames to rnaseq_counts from gene_id as rnaseq_data
rownames(rnaseq_counts) <- rnaseq_data$miRNA_ID

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = rnaseq_counts,
                              colData = rnaseq_samples,
                              design = ~ Sample.Type)

# Filter low count genes
keep <- rowSums(counts(dds)) >= 10                        

# Run DESeq2
dds <- DESeq(dds, betaPrior = TRUE)

# Save the  dds object
saveRDS(dds, "TCGA_data/miR_dds.rds")


# Differential expression analysis
res <- results(dds, contrast = c("Sample.Type", "Primary Tumor", "Solid Tissue Normal"))

# Save the results
write.csv(as.data.frame(res), "TCGA_data/miR_deseq2_results.csv")



# read in methylation data
# Filter out the RNAseq data
rnaseq_samples <- sample_data %>%
  filter(Data.Category == "DNA Methylation") %>%
  # Remove Metastaic samples
    filter(!grepl("Metastatic", Sample.Type))

# Load the RNAseq data using File.ID dir and File.Name as file
rnaseq_files <- file.path("TCGA_data", rnaseq_samples$File.ID, rnaseq_samples$File.Name)

# Add sample name as S1, S2, S3, etc to rnaseq_samples
rnaseq_samples$sample <- paste0("S", 1:nrow(rnaseq_samples))

# Loop through all the files and read the unstranded column into rnaseq_counts data
for (i in 1:length(rnaseq_files)) {
# for (i in 1:10) {
    rnaseq_data <- read.table(rnaseq_files[i], header = FALSE, sep = "\t")
    if (i == 1) {
        rnaseq_counts <- as.data.frame(rnaseq_data[, "V2"])
    } else {
        rnaseq_counts <- cbind(rnaseq_counts, rnaseq_data[, "V2"])
    }

    # print every 100 iterations
    if (i %% 100 == 0) {
        print(i)
    }
    
    # Check dimensions of rnaseq_data is same as rna_seq_counts
    if(dim(rnaseq_data)[1] != dim(rnaseq_counts)[1])
    {   
        print(dim(rnaseq_counts))
        # print message and index
        print(i)
        print("Dimensions do not match")
    }
    
}


# Add sample names to rnaseq_counts
colnames(rnaseq_counts) <- rnaseq_samples$sample

# Add rownames to rnaseq_counts from gene_id as rnaseq_data
rownames(rnaseq_counts) <- rnaseq_data$V1

# Save as tsv
write.table(rnaseq_counts, "TCGA_data/methylation_data.tsv", sep = "\t")
