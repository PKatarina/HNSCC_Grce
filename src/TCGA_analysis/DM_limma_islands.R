# Load the required packages
library(minfi)
library(dplyr)
library(limma)
#library(IlluminaHumanMethylationEPICanno.ilm10b2.hg38)
library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
data(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# Load sample data
sample_data <- read.csv("TCGA_data/gdc_sample_sheet.2024-02-27.tsv", sep = "\t")

#dnam_data_test <- read.csv("TCGA_data/methylation_data_test.tsv", sep = "\t", header = TRUE)

dnam_data <- read.csv("TCGA_data/methylation_data.tsv", sep = "\t", header = TRUE)


methylation_samples <- sample_data %>%
  filter(Data.Category == "DNA Methylation") %>%
    # Remove Metastaic samples
    filter(!grepl("Metastatic", Sample.Type))

methylation_samples_test <- methylation_samples
rownames(methylation_samples_test) <- colnames(dnam_data)
#colData(grs) <- DataFrame(methylation_samples_test)

# Get the annotation data
anno_data <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# Extract the CpG island information
cpg_islands <- anno_data[rownames(dnam_data),]$Islands_Name

# Remove probes that do not map to a CpG island
dnam_data <- dnam_data[cpg_islands != "",]
# remove the elements cpg_islands with empty values
cpg_islands <- cpg_islands[cpg_islands != ""]


# Add the 'cpg_islands' vector as a column to the data frame
dnam_data$cpg_islands <- cpg_islands

# Group by 'cpg_islands' and calculate the mean
mean_data <- aggregate(. ~ cpg_islands, dnam_data, mean)

# Store cpg_islands as rownames
rownames(mean_data) <- mean_data$cpg_islands

# Remove the 'cpg_islands' column
mean_data <- mean_data[, -1]

# Print the first few rows of the mean data
head(mean_data)

# # Createa density plot
# png("test_cgis.png")

# densityPlot(mean_data, sampGroups=methylation_samples_test$Sample.ID,main="Raw", legend=FALSE)
# # Save plot to file
# dev.off()

individual <- factor(methylation_samples_test$Sample.Type)
individual <- relevel(individual, ref = "Solid Tissue Normal")
design <- model.matrix(~individual)

# Fit the linear model
#fit <- lmFit(getBeta(grs), design)
fit <- lmFit(mean_data, design)
fit2 <- eBayes(fit)

top.table <- topTable(fit2, coef=2, number=Inf)
write.csv(top.table, file = "TCGA_data/dnam_limma_out_cgis_topTable.csv")

# write.fit(fit2, file = "TCGA_data/dnam_limma_out_cgis.txt")

# use rownames of mean_data to subset annot_data
anno_data_sub <- anno_data %>%
    as.data.frame() %>%
    dplyr::filter(anno_data$Islands_Name %in% rownames(mean_data))

# Save the annotation data
write.csv(anno_data_sub, "TCGA_data/dnam_limma_annot_cgis.csv")