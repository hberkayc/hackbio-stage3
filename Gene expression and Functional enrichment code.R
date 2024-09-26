if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
library(BiocManager)


#Install Tools:
BiocManager::install("TCGAbiolinks", force = TRUE)
#edgeR and lima for Differential Analysis
BiocManager::install("edgeR")
#EDASeq
BiocManager::install("EDASeq")
#gplots
install.packages('gplots')
#sesameData
BiocManager::install("sesameData", force = TRUE)
installed.packages()["sesameData", ]
library(sesameData)  # if using the sesameData package
library(SummarizedExperiment)

install.packages("httr2", type = "binary")
library(TCGAbiolinks) 

library("sesameData")

#Summarized Experiment

BiocManager::install("SummarizedExperiment", force = TRUE)

# Example to get project summary for TCGA- Lung Adenocarcinoma
query <- GDCquery(project = "TCGA-LUAD")

#Project information
getProjectSummary("TCGA-LUAD")
?GDCquery

#Download and process data
gbmQ <- GDCquery(
project = "TCGA-LUAD",  
data.category = "Transcriptome Profiling", 
data.type = "Gene Expression Quantification"
        )
install.packages("pryr")
library(pryr)

GDCdownload(gbmQ)
gbm.data <- GDCprepare(gbmQ)
str(gbm.data)  
head(gbm.data)
colnames(gbm.data)
row.names(gbm.data)

#Explore metadata information
gbm.data$race
gbm.data$tumor_descriptor
gbm.data$barcode
gbm.data$tissue_type



# Create metadata for our use
simpleMeta <- data.frame(
  barcode = gbm.data$barcode,
  race = gbm.data$race,
  tumor_type = gbm.data$tumor_descriptor,
  tissue_type = gbm.data$tissue_type,
  stringsAsFactors = FALSE  # Prevents automatic conversion of strings to factors
)

# View the first few rows of the new data frame
head(simpleMeta)



#Preprocessing 
#Select unstranded dataset
gbm.raw.data <- assay(gbm.data)

#Transpose the data so that genes are columns and samples are rows
raw.data <- data.frame(t(gbm.raw.data))

raw.data <- data.frame(samples = rownames(raw.data), raw.data)


colnames(simpleMeta)
# Rename the 'barcode' column to 'samples'
colnames(simpleMeta)[colnames(simpleMeta) == "barcode"] <- "samples"
head(simpleMeta)

#Merge metadata with dataset (raw.data) for ML
combined_data <- merge(raw.data, simpleMeta, by = "samples")
View(combined_data)

# Check dimensions directly if it's a data frame or matrix
dim(combined_data)

#Downsize data to 20 tumor and 20 Normal 
#checking the structure of my metadata to ensure the tissue type is clearly labeled
table(simpleMeta$tissue_type)


install.packages("dplyr")
library(dplyr)
# Set seed for reproducibility
set.seed(123)  
tumor_samples <- simpleMeta %>% 
  filter(tissue_type == "Tumor") %>%  # Filter for tumor samples
  sample_n(20)  # Randomly select 20 samples

normal_samples <- simpleMeta %>% 
  filter(tissue_type == "Normal") %>%  # Filter for normal samples
  sample_n(20)  # Randomly select 20 samples

# Combine selected samples
selected_samples <- bind_rows(tumor_samples, normal_samples)


#sample_ids <- selected_samples$samples
#head(selected_samples$samples)

#Restructuring the gene_data for accurate subseting
head(rownames(combined_data))
#gene_data_selected <- combined_data[sample_ids, ]

# Check the first few rows of the gene data
head(combined_data)

# Check the column names to identify sample identifiers
colnames(combined_data)

# Convert the sample identifier column to row names
rownames(combined_data) <- combined_data$samples

# Remove the sample identifier column from the data
combined_data <- combined_data[, -which(names(combined_data) == "samples")]

# Check to ensure row names are now set correctly
head(rownames(combined_data))

# Subset the gene data using the sample IDs from selected_samples
final_data <- combined_data[selected_samples$samples, ]

# Check the dimensions of the new dataset
dim(final_data)

# View the first few rows of the selected gene data
head(final_data)

#save file
write.csv(final_data, "final_data.csv", row.names = TRUE)

#final data with a labeled rowname.(run this code, if you need this file)
#final_data <- data.frame(samples = rownames(final_data), final_data)



