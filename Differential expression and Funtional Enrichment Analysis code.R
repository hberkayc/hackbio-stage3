if (!requireNamespace("BiocManager", quietly = TRUE))
  
install.packages("BiocManager")
# Load BiocManager after installation
library(BiocManager)


#Install the Following Tools:
  
BiocManager::install("TCGAbiolinks")

#edgeR and Limma for Differential Analysis

BiocManager::install("edgeR")

BiocManager::install("edgeR")

##EDASeq

BiocManager::install("EDASeq")

#gplots

install.packages('gplots')

#sesameData

BiocManager::install("sesameData")

installed.packages()["sesameData", ]
library(sesameData)  # if using the sesameData package
library(SummarizedExperiment)

install.packages("httr2", type = "binary")
library(TCGAbiolinks) 

library("sesameData")

#Summarized Experiment

BiocManager::install("SummarizedExperiment")

# to get project summary for TCGA- Lung Adenocarcinoma
query <- GDCquery(project = "TCGA-LUAD")


#Project information
getProjectSummary("TCGA-LUAD")
?GDCquery

# Download and process data
gbmQ <- GDCquery(
project = "TCGA-LUAD",  # Corrected project name
data.category = "Transcriptome Profiling",  # Correct spelling for "Profiling"
data.type = "Gene Expression Quantification"
        )
install.packages("pryr")
library(pryr)


GDCdownload(gbmQ)
gbm.data <- GDCprepare(gbmQ)
str(gbm.data)  
head(gbm.data)

#Explore metadata information
gbm.data$race
gbm.data$tumor_descriptor
gbm.data$barcode
gbm.data$tissue_type


# Create metadata for our use
simpleMeta <- data.frame(
  barcode = gbm.data$barcode,
  race = gbm.data$race,
  tumor_type = gbm.data$tissue_type,
  stringsAsFactors = FALSE  # Prevents automatic conversion of strings to factors
)

# View the first few rows of the new data frame
head(simpleMeta)

#Preprocessing and Normalization

#Select unstranded dataset
gbm.raw.data <- assay(gbm.data)
# Check dimensions directly if it's a data frame or matrix
dim(gbm.raw.data)

#Downsize data to 20 Normal and 20 Tumor 
SelectBarcodes <- c(subset(simpleMeta,tumor_type=="Normal")$barcode[c(1:20)])
SelectBarcodes1 <- c(subset(simpleMeta,tumor_type=="Tumor")$barcode[c(1:20)])

# Combining the two selected barcode vectors into one
CombinedBarcodes <- c(SelectBarcodes, SelectBarcodes1)

SelectData <- gbm.raw.data[,c(CombinedBarcodes)]
dim(SelectData)

#Data Normalization and filtering
normData <- TCGAanalyze_Normalization(tabDF = SelectData,geneInfo = geneInfoHT,method="geneLength")
#then filter
filtData <- TCGAanalyze_Filtering(tabDF = normData,method = "quantile", qnt.cut = 0.25)
dim(filtData)

#Export Data for ML
ml_data <- data.frame(t(filtData))
write.csv(ml_data, file = "ml_data.csv")


#check number of barcodes
length(CombinedBarcodes)


# Differential gene analysis
SelectResults <- TCGAanalyze_DEA(
  mat1 = filtData[, SelectBarcodes[1:10]],  
  mat2 = filtData[, SelectBarcodes1[11:20]], 
  Cond1type = "Tumor",
  Cond2type = "Normal",
  pipeline = "edgeR",
  fdr.cut = 0.01,  # Set false discovery rate cut-off
  logFC.cut = 2    # Set log fold change cut-off
)

# Differential expression analysis with tissue_type
SelectResults.level <- TCGAanalyze_LevelTab(
  SelectResults,
  "Tumor",
  "Normal",
  filtData[, SelectBarcodes[1:10]],  # First 10 samples
  filtData[, SelectBarcodes[11:20]]   # Next 10 samples
)



#Visualise
head(SelectResults.level)
dim(SelectResults.level)
heat.data <- filtData[rownames(SelectResults.level),]

# Color the plot by kind of tumor
cancer.type <- c(rep("Tumor", 10), rep("Normal", 10))
ccodes <- c()

# Assign colors based on cancer type
for (i in cancer.type) {
  if (i == "Tumor") {
    ccodes <- c(ccodes, "red")
  } else {
    ccodes <- c(ccodes, "blue")
  }
}

# Create the heatmap
library(gplots)  # Ensure gplots package is loaded
X11()
# Create the heatmap
heatmap.2(x = as.matrix(heat.data),         
  col = hcl.colors(10, palette = "Blue-Red"),  
  Rowv = FALSE,                     
  Colv = TRUE,                      
  scale = "row",                    
  sepcolor = "black",             
  trace = "none",                 
  key = TRUE,                 
  dendrogram = "column",            
  cexRow = 0.5,                   
  cexCol = 1,                   
  main = "Tumor vs Normal",         
  na.color = "black",         
)


#Functional enrichment analysis

#View volcano
plot(x=SelectResults.level$logFC, y=-log10(SelectResults.level$FDR))

upregulated_genes <- rownames(subset(SelectResults.level,logFC>2))
downregulated_genes <- rownames(subset(SelectResults.level, logFC < -2))

BiocManager::install("biomaRt", force=TRUE)
library(biomaRT)

#Covert ensemble IDs to gene IDs using biomart
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")


upregulated_genes <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'),
                           filters = 'ensembl_gene_id',
                           values = upregulated_genes, 
                           mart = mart)$hgnc_symbol


downregulated_genes <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'),
                             filters = 'ensembl_gene_id',
                             values = downregulated_genes,
                             mart = mart)$hgnc_symbol


#Perform enrichment analysis
up.EA <- TCGAanalyze_EAcomplete(TFname = "upregulated",upregulated_genes)
up.EA <- TCGAanalyze_EAcomplete(TFname = "downregulated",downregulated_genes)


# Visualize enrichment analysis results
TCGAvisualize_EAbarplot(
  tf = rownames(up.EA$ResBP),
  GOBPTab = up.EA$ResBP,
  GOCCTab = up.EA$ResCC,
  GOMFTab = up.EA$ResMF,
  PathTab = up.EA$ResPat,
  nRGTab = upregulated_genes,
  nBar = 5,
  text.size = 2,
  fig.width = 30,
  fig.height = 15
)

TCGAvisualize_EAbarplot(
  tf = rownames(up.EA$ResBP),
  GOBPTab = up.EA$ResBP,
  GOCCTab = up.EA$ResCC,
  GOMFTab = up.EA$ResMF,
  PathTab = up.EA$ResPat,
  nRGTab = downregulated_genes,
  nBar = 5,
  text.size = 2,
  fig.width = 30,
  fig.height = 15
)

#Save data
# Define the path to the Downloads directory
downloads_path <- file.path(Sys.getenv("USERPROFILE"), "Downloads")

# List of data frames to save
data_frames <- list(
  filtData = filtData,
  normData = normData,
  available_marts = available_marts,
  heat.data = heat.data,
  gbm.raw.data = gbm.raw.data,
  SelectData = SelectData,
  SelectResults = SelectResults,
  SelectResults.level = SelectResults.level,
  simpleMeta = simpleMeta
)

# Loop through the list and save each data frame as a CSV file
for (name in names(data_frames)) {
 # Construct the file path
  file_path <- file.path(downloads_path, paste0(name, ".csv"))
  
# Save the data frame
  write.csv(data_frames[[name]], file = file_path, row.names = FALSE)
  
# Print a message confirming the save
  cat("Data has been saved to:", file_path, "\n")
}








