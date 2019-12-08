
# PURPOSE -----------------------------------------------------------------

#' To conduct a differential expression analysis that will be a sample 
#' for conducting the analysis using TCGA data
#' The TCGA RNASeq data is illumina hiseq Level 3 RSEM normalized expression data
start_time <- Sys.time()
cat(paste0("start: ",start_time,"\n"))

#' EXAMPLE RUN
#' -----------
#' Rscript TCGA_differential_expression_analysis.R "" lihc_rnaseq.csv.gz lihc_clinical.csv.gz 100 patient.race asian white
#' Rscript TCGA_differential_expression_analysis.R "" lihc_rnaseq.csv.gz lihc_clinical.csv.gz 100 patient.gender female male
#' Rscript TCGA_differential_expression_analysis.R "" lihc_rnaseq.csv.gz lihc_clinical.csv.gz full patient.ethnicity "hispanic or latino" "not hispanic or latino"
#' Rscript TCGA_differential_expression_analysis.R "" lihc_rnaseq.csv.gz lihc_clinical.csv.gz full patient.vital_status dead alive
#' 
#' 
# install packages/load libraries -----------------------------------------

if (!requireNamespace("BiocManager", quietly = TRUE)){
  suppressWarnings(suppressMessages(install.packages("BiocManager")))}

if(require(DESeq2)){
  suppressWarnings(suppressMessages(library(DESeq2)))
}else{
  suppressWarnings(suppressMessages(BiocManager::install("DESeq2",update=F)))
  suppressWarnings(suppressMessages(library(DESeq2)))
}
if(require(tidyverse)){
  suppressWarnings(suppressMessages(library(tidyverse)))
}else{
  install.packages("tidyverse")  
  suppressWarnings(suppressMessages(library(DESeq2)))
}
if(require(BiocParallel)){
  suppressWarnings(suppressMessages(library(BiocParallel)))
}else{
  BiocManager::install("BiocParallel",update=F)
  suppressWarnings(suppressMessages(library(BiocParallel)))
}
register(MulticoreParam(4))

# set variables -----------------------------------------------------------

set.seed(42)

args = commandArgs(trailingOnly=TRUE)

data_dir <- args[1]
#data_dir <- ""
rnaseq_file <- args[2]
#rnaseq_file <- "lihc_rnaseq.csv.gz"
clinical_file <- args[3]
#clinical_file <- "lihc_clinical.csv.gz"
gene_size <- args[4]
#gene_size <- 1e2

cancer_type <- strsplit(rnaseq_file,"_")[[1]][1]

outfile <- paste0(cancer_type,"_DESeq2_")

response_name <- args[5]
#response_name <- "patient.race"
t <- args[6]
#t <- "white"
c <- args[7]
#c <- "asian"

# load data ---------------------------------------------------------------

rnaseq_data <- read_csv(paste0(data_dir,rnaseq_file))
clinical_data <- read_csv(paste0(data_dir,clinical_file)) %>% 
  data.frame()

# clean rnaseq data for input ----------------------------------------------------

genes <- lapply(
  strsplit(
    colnames(rnaseq_data),
    "\\|"),
  function(x)x[1]
) %>% unlist()
genes_logical <- !is.na(sapply(genes,function(x){nchar(x)>2}))

rnaseq_data_clean <- rnaseq_data[,genes_logical]
colnames(rnaseq_data_clean) <- genes[genes_logical]
rnaseq_data_clean <- data.frame(rnaseq_data_clean)
row.names(rnaseq_data_clean) <- rnaseq_data$bcr_patient_barcode
rnaseq_data_clean <- rnaseq_data_clean[order(row.names(rnaseq_data_clean)),]

# DE analysis -------------------------------------------------------------

rnaseq_data_clean$bcr_patient_barcode <- 
  unlist(
    lapply(
      strsplit(
        str_to_lower(rnaseq_data$bcr_patient_barcode),"-"),
      function(x)paste0(x[1:3],collapse="-")
    )
  )

colData <- 
  inner_join(
    rnaseq_data_clean %>% select(bcr_patient_barcode),
    clinical_data,
    by=c("bcr_patient_barcode"="patient.bcr_patient_barcode")
    )

eligible_responses <- 
  colnames(colData)[
    apply(colData,2,function(x){sum(!is.na(x))>(nrow(colData)*.9)})
    ]

colData_filled <- 
  apply(colData %>% select(eligible_responses),
        2,
        function(x){
          tab <- table(unname(unlist(x)))
          mode <- names(tab)[order(tab)][1]
          x[is.na(x)] <- mode
          factor(x)
          }
        )

rows.to.keep <- 
  rnaseq_data_clean$bcr_patient_barcode %in% colData$bcr_patient_barcode

rsem.in <- data.frame(
  t(
    rnaseq_data_clean[rows.to.keep,colnames(rnaseq_data_clean)!="bcr_patient_barcode"]
    )
  )

if(is.numeric(as.integer(gene_size))){
  rsem.in.sub <- rsem.in[sample(1:nrow(rsem.in),gene_size),]
}
if(gene_size=="full"){
  rsem.in.sub <- rsem.in
}

#converting rsem to whole integers per
#https://www.biostars.org/p/320594/

formula <- as.formula(paste0("~",response_name))

suppressMessages(dds <- DESeqDataSetFromMatrix(ceiling(rsem.in.sub),
                              colData_filled,
                              formula))


res <- DESeq(dds)
results <- results(res,contrast=c(response_name,t,c))
results_df <- data.frame(results)
rownames(results_df) <- 1:nrow(results_df)
results_df$Gene <- row.names(dds)


outfile2 <- paste0(outfile,gene_size,
                   "_sampled_genes_",response_name,"_",
                   t,"_vs_",c,".csv")

results_df %>%  
  write_csv(paste0(data_dir,outfile2))

end_time <- Sys.time()
cat(paste0("end: ",end_time,"\n"))
cat(
  paste0("duration: ",
         round(as.numeric(end_time-start_time,units="mins"),2),
         " minutes\n")
  )


