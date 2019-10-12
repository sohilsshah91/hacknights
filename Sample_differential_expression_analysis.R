
# PURPOSE -----------------------------------------------------------------

#' To conduct a differential expression analysis that will be a sample 
#' for conducting the analysis using TCGA data
#' The TCGA RNASeq data is illumina hiseq Level 3 RSEM normalized expression data

print(paste0("start: ",Sys.time()))
# install packages/load libraries -----------------------------------------

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("DESeq2",update=F)
library(DESeq2)
library(tidyverse)


# set variables -----------------------------------------------------------

data_dir <- "data/"
response_name <- "patient.race"
rnaseq_file <- "lihc_rnaseq.csv.gz"
clinical_file <- "lihc_clinical.csv.gz"

# load data ---------------------------------------------------------------

rnaseq_data <- read_csv(paste0(data_dir,rnaseq_file))
clinical_data <- read_csv(paste0(data_dir,clinical_file))

# clean rnaseq data for input ----------------------------------------------------

ont_map <- read_csv(paste0(data_dir,"gene_ontology_mappings.csv")) %>% 
  distinct(entrezgene_id,hgnc_symbol)

entrez_ids <- lapply(
  strsplit(
    colnames(rnaseq_data),
    "\\|"),
  function(x)x[2]
  ) %>% unlist()

gene_names <- lapply(
  strsplit(
    colnames(rnaseq_data),
    "\\|"),
  function(x)x[1]
) %>% unlist()
rnaseq_data_clean <- rnaseq_data
colnames(rnaseq_data_clean) <- gene_names
rnaseq_data_clean <- rnaseq_data_clean[,
                                       sapply(
                                         colnames(rnaseq_data_clean),
                                         function(x){grepl("[A-z]",x)}
                                         )
                                       ]
rnaseq_data_clean <- data.frame(rnaseq_data[,entrez_ids %in% ont_map$entrezgene_id])
colnames(rnaseq_data_clean) <- ont_map$entrezgene_id[entrez_ids %in% ont_map$entrezgene_id]
rnaseq_data_clean <- rnaseq_data_clean[,na.omit(colnames(rnaseq_data_clean))]
row.names(rnaseq_data_clean) <- rnaseq_data$bcr_patient_barcode
rnaseq_data_clean <- rnaseq_data_clean[order(row.names(rnaseq_data_clean)),]

# clean clinical data for input -------------------------------------------

response_name <- "patient.gender"

response <- clinical_data %>% 
  select(!!response_name) %>% 
  arrange() %>% 
  fill(!!response_name) %>% 
  data.frame()

row.names(response) <- clinical_data$patient.bcr_patient_barcode

# DE analysis -------------------------------------------------------------

rnaseq_rownames <- 
  unlist(
    lapply(
      strsplit(
        str_to_lower(row.names(rnaseq_data_clean)),"-"),
      function(x)paste0(x[1:3],collapse="-")
    )
  )
rnaseq_row_map <- data.frame(
  "old" = row.names(rnaseq_data_clean),
  "new" = rnaseq_rownames
)

clinical_rownames <- row.names(response)

inter_rownames <- intersect(rnaseq_rownames,clinical_rownames)
rows.to.keep <- rnaseq_row_map[rnaseq_row_map$new %in% inter_rownames,][
    !duplicated(rnaseq_row_map$new),"old"]

rsem.in <- data.frame(t(rnaseq_data_clean[rows.to.keep,]))
inter_response <- sapply(
  response[inter_rownames,],
  function(x)str_replace_all(x," ","_")
  ) %>% unname %>% factor

#converting rsem to whole integers per
#https://www.biostars.org/p/320594/
dds <- DESeqDataSetFromMatrix(ceiling(rsem.in),DataFrame(inter_response),~inter_response)

res <- DESeq(dds)
results(res,)
res.vst <- vst(res)
plotPCA(res.vst)

print(paste0("end: ",Sys.time()))


