### SET-UP Codes￣

# Clear workspace
rm(list = ls()); gc()

# get required library
library(Seurat); library(dplyr); library(Matrix)

#	Set output directory 
Merge.outdir = "/media/yschang/SMG1/ITG/Merge/"
# load(paste0(Merge.outdir, 'Baseline.RData'))

# Metadata
Metadata<-read.csv('/media/yschang/SMG1/ITG/Merge/merged_cases.csv', header = T, row.names = 1)

###############################################################################
### Part A. Merging Entire Seurat Object: merging entire datasets ############# ###############################################################################

## SMC dataset
GSE131907_raw.data = read.csv('/media/yschang/SMG1/SMC_Nat_Comm_2020/Data/GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz', header = TRUE, row.names = 1, sep = '\t')
GSE131907.sbj <- CreateSeuratObject(counts = GSE131907_raw.data, project = "GSE131907", min.cells = 3, min.features = 200); rm(GSE131907_raw.data)
# Add Metadata
GSE131907.sbj$dataset.id <- Idents(GSE131907.sbj)
GSE131907.sbj$sample.id <- substr(colnames(GSE131907.sbj), 18, 30)
GSE131907.sbj$orig.ident <- GSE131907.sbj$sample.id
GSE131907.sbj <- SetIdent(GSE131907.sbj, value = 'orig.ident')
# QC with reference (1) percent.mito:<= 20, (2) UMIs count:100 ~ 150,000, (3) gene count:200 ~ 10,000
# Get mito percentage
GSE131907.sbj[["percent.mt"]] <- PercentageFeatureSet(GSE131907.sbj, pattern = "^MT-")
# Get sqj
GSE131907.sqj <- subset(x = GSE131907.sbj, subset = nCount_RNA >= 100 & nCount_RNA <= 150000 & nFeature_RNA >= 200 & nFeature_RNA <= 10000 & percent.mt <= 20)

# HLCA dataset
HLCA_raw.data = read.csv('/media/yschang/SMG1/Lung_ATLAS/krasnow_hlca_10x_UMIs.csv', row.names = 1)
HLCA.sbj <- CreateSeuratObject(counts = HLCA_raw.data, project = "HLCA", min.cells = 3, min.features = 200)

# Add Metadata
HLCA.sbj$dataset.id <- 'HLCA'
HLCA.sbj$sample.id <- substr(colnames(HLCA.sbj), 1, 4)
HLCA.sbj$orig.ident <- HLCA.sbj$sample.id
HLCA.sbj <- SetIdent(HLCA.sbj, value = 'orig.ident')

# QC with reference (1) percent.mito: na (2) UMIs count: >= 1,000, (3) gene count: >= 200
HLCA.sqj <- subset(x = HLCA.sbj, subset = nFeature_RNA >= 500 & nCount_RNA >= 1000)


## EBI_MTAB_6653
matrix_dir = "/media/yschang/SMG1/EBI_MTAB_6653/data/"
barcode.path <- paste0(matrix_dir, "E-MTAB-6653.aggregated_filtered_counts.mtx_cols"); features.path <- paste0(matrix_dir, "genes.tsv"); matrix.path <- paste0(matrix_dir, "E-MTAB-6653.aggregated_filtered_counts.mtx")
EBI_MTAB_6653_raw.data <- readMM(file = matrix.path)
feature.names = read.delim(features.path, header = FALSE,stringsAsFactors = FALSE); barcode.names = read.delim(barcode.path, header = FALSE,stringsAsFactors = FALSE)
colnames(EBI_MTAB_6653_raw.data) = barcode.names$V1; rownames(EBI_MTAB_6653_raw.data) = feature.names$V1
EBI_MTAB_6653.sbj <- CreateSeuratObject(counts = EBI_MTAB_6653_raw.data, project = "EBI_MTAB_6653", min.cells = 3, min.features = 200); rm(EBI_MTAB_6653_raw.data, matrix_dir, matrix.path, barcode.path, barcode.names, features.path, feature.names)

# Metadata barcoding
EBI_MTAB_6653.sbj$dataset.id <- Idents(EBI_MTAB_6653.sbj)
EBI_MTAB_6653.sbj$sample.id <- substr(colnames(EBI_MTAB_6653.sbj), 1, 10)
EBI_MTAB_6653.sbj$orig.ident <- EBI_MTAB_6653.sbj$sample.id
EBI_MTAB_6653.sbj <- SetIdent(EBI_MTAB_6653.sbj, value = 'orig.ident')

# QC with reference (1) percent.mito:<= 10, (2) UMIs count: 200 ~ , (3) gene count: 200 ~ 6,000
# Get mito percentage
EBI_MTAB_6653.sbj[["percent.mt"]] <- PercentageFeatureSet(EBI_MTAB_6653.sbj, pattern = "^MT-")
# Get sqj
EBI_MTAB_6653.sqj <- subset(x = EBI_MTAB_6653.sbj, subset = nCount_RNA >= 200 & nFeature_RNA >= 100 & nFeature_RNA < 6000 & percent.mt <= 10)

## KJA dataset 
# NL
KJA_NL_raw.data <- Read10X(data.dir = "/media/yschang/T/KJA/A_data_NL/filtered_feature_bc_matrix/")
KJA_NL.sbj <- CreateSeuratObject(counts = KJA_NL_raw.data, project = "Yonsei_01_NL", min.cells = 3, min.features = 200); rm(KJA_NL_raw.data)
# Tu
KJA_Tu_raw.data <- Read10X(data.dir = "/media/yschang/T/KJA/A_data_Tu/filtered_feature_bc_matrix/")
KJA_Tu.sbj <- CreateSeuratObject(counts = KJA_Tu_raw.data, project = "Yonsei_01_Tu", min.cells = 3, min.features = 200); rm(KJA_Tu_raw.data)

### LYJ dataset
# NL
LYJ_NL_raw.data <- Read10X(data.dir = "/media/yschang/T/LYJ/A_data_NL/filtered_feature_bc_matrix/")
LYJ_NL.sbj <- CreateSeuratObject(counts = LYJ_NL_raw.data, project = "Yonsei_02_NL", min.cells = 3, min.features = 200); rm(LYJ_NL_raw.data)
# Tu
LYJ_Tu_raw.data <- Read10X(data.dir = "/media/yschang/T/LYJ/A_data_Tu/filtered_feature_bc_matrix/")
LYJ_Tu.sbj <- CreateSeuratObject(counts = LYJ_Tu_raw.data, project = "Yonsei_02_Tu", min.cells = 3, min.features = 200); rm(LYJ_Tu_raw.data)

## NHK dataset 
# NL
NKH_NL_raw.data <- Read10X(data.dir = "/media/yschang/S/NKH/A_data_NL/filtered_feature_bc_matrix/")
NKH_NL.sbj <- CreateSeuratObject(counts = NKH_NL_raw.data, project = "Yonsei_03_NL", min.cells = 3, min.features = 200); rm(NKH_NL_raw.data)
# Tu
NKH_Tu_raw.data <- Read10X(data.dir = "/media/yschang/S/NKH/A_data_Tu/filtered_feature_bc_matrix/")
NKH_Tu.sbj <- CreateSeuratObject(counts = NKH_Tu_raw.data, project = "Yonsei_03_Tu", min.cells = 3, min.features = 200); rm(NKH_Tu_raw.data)

## YJO dataset
# NL
YJO_NL_raw.data <- Read10X(data.dir = "/media/yschang/S/YJO/A_data_NL/filtered_feature_bc_matrix/")
YJO_NL.sbj <- CreateSeuratObject(counts = YJO_NL_raw.data, project = "Yonsei_04_NL", min.cells = 3, min.features = 200); rm(YJO_NL_raw.data)
# Tu
YJO_Tu_raw.data <- Read10X(data.dir = "/media/yschang/S/YJO/A_data_Tu/filtered_feature_bc_matrix/")
YJO_Tu.sbj <- CreateSeuratObject(counts = YJO_Tu_raw.data, project = "Yonsei_04_Tu", min.cells = 3, min.features = 200); rm(YJO_Tu_raw.data)

## LDS dataset
# NL
LDS_NL_raw.data <- Read10X(data.dir = "/media/yschang/S/LDS/A_data_NL/filtered_feature_bc_matrix/")
LDS_NL.sbj <- CreateSeuratObject(counts = LDS_NL_raw.data, project = "Yonsei_05_NL", min.cells = 3, min.features = 200); rm(LDS_NL_raw.data)
# Tu
LDS_Tu_raw.data <- Read10X(data.dir = "/media/yschang/S/LDS/A_data_Tu/filtered_feature_bc_matrix/")
LDS_Tu.sbj <- CreateSeuratObject(counts = LDS_Tu_raw.data, project = "Yonsei_05_Tu", min.cells = 3, min.features = 200); rm(LDS_Tu_raw.data)

## SSS_Tu dataset
SSS_Tu_raw.data <- Read10X(data.dir = "/media/yschang/S/SSS/A_data_Tu/filtered_feature_bc_matrix/")
SSS_Tu.sbj <- CreateSeuratObject(counts = SSS_Tu_raw.data, project = "Yonsei_06_Tu", min.cells = 3, min.features = 200); rm(SSS_Tu_raw.data)

## KBY dataset 
# NL
KBY_NL_raw.data <- Read10X(data.dir = "/media/yschang/T/KBY/A_data_NL/filtered_feature_bc_matrix/")
KBY_NL.sbj <- CreateSeuratObject(counts = KBY_NL_raw.data, project = "Yonsei_07_NL", min.cells = 3, min.features = 200); rm(KBY_NL_raw.data)
# Tu
KBY_Tu_raw.data <- Read10X(data.dir = "/media/yschang/T/KBY/A_data_Tu/filtered_feature_bc_matrix/")
KBY_Tu.sbj <- CreateSeuratObject(counts = KBY_Tu_raw.data, project = "Yonsei_07_Tu", min.cells = 3, min.features = 200); rm(KBY_Tu_raw.data)

# Merge Yonsei dataset
Yonsei.sbj <- merge(KJA_NL.sbj, y = c(KJA_Tu.sbj, LYJ_NL.sbj, LYJ_Tu.sbj, NKH_NL.sbj, NKH_Tu.sbj, YJO_NL.sbj, YJO_Tu.sbj, LDS_NL.sbj, LDS_Tu.sbj, SSS_Tu.sbj,  KBY_NL.sbj, KBY_Tu.sbj))

# Add Metadata
Yonsei.sbj$sample.id <- Idents(Yonsei.sbj)
Yonsei.sbj$orig.ident <- Yonsei.sbj$sample.id
Yonsei.sbj$dataset.id <- 'Yonsei'

rm(KJA_NL.sbj, KJA_Tu.sbj, LYJ_NL.sbj, LYJ_Tu.sbj, KBY_NL.sbj, KBY_Tu.sbj, SSS_Tu.sbj, YJO_NL.sbj, YJO_Tu.sbj, NKH_NL.sbj, NKH_Tu.sbj, LDS_NL.sbj, LDS_Tu.sbj)

# QC with reference (1) percent.mito:<= 20, (2) UMIs count:100 ~ 150,000, (3) gene count:200 ~ 10,000
# Get mito percentage
Yonsei.sbj[["percent.mt"]] <- PercentageFeatureSet(Yonsei.sbj, pattern = "^MT-")
# Get sqj
Yonsei.sqj <- subset(x = Yonsei.sbj, subset = nCount_RNA >= 100 & nCount_RNA <= 150000 & nFeature_RNA >= 200 & nFeature_RNA <= 10000 & percent.mt <= 20)

# Merge by common features
#1. Extract common features
common.features_1 <- intersect(rownames(Yonsei.sqj), rownames(GSE131907.sqj))
common.features_2 <- intersect(common.features_1, rownames(EBI_MTAB_6653.sqj))
common.features <- intersect(common.features_2, rownames(HLCA.sqj))
rm(common.features_1, common.features_2)

#2. Filtering data by common features
EBI_MTAB_6653_common.sqj <- EBI_MTAB_6653.sqj[common.features,]
GSE131907_common.sqj <- GSE131907.sqj[common.features,]
Yonsei_common.sqj <- Yonsei.sqj[common.features,]
HLCA_common.sqj <- HLCA.sqj[common.features,]

#3. merge data
Total.sqj <- merge(Yonsei_common.sqj, c(EBI_MTAB_6653_common.sqj, GSE131907_common.sqj, HLCA_common.sqj), add.cell.ids = c('Yonsei', 'EBI_MTAB_6653', 'GSE131907', 'HLCA')); rm(Yonsei_common.sqj, EBI_MTAB_6653_common.sqj, GSE131907_common.sqj, HLCA_common.sqj)

#4. add metadata
# gender.id
gender <- unique(Metadata$gender.id)
for (idx in gender){
  a <- WhichCells(Total.sqj, idents = row.names(Metadata %>% filter(gender.id == idx))); Idents(Total.sqj, cells = a) <- idx
}
Idents(Total.sqj)
Total.sqj$gender.id <- Idents(Total.sqj)
Idents(Total.sqj) <- Total.sqj$orig.ident

# patient.id
patient <- unique(Metadata$patient.id)
for (idx in patient){
  a <- WhichCells(Total.sqj, idents = row.names(Metadata %>% filter(patient.id == idx))); Idents(Total.sqj, cells = a) <- idx
}
unique(Idents(Total.sqj))
Total.sqj$patient.id <- Idents(Total.sqj)
Idents(Total.sqj) <- Total.sqj$orig.ident

# tissue.id
tissue <- unique(Metadata$tissue.id)
for (idx in tissue){
  a <- WhichCells(Total.sqj, idents = row.names(Metadata %>% filter(tissue.id == idx))); Idents(Total.sqj, cells = a) <- idx
}
unique(Idents(Total.sqj))
Total.sqj$tissue.id <- Idents(Total.sqj)
Idents(Total.sqj) <- Total.sqj$orig.ident

# pathology.id
pathology <- unique(Metadata$pathology.id)
for (idx in pathology){
  a <- WhichCells(Total.sqj, idents = row.names(Metadata %>% filter(pathology.id == idx))); Idents(Total.sqj, cells = a) <- idx
}
unique(Idents(Total.sqj))
Total.sqj$pathology.id <- Idents(Total.sqj)
Idents(Total.sqj) <- Total.sqj$orig.ident

# smoking.id
smoking <- unique(Metadata$smoking.id)
for (idx in smoking){
  a <- WhichCells(Total.sqj, idents = row.names(Metadata %>% filter(smoking.id == idx))); Idents(Total.sqj, cells = a) <- idx
}
unique(Idents(Total.sqj))
Total.sqj$smoking.id <- Idents(Total.sqj)
Idents(Total.sqj) <- Total.sqj$orig.ident

# EGFR.id
EGFR <- unique(Metadata$EGFR.id)
for (idx in EGFR){
  a <- WhichCells(Total.sqj, idents = row.names(Metadata %>% filter(EGFR.id == idx))); Idents(Total.sqj, cells = a) <- idx
}
unique(Idents(Total.sqj))
Total.sqj$EGFR.id <- Idents(Total.sqj)
Idents(Total.sqj) <- Total.sqj$orig.ident

# stage.id
stage <- unique(Metadata$stage.id)
for (idx in stage){
  a <- WhichCells(Total.sqj, idents = row.names(Metadata %>% filter(stage.id == idx))); Idents(Total.sqj, cells = a) <- idx
}
unique(Idents(Total.sqj))
Total.sqj$stage.id <- Idents(Total.sqj)
Idents(Total.sqj) <- Total.sqj$orig.ident

# save
save.image(paste0(Merge.outdir,'Baseline.RData'))
saveRDS(Total.sqj, file = paste0(Merge.outdir,'Total.sqj'), compress = TRUE)

rm(list = ls())
gc()
q()
