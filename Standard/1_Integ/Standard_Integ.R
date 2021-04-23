# Par B. Integration of Dataset 

### SET-UP Codes

# Clear workspace
rm(list = ls()); gc()

# get required library
library(Seurat); library(dplyr); library(Matrix)

#	Set output directory 
ITG.outdir = "/media/yschang/SMG1/ITG/Standard/1_Integ/"

# load previous work
# load('/media/yschang/SMG1/ITG/Merge/Baseline.RData')
Total.sqj <- readRDS('/media/yschang/SMG1/ITG/Merge/Total.sqj')

# get NL and Tu only
Idents(Total.sqj) <- Total.sqj$tissue.id
Merge.sqj <- subset(Total.sxj, idents =c('NL','Tu')); rm(Total.sqj)

# 1. Dataset Preprocessing and normalization
Total_dataset.list <- SplitObject(Merge.sqj, split.by = 'dataset.id')
Total_dataset.list <- lapply(Total_dataset.list, FUN = function(x){
  x<-NormalizeData(x)
  x<-FindVariableFeatures(x, selection.method = 'vst', nfeatures=2000)
})
Total.anchors <- FindIntegrationAnchors(object.list = Total_dataset.list, dims=1:30)

# 2. Generate integrated seurat.obj 
ITG.sbj <- IntegrateData(anchorset=Total.anchors, dims=1:30); rm(Total_dataset.list, Total.anchors)
DefaultAssay(ITG.sbj) <- 'integrated'

# 3. Scaling the data
ITG.ssj<- ScaleData(object = ITG.sbj, features = rownames(ITG.sbj))

# 4. Dimension Reduction
## Linelar Dimension reduction
ITG.spj <- RunPCA(object = ITG.ssj, npcs = 30, verbose = FALSE); rm(object = ITG.ssj)
print(ITG.spj[["pca"]], dims = 1:5, nfeatures=5)
VizDimLoadings(ITG.spj, dims = 1:2, reduction ="pca")
DimPlot(ITG.spj, reduction = 'pca')
DimHeatmap(ITG.spj, dims = 1, cells = 500, balanced =TRUE)
DimHeatmap(ITG.spj, dims = 1:15, cells = 500, balanced =TRUE)

## NON-LINEAR DIMENSION REDUCTION 
ElbowPlot(object = ITG.spj, reduction="pca")
ITG.spj <- JackStraw(object = ITG.spj, num.replicate = 100)
ITG.spj <- ScoreJackStraw(object = ITG.spj, dims = 1:20)
JackStrawPlot(object = ITG.spj, dims = 1:20)

# 5. Clustering analysis: 
ITG.spj <- FindNeighbors(object = ITG.spj, dims = 1:2)
ITG.spj <- FindClusters(object = ITG.spj, resolution = 0.1)
ITG.suj<- RunUMAP(object = ITG.spj, reduction = 'pca',  dims = 1:5) 

# 6. Visualize and save
DimPlot(object = ITG.suj, reduction = "umap", label = T)
DimPlot(object = ITG.suj, group.by= 'orig.ident', label = T)
DimPlot(object = ITG.suj, split.by = 'dataset.id', label = T)

# fig1
jpeg(filename = paste0(ITG.outdir, 'ITG.suj.jpeg'), width = 1200, height = 900, res=600)
DimPlot(object = ITG.suj, group.by= 'orig.ident', label = T)
dev.off()
# fig2
jpeg(filename = paste0(ITG.outdir, 'ITG.suj_by_dataset.jpeg'), width = 3000, height = 1200, res=300)
DimPlot(object = ITG.suj,  split.by = 'dataset.id', label = T)
dev.off()

# save
save.image(paste0(ITG.outdir, 'ITG_Std.RData'))
saveRDS(ITG.suj, paste0(ITG.outdir, 'ITG.suj'), compression =T)
rm(list = ls());gc()
q()
