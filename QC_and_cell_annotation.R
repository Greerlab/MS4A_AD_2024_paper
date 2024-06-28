library(Seurat)
library(DoubletFinder)
library(scCustomize)
library(harmony)
library(patchwork)
library(SingleR)
library(dplyr)
# On cluster, the cellbender seems not able to remove background of the 6CMT4 well,  
# so I screened through the different Soupx setting and 
# plot the most dominate cell type marker to determine the rate for decontamination

samples = scan("../script/sample_name", what = "chracter")
d = rep(0.2, 14)
d[samples=="MAPT5"] = 0.25
d[samples=="6CMT4"] = 0.55

input_files = paste0(samples,"_",d,".rds")

# Doubletfinder and merge
input_files_list = list()
for (i in 1:length(input_files)) {
  obj = readRDS(input_files[i])
  colnames(obj) = gsub("^",paste0(samples[i],"_"),colnames(obj))
  obj = CreateSeuratObject(obj, names.field = 1, names.delim = "_")
  obj[["percent.mt"]] = PercentageFeatureSet(obj, pattern = "^mt-")
  obj[["percent.hbb"]] <- PercentageFeatureSet(obj, pattern = "^Hbb")
  obj <- subset(obj, subset = nFeature_RNA > 300 & nCount_RNA >1000 & percent.mt < 20 & percent.hbb < 1)
  obj = NormalizeData(obj)
  obj = FindVariableFeatures(obj)
  obj = ScaleData(obj)
  obj = RunPCA(obj)
  sweep.res.obj = paramSweep_v3(obj, PCs = 1:30, sct = FALSE)
  sweep.stats_obj <- summarizeSweep(sweep.res.obj, GT = FALSE)
  bcmvn <- find.pK(sweep.stats_obj)
  ggplot(bcmvn, aes(pK, BCmetric, group = 1))+ geom_point()+geom_line()
  x = as.numeric(levels(bcmvn$pK)[which.max(bcmvn$BCmetric)])
  nExp = ncol(obj)*0.016 # estimated doublet rate
  obj_db = doubletFinder_v3(obj, PCs = 1:30, pN = 0.25, pK = x, nExp = nExp, sct =F)
  metadf = obj_db@meta.data
  colnames(metadf)[grep("pANN|DF",colnames(metadf))] = c("db_score", "db_classification")
  obj_db@meta.data = metadf
  obj_db = RunUMAP(obj_db, dims = 1:30)
  DimPlot(obj_db, group.by = "db_classification")
  input_files_list[[i]] = obj_db
}
merged.obj = Merge_Seurat_List(input_files_list)

# filtering
FeatureScatter(merged.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size = 0)
metadf = merged.obj@meta.data
ggplot(metadf, aes(x=nCount_RNA)) + geom_histogram( binwidth = 100)+geom_vline(xintercept = 4000)
summary(metadf)
merged.obj <- subset(merged.obj, 
                     subset = nFeature_RNA >= 800 &
                       nFeature_RNA <= 7000 &
                       nCount_RNA >1200 &
                       nCount_RNA <= 25000 &
                       percent.mt < 20 &
                       percent.hbb < 1)
Idents(merged.obj) = merged.obj$db_classification
merged.obj = subset(merged.obj, idents = "Singlet")

# Dimension reduction and batch correction
merged.obj = NormalizeData(merged.obj)
merged.obj = FindVariableFeatures(merged.obj)
merged.obj = ScaleData(merged.obj)
merged.obj = RunPCA(merged.obj)
ElbowPlot(merged.obj, ndims = 50)
merged.obj = RunUMAP(merged.obj, dims = 1:30)
Idents(merged.obj) = merged.obj$orig.ident
DimPlot(merged.obj, label = T)
merged.obj = RunHarmony(merged.obj, group.by.vars = "orig.ident")
merged.obj = RunUMAP(merged.obj, dims = 1:30, reduction = "harmony")
DimPlot(merged.obj, label = T)
merged.obj = FindNeighbors(merged.obj, dims = 1:2, reduction = "umap")
merged.obj = FindClusters(merged.obj, resolution = 0.15)
DimPlot(merged.obj, label = T)

# Cell identity SingleR
ref = celldex::MouseRNAseqData()

FeaturePlot_scCustom(merged.obj, "percent.mt")
FeaturePlot_scCustom(merged.obj, "db_score")
merged.obj_counts = LayerData(merged.obj,layer = "counts")
pred = SingleR(merged.obj_counts, ref = ref, labels = ref$label.fine)
merged.obj$cell_type_pred = pred$labels[match(colnames(merged.obj), rownames(pred))]
merged.obj$cell_type_pruned = pred$pruned.labels[match(colnames(merged.obj), rownames(pred))]
DimPlot(merged.obj, group.by = "cell_type_pred",label = T)
DimPlot(merged.obj, group.by = "cell_type_pruned",label = T)

#
Genes = c("Pdgfra","Cldn11","Aldh1l1","Cdk1","Sox11","Snap25","Ccdc153",
          "Ttr","Cldn5", "Kcnj8", "Acta2", "Slc47a1", "Tmem119", "Plac8",
          "Pf4", "Cd209a", "S100a9", "Ms4a1", "Cd3g")
Abbreviation = c("OPC","OLG","ASC","NRP","ImmN","mNEUR","EPC","CPC","EC","PC","VSMC","ABC","MC","MNC","MAC","DC","NEUT","BC","TC")
Cell_type = c("Oligodendrocyte precursor cell", "Oligodendrocyte", "Astrocyte", "Neuronal restricted precursor","Immature neurons",
              "Mature neurons","Ependymocytes","Choroid plexus epithelial cells","Endothelial cells","Pericytes","Vascular smooth muscle cells",
              "Arachnoid barrier cells","Microglia","Monocytes","Marcophage","Dendritic cells","Neutrophils","B lymphocyte","T lymphocyte")
celltype = cbind.data.frame(Genes,Abbreviation,Cell_type)

plots = list()
for (i in 1:nrow(celltype)) {
  plots[[i]] = FeaturePlot_scCustom(merged.obj, paste0(celltype$Genes[i]), order = T)+NoLegend()+NoAxes()+labs(title = celltype$Cell_type[i], subtitle = celltype$Genes[i])+
    theme(plot.subtitle = element_text(hjust = 0.5),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(colour = "snow3", size=2))
}
wrap_plots(plots, ncol = 5)
ggsave("plots/merged.obj_cell_types.jpeg", width = 20, height = 16, units = "in", dpi = 300)

# idenify marker genes on hpcc
saveRDS(merged.obj, "../Cell_type_analysis/merged.obj.rds")

markers = list.files("../Cell_type_analysis", ".csv", full.names = T)
marker_list = list()
for (i in 1:length(markers)) {
  tmp = read.csv(markers[i])
  tmp$cluster = gsub("../Cell_type_analysis/marker_|.csv","",markers[i])
  marker_list[[i]] = tmp
}

marker_list = do.call(rbind.data.frame, marker_list)

colnames(marker_list)[1] = "gene"
marker_list[marker_list$cluster==21,]

# 25, 27 Fibroblasts
# 20 Choroid Plexus Cells: (Aqp1, Ttr)
# 10, 13  Astrocyte (reactive and normal)
# 2, 14 Endothelial cells
# 21 doublets of Endothelial cells and Pericyte
# 12 Pericyte
# 19 Vascular smooth muscle cells
# 8, 7 Microglia
# 23 Arachnoid barrier cell
# 24 T cell
# 11, 6, 18, 26, 16 ,1, 3, 4, 9 Neuron
# 22 Ependymocytes
merged.obj = FindSubCluster(merged.obj, cluster = 17, subcluster.name = "Subcluster", graph.name = "RNA_snn")

# subcluster 17
Subcluster = subset(merged.obj, idents = 17) %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA() %>% 
  RunHarmony(group.by.vars = "orig.ident") %>% 
  RunUMAP(reduction = "harmony", dims = 1:30) %>%
  FindNeighbors(reduction = "harmony") %>% FindClusters()

DimPlot_scCustom(Subcluster, group.by = "cell_type_pred", label = T, split.by = "cell_type_pred")

DimPlot_scCustom(Subcluster, label = T)
FeaturePlot_scCustom(Subcluster, "Ms4a1") # 4 B cell
FeaturePlot_scCustom(Subcluster, "Pf4") # 2, 8, 7, 0 Marcophage
FeaturePlot_scCustom(Subcluster, "S1009a") # 5 Neutrophil
FeaturePlot_scCustom(Subcluster, "Plac8") # 9, 10 monocyte
FeaturePlot_scCustom(Subcluster, "Klrd1") #  3 Dendritic cell
FeaturePlot_scCustom(Subcluster, "Tmem119") # 1 microglia
FeaturePlot_scCustom(Subcluster, "percent.mt") # 6 low quality cells with high MT expression and express multiple markers

cells2rmove = c(WhichCells(Subcluster, idents = "6"), WhichCells(merged.obj, idents = 21))
metadf = merged.obj@meta.data
metadf$cell_type = "Neurons"
metadf[metadf$seurat_clusters%in%c(0,5,15),]$cell_type = "Oligodendrocytes"
metadf[metadf$seurat_clusters%in%c(10,13),]$cell_type = "Astrocytes"
metadf[metadf$seurat_clusters%in%c(2,14),]$cell_type = "Endothelial cells"
metadf[metadf$seurat_clusters%in%c(12),]$cell_type = "Pericytes"
metadf[metadf$seurat_clusters%in%c(19),]$cell_type = "Vascular smooth muscle cells"
metadf[metadf$seurat_clusters%in%c(7, 8),]$cell_type = "Microglia"
metadf[metadf$seurat_clusters%in%c(23),]$cell_type = "Arachnoid barrier cell"
metadf[metadf$seurat_clusters%in%c(24),]$cell_type = "T lymphocytes"
metadf[metadf$seurat_clusters%in%c(25,27),]$cell_type = "Fibroblasts"
metadf[metadf$seurat_clusters%in%c(20),]$cell_type = "Choroid plexus cells"
metadf[metadf$seurat_clusters%in%c(22),]$cell_type = "Ependymocytes"
metadf[rownames(metadf)%in%WhichCells(Subcluster,idents = 4),]$cell_type = "B lymphocytes"
metadf[rownames(metadf)%in%WhichCells(Subcluster,idents = c(2, 8, 7, 0)),]$cell_type = "Marcophages"
metadf[rownames(metadf)%in%WhichCells(Subcluster,idents = 3),]$cell_type = "Dendritic cells"
metadf[rownames(metadf)%in%WhichCells(Subcluster,idents = 1),]$cell_type = "Microglia"
metadf[rownames(metadf)%in%WhichCells(Subcluster,idents = 5),]$cell_type = "Neutrophils"
metadf[rownames(metadf)%in%cells2rmove,]$cell_type = "Doublets"

merged.obj@meta.data = metadf
DimPlot_scCustom(merged.obj, group.by = "cell_type", label = T)


# remove doublets and rerun the pipeline
Idents(merged.obj) = merged.obj$cell_type
merged.obj = subset(merged.obj, ident = "Doublets", invert = T)

merged.obj = NormalizeData(merged.obj)
merged.obj = FindVariableFeatures(merged.obj)
merged.obj = ScaleData(merged.obj)
merged.obj = RunPCA(merged.obj)
ElbowPlot(merged.obj, ndims = 50)
merged.obj = RunUMAP(merged.obj, dims = 1:30)
Idents(merged.obj) = merged.obj$orig.ident
DimPlot(merged.obj, label = T)
merged.obj = RunHarmony(merged.obj, group.by.vars = "orig.ident")
merged.obj = RunUMAP(merged.obj, dims = 1:30, reduction = "harmony")
DimPlot(merged.obj, label = T)
merged.obj = FindNeighbors(merged.obj, dims = 1:2, reduction = "umap")
merged.obj = FindClusters(merged.obj, resolution = 0.1)
DimPlot(merged.obj, label = T)

table(merged.obj$cell_type[merged.obj$seurat_clusters==13])
metadf = merged.obj@meta.data
metadf[metadf$seurat_clusters==13&metadf$cell_type=="Microglia",]$cell_type = "Marcophages"
metadf[metadf$seurat_clusters==18,]$cell_type = "Arachnoid barrier cell"
metadf[metadf$seurat_clusters==20,]$cell_type = "Unknown"
metadf[metadf$seurat_clusters==17,]$cell_type = "Ependymocytes"
metadf[metadf$seurat_clusters==16,]$cell_type = "Choroid plexus cells"
metadf[metadf$seurat_clusters==19,]$cell_type = "T lymphocytes"
metadf[metadf$seurat_clusters==21,]$cell_type = "Neurons"
metadf[metadf$seurat_clusters==24,]$cell_type = "Fibroblasts"
metadf[metadf$seurat_clusters==11,]$cell_type = "Microglia"
metadf[metadf$seurat_clusters%in%c(8,22,14,12,6,21,0,4),]$cell_type = "Neurons"
metadf[metadf$seurat_clusters==7,]$cell_type = "Astrocytes"
metadf[metadf$seurat_clusters%in%c(9),]$cell_type = "Pericytes"
metadf[metadf$seurat_clusters%in%c(15),]$cell_type = "Vascular smooth muscle cells"
metadf[metadf$seurat_clusters%in%c(2),]$cell_type = "Endothelial cells"
merged.obj@meta.data = metadf
DimPlot(merged.obj, group.by = "cell_type", label = T)

DimPlot(merged.obj, group.by = "seurat_clusters", label = T)

merged.obj$genotype = gsub("[0-9]$","",merged.obj$orig.ident)
saveRDS(merged.obj, "data/merged.obj.rds")

# Microglia
Idents(merged.obj) = merged.obj$cell_type
MG = subset(merged.obj, idents = "Microglia")
MG = NormalizeData(MG)
MG = FindVariableFeatures(MG)
MG = ScaleData(MG)
MG = RunPCA(MG)
MG = RunHarmony(MG,group.by.vars = "orig.ident")
MG = RunUMAP(MG, reduction = "harmony",  dims = 1:10)
MG = FindNeighbors(MG, reduction = "umap",  dims = 1:2)
MG = FindClusters(MG,resolution = 0.1)
DimPlot(MG, label = T)
FeaturePlot(MG, "Top2a")
prol = WhichCells(MG, expression = Top2a > 1)
DimPlot(MG, cells.highlight = prol, sizes.highlight = 0.1)
cls6 = FindMarkers(MG, ident.1 = 6)

MG = subset(MG, cells = prol, invert = T)
MG = FindVariableFeatures(MG)
MG = ScaleData(MG)
MG = RunPCA(MG)
MG = RunHarmony(MG,group.by.vars = "orig.ident")
MG = RunUMAP(MG, reduction = "harmony",  dims = 1:6)
MG = FindNeighbors(MG, reduction = "umap",  dims = 1:2)
MG = FindClusters(MG,resolution = 0.1)
DimPlot(MG, label = T)
MG$genotype = factor(MG$genotype, levels = c("WT","6CKO","MAPT","6CMT"))
DimPlot(MG, group.by = "genotype", split.by = "genotype", ncol = 2)+NoLegend()
saveRDS(MG, "data/MG.rds")




