library(Seurat)
library(UCell)
library(scCustomize)
library(tidyr)
library(rstatix)
source("script/Seurat_functions.R")
MG = readRDS("data/MG.rds")


DAM_enrich = scan("data/DAM_enriched.txt", "character")
DAM_depleted = scan("data/DAM_depleted.txt", "character")

# score normal and DAM
gene_list = list(DAM_enrich,DAM_depleted)
scores <- ScoreSignatures_UCell(MG@assays$RNA@counts, features=gene_list)
metadf = MG@meta.data
metadf$BC = rownames(metadf)
scores = as.data.frame(scores)
scores$BC = rownames(scores)
colnames(scores) = c("DAM_enrich","DAM_depleted","BC")

metadf = left_join(metadf, scores, by = "BC")
rownames(metadf) = metadf$BC
MG@meta.data = metadf

# score DAM1 and DAM2
DAM1_enrich = scan("data/DAM1_enriched.txt", "character")
DAM2_enrich = scan("data/DAM2_enriched.txt", "character")

gene_list = list(DAM1_enrich,DAM2_enrich)
scores <- ScoreSignatures_UCell(MG@assays$RNA@counts, features=gene_list)
metadf = MG@meta.data
metadf$BC = rownames(metadf)
scores = as.data.frame(scores)
scores$BC = rownames(scores)
colnames(scores) = c("DAM1_enrich","DAM2_enrich","BC")

metadf = left_join(metadf, scores, by = "BC")
rownames(metadf) = metadf$BC
MG@meta.data = metadf

FeaturePlot_scCustom(MG,"DAM_enrich")
FeaturePlot_scCustom(MG,"DAM_depleted")
FeatureScatter(MG, feature1 = "DAM_enrich", feature2 = "DAM_depleted", group.by = "genotype")
FeaturePlot(MG, features = c("DAM_enrich", "DAM_depleted"),blend = T, pt.size = .5,cols = c("grey", "red", "green"), order = T, blend.threshold = 0.5)
#ggsave("plots/MG_pca/MG_umap_DAM_enrich_depleted_blend.jpeg", width = 20, height = 5, units = "in", dpi = 300)

Stacked_VlnPlot(MG, features = c("DAM_enrich","DAM_depleted"),pt.size = 0, group.by = "genotype")
#ggsave("plots/MG_Vln_DAM_enrich_depleted.jpeg", width = 5, height = 12)

out = BlendExpression(MG@meta.data[,c("DAM_enrich", "DAM_depleted")])
bc = rownames(out)
out = apply(out, 2, as.numeric)
rownames(out) = bc
out = as.data.frame(out)
tg = rownames(out)[out$DAM_enrich>=6&out$DAM_depleted<=5]
DimPlot(MG, cells.highlight = tg, sizes.highlight = 0.2)+NoLegend()+ggtitle("DAM")
#ggsave("plots/MG_pca/MG_umap_DAM_enrich_blend_classified.jpeg", width = 5, height = 5, units = "in", dpi = 300)


FeaturePlot_scCustom(MG, "DAM1_enrich")
FeaturePlot_scCustom(MG, "DAM2_enrich")
FeaturePlot(MG, features = c("DAM1_enrich", "DAM2_enrich"),blend = T, pt.size = .5,cols = c("grey", "red", "green"), order = T, blend.threshold = 0.1)
#ggsave("plots/MG_pca/MG_umap_DAM1_DAM2_enrich_blend.jpeg", width = 20, height = 5, units = "in", dpi = 300)

out2 = BlendExpression(MG@meta.data[,c("DAM1_enrich", "DAM2_enrich")])
bc = rownames(out2)
out2 = apply(out2, 2, as.numeric)
rownames(out2) = bc
out2 = as.data.frame(out2)
tg2 = tg[tg%in%rownames(out2)[out2$DAM2_enrich>=6&out2$DAM1_enrich<=5]]
DimPlot(MG, cells.highlight = tg2, sizes.highlight = 0.2)+NoLegend()+ggtitle("DAM2")

MG$classfication = "DAM_depleted"
metadf = MG@meta.data
tg2 = tg2[tg2%in%tg]
metadf[tg,]$classfication = "DAM1"
metadf[tg2,]$classfication = "DAM2"
MG@meta.data = metadf
DimPlot(MG, group.by = "classfication")
#ggsave("plots/MG_pca/MG_umap_DAM_classification.jpeg", width = 7, height = 5, units = "in", dpi = 300)

saveRDS(MG, "data/MG.rds")
