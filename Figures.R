library(Seurat)
library(ggpubr)
library(rstatix)
library(scales)
library(scCustomize)
library(patchwork)
library(harmony)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(topGO)
library(stringr)
library(EnhancedVolcano)
library(openxlsx)
library(reshape2)
#library(ggalt)
library(ggforce)
source("script/Seurat_functions.R")

merged.obj = readRDS("data/merged.obj.rds")
MG = readRDS("data/MG.rds")
MG$classification = MG$classfication
MG$classification[MG$classification%in%c("DAM1","DAM2")]="DAM"
merged.obj$genotype = gsub("WT","Ntg-WT",merged.obj$genotype)
merged.obj$genotype = gsub("6CKO","Ntg-6C KO",merged.obj$genotype)
merged.obj$genotype = gsub("MAPT","Tg-WT",merged.obj$genotype)
merged.obj$genotype = gsub("6CMT","Tg-6C KO",merged.obj$genotype)
merged.obj$genotype = factor(merged.obj$genotype, levels = c("Ntg-WT","Ntg-6C KO","Tg-WT","Tg-6C KO"))

MG$genotype = gsub("WT","Ntg-WT",MG$genotype)
MG$genotype = gsub("6CKO","Ntg-6C KO",MG$genotype)
MG$genotype = gsub("MAPT","Tg-WT",MG$genotype)
MG$genotype = gsub("6CMT","Tg-6C KO",MG$genotype)
MG$genotype = factor(MG$genotype, levels = c("Ntg-WT","Ntg-6C KO","Tg-WT","Tg-6C KO"))

# 3A
p = DimPlot(merged.obj, group.by = "cell_type", label = T)+NoLegend()+ggtitle("")+xlab("UMAP1")+ylab("UMAP2")
ggsave(plot = p, "plots/Figure3A.pdf", width = 7, height = 7, units = "in")
DimPlot_scCustom(merged.obj, group.by = "seurat_clusters", label = T)

# 3B

## remove Ms4a6c
genes = grep("Ms4a6c", rownames(merged.obj), invert = T)
merged.obj2 <- subset(merged.obj, features = genes)
merged.obj2 = NormalizeData(merged.obj2)
merged.obj2 = FindVariableFeatures(merged.obj2)
merged.obj2 = ScaleData(merged.obj2)
merged.obj2 = RunPCA(merged.obj2)
ElbowPlot(merged.obj2, ndims = 50)
merged.obj2 = RunUMAP(merged.obj2, dims = 1:30)
merged.obj2 = RunHarmony(merged.obj2, group.by.vars = "orig.ident")
merged.obj2 = RunUMAP(merged.obj2, dims = 1:30, reduction = "harmony")
merged.obj2 = FindNeighbors(merged.obj2, dims = 1:2, reduction = "umap")
merged.obj2 = FindClusters(merged.obj2, resolution = 0.1)


cell_types = names(table(merged.obj2$cell_type)[table(merged.obj2$cell_type)>1000])
#cell_types = cell_types[1:5]
Idents(merged.obj2) = merged.obj2$cell_type
genotypes = unique(merged.obj2$genotype)
cell_type = c()
samples = c()
value =c()
for (s in 1:100) {
for (i in 1:length(cell_types)) {
  tmp = subset(merged.obj2, idents = cell_types[i])
  mingenotype = min(table(tmp$genotype))
  cell.list = c(sample(colnames(tmp)[tmp$genotype==genotypes[1]],mingenotype),
                sample(colnames(tmp)[tmp$genotype==genotypes[2]],mingenotype),
                sample(colnames(tmp)[tmp$genotype==genotypes[3]],mingenotype),
                sample(colnames(tmp)[tmp$genotype==genotypes[4]],mingenotype))
  tmp = FindNeighbors(tmp, reduction = "harmony", dims = 1:10, return.neighbor = T)
  NN = tmp@neighbors$RNA.nn@nn.idx[,c(1,2)]
  df = cbind.data.frame("Cell_geno" = tmp$genotype[NN[,1]],"Neighbor_geno"= tmp$genotype[NN[,2]], "Sample" = tmp$orig.ident[NN[,1]])
  df$same = ifelse(df$Cell_geno==df$Neighbor_geno,1,0)
  for (x in 1:length(unique(df$Sample))) {
    tmp = df[df$Sample==unique(df$Sample)[x],]
    cell_type = c(cell_type,cell_types[i])
    samples = c(samples,unique(df$Sample)[x])
    value = c(value,sum(tmp$same)/nrow(tmp))
  }
}}
results = cbind.data.frame(cell_type, samples, value)



cell_type = c()
samples = c()
value =c()
for (s in 1:100) {
for (i in 1:length(unique(cell_types))) {
  tmp = subset(merged.obj2, idents = cell_types[i])
  mingenotype = min(table(tmp$genotype))
  cell.list = c(sample(colnames(tmp)[tmp$genotype==genotypes[1]],mingenotype),
                sample(colnames(tmp)[tmp$genotype==genotypes[2]],mingenotype),
                sample(colnames(tmp)[tmp$genotype==genotypes[3]],mingenotype),
                sample(colnames(tmp)[tmp$genotype==genotypes[4]],mingenotype))
  tmp = subset(tmp, cells = cell.list)
  tmp = FindNeighbors(tmp, reduction = "harmony", dims = 1:10, return.neighbor = T)
  NN = tmp@neighbors$RNA.nn@nn.idx[,c(1,2)]
  df = cbind.data.frame("Cell_geno" = sample(tmp$genotype[NN[,1]]),"Neighbor_geno"= tmp$genotype[NN[,2]], "Sample" = tmp$orig.ident[NN[,1]])
  df$same = ifelse(df$Cell_geno==df$Neighbor_geno,1,0)
  for (x in 1:length(unique(df$Sample))) {
    tmp = df[df$Sample==unique(df$Sample)[x],]
    cell_type = c(cell_type,cell_types[i])
    samples = c(samples,unique(df$Sample)[x])
    value = c(value,mean(tmp$same))
  }
}
}
results2 = cbind.data.frame(cell_type, samples, value)

library(ggbeeswarm)
p = ggbarplot(results, x = "cell_type", y = "value", fill = "lightblue",
              title = "Segregation by genotype\n within each cell type", 
              ylab = "Proportion of nearest neighbors in same genotype",
              xlab = F, add = c("mean_se"))+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_y_continuous(limits=c(0,max(results$value)), expand = c(0, 0),breaks = scales::pretty_breaks(n = 10))+
  geom_hline(yintercept = mean(results2$value),linetype = "dashed")


# p = ggbarplot(results, x = "cell_type", y = "value", fill = "lightblue",
#               title = "Segregation by genotype\n within each cell type", 
#               ylab = "Proportion of nearest neighbors in same genotype",
#               xlab = F, add = c("mean_se","jitter"))+ geom_jitter(alpha = 0.2, width = 0.3)+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
#   scale_y_continuous(limits=c(min(results$value)-0.1,max(results$value)+0.05), oob=rescale_none)+
#   geom_hline(yintercept = mean(results2$value),linetype = "dashed")+geom_hline(yintercept = 0.25,linetype = "solid")
ggsave(plot = p, "plots/Figure3B_3.pdf", width = 4, height = 7, units = "in")

## just mapt and 6c;mapt
Idents(merged.obj2)  = merged.obj2$genotype
merged.obj2 <- subset(merged.obj2 , ident = c("Tg-6C KO","Tg-WT"))
merged.obj2$genotype = factor(merged.obj2$genotype, levels = c("Tg-6C KO","Tg-WT"))
merged.obj2 = NormalizeData(merged.obj2)
merged.obj2 = FindVariableFeatures(merged.obj2)
merged.obj2 = ScaleData(merged.obj2)
merged.obj2 = RunPCA(merged.obj2)
ElbowPlot(merged.obj2, ndims = 50)
merged.obj2 = RunUMAP(merged.obj2, dims = 1:30)
merged.obj2 = RunHarmony(merged.obj2, group.by.vars = "orig.ident")
merged.obj2 = RunUMAP(merged.obj2, dims = 1:30, reduction = "harmony")
merged.obj2 = FindNeighbors(merged.obj2, dims = 1:2, reduction = "umap")
merged.obj2 = FindClusters(merged.obj2, resolution = 0.1)

#cell_types = cell_types[1:5]
Idents(merged.obj2) = merged.obj2$cell_type
genotypes = unique(merged.obj2$genotype)
cell_type = c()
samples = c()
value =c()
for (s in 1:100) {
  for (i in 1:length(cell_types)) {
    tmp = subset(merged.obj2, idents = cell_types[i])
    mingenotype = min(table(tmp$genotype))
    cell.list = c(sample(colnames(tmp)[tmp$genotype==genotypes[1]],mingenotype),
                  sample(colnames(tmp)[tmp$genotype==genotypes[2]],mingenotype))
    tmp = FindNeighbors(tmp, reduction = "harmony", dims = 1:10, return.neighbor = T)
    NN = tmp@neighbors$RNA.nn@nn.idx[,c(1,2)]
    df = cbind.data.frame("Cell_geno" = tmp$genotype[NN[,1]],"Neighbor_geno"= tmp$genotype[NN[,2]], "Sample" = tmp$orig.ident[NN[,1]])
    df$same = ifelse(df$Cell_geno==df$Neighbor_geno,1,0)
    for (x in 1:length(unique(df$Sample))) {
      tmp = df[df$Sample==unique(df$Sample)[x],]
      cell_type = c(cell_type,cell_types[i])
      samples = c(samples,unique(df$Sample)[x])
      value = c(value,sum(tmp$same)/nrow(tmp))
    }
  }}
results = cbind.data.frame(cell_type, samples, value)



cell_type = c()
samples = c()
value =c()
for (s in 1:100) {
  for (i in 1:length(unique(cell_types))) {
    tmp = subset(merged.obj2, idents = cell_types[i])
    mingenotype = min(table(tmp$genotype))
    cell.list = c(sample(colnames(tmp)[tmp$genotype==genotypes[1]],mingenotype),
                  sample(colnames(tmp)[tmp$genotype==genotypes[2]],mingenotype))
    tmp = subset(tmp, cells = cell.list)
    tmp = FindNeighbors(tmp, reduction = "harmony", dims = 1:10, return.neighbor = T)
    NN = tmp@neighbors$RNA.nn@nn.idx[,c(1,2)]
    df = cbind.data.frame("Cell_geno" = sample(tmp$genotype[NN[,1]]),"Neighbor_geno"= tmp$genotype[NN[,2]], "Sample" = tmp$orig.ident[NN[,1]])
    df$same = ifelse(df$Cell_geno==df$Neighbor_geno,1,0)
    for (x in 1:length(unique(df$Sample))) {
      tmp = df[df$Sample==unique(df$Sample)[x],]
      cell_type = c(cell_type,cell_types[i])
      samples = c(samples,unique(df$Sample)[x])
      value = c(value,mean(tmp$same))
    }
  }
}
results2 = cbind.data.frame(cell_type, samples, value)

library(ggbeeswarm)
p = ggbarplot(results, x = "cell_type", y = "value", fill = "lightblue",
              title = "Segregation by genotype\n within each cell type", 
              ylab = "Proportion of nearest neighbors in same genotype",
              xlab = F, add = c("mean_se"))+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_y_continuous(limits=c(0,max(results$value)), expand = c(0, 0),breaks = scales::pretty_breaks(n = 10))+
  geom_hline(yintercept = mean(results2$value),linetype = "dashed")


# p = ggbarplot(results, x = "cell_type", y = "value", fill = "lightblue",
#               title = "Segregation by genotype\n within each cell type", 
#               ylab = "Proportion of nearest neighbors in same genotype",
#               xlab = F, add = c("mean_se","jitter"))+ geom_jitter(alpha = 0.2, width = 0.3)+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
#   scale_y_continuous(limits=c(min(results$value)-0.1,max(results$value)+0.05), oob=rescale_none)+
#   geom_hline(yintercept = mean(results2$value),linetype = "dashed")+geom_hline(yintercept = 0.25,linetype = "solid")
ggsave(plot = p, "plots/Figure3B_4.pdf", width = 4, height = 7, units = "in")



# 3C
Ms4as = c("Ms4a1","Ms4a2","Ms4a3","Ms4a4a","Ms4a4b","Ms4a4c","Ms4a6b","Ms4a6c","Ms4a6d","Ms4a7","Ms4a8a","Ms4a15")
FeaturePlot_scCustom(merged.obj, features = Ms4as, num_columns = 6, )+NoAxes()
ggsave(plot = p, "plots/Figure3C.pdf", width = 12, height = 4, units = "in")

plots = list()
for (i in 1:length(Ms4as)) {
  plots[[i]] = FeaturePlot_scCustom(merged.obj, paste0(Ms4as[i]))+NoLegend()+NoAxes()+labs(title = Ms4as[i])
}
wrap_plots(plots, ncol = 6, byrow = T)
ggsave("plots/Figure3C.pdf", width = 12, height = 5, units = "in")

plots = list()
for (i in 1:length(Ms4as)) {
  plots[[i]] = FeaturePlot_scCustom(merged.obj, paste0(Ms4as[i]))+NoLegend()+NoAxes()+labs(title = "")
}
wrap_plots(plots, ncol = 6, byrow = T)
ggsave("plots/Figure3C.jpeg", width = 12, height = 5, units = "in")

FeaturePlot_scCustom(merged.obj, "Ms4a6c")+NoAxes()
ggsave("plots/Figure3C_scale.pdf", width = 5, height = 5, units = "in")


# 3D
# Identify markers of each cluster on hpcc
files = list.files("data/", pattern = "MG_marker")
markers = list()
for (i in 1:length(files)) {
  tmp = read.csv(paste0("data/",files[i]))
  tmp$cluster = gsub("MG_marker_|.csv","",files[i])
  markers[[i]] = tmp
}

markers = do.call(rbind.data.frame,markers)
rownames(markers) = NULL
colnames(markers)[1] = "gene"

markers = markers[grep("mt-|Rpl|Rps",markers$gene, invert = T),]
markers %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC) -> top20
top20 = top20[!duplicated(top20$gene),]
top20$cluster = factor(top20$cluster, levels = c("0","1","6","2","3","4","5"))
top20 = top20[order(top20$cluster),]


ann = top20$cluster
names(ann) = top20$gene
cols = hue_pal()(length(unique(top20$cluster)))
names(cols) = unique(top20$cluster)
cols <- split(unname(cols),names(cols))


ha =  rowAnnotation(cluster = ann, col = list("0" = "#F8766D",
                                              "1" = "#C49A00",
                                              "2" = "#53B400",
                                              "3" = "#00C094",
                                              "4" = "#00B6EB",
                                              "5" = "#A58AFF",
                                              "6" = "#FB61D7"))
#ha = rowAnnotation(foo = anno_block(gp = gpar(fill = 2:4),
#                                    labels = c("group1", "group2", "group3"), 
#                                    labels_gp = gpar(col = "white", fontsize = 10)))

Idents(MG) = MG$genotype
metadf = MG@meta.data
mat<- MG[["RNA"]]@data[top20$gene,] %>% as.matrix()

## scale the rows
mat<- t(scale(t(mat)))
mat[mat>2] = 2
mat[mat<(-2)] = -(2)
mat_1 = mat[,rownames(metadf)[metadf$genotype=="Ntg-WT"]]
mat_2 = mat[,rownames(metadf)[metadf$genotype=="Ntg-6C KO"]]
mat_3 = mat[,rownames(metadf)[metadf$genotype=="Tg-WT"]]
mat_4 = mat[,rownames(metadf)[metadf$genotype=="Tg-6C KO"]]
col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

ht3 = Heatmap(mat_3, column_title = "Tg-WT", name = "Expression",
              show_heatmap_legend = T,
              cluster_columns = T,
              show_column_dend = FALSE,
              cluster_rows = F,
              show_row_dend = FALSE,
              row_names_gp = gpar(fontsize = 4),
              width = unit(4, "cm"),
              show_column_names = FALSE,
              use_raster = TRUE,
              col = col_fun,
              raster_quality = 4)
ht4 = Heatmap(mat_4, column_title = "Tg-6C KO",
              show_heatmap_legend = F,
              cluster_columns = T,
              show_column_dend = FALSE,
              cluster_rows = F,
              show_row_dend = FALSE,
              row_names_gp = gpar(fontsize = 4),
              width = unit(4, "cm"),
              show_column_names = FALSE,
              use_raster = TRUE,
              col = col_fun,
              raster_quality = 4)
ht1 = Heatmap(mat_1, column_title = "Ntg-WT",
              cluster_columns = T,
              show_heatmap_legend = F,
              show_column_dend = FALSE,
              cluster_rows = F,
              show_row_dend = FALSE,
              row_names_gp = gpar(fontsize = 4),
              width = unit(4, "cm"),
              show_column_names = FALSE,
              use_raster = TRUE,
              col = col_fun,
              raster_quality = 4)
ht2 = Heatmap(mat_2, column_title = "Ntg-6C KO",
              cluster_columns = T,
              show_heatmap_legend = F,
              show_column_dend = FALSE,
              cluster_rows = F,
              show_row_dend = FALSE,
              row_names_gp = gpar(fontsize = 4),
              width = unit(4, "cm"),
              show_column_names = FALSE,
              use_raster = TRUE,
              col = col_fun,
              raster_quality = 4)

ha = Heatmap(top20$cluster, col = list("0" = "#F8766D",
                                                          "1" = "#C49A00",
                                                          "2" = "#53B400",
                                                          "3" = "#00C094",
                                                          "4" = "#00B6EB",
                                                          "5" = "#A58AFF",
                                                          "6" = "#FB61D7"), 
             width = unit(0.5, "cm"), name = "cluster")



pdf("plots/Figure3E.pdf", width =10, height = 10)
ha+ht1+ht2+ht3+ht4
dev.off()

# 3D
DimPlot(MG, group.by = "seurat_clusters", label = T)+NoLegend()+ggtitle("")
ggsave("plots/Figure3D_1.pdf", width = 5, height = 5, units = "in")
DimPlot(MG, group.by = "genotype", split.by = "genotype", ncol = 2)+NoLegend()
ggsave("plots/Figure3D_2.pdf", width = 5, height = 5, units = "in")

p = DimPlot(MG, group.by = "seurat_clusters", label = T)
data = p$data
data$circle = NA
data[data$seurat_clusters%in%c(1,6),]$circle = "Y"

p = ggplot(data, aes(umap_1, umap_2, color = seurat_clusters))+geom_point(size =0.3)+theme_classic()+ggtitle("")
p = LabelClusters(plot = p, id = "seurat_clusters", color = "black")+geom_mark_hull(data = data[data$circle=="Y",], aes(group = circle), expand=0.01, colour="black", size = 1, linetype = "dashed")+NoLegend()
ggsave("plots/Figure3D_1_2.pdf", width = 5, height = 5, units = "in")



p2 = DimPlot(MG, group.by = "genotype", split.by = "genotype", ncol = 2)
data = p2$data
data2 = cbind.data.frame("xx" = ggplot_build(p)$data[[3]]$x, "yy" = ggplot_build(p)$data[[3]]$y)
ggplot(data, aes(umap_1, umap_2, color = genotype))+geom_point(size =0.3)+theme_classic()+facet_wrap(~genotype, ncol =2)+
  geom_mark_hull(data = data2,aes(xx, yy), expand=0.01, colour="black", size = 1, linetype = "dashed")+
  theme(strip.background = element_blank())+NoLegend()+ggtitle("genotype") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))+
  theme(strip.text = element_text(size = 12, face = "bold"))
ggsave("plots/Figure3E_2_2.pdf", width = 5, height = 5, units = "in")


# 3F
p2 = DimPlot(MG, group.by = "classfication", split.by = "genotype", ncol = 2)
data = p2$data
data2 = cbind.data.frame("xx" = ggplot_build(p)$data[[3]]$x, "yy" = ggplot_build(p)$data[[3]]$y)
ggplot(data, aes(umap_1, umap_2, color = classfication))+geom_point(size =0.3)+theme_classic()+facet_wrap(~genotype, ncol =2)+
  geom_mark_hull(data = data2,aes(xx, yy), expand=0.01, colour="black", size = 1, linetype = "dashed")+
  theme(strip.background = element_blank())+NoLegend()+ggtitle("genotype") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))+
  theme(strip.text = element_text(size = 12, face = "bold"))
ggsave("plots/Figure3F.pdf", width = 6, height = 5, units = "in")


# 3G
metadf = MG@meta.data
metadf$classfication = factor(metadf$classfication, levels = c("DAM_depleted", "DAM1", "DAM2"))
df = metadf %>% dplyr::group_by(orig.ident) %>% count(classfication, genotype)%>%
  mutate(freq = n / sum(n))
stat.test = df %>% group_by(classfication) %>%
  t_test(freq ~ genotype) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj")
stat.test = stat.test[stat.test$p.adj<0.05,]

stat.test <- stat.test %>%
  add_xy_position(x = "classfication")

ggplot(df, aes(classfication, freq)) +
  geom_boxplot(aes(fill = genotype)) + 
  theme_bw() + 
  stat_pvalue_manual(stat.test,   label = "p.adj.signif", tip.length = 0.005, y.position = c(1.05,1.1,1.15,1.05,1.1,1.15))+
  scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1))+theme_classic2()+ylab("composition (%)")
ggsave("plots/Figure3G.pdf", width = 5, height = 5, units = "in")

# 3H
Stacked_VlnPlot(MG, features = c("Clec7a","Axl","Cst7","Ifi27l2a","Cd74"), group.by = "genotype", pt.size = 0.01)
ggsave("plots/Figure3H.pdf", width = 6, height = 6, units = "in")

# 3I
## DAMs, Mapt vs 6C;Mapt
Idents(MG) = MG$classfication
DAMs = subset(MG, idents = "DAM_depleted", invert = T)
Idents(DAMs) = DAMs$genotype
DEGs1 = FindMarkers(DAMs, ident.1 = "Tg-6C KO",ident.2 = "Tg-WT",test.use = "MAST", latent.vars = "orig.ident")
### 6CMT enrich
geneset_6CMT = rownames(DEGs1)[DEGs1$p_val_adj<0.05&DEGs1$avg_log2FC>0]
GO_6CMT = topGOterms(fg.genes = geneset_6CMT, bg.genes = rownames(MG), organism = "Mouse")
goEnrichment_top = GO_6CMT$res.result
goEnrichment_top$Term = paste(toupper(substr(goEnrichment_top$Term, 1, 1)), substr(goEnrichment_top$Term, 2, nchar(as.character(goEnrichment_top$Term))), sep="")
goEnrichment_top = goEnrichment_top[order(goEnrichment_top$GeneRatio,decreasing = F),]
goEnrichment_top$Term = factor(goEnrichment_top$Term, levels = goEnrichment_top$Term)
ggplot(goEnrichment_top[1:20,], aes(x = Term, y = GeneRatio,
                         size = Significant))+
  geom_point(aes(colour = -log10(padj)))+
  scale_color_gradient(low = "red", high = "blue", name = "Adjusted p-value\n(-log10)")+
  theme_bw()+ggtitle("Biological pathway")+ylab("GeneRatio")+
  xlab("")+
  labs(size="Counts")+scale_size( range = c(3,6))+
  coord_flip()+theme_bw()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave("plots/6CMT_vs_MAPT_GO_in_DAM.pdf", width = 10, height = 5)

# GO analysis between mapt,6c and 6c in all microglia
Idents(MG) = MG$genotype
DEGs = FindMarkers(MG, ident.1 = "Tg-6C KO",ident.2 = "Tg-WT",test.use = "MAST", latent.vars = "orig.ident")

geneset_1 = rownames(DEGs)[DEGs$p_val_adj<0.05&DEGs$avg_log2FC<0]
GO_cls1 = topGOterms(fg.genes = geneset_1, bg.genes = rownames(MG), organism = "Mouse")
goEnrichment_top = GO_cls1$res.result
goEnrichment_top = goEnrichment_top[!is.na(goEnrichment_top$Annotated),]
goEnrichment_top$Term = paste(toupper(substr(goEnrichment_top$Term, 1, 1)), substr(goEnrichment_top$Term, 2, nchar(as.character(goEnrichment_top$Term))), sep="")
goEnrichment_top = goEnrichment_top[order(goEnrichment_top$GeneRatio,decreasing = F),]
goEnrichment_top$Term = factor(goEnrichment_top$Term, levels = goEnrichment_top$Term)

ggplot(goEnrichment_top, aes(x = Term, y = GeneRatio,
                                    size = Significant))+
  geom_point(aes(colour = -log10(padj)))+
  scale_color_gradient(low = "red", high = "blue", name = "Adjusted p-value\n(-log10)")+
  theme_bw()+ggtitle("Biological pathway")+ylab("GeneRatio")+
  xlab("")+
  labs(size="Counts")+scale_size( range = c(3,6))+
  coord_flip()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave("plots/Mapt_vs_Mapt6c_GO_in_MG.pdf", width = 10, height = 5)

geneset_1 = rownames(DEGs)[DEGs$p_val_adj<0.05&DEGs$avg_log2FC>0]
GO_cls1 = topGOterms(fg.genes = geneset_1, bg.genes = rownames(MG), organism = "Mouse")
goEnrichment_top = GO_cls1$res.result
goEnrichment_top = goEnrichment_top[!is.na(goEnrichment_top$Annotated),]
goEnrichment_top$Term = paste(toupper(substr(goEnrichment_top$Term, 1, 1)), substr(goEnrichment_top$Term, 2, nchar(as.character(goEnrichment_top$Term))), sep="")
goEnrichment_top = goEnrichment_top[order(goEnrichment_top$GeneRatio,decreasing = F),]
goEnrichment_top$Term = factor(goEnrichment_top$Term, levels = goEnrichment_top$Term)

ggplot(goEnrichment_top, aes(x = Term, y = GeneRatio,
                             size = Significant))+
  geom_point(aes(colour = -log10(padj)))+
  scale_color_gradient(low = "red", high = "blue", name = "Adjusted p-value\n(-log10)")+
  theme_bw()+ggtitle("Biological pathway")+ylab("GeneRatio")+
  xlab("")+
  labs(size="Counts")+scale_size( range = c(3,6))+
  coord_flip()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave("plots/Mapt6c_vs_Mapt_GO_in_MG.pdf", width = 10, height = 5)

Idents(MG) = MG$genotype
DEGs = FindMarkers(MG, ident.1 = "Tg-WT",ident.2 = "Tg-6C KO",test.use = "MAST", latent.vars = "orig.ident")
DEGs$gene = rownames(DEGs)
write.xlsx(DEGs, "data/Mapt_vs_MT6C_MG_DEGs.xlsx")


# Generate volcano plots comparing Mapt vs 6C;Mapt and WT vs 6C KO using genes from Rangajaru et al. Magenta, Yellow, and Blue and phagocytosis gene lists
gene_table = read.xlsx("data/Rangaraju et al. Magenta_Yellow_Blue and Phagocytosis.xlsx")

Idents(MG) = MG$genotype
Mapt_6CMT = FindMarkers(MG, ident.1 = "Tg-WT",ident.2 = "Tg-6C KO",test.use = "MAST", latent.vars = "orig.ident")
WT_6CKO = FindMarkers(MG, ident.1 = "Ntg-WT",ident.2 = "Ntg-6C KO",test.use = "MAST", latent.vars = "orig.ident")

saveRDS(Mapt_6CMT, "data/DEG_Mapt_6CMT_in_MG.rds")
saveRDS(WT_6CKO, "data/DEG_WT_6CKOin_MG.rds")

Mapt_6CMT = readRDS("data/DEG_Mapt_6CMT_in_MG.rds")


## MAPT vs 6CMT
EnhancedVolcano(Mapt_6CMT[rownames(Mapt_6CMT)%in%gene_table[,1],],
                lab = rownames(Mapt_6CMT[rownames(Mapt_6CMT)%in%gene_table[,1],]),
                x = 'avg_log2FC',
                y = 'p_val_adj', 
                title = "Tg-WT vs Tg-6C KO",
                subtitle = "",
                )
ggsave("plots/Vlo_Mapt_vs_6CMT_Magenta_genes.pdf", width = 7, height = 7)

EnhancedVolcano(Mapt_6CMT[rownames(Mapt_6CMT)%in%gene_table[,2],],
                lab = rownames(Mapt_6CMT[rownames(Mapt_6CMT)%in%gene_table[,2],]),
                x = 'avg_log2FC',
                y = 'p_val_adj', 
                title = "Tg-WT vs Tg-6C KO",
                subtitle = "",
)
ggsave("plots/Vlo_Mapt_vs_6CMT_Yellow_genes.pdf", width = 7, height = 7)

EnhancedVolcano(Mapt_6CMT[rownames(Mapt_6CMT)%in%gene_table[,3],],
                lab = rownames(Mapt_6CMT[rownames(Mapt_6CMT)%in%gene_table[,3],]),
                x = 'avg_log2FC',
                y = 'p_val_adj', 
                title = "Tg-WT vs Tg-6C KO",
                subtitle = "",
)
ggsave("plots/Vlo_Mapt_vs_6CMT_Blue_genes.pdf", width = 7, height = 7)

EnhancedVolcano(Mapt_6CMT[rownames(Mapt_6CMT)%in%gene_table[,4],],
                lab = rownames(Mapt_6CMT[rownames(Mapt_6CMT)%in%gene_table[,4],]),
                x = 'avg_log2FC',
                y = 'p_val_adj', 
                title = "Tg-WT vs Tg-6C KO",
                subtitle = "",
)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave("plots/Vlo_Mapt_vs_6CMT_Phagocytosis_genes.pdf", width = 7, height = 7)

## WT vs 6CKO
EnhancedVolcano(WT_6CKO[rownames(WT_6CKO)%in%gene_table[,1],],
                lab = rownames(WT_6CKO[rownames(WT_6CKO)%in%gene_table[,1],]),
                x = 'avg_log2FC',
                y = 'p_val_adj', 
                title = "Ntg-WT vs Ntg-6C KO",
                subtitle = "",
)
ggsave("plots/Vlo_WT_vs_6CKO_Magenta_genes.pdf", width = 7, height = 7)

EnhancedVolcano(WT_6CKO[rownames(WT_6CKO)%in%gene_table[,2],],
                lab = rownames(WT_6CKO[rownames(WT_6CKO)%in%gene_table[,2],]),
                x = 'avg_log2FC',
                y = 'p_val_adj', 
                title = "Ntg-WT vs Ntg-6C KO",
                subtitle = "",
)
ggsave("plots/Vlo_WT_vs_6CKO_Yellow_genes.pdf", width = 7, height = 7)

EnhancedVolcano(WT_6CKO[rownames(WT_6CKO)%in%gene_table[,3],],
                lab = rownames(WT_6CKO[rownames(WT_6CKO)%in%gene_table[,3],]),
                x = 'avg_log2FC',
                y = 'p_val_adj', 
                title = "Ntg-WT vs Ntg-6C KO",
                subtitle = "",
)
ggsave("plots/Vlo_WT_vs_6CKO_Blue_genes.pdf", width = 7, height = 7)

EnhancedVolcano(WT_6CKO[rownames(WT_6CKO)%in%gene_table[,4],],
                lab = rownames(WT_6CKO[rownames(WT_6CKO)%in%gene_table[,4],]),
                x = 'avg_log2FC',
                y = 'p_val_adj', 
                title = "Ntg-WT vs Ntg-6C KO",
                subtitle = "",
)
ggsave("plots/Vlo_WT_vs_6CKO_Phagocytosis_genes.pdf", width = 7, height = 7)

# UMAP plots of Clec7a, ApoE, Trem2, Itgax, Cst7, P2ry12, and Tmem119

genes = c("Clec7a","Apoe","Trem2","Itgax","Cst7","P2ry12","Tmem119")
FeaturePlot_scCustom(MG, features = genes, num_columns = 6, )+NoAxes()

plots = list()
for (i in 1:length(genes)) {
  plots[[i]] = FeaturePlot_scCustom(MG, paste0(genes[i]))+NoLegend()+NoAxes()+labs(title = genes[i])
}
wrap_plots(plots, ncol = 4, byrow = T)
ggsave("plots/Umap_DAM_genes.pdf", width = 8, height = 5, units = "in")

plots = list()
for (i in 1:length(genes)) {
  plots[[i]] = FeaturePlot_scCustom(MG, paste0(genes[i]))+NoLegend()+NoAxes()+labs(title = "")
}
wrap_plots(plots, ncol = 4, byrow = T)
ggsave("plots/Umap_DAM_genes.jpeg", width = 8, height = 5, units = "in")

# cluster composition

# Proportion / cell number composition per cluster
ggData = data.frame(prop.table(table(MG$genotype, MG$seurat_clusters), margin = 2))
colnames(ggData) = c("genotype", "cluster", "value")

ggplot(ggData, aes(cluster, value, fill = genotype)) +
  geom_col() + xlab("Cluster") + ylab("Proportion of Cells (%)") +
  scale_fill_manual(values = c("#56B4E9", "#009E73",
                               "#F0E442", "#0072B2")) + 
  coord_flip() + theme_classic2()
ggsave("plots/Cluster_proportion.pdf", width = 6, height = 4, units = "in")


ggData = data.frame(prop.table(table(MG$seurat_clusters, MG$genotype), margin = 2))
colnames(ggData) = c("Cluster", "Genotype", "Value")
ggplot(ggData, aes(Genotype, Value, fill = Cluster)) +
  geom_col() + xlab("Genotype") + ylab("Proportion of Cells (%)") + 
  coord_flip() + theme_classic2()
ggsave("plots/Extended_Figure2C.pdf", width = 6, height = 4, units = "in")


# volcano plots for the gene lists attached for Mapt vs 6C;Mapt?
genelist = read.xlsx("data/Gene_lists_TE4T2KO_TFEB_sensome_lysosome.xlsx")


EnhancedVolcano(Mapt_6CMT[rownames(Mapt_6CMT)%in%genelist[,1],],
                lab = rownames(Mapt_6CMT[rownames(Mapt_6CMT)%in%genelist[,1],]),
                x = 'avg_log2FC',
                y = 'p_val_adj', 
                title = "Ntg-WT vs Ntg-6C KO",
                subtitle = "",
)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave("plots/Vlo_Mapt_vs_6CMT_Downregulated_in_TE4-T2KO_genes.pdf", width = 7, height = 7)


EnhancedVolcano(Mapt_6CMT[rownames(Mapt_6CMT)%in%genelist[,2],],
                lab = rownames(Mapt_6CMT[rownames(Mapt_6CMT)%in%genelist[,2],]),
                x = 'avg_log2FC',
                y = 'p_val_adj', 
                title = "Ntg-WT vs Ntg-6C KO",
                subtitle = "",
)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave("plots/Vlo_Mapt_vs_6CMT_Upregulated_in_TE4-T2KO_genes.pdf", width = 7, height = 7)

EnhancedVolcano(Mapt_6CMT[rownames(Mapt_6CMT)%in%genelist[,3],],
                lab = rownames(Mapt_6CMT[rownames(Mapt_6CMT)%in%genelist[,3],]),
                x = 'avg_log2FC',
                y = 'p_val_adj', 
                title = "Ntg-WT vs Ntg-6C KO",
                subtitle = "",
)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave("plots/Vlo_Mapt_vs_6CMT_Downregulated_by_TFEB_genes.pdf", width = 7, height = 7)

EnhancedVolcano(Mapt_6CMT[rownames(Mapt_6CMT)%in%genelist[,4],],
                lab = rownames(Mapt_6CMT[rownames(Mapt_6CMT)%in%genelist[,4],]),
                x = 'avg_log2FC',
                y = 'p_val_adj', 
                title = "Ntg-WT vs Ntg-6C KO",
                subtitle = "",
)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave("plots/Vlo_Mapt_vs_6CMT_Upregulated_by_TFEB_genes.pdf", width = 7, height = 7)

EnhancedVolcano(Mapt_6CMT[rownames(Mapt_6CMT)%in%genelist[,5],],
                lab = rownames(Mapt_6CMT[rownames(Mapt_6CMT)%in%genelist[,5],]),
                x = 'avg_log2FC',
                y = 'p_val_adj', 
                title = "Ntg-WT vs Ntg-6C KO",
                subtitle = "",
)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave("plots/Vlo_Mapt_vs_6CMT_Sensome_genes.pdf", width = 7, height = 7)

EnhancedVolcano(Mapt_6CMT[rownames(Mapt_6CMT)%in%genelist[,6],],
                lab = rownames(Mapt_6CMT[rownames(Mapt_6CMT)%in%genelist[,6],]),
                x = 'avg_log2FC',
                y = 'p_val_adj', 
                title = "Ntg-WT vs Ntg-6C KO",
                subtitle = "",
)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave("plots/Vlo_Mapt_vs_6CMT_Lysosomal_genes.pdf", width = 7, height = 7)

DotPlot_scCustom(MG, features = rownames(MG)[rownames(MG)%in%genelist[,6]], x_lab_rotate = T)+
  scale_colour_gradient2(low = "blue", mid = "white", high = "Red") +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5)

ggsave("plots/Dot_Mapt_vs_6CMT_Lysosomal_genes.pdf", width = 7, height = 7)

#  calculate and plot the % cell type composition by genotype. for neurons, can you break them down by inhibitory and inhibitory populations.
Neurons = subset(merged.obj, ident = "Neurons")
Neurons = FindVariableFeatures(Neurons)
Neurons = ScaleData(Neurons)
Neurons = RunPCA(Neurons)
Neurons = RunHarmony(Neurons, group.by.vars = "orig.ident")
Neurons = RunUMAP(Neurons, reduction = "harmony", dims = 1:10)
Neurons = FindNeighbors(Neurons, reduction = "harmony")
Neurons = FindClusters(Neurons, resolution = 1)
DimPlot_scCustom(Neurons) 

FeaturePlot_scCustom(Neurons, "Gad1") #In #11, 7, 6, 17
FeaturePlot_scCustom(Neurons, "Cck")
FeaturePlot_scCustom(Neurons, "Sst")
FeaturePlot_scCustom(Neurons, "Vip")
FeaturePlot_scCustom(Neurons, "Calb2")
FeaturePlot_scCustom(Neurons, "Dlx5")
FeaturePlot_scCustom(Neurons, "Lhx6")
FeaturePlot_scCustom(Neurons, "Foxg1")

FeaturePlot_scCustom(Neurons, "Neurod6") #EX
FeaturePlot_scCustom(Neurons, "Neurod2")
FeaturePlot_scCustom(Neurons, "Slc17a7") 
FeaturePlot_scCustom(Neurons, "Tbr1")
FeaturePlot_scCustom(Neurons, "Reln")

DotPlot_scCustom(Neurons, c("Neurod6","Neurod2","Slc17a7", "Gad1", "Gad2"))

INN = WhichCells(Neurons, idents = c(6,7,11,17))

merged.obj$cell_type2 = merged.obj$cell_type
merged.obj$cell_type2[merged.obj$cell_type2=="Neurons"] = "Neurons (Ex)"
merged.obj$cell_type2[colnames(merged.obj)%in%INN] = "Neurons (In)"
DimPlot_scCustom(merged.obj, group.by = "cell_type2")

Idents(merged.obj) = merged.obj$genotype
x = Cluster_Stats_All_Samples(merged.obj, group_by_var = "cell_type2")
options(scipen=999)
x2 = melt(x[,c(1,grep("_%",colnames(x)))])
x2 = x2[!x2$Cluster=="Total",]
x2$variable = gsub("_%","",x2$variable)
library(RColorBrewer)
n = length(unique(x2$Cluster))
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col=sample(col_vector, n)
ggplot(x2, aes(fill=Cluster, y=value, x=variable)) + 
  geom_bar(position="stack", stat="identity")+coord_flip()+
  scale_fill_manual(values= col)+xlab("")+ylab("%")
ggsave("plots/Cell_composition_in_each_genotype.pdf", width = 10, height = 6)

Idents(merged.obj) = merged.obj$cell_type2
x = Cluster_Stats_All_Samples(merged.obj, group_by_var = "orig.ident")
x2 = melt(x[,c(1,grep("_%",colnames(x)))])
x2 = x2[!x2$Cluster=="Total",]
x2$variable = gsub("_%","",x2$variable)
x2$genotype = "Tg-WT"
x2$genotype[grep("6CKO",x2$variable)] = "Ntg-6C KO"
x2$genotype[grep("6CMT",x2$variable)] = "Tg-6C KO"
x2$genotype[grep("WT",x2$variable)] = "Ntg-WT"
x2$genotype = factor(x2$genotype, levels = c("Ntg-WT","Ntg-6C KO","Tg-WT","Tg-6C KO"))

ggplot(x2, aes(fill=genotype, y=value, x=Cluster)) + 
  geom_boxplot()+coord_flip()+
  scale_fill_manual(values= col)+xlab("")+ylab("%")

my_comparisons <- list( c("Tg-WT", "Tg-6C KO"), c("Ntg-WT", "Tg-WT") )
ggboxplot(x2, y = "value", fill = "genotype", facet.by = "Cluster", nrow = 3)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  stat_compare_means(aes(group = genotype), label = "p.format")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  stat_compare_means()+ylab("%")
ggsave("plots/Cell_composition_in_each_genotype_by_sample.pdf", width = 10, height = 12)

phg = read.xlsx("data/Phagocytosis_genes.xlsx")
phg = phg$Phagocytosis
phg = phg[phg%in%rownames(MG)]
phg = phg[!phg%in%rownames(MG)[rowSums(MG@assays$RNA@data)==0]]

p = Clustered_DotPlot(MG, features = phg)
x = ComplexHeatmap::row_order(p[[2]])
phg = phg[x]

p = DotPlot_scCustom(MG, features = phg, x_lab_rotate = T)+
  scale_colour_gradient2(low = "blue", mid = "white", high = "Red") +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5)
ggsave("plots/Dot_plot_for_phagocytosis_gene_in_MG_by_genotype.pdf", width = 30, height = 4)

df = p$data
df = reshape(df[,c("features.plot","id","pct.exp")], direction = "wide", idvar = "features.plot", timevar = "id")
 which.max(df)
        
df$max <- apply(df[,grep("pct.exp",colnames(df))], 1, max, na.rm=TRUE)
phg2 = phg[df$max>25]
DotPlot_scCustom(MG, features = phg2, x_lab_rotate = T)+
  scale_colour_gradient2(low = "blue", mid = "white", high = "Red") +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5)
ggsave("plots/Dot_plot_for_phagocytosis_gene_in_MG_by_genotype_expressed_at_least_25_percent.pdf", width = 20, height = 4)

Idents(MG) = MG$genotype
pmarkers = FindAllMarkers(MG, test.use = "MAST", latent.vars = "orig.ident", features = phg2)
pmarkers =pmarkers[pmarkers$p_val_adj<0.05,]
pmarkers = pmarkers[order(pmarkers$avg_log2FC, decreasing = T),]

DotPlot_scCustom(MG, features = pmarkers$gene[1:15], x_lab_rotate = T)+
  scale_colour_gradient2(low = "blue", mid = "white", high = "Red") +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5)
ggsave("plots/Dot_plot_for_top_DEG_phagocytosis_gene_in_MG_by_genotype_expressed_at_least_25_percent.pdf", width = 7, height = 4)

# violin plot to compare ApoE expression levels between genotypes in neurons, astrocytes, and microglia.
p1 = VlnPlot_scCustom(MG, features = "Apoe", group.by = "genotype", pt.size = 0.01, colors_use = scales::hue_pal()(4))+NoLegend()+xlab("")+ggtitle("Microglia")

Idents(merged.obj) = merged.obj$cell_type
Neurons = subset(merged.obj, idents = "Neurons")
Astrocytes = subset(merged.obj, idents = "Astrocytes")

p2 = VlnPlot_scCustom(Neurons, features = "Apoe", group.by = "genotype", pt.size = 0.01, colors_use = scales::hue_pal()(4))+NoLegend()+xlab("")+ggtitle("Neurons")
p3 = VlnPlot_scCustom(Astrocytes, features = "Apoe", group.by = "genotype", pt.size = 0.01, colors_use = scales::hue_pal()(4))+NoLegend()+xlab("")+ggtitle("Astrocytes")

wrap_plots(list(p1,p2,p3), nrow = 1)
ggsave("plots/Vln_plot_Apoe_in_MG_Neu_Ast.pdf", width = 7, height = 4)

# total # of cell types, # of cells and # of microglia included in the analysis

Idents(merged.obj) = merged.obj$cell_type
stats <- Cluster_Stats_All_Samples(seurat_object = merged.obj)
write.xlsx(stats, "data/Cell_numbers.xlsx")

# heatmap of ms4a expression in these cell types (myeloid including microglia and macrophages).
celltype = c("Astrocytes", "Oligodendrocytes", "Neurons", "Endothelial cells", "Marcophages", "Microglia", "Vascular smooth muscle cells","Pericytes")
Idents(merged.obj) = merged.obj$cell_type
sub.obj = subset(merged.obj , idents = celltype)
sub.obj$cell_type[sub.obj$cell_type%in%c("Marcophages", "Microglia")] = "Myeloid"
sub.obj$cell_type[sub.obj$cell_type%in%c("Endothelial cells","Vascular smooth muscle cells","Pericytes")] = "Endothelium"
Idents(sub.obj) = sub.obj$cell_type
sub.obj_ave = AverageExpression(sub.obj, return.seurat = T)
EM = sub.obj_ave@assays$RNA@data
genes = c("Ms4a1","Ms4a2","Ms4a3","Ms4a4a","Ms4a4b","Ms4a4c","Ms4a6b","Ms4a6c","Ms4a6d","Ms4a7","Ms4a8a","Ms4a15")

pdf("plots/cell_type_Ms4a_heatmap.pdf", width = 6, height = 3)
ComplexHeatmap::Heatmap(t(EM[genes,]),col = RColorBrewer::brewer.pal(9, "Reds"),
                        name = "Mean\nExpression", border = T, cluster_rows = F, cluster_columns = F, row_names_side = "left", column_names_rot = 45,
                        )
dev.off()

# Redo 3B with only Tg and Tg-6CKO
# 3B
cell_types = names(table(merged.obj$cell_type)[table(merged.obj$cell_type)>1000])
#cell_types = cell_types[1:5]
Idents(merged.obj) = merged.obj$genotype
genotypes = c("Tg-WT","Tg-6C KO")
merged.obj = subset(merged.obj, idents = genotypes)
merged.obj$genotype = factor(merged.obj$genotype , levels = genotypes)
cell_type = c()
samples = c()
value =c()
Idents(merged.obj) = merged.obj$cell_type
for (s in 1:100) {
  for (i in 1:length(cell_types)) {
    tmp = subset(merged.obj, idents = cell_types[i])
    mingenotype = min(table(tmp$genotype))
    cell.list = c(sample(colnames(tmp)[tmp$genotype==genotypes[1]],mingenotype),
                  sample(colnames(tmp)[tmp$genotype==genotypes[2]],mingenotype))
    tmp = FindNeighbors(tmp, reduction = "harmony", dims = 1:10, return.neighbor = T)
    NN = tmp@neighbors$RNA.nn@nn.idx[,c(1,2)]
    df = cbind.data.frame("Cell_geno" = tmp$genotype[NN[,1]],"Neighbor_geno"= tmp$genotype[NN[,2]], "Sample" = tmp$orig.ident[NN[,1]])
    df$same = ifelse(df$Cell_geno==df$Neighbor_geno,1,0)
    for (x in 1:length(unique(df$Sample))) {
      tmp = df[df$Sample==unique(df$Sample)[x],]
      cell_type = c(cell_type,cell_types[i])
      samples = c(samples,unique(df$Sample)[x])
      value = c(value,sum(tmp$same)/nrow(tmp))
    }
  }}
results = cbind.data.frame(cell_type, samples, value)



cell_type = c()
samples = c()
value =c()
for (s in 1:100) {
  for (i in 1:length(unique(cell_types))) {
    tmp = subset(merged.obj, idents = cell_types[i])
    mingenotype = min(table(tmp$genotype))
    cell.list = c(sample(colnames(tmp)[tmp$genotype==genotypes[1]],mingenotype),
                  sample(colnames(tmp)[tmp$genotype==genotypes[2]],mingenotype))
    tmp = subset(tmp, cells = cell.list)
    tmp = FindNeighbors(tmp, reduction = "harmony", dims = 1:10, return.neighbor = T)
    NN = tmp@neighbors$RNA.nn@nn.idx[,c(1,2)]
    df = cbind.data.frame("Cell_geno" = sample(tmp$genotype[NN[,1]]),"Neighbor_geno"= tmp$genotype[NN[,2]], "Sample" = tmp$orig.ident[NN[,1]])
    df$same = ifelse(df$Cell_geno==df$Neighbor_geno,1,0)
    for (x in 1:length(unique(df$Sample))) {
      tmp = df[df$Sample==unique(df$Sample)[x],]
      cell_type = c(cell_type,cell_types[i])
      samples = c(samples,unique(df$Sample)[x])
      value = c(value,mean(tmp$same))
    }
  }
}
results2 = cbind.data.frame(cell_type, samples, value)

library(ggbeeswarm)
p = ggbarplot(results, x = "cell_type", y = "value", fill = "lightblue",
              title = "Segregation by genotype\n within each cell type", 
              ylab = "Proportion of nearest neighbors in same genotype",
              xlab = F, add = c("mean_se"))+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_y_continuous(limits=c(0,1), expand = c(0, 0),breaks = scales::pretty_breaks(n = 10))+
  geom_hline(yintercept = mean(results2$value),linetype = "dashed")


# p = ggbarplot(results, x = "cell_type", y = "value", fill = "lightblue",
#               title = "Segregation by genotype\n within each cell type", 
#               ylab = "Proportion of nearest neighbors in same genotype",
#               xlab = F, add = c("mean_se","jitter"))+ geom_jitter(alpha = 0.2, width = 0.3)+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
#   scale_y_continuous(limits=c(min(results$value)-0.1,max(results$value)+0.05), oob=rescale_none)+
#   geom_hline(yintercept = mean(results2$value),linetype = "dashed")+geom_hline(yintercept = 0.25,linetype = "solid")
ggsave(plot = p, "plots/Figure3B_3.pdf", width = 4, height = 7, units = "in")


## 3F
p = DimPlot(MG, group.by = "seurat_clusters", label = T)
data = p$data
data$circle = NA
data[data$seurat_clusters%in%c(1,6),]$circle = "Y"

p = ggplot(data, aes(umap_1, umap_2, color = seurat_clusters))+geom_point(size =0.3)+theme_classic()+ggtitle("")
p = LabelClusters(plot = p, id = "seurat_clusters", color = "black")+geom_mark_hull(data = data[data$circle=="Y",], aes(group = circle), expand=0.01, colour="black", size = 1, linetype = "dashed")+NoLegend()



p2 = DimPlot(MG, group.by = "classification", split.by = "genotype", ncol = 2)
data = p2$data
data2 = cbind.data.frame("xx" = ggplot_build(p)$data[[3]]$x, "yy" = ggplot_build(p)$data[[3]]$y)
ggplot(data, aes(umap_1, umap_2, color = classification))+geom_point(size =0.3)+theme_classic()+facet_wrap(~genotype, ncol =2)+
  geom_mark_hull(data = data2,aes(xx, yy), expand=0.01, colour="black", size = 1, linetype = "dashed")+
  theme(strip.background = element_blank())+NoLegend()+ggtitle("genotype") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))+
  theme(strip.text = element_text(size = 12, face = "bold"))
ggsave("plots/Figure3F_2.pdf", width = 5, height = 5, units = "in")

## 3G
metadf = MG@meta.data
metadf$classification = factor(metadf$classification, levels = c("DAM_depleted", "DAM"))
df = metadf %>% dplyr::group_by(orig.ident) %>% count(classification, genotype)%>%
  mutate(freq = n / sum(n))
stat.test = df %>% group_by(classification) %>%
  t_test(freq ~ genotype) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj")
stat.test = stat.test[stat.test$p.adj<0.05,]

stat.test <- stat.test %>%
  add_xy_position(x = "classification")

ggplot(df, aes(classification, freq)) +
  geom_boxplot(aes(fill = genotype)) + 
  theme_bw() + 
  stat_pvalue_manual(stat.test,   label = "p.adj.signif", tip.length = 0.005, y.position = c(1.05,1.1,1.15,1.05,1.1,1.15))+
  scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1,1.5))+theme_classic2()+ylab("composition (%)")
ggsave("plots/Figure3G_2.pdf", width = 5, height = 5, units = "in")


# Volcano plots of Mapt vs WT, 6C;Mapt vs WT, and Mapt vs 6C;Mapt
MAPT_WT = readRDS("data/DEG_Mapt_WT_in_MG.rds")

DEGs = FindMarkers(MG, ident.1 = "Tg-6C KO",ident.2 = "Ntg-WT",test.use = "MAST", latent.vars = "orig.ident")
DEGs$gene = rownames(DEGs)
write.xlsx(DEGs, "data/Mapt6c_vs_WT_MG_DEGs.xlsx")
saveRDS(DEGs, "data/DEG_Mapt6c_WT_in_MG.rds")
MT6C_WT = DEGs

MAPT_MT6C = readRDS("data/DEG_Mapt_6CMT_in_MG.rds")


p = EnhancedVolcano(MAPT_WT,
                lab = rownames(MAPT_WT),
                x = 'avg_log2FC',
                y = 'p_val_adj', 
                title = "Tg-WT vs Ntg-WT",
                subtitle = "",gridlines.major = FALSE,
                gridlines.minor = FALSE)
ggsave(plot = p, filename = "plots/Extended_Figure2A_1.pdf", width = 7, height = 7)


EnhancedVolcano(MAPT_MT6C,
                lab = rownames(MAPT_MT6C),
                x = 'avg_log2FC',
                y = 'p_val_adj', 
                title = "Tg-6C KO vs Ntg-WT",
                subtitle = "",gridlines.major = FALSE,
                gridlines.minor = FALSE,
)
ggsave("plots/Extended_Figure2A_2.pdf", width = 7, height = 7)


