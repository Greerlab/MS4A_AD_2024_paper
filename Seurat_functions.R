require(patchwork)
require(ggplot2)
require(Seurat)
multiFeaturePlot <- function(obj, genelist, ncol) {
  plots = list()
  for (i in 1:length(genelist)) {
    plots[[i]] = FeaturePlot(obj, genelist[i], order = T, cols = c("gray","yellow","orange","red"))+NoLegend()
  }
  return(wrap_plots(plots, ncol = ncol))
}

gene_expression_scoring <- function (object, features1, features2, ctrl = NULL, Names, Identname, 
          ...) 
{
  name <- Identname
  features <- list(G1.Score = features1, G2.Score = features2)
  if (is.null(x = ctrl)) {
    ctrl <- min(vapply(X = features, FUN = length, FUN.VALUE = numeric(length = 1)))
  }
  object.cc <- AddModuleScore(object = object, features = features, 
                              name = name, ctrl = ctrl, ...)
  cc.columns <- grep(pattern = name, x = colnames(x = object.cc[[]]), 
                     value = TRUE)
  cc.scores <- object.cc[[cc.columns]]
  rm(object.cc)
  CheckGC()
  assignments <- apply(X = cc.scores, MARGIN = 1, FUN = function(scores, 
                                                                 first = paste0(Names[1]), second = paste0(Names[2]), null = "Unsure") {
    if (all(scores < 0)) {
      return(null)
    }
    else {
      if (length(which(x = scores == max(scores))) > 1) {
        return("Undecided")
      }
      else {
        return(c(first, second)[which(x = scores == 
                                        max(scores))])
      }
    }
  })
  cc.scores <- merge(x = cc.scores, y = data.frame(assignments), 
                     by = 0)
  colnames(x = cc.scores) <- c("rownames", paste0(Names[1]), paste0(Names[2]), paste0(Identname))
  rownames(x = cc.scores) <- cc.scores$rownames
  cc.scores <- cc.scores[, c(paste0(Names[1]), paste0(Names[2]), paste0(Identname))]
  object[[colnames(x = cc.scores)]] <- cc.scores
  return(object)
}

seurat_heatmap <- function(object, subset, assay = "SCT", genelist, ident, 
                           column_title_rot = 90, row_title_rot = 0,
                           show_row_names = T,row_split = NULL
) {
  if (!is.null(subset)) {
    object = subset(object, ident = subset)
  } 
  Idents(object) = object@meta.data[,ident]
  mat = object[[assay]]@data %>% as.matrix()
  mat = mat[rowSums(mat)>0,]
  mat<- t(scale(t(mat)))
  genelist = unique(genelist)[unique(genelist)%in%rownames(mat)]
  mat = mat[genelist,]
  cluster_anno<- Idents(object)
  if (!is.null(row_split)) {
    row_split = row_split[genelist]
    row_split = factor(names(row_split))
  }
  
  limit = quantile(mat, c(0.1, 0.90))
  col_fun = circlize::colorRamp2(c(limit[1], 0, limit[2]), c("#FF00FF", "black", "#FFFF00"))
  return(Heatmap(mat, name = "Expression",  
                column_split = factor(cluster_anno),
                cluster_columns = T,
                show_column_dend = FALSE,
                cluster_column_slices = TRUE,
                #column_title_gp = gpar(fontsize = 8),
                row_title_rot = row_title_rot,
                column_gap = unit(0.5, "mm"),
                cluster_rows = T,
                show_row_dend = FALSE,
                show_row_names = show_row_names,
                cluster_row_slices = T,
                row_split = row_split,
                col = col_fun,
                row_names_gp = gpar(fontsize = 6),
                column_title_rot = column_title_rot,
                #top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = scales::hue_pal()(9)))),
                show_column_names = FALSE,
                use_raster = TRUE,
                raster_quality = 4))
}


compositionplot <- function(object, ident, composition.ident) {
  require(ggplot2)
  require(dplyr)
  require(tidyr)
  metadf = object@meta.data
  x = metadf %>% group_by(orig.ident) %>% count(get(composition.ident), get(ident))%>%
    mutate(freq = n / sum(n))
  x = ungroup(x) %>% complete(`get(composition.ident)` , nesting(orig.ident, `get(ident)`), fill = list(n = 0, freq = 0))
  ggplot(x, aes(x=`get(composition.ident)` , y=freq, fill=`get(ident)`))+
    geom_dotplot(
      aes(color = `get(ident)`),
      binaxis='y', stackdir='center', dotsize = 0.3,
      position = position_dodge(0.8)
    )+ geom_boxplot( aes(color = `get(ident)`), alpha =0, width = 0.5, size = 0.4,
                     position = position_dodge(0.8)
    )+theme_classic()+ylab("composition (%)")+xlab(paste0(composition.ident))+
    scale_colour_discrete(paste0(ident))+
    scale_fill_discrete(paste0(ident))
}

BlendExpression <- function(data) {
  if (ncol(x = data) != 2) {
    stop("'BlendExpression' only blends two features")
  }
  features <- colnames(x = data)
  data <- as.data.frame(x = apply(
    X = data,
    MARGIN = 2,
    FUN = function(x) {
      return(round(x = 9 * (x - min(x)) / (max(x) - min(x))))
    }
  ))
  data[, 3] <- data[, 1] + data[, 2] * 10
  colnames(x = data) <- c(features, paste(features, collapse = '_'))
  for (i in 1:ncol(x = data)) {
    data[, i] <- factor(x = data[, i])
  }
  return(data)
}

#Credit: https://github.com/karthikshekhar/CellTypeMIMB/blob/master/utilities.R
topGOterms = function( fg.genes = NULL,
                       bg.genes = NULL,
                       organism = "Mouse", 
                       ontology.use = "BP",
                       stats.use = "fisher",
                       algorithm.use = "weight01",
                       num.char =100){
  
  if (is.null(fg.genes) | is.null(bg.genes)){
    stop("Error : Both gene lists are empty")
  }
  
  require(topGO)
  if (organism == "Mouse"){
    mapping.use = "org.Mm.eg.db"
    library(org.Mm.eg.db)
  } else if (organism == "Human"){
    mapping.use = "org.Hs.eg.db"
    library(org.Hs.eg.db)
  } else {
    stop("Error : Organisms other than mouse not supported currently")
  }
  
  n = length(bg.genes)
  geneList = integer(n)
  names(geneList) = bg.genes
  geneList[intersect(names(geneList), fg.genes)]=1
  print(paste0("Total ", length(geneList), " genes. ", sum(geneList), " genes in the foreground"))
  geneList = factor(geneList)
  
  if (ontology.use %in% c("BP", "CC", "MF")){
    print(paste0("Using Ontology : ", ontology.use))
  } else {
    stop("Error: Ontology not available. Should be one of BP, CC or MF")
  }
  # Make GO object
  GOdata <- new("topGOdata",
                description = "GOanalysis",
                ontology = ontology.use,
                allGenes = geneList,
                annot = annFUN.org,
                mapping = mapping.use,
                ID = "SYMBOL",
                nodeSize = 10)
  print(paste0("Using the ", stats.use, " statistic with the ", algorithm.use, " algorithm"))
  res.result <- runTest(GOdata, statistic = stats.use, algorithm = algorithm.use)
  to.return = list()
  to.return$GOdata = GOdata
  
  allGO = usedGO(object = GOdata) 
  goEnrichment = GenTable(GOdata, pval = res.result ,topNodes = length(allGO), numChar = num.char)
  goEnrichment$padj = p.adjust(goEnrichment$pval, method = "BH")
  goEnrichment$GeneRatio = round(goEnrichment$Significant/goEnrichment$Annotated, 3)
  goEnrichment_top = goEnrichment[goEnrichment$padj<0.05,]
  goEnrichment_top = goEnrichment_top[order(goEnrichment_top$padj, decreasing = T),]
  goEnrichment_top$Term = factor(goEnrichment_top$Term, levels = goEnrichment_top$Term)
  
  to.return$res.result <- goEnrichment_top
  return(to.return)
}
