library(Seurat)
library(dplyr)
library(hdf5r)
require(plyr)
library(stringr)
library(ggplot2)
library("ggrepel")
library(pheatmap)
library("genefilter")
library(DataCombine)
library(patchwork)
###devtools::install_github('chris-mcginnis-ucsf/DoubletFinder')
###devtools::install_github('dviraran/SingleR')
##devtools::install_github(repo = "mojaveazure/loomR", ref = "develop")
##devtools::install_github("cellgeni/sceasy")
library(DoubletFinder)
library(scDblFinder)
library(SingleR)
library(tibble)
library(scRNAseq)
library(reshape2)
require(data.table) 
library(RColorBrewer)
library(reticulate)
library(loomR)
# remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
# BiocManager::install("LoomExperiment")
# BiocManager::install("GWASTools")

#########################################################################################################
#########################################################################################################
#####QC#####
dt = Read10X_h5("10_filtered_feature_bc_matrix.h5")
dt <- dt$`Gene Expression`
dt = CreateSeuratObject(counts = dt,min.cells = 3, min.features = 200)

dt2 = Read10X_h5("2_filtered_feature_bc_matrix.h5")
dt2 = CreateSeuratObject(counts = dt2,min.cells = 3, min.features = 200)


# merge into one single seurat object. Add cell ids just in case you have overlapping barcodes between the datasets.


dim(dt@meta.data) #27432  to 33301
dim(dt2@meta.data) #5869

dt[["percent.mt"]] <- PercentageFeatureSet(object = dt, pattern = "^MT-")
dt[["sample"]] <- gsub(".+-", "", rownames(dt@meta.data))
dt[["sample"]] <- mapvalues(as.matrix(dt[["sample"]]), 
                            from=c("1","2","3","4","5","6","7","8","9","10"), 
                            to=c("Her1","Her2","Her3","Her4","Her5","Ctl1","Ctl2","Ctl3","Idio1","Idio3"))
dt[["group"]] = str_extract(as.matrix(dt[["sample"]]), "[aA-zZ]+")

dt2[["percent.mt"]] <- PercentageFeatureSet(object = dt2, pattern = "^MT-")
dt2[["sample"]] <- gsub(".+-", "", rownames(dt2@meta.data))
dt2[["sample"]] <- mapvalues(as.matrix(dt2[["sample"]]), 
                            from=c("1","2"), 
                            to=c("Idio2","Idio4"))
dt2[["group"]] = str_extract(as.matrix(dt2[["sample"]]), "[aA-zZ]+")

dt <- merge(dt, dt2, add.cell.ids=c("10Sample","2Sample"))

VlnPlot(dt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),ncol = 3,pt.size = 0)
VlnPlot(dt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),ncol = 3, group.by = "sample",pt.size = 0)
VlnPlot(dt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),ncol = 3, group.by = "group",pt.size = 0)

######Filter cells#####
dt <- subset(dt, subset = nFeature_RNA > 200 & nFeature_RNA < 50000 & percent.mt < 10)  #33301  to 30860 

# plot1 <- FeatureScatter(dt, feature1 = "nCount_RNA", feature2 = "percent.mt")
# plot2 <- FeatureScatter(dt, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# plot1 + plot2

dt <- NormalizeData(dt)
dt <- FindVariableFeatures(dt)
# top10 <- head(VariableFeatures(dt), 10)
# plot1 <- VariableFeaturePlot(dt)
# plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
# plot2
all.genes <- rownames(dt)
dt <- ScaleData(dt, features = all.genes)
dt <- RunPCA(dt, features = VariableFeatures(object = dt))
# VizDimLoadings(dt, dims = 1:2, reduction = "pca")
# dt <- JackStraw(dt, num.replicate = 100)
# dt <- ScoreJackStraw(dt, dims = 1:20)
# JackStrawPlot(object = dt, dims = 1:15)
# ElbowPlot(dt)
dt <- FindNeighbors(dt, dims = 1:10)
dt <- FindClusters(dt, resolution = 0.5)
#tSNE
dt <- RunTSNE(dt, dims = 1:10)
TSNEPlot(dt)
#uMap
dt <- RunUMAP(dt,dims = 1:10)
UMAPPlot(dt)

Idents(object = dt) <- "sample"
l = unique(dt@meta.data$sample)
l1 = c("Her2","Her3","Her4","Her5","Ctl1","Ctl2","Ctl3","Idio1","Idio3","Idio2","Idio4") 

dt_1 = subset(dt,subset = sample == "Her1")
dt.res <- paramSweep_v3(dt_1, PCs = 1:10, sct = FALSE)
dt.stats <- summarizeSweep(dt.res, GT = FALSE)
bcmvn_kidney <- find.pK(dt.stats)
mpK<-as.numeric(as.vector(bcmvn_kidney$pK[which.max(bcmvn_kidney$BCmetric)]))
annotations <- dt_1@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations) 
nExp_poi <- round(0.075*length(rownames(dt_1@meta.data)))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
dt_1 <- doubletFinder_v3(dt_1, PCs = 1:10, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

for (i in l1) {
  dt_ss = subset(dt,subset = sample == i)
  dt.res <- paramSweep_v3(dt_ss, PCs = 1:10, sct = FALSE)
  dt.stats <- summarizeSweep(dt.res, GT = FALSE)
  bcmvn_kidney <- find.pK(dt.stats)
  mpK<-as.numeric(as.vector(bcmvn_kidney$pK[which.max(bcmvn_kidney$BCmetric)]))
  annotations <- dt_ss@meta.data$ClusteringResults
  homotypic.prop <- modelHomotypic(annotations) 
  nExp_poi <- round(0.075*length(rownames(dt_ss@meta.data)))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  dt_ss <- doubletFinder_v3(dt_ss, PCs = 1:10, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
  
  dt_1 <- merge(dt_1, dt_ss)
}

write.csv(dt_1@meta.data,"finalMetafile_DoubletsAdded.csv")
test = read.csv("finalMetafile_DoubletsAdded.csv")

dt[['Doublet']] = test$DF.classifications_0.25_0.02_110
UMAPPlot(dt, group.by="Doublet")

Idents(object = dt) <- "Doublet"
dt = subset(dt,subset = Doublet == "Singlet")

VlnPlot(dt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),ncol = 3,pt.size = 0)
VlnPlot(dt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),ncol = 3, group.by = "sample",pt.size = 0)
VlnPlot(dt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),ncol = 3, group.by = "group",pt.size = 0)


Idents(object = dt) <- "group"
UMAPPlot(dt,group.by="seurat_clusters")
UMAPPlot(dt, group.by="group")
p2 = DimPlot(dt, reduction = "umap", group.by = "group")

Idents(object = dt) <- "seurat_clusters"

######SIngleR Annotation##########

hpca.se <- HumanPrimaryCellAtlasData()
bpe.se <- BlueprintEncodeData()

test <- as.SingleCellExperiment(dt)
Anno <- SingleR(test=test, assay.type.test=1, ref = list(HP = hpca.se , BP = bpe.se), labels = list(hpca.se$label.main , bpe.se$label.main))

Anno$cluster <- rownames(Anno)
fin <- Anno %>% dplyr::tbl_df() %>% dplyr::select(cluster,labels)
# new.cluster.ids <- fin$labels
# names(new.cluster.ids) <- levels(dt)
# dt <- RenameIdents(dt, new.cluster.ids)

dt[['Anno']] = fin$labels
dt[["Anno"]] <- mapvalues(as.matrix(dt[["Anno"]]), 
                          from=c("B_cell","Endothelial cells","Epithelial cells","Macrophage","Monocyte","NK cells"), 
                          to=c("B-cells","Endothelial_cells","Epithelial_cells","Macrophages","Monocytes","NK_cell"))

UMAPPlot(dt,group.by="Anno")

write.csv(dt@meta.data,"/Users/yanyang/Desktop/Bomi/Aug12_subcluster_Myeloids/UpdatedDoublet_Aug14/final_metafile.csv")
saveRDS(dt, file = "/Users/yanyang/Desktop/Bomi/finalmeta.rds")


for (i in 0:16) {
  cluster1.markers <- FindMarkers(dt, ident.1 = i, min.pct = 0.25)
  name = paste("/Users/yanyang/Desktop/Bomi/Aug12_subcluster_Myeloids/UpdatedDoublet_Aug14/DE/Cluster_",i,"_DEG.csv",sep = "")
  write.csv(cluster1.markers,name, row.names=T)
}



#########
FeaturePlot(dt, features = c("CD3E", "CD8A", "CD68", "CD19", "CD1C", "KIT", "IL1B", "MKI67"))

dt[['count']]  = 1
DT <- data.table(dt@meta.data) 
df = DT[ , .(Totalcount = sum(count)), by = .(sample,seurat_clusters)]
# dff = acast(df, seurat_clusters  ~ sample)
# dff[is.na(dff)] <- 0

colnames(df)[2] <- "Cluster"
DF = setDF(df)
DF$Group = str_extract(DF$sample, "[aA-zZ]+")

# Stacked + percent

colourCount = length(unique(df$Cluster))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
# a <- ifelse(data$category == 0, "red", "blue")
# a = ifelse(DF$Group == 'Ctl', 'red', ifelse(DF$Group == 'Her', "blue", "black"))
a = c("blue","blue","blue","red","red","red","red","red","black","black","black","black")
ggplot(DF, aes(fill=Cluster, y=Totalcount, x=sample)) + 
  geom_bar(position="fill", stat="identity") + scale_fill_manual(values = getPalette(colourCount)) +
  labs(x = 'Samples', y = 'Fraction of total cells') +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.text.x = element_text(size = 12,color = a),
        axis.text.y = element_text(size = 12,color = "black"),
        axis.line = element_line(colour = "black"),
        axis.title=element_text(size=14),
        legend.background = element_rect(colour = NA))

features1<- c("CD3E", "CD8A", "CD4", "CD69", "IL7R", "KLRB1", "CCL5", "CXCR4",
              "GZMA", "CD14", "CD68", "CD163", "ITGAM","ITGAX", "CD1C", 
              "S100A9", "IL1B", "APOE", "MKI67", "KIT","CD19", "CD79A", "IGHA1", "CD34", "PECAM1","COL1A1","LILRA4", "KRT19", "PRSS1")
ts = Embeddings(dt[["tsne"]])
us = Embeddings(dt[["umap"]])
cls = dt@meta.data
ds = merge(cls, ts, by="row.names")
ds = merge(ds, us,by.x = "Row.names", by.y = "row.names")
nord = dt[["RNA"]]@data
mtx.sub = as.data.frame(as.matrix(nord[rownames(nord) %in% features1,]))
mtx.sub11 = as.data.frame(as.matrix(nord))
ds.read = merge(ds, t(mtx.sub), by.x = "Row.names", by.y = "row.names") 
ds.read11 = merge(ds, t(nord), by.x = "Row.names", by.y = "row.names") 
write.csv(ds.read,"/Users/yanyang/Desktop/Bomi/July31/featuresinfo.csv",row.names = F)

###########
table(dt@meta.data$group)

# Ctl   Her  Idio 
# 4169 11786 12592 
#Build subsampled sample to make sure each sample have the same cell # 341
set.seed(66)
ds.gr = {}
n = 4169  #NASH 5881 WT 7359 
for (gr in c("Ctl", "Her","Idio")) {
  sset = ds[ds$group == gr, ]
  sset.sub = sset[sample(nrow(sset), n), ]
  ds.gr = rbind(ds.gr, sset.sub)
}

#Plotting per group
ds.test = ds.gr[ds.gr$group == "Idio",]

ggplot(data=ds.test, aes(x=UMAP_1, y=UMAP_2)) + geom_point(aes(color=as.character(seurat_clusters)), size = 1)  +
  theme_bw() + labs(x="UMAP_1", y="UMAP_2", color="17 Clusters",title="Group Idio") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.text.x = element_text(size = 12,color = "black"),
        axis.text.y = element_text(size = 12,color = "black"),
        axis.line = element_line(colour = "black"),
        axis.title=element_text(size=14),
        legend.title = element_blank(),
        legend.background = element_rect(colour = NA))


#####Calculating the average gene expression within a cluster#####

# dt[["Type"]] <- mapvalues(as.matrix(dt[["seurat_clusters"]]), 
#                             from=c("0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16"), 
#                             to=c("CD8+ T cells","CD4+ T cells","CD8+ cytotoxic","Monocytes","DC","Mac-1 (APOE)","Mono/Mac","Mast cells","B cells","Mac-2 (S100)","Plasma cells","Mac-3 (Ki67)","Endothelial","PSCs","pDC","Epithelial (ductal)","Doublets"))
# dt_test = dt
# Idents(object = dt_test) <- "Type"
# orig.levels <- levels(dt_test)
# Idents(dt_test) <- gsub(pattern = " ", replacement = "_", x = Idents(dt_test))
# orig.levels <- gsub(pattern = " ", replacement = "_", x = orig.levels)
# levels(dt_test) <- orig.levels

levels(dt)
cluster.averages <- AverageExpression(dt, return.seurat = TRUE)

cluster.averages

mapal <- colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(256)
DoHeatmap(cluster.averages, features = features1, size = 3, draw.lines = FALSE) + scale_fill_gradientn(colours = rev(mapal))

p <- ggplot(data=ds1.read, aes(x=tSNE_1, y=tSNE_2)) + geom_point(aes_string(color=colnames(ds.read)[i]), size = 0.5) + 
  theme_bw() + scale_colour_gradient(low = "yellow", high = "blue") + labs(title="Sample IF")

######feature map#########
f1 = c("CD3E", "CD8A", "CD68", "CD19", "CD1C", "KIT", "IL1B", "MKI67")

plot_list = list()
for (i in 1:length(f1)) {
  # filename = paste("/Users/yanyang/Desktop/Bomi/Aug5_subcluster/features/",i,".jpg",sep = "")
  # jpeg(filename, width = 650, height = 550)
  p = FeaturePlot(dt, features = f1[i])
  plot_list[[i]] = p
}

for (i in 1:length(f1)) {
  filename = paste("/Users/yanyang/Desktop/Bomi/Aug12_subcluster_Myeloids/UpdatedDoublet_Aug14/featureplot/",f1[i],".pdf",sep = "")
  pdf(filename)
  print(plot_list[[i]])
  dev.off()
}
dev.off()

####violinplot
for (i in 15:41) {
  print(i)
  p<-ggplot(ds.read, aes_string(x="seurat_clusters", y=colnames(ds.read)[i], fill="seurat_clusters")) +
    geom_violin() +
    xlab(NULL) +
    theme_minimal()
  ggsave(paste("/Users/yanyang/Desktop/Bomi/July31/Violin/", colnames(ds.read)[i], ".pdf", sep = ""), plot = p, width = 8, height = 6)
}

####percentagefigure
ggplot(ds.read, aes(x=as.factor(seurat_clusters), fill=as.factor(group)))+
  geom_bar(aes( y=..count../tapply(..count.., ..x.. ,sum)[..x..]), position="dodge" ) +
  # geom_text(aes( y=..count../tapply(..count.., ..x.. ,sum)[..x..], label=scales::percent(..count../tapply(..count.., ..x.. ,sum)[..x..]) ),
  #           stat="count", position=position_dodge(0.9), vjust=-0.5)+
  ylab('Fraction of cells within each cluster') +
  xlab(NULL)+
  scale_y_continuous(labels = scales::percent) +
  scale_fill_brewer(palette = "Set1")


sdata.loom <- as.loom(x = dt, filename = "/Users/yanyang/Desktop/Bomi/dt_all.loom", verbose = FALSE)
# Always remember to close loom files when done
sdata.loom$close_all()



#####Sep. 16

meta = read.delim("test_meta.txt",sep = "\t")
meta1 = read.delim("cellphonedb_meta.txt",sep = "\t")
count = read.delim("test_counts.txt",sep = "\t")
count1 = read.delim("cellphonedb_count.txt",sep = "\t",check.names = FALSE)

my.symbols = count1$Gene
df <- bitr(my.symbols, fromType = "SYMBOL",
           toType = c("ENSEMBL"),
           OrgDb = org.Hs.eg.db)
count1 = count1[count1$Gene %in% df$SYMBOL,]
df1 = df[!duplicated(df[,'SYMBOL']),]

count11 = count1[match(df1$SYMBOL,count1$Gene),]
count11$Gene = df1$ENSEMBL
write.table(count11,'cellphonedb_count.txt', sep='\t', row.name = FALSE,quote=F)

library(org.Hs.eg.db)
library(clusterProfiler)
dt = readRDS(file = "/Users/yanyang/Desktop/Bomi/finalmeta.rds")
nord = dt[["RNA"]]@data
nord_count = as.data.frame(nord)
nord_count <- cbind(rownames(nord_count), nord_count)
colnames(nord_count)[1] = "Gene"
# count_norm <- apply(nord, 2, function(x) (x/sum(x))*10000)
write.table(nord_count,'cellphonedb_count_all.txt', sep='\t', row.name = FALSE,quote=F)
count1 = read.delim("cellphonedb_count_all.txt",sep = "\t",check.names = FALSE)
my.symbols = count1$Gene
df <- bitr(my.symbols, fromType = "SYMBOL",
           toType = c("ENSEMBL"),
           OrgDb = org.Hs.eg.db)
count1 = count1[count1$Gene %in% df$SYMBOL,]
df1 = df[!duplicated(df[,'SYMBOL']),]

count11 = count1[match(df1$SYMBOL,count1$Gene),]
count11$Gene = df1$ENSEMBL

dt[["A"]] <- mapvalues(as.matrix(dt[["seurat_clusters"]]), 
                            from=c("0","1","2","3","4","5","6","7","8","9",
                                   "10","11","12","13","14","15","16"), 
                            to=c("Tcell","Tcell","Tcell","Myeloids","Myeloids",
                                 "Myeloids","Myeloids","Mastcell","Bcell","Myeloids",
                                 "Bcell","Myeloids","K","K","Bcell","K","Tcell"))

set.seed(1)
dat <- data.frame(col1 = sample(1:4, 10, replace = TRUE))

dat

values <- c("Tcell","Tcell","Tcell","Myeloids","Myeloids",
            "Myeloids","Myeloids","Mastcell","Bcell","Myeloids",
            "Bcell","Myeloids","K","K","Bcell","K","Tcell")
dt[["celltype"]] <- values[as.numeric(dt@meta.data$seurat_clusters)]
saveRDS(dt, file = "/Users/yanyang/Desktop/Bomi/finalmeta.rds")
sdata = dt
iTalk_data <- as.data.frame(t(sdata@assays$RNA@counts))
iTalk_data$cell_type <- sdata@meta.data$celltype
iTalk_data$compare_group <- sdata@meta.data$group

ref =  as.data.frame(dt@meta.data)
colnames(ref)[ncol(ref)] = "AA"
ref = ref[ref$AA != "K",]
ref_test = ref[ref$group == "Ctl",]
count22 = count11[,colnames(count11) %in% c("Gene",rownames(ref_test))]
write.table(count22,'/Users/yanyang/Desktop/Bomi/cellphonedb_count_Idio.txt', sep='\t', row.name = FALSE,quote=F)
meta_data <- cbind(rownames(ref_test), ref_test[,'AA', drop=F])
colnames(meta_data) = c("Cell","group")
write.table(meta_data,'/Users/yanyang/Desktop/Bomi/cellphonedb_meta_Ctl.txt', sep='\t', row.name = FALSE, quote=F)


