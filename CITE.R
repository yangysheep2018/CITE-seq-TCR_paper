setwd("~/Desktop/Bomi")

dt = readRDS(file = "/Users/yanyang/Desktop/Bomi/finalmeta.rds")

daa = Read10X_h5("10_filtered_feature_bc_matrix.h5")
da = CreateSeuratObject(counts = daa$`Gene Expression`)

####
da[["percent.mt"]] <- PercentageFeatureSet(object = da, pattern = "^MT-")
da[["sample"]] <- gsub(".+-", "", rownames(da@meta.data))
da[["sample"]] <- mapvalues(as.matrix(da[["sample"]]), 
                            from=c("1","2","3","4","5","6","7","8","9","10"), 
                            to=c("Her1","Her2","Her3","Her4","Her5","Ctl1","Ctl2","Ctl3","Idio1","Idio3"))
da[["group"]] = str_extract(as.matrix(da[["sample"]]), "[aA-zZ]+")

da[["ADT"]] <- CreateAssayObject(counts = daa$`Antibody Capture`)

test = read.csv("~/Desktop/Bomi/Book2.csv")
da[["check1"]] = rownames(da@meta.data) %in% test$names
da <- subset(da, subset = check1 == "TRUE")

##normalize ADT data
da <- NormalizeData(da, assay = "ADT", normalization.method = "CLR")
da <- ScaleData(da, assay = "ADT")

colnames(ds.read)[16:ncol(ds.read)] = p2  #"--" or "-" cannot imported in ggplot

###########Feature plot############
dt[["barcode"]] <- gsub(".+_", "", rownames(dt@meta.data))
dt[["check"]] <- grepl("10Sample", rownames(dt@meta.data))
dt = subset(dt,subset = check == "TRUE")

us = Embeddings(dt[["umap"]])
cls = dt@meta.data
ds = merge(cls, us, by="row.names")

nord = da@assays$ADT@data
f1 = as.data.frame(t(nord))

ds.read = merge(ds, f1, by.x = "barcode", by.y = "row.names") 
write.csv(ds.read,"~/Desktop/Bomi/CITE-seq/13protein_expression_Metadata.csv")

ds.read = read.csv("~/Desktop/Bomi/CITE-seq/13protein_expression_Metadata.csv")
ds.read = ds.read[ds.read$seurat_clusters %in% c(0,1,2,16),]
colnames(ds.read)[17:ncol(ds.read)] = p2
for (i in 17:ncol(ds.read)) {
  print(i)
  p <- ggplot(data=ds.read, aes(x=UMAP_1, y=UMAP_2)) + geom_point(aes_string(color=colnames(ds.read)[i]), size = 0.5) + 
    scale_colour_gradient(low = "grey80", high = "blue") + 
    labs(title=paste("Protein_",colnames(ds.read)[i],sep = "")) + 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.text.x = element_text(size = 12,color = "black"),
          axis.text.y = element_text(size = 12,color = "black"),
          axis.line = element_line(colour = "black"),
          axis.title=element_text(size=14),
          legend.background = element_rect(colour = NA))
  ggsave(paste("~/Desktop/Bomi/CITE-seq/Tcells/", colnames(ds.read)[i], ".pdf", sep = ""), plot = p, width = 8, height = 6)
  
}
ds.read$seurat_clusters = as.factor(ds.read$seurat_clusters)
p<-ggplot(ds.read, aes_string(x="seurat_clusters", y=colnames(ds.read)[17], fill=colnames(ds.read)[11])) +
  geom_violin(trim=FALSE)

ggsave(paste("~/Desktop/Bomi/CITE-seq/test",  ".pdf", sep = ""), plot = p, width = 5, height = 1)
#######heatmap#####
da[["old_cluster"]] <- dt@meta.data$seurat_clusters

mapal <- colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(256)

Idents(object = da) <- "old_cluster"

da.small <- subset(da, downsample = 300)
adt.markers <- FindAllMarkers(da.small, assay = "ADT", only.pos = TRUE)
DoHeatmap(da.small, features = unique(adt.markers$gene), assay = "ADT", angle = 90) + scale_fill_gradientn(colours = rev(mapal))

cluster.averages <- AverageExpression(da, assay = "ADT",return.seurat = TRUE)
DoHeatmap(cluster.averages, features = p1, size = 3, draw.lines = FALSE) + scale_fill_gradientn(colours = rev(mapal))

# Run PCA, select 13 PCs for tSNE visualization and graph-based clustering
da <- RunPCA(da, verbose = FALSE)
ElbowPlot(da, ndims = 50)

da <- FindNeighbors(da, dims = 1:10)
da <- FindClusters(da, resolution = 0.5)
da <- RunTSNE(da, dims = 1:10)
da <- RunUMAP(da,dims = 1:10)
UMAPPlot(da)
UMAPPlot(da,group.by = "group")

##############

one <- CreateSeuratObject(counts = daa$`Antibody Capture`)

test = read.csv("~/Desktop/Bomi/Book2.csv")
one[["check1"]] = rownames(one@meta.data) %in% test$names
one <- subset(one, subset = check1 == "TRUE")

one <- NormalizeData(one)
one <- FindVariableFeatures(one)
all.genes <- rownames(one)
one <- ScaleData(one, features = all.genes)
one <- RunPCA(one, features = VariableFeatures(object = one))
one <- FindNeighbors(one, dims = 1:10)
one <- FindClusters(one, resolution = 0.5)
#tSNE
one <- RunTSNE(one, dims = 1:10)
TSNEPlot(one)
#uMap
one <- RunUMAP(one,dims = 1:10)
UMAPPlot(one)




##############

p2 = c("CD3_UCHT1_TSC0034",
       "CD19_HIB19_TSC0050",
       "CD25_BC96_TSC0085",
       "HLA_DR_L243_TSC0159",
       "CD45_LCA_2D1_TSC0048",
       "CD4_T4_Sk3_TSC0045",
       "CD8_T8_SK1_TSC0046",
       "CD56_NCAM_5_1H11_TSC0047",
       "CD11c_CR4_S_HCL_3_TSC0053",
       "CD11b_CR3_M1_70_TSC0014",
       "CD206_MMR_15_2_TSC0205",
       "CD197_CCR7_G043H7_TSC0148",
       "CD45RA_GP180_HI100_TSC0063")

p1 = c("CD3--UCHT1-TSC0034",
      "CD19--HIB19-TSC0050",
      "CD25--BC96-TSC0085",
      "HLA-DR--L243-TSC0159",
      "CD45-LCA--2D1-TSC0048",
      "CD4-T4--Sk3-TSC0045",
      "CD8-T8--SK1-TSC0046",
      "CD56-NCAM--5-1H11-TSC0047",
      "CD11c-CR4--S-HCL-3-TSC0053",
      "CD11b-CR3--M1-70-TSC0014",
      "CD206-MMR--15-2-TSC0205",
      "CD197-CCR7--G043H7-TSC0148",
      "CD45RA-GP180--HI100-TSC0063")

################RidgePlot################
plot_list = list()
for (i in 1:length(p1)) {
  # filename = paste("/Users/yanyang/Desktop/Bomi/Aug5_subcluster/features/",i,".jpg",sep = "")
  # jpeg(filename, width = 650, height = 550)
  p = RidgePlot(da, features = p1[i])
  plot_list[[i]] = p
}

for (i in 1:length(p1)) {
  filename = paste("/Users/yanyang/Desktop/Bomi/CITE/ridgeplot/",p1[i],".pdf",sep = "")
  pdf(filename)
  print(plot_list[[i]])
  dev.off()
}
dev.off()




DefaultAssay(da) <- "ADT"
da <- RunPCA(da, features = rownames(da), reduction.name = "pca_adt", reduction.key = "pca_adt_", 
               verbose = FALSE)
DimPlot(da, reduction = "pca_adt")


adt.data <- GetAssayData(da, slot = "data")
adt.dist <- dist(t(adt.data))

# Before we recluster the data on ADT levels, we'll stash the RNA cluster IDs for later
da[["rnaClusterID"]] <- Idents(da)

# Now, we rerun tSNE using our distance matrix defined only on ADT (protein) levels.
da[["tsne_adt"]] <- RunUMAP(adt.dist, assay = "ADT", reduction.key = "adtTSNE_")
da[["adt_snn"]] <- FindNeighbors(adt.dist)$snn
da <- FindClusters(da, resolution = 0.5, graph.name = "adt_snn")

# We can compare the RNA and protein clustering, and use this to annotate the protein clustering
# (we could also of course use FindMarkers)
clustering.table <- table(Idents(cbmc), cbmc$rnaClusterID)
clustering.table

DimPlot(da, reduction = "tsne_adt", pt.size = 0.5) + NoLegend()

us1 = Embeddings(da[["tsne_adt"]])

dr = read.csv("~/Desktop/Bomi/CITE-seq/13protein_expression_Metadata.csv")
ds.read1 = merge(dr, us1, by.x = "barcode", by.y  = "row.names")
write.csv(ds.read1,"~/Desktop/Bomi/CITE-seq/Metadata.csv")
ds.read1 = read.csv("~/Desktop/Bomi/CITE-seq/Metadata.csv")

ggplot() +
  geom_point(aes(x=UMAP_1, y=UMAP_2),data=ds.read1,size=0.1,color = "#00BFC4") +
  geom_point(aes(x=adtTSNE_1, y=adtTSNE_2),data=ds.read1,size=0.1,colour = "#F8766D") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.text.x = element_text(size = 12,color = "black"),
        axis.text.y = element_text(size = 12,color = "black"),
        axis.line = element_line(colour = "black"),
        axis.title=element_text(size=14),
        legend.background = element_rect(colour = NA))


for (i in 18:ncol(ds.read1)-2) {
  print(i)
  p <- ggplot(data=ds.read1, aes(x=adtTSNE_1, y=adtTSNE_2)) + geom_point(aes_string(color=colnames(ds.read1)[i]), size = 0.5) + 
    scale_colour_gradient(low = "grey80", high = "blue") + 
    labs(title=paste("Protein_",colnames(ds.read1)[i],sep = "")) + 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.text.x = element_text(size = 12,color = "black"),
          axis.text.y = element_text(size = 12,color = "black"),
          axis.line = element_line(colour = "black"),
          axis.title=element_text(size=14),
          legend.background = element_rect(colour = NA))
  ggsave(paste("~/Desktop/Bomi/CITE-seq/Updated/Basedonnew/", colnames(ds.read1)[i], ".pdf", sep = ""), plot = p, width = 8, height = 6)
  
}



#########Myeloids cells features plot, heatmap, and violin plot##################
#########Myeloids cells features plot, heatmap, and violin plot##################
#########Myeloids cells features plot, heatmap, and violin plot##################

setwd("~/Desktop/Bomi")

daa = Read10X_h5("10_filtered_feature_bc_matrix.h5")
da = CreateSeuratObject(counts = daa$`Gene Expression`)

da[["percent.mt"]] <- PercentageFeatureSet(object = da, pattern = "^MT-")
da[["sample"]] <- gsub(".+-", "", rownames(da@meta.data))
da[["sample"]] <- mapvalues(as.matrix(da[["sample"]]), 
                            from=c("1","2","3","4","5","6","7","8","9","10"), 
                            to=c("Her1","Her2","Her3","Her4","Her5","Ctl1","Ctl2","Ctl3","Idio1","Idio3"))
da[["group"]] = str_extract(as.matrix(da[["sample"]]), "[aA-zZ]+")

da[["ADT"]] <- CreateAssayObject(counts = daa$`Antibody Capture`)

#Book2 is the way to get rid of low quality cells
test = read.csv("~/Desktop/Bomi/Book2.csv")
da[["check1"]] = rownames(da@meta.data) %in% test$names
da <- subset(da, subset = check1 == "TRUE")

da <- NormalizeData(da, assay = "ADT", normalization.method = "CLR")
da <- ScaleData(da, assay = "ADT")

####modified to remove the "2Sample" and seperate "10Sample"(S col) and "barcode" (B col)
dt = readRDS(file = "/Users/yanyang/Desktop/Bomi/myeloids.rds")
us = Embeddings(dt[["umap"]])
cls = dt@meta.data
ds = merge(cls, us, by="row.names")
write.csv(ds,"~/Desktop/Bomi/CITE-seq/MyeloidCells/dt_Mycells_1.csv")
##Modify as above described
ds = read.csv("~/Desktop/Bomi/CITE-seq/MyeloidCells/dt_Mycells_1.csv")
da[["check2"]] = rownames(da@meta.data) %in% ds$B
da <- subset(da, subset = check2 == "TRUE")

nord = da@assays$ADT@data
f1 = as.data.frame(t(nord))

ds.read = merge(ds, f1, by.x = "B", by.y = "row.names") 
write.csv(ds.read,"~/Desktop/Bomi/CITE-seq/MyeloidCells/my_13protein_expression_Metadata.csv")
ds.read = read.csv("~/Desktop/Bomi/CITE-seq/MyeloidCells/my_13protein_expression_Metadata.csv")

##########featureplot and violinplot
for (i in 18:ncol(ds.read)) {
  print(i)
  p <- ggplot(data=ds.read, aes(x=UMAP_1, y=UMAP_2)) + geom_point(aes_string(color=colnames(ds.read)[i]), size = 0.5) + 
    scale_colour_gradient(low = "grey80", high = "blue") + 
    labs(title=paste("Protein_",colnames(ds.read)[i],sep = "")) + 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.text.x = element_text(size = 12,color = "black"),
          axis.text.y = element_text(size = 12,color = "black"),
          axis.line = element_line(colour = "black"),
          axis.title=element_text(size=14),
          legend.background = element_rect(colour = NA))
  ggsave(paste("~/Desktop/Bomi/CITE-seq/MyeloidCells/featureplot/", colnames(ds.read)[i], ".pdf", sep = ""), plot = p, width = 8, height = 6)
}

ds.read$seurat_clusters = as.factor(ds.read$seurat_clusters)
for (i in 16:ncol(ds.read)) {
  print(i)
  p<-ggplot(ds.read, aes_string(x="seurat_clusters", y=colnames(ds.read)[i], fill="seurat_clusters")) +
    geom_violin(trim=FALSE)+
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.text.x = element_text(size = 12,color = "black"),
          axis.text.y = element_text(size = 12,color = "black"),
          axis.line = element_line(colour = "black"),
          axis.title=element_text(size=14),
          legend.background = element_rect(colour = NA))
  ggsave(paste("~/Desktop/Bomi/CITE-seq/MyeloidCells/violinplot/", colnames(ds.read)[i], ".pdf", sep = ""), plot = p, width = 5, height = 1)
}


##########heatmap
p1 = c("CD3--UCHT1-TSC0034",
       "CD19--HIB19-TSC0050",
       "CD25--BC96-TSC0085",
       "HLA-DR--L243-TSC0159",
       "CD45-LCA--2D1-TSC0048",
       "CD4-T4--Sk3-TSC0045",
       "CD8-T8--SK1-TSC0046",
       "CD56-NCAM--5-1H11-TSC0047",
       "CD11c-CR4--S-HCL-3-TSC0053",
       "CD11b-CR3--M1-70-TSC0014",
       "CD206-MMR--15-2-TSC0205",
       "CD197-CCR7--G043H7-TSC0148",
       "CD45RA-GP180--HI100-TSC0063")

####make sure the barcode order is the same
ds.read_test = ds.read[match(rownames(da@meta.data),ds.read$B),]

da[["old_cluster"]] <- ds.read_test$seurat_clusters
Idents(object = da) <- "old_cluster"

mapal <- colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(256)
cluster.averages <- AverageExpression(da, assay = "ADT",return.seurat = TRUE)
DoHeatmap(cluster.averages, features = p1, size = 3, draw.lines = FALSE) + scale_fill_gradientn(colours = rev(mapal))

#########Mycells features plot, heatmap, and violin plot##################
########## DONE DONE DONE DONE DONE DONE DONE DONE DONE ##################



#########T-cells features plot, heatmap, and violin plot##################
#########T-cells features plot, heatmap, and violin plot##################
#########T-cells features plot, heatmap, and violin plot##################
#########T-cells features plot, heatmap, and violin plot##################

setwd("~/Desktop/Bomi")

daa = Read10X_h5("10_filtered_feature_bc_matrix.h5")
da = CreateSeuratObject(counts = daa$`Gene Expression`)

da[["percent.mt"]] <- PercentageFeatureSet(object = da, pattern = "^MT-")
da[["sample"]] <- gsub(".+-", "", rownames(da@meta.data))
da[["sample"]] <- mapvalues(as.matrix(da[["sample"]]), 
                            from=c("1","2","3","4","5","6","7","8","9","10"), 
                            to=c("Her1","Her2","Her3","Her4","Her5","Ctl1","Ctl2","Ctl3","Idio1","Idio3"))
da[["group"]] = str_extract(as.matrix(da[["sample"]]), "[aA-zZ]+")

da[["ADT"]] <- CreateAssayObject(counts = daa$`Antibody Capture`)

#Book2 is the way to get rid of low quality cells
test = read.csv("~/Desktop/Bomi/Book2.csv")
da[["check1"]] = rownames(da@meta.data) %in% test$names
da <- subset(da, subset = check1 == "TRUE")

da <- NormalizeData(da, assay = "ADT", normalization.method = "CLR")
da <- ScaleData(da, assay = "ADT")

####modified to remove the "2Sample" and seperate "10Sample"(S col) and "barcode" (B col)
dt = readRDS(file = "/Users/yanyang/Desktop/Bomi/part2_Tcell.rds")
us = Embeddings(dt[["umap"]])
cls = dt@meta.data
ds = merge(cls, us, by="row.names")
write.csv(ds,"~/Desktop/Bomi/CITE-seq/Tcells/dt_Tcell_1.csv")

ds = read.csv("~/Desktop/Bomi/CITE-seq/Tcells/dt_Tcell_1.csv")
da[["check2"]] = rownames(da@meta.data) %in% ds$B
da <- subset(da, subset = check2 == "TRUE")

nord = da@assays$ADT@data
f1 = as.data.frame(t(nord))

ds.read = merge(ds, f1, by.x = "B", by.y = "row.names") 
write.csv(ds.read,"~/Desktop/Bomi/CITE-seq/Tcells/CD_13protein_expression_Metadata.csv")
ds.read = read.csv("~/Desktop/Bomi/CITE-seq/Tcells/CD_13protein_expression_Metadata.csv")

##########featureplot and violinplot
for (i in 16:ncol(ds.read)) {
  print(i)
  p <- ggplot(data=ds.read, aes(x=UMAP_1, y=UMAP_2)) + geom_point(aes_string(color=colnames(ds.read)[i]), size = 0.5) + 
    scale_colour_gradient(low = "grey80", high = "blue") + 
    labs(title=paste("Protein_",colnames(ds.read)[i],sep = "")) + 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.text.x = element_text(size = 12,color = "black"),
          axis.text.y = element_text(size = 12,color = "black"),
          axis.line = element_line(colour = "black"),
          axis.title=element_text(size=14),
          legend.background = element_rect(colour = NA))
  ggsave(paste("~/Desktop/Bomi/CITE-seq/Tcells/featureplot/", colnames(ds.read)[i], ".pdf", sep = ""), plot = p, width = 8, height = 6)
}

ds.read$seurat_clusters = as.factor(ds.read$seurat_clusters)
for (i in 16:ncol(ds.read)) {
  print(i)
  p<-ggplot(ds.read, aes_string(x="seurat_clusters", y=colnames(ds.read)[i], fill="seurat_clusters")) +
    geom_violin(trim=FALSE)+
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.text.x = element_text(size = 12,color = "black"),
          axis.text.y = element_text(size = 12,color = "black"),
          axis.line = element_line(colour = "black"),
          axis.title=element_text(size=14),
          legend.background = element_rect(colour = NA))
  ggsave(paste("~/Desktop/Bomi/CITE-seq/Tcells/violinplot/", colnames(ds.read)[i], ".pdf", sep = ""), plot = p, width = 5, height = 1)
}


##########heatmap
p1 = c("CD3--UCHT1-TSC0034",
       "CD19--HIB19-TSC0050",
       "CD25--BC96-TSC0085",
       "HLA-DR--L243-TSC0159",
       "CD45-LCA--2D1-TSC0048",
       "CD4-T4--Sk3-TSC0045",
       "CD8-T8--SK1-TSC0046",
       "CD56-NCAM--5-1H11-TSC0047",
       "CD11c-CR4--S-HCL-3-TSC0053",
       "CD11b-CR3--M1-70-TSC0014",
       "CD206-MMR--15-2-TSC0205",
       "CD197-CCR7--G043H7-TSC0148",
       "CD45RA-GP180--HI100-TSC0063")

ds.read_test = ds.read[match(rownames(da@meta.data),ds.read$B),]
da[["old_cluster"]] <- ds.read_test$seurat_clusters
Idents(object = da) <- "old_cluster"

mapal <- colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(256)
cluster.averages <- AverageExpression(da, assay = "ADT",return.seurat = TRUE)
DoHeatmap(cluster.averages, features = p1, size = 3, draw.lines = FALSE) + scale_fill_gradientn(colours = rev(mapal))

#########T cells features plot, heatmap, and violin plot##################
########## DONE DONE DONE DONE DONE DONE DONE DONE DONE ##################