library("Seurat")
library("patchwork")
library("clustree")
library("magrittr")
library("dplyr")



#these are the cellranger files, processed by medgenome "Neil Patel"
files.h5 = list.files(path="/home/kpradhan/mnt/hpc_home/projects/divij/cc1339/processed/CELLPLEX/R_results_CELLPLEX/per_sample_outs", pattern="sample_filtered_feature_bc_matrix", recursive=T, full=T)
data.10x = lapply(files.h5, Read10X_h5)

f = files.h5[1]
#get the sample name from the folder
dnames = gsub(files.h5, pat=".*/(S[0-9])/.*", rep="\\1")

#get list of gene matrices and adt matricx
data.gex = lapply(data.10x, . %>% '[['(1) %>% CreateSeuratObject)
assay.adt = lapply(data.10x, . %>% '[['(2) %>% CreateAssayObject)

#combine the adt and gex togther
for (i in seq_along(data.gex)){
    data.gex[[i]][['adt']] = assay.adt[[i]]
}

#add sample ID to metadata of seurat objects
for (i in seq_along(dnames)){
    data.gex[[i]]@meta.data$sample = dnames[i]
}

# normalize and identify variable features for each dataset independently
x.list <- lapply(X = data.gex, FUN = function(x) {
    x <- NormalizeData(x)
    x = ScaleData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
    x = RunPCA(x, verbose = F, npcs = 20)
    x = RunUMAP(x, dims = 1:10, verbose = F)
    x
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = x.list)

x.anchors <- FindIntegrationAnchors(object.list = x.list, anchor.features = features)

# this command creates an 'integrated' data assay
x.combined <- IntegrateData(anchorset = x.anchors)

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(x.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
x.combined <- ScaleData(x.combined, verbose = FALSE)
x.combined <- RunPCA(x.combined, npcs = 30, verbose = FALSE)
x.combined <- RunUMAP(x.combined, reduction = "pca", dims = 1:30)
x.combined <- FindNeighbors(x.combined, reduction = "pca", dims = 1:30)
x.combined <- FindClusters(x.combined, resolution = seq(0.1, 1, 0.1))
head(x.combined@meta.data)


#also scale and noramlzie the adt
#x.combined[["adt"]]@data
#x.combined[["adt"]]@counts
DefaultAssay(x.combined) <- "adt"
x.combined <-NormalizeData(x.combined, assay = "adt",  normalization.method = 'CLR', margin = 2)
x.combined <- ScaleData(x.combined, verbose = FALSE)

DefaultAssay(x.combined) <- "integrated"



png("clustree.png", 1000, 1000)
clustree(x.combined)
dev.off()

png("elbow.png")
ElbowPlot(x.combined)
dev.off()


#after looking at the clustree, 0.5 might be a good resolution
Idents(x.combined) = "integrated_snn_res.0.5"


head(x.combined@meta.data)
table(x.combined@meta.data$integrated_snn_res.0.5)
# Visualization
p1 <- DimPlot(x.combined, reduction = "umap", group.by = "sample")
p2 <- DimPlot(x.combined, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2
png("umap_clust_res0.5.png")
p2
dev.off()




###########################################33
#try automatic annotation with celldex and singler
library("celldex")
library("SingleR")
library("SingleCellExperiment")
ref <- MouseRNAseqData()
sce <- as.SingleCellExperiment(DietSeurat(x.combined))
pred <- SingleR(test=sce, ref=ref, labels=ref$label.main)
pred.fine <- SingleR(test=sce, ref=ref, labels=ref$label.fine)

tab1 = table(pred$pruned.labels)
tab2 = table(pred.fine$pruned.labels)
#ignore the cells with less than 100 countso


x.combined@meta.data$pred1 = pred$pruned.labels
x.combined@meta.data$pred2 = pred.fine$pruned.labels


#save the dataset
save(file="x.combined.rdata", x.combined)


x.combined <- SetIdent(x.combined, value = "integrated_snn_res.0.5")
png("umap_res0.5_clusters.png")
DimPlot(x.combined, label = TRUE, repel = TRUE)
dev.off()

x.combined <- SetIdent(x.combined, value = "pred1")
png("umap_celltypes.png")
DimPlot(x.combined, label = TRUE, repel = TRUE)
dev.off()

x.combined <- SetIdent(x.combined, value = "pred2")
png("umap_celltypes2.png", 1000, 1000)
DimPlot(x.combined, label = TRUE, repel = TRUE)
dev.off()

DimPlot(srat, label = T , repel = T, label.size = 3)
t1 = table(pred$labels)
which.max(t1)
head(x.combined@meta.data)

#find the celltype tha appears most in each cluster
sink("cluster0.5_celltypes2.txt")
by( x.combined@meta.data$pred2, 
    x.combined@meta.data$integrated_snn_res.0.3, 
    FUN=function(a){
        t1 = table(a)
        t1 = t1/sum(t1)
        t1[which.max(t1)] 
        data.frame(celltype=names(t1)[which.max(t1)], proportion=max(t1))
        sort(t1)
    })
sink()


x.combined[['adt']]@counts

save(file="x.combined.rdata", "x.combined")

grep(rownames(x.combined[['RNA']]), pat="tet2", ignore.case=T, value=T)
grep(rownames(x.combined[['RNA']]), pat="cd45", ignore.case=T, value=T)
grep(rownames(x.combined[['adt']]), pat="cd45", ignore.case=T, value=T)
grep(rownames(x.combined[['adt']]), pat="tet2", ignore.case=T, value=T)

#make a few feature plots on ADT
p1 <- FeaturePlot(x.combined, "rna_Ptprc", split.by="sample") 
p2 <- FeaturePlot(x.combined, "adt_Ms-CD45-2-TotalA", split.by="sample") 
p1 / p2

pdf("feature_example1.pdf", width=10, height=5)
p1 / p2
dev.off()

