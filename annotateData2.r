library("Seurat")
library("patchwork")
library("clustree")
library("magrittr")
library("dplyr")
library("ggplot2")


library("readxl")

load("x.combined.rdata")

#these are the manually annotated celltypes prepared by Divij
cluster.info = data.frame(read_xlsx("clusterCelltypes_20230810.xlsx"))


#set the cluster ids to the new annotated cell types

Idents(x.combined)
id2ctype = list()
id2ctype = cluster.info$Final.Cluster
names(id2ctype) = as.character(cluster.info$cluster)
cluster.info[1:5,]

x.combined$divij.celltypes = id2ctype[as.character(x.combined$integrated_snn_res.0.5)]



##featureplot on ilrap1
grep(rownames(x.combined[["RNA"]]), pat="il1rap", ignore.case=T, value=T)
grep(rownames(x.combined[["adt"]]), pat="il", ignore.case=T, value=T)
p1 <- DimPlot(x.combined, reduction = "umap", label = TRUE, repel = TRUE)
p2 <- FeaturePlot(x.combined, "rna_Il1rap", pt.size=1.5) 
png("umap_il1rap.png", 1000, 500)
p1 + p2
dev.off()

#make featuresplots for all proteins

#x.combined <- SetIdent(x.combined, value = "integrated_snn_res.0.5")
x.combined <- SetIdent(x.combined, value = "divij.celltypes")
for (adt in rownames(x.combined[["adt"]])){
    print(adt)
    p1 <- DimPlot(x.combined, reduction = "umap", label = TRUE, repel = TRUE)
    p2 <- FeaturePlot(x.combined, paste0("adt_",adt), pt.size=1.5) 

    png(paste0("adtplots/umap_", adt, ".png"), 1000, 500)
    print(p1 + p2)
    dev.off()
}



#for each cluster
#look for gene differences bewteen control and wtc
#for each cluster()
#avg_logFC: log fold-chage of the average expression between the two groups. Positive values indicate that the gene is more highly expressed in the first group 
Idents(x.combined) <- "group"
for (clust in unique(x.combined$divij.celltypes)){
    x1 = subset(x = x.combined, subset = divij.celltypes == clust)
    res1 = FindMarkers(x1, ident.1 = "control", ident.2 = "wtc", logfc.threshold=0, min.pct=0)
    #res1 = FindMarkers(x1, ident.1 = "control", ident.2 = "wtc")
    write.csv(file=paste0("gseaClusters/clust_", clust, "_controlVsWtc.csv"), x=res1)
    #write.csv(file=paste0("deClusters/clust_", clust, "_controlVsWtc.csv"), x=res1)
}

#look for diff in cd45.2 across groups
#, features=c("Ms-CD45-2-TotalA")

DefaultAssay(x.combined) <- "adt"
Idents(x.combined) <- "group"
for (clust in unique(x.combined$divij.celltypes)){
    x1 = subset(x = x.combined, subset = divij.celltypes == clust)
    #res1 = FindMarkers(x1, ident.1 = "control", ident.2 = "wtc", logfc.threshold=0, min.pct=0)
    res1 = FindMarkers(x1, ident.1 = "control", ident.2 = "wtc", logfc.threshold=0)
    write.csv(file=paste0("adtDeClusters/clust_", make.names(clust), "_controlVsWtc.csv"), x=res1)
    #write.csv(file=paste0("gseaClusters/clust_", clust, "_controlVsWtc.csv"), x=res1)
    #write.csv(file=paste0("deClusters/clust_", clust, "_controlVsWtc.csv"), x=res1)
}

    res1 = FindMarkers(x.combined, ident.1 = "control", ident.2 = "wtc", logfc.threshold=0)
#awk '{print FILENAME, $0}' *.csv | awk -F, '$6 < 0.05' | grep -e CD45-2






