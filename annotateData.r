library("Seurat")
library("patchwork")
library("clustree")
library("magrittr")
library("dplyr")
library("ggplot2")


#S1-S2 represents our control group, while S3-S4 refers to the WTC-treated mice. 



#try automatic annotation with celldex and singler
library("celldex")
library("SingleR")
library("SingleCellExperiment")
library("clustifyr")

load("x.combined.rdata")


#make featuresplots for all proteints
x.combined <- SetIdent(x.combined, value = "integrated_snn_res.0.5")
for (adt in rownames(x.combined[["adt"]])){
    print(adt)
    p1 <- DimPlot(x.combined, reduction = "umap", label = TRUE, repel = TRUE)
    p2 <- FeaturePlot(x.combined, paste0("adt_",adt), pt.size=1.5) 

    png(paste0("adtplots/umap_", adt, ".png"), 1000, 500)
    print(p1 + p2)
    dev.off()
}



#run DE on all clusters to help with cell type identification.
DefaultAssay(x.combined) <- "integrated"

x.combined <- SetIdent(x.combined, value = "integrated_snn_res.0.5")
all.markers = FindAllMarkers(x.combined)

write.csv(file="all.markers.csv", x=all.markers)



#S1-S2 represents our control group, while S3-S4 refers to the WTC-treated mice. 
group = NA
group[x.combined$sample %in% c("S1", "S2")] = "control"
group[x.combined$sample %in% c("S3", "S4")] = "wtc"
x.combined$group = group


#for each cluster
#look for gene differences bewteen control and wtc
#for each cluster()
#avg_logFC: log fold-chage of the average expression between the two groups. Positive values indicate that the gene is more highly expressed in the first group 
Idents(x.combined) <- "group"
for (clust in unique(x.combined$integrated_snn_res.0.5)){
    x1 = subset(x = x.combined, subset = integrated_snn_res.0.5 == clust)
    res1 = FindMarkers(x1, ident.1 = "control", ident.2 = "wtc")
    write.csv(file=paste0("deClusters/clust_", clust, "_controlVsWtc.csv"), x=res1)

}

#https://rnabioco.github.io/clustifyrdata/
#ref_immgen

library("clustifyrdatahub")
library("Matrix")

mouse.refs = c(
    "ref_MCA",
    "ref_tabula_muris_drop",
    "ref_tabula_muris_facs",
    "ref_mouse.rnaseq",
    "ref_moca_main",
    "ref_immgen"
)


#integrated features
int.mat = x.combined[["RNA"]]@data[]
int.meta = x.combined@meta.data

#test
ref.immgen = ref_immgen()

runFisherTests <- function(tab1, labels = c("Ctrl", "Wtc")){
    #prepare output matrix
    res.mat = as.data.frame(matrix(NA, nrow=nrow(tab1), ncol=6))
    l1 = paste0("n", labels[1], ".thisCluster")
    l2 = paste0("n", labels[2], ".thisCluster")
    l3 = paste0("n", labels[1], ".otherClusters")
    l4 = paste0("n", labels[2], ".otherClusters")
    colnames(res.mat) = c(l1, l2, l3, l4, "OR", "pval")
    rownames(res.mat) = rownames(tab1)
    for (i in 1:nrow(tab1)){
        a1 = tab1[i,1]
        b1 = tab1[i,2]
        c1 = sum(tab1[-i,1])
        d1 = sum(tab1[-i,2])
        mat1 = matrix(c(a1, b1, c1, d1), ncol=2)
        res = fisher.test(mat1)
        res$estimate
        res$p.value
        res.mat[i,] = c(a1, b1, c1, d1, res$estimate, res$p.value) 
    }
    res.mat$qval = p.adjust(res.mat$pval, method="fdr")
    res.mat
}

reso = "integrated_snn_res.0.5"
int.meta[,reso]
clust.res = data.frame(cluster=sort(as.character(unique(int.meta[,reso]))))
for (ref in mouse.refs){
    f1 = get(ref)
    ref.mat = f1()


    res <- clustify(
      input = int.mat, # matrix of normalized scRNA-seq counts (or SCE/Seurat object)
      metadata = int.meta, # meta.data table containing cell clusters
      cluster_col = reso, # name of column in meta.data containing cell clusters
      ref_mat = ref.mat, # matrix of RNA-seq expression data for each cell type
      query_genes = rownames(x.combined[["integrated"]]) # list of highly varible genes identified with Seurat
    )


    res2 <- data.frame(cor_to_call(
      cor_mat = res                  # matrix correlation coefficients
    ))

    res2 = data.frame(res2)
    colnames(res2)[2] = ref
    clust.res = merge(clust.res, res2[,1:2], by="cluster")


    clust2celltype = list()
    for (i in 1:nrow(res2)){
        clust2celltype[[res2[i,1]]] = res2[i,2]
    }

    #save the reference info in the metadata
    x.combined$celltypes = unlist(clust2celltype[as.character(int.meta[,reso])])
    x.combined@meta.data[,ref] = x.combined$celltypes


    x.combined <- SetIdent(x.combined, value = "celltypes")
    pdf(paste0("celltypes_", ref, ".pdf"), width=15, height=10)
    g = DimPlot(x.combined, label = TRUE, repel = TRUE)
    print(g)
    dev.off()
    pdf(paste0("celltypes2_", ref, ".pdf"), width=25, height=10)
    g = DimPlot(x.combined, split.by="sample", label = TRUE, repel = TRUE)
    print(g)
    dev.off()

    celltypes = x.combined$celltypes
    sample = x.combined$sample
    sample2 = NA
    sample2[x.combined$sample %in% c("S1", "S2")] = "control"
    sample2[x.combined$sample %in% c("S3", "S4")] = "wtc"
    tab1 = table(celltypes, sample)
    tab2 = table(celltypes, sample2)
    fisher.res = runFisherTests(tab2)
    write.csv(file=paste0("fisherCelltypes_", ref, ".csv"), x=fisher.res)

}
write.csv(file="clusterCelltypes.csv", x=clust.res, row.names=F)


#write main data again , this time it's filled with annotatino info

save(file="x.combined.rdata", "x.combined")



