library("Seurat")
library("patchwork")
library("clustree")
library("magrittr")
library("dplyr")
library("ggplot2")

#IL1RAP as well as IL1b expression looks increased in TET2 mutant (CD45.2+) derived granulocytes but it is not significant. Perhaps we can look into the cluster and see if the expression is significantly increased in the sub-cluster within granulocytes? (Originally this cluster had around 3 sub-clusters).
#
#Is there any difference in the expression of IL6 or IL6r in either total cell or CD45.2-derived different cell cluster?
#
#Could you provide the list of Differentially Expressed Genes (DEGs) in CD45.2-derived granulocytes, neutrophils, eosinophils, plasma cells, IgD+ B cells, and IgM+ B cells? Additionally, are there any specific Gene Set Enrichment Analysis (GSEA) results associated with these DEGs.
#
#Please share the list of Differentially Expressed Genes (DEGs) in the entire granulocyte population, neutrophils, eosinophils, plasma cells, IgD+ B cells, and IgM+ B cells population. Are there any specific Gene Set Enrichment Analysis (GSEA) findings related to these DEGs.



#S1-S2 represents our control group, while S3-S4 refers to the WTC-treated mice. 



#try automatic annotation with celldex and singler
library("readxl")

load("x.combined.rdata")
cluster.info = data.frame(read_xlsx("clusterCelltypes_20230810.xlsx"))


#set the cluster ids to the new annotated cell types

Idents(x.combined)
id2ctype = list()
id2ctype = cluster.info$Final.Cluster
names(id2ctype) = as.character(cluster.info$cluster)
cluster.info[1:5,]



#x.combined$divij.celltypes = id2ctype[as.character(x.combined$integrated_snn_res.0.5)]
x.combined@meta.data$divij.celltypes = id2ctype[as.character(x.combined$integrated_snn_res.0.5)]




#run DE on all clusters to help with cell type identification.
DefaultAssay(x.combined) <- "integrated"


#S1-S2 represents our control group, while S3-S4 refers to the WTC-treated mice. 
group = NA
group[x.combined$sample %in% c("S1", "S2")] = "control"
group[x.combined$sample %in% c("S3", "S4")] = "wtc"
x.combined$group = group




#for each cluster
#look for gene differences bewteen control and wtc
#for each cluster()
#avg_logFC: log fold-chage of the average expression between the two groups. Positive values indicate that the gene is more highly expressed in the first group 

DefaultAssay(x.combined) <- "RNA"
Idents(x.combined) <- "group"
for (clust in unique(x.combined$divij.celltypes)){
    x1 = subset(x = x.combined, subset = divij.celltypes == clust)
    #res1 = FindMarkers(x1, ident.1 = "control", ident.2 = "wtc", logfc.threshold=0, min.pct=0)
    res1 = FindMarkers(x1, ident.1 = "control", ident.2 = "wtc", logfc.threshold = 0.1, min.pct=0.0)
    #write.csv(file=paste0("gseaClusters_all/clust_", make.names(clust), "_controlVsWtc.csv"), x=res1)
    write.csv(file=paste0("de_all/clust_", make.names(clust), "_controlVsWtc.csv"), x=res1)
}


#look at cd45-2 across all celltypes
DefaultAssay(x.combined) <- "adt"
Idents(x.combined) <- "group"
res.cd452 = lapply(unique(x.combined$divij.celltypes), function(clust){
    x1 = subset(x = x.combined, subset = divij.celltypes == clust)
    res1 = FindMarkers(x1, features=c("Ms-CD45-2-TotalA"), ident.1 = "control", ident.2 = "wtc", logfc.threshold=0)
    res1
})
names(res.cd452) = unique(x.combined$divij.celltypes)
res.cd452 = do.call(rbind, res.cd452)
res.cd452$celltypes = rownames(res.cd452)
res.cd452$neglog10p = -log10(res.cd452$p_val)
res.cd452$neglog10p = -log10(res.cd452$p_val)
#get the sizes of the celltypes
res.cd452 = merge(res.cd452, data.frame(table(x.combined$divij.celltypes)), by.x="celltypes", by.y="Var1")
res.cd452$celltypeSize = res.cd452$Freq/max(res.cd452$Freq, na.rm=T)*5 + 1
sc = scale_colour_gradientn(colours = c("black", "red"), values=c(0, 1.30103)) 
a = -log10(0.05)
b = -log10(0.001)
p1 = ggplot(res.cd452, aes(x=avg_log2FC, y = celltypes, color=neglog10p, size=celltypeSize)) + 
    geom_point() +
    scale_colour_gradientn(
        colours = c("black", "orange", "red"),
        values = scales::rescale(c(0, a, b)),  # rescale to [0, 1]
        limits = c(0, b),  # or c(min_value, max_value) as needed
        oob = scales::squish  # ensures values beyond b still get red
    )+
    scale_x_reverse() +
    xlim(1.,-1.) +
    ggtitle("CD45-2 ADT expression:  Control vs Wtc")
p1

p2 = ggplot(res.cd452, aes(x=avg_log2FC, y = celltypes, color=p_val_adj<0.05, size=celltypeSize)) + 
    geom_point() +
    ggtitle("CD45-2 ADT expression:  Control vs Wtc")

pdf("de_cd45-2_controlVsWtc_v4.pdf")
print(p1)
dev.off()



#look for differences in il1rap and il1b in just the cd45-2 positive cells
rownames(x.combined[["adt"]]) %in% "Ms-CD45-2-TotalA"
class(x.combined)
class(x.combined[["adt"]])
x.combined[["adt"]]@data[1,]
GetAssayData(x.combined, assay = "adt", slot = "data")

png("density_cd45-2.png")
#plot(density(x.combined[["adt"]]["Ms-CD45-2-TotalA",]))
plot(density(x.combined[["adt"]]@data["Ms-CD45-2-TotalA",]))
dev.off()
#maybe pick > 0.2 for cd45-2 positive cells
DefaultAssay(x.combined) <- "adt"
Idents(x.combined) <- "group"
rownames(x.combined[["adt"]]) %in% "Ms-CD45-2-TotalA"

x.combined@meta.data$cd45.2Pos = GetAssayData(x.combined[["adt"]], slot = "data")["Ms-CD45-2-TotalA",] > 0.2
#x.combined@meta.data$cd45.2Pos = t(x.combined[["adt"]]["Ms-CD45-2-TotalA",] > 0.2)

cd45 = GetAssayData(x.combined[["adt"]], slot = "data")["Ms-CD45-TotalA",]
cd452 = GetAssayData(x.combined[["adt"]], slot = "data")["Ms-CD45-2-TotalA",]
png("scatter_cd45_cd45.2.png")
plot(cd45, cd452, col=rgb(0, 0, 0, 0.25))
dev.off()
ix1 = cd452 < 0.2

png("scatter_cd45_cd45.2(under0.2).png")
plot(cd45[ix1], cd452[ix1]+runif(length(cd452[ix1])), col=rgb(0, 0, 0, 0.25))
dev.off()


x.cd452pos = subset(x = x.combined, subset = cd45.2Pos == T)
x.cd452neg = subset(x = x.combined, subset = cd45.2Pos == F)


GetAssayData(x.cd452neg[["adt"]], slot = "data")["Ms-CD45-1-TotalA",]
grep(rownames(GetAssayData(x.cd452neg[["adt"]], slot = "data")), pat="CD45", value=T)
plot(density(x.cd452neg[["adt"]]@data["Ms-CD45-1-TotalA",]))

x1 = subset(x = x.combined, subset = divij.celltypes == "HSPC")
#p1 <- FeaturePlot(x1, features=c("rna_Il1rap", "adt_Ms-CD45-2-TotalA"), blend=T, split.by="group", pt.size=1, order=T) & theme(legend.position = c(0.1,0.2))
p1 <- FeaturePlot(x1, features=c("rna_Il1rap", "adt_Ms-CD45-2-TotalA"), blend=T, split.by="group", pt.size=1, order=T)
pdf(file="featureplot2_il1rap_cd45.2_hspc.pdf", width=15, height=8)
p1
dev.off()



DefaultAssay(x.cd452pos) <- "RNA"
Idents(x.cd452pos) <- "group"
for (clust in unique(x.cd452pos$divij.celltypes)){
    x1 = subset(x = x.cd452pos, subset = divij.celltypes == clust)
    #res1 = FindMarkers(x1, ident.1 = "control", ident.2 = "wtc", logfc.threshold=0, min.pct=0)
    res1 = FindMarkers(x1, ident.1 = "control", ident.2 = "wtc")
    write.csv(file=paste0("de_cd452pos/clust_", make.names(clust), "_controlVsWtc.csv"), x=res1)
    #write.csv(file=paste0("gseaClusters_cd452/clust_", make.names(clust), "_controlVsWtc.csv"), x=res1)
}

for (clust in unique(x.cd452pos$divij.celltypes)){
    x1 = subset(x = x.cd452pos, subset = divij.celltypes == clust)
    #res1 = FindMarkers(x1, ident.1 = "control", ident.2 = "wtc", logfc.threshold=0, min.pct=0)
    res1 = FindMarkers(x1, ident.1 = "control", ident.2 = "wtc", logfc.threshold=0, min.pct = 0)
    write.csv(file=paste0("de_cd452pos/gsea/clust_", make.names(clust), "_controlVsWtc.csv"), x=res1)
    #write.csv(file=paste0("gseaClusters_cd452/clust_", make.names(clust), "_controlVsWtc.csv"), x=res1)
}

DefaultAssay(x.cd452neg) <- "RNA"
Idents(x.cd452neg) <- "group"
for (clust in unique(x.cd452neg$divij.celltypes)){
    x1 = subset(x = x.cd452neg, subset = divij.celltypes == clust)
    #res1 = FindMarkers(x1, ident.1 = "control", ident.2 = "wtc", logfc.threshold=0, min.pct=0)
    res1 = FindMarkers(x1, ident.1 = "control", ident.2 = "wtc")
    write.csv(file=paste0("de_cd452neg/clust_", make.names(clust), "_controlVsWtc.csv"), x=res1)
    #write.csv(file=paste0("gseaClusters_cd452/clust_", make.names(clust), "_controlVsWtc.csv"), x=res1)
}

for (clust in unique(x.cd452neg$divij.celltypes)){
    x1 = subset(x = x.cd452neg, subset = divij.celltypes == clust)
    #res1 = FindMarkers(x1, ident.1 = "control", ident.2 = "wtc", logfc.threshold=0, min.pct=0)
    res1 = FindMarkers(x1, ident.1 = "control", ident.2 = "wtc", logfc.threshold=0, min.pct = 0)
    write.csv(file=paste0("de_cd452neg/gsea/clust_", make.names(clust), "_controlVsWtc.csv"), x=res1)
    #write.csv(file=paste0("gseaClusters_cd452/clust_", make.names(clust), "_controlVsWtc.csv"), x=res1)
}

#IL6, IL6r, Nlrp3 nek7, gata1, pycard, gsdmd and casp1
my.genes = c("Il6", "Il6ra", "Nlrp3", "Nek7", "Gata1", "Pycard", "Gsdmd", "Casp1")
DefaultAssay(x.cd452pos) <- "RNA"
Idents(x.cd452pos) <- "group"
res2 = FindMarkers(x.cd452pos, features=c("Il1rap", "Il1b"), ident.1 = "control", ident.2 = "wtc", logfc.threshold=0)
#do it for all celltypes
res2.cd452 = lapply(unique(x.cd452pos$divij.celltypes), function(clust){
    x1 = subset(x = x.cd452pos, subset = divij.celltypes == clust)
    res1 = FindMarkers(x1, features=my.genes, ident.1 = "control", ident.2 = "wtc", logfc.threshold=0, min.pc=0)
    res1
})
names(res2.cd452) = unique(x.cd452pos$divij.celltypes)

sink("de_cd45.2_specialGenes.txt")
print(res2.cd452)
sink()


#L1RAP as well as IL1b expression looks increased in TET2 mutant (CD45.2+) derived granulocytes but it is not significant. Perhaps we can look into the cluster and see if the expression is significantly increased in the sub-cluster within granulocytes? (Originally this cluster had around 3 sub-clusters).
gran.clusts = names(id2ctype)[id2ctype == "Granulocytes"]
res3.cd452 = lapply(gran.clusts, function(clust){
    x1 = subset(x = x.cd452pos, subset = integrated_snn_res.0.5 == clust)
    res1 = FindMarkers(x1, features=c("Il1rap", "Il1b"), ident.1 = "control", ident.2 = "wtc", logfc.threshold=0, min.pc=0)
    res1
})
names(res3.cd452) = gran.clusts



sink("de_cd45.2_granulocyteClusters.txt")
print(res3.cd452)
sink()











#look for diffrence in number of cells across groups


cellcounts = (t(table(x.combined$group, x.combined$divij.celltypes)))
cell.counts.res = data.frame(t(sapply(1:nrow(cellcounts), function(i){
    a1 = cellcounts[i,1]
    b1 = cellcounts[i,2]
    c1 = sum(cellcounts[-i,1])
    d1 = sum(cellcounts[-i,2])
    res = fisher.test(matrix(c(a1,b1,c1,d1), nrow=2))
    res2 = c(control.this = a1, wtc.this=b1, control.other=c1, wtc.other=d1, pval=res$p.value, res$estimate, res$conf.int)
    names(res2)[7:8] = c("low95", "high95")
    res2
})))
rownames(cell.counts.res) = rownames(cellcounts)
cell.counts.res$qval = p.adjust(cell.counts.res$pval, method="fdr")
write.csv(file="cellcounts_20250521.csv", cell.counts.res)

cell.counts.res$celltypes = rownames(cell.counts.res)
#show cellcounts results in a forest plot

fp <- ggplot(data=cell.counts.res, aes(x=celltypes, y=odds.ratio, ymin=low95, ymax=high95)) +
        geom_pointrange() + 
        geom_hline(yintercept=1, lty=2) +  # add a dotted line at x=1 after flip
        coord_flip() +  # flip coordinates (puts labels on y axis)
        xlab("celltype") + ylab("Fisher's Test of proportions odds ratio(95% CI)") +
        theme_bw()  # use a white background
png("cellcounts_forest.png")
print(fp)
dev.off()



#look at overall ilraf1 expression across all celltypes
#Can we look if there is any change in the expression of overall 
#IL1RAP, IL1b and CD45.2  between vehicle and WTC group, independent of the cluster?

dim(x.combined[["RNA"]])
x.combined$group
x.combined[["RNA"]]["Il1rap",]
x.combined[["RNA"]]["Il1b",]
x.combined[["adt"]]["Ms-CD45-2-TotalA",]
rownames(x.combined)
DefaultAssay(x.combined) <- "RNA"
res1 = FindMarkers(x.combined, features=c("Il1rap", "Il1b"), ident.1 = "control", ident.2 = "wtc", logfc.threshold=0, min.pct=0)
DefaultAssay(x.combined) <- "adt"
res2 = FindMarkers(x.combined, features=c("Ms-CD45-2-TotalA"), ident.1 = "control", ident.2 = "wtc", logfc.threshold=0)

res3 = rbind(res1, res2)
write.csv(file="deGlobal_3gene_controlVsWtc.csv", res3[,-5])

#show side by side plot of the dots for control/wtc
p1 <- FeaturePlot(x.combined, "rna_Il1rap", split.by="group", pt.size=3) 
p2 <- FeaturePlot(x.combined, "rna_Il1b", split.by="group", pt.size=3) 
p3 <- FeaturePlot(x.combined, "adt_Ms-CD45-2-TotalA", split.by="group", pt.size=3) 
png("plot_3gene_controlVsWtc_v3.png", 1000, 1500)
p1/p2/p3
dev.off()

#maybe show histogram of the significant diff of cd45-2 values?

vals.cd45 = x.combined[["adt"]]@data["Ms-CD45-2-TotalA",]
df1 = dataframe(vals.cd45, group = x.combined$group)

wilcox.test(vals.cd45 ~ group, data=df1)
ggplot(df1, aes(x=vals.cd45, color=group)) +
    geom_density()
    #geom_histogram(position="dodge")

#show a scatterplot of cd45-2 and il1rap
#cd45-2 il1b

p1 = FeatureScatter(object = x.combined, feature2 = 'Il1rap', feature1 = 'Ms-CD45-2-TotalA')
p1
p1 + facet_wrap(vars(celltypes))

library("gridExtra")
ps = lapply(unique(x.combined$divij.celltypes), function(clust){
        x1 = subset(x = x.combined, subset = divij.celltypes == clust)

        p1 = FeatureScatter(object = x1, feature2 = 'Il1rap', feature1 = 'Ms-CD45-2-TotalA', plot.cor=F)
        p1 + ggtitle(clust) 
})
pdf("scatter_cd45-2_il1rap.pdf", 20, 20)
do.call("grid.arrange", c(ps, ncol=4))
dev.off()

ps = lapply(unique(x.combined$divij.celltypes), function(clust){
        x1 = subset(x = x.combined, subset = divij.celltypes == clust)

        p1 = FeatureScatter(object = x1, feature2 = 'rna_Il1b', feature1 = 'Ms-CD45-2-TotalA', plot.cor=F)
        p1 + ggtitle(clust) 
})
pdf("scatter_cd45-2_il1b.pdf", 20, 20)
do.call("grid.arrange", c(ps, ncol=4))
dev.off()








res1
res2
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



#make a function that does this
#plot histogram gene expression of the 4 samples over each celltype.

grep(rownames(x.combined[["RNA"]]), pat="cd45", ignore.case=T, value=T)

DefaultAssay(x.combined) <- "RNA"
gene = "Cxcr4"

p1 <- FeaturePlot(x.combined, "rna_Cxcr4", split.by="sample") 

p1 <- FeaturePlot(x.combined, "adt_Ms-CD45-TotalA", split.by="sample") 
p2 <- FeaturePlot(x.combined, "adt_Ms-CD45-2-TotalA") 
p1+p2


library("ggridges")
cd45.1 = as.numeric(x.combined[['adt']]["Ms-CD45-TotalA",])
cd45.2 = as.numeric(x.combined[['adt']]["Ms-CD45-2-TotalA",])
tet2 = as.numeric(x.combined[['RNA']]["Tet2",])
class(cd45.1)
#plot the distribution of cd451/2 proteins by sample and celltype
df1 = data.frame(sample = x.combined$sample, x.combined$group, cd45.1, cd45.2, tet2)
df1 = cbind(x.combined@meta.data[,mouse.refs], df1)
head(df1)
p1 = ggplot(df1, aes(x=cd45.1, y=sample, fill=group)) + geom_density_ridges() + facet_wrap(. ~ ref_immgen)
p2 = ggplot(df1, aes(x=cd45.2, y=sample, fill=group)) + geom_density_ridges() + facet_wrap(. ~ ref_immgen)
p3 = ggplot(df1, aes(x=tet2, y=sample, fill=group)) + geom_density_ridges() + facet_wrap(. ~ ref_immgen)
pdf("density_adtCD45.1.pdf", width=10, height=10)
p1
dev.off()
ix = grepl(df1$ref_immgen, pat="Macro")
table(ix)
by(df1$cd45.2[ix], df1$sample[ix], FUN=function(a){table(a>1)})
by(df1$cd45.2[ix], df1$sample[ix], FUN=. %>% sum(. > 0))
sum(df1$cd45.2 > 1)

pdf("density_adtCD45.2.pdf", width=10, height=10)
p2
dev.off()

pdf("density_rnaTet2.pdf", width=10, height=10)
p3
dev.off()


dim(df1)
dim(df1)




plot(density(x.combined[['adt']]["Ms-CD45-TotalA",] ))
plot(density(x.combined[['adt']]["Ms-CD45-2-TotalA",] ))

x.combined[['adt']]["Ms-CD45-TotalA",] / x.combined[['adt']]["Ms-CD45-2-TotalA",]
eps = 0.00000001
rat1 = (x.combined[['adt']]["Ms-CD45-2-TotalA",]+eps) / (x.combined[['adt']]["Ms-CD45-TotalA",]+eps)
plot(density(log10(rat1)))
min(rat1, na.rm=T)
max(rat1, na.rm=T)
sort(table(rat1))
p1
grep(rownames(x.combined[["RNA"]]), pat="cxcr4", ignore.case=T, value=T)
grep(rownames(x.combined[["RNA"]]), pat="CD45", ignore.case=T, value=T)

grep(rownames(x.combined[["adt"]]), pat="cd45", ignore.case=T, value=T)
grep(rownames(x.combined[["adt"]]), pat="ly", ignore.case=T, value=T)




#check the counts between samples for each celltype
celltypes = x.combined$celltypes
sample = x.combined$sample
df = data.frame(celltypes, sample)
library("dplyr")
count(df, celltypes, sample)

tab1 = table(celltypes, sample)
mat1 = as.data.frame(tab1)
class(mat1)
dim(tab1)
tab1[,2]


#prepare output matrix
res.mat = as.data.frame(matrix(NA, nrow=nrow(tab1), ncol=6))
colnames(res.mat) = c("nCtrl.thisCluster", "nStat3c.thisCluster", "nCtrl.otherClusters", "nStat3c.otherClusters", "OR", "pval")
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


################################
#set the celltypes based on MCA
ref = mouse.refs[1]


f1 = get(ref)
ref.mat = f1()


res <- clustify(
  input = int.mat, # matrix of normalized scRNA-seq counts (or SCE/Seurat object)
  metadata = int.meta, # meta.data table containing cell clusters
  cluster_col = "integrated_snn_res.0.3", # name of column in meta.data containing cell clusters
  ref_mat = ref.mat, # matrix of RNA-seq expression data for each cell type
  query_genes = rownames(x.combined[["integrated"]]) # list of highly varible genes identified with Seurat
)


res2 <- cor_to_call(
  cor_mat = res                  # matrix correlation coefficients
)

colnames(res2)[2] = ref
clust.res = merge(clust.res, res2[,1:2])


clust2celltype = list()
for (i in 1:nrow(res2)){
    clust2celltype[[unlist(res2[i,1])]] = unlist(res2[i,2])
}

x.combined$celltypes = unlist(clust2celltype[as.character(x.combined@meta.data$integrated_snn_res.0.3)])

x.combined <- SetIdent(x.combined, value = "celltypes")

g = DimPlot(x.combined, split.by="sample", label = TRUE, repel = TRUE)
g
#split by sample

















res <- clustify(
  input = int.mat, # matrix of normalized scRNA-seq counts (or SCE/Seurat object)
  metadata = int.meta, # meta.data table containing cell clusters
  cluster_col = "integrated_snn_res.0.3", # name of column in meta.data containing cell clusters
  ref_mat = ref.mat, # matrix of RNA-seq expression data for each cell type
  query_genes = rownames(x.combined[["integrated"]]) # list of highly varible genes identified with Seurat
)

res[1:5, 1:5]

res2 <- cor_to_call(
  cor_mat = res                  # matrix correlation coefficients
)

clust2celltype = list()
for (i in 1:nrow(res2)){
    clust2celltype[[unlist(res2[i,1])]] = unlist(res2[i,2])
}

x.combined$celltypes = unlist(clust2celltype[x.combined@meta.data$integrated_snn_res.0.3])

x.combined <- SetIdent(x.combined, value = "celltypes")
png(paste0("celltypes_", ref, ".png"))
g = DimPlot(x.combined, label = TRUE, repel = TRUE)
print(g)
dev.off()

#make a plot of annotated clusters


x.combined <- SetIdent(x.combined, value = "integrated_snn_res.0.3")
png("umap_res0.3_clusters.png")

