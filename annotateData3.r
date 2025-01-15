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

library("readxl")

load("x.combined.rdata")


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
#get the sizes of the celltypes
res.cd452 = merge(res.cd452, data.frame(table(x.combined$divij.celltypes)), by.x="celltypes", by.y="Var1")
res.cd452$celltypeSize = res.cd452$Freq/max(res.cd452$Freq, na.rm=T)*5 + 1
sc = scale_colour_gradientn(colours = c("black", "red"), values=c(0, 1.30103)) 
p1 = ggplot(res.cd452, aes(x=avg_log2FC, y = celltypes, color=neglog10p, size=celltypeSize)) + 
    geom_point() +
    scale_colour_gradientn(colours = c("black", "black", "red", "yellow"), values=c(0, 1.30103, 3, 8)/8) + 
    xlim(-0.15, 0.15) +
    ggtitle("CD45-2 ADT expression:  Control vs Wtc")

p2 = ggplot(res.cd452, aes(x=avg_log2FC, y = celltypes, color=p_val_adj<0.05, size=celltypeSize)) + 
    geom_point() +
    ggtitle("CD45-2 ADT expression:  Control vs Wtc")

png("de_cd45-2_controlVsWtc.png")
print(p1)
dev.off()
png("de_cd45-2_controlVsWtc_v2.png")
print(p2)
dev.off()



#look for differences in il1rap and il1b in just the cd45-2 positive cells

png("density_cd45-2.png")
plot(density(x.combined[["adt"]]["Ms-CD45-2-TotalA",]))
dev.off()
#maybe pick > 0.2 for cd45-2 positive cells
x.combined@meta.data$cd45.2Pos = t(x.combined[["adt"]]["Ms-CD45-2-TotalA",] > 0.2)

x.cd452pos = subset(x = x.combined, subset = cd45.2Pos == T)



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
    write.csv(file=paste0("de_cd452/clust_", make.names(clust), "_controlVsWtc.csv"), x=res1)
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
for (i in 1:nrow(cellcounts)){
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
write.csv(file="cellcounts.csv", cell.counts.res)

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



