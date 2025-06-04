library("EnhancedVolcano")


x = read.csv("./de_all/clust_HSPC_controlVsWtc.csv")
head(x)


 #Dynll1, Dnaja1 and Hspa8 along with other top significant genes. 

ix = which(x$X %in% c("Dynll1", "Dnaja1", "Hspa8"))


plot(x$avg_log2FC, -log10(x$p_val_adj))


g1 = EnhancedVolcano(x,
    lab = x$X,
    x = 'avg_log2FC',
    y = 'p_val',
    colAlpha = 1/5,
    pCutoff = 0.0000268504697259559,
    FCcutoff = 0.1,
    boxedLabels=T,
    drawConnectors=T)

pdf("volc_hspc_ctrlVsWtc_v2.pdf")
print(g1)
dev.off()
