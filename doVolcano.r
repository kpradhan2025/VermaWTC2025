library("EnhancedVolcano")



mygenes = c(
    "Il1rap",
    "Hspa8",
    "Hspa5",
    "Cd52",
    "Junb",
    "Itgb2",
    "Lilr4b",
    "Chka",
    "Cxcr2",
    "Gata3",
    "Cd74"
)

#load in the cd45.2+ hspc DE results
x = read.csv("gseaClusters_cd452/clust_HSPC_controlVsWtc.csv")
rownames(x) = x$X

pdf("volcano_cd452_hspc_4.pdf")
EnhancedVolcano(x,
    lab = rownames(x),
    xlim = c(-2, 2),
    ylim = c(0, 7),
    pointSize=1,
    FCcutoff = 0.4,
    pCutoff = NA,
    x = 'avg_log2FC',
    y = 'p_val',
    legendPosition = 'none',
    selectLab = mygenes,
    drawConnectors = TRUE,
    labSize = 6.0,
    labCol = 'black',
    labFace = 'bold',
    colConnectors = 'red',
    boxedLabels = TRUE) + scale_x_reverse()
dev.off()

