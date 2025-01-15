
library("fgsea")
library("ggpubr")


#cfiles = list.files(path="gseaClusters_all", full=T)
cfiles = list.files(path="gseaClusters_cd452", full=T)
f = cfiles[1]
x = read.csv(f)


require("biomaRt")
human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/") 
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x[,1] , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=F)

head(genesV2)
m2h = genesV2[,2]
names(m2h) = genesV2[,1]
mouse.genes = m2h[x[,1]]



gmt.files = list.files(path="/home/kpradhan/Desktop/projects/mdsigdb/", full=T)

gmt.names = gsub(gsub(gmt.files, pat=".*\\/\\/", rep=""), pat=".symbols.gmt", rep="")
gmt.pathways = lapply(gmt.files, gmtPathways)

length(gmt.pathways)

pathways <- gmtPathways("/home/kpradhan/Desktop/projects/mdsigdb/msigdb.v7.4.symbols.gmt")


makeGlist <- function(f){
    x = read.csv(f)
    x$genes = m2h[x[,1]]
    x = x[x$genes != "",]
    head(x)
    glist = x$avg_log2FC
    names(glist) = x$genes

    #remove NA values
    glist = glist[!is.na(glist)]

    glist2 = sapply(unique(names(glist)), function(g){
        vals = glist[names(glist) == g]
        ix = which.max(abs(vals))
        vals[ix]
    })
    names(glist2) = unique(names(glist))

    #when there is more than 1 entry for a gene, retain the one with the highest absolute value


    sort(glist2)
}

g2 = makeGlist(f)



makeGseaTable <- function(glist, paths, fdr.thresh = .20){

    fgRes <- fgsea::fgsea(pathways = paths, 
                               stats = glist,
                               minSize=5,
                               maxSize=600
                               ) %>% 
                      as.data.frame() 
    #sort by nes
    fgRes = fgRes[order(fgRes$NES),]

    #only report rows with sig fdr
    res.sig = fgRes[fgRes$padj < fdr.thresh,]

    #take out null columns
    res.sig = res.sig[!is.na(res.sig$NES),]
    #p1 = as_ggplot(plotGseaTable(paths, 
    #                  glist, 
    #                  res.sig, 
    #                  gseaParam = 0.5,
    #                  render=F))
    #list(res.sig, p1)
    res.sig
}


x = read.csv(de.file, sep="\t")

#get a list of genes up and down regulated
#padj < .05 and l2fc < -2
#padj < .05 and l2fc > 2
saveList <- function(file, x){
    x = x[!is.na(x)]
    x = x[x!=""]
    x
    write.csv(file=file, x=x, quote=F, row.names=F, col.names=F)
}
saveList(file = "sig_l2fc.lt.n2_padj05.csv", x= x$hgnc_symbol[which(x$log2FoldChange < -2 & x$padj < 0.05)])
saveList(file = "sig_l2fc.gt.2_padj05.csv", x= x$hgnc_symbol[which(x$log2FoldChange > 2 & x$padj < 0.05)])
which(res$log2FoldChange > 2 & res$padj < 0.05)

#see what kind of genes (functionally) are differentially expressed in our 
#neutrophil, granulocyte, erythroid, and precursor erythroid cluster between control and WTC?

for (de.file in cfiles){
    for (i in seq_along(gmt.pathways)){
        pname = gmt.names[i]
        paths = gmt.pathways[[i]]
        glist = makeGlist(de.file)
        celltype = make.names(gsub(de.file, pat=".*clust_(.*)_control.*", rep="\\1"))

        print("*****")
        print(de.file)
        print(pname)
        
        outfile = paste0("gsea_cd452/gsea_", pname, "_", celltype, ".txt")

        
        res = makeGseaTable(glist, paths, .05)
        res[,8] = sapply(res[,8], paste0, collapse=";")
        write.table(file = outfile, x=res, row.names=F, sep="\t")
    }
}





#$ grep */*.txt -e inflam -i | grep -e Neu | sed s/.*:// | sed 's/"\([^"]*\).*/\1/' | sort | uniq
#FULCHER_INFLAMMATORY_RESPONSE_LECTIN_VS_LPS_DN
#LAKE_ADULT_KIDNEY_C5_PROXIMAL_TUBULE_EPITHELIAL_CELLS_STRESS_INFLAM
#do GSEA plots for these pathways on the epethelial cells


sig.inf.paths = c(
    "FULCHER_INFLAMMATORY_RESPONSE_LECTIN_VS_LPS_DN",
    "LAKE_ADULT_KIDNEY_C5_PROXIMAL_TUBULE_EPITHELIAL_CELLS_STRESS_INFLAM"
)


#neutrophils
de.file = cfiles[8]
glist = makeGlist(de.file)
celltype = make.names(gsub(de.file, pat=".*clust_(.*)_control.*", rep="\\1"))
for (pname in sig.inf.paths){
    pathway = gmt.pathways[[10]][[pname]]

    gfile = paste0("gsea_neutrophils_", pname, ".png")
    png(gfile)
    p1 = plotEnrichment(pathway, glist) + labs(title=pname) 
    print(p1)
    dev.off()
}

#show gsea for all inflamation pathways

inf.pathways = grep(names(gmt.pathways[[10]]), pat="INFLAM", ignore.case=T, value=T)

fgseaRes <- fgsea(pathways = gmt.pathways[[10]][inf.pathways], stats    = glist)

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))

p2 = ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Inflammatory pathways NES from GSEA (Neutrophils)") + 
  theme_minimal()

pdf("gsea_neutrophils_inflam.pdf", width=20, height=20)
print(p2)
dev.off()

pname = "FULCHER_INFLAMMATORY_RESPONSE_LECTIN_VS_LPS_DN"
pathway = gmt.pathways[[10]][[pname]]
glist = makeGlist(de.file)
celltype = make.names(gsub(de.file, pat=".*clust_(.*)_control.*", rep="\\1"))

fgseaRes <- fgsea(pathways = pathway, stats    = glist)
plotEnrichment(pathway, glist) + labs(title=pname) 

res = makeGseaTable(glist, paths, .05)
        res[,8] = sapply(res[,8], paste0, collapse=";")





