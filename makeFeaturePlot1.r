library("Seurat")
library("readxl")
library("ggplot2")



load("./x.combined.rdata")

head(x.combined@meta.data)



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







p1 <- FeaturePlot(x.combined, features=c("adt_Ms-CD45-2-TotalA"), split.by="group", pt.size=1, order=T)
p1

feature <- "adt_Ms-CD45-2-TotalA"
data <- FetchData(x.combined, vars = c("UMAP_1", "UMAP_2", feature, "group"))
#plot the dots in random order
data=data[sample(1:nrow(data), size=nrow(data)),]


p2 = ggplot(data, aes(x = UMAP_1, y = UMAP_2, alpha=get(feature), color = group)) +
  geom_point(size = 1) +
  theme_minimal() +
  labs(title = feature)
p2

pdf(paste0("feature_", feature, "_ctrlVsWtc_v1.pdf"), 20, 10)
p1
dev.off()

pdf(paste0("feature_", feature, "_ctrlVsWtc_v2.pdf"), 10, 10)
p2
dev.off()



