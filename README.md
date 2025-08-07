The data analysis starts from the output of a cellranger processing pipeline.

"tryEccite.r"
The 4 samples were individually processed using standard Seurat procedures before being integrated together into a combined dataspace.  After integration, a umap was made from the top 30 PCs, and neighbors/clusters were formed based on a resolution of 0.5 resulting in 29 clusters.

The Seurat datastructure can be loaded from the following link:
https://www.dropbox.com/scl/fi/et0dnoxdwdew7y91me7rt/x.combined.rdata?rlkey=mrii169ekg2t0w6tejmmoje7x&st=k02voi9g&dl=0



"annotateData.r"
For each of these 29 clusters, we looked at which genes were most differentially expressed with regard to the remaining clusters.  These gene lists were used to help the manual identification of the clusters.  To aid this process, we also performed automatic cell type identification based on 6 distinct cell type libraries.  The automatic celltypes served as a starting point for the lengthy manual classification procedure.

"annotateData3.r"
Once manual cell type identification was performed, we merged clusters with similar labels resulting in a reduced set of 14 cell types.  
The manual labelling of the X clusters is available from the following link:
https://www.dropbox.com/scl/fi/jnd6zicmbt55fg0598n67/clusterCelltypes_20230810.xlsx?rlkey=oajlnx7x4exu8qlbhn4t5tepi&st=d62ln503&dl=0

In each of these 14 manually curated celltype clusters, we
looked for genes with the differential expression between the 2 control and 2 WTC samples.  
For running Gene set enrichment analysis, we made sure to retain differential expression values for all genes, not just the significantly expressed ones.

DE Results can be found from this following link:
https://www.dropbox.com/scl/fo/lyt6gr6ypm8xq1wpmfksf/AP5jLrCai9J8gZpzcYGsXco?rlkey=pspidkn10eoy9uu9an74i17jd&st=w96b15r7&dl=0

CD45.2 positive cells were determined by looking at the ADT expression of the Ms-CD45-2-TotalA probe.
A histogram showed two distinct peaks and a threshold of 0.2 based on visual examination was used to determine -/+ status.

The differential expression and GSEA procedures were also run using just the CD45.2 Positive cells.


Determining whether the proportion of cells of a specific cell type was different in the two
experimental conditions was done with a fisher test.  Here we looked at the contingency table of
the number of cells of each celltype vs the remaining celltypes in both groups.  The pvalues for the fisher test on the 14 celltypes were reported as false discovery rate qvalues.
Results can be accessed from the following link:
https://www.dropbox.com/scl/fi/upikbpn181q0y5fqaqru9/cellcounts_20250521.csv?rlkey=98dllvw693y4f9t0zvydneez2&st=dtf224dw&dl=0


