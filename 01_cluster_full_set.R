setwd("/Volumes/SHAQ/2020/02_14_20_metaplasia_single_cell/git")

library(monocle3)
library(dplyr)
source(file = "camMonocleHelper.R")

## Loads the individual cellranger output directories and merges them
# d00.cds <- load_cellranger_data(pipestance_path = "../D00_out/")
# d07.cds <- load_cellranger_data(pipestance_path = "../D07_out/")
# d42.cds <- load_cellranger_data(pipestance_path = "../D42_out/")
# cds <- combine_cds(cds_list = list(d00.cds,d07.cds,d42.cds))

## Normalize the count matrix by dividing by a size factor and take the logarithm
## Then compute PCA and use 100 PCs
## Reduce dimension using UMAP algorithm
## Perform clustering using the leiden algorithm
# cds <- recluster.cds(cds = cds,cluster.method = "leiden")

## Plots the clustering results
# print(plot_cells(cds = cds))

## Computes markers of each cluster using Jensen-Shannon specificities
# mtr <- top_markers(cds = cds,verbose = T)

## Identifies top marker genes for the initial clusters
# tm <- filter.marker.results(marker_test_res = mtr,num.to.keep = 4,criterion = "pseudo_R2")
# print(plot_genes_by_group(cds = cds,markers = tm,group_cells_by = "cluster",max.size = 3))

## After manual curation of cluster types, opens the cluster curation file and applies it to the dataset
## Please note that if you are re-running this code from scratch, your cluster numbers may not be the same
## In this case, you will have to re-annotate these clusters
# assignments <- read.csv(file = "clusters_full.csv",header = T)
# cds <- assign.clusters(cds = cds,assignments = assignments)
# print(plot_cells(cds = cds,color_cells_by = "annotation"))
# save(cds,file = "cds_full.rds")

## Grabs just keratinocytes
# barcodes <- colData(cds) %>% data.frame
# cds.ker <- cds[,subset(barcodes,barcodes$annotation=="Keratinocytes") %>% row.names]

## Converts sample ids from a numeric to a factor
# colData(cds.ker)$sample <- as.factor(colData(cds.ker)$sample)

## Re-clusters the keratinocytes to find TZ-derived cells
# cds.ker <- recluster.cds(cds = cds.ker)
# print(plot_cells(cds = cds.ker,show_trajectory_graph = F,group_label_size = 6))
# mtr.ker <- top_markers(cds = cds.ker,verbose = T)
# tm.ker <- filter.marker.results(marker_test_res = mtr.ker,num.to.keep = 4,criterion = "pseudo_R2")
# print(plot_genes_by_group(cds = cds.ker,markers = tm.ker,max.size = 4))

## Filters out likely doublets between keratinocytes and other cell types, and then reclusters
## cluster 18: fibroblast fusions
## cluster 14: RBCs
## cluster 17: RBCs
## cluster 12: IECs
## cluster 13: myeloid cells
# colData(cds.ker)$cell.group <- clusters(cds.ker)
# barcodes.ker <- colData(cds.ker) %>% data.frame
# cds.ker.fil <- cds.ker[,subset(barcodes.ker,!barcodes.ker$cell.group %in% c(18,14,17,12,13)) %>% row.names]
# cds.ker.fil <- recluster.cds(cds = cds.ker.fil)
# mtr.ker.fil <- top_markers(cds = cds.ker.fil,verbose = T)
# tm.ker.fil <- filter.marker.results(marker_test_res = mtr.ker.fil,num.to.keep = 4,criterion = "pseudo_R2")
# print(plot_genes_by_group(cds = cds.ker.fil,markers = tm.ker.fil,max.size = 4))

## Found additional fusion cells, performs additional filtering and reclusters
## cluster 12: IECs
## cluster 13: fibroblasts
# colData(cds.ker.fil)$cell.group <- clusters(cds.ker.fil)
# barcodes.ker.fil <- colData(cds.ker.fil) %>% data.frame
# cds.ker.fil.fil <- cds.ker.fil[,subset(barcodes.ker.fil,!barcodes.ker.fil$cell.group %in% c(12,13)) %>% row.names]
# cds.ker.fil.fil <- recluster.cds(cds = cds.ker.fil.fil)
# mtr.ker.fil.fil <- top_markers(cds = cds.ker.fil.fil,verbose = T)
# tm.ker.fil.fil <- filter.marker.results(marker_test_res = mtr.ker.fil.fil,num.to.keep = 5,criterion = "marker_score")
# print(plot_genes_by_group(cds = cds.ker.fil.fil,markers = tm.ker.fil.fil,max.size = 4))


## Assigns keratinocyte clusters and filters out the unknown cluster
# assignments.ker.fil.fil <- read.csv(file = "clusters_keratinocytes_fil_fil.csv",header = T)
# cds.ker.fil.fil <- assign.clusters(cds = cds.ker.fil.fil,assignments = assignments.ker.fil.fil,new.col.name = "annotation.2")
# barcodes.ker.fil.fil <- colData(cds.ker.fil.fil) %>% data.frame
# cds.ker.fil3 <- cds.ker.fil.fil[,subset(barcodes.ker.fil.fil,!barcodes.ker.fil.fil$annotation.2 %in% c("Unknown")) %>% row.names]
# cds.ker.fil3 <- recluster.cds(cds = cds.ker.fil3)
# save(cds.ker.fil3,file = "cds_keratinocytes.rds")

## Grabs just the exp d 0 keratinocytes
# barcodes.ker.fil3 <- colData(cds.ker.fil3) %>% data.frame
# cds.ker.day00 <- cds.ker.fil3[,subset(barcodes.ker.fil3,barcodes.ker.fil3$sample == 1) %>% row.names]
# cds.ker.day00 <- recluster.cds(cds = cds.ker.day00)
# cds.ker.day00 <- cluster_cells(cds = cds.ker.day00,cluster.method = "leiden",resolution=0.001)
# mtr.ker.day00 <- top_markers(cds = cds.ker.day00,verbose = T)
# tm.ker.day00 <- filter.marker.results(marker_test_res = mtr.ker.day00,num.to.keep = 5,criterion = "pseudo_R2")
# print(plot_genes_by_group(cds = cds.ker.day00,markers = tm.ker.day00,max.size = 4))
# print(plot_cells(cds = cds.ker.day00,group_label_size = 6,show_trajectory_graph = F))
# colData(cds.ker.day00)$cell.group <- clusters(cds.ker.day00)
# save(cds.ker.day00,file = "cds_keratinocytes_day00.rds")

## Grabs just the TZ-derived keratinocytes (from all timepoints)
# barcodes.ker.fil3 <- colData(cds.ker.fil3) %>% data.frame
# Krt7.id <- map.gene.short.names(cds = cds.ker.fil3,gene.short.names = c("Krt7"))
# barcodes.ker.fil3$Krt7.count <- as.vector(exprs(cds.ker.fil3[Krt7.id,]))
# print(summary(subset(barcodes.ker.fil3$Krt7.count,log10(barcodes.ker.fil3$Krt7.count+0.01) > 0)))
# cds.ker.Krt7 <- cds.ker.fil3[,subset(barcodes.ker.fil3,barcodes.ker.fil3$Krt7.count >= 5) %>% row.names] # median expression cut-off
# cds.ker.Krt7 <- recluster.cds(cds = cds.ker.Krt7)
# save(cds.ker.Krt7,file = "cds_Krt7.rds")

## Grabs just the TZ-derived keratinocytes at exp d 42
# barcodes.Krt7 <- colData(cds.ker.Krt7) %>% data.frame
# cds.ker.Krt7.day42 <- cds.ker.Krt7[,subset(barcodes.Krt7,barcodes.Krt7$sample == "3") %>% row.names]
# cds.ker.Krt7.day42 <- recluster.cds(cds = cds.ker.Krt7.day42)
# save(cds.ker.Krt7.day42,file = "cds_Krt7_day42.rds")

## Grabs the colonic epithelial cells (from all timepoints)
# barcodes <- colData(cds) %>% data.frame
# cds.cec <- cds[,subset(barcodes,grepl(pattern = "Colonocytes",x = barcodes$annotation)) %>% row.names]
# cds.cec <- recluster.cds(cds = cds.cec)
# mtr.cec <- top_markers(cds = cds.cec,verbose = T)
# tm.cec <- filter.marker.results(marker_test_res = mtr.cec,num.to.keep = 5,criterion = "pseudo_R2")
# print(plot_genes_by_group(cds = cds.cec,markers = tm.cec,max.size = 5))
# colData(cds.cec)$cell.group <- clusters(cds.cec)

## Filters out colonic epithelial cell doublets
## cluster 7: fibroblasts
## cluster 5: keratinocytes
# barcodes.cec <- colData(cds.cec) %>% data.frame
# cds.cec.fil <- cds.cec[,subset(barcodes.cec,!barcodes.cec$cell.group %in% c(5,7)) %>% row.names]
# cds.cec.fil <- recluster.cds(cds = cds.cec.fil)
# save(cds.cec.fil,file = "cds_cec.rds")
