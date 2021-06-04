setwd("/Volumes/SHAQ/2020/02_14_20_metaplasia_single_cell/git/")

library(monocle3)
library(ggplot2)
library(ggthemes)
library(reshape2)
library(ggrepel)
source("camMonocleHelper.R")

## Loads the relevant keratinocyte and colonocyte monocle3 objects
# load(file = "cds_keratinocytes.rds",verbose = T)
# load(file = "cds_cec.rds",verbose = T)

## Filters keratinocyte data to just look at anal epidermis at baseline conditions
# barcodes.ker <- data.frame(colData(cds.ker.fil3))
# cds.anus <- cds.ker.fil3[,subset(barcodes.ker,barcodes.ker$sample==1 & barcodes.ker$annotation.2 %in% c("IFE - Differentiated","IFE - Progenitor")) %>% row.names]
# colData(cds.anus)$tissue <- "anus"

## Prepares the colonocyte data to just look at baseline conditions
# barcodes.cec <- data.frame(colData(cds.cec.fil))
# cds.colon <- cds.cec.fil[,subset(barcodes.cec,barcodes.cec$sample==1) %>% row.names]
# colData(cds.colon)$tissue <- "colon"

## Performs simple differential expression testing
# cds.diff <- combine_cds(cds_list = list(cds.anus,cds.colon))
# agg <- aggregate.expression.by.factor(cds = cds.diff,grouping = "tissue",do.print = T)
# d <- data.frame(agg)
# d$id <- row.names(d)
# feats <- rowData(cds.diff) %>% data.frame
# d10 <- merge(x = feats,y = d,by = "id",all = F)
# d10$log2ratio <- log2(d10$colon / d10$anus)
# d20 <- filter(d10,colon + anus > 0.1)
# d30 <- filter(d20,abs(log2ratio) > 2)
# write.csv(x = d30,file = "table_colon_anus_markers.csv")