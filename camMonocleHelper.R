library(monocle3)
library(magrittr)
library(dplyr)

assign.clusters <- function(cds,assignments,column = 2,new.col.name="annotation") {
  assigner <- assignments[,column]
  names(assigner) <- as.character(assignments$Cluster)
  colData(cds)[,new.col.name] <- sapply(X = clusters(cds),FUN = function(x) {
    assigner[as.character(x)]
  })
  return(cds)
}

recluster.cds <- function(cds,cluster.method="leiden") {
  cds <- preprocess_cds(cds = cds,num_dim = 100,verbose = T)
  cds <- reduce_dimension(cds = cds,verbose = T)
  cds <- cluster_cells(cds = cds,verbose = T,cluster_method = cluster.method)
  cds <- learn_graph(cds = cds,verbose = T)
  return(cds)
}

filter.marker.results <- function(marker_test_res,num.to.keep=3,criterion="pseudo_R2") {
  spl <- split(x = marker_test_res,f = list(marker_test_res$cell_group))
  top_specific_markers <- do.call("rbind",lapply(X = spl,FUN = function(df) {
    df10 <- df[order(df[,criterion],decreasing = T),]
    df10[1:num.to.keep,]
  }))
  top_specific_markers_ids <- unique(top_specific_markers$gene_id)
  return(top_specific_markers_ids)
}

add.gene.sum.column <- function(cds,genes.of.interest,col.name="gene.sum",use.raw=F) {
  feats <- rowData(x = cds) %>% data.frame %>% filter(gene_short_name %in% genes.of.interest)
  exprs <- exprs(cds[feats$id,])
  if (use.raw) {
    exprs.norm <- exprs
  } else {
    exprs.norm <- t(apply(X = exprs,MARGIN = 1,FUN = function(x) {
      x / max(x)
    }))
  }
  exprs.sum <- apply(X = exprs.norm,MARGIN = 2,FUN = sum)
  colData(cds)[,col.name] <- exprs.sum
  return(cds)
}

find.closest.row <- function(input.matrix,center.point=c(0,0)) {
  ss <- apply(X = input.matrix,MARGIN = 1,FUN = function(x) {
    return((center.point[1]-x[1])^2 + (center.point[2]-x[2])^2)
  })
  return(which.min(ss))
}

aggregate.expression.by.factor <- function(cds,grouping,do.print=F,drop.zeros=F,pct = 1) {
  group.call <- as.factor(colData(cds)[,grouping])
  names(group.call) <- row.names(colData(cds))
  exprs.subset <- exprs(cds)
  aggm <- matrix(nrow = dim(exprs.subset)[1],ncol = length(levels(group.call)))
  for (i in c(1:dim(exprs.subset)[1])) {
    if (do.print) {
      print(i)
    }
    theVec <- as.vector(exprs.subset[i,])
    names(theVec) <- colnames(exprs.subset)

    if (drop.zeros) {
      aggvec <- tapply(X = theVec,INDEX = group.call,FUN = function(x) {
        x.trunc <- subset(x = x,subset = x>0)
        if (length(x.trunc) > pct * 0.01 * length(x)) {
          return(mean(x.trunc))
        } else {
          return(0)
        }
      })
    } else {
      aggvec <- tapply(X = theVec,INDEX = group.call,FUN = mean)
    }
    aggm[i,] <- aggvec
  }
  row.names(aggm) <- row.names(exprs.subset)
  colnames(aggm) <- levels(group.call)
  return(aggm)
}

map.gene.short.names <- function(cds,gene.short.names) {
  feats <- rowData(cds) %>% data.frame
  feats.subset <- filter(feats,gene_short_name %in% gene.short.names)
  return(feats.subset$id)
} 