Do_pHeatmapSamples<-function(data,method='pearson',hmeth='average',
                             title='Clustering Samples',
                             breaks=seq(-1,1,.1),ann.sample,annotation_colors,
                             sample_labels=rownames(data),...){
  cormat<-cor(data, method = method, use = "p")
  pheatmap_msf(cormat,scale="none",cluster_cols=T, cluster_rows=T,
               border_color=FALSE, main=title,
               clustering_method=hmeth,clustering_distance_rows="euclidean",
               clustering_distance_cols ="euclidean",
               color=colorRampPalette(c("navy", "white","firebrick3"))(length(breaks)),breaks=breaks,
               show_colnames=F,show_rownames=T,
               annotation_row=ann.sample,annotation_col=ann.sample,
               annotation_colors=annotation_colors,
               labels_row=sample_labels, labels_col=sample_labels,
               ...)
}
