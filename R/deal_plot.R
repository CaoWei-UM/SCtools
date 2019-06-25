#' A plot Function
#'
#' This function allows you to plot mutations.
#' @param
#' @keywords plot
#' @export
#' @examples
#' PlotCluster(mutation,CTC,'DDX11L5',reduction = 'tsne',threshold=1)
#' PlotTrajectory(mutation,CTC_monocle,'DDX11L5',threshold=1)

PlotCluster<-function(mutation,seurat,gene_name,layer_name='default',reduction = "umap",threshold=3){
  #get reduction point
  if(round(threshold)!=threshold|threshold<1){stop('threshold must be a positive integer')}
  embeddings.use <- GetDimReduction(object = seurat, reduction.type = reduction,
                                    slot = "cell.embeddings")[,1:2]
  if (length(x = embeddings.use) == 0) {
    stop(paste(reduction.use, "has not been run for this object yet."))
  }
  valid_cells<-intersect(mutation$col.attrs$cell_names[],seurat@cell.names)
  embeddings.use<-embeddings.use[valid_cells,]
  #get mutation sites of the gene
  gene_sites<-SitesInGene(mutation,gene_name)
  sifit<-mutation$layers[[layer_name]][,gene_sites,drop=F]
  #remove cells not in intersect
  rownames(sifit)<-mutation$col.attrs$cell_names[]
  sifit<-sifit[valid_cells,,drop=F]
  #remove rare cell sites
  valid_gene_sites<-which(apply(sifit,2,function(x) length(which(x!=3)))>=threshold)
  if (length(valid_gene_sites) == 0) {
    stop(paste(gene_name, "has no valid mutation in mutation list."))
  }
  sifit_all<-data.frame(sifit[,valid_gene_sites,drop=F])
  sifit_all<- replace(sifit_all,sifit_all==3,NA)
  sifit_all<- replace(sifit_all,sifit_all==1,'alert')
  sifit_all<- replace(sifit_all,sifit_all==0,'ref')
  #get site number and site info
  site_info<-mutation$get.attribute.df(MARGIN=1,
             attribute.names = c('chromosome','site','ref_allele','alt_allele'))[gene_sites[valid_gene_sites],]
  site_number=length(valid_gene_sites)
  #start drawing
  plot_list=list()
  for(i in 1:site_number){
    mutation<-sifit_all[,i]
    size<- replace(mutation,mutation %in% c('ref','alert'),'big')
    size<- replace(size,is.na(mutation),'small')
    mutation<-as.data.frame(cbind(embeddings.use,mutation,size))
    colnames(mutation)<-c('Dim1','Dim2','mutation','size')
    mutation$Dim1<-as.numeric(mutation$Dim1)
    mutation$Dim2<-as.numeric(mutation$Dim2)
    mutation$size<-as.factor(mutation$size)
    plot<-ggplot(data=mutation)+
      geom_point(aes(Dim1, Dim2,colour=mutation,size=size))+
      ggtitle(paste(site_info[i,],collapse = '-'))+
      scale_size_manual(values=c(1,0.5))+
      scale_colour_discrete(na.value = "grey80")
    annotation_data<-plot$data
    plot$data<-annotation_data[rev(order(annotation_data$mutation)),]
    plot_list[[i]]<-plot
  }
  pdf('temp.pdf',width=6*round(sqrt(site_number)),height=6*ceiling(site_number/round(sqrt(site_number))))
  multiplot(plotlist=plot_list,cols=round(sqrt(site_number)))
  dev.off()
  # pdf_convert('temp.pdf', format = "png", pages = NULL,
  #         filenames = paste(gene_name,'_mutation_',reduction,'.png',sep=''));
}
#' @export
PlotTrajectory<-function(mutation,monocle,gene_name,layer_name='default',threshold=3){
  #get reduction point
  if(round(threshold)!=threshold|threshold<1){stop('threshold must be a positive integer')}
  valid_cells<-intersect(mutation$col.attrs$cell_names[],rownames(monocle@phenoData@data))
  #get mutation sites of the gene
  gene_sites<-SitesInGene(mutation,gene_name)
  sifit<-mutation$layers[[layer_name]][,gene_sites,drop=F]
  #remove cells not in intersect
  rownames(sifit)<-mutation$col.attrs$cell_names[]
  sifit_all<-as.data.frame(matrix(NA,nrow=length(rownames(monocle@phenoData@data)),
                                  ncol=dim(sifit)[2]))
  rownames(sifit_all)<-rownames(monocle@phenoData@data)
  sifit_all[valid_cells,]<-sifit[valid_cells,,drop=F]
  #remove rare cell sites
  valid_gene_sites<-which(apply(sifit_all,2,function(x) length(which(x!=3)))>=threshold)
  if (length(valid_gene_sites) == 0) {
    stop(paste(gene_name, "has no valid mutation in mutation list."))
  }
  sifit_all<-data.frame(sifit_all[,valid_gene_sites,drop=F])
  sifit_all<- replace(sifit_all,sifit_all==3,NA)
  sifit_all<- replace(sifit_all,sifit_all==1,'alert')
  sifit_all<- replace(sifit_all,sifit_all==0,'ref')
  #get site number and site info
  site_info<-mutation$get.attribute.df(MARGIN=1,
                                       attribute.names = c('chromosome','site','ref_allele','alt_allele'))[gene_sites[valid_gene_sites],]
  site_number=length(valid_gene_sites)
  #start drawing
  plot_list=list()
  for(i in 1:site_number){
    monocle@phenoData@data$mutation<-sifit_all[,i]
    plot<-plot_cell_trajectory(monocle,color_by = "mutation",
                               show_cell_names = FALSE)
    annotation_data<-plot$data
    plot$data<-annotation_data[rev(order(annotation_data$mutation)),]
    plot_list[[i]]<-plot
  }
  pdf('temp.pdf',width=6*round(sqrt(site_number)),height=6*ceiling(site_number/round(sqrt(site_number))))
  multiplot(plotlist=plot_list,cols=round(sqrt(site_number)))
  dev.off()
  #pdf_convert('temp.pdf', format = "png", pages = NULL,
  #            filenames = paste(gene_name,'_mutation_trajectory.png',sep=''));
}
#' @export
plot_gene_feature<-function(gene_list,seurat){
  gene_list=intersect(gene_list,rownames(seurat@data))
  pdf('temp.pdf',width = 10,height = 10)
  FeaturePlot(object = seurat,
              features.plot = gene_list,
              cols.use = c("grey", "red"), reduction.use = "umap")
  dev.off()
  #pdf_convert('temp.pdf', format = "png", pages = NULL,
  #            filenames = 'marker_gene_expression_umap.png');
}
#' @export
plot_gene_density<-function(gene_list,seurat){
  gene_list=intersect(gene_list,rownames(seurat@data))
  seurat_data<-melt(as.data.frame(t(as.matrix(seurat@data)))[,gene_list])
  density_plot<-ggplot(seurat_data, aes(x=value,colour=variable))+
    geom_density() + facet_wrap(~variable,scales = 'free')
  col_num<-round(sqrt(length(gene_list)))
  pdf('temp.pdf',width=6*col_num,height=4*ceiling(length(gene_list)/col_num))
  plot(density_plot)#不知为何不加plot()运行不了，但不在function里可以运行
  dev.off()
#pdf_convert('temp.pdf', format = "png", pages = NULL,
 #             filenames = 'marker_gene_expression_density.png');
}
#' @export
PlotCellAttr<-function(mutation, cell_attr, seurat, reduction = "umap", pdfname='temp'){
  #get reduction point
  embeddings.use <- GetDimReduction(object = seurat, reduction.type = reduction,
                                    slot = "cell.embeddings")[,1:2]
  if (length(x = embeddings.use) == 0) {
    stop(paste(reduction.use, "has not been run for this object yet."))
  }
  valid_cells<-intersect(mutation$col.attrs$cell_names[],seurat@cell.names)
  if(length(valid_cells)==0){stop('No shared cell between mutation and expression dataset')}
  embeddings.use<-embeddings.use[valid_cells,]
  #make attr vector
  attr_vector<-mutation$col.attrs[[cell_attr]][]
  names(attr_vector)<-mutation$col.attrs$cell_names[]
  attr_vector<-attr_vector[valid_cells]
  #check attr type
  if(typeof(attr_vector) %in% c("double","integer")){vector_is_continuous=T}else{vector_is_continuous=F}
  #plot
  dataset<-as.data.frame(cbind(embeddings.use,attr_vector))
  colnames(dataset)<-c('Dim1','Dim2','attr')
  dataset$Dim1<-as.numeric(dataset$Dim1)
  dataset$Dim2<-as.numeric(dataset$Dim2)
  if(vector_is_continuous){
    dataset<-dataset[order(dataset$attr,na.last =F),]
    plot<-ggplot(data=dataset)+
      geom_point(aes(Dim1, Dim2,colour=attr),alpha = 0.8, size=1)+
      ggtitle(cell_attr)+
      scale_color_gradient(low="grey",high = 'red',na.value = 'grey90')
  }else{
    plot<-dataset$attr<-as.factor(dataset$attr)
    ggplot(data=dataset)+
      geom_point(aes(Dim1, Dim2,colour=attr),alpha = 0.8, size=1)+
      ggtitle(cell_attr)+
      scale_colour_discrete(na.value = "grey90")
  }
  pdf(paste(pdfname,'pdf',sep='.'))
  print(plot)
  dev.off()
}
#' @export
PlotSNPHeatmap<-function(mutation, cell_range='', site_range='', layer_name='default'){
  cluater_matrix<-t(mutation$layers[[layer_name]][cell_range,site_range])
  cluater_matrix<- replace(cluater_matrix,cluater_matrix==0,-1)
  cluater_matrix<- replace(cluater_matrix,cluater_matrix==3,0)
  pdf('temp.pdf')
  heatmap_result<-pheatmap(cluater_matrix,  color = colorRampPalette(c("blue","white","red"))(3),
                           show_colnames=F, cluster_cols = T,show_rownames=T,cluster_rows = F,
                           labels_row = 'sites',labels_col = 'cells')
  dev.off()
}





