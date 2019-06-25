#' A expression data deal Function
#'
#' This function allows you to deal expression data.
#' @param
#' @keywords deal expression
#' @export
#' @examples
#' top_expression_cell_list(gene_list,CTC)

top_expression_cell_list<-function(gene_list,seurat,top=3){
  cell_list<-list()
  gene_expression_data<-as.data.frame(t(as.matrix(seurat@data)))[,gene_list]
  for(i in 1:length(gene_list)){
    temp_gene<-as.vector(gene_expression_data[,gene_list[i]])
    names(temp_gene)<-rownames(gene_expression_data)
    top_cells<-sort(temp_gene, decreasing=T)[c(1:top)]
    cell_list[[gene_list[i]]]<-names(top_cells[which(top_cells>0)])
  }
  return(cell_list)
}

AddUMAP<-function(mutation,seurat,force_replace=F){
  cell_name<-mutation$col.attrs$cell_names[]
  UMAP_dim<-as.data.frame(seurat@dr$umap@cell.embeddings[cell_name,])
  AddCellAttr(mutation,'UMAP1',UMAP_dim$UMAP1,force_replace=force_replace)
  AddCellAttr(mutation,'UMAP2',UMAP_dim$UMAP2,force_replace=force_replace)
  return(mutation)
}
