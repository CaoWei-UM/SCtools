#' A mutation data deal Function
#'
#' This function allows you to deal mutation data.
#'
.range_check<-function(vector,mutation_dim){
  vector<-as.integer(vector)
  vector<-vector[!is.na(vector)]
  if(length(vector)==0){return(NULL)}
  if(max(vector)<=mutation_dim&min(vector)>=1){
    return(vector)
  }else{stop('element in range bigger than dimension of mutation')}
}
.margin_check<-function(MARGIN){
  if(MARGIN=='cell'){MARGIN=1}
  if(MARGIN=='site'){MARGIN=2}
  if(!MARGIN %in% c(1,2)){stop('\"MARGIN\" should be 1/cell or 2/site')}
  return(MARGIN)
}

#'
#' @param
#' @keywords deal mutation
#' @export
#' @usage
#' # ReadMutation(mutation_file_path)
#' @examples
#' mutation <- ReadMutation('/home/data/Rworkspace/mutation.txt')
ReadMutation<-function(mutation_file, force_offset=NA, loom_name=NA, force_replace=F){
  #read mutation file
  mutation_list <- read.csv(mutation_file,sep = "\t")
  #split site_annotation and data
  if(is.numeric(force_offset)&force_offset>0){
    mutation_site_annotation<-mutation_list[,1:force_offset]
    mutation_data<-mutation_list[,-c(1:force_offset)]
  }else{
    for(i in 1:dim(mutation_list)[2]){
      one_column<-as.vector(mutation_list[,i])
      if(dim(mutation_list)[1]==length(grep(pattern = "^\\d+:\\d+:\\d+$", x = one_column))){
        break
      }
    }
    mutation_site_annotation<-mutation_list[,1:(i-1)]
    mutation_data<-mutation_list[,-c(1:(i-1))]
  }
  mutation_site_annotation$gene_names<-c(1:dim(mutation_site_annotation)[1])
  # create loom object
  if(is.na(loom_name)){
    loom_name<-strsplit(tail(strsplit(mutation_file,'/')[[1]],n=1),'\\.')[[1]][1]
  }
  loom_link<-loomR::create(paste(loom_name,'loom',sep='.') ,data=mutation_data,
                           gene.attrs=as.list(mutation_site_annotation), overwrite=force_replace)
  return(loom_link)
}

MapMutation<-function(mutation, map_function, layer_name, force_replace=F, replace_default=T){
  #check does layer already exist
  if(layer_name %in% names(mutation$layers)&!force_replace){
    warning(paste(layer_name, 'layer already exist and force_replace is False',sep=''))
    return(mutation)
  }
  #map matrix by map_function
  filtered_matrix<-apply(mutation$matrix[,],1,function(x){
    site_reads_list<-lapply(strsplit(x,':'),as.integer)
    one_site_cells<-lapply(site_reads_list,map_function)
    return(unlist(one_site_cells))
  })
  if(!all(filtered_matrix %in% c(0,1,3))){stop('ERROR: map function return unwanted value')}
  #add layer
  storage.mode(filtered_matrix) <- "integer"
  added_layers<-list()
  added_layers[[layer_name]]<-filtered_matrix
  mutation$add.layer(layers=added_layers, overwrite=force_replace)
  if(replace_default){
    added_layers<-list()
    added_layers[['default']]<-filtered_matrix
    mutation$add.layer(layers=added_layers, overwrite=T)
  }
  return(mutation)
}

#' @export
#' @usage
#' # MapMutationCount(mutation,threshold,min_umi)
#' @examples
#' MapMutationCount(mutation, threshold=2, min_umi=2)
MapMutationCount<-function(mutation, threshold=1, min_umi=1, force_replace=F, replace_default=T){
  #check filter mode
  if(!is.numeric(threshold)){stop("Error: count threshold must be a numerical value")}
  layer_name=paste('min',min_umi,'count',threshold,sep='')
  #make layer
  MapMutation(mutation, layer_name=layer_name, force_replace=force_replace, replace_default=replace_default,
              map_function=function(x){
                if(sum(x)<min_umi){return(3)}
                if(x[2]>=threshold){return(1)}else{return(0)}
                return(3)
              })
  return(mutation)
}

#' @export
#' @usage
#' # MapMutationRatio(mutation,count_threshold,min_umi_counts)
#' @examples
#' MapMutationRatio(mutation, count_threshold=2, min_umi_counts=2)
MapMutationRatio<-function(mutation, ratio_threshold=0.1, min_umi=2, force_replace=F, replace_default=T){
  #check filter mode
  if(!is.numeric(threshold)){stop("Error: ratio threshold must be a numerical value")}
  layer_name=paste('min',min_umi,'ratio',threshold,sep='')
  #make layer
  MapMutation(mutation, layer_name=layer_name, force_replace=force_replace, replace_default=replace_default,
              map_function=function(x){
                if(sum(x)<min_umi){return(3)}
                if(y[2]/sum(x)>=threshold){return(1)}else{return(0)}
                return(3)
              })
  return(mutation)
}

#' @export
#' @usage
#' # UMIcount(mutation)
#' @examples
#' UMIcount(mutation)
UMIcount<-function(mutation, force_replace=F){
  MapMutation(mutation, sum, 'umi_count', force_replace, replace_default=F)
  return(mutation)
}

#' @export
#' @usage
#' # AddSiteAttr(mutation,attr_name,attribute)
#' @examples
#' AddSiteAttr(mutation,'gene',gene_vector)
AddSiteAttr<-function(mutation, attr_name, attribute, force_replace=F){
  if(attr_name %in% names(mutation$row.attrs)&!force_replace){
    warning(paste('\"',attr_name,'\"',' attribute already exist and force_replace is False'))
    return(mutation)
  }
  if(!is.vector(attribute)|length(attribute)!=mutation$shape[1]){
    stop(paste('\"',attr_name,'\"',' must be a vector and length same as mutation\'s site number'))
  }
  add_row_attr<-list()
  add_row_attr[[attr_name]]<-attribute
  mutation$add.row.attribute(add_row_attr,overwrite=force_replace)
  return(mutation)
}

#' @export
#' @usage
#' # AddCellAttr(mutation,attr_name,attribute)
#' @examples
#' AddCellAttr(mutation,'seurat_cluster',cluster_vector)
AddCellAttr<-function(mutation, attr_name, attribute, force_replace=F){
  if(attr_name %in% names(mutation$col.attrs)&!force_replace){
    warning(paste('\"',attr_name,'\"',' attribute already exist and force_replace is False'))
    return(mutation)
  }
  if(!is.vector(attribute)|length(attribute)!=mutation$shape[2]){
    stop(paste('\"',attr_name,'\"',' must be a vector and length same as mutation\'s cell number'))
  }
  add_col_attr<-list()
  add_col_attr[[attr_name]]<-attribute
  mutation$add.col.attribute(add_col_attr,overwrite=force_replace)
  return(mutation)
}

#' @export
FetchAttr<-function(mutation,MARGIN,attr_name,check_function,mode='or',void.value=F,na.value=F){
  #parameter check
  MARGIN<-.margin_check(MARGIN)
  if(MARGIN==1){
    if(attr_name %in% names(mutation$col.attrs)){
      attribute<-mutation$col.attrs[[attr_name]][]
    }else{stop(paste('\"',attr_name,'\"',' attribute not exist in col.attrs'))}
  }
  if(MARGIN==2){
    if(attr_name %in% names(mutation$row.attrs)){
      attribute<-mutation$row.attrs[[attr_name]][]
    }else{stop(paste('\"',attr_name,'\"',' attribute not exist in row.attrs'))}
  }
  if(mode=='or'){op=`|`}else if(mode=='and'){op=`&`}else{stop('mode can be \"or\" or \"and\"')}
  #fetch attr
  result_list<-lapply(attribute,function(x){
    if(length(x)<=0){return(void.value)}
    check_result<-check_function(x)
    check_result[which(is.na(check_result))]<-na.value
    if(typeof(check_result)!='logical'){stop('check_function should return a logical vector')}
    return(Reduce(op,check_result))
  })
  attr_fetched<-which(unlist(result_list))
  return(attr_fetched)
}

#' @export
#' @usage
#' # FetchCell(mutation,attr_name,check_function)
#' @examples
#' TorB_cell<-FetchCell(mutation,'cluster',function(x){x %in% c('T_cell','B_cell')})
FetchCell<-function(mutation,attr_name,check_function,mode='or',void.value=F,na.value=F){
  cell_fetched<-FetchAttr(mutation,1,attr_name,check_function,mode,void.value,na.value)
  return(cell_fetched)
}
#' @export
#' @usage
#' # FetchSite(mutation,attr_name,check_function)
#' @examples
#' mark30<-FetchSite(mutation,'mark',function(x){x>=30})
FetchSite<-function(mutation,attr_name,check_function,mode='or',void.value=F,na.value=F){
  site_fetched<-FetchAttr(mutation,2,attr_name,check_function,mode,void.value,na.value)
  return(site_fetched)
}

#' @export
#' @usage
#' # SitesInGene(mutation,gene_name)
#' @examples
#' site_vector<-SitesInGene(mutation,'EPCAM')
SitesInGene<-function(mutation,gene_name){
  if(!'gene' %in% names(mutation$row.attrs)){stop('You should call \"GeneMutation\" first')}
  gene_site<-FetchSite(mutation,'gene',function(x){x %in% gene_name},mode='or',void.value=F,na.value=F)
  return(gene_site)
}

#' @export
#' @usage
#' # CellInCluster(mutation,cluster_name)
#' @examples
#' cell_vector<-SitesInGene(mutation,'T_cell')
CellInCluster<-function(mutation,cluster_name){
  if(!'cluster' %in% names(mutation$col.attrs)){stop('You should call \"MakeCluster\" first')}
  cell_cluster<-FetchCell(mutation,'cluster',function(x){x %in% cluster_name},mode='or',void.value=F,na.value=F)
  return(cell_cluster)
}


#' @export
#' @usage
#' # GeneMutation(mutation,gene_range_file)
#' @examples
#' GeneMutation(mutation,'/home/data/Rworkspace/gencode_v30_gen_pos.coding.txt')
GeneMutation<-function(mutation, gene_range_file, force_replace=F){
  if(all(c('chromosome','site') %in% names(mutation$row.attrs))){
    chromosome<-mutation$row.attrs$chromosome[]
    site<-mutation$row.attrs$site[]
  }
  gene_range<-read.csv(gene_range_file, sep="\t", header=F)
  unify_chromosome<-function(x){
    nochr_sites<-grep('^chr',x,invert=T)
    x[nochr_sites]<-paste('chr',x[nochr_sites],sep='')
    return(x)
  }
  chromosome<-unify_chromosome(chromosome)
  gene_range$V1<-unify_chromosome(gene_range$V1)
  chromosome_satisfy<-lapply(chromosome,function(x){which(x==gene_range$V1)})
  site_satisfy<-lapply(site,function(x){which(x>=gene_range$V2&x<=gene_range$V3)})
  gene_site<-list()
  for(i in 1:length(chromosome_satisfy)){
    gene_site[[i]]<-gene_range$V4[intersect(chromosome_satisfy[[i]],site_satisfy[[i]])]
  }
  AddSiteAttr(mutation, 'gene', gene_site, force_replace=force_replace)
  return(mutation)
}

MakeCellAttr<-function(mutation, layer_name='default', range_name,
                       attr_function, range=NULL, force_replace=F){
  #parameter check
  range<-.range_check(range,mutation$shape[1])
  if(length(range)==0){stop('range have no element')}
  if(range_name %in% names(mutation$col.attrs)&!force_replace){
    warning(paste('\"',range_name,'\" ',' cell attrs already exist and force_replace is False',sep=''))
    return(mutation)
  }
  #make cell attribute
  attribute<-apply(mutation$layers[[layer_name]][,range,drop=F],1,attr_function)
  AddCellAttr(mutation, range_name, attribute, force_replace=force_replace)
  return(mutation)
}

MakeSiteAttr<-function(mutation, layer_name='default', range_name,
                       attr_function, range=NULL, force_replace=F){
  #parameter check
  range<-.range_check(range,mutation$shape[2])
  if(length(range)==0){stop('range have no element')}
  if(range_name %in% names(mutation$row.attrs)&!force_replace){
    warning(paste('\"',range_name,'\" ',' site attrs already exist and force_replace is False',sep=''))
    return(mutation)
  }
  #make site attribute
  attribute<-apply(mutation$layers[[layer_name]][range,,drop=F],2,attr_function)
  AddSiteAttr(mutation, range_name, attribute, force_replace=force_replace)
  return(mutation)
}

#' @export
#' @usage
#' # CellCoverage(mutation,layer_name,range,range_name)
#' @examples
#' CellCoverage(mutation,range=some_site,range_name='somesite')
CellCoverage<-function(mutation, layer_name='default',
                   range=NULL, range_name=NA, force_replace=F){
  range<-.range_check(range,mutation$shape[1])
  if(is.na(range_name)){
    if(length(range)==0){
      range_name='coverage'
      range=1:mutation$shape[1]
    }else{
      range_name=paste('coverage_range',length(grep('^coverage_range',names(mutation$col.attrs)))+1,sep='_')
    }
  }
  MakeCellAttr(mutation, layer_name=layer_name, range_name=range_name, range=range,
               force_replace=force_replace, attr_function=function(x){length(which(x %in% c(0,1)))})
  return(mutation)
}

#' @export
#' @usage
#' # SiteCoverage(mutation,layer_name,range,range_name)
#' @examples
#' SiteCoverage(mutation,range=some_cell,range_name='somecell')
SiteCoverage<-function(mutation, layer_name='default',
                       range=NULL, range_name=NA, force_replace=F){
  range<-.range_check(range,mutation$shape[2])
  if(is.na(range_name)){
    if(length(range)==0){
      range_name='coverage'
      range=1:mutation$shape[2]
    }else{
      range_name=paste('coverage_range',length(grep('^coverage_range',names(mutation$row.attrs)))+1,sep='_')
    }
  }
  MakeSiteAttr(mutation, layer_name=layer_name, range_name=range_name, range=range,
               force_replace=force_replace, attr_function=function(x){length(which(x %in% c(0,1)))})
  return(mutation)
}

#' @export
#' @usage
#' # CellMutationRatio(mutation,layer_name,range,range_name)
#' @examples
#' CellMutationRatio(mutation, range=some_site,range_name='somesite')
CellMutationRatio<-function(mutation, MARGIN, layer_name='default',
                        min_coverage=5, range=NULL, range_name=NA, force_replace=F){
  range<-.range_check(range,mutation$shape[1])
  if(is.na(range_name)){
    if(length(range)==0){
      range_name='mutation_ratio'
      range=1:mutation$shape[1]
    }else{
      range_name=paste('mutation_ratio_range',
                       length(grep('^mutation_ratio_range',names(mutation$col.attrs)))+1,sep='_')
    }
  }
  MakeCellAttr(mutation, layer_name=layer_name, range_name=range_name, range=range,
               force_replace=force_replace, attr_function=function(x){
                 alert_count<-length(which(x==1))
                 ref_count<-length(which(x==0))
                 if(alert_count+ref_count<min_coverage){return(NA)}
                 return(alert_count/(alert_count+ref_count))
               })
  return(mutation)
}

#' @export
#' @usage
#' # SiteMutationRatio(mutation,layer_name,range,range_name)
#' @examples
#' SiteMutationRatio(mutation, range=some_cell,range_name='somecell')
SiteMutationRatio<-function(mutation, MARGIN, layer_name='default',
                            min_coverage=5, range=NULL, range_name=NA, force_replace=F){
  range<-.range_check(range,mutation$shape[2])
  if(is.na(range_name)){
    if(length(range)==0){
      range_name='mutation_ratio'
      range=1:mutation$shape[2]
    }else{
      range_name=paste('mutation_ratio_range',
                       length(grep('^mutation_ratio_range',names(mutation$row.attrs)))+1,sep='_')
    }
  }
  MakeCellAttr(mutation, layer_name=layer_name, range_name=range_name, range=range,
               force_replace=force_replace, attr_function=function(x){
                 alert_count<-length(which(x==1))
                 ref_count<-length(which(x==0))
                 if(alert_count+ref_count<min_coverage){return(NA)}
                 return(alert_count/(alert_count+ref_count))
               })
  return(mutation)
}

#' @export
#' @usage
#' # DiffSite(mutation,group1,group2,layer_name,range)
#' @examples
#' DiffSite(mutation, group1=some_cell, group2=other_cell)
DiffSite<-function(mutation,group1,group2=NULL,layer_name='default',range=NULL){
  group1<-.range_check(group1,mutation$shape[2])
  group2<-.range_check(group2,mutation$shape[2])
  range<-.range_check(range,mutation$shape[1])
  if(length(group1)==0){stop('\"group1\" have no elements')}
  if(length(group2)==0){group2<-setdiff(c(1:mutation$shape[2]),group1)}
  if(length(range)==0){range=1:mutation$shape[1]}
  stat1<-t(apply(mutation$layers[[layer_name]][group1,range,drop=F],2,function(x){
    c(length(which(x==0)),length(which(x==1)))}))
  stat2<-t(apply(mutation$layers[[layer_name]][group2,range,drop=F],2,function(x){
    c(length(which(x==0)),length(which(x==1)))}))
  stat<-cbind(stat1,stat2)
  colnames(stat)<-c('group1_ref','group1_alert','group2_ref','group2_alert')
  filtered_site<-apply(stat,1,function(x){x[2]+x[4]>=3&x[1]+x[2]>5&x[3]+x[4]>5&x[4]<=1})
  return(range[filtered_site])
}

#' @export
#' @usage
#' # AllDiffSite(mutation,cell_attr_name,layer_name,range)
#' @examples
#' AllDiffSite(mutation, 'cluster')
AllDiffSite<-function(mutation, cell_attr_name, layer_name='default',range=NULL){
  cell_cluster<-mutation$col.attrs[[cell_attr_name]][]
  all_cluster_diff<-list()
  for(i in 1:length(names(table(cell_cluster)))){
    all_cluster_diff[[i]]<-DiffSite(mutation,group1=which(cell_cluster==names(table(cell_cluster))[i]),
                                  layer_name=layer_name,range=range)
  }
  return(all_cluster_diff)
}

