## I. geoseq for main signature analysis ----

#' @description Conducting Gene Oncology test in RNA seq data with potential gene length bias, using Limma model result
#' @param Contr: Contrast to be selected for DEGs
#' @param ResTab: Limma results table, from 'FromEbfit2Table'
#' @param first.ann.cols: annotation attributes 
#' @param annotab: gene annotation table
#' @param CollectionList: a list of gene collection(s) composed of gene symbols
#' @param slope: threshold to determine if gene length adjustment is need, by checking pwf plot
#' @param dir: working directory
#' @param direction_dysregulation: direction of dysregulation:'both','up' or 'down'
#' 

getMeGoSeqfrom_BigTab<-function(Contr,ResTab,first.ann.cols,
                                annotab,CollectionList,slope,dir,direction_dysregulation='both'){
  
  library(goseq)
  library(plyr)
  print(Contr)
  library(BiostatsALL)  ### use my fancy largest 
  
  ## I. get the unique annotation by the rank of FCH from the comparison
  strinContr<-paste0('_',Contr)
  ResTab<-ResTab[,c(grep(strinContr,colnames(ResTab),value=T),first.ann.cols)]
  ResTab<-subset(ResTab,!is.na(external_gene_name))##exclude NAs in genesymbols
  
  status_select=c(1,-1)
  if (direction_dysregulation=='up') status_select=c(1) 
  if (direction_dysregulation=='down') status_select=c(-1)
  
  ResTab<-mutate(ResTab,
                 status=ifelse(ResTab[,grep('Status',colnames(ResTab))]%in%status_select,1,0),
                 LengthGene=end_position-start_position
  )
  FCHlist<-subset(ResTab,select=grep('^FCH',colnames(ResTab),value=T))[,1];names(FCHlist)<-rownames(ResTab)
  FCHlist<-sort(FCHlist*sign(FCHlist),decreasing=T)
  
  geneSymbols<-as.character(subset(annotab,ensembl_gene_id%in%rownames(ResTab))$external_gene_name)
  names(geneSymbols)<-as.character(subset(annotab,ensembl_gene_id%in%rownames(ResTab))$ensembl_gene_id)
  geneSymbols<-geneSymbols[!is.na(geneSymbols)]
  
  symbs_byComprRank<-myfancy_findLargest(probe.ids=rownames(ResTab),
                                         testStat=FCHlist,
                                         gene.names=geneSymbols)###
  
  ResTab.rank.sym<-ResTab[symbs_byComprRank,]## now it's one to one match
  
  genes<-ResTab.rank.sym$status;names(genes)<-ResTab.rank.sym$external_gene_name ##12391 reduced to 12336
  
  ## II. Create Null Distribution
  pwf=nullp(genes,"hg19","symbolGene",bias.data = ResTab.rank.sym$LengthGene)
  pdf(file=paste(dir,'PWF Plot ',Contr,'.pdf'),width=4.5,height=4)
  p<-plotPWF(pwf)+title(paste0('PWF Plot ',Contr))
  dev.off()
  row.names(pwf)<-as.character(ResTab.rank.sym$external_gene_name)
  
  ## III. over presentation analysis
  aux_gene2cat<-lapply(CollectionList, function(L){
    gene2cat_matrix<-L
    gene2cat_string<-sapply(L, function(x){paste0(intersect(x,DEGs), collapse=', ')})
    return(list(gene2cat_matrix=gene2cat_matrix,
                gene2cat_string=data.frame(category=names(gene2cat_string),genes_hit=gene2cat_string)))
  })
  gene2cat_db<-lapply(aux_gene2cat,function(L){L$gene2cat_matrix})
  gene2cat_string<-do.call('rbind',lapply(aux_gene2cat,function(L){L$gene2cat_string}))##ok
  
  
  all.goseq<-lapply(gene2cat_db,function(G){
    Lm<-lm(pwf~bias.data,data=pwf)
    if(Lm$coefficients[2]>slope){
      #all.goseq = goseq(pwf,'hg19','symbolGene',gene2cat = C.All)
      all.goseq = goseq(pwf,genome='hg19',id='geneSymbol',gene2cat = gene2cat_db, method="Wallenius")
    } else {
      print(paste0(Contr,' analyzed without adjusting for selection bias'))
      # all.goseq=goseq(pwf,"hg19","geneSymbol",method="Hypergeometric") ### THIS IS WRONG!
      all.goseq=goseq(pwf,genome="hg19",id="geneSymbol", gene2cat = gene2cat_db, method="Hypergeometric")
    }})
  all.goseq<-do.call('rbind',all.goseq)
  
  ## add hit genes----
  all.goseq.results<-merge(all.goseq, gene2cat_string,by='category', all.x=TRUE,all.y=FALSE)
  all.goseq.results_DE<-subset(all.goseq.results,numDEInCat>1)
  write.csv(all.goseq.results,file=paste0(dir,'Goseq_analysis_',Contr,'_',Sys.Date(),'.csv'))
  return(invisible(list(pwf=pwf,all.goseq=all.goseq.results)))
}
