## I. geoseq for medication signature analysis ----

#' @description Conducting Gene Oncology test in RNA seq data with potential gene length bias, using Limma model result
#' @param Stats: statistics(name) to determine DEG status, can be pvalue or fdr
#' @param ResTab: Model results table
#' @param pcut: cutoffs of to determine DEGs, according to 'Stats'
#' @param rank: statistics to rank the gene. e.g. rank='logFC'
#' @param first.ann.cols: annotation attributes 
#' @param decreasing: The default is TRUE. 
#' @param annotab: gene annotation table
#' @param CollectionList: a list of gene collection(s) composed of gene symbols
#' @param slope: threshold to determine if gene length adjustment is need, by checking pwf plot
#' @param addname: additional name to specify the model, e.g. type of medication
#' @param dir: working directory
#' @param use_genes_without_cat: The default is FALSE. 
#' 

getMeGoSeqfrom_BigTab_Meds<-function(Stats,ResTab,pcut=0.1,rank,first.ann.cols,
                                     decreasing=T,annotab,CollectionList,slope,addname,dir,use_genes_without_cat=FALSE){
  library(goseq)
  library(BiostatsALL)  ### use my fancy largest 
  ## I. get the unique annotation by the rank of FCH from the comparison
  
  ResTab<-ResTab[,c(Stats,rank,first.ann.cols,Stats)]
  ResTab<-subset(ResTab,!is.na(external_gene_name))##exclude NAs in genesymbols
  
  ResTab<-mutate(ResTab,
                 status=ifelse(ResTab[,Stats]>pcut,0,1),
                 LengthGene=end_position-start_position
  )
  FCHlist<-subset(ResTab,select=rank)[,1];names(FCHlist)<-rownames(ResTab)
  FCHlist<-sort(2^(FCHlist)*sign(FCHlist),decreasing=T)
  
  geneSymbols<-as.character(subset(annotab,ensembl_gene_id%in%rownames(ResTab))$external_gene_name)
  names(geneSymbols)<-as.character(subset(annotab,ensembl_gene_id%in%rownames(ResTab))$ensembl_gene_id)
  geneSymbols<-geneSymbols[!is.na(geneSymbols)]
  
  symbs_byComprRank<-myfancy_findLargest(probe.ids=rownames(ResTab),
                                         testStat=FCHlist,
                                         gene.names=geneSymbols)### 
  
  ResTab.rank.sym<-ResTab[symbs_byComprRank,]## now it's one to one match
  
  genes<-ResTab.rank.sym$status;names(genes)<-ResTab.rank.sym$external_gene_name
  
  ## II. Create Null Distribution
  pwf=nullp(genes,"hg19","symbolGene",bias.data = ResTab.rank.sym$LengthGene)
  
  ## sarah, discuss this with Mayte,if
  ## initial point very close to some inequality constraints
  ## then what should we do
  
  pdf(file=paste(dir,'PWF Plot ',addname,'.pdf'),width=4.5,height=4)
  p<-plotPWF(pwf)+title(paste0('PWF Plot '))
  dev.off()
  row.names(pwf)<-as.character(ResTab.rank.sym$external_gene_name)
  
  ## III. over presentation analysis
  
  ## test
  # L<-CollectionList[[1]]
  
  ########
  
  DEGs=names(genes)[genes!=0]
  
  aux_gene2cat<-lapply(CollectionList, function(L){
    # GeneC<-do.call('rbind',L)
    # gene2cat_matrix = cbind.data.frame(GeneC,
    #                 Category=rep(names(L),lapply(L,function(x)length(x))))
    gene2cat_matrix<-L
    gene2cat_string<-sapply(L, function(x){paste0(intersect(x,DEGs), collapse=', ')})
    return(list(gene2cat_matrix=gene2cat_matrix,
                gene2cat_string=data.frame(category=names(gene2cat_string),genes_hit=gene2cat_string)))
  })
  
  #gene2cat_db<-do.call('rbind',lapply(aux_gene2cat,function(L){L$gene2cat_matrix}))##wrong
  #gene2cat_db<-do.call(c,lapply(aux_gene2cat,function(L){L$gene2cat_matrix}))
  gene2cat_db<-lapply(aux_gene2cat,function(L){L$gene2cat_matrix})
  gene2cat_string<-do.call('rbind',lapply(aux_gene2cat,function(L){L$gene2cat_string}))##ok
  
  ## head(supportedGeneIDs())[1;100,c('tablename','GeneID')]
  ## supportedOrganisms()[supportedOrganisms()$Genome=="hg19",]
  
  Lm<-lm(pwf~bias.data,data=pwf)
  
  all.goseq<-lapply(gene2cat_db,function(G){
    if(Lm$coefficients[2]>slope){
      #all.goseq = goseq(pwf,'hg19','symbolGene',gene2cat = C.All)
      goseq = goseq(pwf,genome='hg19',id='geneSymbol',gene2cat = G, method="Wallenius",use_genes_without_cat=FALSE)
    } else {
      print(paste0(' analyzed without adjusting for selection bias'))
      # all.goseq=goseq(pwf,"hg19","geneSymbol",method="Hypergeometric") ### THIS IS WRONG!
      goseq=goseq(pwf,genome="hg19",id="geneSymbol", gene2cat = G, method="Hypergeometric",use_genes_without_cat=FALSE)
    }})
  all.goseq<-do.call('rbind',all.goseq)
  
  ## add hit genes----
  all.goseq.results<-merge(all.goseq, gene2cat_string,by='category', all.x=TRUE,all.y=FALSE)
  all.goseq.results_DE<-subset(all.goseq.results,numDEInCat>1)
  write.csv(all.goseq.results_DE,file=paste0(dir,'Goseq_analysis_',addname,'_',Sys.Date(),'.csv'))
  return(invisible(list(pwf=pwf,all.goseq=all.goseq.results)))
}

