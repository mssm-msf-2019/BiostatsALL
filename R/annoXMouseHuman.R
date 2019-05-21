#' Function for annotating and cross annotating top table of mouse DEGs
#' @description This function does cross annotation of a top table of mouse DEGs to human gene Symbols, making use of biomaRt
#' @param top Top table as derived by calling the function toptable on an ebayes fit object
#' @param mouse The mouse mart to use in the biomaRt function
#' @param human The human mart to use in the biomaRt function
#' @param mouse.attr.in Mouse attributes to use from biomaRt.
#' @param human.attr.in Human attributes to use from biomaRt.
#' @param top.tab TRUE/FALSE if "top" is actually a top table
#' @param AnnPkg Mouse eset annotation package. So which package was used to annotate the mouse expression set.
#' @keywords xannotation, cross annotation, mouse, human, annotation
#' @examples 
#' annoXMouseHuman(top=top.A.flakVSflg)


annoXMouseHuman <- function (top=top.A.flakVSflg,mouse="mmusculus_gene_ensembl",human="hsapiens_gene_ensembl",mouse.attr.in=c("ensembl_gene_id","mgi_symbol"),human.attr.in= c("ensembl_gene_id","hgnc_symbol", "wikigene_description" ), top.tab=TRUE,AnnPkg="mogene21sttranscriptcluster.db") {
  
  require(biomaRt)
  #   require(mogene21sttranscriptcluster.db)
  #   top$PROBEID <- rownames(top)
  #   top$SYMBOL <- unlist(as.list(mogene21sttranscriptclusterSYMBOL[rownames(top)]))
  #   top$ENTREZID <- unlist(as.list(mogene21sttranscriptclusterENTREZID[rownames(top)]))
  #   top$GENENAME <- unlist(as.list(mogene21sttranscriptclusterGENENAME[rownames(top)]))
  #   top$ENSEMBL <- unlist(as.list(mogene21sttranscriptclusterENSEMBL[rownames(top)]))[rownames(top)]
  
  top$PROBEID <- rownames(top)
  top.b <- getAnnotationTable(IDS=top$PROBEID, AnnPkg=AnnPkg, w.annot=c("ENSEMBL","ENTREZID","SYMBOL","GENENAME","CHRLOC","PATH","GO"),w.annot.unique=c("ENSEMBL","ENTREZID","SYMBOL","CHRLOC"))
  top.b <- as.data.frame(top.b)
  top.b <- top.b[top.b$ENSEMBL!="",] ##Getting rid of the non-ENSEMBLE annotated probes
  #     print(top.b);print(top)
  top <- merge(x=top,y=top.b,by="PROBEID")
  
  #     head(top,n=15)
  human = useMart("ensembl", dataset = human)
  mouse = useMart("ensembl", dataset = mouse)
  
  ##Doing cross-annotation
  top.Xanno <- getLDS(attributes = mouse.attr.in,filters = "ensembl_gene_id",values = top$ENSEMBL, mart = mouse , attributesL = human.attr.in, martL = human )
  
  top <- merge(x = top ,y=top.Xanno,by.x = "ENSEMBL",by.y="Ensembl.Gene.ID",all.x=T,all.y=F)
  if(top.tab==TRUE){
    top <- top[with(top,order(-logFC)),]
  }
  
  
  return(top)
}