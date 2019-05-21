#' Function for cross annotating a mouse expression set with human attributes.
#' @description This function takes a mouse expression set as input, and outputs a merged expression set with mouse and human annotations.
#' @param eset Mouse expression set
#' @param AnnPkg.in Annotation package that corresponds to the mouse expression set.
#' @keywords cross annotation, expression set, eset, Xanno, mouse, human
#' @examples 
#' eset.IL23.probes.Xanno.2 <- annoXeset(eset=eset_IL23,AnnPkg.in = "mouse4302.db")
annoXeset <- function(eset,AnnPkg.in){
  ##Annotating Mouse to Human
  
  ##Depends on the following packages:
  require(biomaRt)
  
  ##Getting human and mouse marts
  human = biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  ##Checking mart attributes - just for testing
  # mouse.attr <- attributes(mouse)$attributes$name
  #   mouse.attr [grep(mouse.attr ,pattern="sym")]
  # human.attr <- attributes(human)$attributes$name
  #   human.attr [grep(human.attr ,pattern="descr")]
  
  ##Annotating combined mouse eset
  eset.probes <- data.frame(rownames(exprs(eset)))
  rownames(eset.probes) <- rownames(exprs(eset))
  eset.probes$PROBEID <- rownames(exprs(eset))
  eset.probes.anno <- BiostatsALL::getAnnotationTable(IDS=eset.probes$PROBEID, AnnPkg=AnnPkg.in, w.annot=c("ENSEMBL","ENTREZID","SYMBOL","GENENAME","CHRLOC","PATH","GO"),w.annot.unique=c("ENSEMBL","ENTREZID","SYMBOL","CHRLOC"))
  eset.probes.anno <- as.data.frame(eset.probes.anno)
  eset.probes.anno <- eset.probes.anno[as.data.frame(eset.probes.anno)$ENSEMBL!="",]
  
  ##Doing cross-annotation
  eset.probes.Xanno <- biomaRt::getLDS(attributes = c("ensembl_gene_id","mgi_symbol"),filters = "ensembl_gene_id",values =as.character(eset.probes.anno$ENSEMBL), mart = mouse , attributesL = c("ensembl_gene_id","hgnc_symbol", "wikigene_description" ), martL = human )
  
  ##Merging Xanno and Eset
  eset.probes.Xanno.2 <- merge(x = eset.probes.anno ,y=eset.probes.Xanno,by.x = "ENSEMBL",by.y="Ensembl.Gene.ID",all.x=T,all.y=F)
  
  return(eset.probes.Xanno.2)
  
}