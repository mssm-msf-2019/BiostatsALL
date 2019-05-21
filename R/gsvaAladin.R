#' Internal function of GSVA Aladin and Normalized Aladin Scores 
#' @description It calculates the Aladin and Normalized Aladin scores 
#' @param ESET: expression matrix
#' @param PTWpsets: list of pathways 
#' @param NullHypotesis_index: condition for the selection of the samples 
#' @examples
#' gsvaAladin(ESET=exprs(reset),PTWpsets=JKList,NullHypotesis_index=(pData(reset)$Tissue=="Normal"))


gsvaAladin<-function(ESET, PTWpsets, NullHypotesis_index) {
  nulls.m<-apply(exprs(ESET)[,NullHypotesis_index],1,function(x){c(mean(x,na.rm=T))})
  nulls.sd<-apply(exprs(ESET)[,NullHypotesis_index],1,function(x){sd(x,na.rm=T)})
  
  probes.in<-names(nulls.m)[which(nulls.sd>0)]
  dist2null<-(exprs(ESET)[probes.in,]-nulls.m[probes.in])/nulls.sd[probes.in]
  aladinDF<-lapply(PTWpsets, function(ps.ids,b){
    ps.ids<-intersect(unlist(ps.ids),rownames(b)); 
    if (length(ps.ids)>1) out<-c(colSums(b[ps.ids,]),length(ps.ids)) else out<-c(b[ps.ids,],1)
    return(out)
  }, b=dist2null)
  aladin_mat<-do.call('rbind',lapply(aladinDF,function(x){x[-length(x)]}))
  aladin_mat_length<-do.call('rbind',lapply(aladinDF,function(x){x[-length(x)]/sqrt(x[length(x)])}))
  aladin_scores<-new('ExpressionSet',exprs=aladin_mat, phenoData=phenoData(ESET))
  aladin_scores_norm<-new('ExpressionSet',exprs=aladin_mat_length, phenoData=phenoData(ESET))
  return(list(scores=aladin_scores, normalized.scores=aladin_scores_norm))
}
