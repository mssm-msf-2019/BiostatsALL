#' A function that given common gene names in dermatology, create expression with special character (like greek symbols).
#' @description very useful for title and axis in plots. 
#' @param vs: a gene name
#' @examples 
#' makenames_MSF('IL17')

makenames_MSF<-function(vs){
  #if ((substr(vs,1,2)=='IL')&(!substr(vs,3,3)%in%c('.','_','-')))  {vs<-gsub('IL',replacement='IL-',x=vs,fixed=TRUE)}
  
  #	vs<-gsub('KC',replacement='KC: ',x=vs,fixed=TRUE)
  vs<-gsub('-',".",x=vs,fixed=TRUE)
  
  if (vs=='LEAN') {vs='% LEAN'; return(vs)}
  if (vs=='FAT') {vs='% FAT'; return(vs)}
  if (vs%in%c("TNFaplha","TNF.alpha")) {vs="TNFa"}
  
  if ((vs=='IFNg')|(vs=='IFN.g')|(vs=='IFNG')) {vs=expression("IFN" * gamma); return(vs)}
  if (vs=='IFNa') {vs=expression("IFN" * alpha); return(vs)}
  if (vs=='IFN.a') {vs=expression("IFN" * alpha); return(vs)}
  if ((vs=='FCER1')|(vs=='FcER1')) {vs=expression("Fc" * epsilon*"RI"); return(vs)}
  if (vs=='TSLP.R') {vs= 'TSLPR'; return(vs)}
  if ((vs=='IFNa1')|(vs=='IFNA1')) {vs=expression("IFN"* alpha*"1"); return(vs)}
  if ((vs=='TNFa')|(vs=='TNF')) {vs=expression("TNF"* alpha); return(vs)}
  if (vs=='IL17&TNF') {vs=expression("IL17 & TNF"* alpha); return(vs)}
  if (vs=='TNFnotIL17') {vs=expression("unique TNF"* alpha); return(vs)}
  if (vs=='IL17notTNF') {vs='unique IL17'; return(vs)}
  if (vs=='IL22&TNF') {vs=expression("IL22&TNF"* alpha); return(vs)}
  if (vs=='Syn IL17&TNFa') {vs=expression("Syn. IL17 & TNF"* alpha); return(vs)}
  if (vs=='SynergisticIL17andTNFa') {vs=expression("Syn. IL22&TNF"* alpha); return(vs)}
  if (vs=='AdditiveIL17andTNFa') {vs=expression("Add. IL17&TNF"* alpha); return(vs)}
  if (vs=='Add IL22&TNFa') {vs=expression("Add. IL22 & TNF"* alpha); return(vs)}
  if ((vs=='Elafin')|(vs=='ELAFIN')) {vs="PI3"; return(vs)}
  if ((substr(vs,1,2)=='IL')&(!substr(vs,3,3)%in%c('.','_','-')))  {vs<-gsub('IL',replacement='IL-',x=vs,fixed=TRUE)}
  vs<-gsub('IL.',replacement='IL-',x=vs,fixed=TRUE)
  vs<-gsub('IL17',replacement='IL-17',x=vs,fixed=TRUE)
  if (vs=='IL-1B') {vs=expression("IL-1"* beta); return(vs)}
  if  (any(grep('p19|P19',vs))) {vs="IL-23/p19"}
  if (any(grep('p35|P35',vs))) {vs="IL-12/p35"}
  if (any(grep('p40|P40',vs)))  {vs="IL-12/IL-23p40"}
  if (vs=='IL-1B') {vs=expression("IL-1"* beta); return(vs)}
  if (vs=='SCORAD.Impr.s') {vs='SCORAD Improv'}
  if (vs=='IL13plusCLAplus')  {vs=expression('IL-13'^'+'*'CLA'^'+'); return(vs)}
  if (vs=='IL13plusCLAminus') {vs=expression('IL-13'^'+'*'CLA'^'-'); return(vs)}
  if (vs=='IFNgplusCLAplus') {vs=expression('IFN'*gamma^'+'*'CLA'^'+'); return(vs)}
  if (vs=='IFNgplusCLAminus') {vs=expression('IFN'*gamma^'+'*'CLA'^'-'); return(vs)}
  if (vs=='IL22plusCLAplus') {vs=expression('IL-22'^'+'*'CLA'^'+'); return(vs)}
  if (vs=='IL22plusCLAminus') {vs=expression('IL-22'^'+'*'CLA'^'-'); return(vs)}  
  if (vs=='IL9plusCLAplus') {vs=expression('IL-9'^'+'*'CLA'^'+'); return(vs)}
  if (vs=='IL9plusCLAminus') {vs=expression('IL-9'^'+'*'CLA'^'-'); return(vs)}
  if (vs=='CMICOSplusCLAplus')   {vs=expression('CM ICOS'^'+'*'CLA'^'+'); return(vs)}
  if (vs=='CMICOSplusCLAminus')  {vs=expression('CM ICOS'^'+'*'CLA'^'-'); return(vs)}
  if (vs=='EMICOSplusCLAplus')   {vs=expression('EM ICOS'^'+'*'CLA'^'+'); return(vs)}
  if (vs=='EMICOSplusCLAminus')  {vs=expression('EM ICOS'^'+'*'CLA'^'-'); return(vs)}
  if (vs=='CMHLDRplusCLAplus')   {vs=expression('CM HLA-DR'^'+'*'CLA'^'+'); return(vs)}
  if (vs=='CMHLADRplusCLAminus')  {vs=expression('CM HLA-DR'^'+'*'CLA'^'-'); return(vs)}    
  if (vs=='EMHLADRplusCLAplus')   {vs=expression('EM HLA-DR'^'+'*'CLA'^'+'); return(vs)}
  if (vs=='EMHLADRplusCLAminus')  {vs=expression('EM HLA-DR'^'+'*'CLA'^'-'); return(vs)}        
  vs<-gsub('.',replacement='-',x=vs,fixed=TRUE)
  vs<-gsub("_TH","/E.Th.",vs)
  return(vs)
}