#' A function that creates a table with the results of limma models
#' @description a function that takes ebfit as input and create a table with DEG genes

#' @param ebfit has estimates of contrasts
#' @param AnnPkg package used to annotate data (example: hgu133plus2.db)
#' @param Decision is the matrix with -1 0 1 indicating gene regulation
#' @param coefs is a matrix with the coefficients of the models
#' @param adjmeth is the criterion for adjusting pvalues
#' @param mcut is the cutoff for log2 fold-changes
#' @param fname is the output filename
#' @param pcut is the threshold for pvalues
#' @param html.genetable logical parameter to create html table with results
#' @param exprset is the matrix with gene expressions
#' @param annotTable is a table with annotation about the gens
#' @param which.annot is a selection of which annotations should be used
#' @param which.annot.unique is a selection of which annotations should be unique
#' @param pPath is the pvalue cutoff for pathway analysis
#' @param pGO is the pvalue cutoff for annotation in Gene Ontology Database
#' @param nPath ?????
#' @param pPFAM is the cutoff for annotation in Protein Family DataBase
#' @param html.genetable TRUE if you want to generate html tables
#' @examples
#' CreateOutputsMA(ebfit, AnnPkg, Decision=NULL, coefs=c(1:ncol(ebfit$coef)), adjmeth='fdr', mcut=2,
#' fname="", pcut=0.05, html.genetable =TRUE, exprset=NULL, annotTable=NULL,
#' which.annot= c("ENTREZID","SYMBOL","GENENAME","OMIM","CHRLOC"),
#' which.annot.unique=c("ENTREZID","SYMBOL","OMIM","CHRLOC"),
#' pPath=0.1, pGO=pPath, nPath=2, pPFAM=pPath)

CreateOutputsMA<-function(ebfit, AnnPkg, Decision=NULL, coefs=c(1:ncol(ebfit$coef)), adjmeth='fdr', mcut=2,
                          fname="", pcut=0.05, html.genetable=TRUE, exprset=NULL, annotTable=NULL,
                          which.annot= c("ENTREZID","SYMBOL","GENENAME","OMIM","CHRLOC"),
                          which.annot.unique=c("ENTREZID","SYMBOL","OMIM","CHRLOC"),
                          pPath=0.1, pGO=pPath, nPath=10, pPFAM=pPath){

  library(ReportingTools)
  library(affycoretools)
  library(lattice)
  library(hwriter)
  library(Category)
  addEGIDLink <- function(object, ...){
    object$ENTREZID <- hwrite(as.character(object$ENTREZID),
                              link = paste0("http://www.ncbi.nlm.nih.gov/gene/", as.character(object$ENTREZID)), table = FALSE)
    return(object)
  }
  addOMIMLink <- function(object, ...){
    object$OMIM <- hwrite(as.character(object$OMIM),
                          link = paste0("http://www.ncbi.nlm.nih.gov/omim/", as.character(object$OMIM)), table = FALSE)
    return(object)
  }

  ebfit<-ebfit[,coefs]
  coefname<-colnames(ebfit$coef)
  IDS <-rownames(ebfit$coef)
  ps<-ifelse((adjmeth=='none'),'P','FDR')
  fname=paste(fname,ps,substr(as.character(pcut),3,4),'FCH',as.character(mcut),sep='')
  reportDirectory = paste("./reports_",fname,sep='')

  ### DO DECISION IF NULL!!!!!!!!!!!!!!!!!
  #get annotation
  if (is.null(annotTable)){ annotTable <-getAnnotationTable(IDS, AnnPkg, which.annot, which.annot.unique) }

  ebfit$genes<-annotTable
  TAB<-FromEbfit2Table(ebfit,adj=adjmeth,mcut=mcut,pcut=pcut, annot=FALSE)
  TAB<-cbind.data.frame(TAB, annotTable[rownames(TAB), -grep('Other',colnames(annotTable))])
  if (!is.null(exprset)) TAB<-cbind.data.frame(TAB, exprset[rownames(TAB),])
  TAB<-cbind.data.frame(TAB, annotTable[rownames(TAB), grep('Other',colnames(annotTable))])

 # print(colnames(TAB))

  ### save ALL genes
  ALL.file<-CSVFile(shortName = paste(fname,'AllGenes',sep='_'),reportDirectory = reportDirectory)
  publish(TAB, ALL.file)

  ### SOME DIFF
  StatusTab<-as.matrix(TAB[,grep('Status',colnames(TAB))])
  DEGs<-rownames(TAB)[which(rowSums(abs(StatusTab))>0)]
  DEGs.file<-CSVFile(shortName = paste(fname,'SomeDEGs',sep='_'),reportDirectory = reportDirectory)
  publish(TAB[DEGs,], DEGs.file)

  N.DEG.file<-CSVFile(shortName = paste(fname,'NumberOfDEG.csv',sep='_'),reportDirectory = reportDirectory)
  publish(as.data.frame(t(summary(Decision)[-2,])), N.DEG.file)


  ### save one table per coefficient
  lattice.options<-(default.theme=reporting.theme())
  UnivEntrez<- keys(eval(as.name(AnnPkg)), keytype="ENTREZID")
  for (i in c(1:length(coefs))){
    cni_original<-coefname[i]
    cni<-gsub('_','.',cni_original,fixed=TRUE)  ### This is bc FromEbfit2Table make this changes to names resulting in TAB
    print(paste('Running analysis for contrast:',cni))
    suffix<-c(paste('StatusFCH',mcut,ifelse(adjmeth =='none','P','FDR'),pcut,"_",sep=''),'lgFCH_','FCH_','pvals_','fdrs_')
    colnamesTab<-mgsub(suffix,rep("",length(suffix)),colnames(TAB))

    Table_cni <-TAB[,which(colnamesTab==cni)]
    colnames(Table_cni)<-gsub('_','',gsub(cni,"",colnames(Table_cni)))

    if (any(Table_cni[,grep('Status',colnames(Table_cni))]!=0)){    ## if there is some DEGs

      DEGi<-rownames(Table_cni)[which(Table_cni[,grep('Status',colnames(Table_cni))]!=0) ]      ## get DEG only
      Table_cni<-cbind(TAB[DEGi,c('PROBEID', which.annot[1:4])],
                       Table_cni[DEGi,],
                       TAB[DEGi,setdiff(colnames(annotTable),c('PROBEID', which.annot[1:4]))])
      Table_cni<-Table_cni[order((-Table_cni[,grep('Status',colnames(Table_cni))]),-abs(Table_cni$lgFCH)),]

      Tit=paste(cni_original,' DEG (FCH>',as.character(mcut),', ', ifelse(adjmeth =='none','P','FDR'),'<',as.character(pcut),')',sep='')

      if (length(DEGi)>0){
        #save in csv
        csv.File<-CSVFile(shortName = paste(fname,cni_original,sep='_'), reportDirectory = reportDirectory)
        publish(Table_cni, csv.File)
        print('DEGs published in csv file ..')

        if (html.genetable) {
          ####  save html tables
          colnames(Table_cni)<-gsub(paste('FCH',mcut,ifelse(adjmeth =='none','P','FDR'),pcut,sep=''),"",colnames(Table_cni))
          htab <- HTMLReport(paste(fname,cni_original,sep='_'), reportDirectory = reportDirectory)
          publish(Table_cni, htab, tableTitle =Tit, .modifyDF = list(addEGIDLink, addOMIMLink))
          finish(htab)
          print('DEGs published in html file...')
        }

        if (any(!is.null(c(pPath,pGO,pPFAM)))){
             DEGi_entrez<-uniquenona(AnnotationDbi::select(eval(as.name(AnnPkg)),keys=DEGi,columns='ENTREZID')[,'ENTREZID'])

             if ((!is.null(pGO))&(length(DEGi_entrez)>nPath)){
               ####  go report
               goParams<-new("GOHyperGParams", geneIds = DEGi_entrez, universeGeneIds=UnivEntrez,
                             annotation=eval(as.name(AnnPkg))$packageName, ontology="BP",
                             pvalueCutoff = pGO, conditional=TRUE, testDirection = "over")
               goResults<-hyperGTest(goParams)
               gotab <-HTMLReport(shortName = paste("Gene Ontology BP",cni_original,sep=''),
                                  title='GO analysis', reportDirectory=reportDirectory)
               publish(goResults, gotab, selectedIDs=DEGi_entrez, annotation.db="org.Hs.eg",
                       categorySize=nPath,makePlot=FALSE)
               finish(gotab)
               goParams2<-new("GOHyperGParams", geneIds = DEGi_entrez, universeGeneIds=UnivEntrez,
                              annotation=eval(as.name(AnnPkg))$packageName, ontology="MF",
                              pvalueCutoff = pGO, conditional=TRUE, testDirection = "over")
               goResults2<-hyperGTest(goParams2)
               gotab2 <-HTMLReport(shortName = paste("Gene Ontology MF",cni_original,sep=''), title='GO analysis MF',
                                   reportDirectory=reportDirectory)
               publish(goResults2, gotab2, selectedIDs=DEGi_entrez, annotation.db="org.Hs.eg",
                       categorySize=nPath, makePlot=FALSE)
               finish(gotab2)
               print('GO published in html file...')

             }  	  #### end go report

             if ((!is.null(pPFAM))&(length(DEGi_entrez)>nPath)){
               pfamParams<-new("PFAMHyperGParams", geneIds = DEGi_entrez, universeGeneIds=UnivEntrez,
                               annotation=eval(as.name(AnnPkg))$packageName, pvalueCutoff=pPFAM,
                               testDirection = "over")
               pfamResults<-hyperGTest(pfamParams)
               pfamtab <-HTMLReport(shortName = paste("PFAM Analysis",cni_original,sep=''),
                                    title='PFAM Analysis', reportDirectory=reportDirectory)
               publish(pfamResults, pfamtab, selectedIDs=DEGi_entrez, annotation.db="org.Hs.eg",
                       categorySize=nPath, makePlot=FALSE, pvalueCutoff=pPFAM)
               finish(pfamtab)
               print('PFAM published in html file...')
             }          #### end pfam analysis

             if ((!is.null(pPath))&(length(DEGi_entrez)>nPath)){
               keggParams<-new("KEGGHyperGParams", geneIds = DEGi_entrez, universeGeneIds=UnivEntrez,
                               annotation=eval(as.name(AnnPkg))$packageName, pvalueCutoff = pPath, testDirection = "over")
               keggResults<-hyperGTest(keggParams)
               htmlReport(keggResults, file=file.path(reportDirectory,paste("KEGG Analysis",cni_original,'.html',sep='')), summary.args=list("htmlLinks"=TRUE,  pvalue=pPath,categorySize=20))

               # I communicated with Gabe Becker from genentech and tried his solution but could not make it work
               #                    keggtab <-HTMLReport(shortName = paste("KEGG Analysis",cni_original,sep=''),
               #                                        title='KEGG Analysis', reportDirectory="./reports")
               #
               #                     publish(keggResults, .toDF = getMethod(toReportDF, "PFAMHyperGResult"))
               #                       publish(keggResults, htmlReport=keggtab, .toDF = getMethod(toReportDF, "PFAMHyperGResult"), selectedIDs=DEGi_entrez, annotation.db="org.Hs.eg",categorySize=20, makePlot=FALSE)
               #                       finish(keggtab)
               print('KEGG published in html file...')
             }          #### end pfam a
        }  #end if pathway analisys
      } ## end if no DEGs
      #### save pathway report

      ### create index
      #indx <- HTMLReport("index", "Main page")
      #publish(Link(htab), indx)
      #finish(indx)
    }
  }   #for


}
