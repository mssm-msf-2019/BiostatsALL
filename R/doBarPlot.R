#' A function that creates a barplot with a table for contrasts
#' @description This function creates a barplot with a table with estimated contrasts.
#' @param d data frame with expression data.
#' @param nn name of the gene
#' @param ID name of the variable with subject ID
#' @param Group  name of factor with groups to be compared
#' @keywords barplot , expression vizualisation
#' @examples 
#' doBarPlot(d,nn,ID,Group)


doBarPlot<-function(d,nn,ID,Group){
  
  doLme<-function(d,nn,ID,Group){
    require(lme4)
    d<-d[,c(ID,Group,nn)]
    d<-na.omit(d)
    d$x<-d[,nn]
    
    
    # function that fits a linear model
    f<-lmer(x~-1+Group+(1|ID),data=d) 
    
    # evaluating sample size
    ssize = as.data.frame(table(d$Group))
    names(ssize)<-c('Group','n')
    
    # contrast from lsmeansls
    # this function creates two tables
    ls.res <-lsmeans(f,pairwise~Group,adjust='none')
    
    contrastTable<-as.data.frame(summary(ls.res$contrasts))
    
    TableFCH <- cbind.data.frame(
      Contrast = contrastTable$contrast,
      FCH = round(contrastTable$estimate,4),
    'P-Value' = round(contrastTable$p.value,5)
    )

TableFCH = TableFCH[order(TableFCH$Contrast),]

Means = as.data.frame(summary(ls.res))

return(list(Means = Means,TableFCH=TableFCH))
  }
  
  
  
  require(ggplot2)

  Rx<-doLme(d=d,nn=nn)
  D<-Rx$Means
  Tb<-Rx$TableFCH
  
  D$Group<-factor(D$lsmeans.Group)
  
  
  
  dodge <- position_dodge(width=0.85)
  
  p<-ggplot(D, aes(x=Group, y=lsmeans.lsmean,fill=Group))+ 
  geom_bar(position=dodge,stat="identity")+
  geom_errorbar(aes(ymin=(lsmeans.lsmean-lsmeans.SE),
                    ymax=(lsmeans.lsmean+lsmeans.SE)),
                position=dodge,width=0.25)


 p<-p+ labs(title = nn,y='Y')

tb<-tableGrob(Tb,show.rownames=FALSE)

fn = paste(nn,'barplot.pdf',sep='')
pdf(file=fn,width=15,height=6)
grid.arrange(arrangeGrob(p,tb,widths=c(3,1),ncol=2))
dev.off()

return(list(p=p,tb=tb,D=D))

}

