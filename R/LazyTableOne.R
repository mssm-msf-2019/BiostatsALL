

##LazyTableOne
##before entering this function, we shall make sure all variables are in the right class we want
#' A function that output a summary table with descriptive statistics and group comparison results
#' @description: a function that can produce Table1 and conduct group comparisons with given group variable and covariates of interests
#' @param vars covariates that require statistical analysis (Make sure categorical variables u want is not in integer class)
#' @param strata the grouping variable, could be binary or non-binary
#' @param data per ID dataset, or wide format paired data, as data frame
#' @param paired indicator if data is paired or not
#  LazyTableOne (vars=c('Age','Gender',Ethnicity','IBD_RiskScore','DiseaseSeverity'),strata='DiseaseOnset',data=PD_phase23,paired=F)
#  install package 'tableone' first

LazyTableOne<-function(vars, strata, data, addname,digits=4,paired=FALSE,includeNA = FALSE,test = TRUE,smd = F){

  require(tableone)
  dir.T1<-paste0('./Table1_',addname,'/'); 
  dir.create(dir.T1)
  
  if(length(levels(data[,strata]))==2){
    if (paired){
      testApprox=mcnemar.test
      testExact = mcnemar.test
      testNormal = t.test
      testNonNormal = wilcox.test
      argsApprox=list(correct = TRUE)
      argsExact = list( correct = TRUE)
      argsNormal = list(paired = TRUE, var.equal = FALSE)
      argsNonNormal = list(paired = TRUE, exact = NULL, correct = TRUE)
    } else{
      testApprox=chisq.test
      testExact = fisher.test
      testNormal = t.test
      testNonNormal = wilcox.test
      argsApprox=list(correct = TRUE)
      argsExact = list(workspace = 2*10^5)
      argsNormal = list(paired = FALSE, var.equal = FALSE)
      argsNonNormal = list(paired = FALSE, exact = NULL, correct = TRUE)
    }
  } else if(length(levels(data[,strata]))>2){
    ## under >2 sample comparison 
    testApprox=chisq.test
    testExact = NULL
    testNormal = oneway.test## is it one way anova?
    testNonNormal = kruskal.test
    argsApprox=list(correct = TRUE)
    argsExact = list(NULL)
    argsNormal = list(var.equal = TRUE)##?
    argsNonNormal = list(NULL)
  }
  
  Table1<-CreateTableOne(vars, strata, data, includeNA = FALSE, 
                         testApprox = testApprox, argsApprox =argsApprox, 
                         testExact = testExact, argsExact = argsExact, 
                         testNormal = testNormal, argsNormal =argsNormal, 
                         testNonNormal = testNonNormal, argsNonNormal =argsNonNormal, 
                         smd =F,test = T)
  
  write.csv(print(Table1),file = paste0(dir.T1,'Descriptive_Statistics_',addname,'_',Sys.Date(),'.csv'))
  
  ## to get the both test results that used for each test
  sink(file = paste0(dir.T1,'Descriptive_Statistics_Withbothtests',addname,'_',Sys.Date(),'.txt'))
  summary(Table1$ContTable,digits=digits)
  summary(Table1$CatTable,digits=digits)
  sink()
  ## we can do create 
  return(Table1)
}