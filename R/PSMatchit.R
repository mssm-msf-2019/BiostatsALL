
##PSMatchit

#' A function to perform propensity matching on Binary and Non-Binary group given grouping variable and set of covariates
#' @description:
#' @param dat: dat per id dataset
#' @param id : unique identification
#' @param GrpVar : Binary group variable e.g. Treat vs. Control
#' @param Covars : a set of variables to be matched between the 2 groups.
#' @param Method : nearest is most commonly used, see ?matchit/trimatch--method for other options
#' @param Ratio : ratio of matching, 2 if is 2:1 matching
#' @param addname: name to be added in results folder
#' @param digit:dicimal to be retained in summary table
#' @details This function is based on Matchit from package 'MatchIt';The reference group of your group variable is the one that you want to match for. check for the class of each covariate to ensure it's the right class(categorical or numeric)

# example from MIND study: MINDPso_Psm<-PSMatchit(PD_IBD_pat,'GRID','Mtch_Medications_AntTNFa',Covars,Method='nearest',Ratio=2,digit=3)

PSMatchit<-function(dat,id,GrpVar,Covars,Method,Ratio,addname=GrpVar,digit=3){
  library(MatchIt)
  library(optmatch)
  library(tableone)
  library(TriMatch)
  
  dir.PSM<-paste0('./propensity_score_',addname,'/'); 
  dir.create(dir.PSM)
  
  dat<-dat[,c(id,GrpVar,Covars)]
  dat[,GrpVar]<-as.factor(as.character(dat[,GrpVar]))
  dat[,id]<-as.character(dat[,id])
  
  tab1Unadj <- CreateTableOne(vars = Covars, strata = GrpVar, data = dat)
  write.csv(as.data.frame(print(tab1Unadj)),paste0(dir.PSM,'PSM_Table1_BeforeM',addname,'_',Sys.Date(),'.csv'))
    
## 1. check data 
## ensure no duplicates at id level
if (nrow(dat)==nrow(subset(dat,!duplicated(dat[,id])))){
  print('Check 1:Data is per patient(obs)')
} else {
  stop('Data has duplicates')
}

  ## ensure no missng among covariates
  print(paste0('Original # of cases ',nrow(dat)))
  dat<- dat[complete.cases(dat), ]
  rownames(dat)<-dat[,id]
  print(paste0('Completed # of cases ',nrow(dat)))
  
## ensure group variable a binary variable
if (length(levels(dat[,GrpVar]))==2){
  print('Check 2:Group variable is Binary')

## 2.prepare formula 
fmla <- as.formula(paste(paste0(GrpVar,"~"), paste(Covars, collapse= "+")))

## 3. PSM--for binary outcomes
match.out.all<- matchit(fmla, data = dat, method = Method,  model="logit"
                        ,ratio = Ratio,caliper=0.1,reestatimate=T,replace=F)

## 4. outputs
pdf(file=paste0(dir.PSM,'Plots_PSM_',addname,Sys.Date(),'.pdf'),height=6,width=6,onefile=T)
p1<-plot(match.out.all,type='hist')
p2<-plot(match.out.all,type='jitter') ## press esc to exit..
dev.off()

MatchedID<-match.out.all$match.matrix
colnames(MatchedID)<-paste0('Matchedid_',1:ncol(MatchedID))
write.csv(MatchedID,paste0(dir.PSM,'PSM_Matchedid_',addname,'_',Sys.Date(),'.csv'))

matchedDat<-dat[which(dat[,id]%in%c(as.vector(MatchedID),rownames(MatchedID))),]

tab1Adj <- CreateTableOne(vars = Covars, strata = GrpVar, data = matchedDat)
write.csv(as.data.frame(print(tab1Adj)),paste0(dir.PSM,'PSM_Table1_AfterM',addname,'_',Sys.Date(),'.csv'))

## table one will choose the proper test for you, however you cannot tell which is used...

sink(file = paste0(dir.PSM,'PSM_outputs',addname,'_',Sys.Date(),'.txt'))
print('Model summary')
print(summary(match.out.all$model)$coeff,digit=digit)

print(paste0(levels(dat[,GrpVar])[2],' = ','Treated;',levels(dat[,GrpVar])[1],' = ','Control'))

print('Matching summary')
print(summary(match.out.all,digit=digit))

print('Table 1 Before matching')
print(tab1Unadj, test = T, smd = T)## smd >0.1 is consider imbalanced
print('Table 1 After matching')
print(tab1Adj, test = T, smd = T)## smd >0.1 is consider imbalanced
sink()

} else if(length(levels(dat[,GrpVar]))==3){
  print('Check 2:Group variable is of 3 levels')
  
  fmla_tri <-as.formula(paste("~", paste(Covars, collapse= "+")))
  tpsa <- trips(dat, dat[,GrpVar], fmla_tri)
  tmatch <- trimatch(tpsa,caliper = 0.2, status=FALSE,method=Method)##default='maximumTreat'
  
  pdf(file=paste0(dir.PSM,'Plots_PSM_',addname,Sys.Date(),'.pdf'),height=6,width=6,onefile=T)
  p1<-balance.plot(tmatch,PD_pat[,Covars])
  p2<-multibalance.plot(tpsa, grid = TRUE, Covars)
  dev.off()
  
  MatchedID<-data.frame(tmatch)[,1:3]
  MatchedID<-apply(MatchedID,2,function(x){dat[as.numeric(x),id]})
  write.csv(MatchedID,paste0(dir.PSM,'PSM_Matchedid_trimatch_',addname,'_',Sys.Date(),'.csv'))
  
  matchedDat<-dat[which(dat[,id]%in%as.vector(MatchedID)),]
  
  tab1Aadj_tri <- CreateTableOne(vars = Covars, strata = GrpVar, data = matchedDat)
  write.csv(as.data.frame(print(tab1Aadj_tri)),paste0(dir.PSM,'PSM_Table1_AfterM',addname,'_',Sys.Date(),'.csv'))
  
  sink(file = paste0(dir.PSM,'PSM_matchedtriplets',addname,'_',Sys.Date(),'.txt'))

  print(tmatch)
  print('Table 1 Before matching')
  print(tab1Unadj, test = T, smd = T)## smd >0.1 is consider imbalanced
  print('Table 1 After matching')
  print(tab1Aadj_tri, test = T, smd = T)## 
  print('Number of observation per group under original data')
  print(table(dat[,GrpVar]))
  print('Number of observation per group under matched data')
  print(table(matchedDat[,GrpVar]))
  sink()
  
} else{
  stop ('Group variable is more than 3 levels')
}
  return(list(matchedDat=matchedDat,MatchedID=MatchedID)) ## ouput matched ID and newdata
  
}

