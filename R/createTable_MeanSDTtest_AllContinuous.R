#' A function to create descriptive table with Mean/SD/T-test Pvalue for continuous data
#' TablePeds_ADvsC<-createNiceTable_diagnosis(Z, group.var='Diagnosis')

createTable_MeanSDTtest_AllContinuous<-function(db,group.var){
  callt<-function(x,g,veq){printPVal(t.test(x~g, var.equal = veq)$p.value)}
  printMeanSDN <-function (vec){
    sprintf("%.1f (%.1f) ,n=%.0f", mean(vec, na.rm = TRUE), sd(vec, na.rm = TRUE), length(vec[!is.na(vec)]))
  }
  db<-cbind(db,group=db[,group.var])
  stats=t(ddply(db,.(group),numcolwise(printMeanSDN)))
  p.hom<-numcolwise(callt, g=db$group,veq=T)(db)
  p.het<-numcolwise(callt, g=db$group,veq=F)(db)
  Tab<-merge(cbind.data.frame(Outcome=rownames(a),a),
             cbind.data.frame(Outcome=colnames(p.hom),p.hom=t(p.hom),p.het=t(p.het)),
             by='Outcome')
}


