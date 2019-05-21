normalize.by.patient<-function(x,pat){
  f<-try(lme(x~1,random=~1|pat))
  out<-x;
  if (!inherits(f, "try-error")) {f<-lme(x~1,random=~1|pat); out<-x-f$fitted[,2]}
  return(out)
}
