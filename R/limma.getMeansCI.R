limma.getMeansCI<-function(fit){
  M<-fit$coefficients; colnames(M)<-paste('M', colnames(M),'.')
  CI<-qt(0.975, fit$df.residual + fit$df.prior)*fit$stdev.unscaled *sqrt(fit$s2.post)
  colnames(CI)<-paste('CI', colnames(CI),'.')
  return(cbind(M,CI))
}

