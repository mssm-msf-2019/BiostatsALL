#' Function that runs the original Combat script but returns different parts of the internal model instead of just the adjusted data. It's meant to be used when the design is complex and not fullly modeled,
#'  so one can get the estimated coeficients and play with it. I sed it in the IBM competition
#' @description Inputs and Outputs are the same as in ComBat.R
#' @param exprs: matrix of expression
#' @param g: factor defining the batches
#' #' @examples 
#' abatch=ReadAffy()
#' abatch$Batch<-sapply(protocolData(abatch)$ScanDate,function(x){substr(x,1,10)}
#' eset<-rma(abatch)
#' out<-RunCombat_MSF(exprs(eset),eset$Batch)
#' 

ComBat_MSF<-function(expression_xls, sample_info_file, type='txt', write=T, covariates='all', par.prior=T, filter=F, skip=0, prior.plots=T){
  #debug: expression_xls='exp.txt'; sample_info_file='sam.txt'; type='txt'; write=T; covariates='all'; par.prior=T; filter=F; skip=0; prior.plots=T
  cat('Reading Sample Information File\n')
  saminfo <- read.table(sample_info_file, header=T, sep='\t',comment.char='')
  if(sum(colnames(saminfo)=="Batch")!=1){return('ERROR: Sample Information File does not have a Batch column!')}
  
  cat('Reading Expression Data File\n')
  if(type=='csv'){
    dat <- read.csv(expression_xls,header=T,as.is=T)
    #print(dat[1:2,])
    #	dat <- dat[,trim.dat(dat)]  
    #print(colnames(dat))
    colnames(dat)=scan(expression_xls,what='character',nlines=1,sep=',',quiet=T)[1:ncol(dat)]
    #print(colnames(dat))
  } else{
    dat <- read.table(expression_xls,header=T,comment.char='',fill=T,sep='\t', as.is=T)
    dat <- dat[,trim.dat(dat)]
    colnames(dat)=scan(expression_xls,what='character',nlines=1,sep='\t',quiet=T)[1:ncol(dat)]
  }
  if (skip>0){geneinfo <- as.matrix(dat[,1:skip])
              dat <- dat[,-c(1:skip)]
  } else{geneinfo=NULL}
  #print(geneinfo[1:4])
  #print(dat[1:2,])
  
  if(filter){
    ngenes <- nrow(dat)
    col <- ncol(dat)/2
    #present <- apply(dat, 1, filter.absent, filter)
    SDs<-apply(dat,1,sd)
    present<-(SDs>0.1)
    dat <- dat[present, -(2*(1:col))]
    if (skip>0){geneinfo <- geneinfo[present,]}
    cat('Filtered genes absent in more than',filter,'of samples. Genes remaining:',nrow(dat),'; Genes filtered:',ngenes-nrow(dat),'\n')
  }
  
  if(any(apply(dat,2,mode)!='numeric')){return('ERROR: Array expression columns contain non-numeric values! (Check your .xls file for non-numeric values and if this is not the problem, make a .csv file and use the type=csv option)')}
  
  tmp <- match(colnames(dat),saminfo[,1])
  if(any(is.na(tmp))){return('ERROR: Sample Information File and Data Array Names are not the same!')}
  tmp1 <- match(saminfo[,1],colnames(dat))
  saminfo <- saminfo[tmp1[!is.na(tmp1)],]		
  
  if(any(covariates != 'all')){saminfo <- saminfo[,c(1:2,covariates)]}
  design <- design.mat(saminfo)	
  print(saminfo)
  print(design)
  batches <- list.batch(saminfo)
  n.batch <- length(batches)
  n.batches <- sapply(batches, length)
  n.array <- sum(n.batches)
  
  print(sum(is.na(dat)))
  
  ## Check for missing values
  NAs = any(is.na(dat))
  if(NAs){cat(c('Found',sum(is.na(dat)),'Missing Data Values\n'),sep=' ')}
  #print(dat[1:2,])
  ##Standardize Data across genes
  cat('Standardizing Data across genes\n')
  if (!NAs){B.hat <- solve(t(design)%*%design)%*%t(design)%*%t(as.matrix(dat))}else{B.hat=apply(dat,1,Beta.NA,design)} #Standarization Model
  grand.mean <- t(n.batches/n.array)%*%B.hat[1:n.batch,]
  if (!NAs){var.pooled <- ((dat-t(design%*%B.hat))^2)%*%rep(1/n.array,n.array)}else{var.pooled <- apply(dat-t(design%*%B.hat),1,var,na.rm=T)}
  
  stand.mean <- t(grand.mean)%*%t(rep(1,n.array))
  if(!is.null(design)){tmp <- design;tmp[,c(1:n.batch)] <- 0;stand.mean <- stand.mean+t(tmp%*%B.hat)}	
  s.data <- (dat-stand.mean)/(sqrt(var.pooled)%*%t(rep(1,n.array)))
  
  ##Get regression batch effect parameters
  cat("Fitting L/S model and finding priors\n")
  batch.design <- design[,1:n.batch]
  
  if (!NAs){gamma.hat <- solve(t(batch.design)%*%batch.design)%*%t(batch.design)%*%t(as.matrix(s.data))
  } else {gamma.hat=apply(s.data,1,Beta.NA,batch.design)}
  delta.hat <- NULL
  
  for (i in batches){
    delta.hat <- rbind(delta.hat,apply(s.data[,i], 1, var,na.rm=T))
  }
  
  ##Find Priors
  gamma.bar <- apply(gamma.hat, 1, mean)
  t2 <- apply(gamma.hat, 1, var)
  a.prior <- apply(delta.hat, 1, aprior)
  b.prior <- apply(delta.hat, 1, bprior)
  
  
  ##Plot empirical and parametric priors
  
  if (prior.plots & par.prior){
    par(mfrow=c(2,2))
    tmp <- density(gamma.hat[1,])
    plot(tmp,  type='l', main="Density Plot")
    xx <- seq(min(tmp$x), max(tmp$x), length=100)
    lines(xx,dnorm(xx,gamma.bar[1],sqrt(t2[1])), col=2)
    qqnorm(gamma.hat[1,])	
    qqline(gamma.hat[1,], col=2)	
    
    tmp <- density(delta.hat[1,])
    invgam <- 1/rgamma(ncol(delta.hat),a.prior[1],b.prior[1])
    tmp1 <- density(invgam)
    plot(tmp,  typ='l', main="Density Plot", ylim=c(0,max(tmp$y,tmp1$y)))
    lines(tmp1, col=2)
    qqplot(delta.hat[1,], invgam, xlab="Sample Quantiles", ylab='Theoretical Quantiles')	
    lines(c(0,max(invgam)),c(0,max(invgam)),col=2)	
    title('Q-Q Plot')
  }
  
  ##Find EB batch adjustments
  print('No')
  gamma.star <- delta.star <- NULL
  if(par.prior){
    cat("Finding parametric adjustments\n")
    for (i in 1:n.batch){
      temp <- it.sol(s.data[,batches[[i]]],gamma.hat[i,],delta.hat[i,],gamma.bar[i],t2[i],a.prior[i],b.prior[i])
      gamma.star <- rbind(gamma.star,temp[1,])
      delta.star <- rbind(delta.star,temp[2,])
    }
  }else{
    cat("Finding nonparametric adjustments\n")
    for (i in 1:n.batch){
      temp <- int.eprior(as.matrix(s.data[,batches[[i]]]),gamma.hat[i,],delta.hat[i,])
      gamma.star <- rbind(gamma.star,temp[1,])
      delta.star <- rbind(delta.star,temp[2,])
    }
  }
  
  
  ### Normalize the Data ###
  cat("Adjusting the Data\n")
  
  bayesdata <- s.data
  j <- 1
  print(batch.design)
  for (i in batches){
    bayesdata[,i] <- (bayesdata[,i]-t(batch.design[i,]%*%gamma.star))/(sqrt(delta.star[j,])%*%t(rep(1,n.batches[j])))
    j <- j+1
  }
  
  bayesdata <- (bayesdata*(sqrt(var.pooled)%*%t(rep(1,n.array))))+stand.mean
  return(list(probesets=geneinfo,exprs=bayesdata, s.data=s.data, design.matrix=batch.design, gamma=gamma.star, delta=delta.star, batches=batches, var.pooled=var.pooled,stand.mean=stand.mean))
}