#' A function that produce Swindell plots showing the intersection between the first R ranked human and mouse list of genes (as in SF2014, Plos one)
# with CI. It runs this function for R=1:maxRank and then plots the results over R.#
#' @description The R transcripts with expression most strongly elevated in human psoriasis were identified, along with the 5000 transcripts with expression most strongly decreased in human psoriasis. These transcripts were ranked according to the estimated fold-change expression ratio (psoriasis/control), with lower ranks assigned to transcripts most strongly increased or decreased in human psoriasis. For any rank N, where N = 1, …, 5000, we isolated the top N human transcripts and identified orthologous mouse transcripts, and then determined whether these mouse transcripts overlapped significantly with the top N mouse transcripts increased or decreased in the (A) K5-Tie2 phenotype, (B) IMQ phenotype, (C) K14-AREG phenotype on ear skin, (D) K14-AREG phenotype on tail skin, (E) K5-Stat3C phenotype and (F) K5-TGFβ1 phenotype. In each figure, the red line corresponds to the overlap, at a given rank N, between the top N psoriasiform-increased mouse transcripts and the set of mouse transcripts orthologous to the top N psoriasis-increased human transcripts. Similarly, the green line corresponds to the overlap, at a given rank N, between the top N psoriasiform-decreased mouse transcripts and the set of mouse transcripts orthologous to the top N psoriasis-decreased human transcripts. The grey region outlines the level of overlap expected by chance for any given rank N (i.e., the 95% confidence region of the null hypergeometric distribution). A significant level of overlap is indicated for each psoriasiform phenotype because red and green lines lie above the grey region that spans the null expectation.	
#' @param maxRank, teh maximum rank where the intersectionof the two list will be calculate.
#' @param B: number of siluations to runs for each fixed rank
#' @param rankedList1 ranked List1 
#' @param HomologyMatrix matrix to match genes between List1 and List2 
#' @param B aglomeration strategy 
#' @keywords swindell plots, confidence interval of intersection
#' @examples 
#' Run.swindel.sims(5000,B=100,rankedPsoriasis,rankedMouse, SwindelHomol, 'MAD3', 'IL23Mouse')


Run.swindel.sims<-function(maxRank,B,rankedList1, rankedList2, HomologyMatrix, nameList1, nameList2){
	library(ggplot2)
	TTT<-mclapply(1:maxRank, swindel.sims, rankedList1, rankedList2, HomologyMatrix,B)
	Res<-t(sapply(TTT,function(x){x}))
	bh3<-data.frame(rank=rep(1:nrow(Res),2),overlap=c(Res[,1],Res[,2]),Direction=gl(2,nrow(Res), 
	labels=c('Up-regulated','Down-regulated')),q5=c(Res[,3],Res[,3]),q95=c(Res[,4],Res[,4]))
	doSwindellPlot(bh3,fname=paste('IntersectionPlot', nameList1, nameList2,sep='_'),title=paste(nameList1, nameList2,sep='-'))
	return(bh3)
#
}
