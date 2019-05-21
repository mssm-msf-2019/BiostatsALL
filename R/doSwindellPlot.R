#' A function that produce Swindell plots showing the intersection between the first R ranked human and mouse list of genes (as in SF2014, Plos one)
# with CI. 
#' @description The R transcripts with expression most strongly elevated in human psoriasis were identified, along with the 5000 transcripts with expression most strongly decreased in human psoriasis. These transcripts were ranked according to the estimated fold-change expression ratio (psoriasis/control), with lower ranks assigned to transcripts most strongly increased or decreased in human psoriasis. For any rank N, where N = 1, …, 5000, we isolated the top N human transcripts and identified orthologous mouse transcripts, and then determined whether these mouse transcripts overlapped significantly with the top N mouse transcripts increased or decreased in the (A) K5-Tie2 phenotype, (B) IMQ phenotype, (C) K14-AREG phenotype on ear skin, (D) K14-AREG phenotype on tail skin, (E) K5-Stat3C phenotype and (F) K5-TGFβ1 phenotype. In each figure, the red line corresponds to the overlap, at a given rank N, between the top N psoriasiform-increased mouse transcripts and the set of mouse transcripts orthologous to the top N psoriasis-increased human transcripts. Similarly, the green line corresponds to the overlap, at a given rank N, between the top N psoriasiform-decreased mouse transcripts and the set of mouse transcripts orthologous to the top N psoriasis-decreased human transcripts. The grey region outlines the level of overlap expected by chance for any given rank N (i.e., the 95% confidence region of the null hypergeometric distribution). A significant level of overlap is indicated for each psoriasiform phenotype because red and green lines lie above the grey region that spans the null expectation.	
#' @param Res output of Run.swindel.sims
#' @param fname: file name to save the plot
#' @keywords swindell plots, confidence interval of intersection
#' @examples 
#' doSwindellPlot(Res, 'MAD3_IL23Mouse', title='Swindell Plot')


doSwindellPlot<-function(Res,fname='',title='',width=6.5,height=5){
   p <- ggplot(as.data.frame(Res), aes(x = rank, y = overlap, colour=Direction)) + 
   	geom_point() +
   	theme_bw() + theme(legend.position=c(0,0.8),legend.justification='left') +
    scale_colour_manual(values=c('red', 'blue')) +
  	geom_ribbon(aes(ymin = q5, ymax = q95, fill=Direction),alpha=0.6) +
  	theme_bw() + theme(legend.position=c(0,0.9),legend.justification='left') +
    labs(x="Ranks",y="Overlap",title=title) + ylim(c(0,1800))
  	p
  	ggsave(filename =paste(fname,".pdf",sep=''),width= width,height= height)
}
