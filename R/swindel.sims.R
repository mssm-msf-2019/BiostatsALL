#' A function that runs simulations to produce Swindell plots to show intersection between firsth Rth ranked genes between human and mouse list of genes (as in SF2014, Plos one)
#' @description The R transcripts with expression most strongly elevated in human psoriasis were identified, along with the 5000 transcripts with expression most strongly decreased in human psoriasis. These transcripts were ranked according to the estimated fold-change expression ratio (psoriasis/control), with lower ranks assigned to transcripts most strongly increased or decreased in human psoriasis. For any rank N, where N = 1, …, 5000, we isolated the top N human transcripts and identified orthologous mouse transcripts, and then determined whether these mouse transcripts overlapped significantly with the top N mouse transcripts increased or decreased in the (A) K5-Tie2 phenotype, (B) IMQ phenotype, (C) K14-AREG phenotype on ear skin, (D) K14-AREG phenotype on tail skin, (E) K5-Stat3C phenotype and (F) K5-TGFβ1 phenotype. In each figure, the red line corresponds to the overlap, at a given rank N, between the top N psoriasiform-increased mouse transcripts and the set of mouse transcripts orthologous to the top N psoriasis-increased human transcripts. Similarly, the green line corresponds to the overlap, at a given rank N, between the top N psoriasiform-decreased mouse transcripts and the set of mouse transcripts orthologous to the top N psoriasis-decreased human transcripts. The grey region outlines the level of overlap expected by chance for any given rank N (i.e., the 95% confidence region of the null hypergeometric distribution). A significant level of overlap is indicated for each psoriasiform phenotype because red and green lines lie above the grey region that spans the null expectation.	
#' @param R length of the list, maximum rank.
#' @param hs: Human genes ranged by FCH
#' @param mm  Mouse genes ranged by FCH
#' @param Homol.db factor to separate columns in the heatmaps
#' @param B aglomeration strategy 
#' @keywords swindell plots, confidence interval of intersection
#' @examples 
#' swindel.sims<-function(5000,hs,mm,Homol.db,B=100)

swindel.sims<-function(R,hs,mm,Homol.db,B){
	   ######The 5000 transcripts with expression most strongly elevated in human psoriasis were identified, along with the 5000 transcripts with expression most strongly decreased in human psoriasis. These transcripts were ranked according to the estimated fold-change expression ratio (psoriasis/control), with lower ranks assigned to transcripts most strongly increased or decreased in human psoriasis. For any rank N, where N = 1, …, 5000, we isolated the top N human transcripts and identified orthologous mouse transcripts, and then determined whether these mouse transcripts overlapped significantly with the top N mouse transcripts increased or decreased in the (A) K5-Tie2 phenotype, (B) IMQ phenotype, (C) K14-AREG phenotype on ear skin, (D) K14-AREG phenotype on tail skin, (E) K5-Stat3C phenotype and (F) K5-TGFβ1 phenotype. In each figure, the red line corresponds to the overlap, at a given rank N, between the top N psoriasiform-increased mouse transcripts and the set of mouse transcripts orthologous to the top N psoriasis-increased human transcripts. Similarly, the green line corresponds to the overlap, at a given rank N, between the top N psoriasiform-decreased mouse transcripts and the set of mouse transcripts orthologous to the top N psoriasis-decreased human transcripts. The grey region outlines the level of overlap expected by chance for any given rank N (i.e., the 95% confidence region of the null hypergeometric distribution). A significant level of overlap is indicated for each psoriasiform phenotype because red and green lines lie above the grey region that spans the null expectation.	
	   
	    
	    hs.up.homol <-Homol.db$Ortholog.Probe.Set[which(Homol.db$Probe.Set.ID%in%hs[1:R])]
	    hs.down.homol <-Homol.db$Ortholog.Probe.Set[which(Homol.db$Probe.Set.ID%in%hs[length(hs):(length(hs)-R+1)])]
 
     	hs.homol<-unique(Homol.db$Ortholog.Probe.Set[which(Homol.db$Probe.Set.ID%in%hs)])
     	r<-unlist(mclapply(1:B, function(b, x, y, rR){sum(sample(x,size=rR)%in%sample(y,size=rR))}, mm, hs.homol, R))
     	  	
		oU<-length(intersect(hs.up.homol, mm[1:R]))
		oD<-length(intersect(hs.down.homol,mm[length(mm):(length(mm)-R+1)]))
		return(c(oUp=oU,oDown=oD,quantile(r,c(0.05,0.95))))
		}

