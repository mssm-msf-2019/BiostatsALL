#' A function that prepare the results of a linear model fit
#' @description This function make a Resume of a microarray linear model fit results after eBayes function
#' @param ebfit output of eBayes
#' @param ANNPkg annotation package 
#' @param Decision output of decideTests
#' @param coefs vector of character strings that contains the coeficients to be included
#' @param adjmeth character string specifying p-value adjustment method. 
#'        Possible values are "none", "BH", "fdr" (equivalent to "BH"), "BY" and "holm". 
#' @param method character string specify how probes and contrasts are to be combined in the multiple testing strategy. 
#'        Choices are "separate", "global", "hierarchical", "nestedF" or any partial string.
#' @param mcut numeric value, folchange of cut
#' @param fname
#' @param pvalue critical pvalue
#' @param pPath numeric value, ussually number of pathways
#' @param pGO numeric value
#' @param nPath number of Pathways
#' @param ann_ind logical value equal TRUE by default
#' @param exprset expression set
#' @param which.anncol numeric vector with annotation columns numbers

Resume <- function (ebfit, AnnPkg, Decision = NULL, coefs = c(1:ncol(ebfit$coef)), 
                    adjmeth = "fdr", method = "separate", mcut = 2, fname = "", 
                    pvalue = 0.05, pPath = 0.05, pGO = 0.05, nPath = 2, ann_ind = TRUE, 
                    exprset = NULL, which.anncol = c(1:3, 6:7, 9, 11, 12)) 
{
  IDS <- rownames(ebfit$coef)
  coefname <- colnames(ebfit$coef)[coefs]
  lgFCH <- ebfit$coef[, coefs]
  
  if (length(coefs) == 1) 
    lgFCH <- matrix(ebfit$coef[, coefs], dimnames = list(IDS, coefname))
  colnames(lgFCH) <- paste("logFCH", coefname, sep = "_")
  FCH <- apply(lgFCH, 2, function(x) {sign(x) * 2^(abs(x))
  })
  colnames(FCH) <- paste("FCH", coefname, sep = "_")
  print(coefs)
  if (is.null(Decision)) 
    Decision <- decideTests(ebfit[, coefs], method = method, 
                            adjust.method = adjmeth, p.value = pvalue, lfc = log2(mcut))
  else Decision <- Decision[, coefs]
  if (length(coefs) == 1) {
    Decision <- matrix(Decision, dimnames = list(IDS, coefname))
  }
  Diff <- (Decision != 0)
  Up <- (Decision == 1)
  Down <- (Decision == -1)
  Status <- Decision
  colnames(Status) <- paste("Status", coefname, sep = "_")
  print("Calculate ps")
  reformatps <- function(p) {
    as.numeric(format(p, digit = 3, drop0trailing = TRUE))
  }
  pval <- ebfit$p.value[, coefs]
  if (length(coefs) == 1) {
    pval <- matrix(pval, dimnames = list(IDS, coefname))
  }
  padj <- apply(pval, 2, function(p, adjmeth) {
    reformatps(p.adjust(p, method = adjmeth))
  }, adjmeth)
  if (length(coefs) == 1) {
    padj <- matrix(padj, dimnames = list(IDS, coefname))
  }
  colnames(pval) <- paste("p", coefname, sep = "_")
  dimnames(padj) <- list(IDS, paste("FDR", coefname, sep = "_"))
  ps <- ifelse((adjmeth == "none"), "p", "fdr")
  fname = paste(fname, ps, substr(as.character(pvalue), 3, 
                                  4), "fch", as.character(mcut), sep = "")
  if ((length(coefs) <= 3) & (length(coefs) > 1)) {
    pdf(file = paste(fname, "VennDiagram", ".pdf", sep = ""), 
        width = 10, onefile = TRUE)
    vennDiagram(Decision[, coefs], include = c("up", "down"), 
                counts.col = c("red", "green"))
    dev.off()
  }
  print("hhhh")
  for (i in c(1:length(coefs))) {
    cni <- coefname[i]
    print(cni)
    Table <- cbind(FCH = round(FCH[, i], 3), lgFCH = round(lgFCH[, 
                                                                 i], 3), P = pval[, i], fdr = padj[, i])
    Tit = paste(fname, ":", cni, "differentially expressed genes", 
                "(p<", as.character(pvalue), ")")
    print("j")
    Tablex <- cbind(Table, Status = Decision[, i])
    colnames(Tablex) <- paste(colnames(Tablex), cni, sep = "-")
    if ((ann_ind) & (length(IDS[(Diff[, i] > 0)]))) {
      sig.ps <- IDS[Diff[, i]]
      Annotat.Tab(AnnPkg, probesids = sig.ps, fname = paste(fname, 
                                                            cni, "_Diff", sep = ""), Tit = Tit, cord = 1, 
                  dec = TRUE, adTable = Table[sig.ps, ], which = which.anncol, 
                  savetxt = TRUE, exprset = exprset)
      print("b")
      if (!is.null(pPath)) {
        writefcn(x = AnalyzePathways(IDS[Diff[, i]], 
                                     AnnPkg, p = pPath), fname = paste("Pathways", 
                                                                       fname, cni, "_Diff.txt", sep = ""))
        writefcn(x = AnalyzePathways(IDS[Up[, i]], AnnPkg, 
                                     p = pPath), fname = paste("Pathways", 
                                                               fname, cni, "_UP.txt", sep = ""))
        writefcn(x = AnalyzePathways(IDS[Down[, i]], 
                                     AnnPkg, p = pPath), fname = paste("Pathways", 
                                                                       fname, cni, "_Down.txt", sep = ""))
      }
      if (!is.null(pGO)) {
        write.table(x = AnalyzeGOall(IDS[Diff[, i]], 
                                     AnnPkg, p = pGO, cond = TRUE), file = paste("GOs", 
                                                                                 fname, cni, "_Diff.txt", sep = ""), sep = "\t", 
                    quote = FALSE, row.names = FALSE)
        write.table(x = AnalyzeGOall(IDS[Up[, i]], AnnPkg, 
                                     p = pGO, cond = TRUE), file = paste("GOs", 
                                                                         fname, cni, "_Up.txt", sep = ""), sep = "\t", 
                    quote = FALSE, row.names = FALSE)
        write.table(x = AnalyzeGOall(IDS[Down[, i]], 
                                     AnnPkg, p = pGO, cond = TRUE), file = paste("GOs", 
                                                                                 fname, cni, "_Down.txt", sep = ""), sep = "\t", 
                    quote = FALSE, row.names = FALSE)
      }
    }
  }
  print("ann")
  BigTable <- cbind(lgFCH, FCH, pval, padj, Status)
  if (length(coefs) > 1) {
    nodiff <- (rowSums((Decision == 0)) == length(coef))
    somediff <- (rowSums(Diff) > 0)
    padj.F <- p.adjust(ebfit$F.p.value, method = adjmeth)
    BigTable <- cbind(lgFCH, FCH, pval, padj, Status)
    r <- NULL
    for (i in (1:ncol(lgFCH))) {
      r <- c(r, seq(i, ncol(BigTable), ncol(lgFCH)))
    }
    BigTable <- cbind(BigTable[, r], P_Ftest = ebfit$F.p.value, 
                      FDR_Ftest = padj.F)
    rownames(BigTable) <- IDS
    Annotat.Tab(AnnPkg, probesids = IDS[somediff], fname = paste(fname, 
                                                                 "SomeDiff", sep = ""), Tit = Tit, cord = 3, dec = FALSE, 
                adTable = BigTable[somediff, ], savetxt = TRUE, exprset = exprset, 
                which = which.anncol)
  }
  write.csv(x = cbind(Symbol = getSYMBOL(IDS, AnnPkg), Desc = unlist(lookUp(IDS, 
                                                                            AnnPkg, "GENENAME")), BigTable), file = paste(fname, 
                                                                                                                          "AllGenes.csv", sep = ""))
  return(invisible(Decision))
}