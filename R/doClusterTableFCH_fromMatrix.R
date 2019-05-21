#' A function that creates a heatmap with column indicators on top of each column.
#' @description This function creates a nice heatmap for publication. It puts on top of each column a color indicator of what sample type it is, does clustering, and adds a table in the right margin of the figure with logFoldChanges and respective significance indicators (stars).
#' @param mat numeric matrix to cluster. Usually th raw expression data.
#' @param coefs estimates to be drawn tin the tabl (eg lgFCHs)
#' @param fdrs p values or fdr indicating significance of the estimates; should be the same dim as coefs
#' @param colfc factor to add colored legend to the heatmap
#' @param facsepr factor to separate columns in the heatmaps
#' @param hmeth aglomeration strategy
#' @param dmeth aglomeration distance to cluster
#' @param ColD logical TRUE if cluster by rows and columns
#' @param sn edited names for the rows (if different than the rownames) eg symbols for probesets
#' @param gn gene name
#' @param scl how to scale the matrix (row, column, none, both)
#' @keywords heatmap , expression vizualisation
#' @examples
#' doClusterTableFCH_fromMatrix(mat,ps=rownames(mat),colfac,facsepr=NULL,hmeth='average',dmeth=cor.dist,coefs,fdrs,main="", ColD=FALSE, ss=ps, gn=NULL, breaks=NULL,cexg=1,margins=c(5,20) )

doClusterTableFCH_fromMatrix<-function (mat, ps = rownames(mat), colfac, facsepr = NULL, hmeth = "average",
                                        dmeth = cor.dist, coefs, fdrs, main = "", ColD = FALSE, ss = ps,
                                        gn = NULL, breaks = NULL, cexg = 1, margins = c(5, 20), scl = "row",
                                        p.cuts.stars=c(0.01, 0.05, 0.1), p.symbols.stars= c("**", "*", "+")) {
  #mat: numeric matrix to cluster
  #coefs: estimates to be drawn tin the table (eg lgFCHs)
  #fdrs: p values or fdr indicating significance of the estimates; should be the same dim as coefs
  #colfc: factor to add colored legend to the heatmap
  #facsepr: factor to separate columns in the heatmaps
  #hmeth, dmeth: aglomeration strategy and distance to cluster
  #ColD; logical TRUE if cluster by rows and columns
  #sn edited names for the rows (if different than the rownames) eg symbols for probesets amed vector with probesets as names
  # gn: gene name, named vector with probesets as names

  require(stringr)
  require(weights)

  mat <- mat[ps, ]
  ss <- ss[ps]
  fdrs <- fdrs[ps, ]
  coefs <- coefs[ps, ]
  if (!is.null(gn))
    gn <- gn[ps]
  maxss <- max(sapply(ss, nchar))
  adspace <- function(x, n) {
    paste(x, substr(" ---------------------", 1, n - nchar(x)),
          sep = "")
  }
  ss <- paste(sapply(ss, adspace, maxss + 3), ":  ", sep = "")
  ADDzero <- function(x) {
    if (grepl(pattern = "\\.", x = x) == TRUE) {
      a1 <- unlist(strsplit(x, "\\."))[1]
      a2 <- unlist(strsplit(x, "\\."))[2]
      if (is.na(a2) == T) {
        a2 <- "0"
      }
      b1 <- paste(str_dup(" ", 4 - nchar(a1)), a1, sep = "")
      b1 <- a1
      b2 <- paste(a2, str_dup("0", 2 - nchar(a2)), sep = "")
      out <- paste(b1, b2, sep = ".")
    }
    else {
      out <- paste(paste(x, ".00", sep = ""), sep = "")
    }
    return(out)
  }
  adzero <- function(xv) {
    out <- sapply(xv, ADDzero)
    nmax <- max(sapply(out, nchar))
    sapply(out, function(x, nmx) {
      paste(str_dup(" ", nmx - nchar(x)), x, sep = "")
    }, nmax)
  }
  reformatps <- function(p) {
    as.numeric(format(p, digit = 3, drop0trailing = TRUE))
  }
  transformfch <- function(lgfch) {
    fch <- sign(lgfch) * (2^abs(lgfch))
    return(fch)
  }
  coef1 <- apply(coefs[ps, ], 2, function(x) {
    adzero(as.character(signif(transformfch(x), 3)))
  })
  rownames(coef1) <- ps
  fdrs1 <- apply(fdrs[ps, ], 2, function(x) { ifelse(x < 1e-04, signif(x, 1), round(x, 4)) })  ## MSf deleted this line because uits not needed
  #fdrs2 <- apply(fdrs1, 2, starmaker, p.levels = c(0.01, 0.05, 0.1), symbols = c("**", "*", "+"))
  fdrs2 <- apply(fdrs, 2, starmaker, p.levels =p.cuts.stars, symbols = p.symbols.stars)

  rownames(fdrs2) <- ps
  adSpace <- function(x) {
    out <- paste(x, str_dup(" ", 4 - nchar(x)), sep = "")
  }
  fdrs2 <- apply(fdrs2, 2, adSpace)
  coef2 <- coef1
  print(tail(coef2))
  print(tail(fdrs2))
  a <- DoHeatmap(mat, colfac = colfac, symb = ss, dmeth = dmeth,
                 hmeth = hmeth, cex.genes = cexg, ColD = ColD, main = main,
                 margins = c(5, 10), breaks = breaks, scl = scl)
  Tab <- data.frame(Symbol = ss[a$rowInd])
  if (!is.null(gn))
    Tab <- cbind(Tab, Desc = substr(gn[a$rowInd], 1, 40))
  for (i in c(1:ncol(coef2))) {
    Tab <- cbind(Tab, paste(coef2[a$rowInd, i], fdrs2[a$rowInd,
                                                      i], sep = ""))
  }
  Tab <- print.table(Tab)
  ssTab <- apply(Tab[, ], 1, function(x) {
    stackchar(x, sep = "")
  })
  mat2 <- mat[a$rowInd, ]
  rownames(mat2) <- ssTab
  ssTab_bold <- do.call(expression, sapply(as.character(ssTab),
                                           function(.x) {
                                             substitute(bold(.x), list(.x = .x))
                                           }))
  par(family = "mono")
  a <- DoHeatmap(mat2, colfac = colfac, facsepr = facsepr,
                 symb = ssTab_bold, dmeth = dmeth, hmeth = hmeth, cex.genes = cexg,
                 ColD = ColD, main = main, margins = margins, breaks = breaks,
                 scl = scl)
}

