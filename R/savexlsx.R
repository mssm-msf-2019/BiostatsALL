#' A function that saves data frames in different spreadsheets in excel
#' @description a function that saves
#' @param xs a matrix or a data frame with numeric columns
#' @param ys a vector with the same dimension as the number of rows in x

#' @examples
#' EvaluatePerformance(L)


savexlsx <- function (file, ...)
{
  require(xlsx, quietly = TRUE)
  objects <- list(...)
  fargs <- as.list(match.call(expand.dots = TRUE))
  objnames <- as.character(fargs)[-c(1, 2)]
  nobjects <- length(objects)
  for (i in 1:nobjects) {
    if (i == 1)
      write.xlsx(objects[[i]], file, sheetName = objnames[i])
    else write.xlsx(objects[[i]], file, sheetName = objnames[i],
                    append = TRUE)
  }
  print(paste("Workbook", file, "has", nobjects, "worksheets."))
}
