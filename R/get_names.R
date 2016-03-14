#' Get names of each sample and group names
#' @author Simon J Pelletier
#' @param x Matrix or vector
#' column name has to be separeted
#' @return The name of each groups.
#' @examples
#' EXPRESSION MATRIX
#' expr.matrix <- readRDS("data/expr_matrix_LGVD.rds")
#' get_names(expr.matrix)
#'
#' EXPRESSION SET
#' gset <- getGEO("GSE12654", GSEMatrix =TRUE)
#' exprset <- gset[[1]]
#' names1 <- get_names(sampleNames(exprset))
#' @keywords names names.unique
#' @export
get_names = function(x){
  if (is.matrix(x)) names <- colnames(expr.matrix)
  else names <- x
  names1 = rep("",length(names))
  for (i in 1:length(names)){
    nameSplit = strsplit(names[i],"_")
    for (j in 1:(length(nameSplit[[1]]))){
      if(j==1) names1[i] = nameSplit[[1]][j]
      else if(length(grep("([0-9]+).*$", nameSplit[[1]][j])) == 0) names1[i] = paste0(names1[i],"_",nameSplit[[1]][j])
    }
  }
  return(names1)
}
