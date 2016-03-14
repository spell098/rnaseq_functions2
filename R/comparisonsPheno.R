#' Transform intensity values in a gradient of color from green to red.
#' @author Simon J Pelletier
#' @import GEOquery
#' @param exprset Expression set containing all the information on a dataset
#' @return A vector of colors. The darkest green values are the lowest and the darkest red values are the highest.
#' @keywords comparisons
#' @seealso
#' \code{\link[GEOquery]{getGEO}}
#' \code{\link[Biobase]{ExpressionSet}}
#' @examples
#' gset <- getGEO('GSE54839', GSEMatrix =TRUE) #GSE61276 GSE12654
#' exprset <- gset[[1]]
#' comparisons <- comparisonsPheno(exprset)[[1]]
#' comparisonsTable <- comparisonsPheno(exprset)[[2]]
#' @export
comparisonsPheno = function(exprset){
  comparisons_table <- pData(exprset)
  comparisons = apply(comparisons_table,2,function(x){
    paste(unique(x),collapse="\\")
  })
  comparisons2 = NULL
  j<-1
  for(i in 1:length(comparisons)){
    if(length(strsplit(comparisons[i],"\\\\")[[1]]) > 1 & length(strsplit(comparisons[i],"\\\\")[[1]]) < ncol(exprs(exprset))){
      comparisons2[j] = comparisons[i]
      names(comparisons2)[j] = names(comparisons)[i]
      j<-j+1
    } else {
      comparisons_table = comparisons_table[,-j]
    }
  }
  comparisons2 = comparisons2[!duplicated(comparisons2)]
  comparisons_table = comparisons_table[!duplicated(lapply(comparisons_table,summary))]
  #write.csv(replicates,"testReplicates.txt")
  #if(!is.null(replicates) && replicates != ""){
    #write.csv(replicates,paste0("testReplicates","_notEmpty.txt"))
  #  x <- !duplicated(as.character(comparisons_table[,replicates]))
  #  comparisons_table <- comparisons_table[x,]
  #}
  #write.table(comparisons_table,"test1Comparisons_table.txt")
  return(list(comparisons2,comparisons_table))
}
