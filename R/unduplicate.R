#' Create a new expression matrix with means of technical replicates
#' @examples
#' gset <- getGEO('GSE54839', GSEMatrix =TRUE) #GSE61276 GSE12654
#' exprset <- gset[[1]]
#' expr.matrix <- exprs(exprset)
#' comparisons <- comparisonsPheno(exprset,NULL)[[1]]
#' comparisonsTable <- comparisonsPheno(exprset,NULL)[[2]]
#' selectedReplicates <- "source_name_ch1"
#' @export
unduplicate <- function(comparisonsTable,selectedReplicates,expr.matrix){
  replicates <- as.character(comparisonsTable[,selectedReplicates])
  uniqueReplicates <- unique(replicates)
  if(length(replicates) > length(uniqueReplicates)){
    expr.matrix2 <- matrix(nrow = nrow(expr.matrix),ncol = length(uniqueReplicates))
    comparisonsTable2 <- as.data.frame(matrix(nrow = length(uniqueReplicates),ncol = ncol(comparisonsTable)))
    for (i in 1:length(uniqueReplicates)){
      replicated <- grep(uniqueReplicates[i],replicates)
      expr.matrix2[,i] <- apply(expr.matrix[,replicated],1,mean)
      if(i == 1){
        comparisonsTable2 <- comparisonsTable[replicated[1],]
      } else {
        comparisonsTable2 <- rbind(comparisonsTable2,comparisonsTable[replicated[1],])
      }
    }
    colnames(expr.matrix2) <- uniqueReplicates
    rownames(expr.matrix2) <- rownames(expr.matrix)
    return(list(expr.matrix2,comparisonsTable2))
  } else {
    return(list(expr.matrix,comparisonsTable))
  }
}


