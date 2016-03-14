#' Regroup
#' @examples
#' gset <- getGEO("GSE54839", GSEMatrix =TRUE,getGPL = FALSE) #GSE61276 GSE12654 GSE54839 Methylation: GSE18111
#'
#' exprset <- gset[[1]]
#' expr.matrix <- exprs(exprset)
#' names1 <- get_names(sampleNames(exprset))
#' groupBy <- "group"
#' #expr.matrix <- readRDS("data/expr_matrix_LGVD.rds")
#' noCluster <- 3
#' fileContent <- "Expression Set"
#' names1 <- get_names(expr.matrix)
#' selectedVariables <- "characteristics_ch1"
#' selectedReplicates <- "source_name_ch1"
#' comparisonsTable <- comparisonsPheno(exprset)[[2]]
#' expr.matrix2 <- unduplicate(comparisonsTable,selectedReplicates,expr.matrix)

#' namesTable <- comparisonsTable
#' names2 <- namesSelectedComparisons(namesTable)
#' @export
regroupingSamples <- function(fileContent,groupBy,expr.matrix,selectedVariables, noCluster=NULL,names1,comparisonsTable=NULL,exprset=NULL){
  library(wesanderson)
  names2 <- NULL
  if (groupBy == "fuzzy"){
    clusters1 = fuzzy_clustering(expr.matrix,noCluster)[[1]]
    color = unlist(lapply(clusters1,function(x){wes()[x]}))
    namesColor = paste0("a",clusters1)
    comparisons = get_comparisons(namesColor)
    clusterNames = numberNames(namesColor)

    selectedVariables = comparisons[1]
    names2 = conditionsChoice(selectedVariables,namesColor)
  } else {
    if(fileContent == "Expression Matrix"){
      if(is.null(selectedVariables)){
        names <- names1
      } else {
        names = conditionsChoice(selectedVariables,names1)
      }
      names2 <- names
      namesColor = names1
      clusterNames = names2
      comparisons <- get_comparisons(names1)
    } else if(fileContent == "Expression Set"){
      comparisons <- as.data.frame(comparisonsPheno(exprset)[[1]]) #PROBLEM HERE; not good column name
      #colnames(comparisonsTable) <- comparisons
      namesTable <- as.data.frame(comparisonsTable[,selectedVariables])
      colnames(namesTable) <- selectedVariables
      rownames(namesTable) <- rownames(comparisonsTable)
      names <- namesSelectedComparisons(namesTable)
      for(i in 1:length(names)){
        names2[i] <- strsplit(gsub("[^[:alnum:]_ ]", ".", names[i]), " +")[[1]]
      }
      namesColor = names2
      clusterNames = names2
    }
    write.table(names,"names.txt")
    color <- color_groups(names)
  }
  return(list(names2,comparisons,color,namesColor))
}

