#' Create a vector of the group membership for each sample
#' @author Simon J Pelletier
#' @param namesTable data.frame with a single column containing group names
#' @description Takes a data.frame with a single column containing strings, removes the whitespace
#' of all strings and returns it as a vector
#' @return A vector of the group membership for each sample.The name of each
#' element is the initial name of the sample
#' @examples
#' gset <- getGEO('GSE61276', GSEMatrix =TRUE) #GSE61276 GSE12654
#' exprset <- gset[[1]]
#' comparisons <- comparisonsPheno(exprset)[[1]]
#' comparisonsTable <- comparisonsPheno(exprset)[[2]]
#' namesTable <- comparisonsTable
#' names2 <- namesSelectedComparisons(namesTable)
#' @export
namesSelectedComparisons = function(namesTable){
  namesSelected <- apply(namesTable,1,function(x){
    y<-lapply(x,function(z){
      paste(strsplit(z," ")[[1]],collapse=".")
    })
    y<-lapply(y,function(z){
      paste(strsplit(z,":")[[1]],collapse="")
    })
    paste(y,collapse="_")
  })
  names(namesSelected) = rownames(namesTable)
  return(namesSelected)
}
