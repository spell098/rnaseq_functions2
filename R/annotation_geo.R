#' Get the annotations used in a NCBI GEO dataset
#' @import GEOquery
#' @import annotate
#' @author Simon J Pelletier
#' @param geoAccession The platform ID
#' @return A table with all the annotation for a platform
#' @examples
#'
#' x=annotation_geo('GPL10558') #Illumina HumanHT-12 V4.0 expression beadchip
#' @keywords geo annotation
#' @seealso
#' \code{\link[GEOquery]{getGEO}} ,
#' \code{\link[GEOquery]{Meta}} ,
#' \code{\link[GEOquery]{Table}} ,
#' \code{\link[annotate]{readGEOAnn}}
#' @export
annotation_geo = function(geoAccession){
  annotations = getGEO(geoAccession)
  symbol_column = NULL
  entrez_column = NULL
  iter = 0
  annotation_table <- Table(dataTable(annotations))
  symbol_column <- grep("symbol",colnames(annotation_table),ignore.case=TRUE)
  entrez_column <- grep("entrez",colnames(annotation_table),ignore.case=TRUE)
  iter=iter+1
  if(length(symbol_column) == 0 | length(entrez_column) == 0){
    tmp<-strsplit(Meta(annotations)$relation," ")[[1]]
    geoAccession2<-tmp[grep("GPL",tmp)]
    print("Looking for alternative platform")
    annotation_table<-readGEOAnn(geoAccession2) #Change to make this command only if getGEO return empty table
    symbol_column <- grep("symbol",colnames(annotation_table),ignore.case=TRUE)
    entrez_column <- grep("entrez",colnames(annotation_table),ignore.case=TRUE)
  }
  if(length(symbol_column) == 0 | length(entrez_column) == 0) return(data.frame())
  #symbols <- as.character(annotation_table[,symbol_column])
  colnames(annotation_table)[c(symbol_column,entrez_column)] = c("symbol","entrezgene_id")
  return(annotation_table)
}

