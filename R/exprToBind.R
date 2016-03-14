#' Binds an expression matrix with colors corresponding to the module of the elements
#' in the rows (genes,transcripts...)
#' @import limma
#' @import ggplot2
#' @author Simon J Pelletier
#' @param bnet A list of WGCNA results
#' @param expr.matrix A matrix of expression values. Rows are genes, columns are samples
#' @param typeID Type of IDs used as rows in the expression matrix.Default="ensembl_gene_id"
#' @return
#' An expression matrix additonnal informations with two columns added as the first columns
#'  of the matrix. The columns are, in order:
#' \describe{
#'  \item{module}{Contain the module of the element (gene, transcript...)}
#'  \item{ID}{contain the IDs used afterwards to merge all of the additional
#'            annotations in subsequent operation in the pipeline}
#' }
#' @keywords cytoscape export
#' @seealso \code{\link[limma]{topTable}}
#' @seealso \code{\link[ggplot2]{ggplot}}
#' @examples
#' bnet <- readRDS("data/bnet_LGVD.rds")
#' expr.matrix <- readRDS("data/expr_matrix_LGVD.rds")
#' typeID <- "ensembl_gene_id"
#' expr.toBind <- exprToBind(bnet,expr.matrix,typeID)
#' @export
exprToBind = function(bnet,expr.matrix){
  expr.toBind = cbind(bnet$colors,rownames(expr.matrix),expr.matrix)
  colnames(expr.toBind)[1:2] = c("module","ID")
  return(expr.toBind)
}
