#' Select in an expression matrix only the rows with significant elements (genes,transcripts...)
#' @author Simon J Pelletier
#' @param expr.matrix A matrix of standardized data. Columns = samples, rows = genes,transcripts,CpG...
#' @param resultsSummary A matrix of data from every gene that is significanlty different
#' @return The significant rows of the initial expression matrix
#' between at least two groups. Columns = samples, rows = genes,transcripts,CpG...
#' @examples
#' expr.matrix <- readRDS("data/expr_matrix_LGVD.rds")
#' results <- readRDS("data/results_LGVD.rds")
#' topTable3 <- readRDS("data/topTable3_LGVD.rds")
#' resultsSummary <- results_summary(results,topTable3)
#' @seealso
#' \code{\link[limma]{topTable}}
#' \code{\link{results_summary}}
#' @export
signifRows = function(expr.matrix,selectedElements){
  expr.matrix.signif = matrix(ncol = ncol(expr.matrix),nrow=length(selectedElements))
  expr.matrix.signif[1:length(selectedElements),] = as.matrix(expr.matrix[match(selectedElements,rownames(expr.matrix)),])
  colnames(expr.matrix.signif) = colnames(expr.matrix)
  rownames(expr.matrix.signif) = as.character(selectedElements)
  return(expr.matrix.signif)
}
