#' Find the groups using fuzzy clustering ##FIND AUTOMATICALLY THE NUMBER OF CORES
#' @author Simon J Pelletier
#' @import Mfuzz
#' @import e1071
#' @param expr.matrix A matrix of expression values. Rows are genes, columns are samples
#' @param noCluster The number of clusters to be used. ##COULD BE DETERMINATED AUTOMATICALLY??
#' @return The function return a list of 3 objects
#' \describe{
#'  \item{clusters}{Vector in which every sample is attributed to the most likely cluster}
#'  \item{in_mfuzz}{Expression set}
#'  \item{cl_mfuzz}{Details of fuzzy cmeans}
#' }
#' @keywords fuzzy cluster clustering Hypergeometric
#' @seealso
#' \code{\link[Mfuzz]{standardise}} ,
#' \code{\link{get_names}} ,
#' \code{\link{e1071}}
#' @examples
#' expr.matrix <- readRDS("data/expr_matrix_LGVD.rds")
#' noCluster <- 3
#' fuzzy <- fuzzy_clustering(expr.matrix,noCluster)
#' fuzzy$clusters
#' fuzzy$in_mfuzz
#' fuzzy$cl_mfuzz
#'
#' @export
fuzzy_clustering = function(expr.matrix,noCluster=1){
  mask = !duplicated(rownames(expr.matrix))
  EXPR_MAT_selected= expr.matrix[mask,]
  in_mfuzz = standardise(new("ExpressionSet", exprs = t(as.matrix(EXPR_MAT_selected))))
  cl_mfuzz <- cmeans(t(EXPR_MAT_selected), noCluster, m=1.25)
  if (length(strsplit(colnames(expr.matrix)[1],"_")[[1]]) > 1){
    names = get_names(colnames(expr.matrix))
  } else {
    names = get_names(colnames(expr.matrix))
  }
  clusters=cl_mfuzz$cluster
  names(clusters) = names
  return(list(clusters,in_mfuzz,cl_mfuzz))
}
