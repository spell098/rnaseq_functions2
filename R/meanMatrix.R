#' Mean matrix
#' @description For every row of the initial expression matrix, the mean for each group is calculated.
#' @return Matrix
#' \describe{
#'  \item{columns}{Groups (mean values)}
#'  \item{rows}{Genes}
#' }
#' @examples
#' expr.matrix <- readRDS("data/expr_matrix_LGVD.rds")
#' results <- readRDS("data/results_LGVD.rds")
#' topTable3 <- readRDS("data/results_LGVD.rds")
#' resultsSummary <- results_summary(results,topTable3,adjust="BH")
#' expr.matrix.signif <- signifRows(expr.matrix,resultsSummary)
#' expr.matrix.signif.means <- meanMatrix(expr.matrix.signif)
#' @seealso
#' \code{\link{signifRows}}
#' @export
meanMatrix = function(expr.matrix.signif,names2){
  names.unique <- unique(get_names(names2))
  expr.matrix.signif.means <- matrix(ncol = length(names.unique),nrow=nrow(expr.matrix.signif))
  
  for (i in 1:length(names.unique)){
    mask <- grep(names.unique[i],names2)
    x <- matrix(ncol=length(mask),nrow=nrow(expr.matrix.signif.means))
    x[,1:length(mask)] <- expr.matrix.signif[,mask]
    expr.matrix.signif.means[,i] <- apply(x,1,mean)
  }
  colnames(expr.matrix.signif.means) <- names.unique
  rownames(expr.matrix.signif.means) <- rownames(expr.matrix.signif)
  #heatmap.2(data.matrix(expr.matrix.annotated.signif[,2:12]))
  return(expr.matrix.signif.means)
}
