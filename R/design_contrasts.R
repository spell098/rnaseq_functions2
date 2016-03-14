#' Create a design matrix for linear models analysis
#' @author Simon J Pelletier
#' @param names A vector of the names of all samples
#' @return An experimental design matrix
#' @keywords design linear limma
#' @seealso
#' \code{\link{get_names}}
#' @examples
#' expr.matrix <- readRDS("data/expr_matrix_LGVD.rds")
#' names <- get_names(colnames(expr.matrix))
#' design <- design_contrast(names)
#' @export
design_contrasts = function(expr.matrix){
  names <- get_names(expr.matrix)
  f <- factor(names, levels = unique(names))
  design <- model.matrix(~0+f)
  colnames(design) <- levels(f)
  rownames(design) = names
  return(design)
}
