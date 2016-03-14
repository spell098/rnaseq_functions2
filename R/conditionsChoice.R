#' Transform intensity values in a gradient of color from green to red.
#' @author Simon J Pelletier
#' @import GEOquery
#' @param selectedVariables Selected comparison of 2 groups
#' @return A vector of colors. The darkest green values are the lowest and the darkest red values are the highest.
#' @keywords comparisons
#' @seealso
#' \code{\link[GEOquery]{getGEO}}
#' \code{\link[Biobase]{ExpressionSet}}
#' @examples
#' expr.matrix <- readRDS("data/expr_matrix_LGVD.rds")
#' selectedVariables <- c("HIGH\\LOW","DORSAL\\VENTRAL")
#' names <- get_names(expr.matrix)
#' conditionsChoice(selectedVariables,names)
#' @export
conditionsChoice = function(selectedVariables,names){
  names2=NULL
  if (length(selectedVariables) > 0){
    for (m in length(selectedVariables):1){
      conditionsSelected = strsplit(selectedVariables[m],'\\\\')
      for (j in conditionsSelected[[1]]){
        range=grep(j,names)
        if(m == length(selectedVariables)){
          names2[range] = j
        } else {
          names2[range] = paste(j,names2[range],sep="_")
        }

      }
    }
  }
  return(names2)
}
