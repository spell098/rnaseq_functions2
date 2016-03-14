#' Annotation of a matrix with colors from standardized values
#' @author Simon J Pelletier
#' @param selectedComparison The selected comparison to display
#' @param results A matrix of all the values that correspond to the presence of RNA for each transcript analyzed
#' @return The object resultsContrast annotated with an associated color for the selected comparison
#' @keywords colors
#' @examples
#' results <- expr.matrix <- readRDS("data/results_LGVD.rds")
#' selectedComparison <-
#' resultsContrastRamp <- addColorsMatrix(resultsContrast,selectedComparison)
#' @export
addColorsMatrixResults <- function(results,selectedComparison){
  resultsContrastRamp = data.matrix(results[[selectedComparison]]$logFC)
  colnames(resultsContrastRamp) = selectedComparison
  rownames(resultsContrastRamp) = as.character(results[[selectedComparison]]$symbol)
  resultsContrastRamp[,selectedComparison] = colorNumericValues(resultsContrastRamp)
  return(resultsContrastRamp)
}
