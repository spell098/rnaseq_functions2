#' Annotation of a matrix with colors for contrast values
#' @author Simon J Pelletier
#' @param selectedComparison The selected comparison to display
#' @param resultsContrast A matrix of all the values that correspond to the relative presence of RNA for each transcript analyzed
#' @return The object resultsContrast annotated with an associated color for the selected comparison
#' @keywords colors
#' @seealso \code{\link[rnaseqApp]{colorNumericValues}}
#' @examples
#' resultsContrast <-
#' selectedComparison <-
#' resultsContrastRamp <- addColorsMatrix(resultsContrast,selectedComparison)
#' @export
addColorsMatrix <- function(resultsContrast,selectedComparison){
  resultsContrastRamp = data.matrix(resultsContrast[,selectedComparison])
  colnames(resultsContrastRamp) = selectedComparison
  rownames(resultsContrastRamp) = rownames(resultsContrast)
  resultsContrastRamp[,selectedComparison] = colorNumericValues(resultsContrast[,comparison]) #Transform values in a gradient of colors
  return(resultsContrastRamp)
}
