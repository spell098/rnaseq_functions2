#' Possible Comparisons
#' @description Find all the possible pairs of groups to compare
#' @author Simon J Pelletier
#' @param names A vector of the names of all samples
#' @return A vector of all possible comparisons
#' @keywords comparisons
#' @export
possibleComparisons = function(names){
  possibleComparisons = NULL
  n=1
  for (i in 1:(length(unique(names))-1)){
    for (j in (i+1):length(unique(names))){
      possibleComparisons[n] = paste0(unique(names)[i],"-",unique(names)[j])
      n=n+1
    }
  }
  return(possibleComparisons)
}
