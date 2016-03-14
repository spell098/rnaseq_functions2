#' Take all the name and export all possible comparisons
#' @author Simon J Pelletier
#' @param names Names of all groups compared
#' @return A vector of all the possible comparisons in the format names1//name2 (VERIFY!!)
#' @keywords comparisons
#' @export
get_comparisons = function(names){
  comparisons=NULL
  conditions = matrix(nrow=length(names),ncol = length(strsplit(names[1],"_")[[1]]))
  for (i in 1:length(names)){
    conditions[i,] = strsplit(names[i],"_")[[1]]
  }
  for (j in 1:ncol(conditions)){
    tmp=unique(conditions[,j])
    comparisons[j] = paste(tmp,collapse="\\")
  }
  return(comparisons)
}
