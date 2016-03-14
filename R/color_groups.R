#' Add color to each sample. Samples with the same color belong to the same group
#' @author Simon J Pelletier
#' @param names Names of each sample
#' @param names.unique Name of each group
#' @return A vector of strings corresponding to the color of all samples.
#' @export
color_groups = function(names){
  if(is.null(names)) return(NULL)
  names.unique <- unique(names)
  colors = vector("character",length(names))
  for (i in 1:length(names.unique)){
    if(length(names.unique) < length(wes())){
      colors[(grep(names.unique[i],names))] = wes()[i]
    } else if(length(names.unique) < length(colours(distinct=TRUE))){
      colors[(grep(names.unique[i],names))] = colours(distinct=TRUE)[i]
    } else if(length(names.unique) < length(colours(distinct=FALSE))){
      colors[(grep(names.unique[i],names))] = colours(distinct=FALSE)[i]
    } else { 
      colors[(grep(names.unique[i],names))] = rainbow(length(names.unique))[i]
    }
  }
  return(colors)
}
