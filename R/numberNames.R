#' Add increasing number at the end of all strings. Duplicated names are uniquely named by incrementing
#' the number of each new duplicate found.
#' @author Simon J Pelletier
#' @param names Strings to add numbers
#' @return Vector of the numbered strings
#' @export
numberNames = function(names){
  names=as.character(names)
  names2=NULL
  for (i in unique(names)){
    l=grep(i,names)
    k=0
    for(j in l){
      k=k+1
      names2[j] = paste(names[j],k,sep="_")
    }
  }
  return(names2)
}
