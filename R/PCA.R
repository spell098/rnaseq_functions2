#' Displays a Princpial component analysis graph. USE GGPLOT!!
#' @param expr.matrix Expression set containing all the information on a dataset
#' @param colorGroups Color associated with each sample
#' @seealso
#' \code{\link[vegan]{rda}}
#' \code{\link[vegan]{ggplot}}
#' @export
PCA = function(expr.matrix,names,colorGroups,namesColor,cex.size){

  rs_rda = rda(t(expr.matrix))
  color=NULL
  for (i in 1:length(unique(names))){
    color[grep(unique(names)[i],names)] = wes()[i]
  }
  plot(rs_rda)
  text(rs_rda, display = "sites", labels = colnames(expr.matrix), choices = c(1, 2), cex= cex.size,col=color,adj=0.5,pos=3)
  legend("topright",legend=unique(namesColor),cex=1,col=unique(colorGroups),lwd=1)
  if(length(color)>length(unique(color)) & sum(color != colorGroups)>0){
    legend("topleft",legend=unique(names),cex=1,col=unique(names),lwd=1)
  }
  points(rs_rda,lwd=1,col=color,bg=colorGroups,pch=21)
  #control = grep("control",colnames(expr.matrix))
  #depression= grep("depression",colnames(expr.matrix))
  #bipolar = grep("bipolar",colnames(expr.matrix))
  #schizophernia = grep("schizophrenia",colnames(expr.matrix))
  #col=c("blue","green","orange","red")
  for(i in 1:length(unique(colnames(expr.matrix)))){
    ordiellipse(rs_rda,groups = colorGroups,show.groups = unique(colorGroups)[i],conf = 0.95,col = unique(colorGroups)[i])
  }
}
