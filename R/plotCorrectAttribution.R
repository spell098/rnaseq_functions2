#' Plot the groups correctly attributed
#' @export
plotCorrectAttribution = function(fuzzy){
  correct = NULL
  for(i in 2:length(fuzzy)){
    correct[i] = fuzzy[[i]][[3]]
  }
  plot(correct,type="o",ylim = c(0,1),axes=FALSE)
  axis(side = 1,at=2*(0:10))
  axis(side = 2,at=0.1*(0:10))
  #lines(correct)
  #text(1:length(fuzz), correct, 1:length(fuzz), cex=0.6, pos=3, col="red")
}
