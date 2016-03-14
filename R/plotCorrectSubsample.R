#' Plot the subsamples correctly attributed
#' @export
plotCorrectSubsample = function(fuzzy){
  z=length(fuzzy[[4]][[4]])+1
  correct = matrix(nrow=z,ncol=20)
  for(i in 2:length(fuzzy)){
    for(j in 1:length(fuzzy[[i]][[4]])){
      if(!is.null(fuzzy[[i]][[4]][j])){
        correct[j,i] = fuzzy[[i]][[4]][j]
      } else {
        correct[j,i] = 0
      }
    }
  }
  correct[z,] = apply(correct[-(z),],2,function(x){
    #x[1]
    (15*x[1]+35*x[2])/50 # a corriger
  })
  #rownames(correct) = names(fuzzy[[i]][[4]])
  plot(x=NULL,xlim=c(1, 20), ylim=c(0,1))
  for (k in 1:nrow(correct)){
    lines(correct[k,],type="o",col=k)
  }
  layout(rbind(1,2), heights=c(7,1))
  #rect(0.5, 0.9-(nrow(correct)*0.02), 2.8, 1)
  legend(x=0,y=1,names(fuzzy[[i]][[4]]),lty=c(1,1),col=1:nrow(correct),cex=0.6,text.width=1)
}
