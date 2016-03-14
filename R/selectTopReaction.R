#' Select a reaction from a top reactions table
#' @author Simon J Pelletier
#' @param selectedComparison Selected Comparison
#' @param topReactions Top reactions
#' @return The selected top reactions informations
#' @examples
#' selectedComparison <- "a2-a1"
#' topReactions <- readRDS("data/reactions_mfuzz.rds")
#' selectedTopReactions <- selectTopReaction(selectedComparison,topReactions)
#' head(selectedTopReactions)
#' @export
selectTopReaction = function(selectedComparison,topReactions){
  x = strsplit(selectedComparison,"-")[[1]]
  alt1=paste(x[1],x[2],sep="-")
  alt2=paste(x[2],x[1],sep="-")
  if (!is.null(topReactions[[alt1]])){
    selectedTopReactions = topReactions[[alt1]]
  } else {
    selectedTopReactions = topReactions[[alt2]]
  }
  return(selectedTopReactions)
}
