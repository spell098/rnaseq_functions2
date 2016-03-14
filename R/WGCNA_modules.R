#' Find to which genes belong to the same modules
#' @import WGCNA
#' @author Simon J Pelletier
#' @param expr.matrix  A matrix of expression values. Rows are genes, columns are samples
#' @param ncores Number of cores available for parallel programming (foreach function)
#' @param names.unique Names of all groups of sample
#'
#' @return A list of WGCNA results
#' @keywords cytoscape export
#' @seealso \code{\link{WGCNA}}
#' @examples
#' expr.matrix <- readRDS("data/expr_matrix_LGVD.rds")
#' results <- readRDS("data/results_LGVD.rds")
#'
#' bnet = WGCNA_modules(expr.matrix,pow=3,names=names2)
#' @export
WGCNA_modules = function(expr.matrix,pow,names,deepSplit = 2,minModSize=30,saveTOMs=TRUE,minKMEtoStay=0,mergeCutHeight = 0.25, numericLabels = TRUE,pamRespectsDendro = FALSE){
  print("Looking for consensus modules")
  names.unique <- unique(names)
  
  ncores = detectCores()
  nSets = length(names.unique)
  multiExpr = vector(mode = "list", length = nSets)
  for(i in 1:length(names.unique)){
    multiExpr[[i]] = list(data = as.data.frame(t(expr.matrix[,grep(names.unique[i],names)])))
  }
  
  bnet = blockwiseConsensusModules(
    multiExpr, power = pow, minModuleSize = minModSize, deepSplit = deepSplit,
    pamRespectsDendro = pamRespectsDendro,
    mergeCutHeight = mergeCutHeight, numericLabels = numericLabels,
    minKMEtoStay = minKMEtoStay,
    saveTOMs = saveTOMs, verbose = 5,
    nThreads=ncores)
  print("Modules informations saved as bnet.rdata")
  bnet$colors <- color_groups(bnet$colors)
  return(bnet)
}
