#' Export results in cytoscape format
#' @author Simon J Pelletier
#' @import WGCNA
#' @param expr.matrix A matrix of numeric values. Rows are genes, columns are samples
#' @param results A topTable with only the genes considered significant (according to
#' the p-value and logFC limits provided in the function results_topTable)
#' @param threshold The lowest value of correlation to make the link
#' considered significantly high
#' @param type The type of ID used, e.g.
#' @return The results in cytoscape format
#' @keywords cytoscape export
#' @seealso \code{\link[WGCNA]{adjacency}}
#' @examples
#' expr.matrix <- readRDS("data/expr_matrix_LGVD.rds")
#' results <- readRDS("data/results_LGVD.rds")
#' export2cytoscape <- export2cytoscape(expr.matrix,results[[1]])
#' @export
export2cytoscape = function(expr.matrix,results,threshold=0.9,type="ensembl_gene_id"){
  expr.matrix.selected = expr.matrix[match(as.character(results[,type][!is.na(results[,type])]),rownames(expr.matrix)),]
  expr.adj.selected = adjacency(matrix(as.numeric(t(expr.matrix.selected)),ncol=nrow(expr.matrix.selected),nrow=ncol(expr.matrix.selected)))
  rownames(expr.adj.selected) = rownames(expr.matrix.selected)
  colnames(expr.adj.selected) = rownames(expr.matrix.selected)
  network = exportNetworkToCytoscape(expr.adj.selected, threshold = threshold)
  network$nodeData = merge(network$nodeData,results, by.x = "nodeName", by.y=type,all.y=TRUE)
  network$edgeData = cbind(network$edgeData[,1:4],rep("adjacency(WGCNA)",nrow(network$edgeData)))
  colnames(network$edgeData)[ncol(network$edgeData)] = "type"

  numChildren = rep(0,nrow(network$nodeData))
  #pourrait peut-etre etre changé pour utiliser with()
  for(n in 1:nrow(network$nodeData)){
    if (!is.na(as.character(network$nodeData[,1][n]))){
      numChildren[n] = length(grep(as.character(network$nodeData[,1][n]),c(as.character(network$edgeData[,1]),as.character(network$edgeData[,2]))))
    } else {
      numChildren[n] = 0
    }
  }
  network$nodeData = cbind(symbol=network$nodeData,numChildren)
  #filter="symbol"
  #write.csv(network$nodeData[,-12],file = paste0("node_rnaseq_",name,".csv"))
  # write All result GO
  #htmlReport(go[[i]] , file= paste0("ontology_", colnames(contrast.matrix)[i]), ".html")
  return(network)
}
