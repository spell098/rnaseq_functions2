#' Generates a volcanoplot with ggplot
#' @import limma
#' @import ggplot2
#' @author Simon J Pelletier
#' @param topTable3 topTable results for every gene
#' @param logfc log2 fold-change limit
#' @param pval p-value limit
#' @return A list of WGCNA results
#' @keywords cytoscape export
#' @seealso
#' \code{\link[limma]{topTable}} ,
#' \code{\link[ggplot2]{ggplot}}
#' @seealso \code{\link[limma]{plotMA}}{ : Similar function in the limma package}
#' @examples
#' topTable3 <- readRDS("data/topTable3.rds")
#' @export
volcanoGgPlot <- function(topTable3,logfc=1,pval=0.05){
  data1 = topTable3
  data2 = data1[!duplicated(data1[,1]) & !duplicated(data1$logFC),]
  signifGenes = subset(data2,P.Value < pval)
  signifGenesAdj = subset(data2,adj.P.Val < pval)
  signifGenesPFC = subset(data2,(abs(logFC) > logfc) & data2$P.Value < pval)
  signifGenesFCAdj = subset(data2,(abs(logFC) > logfc) & data2$adj.P.Val < pval)
  g = ggplot() +
    geom_point(data=data2, aes(x=logFC, y=-log2(P.Value)), alpha=1, size=1.75) +
    geom_point(data=signifGenes, aes(x=logFC, y=-log2(P.Value)), alpha=1, size=1.75,colour="grey") +
    geom_point(data=signifGenesPFC, aes(x=logFC, y=-log2(P.Value)), alpha=1, size=1.75,colour="blue") +
    geom_point(data=signifGenesAdj, aes(x=logFC, y=-log2(P.Value)), alpha=1, size=3,colour="red") +
    geom_point(data=signifGenesFCAdj, aes(x=logFC, y=-log2(P.Value)), alpha=1, size=3,colour="orange") +
    geom_text(data = signifGenesFCAdj,mapping=aes(x=logFC, y=-log2(P.Value), label=symbol),size=3, vjust=1.5, hjust=0.5) +
    theme(legend.position = "none") +
    xlab("log2 fold change") + ylab("-log2 p-value")
  g
  #+ geom_hline(aes(yintercept=-log2(pval))) + geom_vline(aes(xintercept=logFC))
}

