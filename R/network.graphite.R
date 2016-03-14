#' Create a graphNEL graph with directed edges from the results for a specified reaction
#' @author Simon J Pelletier
#' @import graphite
#' @param pSymbol Symbol IDs
#' @param expr.matrix Expression set containing all the information on a dataset
#' @param typeID Type of IDs used as rows in the expression matrix.Default="ensembl_gene_id"
#' @param toptableNodes topTable with added column of the reaction nodes it belongs to
#' @param topTable3 topTable results for every gene
#' @param pvalue p-value limit
#' @param selectedReaction Selected reaction from the topReactions
#' @param selectedComparison Comparison studied. Name of the matrices (topTables) in the list "results"
#' @param expr.toBind Expression data.matrix with supplemental colums: module and ID
#' @param threshold Threshold of coexpression
#' @return A graphNEL graph with directed edges
#' @keywords graphite network
#' @seealso
#' \code{\link{export2cytoscape}}
#' \code{\link{graphite}}
#' @examples
#' expr.matrix <- readRDS("data/expr_matrix_LGVD.rds")
#' topReactions <- readRDS("data/reactions_mfuzz.rds")
#' topTable3 <- readRDS("data/topTable3_LGVD.rds")
#' selectedReaction <- "RHO GTPases activate PAKs"
#' ratReactome <- pathways("rnorvegicus", "reactome")
#' p <- ratReactome[[selectedReaction]]
#' pSymbol <- convertIdentifiers(p, "SYMBOL")
#' toptableNodes = cbind(
#'    topTable3[[1]][match(nodes(pSymbol),topTable3[[1]]$symbol),],
#'    Reaction = rep(selectedReaction[[1]],length(nodes(pSymbol)))
#' )
#' toptableNodes1 = toptableNodes[!is.na(toptableNodes[,1]),]
#' typeID = "ensembl_gene_id"
#' threshold = 0.8
#'
#' network.graphite(pSymbol,expr.matrix,toptableNodes,typeID,threshold,selectedReaction,topTable3,1,expr.toBind,pvalue)
#' @export
network.graphite = function(pSymbol,expr.matrix,toptableNodes,typeID,threshold,selectedReaction,toptable3,selectedComparison,expr.toBind,pvalue){
  nodes = cbind(symbol=nodes(pSymbol),Reaction=rep(selectedReaction,length(nodes(pSymbol))))
  colnames(nodes) = c("symbol","reaction")

  edgeColor = NULL
  arrow=NULL
  tail=NULL
  toto=NULL
  networkSelectedReaction = export2cytoscape(expr.matrix,toptableNodes,threshold,typeID)
  networkSelectedReaction$edgeData[,"fromNode"] = as.character(toptableNodes[,"symbol"])[match(networkSelectedReaction$edgeData[,"fromNode"],toptableNodes[,"ensembl_gene_id"])]
  networkSelectedReaction$edgeData[,"toNode"] = as.character(toptableNodes[,"symbol"])[match(networkSelectedReaction$edgeData[,"toNode"],toptableNodes[,"ensembl_gene_id"])]

  edges = merge(networkSelectedReaction$edgeData,edges(pSymbol),by.x=c("fromNode","toNode","type","direction"),by.y=c("src","dest","type","direction"),all.x=TRUE,all.y=TRUE)

  gSymbol <- pathwayGraph(pSymbol,edge.type=NULL)
  g2 <- addEdge(networkSelectedReaction$edgeData[,"fromNode"], networkSelectedReaction$edgeData[,"toNode"], gSymbol, networkSelectedReaction$edgeData[,"weight"])

  g3 <- layoutGraph(g2)

  toptableNodes = merge(toptable3[[selectedComparison]],nodes,by="symbol",all.y=TRUE)
  toptableNodes$logFC[is.na(toptableNodes$logFC)] = 0
  fc=toptableNodes$logFC
  toptableNodes$adj.P.Val[is.na(toptableNodes$adj.P.Val)] = 1
  #modules=rep(1,nrow(toptableNodes))
  if (!is.null(expr.toBind)) {
    toptableNodes$module[is.na(toptableNodes$module)] = 0
    rainbow=rainbow(max(as.numeric(as.character(toptableNodes$module)))+1)
    modules=apply(toptableNodes,1,function(x){
      rainbow[(as.numeric(as.character(x["module"])))+1]
    })
  } else modules = rep("black",nrow(toptableNodes))
  textCol = unlist(modules)
  colorIntensity=apply(toptableNodes,1,function(x){
    if (as.numeric(x["logFC"]) > 0) {col=rgb(red=1,green=0,blue=0,as.numeric(x["logFC"])/(max(toptable3[[selectedComparison]][,"logFC"])+0.000001))
    }else col=rgb(red=0,green=0,blue=1,(abs(as.numeric(x["logFC"])/min(toptable3[[selectedComparison]][,"logFC"])+0.000001)))
  })
  signif=apply(toptableNodes,1,function(x){
    if (as.numeric(x["adj.P.Val"]) < pvalue){ signif="triangle"
    }else signif="circle"
  })

  types.edges = as.character(unique(edges[,"type"]))
  #types.colors = wes[1:length(types.edges)]
  types.colors = c("blue","red","black","cyan","green","dodgerblue","forestgreen","deeppink")
  names(types.colors) = types.edges
  numChildren = networkSelectedReaction$nodeData$numChildren
  names(numChildren) = networkSelectedReaction$nodeData[,4]



  for (i in 1:nrow(edges)){
    a=rownames(edges)[!is.na(match(edges[,2],edges[,1][i]))]
    b=rownames(edges)[!is.na(match(edges[,1],edges[,2][i]))]
    c=match(a,b)
    d=c[!is.na(c)]
    d=a[!is.na(a[c])]

    if (length(d) > 0 && as.character(edges[,"type"][i]) != "binding"){
      if(length(grep("ACTIVATION",as.character(edges[i,"type"]))) > 0 && length(grep("ACTIVATION",as.character(edges[d,"type"]))) > 0){
        arrow[i] = "normal"
        tail[i]="normal"
        edgeColor[i] = "green"
      } else if (length(grep("INHIBITION",as.character(edges[i,"type"]))) > 0 && length(grep("INHIBITION",as.character(edges[d,"type"]))) > 0){
        arrow[i] = "tee"
        tail[i]="tee"
        edgeColor[i] = "green"
      } else if (length(grep("INHIBITION",as.character(edges[i,"type"]))) > 0 && length(grep("ACTIVATION",as.character(edges[d,"type"]))) > 0){
        arrow[i] = "tee"
        tail[i]="normal"
        edgeColor[i] = "green"
      } else if (length(grep("ACTIVATION",as.character(edges[i,"type"]))) > 0 && length(grep("INHIBITION",as.character(edges[d,"type"]))) > 0){
        arrow[i] = "tee"
        tail[i]="normal"
        edgeColor[i] = "green"
      } else {
        tail[i]="normal"
        arrow[i] = "normal"
        edgeColor[i] = "pink"
      }

    }else if(length(grep("ACTIVATION", as.character(edges$type[i]))) == 1){
      arrow[i] = "normal"
      tail[i] = "none"
      edgeColor[i] = "black"
    } else if (length(grep("INHIBITION", as.character(edges$type[i]))) == 1){
      arrow[i] = "none"
      tail[i] = "tee"
      edgeColor[i] = "cyan"
    } else if (edges$type[i] == "binding"){
      arrow[i] = "none"
      tail[i] = "none"
      edgeColor[i] = "red"
    } else if (edges$type[i] == "adjacency(WGCNA)"){
      arrow[i] = "none"
      tail[i] = "none"
      edgeColor[i] = "blue"
    } else {
      arrow[i] = "none"
      tail[i] = "none"
      edgeColor[i] = "yellow"
    }
  }


  edgeWeight=apply(edges,1,function(x){
    if (!is.na(x["weight"])) as.numeric(x["weight"])*3 else 1
  })
  edgeNames = paste0(edges[,"fromNode"],"~",edges[,"toNode"])
  names(edgeWeight) = edgeNames
  names(edgeColor) = edgeNames
  names(tail) = edgeNames
  names(arrow) = edgeNames


  as.numeric(toptable3[[selectedComparison]][,"logFC"])
  names(colorIntensity) = as.character(toptableNodes[,1])
  names(signif) = as.character(toptableNodes[,1])
  names(textCol) = as.character(toptableNodes[,1])

  nodeRenderInfo(g3) <- list(
    fill     = colorIntensity,
    col      = "black",
    lty      = "solid",
    lwd      = (numChildren)+1,
    shape    = signif,
    fontsize = 15,
    textCol=textCol,
    cex=0.9)


  edgeRenderInfo(g3) <-  list(
    col = edgeColor,
    lwd = edgeWeight,
    arrowhead=arrow,
    arrowtail=tail)
  graph.par(list(graph=list(
    main     = selectedComparison,
    sub      = selectedReaction,
    cex.main = 1.8,
    cex.sub  = 1.4,
    col.sub  = "gray")

  )
  )
  return(g3)
}
