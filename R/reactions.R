#' Top reactions
#' @import graphite
#' @description Finds top reactions of each results table
#' @param results Results list (topTables)
#' @param entrezgene_ids All ids, entrez format
#' @examples
#' expr.matrix <- readRDS("data/expr_matrix_LGVD.rds")
#' annotation1 = annotate_ensembl(rownames(expr.matrix))
#' annotations=annotation1[[1]]
#' results <- readRDS("data/results_LGVD.rds")
#' entrezgene_ids <- unique(annotations$entrezgene_id)
#' reactome <- pathways("rnorvegicus", "reactome")
#' prepareSPIA(reactome, "reactome",print.names=TRUE)
#' res <- reactions(results,entrezgene_ids)
#' @seealso
#' \code{\link[graphite]{pathways}}
#' \code{\link[graphite]{prepareSPIA}}
#'
#' @export
reactions = function(results,entrezgene_ids){
  ALL <- as.character(entrezgene_ids[!is.na(entrezgene_ids)])
  print(head(ALL))
  res=vector("list",length(results))
  names(res) = names(results)
  print("Looking for top reactions...")
  for (i in 1:length(results)){
    print(paste("running results list",i,"of",length(results),sep=" "))
    if (nrow(results[[i]]) > 0){
      DE = as.numeric(results[[i]][!is.na(results[[i]][,"entrezgene_id"]),]$logFC)
      names(DE) <- as.character(as.vector(results[[i]][!is.na(results[[i]]$entrezgene),]$entrezgene))
      print(DE)
      x = match(ALL,names(DE))
      DE = DE[x[!is.na(x)]]
      print(DE)
      if(!is.null(DE)){
        res[[i]] <- runSPIA(de=DE, all=ALL, "reactome")
        print(paste0("DONE #",i))
      } else {
        print("NULL")
      }
    } else {
      res[[i]] = NULL
      print(paste0("DONE #",i,"\nNothing Found..."))
    }
  }
  print("All top reactions finished")
  return(res)
}
