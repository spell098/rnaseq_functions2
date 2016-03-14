#' Results by ontology
#' Return the rows from the results topTable
#' @param selectedOntology The ontology term selected
#' @param results The results of one comparison
#' @param ontology An object
#' @examples
#' expr.matrix <- readRDS("data/expr_matrix_LGVD.rds")
#' results <- readRDS("data/results_LGVD.rds")
#' annotation1 = annotate_ensembl(rownames(expr.matrix))
#' go=annotation1[[2]]
#' typeID=annotation1[[4]]
#' result <- results[[1]]
#' ontology <- geneOntology(results,go,typeID,nrow(expr.matrix))
#' selectedOntology <- names(ontology[[1]][[3]])[1]
#' resultsByOntology <- resultsByOntology(selectedOntology,result,ontology[[1]])
#' @seealso
#' \code{\link{}}
#' @export
resultsByOntology = function(selectedOntology,result,ontology){
  if(length(selectedOntology) > 0){
    resultsOntologies = vector("list",length(selectedOntology))
    names(resultsOntologies) = selectedOntology
    for(i in 1:length(selectedOntology)){
      selected <- selectedOntology[i]
      genes <- ontology[[3]][selected]
      resultsOntologies[[i]] <- result[match(genes[[1]],as.character(result$ensembl_gene_id)),]
      resultsOntologies[[i]] <- cbind(resultsOntologies[[i]],rep(selected,nrow(resultsOntologies[[i]])))
      colnames(resultsOntologies[[i]])[ncol(resultsOntologies[[i]])] <- "Ontology"
    }
    resultsOntology = do.call(rbind, resultsOntologies)
    resultsOntology$Ontology = as.character(resultsOntology$Ontology)
    resultsOntology2 = resultsOntology
    for(j in 1:nrow(resultsOntology)){
      duplicates=grep(as.character(resultsOntology[j,]$ensembl_gene_id),as.character(resultsOntology$ensembl_gene_id))
      resultsOntology2[j,]$Ontology = paste(resultsOntology[duplicates,"Ontology"],collapse=",")
    }
    resultsOntology3 = unique(resultsOntology2)
  } else {
    resultsOntology3 <- NULL
  }
  return(resultsOntology3)
}
