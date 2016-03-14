#' Returns the summary for a module
#' @author Simon J Pelletier
#' @param results The results of a comparison between two groups.
#' @param expr.toBind Expression data.matrix with supplemental colums: module and ID
#' @examples
#' results <- readRDS("data/results_LGVD.rds")
#' expr.matrix <- readRDS("data/expr_matrix_LGVD.rds")
#' bnet <- readRDS("data/bnet_LGVD.rds")
#' 
#' expr.toBind <- cbind(module=bnet$colors,ensembl_gene_id=rownames(expr.matrix),expr.matrix)
#' result <- results[[1]]
#' modulesTable <- modules_summary(results[[1]],expr.toBind)
#' 
#' expr.matrix <- readRDS("data/expr_matrix_LGVD.rds")
#' type = "ensembl_gene_id"
#' design <- design_contrasts(get_names(colnames(expr.matrix)))
#' lm2 <- lm2Contrast(expr.matrix,design)
#' lm2.contrast = lm2[[1]]
#' contrasts=lm2[[2]]
#' contrast.matrix=lm2[[3]]
#' names <- get_names(expr.matrix)[[1]]
#' specie = "rnorvegicus"
#' attribute <- "ensembl_gene_id"
#' externalSymbol <- "rgd_symbol"
#' annotationFile <- paste0("annotations/databases/",specie,"_ensembl_",attribute,".csv")
#' #ensemblTable <- annotation_biomart(rownames(expr.matrix),specie,attribute)
#' ensemblTable <- read.csv(annotationFile)
#' genes <- ensemblTable[as.character(ensemblTable[,attribute]) != "",]
#' annotation <- annotate_ensembl(genes[,1],type)
#' annotations <- annotation[[1]]
#' go <- annotation[[2]]
#' genes_annotation_unique <- annotation[[3]]

#' 
#' maskModule <- as.numeric(as.character(results[[1]]$module)) == 
#' resultsModule <- list(resultsModule=results[[1]][maskModule,])
#' ontologyModule <- geneOntology(resultsModule,go,type,length(rownames(voom.matrix)))
#' topModuleGO <- cbind(rownames(ontologyModule[[1]][[2]]),ontologyModule[[1]][[2]])
#' resultsModuleOntology <- resultsByOntology(1,resultsModule[[1]],ontologyModule[[1]])
#' @keywords geo annotation
#' @seealso
#' \code{\link[stats]{phyper}}
#' @export
modules_summary = function(resultTable,expr.toBind){
  ngenes = nrow(expr.toBind)
  saveRDS(resultTable,"resultTableTest.rds")
  saveRDS(expr.toBind,"exprtobind2.rds")
  modules.unique = unique(as.character(resultTable[,"module"]))
  if(length(modules.unique) > 1){
    print(paste0("Number of modules:",length(modules.unique)))
    modules.all = as.character(expr.toBind[,"module"])
    modules = as.character(resultTable[,"module"])
    modules_table = data.frame("Module"=1,"Hits"=1,"total"=1,"pvalue"=1,"qvalue"=1,"ratio"=1)
    totalModules = length(unique(expr.toBind[,"module"]))
    for (i in 1:length(modules.unique)){
      n = sum(modules == modules.unique[i])
      N = sum(modules.all == modules.unique[i])
      modules_table[i,c(1,2,3,6)] = c(modules.unique[i],n,N,n/N)
    }
    modules_table = modules_table[order(modules_table[,"Hits"],decreasing = T),]
    rownames(modules_table) = modules_table[,"Module"]
    for (c in 1:nrow(modules_table)){
      # TO CORRECT P-VALUE
      # # of test = # of annotations tested
      # p-value of GO is calculated using an hypergeometric distribution
      hits = as.numeric(modules_table[,"Hits"][c])
      possibleHits = as.numeric(modules_table[,"total"][c])
      M = totalModules
      N = ngenes-possibleHits #total number of genes
      pvalue = 1-phyper(hits-1,possibleHits,N,nrow(resultTable))
      modules_table[,"pvalue"][c] = pvalue
      
    }
    modules_table = modules_table[order(modules_table[,"pvalue"]),]
    modules_table[,"qvalue"] = round(p.adjust(modules_table[,"pvalue"],method="BH"),digits=4)
    for (i in 1:length(modules_table[,"qvalue"])){
      if (modules_table[,"qvalue"][i] < 0.001){
        modules_table[,"qvalue"][i] = "< 0.001"
      }
    }
    return(modules_table)
  } else {
    print("No module found?")
    return(NULL)
  }
}
