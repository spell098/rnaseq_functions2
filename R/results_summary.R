#' Summary of the results
#' Create a data.frame that displays every gene significant in any comparison and
#' the p-value for each comparison.
#' @author Simon J Pelletier
#' @param expr.matrix A matrix of standardized data. Columns = samples, rows = genes,transcripts,CpG...
#' @param resultsSummary A matrix of data from every gene that is significanlty different
#' @return The significant rows of the initial expression matrix
#' between at least two groups. Columns = samples, rows = genes,transcripts,CpG...
#' @examples
#' #Example with expression matrix
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
#' design <- design_contrasts(get_names(colnames(expr.matrix)))
#' lm2 <- lm2Contrast(expr.matrix,design)
#' lm2.contrast = lm2[[1]]
#' contrasts=lm2[[2]]
#' contrast.matrix=lm2[[3]]
#' names <- get_names(expr.matrix)[[1]]
#' comparisons <- get_comparisons(names)
#' pvalue = 0.05
#' logFC = c(-1.3,1.3)
#' type <- "ensembl_gene_id"
#' results_list <- results_topTable(lm2.contrast,expr.toBind,pvalue,logFC,type,genes_annotation_unique,annotations,"BH",ensemblTable)
#' results = results_list[[1]]
#' topTable3 = results_list[[2]]
#' resultsSummary <- results_summary(results,topTable3,adjust="BH")
#'
#' #Example with online dataset
#' gset <- getGEO("GSE54839", GSEMatrix =TRUE)
#' exprset <- gset[[1]]
#' type <- "ensembl_gene_id"
#' expr.matrix <- exprs(exprset)
#' names1 <- get_names(sampleNames(exprset))
#' comparisons <- comparisonsPheno(exprset)[[1]]
#' comparisonsTable <- comparisonsPheno(exprset)[[2]]
#' specie = "hsapiens"
#' attribute <- "illumina_humanht_12_v3"
#' externalSymbol <- "hgnc_symbol"
#' annotationFile <- paste0("annotations/databases/",specie,"_ensembl_",attribute,".csv")
#' #ensemblTable <- annotation_biomart(rownames(expr.matrix),specie,attribute)
#' ensemblTable <- read.csv(annotationFile)
#' genes <- ensemblTable[as.character(ensemblTable[,attribute]) != "",]
#' annotation <- annotate_ensembl(genes[,1],type)
#' annotations <- annotation[[1]]
#' go <- annotation[[2]]
#' genes_annotation_unique <- annotation[[3]]
#' pvalue = 0.05
#' expr.toBind <- NULL
#' logFC = c(-1,1)
#' selectedVariables = names(comparisons)[1]
#' namesTable <- comparisonsTable
#' names2 <- namesSelectedComparisons(as.data.frame(as.character(namesTable[,2])))
#' design = design_contrasts(names2)
#' bnet = readRDS('data/modules_characteristics_ch1.rds')
#' expr.toBind = exprToBind(bnet,voom.matrix)
#' lm2=lm2Contrast(expr.matrix,design)
#' lm2.contrast = lm2[[1]]
#' contrasts=lm2[[2]]
#' contrast.matrix=lm2[[3]]
#' results_list = results_topTable(lm2.contrast,expr.toBind,pvalue,logFC,type,genes_annotation_unique,annotations,adjust="none",ensemblTable)
#' results = results_list[[1]]
#' topTable3 = results_list[[2]]
#' resultsSummary = results_summary(results,topTable3,"no")
#' 
#' @seealso
#' \code{\link[limma]{topTable}}
#' @export
results_summary = function(results,topTable3,adjust="no"){

  resultsSummary=NULL
  genes = NULL
  if(adjust=="no") sortedBy = "P.Value" else sortedBy = "adj.P.Val"
  for(i in 1:length(results)){
    if(!is.null(results[[i]])){
      genes = unique(c(genes,as.character(results[[i]][,"symbol"])))
      genes = genes[genes != ""]
    }
  }
  if(length(genes) > 0){
    if(grep("ID",colnames(topTable3[[1]]))){
      resultsSummary=topTable3[[1]][match(genes,as.character(topTable3[[1]][,"symbol"])),][,c("ID","symbol","ensembl_gene_id","entrezgene_id","module")]
      numCol <- 5
      rownames(resultsSummary) = resultsSummary[,"ID"]
    } else {
      resultsSummary=topTable3[[1]][match(genes,as.character(topTable3[[1]][,"symbol"])),][,c("symbol","ensembl_gene_id","entrezgene_id","module")]
      numCol <- 4
      rownames(resultsSummary) = resultsSummary[,"ensembl_gene_id"]
    }
    for (i in 1:length(topTable3)){
      resultsSummary[,i+numCol] = topTable3[[i]][match(genes,topTable3[[i]]$symbol),][,sortedBy]
      colnames(resultsSummary)[i+numCol] = names(topTable3)[i]
    }
  } else {
    #resultsSummary=topTable3[[1]][match(genes,as.character(topTable3[[1]][,"symbol"])),][,c("symbol","ensembl_gene_id","entrezgene_id")]
    resultsSummary <- NULL
  }
  return(resultsSummary)
}
