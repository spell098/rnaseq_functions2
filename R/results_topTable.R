#' Calculate results
#' @author Simon J Pelletier
#' @import limma
#' @param lm2.contrast MArrayLM object
#' @param expr.toBind Expression data.matrix with supplemental colums: module and ID
#' @param pvalue p-value limit (alpha)
#' @param logFC log2 fold change limit #VERIFY IT IS REALLY LOG2 (default = 1)
#' @param genes_annotation_unique All the annotations for every ID in the 2nd row in expr.toBind
#' (or rownames of expr.matrix)
#' @param adjust Which correction for multiple analysis to use (default = "no").
#' Note: None is different than no somehow
#' @param annotations Annotation (???)
#' @return
#' \describe{
#'  \item{results}{topTable of only the significant results}
#'  \item{topTable3}{Complete topTable}
#' }
#' @keywords limma linear
#' @seealso
#' \code{\link[limma]{topTable}} ,
#' \code{\link[limma]{MArrayLM}} ,
#' \code{\link{transcript_count}} ,
#' \code{\link{expr.toBind}}
#' @examples
#' expr.matrix <- readRDS("data/expr_matrix_LGVD.rds")
#'
#' type <- "ensembl_gene_id"
#' design <- design_contrasts(get_names(colnames(expr.matrix)))
#' lm2 <- lm2Contrast(expr.matrix,design)
#' lm2.contrast = lm2[[1]]
#' contrasts=lm2[[2]]
#' contrast.matrix=lm2[[3]]
#' names <- get_names(expr.matrix)
#' specie_genome = "hsapiens_gene_ensembl"
#' attribute <- "ensembl_gene_id"
#' externalSymbol <- "hgnc_symbol"
#' ensemblTable <- convert_id(specie_genome,attribute)
#' genes <- ensemblTable[as.character(ensemblTable[,attribute]) != "",]
#' annotation <- annotate_ensembl(unique(as.character(ensemblTable[,type])),type)
#' annotations <- annotation[[1]]
#' go <- annotation[[2]]
#' genes_annotation_unique <- annotation[[3]]
#' pvalue = 0.05
#' logFC = c(-0.9,0.9)
#' type <- "ensembl_gene_id"
#' results_list <- results_topTable(lm2.contrast,expr.toBind,pvalue,logFC,type,genes_annotation_unique,annotations,"no")
#' results = results_list[[1]]
#' topTable3 = results_list[[2]]


#' #Example with online dataset fuzzy
#' gset <- getGEO("GSE25860", GSEMatrix =TRUE) #GSE25860 GSE61276 affy_hg_u133_plus_2
#'
#' exprset <- gset[[1]]
#' type <- "ensembl_gene_id"
#' expr.matrix <- exprs(exprset)
#' names1 <- get_names(sampleNames(exprset))
#' comparisons <- comparisonsPheno(exprset)[[1]]
#' comparisonsTable <- comparisonsPheno(exprset)[[2]]
#' specie_genome = "hsapiens_gene_ensembl"
#' specie <- "hsapiens"
#' attribute <- "illumina_humanht_12_v4"
#' externalSymbol <- "hgnc_symbol"
#' #ensemblTable <- annotation_biomart(rownames(expr.matrix),specie,attribute)
#' annotationFile <- paste0("annotations/databases/",specie,"_ensembl_",attribute,".csv")
#' ensemblTable <- read.csv(annotationFile)
#' genes <- ensemblTable[as.character(ensemblTable[,attribute]) != "",]
#' annotation <- annotate_ensembl(genes[,1],type)
#' annotations <- annotation[[1]]
#' go <- annotation[[2]]
#' genes_annotation_unique <- annotation[[3]]
#' pvalue = 0.05
#' logFC = c(-0.9,0.9)
#' annotationFile <- paste0("annotations/databases/",specie,"_ensembl_",attribute,".csv")
#' ensemblTable <- read.csv(annotationFile)
#'
#' #clusters <- fuzzy_clustering(expr.matrix,3)[[1]]
#' #color = unlist(lapply(clusters,function(x){wes()[x]}))
#' #namesColor = paste0("a",clusters)
#' comparisons = get_comparisons(namesColor)
#' clusterNames = numberNames(namesColor)
#' selectedVariables = comparisons[1]
#' namesTable <- comparisonsTable
#' names2 <- namesSelectedComparisons(namesTable)
#' design = design_contrasts(names2)
#' expr.toBind <- NULL

#' names2 = conditionsChoice(selectedVariables,namesColor)
#' design = design_contrasts(names2)
#' expr.toBind <- cbind(expr.matrix,module = rep("white",nrow(expr.matrix)))
#' lm2=lm2Contrast(expr.matrix,design)
#' lm2.contrast = lm2[[1]]
#' contrasts=lm2[[2]]
#' contrast.matrix=lm2[[3]]
#' results_list = results_topTable(lm2.contrast,expr.toBind,pvalue,logFC,type,genes_annotation_unique,annotations,"no",ensemblTable)
#' results = results_list[[1]]
#' topTable3 = results_list[[2]]
#' @export
results_topTable = function(lm2.contrast,expr.toBind,pvalue=0.05,logFC=1,type,genes_annotation_unique,annotations,adjust="no",ensemblTable){
  #x<-annotation[!duplicated(annotation$ID),]
  print("Running results...")
  results = vector("list", ncol(lm2.contrast$coefficients))
  topTable2 = vector("list", ncol(lm2.contrast$coefficients))
  topTable3 = vector("list", ncol(lm2.contrast$coefficients))
  names(results) = colnames(lm2.contrast$coefficients)
  names(topTable2) = colnames(lm2.contrast$coefficients)
  names(topTable3) = colnames(lm2.contrast$coefficients)
  #saveRDS(genes_annotation_unique, "tests/genes_annotation_unique.csv")
  for (i in 1:ncol(lm2.contrast$coefficients)){
    print(paste0("running ","1/",i))
    topTable1 = topTable(lm2.contrast , coef=i ,  sort.by = "none",
                        n=nrow(lm2.contrast$coefficients),adjust.method=adjust) #ENSG00000124831
    if(length(grep("ID",colnames(topTable1))) < 1){
      print("ID column found")
      if (!is.null(expr.toBind)) {
        print("expr.tobind found")
        tmp <- cbind(module=expr.toBind[,"module"],ID=rownames(topTable1),topTable1)
      } else {
        tmp <- cbind(module=rep("white",nrow(topTable1)),ID=rownames(topTable1),topTable1)
        print("expr.tobind not found")
      }
    } else {
      print("ID column found")
      tmp <- cbind(module=rep("white",nrow(topTable1)),topTable1)
      #print(head(tmp))
    }
    colnames(ensemblTable)[2] <- "ID"
    ensemblTable2 <- ensemblTable[ensemblTable[,2] != "",]
    ensemblTable2 <- ensemblTable2[!duplicated(ensemblTable2[,2]),]
    #print(head(ensemblTable2))
    topTable2 <- merge(tmp,ensemblTable2,by="ID")
    topTable3[[i]] <- merge(annotations,topTable2,by=type)
    topTable3[[i]] <- topTable3[[i]][!duplicated(as.character(topTable3[[i]][,type])),]
    #print(head(topTable2[[i]]))
    #topTable3[[i]] = merge(genes_annotation_unique,topTable2[[i]],by=type)
    #go_results[[i]] = merge(genes_ontology,)
    if(adjust=="no") sortedBy = "P.Value" else sortedBy = "adj.P.Val"

    mask = (topTable3[[i]]$logFC > logFC[2] | topTable3[[i]]$logFC < logFC[1]) & topTable3[[i]][,sortedBy] < pvalue
    results[[i]] = topTable3[[i]][mask,]

    ratio_transcript_signif = vector("numeric",nrow(results[[i]]))
    transcript_total = vector("numeric",nrow(results[[i]]))

    #pas efficace, mais ca marche  # Il y a des replicats de transcripts, pas bon... à fixer
    if(nrow(results[[i]]) >= 1){
      print("counting transcripts...")
      count_transcript = transcript_count(results[[i]],annotations,type)
      results[[i]] = cbind(results[[i]],count_transcript)
      print("counting transcripts done.")
    }
    x <- data.frame(transcript_signif = rep(0,nrow(topTable3[[i]])),ratio_transcript_signif = rep(0,nrow(topTable3[[i]])))
    topTable3[[i]] = cbind(topTable3[[i]],x)

    topTable3[[i]][(match(results[[i]][,type],topTable3[[i]][,type])),] = results[[i]]
    rownames(topTable3[[i]]) <- topTable3[[i]]$ID
  }
  #saveRDS(topTable2,"tests/topTable2.rds")
  #saveRDS(topTable1,"tests/topTable1.rds")

  return(list(results,topTable3))
}
