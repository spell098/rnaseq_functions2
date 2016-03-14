#' Biomart annotation
#' @author Simon J Pelletier
#' @aliases annotation_biomart
#' @examples
#' specieEnsembl <- "hsapiens_gene_ensembl"
#' attribute <- "affy_hg_u95av2"
#' gset <- getGEO("GSE12654")
#' exprset <- gset[[1]]
#' type <- "ensembl_gene_id"
#' expr.matrix <- exprs(exprset)
#' IDs <- rownames(expr.matrix)[1:10]
#' x<-annotation_biomart(IDs,specieEnsembl,attribute)
#' @export
annotation_biomart <- function(IDs,specie,attribute,mart){
  if(mart == "ENSEMBL_MART_ENSEMBL"){
    specieEnsembl <- paste0(specie,"_gene_ensembl") # attribute <- "affy_hg_u95av2"
    attributes1 <- c("ensembl_gene_id",attribute)
    filter <- attribute
  } else if(mart == "ENSEMBL_MART_SNP"){
    specieEnsembl <- paste0(specie,"_snp") # attribute <- "refsnp_id"
    attributes1 <- c("ensembl_gene_stable_id",attribute,"chr_name","chrom_start","chrom_end")
    filter <- "snp_filter"
  }
  ensembl=useMart(mart,specieEnsembl)
  ID.annotated = getBM(attributes=attributes1,
                       filters = filter,
                       values = IDs,
                       mart = ensembl)
  write.csv(ID.annotated,paste0("annotations/databases/",specie,"_ensembl_",attribute,".csv"),quote=FALSE,row.names=FALSE)
}
