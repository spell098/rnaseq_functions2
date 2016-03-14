#' Takes an ensembl_gene_id
#' @author Simon J Pelletier
#' @aliases annotate_ensembl
#' @param IDs All IDs (rat or human); all are ensembl_gene_id
#' @return The object resultsContrast annotated with an associated color for the selected comparison
#' @keywords ensembl IDs
#' @examples
#' expr.matrix <- readRDS("data/expr_matrix_LGVD.rds")
#' IDs <- rownames(expr.matrix)
#' annotation1 <- annotate_ensembl(IDs)
#' annotations <- annotation1[[1]]
#' go <- annotation1[[2]]
#' genes_annotation_unique <- annotation1[[3]]
#' typeID <- annotation1[[4]]
#' 
#' @seealso
#' \code{\link{colorNumericValues}}
#' @export
annotate_ensembl = function(IDs,type){
  typeID <- type
  print("Homemade annotation")
  write.csv(IDs,"annotations/gene_names.csv",quote=FALSE)
  write.table(IDs,"annotations/gene_names.ssv",quote=FALSE,sep=";")
  #write.table(rownames(expr.matrix),"annotations/gene_names.txt",sep="\t",quote=FALSE)
  if (length(grep("ENSRNOG",IDs[1])) > 0){
    system('awk -F "," \'NR==FNR{c[$2]=1;next} c[$1] == 1  && ( $3 || $4 ) {print $1","$3","$4}\' annotations/gene_names.csv annotations/databases/rat_symbols.csv > annotations/genes_annotations.csv')
    system('awk -F ";" \'NR==FNR{c[$2]=1;next} c[$1] == 1  && ( $3 || $4 ) {print $1";"$3";"$4}\' annotations/gene_names.ssv annotations/databases/allGO.ssv > annotations/go.txt')
    annotations=read.csv("annotations/genes_annotations.csv",header=FALSE)
    colnames(annotations) = c(typeID,"symbol","entrezgene_id")
    go=read.table("annotations/go.txt",header=FALSE,sep=";",quote="")
    colnames(go) = c(typeID,"go_term_name","go_domain")
  }else if (length(grep("ENSRNOT",IDs[1])) > 0){
    system('awk -F "," \'NR==FNR{c[$2]=1;next} c[$2] == 1  && ( $3 || $4 ) {print $1","$2","$3","$4}\' annotations/gene_names.csv annotations/databases/rat_symbols.csv > annotations/genes_annotations.csv')
    system('awk -F ";" \'NR==FNR{c[$2]=1;next} c[$2] == 1  && ( $3 || $4 ) {print $1";"$2";"$3";"$4}\' annotations/gene_names.ssv annotations/databases/allGO.ssv > annotations/go.txt')
    annotations=read.csv("annotations/genes_annotations.csv",header=FALSE)
    colnames(annotations) = c("ensembl_gene_id",typeID,"symbol","entrezgene_id")
    go=read.table("annotations/go.txt",header=FALSE,sep=';',quote="")
    colnames(go) = c("ensembl_gene_id",typeID,"go_term_name","go_domain")
  } else if (length(grep("ENSG",IDs[1])) > 0){
    system('awk -F "," \'NR==FNR{c[$2]=1;next} c[$1] == 1  && ( $3 || $4 ) {print $1","$3","$4}\' annotations/gene_names.csv annotations/databases/hgnc.csv > annotations/genes_annotations.csv')
    system('awk -F ";" \'NR==FNR{c[$2]=1;next} c[$1] == 1  && ( $3 || $4 ) {print $1";"$3";"$4}\' annotations/gene_names.ssv annotations/databases/go_human.ssv > annotations/go.txt')
    annotations=read.csv("annotations/genes_annotations.csv",header=FALSE)
    colnames(annotations) = c(typeID,"symbol","entrezgene_id")
    go=read.table("annotations/go.txt",header=FALSE,sep=";",quote="")
    colnames(go) = c(typeID,"go_term_name","go_domain")
  } else if (length(grep("ENST",IDs[1])) > 0){
    system('awk -F "," \'NR==FNR{c[$2]=1;next} c[$2] == 1  && ( $3 || $4 ) {print $1","$2","$3","$4}\' annotations/gene_names.csv annotations/databases/hgnc.csv > annotations/genes_annotations.csv')
    system('awk -F ";" \'NR==FNR{c[$2]=1;next} c[$2] == 1  && ( $3 || $4 ) {print $1","$2","$3","$4}\' annotations/gene_names.ssv annotations/databases/go_human.ssv > annotations/go.txt')
    annotations=read.csv("annotations/genes_annotations.csv",header=FALSE)
    colnames(annotations) =  c("ensembl_gene_id",typeID,"symbol","entrezgene_id")
    go=read.table("annotations/go.txt",header=FALSE,sep=";",quote="")
    colnames(go) = c("ensembl_gene_id",typeID,"go_term_name","go_domain")
  }
  genes_annotation_unique = annotations[!duplicated(as.character(annotations[,typeID])),]
  return(list(annotations,go,genes_annotation_unique))
}
