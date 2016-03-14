#' Convert any type of ID to ensembl ids
#' @author Simon J Pelletier
#' @param specie_genome The ensembl ID of the genome studied e.g. rnorvegicus_gene_ensembl
#' @param attribute The type of ID to convert to ensembl_gene_id and ensembl_transcript_id
#' @return A 3 column data.frame: ensembl_gene_id ensembl_transcript_id attribute(id chosen)
#' @examples
#' specie_genome = "hsapiens_gene_ensembl"
#' attribute <- "affy_hg_u95av2"
#' externalSymbol <- "hgnc_symbol"
#' ensemblID <- conversionID(specie_genome,attribute,externalSymbol)
#' @export
convert_id <- function(specie_genome = "rnorvegicus_gene_ensembl", attribute = "efg_agilent_wholegenome_4x44k_v1", externalSymbol = "rgd_symbol"){
  a<-system(paste("perl biomart-perl/scripts/converter.pl",specie_genome,attribute,externalSymbol,sep=" "),intern=TRUE) #IF ignore.stderr is removed their is often an external 500 error
  write.table(a[-1],"biomart-perl/scripts/affyId.txt",quote=FALSE,col.names=FALSE,row.names=FALSE)
  IDs<-read.table("biomart-perl/scripts/affyId.txt",sep="\t")
  colnames(IDs)<-c("ensembl_gene_id","ID")
  return(IDs)
}
