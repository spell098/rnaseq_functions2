#' Convert affyId to ensemblID
#' @author Simon J Pelletier
#' @param affyId Vector of all affyId to be converted to ensemblID
#' @return The object resultsContrast annotated with an associated color for the selected comparison
#' @keywords convert
#' @export
affy2ensembl = function(affyID){
  write.csv(affyID,"annotations/affyID.csv",quote=FALSE)
  system('awk -F"," \'NR==FNR{c[$2]=1;next} c[$2] == 1 {print $1","$2}\' annotations/affyID.csv annotations/databases/affygene.csv > annotations/affy2ensembl.csv')
  ensemblID = read.csv("annotations/affy2ensembl.csv",header=FALSE)
  colnames(ensemblID) <- c("ensembl_gene_id","ID")
  return(ensemblID)
}
