#' List of informations on Gene Ontology of selected genes
#' @author Simon J Pelletier
#' @param results A topTable with only the genes considered significant (according to
#' the p-value and logFC limits provided in the function results_topTable)
#' @param ngenes A list of WGCNA results
#' @param go An object that contain all gene ontology informations for every gene
#' @param typeID Type of IDs used as rows in the expression matrix.Default="ensembl_gene_id"
#' @return A list of gene ontology informations for all comparisons determinated my the
#' experimental design. For example, if there are 3 groups (a,b,c), 3 comparisons are
#' possible (a-b,a-c,b-c). The object return is a list containing 3 lists named "a-b", etc.
#' Each of these lists contain the following objects:
#' \describe{
#'  \item{geneGO.selected}{For every gene ontology that imply at least 1 significant gene, }
#'  \item{GOarray.percent}{}
#'  \item{GOarray}{}
#'  \item{GOarray_CellComponent}{}
#'  \item{GOarray_MolecularFunction}{}
#'  \item{GOarray_BiologicalProcess}{}
#' }
#' @keywords cytoscape export
#' @seealso
#' \code{\link{annotate_ensembl}}
#' \code{\link{annotate_geo}}
#' \code{\link{results_topTable}}
#' \code{\link[limma]{topTable}}
#' @examples
#' results <- readRDS("data/results_LGVD.rds")
#' go <- read.csv("annotations/go.csv",header=FALSE)
#' typeID <- "ensembl_gene_id"
#' colnames(go) = c(typeID,"go_term_name","go_domain")
#' ngenes <- 26689
#' ontologyInfos <- geneOntology(results,go,typeID,ngenes)
#' @export
geneOntology = function(resultsList,go,type="ensembl_gene_id",ngenes){
  #saveRDS(resultsList,"resultsList.rds")
  print("gene ontology started...")
  ontologyInfos <- vector("list",length(resultsList))
  names(ontologyInfos) <- names(resultsList)
  geneGO <- by(go$go_term_name,
               go[,type],
               function(x) as.character(x))
  GOgene <- by(go[,type],
               go$go_term_name,
               function(x) as.character(x))
  GOgene = GOgene[unique(unlist(geneGO))]
  GOarray = array(dim = length(GOgene))

  for (i in 1:length(resultsList)){
    if (length(resultsList[[i]]) > 0){
      print(paste0("running ","1/",i,": ",names(resultsList)))
      genes=as.character(resultsList[[i]][,type])
      significant_geneGO = geneGO[genes]
      significant_geneGO = significant_geneGO[!is.na(names(significant_geneGO))]
      significant_GO = unique(unlist(significant_geneGO))
      significant_GOlist = GOgene[significant_GO]
      GOarray <- lapply(significant_GO,function(x){
        names(significant_geneGO[grep(x,significant_geneGO,fixed = TRUE)])
      })
      mask <- significant_GO != ""
      GOarray <- GOarray[mask]

      names(GOarray) <- significant_GO[mask]
      GOarray = order_list(GOarray)

      a=match(genes,as.character(go[,type]))
      a=a[!is.na(a)]
      GOarray_CellComponent = GOarray[go[a,]$go_domain == "cellular_component"]
      GOarray_MolecularFunction = GOarray[go[a,]$go_domain == "molecular_function"]
      GOarray_BiologicalProcess = GOarray[go[a,]$go_domain == "biological_process"]

      #Find the most represented GO
      GO.percent = matrix(ncol = 5,nrow = length(GOarray))
      rownames(GO.percent) = names(GOarray)
      colnames(GO.percent) = c("#Hits","#Possible_genes","p-value","q-value","Hit_ratio")

      #Faire un lapply???
      if(length(GOarray) > 0){
        print("Looking for significant ontologies...")
        
        for (z in 1:length(GOarray)){
          # TO CORRECT P-VALUE
          # # of test = # of annotations tested
          # p-value of GO is calculated using an hypergeometric distribution
          hits <- length(GOarray[names(GOarray[z])][[1]])
          if(hits > 0){
            possibleHits <- length(significant_GOlist[names(GOarray[z])][[1]])
            GO.percent[names(GOarray)[z],]["#Hits"] <- hits
            GO.percent[names(GOarray)[z],]["#Possible_genes"] = possibleHits
            GO.percent[names(GOarray)[z],]["Hit_ratio"] = hits/possibleHits
            M = length(significant_GO)
            N = ngenes-possibleHits #total number of genes
            pvalue = 1-phyper(hits-1,possibleHits,N,nrow(resultsList[[i]]))
            if(length(pvalue) == 0) {
              pvalue=1
              print("oops")
            }
            GO.percent[names(GOarray)[z],]["p-value"] = pvalue
          } else if(hits <= 0){
            print("PROBLEM WITH #HITS FOUND; SHOULD'NT BE 0")
            GO.percent[names(GOarray)[z],]["p-value"] = 1
          }
        }
        x <- order(GO.percent[,"p-value"])
        GO.percent[x,] = GO.percent[x,]
        GO.percent[,"q-value"] = p.adjust(GO.percent[,"p-value"],method="BH")
      }
      print("Done.")
      ontologyInfos[[i]] = list(geneGO,GO.percent,GOarray,GOarray_CellComponent,GOarray_MolecularFunction,GOarray_BiologicalProcess)
      names(ontologyInfos[[i]]) = c("geneGO.selected","GOarray.percent","GOarray","GOarray_CellComponent","GOarray_MolecularFunction","GOarray_BiologicalProcess")
    }
  }
  return(ontologyInfos)
}
