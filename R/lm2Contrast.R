#' Linear Model analysis with limma
#' @author Simon J Pelletier
#' @import limma
#' @param expr.matrix A matrix of standardized data. Columns = samples, rows = genes,transcripts,CpG...
#' @param design The design of the experiment (object)
#' @return Linear model contrast ...
#' \describe{
#'  \item{lm2.contrast}{MArrayLM object(limma) result of the eBayes function.}
#'  \item{contrasts}{Vector of contrasts}
#'  \item{contrast.matrix}{Matrix of contrasts}
#' }
#' @keywords limma linear
#' @seealso
#' \code{\link{design_contrasts}}
#' \code{\link{makeContrasts}}
#' \code{\link[limma]{lmFit}}
#' \code{\link[limma]{MArrayLM}}
#' @examples
#' expr.matrix <- readRDS("data/expr_matrix_LGVD.rds")
#' design <- design_contrasts(get_names(colnames(expr.matrix)))
#' results_list <- results_topTable()
#' lm2 <- lm2Contrast(expr.matrix,design)
#' lm2.contrast = lm2[[1]]
#' contrasts=lm2[[2]]
#' contrast.matrix=lm2[[3]]
#' @export
lm2Contrast <- function(expr.matrix,design){
  lm <- lmFit(expr.matrix, design)
  colnames(design) = as.character(colnames(design))
  contrasts = ""
  n = 1
  for (i in 1:(length(colnames(design))-1)){
    for (j in (i+1):length(colnames(design))){
      contrasts[n] = paste0(colnames(design)[i],"-",colnames(design)[j])
      n = n + 1
    }
  }
  contrast.matrix <-makeContrasts(contrasts=contrasts, levels = design)
  lm2 = contrasts.fit(lm, contrast.matrix)
  lm2.contrast <- eBayes(lm2)

  return(list(lm2.contrast,contrasts,contrast.matrix))
}
