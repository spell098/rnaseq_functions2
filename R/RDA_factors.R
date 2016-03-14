#' RDA plots
#' Creates a biplot
#' @param expr.matrix A matrix of standardized data. Columns = samples, rows = genes,transcripts,CpG...
#' @examples
#' expr.matrix <- readRDS('data/expr_matrix_LGVD.rds')
#' #' fileContent <- "Expression Matrix"

#' RDA_factors(expr.matrix[1:100,])
#' @export
RDA_factors <- function(expr.matrix,names){
  colnames(expr.matrix) <- names
  names.unique = unique(names)
  labs = as.data.frame(do.call("rbind", strsplit(names,"_")))
  colnames(labs) = apply(labs,2,function(x){
    paste(unique(x),collapse="_")
  })
  #}
  mod0 <- rda(t(expr.matrix) ~ 1, labs)
  mod1 <- rda(t(expr.matrix) ~ ., labs)
  mod <- ordiR2step(mod0, scope = formula(mod1),direction="forward",R2scope = TRUE)
  #plot(mod)
  #points(mod, display = "sites", cex = 0.8, pch=19, col=color_groups(names))
  #legend("topright",legend=names.unique,cex=0.6,col = color_groups(names.unique),lwd=1.5)
  #for(i in 1:length(unique(colnames(expr.matrix)))){
  #  ordiellipse(mod,groups = names,show.groups = unique(names)[i],conf = 0.95,col = color_groups(unique(names))[i])
  #}

  #text(mod, display = "sites", cex = 0.8, pch=19, col=color_groups(names),pos=3)
  #step.res <- ordiR2step(mod0, mod1, perm.max = 200)
  #step.res$anova  # Summary table
  #plot(step.res)
  return(mod)
}



