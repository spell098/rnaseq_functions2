#' Generate boxplots of every groups for a selected gene. Also display stars for significance between groups
#' @author Simon J Pelletier
#' @param selectedGene Gene selected to display boxplots for each groups compared
#' @param names.unique A vector of the names of all groups
#' @param names A vector of the names of all samples
#' @param expr.matrix A matrix of standardized data. Columns = samples, rows = genes,transcripts,CpG...
#' @param resultsSummary A matrix of data from every gene that is significanlty different
#' between at least two groups. Columns = samples, rows = genes,transcripts,CpG...
#' @param results The results of every comparisons
#' @param color The colors of all boxplot. One color for each group. Could be found in the function, not imported
#' @examples
#' gset <- getGEO("GSE12654", GSEMatrix =TRUE)
#' exprset <- gset[[1]]
#' typeID <- "ID"
#' expr.matrix <- exprs(exprset)
#' names1 <- get_names(sampleNames(exprset))[[1]]
#' comparisons <- comparisonsPheno(exprset)[[1]]
#' comparisonsTable <- comparisonsPheno(exprset)[[2]]
#' annotation1 <- annotation_geo(annotation(exprset))
#' goGPL <- annotation1[,c(colnames(annotation1)[grep("ontology",colnames(annotation1),ignore.case=TRUE)],"symbol")]
#' write.table(goGPL,"goTable.txt",sep=";",quote=FALSE,col.names=FALSE,row.names=FALSE)
#' annotations <- annotation1[,c("ID","symbol","entrezgene_id")]
#' system("bash annotateGPL.sh")
#' go<-read.csv("goTable3.txt")
#' genes_annotation_unique <- annotations[!duplicated(as.character(annotations[,typeID])),]
#' pvalue = 0.05
#' logFC = c(-2,2)
#' clusters <- fuzzy_clustering(expr.matrix,3)[[1]]
#' color = unlist(lapply(clusters,function(x){wes()[x]}))
#' namesColor = paste0("a",clusters)
#' comparisons = get_comparisons(namesColor)
#' clusterNames = numberNames(namesColor)
#' selectedVariables = comparisons[1]
#' names2 = conditionsChoice(selectedVariables,namesColor)
#' design = design_contrasts(names2)
#' expr.toBind <- NULL
#' lm2=lm2Contrast(expr.matrix,design)
#' lm2.contrast = lm2[[1]]
#' contrasts=lm2[[2]]
#' contrast.matrix=lm2[[3]]
#' results_list = results_topTable(lm2.contrast,expr.toBind,pvalue,logFC,typeID,genes_annotation_unique,annotations,"BH")
#' results = results_list[[1]]
#' topTable3 = results_list[[2]]
#' resultsSummary = results_summary(results,topTable3,"BH",typeID)
#' selectedGene <- rownames(resultsSummary[1,])
#' boxplot_element(selectedGene,names2,expr.matrix,resultsSummary,results,color)
#'
#'
#' expr.matrix <- readRDS("data/expr_matrix_LGVD.rds")
#' results <- readRDS("data/expr_matrix_LGVD.rds")
#' selectedGene <- "ENSRNOG00000046319"
#' topTable3 <- readRDS("data/topTable3_LGVD.rds")
#' resultsSummary <- results_summary(results,topTable3)
#' color <- c("green","blue","red")
#' boxplot_element(selectedGene,names2,expr.matrix,resultsSummary,results,color,grouped=FALSE)
#' @keywords geo annotation
#' @seealso
#' \code{\link[GEOquery]{getGEO}} ,
#' \code{\link[GEOquery]{Meta}} ,
#' \code{\link[GEOquery]{Table}} ,
#' \code{\link[annotate]{readGEOAnn}}
#' @note The function should require less input.
#' names and names.unique could be removed by changing the code
#' only the line from the selected gene should be imported, not every object to select it.
#' This would remove expr.matrix,resultsSummary,results,selectedGene
#' colors could also be found in the function. this would leave only a vector from expr.matrix(data) and a vector from resultsSummary(significance)
#' @export
boxplot_element = function(selectedGene,names,expr.matrix,resultsSummary,comparisons,color,grouped=TRUE,nbins=20){
  names.unique <- unique(names)
  if(grouped == TRUE){
    nrows=NULL
    for(i in 1:length(names.unique)){
      nrows[i]=length(grep(names.unique[i],names))
    }
    nrowMax=max(nrows)
    expr=matrix(ncol=length(names.unique),nrow=nrowMax)
    for(i in 1:length(names.unique)){
      expr[1:nrows[i],i] = expr.matrix[selectedGene,grep(names.unique[i],names)]
    }
    colnames(expr) = names.unique
    colors <- color_groups(names.unique)
  } else if(grouped==FALSE){
    expr=matrix(ncol=length(names),nrow=length(selectedGene))
    expr[1:length(selectedGene),] <- expr.matrix[selectedGene,]
    colnames(expr) <- names
    colors <- color_groups(names.unique)
  }
  pvalues = NULL
  maxVal <- max(expr[!is.na(expr)])
  height = maxVal + 0.1*(maxVal)
  if(grouped == TRUE){
    #expr2 <- data.frame(cbind(names=as.character(names),values=as.numeric(expr)))
    #expr2[,"values"] <- as.numeric(as.character(expr2[,"values"]))
    #expr2[,"names"] <- as.character(expr2[,"names"])
    boxplot(expr,ylim = c(min(expr[!is.na(expr)]),height),col=unique(colors),main=selectedGene)
    #g <- ggplot(expr2, aes(x=names,y=values,fill=names)) +
    #  geom_boxplot() 
    #return(g)
  } else {
    expr2 <- data.frame(cbind(names=as.character(names),values=as.numeric(expr)))
    expr2[,"values"] <- as.numeric(expr2[,"values"])
    expr2[,"names"] <- as.character(expr2[,"names"])
    
    g <- ggplot(expr2, aes(x=values, fill=names)) +
      geom_histogram(bins=nbins) +
      geom_vline(aes(xintercept=mean(values, na.rm=T)),   # Ignore NA values for mean
                 color="red", linetype="dashed", size=1)
    
    return(g)
    #x <- vector("list",length(names.unique))
    #names(x) <- names.unique
    #for(i in 1:length(names.unique)){
    #  x[[i]] <- expr[,grep(names.unique[i],names)]
    #}
    #stripchart(x,method = "jitter", pch=20,col = colors, cex = par("cex"),)
  }
  if(grouped==TRUE && !is.na(match(selectedGene,rownames(resultsSummary))) > 0){
    coordinates = matrix(ncol=2,nrow=length(comparisons))
    rownames(coordinates) = names(comparisons)
    colnames(coordinates) = c("a","b")
    for(i in 1:nrow(coordinates)){
      coordinates[i,] = match(strsplit(comparisons[i],"-")[[1]],colnames(expr))
    }
    if(nrow(resultsSummary) > 0){
      for(i in 1:(ncol(resultsSummary)-5)){
        pvalues[i] = resultsSummary[match(selectedGene,rownames(resultsSummary)),i+5]
      }
      for(i in 1:nrow(coordinates)){
        if (pvalues[i] < 0.05){
          first = min(coordinates[i,])
          second=max(coordinates[i,])
          ht = maxVal+0.05*maxVal
          segments(y0=ht,y1=ht,x0=first+0.1,x1=second-0.1,lwd=1)
          if (pvalues[i] < 0.001){
            text(x=mean(c(first,second)), y=maxVal+0.08*maxVal, "***", cex=1.2)
          } else if (pvalues[i] < 0.01){
            text(x=mean(c(first,second)), y=maxVal+0.08*maxVal, "**", cex=1.2)
          } else text(x=mean(c(first,second)), y=maxVal+0.08*maxVal, "*", cex=1.2)
        }
        
      }
    }
  }
}
