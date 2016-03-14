#' Get the annotations used in a NCBI GEO dataset
#' @author Simon J Pelletier
#' @param expr.matrix A matrix of numeric values. Rows are genes, columns are samples
#' @param pval Alpha value (selected p-value)
#' @return A table with all the annotation for a platform
#' @examples
#' expr.matrix <- readRDS("data/expr_matrix_LGVD.rds")
#' pval <- 0.05
#' summaryClust <- clusteringSummary(expr.matrix,pval)
#'
#' gset <- getGEO("GSE44593", GSEMatrix =TRUE)
#' exprset <- gset[[1]]
#' expr.matrix <- exprs(exprset)
#' names1 <- get_names(sampleNames(exprset))
#' groupBy <- "group"
#' noCluster <- 3
#' fileContent <- "Expression Set"
#' names1 <- get_names(expr.matrix)
#' selectedVariables <- NULL
#' comparisonsTable <- comparisonsPheno(exprset)[[2]]
#' namesTable <- comparisonsTable
#' names2 <- namesSelectedComparisons(namesTable)
#' summaryClust <- clusteringSummary(expr.matrix,pval,names2)
#' @keywords fuzzy cluster clustering Hypergeometric
#' @seealso
#' \code{\link{fuzzy_clustering}} ,
#' \code{\link[stats]{phyper}}
#' @export
clusteringSummary = function(expr.matrix,pval,names){
  colnames(expr.matrix) <- numberNames(names)
  if(ncol(expr.matrix) < 20 ) ncluster = ncol(expr.matrix) else ncluster = 20
  fuzzy=vector("list",ncluster)
  for(i in 2:ncluster){
    print(paste0("Running #",i))
    validClusters=NULL
    dominantClusters=NULL
    fuzzy[[i]] = vector("list",4)
    names(fuzzy[[i]]) = c("clusters","clusters_table","correct","correct in subsample")
    clusters = fuzzy_clustering(expr.matrix,i)[[1]]
    clustersUnique = unique(clusters)
    sampleNamesUnique = unique(names(clusters))
    fuzzy[[i]][[1]]=clusters
    fuzzy[[i]][[2]] = as.data.frame(matrix(ncol=length(sampleNamesUnique)*2+4,nrow=length(clustersUnique)+3))
    colnames(fuzzy[[i]][[2]]) = c(sampleNamesUnique,"Dominant samples","% in cluster","% of dominant samples","total in cluster",paste("pvalue",sampleNamesUnique,sep="_"))
    rownames(fuzzy[[i]][[2]]) = c(clustersUnique,"Dominant cluster","% in dominant cluster","total samples")

    # j = rows
    # k & l = column
    for (j in 1:length(clustersUnique)){
      for(k in 1:length(sampleNamesUnique)){
        fuzzy[[i]][[2]][j,k]=length(grep(sampleNamesUnique[k],names(clusters[grep(clustersUnique[j],clusters)])))
      }
    }
    for(l in 1:length(sampleNamesUnique)){
      max = max(as.numeric(fuzzy[[i]][[2]][,l][1:length(clustersUnique)]))
      sum = sum(as.numeric(fuzzy[[i]][[2]][,l][1:length(clustersUnique)]))
      fuzzy[[i]][[2]][length(clustersUnique)+1,l] = rownames(fuzzy[[i]][[2]])[as.numeric(fuzzy[[i]][[2]][,l]) == max][1]
      fuzzy[[i]][[2]][length(clustersUnique)+2,l] = max/sum
      fuzzy[[i]][[2]][length(clustersUnique)+3,l] = sum(as.numeric(fuzzy[[i]][[2]][,l][1:length(clustersUnique)]))
    }
    for(j in 1:length(clustersUnique)){

      noCorrect=matrix(ncol=length(sampleNamesUnique),nrow=length(clustersUnique)) #calculate the ratio of the signif group in cluster
      noCorrect2 = noCorrect #calculate the ratio of signif samples of total samples
      colnames(noCorrect) = sampleNamesUnique
      rownames(noCorrect) = clustersUnique
      samplesValues = as.numeric(fuzzy[[i]][[2]][j,][1:length(sampleNamesUnique)])
      totalSamples = as.numeric(fuzzy[[i]][[2]][length(clustersUnique)+3,][1:length(sampleNamesUnique)])
      ratios = samplesValues/totalSamples
      maxRatio=max(ratios)
      sumValues=sum(samplesValues)

      fuzzy[[i]][[2]][j,length(sampleNamesUnique)+1] = paste(sampleNamesUnique[ratios==maxRatio],collapse=" & ")
      fuzzy[[i]][[2]][j,length(sampleNamesUnique)+2] = max(samplesValues)/sumValues
      fuzzy[[i]][[2]][j,length(sampleNamesUnique)+3] = max(samplesValues)/sumValues
      fuzzy[[i]][[2]][j,length(sampleNamesUnique)+4] = sumValues
      for(l in 1:length(sampleNamesUnique)){
        #N = nsamples-totalSamples #rest
        m = totalSamples[l]
        k = sumValues
        x=samplesValues[l]
        n=sum(totalSamples) - totalSamples[l]
        pvalue = 1-phyper(x,m,n,k)
        fuzzy[[i]][[2]][j,length(sampleNamesUnique)+4+l] = pvalue
        if(pvalue < pval) dominantClusters[l] = sampleNamesUnique[l]
      }
    }
    for (z in 1:length(clustersUnique)){
      for (y in 1:length(sampleNamesUnique)){
        pvalue = fuzzy[[i]][[2]][z,length(sampleNamesUnique)+4+y]
        if (pvalue < pval){
          validClusters[z] = sampleNamesUnique[y]
          noCorrect[z,y] = as.numeric(fuzzy[[i]][[2]][z,y])/as.numeric(fuzzy[[i]][[2]][nrow(fuzzy[[i]][[2]]),y])
          noCorrect2[z,y] = as.numeric(fuzzy[[i]][[2]][z,y])/ncol(expr.matrix)
        } else {
          noCorrect[z,y] = 0
          noCorrect2[z,y] = 0
        }
      }
    }
    if(!is.null(validClusters)){
      correctClusters = fuzzy[[i]][[2]][1:length(validClusters),][!is.na(validClusters),]
      totalSamples = NULL
      sampling = rep(0,length(sampleNamesUnique))
      names(sampling) = sampleNamesUnique
      ratioSamples = NULL
      for (z in 1:length(sampleNamesUnique)){
        totalSamples[z] = sum(as.numeric(correctClusters[,z]))
      }
      for (z in 1:length(sampleNamesUnique)){
        for (y in 1:nrow(correctClusters)){
          if(correctClusters[y,(z+length(sampleNamesUnique)+4)] < pval){
            sampling[z] = sampling[z]+as.numeric(correctClusters[y,z])
          }
        }
      }
      fuzzy[[i]][[3]] = sum(noCorrect2)
      ratioSamples = sampling/totalSamples
      names(ratioSamples) = sampleNamesUnique
      fuzzy[[i]][[4]] = ratioSamples
    } else {
      fuzzy[[i]][[4]] = rep(0,length(sampleNamesUnique))
    }

  }
  return(fuzzy)
}
