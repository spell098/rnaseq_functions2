voom.matrix <-voom(expr.matrix,design,normalize.method = "quantile")$E
bnet = readRDS('data/modules_characteristics_ch1.rds')
expr.toBind = exprToBind(bnet,voom.matrix)

merged <- merge(genes_annotation_unique,ensemblTable,by="ensembl_gene_id")
selectedElements <- sample(1:nrow(voom.matrix),1000)
mod <- rda(t(voom.matrix[selectedElements,]))

v <- vector("list",1)
v2 <- as.data.frame(cbind(as.character(rownames(mod$CA$v)),mod$CA$v)[1:100,])

#
tmp <- as.data.frame(merge(merged,v2,by.x=colnames(merged)[4],by.y=colnames(v2)[1]))
v[[1]] <- as.data.frame(tmp[!duplicated(tmp[,attribute]),])
colnames(v[[1]])[1] <- attribute
rownames(v[[1]]) <- v[[1]][,attribute]
RowSideColors <- expr.toBind[match(as.character(v[[1]][,1]),rownames(expr.toBind)) ,"module"]
expr.matrix.signif <- signifRows(voom.matrix,as.character(v[[1]][,attribute]))
expr.matrix.signif.means <- meanMatrix(expr.matrix.signif,names2)
#resultsContrast <- signifRows(lm2.contrast$coefficients,resultsSummary)
print("Looking for gene ontologies")
selectedComparison <- "disease.state.control-disease.state.cocaine.addiction"
names(v) <- selectedComparison

ontology <- geneOntology(v,go,type,nrow(genes_annotation_unique))
topGO <- cbind(rownames(ontology[[selectedComparison]][[2]]),ontology[[selectedComparison]][[2]])
selectedOntology <- "binding"

if(!is.null(selectedOntology)){
  resultsOntology <- resultsByOntology(selectedOntology,results[[selectedComparison]],ontology[[selectedComparison]])
} else {
  resultsOntology <- NULL
}
print("Gene ontologies: done.")


selectedModule = 1
mask <- match(expr.toBind[rownames(expr.matrix.signif),],rownames(expr.toBind))
topTable1 <- expr.toBind[mask[!is.na(mask)],]
modulesTable1 <- modules_summary(topTable1,expr.toBind)
maskModule1 <- as.numeric(as.character(topTable1[,"module"])) == as.numeric(as.character(selectedModule))
resultsModule1 <- list(resultsModule=topTable1[maskModule,])
resultsModule1[[1]] <- merge(ensemblTable,resultsModule1[[1]],by="ID")
colnames(resultsModule1[[1]]) [2]<- "ensembl_gene_id"
ontologyModule1 <- geneOntology(resultsModule,go,type,length(rownames(voom.matrix)))
topModuleGO1 <- cbind(rownames(ontologyModule[[1]][[2]]),ontologyModule[[1]][[2]])
resultsModuleOntology1 <- resultsByOntology(selectedModule,resultsModule[[1]],ontologyModule[[1]])



modules <- modules_summary(results[[1]],expr.toBind)
maskModule <- match(as.numeric(as.character(1)),as.numeric(as.character(results[[1]]$module)))
resultsModuleDiff <- list(resultsModule=results[[1]][maskModule,])
ontologyModule <- geneOntology(resultsModuleDiff,go,type,length(rownames(voom.matrix)))
topModuleGO <- cbind(rownames(ontologyModule[[1]][[2]]),ontologyModule[[1]][[2]])
resultsModuleOntology <- resultsByOntology(topModuleGO1_rows_selected,resultsModule[[1]],ontologyModule[[1]])

modules$ontologyModule <- geneOntology(modules$resultsModuleDiff,genesInfos$go,input$type,length(rownames(v$voom.matrix)))

