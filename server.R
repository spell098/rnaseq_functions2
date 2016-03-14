# All large files are on ratking and retrieved using the system(wget ...,intern=TRUE) function
library(data.table)
library(rnaseqApp)
library(BiocGenerics)
library(annotate)
library(GEOquery)
library(Mfuzz)
library(vegan)
library(pathview)
library(foreach)
library(DT)
library(parallel)
library(doParallel)
library(doMC)
library(foreign)
library(limma)
library(gdata)
library(gplots)
library(WGCNA)
library(wesanderson)
library(cclust)
library(RCurl)
library(parallel)
library(plyr)
library(graphite)
library(SPIA)
library(org.Rn.eg.db)
library(Rgraphviz)
library(network)
library(shiny)
library(a4Base)
library(graphite)
library(ggplot2)
library(devtools)
library(DT)
library(biomaRt)
#WHEN THE DESIGN IS CHANGED (GROUPS CHANGED) THE SAMPLES NEED TO BE RENORMALIZED,
#THE PIPELINE HAS TO START FROM RAW DATA TO ALLOW THAT


shinyServer(
  function(input, output, session) {
    #ncores = detectCores()
    #cl = makeForkCluster(nnodes = ncores)
    #registerDoParallel(cl)
    #registerDoMC()
    
    ################################################################################
    ###                              RNASEQ                                      ###
    ################################################################################
    print("App Started")
    observe({
      if(input$submitRun > 0){
        ordDiff <- reactiveValues(diffOrd=NULL,maskOrdDiff=NULL,mask1=NULL,resultsOrdDiff=NULL)
        
        maskModule <- reactive({NULL})
        resultsModuleOntology <- reactive({NULL})
        ord <- reactiveValues(v2=NULL,tmp=NULL,mod=NULL, eigenvalues=NULL,RowSideColors=NULL,v=NULL,expr.matrix.signif=NULL,
                              expr.matrix.signif.means=NULL,resultsOntology=NULL,ontology=NULL,topGO=NULL,
                              selectedOntology=NULL)
        v <- reactiveValues(bnet=NULL,expr.toBind=NULL,expr.matrix = NULL,expr.matrix2 = NULL,exprset = NULL,
                            names1 = NULL,comparisonsTable = NULL,unduplicated=NULL, comparisonsTable2 = NULL,
                            voom.matrix=NULL)
        
        groupsInfos <- reactiveValues(names1=NULL,design = NULL, groups = NULL, names2 = NULL, comparisons = NULL, 
                                      color = NULL, namesColor = NULL, replicates = NULL, ensemblTable = NULL,
                                      possibleComp=NULL,possibleComp=NULL)
        
        genesInfos <- reactiveValues(annotationFile = NULL,genes=NULL,annotation=NULL,annotations=NULL,go=NULL,
                                     genes_annotation_unique=NULL,merged=NULL)
        
        analysis <- reactiveValues(resultsSummary=NULL,expr.matrix.signif=NULL,expr.matrix.signif.means=NULL,RowSideColors=NULL,
                                   results_list=NULL,topTable3=NULL,results=NULL,lm2=NULL,lm2.contrast=NULL,
                                   contrasts=NULL,contrast.matrix=NULL,summaryClust=NULL,logFC=NULL)
        
        geneOntology <- reactiveValues(selectedOntology=NULL,ontology=NULL,topGO=NULL,resultsOntology=NULL)
        
        modules <- reactiveValues(modulesTable = NULL,resultsModule=NULL,ontologyModule=NULL,topModuleGO=NULL,
                                  resultsModuleOntology=NULL,maskModule=NULL)
        modules2 <- reactiveValues(modulesTable = NULL,resultsModule=NULL,ontologyModule=NULL,topModuleGO=NULL,
                                   resultsModuleOntology=NULL,maskModule=NULL) 
        
        modules3 <- reactiveValues(modulesTableOrd = NULL,resultsModule=NULL,ontologyModule=NULL,topModuleGO=NULL,
                                     resultsModuleOntology=NULL,maskModule=NULL, topTable=NULL,mask=NULL)
        modules4 <- reactiveValues(resultsModule=NULL,ontologyModule=NULL,topModuleGO=NULL,
                                   resultsModuleOntology=NULL,maskModule=NULL, topTable=NULL)
        network <- reactiveValues(selectedReaction = NULL,p=NULL,pSymbol=NULL,topTableNodes=NULL,
                                  reactome=NULL,networkSelectedOntology = NULL,topReactions=NULL,
                                  selectedTopReactions=NULL,toptableNodes1=NULL,network1=NULL)
        
        col_clu=redgreen(100)
        col_clu = col_clu[100:1]
        topReactions = reactive({NULL})
        modulesTable = reactive({NULL})
        summaryClust = reactive({NULL})
        geo_id = reactive({NULL})
        comparisons <- reactive({NULL})
        clusters1=reactive({NULL})
        observeEvent(input$submitRun,{
          print("Loading Data...")
          if(input$fileContent == "Expression Matrix"){
            v$expr.matrix <- readRDS(input$File1[,"datapath"])
            v$names1 <- get_names(colnames(v$expr.matrix))
          } else if(input$fileContent == "Expression Set"){
            if(length(input$geo_id) > 0 ){
              gset <- getGEO(input$geo_id, GSEMatrix =TRUE, getGPL = FALSE)  #GSE12654 GSE61276 illumina_humanht_12_v4
              v$exprset <- gset[[1]]
            }
            v$names1 <- get_names(sampleNames(v$exprset))
            v$comparisonsTable <- comparisonsPheno(v$exprset)[[2]]
            v$expr.matrix <- exprs(v$exprset)
          }
          print("Data loaded.")
        })
        
        
        observeEvent(input$submitParam,{
          if(!is.null(input$selectedReplicates) && input$selectedReplicates != ""){
            print("Unreplicating the expression matrix")
            v$unduplicated <- unduplicate(v$comparisonsTable,input$selectedReplicates,v$expr.matrix)
            print("Done Unreplicating")
            v$expr.matrix2 <- v$unduplicated[[1]]
            v$names1 <- colnames(v$expr.matrix2)
            v$comparisonsTable2 <- v$unduplicated[[2]]
          } else {
            print("No technical replication used")
            v$comparisonsTable2 <- v$comparisonsTable
            v$expr.matrix2 <- v$expr.matrix
          }
          if(input$randomSamples == TRUE){
            selectedElements <- sample(1:nrow(v$expr.matrix2),input$randomSampling)
            v$expr.matrix2 <- v$expr.matrix2[selectedElements,]
          } 
          
        })
        observeEvent(input$submitParam,{
          print("Retrieving groups informations...")
          if((input$groupBy == "groups")){
            groupsInfos$names1 <- groupsInfos$names2
          } 
          if(input$groupBy == "fuzzy" && is.null(groupsInfos$names1)){
            groupsInfos$names1 <- v$names1
          }
          
          groupsInfos$groups <- regroupingSamples(input$fileContent,input$groupBy,v$expr.matrix2,
                                                  input$selectedVariables,input$noCluster,v$names1,v$comparisonsTable2,
                                                  v$exprset)
          groupsInfos$names2 <- groupsInfos$groups[[1]]
          groupsInfos$comparisons <- groupsInfos$groups[[2]]
          groupsInfos$color <- groupsInfos$groups[[3]]
          groupsInfos$namesColor <- groupsInfos$groups[[4]]
          groupsInfos$possibleComp <- as.character(possibleComparisons(groupsInfos$names2))
          names(groupsInfos$possibleComp) <- groupsInfos$possibleComp
          colnames(v$expr.matrix2) <- groupsInfos$names2
          print("Done.")
        })
        observeEvent(input$submitParam,{
          genesInfos$annotationFile <- paste0("annotations/databases/",input$specie,"_ensembl_",input$attribute,".csv")
          if(file.exists(genesInfos$annotationFile)){
            print("annotation file already exists")
            groupsInfos$ensemblTable <- read.csv(genesInfos$annotationFile)
          } else {
            print("Translating to ensembl gene IDs")
            groupsInfos$ensemblTable <- annotation_biomart(rownames(v$expr.matrix),input$specie,input$attribute,input$mart)
            print("Transclation done.")
          }
          
          genesInfos$genes <- groupsInfos$ensemblTable[as.character(groupsInfos$ensemblTable[,input$attribute]) != "",]
          
          genesInfos$annotation <- annotate_ensembl(genesInfos$genes[,1],input$type)
          genesInfos$annotations <- genesInfos$annotation[[1]]
          genesInfos$go <- genesInfos$annotation[[2]]
          genesInfos$genes_annotation_unique <- genesInfos$annotation[[3]]
          
          groupsInfos$design <- design_contrasts(groupsInfos$names2)
          genesInfos$merged <- merge(genesInfos$genes_annotation_unique,groupsInfos$ensemblTable,by="ensembl_gene_id")
          
          if(input$standardization != "no"){
            v$voom.matrix <-voom(v$expr.matrix2,groupsInfos$design,normalize.method = input$standardization)$E
          } else {
            v$voom.matrix <- v$expr.matrix2
          }
          print("Normalization done.")
          labs = reactive({do.call("rbind", strsplit(groupsInfos$names2,"_"))})
          if (!is.null(input$File2)){
            v$bnet = readRDS(input$File2[,"datapath"])
            v$expr.toBind = exprToBind(v$bnet,v$voom.matrix)
          } else if(input$submitModule < 1){
            print("No modules")
            v$expr.toBind = cbind(module=rep("white",nrow(v$voom.matrix)),v$voom.matrix)
          }
        })
        
        observeEvent(input$submitModule,{
          print("Calculating modules...")
          v$bnet <- WGCNA_modules(v$voom.matrix,input$pow,groupsInfos$names2)
          saveRDS(v$bnet,paste0("data/modules_",paste(input$selectedVariables,collapse = "_"), ".rds"))
          saveRDS(v$voom.matrix,"voommatrix.rds")
          v$expr.toBind <- exprToBind(v$bnet,v$voom.matrix)
          #v$expr.toBind[,"module"] <- color_groups(v$expr.toBind[,"module"])
          saveRDS(v$expr.toBind,"exprtobind1.rds")
          print("Modules calculated.")
        })
        
        observeEvent(input$submitOrdination,{
          print(paste0("Generating ",input$ordination))
          if(input$random == TRUE){
            selectedElements <- sample(1:nrow(v$voom.matrix),input$randomSamplingOrd)
          } else {
            selectedElements <- 1:nrow(v$voom.matrix)
          }
          if(input$ordination == "Principal Component Analysis"){
            
            ord$mod <- rda(t(v$voom.matrix[selectedElements,]))
          } else if (input$ordination == "Redundency analysis"){
            ord$mod <- RDA_factors(v$voom.matrix,groupsInfos$names2,input$fileContent)
          }
          print("Done. Now generating the graphic")
          ord$eigenvalues <- data.frame(t(ord$mod$CA$eig/ord$mod$CA$tot.chi))
          ord$v <- vector("list",1)
          ord$v2 <- as.data.frame(cbind(as.character(rownames(ord$mod$CA$v)),ord$mod$CA$v)[1:100,])
          
          #
          ord$tmp <- as.data.frame(merge(genesInfos$merged,ord$v2,by.x=colnames(genesInfos$merged)[4],by.y=colnames(ord$v2)[1]))
          ord$v[[1]] <- as.data.frame(ord$tmp[!duplicated(ord$tmp[,input$attribute]),])
          colnames(ord$v[[1]])[1] <- input$attribute
          rownames(ord$v[[1]]) <- ord$v[[1]][,input$attribute]
          ord$RowSideColors <- v$expr.toBind[match(as.character(ord$v[[1]][,1]),rownames(v$expr.toBind)) ,"module"]
          ord$expr.matrix.signif <- signifRows(v$voom.matrix,as.character(ord$v[[1]][,input$attribute]))
          ord$expr.matrix.signif.means <- meanMatrix(ord$expr.matrix.signif,groupsInfos$names2)
          #ord$resultsContrast <- signifRows(ord$lm2.contrast$coefficients,ord$resultsSummary)
          print("Looking for gene ontologies")
          names(ord$v) <- input$selectedComparison
          
          ord$ontology <- geneOntology(ord$v,genesInfos$go,input$type,nrow(genesInfos$genes_annotation_unique))
          ord$topGO <- cbind(rownames(ord$ontology[[input$selectedComparison]][[2]]),ord$ontology[[input$selectedComparison]][[2]])
          ord$selectedOntology <- c(input$topGO_ord_rows_selected,strsplit(input$selectedOntology_ord,"\n")[[1]])
          
          if(!is.null(ord$selectedOntology)){
            ord$resultsOntology <- resultsByOntology(ord$selectedOntology,analysis$results[[input$selectedComparison]],geneOntology$ontology[[input$selectedComparison]])
          } else {
            ord$resultsOntology <- NULL
          }
          print("Gene ontologies: done.")
          
        })
        
        observeEvent(input$submit2,{
          analysis$RowSideColors <- rep("white",nrow(v$voom.matrix))
          ord$RowSideColors <- rep("white",nrow(v$voom.matrix))
          print("Preparing differential expressions...")
          analysis$lm2 <- lm2Contrast(v$voom.matrix,groupsInfos$design)
          analysis$lm2.contrast <- analysis$lm2[[1]]
          analysis$contrasts <- analysis$lm2[[2]]
          analysis$contrast.matrix <- analysis$lm2[[3]]
          logfc <- c(input$logFCneg,input$logFCpos)
          print("Calculating differential expressions...")
          
          analysis$results_list <- results_topTable(analysis$lm2.contrast,v$expr.toBind,input$pvalue,
                                                    logFC = logfc,input$type,genesInfos$genes_annotation_unique,
                                                    genesInfos$annotations,input$adjust,groupsInfos$ensemblTable)
          
          print("Differential expressions done.")
          analysis$results <- analysis$results_list[[1]]
          analysis$topTable3 <- analysis$results_list[[2]]
          if(nrow(do.call("rbind", analysis$results)) > 0){
            #print(do.call("rbind", analysis$results))
            #print(head(v$expr.toBind))
            analysis$resultsSummary <- results_summary(analysis$results,analysis$topTable3,input$adjust)
            if(length(analysis$resultsSummary) > 0 && !is.null(analysis$resultsSummary)){
              analysis$expr.matrix.signif <- signifRows(v$voom.matrix,rownames(analysis$resultsSummary))
              analysis$expr.matrix.signif.means <- meanMatrix(analysis$expr.matrix.signif,groupsInfos$names2)
              analysis$resultsContrast <- signifRows(analysis$lm2.contrast$coefficients,rownames(analysis$resultsSummary))
              print("Looking for gene ontologies")
              selectedOntology <- c(input$topGO_rows_selected,strsplit(input$selectedOntology,"\n")[[1]])
              
              geneOntology$ontology <- geneOntology(analysis$results,genesInfos$go,input$type,nrow(genesInfos$genes_annotation_unique))
              geneOntology$topGO <- cbind(rownames(geneOntology$ontology[[input$selectedComparison]][[2]]),geneOntology$ontology[[input$selectedComparison]][[2]])
              if(!is.null(selectedOntology)){
                geneOntology$resultsOntology <- resultsByOntology(selectedOntology,analysis$results[[input$selectedComparison]],geneOntology$ontology[[input$selectedComparison]])
              } else {
                geneOntology$resultsOntology <- NULL
              }
              print("Gene ontologies: done.")
              #networkSelectedComparison = reactive({export2cytoscape(v$voom.matrix,analysis$results[[input$selectedComparison]],input$threshold,input$type)})
            } else {
              output$warnings <- renderText({
                print("No significant genes were found.")
              })
            } 
          } else {
            output$warnings <- renderText({
              print("No significant genes were found.")
            })
          }
        })
        ##########################################################################
        ##                Differential Analysis Ordination Analysis             ##
        ##########################################################################
        
        observe({
          if (length(unique(v$expr.toBind[,"module"])) > 1 && !is.null(analysis$resultsSummary) && length(analysis$resultsSummary) > 0){
            observe({
              print("Looking for differentially expressed genes in top genes that explain the variance")
              ordDiff$mask1 <- match(rownames(ord$expr.matrix.signif),analysis$results[[1]][,"ID"])
              ordDiff$diffOrd <- analysis$results[[input$selectedComparison]][ordDiff$mask1[!is.na(ordDiff$mask1)]]
              print("done.")
            })
              observe({
                print(paste0("Loading differential and top variant gene ",as.character(input$diffOrd_row_last_clicked)))
                #as.numeric(as.character(analysis$results[[input$selectedComparison]]$module)) == as.numeric(as.character(input$modulesTable_rows_selected))
                ordDiff$maskOrdDiff <- as.character(analysis$results[[input$selectedComparison]][,"ID"]) == as.character(input$diffOrd_row_last_clicked)
                ordDiff$resultsOrdDiff <- list(resultsModule=analysis$results[[input$selectedComparison]][ordDiff$maskOrdDiff,])
                #print(ordDiff$resultsOrdDiff[[1]])
                saveRDS(ordDiff$resultsOrdDiff,file="data/resultsModule.rds")
                print("...")
                ordDiff$ontologyOrdDiff <- geneOntology(ordDiff$resultsModule,genesInfos$go,input$type,length(rownames(v$voom.matrix)))
                ordDiff$topOrdDiffGO <- cbind(rownames(ordDiff$ontologyModule[[1]][[2]]),ordDiff$ontologyModule[[1]][[2]])
                print("Done.")
              })  
              observe({
                if(length(input$topOrdDiffGO_rows_selected) > 0){
                  if(nrow(ordDiff$topOrdDiffGO) > 0){
                    ordDiff$resultsOrdDiffOntology <- resultsByOntology(input$topOrdDiffGO2_rows_selected,ordDiff$resultsOrdDiff[[1]],ordDiff$ontologyOrdDiff[[1]])
                  } else print("problem module 1: empty topModuleGO")
                } else print("problem module 1: topModuleGO length 0")
                
              })
              
            
          } else print("No results or ordination")

        })
        
        
        ##########################################################################
        ##                Differential Analysis Module Analysis                 ##
        ##########################################################################
        observe({
          if (length(unique(v$expr.toBind[,"module"])) > 1 && !is.null(analysis$resultsSummary) && length(analysis$resultsSummary) > 0){
              observe({
                print("Looking for modules enriched with differnatially expressed genes")
                modules$modulesTable <- modules_summary(analysis$results[[input$selectedComparison]],v$expr.toBind)
                print("done.")
              })
              if(!is.null(modules$modulesTable)){
                observe({
                  print(paste0("Loading differential modules ",as.numeric(as.character(input$modulesTable_row_last_clicked))," : selected genes..."))
                  #as.numeric(as.character(analysis$results[[input$selectedComparison]]$module)) == as.numeric(as.character(input$modulesTable_rows_selected))
                  modules$maskModule <- as.character(analysis$results[[input$selectedComparison]][,"module"]) == as.character(input$modulesTable_row_last_clicked)
                  modules$resultsModule <- list(resultsModule=analysis$results[[input$selectedComparison]][modules$maskModule,])
                  #print(modules$resultsModule[[1]])
                  saveRDS(modules$resultsModule,file="data/resultsModule.rds")
                  print("...")
                  modules$ontologyModule <- geneOntology(modules$resultsModule,genesInfos$go,input$type,length(rownames(v$voom.matrix)))
                  modules$topModuleGO <- cbind(rownames(modules$ontologyModule[[1]][[2]]),modules$ontologyModule[[1]][[2]])
                  print("Done.")
                })  
                observe({
                  print(paste0("Loading differential modules ", as.numeric(as.character(input$modulesTable_row_last_clicked))))
                  modules2$maskModule <- as.character(analysis$topTable3[[input$selectedComparison]][,"module"]) == as.character(input$modulesTable_row_last_clicked)
                  modules2$resultsModule <- list(resultsModule=analysis$topTable3[[input$selectedComparison]][modules2$maskModule,])
                  modules2$ontologyModule <- geneOntology(modules2$resultsModule,genesInfos$go,input$type,length(rownames(v$voom.matrix)))
                  modules2$topModuleGO <- cbind(rownames(modules2$ontologyModule[[1]][[2]]),modules2$ontologyModule[[1]][[2]])
                  print("Done.")
                })
                observe({
                  if(length(input$topModuleGO1_rows_selected) > 0){
                    if(nrow(modules$topModuleGO) > 0){
                      modules$resultsModuleOntology <- resultsByOntology(input$topModuleGO1_rows_selected,modules$resultsModule[[1]],modules$ontologyModule[[1]])
                    } else print("problem module 1: empty topModuleGO")
                  } else print("problem module 1: topModuleGO length 0")
                  
                })
                observe({
                  if(length(input$topModuleGO2_rows_selected) > 0){
                    if(nrow(modules2$topModuleGO) > 0){
                      modules2$resultsModuleOntology <- resultsByOntology(input$topModuleGO2_rows_selected,modules2$resultsModule[[1]],modules2$ontologyModule[[1]])
                    } else print("problem module 2: empty topModuleGO")
                  } else print("problem module 2: topModuleGO length 0")
                })
            }
          }
        })
        
        
        ##########################################################################
        ##                Ordination Analysis Module Analysis                   ##
        ##########################################################################
        observe({
          if (length(unique(v$expr.toBind[,"module"])) > 1 && !is.null(analysis$topTable3) && !is.null(ord$v)){
            observe({
              print("Looking if top genes from PCA are enriched in modules")
              modules3$mask <- match(v$expr.toBind[rownames(ord$expr.matrix.signif),],analysis$topTable3[[input$selectedComparison]]$ID)
              modules3$topTable <- analysis$topTable3[[input$selectedComparison]][modules3$mask[!is.na(modules3$mask)],]
              modules3$modulesTableOrd <- modules_summary(modules3$topTable,v$expr.toBind)
              print("done") 
            })
            if(!is.null(modules3$modulesTableOrd)){
              observe({
                print(paste0("Loading ordination modules ",as.numeric(as.character(input$modulesTableOrd_row_last_clicked)),": selected genes..."))
                
                modules3$maskModule <- as.character(modules3$topTable[,"module"]) == as.character(input$modulesTableOrd_row_last_clicked)
                modules3$resultsModule <- list(resultsModuleTable=modules3$topTable[modules3$maskModule,])
                modules3$ontologyModule <- geneOntology(modules3$resultsModule,genesInfos$go,input$type,length(rownames(v$voom.matrix)))
                modules3$topModuleGO <- cbind(rownames(modules3$ontologyModule[[1]][[2]]),modules3$ontologyModule[[1]][[2]])
                print("done")
              })
              observe({
                print(paste0("Loading ordination module ",as.numeric(as.character(input$modulesTableOrd_row_last_clicked))))
                
                modules4$maskModule <- as.character(analysis$topTable3[[input$selectedComparison]][,"module"]) == as.character(input$modulesTableOrd_row_last_clicked)
                modules4$resultsModule <- list(resultsModuleTable=analysis$topTable3[[input$selectedComparison]][modules4$maskModule,])
                modules4$ontologyModule <- geneOntology(modules4$resultsModule,genesInfos$go,input$type,length(rownames(v$voom.matrix)))
                modules4$topModuleGO <- cbind(rownames(modules4$ontologyModule[[1]][[2]]),modules4$ontologyModule[[1]][[2]])
                #print(modules4$topModuleGO)
                print("Done.")
              })
              observe({
                if(length(input$topModuleGO3_rows_selected) > 0){
                  if(nrow(modules3$topModuleGO) > 0){
                    modules3$resultsModuleOntology <- resultsByOntology(input$topModuleGO3_rows_selected,modules3$resultsModule[[1]],modules3$ontologyModule[[1]])
                  } else print("problem module 3: empty topModuleGO")
                } else print("problem module 3: topModuleGO length 0")
              })
              observe({
                if(length(input$topModuleGO4_rows_selected) > 0){
                  if(nrow(modules4$topModuleGO) > 0){
                    modules4$resultsModuleOntology <- resultsByOntology(input$topModuleGO4_rows_selected,modules4$resultsModule[[1]],modules4$ontologyModule[[1]])
                  } else print("problem module 4: empty topModuleGO")
                } else print("problem module 4: topModuleGO length 0")
              })
            }
          }
        })
        
        
        
        ########################################################################
        
        observeEvent(input$submitSummaryClust,{
          print("Clustering...")
          analysis$summaryClust <- clusteringSummary(v$voom.matrix,input$pval_clust,groupsInfos$names2)
        })
        
        observeEvent(input$submitNetwork,{
          print("Network analysis started")
          #network$networkSelectedOntology = export2cytoscape(v$voom.matrix,geneOntology$resultsOntology,input$threshold,input$type)
          
          if (!is.null(input$File3)){
            network$topReactions = readRDS(input$File3[,"datapath"])
            network$selectedTopReactions = selectTopReaction(input$selectedComparison,network$topReactions)
          }
          network$reactome <- pathways(input$specie, "reactome")
          observeEvent(input$findTopReactions,{
            print("preparing for SPIA...")
            prepareSPIA(network$reactome, "reactome",print.names=TRUE)
            print("Looking for top reactions...")
            network$topReactions = reactions(analysis$results,unique(genesInfos$annotations$entrezgene_id))
            saveRDS(network$topReactions,paste0("data/reactions_",paste(input$selectedVariables,collapse = "_"), ".rds"))
          })
          
          network$selectedReaction <- c(strsplit(input$selectedReactions,"\n")[[1]],strsplit(network$networkSelectedOntology[input$topReactions_rows_selected,1],"\n")[[1]])
          network$p <- network$reactome[[network$selectedReaction]]
          network$pSymbol <- convertIdentifiers(network$p, "SYMBOL")
          toptableNodes <-  cbind(
            analysis$topTable3[[input$selectedComparison]][match(nodes(network$pSymbol),analysis$topTable3[[input$selectedComparison]]$symbol),],
            Reaction = rep(network$selectedReaction[[1]],length(nodes(network$pSymbol)))
          )
          network$toptableNodes1 = network$toptableNodes[!is.na(network$toptableNodes[,1]),]
          network$network1 <- network.graphite(network$pSymbol,v$voom.matrix,network$toptableNodes1,
                                               input$type,input$threshold,network$selectedReaction[[1]],analysis$topTable3,
                                               input$selectedComparison,NULL,input$pvalue)
          print("Network analysis: done.")
        })
        
        observe({
          updateSelectInput(session, "selectedComparison", choices = groupsInfos$possibleComp)
        })
        observe({
          updateSelectInput(session, "selectedModule", choices = as.numeric(head(modules$modulesTable[,"Module"])))
        })
        observe({
          updateSelectInput(session, "selectedVariables", choices = colnames(v$comparisonsTable), selected = colnames(v$comparisonsTable))
        })
        observe({
          updateSelectInput(session, "selectedReplicates", choices = colnames(v$comparisonsTable), selected = NULL)
        })
        observe({
          updateSelectInput(session, "selectedComponent", choices = colnames(v$comparisonsTable), selected = NULL)
        })
        output$boxplot <- renderPlot({
          boxplot(v$voom.matrix, col= groupsInfos$color, ylab = "expression level", las = 2,
                  cex.axis = input$cex.axis_boxplot,
                  names = colnames(v$voom.matrix),
                  ylim = input$ylim,cex=0.1)
          mtext("Labels", side=1, line=7)
          title("Boxplot")
        })
        ################################################################################
        ###                              OUTPUTS                                     ###
        ################################################################################
        output$text1 <- renderText({
          print(input$v_row_last_clicked)
        })
        ################################################################################
        ###                              PLOTS                                       ###
        ################################################################################
        output$boxplot_element <- renderPlot({
          boxplot_element(input$resultsSummary_row_last_clicked,groupsInfos$names2,
                          v$voom.matrix,analysis$resultsSummary,names(analysis$results),groupsInfos$color,grouped=TRUE)
        })
        output$hist_element <- renderPlot({
          boxplot_element(input$resultsSummary_row_last_clicked,groupsInfos$names2,
                          v$voom.matrix,analysis$resultsSummary,names(analysis$results),groupsInfos$color,grouped=FALSE,input$bins)
        })
        output$boxplot_element_pca <- renderPlot({
          boxplot_element(input$v_row_last_clicked,groupsInfos$names2,
                          v$voom.matrix,analysis$resultsSummary,names(analysis$results),groupsInfos$color,grouped=TRUE)
        })
        output$hist_element_pca <- renderPlot({
          boxplot_element(input$v_row_last_clicked,groupsInfos$names2,
                          v$voom.matrix,analysis$resultsSummary,names(analysis$results),groupsInfos$color,grouped=FALSE,input$bins)
        })
        output$boxplot_element_raw <- renderPlot({
          # Take a dependency on input$submit2
          input$submit2
          x=input$selectedGenes
          boxplot_element(x,groupsInfos$names2,v$voom.matrix,analysis$resultsSummary,analysis$topTable3,groupsInfos$color)
        })
        output$pathway <- renderPlot({
          # Take a dependency on input$submit2
          input$submit2
          
          renderGraph(network$network1)
        },height = 800, width = 1200)
        
        
        output$venn <- renderPlot({
          # Take a dependency on input$submit2
          input$submit2
          vennDiagram(analysis$resultsContrast,cex=input$venn.cex)
        })
        
        
        output$heatmap1 <- renderPlot({
          
          heatmap.2(analysis$expr.matrix.signif,labRow=as.character(analysis$resultsSummary$symbol), col=col_clu, scale="row",
                    main = "Heatmap result", key=TRUE, symkey=FALSE,
                    density.info="none", trace="none", cexRow=input$cexRow_heatmap,cex=1,
                    dendrogram="both", Colv = TRUE, Rowv = TRUE,
                    lhei=c(0.5,input$row_heatmap_height),lwid=c(0.5,2),
                    RowSideColors=color_groups(analysis$resultsSummary$module),
                    ColSideColors=color_groups(groupsInfos$names2))
        },height = 1000, width = 600)
        
        output$heatmap2 <- renderPlot({
          
          heatmap.2(analysis$expr.matrix.signif.means,labRow=analysis$resultsSummary$symbol, col=col_clu, scale="row",
                    main = "Heatmap result", key=TRUE, symkey=FALSE,
                    density.info="none", trace="none", cexRow=input$cexRow_heatmap,cex=0.8,
                    dendrogram=input$dendrogram, Colv = TRUE, Rowv = TRUE,
                    lhei=c(0.5,input$row_heatmap_height),lwid=c(0.5,2),
                    RowSideColors=color_groups(analysis$resultsSummary$module))
        },height = 1000, width = 600)
        
        output$heatmap1_ord <- renderPlot({
          # Take a dependency on input$submit2
          heatmap.2(ord$expr.matrix.signif,labRow=as.character(ord$v[[1]]$symbol), col=col_clu, scale="row",
                    main = "Heatmap result", key=TRUE, symkey=FALSE,
                    density.info="none", trace="none", cexRow=input$cexRow_heatmap_ord,cex=1,
                    dendrogram=input$dendrogram_ord, Colv = TRUE, Rowv = TRUE,
                    lhei=c(0.5,input$row_heatmap_height_ord),lwid=c(0.5,2),
                    RowSideColors=ord$RowSideColors)
        },height = 1000, width = 600)
        
        output$heatmap2_ord <- renderPlot({
          # Take a dependency on input$submit2
          input$submit2
          heatmap.2(ord$expr.matrix.signif.means,labRow=as.character(ord$v[[1]]$symbol), col=col_clu, scale="row",
                    main = "Heatmap result", key=TRUE, symkey=FALSE,
                    density.info="none", trace="none", cexRow=input$cexRow_heatmap_ord,cex=0.8,
                    dendrogram=input$dendrogram_ord, Colv = TRUE, Rowv = TRUE,
                    lhei=c(0.5,input$row_heatmap_height_ord),lwid=c(0.5,2),
                    RowSideColors=ord$RowSideColors)
        },height = 1000, width = 600)
        
        output$MDS <- renderPlot({
          # Take a dependency on input$submit2
          input$submit2
          plotMDS(v$voom.matrix,cex=input$cex.MDS,col=groupsInfos$color)
        })
        
        output$ordination <- renderPlot({
          #print(groupsInfos$names1)
          plot(ord$mod,type="n")
          points(ord$mod, display = "sites", cex = 1.5, pch=21, bg=color_groups(groupsInfos$names2),col=color_groups(groupsInfos$names1),lwd=3)
          legend("topright",legend=unique(groupsInfos$names2),cex=0.6,col = color_groups(unique(groupsInfos$names2)),lwd=1.5)
          if(length(unique(groupsInfos$names1)) < length(groupsInfos$names1)){
            legend("bottomright",legend=unique(groupsInfos$names1),cex=0.6,col = color_groups(unique(groupsInfos$names1)),lwd=1.5)
          }
          for(i in 1:length(unique(colnames(v$voom.matrix)))){
            ordiellipse(ord$mod,groups = groupsInfos$names2,show.groups = unique(groupsInfos$names2)[i],conf = 0.95,col = color_groups(unique(groupsInfos$names2))[i])
          }
        })
        
        output$eigenvalues <- DT::renderDataTable({
          #Eigenvalues for unconstrained axes
          print(ord$eigenvalues)
        })
        output$u <- DT::renderDataTable({
          print(data.frame(ord$mod$CA$u))
        })
        output$v <- DT::renderDataTable({
          ord$v[[1]]
        },selection="single")
        
        output$ordinationAnova <- DT::renderDataTable({
          print(data.frame(ord$mod$anova))
        })
        output$volcanoplot <- renderPlot({
          # Take a dependency on input$submit2
          input$submit2
          volcanoGgPlot(analysis$topTable3[[input$selectedComparison]],input$logFCVolcano,input$pvalueVolcano)
        })
        
        output$MA <- renderPlot({
          # Take a dependency on input$submit2
          input$submit2
          MAPlot(analysis$topTable3[[input$selectedComparison]],input$logFCVolcano,input$pvalueVolcano)
        })
        
        output$heatDiagram <- renderPlot({
          # Take a dependency on input$submit2
          input$submit2
          heatmap.2(analysis$resultsContrast, labRow=analysis$resultsSummary$symbol,col=col_clu, scale="row",
                    main = "Heatmap result", key=TRUE, symkey=FALSE,
                    density.info="none", trace="none", cexRow=input$cexRow_heatmap,cex=0.8,
                    dendrogram=input$dendrogram, Colv = TRUE, Rowv = TRUE,
                    lhei=c(0.5,input$row_heatmap_height),lwid=c(0.5,2),
                    RowSideColors=color_groups(analysis$resultsSummary$module))
        },height = 1000, width = 600)
        
        
        ################################################################################
        ###                              TABLES                                      ###
        ################################################################################
        output$resultsSummary <- DT::renderDataTable({
          # Take a dependency on input$submit2
          input$submit2
          print(analysis$resultsSummary)
        },selection="single")
        output$comparisons_table <- DT::renderDataTable({
          # Take a dependency on input$submit2
          input$submit2
          print(v$comparisonsTable2)
        })
        output$summaryClustering <- DT::renderDataTable({
          # Take a dependency on input$submit2
          input$submit2
          print(analysis$summaryClust[[(input$selectedClustn)]][[2]])
        })
        output$plotCorrectAttribution <- renderPlot({
          # Take a dependency on input$submit2
          input$submit2
          plotCorrectAttribution(analysis$summaryClust)
        })
        output$plotCorrectSubsample <- renderPlot({
          # Take a dependency on input$submit2
          input$submit2
          plotCorrectSubsample(analysis$summaryClust)
        })
        output$topTable <- DT::renderDataTable({
          print(analysis$results[[input$selectedComparison]])
        })
        
        output$pathwayTable <- DT::renderDataTable({
          print(network$toptableNodes1)
        })
        
        output$topGO_ord <- DT::renderDataTable({
          print(ord$topGO)
        })
        
        output$topGO <- DT::renderDataTable({
          print(geneOntology$topGO)
        })
        
        output$modulesInfos <- DT::renderDataTable({
          
        })
        
        
        output$selectedGO <- DT::renderDataTable({
          print(geneOntology$resultsOntology)
        })
        
        output$selectedGO_ord <- DT::renderDataTable({
          print(ord$resultsOntology)
        })
        
        output$modulesTable <- DT::renderDataTable({
          print(modules$modulesTable)
        },selection="single")
        
        ##
        ## 
        
         output$selectedModule1 <- DT::renderDataTable({
          print(modules$resultsModule[[1]])
        },selection="single")
        output$topModuleGO1 <- DT::renderDataTable({
          print(modules$topModuleGO)
        })
        output$resultsModuleOntology1 <- DT::renderDataTable({
          print(modules$resultsModuleOntology)
        })
        
        ##
        
        output$selectedModule2 <- DT::renderDataTable({
          print(modules2$resultsModule[[1]])
        },selection="single")
        output$topModuleGO2 <- DT::renderDataTable({
          print(modules2$topModuleGO)
        })
        output$resultsModuleOntology2 <- DT::renderDataTable({
          print(modules2$resultsModuleOntology)
        })
        
        ##
        
        output$resultsModuleOntology3 <- DT::renderDataTable({
          print(modules3$resultsModuleOntology)
        })
        output$modulesTableOrd <- DT::renderDataTable({
          print(modules3$modulesTableOrd)
        },selection="single")
        output$selectedModule3 <- DT::renderDataTable({
          print(modules3$resultsModule[[1]])
        })
        output$topModuleGO3 <- DT::renderDataTable({
          print(modules3$topModuleGO)
        })
        
        ##
        
        output$diffOrd <- DT::renderDataTable({
          print(ordDiff$diffOrd)
        })
        output$resultsOrdDiff <- DT::renderDataTable({
          print(ordDiff$resultsOrdDiff[[1]])
        })
        output$topOrdDiffGO <- DT::renderDataTable({
          print(ordDiff$topOrdDiffGO)
        })
        output$resultsOrdDiffOntology <- DT::renderDataTable({
          print(ordDiff$resultsOrdDiffOntology)
        })
        
        ## 
        
        output$selectedModule4 <- DT::renderDataTable({
          print(modules4$resultsModule[[1]])
        })
        output$topModuleGO4 <- DT::renderDataTable({
          print(modules4$topModuleGO)
        })
        output$resultsModuleOntology4 <- DT::renderDataTable({
          print(modules4$resultsModuleOntology)
        })
        
        
        ## 

        output$topReactions <- DT::renderDataTable({
          print(network$networkSelectedOntology)
        })
        output$toptableNodes <- DT::renderDataTable({
          # Take a dependency on input$submit2
          input$submit2
          #print(network$toptableNodes1)
        })
        output$rawTopTable <- DT::renderDataTable({
          # Take a dependency on input$submit2
          input$submit2
          print(analysis$topTable3[[input$selectedComparison]][match(input$selectedGenes,as.character(analysis$topTable3[[input$selectedComparison]][,input$type])),])
        })
        
        ################################################################################
        ###                              DOWNLOAD HANDLER                            ###
        ################################################################################
        
        
        output$downloadExprMatrix <- downloadHandler(
          filename = function() { paste0("expr_matrix", '.rds') },
          content = function(file) {
            saveRDS(v$voom.matrix, file)
          }
        )
        
        output$downloadSelectedGO <- downloadHandler(
          filename = function() { paste0("selectedGO", '.rds') },
          content = function(file) {
            saveRDS(geneOntology$resultsOntology, file)
          }
        )
        
        output$downloadSelectedGO_ord <- downloadHandler(
          filename = function() { paste0("selectedGO_ord", '.rds') },
          content = function(file) {
            saveRDS(geneOntology$resultsOntology_ord, file)
          }
        )
        
        output$downloadSummary <- downloadHandler(
          filename = function() { paste0("resultsSummary_",possibleComparisons(groupsInfos$names2), '.rds') },
          content = function(file) {
            saveRDS(analysis$resultsSummary, file)
          }
        )
        
        output$downloadResultsSelectedComparison <- downloadHandler(
          filename = function() { paste0("results_",input$selectedComparison, '.rds') },
          content = function(file) {
            saveRDS(analysis$results[[input$selectedComparison]], file)
          }
        )
        output$downloadAllResults <- downloadHandler(
          filename = function() { paste0("results_", '.rds') },
          content = function(file) {
            saveRDS(analysis$results, file)
          }
        )
        output$downloadNodeDataComparison <- downloadHandler(
          filename = function() { paste("node_rnaseq_",input$selectedComparison, '.csv', sep='') },
          content = function(file) {
            write.csv(networkinput$selectedComparison$nodeData, file)}
        )
        
        output$downloadEdgeDataComparison <- downloadHandler(
          filename = function() { paste("edge_rnaseq_",input$selectedComparison, '.csv', sep='') },
          content = function(file) {
            write.csv(networkinput$selectedComparison$edgeData, file)}
        )
        output$downloadNodeSelectedGO <- downloadHandler(
          filename = function() { paste("node_rnaseq_",input$selectedComparison,input$selectedOntology, '.csv', sep='') },
          content = function(file) {
            write.csv(network$networkSelectedOntology$nodeData, file)}
        )
        
        output$downloadEdgeSelectedGO <- downloadHandler(
          filename = function() { paste("edge_rnaseq_",input$selectedComparison,input$selectedOntology, '.csv', sep='') },
          content = function(file) {
            write.csv(network$networkSelectedOntology$edgeData, file)}
        )
        
        #output$downloadSelectedGO <- downloadHandler(
        #  filename = function() { paste0("selectedGO", '.rds') },
        #  content = function(file) {
        #    saveRDS(geneOntology$resultsOntology, file)}
        #)
        
        #WILL NEED TO BE MODIFIED, or MAKE ANOTHER FUNCTION WITH RESULTS INSTEAD OF CONTRASTS; NOW IT TAKES EVERY GENES
        output$downloadColorMatrix <- downloadHandler(
          filename = function() { paste0("colorMatrix", '.txt') },
          content = function(file) {
            write.table(
              addColorsMatrix(analysis$resultsContrast,input$selectedComparison),
              file,sep="\t\t",quote=FALSE,col.names=FALSE)}
        )
        
        
        output$downloadColorMatrixResults <- downloadHandler(
          filename = function() { paste0("colorMatrixComparison", '.txt') },
          content = function(file) {
            write.table(
              addColorsMatrixResults(analysis$results,input$selectedComparison),
              file,sep="\t\t",quote=FALSE,col.names=FALSE)}
        )
        
        output$downloadTopGO <- downloadHandler(
          filename = function() { paste0("topGO", '.rds') },
          content = function(file) {
            saveRDS(genesInfos$topGO, file)}
        )
        output$downloadModulesTable <- downloadHandler(
          filename = function() { paste0("topGO", '.rds') },
          content = function(file) {
            saveRDS(modules$modulesTable, file)}
        )
        output$downloadtopModuleGO <- downloadHandler(
          filename = function() { paste0("topGO", '.rds') },
          content = function(file) {
            saveRDS(modules$ontologyModule[[1]], file)}
        )
        output$downloadResultsModuleOntology <- downloadHandler(
          filename = function() { paste0("topGO", '.rds') },
          content = function(file) {
            saveRDS(modules$resultsModuleOntology, file)}
        )
      }
    })
    
  })
