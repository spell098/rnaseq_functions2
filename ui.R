
library(shiny)
shinyUI(fluidPage(
  titlePanel("RNAseq analysis"),
  
  tabsetPanel(
    tabPanel("RNA-seq",
             sidebarLayout(
               sidebarPanel(
                 helpText("Options"),
                 fileInput("File1", "Updload RNAseq data file"),
                 selectInput("fileContent","File Content",choices=c("Expression Matrix","Expression Set"), selected="Expression Matrix"),
                 tags$textarea(id = "geo_id", placeholder = 'GEO ID', rows = 1, "GSE54839"),
                 fileInput("File2", "Optional: Upload RNAseq module file"),
                 fileInput("File3", "Optional: Upload RNAseq reactions list"),
                 radioButtons("randomSamples",label=h5("Use a sample of whole dataset"),choices=list("True" = "TRUE","False"="FALSE"),selected=FALSE),
                 numericInput("randomSampling","Select number of random elements to use",value=1000),
                 actionButton("submitRun", "Load"),
                 selectInput("standardization","Standardization method",choices = c("no","none", "scale", "quantile", "Aquantile", "Gquantile", "Rquantile", "Tquantile", "vsn", "cyclicloess"), selected="no"),
                 selectInput("specie","Specie",choices=c("hsapiens","rnorvegicus")),
                 selectInput("mart","Mart type",choices=c("ENSEMBL_MART_ENSEMBL","ENSEMBL_MART_SNP","ENSEMBL_MART_FUNCGEN","ENSEMBL_MART_VEGA","pride")),
                 radioButtons("type",label=h5("transcript or gene results?"),choices=list("Genes" = "ensembl_gene_id","Transcripts"="ensembl_transcript_id")),
                 tags$textarea(id = "attribute", placeholder = 'Type of ID e.g. affy_hg_u95av2', rows = 3, "illumina_humanht_12_v3"),
                 selectInput("groupBy",label="Make group by",choices=c("groups","fuzzy"),selected="groups"),
                 numericInput("noCluster",label= "Number of clusters for fuzzy c-means clustering",value=3),
                 selectInput("selectedVariables",label= "Select the variables to compare",choices=NULL,multiple=TRUE),
                 selectInput("selectedReplicates", label = "Replicates informations", choices = NULL,multiple=TRUE),
                 selectInput("selectedComparison", label = "Select the comparison of interest", choices = NULL),
                 actionButton("submitParam","Update dataset parameters"),
                 actionButton("submitModule","Calculate Modules"),
                 numericInput("pow","Power for WGCNA",value=6),
                 helpText("This action may take a while..."),
                 numericInput("pvalue", label = "Select a p-value",value = 0.01,max=1,min=0),
                 numericInput("logFCneg", label = "Select the negative logFC limit", value = -1.3),
                 numericInput("logFCpos", label = "Select the positive logFC limit", value = 1.3),
                 selectInput("adjust",label="p-value correction method",selected="BH",choices=c("no","none","BH","BY","holm")),
                 selectInput("ordination","Select the type of ordination",choices=c("Principal Component Analysis","Redundency analysis","correspondence analysis","detrended correspondence analysis","non-metric multidimensional scaling")),
                 actionButton("submit2", "Run!"),
                 #actionButton("findSignificantModules","Find Significant Modules"),
                 actionButton("submitNetwork", "Network Analysis")
                 
                 ,witdth=2),
               
               mainPanel(
                 textOutput("warnings"),
                 textOutput("text1"),
                 tabsetPanel(
                   tabPanel("Quality Control",
                            tabsetPanel(
                              tabPanel("Boxplot",
                                       sliderInput("cex.axis_boxplot",label = "Caracters size on the axis",min=0,max=1,value=0.6),
                                       sliderInput("ylim",label = "Limits of the y-axis", min = -30, max = 30, value = c(-15, 15),step=0.1),
                                       plotOutput("boxplot"),
                                       help('Save data with new groups'),
                                       downloadButton("downloadExprMatrix","Download data")
                              ),
                              tabPanel("Comparisons Table",
                                       DT::dataTableOutput("comparisons_table")
                              ),
                              tabPanel("Clustering",
                                       actionButton("submitSummaryClust","Generate clustering summary tables"),
                                       numericInput("selectedClustn",label = "select number of clusters for summary",value=3),
                                       DT::dataTableOutput("summaryClustering"),
                                       numericInput("pval_clust",label="select a p-value for significant clustering group",value=0.2),
                                       plotOutput("plotCorrectAttribution"),
                                       plotOutput("plotCorrectSubsample")
                              )
                            )
                            #N'arrive pas a afficher tous les graphique en loop
                            
                   ),
                   tabPanel("Ordination",
                            tabsetPanel(
                              tabPanel("Ordination",
                                       radioButtons("random",label=h5("Use a sample of whole dataset"),choices=list("True" = "TRUE","False"="FALSE"),selected=FALSE),
                                       numericInput("randomSamplingOrd","Select number of random elements to use",value=1000),
                                       actionButton("submitOrdination","Generate ordination analysis"),
                                       sliderInput("cex.PCA", label = "Labels size of the plot", min=0.5,max = 2, value = 1),
                                       plotOutput("ordination",height="500px"),
                                       sliderInput("cex.MDS", label = "Labels size of the MDS plot", min=0.5,max = 2, value = 1),
                                       plotOutput("MDS")
                                       
                              ),
                              tabPanel("Principal Components top genes",
                                       DT::dataTableOutput("eigenvalues"),
                                       DT::dataTableOutput("u"),
                                       DT::dataTableOutput("v"),
                                       DT::dataTableOutput("ordinationAnova"),
                                       plotOutput("boxplot_element_pca"),
                                       numericInput("bins","Nunmber of columns",value=20),
                                       plotOutput("hist_element_pca"),
                                       selectInput("topOrdinationGenes","Principal Componant of interest",choices=NULL)
                              ),
                              tabPanel("Heatmaps",
                                       selectInput("dendrogram", label = "Dendrograms",choices = c("column","row","both","none"),selected = "none"),
                                       selectInput("scale", label = "Scale", choices = c("column","row"), selected="row"),
                                       sliderInput("cexRow_heatmap_ord",label = "Caracters size on the x-axis",min=0,max=1.5,value=0.5),
                                       sliderInput("row_heatmap_height_ord",label = "Row height",min=0,max=6,value=2,step=0.05),
                                       tabsetPanel(
                                         tabPanel("Individual samples",
                                                  plotOutput("heatmap1_ord",width="100%",height="1000px")
                                         ),
                                         tabPanel("Grouped samples",
                                                  plotOutput("heatmap2_ord",width="100%",height="1000px")
                                         )
                                       )
                                       
                              ),
                              tabPanel("Top gene Ontologies",
                                       DT::dataTableOutput("topGO_ord"),
                                       downloadButton("downloadTopGO_ord","Download"),
                                       tags$textarea(id = "selectedOntology_ord", placeholder = 'ontologies of interest', rows = 8, ""),
                                       
                                       DT::dataTableOutput("selectedGO_ord"),
                                       downloadButton('downloadSelectedGO_ord', 'Download'),
                                       downloadButton('downloadNodeSelectedGO_ord','Download nodes file for cytoscape'),
                                       downloadButton('downloadEdgeSelectedGO_ord','Download edges file for cytoscape')
                              ),
                              tabPanel("Module analysis",
                                       DT::dataTableOutput("modulesTableOrd"),
                                       #downloadButton("downloadModulesTable","Download"),
                                       tabsetPanel(
                                         tabPanel("Genes of interest",
                                                  DT::dataTableOutput("selectedModule3"),
                                                  #downloadButton("downloadSelectedModulePCA","Download"),
                                                  DT::dataTableOutput("topModuleGO3"),
                                                  #downloadButton("downloadtopPCAModuleGO","Download"),
                                                  DT::dataTableOutput("resultsModuleOntology3")
                                                  #downloadButton("downloadResultsModulePCAOntology","Download")
                                         ),
                                         tabPanel("Module analysis",
                                                  DT::dataTableOutput("selectedModule4"),
                                                  #downloadButton("downloadModuleTable","Download"),
                                                  DT::dataTableOutput("topModuleGO4"),
                                                  #downloadButton("downloadtopModuleGO","Download"),
                                                  DT::dataTableOutput("resultsModuleOntology4")
                                                  #downloadButton("downloadResultsModuleTableOntology","Download")
                                         )
                                       )
                              )
                            )
                   ),
                   
                   tabPanel("Differential analysis results",
                            tabsetPanel(
                              tabPanel("Volcano & MA plots",
                                       numericInput("pvalueVolcano",label = "p-val for volcano plot",value = 0.05),
                                       numericInput("logFCVolcano",label = "Absolute logFC for volcano plot", value = 1),
                                       plotOutput("volcanoplot"),
                                       plotOutput("MA")
                                       
                              ),
                              tabPanel("Heatmaps",
                                       selectInput("dendrogram", label = "Dendrograms",choices = c("column","row","both","none"),selected = "both"),
                                       selectInput("scale", label = "Scale", choices = c("column","row"), selected="row"),
                                       sliderInput("cexRow_heatmap",label = "Caracters size on the x-axis",min=0,max=1.5,value=0.5),
                                       sliderInput("row_heatmap_height",label = "Row height",min=0,max=6,value=2,step=0.05),
                                       tabsetPanel(
                                         tabPanel("Individual samples",
                                                  plotOutput("heatmap1",width="100%",height="1000px")
                                         ),
                                         tabPanel("Grouped samples",
                                                  plotOutput("heatmap2",width="100%",height="1000px")
                                         ),
                                         tabPanel("Comparison heatmap",
                                                  sliderInput("heatdiagram_cex",label = "Cex heatdiagram",min=0,max=3,value=0.6,step=0.05),
                                                  plotOutput("heatDiagram",width="100%",height="1000px")
                                         )
                                       )
                                       
                              ) ,
                              tabPanel("Results summary",
                                       help('Download results all comparisons'),
                                       downloadButton('downloadAllResults', 'Download'),
                                       DT::dataTableOutput("resultsSummary"),
                                       downloadButton("downloadSummary","Download"),
                                       plotOutput("boxplot_element"),
                                       plotOutput("hist_element")
                              ),
                              #Mettre ensembl results table, gene ontology et network. Faire apparaitre graph boxplot / chaque table
                              tabPanel("Results table",
                                       DT::dataTableOutput("topTable"),
                                       downloadButton('downloadResultsSelectedComparison', 'Download table'),
                                       downloadButton('downloadNodeDataComparison','Download nodes file for cytoscape'),
                                       downloadButton('downloadEdgeDataComparison','Download edges file for cytoscape'),
                                       downloadButton("downloadColorMatrixResults","Download")
                              ),
                              tabPanel("Venn",
                                       sliderInput("lfc", label = "Select log Fold Change (logFC) limit (absolute)", value = c(1) ,max=5,min=0,step=0.05),
                                       sliderInput("venn.cex", label = "Venn Diagram caracters size", max=2,min=0.5,value=1),
                                       plotOutput("venn"),
                                       downloadButton("downloadVenn","Download")
                              ),
                              tabPanel("Top gene Ontologies",
                                       DT::dataTableOutput("topGO"),
                                       downloadButton("downloadTopGO","Download"),
                                       tags$textarea(id = "selectedOntology", placeholder = 'ontologies of interest', rows = 8, ""),
                                       
                                       DT::dataTableOutput("selectedGO"),
                                       downloadButton('downloadSelectedGO', 'Download'),
                                       downloadButton('downloadNodeSelectedGO','Download nodes file for cytoscape'),
                                       downloadButton('downloadEdgeSelectedGO','Download edges file for cytoscape')
                              ),
                              tabPanel("Top modules",
                                       DT::dataTableOutput("modulesTable"),
                                       #downloadButton("downloadModulesTable","Download"),
                                       tabsetPanel(
                                         tabPanel("Genes of interest",
                                                  DT::dataTableOutput("selectedModule1"),
                                                  #downloadButton("downloadSelectedModule","Download"),
                                                  DT::dataTableOutput("topModuleGO1"),
                                                  #downloadButton("downloadtopModuleGO","Download"),
                                                  DT::dataTableOutput("resultsModuleOntology1")
                                                  #downloadButton("downloadResultsModuleOntology","Download") 
                                         ),
                                         tabPanel("Module analysis",
                                                  DT::dataTableOutput("selectedModule2"),
                                                  #downloadButton("downloadModuleTable","Download"),
                                                  DT::dataTableOutput("topModuleGO2"),
                                                  #downloadButton("downloadtopModuleGO","Download"),
                                                  DT::dataTableOutput("resultsModuleOntology2")
                                                  #downloadButton("downloadResultsModuleTableOntology","Download")
                                         )
                                       )
                              ),
                              tabPanel("top ordination genes",
                                       DT::dataTableOutput("diffOrd"),
                                       DT::dataTableOutput("resultsOrdDiff"),
                                       DT::dataTableOutput("topOrdDiffGO"),
                                       DT::dataTableOutput("resultsOrdDiffOntology")
                              ),
                              tabPanel("Raw data",
                                       #don't display the whole thing, it's useless; make the user enter is genes of interest manually, and display only that
                                       tags$textarea(id = "selectedGenes", placeholder = 'genes of interest', rows = 8, ""),
                                       DT::dataTableOutput("rawTopTable"),
                                       downloadButton("downloadRawTopTable","Download"),
                                       plotOutput("boxplot_element_raw")
                              )
                              
                            )

                   ),
                   tabPanel("Network Analysis",
                            textOutput("confirmFindReactions"),
                            numericInput("threshold","Adjacency threshold for network",value = 0.7),
                            actionButton("findTopReactions","Find top reactions"),
                            DT::dataTableOutput("topReactions"),
                            tags$textarea(id = "selectedReactions", placeholder = 'Reactions of interest', rows = 3, ""),
                            
                            plotOutput("pathway",width="100%",height="800px"),
                            DT::dataTableOutput("pathwayTable")
                   )
                 ))
             ) ),
    tabPanel("Methylation",
             sidebarLayout(
               sidebarPanel(
                 helpText("Option"),
                 fileInput("File5", "Updload metylation table file"),
                 fileInput("File6", "Updload metylation results file"),
                 fileInput("File7", "Optional: Upload methylation module file"),
                 fileInput("File8", "Optional: Upload methylation reactions list"),
                 
                 selectInput("specieEnsembl","specieEnsembl",choices=c("rnorvegicus_gene_ensembl","hsapiens_gene_ensembl")),
                 selectInput("symbol","Type of gene symbols",choices=c("rgd_symbol","hgnc_symbol")),
                 
                 selectInput("selectedComparisonMethyl", label = "Select the comparison of interest", choices = NULL),
                 selectInput("selectedVariablesMethyl",label= "Select the variables to compare",choices=NULL,multiple=TRUE),
                 
                 actionButton("submit21", "Run!"),
                 numericInput("methDiffNeg", label = "Select negative limit",value = -25),
                 numericInput("methDiffPos", label = "Select positive limit",value = 25),
                 
                 numericInput("qvalue", label = "Select a q-value",value = 0.01,max=1,min=0),
                 
                 actionButton("submit22", "Update restriction parameters")
                 
                 ,witdth=2),
               
               mainPanel(
                 textOutput("test2"),
                 textOutput("test3"),
                 
                 actionButton("submit23","Calculate differential analysis!"),
                 tabsetPanel(
                   tabPanel("Quality Control",
                            tabsetPanel(
                              tabPanel("Boxplot",
                                       sliderInput("cex.axis_boxplot2",label = "Caracters size on the axis",min=0,max=1,value=0.6),
                                       sliderInput("ylim2",label = "Limits of the y-axis", min = 0, max = 1, value = c(0, 1),step=0.1),
                                       plotOutput("boxplot2")
                              ),
                              tabPanel("MDS",
                                       sliderInput("cex.MDS", label = "Labels size of the MDS plot", min=0.5,max = 2, value = 1),
                                       plotOutput("MDS2"),
                                       plotOutput("PCA2")
                              ),
                              tabPanel("Volcanoplots",
                                       plotOutput("volcanoplot2")
                              ),
                              tabPanel("Correlation",
                                       plotOutput("correlation")
                              ),
                              tabPanel("Tree",
                                       plotOutput("tree")
                              )
                            )
                            #N'arrive pas a afficher tous les graphique en loop
                            
                   ),
                   
                   tabPanel("Heatmaps",
                            selectInput("dendrogram2", label = "Dendrograms",choices = c("column","row","both","none"),selected = "none"),
                            selectInput("scale2", label = "Scale", choices = c("column","row"), selected="row"),
                            sliderInput("cexRow_heatmap2",label = "Caracters size on the x-axis",min=0,max=1.5,value=0.5),
                            sliderInput("row_heatmap_height2",label = "Row height",min=0,max=6,value=2,step=0.05),
                            tabsetPanel(
                              tabPanel("Individual samples",
                                       plotOutput("heatmap21",width="100%",height="1000px")
                              ),
                              tabPanel("Grouped samples",
                                       plotOutput("heatmap22",width="100%",height="1000px")
                              ),
                              tabPanel("Comparison heatmap",
                                       sliderInput("heatdiagram_cex",label = "Cex heatdiagram",min=0,max=3,value=0.6,step=0.05),
                                       plotOutput("heatDiagram2",width="100%",height="1000px")
                              )
                            )
                   ) ,
                   tabPanel("Results",
                            tabsetPanel(
                              tabPanel("Results table",
                                       DT::dataTableOutput("topTableMethyl"),
                                       downloadButton('downloadResultsMethyl', 'Download'),
                                       downloadButton('downloadNodeDataComparisonMethyl','Download nodes file for cytoscape'),
                                       downloadButton('downloadEdgeDataComparisonMethyl','Download edges file for cytoscape'),
                                       DT::dataTableOutput("topTableMethylOutGene"),
                                       downloadButton('downloadResultsMethylOut', 'Download'),
                                       downloadButton('downloadNodeDataComparisonMethylOutGene','Download nodes file for cytoscape'),
                                       downloadButton('downloadEdgeDataComparisonMethylOutGene','Download edges file for cytoscape'),
                                       
                                       DT::dataTableOutput("cpgByGene"),
                                       downloadButton('downloadResultsByGene', 'Download'),
                                       downloadButton('downloadNodeDataComparisonMethylByGene','Download nodes file for cytoscape'),
                                       downloadButton('downloadEdgeDataComparisonMethylByGene','Download edges file for cytoscape')
                                       
                                       
                              ),
                              tabPanel("Venn",
                                       sliderInput("lfcMethyl", label = "Select log Fold Change (logFC) limit (absolute)", value = c(1) ,max=5,min=0,step=0.05),
                                       sliderInput("venn.cexMethyl", label = "Venn Diagram caracters size", max=2,min=0.5,value=1),
                                       plotOutput("vennMethyl")
                              ),
                              tabPanel("Top gene Ontologies",
                                       DT::dataTableOutput("topGOMethyl")
                              ),
                              tabPanel("Selected ontologies results",
                                       numericInput("headGOMethyl", label = "Number of gene ontologies to display", value = 10),
                                       tags$textarea(id = "selectedOntology", placeholder = 'ontologies of interest', rows = 8, ""),
                                       #selectInput("selectedOntology",label = "Select ontologies of interest",choices = NULL),
                                       numericInput("thresholdMethyl","Adjacency threshold for network",value = 0.7),
                                       
                                       DT::dataTableOutput("selectedGOMethyl"),
                                       downloadButton('downloadSelectedGOMethyl', 'Download'),
                                       downloadButton('downloadNodeSelectedGOMethyl','Download nodes file for cytoscape'),
                                       downloadButton('downloadEdgeSelectedGOMethyl','Download edges file for cytoscape')
                              ),
                              tabPanel("Top modules",
                                       actionButton("submitMethylModule","Calculate Modules"),
                                       helpText("This action may take a while..."),
                                       DT::dataTableOutput("moduleTableMethyl")
                                       
                              ),
                              tabPanel("Selected modules results",
                                       numericInput("headModuleMethyl", label = "Number of Modules to display", value = 20),
                                       selectInput("selectedModuleMethyl", label = "Select the module of interest", choices = NULL, multiple=TRUE),
                                       DT::dataTableOutput("selectedModuleMethyl")
                                       
                              ),
                              tabPanel("Raw CpGs",
                                       tags$textarea(id = "selectedCpgs", placeholder = 'CpGs of interest', rows = 8, ""),
                                       DT::dataTableOutput("rawTopCpG")
                              ),
                              tabPanel("Raw methylated genes",
                                       tags$textarea(id = "selectedCpgGenes", placeholder = 'Genes of interest', rows = 8, ""),
                                       DT::dataTableOutput("rawTopCpgGenes")
                              )
                            )
                   ),
                   tabPanel("Network Analysis",
                            textOutput("confirmFindReactions2"),
                            actionButton("findTopReactions2","Find top reactions"),
                            DT::dataTableOutput("topReactions2"),
                            tags$textarea(id = "selectedReactions2", placeholder = 'Reactions of interest', rows = 3, ""),
                            
                            #selectInput("selectedReactions", label = "Select pathways of interest", choices = NULL),
                            plotOutput("pathway2",width="100%",height="800px"),
                            DT::dataTableOutput("toptableNodes2")
                   )
                 ))
             )
    ),
    tabPanel("RNAseq / Methylation integration",
             tabsetPanel(
               
               tabPanel("CpG/RNA",
                        DT::dataTableOutput("mergeTable")
               ),
               tabPanel("Any CpG/RNA"
               ),
               tabPanel("Network"
               )
             )
    )
  )))
