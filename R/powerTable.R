powerTable <- function(expr.matrix,names){
  names.unique=unique(names)
  print("Running WGCNA to find modules")
  nSets = length(names.unique)
  setLabels = names.unique
  comparison <- get_comparisons(names)
  # Form multi-set expression data: columns starting from 9 contain actual expression data.
  multiExpr = vector(mode = "list", length = nSets)
  
  for(i in 1:length(names.unique)){
    multiExpr[[i]] = list(data = as.data.frame(t(expr.matrix[,grep(names.unique[i],colnames(expr.matrix))])))
  }
  exprSize = checkSets(multiExpr)
  powers = c(seq(4,10,by=1), seq(12,20, by=2));
  powerTables = vector(mode = "list", length = nSets);
  # Call the network topology analysis function for each set in turn
  saveRDS(powerTables,file="powerTables.rds")
  print("Looking for soft thresholds")
  
  for(set in 1:nSets){
    print(paste0("Soft threshold #",set," / ",nSets))
    powerTables[[set]] = pickSoftThreshold(multiExpr[[set]]$data, powerVector=powers,
                                           verbose = 2)[[2]]
    
  }
  saveRDS(powerTables,paste0("powerTables_",comparison,".rds"))
}
