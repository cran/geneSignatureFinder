puttingAllTogether <-
function(workingFile = "",
         nchips = NULL,
         useCpuCluster = FALSE,
         alpha = 0.05,
         saveIndividualSignature = FALSE,
         saveSurvivalCurvesPlot = FALSE,
         saveImportancePlot = FALSE,
         saveRegulationPlot = FALSE) {
         
         
 if(areDataNotLoaded()) return(NULL)
  
  if(useCpuCluster && !is.numeric(nchips))
    stop("Provide the number of cpu's for parallel computation.")
  
  if(workingFile != "")
    message(paste("Using", workingFile, "data."))
  
  n <- nrow(geData)
  m <- ncol(geData)
  message(paste("Found", n, "samples and",m,"gene expression levels."))
  
  genesToUse <- !is.na(geData)
  genesToUse <- apply(genesToUse, 2, sum)/n
  toRemove <-  which(genesToUse <= 0.75)
  genesToUse <- which(genesToUse > 0.75)
  
  if(length(genesToUse) != m) {
    message(paste(m-length(genesToUse), "genes are going to be removed from the data matrix (more than 75% of missing values)."))
    message(paste("Removing:", paste(names(toRemove), collapse = ", ")))
    geData <- geData[, genesToUse]
    m <- ncol(geData)
  }
  
  message(paste("Working on", n, "samples and",m,"gene expression levels."))
  
  message("Starting the scanning of the genes (seeds finder step)")
  if(useCpuCluster) {
    aNCPUS <- NCPUS(nchips = nchips)
    aSeedsFinder <- seedsFinder(cpuCluster = aNCPUS, cutoff = 1.95)
    attr(aSeedsFinder, "Creation date") <- date()
  } else {
    aSeedsFinder <- seedsFinder(cutoff = 1.95)
    attr(aSeedsFinder, "Creation date") <- date()
  }  
  
  fName <- paste(workingFile, "AnsSeedsFinder.RData", sep = "")
  message("Saving the result of the scanning in ", fName, " of the working directory.")
  save(aSeedsFinder, workingFile, file = fName)
  
  bimodalGenes <- aSeedsFinder[,"bic1"] > aSeedsFinder[,"bic2"]
  if(sum(bimodalGenes, na.rm = TRUE) < 1) {
    message("No gene have been found bi-modal")
    bimodalGenes <- rep(TRUE, n)
  } else message(paste("The ", round(sum(bimodalGenes, na.rm = TRUE)/m, 4)*100,"% of the genes have been found bi-modal", sep = ""))
  
  significantGenes <- BHcorrection(aSeedsFinder[, "pValue"]) < alpha
  if(sum(significantGenes) == 0) {
    message("The B & H correction of the p-values returns no significant genes.")
    significantGenes <- aSeedsFinder[, "pValue"] < alpha
	if(sum(significantGenes) == 0) { 
	    message("No significant genes without p-value correction at level.", round(alpha, 4), "\n")
	    message("The minimum of the p-values is ", round(min(aSeedsFinder[, "pValue"]), 4), "\n")
	    message("The procedure stops.")
	    return(NULL)
	}
    message(paste("The analysis proceds on the ", round(sum(significantGenes)/m, 4)*100,"% of the genes found significant (without correction)", sep = ""))
  } else message(paste("The analysis proceds on the ", round(sum(significantGenes)/m, 4)*100,"% of the genes found significant (with correction)", sep = ""))
  
  if(sum(bimodalGenes, na.rm = TRUE) == n) {
    seedGenes <- significantGenes
    message("The analysis proceds on the significant genes (not bi-modal).")
  } else {
    seedGenes <- significantGenes * bimodalGenes
    message(paste("The ", round(sum(seedGenes, na.rm = TRUE)/m, 4)*100,"% of the genes  has been found significant and bimodal.", sep = ""))
  }
  
  seedGenes <- which(seedGenes == 1)
  seedGenesNames <- names(seedGenes)
  message(paste("The seed-genes are:", paste(seedGenesNames, collapse = ", ")))
  
  message("Starting the development of the signature from all seed-genes.")
  
  
  K <- length(seedGenes)
  searchResults <- vector("list", K)
  names(searchResults) <- seedGenesNames
  if(useCpuCluster) 
    aNCPUS <- NCPUS(nchips = nchips)
  
  for(k in 1:K) {
    message(paste("Developing from the seed:", seedGenesNames[k]))
    if(workingFile != "")
      fName <- paste(seedGenesNames[k], "@", workingFile, sep = "") else
    fName <- paste(seedGenesNames[k], sep = "") 
    if(useCpuCluster) 
      aSignatureFinder <- signatureFinder(seedGenes[k], logFilePrefix = fName,
                                          cpuCluster = aNCPUS, stopCpuCluster = FALSE) else
    aSignatureFinder <- signatureFinder(seedGenes[k], logFilePrefix = fName) 
    
    if(saveSurvivalCurvesPlot) {
      sf <- survfit(stData ~ aSignatureFinder$classification)
      plot(sf, main = paste("Survivals for signature starting from:", aSignatureFinder$startingSignature),
           xlab = paste("tValue(Log-Rank test) = ", round(aSignatureFinder$tValue, 3)),
           col = c("green", "red"))  
      dev.copy2pdf(device = x11, file = paste(fName, "SurvivalCurves.pdf", sep = ""))
    }
    
    if(length(aSignatureFinder$signature) > 1) {
      if(useCpuCluster) 
        aSignatureFinder <- importance(aSignatureFinder, cpuCluster = aNCPUS, stopCpuCluster = FALSE) else
      aSignatureFinder <- importance(aSignatureFinder) 
      if(saveImportancePlot) {
        
        barplot(aSignatureFinder$importance,
                main = "Importance based on L1GeneOut",
                sub = paste("Signature starting from:", aSignatureFinder$startingSignature))
        dev.copy2pdf(device = x11,
                     file = paste(fName, "Importance.pdf", sep = ""))
      }
      
      
      if(useCpuCluster) {
        aSignatureFinder <- testGE(aSignatureFinder, cpuCluster = aNCPUS, stopCpuCluster = FALSE)
        if(saveRegulationPlot) {
          barplot2(t(aSignatureFinder$groupMean), beside = TRUE,
                   main = paste("Signature starting from:", aSignatureFinder$startingSignature),
                   legend = paste(colnames(aSignatureFinder$groupMedian), "prognosis group"))
          dev.copy2pdf(device = x11,
                       file = paste(fName, "GenesRegulation.pdf", sep = ""))
        }
      } else message("Tests cannot be sequentially computed")
    } else message("Importance and test are not computed on this signature because of length = 1.")
    aSignatureFinder$workingFile <- workingFile
    if(saveIndividualSignature)
      save(aSignatureFinder, file = paste(fName, ".RData", sep = ""))
    searchResults[[k]] <- aSignatureFinder
    
    if(length(aSignatureFinder$signature) > 1) {
        signatureTable <- signatureSummaryTable(aSignatureFinder)
        if(!is.null(signatureTable)) {
          if(sum(search() == "package:WriteXLS") > 0) 
            WriteXLS("signatureTable", row.names = TRUE, AdjWidth = TRUE, 
                   paste(workingFile, "SummaryTableFor", searchResults[[k]]$startingSignature, ".xls", sep = "")) else 
                     print(signatureTable)
          }      
        }    
    }
  if(useCpuCluster)
    stopCluster(aNCPUS)
  fName <- paste(workingFile, "SearchResults.RData", sep = "")
  message("Saving the the signatures found in ", fName, " of the working directory.")
  save(searchResults, seedGenesNames, K, workingFile, file = fName)

  signaturesTable <- searchResultsSummaryTable(searchResults)
  ensemble <- ensembleTable(searchResults)
  if(sum(search() == "package:WriteXLS") > 0) 
    WriteXLS(c("signaturesTable", "ensemble"), row.names = TRUE, AdjWidth = TRUE, 
             paste(workingFile, "SearchResultsSummaryTable.xls", sep = "")) else {
              message("The WriteXLS library has not found in this system.\nThe tables will be forwarded to the screen.")
              message("\n\nSummary table of the signatures found.\n")
              print(signaturesTable)
	      message("\n\nSummary table of the genes in every signature.\n")
              print(ensemble)
              }

  message("The procedure stops.")
  return(NULL)
}
