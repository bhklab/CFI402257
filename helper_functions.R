#####################################
##### Functions to reproduce     ####
##### the analysis made in the   ####
##### following paper:           ####
##### PMID:                      ####
##### Date: Dec 27, 2017         ####
##### Wail Ba-alawi              ####
#####################################



findDrugResponseAssociationWithGeneList <- function(PSet, cellLines=NULL, drug, geneList, mDataType="rnaseq",sensitivity.measure="auc_recomputed",outputDir){
  
  
  # get molecular data 
  rnaSeq <- t(exprs(summarizeMolecularProfiles(PSet,mDataType = mDataType,fill.missing = F)))
  if(!is.null(cellLines)){
    rnaSeq <- rnaSeq[cellLines,]
  }
  # get sensitivity data 
  AUCs <- summarizeSensitivityProfiles(PSet,sensitivity.measure = sensitivity.measure,fill.missing = F,drugs = c(drug))
  
  commonCellLines <- intersect(rownames(rnaSeq),names(AUCs))
  
  rnaSeq <- rnaSeq[commonCellLines,]
  AUCs <- AUCs[commonCellLines]
  
  
  
  print("Finding univariate associations between each gene in the molecular data and the drug response")
  print("This will take some time")
  
  ## Find univariate association between each gene in the molecular data and the drug response
  drug_assoc <- PharmacoGx:::rankGeneDrugSensitivity(data = rnaSeq,drugpheno = AUCs,single.type = T,nthread = 4)
  print("Done")
  drug_assoc_drug <- as.data.frame(drug_assoc[[1]],stringsAsFactors=F)

  ibx <- match(geneList,PSet@molecularProfiles[[mDataType]]@featureData@data$Symbol)
  geneList_ENS <- rownames(PSet@molecularProfiles[[mDataType]]@featureData@data)[ibx]
  names(geneList_ENS) <- geneList
  
  geneSetInfo <- cbind("g"=geneList_ENS, "s"=rep("geneList",length(geneList_ENS)))
  
  gsc1 <- piano::loadGSC(geneSetInfo)
  
  drug_assoc_drug_NA <- drug_assoc_drug
  drug_assoc_drug_NA[is.na(drug_assoc_drug_NA[,"fdr"]),c("tstat")] <- 0
  
  genelevelstats <- drug_assoc_drug_NA[,"tstat"]
  
  names(genelevelstats) <- rownames(drug_assoc_drug_NA)
  
  print("Running enrichment analysis pipeline")
  print("This might take some time")
  set.seed(12342)
  gsea_out <- piano::runGSA(geneLevelStats=genelevelstats, geneSetStat="fgsea", gsc=gsc1,  nPerm=1000000, ncpus=4, adjMethod="none", verbose=FALSE)
  gseares <- try(piano::GSAsummaryTable(gsea_out))
  D <- length(grep("dist.dir.up",colnames(gseares)))!=0
  pdf(paste(outputDir,"Enrichment_of_geneList_in_genes_associated_with_drugResponse.pdf",sep = "/"),height = 8,width = 10)
  plotTmp <- plotEnrichment(gsc1[[1]][[1]],
                 genelevelstats) + labs(title="GeneList Enrichment Plot [Genes ranked by association with drug response]",  subtitle=paste("Estimate: " ,sprintf("%.3g", gseares[1,"Stat (dist.dir)"]), ", P-value: ",sprintf("%.1E",ifelse(D,gseares[,"p (dist.dir.up)"],gseares[,"p (dist.dir.dn)"])),sep = "" )) +
    theme( panel.grid.minor = element_blank(),plot.title = element_text(hjust = 0.5,face = "bold"), plot.subtitle  = element_text(hjust = 0.5,face = "bold") 
           ,axis.title=element_text(face="bold"),axis.text=element_text(face="bold") )
  print(plotTmp)
  dev.off()
  
  
  print(paste("Enrichment plot created at:", paste(outputDir,"Enrichment_of_geneList_in_genes_associated_with_drugResponse.pdf",sep = "/")))
  
  
  
  
  rnaSeq_geneList <- rnaSeq[,geneList_ENS]
  
  
  leadingEdgeGenes <- genelevelstats[geneList_ENS]
  leadingEdgeGenes <- sort(leadingEdgeGenes,decreasing = T)
  
  medians <- apply(rnaSeq_geneList,2,median)
  
  medians <- medians[order(medians)]
  
  ibx <- match(names(medians),geneList_ENS)
  
  
  testEnrich <- lapply(1:length(leadingEdgeGenes),function(x){
    metaGene <- apply(rnaSeq_geneList[,names(leadingEdgeGenes)[1:x],drop=F],1,mean)
    A <- cor.test(metaGene,AUCs[names(metaGene)],method = "p")
    return(c("cor"=A$estimate,"pval"=A$p.value))
  })
  
  
  corrs <- unlist(lapply(testEnrich,"[[",1))
  pvals <- unlist(lapply(testEnrich,"[[",2))
  fdr <- p.adjust(pvals,method = "fdr")
  
  
  metaGene <- apply(rnaSeq_geneList[,names(leadingEdgeGenes)[1:2],drop=F],1,mean)
  gene1 <- names(geneList_ENS)[which(geneList_ENS==names(leadingEdgeGenes)[1])]
  gene2 <- names(geneList_ENS)[which(geneList_ENS==names(leadingEdgeGenes)[2])]
  pdf(paste(outputDir,"DrugResponse_association_with_top_two_leadingEdgeGenes_in_EnrichmentAnalysis.pdf",sep = "/"))
  plot(metaGene,AUCs[names(metaGene)]*100,ylab = "Area Above the Curve (%)",xlab = "Metagene score of top two leading edge genes in enrichment analysis"
       , main = paste("Drug response association with 2-genes metagene\n[",gene1,", ",gene2,"]\ncor=",sprintf("%.3g",testEnrich[[2]][1]),", P-value=",sprintf("%.1E",testEnrich[[2]][2])))
  
  reg1 <- glm(formula =AUCs[names(metaGene)]*100~metaGene)
  abline(reg1,col="blue")
  dev.off()
  
  
  print(paste("Drug response correlation plot with 2-metaGene score created at:", paste(outputDir,"DrugResponse_association_with_top_two_leadingEdgeGenes_in_EnrichmentAnalysis.pdf",sep = "/")))
  
  
  
  metaGene <- apply(rnaSeq_geneList[,names(leadingEdgeGenes),drop=F],1,mean)
  pdf(paste(outputDir,"DrugResponse_association_with_mean_all_genes_in_geneList.pdf",sep = "/"))
  plot(metaGene,AUCs[names(metaGene)]*100,ylab = "Area Above the Curve (%)",xlab = "MetaGene Score"
       , main = paste("Drug response association with all genes in geneList \nas one metagene score [Mean]\ncor=",sprintf("%.3g",testEnrich[[17]][1]),", P-value=",sprintf("%.1E",testEnrich[[17]][2])))
  
  reg1 <- glm(formula =AUCs[names(metaGene)]*100~metaGene)
  abline(reg1,col="blue")
  dev.off()
  
  print(paste("Drug response correlation plot with all genes metaGene score created at:", paste(outputDir,"DrugResponse_association_with_mean_all_genes_in_geneList.pdf",sep = "/")))
  
  
  
  
  
  
}
