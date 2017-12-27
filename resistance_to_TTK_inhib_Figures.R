require(Biobase)
require(PharmacoGx)
require(ggplot2)
require(gplots)
require(gridExtra)
require(piano)


# set the working directory to CFI402257-master
setwd("~/Projects/SU2C/SU2C_biomarkers_analysis/2257_story/KelsiePaper/resistance_to_TTK_inhib/CFI402257/")


######################################################
######################################################
# Prepare the data for the Breast cancer analysis

# load UHN Breast Cancer cell lines PSet: The molecular profiles are compiled from Marcotte et.al. 2016, Cell. The drug sensitivity data is measured in-house
UHNBreast <- readRDS("../UHNBreast.rda")

rnaSeq_BC <- t(exprs(summarizeMolecularProfiles(UHNBreast,mDataType = "rnaseq",fill.missing = F)))

drug <- "2257"

AUCs <- summarizeSensitivityProfiles(UHNBreast,sensitivity.measure = "AUC",fill.missing = F,drugs = c(drug))

rnaSeq_BC <- rnaSeq_BC[names(AUCs),]

#############################################################
#############################################################
# Univariate analysis


drug_assoc <- PharmacoGx:::rankGeneDrugSensitivity(data = rnaSeq_BC,drugpheno = AUCs,single.type = T,nthread = 4)

drug_assoc_drug <- as.matrix(drug_assoc[[1]])
class(drug_assoc_drug) <- "numeric"
drug_assoc_drug <- drug_assoc_drug[order(drug_assoc_drug[,"fdr"],na.last = T),]

drug_assoc_drug <- as.data.frame(drug_assoc_drug)
drug_assoc_drug <- drug_assoc_drug[
  with(drug_assoc_drug, order(drug_assoc_drug[,"fdr"], -1*drug_assoc_drug[,"estimate"])),
  ]


#############
# APC functional association with 2257

APC_genes <- c("ANAPC1", "ANAPC2", "ANAPC4", "ANAPC5", "ANAPC7", "ANAPC10", "ANAPC11", "ANAPC13", "ANAPC15", "ANAPC16", "CDC16", "CDC20", "CDC23", "CDC26", "CDC27", "UBE2C", "MAD2L1BP")

ibx <- match(APC_genes,UHNBreast@molecularProfiles$rnaseq@featureData@data$Symbol)
APC_genes_ENS <- rownames(UHNBreast@molecularProfiles$rnaseq@featureData@data)[ibx]
names(APC_genes_ENS) <- APC_genes

geneSetInfo <- cbind("g"=APC_genes_ENS, "s"=rep("APC_GeneSet",length(APC_genes_ENS)))

gsc1 <- piano::loadGSC(geneSetInfo)

drug_assoc_drug_NA <- drug_assoc_drug
drug_assoc_drug_NA[is.na(drug_assoc_drug_NA[,"fdr"]),c("tstat")] <- 0

genelevelstats1 <- drug_assoc_drug_NA[,"tstat"]

names(genelevelstats) <- rownames(drug_assoc_drug_NA)

set.seed(12342)
gsea_out <- piano::runGSA(geneLevelStats=genelevelstats, geneSetStat="fgsea", gsc=gsc1,  nPerm=1000000, ncpus=4, adjMethod="none", verbose=FALSE)
gseares <- try(piano::GSAsummaryTable(gsea_out))



library(fgsea)
pdf("Plots/Fig_5C.pdf",height = 8,width = 10)
plotEnrichment(gsc1[[1]][[1]],
               genelevelstats) + labs(title="APC/C gene set [Genes ranked by association with CFI-402257 response on Breast cancer cell lines]",  subtitle=paste("Estimate: " ,sprintf("%.3g", gseares[1,"Stat (dist.dir)"]), ", P-value: ",sprintf("%.1E",gseares[1,"p (dist.dir.up)"]),sep = "" )) +
  theme( panel.grid.minor = element_blank(),plot.title = element_text(hjust = 0.5,face = "bold"), plot.subtitle  = element_text(hjust = 0.5,face = "bold") 
         ,axis.title=element_text(face="bold"),axis.text=element_text(face="bold") )
dev.off()

#############################################################
#############################################################
# Univariate analysis by subtype

rnaSeq_BC_subtypes <- UHNBreast@cell[rownames(rnaSeq_BC), "SCMOD2"]
rnaSeq_BC_subtypes <- unlist(lapply(strsplit(rnaSeq_BC_subtypes," "),"[[",1))
names(rnaSeq_BC_subtypes) <- rownames(rnaSeq_BC)

subtype <- "Basal"
#subtype <- "HER2"
#subtype <- "LumB"

ibxSubtypes <- names(rnaSeq_BC_subtypes)[which(rnaSeq_BC_subtypes==subtype)]

drug_assoc_subtype <- PharmacoGx:::rankGeneDrugSensitivity(data = rnaSeq_BC[ibxSubtypes,],drugpheno = AUCs[ibxSubtypes],single.type = T,nthread = 4)
 
drug_assoc_drug_subtype <- as.matrix(drug_assoc_subtype[[1]])
class(drug_assoc_drug_subtype) <- "numeric"
drug_assoc_drug_subtype <- drug_assoc_drug_subtype[order(drug_assoc_drug_subtype[,"fdr"],na.last = T),]

drug_assoc_drug_subtype <- as.data.frame(drug_assoc_drug_subtype)
drug_assoc_drug_subtype <- drug_assoc_drug_subtype[
 with(drug_assoc_drug_subtype, order(drug_assoc_drug_subtype[,"fdr"], -1*drug_assoc_drug_subtype[,"estimate"])),
]


#############
# APC functional association with 2257 with subtype


drug_assoc_drug_subtype[is.na(drug_assoc_drug_subtype[,"fdr"]),c("tstat")] <- 0
genelevelstats_subtype <- drug_assoc_drug_subtype[,"tstat"]
names(genelevelstats_subtype) <- rownames(drug_assoc_drug_subtype)

gsea_out_subtype <- piano::runGSA(geneLevelStats=genelevelstats_subtype, geneSetStat="fgsea", gsc=gsc1,  nPerm=1000000, ncpus=4, adjMethod="none", verbose=FALSE)
gseares_subtype <- try(piano::GSAsummaryTable(gsea_out_subtype))



library(fgsea)
pdf(paste("Plots/Fig_5D.pdf",sep = ""),height = 8,width = 10)
plotEnrichment(gsc1[[1]][[1]],
               genelevelstats_subtype) + labs(title=paste("APC/C gene set [Genes ranked by association with CFI-402257 response on ",subtype," BC cell lines]",sep=""),  subtitle=paste("Estimate: " ,sprintf("%.3g", gseares_subtype[1,"Stat (dist.dir)"]), ", P-value: ",sprintf("%.1E",ifelse(!is.null(gseares_subtype[1,"p (dist.dir.up)"]),gseares_subtype[1,"p (dist.dir.up)"],gseares_subtype[1,"p (dist.dir.dn)"])),sep = "" )) +
  theme( panel.grid.minor = element_blank(),plot.title = element_text(hjust = 0.5,face = "bold"), plot.subtitle  = element_text(hjust = 0.5,face = "bold") 
         ,axis.title=element_text(face="bold"),axis.text=element_text(face="bold") )
dev.off()


#########################################
#########################################
#########################################
# APC distribution cell lines

rnaSeq_BC_APC <- rnaSeq_BC[,APC_genes_ENS]


leadingEdgeGenes <- genelevelstats[APC_genes_ENS]
leadingEdgeGenes <- sort(leadingEdgeGenes,decreasing = T)

medians <- apply(rnaSeq_BC_APC,2,median)

medians <- medians[order(medians)]

ibx <- match(names(medians),APC_genes_ENS)


testEnrich <- lapply(1:length(leadingEdgeGenes),function(x){
  metaGene <- apply(rnaSeq_BC_APC[,names(leadingEdgeGenes)[1:x],drop=F],1,mean)
  A <- cor.test(metaGene,AUCs[names(metaGene)],method = "p")
  return(c("cor"=A$estimate,"pval"=A$p.value))
})


corrs <- unlist(lapply(testEnrich,"[[",1))
pvals <- unlist(lapply(testEnrich,"[[",2))
fdr <- p.adjust(pvals,method = "fdr")


metaGene <- apply(rnaSeq_BC_APC[,names(leadingEdgeGenes)[1:2],drop=F],1,mean)

pdf("Plots/Fig_5E.pdf")
plot(metaGene,AUCs[names(metaGene)]*100,ylab = "Area Above the Curve (%)",xlab = "ANAPC4/CDC20 Metagene Score"
     , main = paste("CFI-402257 response association with 2-gene metagene\ncor=",sprintf("%.3g",testEnrich[[2]][1]),", P-value=",sprintf("%.1E",testEnrich[[2]][2])))

reg1 <- glm(formula =AUCs[names(metaGene)]*100~metaGene)
abline(reg1,col="blue")
dev.off()

metaGene <- apply(rnaSeq_BC_APC[,names(leadingEdgeGenes),drop=F],1,mean)
pdf("Plots/Fig_S5B.pdf")
plot(metaGene,AUCs[names(metaGene)]*100,ylab = "Area Above the Curve (%)",xlab = "MetaGene Score"
     , main = paste("CFI-402257 response association with 17-gene metagene\ncor=",sprintf("%.3g",testEnrich[[17]][1]),", P-value=",sprintf("%.1E",testEnrich[[17]][2])))

reg1 <- glm(formula =AUCs[names(metaGene)]*100~metaGene)
abline(reg1,col="blue")
dev.off()



#########################################
#########################################
#########################################


#TCGA


samplesInfo <- read.csv("DATA/PANCAN_clinicalMatrix", header=T, row.names = 1, sep="\t",stringsAsFactors = F)

TCGA_tissues <- table(samplesInfo$X_primary_site)

TCGA_tissues_sub <-names(TCGA_tissues)[which(TCGA_tissues>500)]


TCGA_expession_APC_C <- readRDS("DATA/TCGA_expession_APC_C.rda")

medians <- apply(TCGA_expession_APC_C,2,median)

medians <- medians[order(medians)]

ibx <- match(names(medians),APC_genes_ENS)

MetaGene_TwoLeading_mean <- apply(TCGA_expession_APC_C[,names(leadingEdgeGenes)[1:2]],1,mean)



allSamplesIDs <- intersect(rownames(TCGA_expession_APC_C),rownames(samplesInfo))
TCGA_expression.df <- cbind("Tissue"=samplesInfo[allSamplesIDs,"X_primary_site"],"MetaGene_TwoLeading_mean"=MetaGene_TwoLeading_mean[allSamplesIDs])

TCGA_expression.df <- as.data.frame(TCGA_expression.df,stringsAsFactors = F)

TCGA_expression.df_sub <- TCGA_expression.df[TCGA_expression.df$Tissue %in% TCGA_tissues_sub, ]

TCGA_expression.df_sub$Tissue <- as.factor(TCGA_expression.df_sub$Tissue)

class(TCGA_expression.df_sub$MetaGene_TwoLeading_mean) <- "numeric"

pdf("Plots/Fig_5F.pdf",width = 13,height = 8)
ggplot(TCGA_expression.df_sub, aes(x=Tissue, y=MetaGene_TwoLeading_mean)) + 
  geom_violin(trim=F,fill="gray") + ggtitle("ANAPC4/CDC20 Metagene Score Across Various TCGA Tumour Types") +  
  labs(y = "ANAPC4/CDC20 Metagene Score") +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5),panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_boxplot(width=0.03,fill="black",outlier.shape = NA) + 
  stat_summary(geom = "crossbar", width=0.05, fatten=3, color="white",
              fun.data = function(x){c(y=median(x), ymin=median(x), ymax=median(x))}) 
dev.off()


###########################################################
###########################################################
# Lung analysis


load("DATA/lungData.RData")


drug_assoc_Lung <- PharmacoGx:::rankGeneDrugSensitivity(data = rnaSeq_lung,drugpheno = AUCs_Lung[rownames(rnaSeq_lung),"AUC",drop=F],single.type = T,nthread = 4)

drug_assoc_drug_Lung <- as.matrix(drug_assoc_Lung[[1]])
class(drug_assoc_drug_Lung) <- "numeric"
drug_assoc_drug_Lung <- drug_assoc_drug_Lung[order(drug_assoc_drug_Lung[,"fdr"],na.last = T),]

drug_assoc_drug_Lung <- as.data.frame(drug_assoc_drug_Lung)
drug_assoc_drug_Lung <- drug_assoc_drug_Lung[
  with(drug_assoc_drug_Lung, order(drug_assoc_drug_Lung[,"fdr"], -1*drug_assoc_drug_Lung[,"estimate"])),
  ]

drug_assoc_drug_Lung[is.na(drug_assoc_drug_Lung[,"fdr"]),c("tstat")] <- 0
drug_assoc_drug_Lung_NA <- drug_assoc_drug_Lung

genelevelstats_Lung <- drug_assoc_drug_Lung_NA[,"tstat"]
names(genelevelstats_Lung) <- rownames(drug_assoc_drug_Lung_NA)

gsea_out_Lung <- piano::runGSA(geneLevelStats=genelevelstats_Lung, geneSetStat="fgsea", gsc=gsc1,  nPerm=1000000, ncpus=4, adjMethod="none", verbose=FALSE)
gseares_Lung <- try(piano::GSAsummaryTable(gsea_out_Lung))


library(fgsea)
pdf("Plots/Fig_S5A.pdf",height = 8,width = 10)
plotEnrichment(gsc1[[1]][[1]],
               genelevelstats) + labs(title="APC/C gene set [Genes ranked by association with CFI-402257 response on Lung Cancer cell lines]",  subtitle=paste("Enrichment Score: " ,sprintf("%.3g", gseares_Lung[1,"Stat (dist.dir)"]), ", P-value: ",sprintf("%.1E",gseares_Lung[1,"p (dist.dir.up)"]),sep = "" )) +
  theme( panel.grid.minor = element_blank(),plot.title = element_text(hjust = 0.5,face = "bold"), plot.subtitle  = element_text(hjust = 0.5,face = "bold") 
         ,axis.title=element_text(face="bold"),axis.text=element_text(face="bold") )
dev.off()







