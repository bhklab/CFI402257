require(Biobase)
require(PharmacoGx)
require(ggplot2)
require(gplots)
require(gridExtra)
require(piano)
require(fgsea)



# set the working directory to CFI402257-master
setwd("./CFI402257-master")



source("helper_functions.R")

######################################################
######################################################
# Prepare the data for the Breast cancer analysis

# In order to run this pipeline an object of type PSet from PharmacoGx R package is required.
# The example below downloads a PSet created from molecular data produced as part of Marcotte et.al. 2016, Cell. and drug sensitivity data that were produced in-house. 

download.file(url = "https://ndownloader.figshare.com/files/10113153?private_link=075ced04b1e6c245c484",destfile = "UHN_final.rda")

# PSet of interest
UHN_final <- readRDS("UHN_final.rda")

# gene list of interest
APC_genes <- c("ANAPC1", "ANAPC2", "ANAPC4", "ANAPC5", "ANAPC7", "ANAPC10", "ANAPC11", "ANAPC13", "ANAPC15", "ANAPC16", "CDC16", "CDC20", "CDC23", "CDC26", "CDC27", "UBE2C", "MAD2L1BP")






################################################
################################################


# function to find associations between drug response and a geneList of interest
# input is:
# PSet: an object of type PSet built using PharmacoGx R package is required.
# drug: name of drug of interest. Drug name needs to be part of the drugs in the PSet: drugNames(PSet)
# cellLines: a vector of cell line names to subset the molecular data to. CellLines names need to be part of the cell lines in the PSet: cellNames(PSet)
# geneList: a vector of genes (Symbol) of interest that represent a biological entity such as a protein complex or a phenotype such a pathway.
# outputDir: output directory of the plots. This is directory will reflect this path "./"outputDir/"
# mDataType: type of molecular data to extract from the PSet. Default is "rnaseq". For other types, please refer to PharmacoGx package documentation 
# sensitivity.measure: type of sensitivity measure to summarize the drug response curves. Default is "auc_recomputed". For other measures, please refer to PharmacoGx package documentation.

# The example below provides plots to show the association between a geneList of interest (here as an example, we use the APC complex set of genes) 
# and a drug of interest (example: paclitaxel)
findDrugResponseAssociationWithGeneList(PSet = UHN_final,drug = "paclitaxel",geneList = APC_genes,outputDir = "Plots", sensitivity.measure =  )



