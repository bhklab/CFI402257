# CFI402257

Disruption of the anaphase-promoting complex confers resistance to TTK inhibitors in triple-negative breast cancer
=================================================================


Abstract
--------


 

Citation
--------

To cite this work in publication, use



Reproducibility of the Analysis Results
--------------------------------------------

1.  Set up the software environment

2.   Download pSets

3.   Run the R scripts


Set up the software environment
-------------------------------

Analysis of this paper has been done in the following session environment. So to reproduce all the paper results the same session should be prepared by installing all the mentioned packages.

#sessionInfo()
R version 3.4.2 (2017-09-28)
Platform: x86_64-apple-darwin15.6.0 (64-bit)
Running under: OS X El Capitan 10.11.6

locale:
[1] en_CA.UTF-8/en_CA.UTF-8/en_CA.UTF-8/C/en_CA.UTF-8/en_CA.UTF-8

attached base packages:
[1] parallel  stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] fgsea_1.4.0         Rcpp_0.12.14        piano_1.18.0        gridExtra_2.3       gplots_3.0.1        ggplot2_2.2.1       PharmacoGx_1.8.2   
[8] Biobase_2.38.0      BiocGenerics_0.24.0



Download pSets
-------------------------------

These pSets are required to be accessible for the scripts in a directory named data/PSets:

UHNBreast_PSet
CCLE_Pset


Run the R scripts
-------------------------------

#pipeline_to_associate_drugResponse_with_geneList.R: 

Script aims to find associations between a given drug and geneList of interest.

*This script takes some time to be completed as it finds univariate associations between all genes and the drug response. Then, it will measure the enrichemnt of geneList in the list of all genes ranked by their univariate associations. More details of the method are provided in the paper mentioned above.

User needs to modify the script and set the working directory to "./CFI402257-master/"

The script is already set with a general example of using the UHNBreast_PSet to find associations between Paclitaxel and anaphase-promoting complex gene set.

