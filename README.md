# OutDeCo <img src="./vignettes/figures/sticker2.gif" align="right" height = 150/>
*OutDeCo*: Outlier detection through co-expression. The purpose of this package is to assess genes - more specifically differentially expressed genes - with respect to their co-expression properties. 

# Introduction 
This manual provides an overview of the Bioconductor package *OutDeCo* for differential expression outlier analysis using co-expression. 

## What do we mean by "functional outliers"? 
Genes do not act alone. They participate in pathways (genetic interactions) or form complexes (physical interactions), each of which defines their function. An approach to assess their functions is to look to differential expression, where we assay the transcriptome and search for differences between conditions. Since one way to think of gene function is to consider disease, i.e., where systems break down or respond unusually to a perturbation, differential expression is a typical approach. This is usually followed by a gene set enrichment analysis to discover the common theme or function in the set of genes. 

However, this misses an interesting and potentially important counterfactual. What if the genes are no-longer functioning in their respective roles i.e., with their common interacting partners? Rather, it might be that genes acting uncharacteristically are of relevance to the dysfunction. We call these rogue actor genes “functional outliers”. In this scenario, a gene set enrichment analysis will miss these genes. We've created this packaged to search for those outlier genes. 
![principle](./vignettes/figures/fig_outliers_network.png)



## What you can do with this package
The functions implemented in *OutDeCo* can be applied to human gene expression data or a gene list of interest. An outline of the objects and workflow included in this package are in the **schematic**, which we go into more detail in the next sections. In particular, this package allows users to:
  * Run a differential expression (DE) analysis either through the typical case-control approach or a meta-analytic recurrence analysis 
  * Assess your DE results using gene co-expression properties 
  * Report a functional outlier assessment
  * Run a network connectivity analysis of DE results within a gene co-expression network


![principle2](./vignettes/figures/fig_outliers_matrix.png)



## What you need to get started

### Data 
A gene set of interest, typically from a differential expression analysis. Or an expression experiment with cases and controls. Or multiple expression experiments or gene lists to perform a meta-analysis. These genes will be assessed with respect to their co-expression properties from a selection of aggregate co-expression networks. Optionally, you can also provide your own gene network to assess the gene set(s) against. This can be one of protein-protein interaction networks, contact networks, etc.  

### System requirements
Although not necessary, this method runs best on a HPC with 20GB+ RAM. However, datasets of a few hundred samples and up to 30,000 genes can run on smaller CPUs, without the need to parallelise.  


## What is in this user guide
This manual contains a usage guide and descriptions of the package. 


## How to get help
Ask us! Please reach out. 


## Installation
Using devtools: 
```{r , eval=FALSE}
# install.packages("devtools")
devtools::install_github("sarbal/OutDeCo")
```
The co-expression networks needed are large, so before starting you will need to download them. But this only needs to be done once. If the files already exist in your current directory, they will not be downloaded again. 
```{}
download_network_hdf5(network_type="generic") 
download_network_hdf5(network_type="blood") 
download_network_hdf5(network_type="brain") 
```

