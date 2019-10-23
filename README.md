# Genome-wide-mediation-analysis

This repository aims at sharing the scripts and methods used in the paper "Using genetics to disentangle the complex relationship between food choices and health status"

##Installation 

Genome-wide mediation analysis relies heavily on the bGWAs package (https://github.com/n-mounier/bGWAS).
The package has been updated recently and thus we have included the version used in the paper (0.3) in this repository.
The bGWAS package has the following requirements which need to be installed before the package:

Depends:

* R (>= 3.3.0)

Packages:

*  calibrate (>= 1.7.2),
*  data.table (>= 1.10.4.3),
*  fastmatch (>= 1.1.0),
*  ggplot2 (>= 2.2.1),
*  gplots (>= 3.0.1),
*  gtools (>= 3.5.0),
*  magrittr (>= 1.5),
*  qqman (>= 0.1.4),
*  Rcpp (>= 0.12.15),
*  readr (>= 1.1.1),
*  TwoSampleMR (>= 0.3.0)

Also it the Z-matrices files will be needed. 
From inside the data folder type
```
wget https://www.dropbox.com/s/bho8irpf2myavzh/ZMatrices_2018-11-12.tar.gz?dl=1
tar -zxvf ZMatrices_2018-11-12.tar.gz

```

To install the bGWAS package from the terminal in the repository folder:

```
R CMD INSTALL bGWAS_package
```

Once the bGWAS package is installed we need to first load the functions which are contained in the scripts folder.
From R:
```
library(bGWAS)
library(data.table)

source("scripts/Functions_MR_correction.R")

```

## Running a GWMA analysis

First load the example data containing the results from the GWAS of Cheese consumption
```
mygwas=fread("data/Cheese_gwas_results.tsv.gz",data.table=F)

```
The GWAS has the following columns:

*rsid: SNP id 
*a1 : coded allele
*a0 : non-coded allele
*beta1 : effect of the coded allele
*se : standard error of the effect estimate

We are now ready to run the bGWAS run.
To limit the analysis only to the traits used in the paper we first select the traits to be included
```
trait.f=c("cardiogram_gwas_results.txt",
          "ibd_cd.txt",
          "TransEthnic_T2D_GWAS.MegaMeta.2014OCT16.zip",
          "ibd_uc.txt",
          "EDUyears_2016_sumstat.txt",
          "DBP",
          "SBP",
          "All_ancestries_SNP_gwas_mc_merge_nogc.tbl.uniq.gz",
          "jointGwasMc_HDL.txt.gz",
          "jointGwasMc_LDL.txt.gz",
          "jointGwasMc_TG.txt.gz")

mystudies=select_priorGWASs(Z_matrices = "data/ZMatrices_Nicola_2018-11-12",
                            include_files =trait.f)

```

We then proceed to run the bGWAS analysis which will create the prior expected effect.
```
my_bgwas=  bGWAS(name = "Cheese",
               GWAS = res.t,
               Z_matrices = "data/ZMatrices_2018-11-12",
               prior_shrinkage = 1,
               prior_studies = mystudies,
               MR_threshold = 5e-08,
               MR_pruning_dist = 1000,
               MR_pruning_LD = 0.1,
               MR_shrinkage = 1,
               stop_after_prior = T,
               save_files = FALSE,
               verbose=T)

```

The resulting output can be now used for the GWMA analysis

```

corrected.gwas=gwma.data.frame(bGWAS.obj=my_bgwas,         # bGWAS output object
                         results.data.frame=mygwas)        # data.frame of results


```
The output will have the following columns:

* rs : SNP id
* chrm : Chromosome          
* pos  : Position in base pairs    
* alt  : Coded allele
* ref  : Non-coded allele 
* observed_Z : Z-score for the observed result, comes directly from the original gwas
* prior_estimate : prior expected effect in the z-score scale
* prior_std_error : Standard error of the prior expected effect
* posterior_estimate : Posterior estimate (see original bGWAs for details non needed for GWMA ) 
* posterior_std_error : Standard error of posterior estimate (see original bGWAs for details non needed for GWMA ) 
* beta1 : Original effect estimate from GWAS              
* se : Original standard error of the effect estimate from GWAS               
* expected_p : p-value of the prior expected effect         
* z_diff : Z-score difference between the original GWAS and the prior expected effect              
* corr.beta : Corrected effect
* corr.p : p-value of the corrected effect
* p : p-value of the original GWAS               
* corr2raw_ratio : Corrected to Raw ratio value.

It is now possible to plot the miami plot to compare the corrected results to the original ones:
```
gwma.miami(gwma.obj=corrected.gwas,output.file="Cheese_miamiplot.bmp",main="Cheese_conumption")

```
Which results in the following plot.

![](Cheese_miamiplot.bmp)



