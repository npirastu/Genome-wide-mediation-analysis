library(bGWAS)
library(data.table)

source("scripts/Functions_MR_correction.R")

#R --vanilla --slave -q -f ../../scripts/step3_MR_correction.R --args 


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

mygwas=fread("data/Cheese_gwas_results.tsv.gz",data.table=F)


# Run bGWAS analysis
my_bgwas=  bGWAS(name = "Cheese",
               GWAS = res.t,
               Z_matrices = "data/ZMatrices_Nicola_2018-11-12",
               prior_shrinkage = 1,prior_studies = mystudies,
               MR_threshold = 5e-08,
               MR_pruning_dist = 1000,
               MR_pruning_LD = 0.1,
               MR_shrinkage = 1,
               stop_after_prior = T,
               save_files = FALSE,
               verbose=T)

### Coefficient plot



# Run correction analysis
corrected.gwas=gwma.data.frame(bGWAS.obj=my_bgwas,         # bGWAS output object
                         results.data.frame=mygwas) # data.frame of results




# plot miami plot of results
gwma.miami(gwma.obj=corrected.gwas,output.file="Cheese_miamiplot.bmp",main="Cheese_conumption")







