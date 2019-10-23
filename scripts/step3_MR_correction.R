library(bGWAS)
library(data.table)
source("../../../fix_model/scripts/Functions_MR_correction.R")

#R --vanilla --slave -q -f ../../scripts/step3_MR_correction.R --args 

args <- commandArgs(TRUE)
i=args[1]
tratto=gsub("../raw_data/model1/imputed/Nutrition_","",gsub(".eu_s_Imputed.gz","",i))

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

mystudies=select_priorGWASs(Z_matrices = "/opt/working/wilson/projects/prj_063_Taste_genes/fix_model/imputed/package/ZMatrices_Nicola_2018-11-12",
                            include_files =trait.f)

r2.table=c()
all.mr.betas=c()
res.t=fread(paste("zcat",i),data.table=F)
names(res.t)[c(1,5,6,11,12)]=c("rsid","a1","a0","beta1","se")

res.bg=  bGWAS(name = tratto,
               GWAS = res.t,
               Z_matrices = "/opt/working/wilson/projects/prj_063_Taste_genes/fix_model/imputed/package/ZMatrices_Nicola_2018-11-12",
               prior_shrinkage = 1,prior_studies = mystudies,
               MR_threshold = 5e-08,
               MR_pruning_dist = 1000,
               MR_pruning_LD = 0.1,
               MR_shrinkage = 1,
               stop_after_prior = T,
               save_files = FALSE,
               verbose=T)



idx=grep("Out-of-sample R-squared across all chromosomes",res.bg$log_info)
r2.table=c(tratto,res.bg$log_info[idx],res.bg$log_info[idx+1])
write(r2.table,ncol=3,file="r2.log",append=T,sep="\t")

all.mr.betas=cbind(tratto,res.bg$significant_studies)
write.table(all.mr.betas,file="Betas.log",append=T,sep="\t",row.names=F,quote=F,col.names=F)


res.corr=gwma.data.frame(bGWAS.obj=res.bg,chi.df=1
                         ,no.se=T,correct.se=T,rs.lab="rsid"
                       ,z_filter=1,results.data.frame=res.t)
write.table(res.corr,file=paste0(tratto,"_corrected_GWAS.tsv"),row.names=F,quote=F,sep="\t")
system(paste0("gzip ",tratto,"_corrected_GWAS.tsv"))


gwama.plot(gwma.obj=res.corr
           ,output.file=paste0(tratto,".png"))

  





