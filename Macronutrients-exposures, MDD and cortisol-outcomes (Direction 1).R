
####Macronutrients-exposures, depression, cortisol-outcomes
####Loading necessary packages and setting directory

library(remotes)
library(TwoSampleMR)
library(devtools)
library(TwoSampleMR)
library(data.table)
library(gmodels)
library(arrangements)
library(MendelianRandomization)
library(MRPRESSO)
library(MRMix)
library(patchwork)
library(mr.raps)
library(R.utils)
library(ggplot2)
library(ggpubr)
library(LDlinkR)
library(dplyr)
setwd("C:/Users/Matylda/Desktop/MR_files")


####Unpacking macro nutrient data for sugar and fat

R.utils::gunzip(filename = "Diet_Sugar_GWAS_MA_SSGAC_2020_MolPsych.txt.gz")
R.utils::gunzip(filename = "Diet_Fat_GWAS_MA_SSGAC_2020_MolPsych.txt.gz")


########################################################################################################
########################## Sugar-exposure, depression-outcome ############################################

# Getting sugar data
sugar_full <- read.delim("Diet_Sugar_GWAS_MA_SSGAC_2020_MolPsych.txt", header=T)
head(sugar_full)
colnames(sugar_full)
sugar <- data.frame(
  SNP = sugar_full$rsID,
  beta = sugar_full$Beta,
  se = sugar_full$SE,
  effect_allele = sugar_full$A1,
  other_allele = sugar_full$A2,
  pval = sugar_full$Pval,
  eaf = sugar_full$FREQA1_HRC
)
write.csv(sugar, file= "sugar.csv")

sugar_df<-read.csv("sugar.csv")
names(sugar_df)


#Sensitivity
exp_sugar <- format_data(sugar_df, type="exposure",samplesize_col = "samplesize.exposure")
exp_sugar <- exp_sugar[exp_sugar$pval.exposure <=5e-6, ]
exp_sugar_clu<- clump_data(exp_sugar)


#Main
exp_sugar_2<- exp_sugar[exp_sugar$pval.exposure <=5e-8, ]
exp_sugar_2_clu<- clump_data(exp_sugar_2)



###Getting depression data 

dep_full<- read.delim("daner_PGC_MDD_noUKB_no23andMe.txt", header=T, sep="")
head(dep_full)
colnames(dep_full)
dep <- data.frame(
  SNP = dep_full$SNP,
  beta = dep_full$BETA,
  se = dep_full$SE,
  effect_allele = dep_full$A1,
  other_allele = dep_full$A2,
  pval = dep_full$P,
  eaf = dep_full$FRQ_U_95680
)
write.csv(dep, file= "dep.csv")
dep_df<-read.csv("dep.csv")
names(dep_df)

###Prepare outcome data
##Sugar main and depression
out_sugar_2_dep <- read_outcome_data(
  snps = exp_sugar_2_clu$SNP,
  filename = "dep.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "pval",
  eaf_col = "eaf"
)

nrow(out_sugar_2_dep)
nrow(exp_sugar_2_clu)
out_sugar_2_dep_proxy <- out_sugar_2_dep

##Sugar sensitivity and depression
out_sugar_dep <- read_outcome_data(
  snps = exp_sugar_clu$SNP,
  filename = "dep.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "pval",
  eaf_col = "eaf"
)

nrow(out_sugar_dep)
nrow(exp_sugar_clu)
out_sugar_dep_proxy <- out_sugar_dep


### 3. Harmonise data 
dat_sugar_2_dep <- harmonise_data(
  exposure_dat = exp_sugar_2_clu, 
  outcome_dat = out_sugar_2_dep_proxy, 
  action =2)


dat_sugar_dep <- harmonise_data(
  exposure_dat = exp_sugar_clu, 
  outcome_dat = out_sugar_dep_proxy, 
  action =2)

### 4. Drop duplicate exposure-outcome summary sets


## Pruning 
dat_sugar_2_dep$samplesize.outcome<- 138884
dat_sugar_2_dep$samplesize.exposure<- 230648 
dat_sugar_2_dep<-power_prune(dat_sugar_2_dep,method=1,dist.outcome="binary")


dat_sugar_dep$samplesize.outcome<- 138884
dat_sugar_dep$samplesize.exposure<- 230648 
dat_sugar_dep<-power_prune(dat_sugar_dep,method=1,dist.outcome="binary")

## Steiger directionality test 

steiger_filter_sugar_2_dep<-steiger_filtering(dat_sugar_2_dep)
steiger_direct_sugar_2_dep <- directionality_test(dat_sugar_2_dep)


steiger_filter_sugar_dep<-steiger_filtering(dat_sugar_dep)
steiger_direct_sugar_dep <- directionality_test(dat_sugar_dep)


## Edit names of outcome and exposure columns 
dat_sugar_2_dep$exposure <- "Relative sugar intake (p<10e-8)"
dat_sugar_2_dep$outcome<- "MDD"

dat_sugar_dep$exposure <- "Relative sugar intake (p<10e-6)"
dat_sugar_dep$outcome<- "MDD"



### 5. MR analyses 

## Ivw, egger, weighted_median, weighted mode 


results_sugar_2_dep<-mr(dat_sugar_2_dep, method_list=c("mr_ivw_mre",
                                                       "mr_egger_regression",
                                                       "mr_weighted_median", 
                                                       "mr_weighted_mode"))

results_sugar_dep<-mr(dat_sugar_dep, method_list=c("mr_ivw_mre",
                                                         "mr_egger_regression",
                                                         "mr_weighted_median", 
                                                         "mr_weighted_mode"))


results_sugar_dep
results_sugar_2_dep

###Sensitivity analysis
## Raps 
#Main sugar
raps_sugar_2_dep<- mr.raps(dat_sugar_2_dep$beta.exposure, dat_sugar_2_dep$beta.outcome, dat_sugar_2_dep$se.exposure, 
                           dat_sugar_2_dep$se.outcome)

raps_df_sugar_2_dep<- data.frame(id.exposure = results_sugar_2_dep[1,1],
                                 id.outcome = results_sugar_2_dep[1,2], 
                                 outcome = results_sugar_2_dep[1,3],
                                 exposure = results_sugar_2_dep[1,4], 
                                 method = c("mr_raps"), 
                                 nsnp = results_sugar_2_dep[1,6],
                                 b = raps_sugar_2_dep[["beta.hat"]], 
                                 se = raps_sugar_2_dep[["beta.se"]], 
                                 pval = raps_sugar_2_dep[["beta.p.value"]])

results_sugar_2_dep<-rbind(results_sugar_2_dep,raps_df_sugar_2_dep )
results_sugar_2_dep$analysis = "sugar_2_dep"
results_sugar_2_dep <- generate_odds_ratios(results_sugar_2_dep)

hetero_sugar_2_dep<- mr_heterogeneity(dat_sugar_2_dep, method_list=c("mr_ivw_mre", "mr_egger_regression"))
hetero_sugar_2_dep$analysis<- "sugar_2_dep"


#Sensitivity sugar
raps_sugar_dep<- mr.raps(dat_sugar_dep$beta.exposure, dat_sugar_dep$beta.outcome, dat_sugar_dep$se.exposure, 
                            dat_sugar_dep$se.outcome)

raps_df_sugar_dep<- data.frame(id.exposure = results_sugar_dep[1,1],
                                  id.outcome = results_sugar_dep[1,2], 
                                  outcome = results_sugar_dep[1,3],
                                  exposure = results_sugar_dep[1,4], 
                                  method = c("mr_raps"), 
                                  nsnp = results_sugar_dep[1,6],
                                  b = raps_sugar_dep[["beta.hat"]], 
                                  se = raps_sugar_dep[["beta.se"]], 
                                  pval = raps_sugar_dep[["beta.p.value"]])

results_sugar_dep<-rbind(results_sugar_dep,raps_df_sugar_dep )
results_sugar_dep$analysis = "sugar_dep"
results_sugar_dep <- generate_odds_ratios(results_sugar_dep)

hetero_sugar_dep<- mr_heterogeneity(dat_sugar_dep, method_list=c("mr_ivw_mre", "mr_egger_regression"))
hetero_sugar_dep$analysis<- "sugar_dep"



## Check intercept of Egger regression 

#Main sugar

intercept_sugar_2_dep<- mr_pleiotropy_test(dat_sugar_2_dep)
intercept_sugar_2_dep$analysis<- "sugar_2_dep"

#Sensitivity sugar
intercept_sugar_dep<- mr_pleiotropy_test(dat_sugar_dep)
intercept_sugar_dep$analysis<- "sugar_dep"


## MR presso 
#Main sugar
presso_sugar_2_dep<-mr_presso(BetaOutcome = "beta.outcome", 
                              BetaExposure = "beta.exposure", 
                              SdOutcome = "se.outcome", 
                              SdExposure = "se.exposure", 
                              OUTLIERtest = TRUE,
                              DISTORTIONtest = TRUE, 
                              data = dat_sugar_2_dep, 
                              NbDistribution = 1000,  
                              SignifThreshold = 0.05)

presso_sugar_2_dep[["Main MR results"]]
presso_sugar_2_dep$`MR-PRESSO results`$`Global Test`$RSSobs
presso_sugar_2_dep$`MR-PRESSO results`$`Global Test`$Pvalue

presso_sugar_2_dep<- presso_sugar_2_dep[["Main MR results"]]
presso_sugar_2_dep$globaltest<-c(presso_sugar_2_dep$`MR-PRESSO results`$`Global Test`$RSSobs,
                                 presso_sugar_2_dep$`MR-PRESSO results`$`Global Test`$Pvalue)
presso_sugar_2_dep$analysis<- "sugar_2_dep"

#Sensitivity sugar
presso_sugar_dep<-mr_presso(BetaOutcome = "beta.outcome", 
                            BetaExposure = "beta.exposure", 
                            SdOutcome = "se.outcome", 
                            SdExposure = "se.exposure", 
                            OUTLIERtest = TRUE,
                            DISTORTIONtest = TRUE, 
                            data = dat_sugar_dep, 
                            NbDistribution = 1000,  
                            SignifThreshold = 0.05)

presso_sugar_dep[["Main MR results"]]
presso_sugar_dep$`MR-PRESSO results`$`Global Test`$RSSobs
presso_sugar_dep$`MR-PRESSO results`$`Global Test`$Pvalue

presso_sugar_dep<- presso_sugar_dep[["Main MR results"]]
presso_sugar_dep$globaltest<-c(presso_sugar_dep$`MR-PRESSO results`$`Global Test`$RSSobs,
                               presso_sugar_dep$`MR-PRESSO results`$`Global Test`$Pvalue)
presso_sugar_dep$analysis<- "sugar_dep"


##Single SNP analysis 
#Main sugar
res_single_sugar_2_dep <- mr_singlesnp(dat_sugar_2_dep)
res_single_sugar_2_dep_full <- mr_singlesnp(dat_sugar_2_dep, all_method = "mr_two_sample_ml")

res_single_sugar_2_dep <- generate_odds_ratios(res_single_sugar_2_dep)
res_single_sugar_2_dep_full <-generate_odds_ratios(res_single_sugar_2_dep_full)

##Leave-one-out analysis 
res_loo_sugar_2_dep <- mr_leaveoneout(dat_sugar_2_dep)
res_loo_sugar_2_dep <- generate_odds_ratios(res_loo_sugar_2_dep)

#Sensitivity sugar
res_single_sugar_dep<- mr_singlesnp(dat_sugar_dep)
res_single_sugar_dep_full <- mr_singlesnp(dat_sugar_dep, all_method = "mr_two_sample_ml")

res_single_sugar_dep <- generate_odds_ratios(res_single_sugar_dep)
res_single_sugar_dep_full <-generate_odds_ratios(res_single_sugar_dep_full)

##Leave-one-out analysis 
res_loo_sugar_dep <- mr_leaveoneout(dat_sugar_dep)
res_loo_sugar_dep <- generate_odds_ratios(res_loo_sugar_dep)


### 6.Graphs 

##Single SNP
#Main sugar
s_2_d_single <- mr_forest_plot(res_single_sugar_2_dep)
single_s2d_graph <- s_2_d_single[[1]]

s_2_d_single_full <- mr_forest_plot(res_single_sugar_2_dep_full)
s_2_d_single_full[[1]]

#Sensitivity sugar
s_d_single <- mr_forest_plot(res_single_sugar_dep)
single_sd_graph <- s_d_single[[1]]

s_d_single_full <- mr_forest_plot(res_single_sugar_dep_full)
s_d_single_full[[1]]


##Leave one out
#Main sugar
s_2_d_loo <- mr_leaveoneout_plot(res_loo_sugar_2_dep)
loo_s2d_graph <- s_2_d_loo[[1]]

#Sensitivity analysis
s_d_loo <- mr_leaveoneout_plot(res_loo_sugar_dep)
loo_sd_graph <-  s_d_loo[[1]]



########################################################################################################
########################## Fat-exposure,depression-outcome ############################################

# Getting fat data 
fat_full <- read.delim("Diet_Fat_GWAS_MA_SSGAC_2020_MolPsych.txt", header=T)
head(fat_full)
colnames(fat_full)
fat <- data.frame(
  SNP = fat_full$rsID,
  beta = fat_full$Beta,
  se = fat_full$SE,
  effect_allele = fat_full$A1,
  other_allele = fat_full$A2,
  pval = fat_full$Pval,
  eaf = fat_full$FREQA1_HRC
)
write.csv(fat, file= "fat.csv")

fat_df<-read.csv("fat.csv")
names(fat_df)

#Sensitivity
exp_fat <- format_data(fat_df, type="exposure",samplesize_col = "samplesize.exposure")
exp_fat<- exp_fat[exp_fat$pval.exposure <=5e-6, ]
exp_fat_clu<- clump_data(exp_fat)

#Main
exp_fat_2<- exp_fat[exp_fat$pval.exposure <=5e-8, ]
exp_fat_2_clu<- clump_data(exp_fat_2)

###Getting depression data 
dep_df<-read.csv("dep.csv")
names(dep_df)

###Prepare outcome data
##Fat main and depression
out_fat_2_dep <- read_outcome_data(
  snps = exp_fat_2_clu$SNP,
  filename = "dep.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "pval",
  eaf_col = "eaf"
)

nrow(out_fat_2_dep)
nrow(exp_fat_2_clu)
out_fat_2_dep_proxy <- out_fat_2_dep

#Sensitivity fat and depression
out_fat_dep <- read_outcome_data(
  snps = exp_fat_clu$SNP,
  filename = "dep.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "pval",
  eaf_col = "eaf"
)

nrow(out_fat_dep)
nrow(exp_fat_clu)
out_fat_dep_proxy <- out_fat_dep



### 3. Harmonise data 

dat_fat_2_dep <- harmonise_data(
  exposure_dat = exp_fat_2_clu, 
  outcome_dat = out_fat_2_dep_proxy, 
  action =2)

dat_fat_dep <- harmonise_data(
  exposure_dat = exp_fat_clu, 
  outcome_dat = out_fat_dep_proxy, 
  action =2)




### 4. Drop duplicate exposure-outcome summary sets


## Pruning 

dat_fat_2_dep$samplesize.outcome<- 138884
dat_fat_2_dep$samplesize.exposure<- 264181
dat_fat_2_dep<-power_prune(dat_fat_2_dep,method=1,dist.outcome="binary")

dat_fat_dep$samplesize.outcome<- 138884
dat_fat_dep$samplesize.exposure<- 264181
dat_fat_dep<-power_prune(dat_fat_dep,method=1,dist.outcome="binary")



## Steiger directionality test 

steiger_filter_fat_2_dep<-steiger_filtering(dat_fat_2_dep)
steiger_direct_fat_2_dep <- directionality_test(dat_fat_2_dep)

steiger_filter_fat_dep<-steiger_filtering(dat_fat_dep)
steiger_direct_fat_dep <- directionality_test(dat_fat_dep)



## Edit names of outcome and exposure columns 

dat_fat_2_dep$exposure <- "Relative fat intake (p<10e-8)"
dat_fat_2_dep$outcome<- "MDD"

dat_fat_dep$exposure <- "Relative fat intake (p<10e-6)"
dat_fat_dep$outcome<- "MDD"



### 5. MR analyses 

## Ivw, egger, weighted_median, weighted mode 
results_fat_2_dep<-mr(dat_fat_2_dep, method_list=c("mr_ivw_mre",
                                                   "mr_egger_regression",
                                                   "mr_weighted_median", 
                                                   "mr_weighted_mode"))

results_fat_dep<-mr(dat_fat_dep, method_list=c("mr_ivw_mre",
                                               "mr_egger_regression",
                                               "mr_weighted_median", 
                                               "mr_weighted_mode"))



results_fat_2_dep
results_fat_dep


###Sensitivity analysis
## Raps 
#Main fat
raps_fat_2_dep<- mr.raps(dat_fat_2_dep$beta.exposure, dat_fat_2_dep$beta.outcome, dat_fat_2_dep$se.exposure, 
                         dat_fat_2_dep$se.outcome)

raps_df_fat_2_dep<- data.frame(id.exposure = results_fat_2_dep[1,1],
                               id.outcome = results_fat_2_dep[1,2], 
                               outcome = results_fat_2_dep[1,3],
                               exposure = results_fat_2_dep[1,4], 
                               method = c("mr_raps"), 
                               nsnp = results_fat_2_dep[1,6],
                               b = raps_fat_2_dep[["beta.hat"]], 
                               se = raps_fat_2_dep[["beta.se"]], 
                               pval = raps_fat_2_dep[["beta.p.value"]])

results_fat_2_dep<-rbind(results_fat_2_dep,raps_df_fat_2_dep )
results_fat_2_dep$analysis = "fat_2_dep"
results_fat_2_dep <- generate_odds_ratios(results_fat_2_dep)


hetero_fat_2_dep<-mr_heterogeneity(dat_fat_2_dep, method_list=c("mr_ivw_mre", "mr_egger_regression"))
hetero_fat_2_dep$analysis<- "fat_2_dep"

#Sensitivity fat
raps_fat_dep<- mr.raps(dat_fat_dep$beta.exposure, dat_fat_dep$beta.outcome, dat_fat_dep$se.exposure, 
                       dat_fat_dep$se.outcome)

raps_df_fat_dep<- data.frame(id.exposure = results_fat_dep[1,1],
                             id.outcome = results_fat_dep[1,2], 
                             outcome = results_fat_dep[1,3],
                             exposure = results_fat_dep[1,4], 
                             method = c("mr_raps"), 
                             nsnp = results_fat_dep[1,6],
                             b = raps_fat_dep[["beta.hat"]], 
                             se = raps_fat_dep[["beta.se"]], 
                             pval = raps_fat_dep[["beta.p.value"]])

results_fat_dep<-rbind(results_fat_dep,raps_df_fat_dep )
results_fat_dep$analysis = "fat_dep"
results_fat_dep <- generate_odds_ratios(results_fat_dep)

hetero_fat_dep<- mr_heterogeneity(dat_fat_dep, method_list=c("mr_ivw_mre", "mr_egger_regression"))
hetero_fat_dep$analysis<- "fat_dep"



## Check intercept of Egger regression 

#Main fat
intercept_fat_2_dep<- mr_pleiotropy_test(dat_fat_2_dep)
intercept_fat_2_dep$analysis<- "fat_2_dep"

#Sensitivity fat
intercept_fat_dep<- mr_pleiotropy_test(dat_fat_dep)
intercept_fat_dep$analysis<- "fat_dep"


## MR presso 
#Main fat
presso_fat_2_dep<-mr_presso(BetaOutcome = "beta.outcome", 
                            BetaExposure = "beta.exposure", 
                            SdOutcome = "se.outcome", 
                            SdExposure = "se.exposure", 
                            OUTLIERtest = TRUE,
                            DISTORTIONtest = TRUE, 
                            data = dat_fat_2_dep, 
                            NbDistribution = 1000,  
                            SignifThreshold = 0.05)

presso_fat_2_dep[["Main MR results"]]
presso_fat_2_dep$`MR-PRESSO results`$`Global Test`$RSSobs
presso_fat_2_dep$`MR-PRESSO results`$`Global Test`$Pvalue

presso_fat_2_dep<- presso_fat_2_dep[["Main MR results"]]
presso_fat_2_dep$globaltest<-c(presso_fat_2_dep$`MR-PRESSO results`$`Global Test`$RSSobs,
                               presso_fat_2_dep$`MR-PRESSO results`$`Global Test`$Pvalue)
presso_fat_2_dep$analysis<- "fat_2_dep"

#Sensitivity fat
presso_fat_dep<-mr_presso(BetaOutcome = "beta.outcome", 
                          BetaExposure = "beta.exposure", 
                          SdOutcome = "se.outcome", 
                          SdExposure = "se.exposure", 
                          OUTLIERtest = TRUE,
                          DISTORTIONtest = TRUE, 
                          data = dat_fat_dep, 
                          NbDistribution = 1000,  
                          SignifThreshold = 0.05)

presso_fat_dep[["Main MR results"]]
presso_fat_dep$`MR-PRESSO results`$`Global Test`$RSSobs
presso_fat_dep$`MR-PRESSO results`$`Global Test`$Pvalue

presso_fat_dep<- presso_fat_dep[["Main MR results"]]
presso_fat_dep$globaltest<-c(presso_fat_dep$`MR-PRESSO results`$`Global Test`$RSSobs,
                             presso_fat_dep$`MR-PRESSO results`$`Global Test`$Pvalue)
presso_fat_dep$analysis<- "fat_dep"



##Single SNP analysis 
#Main fat
res_single_fat_2_dep <- mr_singlesnp(dat_fat_2_dep)
res_single_fat_2_dep_full <- mr_singlesnp(dat_fat_2_dep, all_method = "mr_two_sample_ml")

res_single_fat_2_dep <- generate_odds_ratios(res_single_fat_2_dep)
res_single_fat_2_dep_full <-generate_odds_ratios(res_single_fat_2_dep_full)


##Leave-one-out analysis 
res_loo_fat_2_dep <- mr_leaveoneout(dat_fat_2_dep)
res_loo_fat_2_dep <- generate_odds_ratios(res_loo_fat_2_dep)

#Sensitivity fat
res_single_fat_dep<- mr_singlesnp(dat_fat_dep)
res_single_fat_dep_full <- mr_singlesnp(dat_fat_dep, all_method = "mr_two_sample_ml")

res_single_fat_dep <- generate_odds_ratios(res_single_fat_dep)
res_single_fat_dep_full <-generate_odds_ratios(res_single_fat_dep_full)


##Leave-one-out analysis 
res_loo_fat_dep <- mr_leaveoneout(dat_fat_dep)
res_loo_fat_dep <- generate_odds_ratios(res_loo_fat_dep)



### 6.Graphs 

##Single SNP
#Main fat
f_2_d_single <- mr_forest_plot(res_single_fat_2_dep)
single_f2d_graph <- f_2_d_single[[1]]

f_2_d_single_full <- mr_forest_plot(res_single_fat_2_dep_full)
f_2_d_single_full[[1]]

#Sensitivity fat
f_d_single <- mr_forest_plot(res_single_fat_dep)
single_fd_graph <- f_d_single[[1]]

f_d_single_full <- mr_forest_plot(res_single_fat_dep_full)
f_d_single_full[[1]]



##Leave one out
#Main fat
f_2_d_loo <- mr_leaveoneout_plot(res_loo_fat_2_dep)
loo_f2d_graph <- f_2_d_loo[[1]]

#Sensitivity fat
f_d_loo <- mr_leaveoneout_plot(res_loo_fat_dep)
loo_fd_graph <- f_d_loo[[1]]



########################################################################################################
########################## Sugar-exposure, Cortisol-outcome ############################################


###Getting sugar data 
sugar_df<-read.csv("sugar.csv")
names(sugar_df)

#Sensitivity
exp_sugar <- format_data(sugar_df, type="exposure",samplesize_col = "samplesize.exposure")
exp_sugar<- exp_sugar[exp_sugar$pval.exposure <=5e-6, ]
exp_sugar_clu<- clump_data(exp_sugar)

#Main
exp_sugar_2<- exp_sugar[exp_sugar$pval.exposure <=5e-8, ]
exp_sugar_2_clu<- clump_data(exp_sugar_2)


### 2. Prepare outcome data (cortisol)
#Sugar main
out_sugar_2_cortisol <- extract_outcome_data(
  snps = exp_sugar_2_clu$SNP,
  outcomes = "ieu-a-1012"
)

# Note: proxy SNPs were automatically added when using code above

nrow(out_sugar_2_cortisol)
nrow(exp_sugar_2_clu)


#Sugar sensitivity
out_sugar_cortisol <- extract_outcome_data(
  snps = exp_sugar_clu$SNP,
  outcomes = "ieu-a-1012"
)

# Note: proxy SNPs were automatically added when using code above

nrow(out_sugar_cortisol)
nrow(exp_sugar_clu)



### 3. Harmonise data 
#Main sugar
dat_sugar_2_cortisol <- harmonise_data(
  exposure_dat = exp_sugar_2_clu, 
  outcome_dat = out_sugar_2_cortisol, 
  action =2)

#Sensitivity sugar
dat_sugar_cortisol <- harmonise_data(
  exposure_dat = exp_sugar_clu, 
  outcome_dat = out_sugar_cortisol, 
  action =2)


### 4. Drop duplicate exposure-outcome summary sets
## Pruning 
dat_sugar_2_cortisol$samplesize.exposure <- 230648 
dat_sugar_2_cortisol$exposure <- "exposure"
dat_sugar_2_cortisol$data_source.exposure <- "textfile"
dat_sugar_2_cortisol$pval_origin.outcome <- "reported"
dat_sugar_2_cortisol<-power_prune(dat_sugar_2_cortisol,method=1,dist.outcome="continuous")


dat_sugar_cortisol$samplesize.exposure <- 230648 
dat_sugar_cortisol$exposure <- "exposure"
dat_sugar_cortisol$data_source.exposure <- "textfile"
dat_sugar_cortisol$pval_origin.outcome <- "reported"
dat_sugar_cortisol<-power_prune(dat_sugar_cortisol,method=1,dist.outcome="continuous")


## Steiger directionality test 
steiger_filter_sugar_2_cortisol<-steiger_filtering(dat_sugar_2_cortisol)
steiger_direct_sugar_2_cortisol <- directionality_test(dat_sugar_2_cortisol)

steiger_filter_sugar_cortisol<-steiger_filtering(dat_sugar_cortisol)
steiger_direct_sugar_cortisol <- directionality_test(dat_sugar_cortisol)


## Edit names of outcome and exposure columns 
dat_sugar_2_cortisol$exposure <- "Relative sugar intake (p<10e-8)"
dat_sugar_2_cortisol$outcome<- "Plasma cortisol"

dat_sugar_cortisol$exposure <- "Relative sugar intake (p<10e-6)"
dat_sugar_cortisol$outcome<- "Plasma cortisol"


### 5. MR analyses 

## Ivw, egger, weighted_median, weighted mode 
results_sugar_2_cortisol<-mr(dat_sugar_2_cortisol, method_list=c("mr_ivw_mre",
                                                                 "mr_egger_regression",
                                                                 "mr_weighted_median", 
                                                                 "mr_weighted_mode"))

results_sugar_cortisol<-mr(dat_sugar_cortisol, method_list=c("mr_ivw_mre",
                                                             "mr_egger_regression",
                                                             "mr_weighted_median", 
                                                             "mr_weighted_mode"))


results_sugar_2_cortisol
results_sugar_cortisol


###Sensitivity analysis
## Raps 
#Main sugar
raps_sugar_2_cortisol<- mr.raps(dat_sugar_2_cortisol$beta.exposure, dat_sugar_2_cortisol$beta.outcome, dat_sugar_2_cortisol$se.exposure, 
                                dat_sugar_2_cortisol$se.outcome)

raps_df_sugar_2_cortisol<- data.frame(id.exposure = results_sugar_2_cortisol[1,1],
                                      id.outcome = results_sugar_2_cortisol[1,2], 
                                      outcome = results_sugar_2_cortisol[1,3],
                                      exposure = results_sugar_2_cortisol[1,4], 
                                      method = c("mr_raps"), 
                                      nsnp = results_sugar_2_cortisol[1,6],
                                      b = raps_sugar_2_cortisol[["beta.hat"]], 
                                      se = raps_sugar_2_cortisol[["beta.se"]], 
                                      pval = raps_sugar_2_cortisol[["beta.p.value"]])

results_sugar_2_cortisol<-rbind(results_sugar_2_cortisol,raps_df_sugar_2_cortisol )
results_sugar_2_cortisol$analysis = "sugar_2_cortisol"
results_sugar_2_cortisol <- generate_odds_ratios(results_sugar_2_cortisol)

hetero_sugar_2_cortisol<- mr_heterogeneity(dat_sugar_2_cortisol, method_list=c("mr_ivw_mre", "mr_egger_regression"))
hetero_sugar_2_cortisol$analysis<- "sugar_2_cortisol"


#Sensitivity sugar
raps_sugar_cortisol<- mr.raps(dat_sugar_cortisol$beta.exposure, dat_sugar_cortisol$beta.outcome, dat_sugar_cortisol$se.exposure, 
                              dat_sugar_cortisol$se.outcome)

raps_df_sugar_cortisol<- data.frame(id.exposure = results_sugar_cortisol[1,1],
                                    id.outcome = results_sugar_cortisol[1,2], 
                                    outcome = results_sugar_cortisol[1,3],
                                    exposure = results_sugar_cortisol[1,4], 
                                    method = c("mr_raps"), 
                                    nsnp = results_sugar_cortisol[1,6],
                                    b = raps_sugar_cortisol[["beta.hat"]], 
                                    se = raps_sugar_cortisol[["beta.se"]], 
                                    pval = raps_sugar_cortisol[["beta.p.value"]])

results_sugar_cortisol<-rbind(results_sugar_cortisol,raps_df_sugar_cortisol )
results_sugar_cortisol$analysis = "sugar_cortisol"
results_sugar_cortisol <- generate_odds_ratios(results_sugar_cortisol)

hetero_sugar_cortisol<- mr_heterogeneity(dat_sugar_cortisol, method_list=c("mr_ivw_mre", "mr_egger_regression"))
hetero_sugar_cortisol$analysis<- "sugar_cortisol"




## Check intercept of Egger regression 

#Main sugar
intercept_sugar_2_cortisol<- mr_pleiotropy_test(dat_sugar_2_cortisol)
intercept_sugar_2_cortisol$analysis<- "sugar_2_cortisol"

#Sensitivity sugar
intercept_sugar_cortisol<- mr_pleiotropy_test(dat_sugar_cortisol)
intercept_sugar_cortisol$analysis<- "sugar_cortisol"


## MR presso 
#Main sugar
presso_sugar_2_cortisol<-mr_presso(BetaOutcome = "beta.outcome", 
                                   BetaExposure = "beta.exposure", 
                                   SdOutcome = "se.outcome", 
                                   SdExposure = "se.exposure", 
                                   OUTLIERtest = TRUE,
                                   DISTORTIONtest = TRUE, 
                                   data = dat_sugar_2_cortisol, 
                                   NbDistribution = 1000,  
                                   SignifThreshold = 0.05)

presso_sugar_2_cortisol[["Main MR results"]]
presso_sugar_2_cortisol$`MR-PRESSO results`$`Global Test`$RSSobs
presso_sugar_2_cortisol$`MR-PRESSO results`$`Global Test`$Pvalue

presso_sugar_2_cortisol<- presso_sugar_2_cortisol[["Main MR results"]]
presso_sugar_2_cortisol$globaltest<-c(presso_sugar_2_cortisol$`MR-PRESSO results`$`Global Test`$RSSobs,
                                      presso_sugar_2_cortisol$`MR-PRESSO results`$`Global Test`$Pvalue)
presso_sugar_2_cortisol$analysis<- "sugar_2_cortisol"


#Sensitivity sugar
presso_sugar_cortisol<-mr_presso(BetaOutcome = "beta.outcome", 
                                 BetaExposure = "beta.exposure", 
                                 SdOutcome = "se.outcome", 
                                 SdExposure = "se.exposure", 
                                 OUTLIERtest = TRUE,
                                 DISTORTIONtest = TRUE, 
                                 data = dat_sugar_cortisol, 
                                 NbDistribution = 1000,  
                                 SignifThreshold = 0.05)

presso_sugar_cortisol[["Main MR results"]]
presso_sugar_cortisol$`MR-PRESSO results`$`Global Test`$RSSobs
presso_sugar_cortisol$`MR-PRESSO results`$`Global Test`$Pvalue

presso_sugar_cortisol<- presso_sugar_cortisol[["Main MR results"]]
presso_sugar_cortisol$globaltest<-c(presso_sugar_cortisol$`MR-PRESSO results`$`Global Test`$RSSobs,
                                    presso_sugar_cortisol$`MR-PRESSO results`$`Global Test`$Pvalue)
presso_sugar_cortisol$analysis<- "sugar_cortisol"



##Single SNP analysis 
#Main sugar
res_single_sugar_2_cortisol <- mr_singlesnp(dat_sugar_2_cortisol)
res_single_sugar_2_cortisol_full <- mr_singlesnp(dat_sugar_2_cortisol, all_method = "mr_two_sample_ml")

res_single_sugar_2_cortisol <- generate_odds_ratios(res_single_sugar_2_cortisol)
res_single_sugar_2_cortisol_full <-generate_odds_ratios(res_single_sugar_2_cortisol_full)

##Leave-one-out analysis 
res_loo_sugar_2_cortisol <- mr_leaveoneout(dat_sugar_2_cortisol)
res_loo_sugar_2_cortisol <- generate_odds_ratios(res_loo_sugar_2_cortisol)


#Sensitivity sugar
res_single_sugar_cortisol<- mr_singlesnp(dat_sugar_cortisol)
res_single_sugar_cortisol_full <- mr_singlesnp(dat_sugar_cortisol, all_method = "mr_two_sample_ml")

res_single_sugar_cortisol <- generate_odds_ratios(res_single_sugar_cortisol)
res_single_sugar_cortisol_full <-generate_odds_ratios(res_single_sugar_cortisol_full)

##Leave-one-out analysis 
res_loo_sugar_cortisol <- mr_leaveoneout(dat_sugar_cortisol)
res_loo_sugar_cortisol <- generate_odds_ratios(res_loo_sugar_cortisol)


### 6.Graphs 


##Single SNP
#Main sugar
s_2_c_single <- mr_forest_plot(res_single_sugar_2_cortisol)
single_s2c_graph <- s_2_c_single[[1]]

s_2_c_single_full <- mr_forest_plot(res_single_sugar_2_cortisol_full)
s_2_c_single_full[[1]]

#Sensitivity sugar
s_c_single <- mr_forest_plot(res_single_sugar_cortisol)
single_sc_graph <- s_c_single[[1]]

s_c_single_full <- mr_forest_plot(res_single_sugar_cortisol_full)
 s_c_single_full[[1]]



##Leave one out
#Main sugar
s_2_c_loo <- mr_leaveoneout_plot(res_loo_sugar_2_cortisol)
loo_s2c_graph <- s_2_c_loo[[1]]

#Sensitivity sugar
s_c_loo <- mr_leaveoneout_plot(res_loo_sugar_cortisol)
loo_sc_graph <- s_c_loo[[1]]



################################################################################################################
####################################Fat-exposure, cortisol-outcome############################################

###Getting fat data 
fat_df<-read.csv("fat.csv")
names(fat_df)

#Sensitivity
exp_fat <- format_data(fat_df, type="exposure",samplesize_col = "samplesize.exposure")
exp_fat<- exp_fat[exp_fat$pval.exposure <=5e-6, ]
exp_fat_clu<- clump_data(exp_fat)

#Main
exp_fat_2<- exp_fat[exp_fat$pval.exposure <=5e-8, ]
exp_fat_2_clu<- clump_data(exp_fat_2)


### 2. Prepare outcome data (cortisol)
#Main fat
out_fat_2_cortisol <- extract_outcome_data(
  snps = exp_fat_2_clu$SNP,
  outcomes = "ieu-a-1012"
)

# Note: proxy SNPs were automatically added when using code above

nrow(out_fat_2_cortisol)
nrow(exp_fat_2_clu)

#Sensitivity fat
out_fat_cortisol <- extract_outcome_data(
  snps = exp_fat_clu$SNP,
  outcomes = "ieu-a-1012"
)

# Note: proxy SNPs were automatically added when using code above

nrow(out_fat_cortisol)
nrow(exp_fat_clu)


### 3. Harmonise data 
#Main fat
dat_fat_2_cortisol <- harmonise_data(
  exposure_dat = exp_fat_2_clu, 
  outcome_dat = out_fat_2_cortisol, 
  action =2)

#Sensitivity fat
dat_fat_cortisol <- harmonise_data(
  exposure_dat = exp_fat_clu, 
  outcome_dat = out_fat_cortisol, 
  action =2)



### 4. Drop duplicate exposure-outcome summary sets
## Pruning 
dat_fat_2_cortisol$samplesize.exposure <- 230648 
dat_fat_2_cortisol<-power_prune(dat_fat_2_cortisol,method=1,dist.outcome="continuous")

dat_fat_cortisol$samplesize.exposure <- 230648 
dat_fat_cortisol<-power_prune(dat_fat_cortisol,method=1,dist.outcome="continuous")


## Steiger directionality test 
steiger_filter_fat_2_cortisol<-steiger_filtering(dat_fat_2_cortisol)
steiger_direct_fat_2_cortisol <- directionality_test(dat_fat_2_cortisol)

steiger_filter_fat_cortisol<-steiger_filtering(dat_fat_cortisol)
steiger_direct_fat_cortisol <- directionality_test(dat_fat_cortisol)


## Edit names of outcome and exposure columns 
dat_fat_2_cortisol$exposure <- "Relative fat intake (p<10e-8)"
dat_fat_2_cortisol$outcome<- "Plasma cortisol"

dat_fat_cortisol$exposure <- "Relative fat intake (p<10e-6)"
dat_fat_cortisol$outcome<- "Plasma cortisol"


### 5. MR analyses 

## Ivw, egger, weighted_median, weighted mode 

results_fat_cortisol<-mr(dat_fat_cortisol, method_list=c("mr_ivw_mre",
                                                         "mr_egger_regression",
                                                         "mr_weighted_median", 
                                                         "mr_weighted_mode"))


results_fat_2_cortisol<-mr(dat_fat_2_cortisol, method_list=c("mr_wald_ratio"))

results_fat_cortisol
results_fat_2_cortisol

###Sensitivity analysis
## Raps 
#Main fat
raps_fat_2_cortisol<- mr.raps(dat_fat_2_cortisol$beta.exposure, dat_fat_2_cortisol$beta.outcome, dat_fat_2_cortisol$se.exposure, 
                              dat_fat_2_cortisol$se.outcome)

raps_df_fat_2_cortisol<- data.frame(id.exposure = results_fat_2_cortisol[1,1],
                                    id.outcome = results_fat_2_cortisol[1,2], 
                                    outcome = results_fat_2_cortisol[1,3],
                                    exposure = results_fat_2_cortisol[1,4], 
                                    method = c("mr_raps"), 
                                    nsnp = results_fat_2_cortisol[1,6],
                                    b = raps_fat_2_cortisol[["beta.hat"]], 
                                    se = raps_fat_2_cortisol[["beta.se"]], 
                                    pval = raps_fat_2_cortisol[["beta.p.value"]])

results_fat_2_cortisol<-rbind(results_fat_2_cortisol,raps_df_fat_2_cortisol )
results_fat_2_cortisol$analysis = "fat_2_cortisol"
results_fat_2_cortisol <- generate_odds_ratios(results_fat_2_cortisol)

#Sensitivity fat
raps_fat_cortisol<- mr.raps(dat_fat_cortisol$beta.exposure, dat_fat_cortisol$beta.outcome, dat_fat_cortisol$se.exposure, 
                            dat_fat_cortisol$se.outcome)

raps_df_fat_cortisol<- data.frame(id.exposure = results_fat_cortisol[1,1],
                                  id.outcome = results_fat_cortisol[1,2], 
                                  outcome = results_fat_cortisol[1,3],
                                  exposure = results_fat_cortisol[1,4], 
                                  method = c("mr_raps"), 
                                  nsnp = results_fat_cortisol[1,6],
                                  b = raps_fat_cortisol[["beta.hat"]], 
                                  se = raps_fat_cortisol[["beta.se"]], 
                                  pval = raps_fat_cortisol[["beta.p.value"]])

results_fat_cortisol<-rbind(results_fat_cortisol,raps_df_fat_cortisol )
results_fat_cortisol$analysis = "fat_cortisol"
results_fat_cortisol <- generate_odds_ratios(results_fat_cortisol)

hetero_fat_cortisol<- mr_heterogeneity(dat_fat_cortisol, method_list=c("mr_ivw_mre", "mr_egger_regression"))
hetero_fat_cortisol$analysis<- "fat_cortisol"


##Single SNP analysis 
#Sensitivity fat
res_single_fat_cortisol<- mr_singlesnp(dat_fat_cortisol)
res_single_fat_cortisol_full <- mr_singlesnp(dat_fat_cortisol, all_method = "mr_two_sample_ml")

res_single_fat_cortisol <- generate_odds_ratios(res_single_fat_cortisol)
res_single_fat_cortisol_full <-generate_odds_ratios(res_single_fat_cortisol_full)

##Leave-one-out analysis 
res_loo_fat_cortisol <- mr_leaveoneout(dat_fat_cortisol)
res_loo_fat_cortisol <- generate_odds_ratios(res_loo_fat_cortisol)


### 6.Graphs 

##Single SNP

#Sensitivity fat
f_c_single <- mr_forest_plot(res_single_fat_cortisol)
single_fc_graph <- f_c_single[[1]]

f_c_single_full <- mr_forest_plot(res_single_fat_cortisol_full)
f_c_single_full[[1]]


###################Figures#########################

single_sensitivity_1<-ggarrange(
  single_fd_graph, single_sd_graph, single_fc_graph, single_sc_graph, labels = c("A", "B", "C", "D"))      

single_main_1<-ggarrange(
  single_f2d_graph, single_s2d_graph, single_s2c_graph, labels = c("A", "B", "C")) 

loo_sensitivity_1 <- ggarrange(
  loo_fd_graph, loo_sd_graph, loo_sc_graph, labels = c("A", "B", "C"))   

loo_main_1 <- ggarrange(
  loo_f2d_graph, loo_s2d_graph, loo_s2c_graph, labels = c("A", "B", "C"))   

#################################################################################################
####################################Combined results#############################################
results_macro_exposures<- rbind(results_fat_2_cortisol,
                        results_fat_cortisol,
                        results_fat_2_dep,
                        results_fat_dep,
                        results_sugar_2_cortisol,
                        results_sugar_cortisol,
                        results_sugar_2_dep,
                        results_sugar_dep)

hetero_macro_exposures<-rbind(hetero_fat_2_cortisol,
                              hetero_fat_cortisol,
                              hetero_fat_dep,
                              hetero_fat_2_dep,
                              hetero_sugar_2_cortisol,
                              hetero_sugar_cortisol,
                              hetero_sugar_2_dep,
                              hetero_sugar_dep)

intercept_macro_exposures<- rbind(intercept_fat_2_cortisol,
                                  intercept_fat_cortisol,
                                  intercept_fat_2_dep,
                         intercept_fat_dep,
                         intercept_sugar_2_cortisol,
                         intercept_sugar_cortisol,
                         intercept_sugar_2_dep,
                         intercept_sugar_dep)

steiger_macro_exposures<-rbind(steiger_direct_fat_2_dep,
                       steiger_direct_fat_2_dep,
                       steiger_direct_fat_2_cortisol,
                       steiger_direct_fat_cortisol,
                       steiger_direct_sugar_2_cortisol,
                       steiger_direct_sugar_cortisol,
                       steiger_direct_sugar_2_dep,
                       steiger_direct_sugar_dep)

presso_macro_exposures<- rbind(presso_fat_2_dep,
                       presso_fat_dep,
                       presso_sugar_2_cortisol,
                       presso_sugar_cortisol,
                       presso_sugar_2_dep,
                       presso_sugar_dep)