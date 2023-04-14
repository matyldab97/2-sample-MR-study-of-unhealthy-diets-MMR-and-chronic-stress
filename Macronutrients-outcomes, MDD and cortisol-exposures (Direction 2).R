####Cortisol and depression-exposures, macronutrients-outcomes
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
########################## Cortisol-exposure, Sugar-outcome ############################################

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


###Getting cortisol data 
#Main 
exp_cortisol_2<- extract_instruments("ieu-a-1012", p1=5e-8)
exp_cortisol_2_clu<- clump_data(exp_cortisol_2)


#Sensitivity
exp_cortisol<- extract_instruments("ieu-a-1012", p1=5e-6)
exp_cortisol_clu<- clump_data(exp_cortisol)



###Prepare outcome data
##Cortisol main and sugar
out_cortisol_2_sugar <- read_outcome_data(
  snps = exp_cortisol_2_clu$SNP,
  filename = "sugar.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "pval",
  eaf_col = "eaf"
)

nrow(out_cortisol_2_sugar)
nrow(exp_cortisol_2_clu)
out_cortisol_2_sugar_proxy <- out_cortisol_2_sugar

#Sensitivity cortisol and sugar
out_cortisol_sugar <- read_outcome_data(
  snps = exp_cortisol_clu$SNP,
  filename = "sugar.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "pval",
  eaf_col = "eaf"
)

nrow(out_cortisol_sugar)
nrow(exp_cortisol_clu)
out_cortisol_sugar_proxy <- out_cortisol_sugar


### 3. Harmonise data 

dat_cortisol_2_sugar <- harmonise_data(
  exposure_dat = exp_cortisol_2_clu, 
  outcome_dat = out_cortisol_2_sugar_proxy, 
  action =2)

dat_cortisol_sugar <- harmonise_data(
  exposure_dat = exp_cortisol_clu, 
  outcome_dat = out_cortisol_sugar_proxy, 
  action =2)

### 4. Drop duplicate exposure-outcome summary sets

## Pruning 
dat_cortisol_2_sugar$samplesize.outcome<- 230648
dat_cortisol_2_sugar<-power_prune(dat_cortisol_2_sugar,method=2,dist.outcome="continuous")

dat_cortisol_sugar$samplesize.outcome<- 230648 
dat_cortisol_sugar<-power_prune(dat_cortisol_sugar,method=2,dist.outcome="continuous")


## Steiger directionality test 
steiger_filter_cortisol_2_sugar<-steiger_filtering(dat_cortisol_2_sugar)
steiger_direct_cortisol_2_sugar <- directionality_test(dat_cortisol_2_sugar)

steiger_filter_cortisol_sugar<-steiger_filtering(dat_cortisol_sugar)
steiger_direct_cortisol_sugar <- directionality_test(dat_cortisol_sugar)


## Edit names of outcome and exposure columns 
dat_cortisol_2_sugar$exposure <- "Plasma cortisol levels (p<10e-8)"
dat_cortisol_2_sugar$outcome <- "Relative sugar intake"


dat_cortisol_sugar$exposure <- "Plasma cortisol levels (p<10e-6)"
dat_cortisol_sugar$outcome <- "Relative sugar intake"



### 5. MR analyses 

## Ivw, egger, weighted_median, weighted mode 
results_cortisol_2_sugar<-mr(dat_cortisol_2_sugar, method_list=c("mr_wald_ratio"))

results_cortisol_sugar<-mr(dat_cortisol_sugar, method_list=c("mr_ivw_mre",
                                                                     "mr_egger_regression",
                                                                     "mr_weighted_median", 
                                                                     "mr_weighted_mode"))



results_cortisol_2_sugar
results_cortisol_sugar


###Sensitivity analysis
## Raps 
#Main cortisol
raps_cortisol_2_sugar<- mr.raps(dat_cortisol_2_sugar$beta.exposure, dat_cortisol_2_sugar$beta.outcome, dat_cortisol_2_sugar$se.exposure, 
                                dat_cortisol_2_sugar$se.outcome)

raps_df_cortisol_2_sugar<- data.frame(id.exposure = results_cortisol_2_sugar[1,1],
                                      id.outcome = results_cortisol_2_sugar[1,2], 
                                      outcome = results_cortisol_2_sugar[1,3],
                                      exposure = results_cortisol_2_sugar[1,4], 
                                      method = c("mr_raps"), 
                                      nsnp = results_cortisol_2_sugar[1,6],
                                      b = raps_cortisol_2_sugar[["beta.hat"]], 
                                      se = raps_cortisol_2_sugar[["beta.se"]], 
                                      pval = raps_cortisol_2_sugar[["beta.p.value"]])

results_cortisol_2_sugar<-rbind(results_cortisol_2_sugar,raps_df_cortisol_2_sugar )
results_cortisol_2_sugar <- generate_odds_ratios(results_cortisol_2_sugar)
results_cortisol_2_sugar$analysis = "cortisol_2_sugar"

#Sensitivity cortisol
raps_cortisol_sugar<- mr.raps(dat_cortisol_sugar$beta.exposure, dat_cortisol_sugar$beta.outcome, dat_cortisol_sugar$se.exposure, 
               dat_cortisol_sugar$se.outcome)

raps_df_cortisol_sugar<- data.frame(id.exposure = results_cortisol_sugar[1,1],
                     id.outcome = results_cortisol_sugar[1,2], 
                     outcome = results_cortisol_sugar[1,3],
                     exposure = results_cortisol_sugar[1,4], 
                     method = c("mr_raps"), 
                     nsnp = results_cortisol_sugar[1,6],
                     b = raps_cortisol_sugar[["beta.hat"]], 
                     se = raps_cortisol_sugar[["beta.se"]], 
                     pval = raps_cortisol_sugar[["beta.p.value"]])

results_cortisol_sugar<-rbind(results_cortisol_sugar,raps_df_cortisol_sugar )
results_cortisol_sugar$analysis = "cortisol_sugar"
results_cortisol_sugar <- generate_odds_ratios(results_cortisol_sugar)

hetero_cortisol_sugar<- mr_heterogeneity(dat_cortisol_sugar, method_list=c("mr_ivw_mre", "mr_egger_regression"))
hetero_cortisol_sugar$analysis<- "cortisol_sugar"


## Check intercept of Egger regression (only  sensitivity-as cortisol main has only 1 SNP)
intercept_cortisol_sugar<- mr_pleiotropy_test(dat_cortisol_sugar)
intercept_cortisol_sugar$analysis<- "cortisol_sugar"


## MR presso (only sensitivity-as cortisol main has only 1 SNP)
presso_cortisol_sugar<-mr_presso(BetaOutcome = "beta.outcome", 
                  BetaExposure = "beta.exposure", 
                  SdOutcome = "se.outcome", 
                  SdExposure = "se.exposure", 
                  OUTLIERtest = TRUE,
                  DISTORTIONtest = TRUE, 
                  data = dat_cortisol_sugar, 
                  NbDistribution = 1000,  
                  SignifThreshold = 0.05)

presso_cortisol_sugar[["Main MR results"]]
presso_cortisol_sugar$`MR-PRESSO results`$`Global Test`$RSSobs
presso_cortisol_sugar$`MR-PRESSO results`$`Global Test`$Pvalue

presso_cortisol_sugar<- presso_cortisol_sugar[["Main MR results"]]
presso_cortisol_sugar$globaltest<-c(presso_cortisol_sugar$`MR-PRESSO results`$`Global Test`$RSSobs,
                                        presso_cortisol_sugar$`MR-PRESSO results`$`Global Test`$Pvalue)
presso_cortisol_sugar$analysis<- "cortisol_sugar"


##Single SNP analysis (only sensitivity-as cortisol main has only 1 SNP)
res_single_cortisol_sugar <- mr_singlesnp(dat_cortisol_sugar)
res_single_cortisol_sugar_full <- mr_singlesnp(dat_cortisol_sugar, all_method = "mr_two_sample_ml")

res_single_cortisol_sugar <- generate_odds_ratios(res_single_cortisol_sugar)
res_single_cortisol_sugar_full <-generate_odds_ratios(res_single_cortisol_sugar_full)


##Leave-one-out analysis (only sensitivity-as cortisol main has only 1 SNP)

res_loo_cortisol_sugar <- mr_leaveoneout(dat_cortisol_sugar)
res_loo_cortisol_sugar <- generate_odds_ratios(res_loo_cortisol_sugar)


### 6.Graphs (only sensitivity-as cortisol main has only 1 SNP)


##Single SNP 
c_s_single <- mr_forest_plot(res_single_cortisol_sugar)
single_cs_graph <- c_s_single[[1]]

c_s_single_full <- mr_forest_plot(res_single_cortisol_sugar_full)
c_s_single_full[[1]]

##Leave one out 
c_s_loo <- mr_leaveoneout_plot(res_loo_cortisol_sugar)
loo_cs_graph <- c_s_loo[[1]]



#######################################################################################################
#######################Cortisol-exposure, fat-outcome##################################################

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


###Getting cortisol data 
#Main
exp_cortisol_2<- extract_instruments("ieu-a-1012", p1=5e-8)
exp_cortisol_2_clu<- clump_data(exp_cortisol_2)

#Sensitivity
exp_cortisol<- extract_instruments("ieu-a-1012", p1=5e-6)
exp_cortisol_clu<- clump_data(exp_cortisol)



###Prepare outcome data
##Cortisol main and fat
out_cortisol_2_fat <- read_outcome_data(
  snps = exp_cortisol_2_clu$SNP,
  filename = "fat.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "pval",
  eaf_col = "eaf"
)

nrow(out_cortisol_2_fat)
nrow(exp_cortisol_2_clu)
out_cortisol_2_fat_proxy <- out_cortisol_2_fat

#Cortisol sensitivity and fat
out_cortisol_fat <- read_outcome_data(
  snps = exp_cortisol_clu$SNP,
  filename = "fat.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "pval",
  eaf_col = "eaf"
)

nrow(out_cortisol_fat)
nrow(exp_cortisol_clu)
out_cortisol_fat_proxy <- out_cortisol_fat



### 3. Harmonise data 

dat_cortisol_2_fat <- harmonise_data(
  exposure_dat = exp_cortisol_2_clu, 
  outcome_dat = out_cortisol_2_fat_proxy, 
  action =2)

dat_cortisol_fat <- harmonise_data(
  exposure_dat = exp_cortisol_clu, 
  outcome_dat = out_cortisol_fat_proxy, 
  action =2)


### 4. Drop duplicate exposure-outcome summary sets
## Pruning 

dat_cortisol_2_fat$samplesize.outcome<- 264181
dat_cortisol_2_fat<-power_prune(dat_cortisol_2_fat,method=2,dist.outcome="continuous")


dat_cortisol_fat$samplesize.outcome<- 264181
dat_cortisol_fat<-power_prune(dat_cortisol_fat,method=2,dist.outcome="continuous")


## Steiger directionality test 

steiger_filter_cortisol_2_fat<-steiger_filtering(dat_cortisol_2_fat)
steiger_direct_cortisol_2_fat <- directionality_test(dat_cortisol_2_fat)


steiger_filter_cortisol_fat<-steiger_filtering(dat_cortisol_fat)
steiger_direct_cortisol_fat <- directionality_test(dat_cortisol_fat)


## Edit names of outcome and exposure columns 

dat_cortisol_2_fat$exposure <- "Plasma cortisol levels (p<10e-8)"
dat_cortisol_2_fat$outcome <- "Relative fat intake"


dat_cortisol_fat$exposure <- "Plasma cortisol levels (p<10e-6)"
dat_cortisol_fat$outcome<- "Relative fat intake"


### 5. MR analyses 

## Ivw, egger, weighted_median, weighted mode 

results_cortisol_2_fat<-mr(dat_cortisol_2_fat, method_list=c("mr_wald_ratio"))


results_cortisol_fat<-mr(dat_cortisol_fat, method_list=c("mr_ivw_mre",
                                                         "mr_egger_regression",
                                                         "mr_weighted_median", 
                                                         "mr_weighted_mode"))


results_cortisol_2_fat
results_cortisol_fat


###Sensitivity analysis
## Raps 
#Main cortisol 
raps_cortisol_2_fat<- mr.raps(dat_cortisol_2_fat$beta.exposure, dat_cortisol_2_fat$beta.outcome, dat_cortisol_2_fat$se.exposure, 
                              dat_cortisol_2_fat$se.outcome)

raps_df_cortisol_2_fat<- data.frame(id.exposure = results_cortisol_2_fat[1,1],
                                    id.outcome = results_cortisol_2_fat[1,2], 
                                    outcome = results_cortisol_2_fat[1,3],
                                    exposure = results_cortisol_2_fat[1,4], 
                                    method = c("mr_raps"), 
                                    nsnp = results_cortisol_2_fat[1,6],
                                    b = raps_cortisol_2_fat[["beta.hat"]], 
                                    se = raps_cortisol_2_fat[["beta.se"]], 
                                    pval = raps_cortisol_2_fat[["beta.p.value"]])

results_cortisol_2_fat<-rbind(results_cortisol_2_fat,raps_df_cortisol_2_fat )
results_cortisol_2_fat$analysis = "cortisol_2_fat"
results_cortisol_2_fat <- generate_odds_ratios(results_cortisol_2_fat)

#Sensitivity cortisol
raps_cortisol_fat<- mr.raps(dat_cortisol_fat$beta.exposure, dat_cortisol_fat$beta.outcome, dat_cortisol_fat$se.exposure, 
                            dat_cortisol_fat$se.outcome)

raps_df_cortisol_fat<- data.frame(id.exposure = results_cortisol_fat[1,1],
                                  id.outcome = results_cortisol_fat[1,2], 
                                  outcome = results_cortisol_fat[1,3],
                                  exposure = results_cortisol_fat[1,4], 
                                  method = c("mr_raps"), 
                                  nsnp = results_cortisol_fat[1,6],
                                  b = raps_cortisol_fat[["beta.hat"]], 
                                  se = raps_cortisol_fat[["beta.se"]], 
                                  pval = raps_cortisol_fat[["beta.p.value"]])

results_cortisol_fat<-rbind(results_cortisol_fat,raps_df_cortisol_fat )
results_cortisol_fat$analysis = "cortisol_fat"
results_cortisol_fat <- generate_odds_ratios(results_cortisol_fat)

hetero_cortisol_fat<- mr_heterogeneity(dat_cortisol_fat, method_list=c("mr_ivw_mre", "mr_egger_regression"))
hetero_cortisol_fat$analysis<- "cortisol_fat"


## Check intercept of Egger regression (only sensitivity-as cortisol main has only 1 SNP)

intercept_cortisol_fat<- mr_pleiotropy_test(dat_cortisol_fat)
intercept_cortisol_fat$analysis<- "cortisol_fat"

## MR presso (only sensitivity-as cortisol main has only 1 SNP)
presso_cortisol_fat<-mr_presso(BetaOutcome = "beta.outcome", 
                               BetaExposure = "beta.exposure", 
                               SdOutcome = "se.outcome", 
                               SdExposure = "se.exposure", 
                               OUTLIERtest = TRUE,
                               DISTORTIONtest = TRUE, 
                               data = dat_cortisol_fat, 
                               NbDistribution = 1000,  
                               SignifThreshold = 0.05)

presso_cortisol_fat[["Main MR results"]]
presso_cortisol_fat$`MR-PRESSO results`$`Global Test`$RSSobs
presso_cortisol_fat$`MR-PRESSO results`$`Global Test`$Pvalue

presso_cortisol_fat<- presso_cortisol_fat[["Main MR results"]]
presso_cortisol_fat$globaltest<-c(presso_cortisol_fat$`MR-PRESSO results`$`Global Test`$RSSobs,
                                  presso_cortisol_fat$`MR-PRESSO results`$`Global Test`$Pvalue)
presso_cortisol_fat$analysis<- "cortisol_fat"


##Single SNP analysis (only sensitivity-as cortisol main has only 1 SNP)
res_single_cortisol_fat <- mr_singlesnp(dat_cortisol_fat)
res_single_cortisol_fat_full <- mr_singlesnp(dat_cortisol_fat, all_method = "mr_two_sample_ml")

res_single_cortisol_fat <- generate_odds_ratios(res_single_cortisol_fat)
res_single_cortisol_fat_full <-generate_odds_ratios(res_single_cortisol_fat_full)


##Leave-one-out analysis (only sensitivity-as cortisol main has only 1 SNP)

res_loo_cortisol_fat <- mr_leaveoneout(dat_cortisol_fat)
res_loo_cortisol_fat <- generate_odds_ratios(res_loo_cortisol_fat)


### 6.Graphs (only sensitivity-as cortisol main has only 1 SNP)


##Single SNP
c_f_single <- mr_forest_plot(res_single_cortisol_fat)
single_cf_graph <- c_f_single[[1]]

c_f_single_full <- mr_forest_plot(res_single_cortisol_fat_full)
c_f_single_full[[1]]

##Leave one out
c_f_loo <- mr_leaveoneout_plot(res_loo_cortisol_fat)
loo_sf_graph <- c_f_loo[[1]]


#######################################################################################################
#######################Depression-exposure, sugar-outcome##################################################


# Getting sugar data
sugar_df<-read.csv("sugar.csv")
names(sugar_df)


###Getting depression data 
#Main
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
exp_dep <- format_data(dep_df, type="exposure",samplesize_col = "samplesize.exposure")

#Main depression
exp_dep_2<- exp_dep[exp_dep$pval.exposure <=5e-8, ]
exp_dep_2_clu<- clump_data(exp_dep_2)

#Sensitivity
exp_dep<- exp_dep[exp_dep$pval.exposure <=5e-6, ]
exp_dep_clu<- clump_data(exp_dep)



###Prepare outcome data
##Depression main and sugar
out_dep_2_sugar <- read_outcome_data(
  snps = exp_dep_2_clu$SNP,
  filename = "sugar.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "pval",
  eaf_col = "eaf"
)

nrow(out_dep_2_sugar)
nrow(exp_dep_2_clu)

## Find proxy SNPs in outcome data 
expSNP_absent_dep_2_sugar<- exp_dep_2_clu[!exp_dep_2_clu$SNP %in% out_dep_2_sugar$SNP, ]

# Loop 
df_2 = data.frame()

for(i in (expSNP_absent_dep_2_sugar$SNP)) {
  trial<-LDproxy(i, pop ="CEU",
                 r2d = "r2", token = "f44d87b4d416")
  df_2 = rbind(df, trial)  
}


# keep only r >=.80
df_2_high<- df_2[df_2$R2 >= 0.80,]

# select only rows for proxy SNPs 

sugar_proxy_2 <- filter(sugar_df, SNP %in% df_2_high$RS_Number) 
nrow(sugar_proxy_2)

# rename column names 


sugar_proxy_2<- sugar_proxy_2 %>%
  mutate(SNP = SNP,
         effect_allele.outcome = effect_allele,
         other_allele.outcome = other_allele,
         beta.outcome = beta,
         se.outcome =  se,
         pval.outcome = pval ,
         outcome = "outcome",
         mr_keep.outcome = "TRUE",
         pval_origin.outcome = "reported" ,
         id.outcome = out_dep_2_sugar[1,11] , 
         eaf.outcome = eaf, 
         data_source.outcome = "textfile") %>%
  select(SNP ,
         effect_allele.outcome,
         other_allele.outcome ,
         beta.outcome ,
         se.outcome ,
         pval.outcome ,
         outcome,
         mr_keep.outcome , 
         pval_origin.outcome , 
         id.outcome ,
         eaf.outcome,
         data_source.outcome)


sugar_proxy_2 <- sugar_proxy_2[!sugar_proxy_2$SNP %in% out_dep_2_sugar$SNP, ]


sugar_proxy_2$effect_allele.outcome[sugar_proxy_2$effect_allele.outcome == 't'] <- 'T'  
sugar_proxy_2$effect_allele.outcome[sugar_proxy_2$effect_allele.outcome == 'a'] <- 'A'  
sugar_proxy_2$effect_allele.outcome[sugar_proxy_2$effect_allele.outcome == 'c'] <- 'C'  
sugar_proxy_2$effect_allele.outcome[sugar_proxy_2$effect_allele.outcome == 'g'] <- 'G' 

sugar_proxy_2$other_allele.outcome[sugar_proxy_2$other_allele.outcome == 't'] <- 'T'  
sugar_proxy_2$other_allele.outcome[sugar_proxy_2$other_allele.outcome == 'a'] <- 'A'  
sugar_proxy_2$other_allele.outcome[sugar_proxy_2$other_allele.outcome == 'c'] <- 'C'  
sugar_proxy_2$other_allele.outcome[sugar_proxy_2$other_allele.outcome == 'g'] <- 'G'  


# combine with outcome data
out_dep_2_sugar_proxy<- rbind(out_dep_2_sugar,sugar_proxy_2)

#Depression sensitivity and sugar
out_dep_sugar <- read_outcome_data(
  snps = exp_dep_clu$SNP,
  filename = "sugar.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "pval",
  eaf_col = "eaf"
)

nrow(out_dep_sugar)
nrow(exp_dep_clu)

## Find proxy SNPs in outcome data 

# Exposure SNPs not available in outcome 
expSNP_absent_dep_sugar<- exp_dep_clu[!exp_dep_clu$SNP %in% out_dep_sugar$SNP, ]

# Loop 
df = data.frame()

for(i in (expSNP_absent_dep_sugar$SNP)) {
  trial<-LDproxy(i, pop ="CEU",
                 r2d = "r2", token = "f44d87b4d416")
  df = rbind(df, trial)  
}



# keep only r >=.80
df_high<- df[df$R2 >= 0.80,]

# select only rows for proxy SNPs 

sugar_proxy <- filter(sugar_df, SNP %in% df_high$RS_Number) 
nrow(sugar_proxy)
colnames(sugar_proxy)
# rename column names 


sugar_proxy<- sugar_proxy %>%
  mutate(SNP = SNP,
         effect_allele.outcome = effect_allele,
         other_allele.outcome = other_allele,
         beta.outcome = beta,
         se.outcome =  se,
         pval.outcome = pval ,
         outcome = "outcome",
         mr_keep.outcome = "TRUE",
         pval_origin.outcome = "reported" ,
         eaf.outcome = eaf,
         id.outcome = out_dep_sugar[1,11] , 
         data_source.outcome = "textfile") %>%
  select(SNP ,
         effect_allele.outcome,
         other_allele.outcome ,
         beta.outcome ,
         se.outcome ,
         pval.outcome ,
         outcome,
         mr_keep.outcome , 
         pval_origin.outcome , 
         eaf.outcome,
         id.outcome ,
         data_source.outcome)

#Remove duplicated
sugar_proxy <- sugar_proxy[!sugar_proxy$SNP %in% out_dep_sugar$SNP, ]


sugar_proxy$effect_allele.outcome[sugar_proxy$effect_allele.outcome == 't'] <- 'T'  
sugar_proxy$effect_allele.outcome[sugar_proxy$effect_allele.outcome == 'a'] <- 'A'  
sugar_proxy$effect_allele.outcome[sugar_proxy$effect_allele.outcome == 'c'] <- 'C'  
sugar_proxy$effect_allele.outcome[sugar_proxy$effect_allele.outcome == 'g'] <- 'G' 

sugar_proxy$other_allele.outcome[sugar_proxy$other_allele.outcome == 't'] <- 'T'  
sugar_proxy$other_allele.outcome[sugar_proxy$other_allele.outcome == 'a'] <- 'A'  
sugar_proxy$other_allele.outcome[sugar_proxy$other_allele.outcome == 'c'] <- 'C'  
sugar_proxy$other_allele.outcome[sugar_proxy$other_allele.outcome == 'g'] <- 'G'  


colnames(out_dep_sugar)
colnames(sugar_proxy)

# combine with outcome data
out_dep_sugar_proxy<- rbind(out_dep_sugar,sugar_proxy)




### 3. Harmonise data 
dat_dep_2_sugar <- harmonise_data(
  exposure_dat = exp_dep_2_clu, 
  outcome_dat = out_dep_2_sugar_proxy, 
  action =2)


dat_dep_sugar <- harmonise_data(
  exposure_dat = exp_dep_clu, 
  outcome_dat = out_dep_sugar_proxy, 
  action =2)



### 4. Drop duplicate exposure-outcome summary sets
## Pruning 
dat_dep_2_sugar$samplesize.outcome<- 230648 
dat_dep_2_sugar$samplesize.exposure<- 138884 
dat_dep_2_sugar<-power_prune(dat_dep_2_sugar,method=2,dist.outcome="continuous")

dat_dep_sugar$samplesize.outcome<- 230648 
dat_dep_sugar$samplesize.exposure<- 138884
dat_dep_sugar<-power_prune(dat_dep_sugar,method=2,dist.outcome="continuous")


## Steiger directionality test 
steiger_filter_dep_2_sugar<-steiger_filtering(dat_dep_2_sugar)
steiger_direct_dep_2_sugar <- directionality_test(dat_dep_2_sugar)

steiger_filter_dep_sugar<-steiger_filtering(dat_dep_sugar)
steiger_direct_dep_sugar<- directionality_test(dat_dep_sugar)


## Edit names of outcome and exposure columns 

dat_dep_2_sugar$exposure <- "MDD (p<10e-8)"
dat_dep_2_sugar$outcome <- "Relative sugar intake"


dat_dep_sugar$exposure <- "MDD (p<10e-6)"
dat_dep_sugar$outcome<- "Relative sugar intake"


### 5. MR analyses 

## Ivw, egger, weighted_median, weighted mode 

results_dep_2_sugar<-mr(dat_dep_2_sugar, method_list=c("mr_ivw_mre"))

results_dep_sugar<-mr(dat_dep_sugar, method_list=c("mr_ivw_mre",
                                                   "mr_egger_regression",
                                                   "mr_weighted_median", 
                                                   "mr_weighted_mode"))


results_dep_2_sugar
results_dep_sugar


###Sensitivity analysis
## Raps 
#Main depression 
raps_dep_2_sugar<- mr.raps(dat_dep_2_sugar$beta.exposure, dat_dep_2_sugar$beta.outcome, dat_dep_2_sugar$se.exposure, 
                           dat_dep_2_sugar$se.outcome)

raps_df_dep_2_sugar<- data.frame(id.exposure = results_dep_2_sugar[1,1],
                                 id.outcome = results_dep_2_sugar[1,2], 
                                 outcome = results_dep_2_sugar[1,3],
                                 exposure = results_dep_2_sugar[1,4], 
                                 method = c("mr_raps"), 
                                 nsnp = results_dep_2_sugar[1,6],
                                 b = raps_dep_2_sugar[["beta.hat"]], 
                                 se = raps_dep_2_sugar[["beta.se"]], 
                                 pval = raps_dep_2_sugar[["beta.p.value"]])

results_dep_2_sugar<-rbind(results_dep_2_sugar,raps_df_dep_2_sugar )
results_dep_2_sugar$analysis = "dep_2_sugar"
results_dep_2_sugar <- generate_odds_ratios(results_dep_2_sugar)

hetero_dep_2_sugar<- mr_heterogeneity(dat_dep_2_sugar, method_list=c("mr_ivw_mre", "mr_egger_regression"))
hetero_dep_2_sugar$analysis<- "dep_2_sugar"

#Sensitivity depression
raps_dep_sugar<- mr.raps(dat_dep_sugar$beta.exposure, dat_dep_sugar$beta.outcome, dat_dep_sugar$se.exposure, 
                         dat_dep_sugar$se.outcome)

raps_df_dep_sugar<- data.frame(id.exposure = results_dep_sugar[1,1],
                               id.outcome = results_dep_sugar[1,2], 
                               outcome = results_dep_sugar[1,3],
                               exposure = results_dep_sugar[1,4], 
                               method = c("mr_raps"), 
                               nsnp = results_dep_sugar[1,6],
                               b = raps_dep_sugar[["beta.hat"]], 
                               se = raps_dep_sugar[["beta.se"]], 
                               pval = raps_dep_sugar[["beta.p.value"]])

results_dep_sugar<-rbind(results_dep_sugar,raps_df_dep_sugar )
results_dep_sugar$analysis = "dep_sugar"
results_dep_sugar <- generate_odds_ratios(results_dep_sugar)

hetero_dep_sugar<- mr_heterogeneity(dat_dep_sugar, method_list=c("mr_ivw_mre", "mr_egger_regression"))
hetero_dep_sugar$analysis<- "dep_sugar"


## Check intercept of Egger regression (only sensitivity-as depression main has only 2 SNPs)
intercept_dep_sugar<- mr_pleiotropy_test(dat_dep_sugar)
intercept_dep_sugar$analysis<- "dep_sugar"


## MR presso (only sensitivity-as depression main has only 2 SNPs)
presso_dep_sugar<-mr_presso(BetaOutcome = "beta.outcome", 
                            BetaExposure = "beta.exposure", 
                            SdOutcome = "se.outcome", 
                            SdExposure = "se.exposure", 
                            OUTLIERtest = TRUE,
                            DISTORTIONtest = TRUE, 
                            data = dat_dep_sugar, 
                            NbDistribution = 1000,  
                            SignifThreshold = 0.05)

presso_dep_sugar[["Main MR results"]]
presso_dep_sugar$`MR-PRESSO results`$`Global Test`$RSSobs
presso_dep_sugar$`MR-PRESSO results`$`Global Test`$Pvalue

presso_dep_sugar<- presso_dep_sugar[["Main MR results"]]
presso_dep_sugar$globaltest<-c(presso_dep_sugar$`MR-PRESSO results`$`Global Test`$RSSobs,
                               presso_dep_sugar$`MR-PRESSO results`$`Global Test`$Pvalue)
presso_dep_sugar$analysis<- "dep_sugar"


##Single SNP analysis (only sensitivity-as depression main has only 2 SNPs)
res_single_dep_sugar <- mr_singlesnp(dat_dep_sugar)
res_single_dep_sugar_full <- mr_singlesnp(dat_dep_sugar, all_method = "mr_two_sample_ml")

res_single_dep_sugar <- generate_odds_ratios(res_single_dep_sugar)
res_single_dep_sugar_full <-generate_odds_ratios(res_single_dep_sugar_full)

##Leave-one-out analysis (only sensitivity-as depression main has only 2 SNPs)
res_loo_dep_sugar <- mr_leaveoneout(dat_dep_sugar)
res_loo_dep_sugar <- generate_odds_ratios(res_loo_dep_sugar)


### 6.Graphs (only sensitivity-as depression main has only 2 SNPs)

##Single SNP
d_s_single <- mr_forest_plot(res_single_dep_sugar)
single_ds_graph <- d_s_single[[1]]

d_s_single_full <- mr_forest_plot(res_single_dep_sugar_full)
d_s_single_full[[1]]

##Leave one out
d_s_loo <- mr_leaveoneout_plot(res_loo_dep_sugar)
loo_ds_graph <- d_s_loo[[1]]


#######################################################################################################
#######################Depression-exposure, fat-outcome##################################################


# Getting fat data
fat_df<-read.csv("fat.csv")
names(fat_df)


###Getting depression data 
#Main
dep_df<-read.csv("dep.csv")
names(dep_df)
exp_dep <- format_data(dep_df, type="exposure",samplesize_col = "samplesize.exposure")

#Main depression
exp_dep_2<- exp_dep[exp_dep$pval.exposure <=5e-8, ]
exp_dep_2_clu<- clump_data(exp_dep_2)

#Sensitivity depression
exp_dep<- exp_dep[exp_dep$pval.exposure <=5e-6, ]
exp_dep_clu<- clump_data(exp_dep)


###Prepare outcome data
##Depression main and fat
out_dep_2_fat <- read_outcome_data(
  snps = exp_dep_2_clu$SNP,
  filename = "fat.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "pval",
  eaf_col = "eaf"  
)

nrow(out_dep_2_fat)
nrow(exp_dep_2_clu)

# Exposure SNPs not available in outcome 
expSNP_absent_dep_2_fat<- exp_dep_2_clu[!exp_dep_2_clu$SNP %in% out_dep_2_fat$SNP, ]

# Loop 
df_2 = data.frame()

for(i in (expSNP_absent_dep_2_fat$SNP)) {
  trial<-LDproxy(i, pop ="CEU",
                 r2d = "r2", token = "f44d87b4d416")
  df_2 = rbind(df_2, trial)  
}



# keep only r >=.80
df_2_high<- df_2[df_2$R2 >= 0.80,]

# select only rows for proxy SNPs 

fat_proxy_2 <- filter(fat_df, SNP %in% df_2_high$RS_Number) 

# rename column names 

fat_proxy_2<- fat_proxy_2 %>%
  mutate(SNP = SNP,
         effect_allele.outcome = effect_allele,
         other_allele.outcome = other_allele,
         beta.outcome = beta,
         se.outcome =  se,
         pval.outcome = pval ,
         outcome = "outcome",
         mr_keep.outcome = "TRUE",
         pval_origin.outcome = "reported" ,
         id.outcome = out_dep_2_fat[1,11] , 
         eaf.outcome = eaf,
         data_source.outcome = "textfile") %>%
  select(SNP ,
         effect_allele.outcome,
         other_allele.outcome ,
         beta.outcome ,
         se.outcome ,
         pval.outcome ,
         outcome,
         mr_keep.outcome , 
         pval_origin.outcome , 
         id.outcome ,
         eaf.outcome,
         data_source.outcome)

#Remove duplicated
fat_proxy_2 <- fat_proxy_2[!fat_proxy_2$SNP %in% out_dep_2_fat$SNP, ]


fat_proxy_2$effect_allele.outcome[fat_proxy_2$effect_allele.outcome == 't'] <- 'T'  
fat_proxy_2$effect_allele.outcome[fat_proxy_2$effect_allele.outcome == 'a'] <- 'A'  
fat_proxy_2$effect_allele.outcome[fat_proxy_2$effect_allele.outcome == 'c'] <- 'C' 
fat_proxy_2$effect_allele.outcome[fat_proxy_2$effect_allele.outcome == 'g'] <- 'G' 

fat_proxy_2$other_allele.outcome[fat_proxy_2$other_allele.outcome == 't'] <- 'T'  
fat_proxy_2$other_allele.outcome[fat_proxy_2$other_allele.outcome == 'a'] <- 'A'  
fat_proxy_2$other_allele.outcome[fat_proxy_2$other_allele.outcome == 'c'] <- 'C'  
fat_proxy_2$other_allele.outcome[fat_proxy_2$other_allele.outcome == 'g'] <- 'G'  


# combine with outcome data
out_dep_2_fat_proxy<- rbind(out_dep_2_fat,fat_proxy_2)
colnames(out_dep_2_fat)
colnames(fat_proxy_2)

##Depression sensitivity and fat
out_dep_fat <- read_outcome_data(
  snps = exp_dep_clu$SNP,
  filename = "fat.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "pval",
  eaf_col = "eaf"
)

nrow(out_dep_fat)
nrow(exp_dep_clu)

## Find proxy SNPs in outcome data 

# Exposure SNPs not available in outcome 
expSNP_absent_dep_fat<- exp_dep_clu[!exp_dep_clu$SNP %in% out_dep_fat$SNP, ]

# Loop 
df = data.frame()

for(i in (expSNP_absent_dep_fat$SNP)) {
  trial<-LDproxy(i, pop ="CEU",
                 r2d = "r2", token = "f44d87b4d416")
  df = rbind(df, trial)  
}



# keep only r >=.80
df_high<- df[df$R2 >= 0.80,]

# select only rows for proxy SNPs 

fat_proxy <- filter(fat_df, SNP %in% df_high$RS_Number) 

# rename column names 

fat_proxy<- fat_proxy %>%
  mutate(SNP = SNP,
         effect_allele.outcome = effect_allele,
         other_allele.outcome = other_allele,
         beta.outcome = beta,
         se.outcome =  se,
         pval.outcome = pval ,
         outcome = "outcome",
         mr_keep.outcome = "TRUE",
         pval_origin.outcome = "reported" ,
         id.outcome = out_dep_fat[1,11] , 
         eaf.outcome = eaf,
         data_source.outcome = "textfile") %>%
  select(SNP ,
         effect_allele.outcome,
         other_allele.outcome ,
         beta.outcome ,
         se.outcome ,
         pval.outcome ,
         outcome,
         mr_keep.outcome , 
         pval_origin.outcome , 
         id.outcome ,
         eaf.outcome,
         data_source.outcome)

#Remove duplicated
fat_proxy <- fat_proxy[!fat_proxy$SNP %in% out_dep_fat$SNP, ]


fat_proxy$effect_allele.outcome[fat_proxy$effect_allele.outcome == 't'] <- 'T'  
fat_proxy$effect_allele.outcome[fat_proxy$effect_allele.outcome == 'a'] <- 'A'  
fat_proxy$effect_allele.outcome[fat_proxy$effect_allele.outcome == 'c'] <- 'C'  
fat_proxy$effect_allele.outcome[fat_proxy$effect_allele.outcome == 'g'] <- 'G' 

fat_proxy$other_allele.outcome[fat_proxy$other_allele.outcome == 't'] <- 'T'  
fat_proxy$other_allele.outcome[fat_proxy$other_allele.outcome == 'a'] <- 'A'  
fat_proxy$other_allele.outcome[fat_proxy$other_allele.outcome == 'c'] <- 'C'  
fat_proxy$other_allele.outcome[fat_proxy$other_allele.outcome == 'g'] <- 'G'  


colnames(out_dep_fat)
colnames(fat_proxy)
colnames(exp_dep_clu)
exp_dep_clu$data_source.exposure <- "textfile"

# combine with outcome data
out_dep_fat_proxy<- rbind(out_dep_fat,fat_proxy)




### 3. Harmonise data 
dat_dep_2_fat <- harmonise_data(
  exposure_dat = exp_dep_2_clu, 
  outcome_dat = out_dep_2_fat_proxy, 
  action =2)


dat_dep_fat <- harmonise_data(
  exposure_dat = exp_dep_clu, 
  outcome_dat = out_dep_fat_proxy, 
  action =2)



### 4. Drop duplicate exposure-outcome summary sets
## Pruning 
dat_dep_2_fat$samplesize.outcome<- 264181
dat_dep_2_fat$samplesize.exposure<- 138884
dat_dep_2_fat<-power_prune(dat_dep_2_fat,method=2,dist.outcome="continuous")


dat_dep_fat$samplesize.outcome<- 264181
dat_dep_fat$samplesize.exposure<- 138884
dat_dep_fat<-power_prune(dat_dep_fat,method=2,dist.outcome="continuous")


## Steiger directionality test 
steiger_filter_dep_2_fat<-steiger_filtering(dat_dep_2_fat)
steiger_direct_dep_2_fat <- directionality_test(dat_dep_2_fat)

steiger_filter_dep_fat<-steiger_filtering(dat_dep_fat)
steiger_direct_dep_fat <- directionality_test(dat_dep_fat)


## Edit names of outcome and exposure columns 
dat_dep_2_fat$exposure <- "MDD (p<10e-8)"
dat_dep_2_fat$outcome <- "Relative fat intake"


dat_dep_fat$exposure <- "MDD (p<10e-6)"
dat_dep_fat$outcome<- "Relative fat intake"


### 5. MR analyses 

## Ivw, egger, weighted_median, weighted mode 

results_dep_2_fat<-mr(dat_dep_2_fat, method_list=c("mr_ivw_mre"))


results_dep_fat<-mr(dat_dep_fat, method_list=c("mr_ivw_mre",
                                               "mr_egger_regression",
                                               "mr_weighted_median", 
                                               "mr_weighted_mode"))


results_dep_2_fat
results_dep_fat


###Sensitivity analysis
## Raps 
#Main depression 
raps_dep_2_fat<- mr.raps(dat_dep_2_fat$beta.exposure, dat_dep_2_fat$beta.outcome, dat_dep_2_fat$se.exposure, 
                         dat_dep_2_fat$se.outcome)

raps_df_dep_2_fat<- data.frame(id.exposure = results_dep_2_fat[1,1],
                               id.outcome = results_dep_2_fat[1,2], 
                               outcome = results_dep_2_fat[1,3],
                               exposure = results_dep_2_fat[1,4], 
                               method = c("mr_raps"), 
                               nsnp = results_dep_2_fat[1,6],
                               b = raps_dep_2_fat[["beta.hat"]], 
                               se = raps_dep_2_fat[["beta.se"]], 
                               pval = raps_dep_2_fat[["beta.p.value"]])

results_dep_2_fat<-rbind(results_dep_2_fat,raps_df_dep_2_fat )
results_dep_2_fat$analysis = "dep_2_fat"
results_dep_2_fat <- generate_odds_ratios(results_dep_2_fat)

hetero_dep_2_fat<- mr_heterogeneity(dat_dep_2_fat, method_list=c("mr_ivw_mre", "mr_egger_regression"))
hetero_dep_2_fat$analysis<- "dep_2_fat"

#Senstivity depression
raps_dep_fat<- mr.raps(dat_dep_fat$beta.exposure, dat_dep_fat$beta.outcome, dat_dep_fat$se.exposure, 
                       dat_dep_fat$se.outcome)

raps_df_dep_fat<- data.frame(id.exposure = results_dep_fat[1,1],
                             id.outcome = results_dep_fat[1,2], 
                             outcome = results_dep_fat[1,3],
                             exposure = results_dep_fat[1,4], 
                             method = c("mr_raps"), 
                             nsnp = results_dep_fat[1,6],
                             b = raps_dep_fat[["beta.hat"]], 
                             se = raps_dep_fat[["beta.se"]], 
                             pval = raps_dep_fat[["beta.p.value"]])

results_dep_fat<-rbind(results_dep_fat,raps_df_dep_fat )
results_dep_fat$analysis = "dep_fat"
results_dep_fat <- generate_odds_ratios(results_dep_fat)

hetero_dep_fat<- mr_heterogeneity(dat_dep_fat, method_list=c("mr_ivw_mre", "mr_egger_regression"))
hetero_dep_fat$analysis<- "dep_fat"


## Check intercept of Egger regression (only sensitivity-as depression main has only 2 SNPs)
intercept_dep_fat<- mr_pleiotropy_test(dat_dep_fat)
intercept_dep_fat$analysis<- "dep_fat"


## MR presso (only sensitivity-as depression main has only 2 SNPs)
presso_dep_fat<-mr_presso(BetaOutcome = "beta.outcome", 
                          BetaExposure = "beta.exposure", 
                          SdOutcome = "se.outcome", 
                          SdExposure = "se.exposure", 
                          OUTLIERtest = TRUE,
                          DISTORTIONtest = TRUE, 
                          data = dat_dep_fat, 
                          NbDistribution = 1000,  
                          SignifThreshold = 0.05)

presso_dep_fat[["Main MR results"]]
presso_dep_fat$`MR-PRESSO results`$`Global Test`$RSSobs
presso_dep_fat$`MR-PRESSO results`$`Global Test`$Pvalue

presso_dep_fat<- presso_dep_fat[["Main MR results"]]
presso_dep_fat$globaltest<-c(presso_dep_fat$`MR-PRESSO results`$`Global Test`$RSSobs,
                             presso_dep_fat$`MR-PRESSO results`$`Global Test`$Pvalue)
presso_dep_fat$analysis<- "dep_fat"




##Single SNP analysis (only sensitivity-as depression main has only 2 SNPs)
res_single_dep_fat <- mr_singlesnp(dat_dep_fat)
res_single_dep_fat_full <- mr_singlesnp(dat_dep_fat, all_method = "mr_two_sample_ml")

res_single_dep_fat <- generate_odds_ratios(res_single_dep_fat)
res_single_dep_fat_full <-generate_odds_ratios(res_single_dep_fat_full)


##Leave-one-out analysis (only sensitivity-as depression main has only 2 SNPs)
res_loo_dep_fat <- mr_leaveoneout(dat_dep_fat)
res_loo_dep_fat <- generate_odds_ratios(res_loo_dep_fat)


### 6.Graphs (only sensitivity-as depression main has only 2 SNPs)

##Single SNP
d_f_single <- mr_forest_plot(res_single_dep_fat)
signle_df_graph <- d_f_single[[1]]

d_f_single_full <- mr_forest_plot(res_single_dep_fat_full)
d_f_single_full[[1]]

##Leave one out
d_f_loo <- mr_leaveoneout_plot(res_loo_dep_fat)
loo_ds_graph <- d_f_loo[[1]]

###################Figures#########################
single_sensitivity_1<-ggarrange(
  single_df_graph, single_ds_graph, single_cf_graph, single_cs_graph, labels = c("A", "B", "C", "D"))      

loo_sensitivity_1 <- ggarrange(
  loo_df_graph, loo_ds_graph, loo_cf_graph, loo_cs_graph, labels = c("A", "B", "C", "D"))   


#################################################################################################
####################################Combined results#############################################
results_macro_outcomes<- rbind(results_cortisol_fat,
                                results_cortisol_2_fat,
                                results_cortisol_sugar,
                                results_cortisol_2_sugar,
                                results_dep_fat,
                                results_dep_2_fat,
                                results_dep_sugar,
                                results_dep_2_sugar)

hetero_macro_outcomes<-rbind(hetero_cortisol_fat,
                              hetero_cortisol_sugar,
                             hetero_dep_2_fat,
                              hetero_dep_fat,
                              hetero_dep_2_sugar,
                              hetero_dep_sugar)

intercept_macro_outcomes<- rbind(intercept_cortisol_fat,
                                  intercept_cortisol_sugar,
                                  intercept_dep_fat,
                                  intercept_dep_sugar)

steiger_macro_outcomes<-rbind(steiger_direct_cortisol_fat,
                               steiger_direct_cortisol_2_fat,
                               steiger_direct_cortisol_sugar,
                               steiger_direct_cortisol_2_sugar,
                               steiger_direct_dep_fat,
                               steiger_direct_dep_2_fat,
                               steiger_direct_dep_sugar,
                               steiger_direct_dep_2_sugar)

presso_macro_outcomes<- rbind(presso_cortisol_fat,
                               presso_cortisol_sugar,
                               presso_dep_fat,
                              presso_dep_sugar)