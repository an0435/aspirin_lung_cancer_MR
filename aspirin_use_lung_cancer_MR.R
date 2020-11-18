#ukb-a-132= aspirin use
#ieu-a-965=lung adenocarcinoma
#ieu-a-966=lung cancer
#ieu-a-967=Squamous cell lung cancer

library(devtools)
library(MRInstruments) 
library(TwoSampleMR)
library(shiny)
library(MRInstruments)
ao<-available_outcomes()
library(plyr)

##  normal analysis  ##
exposure_dat<-extract_instruments(outcomes="ukb-a-132")
outcome_dat<-extract_outcome_data(snps=exposure_dat$SNP, outcomes="ieu-a-967")
dat<-harmonise_data(exposure_dat, outcome_dat)
res<-mr(dat)
res
write.csv(dat, "ukb_a_132_aspirin_ieu_a_967_lung_squamous_dat.csv")
write.csv(res, "ukb_a_132_aspirin_ieu_a_967_lung_squamous_res.csv")

##  transformed analysis  ##
exposure_dat<-extract_instruments(outcomes="ukb-a-132")
#u= cases/(cases+controls)
u<-45012/(45012+292147)
exposure_dat<-rename(exposure_dat, c("beta.exposure"="beta_absolute"))
exposure_dat<-rename(exposure_dat, c("se.exposure"="se_absolute"))
exposure_dat$beta.exposure<-exposure_dat$beta_absolute/(u*(1-u))
exposure_dat$se.exposure<-exposure_dat$se_absolute/(u*(1-u))

outcome_dat<-extract_outcome_data(snps=exposure_dat$SNP, outcomes="ieu-a-967")
dat_OR<-harmonise_data(exposure_dat, outcome_dat)
res_OR<-mr(dat_OR)
res_OR
write.csv(dat_OR, "ukb_a_132_aspirin_ieu_a_967_lung_squamous_transformed_dat.csv")
write.csv(res_OR, "ukb_a_132_aspirin_ieu_a_967_lung_squamous_transformed_res.csv")
