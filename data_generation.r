#!/usr/bin/Rscript

#setwd("@TODO: set path")

#
#
##################################### symptoms ##################################### 
#
#
## Load symptoms-DO.tsv Data
symptoms_DO <- read.table("2_behavior_medical_attribute_values\\Symptoms\\symptoms-DO.tsv", fill = TRUE, sep = '\t', header = TRUE)

#
#
##################################### drugs ##################################### 
#
#
## Load drugs.tsv Data
zz=gzfile('2_behavior_medical_attribute_values\\Drugs\\drugs.tsv.gz','rt')  
drugs <- read.table(zz,sep = "\t", fill = TRUE, header=TRUE, na.strings=c(""))

#
#
##################################### SNPs ##################################### 
#
#
## Load snps.tsv Data
zz=gzfile('2_behavior_medical_attribute_values\\SNPs\\snps.tsv.gz','rt')  
snps <- read.table(zz,sep = "\t", fill = TRUE, header=TRUE, na.strings=c(""))

gens<- snps

#
#
##################################### User details ##################################### 
#
#
#user_data <- read.csv("data\\user_data.csv",  sep = ",", header=TRUE, na.strings=c(""))

zz=gzfile('1_CRM_attribute_values\\user_data.csv.gz','rt')  
dat=read.csv(zz,sep = ",", header=TRUE, na.strings=c(""))#, nrows=100000)

relevant_columns <-  c("GivenName","Surname","Gender","NameSet","Title","StreetAddress","City","State",
                       "StateFull","ZipCode","Country","EmailAddress","Username","Password","TelephoneNumber",
                       "TelephoneCountryCode","Birthday","Age","NationalID","Color","Occupation","Company",
                       "BloodType","Kilograms","Centimeters")
user_data <- dat[, relevant_columns]

#
#
##################################### Dependencies & Settings  ##################################### 
#
#
library(plyr)
library(foreach)
library(iterators)

if (!require("doParallel")) {
  install.packages('doParallel', dependencies=TRUE)
}
library(doParallel)

if (!require("rbenchmark")) {
  install.packages('rbenchmark', dependencies=TRUE)
}
library(rbenchmark)

if (!require("gtools")) {
  install.packages('gtools', dependencies=TRUE)
}
library(gtools)


## parrallel settings
no_cores <- detectCores() - 1 # Calculate the number of cores
registerDoParallel(no_cores)

## vars
max_disease=10
max_drug=10
max_snp=10
max_gen=10
#
#
##################################### Creating Attribute Distribution  ##################################### 
#
#

## Disease Distribution
disease_distribution <- rnorm(nrow(user_data), mean=(nrow(symptoms_DO)/2), sd = (nrow(symptoms_DO)/100)) # Create a sample of nrow(user_data) numbers which are normally distributed.
disease_distribution_rounded <- as.numeric(lapply(disease_distribution,round,0))
png(file = "disease_distribution.png") # Give the chart file a name.
hist(disease_distribution_rounded, main = "Disease Distribution") # Plot the histogram for this sample.
dev.off()

## Drug Distribution
drug_distribution <- rnorm(nrow(user_data), mean=(nrow(drugs)/2), sd = (nrow(drugs)/100)) # Create a sample of nrow(user_data) numbers which are normally distributed.
drug_distribution_rounded <- as.numeric(lapply(drug_distribution,round,0))
png(file = "drug_distribution.png") # Give the chart file a name.
hist(drug_distribution_rounded, main = "Drug Distribution") # Plot the histogram for this sample.
dev.off()

## SNPs Distribution
snp_distribution <- rnorm(nrow(user_data), mean=(nrow(snps)/2), sd = (nrow(snps)/100)) # Create a sample of nrow(user_data) numbers which are normally distributed.
snp_distribution_rounded <- as.numeric(lapply(snp_distribution,round,0))
png(file = "snp_distribution.png") # Give the chart file a name.
hist(snp_distribution_rounded, main = "SNPs Distribution") # Plot the histogram for this sample.
dev.off()

## Gen Distribution
gen_distribution <- rnorm(nrow(user_data), mean=(nrow(gens)/2), sd = (nrow(gens)/100)) # Create a sample of nrow(user_data) numbers which are normally distributed.
gen_distribution_rounded <- as.numeric(lapply(gen_distribution,round,0))
png(file = "gen_distribution.png") # Give the chart file a name.
hist(gen_distribution_rounded, main = "Gen Distribution") # Plot the histogram for this sample.
dev.off()

#
#
##################################### Dependencies & Settings  ##################################### 
#
#

## adding running ids
user_data$id  <-unlist(seq.int(nrow(user_data)))
symptoms_DO$id <- unlist(seq.int(nrow(symptoms_DO)))
drugs$id <- unlist(seq.int(nrow(drugs)))
snps$id <- unlist(seq.int(nrow(snps)))
gens$id <- unlist(seq.int(nrow(gens)))

#foreach(disease_counter=1:max_disease) %dopar% {
disease_counter=0
while(disease_counter < max_disease) {
  attribute_disease <- paste("disease", disease_counter, sep="_")
  
  ## choose randomly
  #symptoms_sample <- symptoms_DO[sample(nrow(symptoms_DO), nrow(user_data), replace = TRUE  ), ]$doid_cod
  #user_data[attribute_disease] <- sample_n(symptoms_DO, nrow(user_data))$doid_code
  
  ## choose with normal distribution
  user_data[attribute_disease] <- subset(symptoms_DO, disease_distribution_rounded %in% symptoms_DO$id)$doid_code
  
  attribute_disease_date <- paste("disease_date", disease_counter, sep="_")
  user_data[attribute_disease_date] <- sample(seq(as.Date("2000-04-14"), Sys.Date(), by="day"), nrow(user_data), replace = TRUE )
  disease_counter <- sum(disease_counter, 1)
}

#foreach(drug_counter=1:max_drug) %dopar% {
drug_counter=0
while(drug_counter < max_drug) {
  attribute_drug <- paste("drug", drug_counter, sep="_")
  
  ## choose randomly
  #user_data[attribute_drug] <- drugs[sample(nrow(drugs), nrow(user_data), replace = TRUE  ), ]$PRODUCTNDC
  
  ## choose with normal distribution
  user_data[attribute_drug] <- subset(drugs, drug_distribution_rounded %in% drugs$id )$PRODUCTNDC
  
  attribute_drug_date <- paste("drug_date", drug_counter, sep="_")
  user_data[attribute_drug_date] <- sample(seq(as.Date("2000-04-14"), Sys.Date(), by="day"), nrow(user_data), replace = TRUE )
  drug_counter <- sum(drug_counter, 1)
}


#foreach(drug_counter=1:max_snp) %dopar% {
snp_counter=0
while(snp_counter < max_snp) {
  attribute_snp <- paste("snp", snp_counter, sep="_")
  
  ## choose with normal distribution
  user_data[attribute_snp] <- subset(snps, snp_distribution_rounded %in% snps$id )$snpId
  snp_counter <- sum(snp_counter, 1)
}

#foreach(drug_counter=1:max_snp) %dopar% {
gen_counter=0
while(gen_counter < max_gen) {
  attribute_gen <- paste("gen", gen_counter, sep="_")
  
  ## choose with normal distribution
  user_data[attribute_gen] <- subset(gens, gen_distribution_rounded %in% gens$id )$geneId
  gen_counter <- sum(gen_counter, 1)
}

stopImplicitCluster()


#
#
##################################### Create User record by merging  ##################################### 
#
#  

## stora as csv.gz
gz1 <- gzfile("extended_user_data.csv.gz", "w")
write.csv(user_data, gz1)
close(gz1)

## stora as csv
#write.csv(user_data, file = "ex_user_data.csv")

summary(user_data)

gc()
