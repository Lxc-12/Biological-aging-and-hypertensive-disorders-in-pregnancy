
# Pre-processing exposure data ----------------------------------
### convert the format of exposure data, for example, vcf to txt etc.
### Change the variable name of exposed data to a unified ones.
library(vcfR)
library(tidyr)
library(dplyr)
library(vroom)
#** Telomere length --------------------------------
path <- ("F:\\research\\aging_PE_GTN\\data\\exposure_gwas\\Telomere_length\\ieu-b-4879.vcf.gz")
TL_vcf <- read.vcfR(path, verbose = FALSE )

fix <- data.frame(TL_vcf@fix)
gt <- data.frame(TL_vcf@gt)

gt_sep <- separate(gt, "ieu.b.4879", into = c("ES","SE","LP","AF","rsID"),
                   sep = ":", remove = TRUE, convert = TRUE)

TL_merge <- data.frame(fix, gt_sep)

TL_merge$exposure <- "Telomere_length"
TL_merge$p_val <- 10^-TL_merge$LP

TL_gwas <- select(TL_merge, exposure,ID,CHROM,POS,ALT,REF,AF,ES,SE,p_val)

colnames(TL_gwas)[c(2:10)] <- c("SNP","chr","position","ea","oa","eaf","b","se","p")

write.table(TL_gwas, "F:\\research\\aging_PE_GTN\\data\\exposure_gwas\\Telomere_length\\Telomere_length.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

#** Facial aging -----------------------------------
path <- ("F:\\research\\aging_PE_GTN\\data\\exposure_gwas\\Facial_aging\\ukb-b-2148.vcf.gz")
FA_vcf <- read.vcfR(path, verbose = FALSE)

fix <- data.frame(FA_vcf@fix)
gt <- data.frame(FA_vcf@gt)

gt_sep <- separate(gt, "UKB.b.2148", into = c("ES","SE","LP","AF","ID"),
                   sep = ":", remove = TRUE, convert = TRUE)

FA_merge <- data.frame(fix, gt_sep)

FA_merge$exposure <- "Facial_aging"
FA_merge$p_val <- 10^-FA_merge$LP

FA_gwas <- select(FA_merge, exposure,ID,CHROM,POS,ALT,REF,AF,ES,SE,p_val)

colnames(FA_gwas)[c(2:10)] <- c("SNP","chr","position","ea","oa","eaf","b","se","p")

write.table(FA_gwas, "F:\\research\\aging_PE_GTN\\data\\exposure_gwas\\Facial_aging\\Facial_aging.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

#** Clonal hematopoiesis ---------------------------
path <- ("F:\\research\\aging_PE_GTN\\data\\exposure_gwas\\CHIP_decode\\ClonalHematopoiesis_All.txt")
CHIP_txt <- vroom(path, col_names = TRUE)

CHIP_txt$exposure <- "Clonal_hematopoiesis"
CHIP_txt$Chr <- as.numeric(gsub("chr", "", CHIP_txt$Chr)) 

## calculate the beta based on OR and convert the effect allele frequency
CHIP_txt$beta <- log(CHIP_txt$OR)
CHIP_txt$eaf <- CHIP_txt$EAFrq/100

## calculate the 'se' based on the 'beta' and 'P value'
CHIP_txt$se <- abs(CHIP_txt$beta/ qnorm(CHIP_txt$P/2))

CHIP_gwas <- select(CHIP_txt, exposure,rsID,Chr,PosB38,EA,OA,eaf,beta,se,P)

colnames(CHIP_gwas)[2:10] <- c("SNP","chr","position","ea","oa","eaf","b","se","p")

write.table(CHIP_gwas, "F:\\research\\aging_PE_GTN\\data\\exposure_gwas\\CHIP_decode\\Clonal_hematopoiesis.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

#** Phenotypic Age Acceleration ---------------------------
path <- ("F:\\research\\aging_PE_GTN\\data\\exposure_gwas\\PhenoAgeAccel\\phenoage.bolt.imputed.qc.txt")
PAA_txt <- vroom(path, col_names = TRUE)

PAA_txt$exposure <- "Pheno_Age_Accel"

PAA_gwas <- select(PAA_txt, exposure,SNP,CHR,BP,ALLELE1,ALLELE0,A1FREQ,BETA,SE,P_BOLT_LMM_INF)

colnames(PAA_gwas)[2:10] <- c("SNP","chr","position","ea","oa","eaf","b","se","p")

write.table(PAA_gwas, "F:\\research\\aging_PE_GTN\\data\\exposure_gwas\\PhenoAgeAccel\\PhenoAgeAccel.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

#** Biological Age Acceleration ---------------------------
path <- ("F:\\research\\aging_PE_GTN\\data\\exposure_gwas\\BioAgeAccel\\bioage.bolt.imputed.qc.txt")
BAA_txt <- vroom(path, col_names = TRUE)

BAA_txt$exposure <- "Bio_Age_Accel"

BAA_gwas <- select(BAA_txt, exposure,SNP,CHR,BP,ALLELE1,ALLELE0,A1FREQ,BETA,SE,P_BOLT_LMM_INF)

colnames(BAA_gwas)[2:10] <- c("SNP","chr","position","ea","oa","eaf","b","se","p")

write.table(BAA_gwas, "F:\\research\\aging_PE_GTN\\data\\exposure_gwas\\BioAgeAccel\\BioAgeAccel.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

#** Extrinsic Epigenetic age acceleration ------------------
path <- ("F:\\research\\aging_PE_GTN\\data\\exposure_gwas\\EpiAgeAccel\\EEAA_metaPlusGS.csv")
EEAA_csv <- vroom(path, col_names = TRUE)

EEAA_csv$exposure <- "Extr_Epi_Age_Accel"

EEAA_gwas <- select(EEAA_csv, exposure,SNP,Chr,bp,metaA1,metaA2,metaFreq,
                    meta_b,meta_se,metaPval)

colnames(EEAA_gwas)[2:10] <- c("SNP","chr","position","ea","oa","eaf","b","se","p")

write.table(EEAA_gwas, "F:\\research\\aging_PE_GTN\\data\\exposure_gwas\\EpiAgeAccel\\Extr_Epi_Age_Accel.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

#** Intrinsic Epigenetic age acceleration ------------------
path <- ("F:\\research\\aging_PE_GTN\\data\\exposure_gwas\\EpiAgeAccel\\IEAA_metaPlusGS.csv")
IEAA_csv <- vroom(path, col_names = TRUE)

IEAA_csv$exposure <- "Intr_Epi_Age_Accel"

IEAA_gwas <- select(IEAA_csv, exposure,SNP,Chr,bp,metaA1,metaA2,metaFreq,
                    meta_b,meta_se,metaPval)

colnames(IEAA_gwas)[2:10] <- c("SNP","chr","position","ea","oa","eaf","b","se","p")

write.table(IEAA_gwas, "F:\\research\\aging_PE_GTN\\data\\exposure_gwas\\EpiAgeAccel\\Ixtr_Epi_Age_Accel.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

#** Frailty index ------------------
path <- ("F:\\research\\aging_PE_GTN\\data\\exposure_gwas\\Frailty_index\\GCST90020053_buildGRCh37.tsv")
FI_tsv <- vroom(path, col_names = TRUE)

FI_tsv$exposure <- "Frailty_index"
colnames(FI_tsv)
FI_gwas <- select(FI_tsv, exposure,variant_id,chromosome,base_pair_location,
                  effect_allele,other_allele,effect_allele_frequency,beta,
                  standard_error,p_value)

colnames(FI_gwas)[2:10] <- c("SNP","chr","position","ea","oa","eaf","b","se","p")

write.table(FI_gwas, "F:\\research\\aging_PE_GTN\\data\\exposure_gwas\\Frailty_index\\Frailty_index.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

# MR main analysis



# Exposure data preparation for MR analysis ----------------
library(TwoSampleMR)
library(MRInstruments)

path_list <- list("F:\\research\\aging_PE_GTN\\data\\exposure_gwas\\BioAgeAccel\\BioAgeAccel.txt",
                  "F:\\research\\aging_PE_GTN\\data\\exposure_gwas\\CHIP_decode\\Clonal_hematopoiesis.txt",
                  "F:\\research\\aging_PE_GTN\\data\\exposure_gwas\\EpiAgeAccel\\Extr_Epi_Age_Accel.txt",
                  "F:\\research\\aging_PE_GTN\\data\\exposure_gwas\\EpiAgeAccel\\Ixtr_Epi_Age_Accel.txt",
                  "F:\\research\\aging_PE_GTN\\data\\exposure_gwas\\Facial_aging\\Facial_aging.txt",
                  "F:\\research\\aging_PE_GTN\\data\\exposure_gwas\\Frailty_index\\Frailty_index.txt",
                  "F:\\research\\aging_PE_GTN\\data\\exposure_gwas\\PhenoAgeAccel\\PhenoAgeAccel.txt",
                  "F:\\research\\aging_PE_GTN\\data\\exposure_gwas\\Telomere_length\\Telomere_length.txt")

name_list <- list("BioAgeAccel",
                  "Clonal_hematopoiesis",
                  "Extr_Epi_Age_Accel",
                  "Ixtr_Epi_Age_Accel",
                  "Facial_aging",
                  "Frailty_index",
                  "PhenoAgeAccel",
                  "Telomere_length")

for (i in 1:8) {gwas <- vroom(paste0(path_list[i]), col_names = TRUE)

## extract SNPs with p<5e-08.
gwas_p <- subset(gwas, p < 5e-08 & eaf > 0.01)

## format the data.
gwas_format <- format_data(gwas_p, type = "exposure",
                          snps = NULL,
                          header = TRUE,
                          phenotype_col = "exposure",
                          snp_col = "SNP",
                          beta_col = "b",
                          se_col = "se",
                          eaf_col = "eaf",
                          effect_allele_col = "ea",
                          other_allele_col = "oa",
                          pval_col = "p",
                          units_col = "units",
                          ncase_col = "ncase",
                          ncontrol_col = "ncontrol",
                          samplesize_col = "N",
                          gene_col = "gene",
                          id_col = "id",
                          min_pval = 0,
                          z_col = "z",
                          info_col = "info",
                          chr_col = "chr",
                          pos_col = "position",
                          log_pval = FALSE)

## LD clumping locally
library(ieugwasr)
library(genetics.binaRies)
##genetics.binaRies::get_plink_binary()
## you should download Plink from
## https://www.cog-genomics.org/plink/
## and replace the Plink in the folder
## "D:/software/R-4.3.1/library/genetics.binaRies/bin/plink.exe"
## by that you downloaded
gwas_ld <- dplyr::tibble(rsid=gwas_format$SNP, pval=gwas_format$pval.exposure, id=gwas_format$id.exposure)

ld_result <- ld_clump(dat =  gwas_ld,
                      clump_kb = 10000,
                      clump_r2 = 0.001,
                      clump_p = 1,
                      pop = "EUR",
                      plink_bin = genetics.binaRies::get_plink_binary(),
                      bfile = "F:/1000_gen_ref_panel/EUR")

assign(toString(name_list[i]), gwas_format[which(gwas_format$SNP %in% ld_result$rsid),])
}

Aging_exposure_data <- rbind(BioAgeAccel, Clonal_hematopoiesis,
                             Extr_Epi_Age_Accel, Ixtr_Epi_Age_Accel,
                             Facial_aging, Frailty_index, PhenoAgeAccel,
                             Telomere_length)

remove(list = setdiff(ls(),"Aging_exposure_data"))

## calculate the strength of instruments (F statistics)
# calculate R2
Aging_exposure_data$R_square <- (get_r_from_bsen(b = Aging_exposure_data$beta.exposure,
                                                 se = Aging_exposure_data$se.exposure,
                                                 n = Aging_exposure_data$Sample.Size))^2

# calculate F statistics
Aging_exposure_data$F.Statistics <- (Aging_exposure_data$R_square/(1-Aging_exposure_data$R_square))*
  ((Aging_exposure_data$Sample.Size-2))

Aging_expusure_data_F <- subset(Aging_exposure_data, F.Statistics > 10)

write.csv(Aging_expusure_data_F, file = "F:\\research\\aging_PE_GTN\\data\\data_for_ana\\Aging_exposure_data_Fstatistics.csv",
          quote = FALSE, row.names = FALSE)


# Mapping the SNPs of exposure to the outcome GWAS data --------------

Aging_exposure_data <- read.csv(file = "F:\\research\\aging_PE_GTN\\data\\data_for_ana\\Aging_exposure_data.csv",
                                header = TRUE, sep = ",")

## extract those SNPs included in exposure data on PE GWAS data.
PE_outcome_data <- read_outcome_data(snps = Aging_exposure_data$SNP,
                                  filename = "F:\\PE_and_geshtn_gwas\\PE_gwas.txt",
                                  sep = "\t",
                                  snp_col = "SNP",
                                  beta_col = "Effect",
                                  se_col = "StdErr",
                                  effect_allele_col = "Allele1",
                                  other_allele_col = "Allele2",
                                  eaf_col = "Freq1",
                                  pval_col = "P_value",
                                  units_col = "Units",
                                  gene_col = "Gene",
                                  samplesize_col = "n")

write.csv(PE_outcome_data, file = "F:\\research\\aging_PE_GTN\\data\\data_for_ana\\PE_outcome_data.csv",
          quote = FALSE, row.names = FALSE)

## extract those SNPs included in exposure data on GTN GWAS data.
GTN_outcome_data <- read_outcome_data(snps = Aging_exposure_data$SNP,
                                     filename = "F:\\PE_and_geshtn_gwas\\geshtn_gwas.txt",
                                     sep = "\t",
                                     snp_col = "SNP",
                                     beta_col = "Effect",
                                     se_col = "StdErr",
                                     effect_allele_col = "Allele1",
                                     other_allele_col = "Allele2",
                                     eaf_col = "Freq1",
                                     pval_col = "P_value",
                                     units_col = "Units",
                                     gene_col = "Gene",
                                     samplesize_col = "n")

write.csv(GTN_outcome_data, file = "F:\\research\\aging_PE_GTN\\data\\data_for_ana\\GTN_outcome_data.csv",
          quote = FALSE, row.names = FALSE)

## LD proxy
## in order to find proxy for the SNPs not mapping to outcome GWAS data
## those SNPs should be extracted from exposure data

# 2 SNPs of exposure data not mapping to PE GWAS data
PE_outcome_data <- read.csv(file = "F:\\research\\aging_PE_GTN\\data\\data_for_ana\\PE_outcome_data.csv",
                            header = TRUE, sep = ",")

GTN_outcome_data <- read.csv(file = "F:\\research\\aging_PE_GTN\\data\\data_for_ana\\GTN_outcome_data.csv",
                             header = TRUE, sep = ",")

SNPs_Diff_to_PE <- Aging_exposure_data[!(Aging_exposure_data$SNP %in% PE_outcome_data$SNP),]

write.csv(SNPs_Diff_to_PE, file = "F:\\research\\aging_PE_GTN\\data\\data_for_ana\\SNPs_of_exposure_not_mapping_to_PE_outcome.csv",
          quote = FALSE, row.names = FALSE)

# 13 SNPs of exposure data not mapping to GTN GWAS data
SNPs_Diff_to_GTN <- Aging_exposure_data[!(Aging_exposure_data$SNP %in% GTN_outcome_data$SNP),]

write.csv(SNPs_Diff_to_GTN, file = "F:\\research\\aging_PE_GTN\\data\\data_for_ana\\SNPs_of_exposure_not_mapping_to_GTN_outcome.csv",
          quote = FALSE, row.names = FALSE)

## We manually searched for proxies for those SNPs that were not mapped to outcome data
## through web-based applications https://ldlink.nci.nih.gov/.
## Finally there was no proxy retrieved for 3 SNPs in PE outcome data, and 4 proxies
## were retrieved for 11 SNPs in GTN outcome data.

## To conducting following MR analysis, we manually add proxies we retrieved into the
## GTN_outcome_data file (.csv) that outputted above steps.


# Main MR analysis for each exposure on PE and GTN ----------------
Aging_exposure_data <- read.csv(file = "F:\\research\\aging_PE_GTN\\data\\data_for_ana\\Aging_exposure_data.csv",
                                header = TRUE, sep = ",")

PE_outcome_data <- read.csv(file = "F:\\research\\aging_PE_GTN\\data\\data_for_ana\\PE_outcome_data.csv",
                            header = TRUE, sep = ",")

GTN_outcome_data <- read.csv(file = "F:\\research\\aging_PE_GTN\\data\\data_for_ana\\GTN_outcome_data.csv",
                             header = TRUE, sep = ",")

# to harmonize data (aging and PE)
Aging_PE_dat <- harmonise_data(exposure_dat = Aging_exposure_data, outcome_dat = PE_outcome_data)

# to harmonize data (aging and GTN)
Aging_GTN_dat <- harmonise_data(exposure_dat = Aging_exposure_data, outcome_dat = GTN_outcome_data)

## MR analysis of aging on PE
## method <- mr_method_list()
Aging_PE_result <- generate_odds_ratios(
      mr_res = mr(Aging_PE_dat,method_list = c("mr_egger_regression",
                                              "mr_weighted_median",
                                              "mr_ivw",
                                              "mr_ivw_fe",
                                              "mr_ivw_mre",
                                              "mr_weighted_mode",
                                              "mr_wald_ratio")))

write.csv(Aging_PE_result, file = "F:\\research\\aging_PE_GTN\\analysis_results\\table\\Aging_PE_result.csv",
          quote = FALSE, row.names = FALSE)

egger_intercept_PE <- mr_pleiotropy_test(Aging_PE_dat)

write.csv(egger_intercept_PE, file = "F:\\research\\aging_PE_GTN\\analysis_results\\table\\egger_intercept_PE.csv",
          quote = FALSE, row.names = FALSE)

## MR analysis of aging on GTN
Aging_GTN_result <- generate_odds_ratios(
  mr_res = mr(Aging_GTN_dat,method_list = c("mr_egger_regression",
                                           "mr_weighted_median",
                                           "mr_ivw",
                                           "mr_ivw_fe",
                                           "mr_ivw_mre",
                                           "mr_weighted_mode",
                                           "mr_wald_ratio")))

write.csv(Aging_GTN_result, file = "F:\\research\\aging_PE_GTN\\analysis_results\\table\\Aging_GTN_result.csv",
          quote = FALSE, row.names = FALSE)

egger_intercept_GTN <- mr_pleiotropy_test(Aging_GTN_dat)

write.csv(egger_intercept_GTN, file = "F:\\research\\aging_PE_GTN\\analysis_results\\table\\egger_intercept_GTN.csv",
          quote = FALSE, row.names = FALSE)



## MR-PRESSO analysis for each aging exposure on PE -----------

data_list <- split(Aging_PE_dat, Aging_PE_dat$exposure)
## Remove Extr_Epi_Age_Accel exposure, 
## because this exposure included only one SNP
data_list <- data_list[- 3] 

PE_presso_output <- list("presso_PE_Bio_Age_Accel",
                         "presso_PE_Clonal_hematopoiesis",
                         "presso_PE_Facial_aging",
                         "presso_PE_Frailty_index",
                         "presso_PE_Intr_Epi_Age_Accel",
                         "presso_PE_Pheno_Age_Accel",
                         "presso_PE_Telomere_length")

for (i in 1:7) {
  test <- subset(data.frame(data_list[[i]]), mr_keep == "TRUE")
  assign(paste0(PE_presso_output[[i]]), 
         run_mr_presso(test, NbDistribution = 3000, SignifThreshold = 0.05))
}

Main_results_presso_PE <- rbind(
    data.frame(presso_PE_Bio_Age_Accel[[1]][["Main MR results"]]),
    data.frame(presso_PE_Clonal_hematopoiesis[[1]][["Main MR results"]]),
    data.frame(presso_PE_Facial_aging[[1]][["Main MR results"]]),
    data.frame(presso_PE_Frailty_index[[1]][["Main MR results"]]),
    data.frame(presso_PE_Intr_Epi_Age_Accel[[1]][["Main MR results"]]),
    data.frame(presso_PE_Pheno_Age_Accel[[1]][["Main MR results"]]),
    data.frame(presso_PE_Telomere_length[[1]][["Main MR results"]]))

write.csv(Main_results_presso_PE, file = "F:\\research\\aging_PE_GTN\\analysis_results\\table\\Main_results_presso_PE.csv",
          quote = FALSE, row.names = FALSE)


Global_test_presso_PE <- rbind(
data.frame(presso_PE_Bio_Age_Accel[[1]][["MR-PRESSO results"]][["Global Test"]]),
data.frame(presso_PE_Clonal_hematopoiesis[[1]][["MR-PRESSO results"]][["Global Test"]]),
data.frame(presso_PE_Facial_aging[[1]][["MR-PRESSO results"]][["Global Test"]]),
data.frame(presso_PE_Frailty_index[[1]][["MR-PRESSO results"]][["Global Test"]]),
data.frame(presso_PE_Intr_Epi_Age_Accel[[1]][["MR-PRESSO results"]][["Global Test"]]),
data.frame(presso_PE_Pheno_Age_Accel[[1]][["MR-PRESSO results"]][["Global Test"]]),
data.frame(presso_PE_Telomere_length[[1]][["MR-PRESSO results"]][["Global Test"]]))

write.csv(Global_test_presso_PE, file = "F:\\research\\aging_PE_GTN\\analysis_results\\table\\Global_test_presso_PE.csv",
          quote = FALSE, row.names = FALSE)


test1 <- data.frame(presso_PE_Bio_Age_Accel[[1]][["MR-PRESSO results"]][["Outlier Test"]])
  test1$exposure <- "Bio_Age_Accel"
test2 <- data.frame(presso_PE_Clonal_hematopoiesis[[1]][["MR-PRESSO results"]][["Outlier Test"]])
  ## No outlier identified for CHIP exposure
test3 <- data.frame(presso_PE_Facial_aging[[1]][["MR-PRESSO results"]][["Outlier Test"]])
  test3$exposure <- "Facial_aging"
test4 <- data.frame(presso_PE_Frailty_index[[1]][["MR-PRESSO results"]][["Outlier Test"]])
  ## No outlier identified for Frailty_index
test5 <- data.frame(presso_PE_Intr_Epi_Age_Accel[[1]][["MR-PRESSO results"]][["Outlier Test"]])
  ## No outlier identified for Intr_Epi_Age_Accel
test6 <- data.frame(presso_PE_Pheno_Age_Accel[[1]][["MR-PRESSO results"]][["Outlier Test"]])
  test6$exposure <- "Pheno_Age_Accel"
test7 <- data.frame(presso_PE_Telomere_length[[1]][["MR-PRESSO results"]][["Outlier Test"]])
  test7$exposure <- "Telomere_length"

presso_PE_outlier <- rbind(test1, test3, test6, test7)

write.csv(presso_PE_outlier, file = "F:\\research\\aging_PE_GTN\\analysis_results\\table\\presso_PE_outlier.csv",
          quote = FALSE, row.names = FALSE)

## MR-PRESSO analysis for each aging exposure on GTN -----------

data_list <- split(Aging_GTN_dat, Aging_GTN_dat$exposure)
## Remove Extr_Epi_Age_Accel exposure, 
## because this exposure included only one SNP
data_list <- data_list[- 3] 

GTN_presso_output <- list("presso_GTN_Bio_Age_Accel",
                         "presso_GTN_Clonal_hematopoiesis",
                         "presso_GTN_Facial_aging",
                         "presso_GTN_Frailty_index",
                         "presso_GTN_Intr_Epi_Age_Accel",
                         "presso_GTN_Pheno_Age_Accel",
                         "presso_GTN_Telomere_length")

for (i in 1:7) {
  test <- subset(data.frame(data_list[[i]]), mr_keep == "TRUE")
  assign(paste0(GTN_presso_output[[i]]), 
         run_mr_presso(test, NbDistribution = 3000, SignifThreshold = 0.05))
}

Main_results_presso_GTN <- rbind(
  data.frame(presso_GTN_Bio_Age_Accel[[1]][["Main MR results"]]),
  data.frame(presso_GTN_Clonal_hematopoiesis[[1]][["Main MR results"]]),
  data.frame(presso_GTN_Facial_aging[[1]][["Main MR results"]]),
  data.frame(presso_GTN_Frailty_index[[1]][["Main MR results"]]),
  data.frame(presso_GTN_Intr_Epi_Age_Accel[[1]][["Main MR results"]]),
  data.frame(presso_GTN_Pheno_Age_Accel[[1]][["Main MR results"]]),
  data.frame(presso_GTN_Telomere_length[[1]][["Main MR results"]]))

write.csv(Main_results_presso_GTN, file = "F:\\research\\aging_PE_GTN\\analysis_results\\table\\Main_results_presso_GTN.csv",
          quote = FALSE, row.names = FALSE)

Global_test_presso_GTN <- rbind(
  data.frame(presso_GTN_Bio_Age_Accel[[1]][["MR-PRESSO results"]][["Global Test"]]),
  data.frame(presso_GTN_Clonal_hematopoiesis[[1]][["MR-PRESSO results"]][["Global Test"]]),
  data.frame(presso_GTN_Facial_aging[[1]][["MR-PRESSO results"]][["Global Test"]]),
  data.frame(presso_GTN_Frailty_index[[1]][["MR-PRESSO results"]][["Global Test"]]),
  data.frame(presso_GTN_Intr_Epi_Age_Accel[[1]][["MR-PRESSO results"]][["Global Test"]]),
  data.frame(presso_GTN_Pheno_Age_Accel[[1]][["MR-PRESSO results"]][["Global Test"]]),
  data.frame(presso_GTN_Telomere_length[[1]][["MR-PRESSO results"]][["Global Test"]]))

write.csv(Global_test_presso_GTN, file = "F:\\research\\aging_PE_GTN\\analysis_results\\table\\Global_test_presso_GTN.csv",
          quote = FALSE, row.names = FALSE)

test1 <- data.frame(presso_GTN_Bio_Age_Accel[[1]][["MR-PRESSO results"]][["Outlier Test"]])
test1$exposure <- "Bio_Age_Accel"
test2 <- data.frame(presso_GTN_Clonal_hematopoiesis[[1]][["MR-PRESSO results"]][["Outlier Test"]])
## No outlier identified for CHIP exposure
test3 <- data.frame(presso_GTN_Facial_aging[[1]][["MR-PRESSO results"]][["Outlier Test"]])
## No outlier identified for Facial_aging
test4 <- data.frame(presso_GTN_Frailty_index[[1]][["MR-PRESSO results"]][["Outlier Test"]])
## No outlier identified for Frailty_index
test5 <- data.frame(presso_GTN_Intr_Epi_Age_Accel[[1]][["MR-PRESSO results"]][["Outlier Test"]])
## No outlier identified for Intr_Epi_Age_Accel
test6 <- data.frame(presso_GTN_Pheno_Age_Accel[[1]][["MR-PRESSO results"]][["Outlier Test"]])
test6$exposure <- "Pheno_Age_Accel"
test7 <- data.frame(presso_GTN_Telomere_length[[1]][["MR-PRESSO results"]][["Outlier Test"]])
## No outlier identified for Telomere_length

presso_GTN_outlier <- rbind(test1, test6)

write.csv(presso_GTN_outlier, file = "F:\\research\\aging_PE_GTN\\analysis_results\\table\\presso_GTN_outlier.csv",
          quote = FALSE, row.names = FALSE)

remove(list = setdiff(ls(),c("Aging_PE_dat","Aging_GTN_dat")))


test1 <- subset(Aging_PE_dat, mr_keep == "TRUE")
test2 <- subset(Aging_GTN_dat, mr_keep == "TRUE")
write.csv()


