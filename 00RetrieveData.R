

## 暴露数据
# deCODE
library(TwoSampleMR)
library(tidyverse)
library(future)
library(future.apply)

plan("multicore", workers = 4) ###compute cores
options(future.globals.maxSize = 30000 * 1024^2)

#ao <- available_outcomes()

exposure <- data.table::fread("data/3168_8_ADAMTS5_ADAMTS_5.txt.gz")%>% 
  subset(Pval < 5e-08) %>% 
  dplyr::filter(!is.na(rsids)) %>% 
  dplyr::mutate(phenotype = "ADAMST5",
                id = "deCODE") %>%
  TwoSampleMR::format_data(.,type = "exposure",phenotype_col = "phenotype",
                           beta_col = "Beta",se_col = "SE",id_col = "id",snp_col = "rsids",
                           pval_col = "Pval",effect_allele_col = "effectAllele",
                           other_allele_col = "otherAllele", chr_col = "Chrom",
                           eaf_col = "ImpMAF",pos_col = "Pos",samplesize_col = "N") %>%
  clump_data(.)

saveRDS(exposure,file = "data/ADAMTS5_deCODE.rds")

## IEU
# 最新数据 
# IEU数据
repeat{
  try({
    exposure_ieu <- extract_instruments(outcomes = "prot-c-3168_8_2")
  })
  if(exists("exposure_ieu")){break}
  Sys.sleep(2)
}

##saveRDS(exposure_ieu,file = "data/ADAMTS5_ieu.rds")

## 结局变量
##############不需要运行
# 在获取结局id
outcomeid <- c("ieu-a-1120","ieu-b-4960", "ieu-b-96","ieu-a-1129","ieu-a-966","ieu-a-822","ieu-b-85",
               "finn-b-C3_COLON_EXALLC","finn-b-C3_RECTUM_EXALLC","finn-b-C3_BLADDER_EXALLC",
               'finn-b-C3_MELANOMA_SKIN',"finn-b-CD2_HODGKIN_LYMPHOMA_EXALLC","finn-b-CD2_NONHODGKIN_NAS_EXALLC",
               "finn-b-C3_STOMACH_EXALLC","finn-b-C3_OESOPHAGUS_EXALLC","finn-b-C3_KIDNEY_NOTRENALPELVIS_EXALLC")

#################################
## 结局
# 在线加本地 
# IEU数据
repeat{
  try({
    outcome <- extract_outcome_data(snps = c(exposure$SNP,exposure_ieu$SNP),proxies = T,rsq = 0.8,
                                    outcomes = c("ieu-a-1120", "ieu-b-96",
                                                 "ieu-a-1129","ieu-a-966",
                                                 "ieu-a-822","ieu-b-85"))
  })
  if(exists("outcome")){break}
  Sys.sleep(2)
}

#saveRDS(outcome,file = "data/outcome_ieu.rds")
## 下载

# colon
fun1 <- paste0("wget -c -nv -P data ", "https://storage.googleapis.com/finngen-public-data-r9/summary_stats/finngen_R9_C3_COLON_EXALLC.gz")
system(fun1)

# rectum
fun2 <- paste0("wget -c -nv -P data ", "https://storage.googleapis.com/finngen-public-data-r9/summary_stats/finngen_R9_C3_RECTUM_EXALLC.gz")
system(fun2)

# 膀胱
fun3 <- paste0("wget -c -nv -P data ", "https://storage.googleapis.com/finngen-public-data-r9/summary_stats/finngen_R9_C3_BLADDER_EXALLC.gz")
system(fun3)

# 肾

fun4 <- paste0("wget -c -nv -P data ", "https://storage.googleapis.com/finngen-public-data-r9/summary_stats/finngen_R9_C3_KIDNEY_NOTRENALPELVIS_EXALLC.gz")
system(fun4)

# 食管
fun5 <- paste0("wget -c -nv -P data ", "https://storage.googleapis.com/finngen-public-data-r9/summary_stats/finngen_R9_C3_OESOPHAGUS_EXALLC.gz")
system(fun5)

# 胃
fun6 <- paste0("wget -c -nv -P data ", "https://storage.googleapis.com/finngen-public-data-r9/summary_stats/finngen_R9_C3_STOMACH_EXALLC.gz")
system(fun6)

# 非霍奇金
fun7 <- paste0("wget -c -nv -P data ", "https://storage.googleapis.com/finngen-public-data-r9/summary_stats/finngen_R9_C3_NONHODGKIN_EXALLC.gz")
system(fun7)

# 霍奇金
fun8 <- paste0("wget -c -nv -P data ", "https://storage.googleapis.com/finngen-public-data-r9/summary_stats/finngen_R9_CD2_HODGKIN_LYMPHOMA_EXALLC.gz")
system(fun8)

# 黑色素
fun9 <- paste0("wget -c -nv -P data ", "https://storage.googleapis.com/finngen-public-data-r9/summary_stats/finngen_R9_C3_MELANOMA_EXALLC.gz")
system(fun9)

# 口腔
fun10 <- paste0("wget -c -nv -P data ", "https://storage.googleapis.com/finngen-public-data-r9/summary_stats/finngen_R9_C3_ORALCAVITY_EXALLC.gz")
system(fun10)

###甲状腺
fun10 <- paste0("wget -c -nv -P data ", "https://storage.googleapis.com/finngen-public-data-r9/summary_stats/finngen_R9_C3_THYROID_GLAND_EXALLC.gz")
system(fun10)

## 读入芬兰数据形成结局文件
files <- list.files("data",pattern = "finngen_R9")

readfinne <- function(file,SNP){
  
  filename <- file.path("data",file)
  manifest <- data.table::fread("data/summary_stats_R9_manifest.tsv",sep = "\t",quote = F)
  
  outcome_finn <- data.table::fread(file.path("data",file)) %>%  
    dplyr::mutate(.,phenotype = strsplit(file,split = "_")[[1]][4],
                  id = sub(".gz","",file),
                  ncase = manifest$num_cases[manifest$phenocode == paste(strsplit(sub(".gz","",file),split = "_")[[1]][-c(1:2)],collapse = "_")],
                  ncontrol = manifest$num_controls[manifest$phenocode == paste(strsplit(sub(".gz","",file),split = "_")[[1]][-c(1:2)],collapse = "_")]) %>%
    dplyr::mutate(N = ncase + ncontrol) %>%
    TwoSampleMR::format_data(.,type = "outcome",snps = SNP,phenotype_col = "phenotype",
                             beta_col = "beta",se_col = "sebeta",id_col = "id",snp_col = "rsids",
                             pval_col = "pval",effect_allele_col = "alt",other_allele_col = "ref",
                             chr_col = "#chrom",eaf_col = "af_alt",ncase_col = "ncase",pos_col = "pos",
                             ncontrol_col = "ncontrol",samplesize_col = "N",gene_col = "nearest_genes")
}


outcome_finn <- lapply(files, function(i){readfinne(i,SNP = c(exposure$SNP,exposure_ieu$SNP))})
outcome_finn <- do.call(rbind,outcome_finn)
###saveRDS(outcome_finn,file = "data/outcome_finn.rds")

outcome <- outcome_finn %>%
  dplyr::mutate(originalname.outcome = outcome,
                outcome.deprecated = paste0(outcome, " || "),
                data_source.outcome = "FinnGen",
                proxy.outcome = NA,
                target_snp.outcome = NA,
                proxy_snp.outcome = NA,
                target_a1.outcome = NA,
                target_a2.outcome = NA,
                proxy_a1.outcome = NA,
                proxy_a2.outcome = NA) %>%
  dplyr::select(SNP,chr.outcome,pos.outcome,beta.outcome,se.outcome,
         samplesize.outcome,pval.outcome,eaf.outcome,
         effect_allele.outcome,other_allele.outcome,outcome,
         id.outcome,originalname.outcome,outcome.deprecated,
         mr_keep.outcome,data_source.outcome,proxy.outcome,
         target_snp.outcome,proxy_snp.outcome,target_a1.outcome,
         target_a2.outcome,proxy_a1.outcome,proxy_a2.outcome) %>%
  dplyr::rename(chr = "chr.outcome",pos = "pos.outcome")  %>%
  rbind(outcome)

saveRDS(outcome,file = "data/canceroutcome.rds")

######
fun12 <- paste0("wget -c -nv -P data ","https://gwas.mrcieu.ac.uk/files/ieu-a-1129/ieu-a-1129.vcf.gz") 
system(fun12)

fun13 <- paste0("wget -c -nv -P data ","https://gwas.mrcieu.ac.uk/files/ieu-a-1120/ieu-a-1120.vcf.gz") 
system(fun13)

fun14 <- paste0("wget -c -nv -P data ","https://gwas.mrcieu.ac.uk/files/ieu-b-85/ieu-b-85.vcf.gz") 
system(fun14)

fun15 <- paste0("wget -c -nv -P data ","https://gwas.mrcieu.ac.uk/files/ieu-b-96/ieu-b-96.vcf.gz") 
system(fun15)

fun16 <- paste0("wget -c -nv -P data ","https://gwas.mrcieu.ac.uk/files/prot-c-3168_8_2/prot-c-3168_8_2.vcf.gz") 
system(fun16)


##ieu暴露

repeat{
  try({
    exposure_ieu <- extract_instruments(outcomes = c("ieu-a-1120", "ieu-b-96","ieu-a-1129",
                                                     "ieu-a-966","ieu-a-822","ieu-b-85")) 
  })
  if(exists("exposure_ieu")){break}
  Sys.sleep(2)
}
exposure_ieu <- exposure_ieu %>%
  dplyr::mutate(.,exposure = unlist(lapply(exposure,function(i){strsplit(i,split = " cancer")[[1]][1]})))

repeat{
  try({
    exposure_ieu_tmp <- extract_instruments(outcomes = "ieu-a-822",p1 = 5e-06)
  })
  if(exists("exposure_ieu_tmp")){break}
  Sys.sleep(2)
}

exposure_ieu <- exposure_ieu %>%
  rbind(.,dplyr::mutate(exposure_ieu_tmp,exposure = "Pancreatic"))


files <- list.files("data",pattern = ".gz")[c(2:13)]

readfinnetoexposure <- function(file){
  
  filename <- file.path("data",file)
  manifest <- data.table::fread("data/summary_stats_R9_manifest.tsv",sep = "\t",quote = F)
  
  exposure_finn <- data.table::fread(file.path("data",file)) %>%  
    dplyr::mutate(.,phenotype = strsplit(file,split = "_")[[1]][4],
                  id = sub(".gz","",file),
                  ncase = manifest$num_cases[manifest$phenocode == paste(strsplit(sub(".gz","",file),split = "_")[[1]][-c(1:2)],collapse = "_")],
                  ncontrol = manifest$num_controls[manifest$phenocode == paste(strsplit(sub(".gz","",file),split = "_")[[1]][-c(1:2)],collapse = "_")]) %>%
    dplyr::mutate(N = ncase + ncontrol) 
  
  if(sum(exposure_finn$pval < 5e-08) == 0){
    message(paste0(file," # no snp < 5e-08"))
    exposure_finn <- subset(exposure_finn,pval < 5e-06) %>%
      TwoSampleMR::format_data(.,type = "exposure",phenotype_col = "phenotype",
                               beta_col = "beta",se_col = "sebeta",id_col = "id",snp_col = "rsids",
                               pval_col = "pval",effect_allele_col = "alt",other_allele_col = "ref",
                               chr_col = "#chrom",eaf_col = "af_alt",ncase_col = "ncase",pos_col = "pos",
                               ncontrol_col = "ncontrol",samplesize_col = "N",gene_col = "nearest_genes") 
  }else{
    exposure_finn <-  subset(exposure_finn,pval < 5e-08) %>%
      TwoSampleMR::format_data(.,type = "exposure",phenotype_col = "phenotype",
                               beta_col = "beta",se_col = "sebeta",id_col = "id",snp_col = "rsids",
                               pval_col = "pval",effect_allele_col = "alt",other_allele_col = "ref",
                               chr_col = "#chrom",eaf_col = "af_alt",ncase_col = "ncase",pos_col = "pos",
                               ncontrol_col = "ncontrol",samplesize_col = "N",gene_col = "nearest_genes") 
  }
  
  repeat{
    try({
      exposure <- TwoSampleMR::clump_data(exposure_finn)
    })
    if(exists("exposure")){break}
    Sys.sleep(2)
  }
  
  return(exposure)
}

exposure_finn <- lapply(files, function(i){readfinnetoexposure(i)})
exposure_finn <- do.call(rbind,exposure_finn)  %>%
  dplyr::select(colnames(exposure_ieu)[-15]) %>%
  rbind(exposure_ieu[,-15]) 

#saveRDS(exposure_finn,file = "data/cancer_exposure.rds")

manifest <- data.table::fread("data/summary_stats_R9_manifest.tsv",sep = "\t",quote = F)
file <- "finngen_R9_C3_KIDNEY_NOTRENALPELVIS_EXALLC.gz"
exposure_finn <- data.table::fread(file.path("data",file)) %>%  
  dplyr::mutate(.,phenotype = strsplit(file,split = "_")[[1]][4],
                id = sub(".gz","",file),
                ncase = manifest$num_cases[manifest$phenocode == paste(strsplit(sub(".gz","",file),split = "_")[[1]][-c(1:2)],collapse = "_")],
                ncontrol = manifest$num_controls[manifest$phenocode == paste(strsplit(sub(".gz","",file),split = "_")[[1]][-c(1:2)],collapse = "_")]) %>%
  dplyr::mutate(N = ncase + ncontrol) %>%
  subset(pval < 5e-06) %>%
  TwoSampleMR::format_data(.,type = "exposure",phenotype_col = "phenotype",
                           beta_col = "beta",se_col = "sebeta",id_col = "id",snp_col = "rsids",
                           pval_col = "pval",effect_allele_col = "alt",other_allele_col = "ref",
                           chr_col = "#chrom",eaf_col = "af_alt",ncase_col = "ncase",pos_col = "pos",
                           ncontrol_col = "ncontrol",samplesize_col = "N",gene_col = "nearest_genes") %>%
  TwoSampleMR::clump_data()


tmp <- readRDS("data/cancer_exposure.rds") %>%
  dplyr::filter(id.exposure != "finngen_R9_C3_KIDNEY_NOTRENALPELVIS_EXALL") %>%
  dplyr::filter(grepl("finng",id.exposure)) %>%
  rbind(.,dplyr::select(exposure_finn,colnames(.))) %>%
  rbind(.,dplyr::select(exposure_ieu,colnames(.)))


saveRDS(tmp,file = "data/cancer_exposurenew.rds")

###结局
ADAMTSoutcome <- data.table::fread("data/3168_8_ADAMTS5_ADAMTS_5.txt.gz") %>%
  dplyr::filter(Pval >= 5e-08) %>%
  dplyr::filter(!is.na(rsids)) %>% 
  dplyr::mutate(phenotype = "ADAMST5",
                id = "deCODE") %>%
  TwoSampleMR::format_data(.,type = "outcome",snp = tmp$SNP,phenotype_col = "phenotype",
                           beta_col = "Beta",se_col = "SE",id_col = "id",snp_col = "rsids",
                           pval_col = "Pval",effect_allele_col = "effectAllele",
                           other_allele_col = "otherAllele", chr_col = "Chrom",
                           eaf_col = "ImpMAF",pos_col = "Pos",samplesize_col = "N")

saveRDS(ADAMTSoutcome,file = "data/ADAMTS5_outcome.rds")




