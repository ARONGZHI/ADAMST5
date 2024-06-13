rm(list = ls())
## coloc 数据
library(vroom)
library(MungeSumstats)
library(tidyverse)
library(locuscomparer)

################数据预处理
leadSNP <- readRDS("data/ADAMTS5_ieu.rds") 
proADAMTS <- VariantAnnotation::readVcf("data/prot-c-3168_8_2.vcf.gz") %>% gwasglue::gwasvcf_to_TwoSampleMR()

data1 <- proADAMTS %>% dplyr::filter(chr.exposure == leadSNP$chr.exposure,
                                     pos.exposure >= leadSNP$pos.exposure - 1000*1000,
                                     pos.exposure <= leadSNP$pos.exposure + 1000*1000)

oropharyngeal <- VariantAnnotation::readVcf("data/ieu-b-96.vcf.gz") %>% 
  gwasglue::gwasvcf_to_TwoSampleMR(type = "outcome")

data2 <- oropharyngeal %>% dplyr::filter(chr.outcome == leadSNP$chr.exposure,
                                         pos.outcome >= leadSNP$pos.exposure - 1000*1000,
                                         pos.outcome <= leadSNP$pos.exposure + 1000*1000)

data <- merge(data1,data2,by = "SNP") %>%
  dplyr::filter(!duplicated(SNP))

#saveRDS(data,file = "data/colocieu.rds")
##########################

### coloc分析
colocieu <- readRDS("data/colocieu.rds")
colocieu <-  colocieu %>% filter((effect_allele.exposure==effect_allele.outcome & other_allele.exposure==other_allele.outcome)|
                                   (effect_allele.exposure==other_allele.outcome & other_allele.exposure==effect_allele.outcome)) %>%
  mutate(.,beta.outcome = ifelse(effect_allele.exposure==effect_allele.outcome,beta.outcome,
                                 -beta.outcome)) %>%
  mutate(.,VAR.x = se.exposure^2,
         VAR.y = se.outcome^2) %>%
  dplyr::filter(VAR.x!=0 & VAR.y!=0 )



exposuretabieu <- colocieu[,c("beta.exposure","VAR.x","SNP","pval.exposure")] %>%
  mutate(.,sdY = 1) %>%
  rename("beta" = beta.exposure,"varbeta" = VAR.x,"snp" =SNP) %>%
  as.list() 

outcometab <-  colocieu[,c("beta.outcome","VAR.y","SNP","pval.outcome")] %>%
  rename("beta" = beta.outcome,"varbeta" = VAR.y,"snp" =SNP) %>%
  as.list()

exposuretabieu$type <- "quant"
outcometab$type = "cc"

resieu <-  coloc::coloc.abf(exposuretabieu,outcometab,p1=1e-4,p2=1e-4,p12=1e-5)
expousreieu_dat <- data.frame(rsid = exposuretabieu$snp,
                              pval = exposuretabieu$pval.exposure) %>%
  .[order(.$rsid),]

outcomeieu_dat <- data.frame(rsid = outcometab$snp,
                             pval = outcometab$pval.outcome) %>%
  .[order(.$rsid),]

pieu <- locuscompare(in_fn1 = outcomeieu_dat, in_fn2 = expousreieu_dat, 
                     title1 = "Oropharyngeal cancer GWAS",title2 = "ADAMTS-5 pQTL",combine = F,legend = F)

##################deCODE
leadSNPdeCODE <- readRDS("data/ADAMTS5_deCODE.rds") %>%
  arrange(pval.exposure)

proADAMTSdeCODE <- data.table::fread("data/3168_8_ADAMTS5_ADAMTS_5.txt.gz")  

data1 <- subset(proADAMTSdeCODE, Chrom == leadSNPdeCODE$chr.exposure[1] &
                  Pos > (leadSNPdeCODE$pos[1] - 1000*1000) & 
                  Pos < (leadSNPdeCODE$pos[1] + 1000*1000)) %>%
  dplyr::filter(!is.na(rsids)) %>% 
  dplyr::mutate(phenotype = "ADAMST5",
                id = "deCODE") %>%
  TwoSampleMR::format_data(.,type = "exposure",phenotype_col = "phenotype",
                           beta_col = "Beta",se_col = "SE",id_col = "id",snp_col = "rsids",
                           pval_col = "Pval",effect_allele_col = "effectAllele",
                           other_allele_col = "otherAllele", chr_col = "Chrom",
                           eaf_col = "ImpMAF",pos_col = "Pos") %>%
  mutate(samplesize.exposure = 35559) %>%
  subset(., !duplicated(SNP)) 

oropharyngeal <- VariantAnnotation::readVcf("data/ieu-b-96.vcf.gz") %>% 
  gwasglue::gwasvcf_to_TwoSampleMR(type = "outcome")

sumstats_dt  <- oropharyngeal %>%
  select("chr.outcome","pos.outcome","other_allele.outcome","effect_allele.outcome",
         "SNP","pval.outcome","beta.outcome","se.outcome","eaf.outcome") %>%
  rename("CHR"= chr.outcome,"BP" = pos.outcome,"A1" = other_allele.outcome,
         "A2" =effect_allele.outcome,"P" = pval.outcome,"BETA" = beta.outcome ,
         "SE" = se.outcome,'FRQ' = eaf.outcome) %>%
  select("SNP","CHR","BP","A1","A2","FRQ",'BETA',"SE",'P') %>%
  dplyr::filter(complete.cases(BP))

# 参考基因组转换，37to38====
data2 <- liftover(sumstats_dt=sumstats_dt, 
                  ref_genome = "hg19",
                  convert_ref_genome="hg38") %>%
  subset(CHR == 21 & BP > (leadSNPdeCODE$pos[1] - 1000*1000) & 
           BP < (leadSNPdeCODE$pos[1] + 1000*1000)) %>% 
  subset(., !duplicated(SNP)) %>%
  dplyr::mutate(phenotype = "oropharyngeal",samplesize = 997) %>%
  TwoSampleMR::format_data(.,type = "outcome",phenotype_col = "phenotype",
                           beta_col = "BETA",se_col = "SE",snp_col = "SNP",
                           pval_col = "P",effect_allele_col = "A2",
                           other_allele_col = "A1", chr_col = "CHR",
                           eaf_col = "FRQ",pos_col = "BP",samplesize_col = "samplesize") 

data <- merge(data1,data2,by = "SNP")  %>%
  dplyr::filter(!duplicated(SNP))

#saveRDS(data,file = "data/colocdeCODE.rds")
###########################
colocdecode <- readRDS("data/colocdeCODE.rds") %>%
  filter((effect_allele.exposure==effect_allele.outcome & other_allele.exposure==other_allele.outcome)|
           (effect_allele.exposure==other_allele.outcome & other_allele.exposure==effect_allele.outcome)) %>%
  mutate(.,beta.outcome = ifelse(effect_allele.exposure==effect_allele.outcome,beta.outcome,
                                 -beta.outcome)) %>%
  mutate(.,VAR.x = se.exposure^2,
         VAR.y = se.outcome^2) %>%
  dplyr::filter(VAR.x!=0 & VAR.y!=0 )

exposuretabdecode <- colocdecode[,c("beta.exposure","VAR.x","SNP","eaf.exposure","samplesize.exposure","pval.exposure")] %>%
  rename("beta" = beta.exposure,"varbeta" = VAR.x,"snp" =SNP,"MAF" = eaf.exposure,N = samplesize.exposure) %>%
  as.list() 

outcometabdecode <-  colocieu[,c("beta.outcome","VAR.y","SNP","pval.outcome")] %>%
  rename("beta" = beta.outcome,"varbeta" = VAR.y,"snp" =SNP) %>%
  as.list()

exposuretabdecode$type <- "quant"
outcometabdecode$type = "cc"

resdecode <-  coloc::coloc.abf(exposuretabdecode,outcometabdecode,p1=1e-4,p2=1e-4,p12=1e-5)


exposuretabdecode_dat  <- data.frame(rsid = exposuretabdecode$snp,
                                     pval = exposuretabdecode$pval.exposure) %>%
  .[order(.$rsid),]

outcometabdecode_dat <- data.frame(rsid = outcometabdecode$snp,
                                   pval = outcometabdecode$pval.outcome) %>%
  .[order(.$rsid),]

pdecode <- locuscompare(in_fn1 = outcometabdecode_dat, in_fn2 = exposuretabdecode_dat, 
                        title1 = "Oropharyngeal cancer GWAS",title2 = "ADAMTS-5 pQTL",combine = F,legend = F)

####Figure4
Figure4 <- gridExtra::grid.arrange(pdecode$locuscompare + ggtitle("A",subtitle = "deCODE database"),
                                   pieu$locuscompare + ggtitle("B",subtitle = 'IEU Open GWAS project'),nrow = 1, ncol = 2)

library(cowplot)
save_plot(Figure4,filename = "image/Figure4.pdf",base_width = 10,base_height = 5)

####TableS5

TableS5 <- as.data.frame(round(resdecode$summary,5)) %>%
  t() %>%
  rbind(.,t(as.data.frame(round(resieu$summary,5)))) 

rownames(TableS5) <- c("deCODE","IEU Open GWAS project")

write.csv(TableS5,file = "result/tableS5.csv",quote =  F,row.names = T)