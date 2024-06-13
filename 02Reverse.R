rm(list = ls())

library(tidyverse)
library(forestplot)
## 暴露
exposure <- readRDS("data/cancer_exposurenew.rds")  %>%
  dplyr::filter(id.exposure != "finngen_R9_C3_COLORECTAL_EXALLC") %>%
  dplyr::mutate(.,#R2 = 2*eaf.exposure*(1 - eaf.exposure) *beta.exposure^2,
                R2 = (beta.exposure^2)/(beta.exposure^2 + samplesize.exposure*se.exposure^2),
                Fstatic =  ((samplesize.exposure - 1 - 1)/1)*(R2/(1 - R2))) %>%
  dplyr::filter(Fstatic > 01) 

##结局
ADAMTSoutcome <- readRDS("data/ADAMTS5_outcome.rds") %>%
  dplyr::filter(pval.outcome > 5e-08)

## Harmonize
dat_ADAMTS <- harmonise_data(exposure_dat = exposure,outcome_dat = ADAMTSoutcome,action = 2)

## sterger 过滤
S2 <- steiger_filtering(dat_ADAMTS) %>%
  dplyr::select(SNP,effect_allele.exposure,other_allele.exposure,effect_allele.outcome,
                other_allele.outcome,exposure,outcome,eaf.exposure,beta.exposure,se.exposure,pval.exposure,
                eaf.outcome,beta.outcome,se.outcome,pval.outcome,R2,Fstatic,steiger_dir)

write.csv(S2,file = "result/S2.csv",quote = F)


dat_ADAMTS  <- steiger_filtering(dat_ADAMTS) %>%
  subset(.,steiger_dir == TRUE)

######MRPRESSO 过滤异常值
library(MRPRESSO)

presso_modified <- function(id,dat){
  tmp <- dplyr::filter(dat,exposure == id )
  presso <-mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", 
                     OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data= tmp, NbDistribution = 1000,  SignifThreshold = 0.05,seed=11)
}


presso_ADAMTS <- lapply(unique(dat_ADAMTS$exposure)[-c(1,2,9,11,14,17)], function(i) {
  presso_modified(i,dat_ADAMTS)})

#saveRDS(presso_ADAMTS,file = "result/ADAMTSpresso.rds")
#presso_ADAMTS <- readRDS("result/ADAMTSpresso.rds")

presso_ADAMTSres <- data.frame(exposure = unique(dat_ADAMTS$exposure)[-c(1,2,9,11,14,17)],
                               pvalue = unlist(lapply(1:length(presso_ADAMTS),function(i){
                                 presso_ADAMTS[[i]]$`MR-PRESSO results`$`Global Test`$Pvalue}))) %>%
  dplyr::mutate(.,`Number of Outlier SNPs` = NA,`Outlier SNPs` = NA,`Outlier-corrected OR (95%CI)`=NA, 
                `Outlier-corrected P-value` = NA, `P-value of distortion test` = NA)
presso_ADAMTSres[9,"Number of Outlier SNPs"] <- 0

## MR分析mr(Mydata, method_list=c("mr_ivw", "mr_ivw_fe", "mr_two_sample_ml","mr_egger_regression", "mr_weighted_median", "mr_penalised_weighted_median","mr_simple_mode", "mr_weighted_mode"))
mr_ADAMTS <- mr(dat_ADAMTS) %>%
  generate_odds_ratios()

breastmr <-  mr(dplyr::filter(dat_ADAMTS,exposure == "Breast"),method_list = c("mr_ivw_mre")) %>%
  generate_odds_ratios()

# figure2
mr_tabledecode <- mr_ADAMTS %>%
  dplyr::filter(method %in% c("Wald ratio","Inverse variance weighted")) %>%
  dplyr::mutate(.,Exposure = Hmisc::capitalize(tolower(unlist(lapply(exposure, 
                                                                     function(i){strsplit(i,split = " ")[[1]][1]})))),
                Outcome = "ADAMTS-5") %>%
  dplyr::mutate(Exposure = ifelse(pval < 0.05, paste0(Exposure,"*"),Exposure),
                "N variants" = nsnp,
                'P Value' = round(as.numeric(pval),3),
                OR = round(or,2),
                lor = round(or_lci95,2),
                uor = round(or_uci95,2)) %>%
  mutate("95%CI" = paste0(lor,
                          " - ",uor),
         " " = paste(rep(" ",12),collapse = " "))


tm <- forest_theme(ci_col = "#008B45FF",refline_col = "#EE0000FF")

#pdf(file = "image/TIMP3vsADAMTS-5.pdf",width = 8,height = 4)
p2 <- forest(mr_tabledecode[,c(15,16,17,19,22,23,18)],est = mr_tabledecode$OR,lower = mr_tabledecode$lor,
             upper = mr_tabledecode$uor,ci_column = 6,xlim = c(0.8,1.2),
             ticks_at = c(0.8,1,1.2),ref_line = 1,theme = tm)

#dev.off()
p2 <- add_border(p2, part = "header", row = c(1,2), where = "top") 
p2 <- edit_plot(p2,row = which(mr_tabledecode$pval <0.05),
                gp = grid::gpar(col = "red", fontface = "bold"))


p2
tiff(filename = "image/ADAMST5.tiff",width = 2000,height = 1600,res = 300)
p2
dev.off()

pdf(file = "image/ADAMST5.pdf",width = 8,height = 10)
p2

dev.off()

###异质性检验
ADAMTS_heter <- mr_heterogeneity(dat_ADAMTS)

###水平性检验
ADAMTS_plei <- mr_pleiotropy_test(dat_ADAMTS)

### 方向性检验 steiger
ADAMTS_steiger <- directionality_test(dat_ADAMTS)

### table S4
tableS4 <- ADAMTS_heter %>%
  dplyr::filter(method == "Inverse variance weighted") %>%
  dplyr::select(exposure,Q,Q_pval) %>%
  dplyr::mutate(.,Q = sprintf("%.2f",Q)) %>%
  left_join(.,dplyr::select(ADAMTS_plei,exposure,pval)) %>%
  left_join(.,dplyr::select(ADAMTS_steiger,exposure,correct_causal_direction,steiger_pval)) %>%
  left_join(.,presso_ADAMTSres) %>%
  dplyr::mutate(.,exposure = Hmisc::capitalize(tolower(unlist(lapply(exposure, 
                                                                     function(i){strsplit(i,split = " ")[[1]][1]}))))) 
colnames(tableS4)[1:7] <- c("Cancer","Q static","P value (Heteroeneity test)","P value (Horizontal pleiotropy test)",
                            "Correct steiger direction","P value (Steiger test)","P value (MR-PRESSO globtal test)")

write.csv(tableS4,file = "result/tableS4.csv",quote = F,row.names = F)

