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
#presso_ADAMTSres[9,"Number of Outlier SNPs"] <- 0

## MR分析mr(Mydata, method_list=c("mr_ivw", "mr_ivw_fe", "mr_two_sample_ml","mr_egger_regression", "mr_weighted_median", "mr_penalised_weighted_median","mr_simple_mode", "mr_weighted_mode"))
mr_ADAMTS <- mr(dat_ADAMTS) %>%
  generate_odds_ratios()

breastmr <-  mr(dplyr::filter(dat_ADAMTS,exposure == "Breast"),method_list = c("mr_ivw_mre")) %>%
  generate_odds_ratios()

# figure2
add_title_for_table <- function(df)
{
  tibble(
    id.exposure=NA,id.outcome = NA,outcome = NA,exposure = NA,
    method = "Method",b = NA, se = NA, lo_ci = NA, up_ci = NA,
    Odds = NA, Low = NA, High =NA, ci="OR (95%CI)",nsnp = "SNP",
    pval="P", pve= "PVE",Cancer = "Cancer") %>% 
    dplyr::select(all_of(names(df))) %>% 
    rbind(df)
}


mr_tableADAMTS <- mr_ADAMTS %>%
  rbind(breastmr) %>%
  dplyr::filter(method %in% c("Wald ratio","Inverse variance weighted")) %>%
  dplyr::mutate(.,Cancer = Hmisc::capitalize(tolower(unlist(lapply(exposure, 
                                                                   function(i){strsplit(i,split = " ")[[1]][1]}))))) %>%
  rename(Odds = or,Low = or_lci95,High = or_uci95) %>% 
  dplyr::mutate(.,ci = sprintf("%.2f (%.2f, %.2f)",Odds,Low,High),
                pval = round(as.numeric(pval),3),
                ci = ifelse(is.na(pval),NA,ci),
                Cancer = ifelse(pval < 0.05, paste0(Cancer,"*"),Cancer))  %>%
  add_title_for_table(.)


tiff(filename = "image/ADAMST5.tiff",width = 800,height = 950,res = 120)
forestplot(labeltext = as.matrix(mr_tableADAMTS[,c("Cancer","nsnp","ci","pval")]),
           mean = mr_tableADAMTS$Odds,lower = mr_tableADAMTS$Low,upper = mr_tableADAMTS$High,
           align = "l",graph.pos = 3, zero = 1,graphwidth = unit(5,'cm'),is.summary = c(T,rep(F,19)),
           hrzl_lines = list("2" = gpar(lty = 1,col = "black")),
           col = fpColors(line = "blue",box = "#008B45FF"),
           txt_gp=fpTxtGp(label=gpar(cex=0.95),ticks=gpar(cex=0.85),xlab=gpar(cex = 0.85)),
           boxsize=0.15,xticks = c(0.9,1,1.1))
dev.off()


pdf(file = "image/ADAMST5.pdf",width = 8,height = 10)
forestplot(labeltext = as.matrix(mr_tableADAMTS[,c("Cancer","nsnp","ci","pval")]),
           mean = mr_tableADAMTS$Odds,lower = mr_tableADAMTS$Low,upper = mr_tableADAMTS$High,
           align = "l",graph.pos = 3, zero = 1,graphwidth = unit(5,'cm'),is.summary = c(T,rep(F,19)),
           hrzl_lines = list("2" = gpar(lty = 1,col = "black")),
           col = fpColors(line = "blue",box = "#008B45FF"),
           txt_gp=fpTxtGp(label=gpar(cex=0.95),ticks=gpar(cex=0.85),xlab=gpar(cex = 0.85)),
           boxsize=0.15,xticks = c(0.9,1,1.1))
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

