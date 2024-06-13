
#############deCOED

## 暴露
rm(list = ls())
#ao <- available_outcomes()

##decode
expoure_decode <- readRDS("data/ADAMTS5_deCODE.rds")

expoure_decode <- expoure_decode %>%
  dplyr::mutate(.,#R2 = 2*eaf.exposure*(1 - eaf.exposure) *beta.exposure^2,
                R2 = (beta.exposure^2)/(beta.exposure^2 + samplesize.exposure*se.exposure^2),
                Fstatic =  ((samplesize.exposure - 1 - 1)/1)*(R2/(1 - R2)))

tableS2 <- expoure_decode %>%
  select("SNP","effect_allele.exposure","other_allele.exposure","beta.exposure",
         "se.exposure", "pval.exposure") %>%
  rename("SNP" = SNP,"Effect Allele" = effect_allele.exposure,"Other Allele" = other_allele.exposure,
         "BETA" = beta.exposure,"SE" = se.exposure,"P" = pval.exposure)

write.csv(tableS2,file = "result/TableS2.csv",quote = F,row.names = F)

##结局
outcome <- readRDS("data/canceroutcome.rds") %>%
  subset(.,pval.outcome > 5e-08) %>%
  subset(.,! id.outcome %in% c("ieu-b-89","ieu-b-94","ukb-b-14956","finngen_R9_C3_COLORECTAL_EXALLC"))

##Harmonise
decode_dat <- harmonise_data(exposure_dat=expoure_decode, outcome_dat=outcome, action=2)
##finn

## sterger 过滤
decode_dat  <- steiger_filtering(decode_dat) %>%
  subset(.,steiger_dir == TRUE)

## Supplement S1  SNPs adopted as genetic instruments in mendelian randomization analysis of ADAMTS-5 on different types of cancer
S1 <- dplyr::select(decode_dat,SNP,effect_allele.exposure,other_allele.exposure,effect_allele.outcome,
                    other_allele.outcome,exposure,outcome,eaf.exposure,beta.exposure,se.exposure,pval.exposure,
                    eaf.outcome,beta.outcome,se.outcome,pval.outcome,proxy.outcome,proxy_snp.outcome,proxy_a1.outcome,proxy_a2.outcome,R2,Fstatic,steiger_dir)

write.csv(S1,file = "result/S1.csv",quote = F)


## MR分析
mr_decode <- mr(decode_dat) %>%
  generate_odds_ratios()

## Figure 1绘制结果森林图
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

library(tibble)
library(tidyverse)
mr_tabledecode <- mr_decode %>%
  dplyr::filter(method %in% c("Wald ratio","Inverse variance weighted")) %>%
  dplyr::mutate(.,Cancer = Hmisc::capitalize(tolower(unlist(lapply(outcome, 
                                                                   function(i){strsplit(i,split = " ")[[1]][1]}))))) %>%
  rename(Odds = or,Low = or_lci95,High = or_uci95) %>% 
  dplyr::mutate(.,ci = sprintf("%.2f (%.2f, %.2f)",Odds,Low,High),
                pval = round(as.numeric(pval),4),
                ci = ifelse(is.na(pval),NA,ci),
                Cancer = ifelse(pval < 0.05, paste0(Cancer,"*"),Cancer))  %>%
  add_title_for_table(.)

library(forestplot)

tiff(filename = "image/forestdeCODE.tiff",width = 800,height = 950,res = 120)
forestplot(labeltext = as.matrix(mr_tabledecode[,c("Cancer","nsnp","ci","pval")]),
           mean = mr_tabledecode$Odds,lower = mr_tabledecode$Low,upper = mr_tabledecode$High,
           align = "l",graph.pos = 3, zero = 1,graphwidth = unit(5,'cm'),is.summary = c(T,rep(F,19)),
           hrzl_lines = list("2" = gpar(lty = 1,col = "black")),
           col = fpColors(line = "blue",box = "#008B45FF"),
           txt_gp=fpTxtGp(label=gpar(cex=0.95),ticks=gpar(cex=0.85),xlab=gpar(cex = 0.85)),
           boxsize=0.15)
dev.off()

pdf(file = "image/forestdeCODE.pdf",width = 8,height = 10)
forestplot(labeltext = as.matrix(mr_tabledecode[,c("Cancer","nsnp","ci","pval")]),
           mean = mr_tabledecode$Odds,lower = mr_tabledecode$Low,upper = mr_tabledecode$High,
           align = "l",graph.pos = 3, zero = 1,graphwidth = unit(5,'cm'),is.summary = c(T,rep(F,19)),
           hrzl_lines = list("2" = gpar(lty = 1,col = "black")),
           col = fpColors(line = "blue",box = "#008B45FF"),
           txt_gp=fpTxtGp(label=gpar(cex=0.95),ticks=gpar(cex=0.85),xlab=gpar(cex = 0.85)),
           boxsize=0.15)

dev.off()

###异质性检验
decode_heter <- mr_heterogeneity(decode_dat)

###水平性检验
decode_plei <- mr_pleiotropy_test(decode_dat)

### 方向性检验 steiger
decode_steiger <- directionality_test(decode_dat)

######MRPRESSO 异常值检验
library(MRPRESSO)

presso_modified <- function(id,dat){
  tmp <- dplyr::filter(dat,id.outcome == id )
  presso <-mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", 
                     OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data= tmp, NbDistribution = 1000,  SignifThreshold = 0.05,seed=121)
}

presso_decode <- lapply(unique(decode_dat$id.outcome)[-14], function(i) {
  presso_modified(i,decode_dat)})

presso_res <- data.frame(id.outcome = unique(decode_dat$id.outcome) [-15],
                         pvalue = unlist(lapply(1:length(presso_decode),function(i){
                           presso_decode[[i]]$`MR-PRESSO results`$`Global Test`$Pvalue})))

#saveRDS(presso_res,file = "result/deCODEpresso.rds")
#presso_res <- readRDS("result/deCODEpresso.rds")
### table S3 
R2 <- decode_dat %>%
  dplyr::select(id.outcome, beta.exposure, se.exposure, samplesize.exposure) %>% 
  dplyr::group_by(id.outcome) %>% 
  dplyr::summarise(R2 = sum((beta.exposure^2)/(beta.exposure^2 + samplesize.exposure*se.exposure^2))) %>%
  dplyr::select(id.outcome,R2) %>%
  dplyr::mutate(.,R2 = paste0(sprintf("%.2f",100*R2),"%")) 

tableS3 <- decode_heter %>%
  dplyr::filter(method == "Inverse variance weighted") %>%
  dplyr::select(id.outcome,outcome,Q,Q_pval) %>%
  inner_join(.,dplyr::select(decode_plei,id.outcome,pval)) %>%
  inner_join(.,dplyr::select(decode_steiger,id.outcome,correct_causal_direction,steiger_pval)) %>%
  dplyr::mutate(.,outcome = Hmisc::capitalize(tolower(unlist(lapply(outcome, 
                                                                    function(i){strsplit(i,split = " ")[[1]][1]}))))) %>%
  left_join(.,presso_res) %>%
  rename("Cancer" = outcome,"Q static" = Q,"P value (Heteroeneity test)" = Q_pval,
         "P value (Horizontal pleiotropy test)" = pval,
         "Correct steiger direction" = correct_causal_direction,
         "P value (Steiger test)" = steiger_pval,
         "P value (MR-PRESSO globtal test)" = pvalue)
tableS3 <- inner_join(tableS3,R2) %>%
  dplyr::select(-id.outcome)


write.csv(tableS3,file = "result/tableS3.csv",quote = F,row.names = F)


res_single <- mr_leaveoneout(dat = dplyr::filter(decode_dat,outcome %in% c("OESOPHAGUS","OROPHARYNGEAL")))

mr_leaveoneout_plot(res_single)