
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
  subset(.,pval.outcome > 5e-08) 

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

library(forestploter)

mr_tabledecode <- mr_decode %>%
  dplyr::filter(method %in% c("Wald ratio","Inverse variance weighted")) %>%
  dplyr::mutate(.,Exposure = "ADAMST-5",
                Outcome = Hmisc::capitalize(tolower(unlist(lapply(outcome, 
                                                                   function(i){strsplit(i,split = " ")[[1]][1]}))))) %>%
  dplyr::mutate(Outcome = ifelse(pval < 0.05, paste0(Outcome,"*"),Outcome),
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
p1 <- forest(mr_tabledecode[,c(15,16,17,19,22,23,18)],est = mr_tabledecode$OR,lower = mr_tabledecode$lor,
       upper = mr_tabledecode$uor,ci_column = 6,xlim = c(0.2,1.8),
       ticks_at = c(0.2,1,1.8),ref_line = 1,theme = tm)

#dev.off()

p1 <- add_border(p1, part = "header", row = c(1,2), where = "top") 


p1 <- edit_plot(p1,row = which(mr_tabledecode$pval <0.05),
          gp = grid::gpar(col = "red", fontface = "bold"))

p1
tiff(filename = "image/forestdeCODE.tiff",width = 2000,height = 1600,res = 300)
p1
dev.off()

pdf(file = "image/forestdeCODE.pdf",width = 8,height = 10)
p1

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