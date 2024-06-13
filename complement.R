rm(list = ls())

exposure <- data.table::fread("data/2480_58_TIMP3_TIMP_3.txt.gz")%>% 
  subset(Pval < 5e-08) %>% 
  dplyr::filter(!is.na(rsids)) %>% 
  dplyr::mutate(phenotype = "TIMP3",
                id = "deCODE") %>%
  TwoSampleMR::format_data(.,type = "exposure",phenotype_col = "phenotype",
                           beta_col = "Beta",se_col = "SE",id_col = "id",snp_col = "rsids",
                           pval_col = "Pval",effect_allele_col = "effectAllele",
                           other_allele_col = "otherAllele", chr_col = "Chrom",
                           eaf_col = "ImpMAF",pos_col = "Pos",samplesize_col = "N") %>%
  clump_data(.)

repeat{
  try({
    outcome <- extract_outcome_data(snps = exposure$SNP,proxies = T,rsq = 0.8,
                                    outcomes = c("ieu-a-1120","ieu-b-96",
                                                 "ieu-a-1129","ieu-a-966",
                                                 "ieu-a-822","ieu-b-85"))
  })
  if(exists("outcome")){break}
  Sys.sleep(2)
}

files <- list.files("data",pattern = "finngen_R")

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


outcome_finn <- lapply(files, function(i){readfinne(i,SNP =exposure$SNP)})
outcome_finn <- do.call(rbind,outcome_finn)


outcome <- outcome_finn %>% as.data.frame() %>%
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




##Harmonise
decode_dat <- harmonise_data(exposure, outcome, action=2) %>%
  subset(.,pval.outcome > 5e-08) %>%
  steiger_filtering() %>%
  subset(.,steiger_dir == TRUE)

## MR分析
mr_decode <- mr(decode_dat) %>%
  generate_odds_ratios()

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
  dplyr::rename(Odds = or,Low = or_lci95,High = or_uci95) %>% 
  dplyr::mutate(.,ci = sprintf("%.2f (%.2f, %.2f)",Odds,Low,High),
                pval = round(as.numeric(pval),4),
                ci = ifelse(is.na(pval),NA,ci),
                Cancer = ifelse(pval < 0.05, paste0(Cancer,"*"),Cancer))  %>%
  add_title_for_table(.)

library(forestplot)

#tiff(filename = "image/TIMP3.tiff",width = 800,height = 950,res = 120)
forestplot(labeltext = as.matrix(mr_tabledecode[,c("Cancer","nsnp","ci","pval")]),
           mean = mr_tabledecode$Odds,lower = mr_tabledecode$Low,upper = mr_tabledecode$High,
           align = "l",graph.pos = 3, zero = 1,graphwidth = unit(5,'cm'),is.summary = c(T,rep(F,19)),
           hrzl_lines = list("2" = gpar(lty = 1,col = "black")),
           col = fpColors(line = "blue",box = "#008B45FF"),
           txt_gp=fpTxtGp(label=gpar(cex=0.95),ticks=gpar(cex=0.85),xlab=gpar(cex = 0.85)),
           boxsize=0.15)
#dev.off()
pdf(file = "image/forestTIMP3.pdf",width = 8,height = 10)
forestplot(labeltext = as.matrix(mr_tabledecode[,c("Cancer","nsnp","ci","pval")]),
           mean = mr_tabledecode$Odds,lower = mr_tabledecode$Low,upper = mr_tabledecode$High,
           align = "l",graph.pos = 3, zero = 1,graphwidth = unit(5,'cm'),is.summary = c(T,rep(F,19)),
           hrzl_lines = list("2" = gpar(lty = 1,col = "black")),
           col = fpColors(line = "blue",box = "#008B45FF"),
           txt_gp=fpTxtGp(label=gpar(cex=0.95),ticks=gpar(cex=0.85),xlab=gpar(cex = 0.85)),
           boxsize=0.15)

dev.off()


###################TIMP3.tiff on risk ADAMTS-5
outcome <- data.table::fread("data/3168_8_ADAMTS5_ADAMTS_5.txt.gz") %>%
  dplyr::filter(Pval >= 5e-08) %>%
  dplyr::filter(!is.na(rsids)) %>% 
  dplyr::mutate(phenotype = "ADAMST-5",
                id = "deCODE") %>%
  TwoSampleMR::format_data(.,type = "outcome",snp = exposure$SNP,phenotype_col = "phenotype",
                           beta_col = "Beta",se_col = "SE",id_col = "id",snp_col = "rsids",
                           pval_col = "Pval",effect_allele_col = "effectAllele",
                           other_allele_col = "otherAllele", chr_col = "Chrom",
                           eaf_col = "ImpMAF",pos_col = "Pos",samplesize_col = "N")

dat <- harmonise_data(exposure, outcome, action=2) %>%
  steiger_filtering() %>%
  subset(.,steiger_dir == TRUE)

## MR analysis
mr_res <- mr(dat) %>%
  generate_odds_ratios()

library(forestploter)
forestdat <- data.frame(Method = mr_res$method,
                        Nsnp = c(12,"","","",""),
                        OR = round(mr_res$or,3),
                        lor = round(mr_res$or_lci95,3),
                        uor = round(mr_res$or_uci95,3)) %>%
  mutate(CI = paste0(OR," (",lor,
                      " - ",uor,")"),
          Index = paste(rep(" ",12),collapse = " "),
          #P = format(mr_res$pval,scientific=T,digits=3)
          P = round(mr_res$pval,3)) %>%
  mutate("P value" = ifelse(P < 0.05,paste0(P,"*"),P))

colnames(forestdat) <- c("Method","N variants",
                         "OR","lor","uor","OR (95%CI)"," ","P","P Value")

tm <- forest_theme(ci_col = "#008B45FF",refline_col = "#EE0000FF")

pdf(file = "image/TIMP3vsADAMTS-5.pdf",width = 8,height = 4)
forest(forestdat[,c(1,2,3,6,7,9)],est = forestdat$OR,lower = forestdat$lor,title = "A",
       upper = forestdat$uor,ci_column = 5,xlim = c(0.8,1.2),
       ticks_at = c(0.8,1,1.2),
       ref_line = 1,theme = tm)

dev.off()

