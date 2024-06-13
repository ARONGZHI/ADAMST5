
rm(list = ls())
##ieu暴露
expoure_ieu <- readRDS("data/ADAMTS5_ieu.rds")

expoure_ieu <- expoure_ieu %>%
  dplyr::mutate(.,R2 = (beta.exposure^2)/(beta.exposure^2 + samplesize.exposure*se.exposure^2),
                Fstatic =  ((samplesize.exposure - 1 - 1)/1)*(R2/(1 - R2)))


##结局
outcome <- readRDS("data/canceroutcome.rds") %>%
  subset(.,pval.outcome > 5e-08) 

###Harmonise          
##########deCODE
ieu_dat <- harmonise_data(exposure_dat=expoure_ieu, outcome_dat=outcome, action=2)


### sterger 过滤
ieu_dat <- steiger_filtering(ieu_dat) %>%
  subset(.,steiger_dir == TRUE)

####MR分析
mr_ieu <- mr(ieu_dat) %>%
  generate_odds_ratios()


# Figure 3
mr_tabledecode <- mr_ieu %>%
  dplyr::filter(method %in% c("Wald ratio","Inverse variance weighted")) %>%
  dplyr::mutate(.,Exposure = "ADAMST-5",
                Outcome = Hmisc::capitalize(tolower(unlist(lapply(outcome, 
                                                                  function(i){strsplit(i,split = " ")[[1]][1]}))))) %>%
  dplyr::mutate(Outcome = ifelse(pval < 0.05, paste0(Outcome,"*"),Outcome),
                "N variants" = nsnp,
                'P Value' = round(as.numeric(pval),4),
                OR = round(or,2),
                lor = round(or_lci95,2),
                uor = round(or_uci95,2)) %>%
  mutate("95%CI" = paste0(lor,
                          " - ",uor),
         " " = paste(rep(" ",12),collapse = " "))


tm <- forest_theme(ci_col = "#008B45FF",refline_col = "#EE0000FF")

p3 <- forest(mr_tabledecode[,c(15,16,17,19,22,23,18)],est = mr_tabledecode$OR,lower = mr_tabledecode$lor,
             upper = mr_tabledecode$uor,ci_column = 6,xlim = c(0.6,1.4),
             ticks_at = c(0.6,1,1.4),ref_line = 1,theme = tm)
p3

p3 <- add_border(p3, part = "header", row = c(1,2), where = "top") 


p3 <- edit_plot(p3,row = which(mr_tabledecode$pval <0.05),
                gp = grid::gpar(col = "red", fontface = "bold"))

p3

tiff(filename = "image/forestdeieu.tiff",width = 2000,height = 1600,res = 300)
p3
dev.off()

pdf(file = "image/forestdeieu.pdf",width = 8,height = 10)
p3
dev.off()
