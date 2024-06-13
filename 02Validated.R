
rm(list = ls())
##ieu暴露
expoure_ieu <- readRDS("data/ADAMTS5_ieu.rds")

expoure_ieu <- expoure_ieu %>%
  dplyr::mutate(.,R2 = (beta.exposure^2)/(beta.exposure^2 + samplesize.exposure*se.exposure^2),
                Fstatic =  ((samplesize.exposure - 1 - 1)/1)*(R2/(1 - R2)))


##结局
outcome <- readRDS("data/canceroutcome.Rda") %>%
  subset(.,pval.outcome > 5e-08) %>%
  subset(.,! id.outcome %in% c("ieu-b-89","ieu-b-94","ukb-b-14956"))

###Harmonise          
##########deCODE
ieu_dat <- harmonise_data(exposure_dat=expoure_ieu, outcome_dat=outcome, action=2)


### sterger 过滤
ieu_dat <- steiger_filtering(ieu_dat) %>%
  subset(.,steiger_dir == TRUE)

####MR分析
mr_ieu <- mr(ieu_dat) %>%
  generate_odds_ratios()

# Figure 1S验证森林图
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

mr_tableieu <- mr_ieu %>%
  dplyr::filter(method %in% c("Wald ratio","Inverse variance weighted")) %>%
  dplyr::mutate(.,Cancer = Hmisc::capitalize(tolower(unlist(lapply(outcome, 
                                                                   function(i){strsplit(i,split = " ")[[1]][1]}))))) %>%
  rename(Odds = or,Low = or_lci95,High = or_uci95) %>% 
  dplyr::mutate(.,ci = sprintf("%.2f (%.2f, %.2f)",Odds,Low,High),
                pval = as.numeric(pval),
                ci = ifelse(is.na(pval),NA,ci),
                Cancer = ifelse(pval < 0.05, paste0(Cancer,"*"),Cancer),
                pval = ifelse(pval<0.001, scales::scientific(pval,digits = 3),sprintf("%.2f",pval)))  %>%
  add_title_for_table(.)

# Figure 1
tiff(filename = "image/forestdeieu.tiff",width = 800,height = 950,res = 120)
forestplot(labeltext = as.matrix(mr_tableieu[,c("Cancer","ci","pval")]),
           mean = mr_tableieu$Odds,lower = mr_tableieu$Low,upper = mr_tableieu$High,
           align = "l",graph.pos = 2, zero = 1,graphwidth = unit(5,'cm'),is.summary = c(T,rep(F,19)),
           hrzl_lines = list("2" = gpar(lty = 1,col = "black")),
           col = fpColors(line = "blue",box = "#008B45FF"),
           txt_gp=fpTxtGp(label=gpar(cex=0.95),ticks=gpar(cex=0.85),xlab=gpar(cex = 0.85)),
           boxsize=0.15,xticks = c(0.7,1,1.3))

dev.off()

pdf(file = "image/forestdeieu.pdf",width = 8,height = 10)
forestplot(labeltext = as.matrix(mr_tableieu[,c("Cancer","ci","pval")]),
           mean = mr_tableieu$Odds,lower = mr_tableieu$Low,upper = mr_tableieu$High,
           align = "l",graph.pos = 2, zero = 1,graphwidth = unit(5,'cm'),is.summary = c(T,rep(F,19)),
           hrzl_lines = list("2" = gpar(lty = 1,col = "black")),
           col = fpColors(line = "blue",box = "#008B45FF"),
           txt_gp=fpTxtGp(label=gpar(cex=0.95),ticks=gpar(cex=0.85),xlab=gpar(cex = 0.85)),
           boxsize=0.15,xticks = c(0.7,1,1.3))

dev.off()
