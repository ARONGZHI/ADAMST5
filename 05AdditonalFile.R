
# TableS1
id <- readRDS("data/cancer_exposure.rds") %>%
  dplyr::filter(id.exposure != "finngen_R9_C3_COLORECTAL_EXALLC") %>%
  dplyr::filter(grepl("finngen",id.exposure))
id <- unique(id$id.exposure) %>% sub("finngen_R9_","",.)

finnmanifest <- data.table::fread("data/summary_stats_R9_manifest.tsv")[,c(1,4,5)] %>%
  dplyr::filter(phenocode %in% id) %>%
  dplyr::mutate(.,"Sample size" = num_cases + num_controls,
                "Source" = "https://www.finngen.fi",
                "Consortium" =  'FinnGen',
                "Ancestry" = "European") %>%
  dplyr::mutate(.,phenocode = Hmisc::capitalize(tolower(unlist(lapply(phenocode, 
                                                                      function(i){strsplit(i,split = "_")[[1]][2]}))))) %>%
  
  dplyr::mutate(.,phenocode = paste0(phenocode, " cancer")) 

colnames(finnmanifest)[1:3] <- c("Traits","Number of cases","Numben of controls")


library(TwoSampleMR)
ao <- available_outcomes()
ieumanifest <- ao[ao$id %in% c("ieu-a-1120", "ieu-b-96","ieu-a-1129","ieu-a-966","ieu-a-822",
                               "ieu-b-85",'prot-c-3168_8_2'),][,c(2,3,7,9,12,10,14)] %>%
  dplyr::mutate(pmid = paste0("PMID ",pmid)) %>%
  rename("Traits" = trait,"Number of cases" = "ncase","Consortium" = consortium,"Numben of controls" = ncontrol,
         "Sample size" = sample_size,"Ancestry" = population,"Source" = pmid)
ieumanifest[3,"Sample size"] <- ieumanifest[3,"Numben of controls"] +ieumanifest[3,"Number of cases"]



manifest <- tibble::tibble("Traits" = "ADAMTS-5 deCDDE","Consortium"='deCODE',`Number of cases`=NA,"Numben of controls" =NA,
                           "Sample size" = 35559,"Ancestry" = "European","Source" = "https://www.decode.com")  %>%
  rbind(.,dplyr::select(ieumanifest,colnames(.))[4,]) %>%
  rbind(.,dplyr::select(ieumanifest,colnames(.))[-4,]) %>%
  rbind(.,dplyr::select(finnmanifest,colnames(.))) 
manifest[2,"Sample size"] <- 997
manifest[2,"Consortium"] <- "IEU Open GWAS project"

write.csv(manifest,file = "result/tableS1.csv")


