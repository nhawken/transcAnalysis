setwd("C:/Users/natha/Documents/Geschwind Rotation 2020/transcAnalysis/")
options(stringsAsFactors = FALSE)

library(tidyverse)
library(RColorBrewer)
library(qqman)
library(patchwork)
library(corrplot)
library(pheatmap)
library(gridExtra)
library(UpSetR)
library(grid)
library(igraph)
library(ggraph)


file_index <- 1
parent_folder  <-  getwd()
# # output_folder  <-  file.path(parent_folder, "LME","dx_day")
output_folder = paste(getwd(), '/analyze_transcLME_OUT_',Sys.Date(),'/', sep = "") # SET AN OUTPUT PATH HERE
dir.create(output_folder, showWarnings = F, recursive = T)
#time0 <- proc.time()

theme_set(theme_classic(base_size = 18))

load("all_transcript_data.Rdata")
load(file.path(parent_folder,"input_for_DE_transc.rdata"))

cond <- "dx_day"

# source("~/Aaron/Scripts/Aaron/CIRM/3_Diagnosis/DE_functions.R")

# read in DE results
de_results <- read_csv(paste0(parent_folder, "/",cond, "_lme_DE.csv"))[,-2] %>% #calls all but the 2nd col of the matrix
  rename("ensembl_transcript_id" = "X1") %>%  #renames unlabeled column
  setNames(gsub("Idiopathic(.?ASD)?","IdiopathicASD", colnames(.))) #rename idiopathic ASD columns 



dxs <- gsub(".*\\.(.*?)_.*","\\1", colnames(de_results)) %>% #gsub one or more of any charater followed by a period, (.*?) is the capture group
  unique() %>%  #gets rid of Beta/p/SE. in front, keeps the next dx, gets rid of the _025/50/late afterwards
  sort() %>%
  .[. != "ensembl_transcript_id"] 

stages <- gsub(".*\\..*?_(.*)","\\1", colnames(de_results)) %>% #get rid of Beta/p/Se. in front then gets rid of Dx_, capture anything after Dx
  unique() %>% 
  sort() %>%
  .[. != "ensembl_transcript_id"] 

dx_stage <- gsub(".*\\.(.*?_.*)","\\1", colnames(de_results)) %>%  #get rid of Beta/p/SE., then capture dx and stages
  unique() %>% 
  sort() %>%
  .[. != "ensembl_transcript_id"] 


# # correct p value and add annotation
de_results_adj <- lapply(dx_stage, function(dx_stg){ #for each dx_stage, create a list with 1 data frame (tibble) containing beta, SE, p, fdr
  p_val_col <- paste0("p.",dx_stg)
  de_results %>% 
    dplyr::select(contains(dx_stg)) %>% 
    mutate(fdr =   !!as.name(p_val_col) %>% p.adjust(., method = "fdr")) %>% #control the false discovery rate, the expected proportion of false discoveries amongst the rejected hypotheses
    # mutate(data fram, new_var = [existing var])
    rename_at(vars(matches("fdr")), list(~paste0("fdr.",dx_stg)))
}) %>% bind_cols() %>%  #combines dataframes, list of dataframs
  mutate(ensembl_transcript_id = de_results$ensembl_transcript_id) %>% #add new column to df with the ensembl tanscript ids
  dplyr::select(ensembl_transcript_id, everything()) #puts all other columns after the ID column



transAnno <-  transcriptAnnoRaw %>%
  filter(ensembl_transcript_id %in% rownames(datTrans_reg_batch)) %>%
  filter(!duplicated(ensembl_transcript_id))

# separate dx_day into seperate dfs
de_dx_day <- lapply(dx_stage, function(dx_stg){
  de_results_adj %>% 
    dplyr::select(ensembl_transcript_id,contains(dx_stg))  %>% #keep only the variables ID, and dx_stage of interest (beta, se, p, fdr)
    setNames(gsub(paste0("\\.",dx_stg), "", colnames(.))) %>% #remove .dx_stage identifier from the columns of beta, se, p, fdr
    left_join(transAnno, by = "ensembl_transcript_id") #add transcript annotation table to the right edge of the table
})  %>% setNames(dx_stage) #rename overall df to be the dx_stage identifier


# Plot number of DEG -------------------------------------------------------------
deg_n <- lapply(dx_stage, function(dx_stg){
  res <- de_dx_day[[dx_stg]] #call df with specific dx_stage id
  stage = gsub(".*?_(.*)","\\1", dx_stg) #pull out stage from id
  dx = gsub("(.*?)_.*","\\1", dx_stg) #pull out dx from id
  
  sigRes1 <- res %>% 
    filter(p < 0.005) #grab only p vals less than 0.005
  sigGenSum1 <- c(stage,dx, paste("p =",0.005),sum(sigRes1$beta > 0),sum(sigRes1$beta < 0)) #number of downreg and upreg values
  
  sigRes2 <- res %>% 
    filter(fdr < 0.05) #grab only fdr vals less than 0.05
  sigGenSum2 <- c(stage,dx, "FDR = 0.05",sum(sigRes2$beta > 0),sum(sigRes2$beta < 0)) #number of downreg and upreg values
  
  sigGenSum <- rbind(sigGenSum1,sigGenSum2) %>% 
    as.data.frame()
  colnames(sigGenSum) <- c("stage","Dx","pvalue", "upRegulated", "downRegulated")
  sigGenSum$upRegulated = as.numeric(sigGenSum$upRegulated)
  sigGenSum$downRegulated = -1*as.numeric(sigGenSum$downRegulated)
  sigGenSum <- sigGenSum %>%  
    gather( key = "variable", value = "value", upRegulated:downRegulated) 
}) %>% 
  do.call(rbind,.)


pdf(file = paste0(output_folder, "/", file_index, "_LME_pvalues.pdf"), height = 7, width = 14)
deg_n_plots <- lapply(unique(deg_n$pvalue), function(p){
  plot_dat <- deg_n %>% 
    filter(pvalue == p) %>% 
    mutate(stage = factor(stage, levels = rev(c("025", "050", "075", "100", "late_stage", "all_timepoints")))) %>% 
    mutate(stg_n = as.numeric(stage)) %>% 
    mutate(stg_n = case_when(stg_n == 1 ~ 1,
                             stg_n == 2 ~ 2, 
                             stg_n %in% c(3:6) ~ stg_n + 0.5
    ))
  limits <- range(deg_n$value) + c(-1,1) * 220
  browser()
  # p1 <- ggplot(plot_dat, aes(x = Dx, y = value, fill = variable)) +
  #   geom_bar(stat = "identity", position = "identity") +
  #   facet_grid(Dx ~ stage, scales = "free") +
  #   scale_fill_manual(values = c("#77dd77", "#ff7878")) +
  #   geom_hline(yintercept = 0, size = 1) +
  #   scale_y_continuous(breaks = NULL,name = "",limits = limits) +
  #   theme_minimal() +
  #   theme(legend.position = "none",
  #         axis.line = element_blank(),
  #         axis.text.y = element_text(face = "bold"),
  #         axis.title = element_blank(),
  #         strip.text = element_text(size = 12, face = "bold"),
  #         strip.background = element_rect(fill = "grey90", color = "white"),
  #         panel.grid.major.y = element_blank(),
  #         strip.text.y = element_blank()) +
  #   coord_flip() +
  #   geom_text(data = plot_dat[plot_dat$variable == "upRegulated",],aes(label = value),hjust = -0.1, size = 3) +
  #   geom_text(data = plot_dat[plot_dat$variable == "downRegulated",],aes(label = abs(value)),hjust = 1.1, size = 3)
  p1 <- ggplot(plot_dat, aes(x = stg_n, y = value, fill = variable)) +
    geom_col(position = "identity") +
    facet_wrap(~Dx)+
    scale_fill_manual(values = c("#77dd77", "#ff7878")) +
    geom_hline(yintercept = 0, size = 1) +
    scale_y_continuous(breaks = NULL,name = "",limits = limits) +
    scale_x_continuous(breaks=sort(unique(plot_dat$stg_n)), labels  = rev(c("025", "050", "075", "100","late_stage", "all_timepoints"))) +
    coord_flip() +   
    theme_minimal() +
    theme(legend.position = "none",
          panel.grid = element_blank(),
          axis.line = element_blank(),
          axis.text.y = element_text(face = "bold"),
          axis.title = element_blank(),
          strip.text = element_text(size = 12, face = "bold"),
          strip.background = element_rect(fill = "grey90", color = "white"),
          panel.grid.major = element_blank(),
          strip.text.y = element_blank()) +

    geom_text(data = plot_dat[plot_dat$variable == "upRegulated",],aes(label = value),hjust = -0.1, size = 3) +
    geom_text(data = plot_dat[plot_dat$variable == "downRegulated",],aes(label = abs(value)),hjust = 1.1, size = 3)

  print(p1)
})

dev.off()

# plot P value histograms
melted_p_values <- de_results_adj %>% 
  select(ensembl_transcript_id, starts_with("p.")) %>% 
  gather("dx_day","pvalue", starts_with("p.")) %>%
  mutate(Dx = gsub("p.(.*?)_(.*)","\\1", dx_day), #grabs dx only
         Stage = gsub("p.(.*?)_(.*)","\\2", dx_day)) #grabs stage only

pdf(file = paste0(output_folder, "/", file_index, ".1_LME_pvalues_histogram.pdf"), height = 7, width = 12)
ggplot(melted_p_values, aes(x = pvalue)) +
  geom_histogram(binwidth = 0.05) +
  facet_grid(Stage ~ Dx) +
  theme_classic(base_size = 10) +
  theme(strip.text = element_text(size = 12)) +
  theme(axis.text.x = element_text(angle = -90, hjust  = 0, vjust = 0.5))

dev.off()

file_index = file_index + 1

# plot logFC around the CNV locus ---------------------------------------------------------------------------------
isoID <- paste(transAnno$hgnc_symbol, str_sub(transAnno$ensembl_transcript_id, -6, -1), sep = "_")
transAnno_id <- mutate(transAnno, hgnc_iso_symbol = isoID)

for (iso in 1:length(transAnno_id$ensembl_transcript_id)) {
  if (transAnno_id$hgnc_symbol[iso] == '') {
    transAnno_id$hgnc_symbol[iso] <- "N/A"
  }
}

loci <- list("15q13" = list("chr" = "15",
                            "surround" = c(28*10^6, 33*10^6),
                            "cnv" = c(30.61*10^6, 32.0*10^6)),
             "16p11del" = list("chr" = "16",
                               "surround" = c(28000000,31000000),
                               "cnv" = c(29.3*10^6,30.4*10^6)),
             "16p11dup" = list("chr" = "16",
                               "surround" = c(28000000,31000000),
                               "cnv" = c(29.3*10^6,30.4*10^6)),
             "22q11del" = list("chr" = "22",
                               "surround" = c(1, 2.3*10^7),
                               "cnv" = c(18.5*10^6,21.5*10^6)),
             "22q13del" = list("chr" = "22",
                               "surround" = c(44*10^6, 55*10^6),
                               "cnv" = c(46.77*10^6, 51.12*10^6)),
             "PCDH19" = list(chr = "X",
                             "surround" = c(90*10^6,101.5*10^6),
                             "gene" = "PCDH19"),
             "SHANK3" = list("chr" = "22",
                             "surround" = c(50*10^6, 55*10^6),
                             "gene" = "SHANK3"),
             "TS" = list("chr" = "12",
                         "surround" = c(1, 3*10^6),
                         "gene" = "CACNA1C"))


get_cnv_logFC <- function(deg,chr, plot_region, mut, type, stage){
  # deg - df from limma
  # chr - which chromosome
  # plot_region - vector of length 2 with start and end of region to be plotted
  # mut - either vector of length 2 with start and end of cnv or gene hgnc name
  # type - single gene or cnv
  if(type == "cnv"){
    dat_plot <- deg %>%
      filter(chromosome_name == chr & start_position >  plot_region[[1]] & end_position < plot_region[[2]]) %>% 
      mutate(inCNV = ifelse(end_position > mut[[1]] & start_position < mut[[2]],T,F),
             significant = ifelse(p < 0.005,T,F), 
             hgnc_iso_symbol = ifelse(hgnc_iso_symbol == "", ensembl_transcript_id, hgnc_iso_symbol)) %>% 
      arrange(start_position) %>% 
      mutate(hgnc_iso_symbol = factor(hgnc_iso_symbol, levels = hgnc_iso_symbol))
    
  } else if (type == "gene"){
    dat_plot <- deg %>%
      filter(chromosome_name == chr & start_position >  plot_region[[1]] & end_position < plot_region[[2]]) %>% 
      mutate(inCNV = ifelse(hgnc_symbol == mut,T,F),
             significant = ifelse(p < 0.005,T,F), 
             hgnc_iso_symbol = ifelse(hgnc_iso_symbol == "", ensembl_transcript_id, hgnc_iso_symbol)) %>% 
      arrange(start_position) %>% 
      mutate(hgnc_iso_symbol = factor(hgnc_iso_symbol, levels = hgnc_iso_symbol))  
  } else {
    stop('type must be one of "cnv"/"gene"')
  }
  dat_plot
}


dx_stage_loci <- dx_stage[!grepl(".*Idiopathic.*", dx_stage )]



defined_forms <-  grep(paste(names(loci),collapse = "|"), dx_stage_loci, value = T) # only run this on dx found in the loci list
cnv_logfc_plots <- lapply(defined_forms,function(dx_stg){
  dx <- gsub("(.*?)_.*","\\1",dx_stg)
  stage <- gsub(".*?_(.*)","\\1",dx_stg)
  type <- names(loci[[dx]])[3]
  de_results_adj %>%
    left_join(transAnno_id, by = "ensembl_transcript_id") %>% 
    select(contains(dx_stg), colnames(transAnno_id))  %>% 
    setNames(gsub(paste0("\\.",dx_stg), "", colnames(.))) %>% 
    get_cnv_logFC(loci[[dx]][["chr"]], loci[[dx]][["surround"]], loci[[dx]][[3]],type, stage) %>%
    mutate(dx_day = dx_stg, 
           Dx = dx, 
           stage = stage) %>% 
    ggplot() +
    geom_hline(yintercept = 0, color = "red") +
    geom_pointrange(aes(x = hgnc_symbol, y = beta, ymin = beta-SE , ymax = beta+SE, shape = significant, color = inCNV), position = position_dodge2(width = 1, preserve = "single", padding = 0.8))  +
    theme_classic(base_size = 12) +
    scale_shape_manual(values = c(1, 8)) +
    scale_color_manual(values = c("grey40", "royalblue")) +
    #facet_wrap(~hgnc_symbol, strip.position = 'bottom', scales = "free_x" ) +   
    #geom_text(position, aes(x=hgnc_symbol, y=0, label = hgnc_symbol)) +
    theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5)) +
    labs(x = "")  +
    ggtitle(dx_stg)
}) %>% setNames(defined_forms) 

pdf(paste0(output_folder, "/", file_index, "_loci_logFC.pdf"), width = 14)
print(cnv_logfc_plots)
dev.off()
file_index <- file_index + 1

pdf(paste0(output_folder, "/", file_index, ".1_loci_logFC.pdf"), width = 20, height = 14)
cnv_logfc_plots_atps <- cnv_logfc_plots[grep("all", names(cnv_logfc_plots))] 
cnv_logfc_plots_atps <- lapply(names(cnv_logfc_plots_atps), function(dx){
  p1 <- cnv_logfc_plots_atps[[dx]] + 
    ggtitle(gsub("_all.*","", dx))
  if(grepl("TS", dx)){
    p1 <- p1 + theme(legend.position = "bottom")
  } else {
    p1 <-  p1 + theme(legend.position = "none")
  }
  return(p1+ theme(plot.title = element_text(size =24)))
})

(cnv_logfc_plots_atps[[4]]) /
  (cnv_logfc_plots_atps[[2]] | cnv_logfc_plots_atps[[3]]) /
  (cnv_logfc_plots_atps[[1]]  | cnv_logfc_plots_atps[[5]]) / 
  (cnv_logfc_plots_atps[[6]] | cnv_logfc_plots_atps[[8]] | cnv_logfc_plots_atps[[7]])
dev.off()
file_index <- file_index + 1

# manhatan plot of DEGs ---------------------------------------------------------------------------------------
dat_manhattan <- lapply(dx_stage, function(dx_stg){
  de_results_adj %>% 
    dplyr::select(ensembl_transcript_id, contains(dx_stg))  %>% 
    left_join(transAnno_id, by = "ensembl_transcript_id") %>% 
    setNames(gsub(paste0("\\.",dx_stg), "", colnames(.)))
  
})  %>% setNames(dx_stage)

pdf(paste0(output_folder,"/", file_index,"_manhattan_plots.pdf"), width = 10, height = 4)#, units = "in",res = 300)
for (dx_stg in dx_stage) {
  manaht_data <- dat_manhattan[[dx_stg]] %>% 
    dplyr::select(chromosome_name, start_position, p) %>% 
    drop_na() %>% 
    mutate(chromosome_name = plyr::mapvalues(chromosome_name,from = c("X","Y"),to = c(23:24))) %>% 
    mutate_all(funs(as.numeric)) %>% 
    mutate(p = if_else(p == 0, 10^-32, p))
  print(manhattan(manaht_data,chr = "chromosome_name",bp = "start_position", p = "p",
                  col = rep(brewer.pal(8, "Spectral"),24), suggestiveline = F, genomewideline = F, 
                  chrlabs = c(1:22,"X","Y"), cex.axis = 1.5, main = paste(dx_stg)))
}
dev.off()

file_index <- file_index + 1

# GO enrichment ---------------------------------------------------------------------------------------------------
# library(clusterProfiler)
dat_clusterProfiler <- de_dx_day 
# save and move to hoffman for analysis
#save(dat_clusterProfiler,  file = paste0(output_folder,"/clusterProfiler_input.rdata"))
#hoffman_folder <- file.path("/u/project/geschwind/aarong/CIRM_iPSC/analysis/rsem/3_Dx/", basename(parent_folder), "LME/dx_day/clusterProfiler" )

#system(paste("ssh aarong@hoffman2.idre.ucla.edu mkdir", hoffman_folder))

#system(paste("scp",file.path(getwd(),output_folder, "clusterProfiler_input.rdata"),  paste0("aarong@hoffman2.idre.ucla.edu:",hoffman_folder)))

# Run clusterProfiler on server ---
# on personal computer, need to increase memory.limit()
# code here for reference: 
library(clusterProfiler)
library(org.Hs.eg.db)

go_results <- lapply(names(de_dx_day), function(dx_day){
  go_enrich <- lapply(c("BP","MF"), function(ontology){
    print(dx_day)
    deg <- de_dx_day[[dx_day]]
    go_enrich_raw <- deg %>%
      dplyr::filter(p < 0.005) %>%
      dplyr::select(ensembl_transcript_id) %>%
      unlist() %>%
      enrichGO(OrgDb = "org.Hs.eg.db",
               universe = deg$ensembl_transcript_id,
               keyType = "ENSEMBLTRANS", #changed from "ENSEMBL"
               ont = ontology,
               readable = T, pvalueCutoff = 1, qvalueCutoff = 1) %>%
      .@result %>%
      mutate(Ontology  = ontology)
  }) %>%
    do.call(rbind,.) %>%
    arrange(pvalue)
})
save(goresultsup_all, file = file.path(parent_folder,"clusterProfiler_results2.rdata"))

go_enrich_raw2 <- de_dx_day[[1]] %>%
  dplyr::filter(p <0.005) %>%
  filter(beta < 0) %>%
  distinct(ensembl_gene_id, .keep_all = TRUE) %>%
  dplyr::select(ensembl_gene_id) %>%
  unlist()
  

xxx <- compareCluster(go_enrich_raw2, fun = "enrichGO",OrgDb = "org.Hs.eg.db",
                      universe = distinct(dat_clusterProfiler[[1]], ensembl_gene_id, .keep_all = TRUE),
                      keyType = "ENSEMBL",
                      readable = T, pvalueCutoff = 0.05)

# Aaron's code for up vs down -------------------------------



go_results_updown44 <- lapply(names(dat_clusterProfiler[44]), function(dx_day){
  print(dx_day)
lapply(c("up_regulated", "down_regulated"),function(direction){
  print(direction)
  if (direction == "up_regulated") {
    go_enrich <- lapply("BP", function(ontology){
      print(ontology)
      deg <- dat_clusterProfiler[[dx_day]]
      go_enrich_raw <- deg %>%
        dplyr::filter(p < 0.005) 
      if(direction == "up_regulated") {
        go_enrich_raw <- go_enrich_raw %>% 
          filter(beta > 0) %>%
          #one instance of each gene id goes into the GO analysis
          distinct(ensembl_gene_id, .keep_all = TRUE)
      } else if (direction == "down_regulated"){
        go_enrich_raw <- go_enrich_raw %>% 
          filter(beta < 0) %>%
          distinct(ensembl_gene_id, .keep_all = TRUE)
      } else {
        print("crash and burn!!")
      }
      go_enrich_raw <- go_enrich_raw %>%
        dplyr::select(ensembl_transcript_id) %>%
        unlist() %>%
        enrichGO(OrgDb = "org.Hs.eg.db",
                 universe = deg$ensembl_transcript_id,
                 keyType = "ENSEMBLTRANS",
                 ont = ontology,
                 readable = T, pvalueCutoff = 1, qvalueCutoff = 1) %>%
        .@result %>%
        mutate(Ontology  = ontology,
               Direction = direction)
    }) %>%
      bind_rows()
  } else if (direction == "down_regulated"){
    go_enrich <- lapply(c("BP","MF"), function(ontology){
      print(ontology)
      deg <- dat_clusterProfiler[[dx_day]]
      go_enrich_raw <- deg %>%
        dplyr::filter(p < 0.005) 
      if(direction == "up_regulated") {
        go_enrich_raw <- go_enrich_raw %>% 
          filter(beta > 0) %>%
          distinct(ensembl_gene_id, .keep_all = TRUE)
      } else if (direction == "down_regulated"){
        go_enrich_raw <- go_enrich_raw %>% 
          filter(beta < 0) %>%
          distinct(ensembl_gene_id, .keep_all = TRUE)
      } else {
        print("crash and burn!!")
      }
      go_enrich_raw <- go_enrich_raw %>%
        dplyr::select(ensembl_transcript_id) %>%
        unlist() %>%
        enrichGO(OrgDb = "org.Hs.eg.db",
                 universe = deg$ensembl_transcript_id,
                 keyType = "ENSEMBLTRANS",
                 ont = ontology,
                 readable = T, pvalueCutoff = 1, qvalueCutoff = 1) %>%
        .@result %>%
        mutate(Ontology  = ontology,
               Direction = direction)
    }) %>%
      bind_rows()
  } 
}) %>%
  bind_rows() %>%
  arrange(pvalue)
}) %>%
  setNames(names(dat_clusterProfiler[44]))



go_results_updown2 <- lapply(names(dat_clusterProfiler[45:54]), function(dx_day){
  print(dx_day)
  lapply(c("up_regulated", "down_regulated"),function(direction){
    go_enrich <- lapply(c("BP","MF"), function(ontology){
      deg <- dat_clusterProfiler[[dx_day]]
      go_enrich_raw <- deg %>%
        dplyr::filter(p < 0.005) 
      if(direction == "up_regulated") {
        go_enrich_raw <- go_enrich_raw %>% 
          filter(beta > 0) %>%
          distinct(ensembl_gene_id, .keep_all = TRUE)
      } else if (direction == "down_regulated"){
        go_enrich_raw <- go_enrich_raw %>% 
          filter(beta < 0) %>%
          distinct(ensembl_gene_id, .keep_all = TRUE)
      } else {
        print("crash and burn!!")
      }
      go_enrich_raw <- go_enrich_raw %>%
        dplyr::select(ensembl_transcript_id) %>%
        unlist() %>%
        enrichGO(OrgDb = "org.Hs.eg.db",
                 universe = deg$ensembl_transcript_id,
                 keyType = "ENSEMBLTRANS",
                 ont = ontology,
                 readable = T, pvalueCutoff = 1, qvalueCutoff = 1) %>%
        .@result %>%
        mutate(Ontology  = ontology,
               Direction = direction)
    }) %>%
      bind_rows()
  }) %>%
    bind_rows() %>%
    arrange(pvalue)
}) %>%
  setNames(names(dat_clusterProfiler[45:54]))

go_results_updown <- lapply(names(dat_clusterProfiler), function(dx_day){
  print(dx_day)
  lapply(c("up_regulated", "down_regulated"),function(direction){
    go_enrich <- lapply(c("BP","MF"), function(ontology){
      deg <- dat_clusterProfiler[[dx_day]]
      go_enrich_raw <- deg %>%
        dplyr::filter(p < 0.005) 
      if(direction == "up_regulated") {
        go_enrich_raw <- go_enrich_raw %>% 
          filter(beta > 0) %>%
          distinct(ensembl_gene_id, .keep_all = TRUE)
        print(dim(go_enrich_raw)[1])
      } else if (direction == "down_regulated"){
        go_enrich_raw <- go_enrich_raw %>% 
          filter(beta < 0) %>%
          distinct(ensembl_gene_id, .keep_all = TRUE)
        print(dim(go_enrich_raw)[1])
      } else {
        print("crash and burn!!")
      }
    })
  })
})


# Determine whether GO cluster is UP or DOWN regulated
names(go_results) <- dx_stage

go_reg <- lapply(dx_stage, function(dx_s){
  go_dx <- go_results[[dx_s]] %>%
    mutate(sep_hgnc = str_split(geneID, "/"))
  go_ensem <- lapply(go_dx$sep_hgnc, function(hgnc_list){
    output <- vector(mode = "list", length = length(hgnc_list))
    for (idx in 1:length(hgnc_list)){
      anno_list <- dplyr::filter(transAnno_id, hgnc_symbol == hgnc_list[idx])
      output[[idx]] <- anno_list$ensembl_transcript_id      
    }
  unlist(output)
  })
  mutate(go_dx, ensemblList = go_ensem)
})

names(go_reg) <- dx_stage



# Analyze cluster profiler results --------------------------------------------------------------------------------
# load(file.path(output_folder,"clusterProfiler_results_up_down.rdata"))
# 
# names(go_results) <- dx_stage

go_plots <- lapply(dxs, function(dx){
  idx <- grep(dx, names(go_results), value = T)
  dx_plots <- lapply(idx, function(dx_day){
    p1 <- goresultsup_all[[dx_day]] %>%  #edited input file
      mutate(log10.p.value = -log10(pvalue), 
             go_size = as.numeric(gsub("/.*","",BgRatio))) %>% 
      filter(pvalue < 0.05 & go_size > 50 & Count > 3) %>% 
      group_by(Direction) %>% 
      arrange(pvalue) %>% 
      dplyr::slice(1:6) %>% 
      arrange(desc(pvalue), .by_group = T) %>% 
      ungroup() %>% 
      mutate(Description = fct_inorder(Description)) %>% 
      ggplot(aes(x = Description , y =  abs(log10.p.value))) +
      geom_col(aes(fill = Direction), position = "dodge") +
      labs(y = "-log10(p value)", x = "") +
      theme(axis.text.y = element_text(hjust = 1, vjust  = 0.3)) +
      scale_x_discrete(labels = function(x) str_wrap(x, width = 40)) +
      scale_fill_manual(values = c("down_regulated" = "#ff7878", "up_regulated" = "#77dd77")) +  #c("grey30", "grey70")
      coord_flip() +
      geom_hline(yintercept = -log10(0.05), color = "blue4") +
      ggtitle(dx_day)
    if(grepl("all_timepoints", dx_day)){
      p1 <- p1 + theme(legend.position = "bottom") 
    } else {
      p1 <- p1 + theme(legend.position = "none") 
    }
    return(p1)
  })
}) %>% setNames(dxs)

pdf(paste0( output_folder,"/", file_index, "_clusterProfilers_results.pdf"), width = 24, height = 14)
lapply(dxs, function(dx){
  idx <- grep(dx, names(goresultsup_all), value = T)
  plot_command <- paste0("go_plots[['",dx,"']][[", 1:length(idx),"]]", collapse = "+") %>% 
    paste0(" + plot_annotation(title = dx) + plot_layout(nrow = 2)")
  eval(parse(text = plot_command))
})
dev.off()  
file_index <- file_index + 1

pdf(paste0( output_folder,"/", file_index, "_clusterProfilers_results_allTPs.pdf"), width = 13.5, height = 7.5)
all_tps_plots <- lapply(go_plots, function(l1){
  l1[[5]] + 
    theme_classic(base_size = 9) + 
    theme(axis.text.y = element_text(hjust = 1, vjust  = 0.3), 
          legend.position = "bottom")
  
})
for (i in c(1:6,8:9)) {
  
  all_tps_plots[[i]] <- all_tps_plots[[i]] + theme(legend.position = "none") 
}

plot_command <- paste0("all_tps_plots[[", 1:length(all_tps_plots),"]]", collapse = "+") %>% 
  paste0(" + plot_layout(nrow = 3)")
eval(parse(text = plot_command))
dev.off()

file_index = file_index + 1

# GSEA ------------------------------------------------------------------------------------------------------------
gsea_dir <- file.path(output_folder, "GSEA")
dir.create(gsea_dir)
####!!!!run gsea on hoffman2
load(file.path(gsea_dir,"gsea_results.rdata"))

plot_gsea_data <- map(fgsea_res, function(res){
  res %>%
    mutate(direction = ifelse(NES > 0 ,"up_regulated","down_regulated"),
           pathway = gsub("GO_", "", pathway)) %>%
    mutate(pathway = str_to_sentence(gsub("_", " ", pathway))) %>%
    filter(padj < 0.05) %>%
    arrange(desc(abs(NES))) %>%
    group_by(direction) %>%
    dplyr::slice(1:7) %>%
    arrange(desc(NES), .by_group = T) %>%
    ungroup() %>%
    mutate(pathway = fct_inorder(pathway))
})

plot_number_go_terms <- function(expr, nm){
  related_go  <- map_df(fgsea_res, function(res){
    res %>%
      mutate(direction = ifelse(NES >0 ,"up_regulated","down_regulated")) %>%
      filter(padj < 0.05) %>%
      filter(grepl(toupper(expr), pathway)) %>% 
      group_by(direction) %>% 
      summarise(n = n())
  },.id = "condition") %>% 
    mutate(condition = gsub("_ASD", "", condition)) %>% 
    separate(condition, into = c("dx", "stage"), sep = "_", remove = FALSE, extra = "merge")  %>% 
    filter(grepl("\\d+", stage)) 
  ggplot(related_go, aes(x = stage, y = n, fill = direction)) +
    geom_col(position = "dodge")+
    facet_wrap(~dx) +
    labs(y = paste("# GO-terms"), title = paste(nm, "related GO-terms"))
}

gsea_p1 <- plot_number_go_terms("synap", "Syanpse")
gsea_p2 <- plot_number_go_terms("mitot", "Mitotic")
pdf(file.path(gsea_dir,"gsea_results_specific.pdf"), w = 14)
wrap_plots(gsea_p1, gsea_p2, guides = "collect")
dev.off()
# tmp <- fgsea_res[["16p11del_025"]] %>% 
#   filter(padj < 0.05) %>% 
#   filter(grepl(toupper("differe"), pathway)) 

gsea_plots <- map(set_names(names(plot_gsea_data)) , function(dx){
  plot_gsea_data[[dx]] %>%
    ggplot(aes(x = pathway , y =  abs(NES))) +
    geom_col(aes(fill = direction), position = "dodge") +
    labs(y = "NES (absolute value)", x = "") +
    theme(axis.text.y = element_text(hjust = 1, vjust  = 0.3)) +
    scale_x_discrete(labels = function(x) str_wrap(x, width = 40)) +
    scale_fill_manual(values = c("grey30", "grey70")) +
    coord_flip() +
    #geom_hline(yintercept = 0, color = "indianred") +
    ggtitle(dx) +
    theme_classic() +
    theme(legend.direction = "horizontal")
})

pdf(file.path(gsea_dir,"gsea_results.pdf"), width = 10, height = 14)
for(dx in dxs){
  dx <- gsub("ASD","", dx)
  print(
    wrap_plots(gsea_plots[grep(dx, names(gsea_plots))], 
               ncol = 2) + 
      plot_layout(guides = 'collect') &
      theme(legend.position = "bottom")
  )
}
dev.off()

gsea_subset_plots <- grep("025|100", names(plot_gsea_data), value = T) %>% 
  set_names() %>% 
  
  map(function(dx){
    plot_gsea_data[[dx]] %>%
      arrange(desc(abs(NES))) %>%
      group_by(direction) %>%
      dplyr::slice(1:2) %>%
      arrange(desc(NES), .by_group = T) %>%    
      ungroup() %>%
      mutate(pathway = fct_inorder(as.character(pathway))) %>%
      ggplot(aes(x = pathway , y =  abs(NES))) +
      geom_col(aes(fill = direction), position = "dodge") +
      labs(y = "absolute NES", x = "") +
      theme(axis.text.y = element_text(hjust = 1, vjust  = 0.3)) +
      scale_x_discrete(labels = function(x) str_wrap(x, width = 40)) +
      #scale_fill_manual(values = c("grey30", "grey70")) +
      coord_flip() +
      #geom_hline(yintercept = 0, color = "indianred") +
      ggtitle(sub("(.*)(_)(.*)?", "\\1\n\\3", dx)) +
      theme_classic() +
      theme(legend.direction = "horizontal", 
            axis.text.y = element_text(size = 12))
  })

pdf(file.path(gsea_dir,"gsea_results_subset.pdf"), width = 20, height = 12)
wrap_plots(gsea_subset_plots, guides = "collect", ncol = 4) &
  theme(legend.position = "bottom")
dev.off()



load(file.path(gsea_dir, "genedat.rdata"))

geneDat_full <- geneDat %>%
  left_join(geneAnno,by = "ensembl_gene_id") %>%
  dplyr::select(ensembl_gene_id, hgnc_symbol, entrezgene_id)

res_for_file <- map(names(fgsea_res),function(f){
  res <- fgsea_res[[f]] %>%
    mutate(pathway = gsub("GO_", "", pathway)) %>%
    mutate(pathway = gsub("_", " ", pathway)) %>%
    filter(padj < 0.05) %>%
    arrange(desc(NES))
  
  
  hgnc_le <- map(res$leadingEdge,function(le){
    data.frame(entrezgene_id = as.numeric(le)) %>%
      left_join(geneDat_full,by = "entrezgene_id") %>%
      dplyr::select(hgnc_symbol) %>%
      unlist() %>%
      paste(collapse = "|")
  })
  res$leadingEdge <- unlist(hgnc_le)
  write_csv(res,
            path = file.path(gsea_dir,paste0("gsea_res_",f,".csv")))
  return(res)
}) %>% setNames(names(fgsea_res))
writexl::write_xlsx(res_for_file, path = file.path(gsea_dir,"gsea_results.xlsx"))



# Correlation of logFC within Dx  --------------------------------------------------------------------------------
dx_logFC_cor <- lapply(dxs, function(dx){
  idx <- grep(paste0(dx,"_\\d+"), names(de_dx_day), value = T)
  
  betas <- lapply(names(de_dx_day[idx]), function(dx_day){
    deg <- de_dx_day[[dx_day]][]
    deg <- deg[order(deg$ensembl_transcript_id),c("beta")] %>% 
      setNames(dx_day)
  })  %>% do.call(cbind,.)
  cormat <- cor(betas, use = "pairwise", method = "spearman")
  cormat[lower.tri(cormat)] <- NA
  diag(cormat) <- NA
  melted_cormat <- cormat %>% 
    as.data.frame() %>% 
    rownames_to_column("var1") %>% 
    gather("var2","cor", -var1) %>% 
    drop_na() %>% 
    mutate_at(vars(var1, var2), funs(gsub(".*_","", .))) %>% 
    mutate(dx = dx)
}) %>% do.call(rbind,.) %>% 
  mutate_at(vars(var1, var2), funs(factor(., levels = union(var1, var2))))

pdf(paste0(output_folder,"/", file_index, "_logFC_cor_within_Dx.pdf"), height = 6, width = 12)
ggplot(dx_logFC_cor, aes(var2, var1, fill = cor)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, limit = c(-1,1), space = "Lab",
                       name = "Spearman's Rho") +
  theme_bw() + 
  coord_fixed() +
  geom_text(aes(var2, var1, label = signif(cor,2)), color = "black", size = 4) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1),
        axis.text.y = element_text(size = 12),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        strip.text = element_text(size = 12, face = "bold"), 
        strip.background = element_rect(fill = "white")) +
  scale_x_discrete(drop = F) +
  scale_y_discrete(drop = F) +
  labs(x = "Differentiation day", y = "Differentiation day") +
  facet_wrap(~dx, nrow = 2)

dev.off()
file_index <- file_index + 1


pdf(paste0(output_folder,"/", file_index, "_logFC_cor_within_Dx_network.pdf"), height = 7, width = 12)
dx_logFC_cor %>% 
  graph_from_data_frame(directed = FALSE) %>% 
  ggraph(layout = 'igraph', algorithm = 'circle')+
  geom_edge_link(aes(edge_alpha = abs(cor), edge_width = abs(cor), color = cor)) +
  guides(edge_alpha = "none", edge_width = "none") +
  scale_edge_colour_gradientn(limits = c(-1, 1), colors = c("dodgerblue2", "white","indianred2")) +
  geom_node_point(color = "grey95", size = 5) +
  geom_node_text(aes(label = name), repel = TRUE) +
  theme_graph(base_family = "") +
  facet_wrap(~dx)
dev.off()
file_index <- file_index + 1


# Correlation of logFC between Dx ---------------------------------------------------------------------------------
stage_logFC_cor <- lapply(stages, function(stg) {
  idx <- grep(stg, names(de_dx_day), value = T)
  betas <- lapply(names(de_dx_day[idx]), function(dx_day){
    deg <- de_dx_day[[dx_day]][]
    deg <- deg[order(deg$ensembl_transcript_id),c("beta")] %>% 
      setNames(gsub("_.*","", dx_day))
  })  %>% do.call(cbind,.) 
  cormat <- cor(betas, use = "pairwise", method = "spearman")
}) %>% setNames(stages)

logfc_heatmaps <- lapply(stages, function(stg) {
  
  curent_plot <- pheatmap(stage_logFC_cor[[stg]], 
                          color = colorRampPalette( c("green","black", "red"))(200),
                          breaks = seq(-1,1,0.01),
                          border_color = NA,  
                          cellwidth = 20, 
                          cellheight = 20, main = gsub("_"," ", stg), silent = T)
  curent_plot[[4]]
})
png(paste0(output_folder,"/", file_index, "_logFC_cor_between_Dx.png"), height = 8, width = 16, units = "in", res = 300)
do.call(grid.arrange, c(logfc_heatmaps, nrow = 2))
dev.off()
file_index <- file_index +1
#plot networks
stage_logFC_cor_tidy <- lapply(stages, function(stg){
  cormat <- stage_logFC_cor[[stg]]
  cormat[lower.tri(cormat)] <- NA
  diag(cormat) <- NA
  melted_cormat <- cormat %>% 
    as.data.frame() %>% 
    rownames_to_column("var1") %>% 
    gather("var2","cor", -var1) %>% 
    drop_na() %>% 
    mutate_at(vars(var1, var2), funs(gsub(".*_","", .))) %>% 
    mutate(stage = stg)
}) %>% do.call(rbind,.)

pdf(paste0(output_folder,"/", file_index, "_logFC_cor_between_Dx_network.pdf"), height = 7, width = 12)
stage_logFC_cor_tidy %>% 
  graph_from_data_frame(directed = FALSE) %>% 
  ggraph(layout = 'igraph', algorithm = 'circle') +
  geom_edge_link(aes(edge_alpha = abs(cor), edge_width = abs(cor), color = cor)) +
  guides(edge_alpha = "none", edge_width = "none") +
  scale_edge_colour_gradientn(limits = c(-1, 1), colors = c("blue", "white","red")) +
  geom_node_point(color = "grey95", size = 5) +
  geom_node_text(aes(label = name), repel = TRUE) +
  theme_graph(base_family = "") +
  facet_wrap(~stage)
dev.off()
file_index <- file_index + 1

# plot all timepoints/Dxs -----------------------------------------

betas <- lapply(names(de_dx_day), function(dx_day){
  deg <- de_dx_day[[dx_day]][]
  deg <- deg[order(deg$ensembl_transcript_id),c("beta")] %>% 
    setNames(dx_day)
})  %>% 
  do.call(cbind,.) %>% 
  .[grepl("\\d{3}",names(.))]
cormat <- cor(betas, use = "pairwise", method = "spearman")
save(cormat, file =  file.path(output_folder,"correlation matrix.rdata"))
# cormat2 <- psych::corr.test(betas, method = "spearman", adjust = "none")
# sum(cormat2$r > 0.3)
# 
# sum(cormat2$p > 0.05)
# dim(cormat2$p)
cormat[lower.tri(cormat)] <- NA
diag(cormat) <- NA
melted_cormat <- cormat %>% 
  as.data.frame() %>% 
  rownames_to_column("var1") %>% 
  gather("var2","cor", -var1) %>% 
  drop_na() 
pdf(paste0(output_folder,"/", file_index, "_logFC_cor_all_points_network.pdf"), height = 12, width = 12)
melted_cormat %>% 
  mutate(Dx = gsub("(.*?)_.*","\\1", var1)) %>% 
  filter(abs(cor) > 0.2) %>% 
  graph_from_data_frame(directed = FALSE) %>% 
  ggraph(layout = 'igraph', algorithm = "kk") +
  geom_edge_link(aes(edge_alpha = abs(cor), edge_width = abs(cor), color = cor)) +
  guides(edge_alpha = "none", edge_width = "none") +
  scale_edge_colour_gradientn(limits = c(-1, 1), colors = c("blue", "white","red")) +
  geom_node_point(aes(color = gsub("(.*?)_.*","\\1", name)), size = 5) +
  scale_color_discrete(name = "Dx") +
  geom_node_text(aes(label = name), repel = TRUE, size = 6) +
  labs(caption = "subset to |cor| > 0.2") + 
  theme_void(base_size = 18)
dev.off()



cormat1 <- cor(betas, use = "pairwise", method = "spearman")

row_anno <- data.frame(Dx = gsub("(.*?)_.*","\\1",rownames(cormat1)), 
                       Day = gsub(".*?_(.*)","\\1",rownames(cormat1)))
rownames(row_anno) <- colnames(cormat1)

col_anno_colors = list(
  # Module = setNames(col_anno$Module, nm = col_anno$Module), 
  # Direction = c("down regulated" = "seagreen3", "up regulated" = "indianred1"),
  Day = setNames(colorRampPalette(c("cadetblue1", "cadetblue4"))(4), nm = c("025", "050", "075", "100"))
  # c("22q11.2" = "white", "15q13.3" = "dodgerblue", "16p11.2" = "dodgerblue4")
)
pdf(paste0(output_folder,"/", file_index, "_logFC_cor_all_points_heatmap.pdf"), height = 8, width = 10)
pheatmap(cormat1,
         clustering_method = "ward.D2",
         color = colorRampPalette( c("darkblue","white", "firebrick3"))(100),
         annotation_row = row_anno,
         annotation_col = row_anno, 
         annotation_colors = col_anno_colors,
         # annotation_names_row = F, 
         breaks = seq(-1,1,0.02),
         cellwidth = 10, 
         cellheight = 10, 
         gaps_row = 4,
         border_color = NA,
         cutree_rows = 6,
         cutree_cols = 6)
dev.off()







file_index = file_index + 1

# overlap of DEGs -------------------------------------------------------------------------------------------------
all_degs <- lapply(de_dx_day, function(x){
  x %>% 
    dplyr::filter(p < 0.005) %>% 
    dplyr::select(ensembl_transcript_id) %>% 
    unlist()
})


pdf(paste0(output_folder,"/", file_index, "_upset_byDx.pdf"),  height = 8, width = 16)#, units = "in", res = 300)
for (dx in dxs) {
  idx <- grep(paste0(dx,"_\\d+"), names(de_dx_day), value = T)
  current_degs <- all_degs[idx] %>% 
    lapply(unname) %>% 
    setNames(gsub(".*_","", names(.)))
  
  print(upset(fromList(current_degs), 
              sets = names(current_degs), 
              keep.order = T, sets.bar.color = "#56B4E9",nintersects=60,
              order.by = "freq", empty.intersections = "off",
              mb.ratio = c(0.65,0.35), text.scale = 1.8))
  grid.text(dx,x = 0.65, y=0.95, gp=gpar(fontsize=16))
}
dev.off()


pdf(paste0(output_folder,"/", file_index, "_upset_byDay.pdf"),  height = 8, width = 16)
for (stg in stages) {
  idx <- grep(stg, names(de_dx_day), value = T)
  current_degs <- all_degs[idx] %>% 
    lapply(unname) %>% 
    setNames(gsub("(.*?)_.*","\\1", names(.)))
  
  print(upset(fromList(current_degs), 
              sets = names(current_degs), 
              keep.order = T, sets.bar.color = "#56B4E9",nintersects=60,
              order.by = "freq", empty.intersections = "off",
              mb.ratio = c(0.65,0.35), text.scale = 1.8))
  grid.text(stg,x = 0.65, y = 0.95, gp = gpar(fontsize=16))
}
dev.off()

#overalp of DEG inside the CNV regions

cnv_overlap <- map_df(set_names(names(de_dx_day)), function(dx_day){
  map_df(set_names(names(loci)), function(lci){
    mut <- loci[[lci]]
    location <- gsub("del|dup", "", lci)
    if(!grepl(location, dx_day)) {
      if(names(mut)[3]=="cnv"){
        de_dx_day[[dx_day]] %>% 
          filter(fdr < 0.05) %>% 
          filter(chromosome_name == mut[["chr"]] &
                   end_position > mut[["cnv"]][1]  &
                   start_position < mut[["cnv"]][2])
      } else if (names(mut)[3] == "gene") {
        de_dx_day[[dx_day]] %>% 
          filter(fdr < 0.05) %>% 
          filter(hgnc_symbol == mut[["gene"]])
      }
    }
  }, .id = "cnv")
}, .id = "dx_day") %>% 
  filter(cnv != "16p11dup")

write.csv(cnv_overlap,file.path(output_folder, "overlap of cnvs and degs.csv"))


# correaltion of expression ---------------------------------------------------------------------------------------
expr_cor <- cor(datExpr_reg_batch, method = "s")

row_anno1 <- datMeta %>% 
  select(SampleID,Dx, Day = Day.grouped ) %>% 
  as.data.frame() %>% 
  `rownames<-`(.$SampleID) %>% 
  .[,-1] %>% 
  .[colnames(expr_cor),]

row_anno2 <- map_df(set_names(unique(datMeta$Dx)), function(dx){
  as.numeric(row_anno1$Dx == dx)
})
row_anno <- cbind(row_anno, row_anno2)


col_anno_colors = list(
  Day = setNames(colorRampPalette(c("grey90", "cadetblue4"))(4), nm = c("025", "050", "075", "100"))
) %>% 
  c(.,map(set_names(colnames(row_anno2)),~setNames(c("white", "black"), nm = c("1", "0"))))
pdf(paste0(output_folder,"/", file_index, "_exprcor_all_heatmap.pdf"), height = 10, width = 15)
pheatmap(expr_cor,
         color = colorRampPalette( c("darkblue","white", "firebrick3"))(100),
         #annotation_row = row_anno,
         annotation_col = row_anno, 
         annotation_colors = col_anno_colors,
         # annotation_names_row = F, 
         #breaks = seq(-1,1,0.02),
         cellwidth = 1, 
         cellheight = 1, 
         gaps_row = 4,
         border_color = NA,
         show_rownames = F,
         show_colnames = F
         # cutree_rows = 6,
         
         #cutree_cols = 6
)
dev.off()


all_unique_degs <- unique(unlist(all_degs)) 

top_var_deg <- apply(datExpr_reg_batch[all_unique_degs,],1, sd) %>% 
  sort(decreasing = T) %>% 
  head(500) %>% 
  names(.)


pdf(paste0(output_folder,"/", file_index, "_expr_all_heatmap.pdf"), height = 20, width = 15)
p1 <- pheatmap(datExpr_reg_batch[top_var_deg,],
               color = colorRampPalette( c("darkblue","white", "firebrick3"))(100),
               #annotation_row = row_anno,
               annotation_col = row_anno, 
               annotation_colors = col_anno_colors,
               # annotation_names_row = F, 
               #breaks = seq(-1,1,0.02),
               cellwidth = 1, 
               cellheight = 1, 
               gaps_row = 4,
               border_color = NA,
               show_rownames = F,
               show_colnames = F,
               # cutree_rows = 6,
               scale = "row"
               #cutree_cols = 6, 
               
)
dev.off()



#EWCE-------
##!!!!! run EWCE on hoffman
ewce_dir <- file.path(output_folder,"EWCE")
datasets <- c("hcs_pasca","PolioudakisNeuralFetal") 
ds <- datasets[1]
for (ds in datasets){
  for(lvl in c(1,2)){  
    fs <- list.files(file.path(ewce_dir,ds), pattern = paste0("_",lvl,".Rdata$"), full.names = T)
    ewce_res <- lapply(fs, function(f){
      load(f)
      ewce <- ls(pattern = "EWCE") %>% 
        get() %>% 
        map(~bind_rows(.x, .id = "direction")) %>% 
        map(~mutate(.x, p = ifelse(p == 0, 1e-6, p))) %>% 
        map(~mutate(.x, fdr = p.adjust(p)))
      return(ewce)
    }) %>% 
      unlist(recursive = F) %>% 
      bind_rows(.id = "dx_day") %>% 
      separate(dx_day, into = c("Dx", "stage"), remove = FALSE, sep = "_", extra = "merge") %>% 
      mutate(fold_change = if_else(fdr < 0.05, fold_change, 0),
             txt_fdr = symnum(fdr, cutpoints = c(0, 0.005,0.01,0.05,1), symbols = c("***","**","*",""))) %>%
      mutate(txt_fc = ifelse(fdr < 0.05,signif(fold_change,2) , ""),
             txt = ifelse(txt_fdr == "", txt_fdr, paste(signif(fold_change,1), txt_fdr, sep = "\n")))
    
    if (ds == "hcs_pasca") {  
      if (lvl == 1) {
        ewce_res <- ewce_res %>% 
          mutate(CellType = factor(CellType, 
                                   levels = c("Radial_Glia", "IMPs", 
                                              "Excitatory_Neurons" , "Interneurons", 
                                              "Astrocytes_ChoroidEndo", "OPCs")))
      } else if (lvl == 2) {
        ewce_res <- ewce_res %>% 
          mutate(CellType = factor(CellType, 
                                   levels = c("MphaseRadial","SphaseRadialGlia","MatureRadialGlia", 
                                              "CyclingIMPs" ,"VGLUT1_IMPs",
                                              "ExtNeuronsL566-like","ExtNeuronsL56-like","ExtNeuronsL556-Like", "ExtNeuronsL456like_Mature", 
                                              "Interneurons(all)" , 
                                              "Astrocytes_ChoroidEndo", "OPCs")))
      }
    } else if (ds == "PolioudakisNeuralFetal") {  
      if (lvl == 1) {
        ewce_res <- ewce_res %>% 
          mutate(CellType = factor(CellType, 
                                   levels = c("CycProg","RG",  
                                              "IP",
                                              "ExcNeu","IntNeu" ,
                                              "Microglia","OPC","Pericyte","Endothelial" )))
      } else if (lvl == 2) {
        ewce_res <- ewce_res %>% 
          mutate(CellType = factor(CellType, 
                                   levels = c("CycProg_G2M","CycProg_S", "oRG","vRG",
                                              "IP",
                                              "ExcNeu_Mat","ExcNeu_DL1","ExcNeu_Mig","ExcNeu_Cal","ExcNeu_DL2", 
                                              "IntNeu_SST","IntNeu_CALB2" ,
                                              "Pericyte", "Microglia","OPC","Endothelial"
                                   )))
      }
    }
    if (lvl == 1) {
      wi = 12
    } else if (lvl == 2) {
      wi = 18
    }
    ewce_res <- ewce_res %>% 
      filter(grepl("\\d",stage)) %>% 
      mutate(fold_change = ifelse(direction == "down_regulated", -fold_change, fold_change), 
             Dx = gsub("ASD","",Dx))
    pdf(file.path(ewce_dir, paste0("ewce_", ds,"_",lvl,".pdf")), w = wi, h = 8)    
    
    print(ggplot(ewce_res, aes(x = stage, fill = direction, y = fold_change)) +
            geom_col() +
            facet_grid(Dx~CellType) +
            theme_bw(base_size = 14) +
            theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5), 
                  strip.background = element_rect(fill = "white"), 
                  panel.grid.minor = element_blank(), 
                  legend.position = "left"))
    
    # print(ggplot(ewce_res, aes(x = CellType, y = direction, fill = fold_change)) +
    #        geom_tile() +
    #        geom_text(aes(label = txt_fc), size = 2) +
    #        facet_grid(stage~Dx, scales = "free") +
    #        theme_bw(base_size = 14) +
    #        theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5)) +
    #        scale_fill_gradient2(low = "white",mid = "white", high = "dodgerblue2", midpoint = 1) +
    #        theme(strip.background = element_rect(fill = "white")))
    dev.off()
  }
}


