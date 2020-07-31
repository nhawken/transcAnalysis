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
output_folder = paste(getwd(), '/compare_DEIvDEG_OUT_',Sys.Date(),'/', sep = "") # SET AN OUTPUT PATH HERE
dir.create(output_folder, showWarnings = F, recursive = T)
#time0 <- proc.time()

theme_set(theme_classic(base_size = 18))

load("all_transcript_data.Rdata")
load(file.path(parent_folder,"input_for_DE_transc.rdata"))

cond <- "dx_day"

# source("~/Aaron/Scripts/Aaron/CIRM/3_Diagnosis/DE_functions.R")

# read in DE results
dei_results <- read_csv(paste0(parent_folder, "/",cond, "_lme_DE.csv"))[,-2] %>% #calls all but the 2nd col of the matrix
  rename("ensembl_transcript_id" = "X1") %>%  #renames unlabeled column
  setNames(gsub("Idiopathic(.?ASD)?","IdiopathicASD", colnames(.))) #rename idiopathic ASD columns 

deg_results <- read_csv(paste0(parent_folder, "/",cond, "_lme_genes.csv"))[,-2] %>% #calls all but the 2nd col of the matrix
  rename("ensembl_gene_id" = "X1") %>%  #renames unlabeled column
  setNames(gsub("Idiopathic(.?ASD)?","IdiopathicASD", colnames(.))) #rename idiopathic ASD columns 


dxs <- gsub(".*\\.(.*?)_.*","\\1", colnames(dei_results)) %>% #gsub one or more of any charater followed by a period, (.*?) is the capture group
  unique() %>%  #gets rid of Beta/p/SE. in front, keeps the next dx, gets rid of the _025/50/late afterwards
  sort() %>%
  .[. != "ensembl_transcript_id"] 

stages <- gsub(".*\\..*?_(.*)","\\1", colnames(dei_results)) %>% #get rid of Beta/p/Se. in front then gets rid of Dx_, capture anything after Dx
  unique() %>% 
  sort() %>%
  .[. != "ensembl_transcript_id"] 

dx_stage <- gsub(".*\\.(.*?_.*)","\\1", colnames(dei_results)) %>%  #get rid of Beta/p/SE., then capture dx and stages
  unique() %>% 
  sort() %>%
  .[. != "ensembl_transcript_id"] 


# # correct p value and add annotation
dei_results_adj <- lapply(dx_stage, function(dx_stg){ #for each dx_stage, create a list with 1 data frame (tibble) containing beta, SE, p, fdr
  p_val_col <- paste0("p.",dx_stg)
  dei_results %>% 
    dplyr::select(contains(dx_stg)) %>% 
    mutate(fdr =   !!as.name(p_val_col) %>% p.adjust(., method = "fdr")) %>% #control the false discovery rate, the expected proportion of false discoveries amongst the rejected hypotheses
    # mutate(data fram, new_var = [existing var])
    rename_at(vars(matches("fdr")), list(~paste0("fdr.",dx_stg)))
}) %>% bind_cols() %>%  #combines dataframes, list of dataframs
  mutate(ensembl_transcript_id = dei_results$ensembl_transcript_id) %>% #add new column to df with the ensembl tanscript ids
  dplyr::select(ensembl_transcript_id, everything()) #puts all other columns after the ID column

deg_results_adj <- lapply(dx_stage, function(dx_stg){ #for each dx_stage, create a list with 1 data frame (tibble) containing beta, SE, p, fdr
  p_val_col <- paste0("p.",dx_stg)
  deg_results %>% 
    dplyr::select(contains(dx_stg)) %>% 
    mutate(fdr =   !!as.name(p_val_col) %>% p.adjust(., method = "fdr")) %>% #control the false discovery rate, the expected proportion of false discoveries amongst the rejected hypotheses
    # mutate(data fram, new_var = [existing var])
    rename_at(vars(matches("fdr")), list(~paste0("fdr.",dx_stg)))
}) %>% bind_cols() %>%  #combines dataframes, list of dataframs
  mutate(ensembl_gene_id = deg_results$ensembl_gene_id) %>% #add new column to df with the ensembl tanscript ids
  dplyr::select(ensembl_gene_id, everything()) #puts all other columns after the ID column


transAnno <-  transcriptAnnoRaw %>%
  filter(ensembl_transcript_id %in% dei_results$ensembl_transcript_id) %>%
  filter(!duplicated(ensembl_transcript_id))

library(biomaRt)
getinfo <- c("ensembl_gene_id","hgnc_symbol","chromosome_name","band","strand","start_position", "end_position","gene_biotype","transcript_length","percentage_gene_gc_content")
humanMart <- useMart(biomart="ensembl",dataset="hsapiens_gene_ensembl")
geneAnnoRaw <- getBM(attributes = getinfo,filters=c("ensembl_gene_id"),values=deg_results$ensembl_gene_id,mart=humanMart)

geneAnno <-  geneAnnoRaw %>%
  filter(ensembl_gene_id %in% deg_results$ensembl_gene_id) %>%
  filter(gene_biotype != "rRNA") %>% 
  filter(!duplicated(ensembl_gene_id))

datExpr <- datExpr[match(geneAnno$ensembl_gene_id,rownames(datExpr)),]


# separate dx_day into seperate dfs
dei_dx_day <- lapply(dx_stage, function(dx_stg){
  dei_results_adj %>% 
    dplyr::select(ensembl_transcript_id,contains(dx_stg))  %>% #keep only the variables ID, and dx_stage of interest (beta, se, p, fdr)
    setNames(gsub(paste0("\\.",dx_stg), "", colnames(.))) %>% #remove .dx_stage identifier from the columns of beta, se, p, fdr
    left_join(transAnno, by = "ensembl_transcript_id") #add transcript annotation table to the right edge of the table
})  %>% setNames(dx_stage) #rename overall df to be the dx_stage identifier

deg_dx_day <- lapply(dx_stage, function(dx_stg){
  deg_results_adj %>% 
    dplyr::select(ensembl_gene_id,contains(dx_stg))  %>% #keep only the variables ID, and dx_stage of interest (beta, se, p, fdr)
    setNames(gsub(paste0("\\.",dx_stg), "", colnames(.))) %>% #remove .dx_stage identifier from the columns of beta, se, p, fdr
    left_join(geneAnno, by = "ensembl_gene_id") #add transcript annotation table to the right edge of the table
})  %>% setNames(dx_stage) #rename overall df to be the dx_stage identifier


# analyze for significant genes/isoforms --------------------------------------

dei_n <- lapply(dx_stage, function(dx_stg){
  res <- dei_dx_day[[dx_stg]] #call df with specific dx_stage id
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
    gather( key = "variable", value = "value", upRegulated:downRegulated) %>%
    mutate(dataset = "transcript") 
}) %>% 
  do.call(rbind,.)

dei_sig <- lapply(dx_stage, function(dx_stg){
  res <- dei_dx_day[[dx_stg]] #call df with specific dx_stage id
  stage = gsub(".*?_(.*)","\\1", dx_stg) #pull out stage from id
  dx = gsub("(.*?)_.*","\\1", dx_stg) #pull out dx from id
  
  sigRes1 <- res %>%
    dplyr::mutate(Direction = if_else(beta > 0, "upRegulated", "downRegulated")) %>%
    filter(p < 0.005) #grab only p vals less than 0.005
}) %>% 
  setNames(names(dei_dx_day))


deg_n <- lapply(dx_stage, function(dx_stg){
  res <- deg_dx_day[[dx_stg]] #call df with specific dx_stage id
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
    gather( key = "variable", value = "value", upRegulated:downRegulated) %>%
    mutate(dataset = "gene")
}) %>% 
  do.call(rbind,.)

deg_sig <- lapply(dx_stage, function(dx_stg){
  res <- deg_dx_day[[dx_stg]] #call df with specific dx_stage id
  stage = gsub(".*?_(.*)","\\1", dx_stg) #pull out stage from id
  dx = gsub("(.*?)_.*","\\1", dx_stg) #pull out dx from id
  
  sigRes1 <- res %>%
    dplyr::mutate(Direction = if_else(beta > 0, "upRegulated", "downRegulated")) %>%
    filter(p < 0.005) #grab only p vals less than 0.005
}) %>% 
  setNames(names(deg_dx_day))

de_all_n <- as_tibble(rbind(dei_n, deg_n)) %>%
  dplyr::group_by(pvalue) %>%
  group_by(Dx)

# plotting counts of differentially expressed genes and transcripts as bar charts
dei_deg_bar_counts <- lapply(unique(de_all_n$pvalue), function(pset){
  de_all_dat <- dplyr::filter(de_all_n, pvalue == pset) %>%
    ggplot( aes( fill = dataset, y = value, x = stage)) +
    geom_bar(position = "dodge2", stat = "identity")+
    facet_wrap(~Dx, ncol = 3, scales = "free_y") +
    ggtitle(paste0("DEI and DEG at ", pset)) +
    theme(axis.text.x = element_text(
      angle = -90,
      hjust = 0,
      vjust = 0.5
    )) +
    theme(legend.text = element_text(size = 10))
}) %>% setNames(unique(de_all_n$pvalue))

pdf(file = paste0(output_folder, "/", file_index, "_DEI_DEG_bar_counts.pdf"), height = 10, width = 14)
print(dei_deg_bar_counts)
dev.off()
file_index <- file_index + 1

dei_genelist <- lapply(dx_stage, function(dx_stg){
  dei_list <- dei_sig[[dx_stg]] %>%
    dplyr::select(ensembl_gene_id, Direction)
}) %>% setNames(dx_stage)

deg_genelist <- lapply(dx_stage, function(dx_stg){
  deg_list <- deg_sig[[dx_stg]] %>%
    dplyr::select(ensembl_gene_id, Direction)
}) %>% setNames(dx_stage)

dei_not_deg_genelist <- lapply(dx_stage, function(dx_stg){
  dei_list <- dei_genelist[[dx_stg]]$ensembl_gene_id
  deg_list <- deg_genelist[[dx_stg]]$ensembl_gene_id
  dei_genelist[[dx_stg]][!(dei_list %in% deg_list),]
}) %>% setNames(dx_stage)

deg_not_dei_genelist <- lapply(dx_stage, function(dx_stg){
  dei_list <- dei_genelist[[dx_stg]]$ensembl_gene_id
  deg_list <- deg_genelist[[dx_stg]]$ensembl_gene_id
  deg_genelist[[dx_stg]][!(deg_list %in% dei_list),]
}) %>% setNames(dx_stage)


# Ensembl Gene id venn diagrams
library(VennDiagram)
library(gridExtra)
dei_deg_venn <- lapply(c("downRegulated", "upRegulated"), function(direc){
  lapply(dx_stage, function(dx_stg){
    dei_list <- dplyr::filter(dei_genelist[[dx_stg]], Direction == direc)
    dei_genes <- dei_list$ensembl_gene_id
    deg_list <- dplyr::filter(deg_genelist[[dx_stg]], Direction == direc)
    deg_genes <- deg_list$ensembl_gene_id
    venn.diagram(x = list(dei_genes, deg_genes),
                 category.names = c("DEI", "DEG"),
                 main = paste0(direc, " Genes at ", dx_stg),
                 filename = NULL)
    }) %>% setNames(dx_stage)
}) %>% setNames(c("downRegulated", "upRegulated")) %>%
  unlist(recursive = FALSE)


pdf(file = paste0(output_folder, "/", file_index, "_DEI_DEG_venn_counts.pdf"), height = 6, width = 16)
for (dx_id in dxs){
    plot_id_down <- paste0("downRegulated", ".", dx_id)
    plots_down <- dei_deg_venn[grepl(plot_id_down, names(dei_deg_venn))]
    plot_id_up <- paste0("upRegulated", ".", dx_id)
    plots_up <- dei_deg_venn[grepl(plot_id_up, names(dei_deg_venn))]
    #grid.arrange(grobs = c(plots_down, plots_up))
    grid.arrange(grobs=lapply(c(plots_down, plots_up), grobTree), ncol=6)
}
dev.off()
file_index <- file_index + 1



library(gt)

# 
# de_all_wide %>%
#   gt(rowname_col = "Dx") %>%
#   tab_header(title = "Differentially Expressed Counts") %>%
#   fmt_number(
#     columns = vars(value),
#     decimals = 0,
#     use_seps = TRUE) %<%
#   tab_spanner(
#     label = "Direction",
#     
#   )
#   
 
# comparing GO results

load("clusterProfiler_results2.rdata")
dei_go <- goresultsup_all
load("gene_clusterProfiler_results_up_down.rdata")
deg_go <- go_results

sig_go_terms<- lapply(dx_stage, function(dx_stg){
  deg_go_set <- deg_go[[dx_stg]] %>%
    mutate(go_size = as.numeric(gsub("/.*","",BgRatio))) %>%
    filter(pvalue < 0.05 & go_size > 50 & Count > 3) %>%
    mutate(dataset = "gene")
  dei_go_set <- dei_go[[dx_stg]] %>%
    mutate(go_size = as.numeric(gsub("/.*","",BgRatio))) %>%
    filter(pvalue < 0.05 & go_size > 50 & Count > 3) %>%
    mutate(dataset = "transcript")
  all_go_set <- rbind(deg_go_set,dei_go_set) %>% 
      as.data.frame()
}) %>% setNames(dx_stage)


go_sig_n <- data.frame(Stage = character(), Dx = character(), direction = character(), dataset = character(), value = numeric())
go_sig_idx = 1
for (dx_stg in dx_stage){
  for (dat_set in c("gene", "transcript")){
    for (direc in c("down_regulated", "up_regulated")){
      stage = gsub(".*?_(.*)","\\1", dx_stg) #pull out stage from id
      dx = gsub("(.*?)_.*","\\1", dx_stg) #pull out dx from id  
      sig_data <- sig_go_terms[[dx_stg]] %>% #call df with specific dx_stage id 
        dplyr::filter(dataset == dat_set) %>%
        dplyr::filter(Direction == direc)
      go_sig_n[go_sig_idx, 1] <- stage
      go_sig_n[go_sig_idx, 2] <- dx
      go_sig_n[go_sig_idx, 3] <- dat_set
      go_sig_n[go_sig_idx, 4] <- direc 
      go_sig_n[go_sig_idx, 5] <- ifelse(direc == "down_regulated", dim(sig_data)[1]*-1,dim(sig_data)[1])
      #sigDatsum <- c(stage,dx, direc, dat_set, ifelse(direc == "down_regulated", dim(sig_data)[1]*-1,dim(sig_data)[1]))
      go_sig_idx = go_sig_idx + 1
    }
  }
}

# GO id venn diagrams
library(VennDiagram)
library(gridExtra)
go_venn <- lapply(c("down_regulated", "up_regulated"), function(direc){
  lapply(dx_stage, function(dx_stg){
    dei_go_list <- dplyr::filter(sig_go_terms[[dx_stg]], Direction == direc & dataset == "transcript")
    dei_go_id <- dei_go_list$ID
    deg_go_list <- dplyr::filter(sig_go_terms[[dx_stg]], Direction == direc & dataset == "gene")
    deg_go_id <- deg_go_list$ID
    if(length(dei_go_id) == 0 | length(deg_go_id) == 0){ plot.new()}
    else{
    venn.diagram(x = list(dei_go_id, deg_go_id),
                 category.names = c("DEI", "DEG"),
                 main = paste0(direc, "\n", " GO Terms at\n", dx_stg),
                 filename = NULL)
    }
  }) %>% setNames(dx_stage)
}) %>% setNames(c("downRegulated", "upRegulated")) %>%
  unlist(recursive = FALSE)

go_partitions <- lapply(c("down_regulated", "up_regulated"), function(direc){
  lapply(dx_stage, function(dx_stg){
    dei_go_list <- dplyr::filter(sig_go_terms[[dx_stg]], Direction == direc & dataset == "transcript")
    dei_go_id <- dei_go_list$ID
    deg_go_list <- dplyr::filter(sig_go_terms[[dx_stg]], Direction == direc & dataset == "gene")
    deg_go_id <- deg_go_list$ID
    if(length(dei_go_id) == 0 | length(deg_go_id) == 0){ plot.new()}
    else{
      get.venn.partitions(x = list(dei_go_id, deg_go_id))
                   
    }
  }) %>% setNames(dx_stage)
}) %>% setNames(c("downRegulated", "upRegulated")) %>%
  unlist(recursive = FALSE)


pdf(file = paste0(output_folder, "/", file_index, "_GO_venn_counts.pdf"), height = 6, width = 16)
for (dx_id in dxs){
  plot_id_down <- paste0("downRegulated", ".", dx_id)
  plots_down <- go_venn[grepl(plot_id_down, names(dei_deg_venn))]
  plot_id_up <- paste0("upRegulated", ".", dx_id)
  plots_up <- go_venn[grepl(plot_id_up, names(dei_deg_venn))]
  #grid.arrange(grobs = c(plots_down, plots_up))
  grid.arrange(grobs=lapply(c(plots_down, plots_up), grobTree), ncol=6)
}
dev.off()
file_index <- file_index + 1
go_dei_only <-lapply(dx_stage, function(dx_stg){
    dei_go_list <- dplyr::filter(sig_go_terms[[dx_stg]], Direction == direc & dataset == "transcript")
    dei_go_id <- dei_go_list$ID