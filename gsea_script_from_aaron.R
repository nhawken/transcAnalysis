options(stringsAsFactors = FALSE)
library(tidyverse)
library(fgsea)
library(biomaRt)
library(patchwork)
library(ggpubr)
setwd("/u/project/geschwind/aarong/CIRM_iPSC")
#setwd("/Volumes/hoffman2/CIRM_iPSC")
wd <- "analysis/rsem/3_Dx/1.2_all_batches_seqPCS_dropQC_source/LME/dx_day/"
output_folder <- "GSEA"
dir.create(file.path(wd, output_folder))
de_results <- read_csv(file.path(wd,"dx_day_lme_DE.csv"))[,-2] %>%
  rename("ensembl_gene_id" = "X1")
load("analysis/rsem/gene_annotation.Rdata")
dx_stage <- gsub(".*\\.(.*?_.*)","\\1", colnames(de_results)) %>%
  unique() %>%
  sort() %>%
  .[. != "ensembl_gene_id"]
de_dx_day <- lapply(dx_stage, function(dx_stg){
  de_results %>%
    dplyr::select(ensembl_gene_id,contains(dx_stg))  %>%
    setNames(gsub(paste0("\\.",dx_stg), "", colnames(.)))
})  %>% setNames(dx_stage)
human_mart = useMart("ensembl",dataset = "hsapiens_gene_ensembl")
getinfo <- c( "ensembl_gene_id", "entrezgene_id")
geneDat <- getBM(attributes = getinfo,filters = "ensembl_gene_id",
                 values = de_results$ensembl_gene_id,mart = human_mart)
save(geneDat, file = file.path(wd, output_folder, "genedat.rdata"))
go_pathways <- gmtPathways("/u/project/geschwind/aarong/DBs/GSEA/v7/c5.all.v7.0.entrez.gmt")
fgsea_res <- map(de_dx_day, function(de){
  de_ranked <- de %>%
    arrange(beta) %>%
    left_join(geneDat, by = "ensembl_gene_id") %>%
    dplyr::select(entrezgene_id,beta) %>%
    drop_na()
  ranked_stat <- de_ranked$beta
  names(ranked_stat)<- de_ranked$entrezgene_id
  fgseaRes <- fgsea(go_pathways, stats = ranked_stat, nperm = 1000000,minSize = 30, maxSize = 500)
})
save(fgsea_res, file = file.path(wd,output_folder,"gsea_results.rdata"))
load(file.path(wd,output_folder,"gsea_results.rdata"))
plot_data <- map(fgsea_res, function(res){
  res %>%
    mutate(direction = ifelse(NES >0 ,"up_regulated","down_regulated"),
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
go_plots <- map(names(plot_data) , function(dx){
  plot_data[[dx]] %>%
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
plot_legend <- get_legend(go_plots[[1]]) %>%
  as_ggplot()
pdf(file.path( wd,output_folder,"gsea_results.pdf"), width = 24, height = 14)
plot_command <- paste0("go_plots[[", 1:length(go_plots),"]] + theme(legend.position = 'none')", collapse = "+") %>%
  paste0("(", .,") + plot_spacer() +plot_legend + plot_layout(nrow = 3, ncol = 3, heights = c(1,1,0.05))")
eval(parse(text = plot_command))
dev.off()