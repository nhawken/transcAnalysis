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
      dplyr::filter(chromosome_name == chr & start_position >  plot_region[[1]] & end_position < plot_region[[2]]) %>% 
      dplyr::mutate(inCNV = ifelse(end_position > mut[[1]] & start_position < mut[[2]],T,F),
             significant = ifelse(p < 0.005,T,F), 
             hgnc_iso_symbol = ifelse(hgnc_iso_symbol == "", ensembl_transcript_id, hgnc_iso_symbol)) %>% 
      arrange(start_position) %>% 
      mutate(hgnc_iso_symbol = factor(hgnc_iso_symbol, levels = hgnc_iso_symbol))
    
  } else if (type == "gene"){
    dat_plot <- deg %>%
      dplyr::filter(chromosome_name == chr & start_position >  plot_region[[1]] & end_position < plot_region[[2]]) %>% 
      dplyr::mutate(inCNV = ifelse(hgnc_symbol == mut,T,F),
             significant = ifelse(p < 0.005,T,F), 
             hgnc_iso_symbol = ifelse(hgnc_iso_symbol == "", ensembl_transcript_id, hgnc_iso_symbol)) %>% 
      arrange(start_position) %>% 
      mutate(hgnc_iso_symbol = factor(hgnc_iso_symbol, levels = hgnc_iso_symbol))  
  } else {
    stop('type must be one of "cnv"/"gene"')
  }
  dat_plot
}


#dx_stage_loci <- dx_stage[!grepl(".*Idiopathic.*", dx_stage )]



defined_forms <-  grep(paste(names(loci),collapse = "|"), dx_stage, value = T) 
cnv_logfc_plots <- lapply(dx_stage,function(dx_stg){
  dx <- gsub("(.*?)_.*","\\1",dx_stg)
  stage <- gsub(".*?_(.*)","\\1",dx_stg)  
  lapply(names(loci), function(locus){
    type <- names(loci[[locus]])[3]
    de_results_adj %>%
      left_join(transAnno_id, by = "ensembl_transcript_id") %>%
      dplyr::select(contains(dx_stg), colnames(transAnno_id))  %>%
      setNames(gsub(paste0("\\.", dx_stg), "", colnames(.))) %>%
      get_cnv_logFC(loci[[locus]][["chr"]], loci[[locus]][["surround"]], loci[[locus]][[3]], type, stage) %>%
      mutate(dx_day = dx_stg,
             Dx = dx,
             stage = stage) %>%
      ggplot() +
      geom_hline(yintercept = 0, color = "red") +
      geom_pointrange(
        aes(
          x = hgnc_symbol,
          y = beta,
          ymin = beta - SE ,
          ymax = beta + SE,
          shape = significant,
          color = inCNV
        ),
        position = position_dodge2(
          width = 1,
          preserve = "single",
          padding = 0.8
        )
      )  +
      theme_classic(base_size = 12) +
      scale_shape_manual(values = c(1, 8)) +
      scale_color_manual(values = c("grey40", "royalblue")) +
      #facet_wrap(~hgnc_symbol, strip.position = 'bottom', scales = "free_x" ) +
      #geom_text(position, aes(x=hgnc_symbol, y=0, label = hgnc_symbol)) +
      theme(axis.text.x = element_text(
        angle = -90,
        hjust = 0,
        vjust = 0.5
      )) +
      labs(x = "")  +
      ggtitle(paste0(dx_stg, " at ", locus))    
  })
}) %>% setNames(dx_stage) 

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


cnv_sig_betas_p <- lapply(dx_stage, function(dx_stg){
  lapply(1:length(loci), function(locus){
    cnv_data <- cnv_logfc_plots[[dx_stg]][[locus]][["data"]] %>%
      dplyr::mutate(direction = ifelse(beta > 0, "up", "down")) %>%
      dplyr::mutate(location = names(loci[locus])) %>%
      dplyr::filter(p < 0.005)
  }) %>% setNames(names(loci))
}) %>% setNames(dx_stage)

deg_cnv <- lapply(1:length(loci), function(locus){
  cnv_transc <- cnv_sig_betas[[1]][[locus]][["ensembl_transcript_id"]]
  lapply(cnv_transc, function(id){
    dplyr::filter(de_results_adj, id %in% ensembl_transcript_id)   
  })
}) %>% setNames(names(loci))

de_results_adj_cnv <- lapply(1:length(loci), function(locus) {
  if (names(loci[[locus]])[3] == "cnv") {
    de_cnv <- de_results_adj_anno %>%
      dplyr::filter(chromosome_name == loci[[locus]][["chr"]] &
                      start_position >  loci[[locus]][["surround"]][[1]] &
                      end_position < loci[[locus]][["surround"]][[2]]) %>%
      dplyr::mutate(inCNV = ifelse(end_position > loci[[locus]][[3]][[1]] &
                                     start_position < loci[[locus]][[3]][[2]], T, F), 
                    hgnc_iso_symbol = ifelse(hgnc_iso_symbol == "",ensembl_transcript_id, hgnc_iso_symbol)) %>%
  arrange(start_position) %>%
  mutate(hgnc_iso_symbol = factor(hgnc_iso_symbol, levels = hgnc_iso_symbol))
  } else if (names(loci[[locus]])[3] == "gene"){
      de_cnv <- de_results_adj_anno %>%
        dplyr::filter(chromosome_name == loci[[locus]][["chr"]] &
                        start_position >  loci[[locus]][["surround"]][[1]] &
                        end_position < loci[[locus]][["surround"]][[2]]) %>%
        dplyr::mutate(inCNV = ifelse(hgnc_symbol == loci[[locus]][[3]],T,F), 
                      hgnc_iso_symbol = ifelse(hgnc_iso_symbol == "", ensembl_transcript_id, hgnc_iso_symbol)) %>% 
        arrange(start_position) %>% 
        mutate(hgnc_iso_symbol = factor(hgnc_iso_symbol, levels = hgnc_iso_symbol)) 
  } else {
    stop('type must be one of "cnv"/"gene"')
  }
}) %>% setNames(names(loci))


#create heatmaps of betas WORKING
de_beta_heat <- lapply(1:length(de_results_adj_cnv), function(locus){
  de_locus <- de_results_adj_cnv[[locus]] %>%
    arrange(start_position) %>% 
    select(starts_with("ensembl_transcript_id")| starts_with("beta.")) %>%
    na.omit() 
  names(de_locus) <- c("ensembl_transcript_id", dx_stage)
  de_pivot <- pivot_longer(de_locus, cols = -c(1), names_to = "dx_stage", values_to = "beta") %>%
      ggplot(aes(x = dx_stage, y = ensembl_transcript_id, fill = beta)) +
      geom_raster() +
      scale_fill_viridis_c() +
      theme(axis.text.x = element_text(angle = -90, hjust  = 0, vjust = 0.5)) +
      ggtitle(names(loci)[locus])
}) %>% setNames(names(loci))

#create heatmaps TESTING HGNC
de_beta_heat <- lapply(1:length(de_results_adj_cnv), function(locus){
  de_locus <- de_results_adj_cnv[[locus]] %>%
    arrange(start_position) %>% 
    select(starts_with("hgnc_iso_symbol")| starts_with("beta.")) %>%
    na.omit() 
  names(de_locus) <- c("hgnc_iso_symbol", dx_stage)
  de_pivot <- pivot_longer(de_locus, cols = -c(1), names_to = "dx_stage", values_to = "beta") %>%
    ggplot(aes(x = dx_stage, y = hgnc_iso_symbol, fill = beta)) +
    geom_raster() +
    scale_fill_viridis_c() +
    theme(axis.text.x = element_text(angle = -90, hjust  = 0, vjust = 0.5)) +
    ggtitle(names(loci)[locus])+
  coord_flip()
}) %>% setNames(names(loci))

pdf(paste0(output_folder, "/", file_index, "_loci_heatmap.pdf"), width = 14)
print(de_beta_heat)
dev.off()
file_index <- file_index + 1