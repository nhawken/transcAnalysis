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


#looking at DGCR8 gene --------------------
de_results_dgcr8 <- dplyr::filter(de_results_only_cnv, hgnc_symbol=="DGCR8")

dgcr8_logfc_plots <- lapply(dx_stage,function(dx_stg){
  dx <- gsub("(.*?)_.*","\\1",dx_stg)
  stage <- gsub(".*?_(.*)","\\1",dx_stg)  
  de_results_dgcr8 %>%
    dplyr::select(contains(dx_stg), colnames(transAnno_id))  %>%
    setNames(gsub(paste0("\\.", dx_stg), "", colnames(.))) %>%
    mutate(dx_day = dx_stg,
             Dx = dx,
             stage = stage, significant = ifelse(p < 0.005,T,F) ) %>%
      ggplot() +
      geom_hline(yintercept = 0, color = "red") +
      geom_pointrange(
        aes(
          x = ensembl_transcript_id,
          y = beta,
          ymin = beta - SE ,
          ymax = beta + SE,
          shape = significant
        ),
        position = position_dodge2(
          width = 1,
          preserve = "single",
          padding = 0.8
        )
      )  +
      theme_classic(base_size = 12) +
      scale_shape_manual(values = c(1, 8)) +
      #facet_wrap(~hgnc_symbol, strip.position = 'bottom', scales = "free_x" ) +
      #geom_text(position, aes(x=hgnc_symbol, y=0, label = hgnc_symbol)) +
      theme(axis.text.x = element_text(
        angle = -90,
        hjust = 0,
        vjust = 0.5
      )) +
      labs(x = "")  +
    ggtitle(paste0(dx_stg, " at ", "DGCR8"))
}) %>% setNames(dx_stage) 

grid.arrange(dgcr8_logfc_plots[[1]], dgcr8_logfc_plots[[2]], dgcr8_logfc_plots[[3]], dgcr8_logfc_plots[[4]] , nrow = 2)

pdf(paste0(output_folder, "/", file_index, "_dgcr8_logFC.pdf"), width = 14)
for (idx in 1:length(dxs)){
grid.arrange(dgcr8_logfc_plots[[(idx-1)*6+1]], dgcr8_logfc_plots[[(idx-1)*6+2]], 
             dgcr8_logfc_plots[[(idx-1)*6+3]], dgcr8_logfc_plots[[(idx-1)*6+4]],
             dgcr8_logfc_plots[[(idx-1)*6+5]], dgcr8_logfc_plots[[(idx-1)*6+6]], nrow = 2)
}

dev.off()
file_index <- file_index + 1


# looking all isoform expression for DGCR8 ------------------

datTrans <- datTrans_reg_batch

dgcr8_ens <- dplyr::filter(transAnno_id, hgnc_symbol == "DGCR8") %>%
  dplyr::select(ensembl_transcript_id) %>%
  flatten()

datTrans_dgcr8 <- lapply(dgcr8_ens, function(ens){ datTrans[ens,]}) %>%
  setNames(dgcr8_ens) %>%
  data.frame() %>%
  t() %>%
  data.frame()

sample_by_dx_day <- lapply(unique(datMeta$Dx), function(dx){ 
  dx_sample <- dplyr::filter(datMeta, Dx == dx)
  lapply(unique(datMeta$Differentiation.day), function(day){
    dx_day_sample <- dplyr::filter(dx_sample, Differentiation.day == day) %>%
      dplyr::select(SampleID)
    dx_day_sample_fix <- lapply(dx_day_sample, function(samp) {gsub("\\-","\\.", samp)}) %>%
      unlist() %>%
      unname()
  }) %>% setNames(unique(datMeta$Differentiation.day))
  })  %>% setNames(unique(datMeta$Dx))

datTrans_dgcr8_df <- lapply(names(sample_by_dx_day), function(dx){
  lapply(names(sample_by_dx_day[[dx]]), function(day){
    dat_d <- dplyr::select(datTrans_dgcr8, sample_by_dx_day[[dx]][[day]])
    # dat_dd <- dplyr::mutate(dat_d, ensembl_transcript_id = as.factor(rownames(dat_d)))
    # dat_sum <- rbind(dat_dd, c(colSums(dat_dd[,-ncol(dat_dd)]), "totExpr"))
  }) %>% setNames(names(sample_by_dx_day[[dx]]))
}) %>% setNames(names(sample_by_dx_day))

datTrans_dgcr8_percent <- lapply(names(sample_by_dx_day), function(dx){
  lapply(names(sample_by_dx_day[[dx]]), function(day){
    dat_d <- datTrans_dgcr8_df[[dx]][[day]]
    per_expr <- dat_d / colSums(dat_d)
  }) %>% setNames(names(sample_by_dx_day[[dx]]))
}) %>% setNames(names(sample_by_dx_day))


datTrans_dgcr8_mean_percent <- lapply(names(sample_by_dx_day), function(dx){
  lapply(names(sample_by_dx_day[[dx]]), function(day){
    dat_d <- datTrans_dgcr8_df[[dx]][[day]] 
    mean_expr <- rowMeans(dat_d)
    mean_sum_expr = mean_expr / sum(mean_expr)
  }) %>% setNames(names(sample_by_dx_day[[dx]]))
}) %>% setNames(names(sample_by_dx_day))
dgcr8_mean_per <- data.frame(unlist(datTrans_dgcr8_mean_percent, recursive = FALSE))

  
dgcr8_mean_per_table <- dplyr::mutate(dgcr8_mean_per, ensembl_transcript_id = as.factor(rownames(dgcr8_mean_per))) %>%
  pivot_longer(-ensembl_transcript_id, names_to = "SampleID", values_to = "PercentExpr") %>%
  separate(col = SampleID, c("Dx", "Day"), sep = "\\.")
dgcr8_mean_per_table$Dx <- gsub("X", "", dgcr8_mean_per_table$Dx)
dgcr8_mean_per_table$Day <- as.factor(dgcr8_mean_per_table$Day) %>%
  fct_relevel("25", "50", "75", "100")


#stacked bar plot based on percent expression of DGCR8
pdf(paste0(output_folder, "/", file_index, "_dgcr8_perexpr_bar.pdf"), width = 14, height = 10)
ggplot(dgcr8_mean_per_table, aes(fill = ensembl_transcript_id,
                                 y = PercentExpr, x = Day)) +
  geom_bar(position = "stack", stat = "identity") +
  facet_wrap(~Dx, ncol = 5) +
  scale_fill_brewer(palette = "Set2")
  theme(legend.text = element_text(size = 10))
dev.off()
file_index <- file_index + 1


#box plot of DGCR8 isoforms
dgcr8_exp_plots <- lapply(names(sample_by_dx_day),function(dx){
  lapply(names(sample_by_dx_day[[dx]]), function(day){
    dgcr8_data <- datTrans_dgcr8_df[[dx]][[day]]
    dplyr::mutate(dgcr8_data, ensembl_transcript_id = as.factor(rownames(dgcr8_data))) %>%
      
      pivot_longer(-ensembl_transcript_id, names_to = "SampleID", values_to = "expr") %>%
      ggplot() +
      geom_hline(yintercept = 0, color = "red") +
      geom_boxplot(
        aes( x = ensembl_transcript_id, y = expr
        # position = position_dodge2(
        #   width = 1,
        #   preserve = "single",
        #   padding = 0.8)
      ))  +
      theme_classic(base_size = 12) +
      #facet_wrap(~hgnc_symbol, strip.position = 'bottom', scales = "free_x" ) +
      #geom_text(position, aes(x=hgnc_symbol, y=0, label = hgnc_symbol)) +
      theme(axis.text.x = element_text(
        angle = -90,
        hjust = 0,
        vjust = 0.5
      )) +
      labs(x = "")  +
      ggtitle(paste0(dx, " at day ", day, " at ", "DGCR8"))
    }) %>% setNames(names(sample_by_dx_day[[dx]]))
  }) %>% setNames(names(sample_by_dx_day))

pdf(paste0(output_folder, "/", file_index, "_dgcr8_expr_box.pdf"), width = 14)
for (idx in names(sample_by_dx_day)){
    grid.arrange(dgcr8_exp_plots[[idx]][["25"]], dgcr8_exp_plots[[idx]][["50"]], dgcr8_exp_plots[[idx]][["75"]], dgcr8_exp_plots[[idx]][["100"]], nrow = 1)
  }
dev.off()
file_index <- file_index + 1

# plotting all timepoint cnvs -----------------

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
    na.omit() %>%
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
    separate(col = dx_stage, sep = "_",into = c("dx", "stage")) %>%
    ggplot(aes(x = stage, y = hgnc_iso_symbol, fill = beta)) +
    geom_raster() +
    scale_fill_gradient2(midpoint=0, low="blue", mid = "white",
                          high="red", space ="Lab" ) +
    theme(axis.text.x = element_text(angle = -90, hjust  = 0, vjust = 0.5)) +
    facet_grid(~ dx) +
    ggtitle(names(loci)[locus]) 
  # coord_flip()
}) %>% setNames(names(loci))

pdf(paste0(output_folder, "/", file_index, "_loci_heatmap.pdf"), width = 14)
print(de_beta_heat)
dev.off()
file_index <- file_index + 1

# determine correlation without the loci


de_results_only_cnv <- lapply(1:length(loci), function(locus) {
  if (names(loci[[locus]])[3] == "cnv") {
    de_cnv <- de_results_adj_anno %>%
      dplyr::filter(chromosome_name == loci[[locus]][["chr"]] &
                      end_position > loci[[locus]][[3]][[1]] &
                      start_position < loci[[locus]][[3]][[2]]) %>%
      dplyr::mutate(hgnc_iso_symbol = ifelse(hgnc_iso_symbol == "", ensembl_transcript_id, hgnc_iso_symbol)) %>%
      arrange(start_position) %>%
      mutate(hgnc_iso_symbol = factor(hgnc_iso_symbol, levels = hgnc_iso_symbol))
  } else if (names(loci[[locus]])[3] == "gene"){
    de_cnv <- de_results_adj_anno %>%
      dplyr::filter(chromosome_name == loci[[locus]][["chr"]] &
                    hgnc_symbol == loci[[locus]][[3]]) %>%
      dplyr::mutate(hgnc_iso_symbol = ifelse(hgnc_iso_symbol == "", ensembl_transcript_id, hgnc_iso_symbol)) %>% 
      arrange(start_position) %>% 
      mutate(hgnc_iso_symbol = factor(hgnc_iso_symbol, levels = hgnc_iso_symbol)) 
  } else {
    stop('type must be one of "cnv"/"gene"')
  }
}) %>% setNames(names(loci))

de_results_only_cnv <- do.call(rbind.data.frame, de_results_only_cnv) %>%
  unique()

cnv_matches <- match(de_results_only_cnv$ensembl_transcript_id, de_results_adj_anno$ensembl_transcript_id)
de_results_without_cnv <- de_results_adj_anno[-cnv_matches, -c(218:228) ]

de_dx_day_without_cnv <- lapply(dx_stage, function(dx_stg){
  de_results_without_cnv %>% 
    dplyr::select(ensembl_transcript_id,contains(dx_stg))  %>% #keep only the variables ID, and dx_stage of interest (beta, se, p, fdr)
    setNames(gsub(paste0("\\.",dx_stg), "", colnames(.))) %>% #remove .dx_stage identifier from the columns of beta, se, p, fdr
    left_join(transAnno_id, by = "ensembl_transcript_id") #add transcript annotation table to the right edge of the table
})  %>% setNames(dx_stage) #rename overall df to be the dx_stage identifier


# Correlation of logFC within Dx  --------------------------------------------------------------------------------
dx_logFC_cor <- lapply(dxs, function(dx){
  idx <- grep(paste0(dx,"_\\d+"), names(de_dx_day_without_cnv), value = T)
  
  betas <- lapply(names(de_dx_day_without_cnv[idx]), function(dx_day){
    deg <- de_dx_day_without_cnv[[dx_day]][]
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

pdf(paste0(output_folder,"/", file_index, "_logFC_cor_within_Dx_woCNV.pdf"), height = 6, width = 12)
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


pdf(paste0(output_folder,"/", file_index, "_logFC_cor_within_Dx_network_woCNV.pdf"), height = 7, width = 12)
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
  idx <- grep(stg, names(de_dx_day_without_cnv), value = T)
  betas <- lapply(names(de_dx_day_without_cnv[idx]), function(dx_day){
    deg <- de_dx_day_without_cnv[[dx_day]][]
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
png(paste0(output_folder,"/", file_index, "_logFC_cor_between_Dx_woCNV.png"), height = 8, width = 16, units = "in", res = 300)
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

pdf(paste0(output_folder,"/", file_index, "_logFC_cor_between_Dx_network_woCNV.pdf"), height = 7, width = 12)
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

betas <- lapply(names(de_dx_day_without_cnv), function(dx_day){
  deg <- de_dx_day_without_cnv[[dx_day]][]
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
pdf(paste0(output_folder,"/", file_index, "_logFC_cor_all_points_network_woCNV.pdf"), height = 12, width = 12)
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
pdf(paste0(output_folder,"/", file_index, "_logFC_cor_all_points_heatmap_woCNV.pdf"), height = 8, width = 10)
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
all_degs <- lapply(de_dx_day_without_cnv, function(x){
  x %>% 
    dplyr::filter(p < 0.005) %>% 
    dplyr::select(ensembl_transcript_id) %>% 
    unlist()
})


pdf(paste0(output_folder,"/", file_index, "_upset_byDx_woCNV.pdf"),  height = 8, width = 16)#, units = "in", res = 300)
for (dx in dxs) {
  idx <- grep(paste0(dx,"_\\d+"), names(de_dx_day_without_cnv), value = T)
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


pdf(paste0(output_folder,"/", file_index, "_upset_byDay_woCNV.pdf"),  height = 8, width = 16)
for (stg in stages) {
  idx <- grep(stg, names(de_dx_day_without_cnv), value = T)
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


