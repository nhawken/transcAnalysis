go_results <- lapply(names(dat_clusterProfiler), function(dx_day){
  lapply(c("up_regulated", "down_regulated"),function(direction){
    go_enrich <- lapply(c("BP","MF"), function(ontology){
      deg <- dat_clusterProfiler[[dx_day]]
      go_enrich_raw <- deg %>%
        dplyr::filter(p < 0.005) 
      if(direction == "up_regulated") {
        go_enrich_raw <- go_enrich_raw %>% 
          filter(beta > 0)
      } else if (direction == "down_regulated"){
        go_enrich_raw <- go_enrich_raw %>% 
          filter(beta < 0)
      } else {
        print("crash and burn!!")
      }
      go_enrich_raw <- go_enrich_raw %>% 
        dplyr::select(ensembl_gene_id) %>%
        unlist() %>%
        enrichGO(OrgDb = "org.Hs.eg.db",
                 universe = deg$ensembl_gene_id,
                 keyType = "ENSEMBL",
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
  setNames(names(dat_clusterProfiler))