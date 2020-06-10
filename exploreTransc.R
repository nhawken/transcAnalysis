setwd( "C:/Users/natha/Documents/Geschwind Rotation 2020/transcAnalysis")
load('all_transcript_data.Rdata')
Meta_row <- nrow(datMeta)
Meta_col <- ncol(datMeta)

# Exploring the transcriptAnnoRaw data
en_tx_id <- transcriptAnnoRaw$ensembl_transcript_id

N_unique_tx_ids <- length(unique(en_tx_id))
N_all_tx_ids <- length(en_tx_id)

dup_en_id <- duplicated(en_tx_id)
dup_names <- transcriptAnnoRaw[dup_en_id,1]
dup_all2 <- sapply(dup_names, function(dup_names) {idx <- transcriptAnnoRaw$ensembl_transcript_id == dup_names
transcriptAnnoRaw[idx, ]})
ddd <- do.call(rbind, dup_all2)
dup_anno_set <- matrix(ddd, nrow = 50, byrow = TRUE)
colnames(dup_anno_set) <- names(transcriptAnnoRaw)

dup_frame2 <- data.frame(matrix(ncol = 11, nrow = 50))
colnames(dup_frame2) <- names(transcriptAnnoRaw)

for(ii in 1:25){
  dup_frame2[(ii*2-1),] <- dup_anno_set[ii,]
  dup_frame2[(ii*2), ] <- dup_anno_set[(ii+25),]
}
