t_dat <- t(datTrans)

zero_dat <- matrix(data = NA, nrow = dim(datTrans_pre)[1], ncol=dim(datTrans_pre)[2])
zero_dat[] <- vapply(datTrans_pre, function(x){x==0}, logical(1))

colnames(zero_dat) <- colnames(datTrans_pre)
rownames(zero_dat) <- rownames(datTrans_pre)

zero_sums <- colSums(zero_dat)

datMeta_zeros <- cbind(datMeta, zero_sums)

idx <- 1
for(bat in unique(datMeta_zeros$batch)){
  datMeta_b <- datMeta_zeros %>% 
    filter(batch == bat)
  datTrans_zero[idx] <- sum(datMeta_b$zero_sums)
  n_samples[idx] <- dim(datMeta_b)[1]
  idx <- idx+1
}
zero_pct <- datTrans_zero / n_samples
names(zero_pct) <- unique(datMeta_zeros$batch)

clrs <- as.factor(unique(datMeta$batch))
df_zero <- data.frame(zero_pct)

ggplot(df_zero, aes(rownames(df_zero),zero_pct, fill = clrs)) + 
  geom_col() + scale_fill_brewer(palette = "Set2") +
  theme(axis.text.x = element_blank(), legend.title = element_blank()) +
  xlab('Batch') + ylab('Zero Counts per Sample')

df_zero <- cbind(df_zero, datTrans_zero)
df_zero <- cbind(df_zero, n_samples)

datTrans_pre <- datTransRaw[idx,]