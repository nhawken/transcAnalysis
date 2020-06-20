# t_dat <- t(datTrans)
#

datTrans <- datTrans_prenorm


zero_dat <- matrix(data = NA, nrow = dim(datTrans)[1], ncol=dim(datTrans)[2])
zero_dat[] <- vapply(datTrans, function(x){x==0}, logical(1))

colnames(zero_dat) <- colnames(datTrans)
rownames(zero_dat) <- rownames(datTrans)

zero_sums <- colSums(zero_dat)

datMeta_zeros <- cbind(datMeta, zero_sums)

# finding out which isoforms are 0 for each batch -------------------
batch_zeroed <- matrix(data = NA, nrow = dim(zero_dat)[1], ncol = length(unique(datMeta_zeros$batch)))
colnames(batch_zeroed) <- unique(datMeta_zeros$batch)
rownames(batch_zeroed) <- rownames(zero_dat)

idx <- 1
for(bat in unique(datMeta_zeros$batch)){
  batch_samples <- row.names(filter(datMeta_zeros, batch == bat))
  zero_dat_b <- zero_dat[ ,batch_samples] 
  zero_sum_b <- rowSums(zero_dat_b) / dim(zero_dat_b)[2]
  batch_zeroed[,idx] <- zero_sum_b
  idx <- idx+1
}


# sum number of zeros per batch -------------------------------------------------

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


# create bar graph of zero count isoforms/number of samples by batch (NOT USED) ------------
clrs <- as.factor(unique(datMeta$batch))
df_zero <- data.frame(zero_pct)

ggplot(df_zero, aes(rownames(df_zero),zero_pct, fill = clrs)) + 
  geom_col() + scale_fill_brewer(palette = "Set2") +
  theme(axis.text.x = element_blank(), legend.title = element_blank()) +
  xlab('Batch') + ylab('Zero Counts per Sample')

df_zero <- cbind(df_zero, datTrans_zero)
df_zero <- cbind(df_zero, n_samples)

# All the box plots -------------------------
#plot zeros by batch as a boxplot -------------------------------
  ggplot(datMeta_zeros, aes(x = batch, y = zero_sums , fill = batch)) +
  geom_boxplot() + 
  theme_classic() +
  theme(axis.text.x = element_blank())
 
  #plot zeros by differentiation day as a boxplot 
  ggplot(datMeta_zeros, aes(x = as.factor(Differentiation.day), y = zero_sums , fill = as.factor(Differentiation.day))) +
    geom_boxplot() + 
    theme_classic() +
    theme(axis.text.x = element_blank())
  
  #plot zeros by diagnosis as a boxplot
  ggplot(datMeta_zeros, aes(x = Dx, y = zero_sums , fill = Dx)) +
    geom_boxplot() + 
    theme_classic() +
    theme(axis.text.x = element_blank())
 


# plotting correlation between the 20 experimental PCs and the number of zero isoforms per sample (colored by batch )  
colnames(topPC.expr) <- paste0("ExprPC",1:pc_n)  
datMeta_zeroPC <- cbind(datMeta_zeros, topPC.expr) 
pdf(file = paste0(outputFolder,outputfile_counter,"_datMeta_PC_batch_zeros_6.pdf"),height = 6,width = 6)
# par(mfrow = c(2,2))
# par(mar = c(5,2,3,2))
for (i in grep("ExprPC", colnames(datMeta_zeroPC), value = T)) {
  A <- lm(datMeta_zeroPC[,i] ~ datMeta_zeroPC$zero_sums)
  p <- summary(A)$coefficient[,"Pr(>|t|)"][2]
  rsq <- summary(A)$adj.r.squared
  names(p) <- NULL
  names(rsq) <- NULL
  pc_0 <- ggplot(datMeta_zeroPC, aes(x = zero_sums, y = get(i) , col = batch )) +
    geom_point() +
    xlab('Zero Count Isoforms') +
    ylab(paste(i, 'Value')) +
    ggtitle(paste(i," p=", signif(p,2), " r^2=",signif(rsq,2)))
  print(pc_0)
  }
  dev.off()
  
#calculate average value of isoform expression per sample after norm and regression of covariates
avg_expr_reg <- colMeans(datTrans_reg)
names(avg_expr_reg) <- colnames(datTrans_reg)

datMeta_zero_pc_mean <- cbind(datMeta_zeroPC, avg_expr_reg)


#plotting correlation between the 20 experimental PCs and the average expression value per sample (colored by batch)
pdf(file = paste0(outputFolder,outputfile_counter,"_datMeta_PC_batch_avg_2.pdf"),height = 6,width = 6)
# par(mfrow = c(2,2))
# par(mar = c(5,2,3,2))
for (i in grep("ExprPC", colnames(datMeta_zero_pc_mean), value = T)) {
  A <- lm(datMeta_zero_pc_mean[,i] ~ datMeta_zero_pc_mean$avg_expr_reg)
  p <- summary(A)$coefficient[,"Pr(>|t|)"][2]
  rsq <- summary(A)$adj.r.squared
  names(p) <- NULL
  names(rsq) <- NULL
  pc_0 <- ggplot(datMeta_zero_pc_mean, aes(x = avg_expr_reg, y = get(i) , col = batch )) +
    geom_point() +
    xlab('Average Isoform Expression') +
    ylab(paste(i, 'Value')) +
    ggtitle(paste(i," p=", signif(p,2), " r^2=",signif(rsq,2)))
  print(pc_0)
}
dev.off()

  
