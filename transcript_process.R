#Aaron Gordon
#Run data normalization and QC on CIRM batches 1-5 RNAseq data


# NATALIE
# you will need to change some of the vaibale names (i.e. datExpr to datTrans)


# set up R and some variables ----------------------------
setwd("C://Users/natha/Documents/Geschwind Rotation 2020/transcAnalysis")
options(stringsAsFactors = FALSE)
library(cowplot)
library(data.table)
library(WGCNA) #had to install via BiocManager
library(plyr)
library(cqn) #had to install via BiocManager
library(tidyverse)
library(stringr)
library(gridExtra)
library(corrplot)
library(wesanderson)
library(patchwork)

outputfile_counter = 1
outputFolder = paste(getwd(), '/process_OUT_20200615/', sep = "") # SET AN OUTPUT PATH HERE
dir.create(outputFolder, showWarnings = F, recursive = T)
#time0 <- proc.time()

# Functions -------------------------------------------------------------------------------------------------------
# View data function edited to change "expr" to "trans" for data variable
viewData <- function(meta_data, trans_data, log_transform = "F", cvrs = c("Dx","batch","Differentiation.day","Sex")){
  clrs <- as.factor(meta_data$Dx) #what is this (turn str into number)
  if (log_transform) {
    trans_data <- log2(trans_data + 1) #why plus one?
  }
  boxplot(trans_data,range = 0, main = paste("Boxplot Counts"),xaxt = "n",
          col = "white", medcol = clrs, whiskcol = clrs, staplecol = clrs, boxcol = clrs)
  axis(1,at=c(1:dim(trans_data)[2]), labels = meta_data$Dx, las = 2, cex.axis = 0.6)
  plot(density(trans_data[,1]),col = as.factor(meta_data$Dx)[1], main = paste("Density Counts Gene"),ylim = c(0, 0.35))
  for (i in 2:dim(trans_data)[2]) {
    lines(density(trans_data[,i]), col = as.factor(meta_data$Dx)[i])
  }
  mdsG = cmdscale(dist(t(trans_data)),eig = TRUE)
  pc1 = mdsG$eig[1]^2 / sum(mdsG$eig^2)  
  pc2 = mdsG$eig[2]^2 / sum(mdsG$eig^2)
  mdsPlots <- cbind(meta_data,as.data.frame(mdsG$points)[match(meta_data$SampleID,rownames(mdsG$points)),])
  colnames(mdsPlots)[c((ncol(mdsPlots)-1):ncol(mdsPlots))] <- c("MDS1","MDS2")
  for (covar in cvrs){
    p1 <-  ggplot(mdsPlots,aes_string(x="MDS1",y="MDS2",color=covar)) +
      geom_point(size = 3) +
      theme_minimal() + 
      theme(legend.position="bottom") +
      ggtitle(covar)
    print(p1)
  }
}

outlierAnalysis <- function(datExpr,datMeta,title){
  sdout <- 2; normadj <- (0.5+0.5*bicor(datExpr, use='pairwise.complete.obs'))^2
  netsummary <- fundamentalNetworkConcepts(normadj); 
  K <- netsummary$Connectivity; Z.K <- (K-mean(K))/sqrt(var(K))
  C <- netsummary$ClusterCoef; Z.C = (C - mean(C))/sqrt(var(C))
  outliers <- (Z.K > mean(Z.K)+sdout*sd(Z.K))|(Z.K < mean(Z.K)-sdout*sd(Z.K))
  # cat(paste("There are ",sum(outliers)," outliers samples based on a bicor distance sample network connectivity standard deviation above ",sdout,sep="")); 
  # cat("\n")
  print(table(outliers))
  # print(datMeta %>% 
  #   filter(SampleID %in% colnames(datExpr)[outliers]) %>% 
  #   select(SampleID, Individual, CellLine,Induction, Day.grouped) 
  # )
  # cat(colnames(datExpr)[outliers]); 
  cat("\n")
  outlierIDs <- (rep(" ", length(outliers))); 
  outlierIDs[outliers] <- levels(datMeta$Dx_day)[as.numeric(datMeta$Dx_day[outliers])]
  plot(Z.K, col = as.numeric(datMeta$Dx), pch=as.numeric(as.factor(datMeta$batch))+14, main=paste("Outlier detection for",title),xaxt="n", ylab="Network connectivity (z score)");text(1:length(Z.K),Z.K,label=outlierIDs,pos=3,cex=0.6);abline(h=-2, lty=2)
  #axis(1, at=1:123, ,labels=datMeta$Day.grouped,las=2)
  return(outliers)
}

outlierAnalysisLoop <- function(datMeta, datExpr, covar){
  outliers <- vector()
  for(cvr in sort(unlist(unique(datMeta[,covar])))){
    cat(cvr)
    cat("\n")
    datMeta_covar <- filter(datMeta,!!as.name(covar)==cvr)
    datExpr_covar <- datExpr[,match(datMeta_covar$SampleID,colnames(datExpr))]
    out <- outlierAnalysis(datExpr_covar,datMeta_covar,cvr)
    outliers <- c(outliers,datMeta_covar$SampleID[out])
  }
  print(datMeta %>% 
          filter(SampleID %in% outliers) %>% 
          select(SampleID, Individual, CellLine,Induction, covar,Dx), n = Inf)
  return(outliers)
}

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) { 
  ## Useful function for comparing multivariate data
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  if (is.numeric(x) & is.numeric(y)) {
    r <- abs(cor(x, y,use = "pairwise.complete.obs",method = "spearman"))
  } else {
    r = abs(summary(lm(y~x))$adj.r.squared)
  }
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if (missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}

lm_vec <- function(matrix1, matrix2){
  lapply(as.data.frame(matrix1), function(v1){
    lapply(as.data.frame(matrix2), function(v2){
      summary(lm(v1~v2))$"adj.r.squared"
    }) %>% do.call(c,.)
  }) %>% do.call(rbind, .)
}

QCheatMap <-  function(datQC,datMeta,covars,outputfolder,outputfile_counter){
  covarAnovasPval <- matrix(NA, nrow = dim(datQC)[2],ncol = length(covars) )
  rownames(covarAnovasPval) <- colnames(datQC)
  colnames(covarAnovasPval) <- covars
  for(covar in covars) {
    for(i in colnames(datQC)) {
      A = anova(lm(as.numeric(datQC[,i]) ~ datMeta[,covar])); p = A$"Pr(>F)"[1]; 
      covarAnovasPval[i,covar] <- p
    }
  }
  # correct for multiple testing
  covarAnovasFDR <- matrix(p.adjust(covarAnovasPval,method = "BH"),nrow=nrow(covarAnovasPval),ncol=ncol(covarAnovasPval))  
  rownames(covarAnovasFDR) <- rownames(covarAnovasPval)
  colnames(covarAnovasFDR) <- colnames(covarAnovasPval)
  
  # prepare data for heat map
  dispMat <- -log10(covarAnovasFDR)
  dispMat[covarAnovasFDR > 0.05] <- 0
  tMat <- matrix(signif(covarAnovasFDR,1),nrow = nrow(dispMat),ncol = ncol(dispMat))
  tMat[covarAnovasFDR > 0.05] <- ""
  covarAnovasFDR[is.na(covarAnovasFDR)] <- 1
  tMat[is.na(covarAnovasFDR)] <- ""
  #plot heatmap
  pdf(paste0(outputfolder,outputfile_counter,"_QCheatmap.pdf"),height = 8, width = 10)
  par(mar = c(5.1,15.1,4.1,2.1))
  labeledHeatmap(Matrix = dispMat,
                 xLabels = colnames(dispMat),
                 xColorLabels = TRUE,
                 xSymbols = colnames(dispMat),
                 yLabels = rownames(dispMat),
                 colors = blueWhiteRed(100)[51:100],
                 cex.lab.x = 1.1,
                 cex.lab.y = 0.8,
                 zlim = c(0, 10),
                 textMatrix = tMat,
                 main = "Anova FDR values",
                 cex.text = 0.4,
                 setStdMargins = FALSE)
  dev.off()
}


#Load Expression, Meta and QC data  and format them ----------------------------------------
#LOAD DATA HERE
load('all_transcript_data.Rdata')


# (2) Filter genes with less than 10 counts in 60% of  samples per ID ------------
if(TRUE){
  gene_filter_list<- lapply(unique(datMeta$Day.grouped), function(d){
    
    filter_conds <- matrix(NA, nrow = 21, ncol = 10)
    datMetaDxd <- datMeta %>% 
      filter(Day.grouped == d)
    datTransDxd <- datTrans[,datMetaDxd$SampleID,drop=F]
    
    for (i in 0:20) {
      for (j in seq(0.1,1,0.1)) {
        
        pres <- apply(datTransDxd > i,1,sum) #1 indicates function is applied over rows
        idx <- length(which(pres >= j*ncol(datTransDxd)))
        filter_conds[i + 1, j *10] <- idx
        
      }
    }
    return(filter_conds)
  })
  
  gene_filter_tidy <- gene_filter_list %>% 
    setNames(unique(datMeta$Day.grouped)) %>% 
    lapply(function(m1){
      m1 %>% as.data.frame() %>% 
        set_names(paste0(seq(10,100,10),"%")) %>% 
        mutate(min_num_expressing = 0:20) %>% #NOT SURE WHAT THIS DOES
        gather("fraction", "gene_num", matches("%")) %>% 
        mutate(fraction = factor(fraction, levels = paste0(seq(0,100,10),"%"))) 
    }) %>% 
    bind_rows(.id = "Day")
  
  pdf(paste0(outputFolder,outputfile_counter,"_genefilter_threshold.pdf"), width = 14, height = 10) #what determines width/height inches!
  ggplot(gene_filter_tidy, aes(x = min_num_expressing, y= gene_num, color = fraction)) +
    geom_point() +
    geom_line(aes(group = fraction)) + 
    facet_wrap(~Day) +
    # scale_y_continuous(breaks = seq(0,28000,2000)) +
    scale_x_continuous(breaks=seq(0,20,2)) +
    theme_bw(base_size = 20)
  dev.off() #what is this?? device off, add tail to end of pdf file
  
}

idx <- vector()
for(day in unique(datMeta$Day.grouped)){
  datMetaDxd <- datMeta %>% 
    filter(Day.grouped == day)
  datTransDxd <- datTransRaw[,datMetaDxd$SampleID,drop=F]
  pres <- apply(datTransDxd>10,1,sum)
  idx <- union(idx, which(pres > 0.3*dim(datTransDxd)[2]))
}
datTrans <- datTransRaw[idx,]
# (3) View Data Pre-Normalization and Pre-QC -------------------------------
## raw statistics

viewData <- function(meta_data, trans_data, log_transform = "F", cvrs = c("Dx","batch","Differentiation.day","Sex")){
  clrs <- as.factor(meta_data$Dx) #what is this (turn str into number)
  if (log_transform) {
    trans_data <- log2(trans_data + 1) #why plus one?
  }
  boxplot(trans_data,range = 0, main = paste("Boxplot Counts"),xaxt = "n",
          col = "white", medcol = clrs, whiskcol = clrs, staplecol = clrs, boxcol = clrs)
  axis(1,at=c(1:dim(trans_data)[2]), labels = meta_data$Dx, las = 2, cex.axis = 0.6)
  plot(density(trans_data[,1]),col = as.factor(meta_data$Dx)[1], main = paste("Density Counts Gene"),ylim = c(0, 0.35))
  for (i in 2:dim(trans_data)[2]) {
    lines(density(trans_data[,i]), col = as.factor(meta_data$Dx)[i])
  }
  mdsG = cmdscale(dist(t(trans_data)),eig = TRUE)
  pc1 = mdsG$eig[1]^2 / sum(mdsG$eig^2)  
  pc2 = mdsG$eig[2]^2 / sum(mdsG$eig^2)
  mdsPlots <- cbind(meta_data,as.data.frame(mdsG$points)[match(meta_data$SampleID,rownames(mdsG$points)),])
  colnames(mdsPlots)[c((ncol(mdsPlots)-1):ncol(mdsPlots))] <- c("MDS1","MDS2")
  for (covar in cvrs){
    p1 <-  ggplot(mdsPlots,aes_string(x="MDS1",y="MDS2",color=covar)) +
      geom_point(size = 3) +
      theme_minimal() + 
      theme(legend.position="bottom") +
      ggtitle(covar)
    print(p1)
  }
}

pdf(file=paste0(outputFolder,outputfile_counter,"_RawStatistics_80pct.pdf"),width=12,height=12)
viewData(datMeta, datTrans, log_transform = T)
dev.off()
outputfile_counter<-outputfile_counter+1



# (4) Get Gene Annotaion --------------------------------------------
# library(biomaRt)
# getinfo <- c("ensembl_gene_id","hgnc_symbol","chromosome_name","band","strand","start_position", "end_position","gene_biotype","transcript_length","percentage_gene_gc_content")
# humanMart <- useMart(biomart="ensembl",dataset="hsapiens_gene_ensembl") #change if working with differnet spieces
# geneAnnoRaw <- getBM(attributes = getinfo,filters=c("ensembl_gene_id"),values=rownames(datExpr),mart=humanMart)
# save(geneAnnoRaw,file="./data/GeneAnnotation_batch4.rda")

transAnno <-  transcriptAnnoRaw %>%
  filter(ensembl_transcript_id %in% rownames(datTrans)) %>%
  filter(!duplicated(ensembl_transcript_id))
  #filter(gene_biotype != "rRNA") %>% 

#parts of the dataset that do not have an annotation
#Ensembl transcript IDs are deprecated 
datTransMissing <- datTrans[!(rownames(datTrans) %in% transAnno$ensembl_transcript_id),]

#all trans data with active annotations
datTrans <- datTrans[match(transAnno$ensembl_transcript_id,rownames(datTrans)),]

stopifnot(all(transAnno$ensembl_transcript_id==rownames(datTrans)))

rm('datMetaDxd','datTransDxd','datTransMissing','datTransRaw','transcriptAnnoRaw','pres')
# (5) Normalize ------------------------------------------------------
#datTrans_preNorm <- datTrans

#why CQN? are there other normalization options? what are the sqn and verbose options?
gc()
datTransCQN <- cqn(datTrans,x=transAnno$percentage_gene_gc_content, lengths= transAnno$transcript_length,sqn=F ,verbose = T)
RPKM.cqn <- datTransCQN$y + datTransCQN$offset #does RPKM mean something? does .htg mean something?
datTrans.htg<-RPKM.cqn

## Check difference in relationship to GC before and after normalization
preNorm <- datTrans
postNorm <- datTrans.htg
keepisoforms <- intersect(rownames(preNorm),rownames(postNorm)) #why would the row names not all intersect?

preNorm <- preNorm[match(keepisoforms,rownames(preNorm)),]
postNorm <- postNorm[match(keepisoforms,rownames(postNorm)),]

#does normalization remove degrees of freedom? no!

qualMat <- matrix(NA,nrow=ncol(preNorm),ncol=4)
colnames(qualMat) <- c("pre.Norm.GC.cor","pre.Norm.Length.cor","post.Norm.GC.cor","post.Norm.Length.cor")

for (i in 1:nrow(qualMat)) {
  qualMat[i,1] <- cor(preNorm[,i],transAnno$percentage_gene_gc_content,method="spearman")
  qualMat[i,2] <- cor(preNorm[,i],transAnno$transcript_length,method="spearman")
  qualMat[i,3] <- cor(postNorm[,i],transAnno$percentage_gene_gc_content,method="spearman")
  qualMat[i,4] <- cor(postNorm[,i],transAnno$transcript_length,method="spearman")
} 

pdf(file=paste0(outputFolder,outputfile_counter,"_GC_noQN_length_correlations_80pct.pdf"),width=8,height=8)
par(mfrow=c(2,2))
hist(qualMat[,1],main=colnames(qualMat)[1],xlim=c(-0.5,0.2),ylim=c(0,80),breaks=seq(-0.5,0.2,by=0.01),xlab="Spearman's rho across samples")
abline(v=0)
hist(qualMat[,3],main=colnames(qualMat)[3],xlim=c(-0.3,0.2),ylim=c(0,80),breaks=seq(-0.3,0.2,by=0.01),xlab="Spearman's rho across samples")
abline(v=0)
hist(qualMat[,2],main=colnames(qualMat)[2],xlim=c(-0.2,0.8),ylim=c(0,60),breaks=seq(-0.1,0.8,by=0.01),xlab="Spearman's rho across samples")
abline(v=0)
hist(qualMat[,4],main=colnames(qualMat)[4],xlim=c(-0.2,0.8),ylim=c(0,60),breaks=seq(-0.1,0.8,by=0.01),xlab="Spearman's rho across samples")
abline(v=0)
dev.off()
outputfile_counter<-outputfile_counter+1

# # (6) Remove Outliers --------------------------------------------------------------------
# sink(file=paste0(outputFolder,outputfile_counter,".1_Outliers.txt"))
# pdf(file=paste0(outputFolder,outputfile_counter,"_SampleNetwork_allSamps.pdf"),height = 8,width=12)
# 
# ## For All Samples
# cat("All samples \n")
# all_days <- outlierAnalysis(datExpr,datMeta,"all samples")
# 
# print(datMeta %>% 
#         filter(all_days) %>% 
#         select(SampleID, Individual, CellLine,Induction, Day.grouped,Dx) )
# 
# #By differentiation Day
# cat ("\n.........................................................................\n\n")
# cat("By differentiation day\n")
# par(mfrow=c(2,2))
# outliers_days <- outlierAnalysisLoop(datMeta, datExpr, "Day.grouped")
# 
# #By Dx
# cat ("\n.........................................................................\n\n")
# cat("By Dx\n")
# par(mfrow=c(3,4))
# outliers_dx <- outlierAnalysisLoop(datMeta, datExpr, "Dx")
# 
# #By batch
# cat ("\n.........................................................................\n\n")
# cat("By batch and Dx\n")
# par(mfrow=c(3,4))
# datMeta1 <- datMeta %>% 
#   mutate(bd = paste(batch, Day.grouped, sep = "_"))
# outliers_batch_dx <- outlierAnalysisLoop(datMeta1, datExpr, "bd")
# sink()
# dev.off()
# outputfile_counter=outputfile_counter+1
# 
# backupDM <- datMeta
# backupDE <- datExpr
# backupDQ <- datQC
# 
# #remove outliers from datMeta, datExpr and datQC
# outliers <- datMeta$SampleID %in% outliers_dx
# datMeta <- datMeta[!outliers, ]
# datExpr.htg <- datExpr.htg[, !outliers]
# datQC <- datQC[!outliers, ]

# (7) View  Data Post-Normalization  ----------------------------------------------------------------
## View basic statistics
datTrans <- datTrans.htg

viewData <- function(meta_data, trans_data, log_transform = "F", cvrs = c("Dx","batch","Differentiation.day","Sex")){
  clrs <- as.factor(meta_data$Dx)
  if (log_transform) {
    trans_data <- log2(trans_data + 1)
  }
  boxplot(trans_data,range = 0, main = paste("Boxplot Counts"),xaxt = "n",
          col = "white", medcol = clrs, whiskcol = clrs, staplecol = clrs, boxcol = clrs)
  axis(1,at=c(1:dim(trans_data)[2]), labels = meta_data$Dx, las = 2, cex.axis = 0.6)
  plot(density(trans_data[,1]),col = as.factor(meta_data$Dx)[1], main = paste("Density Counts Gene"),ylim = c(0, 0.35))
  for (i in 2:dim(trans_data)[2]) {
    lines(density(trans_data[,i]), col = as.factor(meta_data$Dx)[i])
  }
  mdsG = cmdscale(dist(t(trans_data)),eig = TRUE)
  pc1 = mdsG$eig[1]^2 / sum(mdsG$eig^2)  
  pc2 = mdsG$eig[2]^2 / sum(mdsG$eig^2)
  mdsPlots <- cbind(meta_data,as.data.frame(mdsG$points)[match(meta_data$SampleID,rownames(mdsG$points)),])
  colnames(mdsPlots)[c((ncol(mdsPlots)-1):ncol(mdsPlots))] <- c("MDS1","MDS2")
  for (covar in cvrs){
    p1 <-  ggplot(mdsPlots,aes_string(x="MDS1",y="MDS2",color=covar)) +
      geom_point(size = 3) +
      theme_minimal() + 
      theme(legend.position="bottom") +
      ggtitle(covar)
    print(p1)
  }
}

pdf(file=paste0(outputFolder,outputfile_counter,"_ProcessedStatistics_80pct.pdf"))
viewData(datMeta, datTrans.htg)
dev.off()
outputfile_counter=outputfile_counter+1

# (8) Principal Component Analysis per batch ---------------------------------------------------------------
pc_n <- 20

## Get the first 5 PCs in the expr data
thisdat.expr <- t(scale(t(datTrans.htg), scale = F))  #not sure what the transpose, scale, transpose is about
## Centers the mean of all genes - this means the PCA gives us the eigenvectors of the geneXgene covariance matrix, allowing us to assess the proportion of variance each component contributes to the data
PC.expr <- prcomp(thisdat.expr)
topPC.expr <- PC.expr$rotation[,1:pc_n]
varexp <- (PC.expr$sdev)^2 / sum(PC.expr$sdev^2)
topvar <- varexp[1:pc_n]
colnames(topPC.expr) <- paste("expr\n",colnames(topPC.expr)," (",signif(100*topvar[1:pc_n],2),"%)",sep="")

## perform a PCA of the sequencing statistics
# large_qc <- which(apply(datQC, 2, mean) > 10 ^ 4)
# datQC[, large_qc] = log10(datQC[, large_qc])
# PC.datQC <- prcomp(na.omit(t(scale((datQC),scale = T))), center = T)
# 
# 
# varexp <- (PC.datQC$sdev) ^ 2 / sum(PC.datQC$sdev ^ 2)
# 
# topPC.datQC <- PC.datQC$rotation[,1:pc_n]
topPC.datQC <- datMeta[,paste0("SeqPC",1:pc_n)]


colnames(topPC.datQC) <- paste0("SeqPC",1:pc_n)
seq_qc_plot_data <- data.frame(var_exp = varexp) %>% 
  rowid_to_column("SeqPC") %>% 
  mutate(cumulative_varexpr = cumsum(varexp)) 
seq_qc_plot1 <- ggplot(seq_qc_plot_data,aes(x = SeqPC, y = varexp)) +
  geom_point() +
  theme_minimal()

seq_qc_plot2 <-  ggplot(seq_qc_plot_data,aes(x = SeqPC, y = cumulative_varexpr)) +
  geom_point() +
  theme_minimal()

pairsdat <- datMeta[,c("Dx","Day.grouped","Sex","CellLine", "batch","Source.Tissue","ethnicityPC1","ethnicityPC2")]
pairsdat[,c(1:5)] <- lapply(pairsdat[,c(1:5)],as.factor)
pairsdat <-  pairsdat[,sapply(pairsdat, function(covar){length(unique(covar))}) > 1] #remove cols with no variance

# pdf(file=paste0(outputFolder,outputfile_counter,"_SeqStats_Comparison.pdf"),height=20,width=24)
# pairs(cbind(pairsdat,topPC.expr,topPC.datQC),col=as.factor(datMeta$Dx),pch=19,upper.panel = panel.cor,main="Sequencing Statistics, datMeta, and expr Comparison |Spearman's rho|")
# pairs(cbind(pairsdat,topPC.datQC),col=as.factor(datMeta$Dx),pch=19,upper.panel = panel.cor,main="Sequencing Statistics, datMeta, and expr Comparison |Spearman's rho|")
# pairs(cbind(pairsdat,topPC.expr),col=as.factor(datMeta$Dx),pch=19,upper.panel = panel.cor,main="Sequencing Statistics, datMeta, and expr Comparison |Spearman's rho|")
# pairs(cbind(topPC.expr,topPC.datQC),col=as.factor(datMeta$Dx),pch=19,upper.panel = panel.cor,main="PC expr, datQC Comparison |Spearman's rho|")
# dev.off()

colnames(topPC.expr) <- gsub("\n"," ", colnames(topPC.expr))


lm_vec <- function(matrix1, matrix2){
  lapply(as.data.frame(matrix1), function(v1){
    lapply(as.data.frame(matrix2), function(v2){
      summary(lm(v1~v2))$"adj.r.squared"
    }) %>% do.call(c,.)
  }) %>% do.call(rbind, .)
}

qc_pc_rsqr <- lm_vec(topPC.expr,pairsdat)
expr_pc_rsqr <- lm_vec(topPC.datQC,pairsdat)
qc_expr_cor <- cor(x = topPC.expr, y = topPC.datQC, method = "s") 

pdf(file=paste0(outputFolder,outputfile_counter,".1_SeqStats_Comparison_corrplot.pdf"), width = 10)
corrplot(qc_pc_rsqr,method = "number", tl.col	= "black", mar = c(1,1,1,1), cl.align.text = "l",number.cex	= 0.75)
corrplot(expr_pc_rsqr,method = "number", tl.col	= "black", mar = c(1,1,1,1), cl.align.text = "l",number.cex	= 0.75)
corrplot(qc_expr_cor,method = "number", tl.col	= "black", mar = c(1,1,1,1), cl.align.text = "l",number.cex	= 0.75)
# seq_qc_plot1 +seq_qc_plot2
dev.off()

outputfile_counter = outputfile_counter + 1
stopifnot(all(rownames(topPC.datQC) == datMeta$SampleID))
datMeta <- cbind(datMeta,topPC.datQC)

# (9) Covariate Plots - want to look at all datMeta -----------------------------------------------------------
# it is important to graph character covariates differently from numeric covariates - alter the code below to suit your needs


pdf(file = paste0(outputFolder,outputfile_counter,"_datMeta_covariates_dx_tot.pdf"),height = 11,width = 7)
par(mfrow = c(3,2))
par(mar = c(5,2,3,2))

for (i in grep("SeqPC", colnames(datMeta), value = T)) {
  A = anova(lm(datMeta[,i] ~ as.factor(datMeta$Dx)))
  p = A$"Pr(>F)"[1]
  plot(as.numeric(datMeta[,i]) ~ as.factor(datMeta$Dx), col = rainbow(length(unique(datMeta$Dx))),
       main = paste(i," p=", signif(p,2)), ylab = "", xlab = "", las = 3)
}

for (i in grep("expr PC", colnames(topPC.expr), value = T)) {
  A = anova(lm(topPC.expr[,i] ~ as.factor(datMeta$Dx)))
  p = A$"Pr(>F)"[1]
  plot(as.numeric(topPC.expr[,i]) ~ as.factor(datMeta$Dx), col = rainbow(length(unique(datMeta$Dx))),
       main = paste(i," p=", signif(p,2)), ylab = "", xlab = "", las = 3)
}

# #Plot covarience with QC data
# for (i in c(1:dim(datQC)[2])) {
#   A = anova(lm(as.numeric(datQC[,i]) ~ datMeta$Dx))
#   p = A$"Pr(>F)"[1];   
#   plot(as.numeric(as.character(datQC[,i])) ~ as.factor(datMeta$Dx), col = rainbow(length(unique(datMeta$Dx))),
#        main=paste(colnames(datQC)[i]," \np=", signif(p, 2)), ylab = "", xlab = "", las = 3)
# }
dev.off()

covars=c("Dx","Differentiation.day","Sex")
QCheatMap(datQC,datMeta ,covars, outputFolder,paste0(outputfile_counter,".1_QCheatmap"))

outputfile_counter=outputfile_counter+1

# (11) Regress technical and biological covariates  ----------------------------------------------
covars_to_protect <- c("Dx","Day.grouped")
covars_to_remove <- c("Sex", "Source.Tissue","ethnicityPC1","ethnicityPC1", paste0("SeqPC", 1:15))
form2 <- paste("~", paste(c(covars_to_protect,covars_to_remove),collapse = "+"))

X <-   model.matrix(as.formula(form2), data = datMeta)
Y <-  datExpr.htg
beta <-  (solve(t(X) %*% X) %*% t(X)) %*% t(Y)
regress_indx <- grep("Intercept|Dx|Day.grouped",colnames(X), invert = T)
to_regress <-  (as.matrix(X[,regress_indx]) %*% (as.matrix(beta[regress_indx,])))  #SeqPCs + covars
datExpr_reg <-  datExpr.htg - t(to_regress)

pdf(file=paste0(outputFolder,outputfile_counter,"_ProcessedStatistics_reg_biol.pdf"))
viewData(datMeta, datExpr_reg)
dev.off()

outputfile_counter=outputfile_counter+1

# (12) Batch correction  -------------------------------------------------------------------------------------------
if(FALSE){
  #Don't run as it makes results more batchy (batch 4 doesn't have many controls)
  covars_model3 <- model.matrix(~Dx+Differentiation.day, data = datMeta)
  
  datExpr_reg_batch <- sva::ComBat(datExpr_reg, batch = datMeta$batch, mod = covars_model3)
  
  pdf(file=paste0(outputFolder,outputfile_counter,"_ProcessedStatistics_reg_combat.pdf"))
  viewData(datMeta, datExpr_reg_batch, cvrs = c("Dx","batch","Differentiation.day","Sex","Source.Tissue"))
  dev.off()
  
  
  outputfile_counter=outputfile_counter+1
} else {
  datExpr_reg_batch <- datExpr_reg #combat looked like it was pulling out batch4 
}




# (13) Test PCA after corrections-------------------------------------------------------
pc_n <- 20
thisdat.expr <- t(scale(t(datExpr_reg_batch), scale = F)) 
## Centers the mean of all genes - this means the PCA gives us the eigenvectors of the geneXgene covariance matrix, allowing us to assess the proportion of variance each component contributes to the data
PC.expr <- prcomp(thisdat.expr)
topPC.expr <- PC.expr$rotation[,1:pc_n]
varexp <- (PC.expr$sdev)^2 / sum(PC.expr$sdev^2) 
topvar <- varexp[1:pc_n]
colnames(topPC.expr) <- paste("expr ",colnames(topPC.expr)," (",signif(100*topvar[1:pc_n],2),"%)",sep="")

pairsdat <- datMeta[,c("Dx","Day.grouped","Sex", "batch", "Source.Tissue")]
pairsdat[,c(1:4)] <- lapply(pairsdat[,c(1:4)],as.factor)

qc_pc_rsqr <- lm_vec(topPC.expr,pairsdat) %>% t()

pdf(file=paste0(outputFolder,outputfile_counter,"_exprPCs_corrplot.pdf"), width = 10)
corrplot(qc_pc_rsqr,method="ellipse", tl.col	= "black", mar = c(1,1,1,1), cl.align.text = "l")
dev.off()

meta_pc_data <- cbind(datMeta,topPC.expr) %>% 
  gather("PC", "value", starts_with("expr ")) %>% 
  mutate(PC = factor(PC, levels = str_sort(unique(PC), numeric = T)))


pdf(file = paste0(outputFolder,outputfile_counter,".1_exprPCs_boxplot.pdf"), width = 10, height = 5)
for (pcs in split(1:pc_n, ceiling(seq_along(1:pc_n)/(pc_n/2)))) { # split pcs into 2 groups
  p1 <- meta_pc_data %>% 
    filter(PC %in% levels(PC)[pcs]) %>% 
    ggplot(aes(x = Dx, y = value , fill = Dx)) +
    facet_wrap(~PC, nrow = 2) +
    geom_boxplot() + 
    theme_classic() +
    theme(axis.text.x = element_text(angle = -90, hjust  = 0, vjust = 0.5))
  print(p1)
}
dev.off()


pdf(file = paste0(outputFolder,outputfile_counter,".2_exprPCs_Dx_by_sex.pdf"), width = 6.5, height = 5)
pcPlots <- cbind(datMeta,as.data.frame(topPC.expr))
names(pcPlots) <- c(names(datMeta),paste0("PC",1:ncol(topPC.expr)))
pcPlots %>% 
  ggplot(aes(x = PC1, y = PC2,color = Sex)) +
  geom_point(size = 1.5) +
  theme_bw() + 
  theme(legend.position = c(0.75,0.05),
        legend.direction = "horizontal",
        legend.key.size = unit(0.5, "cm")) +
  facet_wrap(~Dx)
dev.off()

pdf(file = paste0(outputFolder,outputfile_counter,".3_exprPCs_Dx_by_day.pdf"), width = 6.5, height = 5)
pcPlots <- cbind(datMeta,as.data.frame(topPC.expr))
names(pcPlots) <- c(names(datMeta),paste0("PC",1:ncol(topPC.expr)))
pcPlots %>% 
  ggplot(aes(x = PC1, y = PC2,color = Differentiation.day)) +
  geom_point(size = 1.5) +
  theme_bw() + 
  theme(legend.position = c(0.75,0.05),
        legend.direction = "horizontal",
        legend.key.size = unit(0.5, "cm")) +
  facet_wrap(~Dx)
dev.off()


#  (14) save and wrap up  -----------------------------------------------------------------------------------------

save(datMeta,
     datQC, 
     datExprRaw,  # raw counts
     datExpr.htg, # normalized counts
     # datExpr_reg_tech, # counts after regressing out seqPC per batch
     datExpr_reg, # counts after regressing out seqPC per batch and biological covars
     datExpr_reg_batch, # counts after regressing out seqPC per batch and biological covars and batch corrected
     file = paste0(outputFolder,"input_for_DE_all.rdata"))




#proc.time() - time0


