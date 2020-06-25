# set up R and some variables ----------------------------
setwd("C://Users/natha/Documents/Geschwind Rotation 2020/transcAnalysis")
options(stringsAsFactors = FALSE)
# library(cowplot)
library(data.table)
library(WGCNA) #had to install via BiocManager
library(plyr)
# library(cqn) #had to install via BiocManager
library(tidyverse)
library(stringr)
# library(gridExtra)
# library(corrplot)
# library(wesanderson)
# library(patchwork)
# library(sva) #missing OG code, had to install via BiocManager
library(dendextend)

outputfile_counter = 1
outputFolder = paste(getwd(), '/cluster_OUT_20200623/', sep = "") # SET AN OUTPUT PATH HERE
dir.create(outputFolder, showWarnings = F, recursive = T)

# (1) Load Raw Expression, Meta and Normalized data  and format them ----------------------------------------
#LOAD DATA HERE
load('input_for_DE_transc.Rdata')

datTrans <- datTrans_reg_batch

# (2) divide data set by time points ----------------------------

samples_accounted <- 0
for(day in unique(datMeta$Day.grouped)){
  datMetaDxd <- datMeta %>% 
    filter(Day.grouped == day)
  datTransDxd <- datTrans[,datMetaDxd$SampleID,drop=F]
  assign(paste0('datTrans',day),datTransDxd)
  
  samples_accounted <- samples_accounted + dim(datTransDxd)[2]
}

stopifnot(samples_accounted == dim(datTrans)[2])
datDays <- ls(pattern = "datTrans[[:digit:]]{3}")


# (3) Calculate median expression for each gene per time point 

datTransDayMed <- matrix(data = NA, nrow = dim(datTrans)[1], ncol = 0)

for (day in datDays) {
  datTransDay <- get(day)
  data.frame(assign(paste0(day, 'median'), apply(datTransDay, 1, median))) 
  datTransDayMed <- cbind(datTransDayMed, get(paste0(day, 'median')))
}
colnames(datTransDayMed) <- gsub("datTrans","", datDays)


# (4) Clustering the days

daydis <- dist(t(datTransDayMed))
dayhc <- hclust(daydis, method = 'ward.D2')

pdf(file = paste0(outputFolder,outputfile_counter,"_median_day_hclust.pdf"), width = 6, height = 6)
plot(dayhc, hang = 0.1,
     main = 'Differentiation Day Median Expression', sub = NULL,
     xlab = "Day", ylab = "Height")

outputfile_counter <- outputfile_counter +1
dev.off()

# (5) Take top 3000 most variable genes ----------------
n_var <- 3000
datTransvar <- apply(datTrans, 1, var)
topTransvar <- order(datTransvar, decreasing = TRUE)[1:n_var]
datVar <- datTrans[topTransvar, ]


vardis <- dist(t(datVar))
varhc <- hclust(vardis, method = 'ward.D2')

# Clustering based on differentiation day 
id_025 <- row.names(filter(datMeta, Day.grouped == "025"))
id_050 <- row.names(filter(datMeta, Day.grouped == "050"))
id_075 <- row.names(filter(datMeta, Day.grouped == "075"))
id_100 <- row.names(filter(datMeta, Day.grouped == "100"))

cluster_cohorts <- c("id_025", "id_050", "id_075", "id_100")
n_cohorts <- length(cluster_cohorts)
pdf(file = paste0(outputFolder,outputfile_counter,"_topvar_day_hclust.pdf"), width = 8, height = 11) 
par(mfrow = c(ceiling(n_cohorts/2),2))

for (coh in cluster_cohorts){
  as.dendrogram(varhc) %>%
    set("by_labels_branches_col", value = get(coh), type = "any") %>%
    set('labels_cex', 0.3) %>%
    plot(main = '3000 Most Variable Isoforms', sub = NULL,
      xlab = coh, ylab = "Height")
}
# 
# as.dendrogram(varhc) %>%
#   set("by_labels_branches_col", value = id_050, type = "any") %>%
#   set('labels_cex', 0.3) %>%
#   plot( main = '3000 Most Variable Isoforms', sub = NULL,
#         xlab = "Day 50 Colored", ylab = "Height")
# as.dendrogram(varhc) %>%
#   set("by_labels_branches_col", value = id_075, type = "any") %>%
#   set('labels_cex', 0.3) %>%
#   plot(main = '3000 Most Variable Isoforms', sub = NULL,
#         xlab = "Day 75 Colored", ylab = "Height")
# as.dendrogram(varhc) %>%
#   set("by_labels_branches_col", value = id_100, type = "any") %>%
#   set('labels_cex', 0.3) %>%
#   plot(main = '3000 Most Variable Isoforms', sub = NULL,
#         xlab = "Day 100 Colored", ylab = "Height")

outputfile_counter <- outputfile_counter +1
dev.off()

# Clustering based on Dx

for (dx in unique(datMeta$Dx)){
  assign(paste0("id_", dx), row.names(filter(datMeta, Dx == dx)))
}

cluster_cohorts <- paste0("id_", unique(datMeta$Dx))
n_cohorts <- length(cluster_cohorts)
pdf(file = paste0(outputFolder,outputfile_counter,"_topvar_Dx_hclust.pdf"), width = 8, height = 11) 
par(mfrow = c(ceiling(n_cohorts/2),2))

for (coh in cluster_cohorts){
  as.dendrogram(varhc) %>%
    set("by_labels_branches_col", value = get(coh), type = "any") %>%
    set('labels_cex', 0.3) %>%
    plot(main = '3000 Most Variable Isoforms', sub = NULL,
         xlab = coh, ylab = "Height")
}
outputfile_counter <- outputfile_counter +1
dev.off()

datDend <- select(datMeta, Day.grouped, Dx)
colDend <- labels2colors(datDend)

plotDendroAndColors(
  dendro = varhc,
  colors = colDend,
  #groupLabels = datDend,
  rowText = datDend,
  dendroLabels = FALSE
)

# (#) find top 3000 var for each time point

n_var <- 3000

for (dday in unique(datMeta$Day.grouped)) {
  datMetaDxd <- datMeta %>% 
    filter(Day.grouped == dday)
  datTransDxd <- datTrans[,datMetaDxd$SampleID,drop=F]
  datTransDxdvar <- apply(datTransDxd, 1, var)
  topTransDxdvar <- order(datTransDxdvar, decreasing = TRUE)[1:n_var]
  data.frame(assign(paste0('topVar', dday), datTransDxd[topTransDxdvar,]))
  dxd_batch = as.factor(datMetaDxd$batch)
  dxd_dx = datMetaDxd$Dx
  assign(paste0('labelVar',dday), data.frame("batch" = dxd_batch, "dx" = dxd_dx, row.names = rownames(datMetaDxd))) 
}

for (dday in unique(datMeta$Day.grouped)) {
  datVarDx <- get(paste0('topVar', dday))
  labVarDx <- get(paste0('labelVar', dday))
  varDistDx <- dist(t(datVarDx))
  varHcDx <- hclust(varDistDx, method = 'ward.D2')
  
  pdf(file = paste0(outputFolder,outputfile_counter, dday,"_topvar__hclust.pdf"), width = 8, height = 11) 
  plotDendroAndColors(
    dendro = varHcDx,
    colors = labels2colors(labVarDx),
    #groupLabels = datDend,
    rowText = labVarDx,
    dendroLabels = FALSE
  )
  dev.off()
  ouputfile_counter <- outputfile_counter + 1
}

