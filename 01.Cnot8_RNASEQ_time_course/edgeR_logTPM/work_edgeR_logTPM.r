##setwd("/home/mjwang/pwd/Quan-Cnot8-RNASEQ_all_mm9/04.analysis_featureCnt/edgeR_DESeq2_logTPM")
##standard procedure to do DGE test with TPM, not raw counts
##output correlation heatmap
##output vocanoplot and simpleDGE scatterplot

library(ComplexHeatmap)
library(viridis)
library(pheatmap)
library(dendsort)
library(dendextend)
library(edgeR)


##1, prepare data
TPM <- read.table("TPM.data",header = TRUE, row.names = 1) #15113    20  


#data <- read.delim("out.ftcount.format.mat", header = T,stringsAsFactors = FALSE) 
#read.delim() is almost the same as read.table(), except the field separator is tab by default.
TPM$gene <- rownames(TPM) #15113    21
TPM <- TPM[,c(21,1:20)]
x <- TPM #edgeR need raw count as input, do not feed log2 since edgeR will transform

test <- c("gene","WT.ESC_rep1", "WT.ESC_rep2", "KO.ESC_rep1", "KO.ESC_rep2", 
                 "WT.6h_rep1", "WT.6h_rep2", "KO.6h_rep1", "KO.6h_rep2",
                 "WT.12h_rep1", "WT.12h_rep2", "KO.12h_rep1","KO.12h_rep2",
                 "WT.24h_rep1", "WT.24h_rep2", "KO.24h_rep1","KO.24h_rep2",
                 "WT.48h_rep1", "WT.48h_rep2", "KO.48h_rep1","KO.48h_rep2")
all.equal(test,colnames(x)) #TRUE


par(mar=c(8,5,1,1))
boxplot(x[,-1],outline=FALSE,las=2,ylab='TPM')




experiment <- read.table("experiment.design",header = T,row.names = 1)
group <- factor(paste(experiment$Treat,experiment$Time,sep="."),
                levels = unique(paste(experiment$Treat,experiment$Time,sep=".")))


y <- DGEList(counts=x[,-1],group=group,genes = x$gene)

design <- model.matrix(~0+group)
colnames(design) <- levels(group)

keep <- filterByExpr(y,design)
table(keep) #10569 kept, total 15113, ###16311 keeped in CPM
y <- y[keep,keep.lib.sizes=FALSE] #filtered

AveLogCPM <- aveLogCPM(y)

hist(AveLogCPM)

#calculate norm factor and disperse factor
y <- calcNormFactors(y) #TMM ass default?
y <- estimateDisp(y,design)
#examine overall dispersion  
#the asymptotic value for the BCV tends to be in range from 0.05 to 0.2 
#for genetically identical mice or cell lines, whereas somewhat larger 
#values (>0.3) are observed for human subjects.
plotBCV(y)
abline(v=0,lty=2,lwd=2);  abline(h=0.3,lty=2,lwd=2); text(12,0.4,"y=0.3")


# ########output edgeR normalized count mat
# count <- y$counts
# #row.names(count) <- as.character(y$genes[,1])
# row.names(count) <- x$gene
# f <- y$samples$norm.factors 
# 
# head(count)
# count.norm <- t(t(count) * f) 
# head(count.norm)
# 
# par(mar=c(8,4,1,1))
# boxplot(count ,las = 3, outline = F, col = rep(c("navy","navy","brown","brown"),5), medcol='white',medlwd=1.2,
#         whisklwd=1.5,outpch = 1, outcol = "black",ylab="")
# title(ylab = "Raw counts Filtered",line = 2,cex.lab=1.5,family='Calibri Light') 
# 
# boxplot(count.norm,las = 3,  outline = F,col = rep(c("navy","navy","brown","brown"),5), medcol='white',medlwd=1.2,
#         whisklwd=1.5,outpch = 1, outcol = "black",ylab="")
# title(ylab = "edgeR normalized counts Filtered",line = 2,cex.lab=1.5,family='Calibri Light') 
# 
# 
# write.table(x = count, row.names = TRUE, quote = FALSE, 
#             sep = '\t', col.names = NA, file = "out.count.edgeR_filter.txt")
# write.table(x = round(count.norm,digits=3), row.names = TRUE, quote = FALSE, 
#             sep = '\t', col.names = NA, file = "out.count.edgeR_filter_norm.txt")
# #####



###calculate CPM or logCPM
logCPM <- cpm(y,prior.count = 1, log=TRUE) #actually log2TPM, different with ../log2TPM.data
#row.names(logCPM) <- x$gene #use all genes without filtering
#row.names(logCPM) <- as.character(y$genes[,1]) #use filtered genes
nrow(logCPM) #10569
nrow(x) #15113
write.table(x = logCPM, row.names = TRUE, quote = FALSE, sep = '\t', col.names = TRUE, file = "out.logTPM.edgeR.txt")
##log2TPM <- read.table("../log2TPM.data",header = TRUE, row.names = 1)
##cor(logCPM["Nuak2",],unlist(log2TPM["Nuak2",]) ) ##0.996, but value not equal

#######output  filtered averaged logCPM

logCPM.WT.ESC.ave <- rowMeans(logCPM[,1:2]) 
logCPM.KO.ESC.ave <- rowMeans(logCPM[,3:4])  

logCPM.WT.6h.ave <- rowMeans(logCPM[,5:6]) 
logCPM.KO.6h.ave <- rowMeans(logCPM[,7:8])

logCPM.WT.12h.ave <- rowMeans(logCPM[,9:10]) 
logCPM.KO.12h.ave <- rowMeans(logCPM[,11:12]) 

logCPM.WT.24h.ave <- rowMeans(logCPM[,13:14]) 
logCPM.KO.24h.ave <- rowMeans(logCPM[,15:16]) 

logCPM.WT.48h.ave <- rowMeans(logCPM[,17:18]) 
logCPM.KO.48h.ave <- rowMeans(logCPM[,19:20]) 

logCPM.ave <- data.frame(
                                       WT.ESC=logCPM.WT.ESC.ave,
                                       WT.6h=logCPM.WT.6h.ave,
                                       WT.12h=logCPM.WT.12h.ave,
                                       WT.24h=logCPM.WT.24h.ave,
                                       WT.48h=logCPM.WT.48h.ave,
                                       KO.ESC=logCPM.KO.ESC.ave,
                                       KO.6h=logCPM.KO.6h.ave,
                                       KO.12h=logCPM.KO.12h.ave,
                                       KO.24h=logCPM.KO.24h.ave,
                                       KO.48h=logCPM.KO.48h.ave
)


write.table(x = round(logCPM.ave,digits=3), row.names = TRUE, quote = FALSE, 
            sep = '\t', col.names = NA, file = "out.logTPM.edgeR.ave.txt")

#ave.log2TPM <- read.table("../ave.log2TPM.combine.all.data",header = TRUE, row.names = 1)
#cor(unlist(logCPM.ave["Nuak2",]),unlist(ave.log2TPM["Nuak2",c(1,3,5,7,9,2,4,6,8,10)]) ) ##0.997864 but value not equal

####



####2 sample clustering, heatmap correlation and PCA
#### sample clustering, use all genes -> plotMDS -> get xyz
points <- c(0,2,5,6,15,16,17,18,23,24)
colors <- c("blue", "darkgreen", "red","orange","brown","purple","steelblue","black","darkgrey","darkmagenta")


##total 10569, use 1000 genes 
pdf("out_PCA_MDS/MDS.log2TPM.useallgenes.pdf",width = 6.8125, height=6.3750,onefile = TRUE)
par(mar=c(4,6,3,1))
res.mds <- limma::plotMDS(logCPM, col=colors[group], pch=points[group],
                          gene.selection = "common",ndim = 3,cex=1,
                          top = 15000,main='Multiple Dimensional Scaling (MDS) with all genes',
                          xlab ="Principal Component 1 (46.9% )",ylab = "Principal Component 2 (17.1%)" )
#res.mds <- limma::plotMDS(logCPM, col=colors[group], pch=points[group],gene.selection = "common",ndim = 3,cex=1.1,top = 100)
legend("bottomright", legend=levels(group), pch=points, col=colors,ncol=2,text.width = 0.2,box.lty=0,bg=NA)
dev.off()
#MDS analysis do not report explained variance percentage,which was calculated from the PCA analysis

mds.xy <- data.frame(x=res.mds$cmdscale.out[,1],y=res.mds$cmdscale.out[,2])
#plot(mds.xy,pch=20,col=rep(c("navy","navy","firebrick","firebrick"),5),cex=1.2)

mds.xyz <- data.frame(x=res.mds$cmdscale.out[,1],y=res.mds$cmdscale.out[,2],z=res.mds$cmdscale.out[,3])
write.table(t(mds.xyz ),row.names = FALSE, col.names = FALSE, sep = ",",quote = FALSE, 
            file = "xyz.data.logTPM.MDS")


##get ave to plot ave for reps
# mds.xy.ave <- data.frame()
# for(i in seq(1,19,by=2)){
#   mds.xy.ave <- rbind(mds.xy.ave ,(mds.xy[i,] + mds.xy[i+1,])/2)
# }
# #points(mds.xy)
# #plot(mds.xy.ave,pch=20,col=rep(c("navy","firebrick"),5),cex=2)

groups <- gsub("_rep[12]","",row.names(mds.xy))
groups <- factor(groups,levels = unique(groups))
mds.xy$groups <- groups

mds.xy.ave = mds.xy %>% dplyr::group_by(groups) %>% dplyr::summarise( meanx = mean(x), meany = mean(y)  )
mds.xy.ave = as.data.frame(mds.xy.ave)

pdf("out_PCA_MDS/ave.MDS.log2TPM.useallgene.pdf",width = 6.8125, height=6.3750,onefile = TRUE)
par(mar=c(4,6,3,1))
x.coord <- mds.xy.ave[,2]
y.coord <- mds.xy.ave[,3]
plot(x.coord,y.coord,type='p',col=colors[mds.xy.ave$groups], pch=points[mds.xy.ave$groups],
     xlab ="Principal Component 1 (46.9% )",ylab = "Principal Component 2 (17.1%)",
     xpd=TRUE,
     main = "Multiple Dimensional Scaling (MDS) with all genes"  )
legend("bottomright", legend=levels(xyz.ave$groups), pch=points, col=colors,ncol=2,text.width = 0.2,box.lty=0,bg=NA)
dev.off()

write.table(t(mds.xy.ave[,-1]) ,row.names = FALSE, col.names = FALSE, sep = ",",quote = FALSE, 
            file = "xy.data.logTPM.MDS.ave")


###customized MDS (multiple dimensional scaling, the same with limma::plotMDS(gene.selection = "common",top = 15000 ))
d <- dist(t(logCPM)) #use all genes
fit <- cmdscale(d,eig=TRUE,k=3)
plot(fit$points[,1], fit$points[,2], xlab="Coordinate 1", ylab="Coordinate 2", 
     main="Metric MDS",col=colors[group], pch=points[group])
legend("bottomright", legend=levels(group), pch=points, col=colors, ncol=2,text.width = 25)

###


################draw correlation mat
res.cor <- cor(logCPM)

# ##use ggcorrplot
# require('ggcorrplot')
# require(ggplot2)
# 
# require(RColorBrewer)
# color<- colorRampPalette(c("red", "white"))(10)
# #color <- brewer.pal(n=8, name="PuOr")
# #barplot(height = rep(2,100),width = rep(2,100), col = color,border = color )
# ggcorrplot(res.cor,outline.col = "white",lab = TRUE, hc.order = FALSE) +
#   #ggcorrplot(res.cor,outline.col = "white",lab = TRUE, hc.order = TRUE, type = "lower") +
#   scale_fill_gradient2(limit = c(0.6,1), low = "blue", high =  "red", mid = "white", midpoint = 0.8)
# 
# #corrplot::corrplot(res.cor,method = 'color',col = rev(color),cl.lim = c(0.8,1))
# 
# require(pheatmap)
# #my.breaks <- c( seq(0.5, 1, by=.1)) 
# #my.colors <- c(colorRampPalette(colors = c("white","red" ))(length(my.breaks)))
# ##my.colors <- c(colorRampPalette(colors = c("blue", "white"))(length(my.breaks)/2), 
# #               colorRampPalette(colors = c("white", "orange", "red", "purple"))(length(my.breaks)/2))
# res.ph <- pheatmap::pheatmap(res.cor,color = c(colorRampPalette(colors = c("white","red" ))
#                                      (20)),cluster_cols = TRUE,cluster_rows = TRUE)


##use correlation as dist
# res_dist <- as.dist( -( res.cor - 1 )/2 )
# res_clus <- hclust(res_dist, method = "complete")
# plot(as.dendrogram(res_clus))


####order the dendgram
dend = hclust(dist(res.cor,method = 'euclidean'),method = "complete")
plot(dend);labels(dend)
#dend = dendsort(dend,isReverse = FALSE,type="min")
#plot(dend);labels(dend)
dend = rotate(dend,c(1:8,20:17,13:16,9:12))
plot(dend);labels(dend)
#dend = dendsort(hclust(dist(res.cor,method = 'euclidean'),method = "complete"),isReverse = FALSE)
#Heatmap(mat, name = "mat", cluster_rows = dend, column_title = "reorder by dendsort")

##dend = dendsort(hclust(dist(res.cor,method = 'euclidean'),method = "complete"),isReverse = FALSE,type="min") #use this
#dend = dendsort(hclust(dist(res.cor,method = 'euclidean'),method = "complete"),isReverse = FALSE,type="average")
pdf("res.cor.log2TPM.pdf",width = 6.8125, height=6.3750,onefile = FALSE)
pheatmap(res.cor,fontsize_col = 16,fontsize_row = 12,cluster_rows = dend, kmeans_k = NA,
         cluster_cols = dend,main = "log2TPM correlation of samples (euclidean,complete)",scale = "none",
         show_rownames = TRUE,show_colnames = TRUE,
         clustering_method = "complete", #complete,ward.D2
         clustering_distance_cols  = "euclidean",#euclidean,correlation 
         clustering_distance_rows  = "euclidean",#euclidean,correlation
         
         #border_color = "black",
         color = colorRampPalette(c("navy","white","firebrick"))(20),
         #color = colorRampPalette(c("grey","white","firebrick"))(20),
         #color = c(colorRampPalette(colors = c("white","red" ))(20)),
         breaks = seq(0.8,1,by=0.2/20)
)
dev.off()

saveRDS(res.cor,'res.cor.rds')

# ##use complexHeatmap and try to sort dendrogram
# #7.739583 6.885417
# 
# hp = Heatmap(res.cor, name = "correlation", cluster_rows = TRUE, 
#              cluster_columns = TRUE, show_row_names = TRUE,show_column_names = TRUE,use_raster = FALSE,#will use raster if >2000 row or cols, however rstudio do not support raster
#              col = circlize::colorRamp2(seq(0.8,1,by=0.2/20), colorRampPalette(c("navy","white","firebrick"))(21) ),
#              rect_gp = gpar(col = "#A9A5A5", lwd = 0.8),
#              row_dend_reorder = TRUE,column_dend_reorder = TRUE,
#              column_title = "log2CPM hclust:complete euclidean",
#              clustering_method_columns = "complete", ##ward.D,complete
#              clustering_distance_columns  = function (m) dist(m,method="euclidean") #euclidean, manhattan,minkowski  #clustering_distance_columns  = "pearson",
#              #heatmap_legend_param = list(color_bar = "continuous")
#              #right_annotation = ha#,heatmap_width=unit(8, "cm")
# ) 
# #width = max(grobWidth(textGrob(labels))))
# draw(hp, heatmap_legend_side = "left",gap = unit(0.1, "cm")) #8.895833 8.843750
# 
# hp = Heatmap(res.cor.log2TPM, name = "correlation", cluster_rows = TRUE, 
#              cluster_columns = TRUE, show_row_names = TRUE,show_column_names = TRUE,use_raster = FALSE,#will use raster if >2000 row or cols, however rstudio do not support raster
#              col = circlize::colorRamp2(seq(0.6,1,by=0.4/20), colorRampPalette(c("navy","white","firebrick"))(21) ),
#              rect_gp = gpar(col = "#A9A5A5", lwd = 0.8),
#              row_dend_reorder = TRUE,column_dend_reorder = TRUE,
#              column_title = "log2TPM hclust:complete euclidean",
#              clustering_method_columns = "complete", ##ward.D,complete
#              clustering_distance_columns  = function (m) dist(m,method="euclidean") #euclidean, manhattan,minkowski  #clustering_distance_columns  = "pearson",
#              #heatmap_legend_param = list(color_bar = "continuous")
#              #right_annotation = ha#,heatmap_width=unit(8, "cm")
# ) 
# #width = max(grobWidth(textGrob(labels))))
# draw(hp, heatmap_legend_side = "left",gap = unit(0.1, "cm")) #8.895833 8.843750
# 

# 
# ####test for manually creat hclust dendragram
#   a <- list()  # initialize empty object
# # define merging pattern: 
# #    negative numbers are leaves, 
# #    positive are merged clusters (defined by row number in $merge)
# a$merge <- matrix(c(-1, -2,
#                     -3, -4,
#                     1,  2), nc=2, byrow=TRUE ) 
# a$height <- c(1, 1.5, 3)    # define merge heights
# a$order <- c(1,2,3,4)              # order of leaves(trivial if hand-entered)
# a$labels <- c("A","B","C","D")   # labels of leaves
# class(a) <- "hclust"        # make it an hclust object
# plot(a)                     # look at the result   
# 
# #convert to a dendrogram object if needed
# ad <- as.dendrogram(a)
# ###########
# 
# 
# res.ph$tree_col$order
# res.ph$tree_row$order
# dd <- as.dendrogram(res.ph$tree_col)
# dd.reorder <- reorder(dd,order = 10:1)
# plot(dd.reorder)
# 
# plot(as.dendrogram(res.ph$tree_col))
# res.ph$tree_col$order <- c(1,  2,  3,  4,  7,  8,  5,  6, 11, 12,  9, 10, 17, 18, 15, 16, 13, 14, 19, 20)
# plot(as.dendrogram(res.ph$tree_col))
# 
# #col.order.new <- c(1,  2,  3,  4,  7,  8,  5,  6, 11, 12,  9, 10, 17, 18, 15, 16, 13, 14, 19, 20)
# #row.order.new <- c(1,  2,  3,  4,  7,  8,  5,  6, 11, 12,  9, 10,17, 18, 15, 16, 13, 14, 19, 20)
# 
# col.order.new <- c(17, 18, 15, 16, 13, 14, 19, 20,11, 12,  9, 10, 7,  8,  5,  6, 1,  2,  3,  4 )
# row.order.new <- c(17, 18, 15, 16, 13, 14, 19, 20,11, 12,  9, 10, 7,  8,  5,  6, 1,  2,  3,  4)
# 
# res.ph.reorder <- pheatmap::pheatmap(res.cor[row.order.new ,col.order.new], color = c(colorRampPalette(colors = c("white","red" ))
#                                      (20)), cluster_cols = FALSE,cluster_rows = FALSE )
# 
# 
# counts.norm <- t(t(y$counts) * y$samples$norm.factors)
# boxplot(y$counts,outline=FALSE,las=2  )
# boxplot(counts.norm,outline=FALSE,las=2  )



##################do PCA analysis
#res.pca <- prcomp(t(log2(x[,-1]+0.001)))
res.pca <- prcomp(t(logCPM)) #with reps, do not filter low genes? 24936
write.table(t(res.pca$x[,1:3]),row.names = FALSE, col.names = FALSE, sep = ",",quote = FALSE, 
            file = "xyz.data.logTPM")

all.equal(res.pca,read.table('xyz.data.logTPM',header = FALSE, sep = ',') )

saveRDS(res.pca,'res.pca.rds')


#######aggregate ave for each two rows####
xyz <- res.pca$x[,1:3]
#class(xyz) #mat
xyz <- as.data.frame(xyz)
groups <- gsub("_rep[12]","",row.names(xyz))
groups <- factor(groups,levels = unique(groups))
xyz$groups <- groups

library(plyr)
xyz.ave <- ddply(xyz,.(groups),summarise,PC1=mean(PC1),PC2=mean(PC2),PC3=mean(PC3))
write.table(round(t(xyz.ave[,-1]),3),row.names = FALSE, col.names = FALSE, sep = ",",quote = FALSE, 
            file = "xyz.data.logTPM.ave")
###

# WT.ESC.ave <- (logCPM[,1]+logCPM[,2])/2
# KO.ESC.ave <- (logCPM[,3]+logCPM[,4])/2
# 
# WT.6h.ave <- (logCPM[,5]+logCPM[,6])/2
# KO.6h.ave <- (logCPM[,7]+logCPM[,8])/2
# 
# WT.12h.ave <- (logCPM[,9]+logCPM[,10])/2
# KO.12h.ave <- (logCPM[,11]+logCPM[,12])/2
# 
# WT.24h.ave <- (logCPM[,13]+logCPM[,14])/2
# KO.24h.ave <- (logCPM[,15]+logCPM[,16])/2
# 
# WT.48h.ave <- (logCPM[,17]+logCPM[,18])/2
# KO.48h.ave <- (logCPM[,19]+count.norm[,20])/2
# 
# logCPM.ave <- data.frame(
#                               WT.ESC=WT.ESC.ave,
#                               KO.ESC=KO.ESC.ave,
#                               WT.6h=WT.6h.ave,
#                               KO.6h=KO.6h.ave,
#                               WT.12h=WT.12h.ave,
#                               KO.12h=KO.12h.ave,
#                               WT.24h=WT.24h.ave,
#                               KO.24h=KO.24h.ave,
#                               WT.48h=WT.48h.ave,
#                               KO.48h=KO.48h.ave
# )
# 
# res.pca.ave <- prcomp(t(logCPM.ave)) #without reps
# write.table(t(res.pca.ave$x[,1:3]),row.names = FALSE, col.names = FALSE, sep = ",",quote = FALSE,
#             file = "xyz.data.logCPM.ave")


##get importance of each PCs
PC_sd <- setNames(res.pca$sdev , paste0("PC",1:length(res.pca$sdev)))
PC_var_explain <- ( PC_sd^2 ) / sum(PC_sd^2) * 100
round(PC_var_explain,3)
#PC1    PC2    PC3    PC4    PC5    PC6    PC7    PC8    PC9   PC10   PC11   PC12   PC13   PC14   PC15   PC16   PC17   PC18   PC19   PC20 
#46.906 17.112  9.660  6.771  4.161  2.782  2.110  1.972  1.762  1.484  0.937  0.788  0.768  0.639  0.620  0.567  0.466  0.267  0.228  0.000 


#############calculate importance of each variable (gene) ############
#http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/118-principal-component-analysis-in-r-prcomp-vs-princomp/

#var.coord = loadings * the component standard deviations
#var.cos2 = var.coord^2
#var.contrib. The contribution of a variable to a given principal component is (in percentage) : (var.cos2 * 100) / (total cos2 of the component)

# Helper function 
#::::::::::::::::::::::::::::::::::::::::
var_coord_func <- function(loadings, comp.sdev){
  loadings*comp.sdev
}
# Compute Coordinates
#::::::::::::::::::::::::::::::::::::::::
loadings <- res.pca$rotation
sdev <- res.pca$sdev
var.coord <- t(apply(loadings, 1, var_coord_func, sdev)) 
head(var.coord[, 1:4])

res.pca$rotation #var.coord is not this

# Compute Cos2
#::::::::::::::::::::::::::::::::::::::::
var.cos2 <- var.coord^2
head(var.cos2[, 1:4])

# Compute variable (gene) contributions
#::::::::::::::::::::::::::::::::::::::::
comp.cos2 <- apply(var.cos2, 2, sum)
contrib <- function(var.cos2, comp.cos2){var.cos2*100/comp.cos2}
var.contrib <- t(apply(var.cos2,1, contrib, comp.cos2))
head(var.contrib[, 1:4])

saveRDS(var.contrib,'var.contrib.rds')

# select the top 500 gene that most contribute to PC1

sum(var.contrib[,1]) #PC1 total 100
sum(var.contrib[,2]) #PC2 100


var.contrib.order.pc1 <- var.contrib[order( var.contrib[,'PC1'] ,decreasing = TRUE ),]
genes.most.contri.pc1 <- rownames(var.contrib.order.pc1 )#[1:500]

write.table(var.contrib.order.pc1,row.names = TRUE, col.names = NA, sep = "\t",quote = FALSE, 
            file = "genes_contribution_to_PC/var.contrib.ordered_by.pc1.txt")
write.table(genes.most.contri.pc1,row.names = FALSE, col.names = FALSE, sep = "\t",quote = FALSE, 
            file = "genes_contribution_to_PC/genes.most.contri.orded_by.pc1.txt")


var.contrib.order.pc2 <- var.contrib[order( var.contrib[,'PC2'] ,decreasing = TRUE ),]
genes.most.contri.pc2 <- rownames(var.contrib.order.pc2 )#[1:500]

write.table(var.contrib.order.pc2,row.names = TRUE, col.names = NA, sep = "\t",quote = FALSE, 
            file = "genes_contribution_to_PC/var.contrib.ordered_by.pc2.txt")
write.table(genes.most.contri.pc2,row.names = FALSE, col.names = FALSE, sep = "\t",quote = FALSE, 
            file = "genes_contribution_to_PC/genes.most.contri.orded_by.pc2.txt")



##plot rank curve

plot(var.contrib.order.pc1[1:500,1:2], pch=19,cex=0.5) #PC1 + PC2
#plot(1:500,var.contrib.order.pc1[1:500,1],pch = 19,cex=0.8)


var.df <- data.frame(
                    importance_of_PC1=var.contrib.order.pc1[1:500,c('PC1')]
                    #importance_of_PC2=var.contrib.order.pc1[1:500,c('PC2')] 
    )

var.df <- data.frame(
                    importance_of_PC2=var.contrib.order.pc2[1:500,c('PC2')]
                    #importance_of_PC2=var.contrib.order.pc2[1:500,c('PC2')] 
    )

library('ggrepel')
options(repr.plot.width=7.5,repr.plot.height=9)
ggplot(data=var.df,aes(x = 1:nrow(var.df),
                       y = importance_of_PC1, 
                       #y = importance_of_PC2, 
                       label = rownames(var.df)
                       #fill = color, 
                       #segment.color = color,
                       #color = color
                      )
      ) + #will be slow if many rows
  geom_point(size = 1,shape=19) +
  #geom_label_repel(color='white') +
  geom_text_repel(point.padding = NA,box.padding = 0.1,size = 5,max.overlaps = 25) +
  #scale_fill_manual(values=c("#999999", "red")) + #must use two at the same time
#   scale_color_manual(values=c("#999999", "red"),
#                      name="Type",
#                      labels=c('other','known TF')
#                     ) +
  xlab("ranking of genes") +
  ylab('Importance to PC1') +
  #ylab('Importance to PC2') +
  theme_bw()+
  theme(axis.text = element_text(size = 2,color='black'),
        axis.title = element_text(size = 25,color='black'),
        legend.title = element_text(size = 25,color='black'),
        legend.text = element_text(size=25,color='black')
       )


ggsave('genes_contribution_to_PC/ranking_of_genes_by_importance_to_PC1.pdf',width = 7.5,height = 9, useDingbats=FALSE)
#ggsave('genes_contribution_to_PC/ranking_of_genes_by_importance_to_PC2.pdf',width = 7.5,height = 9, useDingbats=FALSE)

#########




################plot the PC explained barplot###############
#http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/118-principal-component-analysis-in-r-prcomp-vs-princomp/
library(factoextra)
fviz_eig(res.pca)

##plot for graph of individual (feature gene), individuals with a similar profile are grouped together
fviz_pca_ind(res.pca,         
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping)
)
fviz_pca_var(res.pca)
###########


######customized get gene with high loading of PC1#####
loadings <- res.pca$rotation
#10569 gene x 20 PCs

#pick PC1 gene order by loading value
loadings.order.pc1 <- loadings[order(loadings[,'PC1'],decreasing = TRUE),]
genes.pc1.top100 <- rownames(head(loadings.order.pc1,n = 100) )
#####






##plot 2D, PC1 PC2
points <- c(0,2,5,6,15,16,17,18,23,24)
colors <- c("blue", "darkgreen", "red","orange","brown","purple","steelblue","black","darkgrey","darkmagenta")


pdf("out_PCA_MDS/PCA.log2TPM.pdf",width = 6.8125, height=6.3750,onefile = TRUE)
par(mar=c(4,6,3,1))
x.coord <- res.pca$x[,1]
y.coord <- res.pca$x[,2]*-1 ##flip y axis
plot(x.coord,y.coord,type='p',col=colors[group], pch=points[group],
     xlab ="Principal Component 1 (46.9% )",ylab = "Principal Component 2 (17.1%)",
     xpd=TRUE,
     main = "Principle Component Analysis (PCA) with all 10569 genes"  )
legend("bottomright", legend=levels(group), pch=points, col=colors,ncol=2,text.width = 15,box.lty=0,bg=NA)
dev.off()

##use ave points####
pdf("out_PCA_MDS/ave.PCA.log2TPM.pdf",width = 6.8125, height=6.3750,onefile = TRUE)
par(mar=c(4,6,3,1))
x.coord <- xyz.ave[,2]
y.coord <- xyz.ave[,3]*-1 ##flip y axis
plot(x.coord,y.coord,type='p',col=colors[xyz.ave$groups], pch=points[xyz.ave$groups],
     xlab ="Principal Component 1 (46.9% )",ylab = "Principal Component 2 (17.1%)",
     xpd=TRUE,
     main = "Principle Component Analysis (PCA) with all 10569 genes"  )
legend("bottomright", legend=levels(xyz.ave$groups), pch=points, col=colors,ncol=2,text.width = 15,box.lty=0,bg=NA)
dev.off()

# ##use ggscatter, can automaticly add text labels 
# # Plot and color by groups
# ggpubr::ggscatter(xyz[,1:2], x = "PC1", y = "PC2", 
#           label = rownames(xyz),
#           #color = rep(c(1,2,5,6,15,16,17,18,23,24),each=2),
#           palette = "jco",
#           size = 2, 
#           #ellipse = TRUE,
#           #ellipse.type = "convex",
#           repel = TRUE)



# ##with PC2, PC3
# x <- xyz.ave$PC2
# y <- xyz.ave$PC3
# plot(x.coord,y.coord,type='p',pch=c(0,2,5,6,15,16,17,18,23,24), col=c("steelblue","navy","darkgreen","red","orange",
#                                                           "brown","darkblue",'purple',"pink","black"),
#      xlab="PC1", ylab="PC2",xpd=FALSE,cex = 1.5  )
# text(x.coord,y.coord,col=c(1,2,5,6,15,16,17,18,23,24),
#      #text(res.pca$x[,1]+10,res.pca$x[,2]+200,col=rep(c(1,2,5,6,15,16,17,18,23,24),each=2),
#      labels = c("WT.ESC", "KO.ESC",
#                 "WT.6h", "KO.6h", 
#                 "WT.12h",  "KO.12h",
#                 "WT.24h",  "KO.24h",
#                 "WT.48h",  "KO.48h"
#      )  
# )


##3D scatterplot 
library(scatterplot3d)
# s3d<-scatterplot3d(res.pca$x[,1:3],type="h",color = rep(c("blue", "darkgreen", "red","orange",
#                                                           "brown","purple","steelblue","black",
#                                                           "darkgrey","darkmagenta"),each=2),
#                    pch = rep(c(0,2,5,6,15,16,17,18,23,24),each=2),angle = 150 )

s3d<-scatterplot3d(res.pca$x[,1:3],type="h",lty.hplot = rep(c(1, 1,2, 2),5),
                   color = rep(c("blue", "blue","red", "red"),5),
                   pch = rep(c(0,2,5,6,15),each=4),angle = 210,
                   xlab = "PC1 (46.9%)", ylab = "PC2 (17.1%)", zlab ="PC3 (9.7%)",
                   box = FALSE,main = "Principle Component Analysis", cex.main = 1.5)

legend(s3d$xyz.convert(-120, 15, 0), pch = rep(c(0,2,5,6,15),each=2),
       col= c("blue", "red"), 
       bg=NA, lty=0, lwd=2, yjust=0, 
       legend = levels(groups) , 
       cex = 0.8,box.lty = 0,xpd = TRUE) #box.lty legend frame box

# text(s3d$xyz.convert(res.pca$x[,1:3]),cex= 0.7,col=rep(c("blue", "darkgreen", "red","orange",
#                                                          "brown","purple","steelblue","black",
#                                                          "darkgrey","darkmagenta"),each=2),
#      labels = c("WT.ESC_rep1", "WT.ESC_rep2", "KO.ESC_rep1", "KO.ESC_rep2",
#                 "WT.6h_rep1", "WT.6h_rep2", "KO.6h_rep1", "KO.6h_rep2",
#                 "WT.12h_rep1", "WT.12h_rep2", "KO.12h_rep1","KO.12h_rep2",
#                 "WT.24h_rep1", "WT.24h_rep2", "KO.24h_rep1","KO.24h_rep2",
#                 "WT.48h_rep1", "WT.48h_rep2", "KO.48h_rep1","KO.48h_rep2"
#      )  )

saveRDS(y,"y.Apr6.2020.rds")

######3 DGE analysis for treat pairwise and time pairwise, output annotated scatterplot and volcano plot#########

my.contrasts.treatpair <- makeContrasts(
  kovswt.ESC = KO.ESC-WT.ESC,
  kovswt.6h = KO.6h-WT.6h,
  kovswt.12h = KO.12h-WT.12h,
  kovswt.24h = KO.24h-WT.24h,
  kovswt.48h = KO.48h-WT.48h,
  levels=design
  
)

my.contrasts.timepair <- makeContrasts(
  #ANOVA-like testing, https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
  #then test as a whole, anov <- glmQLFTest(fit, contrast=my.contrasts.timepair)
  wt.6vsESC = WT.6h-WT.ESC,
  wt.12vs6 = WT.12h-WT.6h,
  wt.24vs12 = WT.24h-WT.12h,
  wt.48vs24 = WT.48h-WT.24h,
  
  ko.6vsESC = KO.6h-KO.ESC,
  ko.12vs6 = KO.12h-KO.6h,
  ko.24vs12 = KO.24h-KO.12h,
  ko.48vs24 = KO.48h-KO.24h,
  
  levels=design
  
)


#fit the model quasi-likelihood F-tests:
fit.qlf <- glmQLFit(y,design)
#colnames(fit.qlf)
#head(fit$coefficients)
#plotQLDisp(fit.qlf)

#fit the model likelihood ratio 
fit.lrt <- glmFit(y,design)


##treat pairwise test
qlf.kovswt.ESC <- glmQLFTest(fit.qlf,contrast = my.contrasts.treatpair[,"kovswt.ESC"])
qlf.kovswt.6h <- glmQLFTest(fit.qlf,contrast = my.contrasts.treatpair[,"kovswt.6h"])
qlf.kovswt.12h <- glmQLFTest(fit.qlf,contrast = my.contrasts.treatpair[,"kovswt.12h"])
qlf.kovswt.24h <- glmQLFTest(fit.qlf,contrast = my.contrasts.treatpair[,"kovswt.24h"])
qlf.kovswt.48h <- glmQLFTest(fit.qlf,contrast = my.contrasts.treatpair[,"kovswt.48h"])

lrt.kovswt.ESC <- glmLRT(fit.lrt,contrast = my.contrasts.treatpair[,"kovswt.ESC"])
lrt.kovswt.6h <- glmLRT(fit.lrt,contrast = my.contrasts.treatpair[,"kovswt.6h"])
lrt.kovswt.12h <- glmLRT(fit.lrt,contrast = my.contrasts.treatpair[,"kovswt.12h"])
lrt.kovswt.24h <- glmLRT(fit.lrt,contrast = my.contrasts.treatpair[,"kovswt.24h"])
lrt.kovswt.48h <- glmLRT(fit.lrt,contrast = my.contrasts.treatpair[,"kovswt.48h"])

#quick look
summary(decideTests(qlf.kovswt.ESC))
summary(decideTests(qlf.kovswt.6h))
summary(decideTests(qlf.kovswt.12h))
summary(decideTests(qlf.kovswt.24h))
summary(decideTests(qlf.kovswt.48h))

# sum.qlf.kovswt.ESC = summary(decideTests(qlf.kovswt.ESC))
# sum.qlf.kovswt.6h=summary(decideTests(qlf.kovswt.6h))
# sum.qlf.kovswt.12h=summary(decideTests(qlf.kovswt.12h))
# sum.qlf.kovswt.24h=summary(decideTests(qlf.kovswt.24h))
# sum.qlf.kovswt.48h=summary(decideTests(qlf.kovswt.48h))

# DGE.qlf = data.frame(time = c("ESC","6h","12h","24h","48h"), 
#                       Down = c(sum.qlf.kovswt.ESC[1],sum.qlf.kovswt.6h[1],
#                                sum.qlf.kovswt.12h[1],sum.qlf.kovswt.24h[1],
#                                sum.qlf.kovswt.48h[1]), 
#                       Up =c(sum.qlf.kovswt.ESC[3],sum.qlf.kovswt.6h[3],
#                             sum.qlf.kovswt.12h[3],sum.qlf.kovswt.24h[3],
#                             sum.qlf.kovswt.48h[3])
#                       )

#plot(DGE.qlf$Down,type='b');points(DGE.qlf$Up,type='b')

summary(decideTests(lrt.kovswt.ESC))
summary(decideTests(lrt.kovswt.6h))
summary(decideTests(lrt.kovswt.12h))
summary(decideTests(lrt.kovswt.24h))
summary(decideTests(lrt.kovswt.48h))

# sum.lrt.kovswt.ESC=summary(decideTests(lrt.kovswt.ESC))
# sum.lrt.kovswt.6h=summary(decideTests(lrt.kovswt.6h))
# sum.lrt.kovswt.12h=summary(decideTests(lrt.kovswt.12h))
# sum.lrt.kovswt.24h=summary(decideTests(lrt.kovswt.24h))
# sum.lrt.kovswt.48h=summary(decideTests(lrt.kovswt.48h))
# 
# DGE.lrt = data.frame(time = c("ESC","6h","12h","24h","48h"),
#                       Down = c(sum.lrt.kovswt.ESC[1],sum.lrt.kovswt.6h[1],
#                                sum.lrt.kovswt.12h[1],sum.lrt.kovswt.24h[1],
#                                sum.lrt.kovswt.48h[1]),
#                       Up =c(sum.lrt.kovswt.ESC[3],sum.lrt.kovswt.6h[3],
#                             sum.lrt.kovswt.12h[3],sum.lrt.kovswt.24h[3],
#                             sum.lrt.kovswt.48h[3])
#                       )
# plot(DGE.lrt$Down,type='b');points(DGE.lrt$Up,type='b')


##time course, pairwise test
##time course, all test
anov <- glmQLFTest(fit.qlf, contrast=my.contrasts.timepair)


##############valcano plot DGE results

drawVolcanoDE <- function(data  = NULL,sample = NULL,select.up = NULL, select.down = NULL, 
                          FC.cutoff = NULL, padj.cutoff = NULL){
  #input DGE result table
  plot(data$logFC,-log10(data$padj),pch=20,cex = 0.5,xlim = c(-10,10),col='black',
       main=paste(sample,"\n(FDR < ",padj.cutoff,", FC > ",FC.cutoff,")",sep=""),xlab = expression(Log[2]~(Fold~Change)), ylab=expression(-Log[10]~(FDR)),
       cex.lab=1.5)
       #cex.lab=1.5,family='Calibri Light')
  points(data[select.up,"logFC"],-log10(data[select.up,"padj"]),col='red' ,pch=20,cex = 0.5)
  points(data[select.down,"logFC"],-log10(data[select.down,"padj"]),col='navy' ,pch=20,cex = 0.5)
  abline(v=c(-log2(FC.cutoff),log2(FC.cutoff)),lty=2,col='grey')
  abline(h=-log10(padj.cutoff),lty=2,col='grey')
  y_max <- max(-log10(data$padj))
  text(x = 11, y = y_max,labels = paste("Up =",sum(select.up),sep = ""), col='red', pos = 2) #align right end
  text(x = -11, y = y_max,labels = paste("Down =",sum(select.down),sep = ""), col='navy', pos = 4) #align left end
  text(x=log2(FC.cutoff)+1.2,y=7-0.1,pos=2,srt=90,labels = paste("FC = ",FC.cutoff,sep = " ")) 
  text(x=-log2(FC.cutoff),y=7,pos=2,srt=90,labels = paste("FC = ",-FC.cutoff,sep = " "))
  
  ##add some gene point annotation
  
  # expected.up.new <- read.table("expected.up.all",stringsAsFactors = F)$V1
  # expected.down.new <- read.table("expected.down.all",stringsAsFactors = F)$V1
  # 
  # expected.up.remove <- read.table("expected.up.remove",stringsAsFactors = F)$V1
  # expected.down.remove <- read.table("expected.down.remove",stringsAsFactors = F)$V1
  # 
  # expected.up.new <- expected.up.new[!(expected.up.new %in% expected.up.remove  )]
  # expected.down.new <- expected.down.new[!(expected.down.new %in% expected.down.remove  )]
  # 
  expected.up.new <- read.table("expected.up.new.new",stringsAsFactors = F)$V1
  expected.down.new <- read.table("expected.down.new.new",stringsAsFactors = F)$V1
  
  expected.up.new.data  <- data[match (expected.up.new,rownames(data)), ]
  expected.down.new.data  <- data[match (expected.down.new,rownames(data)), ]
  
  points(expected.up.new.data$logFC,-log10(expected.up.new.data$padj),col="yellow",pch=1,cex=0.65)
  points(expected.down.new.data$logFC,-log10(expected.down.new.data$padj),col="darkgreen",pch=1,cex=0.65)
  
  ##add text info
  ###text(x=14,y=15.8,pos=4,srt=45,labels = paste("FC =",format(2^threshold.FC,digits = 2),sep = " ")) 
  ###text(x=15.3,y=14.3,pos=4,srt=45,labels = paste("FC =",format(-1*2^threshold.FC,digits = 2),sep = " ")) 

  # for(i in 1:nrow(expected.down.new.data)){
  #   lines(c(expected.down.new.data[i,"logFC"],expected.down.new.data[i,"logFC"]+3),c(-log10(expected.down.new.data[i,"padj"]),-log10(expected.down.new.data[i,"padj"])+1),lwd=1,col='grey')
  #   text(x=expected.down.new.data[i,"logFC"]+3,y=-log10(expected.down.new.data[i,"padj"])+1,pos=4,labels=row.names(expected.down.new.data)[i],cex = 1.5 )
  # }
  # for(i in 1:nrow(expected.up.new.data)){
  #   lines(c(expected.up.new.data[i,"logFC"],expected.up.new.data[i,"logFC"]-3),c(-log10(expected.up.new.data[i,"padj"]),-log10(expected.up.new.data[i,"padj"])+3),lwd=1,col='grey')
  #   text(x=expected.up.new.data[i,"logFC"]-2.5,y=-log10(expected.up.new.data[i,"padj"])+3,pos=2,labels=row.names(expected.up.new.data)[i],cex = 1.5 )
  # }

  return(1)
}


#drawVolcanoDE(qlf.kovswt.ESC.table,"qlf.kovswt.ESC",select.up,select.down,FC.cutoff,padj.cutoff )

plotVolcano <- function( res.table = NULL, sample = NULL, FC.cutoff = NULL, 
                         padj.cutoff = NULL, threshold.exp = NULL, plot = FALSE ){
  ##prepare the DGE results table
  #           logFC   logCPM          F      PValue  ##log(2)FC and log(2)TPM, http://78.128.216.61:8081/manual/ngs-dea-edger-RNA.html
  #Nuak2   1.06769429 3.947513 14.1440260 0.016164348
  #### stopifnot( all.equal (row.names(res.table), row.names(data))   )

  res.table.new <- cbind(res.table[,c(-3)], 
  #res.table.new <- cbind(res.table[,c(-2,-3)], 
                         padj = p.adjust(res.table$PValue,method = "BH")  )
  #stopifnot(all.equal (row.names(res.table.new) , row.names(data)) )
  ####stopifnot(all.equal (row.names(res.table.new) , row.names(y$genes ) ) )
  
  #res.table.new <- as.data.frame(res.table.new )
  #qlf.kovswt.ESC.topTags <- topTags(qlf.kovswt.ESC,n = 30000)
  #all.equal(qlf.kovswt.ESC.table["222",],qlf.kovswt.ESC.topTags["222",]) #TRUE
  
  select.up <- res.table.new$padj < padj.cutoff  & res.table.new$logFC > log2(FC.cutoff) & res.table.new$logCPM > 0
  select.down <- res.table.new$padj < padj.cutoff  & res.table.new$logFC < -log2(FC.cutoff) & res.table.new$logCPM > 0
  if(plot == TRUE){
    drawVolcanoDE(res.table.new,sample,select.up,select.down,FC.cutoff,padj.cutoff )
  }
  ##output up and down regulated gene mat
  #logCPM must exist
  select.up.data <- round(cbind(res.table.new[select.up,],logCPM[select.up,]),digits = 3)
  select.down.data <- round(cbind(res.table.new[select.down,],logCPM[select.down,]),digits = 3)
  write.table(select.up.data,row.names = TRUE, col.names = NA, 
              sep = "\t",quote = FALSE, file = paste("out_DGE_qlf/select.up.edgeR.test.log2TPM.",sample,".FC.",FC.cutoff,".mat",sep=""))
  write.table(select.down.data,row.names = TRUE, col.names = NA, 
              sep = "\t",quote = FALSE, file = paste("out_DGE_qlf/select.down.edgeR.test.log2TPM.",sample,".FC.",FC.cutoff,".mat",sep=""))
  
  select <- data.frame(select.up=select.up,select.down=select.down)
  return(select)
}#draw pairwise 
#dev.off()



##simple scatter plot with annotation
#######do the edgeR test fold change visulization and output step
#height=4.166667,width = 4.177083
drawSimpleDE <- function(data  = NULL,sample = NULL,select.up = NULL, select.down = NULL, 
                         padj.cutoff = NULL, threshold.FC = NULL,threshold.exp = NULL,label =NULL){
  r <- format(cor(data[,1],data[,2],method = 'spearman'),digits = 3)
  #plot(data,pch=20,cex=0.65,xlim=c(-1,1000),ylim=c(-1,1000),col='lightgrey',
  #xlab=paste("WT",sample,"RPM",sep=" "),ylab=paste("KO",sample,"RPM",sep=" "))
  plot(data,pch=20,cex=0.65,xlim=c(-1,15),ylim=c(-1,15),
       col='lightgrey',#rgb( col2rgb('lightgrey')[1,], col2rgb('lightgrey')[2,], col2rgb('lightgrey')[3,],alpha = 250,maxColorValue = 255), 
       #pdf use font (ZapfDingbats and will not embed) if totally opaque, set alpha to 0.99 to use real symbol, to be friendly to AI 
       #see https://r.789695.n4.nabble.com/pdf-device-uses-fonts-to-represent-points-data-alteration-td836140.html
       xlab=paste("WT","log2TPM",sample,sep=" "),ylab=paste("KO","log2TPM",sample,sep=" "),
       main=paste(sample,"\n(FDR < ",padj.cutoff,", FC > ",threshold.FC,")",sep=""),
       cex.lab = 1.5)
  points(data[select.up,],pch=20,cex=0.65,col= 'red')#rgb( col2rgb('red')[1,], col2rgb('red')[2,], col2rgb('red')[3,],alpha = 250,maxColorValue = 255) )
  points(data[select.down,],pch=20,cex=0.65,col= 'blue') #rgb( col2rgb('blue')[1,], col2rgb('blue')[2,], col2rgb('blue')[3,],alpha = 250,maxColorValue = 255) )
  text(x = -1, y = 10,labels = paste("Up =",sum(select.up),sep = ""), col='red', pos = 4) #align left end
  text(x = 5.5, y = 1,labels = paste("Down =",sum(select.down),sep = ""), col='blue', pos = 4) 
  #text(x = -1, y = 17.5,labels = paste("R=",r,sep = ""), col='black', pos = 4)
  #text(x = -1, y = 19,labels = paste("N =",nrow(data),sep = ""), col='black', pos = 4)
  text(x = -1, y = 14,labels = paste("FC >=",threshold.FC,",FDR <",padj.cutoff,"\nlog2TPM > ",threshold.exp,sep = " "),col='black', pos = 4)
  #text(x = 0, y = 500,labels = paste("n=",sum(select.up),sep = ""), col='black')
  #text(x = 500, y = 1,labels = paste("n=",sum(select.down),sep = ""), col='black')
  #text(x = -1, y = 18,labels = paste("R=",r,sep = ""), col='black')
  #abline(a=0.5849,b=1,col='black',lty=2,lwd=1.2,untf = FALSE)
  #abline(a=-0.5849,b=1,col='black',lty=2,lwd=1.2,untf = FALSE)
  #abline(a=1,b=1,col='darkgrey',lty=2,lwd=1.2,untf = FALSE)
  #abline(a=-1,b=1,col='darkgrey',lty=2,lwd=1.2,untf = FALSE)
  abline(a=log2(threshold.FC),b=1,col='darkgrey',lty=2,lwd=1.2,untf = FALSE)
  abline(a=-1*log2(threshold.FC),b=1,col='darkgrey',lty=2,lwd=1.2,untf = FALSE)
  ##abline(a=0,b=1,col='darkgrey',lty=1,lwd=1.2)
  #abline(h=0,v=0)
  
  ##add some gene point annotation
  
  # expected.up.new <- read.table("expected.up.all",stringsAsFactors = F)$V1
  # expected.down.new <- read.table("expected.down.all",stringsAsFactors = F)$V1
  # 
  # expected.up.remove <- read.table("expected.up.remove",stringsAsFactors = F)$V1
  # expected.down.remove <- read.table("expected.down.remove",stringsAsFactors = F)$V1
  # 
  # expected.up.new <- expected.up.new[!(expected.up.new %in% expected.up.remove  )]
  # expected.down.new <- expected.down.new[!(expected.down.new %in% expected.down.remove  )]
  # 
  expected.up.new <- read.table("expected.up.new.new",stringsAsFactors = F)$V1
  expected.down.new <- read.table("expected.down.new.new",stringsAsFactors = F)$V1
  
  expected.up.new.data  <- data[match (expected.up.new,rownames(data)), ]
  expected.down.new.data  <- data[match (expected.down.new,rownames(data)), ]
  
  points(expected.up.new.data,col="yellow",pch=1,cex=0.65)
  points(expected.down.new.data,col="darkgreen",pch=1,cex=0.65)
  
  ##add text info
  text(x=12,y=13.5,pos=4,srt=45,labels = paste("FC =",threshold.FC,sep = " ")) 
  text(x=13.3,y=12.3,pos=4,srt=45,labels = paste("FC =",-1*threshold.FC,sep = " ")) 
  for(i in 1:nrow(expected.up.new.data)){
    lines(c(expected.up.new.data[i,1],expected.up.new.data[i,1]+3),c(expected.up.new.data[i,2],expected.up.new.data[i,2]+1),lwd=1,col='black')
    if(label == TRUE){
      text(x=expected.up.new.data[i,1]+3,y=expected.up.new.data[i,2]+1,pos=4,labels=row.names(expected.up.new.data)[i],cex = 0.8 )
    }
  }
  for(i in 1:nrow(expected.down.new.data)){
    lines(c(expected.down.new.data[i,1],expected.down.new.data[i,1]-1),c(expected.down.new.data[i,2],expected.down.new.data[i,2]+3),lwd=1,col='black')
    if(label == TRUE){
      text(x=expected.down.new.data[i,1]-0.5,y=expected.down.new.data[i,2]+3,pos=2,labels=row.names(expected.down.new.data)[i],cex = 0.8 )
    }  
  }

  return(1)
}
###########



#9.500000 7.979167

#threshold.exp <- 1 #log2TPM > 1 : TPM > 2
threshold.exp <- 0 #log2TPM > 0 : TPM > 1

#FC.cutoff <- 2
FC.cutoff <- 1.5
#FC.cutoff <- 1.2

padj.cutoff <- 0.05 #FDR
#padj.cutoff <- 0.1


logCPM.ave.bk = logCPM.ave # TPM -> cpm(y,prior.count = 1, log=TRUE) -> ave by reps #edgeR::cpm() will + 1
#in fact is logTMP.ave
#10569 x 10

ave.log2TPM <- read.table("../ave.log2TPM.combine.all.data",header = TRUE, row.names = 1) 
#this is the log2TPM used in detect_trend, use method: ave.log2.combine.all <- log2(ave.combine.all+0.001)
ave.log2TPM <- ave.log2TPM[,c(1,3,5,7,9,2,4,6,8,10)] #15113    10
#logCPM.ave.bk <-  logCPM.ave#acutally log2TPM, generated by edgeR
#logCPM.ave <- ave.log2TPM #use to test agian
##logCPM.ave <-  logCPM.ave.bk ##recover logCPM.ave

############correlation of logCPM.ave (edgeR::cpm, +1) and  ave.log2TPM (+0.001)##################
shared.geneid = intersect(row.names(logCPM.ave),row.names(ave.log2TPM)  ) #10569
#res.cor.compare = cor(cbind (logCPM.ave[shared.geneid,], ave.log2TPM[shared.geneid,] ) )
for (i in 1:10){
  cor.value = cor(logCPM.ave[shared.geneid,i],ave.log2TPM[shared.geneid,i])
  cat(colnames(logCPM.ave)[i]," vs ", colnames(ave.log2TPM)[i], ":",cor.value,"\n")
}
#WT.ESC  vs  wt.ave.ESC : 0.9697453 
#WT.6h  vs  wt.ave.6h : 0.9809008 
#WT.12h  vs  wt.ave.12h : 0.9903175 
#WT.24h  vs  wt.ave.24h : 0.9895322 
#WT.48h  vs  wt.ave.48h : 0.9779804 
#KO.ESC  vs  ko.ave.ESC : 0.9755488 
#KO.6h  vs  ko.ave.6h : 0.9861166 
#KO.12h  vs  ko.ave.12h : 0.9897031 
#KO.24h  vs  ko.ave.24h : 0.9974456 
#KO.48h  vs  ko.ave.48h : 0.9909722 

# ##good correlation
#############


# ##search for expected.up.new.new against logCPM.ave
# expected.up.new <- read.table("expected.up.new.new",stringsAsFactors = F)$V1 #24
# expected.down.new <- read.table("expected.down.new.new",stringsAsFactors = F)$V1 #10
# 
# expected.down.new[10]
# grep(pattern = expected.down.new[10],rownames(logCPM.ave),ignore.case = TRUE,perl = TRUE,value = TRUE)



###############DGE test for five time points#################

plot_flag = FALSE
label_flag = FALSE
n_plot = ifelse(plot_flag,2,1)

width = 6.635417; height = 6.739583
if(plot_flag){width = 12.635417; height = 6.739583}
#width = 600; height = 600 #for tiff pix

##ESC
pdf("out_DGE_qlf/qlf.kovswt.ESC.pairwise.pdf",width = width,height = height,onefile = TRUE,useDingbats = FALSE )
#tiff("out_DGE_qlf/qlf.kovswt.ESC.pairwise.tif",width = width,height = height )
par(mfrow=c(1,n_plot)) 
sample <- "qlf.kovswt.ESC"
res.table <- qlf.kovswt.ESC$table
#res.table <- lrt.kovswt.ESC$table
res.plotV <- plotVolcano ( res.table = res.table, sample = sample, 
              FC.cutoff = FC.cutoff, padj.cutoff = padj.cutoff, threshold.exp =threshold.exp,
              plot = plot_flag )
drawSimpleDE (data  = logCPM.ave[,c(1,6)],sample = sample,select.up = res.plotV$select.up, select.down = res.plotV$select.down, 
              padj.cutoff = padj.cutoff, threshold.FC = FC.cutoff, threshold.exp =threshold.exp,label=label_flag)

#hist(res.table$PValue)
#hist(p.adjust(res.table$PValue,method = "BH") )
dev.off()

##6h
pdf("out_DGE_qlf/qlf.kovswt.6h.pairwise.pdf",width = width,height = height,onefile = TRUE,useDingbats = FALSE)
#tiff("out_DGE_qlf/qlf.kovswt.6h.pairwise.tif",width = width,height = height,res= 150)
par(mfrow=c(1,n_plot)) 
sample <- "qlf.kovswt.6h"
res.table <- qlf.kovswt.6h$table
res.plotV <- plotVolcano ( res.table = res.table, sample = sample, 
                           FC.cutoff = FC.cutoff, padj.cutoff = padj.cutoff,  threshold.exp =threshold.exp,
                           plot = plot_flag )
drawSimpleDE (data  = logCPM.ave[,c(2,7)],sample = sample,select.up = res.plotV$select.up, select.down = res.plotV$select.down, 
              padj.cutoff = padj.cutoff, threshold.FC = FC.cutoff, threshold.exp =threshold.exp,label=label_flag)
dev.off()


##12h
pdf("out_DGE_qlf/qlf.kovswt.12h.pairwise.pdf",width = width,height = height,onefile = TRUE,useDingbats = FALSE)
#tiff("out_DGE_qlf/qlf.kovswt.12h.pairwise.tif",width = width,height = height,res= 150)
par(mfrow=c(1,n_plot)) 
sample <- "qlf.kovswt.12h"
res.table <- qlf.kovswt.12h$table
res.plotV <- plotVolcano ( res.table = res.table, sample = sample, 
                           FC.cutoff = FC.cutoff, padj.cutoff = padj.cutoff,  threshold.exp =threshold.exp,
                           plot = plot_flag )
drawSimpleDE (data  = logCPM.ave[,c(3,8)],sample = sample,select.up = res.plotV$select.up, select.down = res.plotV$select.down, 
              padj.cutoff = padj.cutoff, threshold.FC = FC.cutoff, threshold.exp =threshold.exp,label=label_flag)
dev.off()


##24h
pdf("out_DGE_qlf/qlf.kovswt.24h.pairwise.pdf",width = width,height = height,onefile = TRUE,useDingbats = FALSE)
#tiff("out_DGE_qlf/qlf.kovswt.24h.pairwise.tif",width = width,height = height,res= 150)
par(mfrow=c(1,n_plot)) 
sample <- "qlf.kovswt.24h"
res.table <- qlf.kovswt.24h$table
res.plotV <- plotVolcano ( res.table = res.table, sample = sample, 
                           FC.cutoff = FC.cutoff, padj.cutoff = padj.cutoff,  threshold.exp =threshold.exp,
                           plot = plot_flag )
drawSimpleDE (data  = logCPM.ave[,c(4,9)],sample = sample,select.up = res.plotV$select.up, select.down = res.plotV$select.down, 
              padj.cutoff = padj.cutoff, threshold.FC = FC.cutoff, threshold.exp =threshold.exp,label=label_flag)
dev.off()

##48h
pdf("out_DGE_qlf/qlf.kovswt.48h.pairwise.pdf",width = width,height = height,onefile = TRUE,useDingbats = FALSE)
#tiff("out_DGE_qlf/qlf.kovswt.48h.pairwise.tif",width = width,height = height,res= 150)
par(mfrow=c(1,n_plot)) 
sample <- "qlf.kovswt.48h"
res.table <- qlf.kovswt.48h$table
res.plotV <- plotVolcano ( res.table = res.table, sample = sample, 
                           FC.cutoff = FC.cutoff, padj.cutoff = padj.cutoff,  threshold.exp =threshold.exp,
                           plot = plot_flag )
drawSimpleDE (data  = logCPM.ave[,c(5,10)],sample = sample,select.up = res.plotV$select.up, select.down = res.plotV$select.down, 
              padj.cutoff = padj.cutoff, threshold.FC = FC.cutoff, threshold.exp =threshold.exp,label=label_flag)

dev.off()

##collect files
dname = paste("out_DGE_qlf/FC_",FC.cutoff,"_",ifelse(plot_flag,"two","one"),"_plots",ifelse(label_flag,"","_nolabels"),sep=""   )
if(dir.exists(dname)){cat('dir exists\n')}else{dir.create(dname,showWarnings =TRUE) }
system(paste("mv -f out_DGE_qlf/qlf*.pdf out_DGE_qlf/select*.mat ",dname) )
cat("results in ",dname,"\n")


##results in /home/mjwang/pwd/Quan-Cnot8-RNASEQ_all_mm9/04.analysis_featureCnt/edgeR_DESeq2_logTPM/out_DGE_qlf


#save.image(file = ".RData.work_edgeR_logTPM")
load(".RData.work_edgeR_logTPM")

################################################GO enrichment (topGO, GOstats)#####################################################
#########https://www.bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf##########
#https://davetang.org/muse/2010/11/10/gene-ontology-enrichment-analysis/
library(topGO)
library(GOstats)
library(GO.db)
require(org.Mm.eg.db) #a sqlite db

############1, use edgeR goana##########
go.48h = edgeR::goana(qlf.kovswt.48h,species='Mm')


#############2, use GOstats############

#i. get universe geneset
entrez_object <- org.Mm.egGO
mapped_genes <- mappedkeys(entrez_object)
#use all the genes with GO terms as the universe
universe <- mapped_genes #24386

#ii, map symbol gene name to Entrezid
##keytypes(org.Mm.eg.db) #see available id types
y$genes$EntrezID <- AnnotationDbi::mapIds(org.Mm.eg.db,keys=as.character(y$genes$genes),keytype = "SYMBOL", column = "ENTREZID" )
#returned 1:many mapping
y$genes$ENSEMBL <- mapIds(org.Mm.eg.db,keys=as.character(y$genes$genes),keytype = "SYMBOL", column = "ENSEMBL")
#returned 1:many mapping

#qlf.DGE.48h = read.table( "out_DGE_qlf/FC_1.5_one_plots/select.down.edgeR.test.log2TPM.qlf.kovswt.48h.FC.1.5.mat", header = TRUE ) #199
qlf.DGE.48h = read.table( "out_DGE_qlf/FC_1.5_one_plots/select.up.edgeR.test.log2TPM.qlf.kovswt.48h.FC.1.5.mat", header = TRUE ) #313

qlf.DGE.48h.genes = rownames(qlf.DGE.48h) 

qlf.DGE.48h.genes.entrez = AnnotationDbi::mapIds(org.Mm.eg.db,keys=qlf.DGE.48h.genes,keytype = "SYMBOL", column = "ENTREZID" )
#returned 1:1 mapping between keys and columns

##iii, do GOstats HyperGTest
params <- new('GOHyperGParams',
              geneIds = qlf.DGE.48h.genes.entrez,
              universeGeneIds = universe,
              ontology = 'BP',
              pvalueCutoff = 0.001,
              conditional = FALSE,
              testDirection = 'over',
              annotation = "org.Mm.eg.db"
)

mmOver <- GOstats::hyperGTest(params)
result <- summary(mmOver)

write.table(result,"GOstats.results.48h.up.FC1.5.txt",col.names = TRUE,row.names = TRUE, quote = FALSE,sep='\t')
##many metabolic process genes? upregulated in cnot8 KO, why? these genes should be normally degraded?


#######GSEA gene set analysis#######


save.image(file = ".RData.work_edgeR_logTPM")
#load(".RData.work_edgeR_logTPM")



#############end of script#########















logCPM <- cpm(y,prior.count = 2, log=TRUE)
# 
ave12 <- rowMeans(logCPM[,1:2]) #WT.ESC
ave34 <- rowMeans(logCPM[,3:4]) #KO.ESC
#qlf <- qlf.kovswt.ESC
lrt <- lrt.kovswt.ESC

ave12 <- rowMeans(logCPM[,5:6]) #WT.6h
ave34 <- rowMeans(logCPM[,7:8]) #KO.6h
#qlf <- qlf.kovswt.6h
lrt <- lrt.kovswt.6h

ave12 <- rowMeans(logCPM[,9:10]) #WT.12h
ave34 <- rowMeans(logCPM[,11:12]) #KO.12h
#qlf <- qlf.kovswt.12h
lrt <- lrt.kovswt.12h

ave12 <- rowMeans(logCPM[,13:14]) #WT.24h
ave34 <- rowMeans(logCPM[,15:16]) #KO.24h
#qlf <- qlf.kovswt.24h
lrt <- lrt.kovswt.24h

ave12 <- rowMeans(logCPM[,17:18]) #WT.48h
ave34 <- rowMeans(logCPM[,19:20]) #KO.48h
#qlf <- qlf.kovswt.48h
lrt <- lrt.kovswt.48h

#select.up <- (qlf$table$PValue <=0.01 & qlf$table$logFC >=1 & ave34 >=2 & !is.na(qlf$genes$Symbol) ) #2 fold change, log2CPM >=2
#select.down <- (qlf$table$PValue <=0.01 & qlf$table$logFC <= -1 & ave12 >=2 & !is.na(qlf$genes$Symbol)  ) #2 fold change, log2CPM >=2

#genes.up <- qlf$genes[select.up,]
#genes.down <- qlf$genes[select.down,]

#select.up <- (lrt$table$PValue <=0.01 & lrt$table$logFC >=1 & ave34 >=2 & !is.na(lrt$genes$Symbol) ) #2 fold change, log2CPM >=2
#select.down <- (lrt$table$PValue <=0.01 & lrt$table$logFC <= -1 & ave12 >=2 & !is.na(lrt$genes$Symbol)  ) #2 fold change, log2CPM >=2

#genes.up <- lrt$genes[select.up,]
#genes.down <- lrt$genes[select.down,]

#select.up <- (qlf$table$PValue <=0.01 & qlf$table$logFC >=1 & ave34 >=1 & !is.na(qlf$genes$Symbol) ) #2 fold change, log2CPM >=1
#select.down <- (qlf$table$PValue <=0.01 & qlf$table$logFC <= -1 & ave12 >=1 & !is.na(qlf$genes$Symbol)  ) #2 fold change, log2CPM >=1

#genes.up <- qlf$genes[select.up,]
#genes.down <- qlf$genes[select.down,]


#select.up <- (lrt$table$PValue <=0.01 & lrt$table$logFC >=1 & ave34 >=1 & !is.na(lrt$genes$Symbol) ) #2 fold change, log2CPM >=1
#select.down <- (lrt$table$PValue <=0.01 & lrt$table$logFC <= -1 & ave12 >=1 & !is.na(lrt$genes$Symbol)  ) #2 fold change, log2CPM >=1

#genes.up <- lrt$genes[select.up,]
#genes.down <- lrt$genes[select.down,]

#select.up <- (qlf$table$PValue <=0.01 & qlf$table$logFC >=1 & ave34 >=1 ) #2 fold change, log2CPM >=1
#select.down <- (qlf$table$PValue <=0.01 & qlf$table$logFC <= -1 & ave12 >=1 ) #2 fold change, log2CPM >=1

#genes.up <- qlf$genes[select.up,]
#genes.down <- qlf$genes[select.down,]

select.up <- (lrt$table$PValue <=0.01 & lrt$table$logFC >=1 & ave34 >=1 ) #2 fold change, log2CPM >=1
select.down <- (lrt$table$PValue <=0.01 & lrt$table$logFC <= -1 & ave12 >=1 ) #2 fold change, log2CPM >=1

genes.up <- lrt$genes[select.up,]
genes.down <- lrt$genes[select.down,]


r <- format(cor(ave12,ave34,method = 'spearman'),digits = 3)
plot(ave12,ave34,pch=20,cex=0.5,col='lightgrey',xlab="WT 48h log2CPM",ylab="KO 48h log2CPM")
points(ave12[select.up],ave34[select.up],pch=20,cex=0.5,col="red" )
points(ave12[select.down],ave34[select.down],pch=20,cex=0.5,col="blue" )
text(x = 0, y = 5,labels = paste("n=",nrow(genes.up),sep = ""), col='black')
text(x = 5, y = 1,labels = paste("n=",nrow(genes.down),sep = ""), col='black')
text(x = -1, y = 10,labels = paste("R=",r,sep = ""), col='black')

#abline(a=2.32,b=1,col='darkgrey',lty=2,lwd=1.2,untf = FALSE)
#abline(a=-2.32,b=1,col='darkgrey',lty=2,lwd=1.2,untf = FALSE)
abline(a=1,b=1,col='darkgrey',lty=2,lwd=1.2,untf = FALSE)
abline(a=-1,b=1,col='darkgrey',lty=2,lwd=1.2,untf = FALSE)
#abline(a=1,b=1,col='darkgrey',lty=2,lwd=1.2,untf = FALSE)
#abline(a=-1,b=1,col='darkgrey',lty=2,lwd=1.2,untf = FALSE)
#abline(h=0,v=0);
abline(a=0,b=1)











