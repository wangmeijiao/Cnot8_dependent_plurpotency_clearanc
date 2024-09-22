#library(Seurat,lib.loc = "/home/mjwang/.conda/envs/myenv/lib/R/library_opt/") #3.2.2

library(Seurat)
library(viridis)
library(magrittr)
library(dplyr)
library(ggplot2)
#library(matrixStats)
library(ComplexHeatmap)
#library(viridis)
#library('ggplot2')

#source("../../04.scRNA_scATAC/greenleaf_archR_doublets/ArchR/R/ColorPalettes.R")


#library(harmony), no need to do this
library(hrbrthemes)
library(patchwork)


library(reticulate)#1.2
library(Matrix)
#library(Seurat)
#library(dplyr)
library(PRROC)
library(pbapply)
#source('utility.R')
set.seed(2020)





############Buen colors########
library(BuenColors)
color_set0 <- jdb_color_maps #17 different colors
names(color_set0) <- NULL
#plot(1:17,1:17,pch = 19, cex = 5,col=jdb_color_maps)

#discrete colors
color_set1 <- jdb_palette("solar_extra") #9 discrete but gradient colors
color_set2 <- jdb_palette("brewer_spectra") #9 discrete but gradient colors
color_set3 <- jdb_palette("flame_light") #9 discrete but gradient colors, good!

color_set3_ext12 <- colorRampPalette(colors = as.character(color_set3))(12)
color_set3_ext17 <- colorRampPalette(colors = as.character(color_set3))(17)

#############ArchR colors############
#hmcols <- colorRamps::blue2green2red(length(bks) ) #colors
color_peak <- ArchR::paletteContinuous(set = 'solarExtra',n=256,reverse=FALSE)  
color_tfdev = ArchR::paletteContinuous(set = 'blueYellow',n=257,reverse=FALSE)                       
#color_ga <- paletteContinuous(set='solarExtra',n=257,reverse=FALSE) 
#color_ga <- paletteContinuous(set='horizon',n=257,reverse=FALSE)                      
#color_ga <- paletteContinuous(set='horizonExtra',n=257,reverse=FALSE)  #good        
#color_ga <- paletteContinuous(set='greenBlue',n=257,reverse=FALSE)
#color_ga <- paletteContinuous(set='blueYellow',n=257,reverse=FALSE)
color_ga <- ArchR::paletteContinuous(set='greyMagma',n=257,reverse=FALSE)


#########customized colors########
color_snap = c('1'='grey','2'='#E31A1C','3'='#FFD700','4'='#771122','5'='#777711','6'='#1F78B4','7'='#68228B','8'='#AAAA44','9'='#60CC52','10'='#771155','11'='#DDDD77','12'='#774411','13'='#AA7744','14'='#AA4455','15'='#117744')
names(color_snap) <- NULL

color_signac = c(
'0'='#E6D55E','1'='#792B8A','2'='#DA703D','3'='#9DC8E5','4'='#BA273C','5'='#C2C184','6'='#7F8084','7'='#65AB53','8'='#D082AF','9'='#496EAB','10'='#DE896D','11'='#491F8B','12'='#E1AD49','13'='#8E1B85','14'='#E7EE77','15'='#7D1A1D','16'='#96B355')
names(color_signac) <- NULL

color_good <- c("#E7D654", "#6F1482" ,"navy", "#9CC8E6", "#BC2439" ,"#C2C281" ,"#7F8084", 
                "#63AC4E", "#D181B0" ,"#476DAD","#DF8969", "#49188D", "#E2AC3B" ,"#8F1287" ,
                "#E7F06F" ,"#7E191A" ,"#95B44F", "#6C3719" ,"#CA362E" ,"#2B3918","#1E1E1E" )

color_gradient_my <- c(
    rgb(5,48,97,maxColorValue = 255),
    rgb(42,113,178,maxColorValue = 255),
    rgb(147,198,222,maxColorValue = 255),
    rgb(239,243,245,maxColorValue = 255),
    rgb(253,219,199,maxColorValue = 255),
    rgb(214,96,77,maxColorValue = 255),
    rgb(121,5,34,maxColorValue = 255)

)

color_pyreds <- c(
"#fff5f0","#fef4ef","#fef3ee","#fef3ed","#fef2ec","#fef1eb","#fef1ea","#fef0e9","#feefe8","#feefe7","#feeee6","#feede5","#feede4","#feece3","#feebe2","#feebe1","#feeae0","#fee9e0","#fee9df","#fee8de"
)


#color_cellranger <- c("#820610",'#C60F1E','#F7592D','#FBB33D','#FDEC4A','#DFFB52','#B4FB68','#8BFB89','#70FBA7','#52FDCE','#3DEEF6','#33CEE5','#2BAED7','#2392CA','#1A70B9','#1354AB','#0A2E94','#061D88','#07177E')


color_cellranger <-c('#820610','#C50F1E','#F42428','#F86D30','#FBB33D','#FCFB4E','#C0FB61','#87FB8F','#41FEF9','#2BAED7','#155CB1','#08238D')

###select one global color set###
color <- color_good



####
sample = "Cnot8_combine"


##eload object

#placenta <- readRDS( "limb.final.rds")




#####

data.dir = c('ESC' = "../01.cellranger/ESC/out_ESC/outs/filtered_feature_bc_matrix/",
        'LC_48h' = "../01.cellranger/LC-48h/out_LC-48h/outs/filtered_feature_bc_matrix/"
       )
#genes identical

##data.dir = '../01.aggregation/limb-aggre_nonorm/filtered_feature_bc_matrix/'
#data.dir = '../01.aggregation/limb-aggre/filtered_feature_bc_matrix/'
#data.dir = "../../placenta_10X/01.data_cellranger_RNASEQ/limb-1/filtered_feature_bc_matrix/"
#data.dir = "../01.data_cellranger_RNASEQ/filtered_feature_bc_matrix/"
placenta.count <- Read10X(data.dir = data.dir,strip.suffix=TRUE)
#32286 x 25611
#29324 x 18970

#modify cellid
cellid <- colnames(placenta.count)
idx1  <- grepl(pattern = '^ESC',x = cellid)
idx2  <- grepl(pattern = '^LC_48h',x = cellid)

unique(idx1 + idx2 ) #1

# sample <- vector()
# sample[idx1] <- 'U'
# sample[idx2] <- 'D'
# sample <- factor(sample,levels = c('U','D'))
# sample
#    U    D 
# 9055 8402 


cellid[idx1] <- paste0(cellid[idx1],'-1',sep='')
cellid[idx2] <- paste0(cellid[idx2],'-2',sep='')

cellid <- gsub(pattern = "^ESC_",replacement = '',x = cellid,)
cellid <- gsub(pattern = "^LC_48h_",replacement = '',x = cellid,)

sum(grepl(pattern = '-1$',x = cellid))
#13184
#11930 #9055

sum(grepl(pattern = '-2$',x = cellid))
#12427
#7040 #8402

colnames(placenta.count) <- cellid

#create seurat obj
placenta <- CreateSeuratObject(counts = placenta.count, project = "Cnot8_combine", min.cells = 3, min.features = 200)
#22643 x 25608

#17641 x 17938

#23028 x 16989

#placenta <- CreateSeuratObject(counts = placenta.count, project = "limb_combine_28d")
#33538 x 17457

##add sample name in meta data

cellid <- colnames(placenta)
cellid <- gsub(pattern = ".+-1$",replacement = 'ESC',x = cellid,)
cellid <- gsub(pattern = ".+-2$",replacement = 'LC_48h',x = cellid,)
#cellid <- gsub(pattern = ".+-3$",replacement = 'TM_D1',x = cellid,)
#cellid <- gsub(pattern = ".+-4$",replacement = 'TM_D2',x = cellid,)




placenta[["sample"]]  <- factor(cellid,levels = c('ESC','LC_48h'))
table(placenta[["sample"]])
 ESC LC_48h 
 13181  12427

#  U     D 
# 11551  6387 

#   D    U 
# 8402 9055 

#   D    U 
# 8282 8707 

grep(pattern = '^mt', x=rownames(placenta) ,value = TRUE )
#'mt-Nd1''mt-Nd2''mt-Co1''mt-Co2''mt-Atp8''mt-Atp6''mt-Co3''mt-Nd3''mt-Nd4l''mt-Nd4''mt-Nd5''mt-Nd6''mt-Cytb'

grep(pattern = '^tdTomato', x=rownames(placenta) ,value = TRUE )
#'tdTomato'

placenta[["percent.mt"]] = PercentageFeatureSet(placenta,pattern = "^mt-")

VlnPlot(placenta,features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),pt.size=0.1 )

plot1 <- FeatureScatter(placenta, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(placenta, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
options(repr.plot.height=3,repr.plot.width=7.5)
CombinePlots(plots = list(plot1, plot2))

placenta.raw <- placenta
#filter
##placenta <- subset(placenta, subset = nFeature_RNA > 2000 & nFeature_RNA < 5000 & percent.mt < 10 ) #use 5000 instead of 5500, 
placenta <- subset(placenta, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 10 )

#placenta <- subset(placenta, subset = nFeature_RNA > 1500 & nFeature_RNA < 5000 & nCount_RNA >3500 & nCount_RNA < 5e4 & percent.mt < 5) #the same condition in scanpy

# ##filter by scanpy result
# cluster_scanpy <- read.table("../02.scanpy_bbknn_combine/cluster_df_add.txt") #16045 x 13
# #U -> 1
# #D -> 2

# table(rownames(cluster_scanpy) %in% colnames(placenta) ) # TRUE 16045 
# placenta <- subset(placenta,cells = rownames(cluster_scanpy))
# #23028 x 16045
# table(placenta[['sample']])
#   U    D 
# 8241 7804 

placenta
2643 features across 20930 samples
#22643 features across 17665


table(placenta[['sample']])
 ESC LC_48h 
 11489   9441

#ESC LC_48h 
#9932   7733

#after filter #24307 x 11560  #5646 #6471
options(repr.plot.height=7.5,repr.plot.width=12.5)
VlnPlot(placenta,features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),pt.size=0.1 )



##normalize
placenta <- NormalizeData(placenta, normalization.method = "LogNormalize", scale.factor = 10000)

##identification of highly variable features
placenta <- FindVariableFeatures(placenta, selection.method = "vst", nfeatures = 2000)

top30 = head(VariableFeatures(placenta),30  )
plot1 = VariableFeaturePlot(placenta   )
plot2 <- LabelPoints(plot = plot1, points = top30, repel = TRUE)
options(repr.plot.height=7.5,repr.plot.width=15)
CombinePlots(plots = list(plot1, plot2))

##scale the data
all.genes <- rownames(placenta)
placenta <- ScaleData(placenta, features = all.genes)

#linear dimensional reduction: PCA
placenta <- RunPCA(placenta, features = VariableFeatures(object = placenta)) #total 50

#quick look the PCs
options(repr.plot.height=7.5,repr.plot.width=7.5)
VizDimLoadings(placenta, dims = 1:5, reduction = "pca")
DimPlot(placenta, reduction = "pca")
DimHeatmap(placenta, dims = 1:15, cells = 500, balanced = TRUE)

#choose PC number
placenta <- JackStraw(placenta, num.replicate = 100,dims = 35) #20 by default
placenta <- ScoreJackStraw(placenta, dims = 1:25)
#JackStrawPlot(placenta, dims = 1:15)
ElbowPlot(placenta,ndims = 40) #use n_pcs = 40

#cluster using KNN graph Louvain method
placenta <- FindNeighbors(placenta, dims = 1:40)  #use 40 PCs from ElbowPlot
#placenta <- FindClusters(placenta, resolution = 1)
placenta <- FindClusters(placenta, resolution = 0.5) #
placenta <- FindClusters(placenta, resolution = 0.3) #use this

######## plot cluster size##
Idents(placenta)
all.equal(Idents(placenta),placenta@meta.data$seurat_clusters) #name id diff
cluster.stat = table(placenta@meta.data$seurat_clusters)
total <- sum(cluster.stat)  #15107 #9461
res.bp <- barplot(cluster.stat,col = color_good,las=2,names.arg = paste0("cluster",names(cluster.stat)),
                  ylab = "cell count", main='cluster cell number' )
text(x=res.bp[,1],y=cluster.stat+55,labels = cluster.stat,srt=90,xpd=TRUE)
text(x=16,y=1500,labels = paste("total = ",total,sep='') )
###

#system("which python")
#Sys.getenv()
#Sys.setenv("PATH"="/home/mjwang/bin:/usr/local/bin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin")
#system("which python") #use /home/mjwang/bin/python
#reticulate::py_config()

##non-linear dimensional reduction (UMAP/tSNE)
placenta = RunUMAP(placenta,dims=1:40,seed.use = 123) 
placenta =  RunTSNE(placenta,dims=1:40,seed.use = 123)




###################doublet removal with python package scrublet in r#################
#https://github.com/xnnba1984/Doublet-Detection-Benchmark/tree/master/real_data_benchmark
#https://escholarship.org/content/qt6gr648np/qt6gr648np_noSplash_32134c0e501d9fdf306702f2a4dcc052.pdf


# read python module
scr <- import('scrublet')
scipy.io <- import('scipy.io')
np <- import('numpy')
os <- import('os')

count <- placenta@assays$RNA@counts
result <- scr$Scrublet(counts_matrix = t(count), random_state = 10L)
results <- result$scrub_doublets(min_counts=2, 
                                 min_cells=3, 
                                 min_gene_variability_pctl=85, 
                                 n_prin_comps=30L)

##
#result$plot_histogram()
#result$set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, n_neighbors=10, min_dist=0.3,random_state=123))
#result$plot_embedding('UMAP', order_points=True)
###

score <- as.vector(results[[1]])
flag_doublet <- as.vector(results[[2]])

table(flag_doublet)
flag_doublet
FALSE  TRUE 
17042   623

placenta <- AddMetaData(object = placenta,metadata = flag_doublet, col.name = 'scrublet')


# # calculate auprc and auroc
# fg <- score[label==1]
# bg <- score[label==0]
# pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T); pr$auc.integral
# roc <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T); roc$auc


FeaturePlot(placenta, features = c('scrublet'), reduction = "umap",cols = c('grey','red') )


placenta <- subset(x = placenta, subset = scrublet == FALSE)

###############





##########quick look dimension reduction result and quanlity control result#######
#saveRDS(placenta, file = "placenta.limb-new.rds",compress = TRUE)
#placenta <- readRDS(file = "placenta.limb-new.rds") #dgCMatrix

options(repr.plot.height=7.5,repr.plot.width=7.5)
#DimPlot(placenta, reduction = "tsne",label=TRUE,cols = color_good,label.size = 8,pt.size = 1)
DimPlot(placenta, reduction = "umap",label=TRUE,cols = color_good,label.size = 8,pt.size = 1)

DimPlot(placenta, reduction = "tsne",label=TRUE,cols = color_good,label.size = 8,pt.size = 1)


#
'tdTomato' %in% VariableFeatures(object = placenta) #FALSE
FeaturePlot(placenta, features = c('tdTomato'), reduction = "umap", pt.size = 2 )


FeaturePlot(placenta, features = c('nCount_RNA'), reduction = "umap",cols = viridis(12,option = 'D') )
FeaturePlot(placenta, features = c('nFeature_RNA'), reduction = "umap",cols = viridis(12,option = 'D') )
FeaturePlot(placenta, features = c('percent.mt'), reduction = "umap",cols = viridis(12,option = 'D') )




# FeaturePlot(placenta, features = c('nCount_RNA'), reduction = "tsne",cols = viridis(12,option = 'D') )
# FeaturePlot(placenta, features = c('nFeature_RNA'), reduction = "tsne",cols = viridis(12,option = 'D') )
# FeaturePlot(placenta, features = c('percent.mt'), reduction = "tsne",cols = viridis(12,option = 'D') )


options(repr.plot.height=7.5,repr.plot.width=15)
FeaturePlot(placenta, split.by  = c('sample'), features = c('nFeature_RNA'), reduction = "umap",cols = viridis(12,option = 'D'))
FeaturePlot(placenta, split.by  = c('sample'), features = c('percent.mt'), reduction = "umap",cols = viridis(12,option = 'D'))

#split by sample
placenta_ESC <- subset(placenta,subset = sample == 'ESC')
placenta_LC_48 <- subset(placenta,subset = sample == 'LC_48h')

options(repr.plot.height=7.5,repr.plot.width=7.5)
DimPlot(placenta_ESC, reduction = "umap",label=TRUE,label.size = 8,pt.size = 1)
DimPlot(placenta_LC_48, reduction = "umap",label=TRUE,label.size = 8,pt.size = 1)



#############



###########marker gene plotting##########


marker.genes.list <- list(
    
# 'Preadipocyte' = c('PDGFRA', 'ITGB', 'CD34', 'CD24', 'EN1', 'MYF5'),  
# 'MSC'= c('PRRX1', 'NT5E', 'ENG'),
 
'MSC'= c('PRRX1', 'NT5E', 'ENG','PDGFRA', 'CD34', 'CD24', 'EN1', 'MYF5'),

'Dermal fibroblast'= c ('COL3A1', 'TWIST2', 'DPT', 'COL1A1', 'LUM', 'PTN', 'LTIH5', 'NFIA', 'MDK', 'TRIL', 'APCDD1', 'PAX1'),
'Muscle'= c ('PAX7', 'PAX3', 'MYOD', 'MYOG', 'MYF5', 'MYH1', 'MYH2', 'MYH3', 'MYH4', 'MYH7'),
'Blood'= c ('HBB', 'CD36', 'TFRC', 'GYPA'),
'Cartilage'= c ('COL10A1', 'COL2A1', 'AGC1', 'SOX9'),
'Tendon'= c ('TNMD', 'SCX', 'MKX', 'SIX2'),
'Macrophage'= c ('CD68', 'C1QB', 'C1QA', 'C1QC', 'CD48', 'HMOX1', 'CCL3', 'CD36', 'CD14', 'CCL4', 'HMOX1', 'CCL2'),
'Endothelial'= c ('PECAM1', 'KDR'),
'Skin'= c ('KRT14', 'KRT15', 'KRT17', 'KRT5', 'PDGFRA' , 'EPAM', 'LHX2', 'BMP7', 'CXCL14', 'KRT1', 'KRT4', 'CDH1', 'EPCAM', 'MEIS2', 'KRT19', 'KRT8', 'KRT18', 'KLF5'),
'Neuron'= c ('NKTR', 'MEG3', 'ATRX', 'CCNL1', 'NSD3', 'RBM25', 'RBM6', 'FGFR1', 'CCNL2', 'PAXBP1') 

)

marker.genes <- unlist(marker.genes.list)

marker.genes <- stringr::str_to_title(marker.genes)

marker.genes <- marker.genes[marker.genes %in% rownames(placenta)]


marker.genes <- c('Col3a1','Krt8','Pdgfra')

options(repr.plot.height=7.5,repr.plot.width=7.5)
for(gene in marker.genes){

    res.p <- FeaturePlot(placenta, features = gene, reduction = "umap",slot = 'scale.data',cols = c("lightgrey","red" ) ,pt.size = 0.5)+ #, min.cutoff = 1, max.cutoff = 2.5)+
        ##p<- FeaturePlot(placenta_D, features = gene, reduction = "umap",slot = 'data',cols = c("lightgrey", "darkred") )+ #c("lightgrey", "#ff0000", "#00ff00") c("lightgrey", "darkred")
          #scale_color_gradientn(colours = color_ga)+
          #scale_color_gradientn(colours = color_peak)+

          #scale_color_gradientn(colours = color_gradient_my)+
          #scale_color_gradientn(colours = color_pyreds)+
          theme(
                legend.position = 'right',
                axis.text=element_blank(), 
                axis.title = element_text(size = 15, face = "bold"),
                axis.ticks = element_blank(),
                panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(),
                panel.background = element_rect(color="black", fill = NA,size=0.8),
        #         panel.background = element_rect(fill = "white", colour = "white", 
        #                 size = rel(1)),
                #panel.border = element_blank(),
                plot.title = element_text(size = 15, face = "bold"),
                #complete = TRUE
                plot.margin = unit(c(1,1,1,1), "lines") #add margin benifit the outside box frame
               ) 
    print(res.p)
}
#########


saveRDS(placenta,'Cnot8_ESC_LC48h_WT_KO_tdTomato.rds')


###########do filtering############



placenta.filter <- subset(placenta, idents = c('0','1','2','3','7','8'))

options(repr.plot.height=7.5,repr.plot.width=7.5)
DimPlot(placenta.filter, reduction = "umap",label=TRUE,cols = color_good,label.size = 8,pt.size = 1)
DimPlot(placenta.filter, reduction = "tsne",label=TRUE,cols = color_good,label.size = 8,pt.size = 1)




#placenta.filter = RunUMAP(placenta.filter,dims=1:40,seed.use = 200,min.dist = 0.5,spread = 0.3) 
#placenta.filter =  RunTSNE(placenta.filter,dims=1:40,seed.use = 2021)

#options(repr.plot.height=7.5,repr.plot.width=7.5)
#DimPlot(placenta.filter, reduction = "umap",label=TRUE,cols = color_good,label.size = 8,pt.size = 1)
#DimPlot(placenta.filter, reduction = "tsne",label=TRUE,cols = color_good,label.size = 8,pt.size = 1)





##########################tuning umap###########################

data.use <- placenta.filter@reductions$pca@cell.embeddings[,1:40]
cluster <- Idents(placenta.filter)

res.umap <- list()
res.p <- list()
for(mdist in seq(from = 0.1,to=1.2,by=0.2) ){
  for(spread in seq(from = 0.8,to=1.5,by=0.2) ){
 
    set.seed(100)
    dr_umap <- uwot::umap(X = data.use, n_components = 2,n_threads = 2,min_dist=mdist, spread=spread)
    #dr_umap <- uwot::umap(X = data.use, n_components = 2,n_threads = 2,min_dist=0.4, spread=1)
    dr_umap.df <- data.frame(cluster=cluster,dr_umap)

    #plot(dr_umap,pch = 19, cex =0.2,ylim=c(-10,10),xlim=c(10,-10))

    p <- ggplot2::ggplot(data = dr_umap.df , aes(x = X1, y = X2, color = cluster) ) + 
      geom_point() +
      scale_color_manual(values = color_good) +
      #xlim(c(10,-10)) +
      #ylim(c(10,-10)) +
      ggtitle(paste0('seed:',100,' mdist:',mdist,' spread:',spread) )+
      theme_bw()

    print(p)
    
    res.p[[paste0('mdist',mdist,'-spread',spread)]] <- p
    res.umap[[paste0('mdist',mdist,'-spread',spread)]] <- dr_umap.df

  }
}


saveRDS(res.umap,'res.umap.tunning.rds')
saveRDS(res.p,'res.plot.tunning.rds')

res.p[['mdist0.9-spread1']]
umap.tunning <-  res.umap[['mdist0.9-spread1']]

#####filtering and rotate######



dr_umap.df <- umap.tunning
colnames(dr_umap.df) <- c('cluster','UMAP_1','UMAP_2')

all.equal(dr_umap.df$cluster, cluster,check.attributes = FALSE) #TRUE


##cluster distribution
options(repr.plot.width=10,repr.plot.height=8)
par(mfrow = c(2,3))
for (i in levels(dr_umap.df$cluster)){

dotDistri(cluster = dr_umap.df, id = i)

}




###============fix final cluster result by dotDist + dotclean () ====================###
centers <- dr_umap.df %>% dplyr::group_by(cluster) %>% dplyr::summarize(x = median(x = UMAP_1), 
        y = median(x = UMAP_2))


dotDist = function (cluster = NULL, id = NULL,center = centers){
    colids = colnames(cluster) #must 'cluster','dim1','dim2' and rowname
    cluster.sel = cluster[ cluster[,1] == id,]
    n_sel = nrow(cluster.sel)
    #cat ('select for cluster ',id,' n = ',n_sel," \n")
    #color = 'red'
    #color = ifelse('cluster_sg' == id,'red','navy')
#     plot(cluster[,2],cluster[,3],pch = 16, type='p',col='grey',cex=0.5,xlab=colids[2],ylab=colids[3],main=paste(" cells cluster ",id," of ",colids[1],"\nn = ",n_sel,sep=''),cex.main = 2.25,xaxt = 'n' ) 
#     points(cluster.sel[,2],cluster.sel[,3],pch = 16, cex=0.5,col=color)
#     points(center[center$cluster == id,2:3],pch = 16, cex=1,col='black')
       
    dx <- cluster.sel[,2] - as.numeric(center[center$cluster == id,2])
    dy <- cluster.sel[,3] - as.numeric(center[center$cluster == id,3])
    d <- sqrt(dx**2 + dy**2)
    
    #options(repr.plot.height=15,repr.plot.width=15)
    hist(d,breaks=100,main=paste("cluster ",id,sep=''),cex.main=2)
    d.quantile <- quantile(d,prob = seq(0,1,0.1))
    abline(v=d.quantile,lty=2,lwd=1,col='red')
    #flag <- d <= as.numeric(d.quantile['80%'])
    #flag <- d <=2
   cat(paste("cluster ",id," ok",sep='') )
    return(d.quantile)
}


dotClean = function (cluster = NULL, id = NULL,center = centers,d.filter = d.filter){
    #return a flag.dist for each cluster
    colids = colnames(cluster) #must 'cluster','dim1','dim2'
    cluster.sel = cluster[ cluster[,1] == id,]
    d.filter.sel = d.filter[id]
    n_sel = nrow(cluster.sel)
    #cat ('select for cluster ',id,' n = ',n_sel," \n")
    color = 'red'
    cat(paste("cluster ",id," filtering dist:",d.filter.sel,sep=''),'\n')
    #color = ifelse('cluster_sg' == id,'red','navy')
#     plot(cluster[,2],cluster[,3],pch = 16, type='p',col='grey',cex=0.5,xlab=colids[2],ylab=colids[3],main=paste(" cells cluster ",id," of ",colids[1],"\nn = ",n_sel,sep=''),cex.main = 2.25,xaxt = 'n' ) 
#     points(cluster.sel[,2],cluster.sel[,3],pch = 16, cex=0.5,col=color)
#     points(center[center$cluster == id,2:3],pch = 16, cex=1,col='black')
       
    dx <- cluster.sel[,2] - as.numeric(center[center$cluster == id,2])
    dy <- cluster.sel[,3] - as.numeric(center[center$cluster == id,3])
    d <- sqrt(dx**2 + dy**2)
    
    flag.dist <- d > d.filter.sel
    flag.dist.cl <- rownames(cluster)  %in% rownames(cluster.sel[flag.dist,])
    
    plot(cluster[,2],cluster[,3],pch = 16, type='p',col='grey',cex=0.5,xlab=colids[2],ylab=colids[3],main=paste(" cells cluster ",id," of ",colids[1],"\nn = ",n_sel,sep=''),cex.main = 2.25,xaxt = 'n' ) 
    points(cluster.sel[,2],cluster.sel[,3],pch = 16, cex=0.5,col='red')
    points(cluster.sel[!flag.dist,2],cluster.sel[!flag.dist,3],pch = 16, cex=0.5,col='blue') 
    points(center[center$cluster == id,2:3],pch = 16, cex=2.5,col='black')
    
    #return(paste("cluster ",id," ok",sep='') )
    return(flag.dist.cl)
    
}

res.d.quantile <- list()
par(mfrow=c(3,2))
options(repr.plot.height=10,repr.plot.width=15)
for(i in levels(dr_umap.df$cluster)){
  res.d.quantile[[i]] <- dotDist(cluster = dr_umap.df[,c('cluster','UMAP_1','UMAP_2')], id = i,center = centers)
}

# par(mfrow=c(2,3))
# options(repr.plot.height=10,repr.plot.width=15)
# for(i in c('12','9','13','15','11','14') ){
#   dotDist(cluster = cluster.df.add[,c('cluster','UMAP_1','UMAP_2')], id = i,center = centers)
# }

##90th percentile with adjustment


d.filter <- sapply(X = res.d.quantile, FUN = function(x){ c(x[10]) } ) #90th percentile
names(d.filter) <- names(res.d.quantile)

d.filter.manual <- c( '0'= 8, #only filtering very far
              '1'=8, 
              '2'=6,
              '3'=6,
              '7'=as.vector(d.filter['7']),
              '8'=as.vector(d.filter['8'] ) )


par(mfrow=c(3,3))
res.flag.dist <- list()
options(repr.plot.height=15,repr.plot.width=15)
for(i in levels(dr_umap.df$cluster) ){
  res.flag.dist[[i]] <- dotClean(cluster = dr_umap.df[,c('cluster','UMAP_1','UMAP_2')], id = i,
           center = centers,d.filter = d.filter.manual)
}


flag.dist.combine <- apply(do.call(cbind,res.flag.dist),1,any) #because filter for n_cluster time and return all flag for each cluster filtering, then get any boolen result
table(flag.dist.combine) #16428

FALSE  TRUE 
16389    39 


##filter distance
all.equal (colnames(placenta.filter) , rownames(dr_umap.df) )

dr_umap.df.filterdist <- dr_umap.df[!flag.dist.combine,] #16389


##cluster distribution agian
options(repr.plot.width=10,repr.plot.height=8)
par(mfrow = c(2,3))
for (i in levels(dr_umap.df.filterdist$cluster)){

dotDistri(cluster = dr_umap.df.filterdist, id = i)

}





############rotate the umap xy coordinate if needed#############
#test <- data.frame(x=c(10,2),y=c(25,45))
#plot(test,xlim = c(-50,50),ylim = c(-50,50),cex=2 );abline(h = 0);abline(v=0)
#test.rotate90 <- data.frame(x=test$x*cos(90)-test$y*sin(90),y=test$y*cos(90)+test$x*sin(90))
#points(test.rotate90,cex=2,pch=19,col='red')

##umap.df <- data.frame(UMAP_1=x.after.sp@umap[,1],UMAP_2=x.after.sp@umap[,2])
##rownames(umap.df) <- rownames(x.after.sp@metaData)

umap.df <- dr_umap.df.filterdist[,c('UMAP_1','UMAP_2')]


##iterative rotate the umap direction
##for(degree in seq(from = 50,300,50) ) { 
for(degree in c(275) ) { 
    rotate_d <- degree*(pi/180) #degree to radian
    cat('degree,rotate_d:',degree,' ',rotate_d)
    plot(umap.df,cex=0.2,col='grey',main= paste('degree:',degree,sep=''));
    abline(h = 0);
    abline(v=0) 
    UMAP.rotate <- data.frame(
        UMAP_1=umap.df$UMAP_1*cos(rotate_d)-umap.df$UMAP_2*sin(rotate_d),
        UMAP_2=umap.df$UMAP_1*sin(rotate_d)+umap.df$UMAP_2*cos(rotate_d),
        row.names=rownames(umap.df)
    )
    #UMAP.rotate <- data.frame(UMAP_1=-1*umap.df$UMAP_2,UMAP_2=umap.df$UMAP_1)
    points(UMAP.rotate,pch=19,col='red',cex=0.2)
}

#fix degree = 275

all.equal(rownames(UMAP.rotate),rownames(dr_umap.df.filterdist))#TRUE

dr_umap.df.filterdist.rotate <- dr_umap.df.filterdist

dr_umap.df.filterdist.rotate[,c('UMAP_1','UMAP_2')] <- UMAP.rotate



table(rownames(dr_umap.df.filterdist.rotate) %in% colnames(placenta.filter) )
 TRUE 
16389

placenta.filter.rotate <- subset(placenta.filter, cells =  rownames(dr_umap.df.filterdist.rotate))

all.equal(rownames(dr_umap.df.filterdist.rotate) ,colnames(placenta.filter.rotate)  ) #TRUE

all.equal(Idents(placenta.filter.rotate),dr_umap.df.filterdist.rotate$cluster,check.attributes = FALSE ) #TRUE

all.equal(rownames(Embeddings(placenta.filter.rotate,reduction = 'umap')), rownames(dr_umap.df.filterdist.rotate) ) #TRUE

placenta.filter.rotate[["umaprotate"]] <- CreateDimReducObject(embeddings = as.matrix(dr_umap.df.filterdist.rotate[,c('UMAP_1','UMAP_2')]), key = "umaprotate_",assay = 'RNA')


###plot and save filtered and rotated object


options(repr.plot.height=7.5,repr.plot.width=7.5)
DimPlot(placenta.filter.rotate, reduction = "umaprotate",label=TRUE,cols = color_good,label.size = 8,pt.size = 1)


saveRDS(placenta.filter.rotate,'Cnot8_ESC_LC48h_WT_KO_tdTomato.filter.rotate.rds')



placenta <- placenta.filter.rotate
rm(placenta.filter)
rm(placenta.filter.rotate)
gc()




# ##########recalculate with harmony remove batch effect###
# V <- placenta@reductions$pca@cell.embeddings
# #11560 x 50
# #15107 x 50

# meta_data <- placenta@meta.data
# #11560 x 6
# #15107 x 6

# #meta_data$sample = ifelse(grepl(pattern = "-1$", x=rownames(meta_data) ), 'D1', 'D2'  )
# #meta_data$sample <- factor(meta_data$sample,levels=c('D1','D2'))
# table(meta_data$sample)

#     U     D 
# 11474  5683

# #  D1   D2 
# #6996 4564 

# #  D1   D2 
# #9461 5646 

# V_harmony <- HarmonyMatrix(V,meta_data, 'sample', do_pca = FALSE)
# # Harmony 1/10

# # Harmony 2/10

# # Harmony 3/10

# # Harmony converged after 3 iterations

# #Harmony converged after 6 iterations
# #11560 x 50

# #do not converged after 10 iterations
# #15107 x 50

# placenta@reductions$pca@cell.embeddings <- V_harmony #overwrite the PCA


# #DefaultAssay(placenta) <- 'peaks'

# #placenta <- RunUMAP(object = placenta, reduction = 'pca', dims = 1:30)
# placenta = RunUMAP(placenta,dims=1:40,seed.use = 123,min.dist = 0.5,spread=1.2) 

# #placenta = RunUMAP(placenta,dims=1:40,seed.use = 123) 
# placenta =  RunTSNE(placenta,dims=1:40,seed.use = 123)



# placenta <- FindNeighbors(object = placenta, reduction = 'pca', dims = 1:40)
# placenta <- FindClusters(placenta, resolution = 0.5) 


# # placenta<- FindClusters(object = placenta, verbose = FALSE, algorithm = 1,resolution=1.2) #ori louvain
# # placenta<- FindClusters(object = placenta, verbose = FALSE, algorithm = 2,resolution=1.2) #enhanced louvain
# # placenta<- FindClusters(object = placenta, verbose = FALSE, algorithm = 4,resolution=1.2) #leiden

# options(repr.plot.height=7.5,repr.plot.width=7.5)
# DimPlot(object = placenta, label = TRUE,cols=color_good,pt.size = 0.5,label.size = 8,reduction = "umap") + NoLegend()
# DimPlot(object = placenta, label = TRUE,cols=color_good,pt.size = 0.5,label.size = 8,reduction = "tsne") + NoLegend()

# # ##quality check for each clusters
# # options(repr.plot.height=7.5,repr.plot.width=12.5)
# # VlnPlot(placenta,features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),pt.size=0.1 )


# # ####annotation
# # options(repr.plot.height=15,repr.plot.width=15)
# # FeaturePlot(placenta, features = c("CGA", "PSG1","DNMT1","PAGE4",'MKI67','ERVFRD-1','LAIR2','PLAC8',"VIM"), reduction = "umap")

# # FeaturePlot(placenta, features = c('CSH1','CSH2',"CSHL1",'PAPPA','GH2','FLT1','INSIG1', "LEP"), reduction = "umap")


# # # FeaturePlot(placenta, features = c('CSH1','CSH2',"CSHL1",'PAPPA','GH2','FLT1','INSIG1', "PRL"), reduction = "umap")


# FeaturePlot(placenta, features = c('nCount_RNA','nFeature_RNA','percent.mt'), reduction = "umap",cols = viridis(12,option = 'D') )#c('black','red') )

# options(repr.plot.height=7.5,repr.plot.width=7.5)
# FeaturePlot(placenta, features = c('nCount_RNA'), reduction = "umap",cols = viridis(12,option = 'D') )
# FeaturePlot(placenta, features = c('nFeature_RNA'), reduction = "umap",cols = viridis(12,option = 'D') )
# FeaturePlot(placenta, features = c('percent.mt'), reduction = "umap",cols = viridis(12,option = 'D') )

# FeaturePlot(placenta, features = c('nCount_RNA'), reduction = "tsne",cols = viridis(12,option = 'D') )
# FeaturePlot(placenta, features = c('nFeature_RNA'), reduction = "tsne",cols = viridis(12,option = 'D') )
# FeaturePlot(placenta, features = c('percent.mt'), reduction = "tsne",cols = viridis(12,option = 'D') )

# FeaturePlot(placenta, features = c('sample'), reduction = "UMAP" )









#####cell cycle scoring#####
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

s.genes <- stringr::str_to_title(s.genes)
g2m.genes <- stringr::str_to_title( g2m.genes)

table(s.genes %in%  rownames(placenta)) 
FALSE  TRUE 
    1    42

s.genes <- s.genes[s.genes %in%  rownames(placenta)]

table(g2m.genes %in%  rownames(placenta)) 
FALSE  TRUE 
    2    52

g2m.genes <- g2m.genes[g2m.genes %in%  rownames(placenta)]

placenta <- CellCycleScoring(placenta, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)




####annotate ko/wt by cluster 


options(repr.plot.height=7.5,repr.plot.width=7.5)
DimPlot(placenta, reduction = "umaprotate",label=TRUE,cols = color_good,label.size = 8,pt.size = 1)

genotype <- Idents(placenta)
table(genotype)

   0    1    2    3    7    8 
5203 3816 3661 3460  163   86 

genotype  <- ifelse( genotype %in% c('2','3'), 'WT','KO'  )
table(genotype)

KO   WT 
9268 7121 




placenta <- AddMetaData(object = placenta,metadata = genotype, col.name = 'genotype')




####plot with customized method####

# ##add sample metadata
##placenta <- AddMetaData(placenta,metadata = meta_data$sample, col.name = 'sample')
#placenta <- AddMetaData(placenta,metadata = rownames(placenta@meta.data), col.name = 'cellid') #for subset object by cell id



#####get cluster.df
cluster <- Idents(placenta)
umap <- placenta@reductions$umaprotate@cell.embeddings
##umap <- placenta@reductions$umap@cell.embeddings
all.equal(names(cluster),rownames(umap)) #TRUE


colnames(umap) <- c('UMAP_1','UMAP_2')

cluster.df <- data.frame(cluster=cluster,umap)

metadata <- placenta@meta.data
all.equal(rownames(cluster.df),rownames(metadata))#TRUE

cluster.df.add <- cbind(cluster.df, metadata)


anno <- paste0(cluster.df.add$sample,'-',cluster.df.add$genotype)
table(anno)
 ESC-KO    ESC-WT LC_48h-KO LC_48h-WT 
 5244      3661      4024      3460

cluster.df.add$anno <- factor(anno,levels = c('ESC-WT','ESC-KO', 'LC_48h-WT', 'LC_48h-KO' ))


placenta <- AddMetaData(object = placenta,metadata = anno, col.name = 'anno')

##colnames(cluster.df.add) 

###filter low depth cells####
table(cluster.df.add$cluster)
   0    1    2    3    7    8 
5203 3816 3661 3460  163   86 

 0    1    2    3    7    8 
5204 3822 3663 3461  182   96 

   0    1    2    3    4    5    6    7    8 
5204 3822 3663 3461  170  244  200  182   96

#    0    1    2    3    4    5    6    7    8    9 
# 4744 3974 3718 3507  589  430  211  210  185   97

#    0    1    2    3    4    5    6    7    8    9   10   11   12   13   14 
# 3220 2566 2341 1732 1531 1490 1351 1314  570  409  262  134  111  103   23

#    0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15 
# 3390 2626 2234 1786 1194 1152  658  654  567  423  323  300  240  144  101   91 
#   16   17 
#   83   79 

#cluster.df.add.filter <- subset(cluster.df.add,! cluster %in% c('10','17')  )
#cluster.df.add <- cluster.df.add.filter

#cluster.df.add$cluster <- droplevels(cluster.df.add$cluster)
#table(cluster.df.add$cluster)
#  0    1    2    3    4    5    6    7    8    9   11   12   13   14   15   16 
#1720 1428 1225 1093 1083  963  893  853  656  365  288  189  173  123   80   74

#11206 cells 

table(cluster.df.add$sample)
  ESC LC_48h 
  8905   7484 

table(cluster.df.add$genotype)
 KO   WT 
9268 7121 

table(cluster.df.add$anno)
 ESC-WT    ESC-KO LC_48h-WT LC_48h-KO 
     3661      5244      3460      4024


saveRDS(cluster.df.add,'cluster.df.add.rds') #use umaprotate
write.table(file='cluster.df.add.txt',x=cluster.df.add,quote = FALSE,sep = '\t',row.names = TRUE,col.names = TRUE)


#####do we need to modify the placenta seurat object?yes####

#placenta@meta.data$seurat_clusters %in% levels(cluster.df.add$cluster)
#WhichCells(placenta,ident=levels(cluster.df.add$cluster) )

###rename
# placenta  <- RenameIdents(placenta,'0'='1', '1'='2','2'='3','3'='4','4'='5',
#            '5'='6','6'='7','7'='8','8'='9','9'='10',
#            '11'='11','12'='12','13'='13','14'='14','15'='15','16'='16' )


###


# placenta@meta.data
# table(colnames(placenta) %in% rownames(cluster.df.add))
# # FALSE  TRUE 
# #   354 11206 

#placenta_filter <- subset(placenta, subset = cellid %in% rownames(cluster.df.add) )
#11206

# placenta_filter <- subset(placenta, cells = rownames(cluster.df.add) )
# dim(placenta_filter)
# #24307 11206


# all.equal (colnames(placenta_filter), rownames(cluster.df.add) )#TRUE
# #all.equal (colnames(placenta_filter), colnames(placenta_filter1) )#TRUE

# #all.equal(placenta_filter,placenta_filter1) #TRUE
# #rm(placenta_filter1)


# ###rename

# table(Idents(placenta_filter))
#  0    1    2    3    4    5    6    7    8    9   11   12   13   14   15   16 
# 1720 1428 1225 1093 1083  963  893  853  656  365  288  189  173  123   80   74



# placenta_filter  <- RenameIdents(placenta_filter,'0'='1', '1'='2','2'='3','3'='4','4'='5',
#            '5'='6','6'='7','7'='8','8'='9','9'='10',
#            '11'='11','12'='12','13'='13','14'='14','15'='15','16'='16' )

# table(Idents(placenta_filter))
#   1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16 
# 1720 1428 1225 1093 1083  963  893  853  656  365  288  189  173  123   80   74 

# ###



########save object###
#saveRDS(placenta, "Cnot8_ESC_LC48h_WT_KO_tdTomato.rds")
#saveRDS(placenta_filter, "limb.final.rds")


# #####get cluster.df again and save
# cluster <- Idents(placenta_filter)
# umap <- placenta_filter@reductions$umap@cell.embeddings
# all.equal(names(cluster),rownames(umap)) #TRUE


# cluster.df <- data.frame(cluster=cluster,umap)

# metadata <- placenta_filter@meta.data
# all.equal(rownames(cluster.df),rownames(metadata))#TRUE

# cluster.df.add <- cbind(cluster.df, metadata)

# saveRDS(cluster.df.add,'cluster.df.add.rds')
# write.table(file='cluster.df.add.txt',x=cluster.df.add,quote = FALSE,sep = '\t',row.names = TRUE,col.names = TRUE)

# table(cluster.df.add$cluster)
#    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16 
# 1720 1428 1225 1093 1083  963  893  853  656  365  288  189  173  123   80   74 

####find out which clusters to filter snRNA donor1 low quality cells

# cluster.donor1 <- read.table("../../placenta_10X/03.seurate/limbSEQ-1/placenta.limb-1.umap.cl.txt",header = TRUE,sep = '\t',row.names = 1)

# cluster.donor1$cluster <- factor(cluster.donor1$cluster,levels=sort(unique(cluster.donor1$cluster)) )

# cluster.donor1.cents = cluster.donor1 %>% dplyr::group_by(cluster) %>%
#                    dplyr::summarise (meanx = mean(UMAP_1),
# 									 meany=mean(UMAP_2))


# options(repr.plot.width=7,repr.plot.height=7)
# ggplot(data=cluster.donor1,aes(x=UMAP_1,y=UMAP_2,col=cluster) ) +
#   geom_point(size=1) +
#   scale_color_manual(values = color_good) +
#   #xlim(-10,10) +
#   #ylim(-10,10) +
#   annotate(geom = 'text',x=cluster.donor1.cents$meanx,
# 		   y=cluster.donor1.cents$meany,
# 		  label=cluster.donor1.cents$cluster,
# 		  size =10,color='black') +
#   theme_classic() +
#   guides( col = guide_legend(override.aes = list(size = 6)) )

# ##filter cluster 1 

# id1 <- rownames(subset( cluster.donor1, cluster == '1'  ))

# cluster.df.add.filter <- cluster.df.add[!rownames(cluster.df.add) %in% id1,]
# #15107

# saveRDS(cluster.df.add.filter,'cluster.df.add.filter.rds')

# cluster.df.add<- cluster.df.add.filter
# rm(cluster.df.add.filter)



# sample <- 'Cnot8'



###############customized way to plot umap-cluster with text halo##########
centers <- cluster.df.add %>% dplyr::group_by(cluster) %>% dplyr::summarize(x = median(x = UMAP_1), 
        y = median(x = UMAP_2))

centers_shift = data.frame()
##add little shift for text x and y, to plot text halo, borrow from snapATAC
theta= seq(0, 2*pi, length.out=50)
r=0.1
strwidth = 0.5 
strheight = 0.5
xo <- r*strwidth # r*strwidth('A')
yo <- r*strheight #r*strheight('A')
for (i in seq_len(nrow(centers))){
  for (j in theta) {
        centers_shift = rbind(centers_shift,
                              data.frame(
                                  cluster=as.character(unlist(centers[i,'cluster'])),
                                  x=centers[i,'x'] + cos(j)*xo, 
                                  y=centers[i,'y'] + sin(j)*yo
                                 )
                       )
      }
}

################the UMAP plot with annotation##############
#label right
options(repr.plot.height=5,repr.plot.width=5.5)
ggplot(cluster.df.add,aes(x=UMAP_1,y=UMAP_2,col=cluster  )) +
  geom_point(size = 0.1,show.legend = TRUE,alpha= 1 ) +
  scale_colour_manual(values = color_good)  +
  #scale_colour_gradientn(colors = viridis(6,option = 'C'))  +
  #theme_classic() +
  #theme_bw() +
  theme(
        legend.position = 'right',
        axis.text=element_blank(), 
        axis.title = element_text(size = 15, face = "bold"),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(color="black", fill = NA,size=1),
#         panel.background = element_rect(fill = "white", colour = "white", 
#                 size = rel(1)),
        #panel.border = element_blank(),
        plot.title = element_text(size = 15, face = "bold"),
        #complete = TRUE
        plot.margin = unit(c(1,1,1,1), "lines") #add margin benifit the outside box frame
       ) +
 #theme(legend.position = 'none',axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 14, face = "bold"),axis.title.y = element_text(size = 14, face = "bold")  ) +
  ggtitle(paste(sample, "two stages, total cells:",nrow(cluster.df.add),  sep=" ") ) +
#   geom_text(data = centers_shift, #the halo
#             mapping = aes(x=x,y=y,label = cluster), 
#             colour = "white", 
#             size = 6.5) +
#   geom_text(data = centers, 
#             mapping = aes(x=x,y=y,label = cluster), 
#             colour = "black", 
#             size = 6) +
  guides(col = guide_legend(override.aes = list(size = 3))) +  ##no effect ??
  #ylim(0,1) + xlim(0,6) +  ##will overwrite scale_x/y_continuous
  #labs(x = "Fragments per barcode", y = "proportion of promoter reads")
  labs(x = "UMAP1", y = "UMAP2")
  #labs(x = "log10(pass_filters_reads)", y = "proportion of enhancer reads")

ggsave(filename = "pdfs/Cnot8-UMAP.pdf",height=5,width=5.5,useDingbats=FALSE)

##label on cluster
options(repr.plot.height=5,repr.plot.width=5)
ggplot(cluster.df.add,aes(x=UMAP_1,y=UMAP_2,col=cluster  )) +
  geom_point(size = 0.1,show.legend = TRUE,alpha= 1 ) +
  scale_colour_manual(values = color_good)  +
  #scale_colour_gradientn(colors = viridis(6,option = 'C'))  +
  #theme_classic() +
  #theme_bw() +
  theme(
        legend.position = 'none',
        axis.text=element_blank(), 
        axis.title = element_text(size = 15, face = "bold"),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(color="black", fill = NA,size=1),
#         panel.background = element_rect(fill = "white", colour = "white", 
#                 size = rel(1)),
        #panel.border = element_blank(),
        plot.title = element_text(size = 15, face = "bold"),
        #complete = TRUE
        plot.margin = unit(c(1,1,1,1), "lines") #add margin benifit the outside box frame
       ) +
 #theme(legend.position = 'none',axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 14, face = "bold"),axis.title.y = element_text(size = 14, face = "bold")  ) +
  ggtitle(paste(sample, "two stages, total cells:",nrow(cluster.df.add),  sep=" ") ) +
  geom_text(data = centers_shift, #the halo
            mapping = aes(x=x,y=y,label = cluster), 
            colour = "white", 
            size = 6.5) +
  geom_text(data = centers, 
            mapping = aes(x=x,y=y,label = cluster), 
            colour = "black", 
            size = 6) +
  ##guides(col = guide_legend(override.aes = list(size = 6))) +  ##no effect ??
  #ylim(0,1) + xlim(0,6) +  ##will overwrite scale_x/y_continuous
  #labs(x = "Fragments per barcode", y = "proportion of promoter reads")
  labs(x = "UMAP1", y = "UMAP2")
  #labs(x = "log10(pass_filters_reads)", y = "proportion of enhancer reads")

ggsave(filename = "pdfs/Cnot8-UMAP.labelon.pdf",height=5,width=5)
###############


##by donors
options(repr.plot.height=5.5,repr.plot.width=5)
#ggplot(cluster.df.add,aes(x=UMAP_1,y=UMAP_2,col= sample  )) +
ggplot(cluster.df.add,aes(x=UMAP_1,y=UMAP_2,col= anno  )) +
  geom_point(size = .1,show.legend = TRUE,alpha= 0.35 ) +
  #scale_colour_manual(values = c('red','navy'))  +
  scale_colour_manual(values = c('pink','red','lightblue','navy'))  +
  #scale_colour_gradientn(colors = viridis(6,option = 'C'))  +
  #theme_classic() +
  #theme(legend.position = 'top',axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 14, face = "bold"),axis.title.y = element_text(size = 14, face = "bold")  ) +
  theme(legend.position = 'top',
        axis.text=element_blank(), 
        axis.title = element_text(size = 15, face = "bold"),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(color="black", fill = NA,size=1),
        plot.title = element_text(size = 15, face = "bold"),
        #complete = TRUE
        plot.margin = unit(c(1,1,1,1), "lines")
       )+
  ggtitle(paste('limb', "Source of sample ",  sep=" ") ) +
  #ggtitle(paste(sample, "nFeature_RNA",  sep=" ") ) +
  guides(col = guide_legend(override.aes = list(size = 5))) +  ##no effect ??
  #ylim(0,1) + xlim(0,6) +  ##will overwrite scale_x/y_continuous
  #labs(x = "Fragments per barcode", y = "proportion of promoter reads")
  labs(x = "UMAP1", y = "UMAP2")
  #labs(x = "log10(pass_filters_reads)", y = "proportion of enhancer reads")
##ggsave(filename = "pdfs/Cnot8-source-or-stages.pdf",height=5.5,width=5)
ggsave(filename = "pdfs/Cnot8-source-or-stages_and_genotype.pdf",height=5.5,width=5)
####


###by depth
options(repr.plot.height=5.5,repr.plot.width=5)
ggplot(cluster.df.add,aes(x=UMAP_1,y=UMAP_2,col= log10(nCount_RNA) )) +
  geom_point(size = .1,show.legend = TRUE,alpha= .5 ) +
  #scale_colour_manual(values = c('red','navy'))  +
  ##scale_colour_gradientn(colors = viridis(6,option = 'C'))  +
  scale_colour_gradientn(colors = rev(color_cellranger) )  +
  #theme_classic() +
  #theme(legend.position = 'top',axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 14, face = "bold"),axis.title.y = element_text(size = 14, face = "bold")  ) +
  theme(legend.position = 'top',
        axis.text=element_blank(), 
        axis.title = element_text(size = 15, face = "bold"),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(color="black", fill = NA,size=1),
        plot.title = element_text(size = 15, face = "bold"),
        #complete = TRUE
        plot.margin = unit(c(1,1,1,1), "lines")
       )+
  ggtitle(paste(sample, "Sequence depth",  sep=" ") ) +
  #ggtitle(paste(sample, "nFeature_RNA",  sep=" ") ) +
  #guides(col = guide_legend(override.aes = list(size = 8))) +  ##no effect ??
  #ylim(0,1) + xlim(0,6) +  ##will overwrite scale_x/y_continuous
  #labs(x = "Fragments per barcode", y = "proportion of promoter reads")
  labs(x = "UMAP1", y = "UMAP2")
  #labs(x = "log10(pass_filters_reads)", y = "proportion of enhancer reads")
ggsave(filename = "pdfs/Cnot8-depth.pdf",height=5.5,width=5)

#by gene captured
options(repr.plot.height=5.5,repr.plot.width=5)
ggplot(cluster.df.add,aes(x=UMAP_1,y=UMAP_2,col= nFeature_RNA )) +
  geom_point(size = 0.5,show.legend = TRUE,alpha= 1 ) +
  #scale_colour_manual(values = c('red','navy'))  +
  scale_colour_gradientn(colors = viridis(6,option = 'C'))  +
  #theme_classic() +
  #theme(legend.position = 'top',axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 14, face = "bold"),axis.title.y = element_text(size = 14, face = "bold")  ) +
  theme(legend.position = 'top',
        axis.text=element_blank(), 
        axis.title = element_text(size = 15, face = "bold"),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(color="black", fill = NA,size=1),
        plot.title = element_text(size = 15, face = "bold"),
        #complete = TRUE
        plot.margin = unit(c(1,1,1,1), "lines")
       )+
  ggtitle(paste(sample, "Gene captured",  sep=" ") ) +
  #ggtitle(paste(sample, "nFeature_RNA",  sep=" ") ) +
  #guides(col = guide_legend(override.aes = list(size = 8))) +  ##no effect ??
  #ylim(0,1) + xlim(0,6) +  ##will overwrite scale_x/y_continuous
  #labs(x = "Fragments per barcode", y = "proportion of promoter reads")
  labs(x = "UMAP1", y = "UMAP2")
  #labs(x = "log10(pass_filters_reads)", y = "proportion of enhancer reads")
ggsave(filename = "pdfs/Cnot8-gene_captured.pdf",height=5.5,width=5)


#by mt perctage
options(repr.plot.height=5.5,repr.plot.width=5)
ggplot(cluster.df.add,aes(x=UMAP_1,y=UMAP_2,col= percent.mt )) +
  geom_point(size = 0.5,show.legend = TRUE,alpha= 1 ) +
  #scale_colour_manual(values = c('red','navy'))  +
  scale_colour_gradientn(colors = viridis(6,option = 'C'))  +
  #theme_classic() +
  #theme(legend.position = 'top',axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 14, face = "bold"),axis.title.y = element_text(size = 14, face = "bold")  ) +
  theme(legend.position = 'top',
        axis.text=element_blank(), 
        axis.title = element_text(size = 15, face = "bold"),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(color="black", fill = NA,size=1),
        plot.title = element_text(size = 15, face = "bold"),
        #complete = TRUE
        plot.margin = unit(c(1,1,1,1), "lines")
       )+
  ggtitle(paste(sample, "Percent of Mt Gene",  sep=" ") ) +
  #ggtitle(paste(sample, "nFeature_RNA",  sep=" ") ) +
  #guides(col = guide_legend(override.aes = list(size = 8))) +  ##no effect ??
  #ylim(0,1) + xlim(0,6) +  ##will overwrite scale_x/y_continuous
  #labs(x = "Fragments per barcode", y = "proportion of promoter reads")
  labs(x = "UMAP1", y = "UMAP2")
  #labs(x = "log10(pass_filters_reads)", y = "proportion of enhancer reads")
ggsave(filename = "pdfs/Cnot8-percent.mt.pdf",height=5.5,width=5)



##by cell cycling stage
options(repr.plot.height=5.5,repr.plot.width=5)
ggplot(cluster.df.add,aes(x=UMAP_1,y=UMAP_2,col= Phase  )) +
  geom_point(size = 0.3,show.legend = TRUE,alpha= 0.6 ) +
  scale_colour_manual(values = c('#282F76','brown','darkgreen'))  +
  #scale_colour_manual(values = c('red','navy','orange'))  +
  #scale_colour_gradientn(colors = viridis(6,option = 'C'))  +
  #theme_classic() +
  #theme(legend.position = 'top',axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 14, face = "bold"),axis.title.y = element_text(size = 14, face = "bold")  ) +
  theme(legend.position = 'top',
        axis.text=element_blank(), 
        axis.title = element_text(size = 15, face = "bold"),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(color="black", fill = NA,size=1),
        plot.title = element_text(size = 15, face = "bold"),
        #complete = TRUE
        plot.margin = unit(c(1,1,1,1), "lines")
       )+
  ggtitle(paste(sample, "cell cycling phase",  sep=" ") ) +
  guides(col = guide_legend(override.aes = list(size = 5))) +  ##no effect ??
  #ylim(0,1) + xlim(0,6) +  ##will overwrite scale_x/y_continuous
  #labs(x = "Fragments per barcode", y = "proportion of promoter reads")
  labs(x = "UMAP1", y = "UMAP2")
  #labs(x = "log10(pass_filters_reads)", y = "proportion of enhancer reads")
ggsave(filename = "pdfs/Cnot8-cellcycling_combined.pdf",height=5.5,width=5)


##by cell cycling stage: plot layers of Phase score
cluster.df.add.G1 <- subset(cluster.df.add,Phase == 'G1') #6246
cluster.df.add.G2M <- subset(cluster.df.add,Phase == 'G2M') #1950
cluster.df.add.S <- subset(cluster.df.add,Phase == 'S') #3010

cluster.df.add.phase <- cluster.df.add.G1
stage_name <- 'G1'
color_stage <- '#282F76'
outpdf <- "pdfs/Cnot8-cellcycling-G1.pdf"

cluster.df.add.phase <- cluster.df.add.G2M
stage_name <- 'G2M'
color_stage <- 'brown'
outpdf <- "pdfs/Cnot8-cellcycling-G2M.pdf"

cluster.df.add.phase <- cluster.df.add.S
stage_name <- 'S'
color_stage <- 'darkgreen'
outpdf <- "pdfs/Cnot8-cellcycling-S.pdf"


options(repr.plot.height=5,repr.plot.width=5)
ggplot(cluster.df.add.phase,aes(x=UMAP_1,y=UMAP_2,col= Phase  )) +
  #geom_point(size = 0.1,alpha= 0.5,col='#282F76' ) +
  ##geom_point(data=cluster.df.add.G1,aes(x=UMAP_1,y=UMAP_2,col= Phase),size = 0.1,alpha= 0.5,col='#282F76' ) +
  geom_point(data=cluster.df.add.phase,aes(x=UMAP_1,y=UMAP_2,col= Phase),size = 0.1,alpha= 0.5,col= color_stage ) +
  ##geom_point(data=cluster.df.add.S,aes(x=UMAP_1,y=UMAP_2,col= Phase),size = 0.1,alpha= 0.5,col='darkgreen' ) +
  #scale_colour_manual(values = c('#282F76','brown','darkgreen'))  +
  #scale_colour_manual(values = c('red','navy','orange'))  +
  #scale_colour_gradientn(colors = viridis(6,option = 'C'))  +
  #theme_classic() +
  #theme(legend.position = 'top',axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 14, face = "bold"),axis.title.y = element_text(size = 14, face = "bold")  ) +
  theme(legend.position = 'top',
        axis.text=element_blank(), 
        axis.title = element_text(size = 15, face = "bold"),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(color="black", fill = NA,size=1),
        plot.title = element_text(size = 15, face = "bold"),
        #complete = TRUE
        plot.margin = unit(c(1,1,1,1), "lines")
       )+
  ggtitle(paste(sample, "cell cycling phase ", stage_name , sep=" ") ) +
  guides(col = guide_legend(override.aes = list(size = 5))) +  ##no effect ??
  #ylim(0,1) + xlim(0,6) +  ##will overwrite scale_x/y_continuous
  #labs(x = "Fragments per barcode", y = "proportion of promoter reads")
  labs(x = "UMAP1", y = "UMAP2")
  #labs(x = "log10(pass_filters_reads)", y = "proportion of enhancer reads")
ggsave(filename = outpdf,height=5.5,width=5)

######plot cluster count barplot, use the same colors as UMAP
# colPanel = SnapATAC:::createColorPanel(length(unique(x.after.sp@cluster)))
# color = scales::alpha(colPanel[factor(unique(x.after.sp@cluster))], 0.8)

color = color_good

options(repr.plot.width=2,repr.plot.height=7)
barplot(table(cluster.df.add$cluster),horiz = TRUE,col = color )

# ##stb percentage
# stb_sum <- sum(table(cluster.df.add$cluster)[c('1','2','3','4','6','7')]) #7322 #5777 #6779
# stb_sum/sum(table(cluster.df.add$cluster)) #11206 #15107
# #65.3% #51.6%  44.9%



#####cluster distribution
dotDistri = function (cluster = NULL, id = NULL){
    #cluster_str = paste(type,'_',id,sep = '')
    cluster.sel = cluster[ (cluster[,'cluster'] == id),]
    #cluster.sel = cluster[(cluster[,'type'] == type) & (cluster[,'cluster_snapatac'] == id),]
    n_sel = nrow(cluster.sel)
    #cat ('select for cluster ',id,' n = ',n_sel," \n")
    color = 'red'
    #color = ifelse(type == 'atac','red','navy')
    
#     if(type == 'w8'){
#       color = ifelse(test=grepl(pattern='_D1$',x=cluster.sel$sample),'red' ,'pink'   )
#     }else if (type == 'TM'){
#         color = ifelse(test=grepl(pattern='_D1$',x=cluster.sel$sample),'navy','lightblue'    ) 
#     }
    
    plot(cluster$UMAP_1,cluster$UMAP_2,pch = 16, type='p',col='grey',cex=0.5,xlab='UMAP1',ylab='UMAP2',main=paste(" cells cluster ",id,"\nn = ",n_sel,sep=''),cex.main=2,cex.axis=2,cex.lab=2 ) 
    points(cluster.sel$UMAP_1,cluster.sel$UMAP_2,pch = 16, cex=0.5,col=color)
    return(paste("cluster ",id," ok",sep='') )
}

options(repr.plot.width=8,repr.plot.height=15)
par(mfrow = c(4,2))
for (i in 0:8){

dotDistri(cluster = cluster.df.add, id = i)

}

######stat d1 d2 and cluster cell number correlation######
res.stat <- table(cluster.df.add$cluster,cluster.df.add$sample)
   
     ESC LC_48h
  0 5194      9
  1    8   3808
  2 3660      1
  3    1   3459
  7    0    163
  8   42     44

#         U    D
#   0  1734 1656
#   1  1418 1208
#   2   945 1289
#   3   794  992
#   4   564  630
#   5   701  451
#   6   455  203
#   7   264  390
#   8   421  146
#   9   251  172
#   10  158  165
#   11  149  151
#   12  127  113
#   13   64   80
#   14   53   48
#   15   46   45
#   16   49   34
#   17   48   31

#       D1  D2
#   1  992 728
#   2  827 601
#   3  712 513
#   4  669 424
#   5  645 438
#   6  529 434
#   7  473 420
#   8  544 309
#   9  406 250
#   10 279  86
#   11 185 103
#   12 111  78
#   13 124  49
#   14  77  46
#   15  53  27
#   16  32  42

#       D1  D2
#   0  992 728
#   1  827 601
#   2  712 513
#   3  669 424
#   4  645 438
#   5  529 434
#   6  473 420
#   7  544 309
#   8  406 250
#   9  279  86
#   11 185 103
#   12 111  78
#   13 124  49
#   14  77  46
#   15  53  27
#   16  32  42
res.cor <- cor(res.stat)#0.978
res.stat[,1] = -1 * res.stat[,1]

res.stat.df <- as.data.frame(res.stat)
#res.stat.df <- reshape2::melt(res.stat)
colnames(res.stat.df) <- c('cluster','sample','count')


##plot horizonal barplot  #  '#94C6DD', '#1273AE'  vs '#F3A585','#C80927'
options(repr.plot.height=5.5,repr.plot.width=4)
ggplot(res.stat.df, aes(fill=sample, y=count, x=cluster )) + 
    #geom_hline(yintercept = c(0,25,50,75),linetype='solid',size=.3,col='black') +
    #geom_hline(yintercept = c(0),linetype='solid',size=1,col='black') +
    #geom_bar(position="stack", stat="identity",alpha=1,width = 0.5) +  
    geom_bar(position=position_stack(reverse=TRUE), stat="identity",width = 0.5) +  
    #xlim(100,0) +
    #scale_x_continuous(breaks = seq(100,0,-25),labels=seq(100,0,-25)) +
    #scale_y_reverse() +
    #annotate('text',x = 16, y = 1000,label = paste('spearman cor=',cor(res.stat)[1,2]) ) +
    coord_flip() +
    #scale_fill_viridis(discrete = T,option = "E") +
    scale_fill_manual(values = c('#F3A585', '#C80927'),labels=c('ESC','LC_48h'),name='sample' ) +
    #ggtitle("cell number count") +
    labs(title = "cell number count", subtitle=paste('spearman cor=',round(res.cor[1,2],digits=2)) ) +
    theme_ipsum(base_family = 'sans') + #to avoid Arial narrow problem in ggsave
    ylab("count")

ggsave(filename="pdfs/stage1.vs.stage2.cluster.cor.pdf",width = 4, height = 5.5)









###violin plot for gene captured and UMI count, split for WT and KO and ESC,LC_48h stage

###violin by depth
options(repr.plot.height=5.5,repr.plot.width=5)
ggplot(cluster.df.add,aes(x = anno,y = log10(nCount_RNA),color=anno,fill=anno )) +
  geom_violin(size = .2,show.legend = TRUE,alpha= 0.5,  ) +
  scale_color_manual(values = c('pink','red','lightblue','navy'))  +
  scale_fill_manual(values = c('pink','red','lightblue','navy'))  +
  ##scale_colour_gradientn(colors = viridis(6,option = 'C'))  +
  ##scale_colour_gradientn(colors = rev(color_cellranger) )  +
  #theme_classic() +
  #theme(legend.position = 'top',axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 14, face = "bold"),axis.title.y = element_text(size = 14, face = "bold")  ) +
  theme(legend.position = 'top',
        #axis.text=element_blank(), 
        axis.title = element_text(size = 15, face = "bold"),
        #axis.ticks = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(color="black", fill = NA,size=1),
        plot.title = element_text(size = 15, face = "bold"),
        #complete = TRUE
        plot.margin = unit(c(1,1,1,1), "lines")
       )+
  ggtitle(paste(sample, "Sequence depth",  sep=" ") ) +
  #ggtitle(paste(sample, "nFeature_RNA",  sep=" ") ) +
  #guides(col = guide_legend(override.aes = list(size = 8))) +  ##no effect ??
  #ylim(0,1) + xlim(0,6) +  ##will overwrite scale_x/y_continuous
  #labs(x = "Fragments per barcode", y = "proportion of promoter reads")
  labs(x = "Genotype", y = "UMI count (log10)")
  #labs(x = "log10(pass_filters_reads)", y = "proportion of enhancer reads")
ggsave(filename = "pdfs/Cnot8-depth.violin.pdf",height=5.5,width=5)



###violin by depth
options(repr.plot.height=5.5,repr.plot.width=5)
ggplot(cluster.df.add,aes(x = anno,y = nFeature_RNA,color=anno,fill=anno )) +
  geom_violin(size = .2,show.legend = TRUE,alpha= 0.5,  ) +
  scale_color_manual(values = c('pink','red','lightblue','navy'))  +
  scale_fill_manual(values = c('pink','red','lightblue','navy'))  +
  ##scale_colour_gradientn(colors = viridis(6,option = 'C'))  +
  ##scale_colour_gradientn(colors = rev(color_cellranger) )  +
  #theme_classic() +
  #theme(legend.position = 'top',axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 14, face = "bold"),axis.title.y = element_text(size = 14, face = "bold")  ) +
  theme(legend.position = 'top',
        #axis.text=element_blank(), 
        axis.title = element_text(size = 15, face = "bold"),
        #axis.ticks = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(color="black", fill = NA,size=1),
        plot.title = element_text(size = 15, face = "bold"),
        #complete = TRUE
        plot.margin = unit(c(1,1,1,1), "lines")
       )+
  ggtitle(paste(sample, "Gene captured",  sep=" ") ) +
  #ggtitle(paste(sample, "nFeature_RNA",  sep=" ") ) +
  #guides(col = guide_legend(override.aes = list(size = 8))) +  ##no effect ??
  #ylim(0,1) + xlim(0,6) +  ##will overwrite scale_x/y_continuous
  #labs(x = "Fragments per barcode", y = "proportion of promoter reads")
  labs(x = "Genotype", y = "Genes per Cell ")
  #labs(x = "log10(pass_filters_reads)", y = "proportion of enhancer reads")
ggsave(filename = "pdfs/Cnot8-gene_captured.violin.pdf",height=5.5,width=5)









####save rds###
#saveRDS(placenta, "limb.final.rds")




##get aggregated matrix##
##exprMat.ave <- AverageExpression(placenta_filter, slot = 'data')[[1]]
exprMat.ave <- AverageExpression(placenta, slot = 'data')[[1]]


##zscore normalize
#exprMat.ave.z <-  AverageExpression(placenta_filter,slot = "scale.data")[[1]]#,return.seurat = TRUE)
exprMat.ave.z <-  AverageExpression(placenta,slot = "scale.data")[[1]]

saveRDS(object = exprMat.ave,'exprMat.ave.rds')
saveRDS(object = exprMat.ave.z,'exprMat.ave.z.rds')

#exprMat.ave <- readRDS('exprMat.ave.rds')
#xprMat.ave.z <- readRDS('exprMat.ave.z.rds')




#use seurat DoHeatmap
##DoHeatmap(exprMat.ave.z, features = unlist(TopFeatures(placenta_filter[["pca"]], balanced = TRUE)), size = 3, draw.lines = FALSE)


# marker.genes = c(
#     "DNMT1", "CDH1", "PAGE4",'ITGA6','MKI67','PCNA',
#     "FLT1", "CSHL1", "PSG8", "PAPPA",'GCM1','CGA','CGB8',
#     "ERVFRD-1", "LAIR2", "PLAC8",'HLA-G','MMP2',
#     'VIM','PECAM1','THY1','DLK1','HIVEP3','CD68','CD14','HLA-A','HLA-DPB1'
#   )

# marker.genes.list <- list(
    
# # 'Preadipocyte' = c('PDGFRA', 'ITGB', 'CD34', 'CD24', 'EN1', 'MYF5'),  
# # 'MSC'= c('PRRX1', 'NT5E', 'ENG'),
 
# 'MSC'= c('PRRX1', 'NT5E', 'ENG','PDGFRA', 'CD34', 'CD24', 'EN1', 'MYF5'),

# 'Dermal fibroblast'= c ('COL3A1', 'TWIST2', 'DPT', 'COL1A1', 'LUM', 'PTN', 'LTIH5', 'NFIA', 'MDK', 'TRIL', 'APCDD1', 'PAX1'),
# 'Neural crest'= c ('PMEL', 'MITF', 'PLP1', 'FABP7', 'SOX10', 'PAX3', 'POU3F1', 'EGR2', 'MAG', 'MBP', 'MPZ'),
# 'Muscle'= c ('PAX7', 'PAX3', 'MYOD', 'MYOG', 'MYF5', 'MYH1', 'MYH2', 'MYH3', 'MYH4', 'MYH7'),
# 'Blood'= c ('HBB', 'CD36', 'TFRC', 'GYPA'),
# 'Cartilage'= c ('COL10A1', 'COL2A1', 'AGC1', 'SOX9'),
# 'Tendon'= c ('TNMD', 'SCX', 'MKX', 'SIX2'),
# 'Macrophage'= c ('CD68', 'C1QB', 'C1QA', 'C1QC', 'CD48', 'HMOX1', 'CCL3', 'CD36', 'CD14', 'CCL4', 'HMOX1', 'CCL2'),
# 'Endothelial'= c ('PECAM1', 'KDR'),
# 'Skin'= c ('KRT14', 'KRT15', 'KRT17', 'KRT5', 'PDGFRA' , 'EPAM', 'LHX2', 'BMP7', 'CXCL14', 'KRT1', 'KRT4', 'CDH1', 'EPCAM', 'MEIS2', 'KRT19', 'KRT8', 'KRT18', 'KLF5'),
# 'Neuron'= c ('NKTR', 'MEG3', 'ATRX', 'CCNL1', 'NSD3', 'RBM25', 'RBM6', 'FGFR1', 'CCNL2', 'PAXBP1'),
# 'UD_diff' = c('TBX5','PITX1','TBX4')
# )


expected.down.all <- read.table('../../Quan-Cnot8-RNASEQ_all_mm9/04.analysis_featureCnt/expected.down.all',stringsAsFactors = FALSE)$V1
expected.up.all <- read.table('../../Quan-Cnot8-RNASEQ_all_mm9/04.analysis_featureCnt/expected.up.all',stringsAsFactors = FALSE)$V1

expected.down.new <- read.table('../../Quan-Cnot8-RNASEQ_all_mm9/04.analysis_featureCnt/expected.down.new',stringsAsFactors = FALSE)$V1
expected.up.new <- read.table('../../Quan-Cnot8-RNASEQ_all_mm9/04.analysis_featureCnt/expected.up.new',stringsAsFactors = FALSE)$V1



expected.Cnot8_dep_formative<- read.table('../../Quan-Cnot8-RNASEQ_all_mm9/04.analysis_featureCnt/expected.Cnot8_dependent.formative',stringsAsFactors = FALSE)$V1
expected.Cnot8_dep_naive <- read.table('../../Quan-Cnot8-RNASEQ_all_mm9/04.analysis_featureCnt/expected.Cnot8_dependent.naive',stringsAsFactors = FALSE)$V1




marker.genes.list <- list( 
    'expected.down.all' = expected.down.all,
    'expected.up.all' = expected.up.all,
     'expected.down.new' = expected.down.new,
     'expected.up.new' = expected.up.new,
    'expected.Cnot8_dep_formative' = expected.Cnot8_dep_formative,
    'expected.Cnot8_dep_naive' = expected.Cnot8_dep_naive,
    'tdTomato' = 'tdTomato'

)


marker.genes <- unlist(marker.genes.list[c('expected.down.all','expected.up.all','tdTomato')],)
table(marker.genes %in% rownames(placenta))
TRUE 
  42


# ##STR markers
# marker.genes = c(
#     "VIM", "PECAM1", "THY1",
#     "NR2F1", "DLK1", "HIVEP3", 
#     "CD68", "CD14", "HLA-A",
#     'HLA-DPA1','HLA-DPB1','MKI67'
#   );


# marker.genes = c("LINC01605", #naive STB denovo spearman correlation detected gene signature
# "KLRD1",
# "BCL6",
# "NABP1",
# "AC108062.1",
# "HLF",
# "TENT5A",
# "CBLB"
# )

# marker.genes = c(#"TNFRSF1", #naive STB gene by liger RNA D1 D2 ATAC D1 D2
# "SH3TC2",
# "EGFR",
# "STS",
# "ELMO1",
# "TNFAIP8",
# #"NAALADL",
# "PREP",
# "DYSF",
# "SMS",
# "DCP2",
# "RYBP")


# marker.genes = c("AJ009632.2",  #LEP, FLT1 stb
# "ARHGEF28",
# "NEK11",
# "PAEP",
# "FAM153C",
# "LEP",
# "TGFB1",
# "PCED1B",
# "SLC7A1",
# "MOCOS",
# "GALNT2",
# "AC100802.1",
# "GRK3",
# "ALOX5",
# "FGD4",
# "PTPRM",
# "OAF",
# "IGF2BP2",
# "FAM167A",
# "FCGR3B",
# "INHBA",
# "SH3PXD2A"
# )

# marker.genes = c( #PAPPA stb
# "LINC01483",
# "AC006378.1",
# "SLC26A7",
# "SLCO2A1",
# "CPS1",
# "KIAA1671",
# "PAPPA",
# "THSD4",
# "KIF6",
# "AC119674.1",
# "ISM1",
# "HOPX",
# "ANGPT2",
# "SLC4A4",
# "LAMA3",
# "POSTN",
# "GNG7",
# "TRIM40",
# "SGSM1",
# "PTPRQ",
# "RNF103-CHMP3",
# "CLIC5"

# )


# marker.genes = c('FLT1','PAPPA','CGA','PSG1','ATF3','ID3','TCF4','INHBA','ANGPT2')

# marker.genes = c('SUMO2','PIAS1','PIAS2','PIAS4','SENP5','SENP6','SENP2') #xiaozy

# marker.genes = c("CSH2",'PSG1','PSG2','PSG3','PSG4','PSG5','PSG6','PSG7','PSG8','PSG9')



flag <- marker.genes %in% rownames(exprMat.ave)
marker.genes <- marker.genes[flag] #42 #79 #91

marker.mat <- t(exprMat.ave[marker.genes,])
marker.mat.z <- t(exprMat.ave.z[marker.genes,])

round(quantile(marker.mat.z),digits = 2)
0%-1.0325%-0.3350%-0.0875%0.22100%2.17

#0%-1.5125%-0.3150%-0.175%0.25100%2.45

#0%-2.36 25%-0.1 50%-0.04 75%0 100%8.84

#0% -1.12 25% -0.13 50% -0.05 75% 0.03 100% 7.14

##use complexHeatmap to visualize 
##plot the marker gene matrix heatmap


# bks <- seq(0, 100, by = 100/255) #breaks
options(repr.plot.height=8,repr.plot.width=10)
hp = Heatmap(marker.mat.z, name = "marker.mat.z", 
             cluster_rows = TRUE, 
             cluster_columns = TRUE, 
             show_row_names = TRUE,
             use_raster = FALSE,
             ##col = circlize::colorRamp2(seq(-1,1,by=2/10), viridis(n = 11,option = "C")),
             column_title = "hclust: ward.D",
             clustering_method_columns = "ward.D", ##ward.D,complete
             #clustering_distance_columns  = function (m) dist(m,method="manhattan") #euclidean, manhattan,minkowski  #clustering_distance_columns  = "pearson",
             #heatmap_legend_param = list(color_bar = "continuous")
             ##right_annotation = ha#,heatmap_width=unit(8, "cm")
) 
#width = max(grobWidth(textGrob(labels))))
draw(hp, heatmap_legend_side = "left",gap = unit(0.1, "cm")) #8.895833 8.843750



saveRDS(hp,file="hp.Rds",compress = TRUE)
##add a outline box frame
#decorate_heatmap_body("mostVar", {grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 2))})
##get row order to avoid repeatly hclust
hp.roworder <- row_order(hp)
hp.colorder <- column_order(hp)

##draw hclust for columns directly
cluster.dist <- dist(t(top50.data),method = "euclidean")
plot(as.dendrogram(hclust(cluster.dist,method = "complete") ),main="hclust: complete euclidean" )
####






###########customized FeaturePlot for marker genes ####
#marker.genes <- c('TNMD','MKX','SCX','SIX2')
#marker.genes <- marker.genes.list$MSC

marker.genes <- marker.genes.list$tdTomato
marker.genes <- marker.genes.list$expected.down.all
marker.genes <- marker.genes.list$expected.up.all

res.marker <- list()
for (gene in marker.genes){
    #options(repr.plot.height=5,repr.plot.width=5.5)
    options(repr.plot.height=12,repr.plot.width=12.5)
    p<- FeaturePlot(placenta, features = gene, reduction = "umaprotate",slot = 'scale.data',cols = c("lightgrey", "red") ,pt.size=0.5,min.cutoff = 0, max.cutoff = 3.5)+ 
    ##p<- FeaturePlot(placenta.filter, features = gene, reduction = "umap",slot = 'scale.data',cols = c("lightgrey","red" ) ,pt.size = 0.5, min.cutoff = 1, max.cutoff = 2.5)+
    ##p<- FeaturePlot(placenta_D, features = gene, reduction = "umap",slot = 'data',cols = c("lightgrey", "darkred") )+ #c("lightgrey", "#ff0000", "#00ff00") c("lightgrey", "darkred")
      #scale_color_gradientn(colours = color_ga)+
      #scale_color_gradientn(colours = color_peak)+
    
      #scale_color_gradientn(colours = color_gradient_my)+
      #scale_color_gradientn(colours = color_pyreds)+
      theme(
            legend.position = 'right',
            axis.text=element_blank(), 
            axis.title = element_text(size = 15, face = "bold"),
            axis.ticks = element_blank(),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_rect(color="black", fill = NA,size=0.8),
    #         panel.background = element_rect(fill = "white", colour = "white", 
    #                 size = rel(1)),
            #panel.border = element_blank(),
            plot.title = element_text(size = 15, face = "bold"),
            #complete = TRUE
            plot.margin = unit(c(1,1,1,1), "lines") #add margin benifit the outside box frame
           ) 
    #print(p)
    ggsave(filename = paste('pdfs/marker.gene.plot/marker_gene.red_blue.',gene,'.pdf',sep=''),height=12,width=12.5,useDingbats=FALSE  )
     res.marker[[gene]] <- p
     print(p)

}


##arrange plot by patchwork
options(repr.plot.height=13,repr.plot.width=13.5)
#options(repr.plot.height=5,repr.plot.width=5.5)
# res.marker[['PSG8']] + res.marker[['DNMT1']] + res.marker[['ERVFRD-1']] + res.marker[['PECAM1']] +
# res.marker[['CSHL1']] + res.marker[['CDH1']] + res.marker[['LAIR2']] + res.marker[['VIM']] +
# res.marker[['FLT1']] + res.marker[['MKI67']] + res.marker[['PLAC8']] + res.marker[['CD14']] +
res.marker[['tdTomato']]
#plot_layout(ncol=4,nrow=3)
#print(res.marker[['FLT1']], vp=viewport(angle=-185))
ggsave(filename = 'pdfs/marker.gene.umap.pdf',height=13.5,width=18,useDingbats=FALSE) #can save grid multiple plots







###marker gene visualization and detect cluster-specific candiate marker gene


#find markers gene for each cluster
cluster.markers <- FindAllMarkers(placenta,only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#12434 total #8091 total
cluster.markers.filter <- cluster.markers[cluster.markers$p_val_adj<0.01,]
#10631 #7186 total


#saveRDS(cluster.markers,"2.DEG/cluster.markers.rds")
#saveRDS(cluster.markers.filter,"2.DEG/cluster.markers.filter.rds" )

saveRDS(cluster.markers,"DEGs/cluster.markers.rds")
saveRDS(cluster.markers.filter,"DEGs/cluster.markers.filter.rds" )


top10 <- cluster.markers %>% group_by(cluster) %>% top_n(n=10,wt=avg_logFC)
top50 <- cluster.markers %>% group_by(cluster) %>% top_n(n=50,wt=avg_logFC)

#write.table(top50,"2.DEG/top50.de.gene.uplimb.txt",col.names = TRUE,row.names = TRUE, quote = FALSE,sep='\t')
write.table(top50,"DEGs/top50.de.gene.Cnot8.txt",col.names = TRUE,row.names = TRUE, quote = FALSE,sep='\t')


for(i in levels(top50$cluster) ){
    cat(i,'\n')
    top50_sel <- subset(top50,cluster==i)$gene
    #write.table(top50_sel,paste("2.DEG/cluster_de_genes/cluster",i,".de.gene.uplimb.txt",sep=''),col.names = FALSE,row.names = FALSE, quote = FALSE,sep='\t')
    write.table(top50_sel,paste("DEGs/cluster",i,".de.gene.Cnot8.txt",sep=''),col.names = FALSE,row.names = FALSE, quote = FALSE,sep='\t')
    
}

options(repr.plot.height=15,repr.plot.width=15)
res.deheatmap.top10 <- DoHeatmap(placenta,features = top10$gene, combine = TRUE) 
ggsave("DEGs/de.top10.heatmap.pdf",height=15,width=15,useDingbats=FALSE)



# options(repr.plot.height=55,repr.plot.width=15)
# DoHeatmap(placenta,features = top10$gene) 


options(repr.plot.height=20,repr.plot.width=15)
res.deheatmap.top50 <- DoHeatmap(placenta,features = top50$gene, combine = TRUE) 
ggsave("DEGs/de.top50.heatmap.pdf",height=15,width=15,useDingbats=FALSE)

# options(repr.plot.height=65,repr.plot.width=15)
# DoHeatmap(placenta,features = top50$gene) 



#saveRDS(res.doheatmap,"2.DEG/res.doheatmap.top50.rds")







saveRDS(placenta, "Cnot8.final.rds")

# saveRDS(placenta, "Cnot8.final.lowfilter.rds")

#placenta <- readRDS( "Cnot8.final.rds")


###############scripts end###############################

























####plot dot-heatmap for sequence depth
options(repr.plot.height=7.5,repr.plot.width=7.5)
ggplot(cluster.df.add ,aes(x=UMAP_1,y=UMAP_2,col=log10(nCount_RNA))  ) +
#ggplot(cluster.df.combine ,aes(x=UMAP_1,y=UMAP_2,col=TSS.enrichment)  ) +
  geom_point(size = 1,show.legend = TRUE,alpha= 1 ) +
  scale_colour_gradientn(colors = viridis(6,option = 'C'))  +
  theme_classic() +
  theme(legend.position = 'bottom',axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 14, face = "bold"),axis.title.y = element_text(size = 14, face = "bold")  ) +
  ggtitle(paste(sample, "nCount_RNA",  sep=" ") ) +
  #guides(shape = guide_legend(override.aes = list(size = 8))) +  ##no effect ??
  #ylim(0,1) + xlim(0,6) +  ##will overwrite scale_x/y_continuous
  #labs(x = "Fragments per barcode", y = "proportion of promoter reads")
  labs(x = "UMAP_1", y = "UMAP_2")
  #labs(x = "log10(pass_filters_reads)", y = "proportion of enhancer reads")


#ggplot(cluster.df.combine ,aes(x=UMAP_1,y=UMAP_2,col=cluster)  ) +
ggplot(cluster.df.combine ,aes(x=UMAP_1,y=UMAP_2,col=sample)  ) +
  geom_point(size = 0.3,show.legend = TRUE,alpha= 0.3 ) +
  #scale_colour_gradientn(colors = viridis(6,option = 'C'))  +
  scale_color_manual(values = c('red','blue')) +
  theme_classic() +
  theme(legend.position = 'bottom',axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 14, face = "bold"),axis.title.y = element_text(size = 14, face = "bold")  ) +
  ggtitle(paste(sample, "Samples",  sep=" ") ) +
  #guides(shape = guide_legend(override.aes = list(size = 8))) +  ##no effect ??
  #ylim(0,1) + xlim(0,6) +  ##will overwrite scale_x/y_continuous
  #labs(x = "Fragments per barcode", y = "proportion of promoter reads")
  labs(x = "UMAP_1", y = "UMAP_2") +
  guides(col = guide_legend(override.aes = list(size = 8)))
  #labs(x = "log10(pass_filters_reads)", y = "proportion of enhancer reads")










####





res.tsne <- placenta@reductions$tsne@cell.embeddings[,c(1,2)]
res.umap <- placenta@reductions$umap@cell.embeddings[,c(1,2)]

cl <- Idents(placenta)

all.equal( names(cl)  , row.names(res.tsne) )  #true
res.tsne.cl <- cbind(cluster=as.numeric(as.character(cl)),res.tsne)
res.tsne.cl <- as.data.frame(res.tsne.cl)
res.tsne.cl$cluster <- as.factor(res.tsne.cl$cluster)
#all.equal(as.character(cl), as.character(res.tsne.cl$cluster) ) ##TRUE
#write.csv(res.tsne,quote = FALSE, file="placenta.FCA7196220.tsne.csv", row.names = TRUE)
write.table(res.tsne.cl,quote = FALSE, file="placenta.limb-new.tsne.cl.txt", row.names = TRUE,col.names = NA, sep='\t')


res.umap.cl <- cbind(cluster=cl,res.umap)
res.umap.cl = as.data.frame(res.umap.cl)
write.table(res.umap.cl,quote = FALSE, file="placenta.limb-new.umap.cl.txt", row.names = TRUE,col.names = NA, sep='\t')



#saveRDS(placenta,'placenta.limb-new.rds')
placenta <- readRDS('placenta.limb-new.rds')

