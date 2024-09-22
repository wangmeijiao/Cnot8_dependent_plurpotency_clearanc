
setwd("/home/mjwang/pwd/Quan-Cnot8-RNASEQ_all_mm9/04.analysis_featureCnt")

###use this scripts to generate (CPM, log2CPM; TPM, log2TPM) ave.log2TPM.combine.all.data for simple DGE scatterplot and following steps (edgeR)
##no statistic test here


data <- read.delim("out.ftcount.format.mat", header = T)
x <- data[,c(1,6:25)]

colnames(x) <- c("gene","WT.ESC_rep1", "WT.ESC_rep2", "KO.ESC_rep1", "KO.ESC_rep2", 
                 "WT.6h_rep1", "WT.6h_rep2", "KO.6h_rep1", "KO.6h_rep2",
                 "WT.12h_rep1", "WT.12h_rep2", "KO.12h_rep1","KO.12h_rep2",
                 "WT.24h_rep1", "WT.24h_rep2", "KO.24h_rep1","KO.24h_rep2",
                 "WT.48h_rep1", "WT.48h_rep2", "KO.48h_rep1","KO.48h_rep2")
row.names(x) <- x$gene
x.bk <- x


##data preprocessing
#https://nbisweden.github.io/course_rnaseq/lab_dge.html
#https://nbisweden.github.io/course_rnaseq/lab_eda.html


cnt <- x[,-1]
#1 remove genes with counts less than 1 CPM across one of the two reps (and too big than?)
par(mar=c(7,5,1,3))
boxplot(cnt,outline=FALSE,las=2, ylab=expression('Raw read counts'))
boxplot(log2(cnt+1),outline=FALSE,las=2,ylab=expression('Log'[2]~'(Read counts+1)'))

##quick look result if filer by cnt > 5
#by columns (samples)
barplot(colSums(cnt > 5),ylim = c(0,nrow(cnt)), ylab = "Number of detected genes",las = 2 )
abline(h=median(colSums(cnt > 5)))
#by row (genes)
barplot(rowSums(cnt > 5),xlab = "Genes", ylab = "Number of samples with >5 count", 
        las = 2, names.arg = "",border = "red", col = "red" )

##filter by cpm > 1 in at least two samples, get rid of empty gene lines
keep <- rowSums (edgeR::cpm(cnt) > 1) >= 2 #at lease two of 20 samples must has CPM > 1 then log2CPM > 0
table(keep)
lost <- cnt[!keep,]
cnt.filter <- cnt[keep,]
nrow(cnt.filter) #15113
nrow(lost)
cnt <- cnt.filter

boxplot(cnt,outline=FALSE,las=2, ylab=expression('Filtered read counts'))
boxplot(log2(cnt+1),las=2,ylab=expression('Log'[2]~'(Filtered read counts+1)'))


saveRDS(cnt,'cnt.rds')
write.table(cnt, col.names = NA, row.names = TRUE, quote = FALSE, sep = "\t", "count.filter.data")


#2 normalize to CPM = 10^6*N/total or TPM ()
#better to supply raw counts (with filter empty genes) to DESeq2 and edgeR, because they have internal process steps
#( No need to normalize if use edgeR, DESeq2 these packages handle the correction and transformation internally)


##CPM
total <- colSums(cnt)

#Cnot8-2-WT-ESC_rep1.accepted_hits.bam   Cnot8-4-WT-ESC_rep2.accepted_hits.bam   Cnot8-5-KO-ESC_rep1.accepted_hits.bam  
#Cnot8-8-KO-ESC_rep2.accepted_hits.bam   Cnot8-2-WT-6h_rep1.accepted_hits.bam    Cnot8-4-WT-6h_rep2.accepted_hits.bam    
#Cnot8-5-KO-6h_rep1.accepted_hits.bam    Cnot8-8-KO-6h_rep2.accepted_hits.bam    Cnot8-2-WT-12h_rep1.accepted_hits.bam  
#Cnot8-4-WT-12h_rep2.accepted_hits.bam   Cnot8-5-KO-12h_rep1.accepted_hits.bam   Cnot8-8-KO-12h_rep2.accepted_hits.bam   
#Cnot8-2-WT-24h_rep1.accepted_hits.bam   Cnot8-4-WT-24h_rep2.accepted_hits.bam   Cnot8-5-KO-24h_rep1.accepted_hits.bam   
#Cnot8-8-KO-24h_rep2.accepted_hits.bam   Cnot8-2-WT-48h_rep1.accepted_hits.bam   Cnot8-4-WT-48h_rep2.accepted_hits.bam  
#Cnot8-5-KO-48h_rep1.accepted_hits.bam   Cnot8-8-KO-48h_rep2.accepted_hits.bam
total_featureCnt <- c(38547838,37301896,35848753,35556234,38928807,37065166,
                      35795240,35600373,38852109,37145900,35313947,35560166,
                      35007442,36865910,43402617,40220587,38895217,40181532,
                      32464837,38437434) #all without unmapped
cor(total,total_featureCnt) #0.679, not use featureCnt total, use colSums(x)

##f <- 10^6/total_featureCnt
##f <- 10^6/total
##CPM <- t(t(cnt) * f)
##log2CPM <- log2(t(t(cnt+1) * f)) # count + 1
CPM <- edgeR::cpm(cnt) ##use this, normalize to colsums
log2CPM <- edgeR::cpm(cnt,prior.count = 1, log=TRUE) #use this, log is log2

#all.equal(CPM, CPM1) TRUE
#all.equal(log2CPM, log2CPM1) ##nearly the same

boxplot(CPM,outline=FALSE,las=2,ylab="CPM")
#boxplot(CPM1,outline=FALSE,las=2,ylab="CPM1") #the same with above CPM
boxplot(log2CPM,las=2,ylab=expression('Log'[2]~'(CPM+1)'))
#boxplot(log2CPM1,las=2,ylab=expression('Log'[2]~'(CPM1+1)')) #the same with above CPM

write.table(CPM, col.names = NA, row.names = TRUE, quote = FALSE, sep = "\t", "CPM.data")
write.table(log2CPM, col.names = NA, row.names = TRUE, quote = FALSE, sep = "\t", "log2CPM.data")

##TPM, normalized to gene length 
#https://nbisweden.github.io/course_rnaseq/lab_dge.html#21_cpmtpm
tpm <- function(counts,len) {##use this
  x <- counts/len
  return(t(t(x)*1e6/colSums(x)))
}

# tpm1 <- function(counts,len,total) { 
#   x <- counts/len
#   return(t(t(x)*1e6/total))
# }

all.equal(row.names(CPM),row.names(cnt)) #true
#match filtered gene length
subset <- match(row.names(CPM), as.character(data$Geneid)) 
# subset1 <- which (as.character(data$Geneid) %in% row.names(CPM) == TRUE) #similar method
# all.equal(subset,subset1) #check: TRUE
len <- data[subset,5]
all.equal(as.character(data[subset,1]), row.names(CPM) ) #check: TRUE

TPM <- tpm(cnt,len) #use this
#TPM <- tpm1(cnt,len,total_featureCnt)
log2TPM <- log2(TPM+0.001)  # TPM + 0.001
boxplot(TPM,outline=FALSE,las=2,ylab="TPM")
boxplot(log2TPM,las=2,ylab=expression('Log'[2]~'(TPM+0.001)'))

write.table(TPM, col.names = NA, row.names = TRUE, quote = FALSE, sep = "\t", "TPM.data")
write.table(log2TPM, col.names = NA, row.names = TRUE, quote = FALSE, sep = "\t", "log2TPM.data")


######output CPM log2CPM; TPM log2TPM done #############






##prepare data matric for pca correlation, DGE simple detection 
#x <- CPM
x <- TPM

data.ESC <- x[,1:4] #0h
sample.ESC <- "0h"

data.6h <- x[,5:8] #6h
sample.6h <- "6h"

data.12h <- x[,9:12] #12h
sample.12h <- "12h"

data.24h <- x[,13:16] #24h
sample.24h <- "24h"

data.48h <- x[,17:20] #48h
sample.48h <- "48h"

#data.x.ave <- log2((data[,1]+data[,2]+2)/2) #log2RPM WT ave use this
#data.y.ave <- log2((data[,3]+data[,4]+2)/2) #log2RPM KO ave use this

data.wt.ave.ESC <- (data.ESC[,1]+data.ESC[,2])/2 #CPM WT ave use this
data.ko.ave.ESC <- (data.ESC[,3]+data.ESC[,4])/2 #CPM KO ave use this

data.wt.ave.6h <- (data.6h[,1]+data.6h[,2])/2 
data.ko.ave.6h <- (data.6h[,3]+data.6h[,4])/2 

data.wt.ave.12h <- (data.12h[,1]+data.12h[,2])/2 
data.ko.ave.12h <- (data.12h[,3]+data.12h[,4])/2 

data.wt.ave.24h <- (data.24h[,1]+data.24h[,2])/2 
data.ko.ave.24h <- (data.24h[,3]+data.24h[,4])/2 

data.wt.ave.48h <- (data.48h[,1]+data.48h[,2])/2 
data.ko.ave.48h <- (data.48h[,3]+data.48h[,4])/2 

ave.combine.all <- data.frame(  wt.ave.ESC = data.wt.ave.ESC, 
                                ko.ave.ESC = data.ko.ave.ESC,
                                wt.ave.6h = data.wt.ave.6h,
                                ko.ave.6h = data.ko.ave.6h,
                                wt.ave.12h = data.wt.ave.12h,
                                ko.ave.12h = data.ko.ave.12h,
                                wt.ave.24h = data.wt.ave.24h,
                                ko.ave.24h = data.ko.ave.24h,
                                wt.ave.48h = data.wt.ave.48h,
                                ko.ave.48h = data.ko.ave.48h
                             )
ave.log2.combine.all <- log2(ave.combine.all+0.001)

#write.table(ave.combine.all, col.names = NA, row.names = TRUE, quote = FALSE, sep = "\t", "ave.CPM.combine.all.data.1")
write.table(ave.combine.all, col.names = NA, row.names = TRUE, quote = FALSE, sep = "\t", "ave.TPM.combine.all.data")
write.table(format(round(ave.log2.combine.all,digits=5),scientific = FALSE,nsmall=5), col.names = NA, row.names = TRUE, quote = FALSE, sep = "\t", "ave.log2TPM.combine.all.data")
##use this for select.up.down.mat dir analysis
#note col.names = FALSE will not add header line whie col.names = NA will add with a proper space



######output ave TPM log2TPM done #############









#data.x <- data.x.ave
#data.y <- data.y.ave
#data.combine <- data.frame(data.x=data.x,data.y=data.y)

#boxplot(data.combine, outline=FALSE)
#plot(1:length(data.x),sort(data.x),pch=20,cex=0.2)
#plot(1:length(data.y),sort(data.y),pch=20,cex=0.2)
#abline(h=seq(0,5,0.5))

##plot ecdf cumulative curve

##log2FC
log2FC.wt6vswt0 <- log2((data.wt.ave.6h+1)/(data.wt.ave.ESC+1))
log2FC.ko6vsko0 <- log2((data.ko.ave.6h+1)/(data.ko.ave.ESC+1))

log2FC.wt12vswt0 <- log2((data.wt.ave.12h+1)/(data.wt.ave.ESC+1))
log2FC.ko12vsko0 <- log2((data.ko.ave.12h+1)/(data.ko.ave.ESC+1))

log2FC.wt24vswt0 <- log2((data.wt.ave.24h+1)/(data.wt.ave.ESC+1))
log2FC.ko24vsko0 <- log2((data.ko.ave.24h+1)/(data.ko.ave.ESC+1))

log2FC.wt48vswt0 <- log2((data.wt.ave.48h+1)/(data.wt.ave.ESC+1))
log2FC.ko48vsko0 <- log2((data.ko.ave.48h+1)/(data.ko.ave.ESC+1))

log2FC.wt12vswt6 <- log2((data.wt.ave.12h+1)/(data.wt.ave.6h+1))
log2FC.ko12vsko6 <- log2((data.ko.ave.12h+1)/(data.ko.ave.6h+1))

log2FC.wt24vswt12 <- log2((data.wt.ave.24h+1)/(data.wt.ave.12h+1))
log2FC.ko24vsko12 <- log2((data.ko.ave.24h+1)/(data.ko.ave.12h+1))

log2FC.wt48vswt24 <- log2((data.wt.ave.48h+1)/(data.wt.ave.24h+1))
log2FC.ko48vsko24 <- log2((data.ko.ave.48h+1)/(data.ko.ave.24h+1))

log2FC.combine.all.stepwise <- data.frame(                                  
                                           log2FC.wt6vswt0=log2FC.wt6vswt0,
                                           log2FC.wt12vswt6=log2FC.wt12vswt6,
                                           log2FC.wt24vswt12=log2FC.wt24vswt12,
                                           log2FC.wt48vswt24=log2FC.wt48vswt24,
                                           log2FC.ko6vsko0=log2FC.ko6vsko0,
                                           log2FC.ko12vsko6=log2FC.ko12vsko6,
                                           log2FC.ko24vsko12=log2FC.ko24vsko12,
                                           log2FC.ko48vsko24=log2FC.ko48vsko24
                                           
                                          )
par(mar=c(8,3,3,3),mgp=c(2,0.5,0))
boxplot(log2FC.combine.all.stepwise,outline=FALSE,las=2);

write.table(log2FC.combine.all.stepwise, col.names = NA, row.names = TRUE, quote = FALSE, sep = "\t", "log2FC.combine.all.stepwise.data")
write.table(log2FC.combine.all.stepwise, col.names = NA, row.names = TRUE, quote = FALSE, sep = "\t", "log2FC.tpm.combine.all.stepwise.data")
#note col.names = FALSE will not add header line whie col.names = NA will add with a proper space

#log2FC.combine.all <- data.frame(
#  log2FC.wt6vswt0=log2FC.wt6vswt0,
#  log2FC.ko6vsko0=log2FC.ko6vsko0,
#  log2FC.wt12vswt0=log2FC.wt12vswt0,
#  log2FC.ko12vsko0=log2FC.ko12vsko0,
#  log2FC.wt24vswt0=log2FC.wt24vswt0,
#  log2FC.ko24vsko0=log2FC.ko24vsko0,
#  log2FC.wt48vswt0=log2FC.wt48vswt0,
#  log2FC.ko48vsko0=log2FC.ko48vsko0
  
#)
#par(mar=c(8,3,3,3),mgp=c(2,0.5,0))
#boxplot(log2FC.combine.all,outline=FALSE,las=2);

log2FC.combine.all.time <- data.frame(
  log2FC.wt6vswt0=log2FC.wt6vswt0,
  log2FC.wt12vswt0=log2FC.wt12vswt0,
  log2FC.wt24vswt0=log2FC.wt24vswt0,
  log2FC.wt48vswt0=log2FC.wt48vswt0,
  log2FC.ko6vsko0=log2FC.ko6vsko0, 
  log2FC.ko12vsko0=log2FC.ko12vsko0,
  log2FC.ko24vsko0=log2FC.ko24vsko0, 
  log2FC.ko48vsko0=log2FC.ko48vsko0
)


write.table(log2FC.combine.all.time, col.names = NA, row.names = TRUE, quote = FALSE, sep = "\t", "log2FC.combine.all.time.data")
write.table(log2FC.combine.all.time, col.names = NA, row.names = TRUE, quote = FALSE, sep = "\t", "log2FC.tpm.combine.all.time.data")
#note col.names = FALSE will not add header line whie col.names = NA will add with a proper space


par(mar=c(2,5,3,3),mgp=c(2,0.5,0))
boxplot(log2FC.wt24vswt12 , outline=FALSE)
boxplot(log2FC.ko24vsko12 , outline=FALSE)

#pdf(file = "log2FC.all-vs-0h.pdf",height = 5.86,width = 5.60 )
par(mar=c(8,5,3,3),mgp=c(3,0.5,0))
par(xpd=FALSE)
boxplot(log2FC.combine.all.time,outline=FALSE,las=2,border=c(rep("blue",4),rep("red",4)),ylab="log2FC");
abline(h=0.467,lty=2,col='grey'); #1.38 fold
abline(h=-0.224,lty=2,col='grey') #0.856 fold
text(x=1.5,y=0.5,labels = "log2FC = 0.467",cex=0.8)
text(x=1.5,y=-0.35,labels = "log2FC = -0.224",cex=0.8)
text(x=7.5,y=2,labels = paste("genes =",nrow(log2FC.combine.all.time),sep=" "), cex=0.8 )
#dev.off()

#do some adventure
log2FC.wt <- log2FC.wt6vswt0
log2FC.ko <-  log2FC.ko6vsko0
log.sample <- expression("Log"[2]~"(6h/ESC)")

log2FC.wt <- log2FC.wt12vswt0
log2FC.ko <-  log2FC.ko12vsko0
log.sample <- expression("Log"[2]~"(12h/ESC)")

log2FC.wt <- log2FC.wt24vswt0
log2FC.ko <-  log2FC.ko24vsko0
log.sample <- expression("Log"[2]~"(24h/ESC)")

log2FC.wt <- log2FC.wt48vswt0
log2FC.ko <-  log2FC.ko48vsko0
log.sample <- expression("Log"[2]~"(48h/ESC)")

log2FC.wt <- log2FC.wt12vswt6
log2FC.ko <-  log2FC.ko12vsko6
log.sample <- expression("Log"[2]~"(12h/6h)")

log2FC.wt <- log2FC.wt24vswt12
log2FC.ko <-  log2FC.ko24vsko12
log.sample <- expression("Log"[2]~"(24h/12h)")

log2FC.wt <- log2FC.wt48vswt24
log2FC.ko <-  log2FC.ko48vsko24
log.sample <- expression("Log"[2]~"(48h/24h)")

par(mar=c(3,3,5,3))
hist(log2FC.ko,breaks = 60,col='#d7134548',xlab=log.sample,main="") #transparent 50% = #80 30%=48
hist(log2FC.wt,breaks = 60,add=TRUE,col="#102b6a80");
legend(x = 2,y=1500,legend=c("KO","WT"),fill=c('#d7134548',"#102b6a")  )

log2FC.combine <- data.frame(wt=log2FC.wt,ko=log2FC.ko)
boxplot(log2FC.combine,outline=FALSE,main=log.sample);
stripchart(log2FC.combine,method="jitter",vertical=TRUE,pch=16,cex=0.5, at = c(1,2),add=TRUE)

plot(log2FC.wt,log2FC.ko,pch=20,cex=0.5,las=1,main=log.sample);abline(h=0);abline(v=0)
#library(LSD);heatscatter(log2FC.wt,log2FC.ko,colpal="bl2gr2rd",add.contour=T, ncol=8);
#library(MASS);d <- kde2d(log2FC.wt,log2FC.ko);
#ColorLevels <- round(seq(min(d$z), max(d$z), length=8),2);
#legend("bottomright", legend=ColorLevels, pch=19, col=bl2gr2rd(8))

#pdf(file = "ecdf.log2FC.48h-vs-24h.pdf",height = 4.29,width = 4.38 )
par(xpd=FALSE)
plot(ecdf(log2FC.wt),verticals=TRUE, do.points=FALSE,col='blue',lwd=1,xlim=c(-4,4), las=1,mgp=c(2,0.5,0),#yaxs="i": space y between x
      main="",cex.lab=1.4, cex.axis=1.2,   #panel.first = grid() ,main=paste(id.log2FC.wt,id.log2FC.ko,sep=' VS '),
     #xlab=expression("Log"[2]~"Fold change (48h/24h)"), ylab="Cumulative fraction"  )
     xlab=log.sample, ylab="Cumulative fraction"  )
lines(ecdf(log2FC.ko), do.points=FALSE,col='red',lwd=1,lty=1,xlim=c(-2,2))
box(lwd=2)
par(xpd=TRUE) #can draw outside the box
legend(-1,1.25,legend=c("WT","KO"), lty = c(1,1),lwd=2,col = c('blue','red')    )
#dev.off()

#another method ecdfplot
library(lattice)
library(latticeExtra)

vals <- log2FC.combine.all.time
ecdfplot(~ log2FC.wt6vswt0 + log2FC.ko6vsko0 , data=vals, prepanel = "prepanel.ecdfplot", panel = "panel.ecdfplot", ref = FALSE, lwd = 1.2)
ecdfplot(~ log2FC.wt12vswt0 + log2FC.ko12vsko0 , data=vals, prepanel = "prepanel.ecdfplot", panel = "panel.ecdfplot", ref = FALSE, lwd = 1.2)
ecdfplot(~ log2FC.wt24vswt0 + log2FC.ko24vsko0 , data=vals, prepanel = "prepanel.ecdfplot", panel = "panel.ecdfplot", ref = TRUE, lwd = 1.2)
ecdfplot(~ log2FC.wt48vswt0 + log2FC.ko48vsko0 , data=vals, prepanel = "prepanel.ecdfplot", panel = "panel.ecdfplot", ref = FALSE, lwd = 1.2)
###


##FC
FC.wt6vswt0 <- ((data.wt.ave.6h+1)/(data.wt.ave.ESC+1))
FC.ko6vsko0 <- ((data.ko.ave.6h+1)/(data.ko.ave.ESC+1))

FC.wt12vswt6 <- ((data.wt.ave.12h+1)/(data.wt.ave.6h+1))
FC.ko12vsko6 <- ((data.ko.ave.12h+1)/(data.ko.ave.6h+1))

FC.wt12vswt0 <- ((data.wt.ave.12h+1)/(data.wt.ave.ESC+1))
FC.ko12vsko0 <- ((data.ko.ave.12h+1)/(data.ko.ave.ESC+1))

FC.wt24vswt0 <- ((data.wt.ave.24h+1)/(data.wt.ave.ESC+1))
FC.ko24vsko0 <- ((data.ko.ave.24h+1)/(data.ko.ave.ESC+1))

FC.wt48vswt0 <- ((data.wt.ave.48h+1)/(data.wt.ave.ESC+1))
FC.ko48vsko0 <- ((data.ko.ave.48h+1)/(data.ko.ave.ESC+1))

FC.wt24vswt12 <- ((data.wt.ave.24h+1)/(data.wt.ave.12h+1))
FC.ko24vsko12 <- ((data.ko.ave.24h+1)/(data.ko.ave.12h+1))

FC.wt48vswt24 <- ((data.wt.ave.48h+1)/(data.wt.ave.24h+1))
FC.ko48vsko24 <- ((data.ko.ave.48h+1)/(data.ko.ave.24h+1))

FC.wt <- FC.wt6vswt0
FC.ko <-  FC.ko6vsko0
sample <- "(6h/ESC)"

FC.wt <- FC.wt12vswt6
FC.ko <-  FC.ko12vsko6
sample <- "(12h/6h)"

FC.wt <- FC.wt12vswt0
FC.ko <-  FC.ko12vsko0
sample <- "(12h/ESC)"

FC.wt <- FC.wt24vswt0
FC.ko <-  FC.ko24vsko0
sample <- "(24h/ESC)"

FC.wt <- FC.wt48vswt0
FC.ko <-  FC.ko48vsko0
sample <- "(48h/ESC)"

FC.wt <- FC.wt24vswt12
FC.ko <-  FC.ko24vsko12
sample <- "(24h/12h)"

FC.wt <- FC.wt48vswt24
FC.ko <-  FC.ko48vsko24
sample <- "(48h/24h)"

#pdf(file = "ecdf.log2FC.48h-vs-24h.pdf",height = 4.29,width = 4.38 )
par(xpd=FALSE)
plot(ecdf(FC.wt),verticals=TRUE, do.points=FALSE,col='blue',lwd=1,xlim=c(0,10), las=1,mgp=c(2,0.5,0),#yaxs="i": space y between x
     main="",cex.lab=1.4, cex.axis=1.2,   #panel.first = grid() ,main=paste(id.log2FC.wt,id.log2FC.ko,sep=' VS '),
     #xlab=expression("Log"[2]~"Fold change (48h/24h)"), ylab="Cumulative fraction"  )
     xlab=paste("Fold change",sample,sep=' '), ylab="Cumulative fraction")
lines(ecdf(FC.ko), do.points=FALSE,col='red',lwd=1,lty=1,xlim=c(-2,2))
box(lwd=2)
par(xpd=TRUE) #can draw outside the box
legend(1.5,1.3,legend=c("WT","KO"), lty = c(1,1),lwd=2,col = c('blue','red')    )
#dev.off()


##look at disperse: mean vs var
CPM.mean <- rowMeans(CPM)
CPM.var <- apply(CPM, 1, var )
plot(CPM.mean,CPM.var,cex=0.1)

log2CPM.mean <- rowMeans(log2CPM)
log2CPM.var <- apply(log2CPM, 1, var )
plot(log2CPM.mean,log2CPM.var,cex=0.1)

TPM.mean <- rowMeans(TPM)
TPM.var <- apply(TPM, 1, var )
plot(TPM.mean,TPM.var,cex=0.1)

log2TPM.mean <- rowMeans(log2TPM)
log2TPM.var <- apply(log2TPM, 1, var )
plot(log2TPM.mean,log2TPM.var,cex=0.1)


##looks correlation between reps, dist (dendgram) between samples



##plot PCA




##########################do the simple pairwise fold change visulization and output step
#height=4.166667,width = 4.177083
drawSimpleDE <- function(data  = NULL,sample = NULL,select.up = NULL, select.down = NULL, threshold.exp = NULL, threshold.FC = NULL){
  r <- format(cor(data[,1],data[,2],method = 'spearman'),digits = 3)
  #plot(data,pch=20,cex=0.65,xlim=c(-1,1000),ylim=c(-1,1000),col='lightgrey',
      #xlab=paste("WT",sample,"RPM",sep=" "),ylab=paste("KO",sample,"RPM",sep=" "))
  plot(data,pch=20,cex=0.65,xlim=c(-1,20),ylim=c(-1,20),col='lightgrey',
       xlab=paste("WT",sample,"log2TPM",sep=" "),ylab=paste("KO",sample,"log2TPM",sep=" "))
  points(data[select.up,],pch=20,cex=0.65,col="red" )
  points(data[select.down,],pch=20,cex=0.65,col="blue" )
  text(x = -1, y = 12,labels = paste("Up =",sum(select.up),sep = ""), col='black', pos = 4) #align left end
  text(x = 5.5, y = 1,labels = paste("Down =",sum(select.down),sep = ""), col='black', pos = 4) 
  text(x = -1, y = 17.5,labels = paste("R=",r,sep = ""), col='black', pos = 4)
  text(x = -1, y = 19,labels = paste("N =",nrow(data),sep = ""), col='black', pos = 4)
  text(x = -1, y = 16,labels = paste("Log2TPM >=",threshold.exp, "FC >=",format(2^threshold.FC,digits = 2),sep = " "),col='black', pos = 4)
  #text(x = 0, y = 500,labels = paste("n=",sum(select.up),sep = ""), col='black')
  #text(x = 500, y = 1,labels = paste("n=",sum(select.down),sep = ""), col='black')
  #text(x = -1, y = 18,labels = paste("R=",r,sep = ""), col='black')
  #abline(a=0.5849,b=1,col='black',lty=2,lwd=1.2,untf = FALSE)
  #abline(a=-0.5849,b=1,col='black',lty=2,lwd=1.2,untf = FALSE)
  #abline(a=1,b=1,col='darkgrey',lty=2,lwd=1.2,untf = FALSE)
  #abline(a=-1,b=1,col='darkgrey',lty=2,lwd=1.2,untf = FALSE)
  abline(a=threshold.FC,b=1,col='black',lty=2,lwd=1.2,untf = FALSE)
  abline(a=-1*threshold.FC,b=1,col='black',lty=2,lwd=1.2,untf = FALSE)
  ##abline(a=0,b=1,col='darkgrey',lty=1,lwd=1.2)
  #abline(h=0,v=0)
  
  ##add some gene point annotation
  
  expected.up.new <- read.table("expected.up.all",stringsAsFactors = F)$V1
  expected.down.new <- read.table("expected.down.all",stringsAsFactors = F)$V1
  
  expected.up.remove <- read.table("expected.up.remove",stringsAsFactors = F)$V1
  expected.down.remove <- read.table("expected.down.remove",stringsAsFactors = F)$V1
 
  expected.up.new <- expected.up.new[!(expected.up.new %in% expected.up.remove  )]
  expected.down.new <- expected.down.new[!(expected.down.new %in% expected.down.remove  )]
  
  #expected.up.new <- read.table("expected.up.new",stringsAsFactors = F)$V1
  #expected.down.new <- read.table("expected.down.new",stringsAsFactors = F)$V1
  
  expected.up.new.data  <- data[match (expected.up.new,rownames(data)), ]
  expected.down.new.data  <- data[match (expected.down.new,rownames(data)), ]
  
  points(expected.up.new.data,col="yellow",pch=1,cex=0.65)
  points(expected.down.new.data,col="darkgreen",pch=1,cex=0.65)
  
  ##add text info
  text(x=14,y=15.8,pos=4,srt=45,labels = paste("FC =",format(2^threshold.FC,digits = 2),sep = " ")) 
  text(x=15.3,y=14.3,pos=4,srt=45,labels = paste("FC =",format(-1*2^threshold.FC,digits = 2),sep = " ")) 
#   for(i in 1:nrow(expected.up.new.data)){
#     lines(c(expected.up.new.data[i,1],expected.up.new.data[i,1]+3),c(expected.up.new.data[i,2],expected.up.new.data[i,2]+1),lwd=1,col='black')
#     text(x=expected.up.new.data[i,1]+3,y=expected.up.new.data[i,2]+1,pos=4,labels=row.names(expected.up.new.data)[i],cex = 0.8 )
#   }
#   for(i in 1:nrow(expected.down.new.data)){
#     lines(c(expected.down.new.data[i,1],expected.down.new.data[i,1]-1),c(expected.down.new.data[i,2],expected.down.new.data[i,2]+3),lwd=1,col='black')
#     text(x=expected.down.new.data[i,1]-0.5,y=expected.down.new.data[i,2]+3,pos=2,labels=row.names(expected.down.new.data)[i],cex = 0.8 )
#   }
  
  return(1)
}
###########

#wt.ave.ESC ko.ave.ESC  wt.ave.6h  ko.ave.6h wt.ave.12h ko.ave.12h wt.ave.24h ko.ave.24h wt.ave.48h ko.ave.48h
ave.log2.combine.all <- log2(ave.combine.all+0.001)
samples <- c("ESC","6h","12h","24h","48h")

data.combine <- log2(ave.combine.all[,1:2]+0.001)
sample <- "ESC"

data.combine <- log2(ave.combine.all[,3:4]+0.001)
sample <- "6h"

data.combine <- log2(ave.combine.all[,5:6]+0.001)
sample <- "12h"

data.combine <- log2(ave.combine.all[,7:8]+0.001)
sample <- "24h"

data.combine <- log2(ave.combine.all[,9:10]+0.001)
sample <- "48h"

#use wt as reference, up means KO > WT, down means KO < WT

#select.up <- ( (data.combine[,2]+1)/(data.combine[,1]+1) >= 2 ) #FC >= 2, KO/WT
#select.down <- (   (data.combine[,2]+1)/(data.combine[,1]+1) <=0.5  ) #FC <= 1/2, KO/WT

#select.up <- (  data.combine[,2]-data.combine[,1] > 1 ) #log2FC > 1, 2 fold change
#select.down <- ( data.combine[,2]-data.combine[,1] < -1 ) # log2FC < -1, 2 fold change

#select.up <- ( data.combine[,1] >= 1 & data.combine[,2] >= 1  & data.combine[,2]-data.combine[,1] >= 1 ) #log2TPM >=1, log2FC >= 1, 2 fold change
#select.down <- (  data.combine[,1] >= 1 & data.combine[,2] >= 1 & data.combine[,2]-data.combine[,1] <= -1 ) ##log2TPM >=1, log2FC >= 1, 2 fold change

#select.up <- ( data.combine[,1] >= 1 & data.combine[,2] >= 1  & data.combine[,2]-data.combine[,1] >= 0.5849 ) #log2TPM >=1, log2FC >= 0.5849, 1.5 fold change
#select.down <- (  data.combine[,1] >= 1 & data.combine[,2] >= 1 & data.combine[,2]-data.combine[,1] <= -0.5849 ) ##log2TPM >=1, log2FC >= 0.5849, 1.5 fold change

log2thresh_exp <- 1 # log2TPM >= 1, TPM >= 2
#log2thresh_exp <- 2.32 # log2TPM >= 2.32, TPM >= 5
#logFC_thresh <- 0.678 # FC >= 1.6
logFC_thresh <- 0.5849 # FC >= 1.5
##logFC_thresh <- 1 # FC >= 2
#logFC_thresh <- 0.848 #FC 1.8


#15.062500  9.635417
##pdf(file = "select.up.down.mat_FC1.5/out.pairwise.scatter.all.FC1.5.pdf",height = 9.6,width = 15,useDingbats = FALSE  )
##pdf(file = "select.up.down.mat_FC1.5/out.pairwise.scatter.all.FC1.5.nolabel.pdf",height = 9.6,width = 15,useDingbats = FALSE  )

##pdf(file = "select.up.down.mat_FC1.5/out.pairwise.scatter.all.FC2.0.pdf",height = 9.6,width = 15,useDingbats = FALSE )
#pdf(file = "select.up.down.mat_FC1.5/out.pairwise.scatter.all.FC2.0.nolabel.pdf",height = 9.6,width = 15,useDingbats = FALSE )

#options(repr.plot.width=15,repr.plot.height=10)
##par(mar=c(5,5,3,3),mfrow=c(2,3))

for(i in seq(1,10,2) ){
  sample <- samples[(i+1)/2]
  data.combine <- ave.log2.combine.all[,i:(i+1)]

  select.up <- ( data.combine[,1] >= log2thresh_exp & data.combine[,2] >= log2thresh_exp  & 
                 data.combine[,2]-data.combine[,1] >= logFC_thresh) #log2TPM >=1, log2FC >= 0.5849, 1.5 fold change
  select.down <- (  data.combine[,1] >= log2thresh_exp & data.combine[,2] >= log2thresh_exp & 
                    data.combine[,2]-data.combine[,1] <= -1*logFC_thresh) ##log2TPM >=1, log2FC >= 0.5849, 1.5 fold change

#select.up <- ( data.combine[,1] >= 5 & data.combine[,2] >= 5  & data.combine[,2]-data.combine[,1] >= 1 ) #log2TPM >=1, log2FC >= 0.5849, 1.5 fold change
#select.down <- (  data.combine[,1] >= 5 & data.combine[,2] >= 5 & data.combine[,2]-data.combine[,1] <= -1 ) ##log2TPM >=1, log2FC >= 0.5849, 1.5 fold change


#select.up <- ( data.combine[,2] >= 1  & data.combine[,2]-data.combine[,1] >= 1.3219 ) #log2FPKM >=1, 2.5 fold change
#select.down <- (  data.combine[,1] >= 1 & data.combine[,2]-data.combine[,1] <= -1.3219 ) ##log2FPKM >=1, 2.5 fold change

#select.up <- ( data.combine[,2] >= 1  & data.combine[,2]-data.combine[,1] >= 2.32 ) #log2FPKM >=1, 5 fold change
#select.down <- (  data.combine[,1] >= 1 & data.combine[,2]-data.combine[,1] <= -2.32 ) ##log2FPKM >=1, 5 fold change


#select.up <- ( data.combine[,2] >= 2  & data.combine[,2]-data.combine[,1] >= 1 ) #log2FPKM >=2, 2 fold change
#select.down <- (  data.combine[,1] >= 2 & data.combine[,2]-data.combine[,1] <= -1 ) ##log2FPKM >=2, 2 fold change

#select.up <- ( data.combine[,2] >= 2  & data.combine[,2]-data.combine[,1] >= 1.3219 ) #log2FPKM >=2, 2.5 fold change
#select.down <- (  data.combine[,1] >= 2 & data.combine[,2]-data.combine[,1] <= -1.3219 ) ##log2FPKM >=2, 2.5 fold change

#select.up <- ( data.combine[,2] >= 2  & data.combine[,2]-data.combine[,1] >= 2.32 ) #log2FPKM >=2, 5 fold change
#select.down <- (  data.combine[,1] >= 2 & data.combine[,2]-data.combine[,1] <= -2.32 ) ##log2FPKM >=2, 5 fold change

  ##outfile = paste("select.up.down.mat_FC1.5/out.pairwise.scatter.all.FC1.5.",sample,".pdf",sep='')
  outfile = paste("select.up.down.mat_FC1.5/out.pairwise.scatter.all.FC1.5.nolabel.",sample,".pdf",sep='')
  pdf(file = outfile,height = 7.5,width = 7.5,useDingbats = FALSE  )
  drawSimpleDE(data.combine,sample,select.up = select.up,select.down = select.down, 
               threshold.exp = log2thresh_exp, threshold.FC = logFC_thresh)

  dev.off()
    
##output up and down regulated gene mat
  write.table(format(ave.combine.all[select.up,],digits = 3),row.names = TRUE, col.names = NA, 
            sep = "\t",quote = FALSE, file = paste("select.up.down.mat_FC1.5/select.up.TPM.",sample,".FC.",format(2^logFC_thresh,digits = 2),".mat",sep=""))
  write.table(format(ave.combine.all[select.down,],digits = 3),row.names = TRUE, col.names = NA, 
            sep = "\t",quote = FALSE, file = paste("select.up.down.mat_FC1.5/select.down.TPM.",sample,".FC.",format(2^logFC_thresh,digits = 2),".mat",sep=""))

  write.table(format(ave.log2.combine.all[select.up,],digits = 3),row.names = TRUE, col.names = NA, 
            sep = "\t",quote = FALSE, file = paste("select.up.down.mat_FC1.5/select.up.log2TPM.",sample,".FC.",format(2^logFC_thresh,digits = 2),".mat",sep=""))
  write.table(format(ave.log2.combine.all[select.down,],digits = 3),row.names = TRUE, col.names = NA, 
            sep = "\t",quote = FALSE, file = paste("select.up.down.mat_FC1.5/select.down.log2TPM.",sample,".FC.",format(2^logFC_thresh,digits = 2),".mat",sep=""))

}#draw pairwise 
#dev.off()







##do clustering and heatmap
par(mar=c(8,5,3,3))



##draw expected gene list trend
expected.up.list <- read.table("cluster3/expected.up.list",stringsAsFactors = FALSE)$V1
expected.down.list <- read.table("cluster3/expected.down.list",stringsAsFactors = FALSE)$V1

expected.up.list[6] <- "Dnmt3l" # "Dnmt3L" -> "Dnmt3l"
expected.up.list[14] <- "Pou3f1" #Oct6 -> Pou3f1
expected.down.list[1] <- "Pou5f1"  #Oct4 -> Pou5f1
expected.down.list[16] <- "Stat3" #stat3 -> Stat3

write.table(expected.up.list,sep = "",quote = F, row.names = F, col.names = F,file = "expected.up.new")
write.table(expected.down.list,sep = "",quote = F, row.names = F, col.names = F,file = "expected.down.new")

#expected.up.data  <- data.combine[match (expected.up.list,row.names(data.combine)), ]
#expected.down.data  <- data.combine[match (expected.down.list,row.names(data.combine)),]
#pheatmap::pheatmap(expected.up.data,cluster_rows = FALSE,cluster_cols = FALSE,main = paste(sample,"expected.up.data",sep = " ") )
#pheatmap::pheatmap(expected.down.data,cluster_rows = FALSE, cluster_cols = FALSE,main = paste(sample,"expected.down.data",sep = " ") )

#all together
expected.up.data  <- ave.log2.combine.all[match (expected.up.list,row.names(ave.log2.combine.all)), ]
expected.down.data  <- ave.log2.combine.all[match (expected.down.list,row.names(ave.log2.combine.all)),]

expected.up.data.z <- t(apply(expected.up.data, 1, function(x) scale(x, center = T, scale = T) ) )  #zscore
expected.down.data.z <- t(apply(expected.down.data, 1, function(x) scale(x, center = T, scale = T) ) )  
colnames(expected.up.data.z) <- colnames(expected.up.data)
colnames(expected.down.data.z) <- colnames(expected.down.data)

#pairwise
pheatmap::pheatmap(expected.up.data,fontsize_col = 16,fontsize_row = 12,cluster_rows = TRUE,
                   cluster_cols = FALSE,main = paste("expected.up",sep = " "),scale = "row",
                   color = colorRampPalette(c("navy","white","firebrick"))(20) )
pheatmap::pheatmap(expected.down.data,fontsize_col = 16,fontsize_row = 12,cluster_rows = TRUE, 
                   cluster_cols = FALSE,main = paste("expected.down",sep = " "),scale = "row",
                   color = colorRampPalette(c("navy","white","firebrick"))(20)  )
#time course
pheatmap::pheatmap(expected.up.data[,c(1,3,5,7,9,2,4,6,8,10)],fontsize_col = 16,fontsize_row = 12,cluster_rows = TRUE,
                   cluster_cols = FALSE,main = paste("expected.up",sep = " "),scale = "row",
                   color = colorRampPalette(c("navy","white","firebrick"))(20) )
pheatmap::pheatmap(expected.down.data[,c(1,3,5,7,9,2,4,6,8,10)],fontsize_col = 16,fontsize_row = 12,cluster_rows = TRUE, 
                   cluster_cols = FALSE,main = paste("expected.down",sep = " "),scale = "row",
                   color = colorRampPalette(c("navy","white","firebrick"))(20) ) 

pheatmap::pheatmap(expected.up.data.z[,c(1,3,5,7,9,2,4,6,8,10)],fontsize_col = 16,fontsize_row = 12,cluster_rows = TRUE,
                   cluster_cols = FALSE,main = paste("expected.up",sep = " "),
                   color = colorRampPalette(c("navy","white","firebrick"))(20) )
pheatmap::pheatmap(expected.down.data.z[,c(1,3,5,7,9,2,4,6,8,10)],fontsize_col = 16,fontsize_row = 12,cluster_rows = TRUE, 
                   cluster_cols = FALSE,main = paste("expected.down",sep = " "),
                   color = colorRampPalette(c("navy","white","firebrick"))(20) ) 

#######################draw gene expression trend line 
##simple plot (transparent set to 25 (max 255), lwd=2)

drawTrendLine <- function(data=NULL,i=NULL){
  n <- ncol(data)
  plot(1:n,data[i,1:n],type='n',xaxt='n',xlab="",ylab="log2TPM",
       ylim=c(-2,2))
  title(row.names(data)[i],line = 0, font.main = 2, cex.main = 2, col.main = "black")
  #rect(xleft = 1,ybottom = 7.5, xright = 6.5, ytop = 9, col = rgb(100, 0, 0, 80, maxColorValue=255))
  #rect(xleft = 6.5,ybottom = 7.5, xright = 12.5, ytop = 9, col = rgb(0, 100, 0, 80, maxColorValue=255))
  lines(x = 1:5,y = data[i,1:5],type='b',lty=1,
        col=rgb(0, 200, 0, 90, maxColorValue=255),lwd=2.5)
  lines(x = 6:10,y = data[i,6:10],type='b',lty=1,
        col=rgb(200, 0, 0, 80, maxColorValue=255),lwd=2.5)
  abline(v=5.5,lty=2,lwd=1);
  axis(side = 1,at = 1:10,labels = c("WT.ESC","WT.6h",  "WT.12h",       "WT.24h",       "WT.48h",
                                     "KO.ESC","KO.6h","KO.12h","KO.24h","KO.48h"),
       las=2)
  box <- par('usr') #plot box range coordinates, x1,x2,y1,y2
  cat(box,row.names(data)[i],"\n")
  text(x=box[1]+1,y=box[4]-(box[4]-box[3])*0.1,labels = "WT")
  text(x=box[2]-1,y=box[4]-(box[4]-box[3])*0.1,labels = "KO")
  #legend("topright",legend = c("WT","KO"), col = c(rgb(0, 200, 0, 90, maxColorValue=255),
  #                                                 rgb(200, 0, 0, 80, maxColorValue=255)),lwd = 1
  #       )
  return (1)
}
######

par(mfrow=c(5,3),mar=c(3,5,3,3),oma=c(1,1,1,1))
#par(mfrow=c(6,3),mar=c(3,5,3,3),oma=c(1,1,1,1))
data <- expected.up.data.z 
#data <- expected.down.data.z
#data <- expected.down.data, the same with zscored, only y axis scale
for(i in 1:nrow(data)){
  #drawTrendLine(data=expected.down.logCPM,i)
  drawTrendLine(data=data,i)
}




##draw select.up select.down gene set
log2.select.up <- ave.log2.combine.all[select.up,c(1,3,5,7,9,2,4,6,8,10)]
log2.select.down <- ave.log2.combine.all[select.down,c(1,3,5,7,9,2,4,6,8,10)]
boxplot(log2.select.up,las=2)
boxplot(log2.select.down,las=2)
pheatmap::pheatmap(log2.select.up,cluster_cols = FALSE,main = paste(sample,"log2.select.up",sep = " ") )
pheatmap::pheatmap(log2.select.down,cluster_cols = FALSE,main = paste(sample,"log2.select.down",sep = " ") )


