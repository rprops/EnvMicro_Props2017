library("ggplot2")
library("dplyr")
library("plyr")
library("grid")
library("gridExtra")

div.FCM <- read.csv2("fcm.diversity_total.csv")
div.16S <- read.csv2("otu.diversity16S_F.csv")
metadata <- read.csv2("metadata.csv")
metadata <- metadata[metadata$Platform == "Accuri",]
metadata$Sample_fcm <- gsub(metadata$Sample_fcm, pattern="_rep.*", replacement="")
metadata <- do.call(rbind,by(metadata, INDICES = factor(metadata$Sample_fcm), 
                        FUN = unique))

### Calculate means + errors for FCM data
groupLevels <- factor(gsub(div.FCM$Sample_names, pattern="_rep.*", replacement="", fixed=FALSE))
means <- do.call(rbind,by(div.FCM[,c(3:5,9:11)], INDICES = groupLevels, 
                          FUN = colMeans))
errors1 <- do.call(rbind,by(div.FCM[,6:8], INDICES = groupLevels, 
                           FUN = function(x) sqrt(colSums(x^2))/ncol(x)))
errors2 <- do.call(rbind,by(div.FCM[,9:10], INDICES = groupLevels, 
                            FUN = function(x) apply(x,2,sd)))
colnames(errors2) <- c("counts.sd", "volume.sd")
div.FCM.merged <- data.frame(sample_fcm = rownames(means), cbind(means[,1:3], errors1, means[,4:6], errors2))
### Rename colnames
colnames(div.FCM.merged)[1:7] <- c("Sample_fcm","D0.fcm","D1.fcm","D2.fcm","sd.D0.fcm","sd.D1.fcm","sd.D2.fcm")

### Replace .renamed for MI13 data
div.16S$Sample <- gsub(div.16S$Sample, pattern=".renamed",replacement="",fixed=TRUE)

### Replace "-" from div.16S
div.16S$Sample <- gsub(div.16S$Sample, pattern="-", replacement="")

### Remove RNA samples
div.16S <- div.16S[sapply(as.character(div.16S$Sample), function(x) substr(x, nchar(x), nchar(x))) != "R",]

### Replace "D" of DNA samples by ""

### Merge data
data.16s <- inner_join(div.16S, metadata, by=c("Sample"="Sample_16S"))

### Join all data in one dataframe 
data.16s$Sample_fcm <- gsub(data.16s$Sample_fcm, pattern="_rep.*", replacement="")
data.total <- inner_join(div.FCM.merged, data.16s, by="Sample_fcm")
  
### Select samples with > 10000 counts
# data.total <- data.total[data.total$counts>10000,]

### remove outlier
# data.total <- data.total[data.total$D1.y>10,]

### Add data from previous publication (cooling water)
div.ref <- read.csv2("div.ref.merged.csv")

data.total.final <- rbind.fill(data.total, div.ref[, c(2,4:9)])
data.total.final[(nrow(data.total)+1):nrow(data.total.final), 14:23] <- div.ref[, c(10:19)]
data.total.final$Lake <- as.character(data.total.final$Lake)
data.total.final$Lake[is.na(data.total.final$Lake)] <- "Cooling water"
data.total.final$Lake <- factor(data.total.final$Lake, levels=c("Michigan","Muskegon","Cooling water"))
data.total.final$Lake <- revalue(data.total.final$Lake, c("Michigan"="Lake Michigan", "Muskegon"="Muskegon lake"))

### Get R squared
lm.F <- lm(log2(D2)~log2(D2.fcm), data=data.total.final)
summary(lm.F)$r.squared

### Prepare to plot r squared / pearson's correlation
my_grob = grobTree(textGrob(bquote(r^2 == .(round(summary(lm.F)$r.squared, 2))), x=0.8,  y=0.16, hjust=0,
                            gp=gpar(col="black", fontsize=20, fontface="italic")))
my_grob2 = grobTree(textGrob(bquote(r[p] == .(round(cor(y=log2(data.total.final$D2), x=log2(data.total.final$D2.fcm)), 2))), x=0.8,  y=0.08, hjust=0,
                            gp=gpar(col="black", fontsize=20, fontface="italic")))

### Plot D2
png("alpha-div_D2_FREE_log.png",width=7*1.65,height=5*1.5,res=500,units="in")
p5 <- ggplot(data=data.total.final,aes(x=D2.fcm,y=D2, fill=Lake))+ scale_fill_brewer(palette="Paired") +geom_point(shape=21,size=6,alpha=0.6,aes(fill=Lake))+
  theme_bw()+labs(y=expression('Taxonomic diversity - D'[2]),x=expression('Phenotypic diversity - D'[2]), fill="Environment")+
  theme(axis.text=element_text(size=15),axis.title=element_text(size=20,face="bold"),legend.text=element_text(size=15),legend.title=element_text(size=16),strip.text.x = element_text(size = 22))+
  scale_y_continuous(trans='log2', breaks = seq(5,60, 10),minor_breaks =NULL) +
  scale_x_continuous(trans='log2', breaks = seq(1000,4250,250),minor_breaks =NULL) +
  geom_smooth(method="lm",color="black",fill="lightblue",formula=y~x)+
  annotation_custom(my_grob)+
  annotation_custom(my_grob2)
print(p5)
dev.off()

### Get R squared
lm.F <- lm(log2(D1)~log2(D1.fcm), data=data.total.final)
summary(lm.F)$r.squared


### Prepare to plot r squared / pearson's correlation
my_grob = grobTree(textGrob(bquote(r^2 == .(round(summary(lm.F)$r.squared, 2))), x=0.8,  y=0.16, hjust=0,
                            gp=gpar(col="black", fontsize=20, fontface="italic")))
my_grob2 = grobTree(textGrob(bquote(r[p] == .(round(cor(y=log2(data.total.final$D1), x=log2(data.total.final$D1.fcm)), 2))), x=0.8,  y=0.08, hjust=0,
                             gp=gpar(col="black", fontsize=20, fontface="italic")))

### Plot D1
png("alpha-div_D1_FREE_log.png",width=7*1.65,height=5*1.5,res=500,units="in")
p6 <- ggplot(data=data.total.final,aes(x=D1.fcm,y=D1, fill=Lake))+ scale_fill_brewer(palette="Paired") +geom_point(shape=21,size=6,alpha=0.6,aes(fill=Lake))+
  theme_bw()+labs(y=expression('Taxonomic diversity - D'[1]),x=expression('Phenotypic diversity - D'[1]), fill="Environment")+
  theme(axis.text=element_text(size=15),axis.title=element_text(size=18,face="bold"),legend.text=element_text(size=15),legend.title=element_text(size=16),strip.text.x = element_text(size = 22))+
  scale_y_continuous(trans='log2', breaks = seq(20,200, 40),minor_breaks =NULL) +
  scale_x_continuous(trans='log2', breaks = seq(1500,4250,500),minor_breaks =NULL) +
  geom_smooth(method="lm",color="black",fill="lightblue",formula=y~x)+
  annotation_custom(my_grob)+
  annotation_custom(my_grob2)
print(p6)
dev.off()

### Get R squared
lm.F <- lm(log2(D0)~log2(D0.fcm), data=data.total.final)
summary(lm.F)$r.squared

### Prepare to plot r squared / pearson's correlation
my_grob = grobTree(textGrob(bquote(r^2 == .(round(summary(lm.F)$r.squared, 2))), x=0.8,  y=0.16, hjust=0,
                            gp=gpar(col="black", fontsize=20, fontface="italic")))
my_grob2 = grobTree(textGrob(bquote(r[p] == .(round(cor(y=log2(data.total.final$D0), x=log2(data.total.final$D0.fcm)), 2))), x=0.8,  y=0.08, hjust=0,
                             gp=gpar(col="black", fontsize=20, fontface="italic")))

### Plot D0
png("alpha-div_D0_FREE_log.png",width=7*1.65,height=5*1.5,res=500,units="in")
p7 <- ggplot(data=data.total.final,aes(x=D0.fcm,y=D0,fill=Lake))+ scale_fill_brewer(palette="Paired") +geom_point(shape=21,size=6,alpha=0.6,aes(fill=Lake))+
  theme_bw()+labs(y=expression('Taxonomic diversity - D'[0]),x=expression('Phenotypic diversity - D'[0]), fill="Environment")+
  theme(axis.text=element_text(size=15),axis.title=element_text(size=18,face="bold"),legend.text=element_text(size=15),legend.title=element_text(size=16),strip.text.x = element_text(size = 22))+
  scale_y_continuous(trans='log2', breaks = seq(0,4000, 500),minor_breaks =NULL) +
  scale_x_continuous(trans='log2', breaks = seq(1500,20000,2500),minor_breaks =NULL) +
  geom_smooth(method="lm",color="black",fill="lightblue",formula=y~x)+
  annotation_custom(my_grob)+
  annotation_custom(my_grob2)
print(p7)
dev.off()

### All together
png("alpha-div_log_D1D2_together.png",width=1.6*7*1.65,height=5*1.5,res=500,units="in")
grid.arrange(p7,p6, ncol=2)
dev.off()

##################################################################################################
### The same figures but with error bars
##################################################################################################
### Get R squared

lm.F <- lm(log2(D0)~log2(D0.fcm), data=data.total.final)
summary(lm.F)$r.squared

### Prepare to plot r squared / pearson's correlation
my_grob = grobTree(textGrob(bquote(r^2 == .(round(summary(lm.F)$r.squared, 2))), x=0.8,  y=0.16, hjust=0,
                            gp=gpar(col="black", fontsize=20, fontface="italic")))
my_grob2 = grobTree(textGrob(bquote(r[p] == .(round(cor(y=log2(data.total.final$D0), x=log2(data.total.final$D0.fcm)), 2))), x=0.8,  y=0.08, hjust=0,
                             gp=gpar(col="black", fontsize=20, fontface="italic")))

p8 <- ggplot(data=data.total.final,aes(x=D0.fcm,y=D0,fill=Lake))+ scale_fill_brewer(palette="Paired") +geom_point(shape=21,size=6,alpha=0.6,aes(fill=Lake))+
  theme_bw()+labs(y="Taxonomic diversity (D0) - 16S",x="Phenotypic diversity (D0) - FCM", fill="Environment")+
  theme(axis.text=element_text(size=15),axis.title=element_text(size=18,face="bold"),legend.text=element_text(size=15),legend.title=element_text(size=16),strip.text.x = element_text(size = 22))+
  scale_y_continuous(trans='log2', breaks = seq(0,4000, 500),minor_breaks =NULL) +
  scale_x_continuous(trans='log2', breaks = seq(1500,20000,2500),minor_breaks =NULL) +
  geom_smooth(method="lm",color="black",fill="lightblue",formula=y~x)+
  annotation_custom(my_grob)+
  annotation_custom(my_grob2)+
  geom_errorbar(aes(ymin = D0 - sd.D0, ymax = D0 + 
                               sd.D0), width = 0.02, color="black")+
  geom_errorbarh(aes(xmin = D0.fcm - sd.D0.fcm, xmax = D0.fcm + 
                       sd.D0.fcm), height = 0.02, color="black")
print(p8)


## Get R squared
lm.F <- lm(log2(D1)~log2(D1.fcm), data=data.total.final)
summary(lm.F)$r.squared


### Prepare to plot r squared / pearson's correlation
my_grob = grobTree(textGrob(bquote(r^2 == .(round(summary(lm.F)$r.squared, 2))), x=0.8,  y=0.16, hjust=0,
                            gp=gpar(col="black", fontsize=20, fontface="italic")))
my_grob2 = grobTree(textGrob(bquote(r[p] == .(round(cor(y=log2(data.total.final$D1), x=log2(data.total.final$D1.fcm)), 2))), x=0.8,  y=0.08, hjust=0,
                             gp=gpar(col="black", fontsize=20, fontface="italic")))


p6 <- ggplot(data=data.total.final,aes(x=D1.fcm,y=D1, fill=Lake))+ scale_fill_brewer(palette="Paired") +geom_point(shape=21,size=6,alpha=0.6,aes(fill=Lake))+
  theme_bw()+labs(y="Taxonomic diversity (D1) - 16S",x="Phenotypic diversity (D1) - FCM", fill="Environment")+
  theme(axis.text=element_text(size=15),axis.title=element_text(size=18,face="bold"),legend.text=element_text(size=15),legend.title=element_text(size=16),strip.text.x = element_text(size = 22))+
  scale_y_continuous(trans='log2', breaks = seq(20,200, 40),minor_breaks =NULL) +
  scale_x_continuous(trans='log2', breaks = seq(1500,4250,500),minor_breaks =NULL) +
  geom_smooth(method="lm",color="black",fill="lightblue",formula=y~x)+
  annotation_custom(my_grob)+
  annotation_custom(my_grob2)+
  geom_errorbar(aes(ymin = D2 - sd.D2, ymax = D2 + 
                    sd.D2), width = 0.02, color="black")+
  geom_errorbarh(aes(xmin = D2.fcm - sd.D2.fcm, xmax = D2.fcm + 
                       sd.D2.fcm), height = 0.02, color="black")
print(p9)

### Get R squared
lm.F <- lm(log2(D2)~log2(D2.fcm), data=data.total.final)
summary(lm.F)$r.squared

### Prepare to plot r squared / pearson's correlation
my_grob = grobTree(textGrob(bquote(r^2 == .(round(summary(lm.F)$r.squared, 2))), x=0.8,  y=0.16, hjust=0,
                            gp=gpar(col="black", fontsize=20, fontface="italic")))
my_grob2 = grobTree(textGrob(bquote(r[p] == .(round(cor(y=log2(data.total.final$D2), x=log2(data.total.final$D2.fcm)), 2))), x=0.8,  y=0.08, hjust=0,
                             gp=gpar(col="black", fontsize=20, fontface="italic")))

p10 <- ggplot(data=data.total.final,aes(x=D2.fcm,y=D2, fill=Lake))+ scale_fill_brewer(palette="Paired") +geom_point(shape=21,size=6,alpha=0.6,aes(fill=Lake))+
  theme_bw()+labs(y="Taxonomic diversity (D2) - 16S",x="Phenotypic diversity (D2) - FCM", fill="Environment")+
  theme(axis.text=element_text(size=15),axis.title=element_text(size=18,face="bold"),legend.text=element_text(size=15),legend.title=element_text(size=16),strip.text.x = element_text(size = 22))+
  scale_y_continuous(trans='log2', breaks = seq(5,60, 10),minor_breaks =NULL) +
  scale_x_continuous(trans='log2', breaks = seq(1000,4250,250),minor_breaks =NULL) +
  geom_smooth(method="lm",color="black",fill="lightblue",formula=y~x)+
  annotation_custom(my_grob)+
  annotation_custom(my_grob2)+
  geom_errorbar(aes(ymin = D2 - sd.D2, ymax = D2 + 
                      sd.D2), width = 0.02, color="black")+
  geom_errorbarh(aes(xmin = D2.fcm - sd.D2.fcm, xmax = D2.fcm + 
                       sd.D2.fcm), height = 0.02, color="black")
print(p10)


### All together with error
png("alpha-div_log_all_together_errorbars.png",width=3*7*1.65,height=5*1.5,res=500,units="in")
grid.arrange(p8,p9,p10, ncol=3)
dev.off()
