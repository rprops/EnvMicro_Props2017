library("ggplot2")
library("dplyr")
library("plyr")

div.FCM <- read.csv2("fcm.diversity_total.csv")
div.16S <- read.csv2("otu.diversity16S_P.csv")
metadata <- read.csv2("metadata.csv")
metadata <- metadata[metadata$Platform == "Accuri",]
metadata$Sample_fcm <- gsub(metadata$Sample_fcm, pattern="_rep.*", replacement="")
metadata <- do.call(rbind,by(metadata, INDICES = factor(metadata$Sample_fcm), 
                             FUN = unique))

### Discard MI2013 data and inland lakes
metadata <- metadata[metadata$Lake != "Inland" & metadata$Year != 2013,]

### Replace F by Particle
metadata$Sample_16S <- gsub(metadata$Sample_16S, pattern="F", replace="P")


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
data.total.final$Lake <- factor(data.total.final$Lake)

### Explorative plot
png("alpha-div_PARTICLE_log.png",width=7*1.65,height=5*1.5,res=500,units="in")
p5 <- ggplot(data=data.total.final,aes(x=D2,y=D2.fcm))+ scale_fill_brewer(palette="Paired") +geom_point(shape=21,size=6,alpha=0.6,aes(fill=Lake))+
  theme_bw()+labs(x="Taxonomic diversity (D1) - 16S",y="Phenotypic diversity (D1) - FCM")+
  theme(axis.text=element_text(size=13),axis.title=element_text(size=18,face="bold"),legend.text=element_text(size=13),legend.title=element_text(size=13),strip.text.x = element_text(size = 22))+
  scale_x_continuous(trans='log2', breaks = seq(0,60, 10),minor_breaks =NULL) +
  scale_y_continuous(trans='log2', breaks = seq(1500,4250,250),minor_breaks =NULL) +
  geom_smooth(method="lm",color="black",fill="lightblue",formula=y~x)
print(p5)
dev.off()

### Get R squared
lm.P <- lm(log2(D2)~log2(D2.fcm), data=data.total.final)
summary(lm.P)$r.squared

cor.test(y=log2(data.total.final$D2), x=log2(data.total.final$D2.fcm))
