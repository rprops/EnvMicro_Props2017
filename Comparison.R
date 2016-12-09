library("ggplot2")
library("dplyr")

div.FCM <- read.csv2("fcm.diversity.csv")
div.16S <- read.csv2("otu.diversity16S.csv")
metadata <- read.csv2("metadata_reference.csv")

### Replace .renamed for MI13 data
div.16S$Sample <- gsub(div.16S$Sample, pattern=".renamed",replacement="",fixed=TRUE)
### Replace "-" from div.FCM
div.16S$Sample <- gsub(div.16S$Sample, pattern="-", replacement="")

### Merge data
data.16s <- inner_join(div.16S,metadata, by=c("Sample"="Sample_16S"))


### Join all data in one dataframe
data.total <- inner_join(div.FCM,data.16s, by=c("Sample_name"="Sample_fcm"))
  
### Select samples with > 1000 counts
data.total <- data.total[data.total$counts>10000,]

### remove outlier
data.total <- data.total[data.total$D1.y>10,]

df3 <- aggregate(D0.fcm~Sample,data.count,mean)
df3 <- left_join(df3,data.count,by="Sample")
### remove NA values
df3 <- df3[!is.na(df3$D0),]
df3 <- droplevels(df3)

### Explorative plot
png("alpha-div_combined_log.png",width=7*1.65,height=5*1.5,res=500,units="in")
p5 <- ggplot(data=data.total,aes(x=D0.y,y=D0.x))+ scale_fill_brewer(palette="Paired") +geom_point(shape=21,size=6,alpha=0.6,aes(fill=Lake))+
  theme_bw()+labs(x="Taxonomic diversity (D1) - 16S",y="Phenotypic diversity (D1) - FCM")+
  theme(axis.text=element_text(size=13),axis.title=element_text(size=18,face="bold"),legend.text=element_text(size=13),legend.title=element_text(size=13),strip.text.x = element_text(size = 22))+
  scale_x_continuous(trans='log2', breaks = seq(0,200, 40),minor_breaks =NULL) +
  scale_y_continuous(trans='log2', breaks = seq(1500,4250,250),minor_breaks =NULL) +
  geom_smooth(method="lm",color="black",fill="lightblue",formula=y~x)
print(p5)
dev.off()
