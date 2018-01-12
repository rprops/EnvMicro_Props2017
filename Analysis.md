---
title: "Analysis"
author: "Ruben Props"
date: "14 januari 2017"
output:
  html_document:
    code_folding: show
    highlight: haddock
    keep_md: yes
    theme: united
    toc: yes
    toc_float:
      collapsed: no
      smooth_scroll: yes
      toc_depth: 2
    css: report_styles.css
editor_options: 
  chunk_output_type: console
---




The full analysis of the submitted paper: *Invasive dreissenid mussels induce shifts in bacterioplankton diversity through selective feeding on high nucleic acid bacteria*. 

Before starting the analysis please unzip *16S.zip* and download the FCM data for the mussel experiment from [here](https://flowrepository.org/experiments/1034) and store them in a directory named *data_mussel*. The reference FCM data (approx. 1GB) for the cooling water system and Lake Michigan/Muskegon Lake are available from [here](https://flowrepository.org/experiments/746) and [here](https://flowrepository.org/experiments/1047). These data should be stored in a directory called *data_reference* with subdirectories *FCM_CW* for the cooling water data and *FCM_MI* for the Lake Michigan and Muskegon Lake data.

The phenotypic and taxonomic diversities (`REF_diversity.csv`, `Lakes_diversity16S_F.csv` and `Lakes_diversityFCM_F.csv`) were calculated by means of the external scripts in `/extra_scripts` as the computation can take some time.

**Notice:** The bootstrap repeats are set rather low (`R = 3`) to facilitate faster run times. For exactly reproducing the output described in the manuscript please adjust the parameters to those specified in the methods section (`R = 100`). The output from this markdown file is provided in html format.

# Accuracy
Adjust this parameter (`R` = number of bootstraps) to increase the accuracy of the phenotypic diversity estimation. For exactly reproducing the results from the manuscript set `R.main = 100`. Be aware that this will significantly increase the run time.


```r
R.main <- 3
```

# Load libraries


```r
library("Phenoflow")
library("plyr")
library("dplyr")
library("gridExtra")
library("ggplot2")
library("formatR")
library("gridExtra")
library("cowplot")
library("RColorBrewer")
library("vegan")
library("sandwich")
library("lmtest")
library("grid")
library("car")
library("egg")
library("phyloseq")
library("splines")
library("caret") # for cross validation
library("foreach") # for parallelization
source("functions.R")
my.settings <- list(
  strip.background=list(col="transparent"),
  strip.border=list(col="transparent", cex=5),
  gate=list(col="black", fill="lightblue", alpha=0.2,border=NA,lwd=2),
  panel.background=list(col="lightgray"),
  background=list(col="white"))
```

# Part 1: Validation of alpha diversity  


```r
div.FCM <- read.csv2("files/Lakes_diversityFCM_F.csv")
div.16S <- read.csv2("files/Lakes_diversity16S_F.csv")
metadata <- read.csv2("files/Lakes_metadata.csv")
metadata <- metadata[metadata$Platform == "Accuri",]
metadata$Sample_fcm <- gsub(metadata$Sample_fcm, pattern="_rep.*", replacement="")
metadata <- do.call(rbind,by(metadata, INDICES = factor(metadata$Sample_fcm), 
                        FUN = unique))

# Calculate means + errors for FCM data
groupLevels <- factor(gsub(div.FCM$Sample_names, pattern="_rep.*", replacement="", fixed=FALSE))
means <- do.call(rbind,by(div.FCM[,c(3:5,9:12)], INDICES = groupLevels, 
                          FUN = colMeans))
errors1 <- do.call(rbind,by(div.FCM[,6:8], INDICES = groupLevels, 
                           FUN = function(x) sqrt(colSums(x^2))/ncol(x)))
errors2 <- do.call(rbind,by(div.FCM[,c(9,10,12)], INDICES = groupLevels, 
                            FUN = function(x) apply(x,2,sd)))
colnames(errors2) <- c("counts.sd", "volume.sd", "HNA_counts.sd")
div.FCM.merged <- data.frame(sample_fcm = rownames(means), cbind(means[,1:3], errors1, means[,4:7], errors2))

# Rename colnames
colnames(div.FCM.merged)[1:7] <- c("Sample_fcm","D0.fcm","D1.fcm","D2.fcm","sd.D0.fcm","sd.D1.fcm","sd.D2.fcm")

# Replace .renamed for MI13 data
div.16S$Sample <- gsub(div.16S$Sample, pattern=".renamed",replacement="",fixed=TRUE)

# Replace "-" from div.16S
div.16S$Sample <- gsub(div.16S$Sample, pattern="-", replacement="")

### Remove RNA samples
div.16S <- div.16S[sapply(as.character(div.16S$Sample), function(x) substr(x, nchar(x), nchar(x))) != "R",]

# Merge data
data.16s <- inner_join(div.16S, metadata, by=c("Sample"="Sample_16S"))
```

```
## Warning: Column `Sample`/`Sample_16S` joining character vector and factor,
## coercing into character vector
```

```r
# Join all data in one dataframe 
data.16s$Sample_fcm <- gsub(data.16s$Sample_fcm, pattern="_rep.*", replacement="")
data.total <- inner_join(div.FCM.merged, data.16s, by="Sample_fcm")
```

```
## Warning: Column `Sample_fcm` joining factor and character vector, coercing
## into character vector
```

```r
# Select samples with > 10000 counts
data.total <- data.total[data.total$counts>10000,]

# Add fcm/16S data from previous publication (cooling water)
div.ref <- read.csv2("files/REF_diversity.csv")

data.total.final <- rbind.fill(data.total, div.ref[, c(2,4:9)])
data.total.final[(nrow(data.total)+1):nrow(data.total.final), 16:25] <- div.ref[, c(10:19)]
data.total.final$Lake <- as.character(data.total.final$Lake)
data.total.final$Lake[is.na(data.total.final$Lake)] <- "Cooling water"
data.total.final$Lake <- factor(data.total.final$Lake, levels=c("Michigan","Muskegon","Cooling water"))
data.total.final$Lake <- revalue(data.total.final$Lake, c("Michigan"="Lake Michigan", "Muskegon"="Muskegon Lake"))
lb <- read.csv2("files/REF_labels.csv")
tmp <- left_join(data.total.final,lb, by=c("Sample_fcm"="Sample_fcm"))
```

```
## Warning: Column `Sample_fcm` joining character vector and factor, coercing
## into character vector
```

```r
data.total.final$Sample[data.total.final$Lake == "Cooling water"] <- tmp$Sample_seq[!is.na(tmp$Sample_seq)]; remove(tmp)
data.total.final$Lake[1:87][data.total.final$Site[1:87] == "MLB"] <- "Muskegon Lake"
data.total.final$Season[1:87][data.total.final$Season[1:87]=="Winter"] <- "Fall"

# Chlorophyl data
Chl <- read.csv("files/Lakes_metadata_Chl.csv")
Chl <- Chl[Chl$Sample_16S %in% data.total.final$Sample[data.total.final$Lake == "Lake Michigan"],]
mean(aggregate(Chl~Sample_16S, data=Chl, mean)$Chl)
```

```
## [1] 1.459744
```

```r
sd(aggregate(Chl~Sample_16S, data=Chl, mean)$Chl)
```

```
## [1] 1.191332
```

```r
# Get average HNA percentage and error in lake Mi
mean(100*data.total.final$HNA_counts[data.total.final$Lake=="Lake Michigan"]/data.total.final$counts[data.total.final$Lake=="Lake Michigan"])
```

```
## [1] 29.55895
```

```r
100*sqrt(sum((data.total.final$HNA_counts.sd[data.total.final$Lake=="Lake Michigan"]/data.total.final$HNA_counts[data.total.final$Lake=="Lake Michigan"])^2 + (data.total.final$counts.sd[data.total.final$Lake=="Lake Michigan"]/data.total.final$counts[data.total.final$Lake=="Lake Michigan"])^2)/length(data.total.final$counts[data.total.final$Lake=="Lake Michigan"]))
```

```
## [1] 4.240723
```

```r
length(data.total.final$counts[data.total.final$Lake=="Lake Michigan"])
```

```
## [1] 30
```

```r
# Get R squared for all diversity metrics
lm.F.D2 <- lm(log2(D2)~log2(D2.fcm), data=data.total.final)
summary(lm.F.D2)$r.squared
```

```
## [1] 0.8793401
```

```r
lm.F.D1 <- lm(log2(D1)~log2(D1.fcm), data=data.total.final)
summary(lm.F.D1)$r.squared
```

```
## [1] 0.8840665
```

```r
lm.F.D0 <- lm(log2(D0)~log2(D0.fcm), data=data.total.final)
summary(lm.F.D0)$r.squared
```

```
## [1] 0.2868273
```

```r
# Calculate unbiased R squared by tenfold cross validation
lmGrid <- expand.grid(intercept = TRUE)
R.cv.D2 <- train(log2(D2)~log2(D2.fcm), data=data.total.final, method ='lm', trControl = trainControl(method ="repeatedcv", repeats = 100), tuneGrid = lmGrid)
R.cv.D2$results
```

```
##   intercept      RMSE  Rsquared       MAE     RMSESD RsquaredSD     MAESD
## 1      TRUE 0.5694714 0.8864893 0.4790899 0.07833065 0.03649413 0.0745586
```

```r
R.cv.D1 <- train(log2(D1)~log2(D1.fcm), data=data.total.final, method ='lm', trControl = trainControl(method ="repeatedcv", repeats = 100), tuneGrid = lmGrid)
R.cv.D1$results
```

```
##   intercept      RMSE  Rsquared      MAE     RMSESD RsquaredSD      MAESD
## 1      TRUE 0.6091693 0.8938301 0.513774 0.09651116  0.0353859 0.08251614
```

```r
R.cv.D0 <- train(log2(D0)~log2(D0.fcm), data=data.total.final, method ='lm', trControl = trainControl(method ="repeatedcv", repeats = 100), tuneGrid = lmGrid)
R.cv.D0$results
```

```
##   intercept      RMSE  Rsquared       MAE    RMSESD RsquaredSD     MAESD
## 1      TRUE 0.7057909 0.3209198 0.5490113 0.1314418  0.1608102 0.1003106
```

```r
# Dynamic range of D2
max(data.total.final$D2)/min(data.total.final$D2)
```

```
## [1] 42.47004
```

```r
max(data.total.final$D1)/min(data.total.final$D1)
```

```
## [1] 88.71636
```

```r
# Dynamic range of D2 in previous study
max(data.total.final$D2[data.total.final$Lake=="Cooling water"])/min(data.total.final$D2[data.total.final$Lake=="Cooling water"])
```

```
## [1] 10.31645
```


## Figure S1: Check model assumptions  


```r
# D2 and D1 are highly correlated so we only use D2 in the feeding experiment
cor(data.total.final$D1.fcm,data.total.final$D2.fcm)
```

```
## [1] 0.9946115
```

```r
# Residual analysis of the D2 model
par(mfrow=c(1,3))
qqPlot(lm.F.D2, col="blue", reps=10000, ylab="Studentized residuals", xlab="Theoretical quantiles (t-distribution)",
       cex=1.5, las=1)
plot(residuals(lm.F.D2,"pearson"), x=predict(lm.F.D2), col="blue", las=1,
     ylab="Pearson residuals",xlab="Predicted values", cex=1.5)
lines(x=c(0,10), y=c(0,0), lty=2)
plot(y=log2(data.total.final$D2), x=predict(lm.F.D2), col="blue",
     ylab="Observed values",xlab="Predicted values", cex=1.5,
     las=1)
```

<img src="Figures/cached/Plot D2 residuals-1.png" style="display: block; margin: auto;" />

## Figure 1: Regression analysis  


```r
# Prepare to plot r squared / pearson's correlation
my_grob = grobTree(textGrob(bquote(r^2 == .(paste(round(R.cv.D2$results$Rsquared, 2)))), x=0.8,  y=0.16, hjust=0,
                            gp=gpar(col="black", fontsize=20, fontface="italic")))
my_grob2 = grobTree(textGrob(bquote(r[p] == .(round(cor(y=log2(data.total.final$D2), x=log2(data.total.final$D2.fcm)), 2))), x=0.8,  y=0.08, hjust=0,
                            gp=gpar(col="black", fontsize=20, fontface="italic")))

# Plot D2
p5 <- ggplot(data=data.total.final,aes(x=D2.fcm,y=D2, fill=Lake))+ scale_fill_manual(values=c("#88419d","#a6cee3","#fc8d62")) +
  geom_point(shape=21,size=6,alpha=0.6,aes(fill=Lake))+
  theme_bw()+labs(y=expression('Taxonomic diversity - D'[2]),x=expression('Phenotypic diversity - D'[2]), fill="Environment")+
  theme(axis.text=element_text(size=12.5),axis.title=element_text(size=20,face="bold"),legend.text=element_text(size=15),
        legend.title=element_text(size=16),strip.text.x = element_text(size = 22),
         legend.position = c(0.05, .97), legend.justification = c(0, 1),
        legend.background = element_rect(fill="transparent"))+
  scale_y_continuous(trans='log2', breaks = seq(5,60, 10),minor_breaks =NULL) +
  scale_x_continuous(trans='log2', breaks = seq(1000,4250,250),minor_breaks =NULL, limits = c(NA, 2700)) +
  geom_smooth(method="lm",color="black", fill ="lightblue", formula=y~x)+
  annotation_custom(my_grob)+
  annotation_custom(my_grob2)

# pdf("Fig1.pdf", width = 7, height = 6)
print(p5)
```

<img src="Figures/cached/Plot D2 regression-1.png" style="display: block; margin: auto;" />

```r
# dev.off()
```

## Figure S1: Compare slopes of individual regressions  


```r
# Supplementary figure for showing that the slopes for individual regressions of each lake are not significantly different
# from each other. This allows us to make one global regression for further analysis.

p5.si <- ggplot(data=data.total.final,aes(x=D2.fcm,y=D2, fill=Lake))+ scale_fill_manual(values=c("#88419d","#a6cee3","#fc8d62")) +
  geom_point(shape=21,size=6,alpha=0.6,aes(fill=Lake))+
  theme_bw()+labs(y=expression('Taxonomic diversity - D'[2]),x=expression('Phenotypic diversity - D'[2]), fill="Environment")+
  theme(axis.text=element_text(size=15),axis.title=element_text(size=20,face="bold"),legend.text=element_text(size=15),
        legend.title=element_text(size=16),strip.text.x = element_text(size = 22))+
  scale_y_continuous(trans='log2', breaks = seq(5,60, 10),minor_breaks =NULL) +
  scale_x_continuous(trans='log2', breaks = seq(1000,4250,250),minor_breaks =NULL, limits = c(NA, 2700)) +
  geom_smooth(method="lm",color="black",formula=y~x)


lm_comp <- lm(log2(D2)~log2(D2.fcm)*Lake, data=data.total.final)
contrast_matrix <- rbind(
  "Cooling water vs Muskegon Lake" = c(0,0,0,0,-1,1),
  "Lake Michigan vs Cooling water" = c(0,0,0,0,0,1),
  "Lake Michigan vs Muskegon Lake" = c(0,0,0,0,1,0)
)

# Adjusted R squared of individual models
summary(lm_comp)$adj.r.squared
```

```
## [1] 0.9269932
```

```r
# Adjusted R squared of single model
summary(lm(log2(D2)~log2(D2.fcm), data=data.total.final))$adj.r.squared
```

```
## [1] 0.8784529
```

```r
# Adjusted p-values using Benjamini & Hochberg (1995) correction
summary(glht(model = lm_comp, contrast_matrix), test = adjusted("BH"))
```

```
## 
## 	 Simultaneous Tests for General Linear Hypotheses
## 
## Fit: lm(formula = log2(D2) ~ log2(D2.fcm) * Lake, data = data.total.final)
## 
## Linear Hypotheses:
##                                     Estimate Std. Error t value Pr(>|t|)  
## Cooling water vs Muskegon Lake == 0   1.3229     0.5908   2.239   0.0647 .
## Lake Michigan vs Cooling water == 0   1.4795     0.7245   2.042   0.0647 .
## Lake Michigan vs Muskegon Lake == 0   0.1566     0.8692   0.180   0.8573  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## (Adjusted p values reported -- BH method)
```

```r
# pdf("Fig1_SI.pdf", width = 9, height = 6)
# png("Fig1_SI.png", width = 9, height = 6, res = 500, units = "in")
print(p5.si)
```

<img src="Figures/cached/Compare slopes of individual regressions-1.png" style="display: block; margin: auto;" />

```r
# dev.off()
```

## Figure S2: Regression of D0/D1  


```r
### Prepare to plot r squared / pearson's correlation
my_grob = grobTree(textGrob(bquote(r^2 == .(round(R.cv.D1$results$Rsquared, 2))), x=0.8,  y=0.16, hjust=0,
                            gp=gpar(col="black", fontsize=20, fontface="italic")))
my_grob2 = grobTree(textGrob(bquote(r[p] == .(round(cor(y=log2(data.total.final$D1), x=log2(data.total.final$D1.fcm)), 2))), x=0.8,  y=0.08, hjust=0,
                             gp=gpar(col="black", fontsize=20, fontface="italic")))

### Plot D1
p6 <- ggplot(data=data.total.final,aes(x=D1.fcm,y=D1, fill=Lake))+ scale_fill_manual(values=c("#88419d","#a6cee3","#fc8d62")) +
  geom_point(shape=21,size=6,alpha=0.6,aes(fill=Lake))+
  theme_bw()+labs(y=expression('Taxonomic diversity - D'[1]),x=expression('Phenotypic diversity - D'[1]), fill="Environment")+
  theme(axis.text=element_text(size=15),axis.title=element_text(size=18,face="bold"),legend.text=element_text(size=15), legend.title=element_text(size=16),strip.text.x = element_text(size = 22))+
  scale_y_continuous(trans='log2', breaks = seq(20,200, 40),minor_breaks =NULL) +
  scale_x_continuous(trans='log2', breaks = seq(1500,4250,500),minor_breaks =NULL) +
  geom_smooth(method="lm",color="black",fill="lightblue",formula=y~x)+
  annotation_custom(my_grob)+
  annotation_custom(my_grob2)

### Prepare to plot r squared / pearson's correlation
my_grob = grobTree(textGrob(bquote(r^2 == .(round(R.cv.D0$results$Rsquared, 2))), x=0.8,  y=0.16, hjust=0,
                            gp=gpar(col="black", fontsize=20, fontface="italic")))
my_grob2 = grobTree(textGrob(bquote(r[p] == .(round(cor(y=log2(data.total.final$D0), x=log2(data.total.final$D0.fcm)), 2))), x=0.8,  y=0.08, hjust=0,
                             gp=gpar(col="black", fontsize=20, fontface="italic")))

### Plot D0
p7 <- ggplot(data=data.total.final,aes(x=D0.fcm,y=D0,fill=Lake))+ scale_fill_manual(values=c("#88419d","#a6cee3","#fc8d62")) +
  geom_point(shape=21,size=6,alpha=0.6,aes(fill=Lake))+
  theme_bw()+labs(y=expression('Taxonomic diversity - D'[0]),x=expression('Phenotypic diversity - D'[0]),
                  fill="Environment")+
  theme(axis.text=element_text(size=15),axis.title=element_text(size=18,face="bold"),legend.text=element_text(size=15),
        legend.title=element_text(size=16),strip.text.x = element_text(size = 22))+
  scale_y_continuous(trans='log2', breaks = seq(0,4000, 500),minor_breaks =NULL) +
  scale_x_continuous(trans='log2', breaks = seq(1500,20000,2500),minor_breaks =NULL) +
  geom_smooth(method="lm",color="black",fill="lightblue",formula=y~x)+
  annotation_custom(my_grob)+
  annotation_custom(my_grob2)

### All together
grid.arrange(p7, p6, ncol=2)
```

<img src="Figures/cached/Plot D0 & D1 regression-1.png" style="display: block; margin: auto;" />

## Figure S4: Effect of sample size  


```r
# Load data
path = "data_mussel"
flowData_sc <- read.flowSet(path = path, transformation = FALSE, pattern=".fcs")

# Preprocess data according to standard protocol
flowData_sc <- transform(flowData_sc, `FL1-H` = asinh(`FL1-H`), `SSC-H` = asinh(`SSC-H`), 
                                  `FL3-H` = asinh(`FL3-H`), `FSC-H` = asinh(`FSC-H`))
param = c("FL1-H", "FL3-H", "SSC-H", "FSC-H")
flowData_sc = flowData_sc[, param]

# Test this for only one sample
flowData_sc <- flowData_sc[10]

# Create a PolygonGate for denoising the dataset Define coordinates for
# gate in sqrcut1 in format: c(x,x,x,x,y,y,y,y)
sqrcut1 <- matrix(c(8.5,8.5,15,15,3,8,14,3),ncol=2, nrow=4)
colnames(sqrcut1) <- c("FL1-H", "FL3-H")
polyGate1 <- polygonGate(.gate = sqrcut1, filterId = "Total Cells")

# Isolate only the cellular information based on the polyGate1
flowData_sc<- Subset(flowData_sc, polyGate1)

summary <- fsApply(x = flowData_sc, FUN = function(x) apply(x, 2, max), use.exprs = TRUE)
max = max(summary[, 1])
mytrans <- function(x) x/max
flowData_sc <- transform(flowData_sc, `FL1-H` = mytrans(`FL1-H`), 
                                  `FL3-H` = mytrans(`FL3-H`), `SSC-H` = mytrans(`SSC-H`), `FSC-H` = mytrans(`FSC-H`))

# Subsample at various depths and calculate diversity metrics with 100
# bootstraps Notice: this will use some CPU/RAM
for (i in c(10, 100, 200, 300, 400, 500, 750, 1000, 1250, 1500, 2000, 2500, 
            3000, 5000, 10000, 15000, 20000, 30000)) {
  for (j in 1:100) {
    fs1 <- FCS_resample(flowData_sc, replace = TRUE, sample = i, progress = FALSE)
    fp <- flowBasis(fs1, param, nbin = 128, bw = 0.01, normalize = function(x) x)
    div.tmp <- Diversity(fp, d = 3, R = 100, progress = FALSE)
    div.tmp <- cbind(div.tmp, size = i)
    if (j == 1) 
      results <- div.tmp else results <- rbind(results, div.tmp)
  }
  if (i == 10) 
    results.tot <- results else results.tot <- rbind(results.tot, results)
}

# Create plots
D0 <- ggplot(data = results.tot, aes(x = factor(size), y = D0)) + # geom_jitter(alpha=0.7, size=1)+
  geom_boxplot(alpha = 0.2, color = "blue", fill = "blue", size = 1) + 
  labs(x = "Sample size (nr. of cells)", y = expression('Phenotypic diversity - D'[0])) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45,hjust = 1),
                     axis.text=element_text(size=11),axis.title=element_text(size=15))

D1 <- ggplot(data = results.tot, aes(x = factor(size), y = D1)) + # geom_jitter(alpha=0.7, size=1)+
  geom_boxplot(alpha = 0.2, color = "blue", fill = "blue", size = 1) + 
  labs(x = "Sample size (nr. of cells)", y = expression('Phenotypic diversity - D'[1])) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.text=element_text(size=11), axis.title=element_text(size=15))

D2 <- ggplot(data = results.tot, aes(x = factor(size), y = D2)) + # geom_jitter(alpha=0.7, size=1)+
  geom_boxplot(alpha = 0.2, color = "blue", fill = "blue", size = 1) + 
  labs(x = "Sample size (nr. of cells)", y = expression('Phenotypic diversity - D'[2])) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.text=element_text(size=11), axis.title=element_text(size=15))

# Sample used for this assessment
grid.arrange(D0, D1, D2, ncol = 3, top = textGrob(paste("Sample used:",flowCore::sampleNames(flowData_sc)),                                                   gp = gpar(fontsize = 20, font = 3)))
```

<img src="Figures/cached/sample-size-1.png" style="display: block; margin: auto;" />

# Part 2: Validation of beta diversity  


```r
myColours2 <- brewer.pal(n=12,"Paired"); myColours2 <- myColours2[c(2,6,7)]

# Import otu data
physeq.otu <- import_mothur(mothur_shared_file ="16S/stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.shared" ,
                            mothur_constaxonomy_file = "16S/stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.0.03.cons.taxonomy")
otu_table(physeq.otu) <- t(otu_table(physeq.otu))
physeq.otu <- prune_samples(sample_sums(physeq.otu)>10000, physeq.otu)
physeq.otu <- prune_taxa(taxa_sums(physeq.otu)>0, physeq.otu)

# Select samples for which you have FCM data
sample_names(physeq.otu) <- gsub(sample_names(physeq.otu), pattern=".renamed", replacement="")
sample_names(physeq.otu) <- gsub(sample_names(physeq.otu), pattern=".renamed", replacement="")
meta.seq <- read.csv2("files/Lakes_metadata.csv")
physeq.otu <- prune_samples(sample_names(physeq.otu) %in% data.total.final$Sample, physeq.otu)

# Annotate with metadata
rownames(data.total.final) <- data.total.final$Sample
sample_data(physeq.otu) <- data.total.final

# Rescale
sample_sums(physeq.otu)
```

```
##          110D2F515          110D2F715          110D2F815 
##              15990              16536              21815 
##          110M2F115          110M2F415          110M2F515 
##              13836              13875              14547 
##          110M2F715          110S2F415          110S2F515 
##              14070              15451              13655 
##   Fa13.BD.MLB.DN.1   Fa13.BD.MLB.SN.1 Fa13.BD.MM110.SD.1 
##              34307              33378              33727 
## Fa13.BD.MM110.SN.1  Fa13.BD.MM15.DN.1  Fa13.BD.MM15.SD.1 
##              35194              34979              32169 
##  Fa13.BD.MM15.SN.1          M15S2F115          M15S2F515 
##              31872              15229              14730 
##          M15S2F815          M45D2F515          M45D2F715 
##              12456              15392              15393 
##          M45D2F815          M45S2F415          M45S2F515 
##              11619              19377              12282 
##          M45S2F715          M45S2F815          MBRE1F714 
##              16160              16068              33085 
##          MBRE1F715          MBRE1F914          MBRE2F515 
##              17553              24648              27886 
##          MBRE2F915          MBRH1F714          MBRH1F715 
##              21968              33640              17222 
##          MBRH1F914          MBRH2F515          MBRH2F915 
##              27095              26074              19342 
##          MDPE1F514          MDPE1F714          MDPE1F715 
##             106167              24136              18337 
##          MDPE1F914          MDPE2F915          MDPH1F714 
##              26251              21500              28656 
##          MDPH1F715          MDPH1F914          MDPH2F514 
##              20006              30768              21473 
##          MDPH2F515          MDPH2F915          MINE1F514 
##              21782              22252              54347 
##          MINE1F714          MINE1F715          MINE1F914 
##              11147              28549              24892 
##          MINE2F515          MINE2F915          MINH1F514 
##              11224              18584              39119 
##          MINH1F714          MINH1F715          MINH1F914 
##              30650              23635              26647 
##          MINH2F515          MINH2F915          MLBD2F115 
##              13599              22262              13959 
##          MLBD2F515          MLBD2F715          MLBD2F915 
##              13971              15557              24208 
##          MLBS2F115          MLBS2F815          MLBS2F915 
##              11569              11834              24681 
##          MOTE1F514          MOTE1F714          MOTE1F715 
##              41401              19672              22107 
##          MOTE1F914          MOTE2F515          MOTE2F915 
##              26792              16698              20850 
##          MOTH1F714          MOTH1F715          MOTH1F914 
##              30164              17235              25123 
##          MOTH1F915          MOTH2F514          MOTH2F515 
##              23287              31548              22608 
##   Sp13.BD.MLB.SN.1  Sp13.BD.MM15.DD.1  Sp13.BD.MM15.SD.1 
##              37613              36580              34122 
##  Sp13.BD.MM15.SN.1   Su13.BD.MLB.DD.1   Su13.BD.MLB.SD.1 
##              36650              34605              33366 
##  Su13.BD.MM15.DN.1  Su13.BD.MM15.SD.1  Su13.BD.MM15.SN.1 
##              31872              32172              35076
```

```r
physeq.otu <- scale_reads(physeq.otu)
sample_sums(physeq.otu)
```

```
##          110D2F515          110D2F715          110D2F815 
##              11027              11005              10996 
##          110M2F115          110M2F415          110M2F515 
##              11008              11007              11013 
##          110M2F715          110S2F415          110S2F515 
##              11009              11013              11027 
##   Fa13.BD.MLB.DN.1   Fa13.BD.MLB.SN.1 Fa13.BD.MM110.SD.1 
##              10502              10695              10925 
## Fa13.BD.MM110.SN.1  Fa13.BD.MM15.DN.1  Fa13.BD.MM15.SD.1 
##              10935              10820              10902 
##  Fa13.BD.MM15.SN.1          M15S2F115          M15S2F515 
##              10830              10972              10951 
##          M15S2F815          M45D2F515          M45D2F715 
##              10952              10982              10946 
##          M45D2F815          M45S2F415          M45S2F515 
##              10936              11022              10989 
##          M45S2F715          M45S2F815          MBRE1F714 
##              11016              10987              10897 
##          MBRE1F715          MBRE1F914          MBRE2F515 
##              10892              10785              10753 
##          MBRE2F915          MBRH1F714          MBRH1F715 
##              10877              10599              10662 
##          MBRH1F914          MBRH2F515          MBRH2F915 
##              10347              10625              10833 
##          MDPE1F514          MDPE1F714          MDPE1F715 
##              10503              10827              10866 
##          MDPE1F914          MDPE2F915          MDPH1F714 
##              10802              10841              10652 
##          MDPH1F715          MDPH1F914          MDPH2F514 
##              10653              10268              10851 
##          MDPH2F515          MDPH2F915          MINE1F514 
##              10783              10680               9705 
##          MINE1F714          MINE1F715          MINE1F914 
##              11147              10493              10459 
##          MINE2F515          MINE2F915          MINH1F514 
##              10270              10470               9970 
##          MINH1F714          MINH1F715          MINH1F914 
##              10611              10251              10485 
##          MINH2F515          MINH2F915          MLBD2F115 
##              10451              10646              10675 
##          MLBD2F515          MLBD2F715          MLBD2F915 
##              10761              10823              10634 
##          MLBS2F115          MLBS2F815          MLBS2F915 
##              10612              10722              10729 
##          MOTE1F514          MOTE1F714          MOTE1F715 
##              10525              10920              10960 
##          MOTE1F914          MOTE2F515          MOTE2F915 
##              10823              10888              10838 
##          MOTH1F714          MOTH1F715          MOTH1F914 
##              10848              10865              10819 
##          MOTH1F915          MOTH2F514          MOTH2F515 
##              10904              10822              10745 
##   Sp13.BD.MLB.SN.1  Sp13.BD.MM15.DD.1  Sp13.BD.MM15.SD.1 
##              10054              10325              10617 
##  Sp13.BD.MM15.SN.1   Su13.BD.MLB.DD.1   Su13.BD.MLB.SD.1 
##              10432              10718              10816 
##  Su13.BD.MM15.DN.1  Su13.BD.MM15.SD.1  Su13.BD.MM15.SN.1 
##              10913              10930              10892
```

```r
physeq.otu <- transform_sample_counts(physeq.otu, function(x) x/sum(x))

# Run beta diversity analysis on 16s data
pcoa <- ordinate(
  physeq = physeq.otu, 
  method = "PCoA", 
  distance = "bray",
  correction = "lingoes",
  k=2
)

pcoa.df <- data.frame(pcoa$vectors, sample_data(physeq.otu))
var <- round(pcoa$values$Eigenvalues/sum(pcoa$values$Eigenvalues)*100,1)

# Start beta diversity analysis on FCM data
path = "data_reference/FCM_MI"
flowData_transformed  <- read.flowSet(path = path, transformation = FALSE, pattern=".fcs")

flowData_transformed <- transform(flowData_transformed,`FL1-H`=asinh(`FL1-H`), 
                                  `SSC-H`=asinh(`SSC-H`), 
                                  `FL3-H`=asinh(`FL3-H`), 
                                  `FSC-H`=asinh(`FSC-H`))
param=c("FL1-H", "FL3-H","SSC-H","FSC-H")
flowData_transformed = flowData_transformed[,param]

# Create a PolygonGate for denoising the dataset
# Define coordinates for gate in sqrcut1 in format: c(x,x,x,x,y,y,y,y)
sqrcut1 <- matrix(c(8.5,8.5,15,15,3,8,14,3),ncol=2, nrow=4)
colnames(sqrcut1) <- c("FL1-H","FL3-H")
polyGate1 <- polygonGate(.gate=sqrcut1, filterId = "Total Cells")

#  Gating quality check
xyplot(`FL3-H` ~ `FL1-H`, data=flowData_transformed[1], filter=polyGate1,
       scales=list(y=list(limits=c(0,14)),
                   x=list(limits=c(6,16))),
       par.settings = my.settings,
       axis = axis.default, nbin=125, 
       par.strip.text=list(col="white", font=2, cex=2), smooth=FALSE)
```

<img src="Figures/cached/beta-diversity analysis-1.png" style="display: block; margin: auto;" />

```r
# Isolate only the cellular information based on the polyGate1
flowData_transformed <- Subset(flowData_transformed, polyGate1)

summary <- fsApply(x=flowData_transformed,FUN=function(x) apply(x,2,max),use.exprs=TRUE)
max = max(summary[,1])
mytrans <- function(x) x/max
flowData_transformed <- transform(flowData_transformed,`FL1-H`=mytrans(`FL1-H`),
                                  `FL3-H`=mytrans(`FL3-H`), 
                                  `SSC-H`=mytrans(`SSC-H`),
                                  `FSC-H`=mytrans(`FSC-H`))

# Calculate fingerprint with bw = 0.01
fbasis <- flowBasis(flowData_transformed, param, nbin=128, 
                    bw=0.01,normalize=function(x) x)

# Run beta diversity
pos <- gsub(rownames(fbasis@basis),pattern="_rep.*", replacement="")
beta.div <- beta_div_fcm(fbasis, INDICES=pos, ord.type="PCoA")
df.beta.fcm <- data.frame(Sample = rownames(beta.div$points), beta.div$points)
df.beta.fcm$Sample <- gsub(df.beta.fcm$Sample, pattern="MI5", replacement = "M15")
df.beta.fcm <- inner_join(df.beta.fcm, data.total.final, by=c("Sample"="Sample_fcm"))
var2 <- round(beta.div$eig/sum(beta.div$eig)*100, 1)
df.beta.fcm <- droplevels(df.beta.fcm)

# Procrustes analysis
fbasis1 <- fbasis
fbasis1@basis <- fbasis1@basis/apply(fbasis1@basis, 1, max)
fbasis1@basis <- round(fbasis1@basis, 3)
x <- by(fbasis1@basis, INDICES = pos, FUN = colMeans)
x <- do.call(rbind, x)
rownames(x) <- gsub(rownames(x), pattern = "MI5", replacement = "M15")
x <- x[(rownames(x) %in% data.total.final$Sample_fcm), ]

# Rename and order rows of fcm/seq data
tmp <- data.frame(sample_fcm=rownames(x))
tmp <- left_join(tmp, data.total.final, by=c("sample_fcm"="Sample_fcm"))
```

```
## Warning: Column `sample_fcm`/`Sample_fcm` joining factor and character
## vector, coercing into character vector
```

```r
rownames(x) <- tmp$Sample; remove(tmp)
x <- x[order(rownames(x)),]
x.data <- data.total.final[data.total.final$Sample %in% rownames(x),]
x.data <- x.data[order(x.data$Sample),]
otu_table(physeq.otu) <- otu_table(physeq.otu)[order(rownames(otu_table(physeq.otu))),]

# Run PcoA
dist.fcm <- vegan::vegdist(x)
pcoa.fcm <- cmdscale(dist.fcm)
dist.seq <- vegan::vegdist(otu_table(physeq.otu))
pcoa.seq <- cmdscale(dist.seq)
```

## Procrustes analysis  


```r
# Run procrustes + permutation
# Permutations are constrained within each sampling year 
perm <- how(nperm = 999)
setBlocks(perm) <- with(x.data, Year)
dist.prot <- vegan::protest(pcoa.seq,pcoa.fcm, permutations = perm)
summary(dist.prot)
```

```
## 
## Call:
## vegan::protest(X = pcoa.seq, Y = pcoa.fcm, permutations = perm) 
## 
## Number of objects: 87    Number of dimensions: 2 
## 
## Procrustes sum of squares:  
##  0.5749417 
## Procrustes root mean squared error: 
##  0.08129284 
## Quantiles of Procrustes errors:
##         Min          1Q      Median          3Q         Max 
## 0.006092608 0.050629245 0.070314993 0.093284024 0.223174785 
## 
## Rotation matrix:
##           [,1]       [,2]
## [1,] 0.9926241  0.1212331
## [2,] 0.1212331 -0.9926241
## 
## Translation of averages:
##               [,1]          [,2]
## [1,] -1.368612e-18 -1.536175e-18
## 
## Scaling of target:
## [1] 0.651965
```

```r
plot(dist.prot)
```

<img src="Figures/cached/procrustes-analysis-1.png" style="display: block; margin: auto;" />

## Figure 2: Beta diversity analysis  


```r
# Run PERMANOVA to evaluate if conclusions are similar between FCM/seq
# Similar variances across seasons
disper.fcm <- betadisper(dist.fcm, group=x.data$Season)
print(disper.fcm)
```

```
## 
## 	Homogeneity of multivariate dispersions
## 
## Call: betadisper(d = dist.fcm, group = x.data$Season)
## 
## No. of Positive Eigenvalues: 49
## No. of Negative Eigenvalues: 37
## 
## Average distance to median:
##   Fall Spring Summer 
## 0.1561 0.1137 0.1163 
## 
## Eigenvalues for PCoA axes:
##  PCoA1  PCoA2  PCoA3  PCoA4  PCoA5  PCoA6  PCoA7  PCoA8 
## 1.3346 0.3899 0.1224 0.1047 0.0752 0.0671 0.0580 0.0345
```

```r
disper.seq <- betadisper(dist.seq, group=x.data$Season)
print(disper.seq)
```

```
## 
## 	Homogeneity of multivariate dispersions
## 
## Call: betadisper(d = dist.seq, group = x.data$Season)
## 
## No. of Positive Eigenvalues: 60
## No. of Negative Eigenvalues: 26
## 
## Average distance to median:
##   Fall Spring Summer 
## 0.3110 0.3192 0.3248 
## 
## Eigenvalues for PCoA axes:
##  PCoA1  PCoA2  PCoA3  PCoA4  PCoA5  PCoA6  PCoA7  PCoA8 
## 2.6076 2.1516 1.1698 1.0706 0.5386 0.4867 0.4393 0.3406
```

```r
# Permutations are constrained within each sampling year 
perm <- how(nperm = 999)
setBlocks(perm) <- with(x.data, Year)
permanova.fcm <- adonis(dist.fcm ~ Season * Lake, data = x.data, 
    permutations = perm)
permanova.seq <- adonis(dist.seq ~ Season * Lake, data = data.frame(sample_data(physeq.otu)), 
    permutations = perm)

# Plot FCM beta diversity
my_grob = grobTree(textGrob(bquote(paste(r[Season]^2 == 
    .(round(100 * permanova.fcm$aov.tab[1, 5], 1)), 
    "%")), x = 0.7, y = 0.95, hjust = 0, gp = gpar(col = "black", 
    fontsize = 14, fontface = "italic")))
my_grob2 = grobTree(textGrob(bquote(paste(r[Lake]^2 == 
    .(format(round(100 * permanova.fcm$aov.tab[2, 5], 
        1), nsmall = 1)), "%")), x = 0.7, y = 0.87, 
    hjust = 0, gp = gpar(col = "black", fontsize = 14, 
        fontface = "italic")))
my_grob3 = grobTree(textGrob(bquote(paste(r[Season:Lake]^2 == 
    .(round(100 * permanova.fcm$aov.tab[3, 5], 1)), 
    "%")), x = 0.7, y = 0.79, hjust = 0, gp = gpar(col = "black", 
    fontsize = 14, fontface = "italic")))

df.beta.fcm$Season <- factor(df.beta.fcm$Season, levels = c("Spring", "Summer", "Fall"))
pcoa.df$Season <- factor(pcoa.df$Season, levels = c("Spring", "Summer", "Fall"))

beta.pcoa.fcm <- ggplot(data=df.beta.fcm, aes(x=X1, y=-X2, shape=Lake))+
  scale_shape_manual(values=c(21,24))+
  geom_point(size=7, aes(fill = Season), alpha=0.7)+
  theme_bw()+
  scale_fill_manual(values=myColours2)+
  scale_colour_manual(values=myColours2)+
  labs(x = paste0("PCoA axis 1 (",var2[1], "%)"), y = paste0("PCoA axis 2 (",var2[2], "%)"), fill="", shape = "", colour="",
       title="B")+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=16))+ 
  guides(fill = guide_legend(override.aes = list(shape = 22)))+
  annotation_custom(my_grob)+
  annotation_custom(my_grob2)+
  annotation_custom(my_grob3)

# Plot 16S beta diversity
my_grob = grobTree(textGrob(bquote(paste(r[Season]^2 == 
    .(round(100 * permanova.seq$aov.tab[1, 5], 1)), 
    "%")), x = 0.7, y = 0.95, hjust = 0, gp = gpar(col = "black", 
    fontsize = 14, fontface = "italic")))
my_grob2 = grobTree(textGrob(bquote(paste(r[Lake]^2 == 
    .(round(100 * permanova.seq$aov.tab[2, 5], 1)), 
    "%")), x = 0.7, y = 0.87, hjust = 0, gp = gpar(col = "black", 
    fontsize = 14, fontface = "italic")))
my_grob3 = grobTree(textGrob(bquote(paste(r[Season:Lake]^2 == 
    .(round(100 * permanova.seq$aov.tab[3, 5], 1)), 
    "%")), x = 0.7, y = 0.79, hjust = 0, gp = gpar(col = "black", 
    fontsize = 14, fontface = "italic")))

beta.pcoa <- ggplot(data=pcoa.df, aes(x=Axis.1, y=Axis.2, shape=Lake))+
  geom_point(alpha=0.7, size=7, aes(fill=Season))+
  scale_shape_manual(values=c(21,24))+
  theme_bw()+
  scale_fill_manual(values=myColours2)+
  scale_colour_manual(values=myColours2)+
  labs(x = paste0("PCoA axis 1 (",var[1], "%)"), y = paste0("PCoA axis 2 (",var[2], "%)"), fill="", shape="", colour="",
       title="A")+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=16))+
  guides(fill = guide_legend(override.aes = list(shape = 22)))+
  annotation_custom(my_grob)+
  annotation_custom(my_grob2)+
  annotation_custom(my_grob3)+
  ylim(-0.4,0.3)

# Both beta diversity plots together
grid_arrange_shared_legend(beta.pcoa, beta.pcoa.fcm, ncol=2)
```

<img src="Figures/cached/permanova-analysis-1.png" style="display: block; margin: auto;" />

# Part 3: Mussel experiment


```r
# Set seed for reproducible analysis
set.seed(777)

myColours <- brewer.pal("Accent",n=3)

# Samples were diluted 2x
dilution <- 2

path = "data_mussel"
flowData <- read.flowSet(path = path, transformation = FALSE, pattern=".fcs")

flowData_transformed <- transform(flowData,`FL1-H`=asinh(`FL1-H`), 
                                  `SSC-H`=asinh(`SSC-H`), 
                                  `FL3-H`=asinh(`FL3-H`), 
                                  `FSC-H`=asinh(`FSC-H`))
param=c("FL1-H", "FL3-H","SSC-H","FSC-H")
flowData_transformed = flowData_transformed[,param]
remove(flowData)

# Create a PolygonGate for denoising the dataset
# Define coordinates for gate in sqrcut1 in format: c(x,x,x,x,y,y,y,y)
sqrcut1 <- matrix(c(8.5,8.5,15,15,3,8,14,3),ncol=2, nrow=4)
colnames(sqrcut1) <- c("FL1-H","FL3-H")
polyGate1 <- polygonGate(.gate=sqrcut1, filterId = "Total Cells")

#  Gating quality check
xyplot(`FL3-H` ~ `FL1-H`, data=flowData_transformed[1], filter=polyGate1,
       scales=list(y=list(limits=c(0,14)),
                   x=list(limits=c(6,16))),
       par.settings = my.settings,
       axis = axis.default, nbin=125, 
       par.strip.text=list(col="white", font=2, cex=2), smooth=FALSE)
```

<img src="Figures/cached/Mussel-fcm-analysis-1.png" style="display: block; margin: auto;" />

```r
# Isolate only the cellular information based on the polyGate1
flowData_transformed <- Subset(flowData_transformed, polyGate1)

summary <- fsApply(x=flowData_transformed,FUN=function(x) apply(x,2,max),use.exprs=TRUE)
max = max(summary[,1])
mytrans <- function(x) x/max
flowData_transformed <- transform(flowData_transformed,`FL1-H`=mytrans(`FL1-H`),
                                  `FL3-H`=mytrans(`FL3-H`), 
                                  `SSC-H`=mytrans(`SSC-H`),
                                  `FSC-H`=mytrans(`FSC-H`))


# Calculate fingerprint with bw = 0.01
fbasis <- flowBasis(flowData_transformed, param, nbin=128, 
                    bw=0.01,normalize=function(x) x)

# Calculate ecological parameters from normalized fingerprint 
# Densities will be normalized to the interval [0,1]
# d = rounding factor
# Diversity.fbasis <- Diversity(fbasis, d = 3, plot = TRUE, R = 999)
Diversity.fbasis <- Diversity_rf(flowData_transformed, d=3, param = param, R = R.main,
                                 parallel = TRUE, ncores = 10)
```

```
## --- parameters are already normalized at: 1
## Fri Jan 12 13:37:51 2018 --- Using 10 cores for calculations
## Fri Jan 12 13:39:32 2018 --- Closing workers
## Fri Jan 12 13:39:32 2018 --- Alpha diversity metrics (D0,D1,D2) have been computed after 3 bootstraps
## -----------------------------------------------------------------------------------------------------
## 
```

```r
# make metadata table
tmp <- strsplit(rownames(Diversity.fbasis)[1:42], "_")
meta.div <- cbind(do.call(rbind, lapply(tmp, rbind))[,1],rep("C",42), do.call(rbind, lapply(tmp, rbind))[,2:3])
tmp <- strsplit(rownames(Diversity.fbasis)[43:nrow(Diversity.fbasis)], "_")
meta.div <- data.frame(rbind(meta.div, do.call(rbind, lapply(tmp, rbind))))
colnames(meta.div) <- c("Sample","Treatment","Time","Replicate")
meta.div$Replicate <- gsub(meta.div$Replicate,pattern=".fcs",replacement="")
meta.div$Time <- as.numeric(gsub(meta.div$Time,pattern="t",replacement=""))


# Merge with Diversity.fbasis
Diversity.fbasis <- cbind(Diversity.fbasis, meta.div)

# Remove outliers due to human errors (forgot staining)
fbasis@basis <- fbasis@basis[Diversity.fbasis$D2>1600,]
flowData_transformed <- flowData_transformed[Diversity.fbasis$D2>1600]

sqrcut1 <- matrix(c(asinh(12500),asinh(12500),15,15,3,9.55,14,3)/max,ncol=2, nrow=4)
colnames(sqrcut1) <- c("FL1-H","FL3-H")
rGate_HNA <- polygonGate(.gate=sqrcut1, filterId = "HNA")
sqrcut1 <- matrix(c(8.5,8.5,asinh(12500),asinh(12500),3,8,9.55,3)/max,ncol=2, nrow=4)
colnames(sqrcut1) <- c("FL1-H","FL3-H")
rGate_LNA <- polygonGate(.gate=sqrcut1, filterId = "LNA")

# Diversities of HNA/LNA populations
flowData_HNA <- split(flowData_transformed, rGate_HNA)$`HNA+`
flowData_LNA <- split(flowData_transformed, rGate_LNA)$`LNA+`

Diversity.HNA <- Diversity_rf(flowData_HNA, d=3, param = param, R = R.main,
                                 parallel = TRUE, ncores = 10)
```

```
## --- parameters are already normalized at: 1
## Fri Jan 12 13:39:51 2018 --- Using 10 cores for calculations
## Fri Jan 12 13:41:04 2018 --- Closing workers
## Fri Jan 12 13:41:04 2018 --- Alpha diversity metrics (D0,D1,D2) have been computed after 3 bootstraps
## -----------------------------------------------------------------------------------------------------
## 
```

```r
Diversity.LNA <- Diversity_rf(flowData_LNA, d=3, param = param, R = R.main,
                                 parallel = TRUE, ncores = 10)
```

```
## --- parameters are already normalized at: 0.675405293390532
## Fri Jan 12 13:41:11 2018 --- Using 10 cores for calculations
## Fri Jan 12 13:42:33 2018 --- Closing workers
## Fri Jan 12 13:42:33 2018 --- Alpha diversity metrics (D0,D1,D2) have been computed after 3 bootstraps
## -----------------------------------------------------------------------------------------------------
## 
```

```r
# Counts
# Normalize total cell gate
sqrcut1 <- matrix(c(8.5,8.5,15,15,3,8,14,3)/max,ncol=2, nrow=4)
colnames(sqrcut1) <- c("FL1-H","FL3-H")
polyGate1 <- polygonGate(.gate=sqrcut1, filterId = "Total Cells")

# Check if rectangle gate is correct, if not, adjust rGate_HNA
xyplot(`FL3-H` ~ `FL1-H`, data=flowData_transformed[1], filter=rGate_HNA,
       scales=list(y=list(limits=c(0,1)),
                   x=list(limits=c(0.4,1))),
       par.settings = my.settings,
       axis = axis.default, nbin=125, par.strip.text=list(col="white", font=2, 
                                                          cex=2), smooth=FALSE)
```

<img src="Figures/cached/Mussel-fcm-analysis-2.png" style="display: block; margin: auto;" />

```r
# Extract the cell counts
a <- flowCore::filter(flowData_transformed, rGate_HNA) 
HNACount <- flowCore::summary(a);HNACount <- toTable(HNACount)
```

```
## filter summary for frame 'C1_t0.5_rep1.fcs'
##  HNA+: 10977 of 39921 events (27.50%)
## 
## filter summary for frame 'C1_t0.5_rep2.fcs'
##  HNA+: 9725 of 36349 events (26.75%)
## 
## filter summary for frame 'C1_t0.5_rep3.fcs'
##  HNA+: 11115 of 39708 events (27.99%)
## 
## filter summary for frame 'C1_t0_rep1.fcs'
##  HNA+: 12719 of 43774 events (29.06%)
## 
## filter summary for frame 'C1_t0_rep2.fcs'
##  HNA+: 12661 of 44090 events (28.72%)
## 
## filter summary for frame 'C1_t0_rep3.fcs'
##  HNA+: 12771 of 44245 events (28.86%)
## 
## filter summary for frame 'C1_t1.5_rep1.fcs'
##  HNA+: 11536 of 40154 events (28.73%)
## 
## filter summary for frame 'C1_t1.5_rep2.fcs'
##  HNA+: 11433 of 40051 events (28.55%)
## 
## filter summary for frame 'C1_t1.5_rep3.fcs'
##  HNA+: 11106 of 39771 events (27.92%)
## 
## filter summary for frame 'C1_t1_rep1.fcs'
##  HNA+: 10811 of 38171 events (28.32%)
## 
## filter summary for frame 'C1_t1_rep2.fcs'
##  HNA+: 10494 of 37386 events (28.07%)
## 
## filter summary for frame 'C1_t1_rep3.fcs'
##  HNA+: 10771 of 38103 events (28.27%)
## 
## filter summary for frame 'C1_t2.5_rep1.fcs'
##  HNA+: 11737 of 39868 events (29.44%)
## 
## filter summary for frame 'C1_t2.5_rep2.fcs'
##  HNA+: 10859 of 37900 events (28.65%)
## 
## filter summary for frame 'C1_t2.5_rep3.fcs'
##  HNA+: 11482 of 38643 events (29.71%)
## 
## filter summary for frame 'C1_t2_rep1.fcs'
##  HNA+: 10917 of 37341 events (29.24%)
## 
## filter summary for frame 'C1_t2_rep2.fcs'
##  HNA+: 11163 of 37141 events (30.06%)
## 
## filter summary for frame 'C1_t2_rep3.fcs'
##  HNA+: 11063 of 36991 events (29.91%)
## 
## filter summary for frame 'C1_t3_rep1.fcs'
##  HNA+: 11535 of 37261 events (30.96%)
## 
## filter summary for frame 'C1_t3_rep2.fcs'
##  HNA+: 11542 of 37972 events (30.40%)
## 
## filter summary for frame 'C1_t3_rep3.fcs'
##  HNA+: 10985 of 36724 events (29.91%)
## 
## filter summary for frame 'C2_t0.5_rep1.fcs'
##  HNA+: 11085 of 39496 events (28.07%)
## 
## filter summary for frame 'C2_t0.5_rep2.fcs'
##  HNA+: 11129 of 39259 events (28.35%)
## 
## filter summary for frame 'C2_t0.5_rep3.fcs'
##  HNA+: 10915 of 38377 events (28.44%)
## 
## filter summary for frame 'C2_t0_rep1.fcs'
##  HNA+: 11929 of 42475 events (28.08%)
## 
## filter summary for frame 'C2_t0_rep2.fcs'
##  HNA+: 12176 of 42977 events (28.33%)
## 
## filter summary for frame 'C2_t0_rep3.fcs'
##  HNA+: 11914 of 42520 events (28.02%)
## 
## filter summary for frame 'C2_t1.5_rep1.fcs'
##  HNA+: 10945 of 39550 events (27.67%)
## 
## filter summary for frame 'C2_t1.5_rep2.fcs'
##  HNA+: 10426 of 37784 events (27.59%)
## 
## filter summary for frame 'C2_t1.5_rep3.fcs'
##  HNA+: 11135 of 39448 events (28.23%)
## 
## filter summary for frame 'C2_t1_rep1.fcs'
##  HNA+: 10572 of 38337 events (27.58%)
## 
## filter summary for frame 'C2_t1_rep2.fcs'
##  HNA+: 10710 of 38460 events (27.85%)
## 
## filter summary for frame 'C2_t1_rep3.fcs'
##  HNA+: 10891 of 39325 events (27.69%)
## 
## filter summary for frame 'C2_t2.5_rep1.fcs'
##  HNA+: 11814 of 40801 events (28.96%)
## 
## filter summary for frame 'C2_t2.5_rep2.fcs'
##  HNA+: 11141 of 39612 events (28.13%)
## 
## filter summary for frame 'C2_t2.5_rep3.fcs'
##  HNA+: 11583 of 40514 events (28.59%)
## 
## filter summary for frame 'C2_t2_rep1.fcs'
##  HNA+: 10730 of 37025 events (28.98%)
## 
## filter summary for frame 'C2_t2_rep2.fcs'
##  HNA+: 10951 of 37068 events (29.54%)
## 
## filter summary for frame 'C2_t2_rep3.fcs'
##  HNA+: 10885 of 37070 events (29.36%)
## 
## filter summary for frame 'C2_t3_rep1.fcs'
##  HNA+: 11452 of 37833 events (30.27%)
## 
## filter summary for frame 'C2_t3_rep2.fcs'
##  HNA+: 11138 of 37492 events (29.71%)
## 
## filter summary for frame 'C2_t3_rep3.fcs'
##  HNA+: 10909 of 36336 events (30.02%)
## 
## filter summary for frame 'Q1_S_t0.5_rep1.fcs'
##  HNA+: 10313 of 39437 events (26.15%)
## 
## filter summary for frame 'Q1_S_t0.5_rep2.fcs'
##  HNA+: 9616 of 36845 events (26.10%)
## 
## filter summary for frame 'Q1_S_t0.5_rep3.fcs'
##  HNA+: 10023 of 37990 events (26.38%)
## 
## filter summary for frame 'Q1_S_t0_rep1.fcs'
##  HNA+: 10698 of 39258 events (27.25%)
## 
## filter summary for frame 'Q1_S_t0_rep2.fcs'
##  HNA+: 10612 of 38206 events (27.78%)
## 
## filter summary for frame 'Q1_S_t0_rep3.fcs'
##  HNA+: 10619 of 38614 events (27.50%)
## 
## filter summary for frame 'Q1_S_t1.5_rep1.fcs'
##  HNA+: 8525 of 36130 events (23.60%)
## 
## filter summary for frame 'Q1_S_t1.5_rep2.fcs'
##  HNA+: 8274 of 34500 events (23.98%)
## 
## filter summary for frame 'Q1_S_t1.5_rep3.fcs'
##  HNA+: 8630 of 36188 events (23.85%)
## 
## filter summary for frame 'Q1_S_t1_rep1.fcs'
##  HNA+: 9790 of 39855 events (24.56%)
## 
## filter summary for frame 'Q1_S_t1_rep2.fcs'
##  HNA+: 9790 of 39920 events (24.52%)
## 
## filter summary for frame 'Q1_S_t1_rep3.fcs'
##  HNA+: 9789 of 39579 events (24.73%)
## 
## filter summary for frame 'Q1_S_t2.5_rep1.fcs'
##  HNA+: 8238 of 36594 events (22.51%)
## 
## filter summary for frame 'Q1_S_t2.5_rep2.fcs'
##  HNA+: 8056 of 35549 events (22.66%)
## 
## filter summary for frame 'Q1_S_t2.5_rep3.fcs'
##  HNA+: 7968 of 34987 events (22.77%)
## 
## filter summary for frame 'Q1_S_t2_rep1.fcs'
##  HNA+: 8657 of 37802 events (22.90%)
## 
## filter summary for frame 'Q1_S_t2_rep2.fcs'
##  HNA+: 8580 of 37597 events (22.82%)
## 
## filter summary for frame 'Q1_S_t2_rep3.fcs'
##  HNA+: 8219 of 36221 events (22.69%)
## 
## filter summary for frame 'Q1_S_t3_rep1.fcs'
##  HNA+: 7391 of 33469 events (22.08%)
## 
## filter summary for frame 'Q1_S_t3_rep2.fcs'
##  HNA+: 7469 of 33167 events (22.52%)
## 
## filter summary for frame 'Q1_S_t3_rep3.fcs'
##  HNA+: 7596 of 33781 events (22.49%)
## 
## filter summary for frame 'Q1_T_t0.5_rep1.fcs'
##  HNA+: 10223 of 38112 events (26.82%)
## 
## filter summary for frame 'Q1_T_t0.5_rep2.fcs'
##  HNA+: 10000 of 37470 events (26.69%)
## 
## filter summary for frame 'Q1_T_t0.5_rep3.fcs'
##  HNA+: 10051 of 37234 events (26.99%)
## 
## filter summary for frame 'Q1_T_t0_rep1.fcs'
##  HNA+: 11448 of 41808 events (27.38%)
## 
## filter summary for frame 'Q1_T_t0_rep2.fcs'
##  HNA+: 11378 of 41238 events (27.59%)
## 
## filter summary for frame 'Q1_T_t0_rep3.fcs'
##  HNA+: 11631 of 42362 events (27.46%)
## 
## filter summary for frame 'Q1_T_t1.5_rep1.fcs'
##  HNA+: 9056 of 36969 events (24.50%)
## 
## filter summary for frame 'Q1_T_t1.5_rep2.fcs'
##  HNA+: 8901 of 35827 events (24.84%)
## 
## filter summary for frame 'Q1_T_t1.5_rep3.fcs'
##  HNA+: 9429 of 37267 events (25.30%)
## 
## filter summary for frame 'Q1_T_t1_rep1.fcs'
##  HNA+: 9802 of 38324 events (25.58%)
## 
## filter summary for frame 'Q1_T_t1_rep2.fcs'
##  HNA+: 9269 of 37267 events (24.87%)
## 
## filter summary for frame 'Q1_T_t1_rep3.fcs'
##  HNA+: 9849 of 37978 events (25.93%)
## 
## filter summary for frame 'Q1_T_t2.5_rep1.fcs'
##  HNA+: 12806 of 40318 events (31.76%)
## 
## filter summary for frame 'Q1_T_t2.5_rep2.fcs'
##  HNA+: 12254 of 39291 events (31.19%)
## 
## filter summary for frame 'Q1_T_t2.5_rep3.fcs'
##  HNA+: 12212 of 40051 events (30.49%)
## 
## filter summary for frame 'Q1_T_t2_rep1.fcs'
##  HNA+: 8589 of 35543 events (24.17%)
## 
## filter summary for frame 'Q1_T_t2_rep2.fcs'
##  HNA+: 8391 of 35093 events (23.91%)
## 
## filter summary for frame 'Q1_T_t2_rep3.fcs'
##  HNA+: 8410 of 34497 events (24.38%)
## 
## filter summary for frame 'Q1_T_t3_rep1.fcs'
##  HNA+: 8239 of 34647 events (23.78%)
## 
## filter summary for frame 'Q1_T_t3_rep2.fcs'
##  HNA+: 6960 of 31325 events (22.22%)
## 
## filter summary for frame 'Q2_S_t0.5_rep1.fcs'
##  HNA+: 10986 of 42808 events (25.66%)
## 
## filter summary for frame 'Q2_S_t0.5_rep2.fcs'
##  HNA+: 10577 of 41557 events (25.45%)
## 
## filter summary for frame 'Q2_S_t0.5_rep3.fcs'
##  HNA+: 11017 of 43172 events (25.52%)
## 
## filter summary for frame 'Q2_S_t0_rep1.fcs'
##  HNA+: 11151 of 39581 events (28.17%)
## 
## filter summary for frame 'Q2_S_t0_rep2.fcs'
##  HNA+: 10828 of 38881 events (27.85%)
## 
## filter summary for frame 'Q2_S_t0_rep3.fcs'
##  HNA+: 10704 of 38546 events (27.77%)
## 
## filter summary for frame 'Q2_S_t1.5_rep1.fcs'
##  HNA+: 8645 of 35695 events (24.22%)
## 
## filter summary for frame 'Q2_S_t1.5_rep2.fcs'
##  HNA+: 8724 of 36598 events (23.84%)
## 
## filter summary for frame 'Q2_S_t1.5_rep3.fcs'
##  HNA+: 8283 of 35348 events (23.43%)
## 
## filter summary for frame 'Q2_S_t1_rep1.fcs'
##  HNA+: 9399 of 38506 events (24.41%)
## 
## filter summary for frame 'Q2_S_t1_rep2.fcs'
##  HNA+: 9678 of 38967 events (24.84%)
## 
## filter summary for frame 'Q2_S_t1_rep3.fcs'
##  HNA+: 9376 of 38451 events (24.38%)
## 
## filter summary for frame 'Q2_S_t2.5_rep1.fcs'
##  HNA+: 8611 of 36921 events (23.32%)
## 
## filter summary for frame 'Q2_S_t2.5_rep2.fcs'
##  HNA+: 8768 of 37198 events (23.57%)
## 
## filter summary for frame 'Q2_S_t2.5_rep3.fcs'
##  HNA+: 8268 of 35876 events (23.05%)
## 
## filter summary for frame 'Q2_S_t2_rep1.fcs'
##  HNA+: 9008 of 38345 events (23.49%)
## 
## filter summary for frame 'Q2_S_t2_rep2.fcs'
##  HNA+: 8692 of 37232 events (23.35%)
## 
## filter summary for frame 'Q2_S_t2_rep3.fcs'
##  HNA+: 8228 of 35964 events (22.88%)
## 
## filter summary for frame 'Q2_S_t3_rep1.fcs'
##  HNA+: 7906 of 34635 events (22.83%)
## 
## filter summary for frame 'Q2_S_t3_rep2.fcs'
##  HNA+: 7607 of 33967 events (22.40%)
## 
## filter summary for frame 'Q2_S_t3_rep3.fcs'
##  HNA+: 7629 of 34083 events (22.38%)
## 
## filter summary for frame 'Q2_T_t0.5_rep1.fcs'
##  HNA+: 11049 of 41673 events (26.51%)
## 
## filter summary for frame 'Q2_T_t0.5_rep2.fcs'
##  HNA+: 11004 of 41576 events (26.47%)
## 
## filter summary for frame 'Q2_T_t0.5_rep3.fcs'
##  HNA+: 10966 of 41430 events (26.47%)
## 
## filter summary for frame 'Q2_T_t0_rep1.fcs'
##  HNA+: 11726 of 43282 events (27.09%)
## 
## filter summary for frame 'Q2_T_t0_rep2.fcs'
##  HNA+: 11147 of 40883 events (27.27%)
## 
## filter summary for frame 'Q2_T_t0_rep3.fcs'
##  HNA+: 12050 of 44044 events (27.36%)
## 
## filter summary for frame 'Q2_T_t1.5_rep1.fcs'
##  HNA+: 8759 of 36219 events (24.18%)
## 
## filter summary for frame 'Q2_T_t1.5_rep2.fcs'
##  HNA+: 8469 of 35309 events (23.99%)
## 
## filter summary for frame 'Q2_T_t1.5_rep3.fcs'
##  HNA+: 8591 of 35596 events (24.13%)
## 
## filter summary for frame 'Q2_T_t1_rep1.fcs'
##  HNA+: 9739 of 38373 events (25.38%)
## 
## filter summary for frame 'Q2_T_t1_rep2.fcs'
##  HNA+: 9692 of 38487 events (25.18%)
## 
## filter summary for frame 'Q2_T_t1_rep3.fcs'
##  HNA+: 9571 of 38141 events (25.09%)
## 
## filter summary for frame 'Q2_T_t2.5_rep1.fcs'
##  HNA+: 8629 of 36953 events (23.35%)
## 
## filter summary for frame 'Q2_T_t2.5_rep3.fcs'
##  HNA+: 8630 of 37502 events (23.01%)
## 
## filter summary for frame 'Q2_T_t2_rep1.fcs'
##  HNA+: 9067 of 38297 events (23.68%)
## 
## filter summary for frame 'Q2_T_t2_rep2.fcs'
##  HNA+: 8422 of 36492 events (23.08%)
## 
## filter summary for frame 'Q2_T_t2_rep3.fcs'
##  HNA+: 8701 of 37242 events (23.36%)
## 
## filter summary for frame 'Q2_T_t3_rep1.fcs'
##  HNA+: 7900 of 33150 events (23.83%)
## 
## filter summary for frame 'Q2_T_t3_rep2.fcs'
##  HNA+: 7984 of 33737 events (23.67%)
## 
## filter summary for frame 'Q2_T_t3_rep3.fcs'
##  HNA+: 7810 of 33290 events (23.46%)
## 
## filter summary for frame 'Q3_S_t0.5_rep1.fcs'
##  HNA+: 10953 of 41447 events (26.43%)
## 
## filter summary for frame 'Q3_S_t0.5_rep2.fcs'
##  HNA+: 10328 of 40278 events (25.64%)
## 
## filter summary for frame 'Q3_S_t0.5_rep3.fcs'
##  HNA+: 10676 of 40841 events (26.14%)
## 
## filter summary for frame 'Q3_S_t0_rep1.fcs'
##  HNA+: 10668 of 38163 events (27.95%)
## 
## filter summary for frame 'Q3_S_t0_rep3.fcs'
##  HNA+: 10986 of 38368 events (28.63%)
## 
## filter summary for frame 'Q3_S_t1.5_rep1.fcs'
##  HNA+: 8835 of 35328 events (25.01%)
## 
## filter summary for frame 'Q3_S_t1.5_rep2.fcs'
##  HNA+: 8421 of 34111 events (24.69%)
## 
## filter summary for frame 'Q3_S_t1.5_rep3.fcs'
##  HNA+: 8731 of 35854 events (24.35%)
## 
## filter summary for frame 'Q3_S_t1_rep1.fcs'
##  HNA+: 10098 of 39785 events (25.38%)
## 
## filter summary for frame 'Q3_S_t1_rep2.fcs'
##  HNA+: 9883 of 39500 events (25.02%)
## 
## filter summary for frame 'Q3_S_t1_rep3.fcs'
##  HNA+: 10031 of 40081 events (25.03%)
## 
## filter summary for frame 'Q3_S_t2.5_rep1.fcs'
##  HNA+: 7909 of 34110 events (23.19%)
## 
## filter summary for frame 'Q3_S_t2.5_rep2.fcs'
##  HNA+: 7943 of 34462 events (23.05%)
## 
## filter summary for frame 'Q3_S_t2.5_rep3.fcs'
##  HNA+: 7851 of 33877 events (23.18%)
## 
## filter summary for frame 'Q3_S_t2_rep1.fcs'
##  HNA+: 8455 of 36134 events (23.40%)
## 
## filter summary for frame 'Q3_S_t2_rep2.fcs'
##  HNA+: 8379 of 35622 events (23.52%)
## 
## filter summary for frame 'Q3_S_t2_rep3.fcs'
##  HNA+: 8381 of 35652 events (23.51%)
## 
## filter summary for frame 'Q3_S_t3_rep1.fcs'
##  HNA+: 7947 of 35052 events (22.67%)
## 
## filter summary for frame 'Q3_S_t3_rep2.fcs'
##  HNA+: 7494 of 34025 events (22.02%)
## 
## filter summary for frame 'Q3_S_t3_rep3.fcs'
##  HNA+: 7885 of 34806 events (22.65%)
## 
## filter summary for frame 'Q3_T_t0.5_rep1.fcs'
##  HNA+: 10602 of 39768 events (26.66%)
## 
## filter summary for frame 'Q3_T_t0.5_rep2.fcs'
##  HNA+: 10413 of 40221 events (25.89%)
## 
## filter summary for frame 'Q3_T_t0.5_rep3.fcs'
##  HNA+: 10517 of 40525 events (25.95%)
## 
## filter summary for frame 'Q3_T_t0_rep1.fcs'
##  HNA+: 11141 of 40427 events (27.56%)
## 
## filter summary for frame 'Q3_T_t0_rep2.fcs'
##  HNA+: 10778 of 39238 events (27.47%)
## 
## filter summary for frame 'Q3_T_t0_rep3.fcs'
##  HNA+: 11057 of 39677 events (27.87%)
## 
## filter summary for frame 'Q3_T_t1.5_rep1.fcs'
##  HNA+: 8857 of 36149 events (24.50%)
## 
## filter summary for frame 'Q3_T_t1.5_rep2.fcs'
##  HNA+: 8584 of 35048 events (24.49%)
## 
## filter summary for frame 'Q3_T_t1.5_rep3.fcs'
##  HNA+: 8534 of 34974 events (24.40%)
## 
## filter summary for frame 'Q3_T_t1_rep1.fcs'
##  HNA+: 10005 of 39922 events (25.06%)
## 
## filter summary for frame 'Q3_T_t1_rep2.fcs'
##  HNA+: 9380 of 38560 events (24.33%)
## 
## filter summary for frame 'Q3_T_t1_rep3.fcs'
##  HNA+: 9926 of 39378 events (25.21%)
## 
## filter summary for frame 'Q3_T_t2.5_rep1.fcs'
##  HNA+: 8073 of 35912 events (22.48%)
## 
## filter summary for frame 'Q3_T_t2.5_rep2.fcs'
##  HNA+: 6763 of 32456 events (20.84%)
## 
## filter summary for frame 'Q3_T_t2.5_rep3.fcs'
##  HNA+: 7589 of 34403 events (22.06%)
## 
## filter summary for frame 'Q3_T_t2_rep1.fcs'
##  HNA+: 8605 of 38700 events (22.24%)
## 
## filter summary for frame 'Q3_T_t2_rep2.fcs'
##  HNA+: 8126 of 35818 events (22.69%)
## 
## filter summary for frame 'Q3_T_t2_rep3.fcs'
##  HNA+: 8454 of 38120 events (22.18%)
## 
## filter summary for frame 'Q3_T_t3_rep1.fcs'
##  HNA+: 7265 of 32594 events (22.29%)
## 
## filter summary for frame 'Q3_T_t3_rep2.fcs'
##  HNA+: 7427 of 32993 events (22.51%)
## 
## filter summary for frame 'Q3_T_t3_rep3.fcs'
##  HNA+: 7486 of 33277 events (22.50%)
```

```r
s <- flowCore::filter(flowData_transformed, polyGate1)
TotalCount <- flowCore::summary(s);TotalCount <- flowCore::toTable(TotalCount)
```

```
## filter summary for frame 'C1_t0.5_rep1.fcs'
##  Total Cells+: 39921 of 39921 events (100.00%)
## 
## filter summary for frame 'C1_t0.5_rep2.fcs'
##  Total Cells+: 36349 of 36349 events (100.00%)
## 
## filter summary for frame 'C1_t0.5_rep3.fcs'
##  Total Cells+: 39708 of 39708 events (100.00%)
## 
## filter summary for frame 'C1_t0_rep1.fcs'
##  Total Cells+: 43774 of 43774 events (100.00%)
## 
## filter summary for frame 'C1_t0_rep2.fcs'
##  Total Cells+: 44090 of 44090 events (100.00%)
## 
## filter summary for frame 'C1_t0_rep3.fcs'
##  Total Cells+: 44245 of 44245 events (100.00%)
## 
## filter summary for frame 'C1_t1.5_rep1.fcs'
##  Total Cells+: 40154 of 40154 events (100.00%)
## 
## filter summary for frame 'C1_t1.5_rep2.fcs'
##  Total Cells+: 40051 of 40051 events (100.00%)
## 
## filter summary for frame 'C1_t1.5_rep3.fcs'
##  Total Cells+: 39771 of 39771 events (100.00%)
## 
## filter summary for frame 'C1_t1_rep1.fcs'
##  Total Cells+: 38171 of 38171 events (100.00%)
## 
## filter summary for frame 'C1_t1_rep2.fcs'
##  Total Cells+: 37386 of 37386 events (100.00%)
## 
## filter summary for frame 'C1_t1_rep3.fcs'
##  Total Cells+: 38103 of 38103 events (100.00%)
## 
## filter summary for frame 'C1_t2.5_rep1.fcs'
##  Total Cells+: 39868 of 39868 events (100.00%)
## 
## filter summary for frame 'C1_t2.5_rep2.fcs'
##  Total Cells+: 37900 of 37900 events (100.00%)
## 
## filter summary for frame 'C1_t2.5_rep3.fcs'
##  Total Cells+: 38643 of 38643 events (100.00%)
## 
## filter summary for frame 'C1_t2_rep1.fcs'
##  Total Cells+: 37341 of 37341 events (100.00%)
## 
## filter summary for frame 'C1_t2_rep2.fcs'
##  Total Cells+: 37141 of 37141 events (100.00%)
## 
## filter summary for frame 'C1_t2_rep3.fcs'
##  Total Cells+: 36991 of 36991 events (100.00%)
## 
## filter summary for frame 'C1_t3_rep1.fcs'
##  Total Cells+: 37261 of 37261 events (100.00%)
## 
## filter summary for frame 'C1_t3_rep2.fcs'
##  Total Cells+: 37972 of 37972 events (100.00%)
## 
## filter summary for frame 'C1_t3_rep3.fcs'
##  Total Cells+: 36724 of 36724 events (100.00%)
## 
## filter summary for frame 'C2_t0.5_rep1.fcs'
##  Total Cells+: 39496 of 39496 events (100.00%)
## 
## filter summary for frame 'C2_t0.5_rep2.fcs'
##  Total Cells+: 39259 of 39259 events (100.00%)
## 
## filter summary for frame 'C2_t0.5_rep3.fcs'
##  Total Cells+: 38377 of 38377 events (100.00%)
## 
## filter summary for frame 'C2_t0_rep1.fcs'
##  Total Cells+: 42475 of 42475 events (100.00%)
## 
## filter summary for frame 'C2_t0_rep2.fcs'
##  Total Cells+: 42977 of 42977 events (100.00%)
## 
## filter summary for frame 'C2_t0_rep3.fcs'
##  Total Cells+: 42520 of 42520 events (100.00%)
## 
## filter summary for frame 'C2_t1.5_rep1.fcs'
##  Total Cells+: 39550 of 39550 events (100.00%)
## 
## filter summary for frame 'C2_t1.5_rep2.fcs'
##  Total Cells+: 37784 of 37784 events (100.00%)
## 
## filter summary for frame 'C2_t1.5_rep3.fcs'
##  Total Cells+: 39448 of 39448 events (100.00%)
## 
## filter summary for frame 'C2_t1_rep1.fcs'
##  Total Cells+: 38337 of 38337 events (100.00%)
## 
## filter summary for frame 'C2_t1_rep2.fcs'
##  Total Cells+: 38460 of 38460 events (100.00%)
## 
## filter summary for frame 'C2_t1_rep3.fcs'
##  Total Cells+: 39325 of 39325 events (100.00%)
## 
## filter summary for frame 'C2_t2.5_rep1.fcs'
##  Total Cells+: 40801 of 40801 events (100.00%)
## 
## filter summary for frame 'C2_t2.5_rep2.fcs'
##  Total Cells+: 39612 of 39612 events (100.00%)
## 
## filter summary for frame 'C2_t2.5_rep3.fcs'
##  Total Cells+: 40514 of 40514 events (100.00%)
## 
## filter summary for frame 'C2_t2_rep1.fcs'
##  Total Cells+: 37025 of 37025 events (100.00%)
## 
## filter summary for frame 'C2_t2_rep2.fcs'
##  Total Cells+: 37068 of 37068 events (100.00%)
## 
## filter summary for frame 'C2_t2_rep3.fcs'
##  Total Cells+: 37070 of 37070 events (100.00%)
## 
## filter summary for frame 'C2_t3_rep1.fcs'
##  Total Cells+: 37833 of 37833 events (100.00%)
## 
## filter summary for frame 'C2_t3_rep2.fcs'
##  Total Cells+: 37492 of 37492 events (100.00%)
## 
## filter summary for frame 'C2_t3_rep3.fcs'
##  Total Cells+: 36336 of 36336 events (100.00%)
## 
## filter summary for frame 'Q1_S_t0.5_rep1.fcs'
##  Total Cells+: 39437 of 39437 events (100.00%)
## 
## filter summary for frame 'Q1_S_t0.5_rep2.fcs'
##  Total Cells+: 36845 of 36845 events (100.00%)
## 
## filter summary for frame 'Q1_S_t0.5_rep3.fcs'
##  Total Cells+: 37990 of 37990 events (100.00%)
## 
## filter summary for frame 'Q1_S_t0_rep1.fcs'
##  Total Cells+: 39258 of 39258 events (100.00%)
## 
## filter summary for frame 'Q1_S_t0_rep2.fcs'
##  Total Cells+: 38206 of 38206 events (100.00%)
## 
## filter summary for frame 'Q1_S_t0_rep3.fcs'
##  Total Cells+: 38614 of 38614 events (100.00%)
## 
## filter summary for frame 'Q1_S_t1.5_rep1.fcs'
##  Total Cells+: 36130 of 36130 events (100.00%)
## 
## filter summary for frame 'Q1_S_t1.5_rep2.fcs'
##  Total Cells+: 34500 of 34500 events (100.00%)
## 
## filter summary for frame 'Q1_S_t1.5_rep3.fcs'
##  Total Cells+: 36188 of 36188 events (100.00%)
## 
## filter summary for frame 'Q1_S_t1_rep1.fcs'
##  Total Cells+: 39855 of 39855 events (100.00%)
## 
## filter summary for frame 'Q1_S_t1_rep2.fcs'
##  Total Cells+: 39920 of 39920 events (100.00%)
## 
## filter summary for frame 'Q1_S_t1_rep3.fcs'
##  Total Cells+: 39579 of 39579 events (100.00%)
## 
## filter summary for frame 'Q1_S_t2.5_rep1.fcs'
##  Total Cells+: 36594 of 36594 events (100.00%)
## 
## filter summary for frame 'Q1_S_t2.5_rep2.fcs'
##  Total Cells+: 35549 of 35549 events (100.00%)
## 
## filter summary for frame 'Q1_S_t2.5_rep3.fcs'
##  Total Cells+: 34987 of 34987 events (100.00%)
## 
## filter summary for frame 'Q1_S_t2_rep1.fcs'
##  Total Cells+: 37802 of 37802 events (100.00%)
## 
## filter summary for frame 'Q1_S_t2_rep2.fcs'
##  Total Cells+: 37597 of 37597 events (100.00%)
## 
## filter summary for frame 'Q1_S_t2_rep3.fcs'
##  Total Cells+: 36221 of 36221 events (100.00%)
## 
## filter summary for frame 'Q1_S_t3_rep1.fcs'
##  Total Cells+: 33469 of 33469 events (100.00%)
## 
## filter summary for frame 'Q1_S_t3_rep2.fcs'
##  Total Cells+: 33167 of 33167 events (100.00%)
## 
## filter summary for frame 'Q1_S_t3_rep3.fcs'
##  Total Cells+: 33781 of 33781 events (100.00%)
## 
## filter summary for frame 'Q1_T_t0.5_rep1.fcs'
##  Total Cells+: 38112 of 38112 events (100.00%)
## 
## filter summary for frame 'Q1_T_t0.5_rep2.fcs'
##  Total Cells+: 37470 of 37470 events (100.00%)
## 
## filter summary for frame 'Q1_T_t0.5_rep3.fcs'
##  Total Cells+: 37234 of 37234 events (100.00%)
## 
## filter summary for frame 'Q1_T_t0_rep1.fcs'
##  Total Cells+: 41808 of 41808 events (100.00%)
## 
## filter summary for frame 'Q1_T_t0_rep2.fcs'
##  Total Cells+: 41238 of 41238 events (100.00%)
## 
## filter summary for frame 'Q1_T_t0_rep3.fcs'
##  Total Cells+: 42362 of 42362 events (100.00%)
## 
## filter summary for frame 'Q1_T_t1.5_rep1.fcs'
##  Total Cells+: 36969 of 36969 events (100.00%)
## 
## filter summary for frame 'Q1_T_t1.5_rep2.fcs'
##  Total Cells+: 35827 of 35827 events (100.00%)
## 
## filter summary for frame 'Q1_T_t1.5_rep3.fcs'
##  Total Cells+: 37267 of 37267 events (100.00%)
## 
## filter summary for frame 'Q1_T_t1_rep1.fcs'
##  Total Cells+: 38324 of 38324 events (100.00%)
## 
## filter summary for frame 'Q1_T_t1_rep2.fcs'
##  Total Cells+: 37267 of 37267 events (100.00%)
## 
## filter summary for frame 'Q1_T_t1_rep3.fcs'
##  Total Cells+: 37978 of 37978 events (100.00%)
## 
## filter summary for frame 'Q1_T_t2.5_rep1.fcs'
##  Total Cells+: 40318 of 40318 events (100.00%)
## 
## filter summary for frame 'Q1_T_t2.5_rep2.fcs'
##  Total Cells+: 39291 of 39291 events (100.00%)
## 
## filter summary for frame 'Q1_T_t2.5_rep3.fcs'
##  Total Cells+: 40051 of 40051 events (100.00%)
## 
## filter summary for frame 'Q1_T_t2_rep1.fcs'
##  Total Cells+: 35543 of 35543 events (100.00%)
## 
## filter summary for frame 'Q1_T_t2_rep2.fcs'
##  Total Cells+: 35093 of 35093 events (100.00%)
## 
## filter summary for frame 'Q1_T_t2_rep3.fcs'
##  Total Cells+: 34497 of 34497 events (100.00%)
## 
## filter summary for frame 'Q1_T_t3_rep1.fcs'
##  Total Cells+: 34647 of 34647 events (100.00%)
## 
## filter summary for frame 'Q1_T_t3_rep2.fcs'
##  Total Cells+: 31325 of 31325 events (100.00%)
## 
## filter summary for frame 'Q2_S_t0.5_rep1.fcs'
##  Total Cells+: 42808 of 42808 events (100.00%)
## 
## filter summary for frame 'Q2_S_t0.5_rep2.fcs'
##  Total Cells+: 41557 of 41557 events (100.00%)
## 
## filter summary for frame 'Q2_S_t0.5_rep3.fcs'
##  Total Cells+: 43172 of 43172 events (100.00%)
## 
## filter summary for frame 'Q2_S_t0_rep1.fcs'
##  Total Cells+: 39581 of 39581 events (100.00%)
## 
## filter summary for frame 'Q2_S_t0_rep2.fcs'
##  Total Cells+: 38881 of 38881 events (100.00%)
## 
## filter summary for frame 'Q2_S_t0_rep3.fcs'
##  Total Cells+: 38546 of 38546 events (100.00%)
## 
## filter summary for frame 'Q2_S_t1.5_rep1.fcs'
##  Total Cells+: 35695 of 35695 events (100.00%)
## 
## filter summary for frame 'Q2_S_t1.5_rep2.fcs'
##  Total Cells+: 36598 of 36598 events (100.00%)
## 
## filter summary for frame 'Q2_S_t1.5_rep3.fcs'
##  Total Cells+: 35348 of 35348 events (100.00%)
## 
## filter summary for frame 'Q2_S_t1_rep1.fcs'
##  Total Cells+: 38506 of 38506 events (100.00%)
## 
## filter summary for frame 'Q2_S_t1_rep2.fcs'
##  Total Cells+: 38967 of 38967 events (100.00%)
## 
## filter summary for frame 'Q2_S_t1_rep3.fcs'
##  Total Cells+: 38451 of 38451 events (100.00%)
## 
## filter summary for frame 'Q2_S_t2.5_rep1.fcs'
##  Total Cells+: 36921 of 36921 events (100.00%)
## 
## filter summary for frame 'Q2_S_t2.5_rep2.fcs'
##  Total Cells+: 37198 of 37198 events (100.00%)
## 
## filter summary for frame 'Q2_S_t2.5_rep3.fcs'
##  Total Cells+: 35876 of 35876 events (100.00%)
## 
## filter summary for frame 'Q2_S_t2_rep1.fcs'
##  Total Cells+: 38345 of 38345 events (100.00%)
## 
## filter summary for frame 'Q2_S_t2_rep2.fcs'
##  Total Cells+: 37232 of 37232 events (100.00%)
## 
## filter summary for frame 'Q2_S_t2_rep3.fcs'
##  Total Cells+: 35964 of 35964 events (100.00%)
## 
## filter summary for frame 'Q2_S_t3_rep1.fcs'
##  Total Cells+: 34635 of 34635 events (100.00%)
## 
## filter summary for frame 'Q2_S_t3_rep2.fcs'
##  Total Cells+: 33967 of 33967 events (100.00%)
## 
## filter summary for frame 'Q2_S_t3_rep3.fcs'
##  Total Cells+: 34083 of 34083 events (100.00%)
## 
## filter summary for frame 'Q2_T_t0.5_rep1.fcs'
##  Total Cells+: 41673 of 41673 events (100.00%)
## 
## filter summary for frame 'Q2_T_t0.5_rep2.fcs'
##  Total Cells+: 41576 of 41576 events (100.00%)
## 
## filter summary for frame 'Q2_T_t0.5_rep3.fcs'
##  Total Cells+: 41430 of 41430 events (100.00%)
## 
## filter summary for frame 'Q2_T_t0_rep1.fcs'
##  Total Cells+: 43282 of 43282 events (100.00%)
## 
## filter summary for frame 'Q2_T_t0_rep2.fcs'
##  Total Cells+: 40883 of 40883 events (100.00%)
## 
## filter summary for frame 'Q2_T_t0_rep3.fcs'
##  Total Cells+: 44044 of 44044 events (100.00%)
## 
## filter summary for frame 'Q2_T_t1.5_rep1.fcs'
##  Total Cells+: 36219 of 36219 events (100.00%)
## 
## filter summary for frame 'Q2_T_t1.5_rep2.fcs'
##  Total Cells+: 35309 of 35309 events (100.00%)
## 
## filter summary for frame 'Q2_T_t1.5_rep3.fcs'
##  Total Cells+: 35596 of 35596 events (100.00%)
## 
## filter summary for frame 'Q2_T_t1_rep1.fcs'
##  Total Cells+: 38373 of 38373 events (100.00%)
## 
## filter summary for frame 'Q2_T_t1_rep2.fcs'
##  Total Cells+: 38487 of 38487 events (100.00%)
## 
## filter summary for frame 'Q2_T_t1_rep3.fcs'
##  Total Cells+: 38141 of 38141 events (100.00%)
## 
## filter summary for frame 'Q2_T_t2.5_rep1.fcs'
##  Total Cells+: 36953 of 36953 events (100.00%)
## 
## filter summary for frame 'Q2_T_t2.5_rep3.fcs'
##  Total Cells+: 37502 of 37502 events (100.00%)
## 
## filter summary for frame 'Q2_T_t2_rep1.fcs'
##  Total Cells+: 38297 of 38297 events (100.00%)
## 
## filter summary for frame 'Q2_T_t2_rep2.fcs'
##  Total Cells+: 36492 of 36492 events (100.00%)
## 
## filter summary for frame 'Q2_T_t2_rep3.fcs'
##  Total Cells+: 37242 of 37242 events (100.00%)
## 
## filter summary for frame 'Q2_T_t3_rep1.fcs'
##  Total Cells+: 33150 of 33150 events (100.00%)
## 
## filter summary for frame 'Q2_T_t3_rep2.fcs'
##  Total Cells+: 33737 of 33737 events (100.00%)
## 
## filter summary for frame 'Q2_T_t3_rep3.fcs'
##  Total Cells+: 33290 of 33290 events (100.00%)
## 
## filter summary for frame 'Q3_S_t0.5_rep1.fcs'
##  Total Cells+: 41447 of 41447 events (100.00%)
## 
## filter summary for frame 'Q3_S_t0.5_rep2.fcs'
##  Total Cells+: 40278 of 40278 events (100.00%)
## 
## filter summary for frame 'Q3_S_t0.5_rep3.fcs'
##  Total Cells+: 40841 of 40841 events (100.00%)
## 
## filter summary for frame 'Q3_S_t0_rep1.fcs'
##  Total Cells+: 38163 of 38163 events (100.00%)
## 
## filter summary for frame 'Q3_S_t0_rep3.fcs'
##  Total Cells+: 38368 of 38368 events (100.00%)
## 
## filter summary for frame 'Q3_S_t1.5_rep1.fcs'
##  Total Cells+: 35328 of 35328 events (100.00%)
## 
## filter summary for frame 'Q3_S_t1.5_rep2.fcs'
##  Total Cells+: 34111 of 34111 events (100.00%)
## 
## filter summary for frame 'Q3_S_t1.5_rep3.fcs'
##  Total Cells+: 35854 of 35854 events (100.00%)
## 
## filter summary for frame 'Q3_S_t1_rep1.fcs'
##  Total Cells+: 39785 of 39785 events (100.00%)
## 
## filter summary for frame 'Q3_S_t1_rep2.fcs'
##  Total Cells+: 39500 of 39500 events (100.00%)
## 
## filter summary for frame 'Q3_S_t1_rep3.fcs'
##  Total Cells+: 40081 of 40081 events (100.00%)
## 
## filter summary for frame 'Q3_S_t2.5_rep1.fcs'
##  Total Cells+: 34110 of 34110 events (100.00%)
## 
## filter summary for frame 'Q3_S_t2.5_rep2.fcs'
##  Total Cells+: 34462 of 34462 events (100.00%)
## 
## filter summary for frame 'Q3_S_t2.5_rep3.fcs'
##  Total Cells+: 33877 of 33877 events (100.00%)
## 
## filter summary for frame 'Q3_S_t2_rep1.fcs'
##  Total Cells+: 36134 of 36134 events (100.00%)
## 
## filter summary for frame 'Q3_S_t2_rep2.fcs'
##  Total Cells+: 35622 of 35622 events (100.00%)
## 
## filter summary for frame 'Q3_S_t2_rep3.fcs'
##  Total Cells+: 35652 of 35652 events (100.00%)
## 
## filter summary for frame 'Q3_S_t3_rep1.fcs'
##  Total Cells+: 35052 of 35052 events (100.00%)
## 
## filter summary for frame 'Q3_S_t3_rep2.fcs'
##  Total Cells+: 34025 of 34025 events (100.00%)
## 
## filter summary for frame 'Q3_S_t3_rep3.fcs'
##  Total Cells+: 34806 of 34806 events (100.00%)
## 
## filter summary for frame 'Q3_T_t0.5_rep1.fcs'
##  Total Cells+: 39768 of 39768 events (100.00%)
## 
## filter summary for frame 'Q3_T_t0.5_rep2.fcs'
##  Total Cells+: 40221 of 40221 events (100.00%)
## 
## filter summary for frame 'Q3_T_t0.5_rep3.fcs'
##  Total Cells+: 40525 of 40525 events (100.00%)
## 
## filter summary for frame 'Q3_T_t0_rep1.fcs'
##  Total Cells+: 40427 of 40427 events (100.00%)
## 
## filter summary for frame 'Q3_T_t0_rep2.fcs'
##  Total Cells+: 39238 of 39238 events (100.00%)
## 
## filter summary for frame 'Q3_T_t0_rep3.fcs'
##  Total Cells+: 39677 of 39677 events (100.00%)
## 
## filter summary for frame 'Q3_T_t1.5_rep1.fcs'
##  Total Cells+: 36149 of 36149 events (100.00%)
## 
## filter summary for frame 'Q3_T_t1.5_rep2.fcs'
##  Total Cells+: 35048 of 35048 events (100.00%)
## 
## filter summary for frame 'Q3_T_t1.5_rep3.fcs'
##  Total Cells+: 34974 of 34974 events (100.00%)
## 
## filter summary for frame 'Q3_T_t1_rep1.fcs'
##  Total Cells+: 39922 of 39922 events (100.00%)
## 
## filter summary for frame 'Q3_T_t1_rep2.fcs'
##  Total Cells+: 38560 of 38560 events (100.00%)
## 
## filter summary for frame 'Q3_T_t1_rep3.fcs'
##  Total Cells+: 39378 of 39378 events (100.00%)
## 
## filter summary for frame 'Q3_T_t2.5_rep1.fcs'
##  Total Cells+: 35912 of 35912 events (100.00%)
## 
## filter summary for frame 'Q3_T_t2.5_rep2.fcs'
##  Total Cells+: 32456 of 32456 events (100.00%)
## 
## filter summary for frame 'Q3_T_t2.5_rep3.fcs'
##  Total Cells+: 34403 of 34403 events (100.00%)
## 
## filter summary for frame 'Q3_T_t2_rep1.fcs'
##  Total Cells+: 38700 of 38700 events (100.00%)
## 
## filter summary for frame 'Q3_T_t2_rep2.fcs'
##  Total Cells+: 35818 of 35818 events (100.00%)
## 
## filter summary for frame 'Q3_T_t2_rep3.fcs'
##  Total Cells+: 38120 of 38120 events (100.00%)
## 
## filter summary for frame 'Q3_T_t3_rep1.fcs'
##  Total Cells+: 32594 of 32594 events (100.00%)
## 
## filter summary for frame 'Q3_T_t3_rep2.fcs'
##  Total Cells+: 32993 of 32993 events (100.00%)
## 
## filter summary for frame 'Q3_T_t3_rep3.fcs'
##  Total Cells+: 33277 of 33277 events (100.00%)
```

```r
# Extract the volume
vol <- c()
for(i in 1:length(flowData_transformed)){
  vol[i] <- as.numeric(flowData_transformed[[i]]@description$`$VOL`)/1000
}

# Save counts
counts <- data.frame(Samples=flowCore::sampleNames(flowData_transformed), 
                     Total.cells = dilution*TotalCount$true/vol, HNA.cells = dilution*HNACount$true/vol,
                     LNA.cells = dilution*(TotalCount$true-HNACount$true)/vol)

# Pool results
Diversity.fbasis <- Diversity.fbasis[Diversity.fbasis$D2>1600,]
Storage <- c(); Storage[Diversity.fbasis$Treatment=="T"] <- "Covered"; Storage[Diversity.fbasis$Treatment=="S"] <- "Submerged"
results <- cbind(Diversity.fbasis,counts, Storage, Diversity.HNA, Diversity.LNA)

# Add an extra column for biological replicates
bio_rep <- c(rep(1,21),rep(2,21),rep(1,41),rep(2,41),rep(3,41))
results <- cbind(results,bio_rep=factor(bio_rep))

# Select only S experiment for beta-div and further analysis

levels(results$Treatment)[levels(results$Treatment)=="C"] <- "Control"
levels(results$Treatment)[levels(results$Treatment)=="T"] <- "Feeding - T"
levels(results$Treatment)[levels(results$Treatment)=="S"] <- "Feeding - S"

# Beta diversity analysis
pos <- gsub(rownames(fbasis@basis),pattern="_rep1", replacement="")
pos <- gsub(pos,pattern="_rep2", replacement="")
pos <- gsub(pos,pattern="_rep3", replacement="")
pos <- factor(pos)
tmp <- results[,c(8,9,10,16,31)]
tmp <- by(tmp, INDICES=pos, FUN=unique)
tmp <- do.call(rbind, tmp)

# All samples
beta.div <- beta_div_fcm(fbasis, INDICES=pos, ord.type="PCoA")
fbasis1<-fbasis; fbasis1@basis <- fbasis@basis[!Diversity.fbasis$Treatment=="T",]
fbasis2<-fbasis; fbasis2@basis <- fbasis@basis[!Diversity.fbasis$Treatment=="S",]

# Without T and S versus control separately
beta.div.S <- beta_div_fcm(fbasis1, INDICES=pos[!Diversity.fbasis$Treatment=="T"], ord.type="PCoA")
# Take average fingerprint for technical replicates
fbasis1@basis <- fbasis1@basis/apply(fbasis1@basis, 1, max)
fbasis1@basis <- round(fbasis1@basis, 3)

# Remove shared zeroes for Bray-curtis
fbasis1@basis <- fbasis1@basis[, !apply(fbasis1@basis==0, 2, all)]
x <- by(fbasis1@basis, INDICES = pos[!Diversity.fbasis$Treatment=="T"], FUN = colMeans)
x <- do.call(rbind, x)


# Calculate distance matrix
dist.S <- vegdist(x, method="bray")

int_Treat_Time <- interaction(Diversity.fbasis$Treatment,Diversity.fbasis$Time)
dist.S.formatted <- as.matrix(dist.S)
dist.S_C_T0 <- dist.S.formatted[colnames(dist.S.formatted) == "C1_t0.fcs" | colnames(dist.S.formatted) == "C2_t0.fcs", colnames(dist.S.formatted) == "Q2_S_t0.fcs" | colnames(dist.S.formatted) == "Q1_S_t0.fcs"| colnames(dist.S.formatted) == "Q3_S_t3.fcs"]
mean(dist.S_C_T0)
```

```
## [1] 0.04973753
```

```r
sd(dist.S_C_T0)
```

```
## [1] 0.008512372
```

```r
dist.S_C_T3 <- dist.S.formatted[colnames(dist.S.formatted) == "C1_t3.fcs" | colnames(dist.S.formatted) == "C2_t3.fcs", colnames(dist.S.formatted) == "Q2_S_t3.fcs" | colnames(dist.S.formatted) == "Q1_S_t3.fcs"| colnames(dist.S.formatted) == "Q3_S_t3.fcs"]
mean(dist.S_C_T3)
```

```
## [1] 0.07092215
```

```r
sd(dist.S_C_T3)
```

```
## [1] 0.005427802
```

```r
beta.div.T <- beta_div_fcm(fbasis2, INDICES=pos[!Diversity.fbasis$Treatment=="S"], ord.type="PCoA")

results <- data.frame(cbind(Sample_names=levels(pos), do.call(rbind,by(results[,c(2,3,4,13,14,15,18,19,20,25,26,27)], INDICES=pos, FUN=colMeans)), 
                   do.call(rbind,by(results[,c(5,6,7,21,22,23,28,29,30)], INDICES=pos, FUN = function(x) sqrt(colSums(x^2))/nrow(x))),
                   do.call(rbind,by(results[,c(13,14,15)], INDICES=pos, FUN = function(x) apply(x,2,sd))),
                   tmp))

colnames(results)[colnames(results)=="Total.cells.1"|colnames(results)=="HNA.cells.1"|colnames(results)=="LNA.cells.1"] <- 
  c("sd.Total.cells","sd.HNA.cells","sd.LNA.cells")

colnames(results)[colnames(results)=="sd.D0.1"|colnames(results)=="sd.D1.1"|colnames(results)=="sd.D2.1"] <- 
  c("sd.D0.HNA","sd.D1.HNA","sd.D2.HNA")
colnames(results)[colnames(results)=="D0.1"|colnames(results)=="D1.1"|colnames(results)=="D2.1"] <- 
  c("D0.HNA","D1.HNA","D2.HNA")

colnames(results)[colnames(results)=="sd.D0.2"|colnames(results)=="sd.D1.2"|colnames(results)=="sd.D2.2"] <- 
  c("sd.D0.LNA","sd.D1.LNA","sd.D2.LNA")
colnames(results)[colnames(results)=="D0.2"|colnames(results)=="D1.2"|colnames(results)=="D2.2"] <- 
  c("D0.LNA","D1.LNA","D2.LNA")

# Import metadata of the mussels for M&M
meta.mus <- read.csv2("files/Experiment_metadata_mussels.csv")
meta.mus <- meta.mus[meta.mus$Treatment=="S",]
mean(meta.mus$size_mm)
```

```
## [1] 22.74111
```

```r
sd(meta.mus$size_mm)
```

```
## [1] 2.342459
```

```r
# No significant difference between mesocosms in mussel size distribution (p=0.081).
kruskal.test(size_mm ~ Sample, data=meta.mus)
```

```
## 
## 	Kruskal-Wallis rank sum test
## 
## data:  size_mm by Sample
## Kruskal-Wallis chi-squared = 5.0157, df = 2, p-value = 0.08145
```

```r
# Subset results for mussels transported by the best method (submerged)
result.tmp <- results
results <- result.tmp[!result.tmp$Treatment=="Feeding - T", ]
results$Treatment <- droplevels(plyr::revalue(results$Treatment, c("Feeding - S"="Feeding")))
```

## Diversity dynamics during IDM feeding


```r
# Calculate the dynamics of D2 in time for treatment/control
# We can't use a linear regression here (for obvious reasons)
# Hence we will try to fit robust smoothings splines and see if we can make inference this way
sp_T <- rlm(D2~splines::ns(Time, df=3), data=results[results$Treatment=="Feeding",])
sp_C <- rlm(D2~splines::ns(Time, df=3), data=results[results$Treatment=="Control",])

# Order studentized residuals according to timepoint
res.T <- residuals(sp_T, "pearson")[order(results[results$Treatment=="Feeding",]$Time)]
res.C <- residuals(sp_C, "pearson")[order(results[results$Treatment=="Control",]$Time)]

# Location of knots
attr(splines::ns(results$Time, df=3), "knots")
```

```
## 33.33333% 66.66667% 
##         1         2
```

## Figure S5: Check for autocorrelation  


```r
# Check for temporal autocorrelation in model residuals
par(mfrow=c(1,2))

acf(res.C, main="Treatment: control", las=1)
acf(res.T, main="Treatment: feeding", las=1)
```

<img src="Figures/cached/evaluate auto-correlation-1.png" style="display: block; margin: auto;" />

```r
# Perform statistical inference on splines
car::Anova(sp_C, vcov = vcovHAC(sp_C)) # p = 0.03795
```

```
## Analysis of Deviance Table (Type II tests)
## 
## Response: D2
##                           Df      F  Pr(>F)  
## splines::ns(Time, df = 3)  3 2.8069 0.09423 .
## Residuals                 10                 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
car::Anova(sp_T, vcov = vcovHAC(sp_T)) # p = 5.599e-06
```

```
## Analysis of Deviance Table (Type II tests)
## 
## Response: D2
##                           Df      F    Pr(>F)    
## splines::ns(Time, df = 3)  3 12.634 0.0001372 ***
## Residuals                 17                     
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
# Compare errors before and after correction
  # Treatment
summary(sp_T)$coefficients[,2]
```

```
##                (Intercept) splines::ns(Time, df = 3)1 
##                   13.57531                   22.46945 
## splines::ns(Time, df = 3)2 splines::ns(Time, df = 3)3 
##                   34.13373                   15.79189
```

```r
sqrt(diag(vcovHAC(sp_T)))
```

```
## [1] 10.92364 24.72847 29.67099 16.69775
```

```r
  # Control
summary(sp_C)$coefficients[,2]
```

```
##                (Intercept) splines::ns(Time, df = 3)1 
##                   12.05829                   19.95853 
## splines::ns(Time, df = 3)2 splines::ns(Time, df = 3)3 
##                   30.31935                   14.02717
```

```r
sqrt(diag(vcovHAC(sp_C)))
```

```
## [1]  8.001633 16.487676 25.016010 21.926641
```

## Figure 3: Diversity analysis during feeding  


```r
# Prepare plots
p.alpha <- ggplot(data=results, aes(x=Time, y=D2, fill=Treatment)) + 
  geom_point(shape=21, size=7,alpha=0.9)+
  scale_fill_manual(values=myColours[c(1,2)])+
  theme_bw()+
  labs(y=expression('Phenotypic diversity - D'[2]), x="Time (h)", title="A",
       fill="")+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=16),
        legend.direction = "horizontal",legend.position = "bottom")+ 
  geom_smooth(method="rlm", color="black", alpha=0.2, formula = y ~ splines::ns(x,df=3))+
  ylim(1990,2150)

# Separate for submerged data
# PERMANOVA on beta diversity analysis
disper.test <- betadisper(dist.S, group=results$Treatment)
disper.test # average distance to mean 0.03 for both groups
```

```
## 
## 	Homogeneity of multivariate dispersions
## 
## Call: betadisper(d = dist.S, group = results$Treatment)
## 
## No. of Positive Eigenvalues: 27
## No. of Negative Eigenvalues: 7
## 
## Average distance to median:
## Control Feeding 
## 0.03113 0.03057 
## 
## Eigenvalues for PCoA axes:
##  PCoA1  PCoA2  PCoA3  PCoA4  PCoA5  PCoA6  PCoA7  PCoA8 
## 0.0206 0.0093 0.0086 0.0037 0.0026 0.0013 0.0012 0.0008
```

```r
anova(disper.test) # P = 0.892
```

```
## Analysis of Variance Table
## 
## Response: Distances
##           Df    Sum Sq    Mean Sq F value Pr(>F)
## Groups     1 0.0000026 2.5890e-06  0.0189 0.8914
## Residuals 33 0.0045099 1.3666e-04
```

```r
perma.beta <- adonis(dist.S~Treatment*Time, data=results)
my_grob = grobTree(textGrob(bquote(paste(r[Feeding]^2 == .(round(100 * perma.beta$aov.tab[1, 
    5], 1)), "%")), x = 0.67, y = 0.95, hjust = 0, gp = gpar(col = "black", fontsize = 14, 
    fontface = "italic")))
my_grob2 = grobTree(textGrob(bquote(paste(r[Time]^2 == .(format(round(100 * perma.beta$aov.tab[2, 
    5], 1), nsmall = 1)), "%")), x = 0.67, y = 0.87, hjust = 0, gp = gpar(col = "black", 
    fontsize = 14, fontface = "italic")))
my_grob3 = grobTree(textGrob(bquote(paste(r[Feeding:Time]^2 == .(round(100 * perma.beta$aov.tab[3, 
    5], 1)), "%")), x = 0.67, y = 0.79, hjust = 0, gp = gpar(col = "black", fontsize = 14, 
    fontface = "italic")))

beta.div.data.S <- data.frame(beta.div.S$points, tmp[!tmp$Treatment=="Feeding - T",])
beta.div.data.S <- droplevels(beta.div.data.S)
beta.div.data.S$Treatment <- plyr::revalue(beta.div.data.S$Treatment, c("Feeding - S"="Feeding"))
var <- round(vegan::eigenvals(beta.div.S)/sum(vegan::eigenvals(beta.div.S))*100,1)

p.beta.S <- ggplot(data=beta.div.data.S, aes(x=X1, y=X2, fill=Treatment, size=Time))+
  geom_point(shape=21, alpha=1)+
  scale_size(range=c(4,10), breaks=c(0,0.5,1,1.5,2,2.5,3))+ 
  guides(fill = FALSE, size = guide_legend(nrow = 2))+
  theme_bw()+
  scale_fill_manual(values=myColours[c(1,2)])+
  labs(x = paste0("PCoA axis 1 (",var[1], "%)"), y = paste0("PCoA axis 2 (",var[2], "%)"), title="B",
       fill="", size = "Time (h)")+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=16),
        legend.direction = "horizontal",
        legend.position = c(0.15, -0.13), 
        legend.justification = c(0, 1),
        legend.background = element_rect(fill="transparent"))+
  annotation_custom(my_grob)+
  annotation_custom(my_grob2)+
  annotation_custom(my_grob3)+
  ylim(-0.035,0.045)

### Plot diversity dynamics
p.beta.S.plot <-  p.beta.S + theme(legend.position = "none")
p.alpha.plot <-  p.alpha + theme(legend.position = "none")

# pdf("Fig3_revised.pdf", width = 12, height = 6)
grid.arrange(p.alpha.plot, p.beta.S, g_legend(p.alpha), layout_matrix = rbind(c(1,2),
                                                                               c(3,3)),
             heights=c(6, 1))
```

<img src="Figures/cached/plot-diversity-1.png" style="display: block; margin: auto;" />

```r
# dev.off()
```

## Infer taxonomic diversity dynamics  


```r
# D2 regression from part I
lm.F <- lm(log2(D2)~log2(D2.fcm), data=data.total.final)

# Calculate mean predicted decrease in D2 at 0 and 3 hours
df <- data.frame(D2.fcm = results$D2[results$Treatment=="Feeding"][results[results$Treatment=="Feeding",]$Time==3 | results[results$Treatment=="Feeding",]$Time==0],
                 time=results$Time[results$Treatment=="Feeding"][results[results$Treatment=="Feeding",]$Time==3 | results[results$Treatment=="Feeding",]$Time==0])
df.pred <- data.frame(D2.16S = predict(lm.F, df, se=TRUE)$fit, D2.16S.error = predict(lm.F, df, se=TRUE)$se.fit, time = df$time)
df.pred <- data.frame(aggregate(D2.16S~time, df.pred, mean), aggregate(D2.16S.error~time, df.pred, FUN = function(x) sqrt(sum(x^2))/length(x)))

# Note: After transformation back to linear scale, errors are not necessarily normal distributed anymore, so calculate minimum and maximum and report as such.
df.pred.res <- df.pred 
# Maximum decrease in D2 (3.63 units or 15.6 %)
2^(df.pred$D2.16S[df.pred.res$time==0] + df.pred$D2.16S.error[df.pred.res$time==0]) -  2^(df.pred$D2.16S[df.pred.res$time==3] - df.pred$D2.16S.error[df.pred.res$time==3])
```

```
## [1] 3.621665
```

```r
(2^(df.pred$D2.16S[df.pred.res$time==0] + df.pred$D2.16S.error[df.pred.res$time==3]) -  2^(df.pred$D2.16S[df.pred.res$time==3] - df.pred$D2.16S.error[df.pred.res$time==3])) / 2^(df.pred$D2.16S[df.pred.res$time==0] + df.pred$D2.16S.error[df.pred.res$time==0])
```

```
## [1] 0.1546826
```

```r
# Minimum decrease in D2 (1.65 units or 7.5 %)
2^(df.pred$D2.16S[df.pred.res$time==0] - df.pred$D2.16S.error[df.pred.res$time==0]) -  2^(df.pred$D2.16S[df.pred.res$time==3] + df.pred$D2.16S.error[df.pred.res$time==3])
```

```
## [1] 1.645414
```

```r
(2^(df.pred$D2.16S[df.pred.res$time==0] - df.pred$D2.16S.error[df.pred.res$time==3]) -  2^(df.pred$D2.16S[df.pred.res$time==3] + df.pred$D2.16S.error[df.pred.res$time==3])) / 2^(df.pred$D2.16S[df.pred.res$time==0] - df.pred$D2.16S.error[df.pred.res$time==0])
```

```
## [1] 0.07544289
```

```r
# Average decrease in D2 (2.64 units or 11.6 %)
2^df.pred$D2.16S[df.pred.res$time==0] -  2^df.pred$D2.16S[df.pred.res$time==3]
```

```
## [1] 2.63225
```

```r
(2^df.pred$D2.16S[df.pred.res$time==0] -  2^df.pred$D2.16S[df.pred.res$time==3]) / 2^df.pred$D2.16S[df.pred.res$time==0] 
```

```
## [1] 0.1160036
```

```r
# Note: Errors seem (approx.) randomly distributed even on linear scale, so we can just summarize by +/- the average on the mean estimate.

# Calculate mean predicted decrease in D2 at 0 and 1 hour
df <- data.frame(D2.fcm = results$D2[results$Treatment=="Feeding"][results[results$Treatment=="Feeding",]$Time==1 | results[results$Treatment=="Feeding",]$Time==0],
                 time=results$Time[results$Treatment=="Feeding"][results[results$Treatment=="Feeding",]$Time==1 | results[results$Treatment=="Feeding",]$Time==0])
df.pred <- data.frame(D2.16S = predict(lm.F, df, se=TRUE)$fit, D2.16S.error = predict(lm.F, df, se=TRUE)$se.fit, time = df$time)
df.pred <- data.frame(aggregate(D2.16S~time, df.pred, mean), aggregate(D2.16S.error~time, df.pred, FUN = function(x) sqrt(sum(x^2))/length(x)))

# Note: After transformation back to linear scale, errors are not necessarily normal distributed anymore, so calculate minimum and maximum and report as such.

df.pred.res <- df.pred 
# Maximum decrease in D2 (2.82 units or 12.05 %)
2^(df.pred$D2.16S[df.pred.res$time==0] + df.pred$D2.16S.error[df.pred.res$time==0]) -  2^(df.pred$D2.16S[df.pred.res$time==1] - df.pred$D2.16S.error[df.pred.res$time==1])
```

```
## [1] 2.838505
```

```r
(2^(df.pred$D2.16S[df.pred.res$time==0] + df.pred$D2.16S.error[df.pred.res$time==1]) -  2^(df.pred$D2.16S[df.pred.res$time==1] - df.pred$D2.16S.error[df.pred.res$time==1])) / 2^(df.pred$D2.16S[df.pred.res$time==0] + df.pred$D2.16S.error[df.pred.res$time==0])
```

```
## [1] 0.1213341
```

```r
# Minimum decrease in D2 (0.79 units or 3.66 %)
2^(df.pred$D2.16S[df.pred.res$time==0] - df.pred$D2.16S.error[df.pred.res$time==0]) -  2^(df.pred$D2.16S[df.pred.res$time==1] + df.pred$D2.16S.error[df.pred.res$time==1])
```

```
## [1] 0.8109439
```

```r
(2^(df.pred$D2.16S[df.pred.res$time==0] - df.pred$D2.16S.error[df.pred.res$time==1]) -  2^(df.pred$D2.16S[df.pred.res$time==1] + df.pred$D2.16S.error[df.pred.res$time==1])) / 2^(df.pred$D2.16S[df.pred.res$time==0] - df.pred$D2.16S.error[df.pred.res$time==0])
```

```
## [1] 0.03742805
```

```r
# Average decrease in D2 (1.8 units or 7.96 %)
2^df.pred$D2.16S[df.pred.res$time==0] -  2^df.pred$D2.16S[df.pred.res$time==1]
```

```
## [1] 1.823809
```

```r
(2^df.pred$D2.16S[df.pred.res$time==0] -  2^df.pred$D2.16S[df.pred.res$time==1]) / 2^df.pred$D2.16S[df.pred.res$time==0] 
```

```
## [1] 0.08037551
```

```r
# Calculate delta in D2.fcm for an increase of 10 taxonomic units (5->10 and 35->45)
df2 <- data.frame(D2.fcm=c(1000,1500,2000,2500))
df2.res <- predict(lm.F,df2)
# at low diversities
2^(df2.res[2] - df2.res[1])
```

```
##       2 
## 3.35925
```

```r
# at higher diversities
2^(df2.res[4] - df2.res[3])
```

```
##        4 
## 1.948104
```

## Compare diversity dynamics to literature (DOI: 10.1128/mSphere.00189-17)  


```r
physeq_msphere <- import_mothur(mothur_shared_file = "./data_published/QMexp_IL_LM.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.shared", 
                                mothur_constaxonomy_file = "./data_published/QMexp_IL_LM.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.0.03.cons.taxonomy")
meta_msphere <- read.table("./data_published/QMexp_IL_LM_womockblanks.tsv", header = TRUE)
rownames(meta_msphere) <- meta_msphere$Sample_ID
sample_data(physeq_msphere) <- sample_data(meta_msphere)
# Subset for lab feeding experiment data
physeq_msphere <- subset_samples(physeq_msphere, Group == "Experiment")

# Calculate diversity before and after rescaling
div_msphere <- Diversity_16S(physeq_msphere, R=100, parallel = TRUE, ncores = 10)
```

```
## 	**WARNING** this functions assumes that rows are samples and columns
##       	are taxa in your phyloseq object, please verify.
## Fri Jan 12 13:48:19 2018 	Using 10 cores for calculations
## Fri Jan 12 13:48:19 2018	Calculating diversity for sample 1/48 --- C1T0.1214
## Fri Jan 12 13:48:30 2018	Done with sample C1T0.1214
## Fri Jan 12 13:48:30 2018	Calculating diversity for sample 2/48 --- C1T0.813
## Fri Jan 12 13:48:35 2018	Done with sample C1T0.813
## Fri Jan 12 13:48:35 2018	Calculating diversity for sample 3/48 --- C1T0.814
## Fri Jan 12 13:48:39 2018	Done with sample C1T0.814
## Fri Jan 12 13:48:39 2018	Calculating diversity for sample 4/48 --- C1T3.1214
## Fri Jan 12 13:48:42 2018	Done with sample C1T3.1214
## Fri Jan 12 13:48:42 2018	Calculating diversity for sample 5/48 --- C1T3.813
## Fri Jan 12 13:48:45 2018	Done with sample C1T3.813
## Fri Jan 12 13:48:45 2018	Calculating diversity for sample 6/48 --- C1T3.814
## Fri Jan 12 13:48:49 2018	Done with sample C1T3.814
## Fri Jan 12 13:48:49 2018	Calculating diversity for sample 7/48 --- C2T0.1214
## Fri Jan 12 13:48:52 2018	Done with sample C2T0.1214
## Fri Jan 12 13:48:52 2018	Calculating diversity for sample 8/48 --- C2T0.813
## Fri Jan 12 13:48:55 2018	Done with sample C2T0.813
## Fri Jan 12 13:48:55 2018	Calculating diversity for sample 9/48 --- C2T0.814
## Fri Jan 12 13:48:57 2018	Done with sample C2T0.814
## Fri Jan 12 13:48:57 2018	Calculating diversity for sample 10/48 --- C2T3.1214
## Fri Jan 12 13:49:01 2018	Done with sample C2T3.1214
## Fri Jan 12 13:49:01 2018	Calculating diversity for sample 11/48 --- C2T3.813
## Fri Jan 12 13:49:04 2018	Done with sample C2T3.813
## Fri Jan 12 13:49:04 2018	Calculating diversity for sample 12/48 --- C2T3.814
## Fri Jan 12 13:49:10 2018	Done with sample C2T3.814
## Fri Jan 12 13:49:10 2018	Calculating diversity for sample 13/48 --- C3T0.1214
## Fri Jan 12 13:49:15 2018	Done with sample C3T0.1214
## Fri Jan 12 13:49:15 2018	Calculating diversity for sample 14/48 --- C3T0.813
## Fri Jan 12 13:49:18 2018	Done with sample C3T0.813
## Fri Jan 12 13:49:18 2018	Calculating diversity for sample 15/48 --- C3T0.814
## Fri Jan 12 13:49:21 2018	Done with sample C3T0.814
## Fri Jan 12 13:49:21 2018	Calculating diversity for sample 16/48 --- C3T3.1214
## Fri Jan 12 13:49:24 2018	Done with sample C3T3.1214
## Fri Jan 12 13:49:24 2018	Calculating diversity for sample 17/48 --- C3T3.813
## Fri Jan 12 13:49:27 2018	Done with sample C3T3.813
## Fri Jan 12 13:49:27 2018	Calculating diversity for sample 18/48 --- C3T3.814
## Fri Jan 12 13:49:30 2018	Done with sample C3T3.814
## Fri Jan 12 13:49:30 2018	Calculating diversity for sample 19/48 --- E1T0.1214
## Fri Jan 12 13:49:36 2018	Done with sample E1T0.1214
## Fri Jan 12 13:49:36 2018	Calculating diversity for sample 20/48 --- E1T0.813
## Fri Jan 12 13:49:41 2018	Done with sample E1T0.813
## Fri Jan 12 13:49:41 2018	Calculating diversity for sample 21/48 --- E1T0.814
## Fri Jan 12 13:49:45 2018	Done with sample E1T0.814
## Fri Jan 12 13:49:45 2018	Calculating diversity for sample 22/48 --- E1T3.1214
## Fri Jan 12 13:49:50 2018	Done with sample E1T3.1214
## Fri Jan 12 13:49:50 2018	Calculating diversity for sample 23/48 --- E1T3.813
## Fri Jan 12 13:49:53 2018	Done with sample E1T3.813
## Fri Jan 12 13:49:53 2018	Calculating diversity for sample 24/48 --- E1T3.814
## Fri Jan 12 13:49:57 2018	Done with sample E1T3.814
## Fri Jan 12 13:49:57 2018	Calculating diversity for sample 25/48 --- E2T0.1214
## Fri Jan 12 13:50:00 2018	Done with sample E2T0.1214
## Fri Jan 12 13:50:00 2018	Calculating diversity for sample 26/48 --- E2T0.813
## Fri Jan 12 13:50:03 2018	Done with sample E2T0.813
## Fri Jan 12 13:50:03 2018	Calculating diversity for sample 27/48 --- E2T0.814
## Fri Jan 12 13:50:06 2018	Done with sample E2T0.814
## Fri Jan 12 13:50:06 2018	Calculating diversity for sample 28/48 --- E2T3.1214
## Fri Jan 12 13:50:10 2018	Done with sample E2T3.1214
## Fri Jan 12 13:50:10 2018	Calculating diversity for sample 29/48 --- E2T3.813
## Fri Jan 12 13:50:13 2018	Done with sample E2T3.813
## Fri Jan 12 13:50:13 2018	Calculating diversity for sample 30/48 --- E2T3.814
## Fri Jan 12 13:50:16 2018	Done with sample E2T3.814
## Fri Jan 12 13:50:16 2018	Calculating diversity for sample 31/48 --- E3T0.1214
## Fri Jan 12 13:50:19 2018	Done with sample E3T0.1214
## Fri Jan 12 13:50:19 2018	Calculating diversity for sample 32/48 --- E3T0.813
## Fri Jan 12 13:50:22 2018	Done with sample E3T0.813
## Fri Jan 12 13:50:22 2018	Calculating diversity for sample 33/48 --- E3T0.814
## Fri Jan 12 13:50:25 2018	Done with sample E3T0.814
## Fri Jan 12 13:50:25 2018	Calculating diversity for sample 34/48 --- E3T3.1214
## Fri Jan 12 13:50:28 2018	Done with sample E3T3.1214
## Fri Jan 12 13:50:28 2018	Calculating diversity for sample 35/48 --- E3T3.813
## Fri Jan 12 13:50:31 2018	Done with sample E3T3.813
## Fri Jan 12 13:50:31 2018	Calculating diversity for sample 36/48 --- E3T3.814
## Fri Jan 12 13:50:35 2018	Done with sample E3T3.814
## Fri Jan 12 13:50:35 2018	Calculating diversity for sample 37/48 --- E4T0.1214
## Fri Jan 12 13:50:39 2018	Done with sample E4T0.1214
## Fri Jan 12 13:50:39 2018	Calculating diversity for sample 38/48 --- E4T0.813
## Fri Jan 12 13:50:41 2018	Done with sample E4T0.813
## Fri Jan 12 13:50:41 2018	Calculating diversity for sample 39/48 --- E4T0.814
## Fri Jan 12 13:50:44 2018	Done with sample E4T0.814
## Fri Jan 12 13:50:44 2018	Calculating diversity for sample 40/48 --- E4T3.1214
## Fri Jan 12 13:50:47 2018	Done with sample E4T3.1214
## Fri Jan 12 13:50:47 2018	Calculating diversity for sample 41/48 --- E4T3.813
## Fri Jan 12 13:50:50 2018	Done with sample E4T3.813
## Fri Jan 12 13:50:50 2018	Calculating diversity for sample 42/48 --- E4T3.814
## Fri Jan 12 13:50:53 2018	Done with sample E4T3.814
## Fri Jan 12 13:50:53 2018	Calculating diversity for sample 43/48 --- W1FR.814
## Fri Jan 12 13:50:57 2018	Done with sample W1FR.814
## Fri Jan 12 13:50:57 2018	Calculating diversity for sample 44/48 --- W1WH.814
## Fri Jan 12 13:51:00 2018	Done with sample W1WH.814
## Fri Jan 12 13:51:00 2018	Calculating diversity for sample 45/48 --- W2FR.814
## Fri Jan 12 13:51:03 2018	Done with sample W2FR.814
## Fri Jan 12 13:51:03 2018	Calculating diversity for sample 46/48 --- W2WH.814
## Fri Jan 12 13:51:06 2018	Done with sample W2WH.814
## Fri Jan 12 13:51:06 2018	Calculating diversity for sample 47/48 --- W3FR.814
## Fri Jan 12 13:51:09 2018	Done with sample W3FR.814
## Fri Jan 12 13:51:09 2018	Calculating diversity for sample 48/48 --- W3WH.814
## Fri Jan 12 13:51:12 2018	Done with sample W3WH.814
## Fri Jan 12 13:51:12 2018 	Closing workers
## Fri Jan 12 13:51:12 2018 	Done with all 48 samples
```

```r
div_msphere_scaled <- Diversity_16S(scale_reads(physeq_msphere), R=100, parallel = TRUE, ncores = 10)
```

```
## 	**WARNING** this functions assumes that rows are samples and columns
##       	are taxa in your phyloseq object, please verify.
## Fri Jan 12 13:51:16 2018 	Using 10 cores for calculations
## Fri Jan 12 13:51:16 2018	Calculating diversity for sample 1/48 --- C1T0.1214
## Fri Jan 12 13:51:25 2018	Done with sample C1T0.1214
## Fri Jan 12 13:51:25 2018	Calculating diversity for sample 2/48 --- C1T0.813
## Fri Jan 12 13:51:26 2018	Done with sample C1T0.813
## Fri Jan 12 13:51:26 2018	Calculating diversity for sample 3/48 --- C1T0.814
## Fri Jan 12 13:51:28 2018	Done with sample C1T0.814
## Fri Jan 12 13:51:28 2018	Calculating diversity for sample 4/48 --- C1T3.1214
## Fri Jan 12 13:51:29 2018	Done with sample C1T3.1214
## Fri Jan 12 13:51:29 2018	Calculating diversity for sample 5/48 --- C1T3.813
## Fri Jan 12 13:51:30 2018	Done with sample C1T3.813
## Fri Jan 12 13:51:30 2018	Calculating diversity for sample 6/48 --- C1T3.814
## Fri Jan 12 13:51:31 2018	Done with sample C1T3.814
## Fri Jan 12 13:51:31 2018	Calculating diversity for sample 7/48 --- C2T0.1214
## Fri Jan 12 13:51:33 2018	Done with sample C2T0.1214
## Fri Jan 12 13:51:33 2018	Calculating diversity for sample 8/48 --- C2T0.813
## Fri Jan 12 13:51:34 2018	Done with sample C2T0.813
## Fri Jan 12 13:51:34 2018	Calculating diversity for sample 9/48 --- C2T0.814
## Fri Jan 12 13:51:35 2018	Done with sample C2T0.814
## Fri Jan 12 13:51:35 2018	Calculating diversity for sample 10/48 --- C2T3.1214
## Fri Jan 12 13:51:37 2018	Done with sample C2T3.1214
## Fri Jan 12 13:51:37 2018	Calculating diversity for sample 11/48 --- C2T3.813
## Fri Jan 12 13:51:38 2018	Done with sample C2T3.813
## Fri Jan 12 13:51:38 2018	Calculating diversity for sample 12/48 --- C2T3.814
## Fri Jan 12 13:51:41 2018	Done with sample C2T3.814
## Fri Jan 12 13:51:41 2018	Calculating diversity for sample 13/48 --- C3T0.1214
## Fri Jan 12 13:51:42 2018	Done with sample C3T0.1214
## Fri Jan 12 13:51:43 2018	Calculating diversity for sample 14/48 --- C3T0.813
## Fri Jan 12 13:51:44 2018	Done with sample C3T0.813
## Fri Jan 12 13:51:44 2018	Calculating diversity for sample 15/48 --- C3T0.814
## Fri Jan 12 13:51:45 2018	Done with sample C3T0.814
## Fri Jan 12 13:51:45 2018	Calculating diversity for sample 16/48 --- C3T3.1214
## Fri Jan 12 13:51:48 2018	Done with sample C3T3.1214
## Fri Jan 12 13:51:48 2018	Calculating diversity for sample 17/48 --- C3T3.813
## Fri Jan 12 13:51:49 2018	Done with sample C3T3.813
## Fri Jan 12 13:51:49 2018	Calculating diversity for sample 18/48 --- C3T3.814
## Fri Jan 12 13:51:51 2018	Done with sample C3T3.814
## Fri Jan 12 13:51:51 2018	Calculating diversity for sample 19/48 --- E1T0.1214
## Fri Jan 12 13:51:54 2018	Done with sample E1T0.1214
## Fri Jan 12 13:51:54 2018	Calculating diversity for sample 20/48 --- E1T0.813
## Fri Jan 12 13:51:55 2018	Done with sample E1T0.813
## Fri Jan 12 13:51:55 2018	Calculating diversity for sample 21/48 --- E1T0.814
## Fri Jan 12 13:51:56 2018	Done with sample E1T0.814
## Fri Jan 12 13:51:56 2018	Calculating diversity for sample 22/48 --- E1T3.1214
## Fri Jan 12 13:51:58 2018	Done with sample E1T3.1214
## Fri Jan 12 13:51:58 2018	Calculating diversity for sample 23/48 --- E1T3.813
## Fri Jan 12 13:51:59 2018	Done with sample E1T3.813
## Fri Jan 12 13:51:59 2018	Calculating diversity for sample 24/48 --- E1T3.814
## Fri Jan 12 13:52:00 2018	Done with sample E1T3.814
## Fri Jan 12 13:52:00 2018	Calculating diversity for sample 25/48 --- E2T0.1214
## Fri Jan 12 13:52:01 2018	Done with sample E2T0.1214
## Fri Jan 12 13:52:01 2018	Calculating diversity for sample 26/48 --- E2T0.813
## Fri Jan 12 13:52:02 2018	Done with sample E2T0.813
## Fri Jan 12 13:52:02 2018	Calculating diversity for sample 27/48 --- E2T0.814
## Fri Jan 12 13:52:04 2018	Done with sample E2T0.814
## Fri Jan 12 13:52:04 2018	Calculating diversity for sample 28/48 --- E2T3.1214
## Fri Jan 12 13:52:06 2018	Done with sample E2T3.1214
## Fri Jan 12 13:52:06 2018	Calculating diversity for sample 29/48 --- E2T3.813
## Fri Jan 12 13:52:07 2018	Done with sample E2T3.813
## Fri Jan 12 13:52:07 2018	Calculating diversity for sample 30/48 --- E2T3.814
## Fri Jan 12 13:52:10 2018	Done with sample E2T3.814
## Fri Jan 12 13:52:10 2018	Calculating diversity for sample 31/48 --- E3T0.1214
## Fri Jan 12 13:52:12 2018	Done with sample E3T0.1214
## Fri Jan 12 13:52:12 2018	Calculating diversity for sample 32/48 --- E3T0.813
## Fri Jan 12 13:52:13 2018	Done with sample E3T0.813
## Fri Jan 12 13:52:13 2018	Calculating diversity for sample 33/48 --- E3T0.814
## Fri Jan 12 13:52:14 2018	Done with sample E3T0.814
## Fri Jan 12 13:52:14 2018	Calculating diversity for sample 34/48 --- E3T3.1214
## Fri Jan 12 13:52:15 2018	Done with sample E3T3.1214
## Fri Jan 12 13:52:16 2018	Calculating diversity for sample 35/48 --- E3T3.813
## Fri Jan 12 13:52:17 2018	Done with sample E3T3.813
## Fri Jan 12 13:52:17 2018	Calculating diversity for sample 36/48 --- E3T3.814
## Fri Jan 12 13:52:18 2018	Done with sample E3T3.814
## Fri Jan 12 13:52:18 2018	Calculating diversity for sample 37/48 --- E4T0.1214
## Fri Jan 12 13:52:19 2018	Done with sample E4T0.1214
## Fri Jan 12 13:52:19 2018	Calculating diversity for sample 38/48 --- E4T0.813
## Fri Jan 12 13:52:20 2018	Done with sample E4T0.813
## Fri Jan 12 13:52:20 2018	Calculating diversity for sample 39/48 --- E4T0.814
## Fri Jan 12 13:52:22 2018	Done with sample E4T0.814
## Fri Jan 12 13:52:22 2018	Calculating diversity for sample 40/48 --- E4T3.1214
## Fri Jan 12 13:52:23 2018	Done with sample E4T3.1214
## Fri Jan 12 13:52:23 2018	Calculating diversity for sample 41/48 --- E4T3.813
## Fri Jan 12 13:52:24 2018	Done with sample E4T3.813
## Fri Jan 12 13:52:24 2018	Calculating diversity for sample 42/48 --- E4T3.814
## Fri Jan 12 13:52:27 2018	Done with sample E4T3.814
## Fri Jan 12 13:52:27 2018	Calculating diversity for sample 43/48 --- W1FR.814
## Fri Jan 12 13:52:28 2018	Done with sample W1FR.814
## Fri Jan 12 13:52:28 2018	Calculating diversity for sample 44/48 --- W1WH.814
## Fri Jan 12 13:52:32 2018	Done with sample W1WH.814
## Fri Jan 12 13:52:32 2018	Calculating diversity for sample 45/48 --- W2FR.814
## Fri Jan 12 13:52:34 2018	Done with sample W2FR.814
## Fri Jan 12 13:52:34 2018	Calculating diversity for sample 46/48 --- W2WH.814
## Fri Jan 12 13:52:36 2018	Done with sample W2WH.814
## Fri Jan 12 13:52:36 2018	Calculating diversity for sample 47/48 --- W3FR.814
## Fri Jan 12 13:52:38 2018	Done with sample W3FR.814
## Fri Jan 12 13:52:38 2018	Calculating diversity for sample 48/48 --- W3WH.814
## Fri Jan 12 13:52:39 2018	Done with sample W3WH.814
## Fri Jan 12 13:52:39 2018 	Closing workers
## Fri Jan 12 13:52:39 2018 	Done with all 48 samples
```

```r
div_msphere_scaled <- data.frame(div_msphere_scaled, Sample_ID = rownames(div_msphere_scaled))
div_msphere <- data.frame(div_msphere, Sample_ID = rownames(div_msphere))

# Compare pre- and postfeeding diversities
results_msphere <- left_join(div_msphere, meta_msphere, by = "Sample_ID")
```

```
## Warning: Column `Sample_ID` joining factors with different levels, coercing
## to character vector
```

```r
results_msphere_scale <- left_join(div_msphere_scaled, meta_msphere, by = "Sample_ID")
```

```
## Warning: Column `Sample_ID` joining factors with different levels, coercing
## to character vector
```

```r
results_msphere <- results_msphere[results_msphere$Treatment == "Control" | results_msphere$Treatment == "Mussel", ]
results_msphere$Treatment <- revalue(results_msphere$Treatment, c("Mussel" = "Feeding"))
results_msphere <- droplevels(results_msphere)
results_msphere$Date <- factor(results_msphere$Date,levels(results_msphere$Date)[c(3,1,2)])
results_msphere_scale <- results_msphere_scale[results_msphere_scale$Treatment == "Control" | results_msphere_scale$Treatment == "Mussel", ]
results_msphere_scale$Treatment <- revalue(results_msphere_scale$Treatment, c("Mussel" = "Feeding"))
results_msphere_scale <- droplevels(results_msphere_scale)
results_msphere_scale$Date <- factor(results_msphere_scale$Date,levels(results_msphere_scale$Date)[c(3,1,2)])

# Plot diversity dynamics from literature study
p_msphere <- ggplot(data = results_msphere, aes(x = Time_point, y = D2, fill = Treatment))+
  geom_boxplot(alpha = 0.5)+
  geom_point(size = 4, shape = 21, position=position_dodge(width=0.75))+
  theme_bw()+
  facet_grid(.~Date)+
  scale_fill_manual(values=myColours[c(1:3)])+
  labs(y = expression('Taxonomic diversity - D'[2]), x = "Timepoint")+
  theme(axis.title=element_text(size=16), strip.text.x=element_text(size=16),
        legend.title=element_text(size=15),legend.text=element_text(size=14),
        axis.text = element_text(size=14),title=element_text(size=20),
        strip.background=element_rect(fill=adjustcolor("lightgray",0.2))
        #,panel.grid.major = element_blank(), panel.grid.minor = element_blank()
        )
# png(file="msphere_data.png",width=12,height=6,res=500,units="in", pointsize=12)
p_msphere
```

<img src="Figures/cached/compare-diversity-dynamics2MSphere-1.png" style="display: block; margin: auto;" />

```r
# dev.off()

p_msphere_scale <- ggplot(data = results_msphere_scale, aes(x = Time_point, y = D2, fill = Treatment))+
  geom_boxplot(alpha = 0.5)+
  geom_point(size = 4, shape = 21, position=position_dodge(width=0.75))+
  theme_bw()+
  facet_grid(.~Date)+
  scale_fill_manual(values=myColours[c(1:3)])+
  labs(y = expression('Taxonomic diversity - D'[2]), x = "Timepoint")+
  theme(axis.title=element_text(size=16), strip.text.x=element_text(size=16),
        legend.title=element_text(size=15),legend.text=element_text(size=14),
        axis.text = element_text(size=14),title=element_text(size=20),
        strip.background=element_rect(fill=adjustcolor("lightgray",0.2))
        #,panel.grid.major = element_blank(), panel.grid.minor = element_blank()
        )
p_msphere_scale
```

<img src="Figures/cached/compare-diversity-dynamics2MSphere-2.png" style="display: block; margin: auto;" />

```r
# Continuing with the non-scaled diversity we calculate the average diversity reduction over all experiments
## Mean reduction
mean(results_msphere$D2[results_msphere$Treatment == "Feeding" & results_msphere$Time_point == "T0"] - results_msphere$D2[results_msphere$Treatment == "Feeding" & results_msphere$Time_point == "T3"]) # 5.32
```

```
## [1] 5.314181
```

```r
## Standard deviation
sd(results_msphere$D2[results_msphere$Treatment == "Feeding" & results_msphere$Time_point == "T0"] - results_msphere$D2[results_msphere$Treatment == "Feeding" & results_msphere$Time_point == "T3"]) # 4.65
```

```
## [1] 4.684461
```

```r
# Same calculation but without the December experiment
## Mean reduction
mean(results_msphere$D2[results_msphere$Treatment == "Feeding" & results_msphere$Time_point == "T0" & results_msphere$Date != "Dec14"] - results_msphere$D2[results_msphere$Treatment == "Feeding" & results_msphere$Time_point == "T3" & results_msphere$Date != "Dec14"]) # 7.46
```

```
## [1] 7.467301
```

```r
## Standard deviation
sd(results_msphere$D2[results_msphere$Treatment == "Feeding" & results_msphere$Time_point == "T0" & results_msphere$Date != "Dec14"] - results_msphere$D2[results_msphere$Treatment == "Feeding" & results_msphere$Time_point == "T3" & results_msphere$Date != "Dec14"]) # 4.13
```

```
## [1] 4.170674
```


## Figure 4: Predict diversity and put in juxtaposition of literature (DOI: 10.1128/mSphere.00189-17)  


```r
# Predict the shift in taxonomic diversity based on FCM for the T0 and T3 samples
df_mp <- data.frame(D2.fcm = results$D2[results$Time==3 | results$Time==0],
                 time= results$Time[results$Time==3 | results$Time==0],
                 Treatment = results$Treatment[results$Time==3 | results$Time==0])
df_mp_pred <- data.frame(D2.16S = predict(lm.F, df_mp, se=TRUE)$fit, D2.16S.error = predict(lm.F, df_mp, se=TRUE)$se.fit, time = df_mp$time, treatment = df_mp$Treatment)
df_mp_pred$D2.16S.up <- df_mp_pred$D2.16S + df_mp_pred$D2.16S.error
df_mp_pred$D2.16S.lo <- df_mp_pred$D2.16S - df_mp_pred$D2.16S.error
df_mp_pred$D2.16S <- 2^df_mp_pred$D2.16S
df_mp_pred$D2.16S.up <- 2^df_mp_pred$D2.16S.up
df_mp_pred$D2.16S.lo <- 2^df_mp_pred$D2.16S.lo
df_mp_pred$time[df_mp_pred$time == 0] <- "T0"
df_mp_pred$time[df_mp_pred$time == 3] <- "T3"
df_mp_pred$Date <- "November '16"
levels(results_msphere$Date) <- c("July '13", "August '14", "December '14")
  
# Merge these predicted values into the data set from (DOI: 10.1128/mSphere.00189-17)
results_prediction <- data.frame(D2 = c(results_msphere$D2, df_mp_pred$D2.16S),
                                 Time_point = c(as.character(results_msphere$Time_point), df_mp_pred$time),
                                 Treatment = c(as.character(results_msphere$Treatment), as.character(df_mp_pred$treatment)),
                                 Measurement = as.factor(c(rep("Measured",42), rep("Predicted", 10))),
                                 Date = c(as.character(results_msphere$Date), df_mp_pred$Date),
                                 error.up = c(rep(NA, nrow(results_msphere)), df_mp_pred$D2.16S.up),
                                 error.lo = c(rep(NA, nrow(results_msphere)), df_mp_pred$D2.16S.lo)
                                 
                                 )
results_prediction$Date <- factor(results_prediction$Date, levels = c("July '13", "August '14", "November '16", "December '14"))

# Create a data frame with the faceting variables
# and some dummy data (that will be overwritten)
tp <- results_prediction

# Make Figure 5: compare predictions to literature data
p_msphere_predict <- ggplot(data = results_prediction, aes(x = Time_point, y = D2))+
  geom_rect(data = tp, aes(fill = Measurement), xmin = -Inf, xmax = Inf,
            ymin = -Inf,ymax = Inf, alpha = 0.005, show.legend =FALSE, inherit.aes = FALSE)+
  geom_boxplot(alpha = 0.5, aes(fill = Treatment), position = position_dodge(width = .8), width =0.75)+
  geom_point(size = 4, shape = 21, position = position_dodge(width=0.8), aes(fill = Treatment))+
  geom_errorbar(aes(ymin = error.lo, ymax = error.up, group = Treatment), width = 0.1,
                position = position_dodge(width = 0.8))+
  theme_bw()+
  facet_grid(.~Date)+
  scale_fill_manual(breaks = c("Control","Feeding"), values=c(myColours[c(1:2)], "blue","red"))+
  labs(y = expression('Taxonomic diversity - D'[2]), x = "Timepoint", fill = "Treatment")+
  theme(axis.title=element_text(size=16), strip.text.x=element_text(size=16),
        legend.title=element_text(size=15),legend.text=element_text(size=14),
        axis.text = element_text(size=14),title=element_text(size=20),
        strip.background=element_rect(fill=adjustcolor("lightgray",0.2))
        #,panel.grid.major = element_blank(), panel.grid.minor = element_blank()
        )+ 
  geom_text(data = tp, aes(x=1.5,y=11,label=Measurement), size = 8, inherit.aes = FALSE,
            hjust = 0.5)

png(file="Fig5_PROPS.png",width=12,height=6,res=500,units="in", pointsize=12)
pdf(file="Fig5_PROPS.pdf",width=12,height=6)

# png(file="Fig5_PROPS.png",width=12,height=6,res=500,units="in", pointsize=12)
# pdf(file="Fig5_PROPS.pdf",width=12,height=6)
p_msphere_predict
```

```
## Warning: Removed 42 rows containing missing values (geom_errorbar).
```

```r
# dev.off()
```


## Figure S3: Gating strategy  


```r
path = "data_mussel"
flowData <- read.flowSet(path = path, transformation = FALSE, pattern=".fcs")

flowData_transformed <- transform(flowData,`FL1-H`=asinh(`FL1-H`), 
                                  `SSC-H`=asinh(`SSC-H`), 
                                  `FL3-H`=asinh(`FL3-H`), 
                                  `FSC-H`=asinh(`FSC-H`))
param=c("FL1-H", "FL3-H","SSC-H","FSC-H")
flowData_transformed = flowData_transformed[,param]
remove(flowData)

# Create a PolygonGate for denoising the dataset
# Define coordinates for gate in sqrcut1 in format: c(x,x,x,x,y,y,y,y)
sqrcut1 <- matrix(c(8.5,8.5,15,15,3,8,14,3),ncol=2, nrow=4)
colnames(sqrcut1) <- c("FL1-H","FL3-H")
polyGate1 <- polygonGate(.gate=sqrcut1, filterId = "Total Cells")

summary <- fsApply(x=flowData_transformed,FUN=function(x) apply(x,2,max),use.exprs=TRUE)
max = max(summary[,1])
max <- 14.99341
mytrans <- function(x) x/max
flowData_transformed <- transform(flowData_transformed,`FL1-H`=mytrans(`FL1-H`),
                                  `FL3-H`=mytrans(`FL3-H`), 
                                  `SSC-H`=mytrans(`SSC-H`),
                                  `FSC-H`=mytrans(`FSC-H`))

flowData_transformed <- flowData_transformed[which(flowCore::sampleNames(flowData_transformed) 
                                                   %in% c("Q3_T_t0_rep1.fcs","Q3_T_t1.5_rep1.fcs","Q3_T_t3_rep1.fcs",
                                                         "C1_t0_rep1.fcs","C1_t1.5_rep1.fcs","C1_t3_rep1.fcs"))]

MyText <- c("Control - 0h", "Control - 1.5h", "Control - 3h",
            "Feeding - 0h", "Feeding - 1.5h", "Feeding - 3h")

sqrcut1 <- matrix(c(asinh(12500),asinh(12500),15,15,3,9.55,14,3)/max,ncol=2, nrow=4)
colnames(sqrcut1) <- c("FL1-H","FL3-H")
rGate_HNA <- polygonGate(.gate=sqrcut1, filterId = "HNA")
sqrcut1 <- matrix(c(8.5,8.5,asinh(12500),asinh(12500),3,8,9.55,3)/max,ncol=2, nrow=4)
colnames(sqrcut1) <- c("FL1-H","FL3-H")
rGate_LNA <- polygonGate(.gate=sqrcut1, filterId = "LNA")

filters <- filters(list(rGate_LNA, rGate_HNA))
flist <- list(filters , filters, filters, filters, filters, filters)
names(flist) <- flowCore::sampleNames(flowData_transformed)

print(xyplot(`FL3-H`~`FL1-H`, data=flowData_transformed,index.cond=list(c(1:6)),
             filter=flist,
             xbins=200,nbin=128, par.strip.text=list(col="black", font=3,cex=1.85), 
             smooth=FALSE, xlim=c(0.5,1),ylim=c(0.1,1),xlab=list(label="Green fluorescence intensity (FL1-H)",cex=2),ylab=list(label="Red fluorescence intensity (FL3-H)",cex=2),
             par.settings=my.settings,
             scales=list(x=list(at=seq(from=0, to=1, by=.1),cex=1),
                         y=list(at=seq(from=0, to=1, by=.2),cex=1)), layout=c(3,2),
             strip=strip.custom(factor.levels=MyText),
             margin=TRUE,
             binTrans="log"
      )
)
```

<img src="Figures/cached/Gating-1.png" style="display: block; margin: auto;" />

## Figure 5: Identify physiological populations  


```r
comp_t0 <- fp_contrasts(x=fbasis, comp2=(pos=="C1_t0.fcs" | pos== "C2_t0.fcs"),
                        comp1=(pos=="Q1_T_t0.fcs" | pos== "Q2_T_t0.fcs" | pos== "Q3_T_t0.fcs"), 
                        thresh=0.04)
```

```
## 	Region used for contrasts 1 16384
## 	Returning contrasts for Q1_T_t0_rep1.fcs C1_t0_rep1.fcs
##  	Returning contrasts for Q1_T_t0_rep2.fcs C1_t0_rep2.fcs
##  	Returning contrasts for Q1_T_t0_rep3.fcs C1_t0_rep3.fcs
##  	Returning contrasts for Q2_T_t0_rep1.fcs C2_t0_rep1.fcs
##  	Returning contrasts for Q2_T_t0_rep2.fcs C2_t0_rep2.fcs
##  	Returning contrasts for Q2_T_t0_rep3.fcs C2_t0_rep3.fcs
##  	Returning contrasts for Q3_T_t0_rep1.fcs C1_t0_rep1.fcs
##  	Returning contrasts for Q3_T_t0_rep2.fcs C1_t0_rep2.fcs
##  	Returning contrasts for Q3_T_t0_rep3.fcs C1_t0_rep3.fcs
```

```r
comp_t1.5 <- fp_contrasts(x=fbasis, comp2=(pos=="C1_t1.5.fcs" | pos== "C2_t1.5.fcs"),
                          comp1=(pos=="Q1_T_t1.5.fcs" | pos== "Q2_T_t1.5.fcs" | pos== "Q3_T_t1.5.fcs"), 
                          thresh=0.04)
```

```
## 	Region used for contrasts 1 16384
## 	Returning contrasts for Q1_T_t1.5_rep1.fcs C1_t1.5_rep1.fcs
##  	Returning contrasts for Q1_T_t1.5_rep2.fcs C1_t1.5_rep2.fcs
##  	Returning contrasts for Q1_T_t1.5_rep3.fcs C1_t1.5_rep3.fcs
##  	Returning contrasts for Q2_T_t1.5_rep1.fcs C2_t1.5_rep1.fcs
##  	Returning contrasts for Q2_T_t1.5_rep2.fcs C2_t1.5_rep2.fcs
##  	Returning contrasts for Q2_T_t1.5_rep3.fcs C2_t1.5_rep3.fcs
##  	Returning contrasts for Q3_T_t1.5_rep1.fcs C1_t1.5_rep1.fcs
##  	Returning contrasts for Q3_T_t1.5_rep2.fcs C1_t1.5_rep2.fcs
##  	Returning contrasts for Q3_T_t1.5_rep3.fcs C1_t1.5_rep3.fcs
```

```r
comp_t3 <- fp_contrasts(x=fbasis, comp2=(pos=="C1_t3.fcs" | pos== "C2_t3.fcs"),
                        comp1=(pos=="Q1_T_t3.fcs" | pos== "Q2_T_t3.fcs" | pos== "Q3_T_t3.fcs"), 
                        thresh=0.04)
```

```
## 	Region used for contrasts 1 16384
## 	Returning contrasts for Q1_T_t3_rep1.fcs C1_t3_rep1.fcs
##  	Returning contrasts for Q1_T_t3_rep2.fcs C1_t3_rep2.fcs
##  	Returning contrasts for Q2_T_t3_rep1.fcs C1_t3_rep3.fcs
##  	Returning contrasts for Q2_T_t3_rep2.fcs C2_t3_rep1.fcs
##  	Returning contrasts for Q2_T_t3_rep3.fcs C2_t3_rep2.fcs
##  	Returning contrasts for Q3_T_t3_rep1.fcs C2_t3_rep3.fcs
##  	Returning contrasts for Q3_T_t3_rep2.fcs C1_t3_rep1.fcs
##  	Returning contrasts for Q3_T_t3_rep3.fcs C1_t3_rep2.fcs
```

```r
# Merge data frames
comp_total <- rbind(comp_t0, comp_t1.5, comp_t3)
comp_total <- data.frame(comp_total, Timepoint=as.factor(c(rep("0h Feeding", nrow(comp_t0)), 
                                                 rep("1.5h Feeding", nrow(comp_t1.5)),
                                                 rep("3h Feeding", nrow(comp_t3)))))

p12b <- ggplot(data=results, aes(x=Time, y=HNA.cells, fill=Treatment)) + 
  # geom_boxplot(alpha=0.9)+
  geom_point(shape=21, size=5,alpha=0.9)+
  scale_fill_manual(values=myColours[c(1,2)])+
  theme_bw()+
  labs(y=expression("HNA cells µL"^"-1"), x="Time (h)", title="C", fill="")+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=16),
        title=element_text(size=20), legend.text=element_text(size=14),
        legend.title=element_text(size=15),
        # panel.grid.major = element_blank(), panel.grid.minor = element_blank()
        legend.direction = "horizontal",legend.position = "bottom",
        panel.grid.minor = element_blank()
  )+ 
  # guides(fill=FALSE)+
  geom_errorbar(aes(ymin=HNA.cells-sd.HNA.cells, ymax=HNA.cells+sd.HNA.cells), width=0.075)+
  ylim(200,525)+ 
  geom_smooth(method="rlm",color="black", alpha=0.2)+
  scale_x_continuous(breaks=c(0,0.5,1,1.5,2,2.5,3))

p13 <- ggplot(data=results, aes(x=Time, y=LNA.cells, fill=Treatment)) + 
  # geom_boxplot(alpha=0.9)+
  geom_point(shape=21, size=5,alpha=0.9)+
  scale_fill_manual(values=myColours[c(1,2)])+
  # geom_smooth(formula=y ~ x, color="black")+
  # geom_boxplot(mapping=factor(Time),alpha=0.4,outlier.shape=NA)+
  theme_bw()+
  labs(y=expression("LNA cells µL"^"-1"), x="Time (h)", title="B", fill="")+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=16),
        title=element_text(size=20), legend.text=element_text(size=14),
        legend.title=element_text(size=15),
        legend.direction = "horizontal",legend.position = "bottom",
        panel.grid.minor = element_blank())+ 
  # guides(fill=FALSE)+
  geom_errorbar(aes(ymin=LNA.cells-sd.LNA.cells, ymax=LNA.cells+sd.LNA.cells), width=0.075)+
  ylim(1000,1325)+
  scale_x_continuous(breaks=c(0,0.5,1,1.5,2,2.5,3))

# Contrast plot
vtot <- ggplot(comp_total, aes(`FL1.H`, `FL3.H`, z = Density))+
  geom_tile(aes(fill=Density)) + 
  geom_point(colour="gray", alpha=0.4)+
  scale_fill_distiller(palette="RdBu", na.value="white") + 
  stat_contour(aes(fill=..level..), geom="polygon", binwidth=0.1)+
  theme_bw()+
  geom_contour(color = "white", alpha = 1)+
  facet_grid(~Timepoint)+
  labs(x="Green fluorescence intensity (a.u.)", y="Red fluorescence intensity (a.u.)",title="A")+
  theme(axis.title=element_text(size=16), strip.text.x=element_text(size=16),
        legend.title=element_text(size=15),legend.text=element_text(size=14),
        axis.text = element_text(size=14),title=element_text(size=20),
        strip.background=element_rect(fill=adjustcolor("lightgray",0.2))
        #,panel.grid.major = element_blank(), panel.grid.minor = element_blank()
        )

# Make the plot
g1 <- ggplotGrob(p13);
g2 <- ggplotGrob(p12b);
gv <- ggplotGrob(vtot);
fg1 <- gtable_frame(g1)
fg2 <- gtable_frame(g2)
fg12 <- gtable_cbind(fg1,fg2)
fg3 <- gtable_frame(gv)
grid.newpage()
# combined <- rbind(fg3, fg12)
# png(file = "contrasts_0.02.png", width = 12, height = 12, res = 500, 
    # units = "in", pointsize = 10)
grid.draw(fg3)
```

<img src="Figures/cached/Contrasts-1.png" style="display: block; margin: auto;" />

```r
grid.newpage()
grid.draw(fg12)
```

<img src="Figures/cached/Contrasts-2.png" style="display: block; margin: auto;" />

```r
# dev.off()
```

## Calculate removal statistics  


```r
# Calculate coefficient of variation within LNA and DNA of control/treatment
# CV for control in LNA population
100*sd(results$LNA[results$Treatment=="Control"])/mean(results$LNA[results$Treatment=="Control"])
```

```
## [1] 5.820442
```

```r
# CV for treatment in LNA population
100*sd(results$LNA[results$Treatment=="Feeding"])/mean(results$LNA[results$Treatment=="Feeding"])
```

```
## [1] 5.231828
```

```r
# CV for control in HNA population
100*sd(results$HNA[results$Treatment=="Control"])/mean(results$HNA[results$Treatment=="Control"])
```

```
## [1] 5.078665
```

```r
# CV for treatment in HNA population
100*sd(results$HNA[results$Treatment=="Feeding"])/mean(results$HNA[results$Treatment=="Feeding"])
```

```
## [1] 12.79566
```

```r
# Calculate the removal rate of HNA bacteria
# We use robust regression because of the two (suspected) outliers at t0.
lm.HNA <- rlm(HNA.cells~Time, data=results[results$Treatment=="Feeding",])
lm.HNA_C <- rlm(HNA.cells~Time, data=results[results$Treatment=="Control",])
car::Anova(lm.HNA_C) # p = 0.9826
```

```
## Analysis of Deviance Table (Type II tests)
## 
## Response: HNA.cells
##           Df     F Pr(>F)
## Time       1 5e-04 0.9826
## Residuals 12
```

```r
car::Anova(lm.HNA) # p = 9.02e-11
```

```
## Analysis of Deviance Table (Type II tests)
## 
## Response: HNA.cells
##           Df      F   Pr(>F)    
## Time       1 163.03 9.02e-11 ***
## Residuals 19                    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
# To get std. error on parameter estimation
summary(lm.HNA)
```

```
## 
## Call: rlm(formula = HNA.cells ~ Time, data = results[results$Treatment == 
##     "Feeding", ])
## Residuals:
##      Min       1Q   Median       3Q      Max 
## -26.7479  -7.8385   0.8337   7.6519  26.0201 
## 
## Coefficients:
##             Value    Std. Error t value 
## (Intercept) 429.7656   6.0957    70.5027
## Time        -43.1736   3.3813   -12.7683
## 
## Residual standard error: 11.62 on 19 degrees of freedom
```

```r
# Average HNA removal relative to HNA pool
t0 <- mean(results$HNA.cells[results$Treatment=="Feeding"][results$Time[results$Treatment=="Feeding"]=="0"])
t3 <- mean(results$HNA.cells[results$Treatment=="Feeding"][results$Time[results$Treatment=="Feeding"]=="3"])
100*(t0-t3)/t0
```

```
## [1] 28.98963
```

```r
# Error on average HNA densities
e0 <- sqrt(sum(results$sd.HNA.cells[results$Treatment=="Feeding"][results$Time[results$Treatment=="Feeding"]=="0"]^2))/3
e3 <- sqrt(sum(results$sd.HNA.cells[results$Treatment=="Feeding"][results$Time[results$Treatment=="Feeding"]=="3"]^2))/3

# Error on HNA removal ratio
100*sqrt((sqrt(e0^2 + e3^2)/(t0-t3))^2 + (e0/t0)^2)
```

```
## [1] 4.952725
```

```r
# Average HNA removal relative to total cell pool
t0 <- mean(results$HNA.cells[results$Treatment=="Feeding"][results$Time[results$Treatment=="Feeding"]=="0"])
t3 <- mean(results$HNA.cells[results$Treatment=="Feeding"][results$Time[results$Treatment=="Feeding"]=="3"])
t.C <- mean(results$Total.cells[results$Treatment=="Feeding"][results$Time[results$Treatment=="Feeding"]=="3"])
100*(t0-t3)/t.C
```

```
## [1] 9.165898
```

```r
# Error on average HNA densities
e0 <- sqrt(sum(results$sd.HNA.cells[results$Treatment=="Feeding"][results$Time[results$Treatment=="Feeding"]=="0"]^2))/3
e3 <- sqrt(sum(results$sd.HNA.cells[results$Treatment=="Feeding"][results$Time[results$Treatment=="Feeding"]=="3"]^2))/3
e.C <- sqrt(sum(results$sd.Total.cells[results$Treatment=="Feeding"][results$Time[results$Treatment=="Feeding"]=="3"]^2))/3

# Error on HNA removal ratio
100*sqrt((sqrt(e0^2 + e3^2)/(t0-t3))^2 + (e.C/t.C)^2)
```

```
## [1] 4.900268
```

```r
# Rate per mg DW of IDM and associated error
1000*lm.HNA$coefficients[2]/mean(1000*meta.mus$weight_g)
```

```
##      Time 
## -178.1615
```

```r
sqrt((sd(meta.mus$weight_g)/mean(meta.mus$weight_g))^2 + (summary(lm.HNA)$coefficients[2,2]/lm.HNA$coefficients[2])^2)
```

```
##      Time 
## 0.1072376
```

```r
# Calculage clearance rate (linear regression used to get dC/dt)
# Volume = 11L for each bucket
11*abs(lm.HNA$coefficients[2])/lm.HNA$coefficients[1]/mean(meta.mus$weight_g)
```

```
##     Time 
## 4.560106
```

```r
# error on first quotient
eq1 <- sqrt((summary(lm.HNA)$coefficients[1,2]/lm.HNA$coefficients[1])^2 + (summary(lm.HNA)$coefficients[2,2]/lm.HNA$coefficients[2])^2)
# error on final clearance rate
11*sqrt(((eq1/abs(lm.HNA$coefficients[2])/lm.HNA$coefficients[1])^2) + (sd(meta.mus$weight_g)/mean(meta.mus$weight_g))^2)
```

```
## (Intercept) 
##   0.8057903
```
