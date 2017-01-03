library("Phenoflow")
library("dplyr")
library("gridExtra")
library("ggplot2")
library("formatR")
library("gridExtra")
library("cowplot")
library("RColorBrewer")
library("vegan")
source("functions.R")
myColours <- brewer.pal("Accent",n=3)
# For file renaming
# path <- "data/Mosselstaaltjes_deel2"
# setwd(path)
# baseFolder = path
# a = list.files(pattern="")
# file.rename(a, paste0(gsub(",",".", a)))

### Samples were diluted 2x
dilution <- 2

path = "data"
flowData <- read.flowSet(path = path, transformation = FALSE, pattern=".fcs")

flowData_transformed <- transform(flowData,`FL1-H`=asinh(`FL1-H`), 
                                  `SSC-H`=asinh(`SSC-H`), 
                                  `FL3-H`=asinh(`FL3-H`), 
                                  `FSC-H`=asinh(`FSC-H`))
param=c("FL1-H", "FL3-H","SSC-H","FSC-H")
flowData_transformed = flowData_transformed[,param]
remove(flowData)

### Create a PolygonGate for denoising the dataset
### Define coordinates for gate in sqrcut1 in format: c(x,x,x,x,y,y,y,y)
sqrcut1 <- matrix(c(8.5,8.5,15,15,3,8,14,3),ncol=2, nrow=4)
colnames(sqrcut1) <- c("FL1-H","FL3-H")
polyGate1 <- polygonGate(.gate=sqrcut1, filterId = "Total Cells")

###  Gating quality check
xyplot(`FL3-H` ~ `FL1-H`, data=flowData_transformed[1], filter=polyGate1,
       scales=list(y=list(limits=c(0,14)),
                   x=list(limits=c(6,16))),
       axis = axis.default, nbin=125, 
       par.strip.text=list(col="white", font=2, cex=2), smooth=FALSE)

### Isolate only the cellular information based on the polyGate1
flowData_transformed <- Subset(flowData_transformed, polyGate1)

summary <- fsApply(x=flowData_transformed,FUN=function(x) apply(x,2,max),use.exprs=TRUE)
max = max(summary[,1])
mytrans <- function(x) x/max
flowData_transformed <- transform(flowData_transformed,`FL1-H`=mytrans(`FL1-H`),
                                  `FL3-H`=mytrans(`FL3-H`), 
                                  `SSC-H`=mytrans(`SSC-H`),
                                  `FSC-H`=mytrans(`FSC-H`))


### Calculate fingerprint with bw = 0.01
fbasis <- flowBasis(flowData_transformed, param, nbin=128, 
                    bw=0.01,normalize=function(x) x)

### Calculate ecological parameters from normalized fingerprint 
### Densities will be normalized to the interval [0,1]
### d = rounding factor
# Diversity.fbasis <- Diversity(fbasis, d = 3, plot = TRUE, R = 999)
Diversity.fbasis <- Diversity_rf(flowData_transformed, d=3, param = param, R = 3)


### make metadata table
tmp <- strsplit(rownames(Diversity.fbasis)[1:42], "_")
meta.div <- cbind(do.call(rbind, lapply(tmp, rbind))[,1],rep("C",42), do.call(rbind, lapply(tmp, rbind))[,2:3])
tmp <- strsplit(rownames(Diversity.fbasis)[43:nrow(Diversity.fbasis)], "_")
meta.div <- data.frame(rbind(meta.div, do.call(rbind, lapply(tmp, rbind))))
colnames(meta.div) <- c("Sample","Treatment","Time","Replicate")
meta.div$Replicate <- gsub(meta.div$Replicate,pattern=".fcs",replacement="")
meta.div$Time <- as.numeric(gsub(meta.div$Time,pattern="t",replacement=""))


### Merge with Diversity.fbasis
Diversity.fbasis <- cbind(Diversity.fbasis, meta.div)

### Remove outliers due to human errors (forgot staining)
fbasis@basis <- fbasis@basis[Diversity.fbasis$D2>1600,]
flowData_transformed <- flowData_transformed[Diversity.fbasis$D2>1600]

sqrcut1 <- matrix(c(asinh(12500),asinh(12500),15,15,3,9.55,14,3)/max,ncol=2, nrow=4)
colnames(sqrcut1) <- c("FL1-H","FL3-H")
rGate_HNA <- polygonGate(.gate=sqrcut1, filterId = "HNA")
sqrcut1 <- matrix(c(8.5,8.5,asinh(12500),asinh(12500),3,8,9.55,3)/max,ncol=2, nrow=4)
colnames(sqrcut1) <- c("FL1-H","FL3-H")
rGate_LNA <- polygonGate(.gate=sqrcut1, filterId = "LNA")

### Diversities of HNA/LNA populations
flowData_HNA <- split(flowData_transformed, rGate_HNA)$`HNA+`
flowData_LNA <- split(flowData_transformed, rGate_LNA)$`LNA+`

Diversity.HNA <- Diversity_rf(flowData_HNA, d=3, param = param, R = 3)
Diversity.LNA <- Diversity_rf(flowData_LNA, d=3, param = param, R = 3)

### Counts
### Normalize total cell gate
sqrcut1 <- matrix(c(8.5,8.5,15,15,3,8,14,3)/max,ncol=2, nrow=4)
colnames(sqrcut1) <- c("FL1-H","FL3-H")
polyGate1 <- polygonGate(.gate=sqrcut1, filterId = "Total Cells")

### Check if rectangle gate is correct, if not, adjust rGate_HNA
xyplot(`FL3-H` ~ `FL1-H`, data=flowData_transformed[1], filter=rGate_HNA,
       scales=list(y=list(limits=c(0,1)),
                   x=list(limits=c(0.4,1))),
       axis = axis.default, nbin=125, par.strip.text=list(col="white", font=2, 
                                                          cex=2), smooth=FALSE)
### Extract the cell counts
a <- flowCore::filter(flowData_transformed, rGate_HNA) 
HNACount <- flowCore::summary(a);HNACount <- toTable(HNACount)
s <- flowCore::filter(flowData_transformed, polyGate1)
TotalCount <- flowCore::summary(s);TotalCount <- flowCore::toTable(TotalCount)

### Extract the volume
vol <- c()
for(i in 1:length(flowData_transformed)){
  vol[i] <- as.numeric(flowData_transformed[[i]]@description$`$VOL`)/1000
}

### Save counts
counts <- data.frame(Samples=flowCore::sampleNames(flowData_transformed), 
                     Total.cells = dilution*TotalCount$true/vol, HNA.cells = dilution*HNACount$true/vol,
                     LNA.cells = dilution*(TotalCount$true-HNACount$true)/vol)

### Pool results
Diversity.fbasis <- Diversity.fbasis[Diversity.fbasis$D2>1600,]
Storage <- c(); Storage[Diversity.fbasis$Treatment=="T"] <- "Covered"; Storage[Diversity.fbasis$Treatment=="S"] <- "Submerged"
results <- cbind(Diversity.fbasis,counts, Storage, Diversity.HNA, Diversity.LNA)

### Add an extra column for biological replicates
bio_rep <- c(rep(1,21),rep(2,21),rep(1,41),rep(2,41),rep(3,41))
results <- cbind(results,bio_rep=factor(bio_rep))

### Select only S experiment for beta-div and further analysis

levels(results$Treatment)[levels(results$Treatment)=="C"] <- "Control"
levels(results$Treatment)[levels(results$Treatment)=="T"] <- "Feeding - T"
levels(results$Treatment)[levels(results$Treatment)=="S"] <- "Feeding - S"

### Beta diversity analysis
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
fbasis1@basis <- round(fbasis1@basis, 4)
x <- by(fbasis1@basis, INDICES = pos[!Diversity.fbasis$Treatment=="T"], FUN = colMeans)
x <- do.call(rbind, x)
# Calculate distance matrix
dist.S <- vegdist(x, method="bray")

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


##############################################################################
### Make some plots
##############################################################################
result.tmp <- results
results <- result.tmp[!result.tmp$Treatment=="Feeding - T", ]
results$Treatment <- droplevels(plyr::revalue(results$Treatment, c("Feeding - S"="Feeding")))
p1 <- ggplot(data=results, aes(x=factor(Time), y=D2, fill=Treatment)) + 
  # geom_boxplot(alpha=0.9)+
  geom_point(shape=21, size=7,alpha=0.9)+
  scale_fill_manual(values=myColours[c(1,2)])+
  # geom_smooth(formula=y ~ x, color="black")+
  # geom_boxplot(mapping=factor(Time),alpha=0.4,outlier.shape=NA)+
  theme_bw()+
  labs(y=expression('Phenotypic diversity - D'[2]), x="Time (h)", title="A",
       fill="")+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=16),
        # panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.direction = "horizontal",legend.position = "bottom")+ 
  geom_errorbar(aes(ymin=D2-sd.D2, ymax=D2+sd.D2), width=0.075)


# p2 <- ggplot(data=results, aes(x=factor(Time), y=Total.cells, fill=Treatment)) + 
#   # geom_boxplot(alpha=0.9)+
#   geom_point(shape=21, size=5,alpha=0.9)+
#   geom_line()+
#   scale_fill_manual(values=myColours[c(1,2)])+
#   # geom_smooth(formula=y ~ x, color="black")+
#   theme_bw()+
#   labs(y="Cells/µL", x="Time (h)", title="A. Cell density")+
#   theme(axis.text=element_text(size=14), axis.title=element_text(size=16,face="bold"),
#         title=element_text(size=20), legend.text=element_text(size=14),
#         panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ 
#   guides(fill=FALSE)+
#   geom_errorbar(aes(ymin=Total.cells-sd.Total.cells, ymax=Total.cells+sd.Total.cells), width=0.075)

### Beta diversity
beta.div.data <- data.frame(beta.div$points, tmp)
beta.div.data <- droplevels(beta.div.data)
var <- round(vegan::eigenvals(beta.div)/sum(vegan::eigenvals(beta.div))*100,1)

p.beta <- ggplot(data=beta.div.data, aes(x=X1, y=X2, fill=Treatment, size=Time))+
  geom_point(shape=21, alpha=0.9, size=7)+
  # scale_size(range=c(4,10), breaks=c(0,0.5,1,1.5,2,2.5,3))+ 
  # guides(fill = guide_legend(override.aes = list(size=5)))+
  theme_bw()+
  scale_fill_brewer(palette="Accent")+
  labs(x = paste0("Axis1 (",var[1], "%)"), y = paste0("Axis2 (",var[2], "%)"), title="C. Phenotypic beta diversity")+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=16,face="bold"),
        title=element_text(size=20), legend.text=element_text(size=14))
print(p.beta)
  
png(file="Submerged.fcm_pooledD.png",width=12,height=6,res=500,units="in", pointsize=12)
# grid.arrange(arrangeGrob(p2,p1, ncol=2), p.beta.S, heights=c(4/4, 4/4), ncol=1)
# grid.arrange(p1,p.beta.S, ncol=2)
grid_arrange_shared_legend(p1b,p.beta.S, ncol=2)
dev.off()

### Separate S
beta.div.data.S <- data.frame(beta.div.S$points, tmp[!tmp$Treatment=="Feeding - T",])
beta.div.data.S <- droplevels(beta.div.data.S)
beta.div.data.S$Treatment <- plyr::revalue(beta.div.data.S$Treatment, c("Feeding - S"="Feeding"))
var <- round(vegan::eigenvals(beta.div.S)/sum(vegan::eigenvals(beta.div.S))*100,1)
p.beta.S <- ggplot(data=beta.div.data.S, aes(x=X1, y=X2, fill=Treatment, size=Time))+
  geom_point(shape=21, alpha=1)+
  scale_size(range=c(4,10), breaks=c(0,0.5,1,1.5,2,2.5,3))+ 
  guides(fill = guide_legend(override.aes = list(size=5)), size=FALSE)+
  theme_bw()+
  scale_fill_manual(values=myColours[c(1,2)])+
  labs(x = paste0("PCoA axis 1 (",var[1], "%)"), y = paste0("PCoA axis 2 (",var[2], "%)"), title="B",
       fill="")+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=16),
        # panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.direction = "horizontal",legend.position = "bottom")
print(p.beta.S)

### Separate T
beta.div.data.T <- data.frame(beta.div.T$points, tmp[!tmp$Treatment=="Feeding - S",])
beta.div.data.T <- droplevels(beta.div.data.T)
var <- round(vegan::eigenvals(beta.div.T)/sum(vegan::eigenvals(beta.div.T))*100,1)
p.beta.T <- ggplot(data=beta.div.data.T, aes(x=X1, y=X2, fill=Treatment, size=Time))+
  geom_point(shape=21, alpha=1)+
  scale_size(range=c(4,10), breaks=c(0,0.5,1,1.5,2,2.5,3))+ 
  guides(fill = guide_legend(override.aes = list(size=5)))+
  theme_bw()+
  scale_fill_manual(values=myColours[c(1,3)])+
  labs(x = paste0("Axis1 (",var[1], "%)"), y = paste0("Axis2 (",var[2], "%)"), title="C. Phenotypic beta diversity")+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=16,face="bold"),
        title=element_text(size=20), legend.text=element_text(size=14))
print(p.beta.T)

##############################################################################
### Additional statistics
##############################################################################

### Calculate coefficient of variation within LNA and DNA of control/treatment
# CV for control in LNA population
100*sd(results$LNA[results$Treatment=="Control"])/mean(results$LNA[results$Treatment=="Control"])
# CV for treatment in LNA population
100*sd(results$LNA[results$Treatment=="Feeding"])/mean(results$LNA[results$Treatment=="Feeding"])
# CV for control in HNA population
100*sd(results$HNA[results$Treatment=="Control"])/mean(results$HNA[results$Treatment=="Control"])
# CV for treatment in HNA population
100*sd(results$HNA[results$Treatment=="Feeding"])/mean(results$HNA[results$Treatment=="Feeding"])

### Calculate the removal rate of HNA bacteria
### We use robust regression because of the two (suspected) outliers at t0.
lm.HNA <- rlm(HNA.cells~Time, data=results[results$Treatment=="Feeding",])
lm.HNA_C <- rlm(HNA.cells~Time, data=results[results$Treatment=="Control",])
car::Anova(lm.HNA_C) # p = 0.9826
car::Anova(lm.HNA) # p = 9.02e-11

### Plot these regressions in a new figure
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
        legend.direction = "horizontal",legend.position = "bottom"
  )+ 
  # guides(fill=FALSE)+
  geom_errorbar(aes(ymin=HNA.cells-sd.HNA.cells, ymax=HNA.cells+sd.HNA.cells), width=0.075)+
  ylim(200,525)+ 
  geom_smooth(method="rlm",color="black", alpha=0.2)

### Calculate the dynamics of D2 in time for treatment/control
### We can't use a linear regression here (for obvious reasons)
### Hence we will try to fit a spline and see if we can make inference
### this way
sp_T <- lm(D2~ns(Time, df=3), data=results[results$Treatment=="Feeding",])
sp_C <- lm(D2~ns(Time, df=3), data=results[results$Treatment=="Control",])

### Order studentized residuals according to timepoint
res.T <- studres(sp_T)[order(results[results$Treatment=="Feeding",]$Time)]
res.C <- studres(sp_C)[order(results[results$Treatment=="Control",]$Time)]

## Location of knots
attr(ns(results$Time, df=3), "knots")

### Check for temporal autocorrelation in model residuals
png(file="acf_D2.png",width=10,height=5,res=500,units="in", pointsize=12)
par(mfrow=c(1,2))
acf(res.C, main="Treatment: control", las=1)
acf(res.T, main="Treatment: feeding", las=1)
dev.off()

### Perform statistical inference on splines
waldtest(sp_C) # p = 0.168
waldtest(sp_T, vcov = vcovHAC(sp_T)) # p = 1.292e-05

p1b <- ggplot(data=results, aes(x=Time, y=D2, fill=Treatment)) + 
  # geom_boxplot(alpha=0.9)+
  geom_point(shape=21, size=7,alpha=0.9)+
  scale_fill_manual(values=myColours[c(1,2)])+
  # geom_smooth(formula=y ~ x, color="black")+
  # geom_boxplot(mapping=factor(Time),alpha=0.4,outlier.shape=NA)+
  theme_bw()+
  labs(y=expression('Phenotypic diversity - D'[2]), x="Time (h)", title="A",
       fill="")+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=16),
        # panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.direction = "horizontal",legend.position = "bottom")+ 
  geom_errorbar(aes(ymin=D2-sd.D2, ymax=D2+sd.D2), width=0.075)+
  geom_smooth(method="lm", color="black", alpha=0.2, formula = y ~ splines::ns(x,df=3))


### PERMANOVA on beta diversity analysis
disper.test <- betadisper(dist.S, group=results$Treatment)
disper.test # average distance to mean 0.03 for both groups
anova(disper.test) # P = 0.891
adonis(dist.S~Time*Treatment, data=results) 

### No effect on D2 within populations
# p14 <- ggplot(data=results, aes(x=factor(Time), y=D2.HNA, fill=Treatment)) + 
#   # geom_boxplot(alpha=0.9)+
#   geom_point(shape=21, size=5,alpha=0.9)+
#   scale_fill_manual(values=myColours[c(1,2)])+
#   # geom_smooth(formula=y ~ x, color="black")+
#   # geom_boxplot(mapping=factor(Time),alpha=0.4,outlier.shape=NA)+
#   theme_bw()+
#   labs(y="Phenotypic diversity D2", x="Time (h)", title="HNA population")+
#   theme(axis.text=element_text(size=14), axis.title=element_text(size=16,face="bold"),
#         title=element_text(size=20), legend.text=element_text(size=14))+ 
#   guides(fill=FALSE)+
#   geom_errorbar(aes(ymin=D2.HNA -sd.D2.HNA , ymax=D2.HNA +sd.D2.HNA ), width=0.075)+
#   ylim(1000,1350)
# 
# p15 <- ggplot(data=results, aes(x=factor(Time), y=LNA.cells, fill=Treatment)) + 
#   # geom_boxplot(alpha=0.9)+
#   geom_point(shape=21, size=5,alpha=0.9)+
#   scale_fill_manual(values=myColours[c(1,2)])+
#   # geom_smooth(formula=y ~ x, color="black")+
#   # geom_boxplot(mapping=factor(Time),alpha=0.4,outlier.shape=NA)+
#   theme_bw()+
#   labs(y="Phenotypic diversity D2", x="Time (h)", title="LNA population")+
#   theme(axis.text=element_text(size=14), axis.title=element_text(size=16,face="bold"),
#         title=element_text(size=20), legend.text=element_text(size=14))+ 
#   guides(fill=FALSE)+
#   geom_errorbar(aes(ymin=LNA.cells-sd.LNA.cells, ymax=LNA.cells+sd.LNA.cells), width=0.075)+
#   ylim(1000,1350)
# 
# png(file="HNA.LNA.png",width=12,height=5,res=500,units="in", pointsize=12)
# grid.arrange(p14, p15, ncol=2)
# dev.off()
# 
# 
