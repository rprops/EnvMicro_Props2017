library("Phenoflow")
library("dplyr")
library("gridExtra")
library("ggplot2")
library("formatR")
library("gridExtra")
# For file renaming
# path <- "data"
# setwd(path)
# baseFolder = path
# a = list.files(pattern=".")
# file.rename(a, paste0(gsub("l","", a)))


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
Diversity.fbasis <- Diversity(fbasis,d=3,plot=TRUE, R=999)

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

### Counts
### Creating a rectangle gate for counting HNA and LNA cells
rGate_HNA <- rectangleGate("FL1-H"=c(asinh(13500), 20)/max,"FL3-H"=c(0,20)/max, 
                           filterId = "HNA bacteria")
### Normalize total cell gate
sqrcut1 <- matrix(c(8.5,8.5,14,14,3,7.5,14,3)/max,ncol=2, nrow=4)
colnames(sqrcut1) <- c("FL1-H","FL3-H")
polyGate1 <- polygonGate(.gate=sqrcut1, filterId = "Total Cells")

### Check if rectangle gate is correct, if not, adjust rGate_HNA
xyplot(`FL3-H` ~ `FL1-H`, data=flowData_transformed[1], filter=rGate_HNA,
       scales=list(y=list(limits=c(0,1)),
                   x=list(limits=c(0.4,1))),
       axis = axis.default, nbin=125, par.strip.text=list(col="white", font=2, 
                                                          cex=2), smooth=FALSE)
### Extract the cell counts
a <- filter(flowData_transformed, rGate_HNA) 
HNACount <- summary(a);HNACount <- toTable(HNACount)
s <- filter(flowData_transformed, polyGate1)
TotalCount <- summary(s);TotalCount <- toTable(TotalCount)

### Extract the volume
vol <- c()
for(i in 1:length(flowData_transformed)){
  vol[i] <- as.numeric(flowData_transformed[[i]]@description$`$VOL`)/1000
}

### Save counts
dilution <- 2
counts <- data.frame(Samples=flowCore::sampleNames(flowData_transformed), 
                     Total.cells = dilution*TotalCount$true/vol, HNA.cells = dilution*HNACount$true/vol)

### Pool results
results <- cbind(Diversity.fbasis[Diversity.fbasis$D2>1600,],counts[Diversity.fbasis$D2>1600,])

### Add an extra column for biological replicates
bio_rep <- c(rep(1,21),rep(2,21),rep(1,41),rep(2,41),rep(3,41))
results <- cbind(results,bio_rep=factor(bio_rep))

### Select only S experiment for beta-div and further analysis
fbasis@basis <- fbasis@basis[!results$Treatment=="T",]
results <- results[!results$Treatment=="T",]  
results <- droplevels(results)
levels(results$Treatment) <- c("Control","Feeding")

### Beta diversity analysis
beta.div <- beta_div_fcm(fbasis,n=1,ord.type="PCoA")

##############################################################################
### Make some plots
##############################################################################
p1 <- ggplot(data=results, aes(x=Time, y=D2, fill=Treatment)) + 
  geom_point(shape=21, size=4,alpha=0.9)+
  scale_fill_brewer(palette="Accent")+
  geom_smooth(formula=y ~ x, color="black")+
  # geom_boxplot(mapping=factor(Time),alpha=0.4,outlier.shape=NA)+
  theme_bw()+
  labs(y="Phenotypic diversity - D2", x="Time (h)")
print(p1)

# p2 <- ggplot(data=results, aes(x=Time, y=D1, fill=bio_rep)) + 
#   geom_point(width = 0.1, shape=21, size=4,alpha=0.9)+
#   scale_fill_brewer(palette="Accent")+
#   # geom_boxplot(alpha=0.4,outlier.shape=NA)+
#   theme_bw()+
#   facet_grid(~Treatment)+
#   labs(y="Phenotypic diversity - D1", x="Time (h)")
# 
# p3 <- ggplot(data=results, aes(x=Time, y=D0, fill=bio_rep)) + 
#   geom_point(width = 0.1, shape=21, size=4,alpha=0.9)+
#   scale_fill_brewer(palette="Accent")+
#   # geom_boxplot(alpha=0.4,outlier.shape=NA)+
#   theme_bw()+
#   facet_grid(~Treatment)+
#   labs(y="Phenotypic diversity - D0", x="Time (h)")
# 
# grid.arrange(p3,p2,p1,nrow=3)

p1 <- ggplot(data=results, aes(x=factor(Time), y=D2, fill=Treatment)) + 
  geom_boxplot(alpha=0.9)+
  # geom_point(shape=21, size=5,alpha=0.9)+
  scale_fill_brewer(palette="Accent")+
  # geom_smooth(formula=y ~ x, color="black")+
  # geom_boxplot(mapping=factor(Time),alpha=0.4,outlier.shape=NA)+
  theme_bw()+
  labs(y="Phenotypic diversity - D2", x="Time (h)", title="Phenotypic alpha diversity")+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=16,face="bold"),
        title=element_text(size=20), legend.text=element_text(size=14))+ 
  guides(fill=FALSE)

p2 <- ggplot(data=results, aes(x=factor(Time), y=Total.cells, fill=Treatment)) + 
  geom_boxplot(alpha=0.9)+
  # geom_point(shape=21, size=5,alpha=0.9)+
  scale_fill_brewer(palette="Accent")+
  # geom_smooth(formula=y ~ x, color="black")+
  theme_bw()+
  labs(y="Cells/µL", x="Time (h)", title="Cell density")+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=16,face="bold"),
        title=element_text(size=20), legend.text=element_text(size=14))+ 
  guides(fill=FALSE)

### Beta diversity
beta.div.data <- data.frame(beta.div$points, Treatment=results$Treatment)
beta.div.data <- droplevels(beta.div.data)
levels(beta.div.data$Treatment) <- c("Control","Feeding")
var <- round(vegan::eigenvals(beta.div)/sum(vegan::eigenvals(beta.div))*100,1)

p.beta <- ggplot(data=beta.div.data, aes(x=X1, y=X2, fill=Treatment))+
  geom_point(shape=21, size=6, alpha=0.9)+
  theme_bw()+
  scale_fill_brewer(palette="Accent")+
  labs(x = paste0("Axis1 (",var[1], "%)"), y = paste0("Axis2 (",var[2], "%)"), title="Phenotypic beta diversity")+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=16,face="bold"),
        title=element_text(size=20), legend.text=element_text(size=14))
print(p.beta)
  
png(file="Combined.fcm.png",width=12,height=12,res=500,units="in", pointsize=12)
grid.arrange(arrangeGrob(p1,p2, ncol=2), p.beta, heights=c(4/4, 4/4), ncol=1)
dev.off()

