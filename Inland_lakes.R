library("Phenoflow")
library("dplyr")
library("gridExtra")
library("ggplot2")
library("formatR")
library("gridExtra")
library("cowplot")
library("RColorBrewer")

### Samples were diluted 2x
dilution <- 10

meta_inland <- read.csv2("/Users/rprops/Documents/FCM.Lakes/meta.inland2.csv")
meta_inland <- meta_inland[meta_inland$Fraction == "Free",]
meta_inland$Sampnum <- gsub(meta_inland$Sampnum, pattern=".", replacement="-", fixed=TRUE)
meta_inland$Sampnum <- paste(substring(meta_inland$Sampnum,1,3),seq(1:80), sep="-")
meta_inland$Sampnum <- gsub(meta_inland$Sampnum, pattern="Z15", replacement="Z14")
  
path = "/Users/rprops/Documents/FCM.Lakes/Accuri/Inland_Lakes"

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

Diversity.inland <- Diversity_rf(flowData_transformed, d=3, param = param, R = 3)

sqrcut1 <- matrix(c(asinh(12500),asinh(12500),15,15,3,9.55,14,3)/max,ncol=2, nrow=4)
colnames(sqrcut1) <- c("FL1-H","FL3-H")
rGate_HNA <- polygonGate(.gate=sqrcut1, filterId = "HNA")
sqrcut1 <- matrix(c(8.5,8.5,asinh(12500),asinh(12500),3,8,9.55,3)/max,ncol=2, nrow=4)
colnames(sqrcut1) <- c("FL1-H","FL3-H")
rGate_LNA <- polygonGate(.gate=sqrcut1, filterId = "LNA")

### Normalize total cell gate
sqrcut1 <- matrix(c(8.5,8.5,15,15,3,8,14,3)/max,ncol=2, nrow=4)
colnames(sqrcut1) <- c("FL1-H","FL3-H")
polyGate1 <- polygonGate(.gate=sqrcut1, filterId = "Total Cells")

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

counts

### Calculate means + errors for counts
grouplevels <- factor(gsub(counts$Samples, pattern="-1_rep.*", replacement="", fixed=FALSE))

means <- do.call(rbind,by(counts[,c(2,3,4)], INDICES = grouplevels, 
                          FUN = colMeans))

errors <- do.call(rbind,by(counts[,c(2,3,4)], INDICES = grouplevels, 
                            FUN = function(x) apply(x,2,sd)))

colnames(errors) <- c("Total.count.sd", "HNA.sd", "LNA.sd")

results.inland <- data.frame(samples = rownames(means), means, errors)

### Merge with metadata
results.inland.total <- inner_join(results.inland, meta_inland , by = c("samples"="Sampnum"))
results.inland.total <- results.inland.total[results.inland.total$Season=="Summer",]
results.inland.total$Depth
ggplot(data=results.inland.total, aes(x=Lake, y=HNA.cells/LNA.cells, fill=Depth))+
  geom_point(shape=21, size=8, alpha=0.4)

ggplot(data=results.inland.total, aes(x=Lake, y=HNA.cells/Total.cells, fill=Depth))+
  geom_point(shape=21, size=8, alpha=0.4)+
  facet_grid(~Invaded)
