library("Phenoflow")
library("dplyr")

### Output files will be stored in this directory
path = c("/data_reference/FCM_MI")

### Import .fcs data
### Samples are automatically sorted according to name...
flowData <- read.flowSet(path = path, 
                         transformation = FALSE, pattern=".fcs")

### Select parameters (standard: two scatters and two FL) and 
### Transform data using the inverse hyperbolic sine
flowData_transformed <- transform(flowData,`FL1-H`=asinh(`FL1-H`), 
                                  `SSC-H`=asinh(`SSC-H`), 
                                  `FL3-H`=asinh(`FL3-H`), 
                                  `FSC-H`=asinh(`FSC-H`))
param=c("FL1-H", "FL3-H","SSC-H","FSC-H")
flowData_transformed = flowData_transformed[,param]
remove(flowData)

### Create a PolygonGate for extracting the single-cell information
### Input coordinates for gate in sqrcut1 in format: c(x,x,x,x,y,y,y,y)
sqrcut1 <- matrix(c(8.5,8.5,15,15,3,8,14,3),ncol=2, nrow=4)
colnames(sqrcut1) <- c("FL1-H","FL3-H")
polyGate1 <- polygonGate(.gate=sqrcut1, filterId = "Total Cells")

###  Gating quality check
xyplot(`FL3-H` ~ `FL1-H`, data=flowData_transformed[100], filter=polyGate1,
       scales=list(y=list(limits=c(0,15)),
                   x=list(limits=c(6,15))),
       axis = axis.default, nbin=125, 
       par.strip.text=list(col="white", font=2, cex=2), smooth=FALSE,xbins=750)

### Isolate only the cellular information based on the polyGate1
flowData_transformed <- flowCore::Subset(flowData_transformed, polyGate1)

### Normalize data between [0,1] on average, 
### this is required for using the bw=0.01 in the fingerprint calculation
summary <- fsApply(x=flowData_transformed,FUN=function(x) apply(x,2,max),use.exprs=TRUE)
max = mean(summary[,1])
mytrans <- function(x) x/max
flowData_transformed <- transform(flowData_transformed,`FL1-H`=mytrans(`FL1-H`),
                                  `FL3-H`=mytrans(`FL3-H`), 
                                  `SSC-H`=mytrans(`SSC-H`),
                                  `FSC-H`=mytrans(`FSC-H`))


### optional resample
### Calculate phenotypic diversity
Diversity.Accuri <- Phenoflow::Diversity_rf(flowData_transformed, d=3, R=100,
                                             param=param)

### Count nr of cells
sqrcut1 <- matrix(c(8.5,8.5,15,15,3,8,14,3)/max,ncol=2, nrow=4)
colnames(sqrcut1) <- c("FL1-H","FL3-H")
polyGate1 <- polygonGate(.gate=sqrcut1, filterId = "Total Cells")

s <- flowCore::filter(flowData_transformed, polyGate1)
TotalCount <- summary(s);TotalCount <- toTable(TotalCount)

### Extract the volumes
vol.temp<-c()
for(i in 1:length(flowData_transformed)){
  vol.temp[i] <- as.numeric(flowData_transformed[[i]]@description$`$VOL`)/1000
}

### Make count dataframe
### Counts in cells per µL
Counts.Accuri <- data.frame(Samples=flowData_transformed@phenoData@data$name, 
                                 counts = TotalCount$true, volume=vol.temp)

### Merge counts/diversity
tmp <- inner_join(Diversity.Accuri, Counts.Accuri, by=c("Sample_names"="Samples"))

### Write to file
write.csv2(results,file="Lakes_diversityFCM_F.csv")
