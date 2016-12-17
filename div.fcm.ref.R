library("Phenoflow")
library("dplyr")

### Calculate 16S diversity for reference samples
phy.ref <- import_mothur(mothur_shared_file = "16S/stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list_CYANO.shared")
tax.ref <- data.frame(fread("16S/stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.0.03.cons_CYANO.taxonomy"))
tax.ref <- genX(tax.ref, "(", ")"); tax.ref  <- data.frame(tax.ref,stringsAsFactors =FALSE)
tax.ref <- tax_table(as.matrix(tax.ref)); rownames(tax.ref) <- tax.ref[,1]
phy.ref <- phyloseq(phy.ref, tax.ref)

### Only consider samples with > 10000 reads
phy.ref <- prune_samples(sample_sums(phy.ref)>10000, phy.ref)

div.ref <- Diversity_16S(phy.ref, R=100, brea=FALSE, thresh=500)
div.ref <- data.frame(div.ref)

write.csv2(div.ref, "otu.diversity16S.ref.csv")

### Calculate FCM diversity for reference samples

path = "/media/projects1/Ruben.FCM/Michigan/Accuri/Reference"

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
sqrcut1 <- matrix(c(8.75,8.75,14,14,3,7.5,14,3),ncol=2, nrow=4)
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

### Randomly resample to the lowest sample size
# flowData_transformed <- FCS_resample(flowData_transformed, replace=FALSE)

### Calculate fingerprint with bw = 0.01
# fbasis <- flowBasis(flowData_transformed, param, nbin=128, 
#                     bw=0.01,normalize=function(x) x)

### Calculate ecological parameters from normalized fingerprint 
### Densities will be normalized to the interval [0,1]
### n = number of replicates
### d = rounding factor
# Diversity.fbasis <- Diversity(fbasis,d=3,plot=FALSE, R=999)
Diversity.ref <- Diversity_rf(flowData_transformed, d=3, param=param)

### Rename for further processing
div.fcm.ref <- Diversity.ref 

### Now select the ones for which parallel sequencing data is available
meta.seq <- read.csv2("metadata_ref.csv")
meta.seq <- meta.seq[meta.seq$Sequenced=="yes",]

div.fcm.ref <- div.fcm.ref[grep(pattern = paste(meta.seq$Sample.name,collapse="|"),  x = as.character(div.fcm.ref$Sample_names)),]

### Now calculate the averages of these technical replicates
groupLevels <- factor(gsub(div.fcm.ref$Sample_names, pattern="_rep.*", replacement="", fixed=FALSE))
means <- do.call(rbind,by(div.fcm.ref[,3:5], INDICES = groupLevels, 
                FUN = colMeans))
errors <- do.call(rbind,by(div.fcm.ref[,6:8], INDICES = groupLevels, 
             FUN = function(x) sqrt(colSums(x^2))/ncol(x)))
div.fcm.ref.merged <- data.frame(cbind(means,errors), sample_fcm = rownames(means))

### And match them to the corresponding sequencing data
div.seq.ref <- read.csv2("otu.diversity16S.ref.csv")

lb <- read.csv2("labels_ref.csv")

div.fcm.ref.merged <- inner_join(div.fcm.ref.merged, lb, by=c("sample_fcm"="Sample_fcm"))

### Adjust colnames
colnames(div.fcm.ref.merged) <- c("D0.fcm","D1.fcm","D2.fcm","sd.D0.fcm","sd.D1.fcm","sd.D2.fcm","Sample_fcm","Sample_seq")

### Merge
div.total.ref <- inner_join(div.fcm.ref.merged, div.seq.ref, by=c("Sample_seq"="Sample"))

### Tidy up a bit

div.total.ref <- data.frame(div.total.ref[,7:8], div.total.ref[,1:6], div.total.ref[,9:18])

### Write csv

write.csv2(div.total.ref, "div.ref.merged.csv")
